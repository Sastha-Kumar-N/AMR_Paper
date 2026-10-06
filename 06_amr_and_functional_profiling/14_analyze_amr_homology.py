#!/usr/bin/env python3
"""Summarize coordinate-supported Tier-A AMR protein homology clusters."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict

import pandas as pd

from analysis_common import ROOT, atomic_write_tsv, prevalence_band


N_EVALUABLE = 1350


def read_fasta(path):
    records: dict[str, str] = {}
    header = ""
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records[header] = "".join(sequence)
                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
        if header:
            records[header] = "".join(sequence)
    return records


def split_semicolon(value: str) -> list[str]:
    return [x.strip() for x in value.split(";") if x.strip()]


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    sequence_map = pd.read_csv(ROOT / "02_cleaned/AMR_sequence_id_map.tsv", sep="\t", dtype=str, keep_default_na=False)
    metadata = sequence_map.set_index("sequence_id").to_dict("index")
    sequences = read_fasta(ROOT / "02_cleaned/AMR_mapped_sequences.faa")

    amr = pd.read_csv(ROOT / "02_cleaned/AMR_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    primary = amr[amr["primary_analysis"] == "yes"]
    family_classes: dict[str, set[str]] = defaultdict(set)
    family_mechanisms: dict[str, set[str]] = defaultdict(set)
    for row in primary.itertuples(index=False):
        family_classes[row.gene_family].update(split_semicolon(row.normalized_drug_classes))
        family_mechanisms[row.gene_family].update(split_semicolon(row.resistance_mechanism))

    locus_to_panaroo: dict[tuple[str, str], str] = {}
    for chunk in pd.read_csv(
        ROOT / "02_cleaned/funcscan_to_panaroo_mapping.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=["genome_id", "funcscan_category", "prokka_gene_id", "panaroo_cluster"],
        chunksize=250000,
    ):
        mapped = chunk[(chunk["funcscan_category"] == "AMR") & (chunk["prokka_gene_id"] != "") & (chunk["panaroo_cluster"] != "")]
        for row in mapped.itertuples(index=False):
            locus_to_panaroo[(row.genome_id, row.prokka_gene_id)] = row.panaroo_cluster

    membership_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    representative_rows: list[dict[str, object]] = []
    primary_cluster_by_sequence: dict[str, str] = {}
    for threshold in (90, 70, 50):
        directory = ROOT / f"03_matrices/mmseqs_amr_{threshold}"
        cluster_path = directory / f"amr_{threshold}_cluster.tsv"
        groups: dict[str, list[str]] = defaultdict(list)
        with cluster_path.open(encoding="utf-8") as handle:
            for line in handle:
                representative, member = line.rstrip("\n").split("\t")[:2]
                groups[representative].append(member)
        ordered = sorted(groups, key=lambda rep: (-len(groups[rep]), rep))
        stable_ids = {rep: f"AMR{threshold}_C{index:04d}" for index, rep in enumerate(ordered, start=1)}

        identities: dict[tuple[str, str], float] = {}
        align_path = directory / f"amr_{threshold}_representative_alignments.tsv"
        with align_path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    identities[(fields[0], fields[1])] = float(fields[2])

        for representative in ordered:
            cluster_id = stable_ids[representative]
            members = groups[representative]
            genomes = {metadata[member]["genome_id"] for member in members}
            families = sorted(
                {family for member in members for family in split_semicolon(metadata[member]["gene_families"])}
            )
            annotations = Counter(metadata[member]["annotation"] for member in members if metadata[member]["annotation"])
            panaroo_clusters = sorted(
                {
                    locus_to_panaroo.get((metadata[member]["genome_id"], metadata[member]["prokka_gene_id"]), "")
                    for member in members
                }
                - {""}
            )
            drug_classes = sorted({term for family in families for term in family_classes.get(family, set())})
            mechanisms = sorted({term for family in families for term in family_mechanisms.get(family, set())})
            identity_values = [identities[(representative, member)] for member in members if (representative, member) in identities]
            row = {
                "cluster": cluster_id,
                "identity_threshold_percent": threshold,
                "alignment_coverage_threshold": 0.8,
                "representative": representative,
                "representative_sequence": sequences.get(representative, ""),
                "gene_annotations": ";".join(families),
                "representative_Prokka_annotation": annotations.most_common(1)[0][0] if annotations else "",
                "number_sequences": len(members),
                "number_genomes": len(genomes),
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * len(genomes) / N_EVALUABLE,
                "mean_identity_to_representative_percent": sum(identity_values) / len(identity_values) if identity_values else "",
                "minimum_identity_to_representative_percent": min(identity_values) if identity_values else "",
                "identity_observations": len(identity_values),
                "drug_classes": ";".join(drug_classes),
                "mechanisms": ";".join(mechanisms),
                "Panaroo_clusters": ";".join(panaroo_clusters),
                "Panaroo_cluster_count": len(panaroo_clusters),
                "annotation_family_count": len(families),
                "prevalence_band": prevalence_band(100 * len(genomes) / N_EVALUABLE),
                "evidence_scope": "MMseqs2 protein homology cluster; not definitive orthology",
            }
            summary_rows.append(row)
            representative_rows.append(
                {
                    "identity_threshold_percent": threshold,
                    "cluster": cluster_id,
                    "representative_sequence_id": representative,
                    "representative_sequence": sequences.get(representative, ""),
                    "representative_genome": metadata[representative]["genome_id"],
                    "representative_Prokka_gene_id": metadata[representative]["prokka_gene_id"],
                    "gene_annotations": ";".join(families),
                    "number_sequences": len(members),
                    "number_genomes": len(genomes),
                }
            )
            for member in members:
                member_identity = identities.get((representative, member), "")
                membership_rows.append(
                    {
                        "identity_threshold_percent": threshold,
                        "alignment_coverage_threshold": 0.8,
                        "cluster": cluster_id,
                        "representative_sequence_id": representative,
                        "sequence_id": member,
                        "genome_id": metadata[member]["genome_id"],
                        "prokka_gene_id": metadata[member]["prokka_gene_id"],
                        "gene_families": metadata[member]["gene_families"],
                        "identity_to_representative_percent": member_identity,
                    }
                )
                if threshold == 50:
                    primary_cluster_by_sequence[member] = cluster_id

    atomic_write_tsv(ROOT / "03_matrices/AMR_sequence_clusters.tsv", membership_rows, list(membership_rows[0]))
    primary_summary = [row for row in summary_rows if row["identity_threshold_percent"] == 50]
    atomic_write_tsv(ROOT / "06_tables/AMR_homology_clusters.tsv", primary_summary, list(primary_summary[0]))
    atomic_write_tsv(ROOT / "06_tables/AMR_homology_clusters_all_thresholds.tsv", summary_rows, list(summary_rows[0]))
    atomic_write_tsv(ROOT / "06_tables/AMR_cluster_representatives.tsv", representative_rows, list(representative_rows[0]))

    # Homology-cluster multiplicity is an independent, explicitly putative
    # duplication screen complementary to Panaroo family copy counts.
    duplicate_rows: list[dict[str, object]] = []
    primary_members = [row for row in membership_rows if row["identity_threshold_percent"] == 50]
    by_genome_cluster: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in primary_members:
        by_genome_cluster[(str(row["genome_id"]), str(row["cluster"]))].append(row)
    for (genome, cluster), rows in by_genome_cluster.items():
        loci = sorted({str(row["prokka_gene_id"]) for row in rows})
        if len(loci) <= 1:
            continue
        duplicate_rows.append(
            {
                "genome": genome,
                "homolog_cluster": cluster,
                "copy_number": len(loci),
                "member_ids": ";".join(loci),
                "annotations": ";".join(sorted({str(row["gene_families"]) for row in rows})),
                "functional_category": "AMR",
                "interpretation": "putative paralogous or duplicated homolog-family members; no gene-tree validation",
            }
        )
    duplicate_rows.sort(key=lambda row: (-int(row["copy_number"]), str(row["genome"]), str(row["homolog_cluster"])))
    fields = ["genome", "homolog_cluster", "copy_number", "member_ids", "annotations", "functional_category", "interpretation"]
    atomic_write_tsv(ROOT / "06_tables/AMR_putative_paralogous_homologs.tsv", duplicate_rows, fields)

    print(
        f"AMR homology: {len(primary_summary)} clusters at 50% identity; "
        f"largest={primary_summary[0]['cluster']} ({primary_summary[0]['number_sequences']} sequences)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
