#!/usr/bin/env python3
"""Map FuncScan features to Prokka CDSs and the existing clinical Panaroo run."""

from __future__ import annotations

import argparse
import bisect
import csv
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
from scipy.stats import fisher_exact

from analysis_common import (
    PANAROO_CSV,
    PANAROO_RTAB,
    PROKKA,
    ROOT,
    atomic_write_tsv,
    bh_fdr,
    clean_missing,
    load_clinical_metadata,
    parse_gff_attributes,
)


N_EVALUABLE = 1350
_INTERVAL_INDEX: dict[int, tuple[list[dict[str, object]], list[int], list[int]]] = {}


def normalized_contig_key(value: str) -> str:
    """Return a conservative separator-insensitive contig comparison key.

    The three non-accession FuncScan sample names retain NCBI-style contig
    labels such as ``refseq|NZ_CP081477.2|chromosome``, whereas their Prokka
    GFF files use ``refseq_NZ_CP081477.2_chromosome``.  Normalization changes
    separators only; a mapping is accepted below only when the key uniquely
    identifies one Prokka contig.
    """
    return re.sub(r"[^A-Za-z0-9.]+", "_", clean_missing(value)).strip("_")


def resolve_contig_id(genes: dict[str, list[dict[str, object]]], contig: str) -> str:
    if contig in genes:
        return contig
    key = normalized_contig_key(contig)
    matches = [candidate for candidate in genes if normalized_contig_key(candidate) == key]
    return matches[0] if len(matches) == 1 else ""


def pangenome_category(n_genomes: int) -> str:
    percent = 100 * n_genomes / N_EVALUABLE
    if percent >= 99:
        return "core"
    if percent >= 95:
        return "soft-core"
    if percent >= 15:
        return "shell"
    return "cloud"


def parse_gff(path: Path) -> tuple[dict[str, list[dict[str, object]]], dict[str, int]]:
    genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    contig_lengths: dict[str, int] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.rstrip().split()
                if len(parts) >= 4:
                    contig_lengths[parts[1]] = int(parts[3]) - int(parts[2]) + 1
                continue
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = parse_gff_attributes(parts[8])
            locus = attrs.get("locus_tag") or attrs.get("ID", "")
            if not locus:
                continue
            genes[parts[0]].append(
                {
                    "contig_id": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "prokka_gene_id": locus,
                    "gene_name": attrs.get("gene") or attrs.get("Name", ""),
                    "product": attrs.get("product", ""),
                }
            )
    for contig in genes:
        genes[contig].sort(key=lambda g: (int(g["start"]), int(g["end"]), str(g["prokka_gene_id"])))
        if contig not in contig_lengths:
            contig_lengths[contig] = max(int(g["end"]) for g in genes[contig])
    return genes, contig_lengths


def parse_ordinal(sequence_id: str) -> tuple[str, int] | None:
    match = re.fullmatch(r"(.+)_([0-9]+)", clean_missing(sequence_id))
    if not match:
        return None
    return match.group(1), int(match.group(2))


def ordinal_candidate(genes: dict[str, list[dict[str, object]]], sequence_id: str) -> dict[str, object] | None:
    parsed = parse_ordinal(sequence_id)
    if not parsed:
        return None
    contig, ordinal = parsed
    options = genes.get(resolve_contig_id(genes, contig), [])
    if 1 <= ordinal <= len(options):
        return options[ordinal - 1]
    return None


def interval_candidates(
    genes: dict[str, list[dict[str, object]]], contig: str, start: int, end: int
) -> list[tuple[dict[str, object], int, float, float]]:
    if end < start:
        start, end = end, start
    options = genes.get(resolve_contig_id(genes, contig), [])
    if not options:
        return []
    cache_key = id(options)
    cached = _INTERVAL_INDEX.get(cache_key)
    if cached is None or cached[0] is not options:
        starts = [int(g["start"]) for g in options]
        prefix_max_end: list[int] = []
        running_max = 0
        for gene in options:
            running_max = max(running_max, int(gene["end"]))
            prefix_max_end.append(running_max)
        _INTERVAL_INDEX[cache_key] = (options, starts, prefix_max_end)
    else:
        _, starts, prefix_max_end = cached
    upper = bisect.bisect_right(starts, end)
    output: list[tuple[dict[str, object], int, float, float]] = []
    query_len = end - start + 1
    index = upper - 1
    while index >= 0:
        # The prefix maximum makes this an exact interval query: once every
        # preceding feature ends before the query, no earlier overlap is possible.
        if prefix_max_end[index] < start:
            break
        gene = options[index]
        gene_start, gene_end = int(gene["start"]), int(gene["end"])
        if gene_end >= start:
            intersection = max(0, min(end, gene_end) - max(start, gene_start) + 1)
            if intersection:
                gene_len = gene_end - gene_start + 1
                output.append((gene, intersection, intersection / query_len, intersection / gene_len))
        index -= 1
    return output


def best_gene_mapping(
    genes: dict[str, list[dict[str, object]]],
    contig: str,
    start: str,
    end: str,
    sequence_id: str,
) -> tuple[dict[str, object] | None, str, str, str, bool | None]:
    ordinal = ordinal_candidate(genes, sequence_id)
    coordinate = None
    coverage = ""
    if clean_missing(contig) and clean_missing(start) and clean_missing(end):
        candidates = interval_candidates(genes, contig, int(float(start)), int(float(end)))
        if candidates:
            coordinate, _, query_cov, gene_cov = max(candidates, key=lambda x: (min(x[2], x[3]), x[1]))
            coverage = f"{min(query_cov, gene_cov):.6f}"
    concordant = None
    if coordinate is not None and ordinal is not None:
        concordant = coordinate["prokka_gene_id"] == ordinal["prokka_gene_id"]
    if coordinate is not None:
        score = float(coverage)
        method = "coordinate_overlap"
        if concordant is True:
            method = "coordinate_and_ordinal_concordant"
        confidence = "high" if score >= 0.8 else "moderate" if score >= 0.5 else "low"
        if score < 0.5:
            return None, "coordinate_overlap_below_threshold", coverage, "unresolved", concordant
        return coordinate, method, coverage, confidence, concordant
    if ordinal is not None:
        return None, "pyrodigal_ordinal_only_not_mapped", "", "unresolved", concordant
    return None, "unresolved", "", "unresolved", concordant


def iter_fasta(path: Path) -> Iterable[tuple[str, str, str]]:
    header = ""
    sequence: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    parts = header.split(maxsplit=1)
                    yield parts[0], parts[1] if len(parts) > 1 else "", "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            parts = header.split(maxsplit=1)
            yield parts[0], parts[1] if len(parts) > 1 else "", "".join(sequence)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    clinical, _, canonical_to_raw = load_clinical_metadata()
    genomes = sorted(canonical_to_raw)
    amp = pd.read_csv(ROOT / "02_cleaned/AMP_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    amr = pd.read_csv(ROOT / "02_cleaned/AMR_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    cazyme = pd.read_csv(ROOT / "02_cleaned/CAZyme_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    bgc = pd.read_csv(ROOT / "02_cleaned/BGC_consensus_loci.tsv", sep="\t", dtype=str, keep_default_na=False)
    amp_groups = amp.groupby("genome_id", sort=False)
    amr_groups = amr.groupby("genome_id", sort=False)
    cazyme_groups = cazyme.groupby("genome_id", sort=False)
    bgc_groups = bgc.groupby("genome_id", sort=False)

    prelim_fields = [
        "genome_id",
        "funcscan_feature",
        "funcscan_category",
        "functional_subfeature",
        "evidence_tier",
        "sequence_id",
        "prokka_gene_id",
        "prokka_gene_name",
        "prokka_product",
        "contig_id",
        "gene_start",
        "gene_end",
        "match_method",
        "match_identity",
        "match_coverage",
        "mapping_confidence",
        "source_feature_coordinates",
    ]
    fd, prelim_name = tempfile.mkstemp(prefix=".funcscan_prokka_mapping.", dir=ROOT / "02_cleaned")
    target_loci: set[str] = set()
    genome_stats: list[dict[str, object]] = []
    amp_validation = Counter()
    bgc_gene_ids_by_genome: dict[str, set[str]] = defaultdict(set)

    amr_fasta_fd, amr_fasta_tmp = tempfile.mkstemp(prefix=".AMR_mapped_sequences.", dir=ROOT / "02_cleaned")
    amr_map_fd, amr_map_tmp = tempfile.mkstemp(prefix=".AMR_sequence_id_map.", dir=ROOT / "02_cleaned")
    amr_sequence_counter = 0
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as prelim_handle, os.fdopen(
        amr_fasta_fd, "w", encoding="utf-8", newline=""
    ) as amr_fasta_handle, os.fdopen(amr_map_fd, "w", encoding="utf-8", newline="") as amr_map_handle:
        writer = csv.DictWriter(prelim_handle, fieldnames=prelim_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        amr_map_writer = csv.writer(amr_map_handle, delimiter="\t", lineterminator="\n")
        amr_map_writer.writerow(["sequence_id", "genome_id", "prokka_gene_id", "gene_families", "annotation"])

        def write_mapping(
            genome: str,
            feature: str,
            category: str,
            subfeature: str,
            evidence: str,
            sequence_id: str,
            gene: dict[str, object] | None,
            method: str,
            coverage: str,
            confidence: str,
            source_coordinates: str,
            identity: str = "",
        ) -> None:
            if gene:
                target_loci.add(str(gene["prokka_gene_id"]))
            writer.writerow(
                {
                    "genome_id": genome,
                    "funcscan_feature": feature,
                    "funcscan_category": category,
                    "functional_subfeature": subfeature,
                    "evidence_tier": evidence,
                    "sequence_id": sequence_id,
                    "prokka_gene_id": gene["prokka_gene_id"] if gene else "",
                    "prokka_gene_name": gene["gene_name"] if gene else "",
                    "prokka_product": gene["product"] if gene else "",
                    "contig_id": gene["contig_id"] if gene else "",
                    "gene_start": gene["start"] if gene else "",
                    "gene_end": gene["end"] if gene else "",
                    "match_method": method,
                    "match_identity": identity,
                    "match_coverage": coverage,
                    "mapping_confidence": confidence,
                    "source_feature_coordinates": source_coordinates,
                }
            )

        for position, genome in enumerate(genomes, start=1):
            _INTERVAL_INDEX.clear()
            raw_id = canonical_to_raw[genome]
            gff = PROKKA / raw_id / f"{raw_id}.gff"
            faa = PROKKA / raw_id / f"{raw_id}.faa"
            if not gff.exists():
                raise FileNotFoundError(gff)
            if not faa.exists():
                raise FileNotFoundError(faa)
            genes, contig_lengths = parse_gff(gff)
            gene_by_locus = {
                str(gene["prokka_gene_id"]): gene
                for contig_genes in genes.values()
                for gene in contig_genes
            }
            faa_records: dict[str, tuple[str, str]] = {}
            loci_by_protein_sequence: dict[str, list[str]] = defaultdict(list)
            for locus, description, sequence in iter_fasta(faa):
                faa_records[locus] = (description, sequence)
                loci_by_protein_sequence[sequence].append(locus)
            total_cds = sum(len(v) for v in genes.values())
            genome_stats.append(
                {
                    "genome_id": genome,
                    "raw_genome_id": raw_id,
                    "genome_size_bp": sum(contig_lengths.values()),
                    "total_CDS": total_cds,
                    "contig_count": len(contig_lengths),
                    "source_gff": str(gff),
                }
            )

            if genome in amp_groups.groups:
                for row in amp_groups.get_group(genome).to_dict("records"):
                    gene, method, coverage, confidence, concordant = best_gene_mapping(
                        genes,
                        row["contig_id"],
                        row["CDS_start"],
                        row["CDS_end"],
                        row["CDS_id"],
                    )
                    identity = ""
                    exact_loci = loci_by_protein_sequence.get(row["aa_sequence"], [])
                    if gene is not None and str(gene["prokka_gene_id"]) in exact_loci:
                        method = "coordinate_and_exact_protein_sequence"
                        confidence = "high"
                        identity = "100"
                        amp_validation["coordinate_exact_sequence_concordant"] += 1
                    elif len(exact_loci) == 1:
                        if gene is not None:
                            amp_validation["coordinate_exact_sequence_discordant_exact_sequence_used"] += 1
                        else:
                            amp_validation["exact_sequence_only_unique"] += 1
                        gene = gene_by_locus.get(exact_loci[0])
                        method = "exact_protein_sequence_unique"
                        coverage = "1.000000"
                        confidence = "high"
                        identity = "100"
                    elif gene is not None:
                        amp_validation["coordinate_mapped_no_exact_Prokka_protein_match"] += 1
                    elif len(exact_loci) > 1:
                        amp_validation["exact_sequence_nonunique_unresolved"] += 1
                    else:
                        amp_validation["unresolved"] += 1
                    write_mapping(
                        genome,
                        row["feature_id"],
                        "AMP",
                        row["normalized_feature_name"],
                        row["evidence_tier"],
                        row["CDS_id"],
                        gene,
                        method,
                        coverage,
                        confidence,
                        row["coordinates_if_available"],
                        identity,
                    )

            genome_amr_loci: dict[str, set[str]] = defaultdict(set)
            if genome in amr_groups.groups:
                for row in amr_groups.get_group(genome).to_dict("records"):
                    gene, method, coverage, confidence, _ = best_gene_mapping(
                        genes,
                        row["contig_id"],
                        row["coordinate_start"],
                        row["coordinate_end"],
                        row["sequence_id"],
                    )
                    write_mapping(
                        genome,
                        row["feature_id"],
                        "AMR",
                        row["gene_family"],
                        row["evidence_tier"],
                        row["sequence_id"],
                        gene,
                        method,
                        coverage,
                        confidence,
                        row["coordinates_if_available"],
                    )
                    if gene and row["primary_analysis"] == "yes":
                        genome_amr_loci[str(gene["prokka_gene_id"])].add(row["gene_family"])

            if genome in cazyme_groups.groups:
                for row in cazyme_groups.get_group(genome).to_dict("records"):
                    sequence_id = row["sequence_id"]
                    write_mapping(
                        genome,
                        row["feature_id"],
                        "CAZyme",
                        row["CAZyme_family"],
                        row["analysis_tier"],
                        sequence_id,
                        None,
                        "pyrodigal_ordinal_only_not_mapped",
                        "",
                        "unresolved",
                        "",
                    )

            if genome in bgc_groups.groups:
                for row in bgc_groups.get_group(genome).to_dict("records"):
                    candidates = interval_candidates(
                        genes,
                        row["contig_id"],
                        int(row["BGC_start"]),
                        int(row["BGC_end"]),
                    )
                    for gene, _, _, gene_cov in candidates:
                        if gene_cov < 0.5:
                            continue
                        bgc_gene_ids_by_genome[genome].add(str(gene["prokka_gene_id"]))
                        write_mapping(
                            genome,
                            row["consensus_BGC_id"],
                            "BGC",
                            row["normalized_BGC_types"],
                            row["confidence_class"],
                            "",
                            gene,
                            "BGC_region_coordinate_overlap",
                            f"{gene_cov:.6f}",
                            "high",
                            f"{row['contig_id']}:{row['BGC_start']}-{row['BGC_end']}",
                        )

            if genome_amr_loci:
                for locus in sorted(genome_amr_loci):
                    if locus not in faa_records:
                        continue
                    description, sequence = faa_records[locus]
                    amr_sequence_counter += 1
                    seq_id = f"AMRSEQ{amr_sequence_counter:08d}"
                    families = ";".join(sorted(genome_amr_loci[locus]))
                    amr_fasta_handle.write(f">{seq_id}\n{sequence}\n")
                    amr_map_writer.writerow([seq_id, genome, locus, families, description])

    atomic_write_tsv(
        ROOT / "02_cleaned/genome_size_cds.tsv",
        genome_stats,
        list(genome_stats[0]),
    )
    atomic_write_tsv(
        ROOT / "01_audit/AMP_mapping_validation.tsv",
        ({"status": key, "count": value} for key, value in sorted(amp_validation.items())),
        ["status", "count"],
    )
    os.replace(amr_fasta_tmp, ROOT / "02_cleaned/AMR_mapped_sequences.faa")
    os.replace(amr_map_tmp, ROOT / "02_cleaned/AMR_sequence_id_map.tsv")

    # Map only relevant Prokka loci to clinical Panaroo clusters, while also
    # recording cluster copy structure needed for putative-paralog statements.
    csv.field_size_limit(sys.maxsize)
    locus_to_cluster: dict[str, str] = {}
    cluster_annotation: dict[str, tuple[str, str]] = {}
    cluster_max_copy: dict[str, int] = {}
    cluster_multicopy_genomes: dict[str, int] = {}
    with PANAROO_CSV.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            cluster = row[0]
            relevant = False
            max_copy = 0
            multicopy = 0
            for cell in row[3:]:
                ids = [x.strip() for x in cell.split(";") if x.strip()]
                max_copy = max(max_copy, len(ids))
                if len(ids) > 1:
                    multicopy += 1
                for locus in ids:
                    if locus in target_loci:
                        locus_to_cluster[locus] = cluster
                        relevant = True
            if relevant:
                cluster_annotation[cluster] = (row[1], row[2])
                cluster_max_copy[cluster] = max_copy
                cluster_multicopy_genomes[cluster] = multicopy

    cluster_prevalence: dict[str, tuple[int, float, str]] = {}
    background_counts = Counter()
    with PANAROO_RTAB.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) - 1 != N_EVALUABLE:
            raise ValueError(f"Panaroo Rtab has {len(header)-1} genomes, expected {N_EVALUABLE}")
        for row in reader:
            n = sum(int(x) for x in row[1:])
            category = pangenome_category(n)
            background_counts[category] += 1
            cluster_prevalence[row[0]] = (n, 100 * n / N_EVALUABLE, category)

    final_path = ROOT / "02_cleaned/funcscan_to_panaroo_mapping.tsv"
    final_fields = prelim_fields + [
        "panaroo_cluster",
        "panaroo_nonunique_gene_name",
        "panaroo_annotation",
        "panaroo_n_genomes",
        "panaroo_prevalence_percent",
        "panaroo_frequency_category",
        "panaroo_max_copies_per_genome",
        "panaroo_multicopy_genomes",
    ]
    final_fd, final_tmp = tempfile.mkstemp(prefix=f".{final_path.name}.", dir=final_path.parent)
    mapping_counts = Counter()
    with open(prelim_name, encoding="utf-8", newline="") as source, os.fdopen(
        final_fd, "w", encoding="utf-8", newline=""
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(destination, fieldnames=final_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            cluster = locus_to_cluster.get(row["prokka_gene_id"], "")
            nonunique, annotation = cluster_annotation.get(cluster, ("", ""))
            n, prevalence, category = cluster_prevalence.get(cluster, (0, 0.0, ""))
            row.update(
                {
                    "panaroo_cluster": cluster,
                    "panaroo_nonunique_gene_name": nonunique,
                    "panaroo_annotation": annotation,
                    "panaroo_n_genomes": n if cluster else "",
                    "panaroo_prevalence_percent": prevalence if cluster else "",
                    "panaroo_frequency_category": category,
                    "panaroo_max_copies_per_genome": cluster_max_copy.get(cluster, ""),
                    "panaroo_multicopy_genomes": cluster_multicopy_genomes.get(cluster, ""),
                }
            )
            mapping_counts[(row["funcscan_category"], "total")] += 1
            if row["prokka_gene_id"]:
                mapping_counts[(row["funcscan_category"], "prokka_mapped")] += 1
            if cluster:
                mapping_counts[(row["funcscan_category"], "panaroo_mapped")] += 1
            writer.writerow(row)
    os.replace(final_tmp, final_path)
    os.unlink(prelim_name)

    mapping_summary = []
    for category in ("AMP", "AMR", "BGC", "CAZyme"):
        total = mapping_counts[(category, "total")]
        prokka = mapping_counts[(category, "prokka_mapped")]
        panaroo = mapping_counts[(category, "panaroo_mapped")]
        mapping_summary.append(
            {
                "functional_category": category,
                "total_feature_mapping_rows": total,
                "prokka_mapped_rows": prokka,
                "prokka_mapping_percent": 100 * prokka / total if total else 0,
                "panaroo_mapped_rows": panaroo,
                "panaroo_mapping_percent": 100 * panaroo / total if total else 0,
            }
        )
    atomic_write_tsv(
        ROOT / "01_audit/funcscan_panaroo_mapping_summary.tsv",
        mapping_summary,
        list(mapping_summary[0]),
    )

    # Chunked aggregation prevents the large BGC gene mapping from becoming a
    # memory bottleneck while retaining exact locus-level deduplication.
    unique_loci: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    annotations_by_cluster: dict[str, Counter[str]] = defaultdict(Counter)
    subfeatures_by_cluster: dict[str, set[str]] = defaultdict(set)
    categories_by_cluster: dict[str, set[str]] = defaultdict(set)
    paralog_members: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    bgc_gene_aggregate: dict[str, dict[str, object]] = {}
    for chunk in pd.read_csv(final_path, sep="\t", dtype=str, keep_default_na=False, chunksize=200000):
        mapped = chunk[(chunk["prokka_gene_id"] != "") & (chunk["panaroo_cluster"] != "")]
        for row in mapped.itertuples(index=False):
            category = row.funcscan_category
            unique_loci[category].add((row.genome_id, row.prokka_gene_id, row.panaroo_cluster))
            annotation = row.prokka_gene_name or row.prokka_product or row.panaroo_annotation
            if annotation:
                annotations_by_cluster[row.panaroo_cluster][annotation] += 1
            subfeatures_by_cluster[row.panaroo_cluster].add(row.functional_subfeature)
            categories_by_cluster[row.panaroo_cluster].add(category)
            paralog_members[(row.genome_id, row.panaroo_cluster, category)].add(row.prokka_gene_id)
            if category == "BGC":
                agg = bgc_gene_aggregate.setdefault(
                    row.panaroo_cluster,
                    {
                        "genomes": set(),
                        "regions": set(),
                        "types": set(),
                        "genes": Counter(),
                        "products": Counter(),
                    },
                )
                agg["genomes"].add(row.genome_id)
                agg["regions"].add(row.funcscan_feature)
                agg["types"].update(x for x in row.functional_subfeature.split(";") if x)
                if row.prokka_gene_name:
                    agg["genes"][row.prokka_gene_name] += 1
                if row.prokka_product:
                    agg["products"][row.prokka_product] += 1

    distribution_paths = {
        "AMP": "AMP_pangenome_distribution.tsv",
        "AMR": "AMR_pangenome_distribution.tsv",
        "BGC": "BGC_gene_pangenome_distribution.tsv",
        "CAZyme": "CAZyme_pangenome_distribution.tsv",
    }
    enrichment_rows: list[dict[str, object]] = []
    for category in distribution_paths:
        loci = unique_loci.get(category, set())
        by_frequency = Counter(cluster_prevalence[cluster][2] for _, _, cluster in loci)
        cluster_sets = {freq: {cluster for _, _, cluster in loci if cluster_prevalence[cluster][2] == freq} for freq in background_counts}
        rows = []
        all_clusters = {cluster for _, _, cluster in loci}
        for freq in ("core", "soft-core", "shell", "cloud"):
            clusters = cluster_sets.get(freq, set())
            locus_count = sum(1 for _, _, cluster in loci if cluster in clusters)
            rows.append(
                {
                    "functional_category": category,
                    "pangenome_category": freq,
                    "unique_functional_gene_occurrences": locus_count,
                    "unique_panaroo_clusters": len(clusters),
                    "percent_of_mapped_functional_clusters": 100 * len(clusters) / len(all_clusters) if all_clusters else 0,
                    "background_panaroo_clusters": background_counts[freq],
                    "clinical_pangenome_clusters_total": sum(background_counts.values()),
                    "category_threshold": {"core": ">=99%", "soft-core": "95-<99%", "shell": "15-<95%", "cloud": "<15%"}[freq],
                }
            )
            if all_clusters:
                a = len(clusters)
                b = len(all_clusters) - a
                c = background_counts[freq] - a
                d = sum(background_counts.values()) - background_counts[freq] - b
                odds, pvalue = fisher_exact([[a, b], [c, d]])
                enrichment_rows.append(
                    {
                        "functional_category": category,
                        "pangenome_category": freq,
                        "null_hypothesis": "functional Panaroo clusters are not enriched in this frequency category",
                        "functional_clusters_in_category": a,
                        "functional_clusters_outside_category": b,
                        "background_clusters_in_category": c,
                        "background_clusters_outside_category": d,
                        "odds_ratio": odds,
                        "raw_p": pvalue,
                    }
                )
        atomic_write_tsv(ROOT / f"06_tables/{distribution_paths[category]}", rows, list(rows[0]))
    qvalues = bh_fdr([float(row["raw_p"]) for row in enrichment_rows])
    for row, qvalue in zip(enrichment_rows, qvalues):
        row["q_value"] = qvalue
    atomic_write_tsv(
        ROOT / "04_statistics/functional_pangenome_category_enrichment.tsv",
        enrichment_rows,
        list(enrichment_rows[0]),
    )

    orthogroup_rows: list[dict[str, object]] = []
    for cluster in sorted(categories_by_cluster):
        n, prevalence, freq = cluster_prevalence[cluster]
        annotation = annotations_by_cluster[cluster].most_common(1)[0][0] if annotations_by_cluster[cluster] else cluster_annotation.get(cluster, ("", ""))[1]
        orthogroup_rows.append(
            {
                "panaroo_cluster": cluster,
                "annotation": annotation,
                "functional_categories": ";".join(sorted(categories_by_cluster[cluster])),
                "functional_subfeatures": ";".join(sorted(subfeatures_by_cluster[cluster])),
                "n_genomes": n,
                "prevalence_percent": prevalence,
                "pangenome_category": freq,
                "max_copies_per_genome": cluster_max_copy.get(cluster, ""),
                "multicopy_genomes": cluster_multicopy_genomes.get(cluster, ""),
                "single_copy_across_clinical_pangenome": "yes" if cluster_max_copy.get(cluster, 2) <= 1 else "no",
                "evidence_scope": "Panaroo gene-family grouping; orthology is not independently gene-tree validated",
            }
        )
    orthogroup_rows.sort(key=lambda r: (-int(r["n_genomes"]), str(r["panaroo_cluster"])))
    atomic_write_tsv(
        ROOT / "06_tables/functional_panaroo_orthogroups.tsv",
        orthogroup_rows,
        list(orthogroup_rows[0]),
    )
    single_copy_rows = [row for row in orthogroup_rows if row["single_copy_across_clinical_pangenome"] == "yes"]
    atomic_write_tsv(
        ROOT / "06_tables/single_copy_orthologous_functional_genes.tsv",
        single_copy_rows,
        list(orthogroup_rows[0]),
    )

    paralog_rows: list[dict[str, object]] = []
    for (genome, cluster, category), members in paralog_members.items():
        if len(members) <= 1:
            continue
        annotation = annotations_by_cluster[cluster].most_common(1)[0][0] if annotations_by_cluster[cluster] else ""
        paralog_rows.append(
            {
                "genome": genome,
                "orthogroup_or_homolog_cluster": cluster,
                "copy_number": len(members),
                "member_ids": ";".join(sorted(members)),
                "annotations": annotation,
                "functional_category": category,
                "interpretation": "putative paralogs or duplicated family members; no gene-tree validation",
            }
        )
    paralog_rows.sort(key=lambda r: (-int(r["copy_number"]), str(r["genome"]), str(r["orthogroup_or_homolog_cluster"])))
    atomic_write_tsv(
        ROOT / "06_tables/putative_paralogous_features.tsv",
        paralog_rows,
        list(paralog_rows[0]) if paralog_rows else ["genome", "orthogroup_or_homolog_cluster", "copy_number", "member_ids", "annotations", "functional_category", "interpretation"],
    )

    bgc_gene_rows: list[dict[str, object]] = []
    for cluster, agg in bgc_gene_aggregate.items():
        genomes_set = agg["genomes"]
        gene_name = agg["genes"].most_common(1)[0][0] if agg["genes"] else cluster
        product = agg["products"].most_common(1)[0][0] if agg["products"] else cluster_annotation.get(cluster, ("", ""))[1]
        bgc_gene_rows.append(
            {
                "gene": gene_name,
                "panaroo_cluster": cluster,
                "annotation": product,
                "BGC_type": ";".join(sorted(agg["types"])),
                "number_genomes": len(genomes_set),
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * len(genomes_set) / N_EVALUABLE,
                "number_regions": len(agg["regions"]),
                "functional_role": "Prokka annotation for a gene located within a predicted BGC region; biosynthetic role not established",
            }
        )
    bgc_gene_rows.sort(key=lambda r: (-int(r["number_genomes"]), str(r["panaroo_cluster"])))
    atomic_write_tsv(ROOT / "06_tables/BGC_gene_prevalence.tsv", bgc_gene_rows, list(bgc_gene_rows[0]))

    burden = pd.read_csv(ROOT / "03_matrices/BGC_burden_by_genome.tsv", sep="\t", dtype=str, keep_default_na=False)
    burden["BGC_region_gene_count"] = burden["genome_id"].map(lambda g: len(bgc_gene_ids_by_genome.get(g, set())))
    burden["biosynthetic_gene_count_if_available"] = ""
    atomic_write_tsv(ROOT / "03_matrices/BGC_burden_by_genome.tsv", burden.to_dict("records"), list(burden.columns))

    print(
        f"Feature mapping: {len(target_loci)} Prokka loci targeted, {len(locus_to_cluster)} mapped to Panaroo; "
        f"{amr_sequence_counter} Tier-A AMR-associated Prokka sequences exported"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
