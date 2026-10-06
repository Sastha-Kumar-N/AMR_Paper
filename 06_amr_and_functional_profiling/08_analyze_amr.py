#!/usr/bin/env python3
"""Analyze primary (Tier-A) AMR calls without equating genotype with phenotype."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter

import pandas as pd

from analysis_common import ROOT, atomic_write_tsv, describe, load_clinical_metadata, prevalence_band


N_EVALUABLE = 1350


def split_terms(value: str) -> list[str]:
    return [x.strip() for x in re.split(r";\s*", value or "") if x.strip()]


def consolidate_evidence_streams(records: pd.DataFrame, feature_col: str) -> pd.DataFrame:
    """Avoid adding copy counts across callers whose locus IDs are incomparable.

    Within each genome/feature, retain the evidence stream with the largest
    number of distinct loci; ties favor coordinate-supported Prokka mappings.
    This estimates copies conservatively without pretending ordinal IDs can be
    joined to coordinate IDs.
    """
    priority = {"prokka_coordinate": 0, "coordinate_only": 1, "prediction_sequence": 2, "unresolved": 3}
    frame = records.copy()
    frame["evidence_stream"] = frame["occurrence_id"].map(
        lambda value: (
            "prokka_coordinate"
            if value.startswith("prokka:")
            else "coordinate_only"
            if value.startswith("coordinate:")
            else "prediction_sequence"
            if value.startswith("prediction_sequence:")
            else "unresolved"
        )
    )
    output: list[dict[str, object]] = []
    for (genome, feature), group in frame.groupby(["genome_id", feature_col], sort=False):
        stream_counts = group.groupby("evidence_stream")["occurrence_id"].nunique().to_dict()
        selected = sorted(stream_counts, key=lambda stream: (-stream_counts[stream], priority[stream]))[0]
        selected_rows = (
            group[group["evidence_stream"] == selected]
            .sort_values(["occurrence_id", "gene_family_tool_support_count", "prediction_tool"], ascending=[True, False, True])
            .drop_duplicates("occurrence_id")
        )
        for index, row in enumerate(selected_rows.to_dict("records"), start=1):
            row["source_occurrence_id"] = row["occurrence_id"]
            row["occurrence_id"] = f"conservative:{genome}:{feature}:{index}"
            row["copy_number_evidence_stream"] = selected
            row["copy_number_rule"] = "maximum distinct loci in any single evidence stream; streams are not added"
            output.append(row)
    return pd.DataFrame(output)


def prevalence_table(records: pd.DataFrame, feature_col: str, label: str) -> pd.DataFrame:
    dedup = records.drop_duplicates(["genome_id", feature_col, "occurrence_id"])
    rows: list[dict[str, object]] = []
    for feature, group in dedup.groupby(feature_col, sort=False):
        if not feature:
            continue
        counts = group.groupby("genome_id")["occurrence_id"].nunique()
        annotations = Counter(x for x in group["functional_annotation"] if x)
        classes = sorted({x for value in group["normalized_drug_classes"] for x in split_terms(value)})
        mechanisms = sorted({x for value in group["resistance_mechanism"] for x in split_terms(value)})
        tools = sorted({x for x in group["prediction_tool"] if x})
        n_genomes = int(counts.size)
        total = int(counts.sum())
        rows.append(
            {
                label: feature,
                "annotation": annotations.most_common(1)[0][0] if annotations else "",
                "n_genomes": n_genomes,
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * n_genomes / N_EVALUABLE,
                "total_occurrences": total,
                "mean_copies_per_positive_genome": float(counts.mean()),
                "median_copies_per_positive_genome": float(counts.median()),
                "max_copies_per_genome": int(counts.max()),
                "drug_class": ";".join(classes),
                "mechanism": ";".join(mechanisms),
                "prediction_tools": ";".join(tools),
                "prevalence_band": prevalence_band(100 * n_genomes / N_EVALUABLE),
                "interpretation_scope": "genotypic prediction only; phenotype is not established",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["n_genomes", "total_occurrences", label], ascending=[False, False, True]
    )


def term_prevalence(records: pd.DataFrame, source_col: str, label: str) -> pd.DataFrame:
    exploded: list[dict[str, str]] = []
    for row in records[["genome_id", "occurrence_id", source_col]].drop_duplicates().itertuples(index=False):
        for term in split_terms(getattr(row, source_col)):
            exploded.append({"genome_id": row.genome_id, "occurrence_id": row.occurrence_id, label: term})
    if not exploded:
        return pd.DataFrame(columns=[label, "n_genomes", "n_evaluable", "prevalence_percent", "total_occurrences"])
    frame = pd.DataFrame(exploded).drop_duplicates()
    rows = []
    for term, group in frame.groupby(label, sort=False):
        counts = group.groupby("genome_id")["occurrence_id"].nunique()
        n = int(counts.size)
        rows.append(
            {
                label: term,
                "n_genomes": n,
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * n / N_EVALUABLE,
                "total_occurrences": int(counts.sum()),
                "mean_occurrences_per_positive_genome": float(counts.mean()),
                "median_occurrences_per_positive_genome": float(counts.median()),
                "max_occurrences_per_genome": int(counts.max()),
                "prevalence_band": prevalence_band(100 * n / N_EVALUABLE),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["n_genomes", "total_occurrences", label], ascending=[False, False, True]
    )


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    clinical, _, canonical_to_raw = load_clinical_metadata()
    genomes = sorted(canonical_to_raw)
    amr = pd.read_csv(ROOT / "02_cleaned/AMR_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    primary = amr[amr["primary_analysis"] == "yes"].copy()

    feature_to_locus: dict[str, str] = {}
    mapping_method: dict[str, str] = {}
    usecols = ["funcscan_feature", "funcscan_category", "prokka_gene_id", "match_method"]
    for chunk in pd.read_csv(
        ROOT / "02_cleaned/funcscan_to_panaroo_mapping.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=usecols,
        chunksize=250000,
    ):
        for row in chunk[chunk["funcscan_category"] == "AMR"].itertuples(index=False):
            if row.prokka_gene_id:
                feature_to_locus[row.funcscan_feature] = f"prokka:{row.prokka_gene_id}"
                mapping_method[row.funcscan_feature] = row.match_method

    def fallback_locus(row: pd.Series) -> str:
        mapped = feature_to_locus.get(row["feature_id"])
        if mapped:
            return mapped
        if row["contig_id"] and row["coordinate_start"] and row["coordinate_end"]:
            return f"coordinate:{row['contig_id']}:{row['coordinate_start']}-{row['coordinate_end']}"
        if row["sequence_id"]:
            return f"prediction_sequence:{row['sequence_id']}"
        return f"unresolved_call:{row['feature_id']}"

    primary["occurrence_id"] = primary.apply(fallback_locus, axis=1)
    primary["coordinate_supported_mapping"] = primary["feature_id"].map(mapping_method).fillna("")

    # One row per genome, family and inferred locus prevents multi-tool support
    # from being miscounted as biological copy number.
    family_occ = consolidate_evidence_streams(primary, "gene_family")
    gene_occ = consolidate_evidence_streams(primary, "normalized_feature_name")

    gene_prev = prevalence_table(gene_occ, "normalized_feature_name", "gene")
    family_prev = prevalence_table(family_occ, "gene_family", "gene_family")
    atomic_write_tsv(ROOT / "06_tables/AMR_gene_prevalence_all.tsv", gene_prev.to_dict("records"), list(gene_prev.columns))
    atomic_write_tsv(ROOT / "06_tables/AMR_top50_genes.tsv", gene_prev.head(50).to_dict("records"), list(gene_prev.columns))
    atomic_write_tsv(ROOT / "06_tables/AMR_family_prevalence_all.tsv", family_prev.to_dict("records"), list(family_prev.columns))

    class_prev = term_prevalence(family_occ, "normalized_drug_classes", "drug_class")
    mechanism_prev = term_prevalence(family_occ, "resistance_mechanism", "mechanism")
    atomic_write_tsv(ROOT / "06_tables/AMR_drug_class_prevalence.tsv", class_prev.to_dict("records"), list(class_prev.columns))
    atomic_write_tsv(ROOT / "06_tables/AMR_mechanism_prevalence.tsv", mechanism_prev.to_dict("records"), list(mechanism_prev.columns))

    presence = (
        family_occ.assign(value=1)
        .pivot_table(index="genome_id", columns="gene_family", values="value", aggfunc="max", fill_value=0)
        .reindex(genomes, fill_value=0)
        .astype(int)
    )
    presence.index.name = "genome_id"
    atomic_write_tsv(
        ROOT / "03_matrices/AMR_presence_absence.tsv",
        ({"genome_id": genome, **presence.loc[genome].to_dict()} for genome in presence.index),
        ["genome_id", *presence.columns],
    )
    prevalence_n = presence.sum(axis=0)
    filtered_features = sorted(prevalence_n[(prevalence_n >= 14) & (prevalence_n <= N_EVALUABLE - 14)].index)
    atomic_write_tsv(
        ROOT / "04_statistics/AMR_association_matrix_filter.tsv",
        (
            {
                "gene_family": feature,
                "n_genomes": int(prevalence_n[feature]),
                "prevalence_percent": 100 * int(prevalence_n[feature]) / N_EVALUABLE,
                "included": "yes" if feature in filtered_features else "no",
                "threshold": "14-1336 genomes (approximately 1-99%)",
            }
            for feature in presence.columns
        ),
        ["gene_family", "n_genomes", "prevalence_percent", "included", "threshold"],
    )

    family_sets = family_occ.groupby("genome_id")["gene_family"].agg(set).to_dict()
    exact_sets = gene_occ.groupby("genome_id")["normalized_feature_name"].agg(set).to_dict()
    occurrence_counts = family_occ.groupby("genome_id").size().to_dict()
    class_sets: dict[str, set[str]] = {g: set() for g in genomes}
    mechanism_sets: dict[str, set[str]] = {g: set() for g in genomes}
    for row in family_occ.itertuples(index=False):
        class_sets[row.genome_id].update(split_terms(row.normalized_drug_classes))
        mechanism_sets[row.genome_id].update(split_terms(row.resistance_mechanism))
    raw_counts = primary.groupby("genome_id").size().to_dict()
    primary["coordinate_locus"] = primary["feature_id"].map(feature_to_locus).fillna("")
    mapped_counts = primary[primary["coordinate_locus"] != ""].groupby("genome_id")["coordinate_locus"].nunique().to_dict()
    burden_rows = []
    for genome in genomes:
        burden_rows.append(
            {
                "genome_id": genome,
                "total_ARG_count": int(occurrence_counts.get(genome, 0)),
                "unique_ARG_count": len(family_sets.get(genome, set())),
                "unique_specific_ARG_names": len(exact_sets.get(genome, set())),
                "drug_class_count": len(class_sets[genome]),
                "resistance_mechanism_count": len(mechanism_sets[genome]),
                "coordinate_supported_ARG_loci": mapped_counts.get(genome, 0),
                "raw_Tier_A_tool_calls": raw_counts.get(genome, 0),
                "scope": "Tier-A family occurrences consolidated by the maximum count in one evidence stream; not an MDR or phenotype classification",
            }
        )
    atomic_write_tsv(ROOT / "03_matrices/AMR_burden_by_genome.tsv", burden_rows, list(burden_rows[0]))
    summary_rows = []
    for column in (
        "total_ARG_count",
        "unique_ARG_count",
        "unique_specific_ARG_names",
        "drug_class_count",
        "resistance_mechanism_count",
        "coordinate_supported_ARG_loci",
        "raw_Tier_A_tool_calls",
    ):
        stats = describe([float(row[column]) for row in burden_rows])
        summary_rows.append({"metric": column, **stats})
    atomic_write_tsv(ROOT / "06_tables/AMR_burden_summary.tsv", summary_rows, list(summary_rows[0]))

    functional_rows: list[dict[str, object]] = []
    for definition, frame, feature_col in (
        ("gene_family", family_occ, "gene_family"),
        ("drug_class", None, "drug_class"),
        ("mechanism", None, "mechanism"),
    ):
        if definition == "gene_family":
            iterator = ((name, group) for name, group in frame.groupby(feature_col, sort=False))
        else:
            source = "normalized_drug_classes" if definition == "drug_class" else "resistance_mechanism"
            expanded = []
            for row in family_occ.itertuples(index=False):
                for term in split_terms(getattr(row, source)):
                    expanded.append(
                        {
                            "name": term,
                            "genome_id": row.genome_id,
                            "occurrence_id": row.occurrence_id,
                            "gene_family": row.gene_family,
                            "functional_annotation": row.functional_annotation,
                        }
                    )
            expanded_frame = pd.DataFrame(expanded)
            iterator = ((name, group) for name, group in expanded_frame.groupby("name", sort=False))
        for name, group in iterator:
            genomes_set = set(group["genome_id"])
            loci = set(zip(group["genome_id"], group["occurrence_id"]))
            members = sorted(set(group["gene_family"]))
            annotations = Counter(x for x in group["functional_annotation"] if x)
            functional_rows.append(
                {
                    "cluster_definition": definition,
                    "functional_cluster": name,
                    "member_gene_families": ";".join(members),
                    "representative_annotation": annotations.most_common(1)[0][0] if annotations else "",
                    "number_occurrences": len(loci),
                    "number_genomes": len(genomes_set),
                    "n_evaluable": N_EVALUABLE,
                    "prevalence_percent": 100 * len(genomes_set) / N_EVALUABLE,
                }
            )
    functional_rows.sort(key=lambda r: (-int(r["number_genomes"]), str(r["cluster_definition"]), str(r["functional_cluster"])))
    atomic_write_tsv(ROOT / "06_tables/AMR_functional_clusters.tsv", functional_rows, list(functional_rows[0]))

    print(
        f"AMR primary analysis: {len(primary)} Tier-A tool calls, {len(family_prev)} families, "
        f"{len(gene_prev)} specific names, {len(filtered_features)} association-evaluable families"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
