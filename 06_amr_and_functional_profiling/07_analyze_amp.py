#!/usr/bin/env python3
"""Analyze AMP sequence families, prevalence, burden, and molecular properties."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from analysis_common import (
    ROOT,
    atomic_write_tsv,
    bh_fdr,
    describe,
    prevalence_band,
)


THRESHOLDS = (90, 70, 50)
N_EVALUABLE = 1350


def read_cluster_pairs(threshold: int) -> list[tuple[str, str]]:
    path = ROOT / f"03_matrices/mmseqs_amp_{threshold}/amp_{threshold}_cluster.tsv"
    pairs: list[tuple[str, str]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                pairs.append((row[0], row[1]))
    return pairs


def read_hit_identities(threshold: int) -> dict[tuple[str, str], float]:
    path = ROOT / f"03_matrices/mmseqs_amp_{threshold}/amp_{threshold}_representative_hits.tsv"
    result: dict[tuple[str, str], float] = {}
    if not path.exists():
        return result
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 3:
                key = (row[0], row[1])
                value = float(row[2])
                result[key] = max(value, result.get(key, 0.0))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    master = pd.read_csv(ROOT / "02_cleaned/AMP_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    master_by_sequence = master.set_index("sequence_id", drop=False).to_dict("index")
    mapped_gene_context: dict[str, dict[str, str]] = {}
    mapping_columns = [
        "funcscan_feature",
        "funcscan_category",
        "prokka_gene_name",
        "prokka_product",
        "match_method",
        "panaroo_cluster",
        "panaroo_annotation",
    ]
    for chunk in pd.read_csv(
        ROOT / "02_cleaned/funcscan_to_panaroo_mapping.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=mapping_columns,
        chunksize=250000,
    ):
        for row in chunk[chunk["funcscan_category"] == "AMP"].itertuples(index=False):
            mapped_gene_context[row.funcscan_feature] = {
                "gene_context": row.prokka_gene_name or row.prokka_product or row.panaroo_annotation,
                "panaroo_cluster": row.panaroo_cluster,
                "mapping_method": row.match_method,
            }
    assignment_rows: list[dict[str, object]] = []
    cluster_rows: list[dict[str, object]] = []
    representative_rows: list[dict[str, object]] = []
    primary_cluster_for_sequence: dict[str, str] = {}

    for threshold in THRESHOLDS:
        pairs = read_cluster_pairs(threshold)
        raw_groups: dict[str, list[str]] = defaultdict(list)
        for representative, member in pairs:
            raw_groups[representative].append(member)
        ordered = sorted(
            raw_groups,
            key=lambda rep: (
                -len({master_by_sequence[x]["genome_id"] for x in raw_groups[rep]}),
                -len(raw_groups[rep]),
                rep,
            ),
        )
        stable = {rep: f"AMP{threshold}_C{i:04d}" for i, rep in enumerate(ordered, start=1)}
        identities = read_hit_identities(threshold)
        for representative in ordered:
            cluster_id = stable[representative]
            members = raw_groups[representative]
            records = [master_by_sequence[x] for x in members]
            genomes = {r["genome_id"] for r in records}
            annotations = [r["functional_annotation"] for r in records if r["functional_annotation"] != "unannotated predicted AMP"]
            contexts = Counter(
                mapped_gene_context.get(r["feature_id"], {}).get("gene_context", "")
                for r in records
                if mapped_gene_context.get(r["feature_id"], {}).get("gene_context", "")
            )
            context_annotation = contexts.most_common(1)[0][0] if contexts else ""
            annotation = Counter(annotations).most_common(1)[0][0] if annotations else "unannotated predicted AMP family"
            if context_annotation and not annotations:
                annotation = f"unannotated predicted AMP family; genomic gene context: {context_annotation}"
            panaroo_clusters = sorted(
                {
                    mapped_gene_context.get(r["feature_id"], {}).get("panaroo_cluster", "")
                    for r in records
                }
                - {""}
            )
            identity_values = [identities[(member, representative)] for member in members if (member, representative) in identities]
            high_conf = sum(r["high_confidence"] == "yes" for r in records)
            database = sum(r["database_match"] == "yes" for r in records)
            predictor_combinations = Counter(r["predictor_combination"] for r in records)
            representative_record = master_by_sequence[representative]
            summary = {
                "identity_threshold_percent": threshold,
                "alignment_coverage_threshold": 0.8,
                "cluster_id": cluster_id,
                "representative_sequence_id": representative,
                "representative_feature_id": representative_record["feature_id"],
                "representative_sequence": representative_record["aa_sequence"],
                "number_of_sequences": len(members),
                "number_unique_sequences": len({r["sequence_sha256"] for r in records}),
                "number_of_genomes": len(genomes),
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * len(genomes) / N_EVALUABLE,
                "mean_identity_to_representative_percent": np.mean(identity_values) if identity_values else "",
                "identity_observations": len(identity_values),
                "annotation": annotation,
                "genomic_gene_context_annotation": context_annotation,
                "Panaroo_clusters": ";".join(panaroo_clusters),
                "Panaroo_cluster_count": len(panaroo_clusters),
                "Prokka_Panaroo_mapped_predictions": sum(
                    bool(mapped_gene_context.get(r["feature_id"], {}).get("panaroo_cluster", "")) for r in records
                ),
                "DRAMP_annotation": ";".join(sorted({r["DRAMP_ID"] for r in records if r["DRAMP_ID"]})),
                "SwissProt_annotation": ";".join(sorted({r["Swiss_Prot_Entry"] for r in records if r["Swiss_Prot_Entry"]})),
                "predictor_support": ";".join(f"{k}:{v}" for k, v in sorted(predictor_combinations.items())),
                "high_confidence_predictions": high_conf,
                "database_supported_predictions": database,
                "prevalence_band": prevalence_band(100 * len(genomes) / N_EVALUABLE),
            }
            cluster_rows.append(summary)
            representative_rows.append(
                {
                    "identity_threshold_percent": threshold,
                    "cluster_id": cluster_id,
                    "representative_sequence_id": representative,
                    "representative_feature_id": representative_record["feature_id"],
                    "representative_sequence": representative_record["aa_sequence"],
                    "number_of_sequences": len(members),
                    "number_of_genomes": len(genomes),
                    "annotation": annotation,
                }
            )
            for member in members:
                record = master_by_sequence[member]
                assignment_rows.append(
                    {
                        "identity_threshold_percent": threshold,
                        "alignment_coverage_threshold": 0.8,
                        "cluster_id": cluster_id,
                        "representative_sequence_id": representative,
                        "sequence_id": member,
                        "feature_id": record["feature_id"],
                        "genome_id": record["genome_id"],
                        "sequence_sha256": record["sequence_sha256"],
                        "identity_to_representative_percent": identities.get((member, representative), ""),
                    }
                )
                if threshold == 50:
                    primary_cluster_for_sequence[member] = cluster_id

    atomic_write_tsv(
        ROOT / "03_matrices/AMP_sequence_clusters.tsv",
        assignment_rows,
        list(assignment_rows[0]),
    )
    atomic_write_tsv(ROOT / "06_tables/AMP_cluster_summary.tsv", cluster_rows, list(cluster_rows[0]))
    atomic_write_tsv(
        ROOT / "06_tables/AMP_cluster_representatives.tsv",
        representative_rows,
        list(representative_rows[0]),
    )

    master["AMP_cluster_50"] = master["sequence_id"].map(primary_cluster_for_sequence)
    primary_summary = {
        row["cluster_id"]: row
        for row in cluster_rows
        if int(row["identity_threshold_percent"]) == 50
    }
    prevalence: list[dict[str, object]] = []
    high_prevalence: list[dict[str, object]] = []
    for cluster_id, group in master.groupby("AMP_cluster_50", sort=False):
        summary = primary_summary[cluster_id]
        genomes = group["genome_id"].nunique()
        positives = group.groupby("genome_id").size()
        high_group = group[group["high_confidence"] == "yes"]
        record = {
            "AMP": cluster_id,
            "annotation": summary["annotation"],
            "n_genomes": genomes,
            "n_evaluable": N_EVALUABLE,
            "prevalence_percent": 100 * genomes / N_EVALUABLE,
            "total_occurrences": len(group),
            "mean_copies_per_positive_genome": positives.mean(),
            "median_copies_per_positive_genome": positives.median(),
            "prediction_support": f"high_confidence={len(high_group)}/{len(group)}",
            "high_confidence_occurrences": len(high_group),
            "high_confidence_genomes": high_group["genome_id"].nunique(),
            "database_match": "yes" if (group["database_match"] == "yes").any() else "no",
            "cluster_id": cluster_id,
            "identity_threshold_percent": 50,
            "alignment_coverage_threshold": 0.8,
            "prevalence_band": prevalence_band(100 * genomes / N_EVALUABLE),
        }
        prevalence.append(record)
        if len(high_group):
            high_record = dict(record)
            high_record["n_genomes"] = high_group["genome_id"].nunique()
            high_record["prevalence_percent"] = 100 * high_record["n_genomes"] / N_EVALUABLE
            high_record["total_occurrences"] = len(high_group)
            high_record["mean_copies_per_positive_genome"] = high_group.groupby("genome_id").size().mean()
            high_record["median_copies_per_positive_genome"] = high_group.groupby("genome_id").size().median()
            high_prevalence.append(high_record)
    prevalence.sort(key=lambda r: (-int(r["n_genomes"]), -int(r["total_occurrences"]), str(r["AMP"])))
    high_prevalence.sort(key=lambda r: (-int(r["n_genomes"]), -int(r["total_occurrences"]), str(r["AMP"])))
    atomic_write_tsv(ROOT / "06_tables/AMP_prevalence_all.tsv", prevalence, list(prevalence[0]))
    atomic_write_tsv(ROOT / "06_tables/AMP_top50_prevalent.tsv", prevalence[:50], list(prevalence[0]))
    atomic_write_tsv(
        ROOT / "06_tables/AMP_high_confidence_prevalence.tsv",
        high_prevalence,
        list(prevalence[0]),
    )

    genomes = sorted(master["genome_id"].unique())
    clusters = [r["AMP"] for r in prevalence]
    present = set(zip(master["genome_id"], master["AMP_cluster_50"]))
    matrix_rows = []
    for genome in genomes:
        row = {"genome_id": genome}
        row.update({cluster: int((genome, cluster) in present) for cluster in clusters})
        matrix_rows.append(row)
    atomic_write_tsv(
        ROOT / "03_matrices/AMP_presence_absence.tsv",
        matrix_rows,
        ["genome_id", *clusters],
    )

    burden_rows = []
    for genome in genomes:
        group = master[master["genome_id"] == genome]
        burden_rows.append(
            {
                "genome_id": genome,
                "AMP_total": len(group),
                "AMP_unique": group["AMP_cluster_50"].nunique(),
                "AMP_high_confidence": (group["high_confidence"] == "yes").sum(),
                "AMP_high_confidence_unique": group.loc[group["high_confidence"] == "yes", "AMP_cluster_50"].nunique(),
                "AMP_predictor_consensus_fraction": (group["high_confidence"] == "yes").mean(),
            }
        )
    atomic_write_tsv(ROOT / "03_matrices/AMP_burden_by_genome.tsv", burden_rows, list(burden_rows[0]))

    cluster_prevalence_map = {r["AMP"]: float(r["prevalence_percent"]) for r in prevalence}
    master["cluster_prevalence"] = master["AMP_cluster_50"].map(cluster_prevalence_map)
    numeric_columns = [
        "sequence_length_aa",
        "molecular_weight",
        "hydrophobicity",
        "isoelectric_point",
        "helix_fraction",
        "turn_fraction",
        "sheet_fraction",
    ]
    for column in numeric_columns:
        master[column] = pd.to_numeric(master[column], errors="coerce")
    subsets = {
        "all_predictions": master,
        "high_confidence": master[master["high_confidence"] == "yes"],
        "tier_1": master[master["evidence_tier"].str.startswith("Tier 1")],
        "tier_2": master[master["evidence_tier"].str.startswith("Tier 2")],
        "tier_3": master[master["evidence_tier"].str.startswith("Tier 3")],
        "high_prevalence_clusters_ge50pct": master[master["cluster_prevalence"] >= 50],
        "rare_clusters_lt10pct": master[master["cluster_prevalence"] < 10],
    }
    property_rows: list[dict[str, object]] = []
    for subset_name, subset in subsets.items():
        for column in numeric_columns:
            vals = subset[column].dropna().astype(float).tolist()
            stats = describe(vals)
            property_rows.append({"subset": subset_name, "property": column, **stats})
        property_rows.append(
            {
                "subset": subset_name,
                "property": "transporter_protein_yes_fraction",
                "n": len(subset),
                "mean": (subset["transporter_protein"].str.lower() == "yes").mean() if len(subset) else "",
                "median": "",
                "sd": "",
                "q1": "",
                "q3": "",
                "iqr": "",
                "min": "",
                "max": "",
                "p05": "",
                "p95": "",
                "p99": "",
            }
        )
    atomic_write_tsv(
        ROOT / "04_statistics/AMP_physicochemical_summary.tsv",
        property_rows,
        list(property_rows[0]),
    )
    tests: list[dict[str, object]] = []
    comparisons = [
        ("high_confidence", subsets["high_confidence"], "tier_3", subsets["tier_3"]),
        (
            "high_prevalence_clusters_ge50pct",
            subsets["high_prevalence_clusters_ge50pct"],
            "rare_clusters_lt10pct",
            subsets["rare_clusters_lt10pct"],
        ),
    ]
    for left_name, left, right_name, right in comparisons:
        for column in numeric_columns:
            x = left[column].dropna().astype(float)
            y = right[column].dropna().astype(float)
            if len(x) and len(y):
                result = mannwhitneyu(x, y, alternative="two-sided")
                tests.append(
                    {
                        "comparison": f"{left_name}_vs_{right_name}",
                        "property": column,
                        "null_hypothesis": "the two groups have the same distribution",
                        "n_group_1": len(x),
                        "n_group_2": len(y),
                        "median_group_1": x.median(),
                        "median_group_2": y.median(),
                        "median_difference": x.median() - y.median(),
                        "U_statistic": result.statistic,
                        "raw_p": result.pvalue,
                    }
                )
    qvalues = bh_fdr([float(r["raw_p"]) for r in tests])
    for row, qvalue in zip(tests, qvalues):
        row["q_value"] = qvalue
    atomic_write_tsv(ROOT / "04_statistics/AMP_physicochemical_tests.tsv", tests, list(tests[0]))
    print(f"AMP analysis: {len(prevalence)} primary 50%-identity clusters; {len(high_prevalence)} clusters with high-confidence members")
    return 0


if __name__ == "__main__":
    sys.exit(main())
