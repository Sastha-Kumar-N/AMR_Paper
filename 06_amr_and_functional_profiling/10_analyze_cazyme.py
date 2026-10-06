#!/usr/bin/env python3
"""Compute CAZyme family/class prevalence, abundance, burdens, and matrices."""

from __future__ import annotations

import argparse
import sys
from collections import Counter

import pandas as pd

from analysis_common import ROOT, atomic_write_tsv, prevalence_band


N_EVALUABLE = 1350


def prevalence_table(df: pd.DataFrame, column: str, label: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for value, group in df.groupby(column):
        genomes = group["genome_id"].nunique()
        per_genome = group.groupby("genome_id").size()
        rows.append(
            {
                label: value,
                "n_genomes": genomes,
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * genomes / N_EVALUABLE,
                "total_occurrences": len(group),
                "mean_copies_per_positive_genome": per_genome.mean(),
                "median_copies_per_positive_genome": per_genome.median(),
                "prevalence_band": prevalence_band(100 * genomes / N_EVALUABLE),
            }
        )
    rows.sort(key=lambda r: (-int(r["n_genomes"]), -int(r["total_occurrences"]), str(r[label])))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    master = pd.read_csv(ROOT / "02_cleaned/CAZyme_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    primary = master[master["analysis_tier"] == "primary_consensus"].copy()
    family_rows = prevalence_table(primary, "CAZyme_family", "CAZyme_family")
    class_rows = prevalence_table(primary, "CAZyme_class", "CAZyme_class")
    subfamily = primary[primary["CAZyme_subfamily"] != ""]
    subfamily_rows = prevalence_table(subfamily, "CAZyme_subfamily", "CAZyme_subfamily") if len(subfamily) else []
    atomic_write_tsv(ROOT / "06_tables/CAZyme_family_prevalence.tsv", family_rows, list(family_rows[0]))
    atomic_write_tsv(ROOT / "06_tables/CAZyme_class_prevalence.tsv", class_rows, list(class_rows[0]))
    if subfamily_rows:
        atomic_write_tsv(ROOT / "06_tables/CAZyme_subfamily_prevalence.tsv", subfamily_rows, list(subfamily_rows[0]))

    cgc = pd.read_csv(ROOT / "02_cleaned/CAZyme_CGC_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    genomes = sorted(master["genome_id"].unique())
    burden_rows = []
    for genome in genomes:
        group = primary[primary["genome_id"] == genome]
        cgc_group = cgc[cgc["genome_id"] == genome]
        burden_rows.append(
            {
                "genome_id": genome,
                "CAZyme_total": len(group),
                "CAZyme_gene_count": group["sequence_id"].nunique(),
                "CAZyme_families": group["CAZyme_family"].nunique(),
                "CAZyme_subfamilies": group.loc[group["CAZyme_subfamily"] != "", "CAZyme_subfamily"].nunique(),
                "CAZyme_classes": group["CAZyme_class"].nunique(),
                "CAZyme_CGC_count": len(cgc_group),
                "CAZyme_single_tool_only_records": len(master[(master["genome_id"] == genome) & (master["analysis_tier"] == "single_tool_only")]),
            }
        )
    atomic_write_tsv(ROOT / "03_matrices/CAZyme_burden_by_genome.tsv", burden_rows, list(burden_rows[0]))

    family_names = [row["CAZyme_family"] for row in family_rows]
    present = set(zip(primary["genome_id"], primary["CAZyme_family"]))
    matrix = []
    for genome in genomes:
        row = {"genome_id": genome}
        row.update({family: int((genome, family) in present) for family in family_names})
        matrix.append(row)
    atomic_write_tsv(ROOT / "03_matrices/CAZyme_presence_absence.tsv", matrix, ["genome_id", *family_names])
    print(
        f"CAZyme analysis: {len(primary)} consensus family occurrences; "
        f"{len(family_rows)} families, {len(subfamily_rows)} subfamilies, {len(cgc)} CGCs"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
