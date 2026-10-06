#!/usr/bin/env python3
"""Feature-level AMP/AMR/BGC co-occurrence tests with BH-FDR correction."""

from __future__ import annotations

import argparse
import math
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from analysis_common import ROOT, atomic_write_tsv, bh_fdr


N = 1350
MIN_PRESENT = 14
MAX_PRESENT = N - MIN_PRESENT


def read_binary(path, prefix=""):
    frame = pd.read_csv(path, sep="\t", index_col="genome_id")
    frame = frame.sort_index().astype(np.int64)
    if prefix:
        frame.columns = [f"{prefix}{column}" for column in frame.columns]
    return frame


def association_rows(
    left: pd.DataFrame,
    right: pd.DataFrame,
    left_category: str,
    right_category: str,
    annotations: dict[str, str],
    symmetric: bool = False,
) -> list[dict[str, object]]:
    if not left.index.equals(right.index):
        raise ValueError(f"Genome order mismatch for {left_category}-{right_category}")
    x = left.to_numpy(dtype=np.int64)
    y = right.to_numpy(dtype=np.int64)
    intersections = x.T @ y
    left_n = x.sum(axis=0)
    right_n = y.sum(axis=0)
    pairs = combinations(range(left.shape[1]), 2) if symmetric else (
        (i, j) for i in range(left.shape[1]) for j in range(right.shape[1])
    )
    rows: list[dict[str, object]] = []
    for i, j in pairs:
        a = int(intersections[i, j])
        b = int(left_n[i] - a)
        c = int(right_n[j] - a)
        d = N - a - b - c
        denominator = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
        phi = (a * d - b * c) / denominator if denominator else float("nan")
        union = a + b + c
        jaccard = a / union if union else float("nan")
        fisher_or, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        corrected_or = aa * dd / (bb * cc)
        se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
        ci_low = math.exp(math.log(corrected_or) - 1.96 * se)
        ci_high = math.exp(math.log(corrected_or) + 1.96 * se)
        feature_1 = str(left.columns[i])
        feature_2 = str(right.columns[j])
        rows.append(
            {
                "feature_1": feature_1,
                "category_1": left_category,
                "annotation_1": annotations.get(feature_1, ""),
                "n_genomes_1": int(left_n[i]),
                "prevalence_percent_1": 100 * int(left_n[i]) / N,
                "feature_2": feature_2,
                "category_2": right_category,
                "annotation_2": annotations.get(feature_2, ""),
                "n_genomes_2": int(right_n[j]),
                "prevalence_percent_2": 100 * int(right_n[j]) / N,
                "n_both": a,
                "n_feature_1_only": b,
                "n_feature_2_only": c,
                "n_neither": d,
                "Jaccard_index": jaccard,
                "Phi_coefficient": phi,
                "Fisher_odds_ratio": fisher_or,
                "odds_ratio_Haldane_Anscombe": corrected_or,
                "OR_ci_95_low_Haldane_Anscombe": ci_low,
                "OR_ci_95_high_Haldane_Anscombe": ci_high,
                "raw_p": float(pvalue),
                "q_value": "",
                "relationship": "",
                "network_eligible": "",
                "null_hypothesis": "feature presences are independent across evaluable genomes",
                "n_evaluable": N,
                "minimum_prevalence_filter": f"{MIN_PRESENT} genomes ({100*MIN_PRESENT/N:.3f}%)",
                "maximum_prevalence_filter": f"{MAX_PRESENT} genomes ({100*MAX_PRESENT/N:.3f}%)",
                "interpretation_scope": "statistical association only; not causality, physical linkage, or interaction",
            }
        )
    qvalues = bh_fdr([float(row["raw_p"]) for row in rows])
    for row, qvalue in zip(rows, qvalues):
        row["q_value"] = qvalue
        phi = float(row["Phi_coefficient"])
        row["relationship"] = (
            "positive association"
            if qvalue < 0.05 and phi > 0
            else "negative association"
            if qvalue < 0.05 and phi < 0
            else "no statistically supported association"
        )
        row["network_eligible"] = "yes" if qvalue < 0.05 and abs(phi) >= 0.20 else "no"
    rows.sort(key=lambda row: (float(row["q_value"]), -abs(float(row["Phi_coefficient"])), str(row["feature_1"]), str(row["feature_2"])))
    return rows


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    amp_all = read_binary(ROOT / "03_matrices/AMP_presence_absence.tsv")
    amr_all = read_binary(ROOT / "03_matrices/AMR_presence_absence.tsv")
    bgc_family_all = read_binary(ROOT / "03_matrices/BGC_family_presence_absence.tsv", "BGC_FAMILY::")
    bgc_type_all = read_binary(ROOT / "03_matrices/BGC_type_presence_absence.tsv", "BGC_TYPE::")
    if not (amp_all.index.equals(amr_all.index) and amp_all.index.equals(bgc_family_all.index) and amp_all.index.equals(bgc_type_all.index)):
        raise ValueError("Genome rows differ among feature matrices")

    def filter_frame(frame: pd.DataFrame) -> pd.DataFrame:
        counts = frame.sum(axis=0)
        return frame.loc[:, (counts >= MIN_PRESENT) & (counts <= MAX_PRESENT)]

    amp = filter_frame(amp_all)
    amr = filter_frame(amr_all)
    bgc = filter_frame(pd.concat([bgc_family_all, bgc_type_all], axis=1))

    annotations: dict[str, str] = {}
    amp_meta = pd.read_csv(ROOT / "06_tables/AMP_prevalence_all.tsv", sep="\t", dtype=str, keep_default_na=False)
    annotations.update(dict(zip(amp_meta["AMP"], amp_meta["annotation"])))
    amr_meta = pd.read_csv(ROOT / "06_tables/AMR_family_prevalence_all.tsv", sep="\t", dtype=str, keep_default_na=False)
    annotations.update(dict(zip(amr_meta["gene_family"], amr_meta["annotation"])))
    family_meta = pd.read_csv(ROOT / "06_tables/BGC_gene_cluster_families.tsv", sep="\t", dtype=str, keep_default_na=False)
    annotations.update(
        {f"BGC_FAMILY::{row.family_id}": row.predicted_product for row in family_meta.itertuples(index=False)}
    )
    type_meta = pd.read_csv(ROOT / "06_tables/BGC_type_prevalence.tsv", sep="\t", dtype=str, keep_default_na=False)
    annotations.update({f"BGC_TYPE::{row.BGC_type}": row.BGC_type for row in type_meta.itertuples(index=False)})

    filter_rows = []
    for category, frame_all, frame_filtered in (
        ("AMP", amp_all, amp),
        ("AMR", amr_all, amr),
        ("BGC", pd.concat([bgc_family_all, bgc_type_all], axis=1), bgc),
    ):
        counts = frame_all.sum(axis=0)
        for feature in frame_all.columns:
            filter_rows.append(
                {
                    "feature": feature,
                    "functional_category": category,
                    "n_genomes": int(counts[feature]),
                    "prevalence_percent": 100 * int(counts[feature]) / N,
                    "included": "yes" if feature in frame_filtered.columns else "no",
                    "filter": f"{MIN_PRESENT}-{MAX_PRESENT} of {N} genomes",
                }
            )
    atomic_write_tsv(ROOT / "04_statistics/feature_association_filter.tsv", filter_rows, list(filter_rows[0]))

    analyses = [
        ("AMP_ARG_associations.tsv", amp, amr, "AMP", "AMR", False),
        ("AMP_BGC_associations.tsv", amp, bgc, "AMP", "BGC", False),
        ("ARG_BGC_associations.tsv", amr, bgc, "AMR", "BGC", False),
        ("AMR_gene_cooccurrence.tsv", amr, amr, "AMR", "AMR", True),
    ]
    summary = []
    for filename, left, right, left_category, right_category, symmetric in analyses:
        rows = association_rows(left, right, left_category, right_category, annotations, symmetric=symmetric)
        atomic_write_tsv(ROOT / f"04_statistics/{filename}", rows, list(rows[0]))
        summary.append(
            {
                "analysis": filename.removesuffix(".tsv"),
                "features_category_1": left.shape[1],
                "features_category_2": right.shape[1],
                "tests": len(rows),
                "FDR_significant": sum(float(row["q_value"]) < 0.05 for row in rows),
                "positive_FDR_significant": sum(row["relationship"] == "positive association" for row in rows),
                "negative_FDR_significant": sum(row["relationship"] == "negative association" for row in rows),
                "network_eligible_q_lt_0.05_abs_phi_ge_0.20": sum(row["network_eligible"] == "yes" for row in rows),
                "FDR_scope": "BH correction within this analysis family",
            }
        )
    atomic_write_tsv(ROOT / "04_statistics/feature_association_summary.tsv", summary, list(summary[0]))
    print(
        f"Feature associations: AMP={amp.shape[1]}, AMR={amr.shape[1]}, BGC={bgc.shape[1]} filtered features; "
        f"tests={sum(row['tests'] for row in summary)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
