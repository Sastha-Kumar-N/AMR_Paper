#!/usr/bin/env python3
"""Genome-level burden correlations and CDS-adjusted count regressions."""

from __future__ import annotations

import argparse
import math
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr, rankdata, spearmanr

from analysis_common import ROOT, atomic_write_tsv, bh_fdr


BOOTSTRAPS = 1000
SEED = 13502026


def effect_label(rho: float) -> str:
    value = abs(rho)
    if value < 0.10:
        return "negligible"
    if value < 0.30:
        return "weak"
    if value < 0.50:
        return "moderate"
    return "strong"


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float, int]:
    values = []
    n = len(x)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(BOOTSTRAPS):
            indices = rng.integers(0, n, n)
            rho = spearmanr(x[indices], y[indices]).statistic
            if np.isfinite(rho):
                values.append(float(rho))
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)), len(values)


def residualized_rank(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    response = rankdata(values, method="average")
    design = sm.add_constant(rankdata(covariate, method="average"))
    return np.asarray(sm.OLS(response, design).fit().resid)


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    data = pd.read_csv(ROOT / "03_matrices/genome_functional_burden_matrix.tsv", sep="\t")
    variables = [
        "AMP_total",
        "AMP_high_confidence",
        "ARG_total",
        "ARG_unique",
        "AMR_drug_classes",
        "BGC_total",
        "BGC_types",
        "CAZyme_total",
    ]
    rng = np.random.default_rng(SEED)
    correlation_rows: list[dict[str, object]] = []
    partial_rows: list[dict[str, object]] = []
    for left_index, left in enumerate(variables):
        for right in variables[left_index + 1 :]:
            subset = data[[left, right, "total_CDS"]].dropna()
            x = subset[left].to_numpy(dtype=float)
            y = subset[right].to_numpy(dtype=float)
            result = spearmanr(x, y)
            ci_low, ci_high, valid = bootstrap_spearman(x, y, rng)
            correlation_rows.append(
                {
                    "feature_1": left,
                    "feature_2": right,
                    "null_hypothesis": "Spearman rho equals zero",
                    "n_genomes": len(subset),
                    "rho": float(result.statistic),
                    "ci_95_low": ci_low,
                    "ci_95_high": ci_high,
                    "bootstrap_replicates": BOOTSTRAPS,
                    "valid_bootstrap_replicates": valid,
                    "raw_p": float(result.pvalue),
                    "effect_label": effect_label(float(result.statistic)),
                }
            )

            x_resid = residualized_rank(x, subset["total_CDS"].to_numpy(dtype=float))
            y_resid = residualized_rank(y, subset["total_CDS"].to_numpy(dtype=float))
            partial = pearsonr(x_resid, y_resid)
            partial_rows.append(
                {
                    "feature_1": left,
                    "feature_2": right,
                    "adjustment": "rank residualization on total_CDS",
                    "null_hypothesis": "partial rank correlation equals zero after total_CDS adjustment",
                    "n_genomes": len(subset),
                    "partial_rho": float(partial.statistic),
                    "raw_p": float(partial.pvalue),
                    "effect_label": effect_label(float(partial.statistic)),
                    "interpretation_scope": "association; residual adjustment does not establish causality",
                }
            )
    for rows in (correlation_rows, partial_rows):
        qvalues = bh_fdr([float(row["raw_p"]) for row in rows])
        for row, qvalue in zip(rows, qvalues):
            row["q_value"] = qvalue
    correlation_rows.sort(key=lambda row: (-abs(float(row["rho"])), str(row["feature_1"]), str(row["feature_2"])))
    partial_rows.sort(key=lambda row: (-abs(float(row["partial_rho"])), str(row["feature_1"]), str(row["feature_2"])))
    atomic_write_tsv(ROOT / "04_statistics/functional_burden_spearman.tsv", correlation_rows, list(correlation_rows[0]))
    atomic_write_tsv(ROOT / "04_statistics/functional_burden_partial_spearman_CDS.tsv", partial_rows, list(partial_rows[0]))

    model_specs = [
        ("ARG_total", ["BGC_total", "total_CDS"]),
        ("AMP_total", ["ARG_total", "total_CDS"]),
        ("BGC_total", ["CAZyme_total", "total_CDS"]),
        ("CAZyme_total", ["BGC_total", "total_CDS"]),
    ]
    regression_rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for outcome, predictors in model_specs:
        subset = data[[outcome, *predictors]].dropna().copy()
        standardized = subset[predictors].apply(lambda column: (column - column.mean()) / column.std(ddof=1))
        design = sm.add_constant(standardized, has_constant="add")
        poisson = sm.GLM(subset[outcome], design, family=sm.families.Poisson()).fit()
        dispersion = float(np.sum(poisson.resid_pearson**2) / poisson.df_resid)
        if dispersion > 1.5:
            fitted = sm.GLM(subset[outcome], design, family=sm.families.Poisson()).fit(scale="X2")
            family = "quasi-Poisson (Pearson-scaled covariance)"
        else:
            fitted = poisson
            family = "Poisson"
        confidence = fitted.conf_int()
        diagnostics.append(
            {
                "model": f"{outcome} ~ {' + '.join(predictors)}",
                "outcome": outcome,
                "n_genomes": len(subset),
                "selected_family": family,
                "Poisson_Pearson_dispersion": dispersion,
                "overdispersion_rule": "quasi-Poisson selected when Pearson dispersion >1.5",
                "AIC_if_Poisson": float(poisson.aic),
                "converged": bool(fitted.converged),
            }
        )
        for term in fitted.params.index:
            coefficient = float(fitted.params[term])
            low, high = map(float, confidence.loc[term])
            regression_rows.append(
                {
                    "model": f"{outcome} ~ {' + '.join(predictors)}",
                    "outcome": outcome,
                    "predictor": term,
                    "null_hypothesis": "log-count coefficient equals zero",
                    "n_genomes": len(subset),
                    "model_family": family,
                    "Poisson_Pearson_dispersion": dispersion,
                    "predictor_scaling": "one standard deviation" if term != "const" else "intercept",
                    "coefficient_log_scale": coefficient,
                    "standard_error": float(fitted.bse[term]),
                    "incidence_rate_ratio": math.exp(coefficient),
                    "IRR_ci_95_low": math.exp(low),
                    "IRR_ci_95_high": math.exp(high),
                    "raw_p": float(fitted.pvalues[term]),
                    "q_value": "",
                    "interpretation_scope": "conditional count association; not causal",
                }
            )
    tested = [index for index, row in enumerate(regression_rows) if row["predictor"] != "const"]
    qvalues = bh_fdr([float(regression_rows[index]["raw_p"]) for index in tested])
    for index, qvalue in zip(tested, qvalues):
        regression_rows[index]["q_value"] = qvalue
    atomic_write_tsv(ROOT / "04_statistics/functional_count_regression.tsv", regression_rows, list(regression_rows[0]))
    atomic_write_tsv(ROOT / "04_statistics/functional_count_regression_diagnostics.tsv", diagnostics, list(diagnostics[0]))

    strongest = correlation_rows[0]
    print(
        f"Burden statistics: {len(correlation_rows)} Spearman tests ({BOOTSTRAPS} bootstraps each); "
        f"strongest={strongest['feature_1']}~{strongest['feature_2']} rho={float(strongest['rho']):.3f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
