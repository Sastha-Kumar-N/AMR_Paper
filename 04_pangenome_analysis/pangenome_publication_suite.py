#!/usr/bin/env python3
"""
Publication-grade pangenome interpretation suite.

Purpose
-------
This script reads Panaroo/Roary-like pangenome outputs and generates validated
summary tables, publication-quality figures, comparison plots, and manuscript-ready
interpretation text for:
    1. whole-species pangenome
    2. clinical pangenome
    3. environmental pangenome
    4. clinical-vs-environmental gene prevalence comparison using a unified matrix

Primary input expected
----------------------
- gene_presence_absence.Rtab from each pangenome run.
- Optional gene_presence_absence.csv from each pangenome run for annotation support.
- Optional metadata TSV for unified clinical/environmental comparison.

Important concept
-----------------
Separate pangenome runs do not always produce directly comparable gene-cluster IDs.
Therefore, aggregate comparisons across whole/clinical/environmental pangenomes are
valid from separate runs, but gene-level clinical-vs-environmental comparison should
be performed using ONE unified whole-pangenome matrix plus metadata classification.

Example
-------
python pangenome_publication_suite.py \
  --whole-rtab /path/to/whole/gene_presence_absence.Rtab \
  --whole-csv /path/to/whole/gene_presence_absence.csv \
  --clinical-rtab /path/to/clinical/gene_presence_absence.Rtab \
  --clinical-csv /path/to/clinical/gene_presence_absence.csv \
  --environmental-rtab /path/to/environmental/gene_presence_absence.Rtab \
  --environmental-csv /path/to/environmental/gene_presence_absence.csv \
  --metadata /mnt/phoenix_2/11.metadata/metadata/final_genome_metadata_classified.tsv \
  --metadata-id-col assembly_accession \
  --metadata-class-col final_classification \
  --outdir /mnt/phoenix_2/16.publication_pangenome_interpretation \
  --curve-reps 50 --curve-step 10

Author: generated for Sastha Kumar's P. aeruginosa pangenome workflow
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import re
import sys
import textwrap
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# matplotlib is imported after backend selection for headless Linux/HPC use.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

try:
    from scipy.optimize import curve_fit
    SCIPY_OPT_AVAILABLE = True
except Exception:
    SCIPY_OPT_AVAILABLE = False

try:
    from scipy.stats import fisher_exact
    SCIPY_STATS_AVAILABLE = True
except Exception:
    SCIPY_STATS_AVAILABLE = False


# -----------------------------
# Global publication aesthetics
# -----------------------------
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


@dataclass
class DatasetResult:
    name: str
    rtab_path: Path
    csv_path: Optional[Path]
    matrix: pd.DataFrame          # rows = gene clusters, columns = genomes, binary 0/1
    annotations: pd.DataFrame     # rows = gene clusters, columns annotation/product/etc if available
    metrics: Dict[str, float]
    per_gene: pd.DataFrame
    per_genome: pd.DataFrame
    curves: Optional[pd.DataFrame]
    validation_messages: List[str]


# -----------------------------
# Generic helpers
# -----------------------------
def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def clean_label(x: str) -> str:
    x = str(x)
    x = x.strip()
    # Common genome labels may contain paths; reduce only if exact matching later fails.
    return x


def stem_variants(label: str) -> List[str]:
    """Return possible comparable forms for genome identifiers."""
    s = str(label).strip()
    base = os.path.basename(s)
    variants = [s, base]
    for ext in [".gff", ".gff3", ".fna", ".fa", ".fasta", ".gz"]:
        if base.endswith(ext):
            variants.append(base[:-len(ext)])
    # Accession-like extraction.
    m = re.search(r"(GC[AF]_[0-9]+\.[0-9]+)", base)
    if m:
        variants.append(m.group(1))
    return list(dict.fromkeys(variants))


def infer_separator(path: Path) -> str:
    if path.suffix.lower() in [".csv"]:
        return ","
    return "\t"


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def savefig(fig: plt.Figure, outbase: Path) -> None:
    """Save a figure as PNG, PDF, and SVG."""
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outbase.with_suffix(f".{ext}"), bbox_inches="tight")
    plt.close(fig)


def bh_fdr(pvalues: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment without external dependency."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    out = np.full(n, np.nan)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * len(ranked) / (np.arange(1, len(ranked) + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    out[np.where(valid)[0]] = q[inv]
    return out


def prevalence_class(freq: float) -> str:
    """Standard pangenome prevalence categories."""
    if freq >= 0.99:
        return "core_99_100"
    if freq >= 0.95:
        return "soft_core_95_99"
    if freq >= 0.15:
        return "shell_15_95"
    return "cloud_0_15"


def category_order() -> List[str]:
    return ["core_99_100", "soft_core_95_99", "shell_15_95", "cloud_0_15"]


def pretty_category(cat: str) -> str:
    return {
        "core_99_100": "Core (≥99%)",
        "soft_core_95_99": "Soft core (95–<99%)",
        "shell_15_95": "Shell (15–<95%)",
        "cloud_0_15": "Cloud (<15%)",
    }.get(cat, cat)


# -----------------------------
# Input parsing and validation
# -----------------------------
def read_rtab(path: Path) -> pd.DataFrame:
    """Read Panaroo/Roary gene_presence_absence.Rtab into binary matrix."""
    if not path.exists():
        raise FileNotFoundError(f"RTAB file not found: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df.shape[1] < 2:
        raise ValueError(f"RTAB has fewer than two columns: {path}")

    first_col = df.columns[0]
    df = df.rename(columns={first_col: "gene_cluster"})
    df["gene_cluster"] = df["gene_cluster"].astype(str)

    # Convert all sample columns to numeric binary values.
    sample_cols = [c for c in df.columns if c != "gene_cluster"]
    mat = df.set_index("gene_cluster")[sample_cols]
    for c in mat.columns:
        mat[c] = pd.to_numeric(mat[c], errors="coerce")

    # Panaroo/Roary Rtab should be 0/1. Missing values become 0 only after warning.
    if mat.isna().any().any():
        n_missing = int(mat.isna().sum().sum())
        warnings.warn(f"{path}: {n_missing} non-numeric or missing matrix values converted to 0")
        mat = mat.fillna(0)
    mat = (mat > 0).astype(np.uint8)
    mat.index.name = "gene_cluster"
    mat.columns = [clean_label(c) for c in mat.columns]
    return mat


def read_annotations(csv_path: Optional[Path], matrix_index: pd.Index) -> pd.DataFrame:
    """Read optional gene_presence_absence.csv for annotation/product columns."""
    if csv_path is None or not csv_path.exists():
        return pd.DataFrame(index=matrix_index)

    df = pd.read_csv(csv_path, low_memory=False)
    if df.empty:
        return pd.DataFrame(index=matrix_index)

    # Panaroo/Roary usually uses Gene as the cluster column.
    possible_id_cols = ["Gene", "gene", "gene_cluster", "Cluster", "cluster"]
    id_col = next((c for c in possible_id_cols if c in df.columns), df.columns[0])
    df[id_col] = df[id_col].astype(str)
    ann = df.set_index(id_col)

    useful_cols = [
        c for c in [
            "Annotation", "annotation", "Product", "product", "Non-unique Gene name",
            "No. isolates", "No. sequences", "Avg sequences per isolate", "Genome Fragment",
            "Order within Fragment", "Accessory Fragment", "QC", "Min group size nuc",
            "Max group size nuc", "Avg group size nuc"
        ] if c in ann.columns
    ]
    if not useful_cols:
        # Keep a small number of non-genome columns if recognizable annotation columns do not exist.
        useful_cols = [c for c in ann.columns[:8]]
    ann = ann[useful_cols].copy()
    ann = ann.reindex(matrix_index)
    return ann


def validate_matrix(name: str, mat: pd.DataFrame, expected_n: Optional[int]) -> List[str]:
    messages = []
    n_genes, n_genomes = mat.shape
    messages.append(f"{name}: loaded {n_genes:,} gene clusters across {n_genomes:,} genomes.")

    if expected_n is not None:
        if n_genomes == expected_n:
            messages.append(f"{name}: genome count matches expected value ({expected_n:,}).")
        else:
            messages.append(
                f"WARNING - {name}: genome count {n_genomes:,} does not match expected value {expected_n:,}."
            )

    if mat.index.duplicated().any():
        dup = int(mat.index.duplicated().sum())
        messages.append(f"WARNING - {name}: {dup:,} duplicated gene cluster IDs detected.")
    else:
        messages.append(f"{name}: gene cluster IDs are unique.")

    if pd.Index(mat.columns).duplicated().any():
        dup = int(pd.Index(mat.columns).duplicated().sum())
        messages.append(f"WARNING - {name}: {dup:,} duplicated genome/sample columns detected.")
    else:
        messages.append(f"{name}: genome/sample columns are unique.")

    values = pd.unique(mat.to_numpy().ravel())
    if set(values).issubset({0, 1}):
        messages.append(f"{name}: matrix is binary after coercion.")
    else:
        messages.append(f"WARNING - {name}: matrix contains non-binary values after coercion: {values[:10]}")

    empty_genes = int((mat.sum(axis=1) == 0).sum())
    empty_genomes = int((mat.sum(axis=0) == 0).sum())
    if empty_genes:
        messages.append(f"WARNING - {name}: {empty_genes:,} gene clusters are absent from all genomes.")
    if empty_genomes:
        messages.append(f"WARNING - {name}: {empty_genomes:,} genomes have zero gene clusters.")
    if not empty_genes and not empty_genomes:
        messages.append(f"{name}: no all-zero gene clusters or all-zero genomes detected.")

    return messages


# -----------------------------
# Per-dataset analysis
# -----------------------------
def compute_per_gene(mat: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    n_genomes = mat.shape[1]
    counts = mat.sum(axis=1).astype(int)
    freq = counts / n_genomes
    out = pd.DataFrame({
        "gene_cluster": mat.index.astype(str),
        "present_genomes": counts.values,
        "total_genomes": n_genomes,
        "prevalence": freq.values,
        "prevalence_percent": (freq.values * 100),
        "pangenome_category": [prevalence_class(x) for x in freq.values],
    })
    if annotations is not None and not annotations.empty:
        ann = annotations.copy()
        ann = ann.reset_index().rename(columns={ann.index.name or "index": "gene_cluster"})
        out = out.merge(ann, on="gene_cluster", how="left")
    return out


def compute_per_genome(mat: pd.DataFrame) -> pd.DataFrame:
    gene_counts = mat.sum(axis=0).astype(int)
    return pd.DataFrame({
        "genome": gene_counts.index.astype(str),
        "gene_clusters_present": gene_counts.values,
    })


def compute_metrics(name: str, mat: pd.DataFrame, per_gene: pd.DataFrame, per_genome: pd.DataFrame) -> Dict[str, float]:
    n_genomes = int(mat.shape[1])
    n_gene_clusters = int(mat.shape[0])
    cat_counts = per_gene["pangenome_category"].value_counts().to_dict()
    metrics = {
        "dataset": name,
        "n_genomes": n_genomes,
        "n_gene_clusters_total_pan": n_gene_clusters,
        "core_99_100": int(cat_counts.get("core_99_100", 0)),
        "soft_core_95_99": int(cat_counts.get("soft_core_95_99", 0)),
        "shell_15_95": int(cat_counts.get("shell_15_95", 0)),
        "cloud_0_15": int(cat_counts.get("cloud_0_15", 0)),
        "core_plus_soft_core_95_100": int(cat_counts.get("core_99_100", 0) + cat_counts.get("soft_core_95_99", 0)),
        "accessory_0_95": int(cat_counts.get("shell_15_95", 0) + cat_counts.get("cloud_0_15", 0)),
        "mean_gene_clusters_per_genome": float(per_genome["gene_clusters_present"].mean()),
        "median_gene_clusters_per_genome": float(per_genome["gene_clusters_present"].median()),
        "min_gene_clusters_per_genome": int(per_genome["gene_clusters_present"].min()),
        "max_gene_clusters_per_genome": int(per_genome["gene_clusters_present"].max()),
        "sd_gene_clusters_per_genome": float(per_genome["gene_clusters_present"].std(ddof=1)) if n_genomes > 1 else 0.0,
    }
    for cat in category_order():
        metrics[f"{cat}_percent_of_pan"] = 100.0 * metrics[cat] / max(n_gene_clusters, 1)
    return metrics


def pan_model(n: np.ndarray, k: float, gamma: float) -> np.ndarray:
    return k * np.power(n, gamma)


def core_model(n: np.ndarray, c: float, k: float, lam: float) -> np.ndarray:
    return c + k * np.exp(-lam * n)


def compute_rarefaction_curves(
    mat: pd.DataFrame,
    reps: int = 50,
    step: int = 10,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Randomized pan/core accumulation curves.

    Uses random genome order permutations. Values are recorded at interval step
    and at final N. The calculation is intentionally matrix-based and avoids
    repeatedly re-reading files.
    """
    rng = np.random.default_rng(random_seed)
    arr = mat.to_numpy(dtype=bool)  # genes x genomes
    n_genes, n_genomes = arr.shape
    if n_genomes == 0 or n_genes == 0:
        return pd.DataFrame()

    points = list(range(1, n_genomes + 1, max(1, step)))
    if points[-1] != n_genomes:
        points.append(n_genomes)
    points_set = set(points)

    rows = []
    for rep in range(1, reps + 1):
        perm = rng.permutation(n_genomes)
        pan_seen = np.zeros(n_genes, dtype=bool)
        core_seen = np.ones(n_genes, dtype=bool)
        for i, genome_idx in enumerate(perm, start=1):
            col = arr[:, genome_idx]
            pan_seen |= col
            core_seen &= col
            if i in points_set:
                rows.append({
                    "replicate": rep,
                    "n_genomes": i,
                    "pan_gene_clusters": int(pan_seen.sum()),
                    "core_gene_clusters": int(core_seen.sum()),
                })
    curve = pd.DataFrame(rows)
    if curve.empty:
        return curve

    grouped = curve.groupby("n_genomes", as_index=False).agg(
        pan_mean=("pan_gene_clusters", "mean"),
        pan_sd=("pan_gene_clusters", "std"),
        pan_q025=("pan_gene_clusters", lambda x: np.quantile(x, 0.025)),
        pan_q975=("pan_gene_clusters", lambda x: np.quantile(x, 0.975)),
        core_mean=("core_gene_clusters", "mean"),
        core_sd=("core_gene_clusters", "std"),
        core_q025=("core_gene_clusters", lambda x: np.quantile(x, 0.025)),
        core_q975=("core_gene_clusters", lambda x: np.quantile(x, 0.975)),
    )
    grouped["pan_sd"] = grouped["pan_sd"].fillna(0)
    grouped["core_sd"] = grouped["core_sd"].fillna(0)

    # Optional model fits for interpretability.
    grouped["pan_heaps_k"] = np.nan
    grouped["pan_heaps_gamma"] = np.nan
    grouped["core_decay_asymptote"] = np.nan
    grouped["core_decay_lambda"] = np.nan
    if SCIPY_OPT_AVAILABLE and len(grouped) >= 4:
        n = grouped["n_genomes"].to_numpy(dtype=float)
        pan_y = grouped["pan_mean"].to_numpy(dtype=float)
        core_y = grouped["core_mean"].to_numpy(dtype=float)
        try:
            p0 = [max(pan_y[0], 1.0), 0.2]
            popt, _ = curve_fit(pan_model, n, pan_y, p0=p0, maxfev=10000, bounds=([0, -5], [np.inf, 5]))
            grouped["pan_heaps_k"] = popt[0]
            grouped["pan_heaps_gamma"] = popt[1]
        except Exception as exc:
            warnings.warn(f"Pan-genome Heaps fit failed: {exc}")
        try:
            p0 = [max(core_y[-1], 1.0), max(core_y[0] - core_y[-1], 1.0), 0.01]
            popt, _ = curve_fit(core_model, n, core_y, p0=p0, maxfev=10000, bounds=([0, 0, 0], [np.inf, np.inf, 10]))
            grouped["core_decay_asymptote"] = popt[0]
            grouped["core_decay_lambda"] = popt[2]
        except Exception as exc:
            warnings.warn(f"Core-genome decay fit failed: {exc}")
    return grouped


def analyze_dataset(
    name: str,
    rtab_path: Path,
    csv_path: Optional[Path],
    expected_n: Optional[int],
    outdir: Path,
    curve_reps: int,
    curve_step: int,
    random_seed: int,
    skip_curves: bool,
) -> DatasetResult:
    mat = read_rtab(rtab_path)
    annotations = read_annotations(csv_path, mat.index)
    validation = validate_matrix(name, mat, expected_n)
    per_gene = compute_per_gene(mat, annotations)
    per_genome = compute_per_genome(mat)
    metrics = compute_metrics(name, mat, per_gene, per_genome)
    curves = None if skip_curves else compute_rarefaction_curves(mat, curve_reps, curve_step, random_seed)

    tables_dir = outdir / "tables" / name
    mkdir(tables_dir)
    write_tsv(per_gene, tables_dir / f"{name}_gene_prevalence_categories.tsv")
    write_tsv(per_genome, tables_dir / f"{name}_per_genome_gene_cluster_counts.tsv")
    pd.DataFrame([metrics]).to_csv(tables_dir / f"{name}_summary_metrics.tsv", sep="\t", index=False)
    if curves is not None and not curves.empty:
        curves.to_csv(tables_dir / f"{name}_rarefaction_curves.tsv", sep="\t", index=False)

    return DatasetResult(
        name=name,
        rtab_path=rtab_path,
        csv_path=csv_path,
        matrix=mat,
        annotations=annotations,
        metrics=metrics,
        per_gene=per_gene,
        per_genome=per_genome,
        curves=curves,
        validation_messages=validation,
    )


# -----------------------------
# Figures: per dataset
# -----------------------------
def plot_category_bar(ds: DatasetResult, figdir: Path) -> None:
    cats = category_order()
    counts = [ds.metrics[c] for c in cats]
    labels = [pretty_category(c) for c in cats]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    bars = ax.bar(labels, counts)
    ax.set_ylabel("Gene clusters")
    ax.set_title(f"{ds.name}: pangenome frequency classes")
    ax.tick_params(axis="x", rotation=30)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    for b, v in zip(bars, counts):
        ax.text(b.get_x() + b.get_width()/2, b.get_height(), f"{v:,}", ha="center", va="bottom", fontsize=8)
    savefig(fig, figdir / f"{ds.name}_frequency_class_bar")


def plot_prevalence_histogram(ds: DatasetResult, figdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.hist(ds.per_gene["prevalence_percent"], bins=np.linspace(0, 100, 51), edgecolor="black", linewidth=0.25)
    ax.axvline(15, linestyle="--", linewidth=1)
    ax.axvline(95, linestyle="--", linewidth=1)
    ax.axvline(99, linestyle="--", linewidth=1)
    ax.set_xlabel("Prevalence across genomes (%)")
    ax.set_ylabel("Gene clusters")
    ax.set_title(f"{ds.name}: gene-cluster prevalence distribution")
    savefig(fig, figdir / f"{ds.name}_gene_prevalence_histogram")


def plot_per_genome_distribution(ds: DatasetResult, figdir: Path) -> None:
    data = ds.per_genome["gene_clusters_present"].to_numpy()
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    # Matplotlib 3.11+ deprecated vert=; fallback supports older versions.
    try:
        ax.boxplot(data, orientation="vertical", showfliers=True)
    except TypeError:
        ax.boxplot(data, vert=True, showfliers=True)
    jitter = np.random.default_rng(42).normal(1, 0.025, size=len(data))
    # If many genomes, subsample points for readability.
    if len(data) > 500:
        idx = np.random.default_rng(42).choice(len(data), size=500, replace=False)
        ax.scatter(jitter[idx], data[idx], s=6, alpha=0.35)
    else:
        ax.scatter(jitter, data, s=8, alpha=0.45)
    ax.set_xticks([1])
    ax.set_xticklabels([ds.name])
    ax.set_ylabel("Gene clusters per genome")
    ax.set_title(f"{ds.name}: per-genome gene-cluster burden")
    savefig(fig, figdir / f"{ds.name}_per_genome_gene_cluster_boxplot")


def plot_rarefaction(ds: DatasetResult, figdir: Path) -> None:
    if ds.curves is None or ds.curves.empty:
        return
    curve = ds.curves
    x = curve["n_genomes"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(x, curve["pan_mean"], label="Pan-genome")
    ax.fill_between(x, curve["pan_q025"], curve["pan_q975"], alpha=0.18)
    ax.plot(x, curve["core_mean"], label="Core genome")
    ax.fill_between(x, curve["core_q025"], curve["core_q975"], alpha=0.18)
    ax.set_xlabel("Number of genomes sampled")
    ax.set_ylabel("Gene clusters")
    ax.set_title(f"{ds.name}: randomized pan/core accumulation")
    ax.legend(frameon=False)
    savefig(fig, figdir / f"{ds.name}_pan_core_rarefaction")


def plot_presence_absence_heatmap(ds: DatasetResult, figdir: Path, max_genes: int = 600, max_genomes: int = 300) -> None:
    """Lightweight heatmap of the most variable gene clusters.

    This is not a full phylogenetic heatmap. It is designed as a compact QC/publication
    panel showing accessory-gene structure without requiring SciPy clustering.
    """
    mat = ds.matrix
    if mat.empty:
        return
    freq = mat.mean(axis=1)
    # Select most variable gene clusters, excluding invariant core/absent clusters.
    variance = freq * (1 - freq)
    selected_genes = variance.sort_values(ascending=False).head(min(max_genes, len(variance))).index
    sub = mat.loc[selected_genes]
    # Sort genomes by total selected accessory burden.
    genome_order = sub.sum(axis=0).sort_values(ascending=False).index
    if len(genome_order) > max_genomes:
        genome_order = genome_order[:max_genomes]
    sub = sub[genome_order]
    fig_w = max(6.2, min(12, sub.shape[1] / 40))
    fig_h = max(4.5, min(10, sub.shape[0] / 80))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(sub.to_numpy(), aspect="auto", interpolation="nearest")
    ax.set_xlabel(f"Genomes (top {sub.shape[1]} shown by accessory burden)")
    ax.set_ylabel(f"Variable gene clusters (top {sub.shape[0]})")
    ax.set_title(f"{ds.name}: variable gene presence/absence structure")
    ax.set_xticks([])
    ax.set_yticks([])
    savefig(fig, figdir / f"{ds.name}_variable_gene_presence_absence_heatmap")


def make_per_dataset_figures(ds: DatasetResult, outdir: Path) -> None:
    figdir = outdir / "figures" / ds.name
    mkdir(figdir)
    plot_category_bar(ds, figdir)
    plot_prevalence_histogram(ds, figdir)
    plot_per_genome_distribution(ds, figdir)
    plot_rarefaction(ds, figdir)
    plot_presence_absence_heatmap(ds, figdir)


# -----------------------------
# Cross-dataset aggregate comparison
# -----------------------------
def make_comparison_tables(results: Dict[str, DatasetResult], outdir: Path) -> pd.DataFrame:
    rows = []
    for name, ds in results.items():
        rows.append(ds.metrics)
    comp = pd.DataFrame(rows)
    mkdir(outdir / "tables")
    comp.to_csv(outdir / "tables" / "pangenome_three_group_summary_metrics.tsv", sep="\t", index=False)

    long_rows = []
    for _, row in comp.iterrows():
        for cat in category_order():
            long_rows.append({
                "dataset": row["dataset"],
                "category": cat,
                "category_label": pretty_category(cat),
                "gene_clusters": row[cat],
                "percent_of_pan": row[f"{cat}_percent_of_pan"],
                "n_genomes": row["n_genomes"],
                "total_pan_gene_clusters": row["n_gene_clusters_total_pan"],
            })
    comp_long = pd.DataFrame(long_rows)
    comp_long.to_csv(outdir / "tables" / "pangenome_frequency_class_comparison_long.tsv", sep="\t", index=False)
    return comp


def plot_comparison_stacked_percent(results: Dict[str, DatasetResult], outdir: Path) -> None:
    names = list(results.keys())
    cats = category_order()
    data = np.array([[results[n].metrics[f"{c}_percent_of_pan"] for c in cats] for n in names])
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    bottom = np.zeros(len(names))
    for j, cat in enumerate(cats):
        ax.bar(names, data[:, j], bottom=bottom, label=pretty_category(cat))
        bottom += data[:, j]
    ax.set_ylabel("Percentage of pangenome gene clusters")
    ax.set_ylim(0, 100)
    ax.set_title("Pangenome composition across whole, clinical, and environmental datasets")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, outdir / "figures" / "comparison_frequency_classes_stacked_percent")


def plot_comparison_absolute_categories(results: Dict[str, DatasetResult], outdir: Path) -> None:
    names = list(results.keys())
    cats = category_order()
    x = np.arange(len(names))
    width = 0.18
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    for j, cat in enumerate(cats):
        vals = [results[n].metrics[cat] for n in names]
        ax.bar(x + (j - 1.5) * width, vals, width=width, label=pretty_category(cat))
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Gene clusters")
    ax.set_title("Absolute pangenome frequency classes")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, outdir / "figures" / "comparison_frequency_classes_absolute")


def plot_comparison_per_genome_boxplot(results: Dict[str, DatasetResult], outdir: Path) -> None:
    names = list(results.keys())
    data = [results[n].per_genome["gene_clusters_present"].to_numpy() for n in names]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    # Matplotlib 3.11+ renamed labels= to tick_labels=; fallback supports older versions.
    try:
        ax.boxplot(data, tick_labels=names, showfliers=True)
    except TypeError:
        ax.boxplot(data, labels=names, showfliers=True)
    ax.set_ylabel("Gene clusters per genome")
    ax.set_title("Per-genome gene-cluster burden across datasets")
    savefig(fig, outdir / "figures" / "comparison_per_genome_gene_cluster_boxplot")


def plot_comparison_rarefaction_overlay(results: Dict[str, DatasetResult], outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    plotted = False
    for name, ds in results.items():
        if ds.curves is None or ds.curves.empty:
            continue
        ax.plot(ds.curves["n_genomes"], ds.curves["pan_mean"], label=f"{name} pan")
        ax.plot(ds.curves["n_genomes"], ds.curves["core_mean"], linestyle="--", label=f"{name} core")
        plotted = True
    if not plotted:
        plt.close(fig)
        return
    ax.set_xlabel("Number of genomes sampled")
    ax.set_ylabel("Gene clusters")
    ax.set_title("Randomized pan/core accumulation across datasets")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    savefig(fig, outdir / "figures" / "comparison_pan_core_rarefaction_overlay")


def plot_pan_size_vs_genomes(results: Dict[str, DatasetResult], outdir: Path) -> None:
    names = list(results.keys())
    x = [results[n].metrics["n_genomes"] for n in names]
    y = [results[n].metrics["n_gene_clusters_total_pan"] for n in names]
    labels = names
    fig, ax = plt.subplots(figsize=(5.8, 4.5))
    ax.scatter(x, y, s=55)
    for xi, yi, lab in zip(x, y, labels):
        ax.text(xi, yi, f" {lab}", va="center", fontsize=9)
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("Total pangenome gene clusters")
    ax.set_title("Pangenome size relative to sample size")
    savefig(fig, outdir / "figures" / "comparison_pan_size_vs_genome_count")


def compare_gene_cluster_id_overlap(results: Dict[str, DatasetResult], outdir: Path) -> None:
    """Only meaningful if pangenome IDs are harmonized. Output warns user."""
    names = list(results.keys())
    sets = {n: set(results[n].matrix.index.astype(str)) for n in names}
    rows = []
    for a, b in itertools.product(names, names):
        inter = len(sets[a] & sets[b])
        union = len(sets[a] | sets[b])
        rows.append({
            "dataset_a": a,
            "dataset_b": b,
            "shared_gene_cluster_ids": inter,
            "union_gene_cluster_ids": union,
            "jaccard_id_overlap": inter / union if union else np.nan,
            "warning": "Interpret only if gene-cluster IDs are harmonized across runs. Separate Panaroo/Roary runs may not be ID-comparable.",
        })
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "tables" / "gene_cluster_id_overlap_warning_table.tsv", sep="\t", index=False)

    pivot = df.pivot(index="dataset_a", columns="dataset_b", values="jaccard_id_overlap").loc[names, names]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(pivot.to_numpy(), vmin=0, vmax=1)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{pivot.iloc[i, j]:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title("Gene-cluster ID overlap\n(use only for harmonized IDs)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Jaccard overlap")
    savefig(fig, outdir / "figures" / "comparison_gene_cluster_id_overlap_warning_heatmap")


def make_comparison_figures(results: Dict[str, DatasetResult], outdir: Path) -> None:
    mkdir(outdir / "figures")
    plot_comparison_stacked_percent(results, outdir)
    plot_comparison_absolute_categories(results, outdir)
    plot_comparison_per_genome_boxplot(results, outdir)
    plot_comparison_rarefaction_overlay(results, outdir)
    plot_pan_size_vs_genomes(results, outdir)
    compare_gene_cluster_id_overlap(results, outdir)


# -----------------------------
# Unified clinical-vs-environmental analysis
# -----------------------------
def normalize_class_label(x: str) -> str:
    s = str(x).strip().lower().replace("-", "_").replace(" ", "_")
    if "clinical" in s or "human" in s or "infection" in s:
        return "clinical"
    if "environment" in s or "soil" in s or "water" in s or "plant" in s or "industrial" in s or "natural" in s:
        return "environmental"
    return "unknown_other"


def match_metadata_to_matrix_columns(
    matrix_cols: Sequence[str],
    metadata: pd.DataFrame,
    metadata_id_col: str,
) -> pd.DataFrame:
    """Map metadata rows to pangenome matrix columns using exact, basename, stem, and accession variants."""
    if metadata_id_col not in metadata.columns:
        raise ValueError(f"metadata-id-col '{metadata_id_col}' not found in metadata columns: {list(metadata.columns)}")

    lookup = {}
    for _, row in metadata.iterrows():
        raw = str(row[metadata_id_col]).strip()
        for v in stem_variants(raw):
            lookup.setdefault(v, row)

    mapped_rows = []
    unmatched = []
    for col in matrix_cols:
        match = None
        for v in stem_variants(col):
            if v in lookup:
                match = lookup[v]
                break
        if match is None:
            unmatched.append(col)
        else:
            d = match.to_dict()
            d["matrix_column"] = col
            mapped_rows.append(d)
    mapped = pd.DataFrame(mapped_rows)
    mapped.attrs["unmatched_columns"] = unmatched
    return mapped


def run_unified_clinical_environmental_analysis(
    unified_rtab: Path,
    unified_csv: Optional[Path],
    metadata_path: Path,
    metadata_id_col: str,
    metadata_class_col: str,
    outdir: Path,
    min_abs_delta: float,
    fdr_cutoff: float,
) -> Optional[pd.DataFrame]:
    if metadata_path is None or not metadata_path.exists():
        return None
    if unified_rtab is None or not unified_rtab.exists():
        return None

    mat = read_rtab(unified_rtab)
    annotations = read_annotations(unified_csv, mat.index)
    metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    if metadata_class_col not in metadata.columns:
        # Try common alternatives.
        candidates = [c for c in metadata.columns if c.lower() in [
            "final_classification", "classification", "source_class", "class", "isolate_class"
        ]]
        if candidates:
            metadata_class_col = candidates[0]
        else:
            raise ValueError(f"metadata-class-col not found: {metadata_class_col}")

    mapped = match_metadata_to_matrix_columns(mat.columns, metadata, metadata_id_col)
    unmatched = mapped.attrs.get("unmatched_columns", [])
    validation_lines = []
    validation_lines.append(f"Unified matrix: {mat.shape[0]:,} gene clusters x {mat.shape[1]:,} genome columns.")
    validation_lines.append(f"Metadata rows: {metadata.shape[0]:,}.")
    validation_lines.append(f"Matched matrix columns to metadata: {mapped.shape[0]:,}/{mat.shape[1]:,}.")
    if unmatched:
        validation_lines.append(f"WARNING - Unmatched matrix columns: {len(unmatched):,}. First 20: {unmatched[:20]}")

    mapped["normalized_class"] = mapped[metadata_class_col].map(normalize_class_label)
    class_counts = mapped["normalized_class"].value_counts().to_dict()
    validation_lines.append(f"Matched class counts: {class_counts}")

    clinical_cols = mapped.loc[mapped["normalized_class"] == "clinical", "matrix_column"].tolist()
    env_cols = mapped.loc[mapped["normalized_class"] == "environmental", "matrix_column"].tolist()
    if len(clinical_cols) == 0 or len(env_cols) == 0:
        validation_lines.append("ERROR - clinical or environmental group has zero matched genomes; differential analysis skipped.")
        mkdir(outdir / "validation")
        (outdir / "validation" / "unified_clinical_environmental_validation.md").write_text("\n".join(validation_lines))
        return None

    clinical_present = mat[clinical_cols].sum(axis=1).astype(int)
    env_present = mat[env_cols].sum(axis=1).astype(int)
    n_clin = len(clinical_cols)
    n_env = len(env_cols)
    clinical_freq = clinical_present / n_clin
    env_freq = env_present / n_env
    delta = clinical_freq - env_freq

    if not SCIPY_STATS_AVAILABLE:
        raise ImportError("scipy is required for Fisher exact testing. Install with: pip install scipy")

    pvals = []
    odds = []
    for cp, ep in zip(clinical_present.values, env_present.values):
        table = [[int(cp), int(n_clin - cp)], [int(ep), int(n_env - ep)]]
        try:
            od, p = fisher_exact(table, alternative="two-sided")
        except Exception:
            od, p = np.nan, np.nan
        odds.append(od)
        pvals.append(p)
    fdr = bh_fdr(pvals)

    diff = pd.DataFrame({
        "gene_cluster": mat.index.astype(str),
        "clinical_present": clinical_present.values,
        "clinical_total": n_clin,
        "clinical_prevalence": clinical_freq.values,
        "clinical_prevalence_percent": clinical_freq.values * 100,
        "environmental_present": env_present.values,
        "environmental_total": n_env,
        "environmental_prevalence": env_freq.values,
        "environmental_prevalence_percent": env_freq.values * 100,
        "clinical_minus_environmental_prevalence": delta.values,
        "clinical_minus_environmental_prevalence_percent": delta.values * 100,
        "odds_ratio": odds,
        "fisher_pvalue": pvals,
        "fdr_bh": fdr,
        "direction": np.where(delta.values > 0, "clinical_higher", np.where(delta.values < 0, "environmental_higher", "equal")),
    })
    diff["passes_fdr"] = diff["fdr_bh"] <= fdr_cutoff
    diff["passes_abs_delta"] = diff["clinical_minus_environmental_prevalence"].abs() >= min_abs_delta
    diff["candidate"] = diff["passes_fdr"] & diff["passes_abs_delta"]

    if annotations is not None and not annotations.empty:
        ann = annotations.reset_index().rename(columns={annotations.index.name or "index": "gene_cluster"})
        diff = diff.merge(ann, on="gene_cluster", how="left")

    diff = diff.sort_values(["candidate", "fdr_bh", "clinical_minus_environmental_prevalence"], ascending=[False, True, False])

    tables_dir = outdir / "tables" / "unified_clinical_environmental"
    figdir = outdir / "figures" / "unified_clinical_environmental"
    valdir = outdir / "validation"
    mkdir(tables_dir); mkdir(figdir); mkdir(valdir)
    diff.to_csv(tables_dir / "unified_clinical_environmental_differential_prevalence_all_genes.tsv", sep="\t", index=False)
    diff.loc[diff["candidate"]].to_csv(tables_dir / "unified_clinical_environmental_candidate_genes_fdr_delta.tsv", sep="\t", index=False)
    mapped.to_csv(tables_dir / "unified_matrix_metadata_matched_columns.tsv", sep="\t", index=False)

    validation_lines.append(f"Clinical genomes used: {n_clin:,}.")
    validation_lines.append(f"Environmental genomes used: {n_env:,}.")
    validation_lines.append(f"FDR cutoff: {fdr_cutoff}; absolute prevalence delta cutoff: {min_abs_delta}.")
    validation_lines.append(f"Candidate genes passing both cutoffs: {int(diff['candidate'].sum()):,}.")
    (valdir / "unified_clinical_environmental_validation.md").write_text("\n".join(validation_lines) + "\n")

    plot_unified_prevalence_scatter(diff, figdir, fdr_cutoff, min_abs_delta)
    plot_unified_volcano(diff, figdir, fdr_cutoff, min_abs_delta)
    plot_unified_candidate_bar(diff, figdir)
    plot_unified_candidate_heatmap(mat, diff, clinical_cols, env_cols, figdir)
    return diff


def plot_unified_prevalence_scatter(diff: pd.DataFrame, figdir: Path, fdr_cutoff: float, min_abs_delta: float) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    non = ~diff["candidate"]
    cand = diff["candidate"]
    ax.scatter(diff.loc[non, "environmental_prevalence_percent"], diff.loc[non, "clinical_prevalence_percent"], s=8, alpha=0.35)
    ax.scatter(diff.loc[cand, "environmental_prevalence_percent"], diff.loc[cand, "clinical_prevalence_percent"], s=18, alpha=0.9)
    ax.plot([0, 100], [0, 100], linestyle="--", linewidth=1)
    ax.set_xlabel("Environmental prevalence (%)")
    ax.set_ylabel("Clinical prevalence (%)")
    ax.set_title("Unified pangenome: clinical vs environmental prevalence")
    txt = f"Candidate: FDR ≤ {fdr_cutoff}, |Δ| ≥ {min_abs_delta*100:.0f}%"
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, ha="left", va="top", fontsize=8)
    savefig(fig, figdir / "unified_clinical_environmental_prevalence_scatter")


def plot_unified_volcano(diff: pd.DataFrame, figdir: Path, fdr_cutoff: float, min_abs_delta: float) -> None:
    plot_df = diff.copy()
    safe_fdr = plot_df["fdr_bh"].replace(0, np.nextafter(0, 1))
    plot_df["minus_log10_fdr"] = -np.log10(safe_fdr)
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    non = ~plot_df["candidate"]
    cand = plot_df["candidate"]
    ax.scatter(plot_df.loc[non, "clinical_minus_environmental_prevalence_percent"], plot_df.loc[non, "minus_log10_fdr"], s=8, alpha=0.35)
    ax.scatter(plot_df.loc[cand, "clinical_minus_environmental_prevalence_percent"], plot_df.loc[cand, "minus_log10_fdr"], s=20, alpha=0.9)
    ax.axvline(min_abs_delta * 100, linestyle="--", linewidth=1)
    ax.axvline(-min_abs_delta * 100, linestyle="--", linewidth=1)
    ax.axhline(-math.log10(fdr_cutoff), linestyle="--", linewidth=1)
    ax.set_xlabel("Clinical − environmental prevalence difference (%)")
    ax.set_ylabel("−log10(FDR)")
    ax.set_title("Differential gene prevalence volcano plot")
    savefig(fig, figdir / "unified_clinical_environmental_volcano")


def best_annotation_column(df: pd.DataFrame) -> Optional[str]:
    for c in ["Annotation", "annotation", "Product", "product", "Non-unique Gene name"]:
        if c in df.columns:
            return c
    return None


def plot_unified_candidate_bar(diff: pd.DataFrame, figdir: Path, top_n: int = 30) -> None:
    cand = diff.loc[diff["candidate"]].copy()
    if cand.empty:
        return
    cand["abs_delta"] = cand["clinical_minus_environmental_prevalence_percent"].abs()
    cand = cand.sort_values("abs_delta", ascending=False).head(top_n)
    ann_col = best_annotation_column(cand)
    labels = cand["gene_cluster"].astype(str)
    if ann_col:
        labels = cand[ann_col].fillna(cand["gene_cluster"]).astype(str)
        labels = labels.str.slice(0, 80)
    y = np.arange(len(cand))
    fig_h = max(4.5, 0.28 * len(cand) + 1.5)
    fig, ax = plt.subplots(figsize=(7.0, fig_h))
    ax.barh(y, cand["clinical_minus_environmental_prevalence_percent"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.axvline(0, linewidth=1)
    ax.set_xlabel("Clinical − environmental prevalence difference (%)")
    ax.set_title("Top differential clinical/environmental candidate genes")
    ax.invert_yaxis()
    savefig(fig, figdir / "unified_clinical_environmental_candidate_delta_bar")


def plot_unified_candidate_heatmap(
    mat: pd.DataFrame,
    diff: pd.DataFrame,
    clinical_cols: List[str],
    env_cols: List[str],
    figdir: Path,
    top_n: int = 80,
    max_per_group: int = 150,
) -> None:
    cand = diff.loc[diff["candidate"]].copy()
    if cand.empty:
        return
    cand["abs_delta"] = cand["clinical_minus_environmental_prevalence"].abs()
    genes = cand.sort_values("abs_delta", ascending=False).head(top_n)["gene_cluster"].tolist()
    rng = np.random.default_rng(42)
    clin = clinical_cols if len(clinical_cols) <= max_per_group else list(rng.choice(clinical_cols, max_per_group, replace=False))
    env = env_cols if len(env_cols) <= max_per_group else list(rng.choice(env_cols, max_per_group, replace=False))
    cols = clin + env
    sub = mat.loc[genes, cols]

    fig_w = max(7, min(12, len(cols) / 25))
    fig_h = max(4.5, min(10, len(genes) / 8))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(sub.to_numpy(), aspect="auto", interpolation="nearest")
    ax.axvline(len(clin) - 0.5, linewidth=1.2)
    ax.set_xticks([len(clin)/2, len(clin) + len(env)/2])
    ax.set_xticklabels([f"Clinical (n={len(clin)})", f"Environmental (n={len(env)})"])
    ax.set_yticks([])
    ax.set_ylabel(f"Candidate genes (top {len(genes)})")
    ax.set_title("Presence/absence pattern of differential candidate genes")
    savefig(fig, figdir / "unified_clinical_environmental_candidate_presence_absence_heatmap")


# -----------------------------
# Validation and manuscript text
# -----------------------------
def write_validation_report(results: Dict[str, DatasetResult], outdir: Path) -> None:
    mkdir(outdir / "validation")
    lines = ["# Pangenome input validation report", ""]
    for name, ds in results.items():
        lines.append(f"## {name}")
        lines.extend([f"- {m}" for m in ds.validation_messages])
        lines.append("")
    lines.append("## Cross-run comparison caution")
    lines.append("- Aggregate quantities such as pan-genome size, core size, shell/cloud fraction, and per-genome gene-cluster counts can be compared across independently generated pangenomes.")
    lines.append("- Gene-level clinical-vs-environmental interpretation should use the unified whole-pangenome matrix plus metadata classification, because independent pangenome runs may assign non-identical cluster IDs.")
    (outdir / "validation" / "pangenome_validation_report.md").write_text("\n".join(lines) + "\n")


def interpretation_sentence(ds: DatasetResult) -> str:
    m = ds.metrics
    gamma_txt = "not estimated"
    asym_txt = "not estimated"
    if ds.curves is not None and not ds.curves.empty:
        gamma = ds.curves["pan_heaps_gamma"].dropna()
        asym = ds.curves["core_decay_asymptote"].dropna()
        if not gamma.empty:
            gamma_txt = f"{gamma.iloc[0]:.3f}"
        if not asym.empty:
            asym_txt = f"{asym.iloc[0]:.0f} gene clusters"
    return (
        f"The {ds.name} pangenome comprised {m['n_gene_clusters_total_pan']:,} gene clusters across "
        f"{m['n_genomes']:,} genomes. The strict core genome represented {m['core_99_100']:,} clusters "
        f"({m['core_99_100_percent_of_pan']:.2f}% of the pangenome), while the soft-core, shell, and cloud "
        f"fractions contained {m['soft_core_95_99']:,}, {m['shell_15_95']:,}, and {m['cloud_0_15']:,} clusters, respectively. "
        f"Individual genomes encoded a median of {m['median_gene_clusters_per_genome']:.0f} gene clusters "
        f"(range {m['min_gene_clusters_per_genome']:.0f}–{m['max_gene_clusters_per_genome']:.0f}). "
        f"Randomized accumulation analysis estimated a Heaps-law pan-genome exponent of {gamma_txt} and a fitted "
        f"core-genome asymptote of {asym_txt}."
    )


def write_manuscript_interpretation(results: Dict[str, DatasetResult], diff: Optional[pd.DataFrame], outdir: Path) -> None:
    lines = []
    lines.append("# Manuscript-ready pangenome result interpretation")
    lines.append("")
    lines.append("## Individual pangenome structure")
    for name, ds in results.items():
        lines.append(f"### {name.capitalize()} pangenome")
        lines.append(interpretation_sentence(ds))
        lines.append("")

    lines.append("## Comparative interpretation")
    if len(results) >= 2:
        comp = pd.DataFrame([ds.metrics for ds in results.values()])
        largest_pan = comp.sort_values("n_gene_clusters_total_pan", ascending=False).iloc[0]
        highest_cloud = comp.sort_values("cloud_0_15_percent_of_pan", ascending=False).iloc[0]
        highest_core_frac = comp.sort_values("core_99_100_percent_of_pan", ascending=False).iloc[0]
        lines.append(
            f"Across the analyzed pangenomes, {largest_pan['dataset']} had the largest total pangenome size "
            f"({int(largest_pan['n_gene_clusters_total_pan']):,} clusters), whereas {highest_cloud['dataset']} showed the "
            f"highest proportional cloud-genome component ({highest_cloud['cloud_0_15_percent_of_pan']:.2f}%). "
            f"The greatest proportional strict-core fraction was observed in {highest_core_frac['dataset']} "
            f"({highest_core_frac['core_99_100_percent_of_pan']:.2f}%). These contrasts should be interpreted together "
            f"with genome sample size and source diversity, because pan-genome breadth is sensitive to both sampling depth and ecological heterogeneity."
        )
    lines.append("")

    lines.append("## Clinical versus environmental gene prevalence from unified matrix")
    if diff is not None and not diff.empty:
        cand = diff.loc[diff["candidate"]].copy()
        n_cand = cand.shape[0]
        clin_higher = int((cand["direction"] == "clinical_higher").sum())
        env_higher = int((cand["direction"] == "environmental_higher").sum())
        lines.append(
            f"Using the unified pangenome matrix, {n_cand:,} gene clusters passed the combined differential-prevalence "
            f"criteria. Of these, {clin_higher:,} were enriched in clinical genomes and {env_higher:,} were enriched "
            f"in environmental genomes. Candidate loci should be prioritized by effect size, FDR, annotation quality, "
            f"and biological relevance to antimicrobial resistance, mobile genetic elements, virulence, persistence, "
            f"surface exposure, conservation, and vaccine/drug targetability."
        )
        if n_cand:
            ann_col = best_annotation_column(cand)
            top = cand.sort_values("fdr_bh").head(10)
            lines.append("")
            lines.append("Top candidate loci by FDR:")
            for _, r in top.iterrows():
                label = r["gene_cluster"]
                if ann_col and pd.notna(r.get(ann_col, np.nan)):
                    label = f"{label} ({r[ann_col]})"
                lines.append(
                    f"- {label}: clinical {r['clinical_present']:.0f}/{r['clinical_total']:.0f} "
                    f"({r['clinical_prevalence_percent']:.1f}%), environmental {r['environmental_present']:.0f}/"
                    f"{r['environmental_total']:.0f} ({r['environmental_prevalence_percent']:.1f}%), "
                    f"Δ={r['clinical_minus_environmental_prevalence_percent']:.1f}%, FDR={r['fdr_bh']:.3g}."
                )
    else:
        lines.append(
            "Unified clinical/environmental differential-prevalence analysis was not run. Provide --metadata and use the "
            "whole unified gene_presence_absence.Rtab as --unified-rtab or --whole-rtab."
        )
    lines.append("")

    lines.append("## Suggested figure set for a high-impact manuscript")
    lines.append("- Figure 1: Pangenome composition across whole, clinical, and environmental datasets as stacked proportional bars.")
    lines.append("- Figure 2: Randomized pan/core accumulation curves for each dataset, with confidence intervals.")
    lines.append("- Figure 3: Gene-cluster prevalence histograms for whole, clinical, and environmental pangenomes.")
    lines.append("- Figure 4: Clinical-vs-environmental unified-matrix prevalence scatter and volcano plot.")
    lines.append("- Figure 5: Presence/absence heatmap of differential candidate genes.")
    lines.append("- Supplementary Figure S1: Per-genome gene-cluster burden distributions.")
    lines.append("- Supplementary Table S1: Complete per-gene prevalence classification for each pangenome.")
    lines.append("- Supplementary Table S2: Complete unified clinical/environmental differential-prevalence statistics.")

    (outdir / "manuscript_interpretation.md").write_text("\n".join(lines) + "\n")


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate publication-grade pangenome interpretation tables, validation reports, and figures."
    )
    p.add_argument("--whole-rtab", type=Path, required=False, help="Whole pangenome gene_presence_absence.Rtab")
    p.add_argument("--whole-csv", type=Path, required=False, help="Whole pangenome gene_presence_absence.csv")
    p.add_argument("--clinical-rtab", type=Path, required=False, help="Clinical pangenome gene_presence_absence.Rtab")
    p.add_argument("--clinical-csv", type=Path, required=False, help="Clinical pangenome gene_presence_absence.csv")
    p.add_argument("--environmental-rtab", type=Path, required=False, help="Environmental pangenome gene_presence_absence.Rtab")
    p.add_argument("--environmental-csv", type=Path, required=False, help="Environmental pangenome gene_presence_absence.csv")
    p.add_argument("--unified-rtab", type=Path, required=False, help="Unified whole-pangenome matrix for clinical/environmental comparison. Defaults to --whole-rtab.")
    p.add_argument("--unified-csv", type=Path, required=False, help="Unified whole-pangenome CSV annotation file. Defaults to --whole-csv.")
    p.add_argument("--metadata", type=Path, required=False, help="Metadata TSV containing genome IDs and clinical/environmental classification")
    p.add_argument("--metadata-id-col", default="assembly_accession", help="Metadata column used to match matrix columns")
    p.add_argument("--metadata-class-col", default="final_classification", help="Metadata column containing Clinical/Environmental/Unknown classification")
    p.add_argument("--outdir", type=Path, required=True, help="Output directory")
    p.add_argument("--expected-whole-n", type=int, default=1595, help="Expected number of whole genomes; use 0 to disable")
    p.add_argument("--expected-clinical-n", type=int, default=1350, help="Expected number of clinical genomes; use 0 to disable")
    p.add_argument("--expected-environmental-n", type=int, default=99, help="Expected number of environmental genomes; use 0 to disable")
    p.add_argument("--curve-reps", type=int, default=50, help="Random permutations for pan/core accumulation curves")
    p.add_argument("--curve-step", type=int, default=10, help="Record rarefaction every N genomes")
    p.add_argument("--random-seed", type=int, default=42, help="Random seed")
    p.add_argument("--skip-curves", action="store_true", help="Skip randomized pan/core rarefaction curves")
    p.add_argument("--min-abs-delta", type=float, default=0.10, help="Minimum absolute clinical-environmental prevalence difference for candidates")
    p.add_argument("--fdr-cutoff", type=float, default=0.05, help="FDR cutoff for candidates")
    return p.parse_args()


def maybe_expected(x: int) -> Optional[int]:
    return None if x == 0 else x


def main() -> None:
    args = parse_args()
    mkdir(args.outdir)
    mkdir(args.outdir / "tables")
    mkdir(args.outdir / "figures")
    mkdir(args.outdir / "validation")

    datasets_spec = []
    if args.whole_rtab:
        datasets_spec.append(("whole", args.whole_rtab, args.whole_csv, maybe_expected(args.expected_whole_n)))
    if args.clinical_rtab:
        datasets_spec.append(("clinical", args.clinical_rtab, args.clinical_csv, maybe_expected(args.expected_clinical_n)))
    if args.environmental_rtab:
        datasets_spec.append(("environmental", args.environmental_rtab, args.environmental_csv, maybe_expected(args.expected_environmental_n)))

    if not datasets_spec:
        raise SystemExit("No pangenome RTAB files supplied. Provide at least --whole-rtab, --clinical-rtab, or --environmental-rtab.")

    results: Dict[str, DatasetResult] = {}
    for name, rtab, csv, expected_n in datasets_spec:
        print(f"[INFO] Analyzing {name}: {rtab}", file=sys.stderr)
        ds = analyze_dataset(
            name=name,
            rtab_path=rtab,
            csv_path=csv,
            expected_n=expected_n,
            outdir=args.outdir,
            curve_reps=args.curve_reps,
            curve_step=args.curve_step,
            random_seed=args.random_seed,
            skip_curves=args.skip_curves,
        )
        results[name] = ds
        make_per_dataset_figures(ds, args.outdir)

    make_comparison_tables(results, args.outdir)
    make_comparison_figures(results, args.outdir)
    write_validation_report(results, args.outdir)

    unified_rtab = args.unified_rtab or args.whole_rtab
    unified_csv = args.unified_csv or args.whole_csv
    diff = None
    if args.metadata and unified_rtab:
        print("[INFO] Running unified clinical-vs-environmental differential prevalence analysis", file=sys.stderr)
        diff = run_unified_clinical_environmental_analysis(
            unified_rtab=unified_rtab,
            unified_csv=unified_csv,
            metadata_path=args.metadata,
            metadata_id_col=args.metadata_id_col,
            metadata_class_col=args.metadata_class_col,
            outdir=args.outdir,
            min_abs_delta=args.min_abs_delta,
            fdr_cutoff=args.fdr_cutoff,
        )
    else:
        print("[INFO] Unified differential analysis skipped; provide --metadata and --whole-rtab/--unified-rtab to enable it.", file=sys.stderr)

    write_manuscript_interpretation(results, diff, args.outdir)

    print("[DONE] Publication pangenome analysis completed.", file=sys.stderr)
    print(f"[DONE] Results written to: {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
