#!/usr/bin/env python3
"""
Create ONE publication-ready composite QUAST figure with non-overlapping labels.

Designed for QUAST transposed_report.tsv using these fields:
Assembly
# contigs
Largest contig
Total length
GC (%)
N50
N90
auN
L50
L90
# N's per 100 kbp

Output:
  publication_ready_quast_figure.pdf
  publication_ready_quast_figure.svg
  publication_ready_quast_figure.png
  publication_ready_quast_figure.tiff
  quast_cleaned_metrics.tsv
  quast_summary_stats.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLUMN_MAP = {
    "assembly": ["Assembly"],
    "contigs": ["# contigs"],
    "largest_contig": ["Largest contig"],
    "total_length": ["Total length"],
    "gc": ["GC (%)"],
    "n50": ["N50"],
    "n90": ["N90"],
    "aun": ["auN"],
    "l50": ["L50"],
    "l90": ["L90"],
    "ns_per_100kbp": ["# N's per 100 kbp", "N's per 100 kbp"],
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a single publication-ready QUAST composite figure."
    )
    p.add_argument(
        "--report",
        required=True,
        type=Path,
        help="Path to QUAST transposed_report.tsv"
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for outputs"
    )
    return p.parse_args()


def clean_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(
        series.astype(str)
        .str.replace(",", "", regex=False)
        .str.replace(" ", "", regex=False)
        .str.strip()
        .replace({"": np.nan, "-": np.nan, "NA": np.nan, "N/A": np.nan, "nan": np.nan}),
        errors="coerce",
    )


def first_existing(df: pd.DataFrame, cols: list[str]) -> str | None:
    for c in cols:
        if c in df.columns:
            return c
    return None


def load_data(report: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(report, sep="\t", dtype=str)
    if df.empty:
        raise ValueError(f"Empty report: {report}")

    out = pd.DataFrame()
    missing = []

    for key, candidates in COLUMN_MAP.items():
        col = first_existing(df, candidates)
        if col is None:
            missing.append(candidates[0])
            out[key] = np.nan
        elif key == "assembly":
            out[key] = df[col].astype(str).str.strip()
        else:
            out[key] = clean_numeric(df[col])

    out = out.drop_duplicates(subset="assembly", keep="first").reset_index(drop=True)
    return out, missing


def summary_stats(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in df.columns:
        if c == "assembly":
            continue
        vals = df[c].dropna()
        if vals.empty:
            continue
        rows.append({
            "metric": c,
            "n": int(vals.size),
            "median": float(vals.median()),
            "q1": float(vals.quantile(0.25)),
            "q3": float(vals.quantile(0.75)),
            "mean": float(vals.mean()),
            "sd": float(vals.std(ddof=1)),
            "minimum": float(vals.min()),
            "maximum": float(vals.max()),
        })
    return pd.DataFrame(rows)


def median_iqr_text(values: pd.Series, scale: float = 1.0, decimals: int = 2) -> str:
    values = values.dropna() / scale
    if values.empty:
        return "NA"
    med = values.median()
    q1 = values.quantile(0.25)
    q3 = values.quantile(0.75)
    return f"{med:.{decimals}f} [{q1:.{decimals}f}–{q3:.{decimals}f}]"


def safe_bins(values: pd.Series):
    values = values.dropna()
    if len(values) < 2 or values.nunique() < 2:
        return 1
    return "fd"


def format_panel(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def add_hist(ax, values: pd.Series, xlabel: str, panel: str, scale: float = 1.0, note: str | None = None):
    vals = values.dropna() / scale
    if vals.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=9)
        ax.set_axis_off()
        return

    ax.hist(vals, bins=safe_bins(vals), edgecolor="black", linewidth=0.5)
    med = vals.median()
    ax.axvline(med, linestyle="--", linewidth=1.0)

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Genomes", fontsize=9)
    ax.set_title(panel, loc="left", fontweight="bold", fontsize=11)

    annotation = f"n={len(vals):,}\nMedian={med:.2f}"
    if note:
        annotation += f"\n{note}"
    ax.text(
        0.98, 0.98, annotation,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.75, linewidth=0.4)
    )
    format_panel(ax)


def add_scatter(ax, df: pd.DataFrame, panel: str):
    data = df[["contigs", "n50"]].dropna()
    data = data[(data["contigs"] > 0) & (data["n50"] > 0)]

    if data.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=9)
        ax.set_axis_off()
        return

    rho = data["contigs"].corr(data["n50"], method="spearman")

    ax.scatter(
        data["contigs"],
        data["n50"] / 1000.0,
        s=14,
        alpha=0.45,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Contigs (log scale)", fontsize=9)
    ax.set_ylabel("N50 (kb, log scale)", fontsize=9)
    ax.set_title(panel, loc="left", fontweight="bold", fontsize=11)
    ax.text(
        0.98, 0.98,
        f"n={len(data):,}\nSpearman ρ={rho:.3f}",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.75, linewidth=0.4)
    )
    format_panel(ax)


def save_all(fig: plt.Figure, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.titlesize": 13,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    df, missing = load_data(args.report)

    if missing:
        print("Warning: missing expected columns:", file=sys.stderr)
        for col in missing:
            print(f"  - {col}", file=sys.stderr)

    df.to_csv(args.outdir / "quast_cleaned_metrics.tsv", sep="\t", index=False)
    summary_stats(df).to_csv(args.outdir / "quast_summary_stats.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2), constrained_layout=True)

    add_hist(
        axes[0, 0], df["total_length"],
        "Total assembly length (Mb)", "A",
        scale=1_000_000
    )
    add_hist(
        axes[0, 1], df["n50"],
        "N50 (kb)", "B",
        scale=1_000
    )
    add_hist(
        axes[0, 2], df["contigs"],
        "Number of contigs (≥500 bp)", "C"
    )
    add_hist(
        axes[1, 0], df["gc"],
        "GC content (%)", "D"
    )
    add_hist(
        axes[1, 1], df["largest_contig"],
        "Largest contig (Mb)", "E",
        scale=1_000_000
    )
    add_scatter(
        axes[1, 2], df,
        "F"
    )

    n = len(df)
    title = "QUAST assembly-quality profile of the genome collection"
    subtitle = (
        f"n = {n:,} genomes | Panels show core assembly-quality and contiguity metrics"
    )
    fig.suptitle(title, y=1.02, fontweight="bold")
    fig.text(0.5, 0.99, subtitle, ha="center", va="top", fontsize=9)

    caption = (
        "A, total assembly length; B, N50; C, contig count; D, GC content; "
        "E, largest contig; F, relationship between contig count and N50."
    )
    fig.text(0.5, -0.01, caption, ha="center", va="top", fontsize=8)

    save_all(fig, args.outdir / "publication_ready_quast_figure")
    plt.close(fig)

    print(f"Processed {n:,} assemblies.")
    print(f"Output written to: {args.outdir.resolve()}")
    print("Main figure: publication_ready_quast_figure.[pdf|svg|png|tiff]")


if __name__ == "__main__":
    raise SystemExit(main())
