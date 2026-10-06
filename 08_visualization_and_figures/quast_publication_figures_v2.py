#!/usr/bin/env python3
"""
Generate publication-ready figures, tables, and manuscript-ready statistics
from a QUAST transposed_report.tsv file.

Exact QUAST fields used:
Assembly; contig counts at >=0/1k/5k/10k/25k/50k bp; retained lengths at
those thresholds; # contigs; Largest contig; Total length; GC (%); N50; N90;
auN; L50; L90; and # N's per 100 kbp.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLUMN_MAP = {
    "assembly": ["Assembly"],
    "contigs_ge_0": ["# contigs (>= 0 bp)"],
    "contigs_ge_1000": ["# contigs (>= 1000 bp)"],
    "contigs_ge_5000": ["# contigs (>= 5000 bp)"],
    "contigs_ge_10000": ["# contigs (>= 10000 bp)"],
    "contigs_ge_25000": ["# contigs (>= 25000 bp)"],
    "contigs_ge_50000": ["# contigs (>= 50000 bp)"],
    "length_ge_0": ["Total length (>= 0 bp)"],
    "length_ge_1000": ["Total length (>= 1000 bp)"],
    "length_ge_5000": ["Total length (>= 5000 bp)"],
    "length_ge_10000": ["Total length (>= 10000 bp)"],
    "length_ge_25000": ["Total length (>= 25000 bp)"],
    "length_ge_50000": ["Total length (>= 50000 bp)"],
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

THRESHOLD_SERIES = [
    ("0 bp", "contigs_ge_0", "length_ge_0"),
    ("1 kb", "contigs_ge_1000", "length_ge_1000"),
    ("5 kb", "contigs_ge_5000", "length_ge_5000"),
    ("10 kb", "contigs_ge_10000", "length_ge_10000"),
    ("25 kb", "contigs_ge_25000", "length_ge_25000"),
    ("50 kb", "contigs_ge_50000", "length_ge_50000"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create publication-ready QUAST figures and statistics."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--report", type=Path,
                        help="Exact QUAST transposed_report.tsv file.")
    source.add_argument("--quast-dir", type=Path,
                        help="QUAST output directory.")
    parser.add_argument("--outdir", required=True, type=Path)

    parser.add_argument("--max-contigs", type=float, default=None)
    parser.add_argument("--min-n50", type=float, default=None)
    parser.add_argument("--min-length", type=float, default=None)
    parser.add_argument("--max-length", type=float, default=None)
    parser.add_argument("--min-gc", type=float, default=None)
    parser.add_argument("--max-gc", type=float, default=None)
    parser.add_argument("--max-ns-per-100kbp", type=float, default=None)
    return parser.parse_args()


def locate_report(args: argparse.Namespace) -> Path:
    if args.report is not None:
        report = args.report.resolve()
        if not report.is_file():
            raise FileNotFoundError(f"QUAST report not found: {report}")
        return report

    quast_dir = args.quast_dir.resolve()
    direct = quast_dir / "transposed_report.tsv"
    if direct.is_file():
        return direct

    matches = sorted(quast_dir.rglob("transposed_report.tsv"))
    matches = [p for p in matches if ".quast_figures_env" not in p.parts]
    if not matches:
        raise FileNotFoundError(f"No transposed_report.tsv found under: {quast_dir}")
    if len(matches) > 1:
        listing = "\n".join(f"  - {p}" for p in matches[:20])
        raise RuntimeError(
            "Multiple transposed_report.tsv files were found. Use --report:\n"
            f"{listing}"
        )
    return matches[0]


def clean_numeric(series: pd.Series) -> pd.Series:
    cleaned = (
        series.astype(str)
        .str.replace(",", "", regex=False)
        .str.replace(" ", "", regex=False)
        .str.strip()
        .replace({"": np.nan, "-": np.nan, "NA": np.nan,
                  "N/A": np.nan, "nan": np.nan, "None": np.nan})
    )
    return pd.to_numeric(cleaned, errors="coerce")


def load_metrics(report: Path) -> tuple[pd.DataFrame, list[str]]:
    raw = pd.read_csv(report, sep="\t", dtype=str)
    if raw.empty:
        raise ValueError(f"The QUAST report is empty: {report}")

    cleaned = pd.DataFrame()
    missing: list[str] = []
    for output_name, candidates in COLUMN_MAP.items():
        column = next((c for c in candidates if c in raw.columns), None)
        if column is None:
            cleaned[output_name] = np.nan
            missing.append(candidates[0])
        elif output_name == "assembly":
            cleaned[output_name] = raw[column].astype(str).str.strip()
        else:
            cleaned[output_name] = clean_numeric(raw[column])

    duplicate_count = int(cleaned["assembly"].duplicated().sum())
    if duplicate_count:
        print(
            f"Warning: {duplicate_count:,} duplicate Assembly IDs found; "
            "retaining the first occurrence.", file=sys.stderr
        )
        cleaned = cleaned.drop_duplicates("assembly", keep="first")

    return cleaned.reset_index(drop=True), missing


def make_summary_table(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for column in metrics.columns:
        if column == "assembly":
            continue
        values = metrics[column].dropna()
        if values.empty:
            continue
        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        rows.append({
            "metric": column,
            "n": int(values.size),
            "missing": int(metrics[column].isna().sum()),
            "mean": values.mean(),
            "sd": values.std(ddof=1),
            "median": values.median(),
            "q1": q1,
            "q3": q3,
            "iqr": q3 - q1,
            "minimum": values.min(),
            "maximum": values.max(),
        })
    return pd.DataFrame(rows)


def safe_hist_bins(values: pd.Series):
    values = values.dropna()
    return 1 if values.size < 2 or values.nunique() < 2 else "fd"


def add_histogram(ax, values, xlabel, panel_label, scale=1.0):
    plotted = values.dropna() / scale
    if plotted.empty:
        ax.text(0.5, 0.5, "Metric unavailable", ha="center", va="center")
        ax.set_axis_off()
        return

    ax.hist(plotted, bins=safe_hist_bins(plotted),
            edgecolor="black", linewidth=0.45)
    median = plotted.median()
    ax.axvline(median, linestyle="--", linewidth=1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of genomes")
    ax.text(0.02, 0.96, panel_label, transform=ax.transAxes,
            ha="left", va="top", fontweight="bold", fontsize=11)
    ax.text(0.98, 0.96,
            f"n = {len(plotted):,}\nMedian = {median:,.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})


def plot_core_quality(metrics, outdir):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
    add_histogram(axes[0, 0], metrics["total_length"],
                  "Total assembly length (Mb)", "A", 1_000_000)
    add_histogram(axes[0, 1], metrics["n50"], "N50 (kb)", "B", 1_000)
    add_histogram(axes[1, 0], metrics["contigs"],
                  "Number of contigs (≥500 bp)", "C")
    add_histogram(axes[1, 1], metrics["gc"], "GC content (%)", "D")
    fig.suptitle("Assembly-quality characteristics of the genome collection", y=1.01)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_1_QUAST_core_quality")
    plt.close(fig)


def plot_contiguity_metrics(metrics, outdir):
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 6.2))
    add_histogram(axes[0, 0], metrics["largest_contig"],
                  "Largest contig (Mb)", "A", 1_000_000)
    add_histogram(axes[0, 1], metrics["n90"], "N90 (kb)", "B", 1_000)
    add_histogram(axes[0, 2], metrics["aun"], "auN (kb)", "C", 1_000)
    add_histogram(axes[1, 0], metrics["l50"], "L50", "D")
    add_histogram(axes[1, 1], metrics["l90"], "L90", "E")
    add_histogram(axes[1, 2], metrics["ns_per_100kbp"],
                  "Ambiguous bases per 100 kbp", "F")
    fig.suptitle("Additional QUAST contiguity and sequence-quality metrics", y=1.01)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_2_QUAST_contiguity_metrics")
    plt.close(fig)


def plot_threshold_profile(metrics, outdir):
    labels, cm, cq1, cq3, lm, lq1, lq3 = [], [], [], [], [], [], []
    for label, ccol, lcol in THRESHOLD_SERIES:
        labels.append(label)
        c = metrics[ccol].dropna()
        l = metrics[lcol].dropna() / 1_000_000
        cm.append(c.median()); cq1.append(c.quantile(0.25)); cq3.append(c.quantile(0.75))
        lm.append(l.median()); lq1.append(l.quantile(0.25)); lq3.append(l.quantile(0.75))

    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))

    axes[0].plot(x, cm, marker="o", linewidth=1.5)
    axes[0].fill_between(x, cq1, cq3, alpha=0.2)
    axes[0].set_xticks(x, labels)
    axes[0].set_xlabel("Minimum contig length threshold")
    axes[0].set_ylabel("Median number of contigs")
    axes[0].set_title("A  Contig-count profile", loc="left", fontweight="bold")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    axes[1].plot(x, lm, marker="o", linewidth=1.5)
    axes[1].fill_between(x, lq1, lq3, alpha=0.2)
    axes[1].set_xticks(x, labels)
    axes[1].set_xlabel("Minimum contig length threshold")
    axes[1].set_ylabel("Median retained assembly length (Mb)")
    axes[1].set_title("B  Retained-length profile", loc="left", fontweight="bold")
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    fig.suptitle("Effect of minimum contig-size thresholds on assembly statistics", y=1.03)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_3_QUAST_size_threshold_profile")
    plt.close(fig)


def plot_n50_vs_contigs(metrics, outdir):
    data = metrics[["contigs", "n50"]].dropna()
    data = data[(data["contigs"] > 0) & (data["n50"] > 0)]
    if data.empty:
        return float("nan")

    rho = float(data["contigs"].corr(data["n50"], method="spearman"))
    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    ax.scatter(data["contigs"], data["n50"] / 1_000,
               s=14, alpha=0.45, linewidths=0, rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of contigs (log scale)")
    ax.set_ylabel("N50 (kb; log scale)")
    ax.set_title("Assembly fragmentation versus contiguity")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(0.03, 0.97,
            f"n = {len(data):,}\nSpearman ρ = {rho:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=9)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_4_QUAST_N50_vs_contigs")
    plt.close(fig)
    return rho


def classify_qc(metrics, args):
    values = [args.max_contigs, args.min_n50, args.min_length,
              args.max_length, args.min_gc, args.max_gc,
              args.max_ns_per_100kbp]
    if all(v is None for v in values):
        return None

    qc = metrics.copy()
    qc["qc_status"] = "PASS"
    qc["qc_reason"] = ""
    checks = []
    if args.max_contigs is not None:
        checks.append(("contigs", qc["contigs"] <= args.max_contigs,
                       f"contigs>{args.max_contigs:g}"))
    if args.min_n50 is not None:
        checks.append(("n50", qc["n50"] >= args.min_n50,
                       f"N50<{args.min_n50:g}"))
    if args.min_length is not None:
        checks.append(("total_length", qc["total_length"] >= args.min_length,
                       f"length<{args.min_length:g}"))
    if args.max_length is not None:
        checks.append(("total_length", qc["total_length"] <= args.max_length,
                       f"length>{args.max_length:g}"))
    if args.min_gc is not None:
        checks.append(("gc", qc["gc"] >= args.min_gc,
                       f"GC<{args.min_gc:g}"))
    if args.max_gc is not None:
        checks.append(("gc", qc["gc"] <= args.max_gc,
                       f"GC>{args.max_gc:g}"))
    if args.max_ns_per_100kbp is not None:
        checks.append(("ns_per_100kbp",
                       qc["ns_per_100kbp"] <= args.max_ns_per_100kbp,
                       f"Ns_per_100kbp>{args.max_ns_per_100kbp:g}"))

    for metric, passed, reason in checks:
        missing = qc[metric].isna()
        failed = (~passed) & (~missing)
        qc.loc[missing & (qc["qc_status"] == "PASS"), "qc_status"] = "REVIEW"
        qc.loc[missing, "qc_reason"] += f"missing_{metric};"
        qc.loc[failed, "qc_status"] = "FAIL"
        qc.loc[failed, "qc_reason"] += f"{reason};"

    qc["qc_reason"] = qc["qc_reason"].str.rstrip(";")
    return qc


def plot_qc_status(qc, outdir):
    order = ["PASS", "REVIEW", "FAIL"]
    counts = qc["qc_status"].value_counts().reindex(order, fill_value=0)
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    bars = ax.bar(counts.index, counts.values,
                  edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Number of genomes")
    ax.set_title("QUAST quality-control classification")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    total = int(counts.sum())
    for bar, value in zip(bars, counts.values):
        pct = 100 * value / total if total else 0
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f"{value:,}\n({pct:.1f}%)",
                ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_5_QUAST_QC_status")
    plt.close(fig)


def median_iqr(values, scale=1.0):
    values = values.dropna() / scale
    if values.empty:
        return "unavailable"
    return (f"{values.median():,.2f} "
            f"[{values.quantile(0.25):,.2f}–{values.quantile(0.75):,.2f}]")


def write_manuscript_summary(metrics, report, rho, outpath):
    n = len(metrics)
    complete = int((metrics["contigs"] == 1).sum())
    complete_pct = 100 * complete / n if n else 0
    lines = [
        "QUAST PUBLICATION SUMMARY",
        "=" * 72,
        f"Source report: {report}",
        f"Unique assemblies analysed: {n:,}",
        "",
        "Recommended descriptive statistics: median [Q1–Q3]",
        f"Total assembly length (Mb): {median_iqr(metrics['total_length'], 1_000_000)}",
        f"GC content (%): {median_iqr(metrics['gc'])}",
        f"Number of contigs: {median_iqr(metrics['contigs'])}",
        f"Largest contig (Mb): {median_iqr(metrics['largest_contig'], 1_000_000)}",
        f"N50 (kb): {median_iqr(metrics['n50'], 1_000)}",
        f"N90 (kb): {median_iqr(metrics['n90'], 1_000)}",
        f"auN (kb): {median_iqr(metrics['aun'], 1_000)}",
        f"L50: {median_iqr(metrics['l50'])}",
        f"L90: {median_iqr(metrics['l90'])}",
        f"Ambiguous bases per 100 kbp: {median_iqr(metrics['ns_per_100kbp'])}",
        "",
        f"Single-contig assemblies: {complete:,}/{n:,} ({complete_pct:.2f}%)",
        f"Spearman correlation, contig count versus N50: {rho:.4f}",
        "",
        "Manuscript-ready Results template:",
        (
            f"QUAST analysis was performed on {n:,} assemblies. The median total "
            f"assembly length was {median_iqr(metrics['total_length'], 1_000_000)} Mb, "
            f"with a median GC content of {median_iqr(metrics['gc'])}%. Assemblies "
            f"contained a median of {median_iqr(metrics['contigs'])} contigs, and the "
            f"median N50 was {median_iqr(metrics['n50'], 1_000)} kb. "
            f"{complete:,} assemblies ({complete_pct:.2f}%) consisted of a single "
            f"contig. Contig count and N50 showed a Spearman correlation of "
            f"ρ = {rho:.3f}."
        ),
        "",
        "Verify all values and adapt wording before manuscript submission.",
    ]
    outpath.write_text("\n".join(lines) + "\n")


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    report = locate_report(args)
    print(f"Using QUAST report: {report}")

    metrics, missing_columns = load_metrics(report)
    if missing_columns:
        print("Warning: expected QUAST columns not found:", file=sys.stderr)
        for column in missing_columns:
            print(f"  - {column}", file=sys.stderr)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    metrics.to_csv(args.outdir / "quast_metrics_cleaned.tsv",
                   sep="\t", index=False)
    make_summary_table(metrics).to_csv(
        args.outdir / "quast_summary_statistics.tsv", sep="\t", index=False)

    plot_core_quality(metrics, args.outdir)
    plot_contiguity_metrics(metrics, args.outdir)
    plot_threshold_profile(metrics, args.outdir)
    rho = plot_n50_vs_contigs(metrics, args.outdir)

    qc = classify_qc(metrics, args)
    if qc is not None:
        qc.to_csv(args.outdir / "quast_qc_classification.tsv",
                  sep="\t", index=False)
        plot_qc_status(qc, args.outdir)

    write_manuscript_summary(
        metrics, report, rho,
        args.outdir / "quast_manuscript_summary.txt")

    print(f"Processed {len(metrics):,} unique assemblies.")
    print(f"Results written to: {args.outdir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
