#!/usr/bin/env python3

"""
Publishable CheckM2 quality visualization for genome assemblies.

Input:
    checkm2_pass_genomes.csv

Outputs:
    1. checkm2_pass_quality_multipanel.pdf
    2. checkm2_pass_quality_multipanel.png
    3. checkm2_pass_quality_multipanel.svg
    4. checkm2_pass_quality_summary.csv

Recommended use:
    python3 plot_checkm2_pass_quality.py \
        --input /mnt/phoenix/v4/1_raw_files/8.checkm2_filtered_genomes/reports/checkm2_pass_genomes.csv \
        --outdir /mnt/phoenix/v4/1_raw_files/8.checkm2_filtered_genomes/figures
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def detect_column(df, candidates, required=True):
    """
    Detect a column from a list of possible names.
    """
    lower_map = {c.lower(): c for c in df.columns}

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    if required:
        raise ValueError(
            f"Could not find any of these columns: {candidates}\n"
            f"Available columns: {list(df.columns)}"
        )

    return None


def clean_numeric(series):
    """
    Convert numeric-like columns to float safely.
    """
    return pd.to_numeric(series, errors="coerce")


def main():
    parser = argparse.ArgumentParser(
        description="Create publishable CheckM2 quality figures from checkm2_pass_genomes.csv"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to checkm2_pass_genomes.csv"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for figures and summary table"
    )

    parser.add_argument(
        "--prefix",
        default="checkm2_pass_quality",
        help="Output filename prefix"
    )

    parser.add_argument(
        "--completeness-threshold",
        type=float,
        default=95.0,
        help="Completeness threshold used for standard CheckM2 PASS"
    )

    parser.add_argument(
        "--contamination-threshold",
        type=float,
        default=5.0,
        help="Contamination threshold used for standard CheckM2 PASS"
    )

    parser.add_argument(
        "--strict-completeness",
        type=float,
        default=98.0,
        help="Completeness threshold for strict high-quality subset"
    )

    parser.add_argument(
        "--strict-contamination",
        type=float,
        default=2.0,
        help="Contamination threshold for strict high-quality subset"
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    df = pd.read_csv(input_file)

    # Detect required columns
    genome_col = detect_column(
        df,
        ["Name", "Bin Id", "Bin_Id", "Genome", "genome", "sample", "Sample"],
        required=False
    )

    completeness_col = detect_column(
        df,
        ["Completeness", "completeness"],
        required=True
    )

    contamination_col = detect_column(
        df,
        ["Contamination", "contamination"],
        required=True
    )

    df[completeness_col] = clean_numeric(df[completeness_col])
    df[contamination_col] = clean_numeric(df[contamination_col])

    df = df.dropna(subset=[completeness_col, contamination_col]).copy()

    if df.empty:
        raise ValueError("No valid rows after filtering missing completeness/contamination values.")

    total = len(df)

    standard_pass = (
        (df[completeness_col] >= args.completeness_threshold) &
        (df[contamination_col] <= args.contamination_threshold)
    )

    strict_pass = (
        (df[completeness_col] >= args.strict_completeness) &
        (df[contamination_col] <= args.strict_contamination)
    )

    df["quality_group"] = np.where(
        strict_pass,
        "Strict high-quality",
        np.where(standard_pass, "Standard pass only", "Outside standard pass")
    )

    strict_count = int((df["quality_group"] == "Strict high-quality").sum())
    standard_only_count = int((df["quality_group"] == "Standard pass only").sum())
    outside_count = int((df["quality_group"] == "Outside standard pass").sum())

    summary = pd.DataFrame({
        "metric": [
            "Total genomes",
            "Mean completeness",
            "Median completeness",
            "Minimum completeness",
            "Maximum completeness",
            "Mean contamination",
            "Median contamination",
            "Minimum contamination",
            "Maximum contamination",
            f"Standard pass: completeness >= {args.completeness_threshold} and contamination <= {args.contamination_threshold}",
            f"Strict high-quality: completeness >= {args.strict_completeness} and contamination <= {args.strict_contamination}",
            "Standard pass only",
            "Outside standard pass"
        ],
        "value": [
            total,
            round(df[completeness_col].mean(), 4),
            round(df[completeness_col].median(), 4),
            round(df[completeness_col].min(), 4),
            round(df[completeness_col].max(), 4),
            round(df[contamination_col].mean(), 4),
            round(df[contamination_col].median(), 4),
            round(df[contamination_col].min(), 4),
            round(df[contamination_col].max(), 4),
            int(standard_pass.sum()),
            strict_count,
            standard_only_count,
            outside_count
        ]
    })

    summary_file = outdir / f"{args.prefix}_summary.csv"
    summary.to_csv(summary_file, index=False)

    # Publication-friendly Matplotlib settings
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "figure.titlesize": 12,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none"
    })

    # Colorblind-aware palette
    color_main = "#0072B2"
    color_secondary = "#D55E00"
    color_strict = "#009E73"
    color_standard = "#E69F00"
    color_outside = "#CC79A7"
    color_threshold = "#333333"

    fig, axes = plt.subplots(
        2, 2,
        figsize=(7.4, 6.2),
        constrained_layout=True
    )

    ax1, ax2, ax3, ax4 = axes.flatten()

    # Panel A: Completeness histogram
    ax1.hist(
        df[completeness_col],
        bins=30,
        color=color_main,
        edgecolor="white",
        linewidth=0.5
    )
    ax1.axvline(
        args.completeness_threshold,
        color=color_threshold,
        linestyle="--",
        linewidth=1.1,
        label=f"{args.completeness_threshold:.0f}% threshold"
    )
    ax1.axvline(
        df[completeness_col].median(),
        color=color_secondary,
        linestyle="-",
        linewidth=1.1,
        label=f"Median = {df[completeness_col].median():.2f}%"
    )
    ax1.set_xlabel("Completeness (%)")
    ax1.set_ylabel("Genome count")
    ax1.set_title("A. Completeness distribution")
    ax1.legend(frameon=False)
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Panel B: Contamination histogram
    ax2.hist(
        df[contamination_col],
        bins=30,
        color=color_secondary,
        edgecolor="white",
        linewidth=0.5
    )
    ax2.axvline(
        args.contamination_threshold,
        color=color_threshold,
        linestyle="--",
        linewidth=1.1,
        label=f"{args.contamination_threshold:.0f}% threshold"
    )
    ax2.axvline(
        df[contamination_col].median(),
        color=color_main,
        linestyle="-",
        linewidth=1.1,
        label=f"Median = {df[contamination_col].median():.2f}%"
    )
    ax2.set_xlabel("Contamination (%)")
    ax2.set_ylabel("Genome count")
    ax2.set_title("B. Contamination distribution")
    ax2.legend(frameon=False)
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Panel C: Completeness vs contamination scatter
    group_order = [
        "Strict high-quality",
        "Standard pass only",
        "Outside standard pass"
    ]

    group_colors = {
        "Strict high-quality": color_strict,
        "Standard pass only": color_standard,
        "Outside standard pass": color_outside
    }

    for group in group_order:
        sub = df[df["quality_group"] == group]
        if sub.empty:
            continue
        ax3.scatter(
            sub[contamination_col],
            sub[completeness_col],
            s=12,
            alpha=0.70,
            linewidths=0,
            color=group_colors[group],
            label=f"{group} (n={len(sub)})"
        )

    ax3.axhline(
        args.completeness_threshold,
        color=color_threshold,
        linestyle="--",
        linewidth=1.0
    )
    ax3.axvline(
        args.contamination_threshold,
        color=color_threshold,
        linestyle="--",
        linewidth=1.0
    )
    ax3.set_xlabel("Contamination (%)")
    ax3.set_ylabel("Completeness (%)")
    ax3.set_title("C. Completeness-contamination profile")
    ax3.legend(frameon=False, loc="lower left")

    # Panel D: Quality class count
    class_counts = pd.Series({
        "Strict\nhigh-quality": strict_count,
        "Standard\npass only": standard_only_count,
        "Outside\nstandard pass": outside_count
    })

    bar_colors = [color_strict, color_standard, color_outside]

    bars = ax4.bar(
        class_counts.index,
        class_counts.values,
        color=bar_colors,
        edgecolor="black",
        linewidth=0.5
    )

    ax4.set_ylabel("Genome count")
    ax4.set_title("D. CheckM2 quality classification")
    ax4.yaxis.set_major_locator(MaxNLocator(integer=True))

    for bar in bars:
        height = bar.get_height()
        ax4.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f"{int(height)}",
            ha="center",
            va="bottom",
            fontsize=9
        )

    # Clean journal-style axes
    for ax in [ax1, ax2, ax3, ax4]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

    fig.suptitle(
        f"CheckM2 quality assessment of QUAST-pass genomes (n = {total})",
        y=1.02,
        fontweight="bold"
    )

    pdf_file = outdir / f"{args.prefix}_multipanel.pdf"
    png_file = outdir / f"{args.prefix}_multipanel.png"
    svg_file = outdir / f"{args.prefix}_multipanel.svg"

    fig.savefig(pdf_file, bbox_inches="tight")
    fig.savefig(svg_file, bbox_inches="tight")
    fig.savefig(png_file, dpi=600, bbox_inches="tight")

    plt.close(fig)

    print("CheckM2 visualization completed.")
    print(f"Input file: {input_file}")
    print(f"Total genomes plotted: {total}")
    print(f"Strict high-quality genomes: {strict_count}")
    print(f"Standard pass only genomes: {standard_only_count}")
    print(f"Outside standard pass genomes: {outside_count}")
    print()
    print("Output files:")
    print(f"PDF: {pdf_file}")
    print(f"PNG: {png_file}")
    print(f"SVG: {svg_file}")
    print(f"Summary CSV: {summary_file}")


if __name__ == "__main__":
    main()
