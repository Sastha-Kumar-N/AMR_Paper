#!/usr/bin/env python3

"""
Create one publication-ready composite figure from a GTDB-Tk summary table.

Panels
------
A. GTDB-Tk species-classification outcome
B. Closest-reference ANI versus alignment fraction
C. ANI margin relative to the species-specific GTDB boundary
D. Marker-gene alignment coverage and warning summary

The script produces one figure in PNG, TIFF, PDF, and SVG formats.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def extract_taxon(classification, prefix):
    """Extract a taxonomic rank from a GTDB taxonomy string."""
    for taxon in str(classification).split(";"):
        taxon = taxon.strip()

        if taxon.startswith(prefix):
            value = taxon[len(prefix):].strip()
            return value if value else "Unresolved"

    return "Unresolved"


def numeric_column(dataframe, column):
    """Safely convert a GTDB-Tk column to numeric values."""
    if column not in dataframe.columns:
        return pd.Series(
            np.nan,
            index=dataframe.index,
            dtype=float,
        )

    missing_values = {
        "": np.nan,
        "N/A": np.nan,
        "NA": np.nan,
        "None": np.nan,
        "none": np.nan,
        "nan": np.nan,
    }

    return pd.to_numeric(
        dataframe[column].replace(missing_values),
        errors="coerce",
    )


def abbreviated_species_name(species):
    """Convert Pseudomonas aeruginosa to P. aeruginosa."""
    parts = str(species).split()

    if len(parts) >= 2:
        return f"{parts[0][0]}. {' '.join(parts[1:])}"

    return species


def add_panel_label(axis, label):
    """Add a publication-style panel label."""
    axis.text(
        -0.12,
        1.07,
        label,
        transform=axis.transAxes,
        fontsize=17,
        fontweight="bold",
        va="top",
        ha="left",
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate one publication-ready GTDB-Tk composite figure."
        )
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Path to gtdbtk.bac120.summary.tsv",
    )

    parser.add_argument(
        "--target-species",
        default="Pseudomonas aeruginosa",
        help="Expected target species without the s__ prefix",
    )

    parser.add_argument(
        "--outdir",
        default="gtdbtk_composite_figure",
        help="Output directory",
    )

    parser.add_argument(
        "--prefix",
        default="GTDBTk_taxonomic_confirmation",
        help="Output filename prefix",
    )

    args = parser.parse_args()

    summary_file = Path(args.summary)
    output_directory = Path(args.outdir)
    output_directory.mkdir(parents=True, exist_ok=True)

    if not summary_file.exists():
        raise SystemExit(
            f"ERROR: Summary file not found: {summary_file}"
        )

    dataframe = pd.read_csv(
        summary_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    if "classification" not in dataframe.columns:
        raise SystemExit(
            "ERROR: The GTDB-Tk table does not contain "
            "a 'classification' column."
        )

    # ----------------------------------------------------------
    # Extract species assignments
    # ----------------------------------------------------------
    dataframe["species"] = dataframe["classification"].apply(
        lambda value: extract_taxon(value, "s__")
    )

    target_species = args.target_species
    target_short = abbreviated_species_name(target_species)

    dataframe["is_target"] = (
        dataframe["species"] == target_species
    )

    dataframe["classification_outcome"] = np.select(
        [
            dataframe["species"].eq(target_species),
            dataframe["species"].eq("Unresolved"),
        ],
        [
            target_short,
            "Species unresolved",
        ],
        default="Other species",
    )

    # ----------------------------------------------------------
    # Numeric GTDB-Tk fields
    # ----------------------------------------------------------
    numeric_fields = [
        "closest_genome_ani",
        "closest_genome_af",
        "closest_genome_reference_radius",
        "msa_percent",
        "red_value",
    ]

    for field in numeric_fields:
        dataframe[field] = numeric_column(dataframe, field)

    dataframe["ani_margin"] = (
        dataframe["closest_genome_ani"]
        - dataframe["closest_genome_reference_radius"]
    )

    # ----------------------------------------------------------
    # Warning status
    # ----------------------------------------------------------
    if "warnings" in dataframe.columns:
        no_warning_values = {
            "",
            "N/A",
            "NA",
            "None",
            "none",
            "nan",
        }

        dataframe["has_warning"] = ~dataframe[
            "warnings"
        ].str.strip().isin(no_warning_values)

    else:
        dataframe["has_warning"] = False

    # ----------------------------------------------------------
    # Publication style
    # ----------------------------------------------------------
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.9,
            "savefig.dpi": 600,
        }
    )

    target_colour = "#2369A1"
    other_colour = "#C9473D"
    unresolved_colour = "#D5A229"
    neutral_colour = "#7A7A7A"
    light_neutral = "#B8B8B8"

    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14.2, 10.4),
        constrained_layout=True,
    )

    ax_a = axes[0, 0]
    ax_b = axes[0, 1]
    ax_c = axes[1, 0]
    ax_d = axes[1, 1]

    # ==========================================================
    # Panel A: classification outcome
    # ==========================================================
    outcome_order = [
        target_short,
        "Other species",
        "Species unresolved",
    ]

    outcome_counts = (
        dataframe["classification_outcome"]
        .value_counts()
        .reindex(outcome_order, fill_value=0)
    )

    outcome_counts = outcome_counts[
        outcome_counts > 0
    ]

    total_genomes = len(dataframe)

    outcome_colours = {
        target_short: target_colour,
        "Other species": other_colour,
        "Species unresolved": unresolved_colour,
    }

    bars = ax_a.barh(
        outcome_counts.index,
        outcome_counts.values,
        color=[
            outcome_colours.get(label, neutral_colour)
            for label in outcome_counts.index
        ],
        edgecolor="black",
        linewidth=0.5,
        height=0.62,
    )

    # Log scale is used because target and excluded counts may
    # differ by several orders of magnitude.
    ax_a.set_xscale("log")

    ax_a.set_xlabel("Number of genomes, logarithmic scale")
    ax_a.set_title(
        "GTDB-Tk species-classification outcome",
        loc="left",
        fontweight="bold",
    )

    ax_a.invert_yaxis()
    ax_a.grid(
        axis="x",
        linestyle=":",
        linewidth=0.7,
        alpha=0.55,
    )

    maximum_count = outcome_counts.max()

    for bar, count in zip(bars, outcome_counts.values):
        percentage = (
            count / total_genomes * 100
            if total_genomes > 0
            else 0
        )

        ax_a.text(
            count * 1.08,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,} ({percentage:.2f}%)",
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold",
        )

    ax_a.set_xlim(
        left=0.8,
        right=max(maximum_count * 3.2, 10),
    )

    # Add a compact summary of classification methods.
    if "classification_method" in dataframe.columns:
        method_series = dataframe[
            "classification_method"
        ].replace("", "Not reported")

        method_counts = method_series.value_counts()

        method_lines = ["Classification evidence"]

        for method, count in method_counts.head(4).items():
            shortened_method = str(method)

            if len(shortened_method) > 34:
                shortened_method = (
                    shortened_method[:31] + "..."
                )

            method_lines.append(
                f"{shortened_method}: {count:,}"
            )

        ax_a.text(
            0.98,
            0.05,
            "\n".join(method_lines),
            transform=ax_a.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.5,
            bbox={
                "boxstyle": "round,pad=0.45",
                "facecolor": "white",
                "edgecolor": "#9A9A9A",
                "alpha": 0.93,
            },
        )

    add_panel_label(ax_a, "A")

    # ==========================================================
    # Panel B: ANI versus alignment fraction
    # ==========================================================
    ani_af = dataframe.dropna(
        subset=[
            "closest_genome_ani",
            "closest_genome_af",
        ]
    ).copy()

    target_points = ani_af[
        ani_af["is_target"]
    ]

    other_points = ani_af[
        ~ani_af["is_target"]
    ]

    if not target_points.empty:
        ax_b.scatter(
            target_points["closest_genome_af"],
            target_points["closest_genome_ani"],
            s=24,
            alpha=0.35,
            color=target_colour,
            edgecolors="none",
            rasterized=True,
            label=(
                f"{target_short} "
                f"(n={len(target_points):,})"
            ),
        )

    if not other_points.empty:
        ax_b.scatter(
            other_points["closest_genome_af"],
            other_points["closest_genome_ani"],
            s=52,
            alpha=0.90,
            color=other_colour,
            edgecolors="black",
            linewidths=0.55,
            rasterized=True,
            label=(
                "Other or unresolved "
                f"(n={len(other_points):,})"
            ),
            zorder=4,
        )

    reference_radii = target_points[
        "closest_genome_reference_radius"
    ].dropna()

    if not reference_radii.empty:
        median_radius = reference_radii.median()

        ax_b.axhline(
            median_radius,
            color="black",
            linestyle="--",
            linewidth=1.2,
            label=(
                "Median species-specific ANI radius "
                f"({median_radius:.2f}%)"
            ),
        )

    ax_b.set_xlabel("Alignment fraction to closest GTDB reference")
    ax_b.set_ylabel("ANI to closest GTDB reference (%)")

    ax_b.set_title(
        "Reference-genome similarity supporting classification",
        loc="left",
        fontweight="bold",
    )

    if not ani_af.empty:
        af_maximum = ani_af["closest_genome_af"].max()

        if af_maximum <= 1.1:
            ax_b.set_xlim(
                max(0, ani_af["closest_genome_af"].min() - 0.03),
                1.02,
            )

    ax_b.grid(
        linestyle=":",
        linewidth=0.65,
        alpha=0.45,
    )

    ax_b.legend(
        frameon=False,
        loc="lower right",
    )

    add_panel_label(ax_b, "B")

    # ==========================================================
    # Panel C: ANI margin from species-specific radius
    # ==========================================================
    margin_values = dataframe[
        "ani_margin"
    ].dropna()

    if not margin_values.empty:
        number_of_bins = min(
            45,
            max(15, int(np.sqrt(len(margin_values)))),
        )

        ax_c.hist(
            margin_values,
            bins=number_of_bins,
            color=target_colour,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.45,
        )

        ax_c.axvline(
            0,
            color=other_colour,
            linestyle="--",
            linewidth=1.7,
            label="Species-specific ANI boundary",
        )

        median_margin = margin_values.median()
        below_boundary = int(
            (margin_values < 0).sum()
        )

        ax_c.axvline(
            median_margin,
            color="black",
            linestyle=":",
            linewidth=1.4,
            label=(
                f"Median ANI margin: "
                f"{median_margin:.2f} percentage points"
            ),
        )

        ax_c.text(
            0.97,
            0.94,
            (
                f"Genomes below ANI boundary: "
                f"{below_boundary:,}\n"
                f"Genomes with ANI-radius data: "
                f"{len(margin_values):,}"
            ),
            transform=ax_c.transAxes,
            ha="right",
            va="top",
            fontsize=9.2,
            bbox={
                "boxstyle": "round,pad=0.4",
                "facecolor": "white",
                "edgecolor": "#A0A0A0",
                "alpha": 0.93,
            },
        )

        ax_c.legend(
            frameon=False,
            loc="upper left",
        )

    else:
        ax_c.text(
            0.5,
            0.5,
            "ANI-radius information was not available",
            transform=ax_c.transAxes,
            ha="center",
            va="center",
        )

    ax_c.set_xlabel(
        "ANI minus species-specific reference radius "
        "(percentage points)"
    )

    ax_c.set_ylabel("Number of genomes")

    ax_c.set_title(
        "ANI support relative to the GTDB species boundary",
        loc="left",
        fontweight="bold",
    )

    ax_c.grid(
        axis="y",
        linestyle=":",
        linewidth=0.65,
        alpha=0.45,
    )

    add_panel_label(ax_c, "C")

    # ==========================================================
    # Panel D: marker MSA coverage and warnings
    # ==========================================================
    msa_values = dataframe[
        "msa_percent"
    ].dropna()

    if not msa_values.empty:
        number_of_bins = min(
            40,
            max(15, int(np.sqrt(len(msa_values)))),
        )

        histogram_values, bin_edges, patches = ax_d.hist(
            msa_values,
            bins=number_of_bins,
            color=neutral_colour,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.45,
        )

        median_msa = msa_values.median()

        ax_d.axvline(
            median_msa,
            color=target_colour,
            linestyle="--",
            linewidth=1.7,
            label=f"Median MSA coverage: {median_msa:.2f}%",
        )

        warning_msa = dataframe.loc[
            dataframe["has_warning"],
            "msa_percent",
        ].dropna()

        if (
            not warning_msa.empty
            and len(histogram_values) > 0
        ):
            rug_height = max(histogram_values) * 0.035

            ax_d.vlines(
                warning_msa,
                ymin=0,
                ymax=rug_height,
                color=other_colour,
                linewidth=0.8,
                alpha=0.75,
                label="Genomes with GTDB-Tk warnings",
            )

        warning_count = int(
            dataframe["has_warning"].sum()
        )

        warning_percentage = (
            warning_count / total_genomes * 100
            if total_genomes > 0
            else 0
        )

        ax_d.text(
            0.97,
            0.94,
            (
                f"Median marker coverage: {median_msa:.2f}%\n"
                f"Warnings: {warning_count:,} "
                f"({warning_percentage:.2f}%)"
            ),
            transform=ax_d.transAxes,
            ha="right",
            va="top",
            fontsize=9.3,
            bbox={
                "boxstyle": "round,pad=0.4",
                "facecolor": "white",
                "edgecolor": "#A0A0A0",
                "alpha": 0.93,
            },
        )

        ax_d.legend(
            frameon=False,
            loc="upper left",
        )

    else:
        ax_d.text(
            0.5,
            0.5,
            "Marker-gene MSA coverage was not available",
            transform=ax_d.transAxes,
            ha="center",
            va="center",
        )

    ax_d.set_xlabel("GTDB marker-gene MSA coverage (%)")
    ax_d.set_ylabel("Number of genomes")

    ax_d.set_title(
        "Marker alignment coverage and warning status",
        loc="left",
        fontweight="bold",
    )

    ax_d.grid(
        axis="y",
        linestyle=":",
        linewidth=0.65,
        alpha=0.45,
    )

    add_panel_label(ax_d, "D")

    # ----------------------------------------------------------
    # Overall title and footer
    # ----------------------------------------------------------
    target_count = int(
        dataframe["is_target"].sum()
    )

    target_percentage = (
        target_count / total_genomes * 100
        if total_genomes > 0
        else 0
    )

    figure.suptitle(
        (
            "GTDB-Tk taxonomic confirmation of the curated "
            r"$\it{Pseudomonas\ aeruginosa}$ genome collection"
        ),
        fontsize=17,
        fontweight="bold",
    )

    figure.text(
        0.5,
        -0.012,
        (
            f"GTDB-Tk bacterial marker analysis: "
            f"{target_count:,} of {total_genomes:,} genomes "
            f"({target_percentage:.2f}%) were classified as "
            f"{target_species}."
        ),
        ha="center",
        va="bottom",
        fontsize=10.5,
    )

    # ----------------------------------------------------------
    # Export annotated table and one composite figure
    # ----------------------------------------------------------
    dataframe.to_csv(
        output_directory / "gtdbtk_plotting_data.tsv",
        sep="\t",
        index=False,
    )

    output_base = output_directory / args.prefix

    figure.savefig(
        f"{output_base}.png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )

    figure.savefig(
        f"{output_base}.tiff",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )

    figure.savefig(
        f"{output_base}.pdf",
        bbox_inches="tight",
        facecolor="white",
    )

    figure.savefig(
        f"{output_base}.svg",
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(figure)

    # ----------------------------------------------------------
    # Export statistics used in the figure
    # ----------------------------------------------------------
    summary_metrics = pd.DataFrame(
        [
            ["Total GTDB-Tk genomes", total_genomes],
            [
                f"Classified as {target_species}",
                target_count,
            ],
            [
                "Target-species percentage",
                target_percentage,
            ],
            [
                "Other or unresolved genomes",
                total_genomes - target_count,
            ],
            [
                "Median closest-reference ANI",
                dataframe["closest_genome_ani"].median(),
            ],
            [
                "Median alignment fraction",
                dataframe["closest_genome_af"].median(),
            ],
            [
                "Median ANI margin",
                dataframe["ani_margin"].median(),
            ],
            [
                "Median marker MSA coverage",
                dataframe["msa_percent"].median(),
            ],
            [
                "Genomes with GTDB-Tk warnings",
                int(dataframe["has_warning"].sum()),
            ],
        ],
        columns=["metric", "value"],
    )

    summary_metrics.to_csv(
        output_directory / "gtdbtk_figure_summary.tsv",
        sep="\t",
        index=False,
    )

    caption = (
        "Figure X. GTDB-Tk-based taxonomic confirmation of the "
        f"curated {target_species} genome collection. "
        "(A) Distribution of species-level classification outcomes. "
        "(B) Average nucleotide identity (ANI) and alignment fraction "
        "relative to the closest GTDB reference genome, with the median "
        "species-specific ANI radius indicated. "
        "(C) Distribution of the ANI margin calculated as closest-"
        "reference ANI minus the corresponding species-specific GTDB "
        "reference radius; values above zero support assignment within "
        "the reference species boundary. "
        "(D) Distribution of bacterial marker-gene multiple-sequence "
        "alignment coverage, with genomes carrying GTDB-Tk warnings "
        "indicated at the base of the distribution."
    )

    with open(
        output_directory / "gtdbtk_figure_caption.txt",
        "w",
        encoding="utf-8",
    ) as caption_file:
        caption_file.write(caption + "\n")

    print()
    print("GTDB-Tk composite publication figure completed")
    print("================================================")
    print(f"Total genomes:                 {total_genomes:,}")
    print(
        f"{target_species}: "
        f"{target_count:,} ({target_percentage:.2f}%)"
    )
    print(
        "Other or unresolved:          "
        f"{total_genomes - target_count:,}"
    )
    print(
        "Genomes with warnings:        "
        f"{int(dataframe['has_warning'].sum()):,}"
    )
    print()
    print("Composite figure files:")
    print(f"  {output_base}.png")
    print(f"  {output_base}.tiff")
    print(f"  {output_base}.pdf")
    print(f"  {output_base}.svg")


if __name__ == "__main__":
    main()
