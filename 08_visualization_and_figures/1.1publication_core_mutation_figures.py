#!/usr/bin/env python3

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt


STATUS_COLORS = {
    "Conserved-dominant": "#2ca02c",
    "Mixed": "#ffbf00",
    "Core-mutation hotspot": "#d62728",
}


STATUS_ORDER = [
    "Conserved-dominant",
    "Mixed",
    "Core-mutation hotspot",
]


def safe_read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str)


def as_numeric(df, cols):
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def clean_text(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"nan", "none", "unknown", "not available", "not collected", ""}:
        return ""
    return value


def shorten_label(label, max_len=45):
    label = clean_text(label)
    if len(label) <= max_len:
        return label
    return label[: max_len - 3] + "..."


def make_location_label(row):
    candidates = []

    for col in ["locations", "countries", "country_clean", "geo_loc_name"]:
        if col in row.index:
            val = clean_text(row[col])
            if val:
                candidates.append(val)

    if candidates:
        label = candidates[0]
    else:
        label = f"{row.get('latitude', '')}, {row.get('longitude', '')}"

    return shorten_label(label, 55)


def infer_score_column(location_df):
    candidates = [
        "mean_core_variant_count",
        "median_core_variant_count",
        "max_core_variant_count",
        "mean_core_variant_rate",
        "high_core_variant_fraction",
    ]

    for col in candidates:
        if col in location_df.columns:
            return col

    raise ValueError(
        "Could not find a usable mutation score column in location summary."
    )


def prepare_location_table(location_df):
    numeric_cols = [
        "latitude",
        "longitude",
        "genome_count",
        "conserved_count",
        "high_core_variant_count",
        "high_core_variant_fraction",
        "mean_core_variant_count",
        "median_core_variant_count",
        "max_core_variant_count",
        "mean_core_variant_rate",
    ]

    location_df = as_numeric(location_df.copy(), numeric_cols)

    location_df = location_df.dropna(subset=["latitude", "longitude"]).copy()

    if "region_core_status" not in location_df.columns:
        location_df["region_core_status"] = "Conserved-dominant"

    if "location_label" not in location_df.columns:
        location_df["location_label"] = location_df.apply(make_location_label, axis=1)

    return location_df


def prepare_merged_table(merged_df):
    numeric_cols = [
        "latitude",
        "longitude",
        "core_variant_count_vs_consensus",
        "core_variant_rate_vs_consensus",
        "core_compared_sites",
        "core_missing_sites",
    ]

    merged_df = as_numeric(merged_df.copy(), numeric_cols)

    if "country_clean" not in merged_df.columns:
        merged_df["country_clean"] = "Unknown"

    if "core_mutation_status" not in merged_df.columns:
        merged_df["core_mutation_status"] = "No mutation data"

    return merged_df


def compute_figure_metrics(location_df, merged_df):
    mapped_isolates = merged_df.dropna(subset=["latitude", "longitude"]).shape[0]

    mutation_mapped = merged_df[
        merged_df["core_variant_count_vs_consensus"].notna()
        & merged_df["latitude"].notna()
        & merged_df["longitude"].notna()
    ].shape[0]

    high_isolates = (
        merged_df["core_mutation_status"]
        .astype(str)
        .eq("High core-variant burden")
        .sum()
    )

    n_locations = location_df.shape[0]

    hotspot_locations = (
        location_df["region_core_status"]
        .astype(str)
        .eq("Core-mutation hotspot")
        .sum()
    )

    mixed_locations = (
        location_df["region_core_status"]
        .astype(str)
        .eq("Mixed")
        .sum()
    )

    conserved_locations = (
        location_df["region_core_status"]
        .astype(str)
        .eq("Conserved-dominant")
        .sum()
    )

    metrics = {
        "total_metadata_rows": merged_df.shape[0],
        "coordinate_resolved_isolates": mapped_isolates,
        "coordinate_and_mutation_resolved_isolates": mutation_mapped,
        "high_core_variant_isolates": int(high_isolates),
        "distinct_bubble_locations": int(n_locations),
        "conserved_dominant_locations": int(conserved_locations),
        "mixed_locations": int(mixed_locations),
        "hotspot_locations": int(hotspot_locations),
    }

    return metrics


def create_publication_bubble_map(
    location_df,
    merged_df,
    outdir,
    score_col,
    top_label_n=10,
    width=2200,
    height=1300,
):
    figdir = outdir / "figures"
    tabledir = outdir / "tables"

    metrics = compute_figure_metrics(location_df, merged_df)

    location_df = location_df.copy()
    location_df["plot_score"] = pd.to_numeric(location_df[score_col], errors="coerce")

    max_size = location_df["genome_count"].max()
    if pd.isna(max_size) or max_size <= 0:
        max_size = 1

    desired_max_marker = 55
    sizeref = 2.0 * max_size / (desired_max_marker ** 2)

    fig = go.Figure()

    for status in STATUS_ORDER:
        subset = location_df[location_df["region_core_status"] == status].copy()

        if subset.empty:
            continue

        fig.add_trace(
            go.Scattergeo(
                lon=subset["longitude"],
                lat=subset["latitude"],
                mode="markers",
                name=status,
                marker=dict(
                    size=subset["genome_count"],
                    sizemode="area",
                    sizeref=sizeref,
                    sizemin=5,
                    color=subset["plot_score"],
                    colorscale="Turbo",
                    cmin=location_df["plot_score"].quantile(0.02),
                    cmax=location_df["plot_score"].quantile(0.98),
                    colorbar=dict(
                        title="Mean core<br>variant burden",
                        len=0.65,
                        thickness=22,
                    ) if status == STATUS_ORDER[0] else None,
                    showscale=True if status == STATUS_ORDER[0] else False,
                    opacity=0.78,
                    line=dict(
                        width=2.2 if status == "Core-mutation hotspot" else 1.0,
                        color=STATUS_COLORS.get(status, "black"),
                    ),
                ),
                customdata=np.stack(
                    [
                        subset["location_label"],
                        subset["genome_count"],
                        subset.get("high_core_variant_count", pd.Series([""] * len(subset))),
                        subset.get("high_core_variant_fraction", pd.Series([""] * len(subset))),
                        subset.get("mean_core_variant_count", pd.Series([""] * len(subset))),
                        subset.get("median_core_variant_count", pd.Series([""] * len(subset))),
                        subset.get("max_core_variant_count", pd.Series([""] * len(subset))),
                        subset.get("countries", pd.Series([""] * len(subset))),
                    ],
                    axis=-1,
                ),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Isolates: %{customdata[1]}<br>"
                    "High-burden isolates: %{customdata[2]}<br>"
                    "High-burden fraction: %{customdata[3]}<br>"
                    "Mean variant count: %{customdata[4]}<br>"
                    "Median variant count: %{customdata[5]}<br>"
                    "Max variant count: %{customdata[6]}<br>"
                    "Countries: %{customdata[7]}<br>"
                    "<extra></extra>"
                ),
            )
        )

    # Label top hotspots / high-burden locations
    label_df = location_df.copy()
    label_df["is_hotspot"] = (
        label_df["region_core_status"] == "Core-mutation hotspot"
    ).astype(int)

    label_df = label_df.sort_values(
        [
            "is_hotspot",
            "high_core_variant_count",
            "high_core_variant_fraction",
            "genome_count",
            "plot_score",
        ],
        ascending=[False, False, False, False, False],
    ).head(top_label_n)

    fig.add_trace(
        go.Scattergeo(
            lon=label_df["longitude"],
            lat=label_df["latitude"],
            mode="text",
            text=label_df["location_label"].apply(lambda x: shorten_label(x, 22)),
            textfont=dict(size=16, color="black"),
            textposition="top center",
            showlegend=False,
            hoverinfo="skip",
        )
    )

    annotation_text = (
        f"<b>Clinical core-genome mutation map</b><br>"
        f"Coordinate-resolved isolates: {metrics['coordinate_and_mutation_resolved_isolates']}<br>"
        f"Distinct locations: {metrics['distinct_bubble_locations']}<br>"
        f"High-burden isolates: {metrics['high_core_variant_isolates']}<br>"
        f"Hotspot locations: {metrics['hotspot_locations']}<br>"
        f"Bubble size = isolates per location<br>"
        f"Bubble fill = mean core-variant burden<br>"
        f"Bubble outline = location status"
    )

    fig.update_layout(
        title=dict(
            text=(
                "Geographic Distribution of Clinical Isolates and "
                "Core-Genome Variant Burden"
            ),
            x=0.5,
            y=0.96,
            font=dict(size=30),
        ),
        width=width,
        height=height,
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend=dict(
            title="Location classification",
            orientation="h",
            x=0.5,
            y=0.03,
            xanchor="center",
            yanchor="bottom",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="rgba(0,0,0,0.25)",
            borderwidth=1,
            font=dict(size=18),
        ),
        margin=dict(l=20, r=20, t=90, b=80),
        annotations=[
            dict(
                x=0.02,
                y=0.07,
                xref="paper",
                yref="paper",
                text=annotation_text,
                align="left",
                showarrow=False,
                bgcolor="rgba(255,255,255,0.88)",
                bordercolor="rgba(0,0,0,0.35)",
                borderwidth=1,
                font=dict(size=16),
            )
        ],
    )

    fig.update_geos(
        projection_type="natural earth",
        showland=True,
        landcolor="rgb(240,240,240)",
        showocean=True,
        oceancolor="rgb(250,250,250)",
        showcountries=True,
        countrycolor="rgb(150,150,150)",
        countrywidth=0.7,
        showcoastlines=True,
        coastlinecolor="rgb(120,120,120)",
        coastlinewidth=0.8,
        showframe=False,
    )

    html_file = figdir / "Fig1_world_core_mutation_bubble_map.html"
    png_file = figdir / "Fig1_world_core_mutation_bubble_map.png"
    svg_file = figdir / "Fig1_world_core_mutation_bubble_map.svg"
    pdf_file = figdir / "Fig1_world_core_mutation_bubble_map.pdf"

    fig.write_html(html_file)

    for path in [png_file, svg_file, pdf_file]:
        try:
            fig.write_image(str(path), width=width, height=height, scale=3)
            print(f"[INFO] Saved: {path}")
        except Exception as e:
            print(f"[WARNING] Could not save {path.name}: {e}")

    label_table = figdir / "Fig1_labeled_locations.tsv"
    label_df.to_csv(label_table, sep="\t", index=False)

    metrics_table = tabledir / "figure_metrics_summary.tsv"
    pd.DataFrame([metrics]).to_csv(metrics_table, sep="\t", index=False)

    print(f"[INFO] Saved HTML: {html_file}")
    print(f"[INFO] Saved label table: {label_table}")
    print(f"[INFO] Saved metrics table: {metrics_table}")


def top_category_summary(series, top_n=5):
    if series is None:
        return pd.Series(dtype=int)
    cleaned = series.dropna().astype(str)
    cleaned = cleaned[~cleaned.str.lower().isin(["nan", "none", "unknown", "not available", ""])]
    return cleaned.value_counts().head(top_n)


def create_summary_panels(location_df, merged_df, outdir, top_n=12):
    figdir = outdir / "figures"
    tabledir = outdir / "tables"

    loc = location_df.copy()

    loc["location_label_short"] = loc["location_label"].apply(lambda x: shorten_label(x, 30))

    hotspot_rank = loc.sort_values(
        [
            "high_core_variant_count",
            "high_core_variant_fraction",
            "genome_count",
            "mean_core_variant_count",
        ],
        ascending=False,
    ).head(top_n)

    hotspot_rank.to_csv(
        tabledir / "top_hotspot_locations_for_publication.tsv",
        sep="\t",
        index=False,
    )

    country = create_country_summary(merged_df, tabledir)

    top_country = country.head(top_n).copy()

    iso = merged_df.copy()
    iso = iso[pd.to_numeric(iso["core_variant_count_vs_consensus"], errors="coerce").notna()].copy()

    high_iso = iso[iso["core_mutation_status"] == "High core-variant burden"].copy()

    source_counts = None
    for col in ["clinical_source", "sample_type", "isolation_source", "source_type"]:
        if col in high_iso.columns:
            source_counts = top_category_summary(high_iso[col], top_n=8)
            source_col_used = col
            if not source_counts.empty:
                break

    fig = plt.figure(figsize=(18, 14), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[0, 0])
    y = np.arange(len(hotspot_rank))[::-1]
    ax1.barh(y, hotspot_rank["high_core_variant_count"].values[::-1])
    ax1.set_yticks(y)
    ax1.set_yticklabels(hotspot_rank["location_label_short"].values[::-1], fontsize=9)
    ax1.set_xlabel("High core-variant isolates")
    ax1.set_title("A. Top locations by high-burden isolate count", loc="left", fontweight="bold")
    ax1.grid(axis="x", alpha=0.25)

    ax2 = fig.add_subplot(gs[0, 1])
    y2 = np.arange(len(top_country))[::-1]
    ax2.barh(y2, top_country["genome_count"].values[::-1])
    ax2.set_yticks(y2)
    ax2.set_yticklabels(top_country["country_clean"].values[::-1], fontsize=9)
    ax2.set_xlabel("Clinical isolates with mutation data")
    ax2.set_title("B. Top countries by mapped isolate count", loc="left", fontweight="bold")
    ax2.grid(axis="x", alpha=0.25)

    ax3 = fig.add_subplot(gs[1, 0])

    status_values = [
        "Conserved / lower core-variant burden",
        "High core-variant burden",
    ]

    for status in status_values:
        vals = iso.loc[
            iso["core_mutation_status"] == status,
            "core_variant_count_vs_consensus"
        ].dropna()

        if len(vals) > 0:
            ax3.hist(vals, bins=40, alpha=0.65, label=status)

    ax3.set_xlabel("Core variant count vs consensus")
    ax3.set_ylabel("Number of isolates")
    ax3.set_title("C. Per-isolate core variant burden distribution", loc="left", fontweight="bold")
    ax3.legend(fontsize=8)
    ax3.grid(axis="y", alpha=0.25)

    ax4 = fig.add_subplot(gs[1, 1])

    if source_counts is not None and not source_counts.empty:
        labels = [shorten_label(x, 28) for x in source_counts.index.tolist()]
        y4 = np.arange(len(labels))[::-1]
        ax4.barh(y4, source_counts.values[::-1])
        ax4.set_yticks(y4)
        ax4.set_yticklabels(labels[::-1], fontsize=9)
        ax4.set_xlabel("High-burden isolates")
        ax4.set_title(
            f"D. High-burden isolates by {source_col_used}",
            loc="left",
            fontweight="bold",
        )
        ax4.grid(axis="x", alpha=0.25)
    else:
        ax4.text(
            0.5,
            0.5,
            "No usable clinical/source category found",
            ha="center",
            va="center",
            fontsize=12,
        )
        ax4.set_axis_off()
        ax4.set_title("D. Source composition", loc="left", fontweight="bold")

    fig.suptitle(
        "Summary of Geographic Core-Genome Variant Burden in Clinical Isolates",
        fontsize=18,
        fontweight="bold",
        y=0.98,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    for ext in ["png", "svg", "pdf"]:
        path = figdir / f"Fig2_core_mutation_summary_panels.{ext}"
        fig.savefig(path, bbox_inches="tight")
        print(f"[INFO] Saved: {path}")

    plt.close(fig)


def create_country_summary(merged_df, tabledir):
    df = merged_df.copy()

    if "country_clean" not in df.columns:
        df["country_clean"] = "Unknown"

    df["country_clean"] = df["country_clean"].fillna("Unknown")

    df = df[
        df["core_variant_count_vs_consensus"].notna()
        & df["latitude"].notna()
        & df["longitude"].notna()
    ].copy()

    df["is_high_core_variant"] = (
        df["core_mutation_status"] == "High core-variant burden"
    )

    grouped = (
        df.groupby("country_clean", dropna=False)
        .agg(
            genome_count=("country_clean", "size"),
            high_core_variant_count=("is_high_core_variant", "sum"),
            mean_core_variant_count=("core_variant_count_vs_consensus", "mean"),
            median_core_variant_count=("core_variant_count_vs_consensus", "median"),
            max_core_variant_count=("core_variant_count_vs_consensus", "max"),
            mean_core_variant_rate=("core_variant_rate_vs_consensus", "mean"),
        )
        .reset_index()
    )

    grouped["high_core_variant_fraction"] = (
        grouped["high_core_variant_count"] / grouped["genome_count"]
    )

    grouped = grouped.sort_values(
        ["genome_count", "high_core_variant_count"],
        ascending=False,
    )

    outfile = tabledir / "country_level_core_mutation_summary.tsv"
    grouped.to_csv(outfile, sep="\t", index=False)
    print(f"[INFO] Saved: {outfile}")

    return grouped


def create_country_level_figure(merged_df, outdir, top_n=20):
    figdir = outdir / "figures"
    tabledir = outdir / "tables"

    country = create_country_summary(merged_df, tabledir)
    top = country.head(top_n).copy()

    fig, ax1 = plt.subplots(figsize=(16, 8), dpi=300)

    x = np.arange(len(top))
    bars = ax1.bar(x, top["genome_count"])
    ax1.set_ylabel("Clinical isolates with coordinate + mutation data")
    ax1.set_xticks(x)
    ax1.set_xticklabels(top["country_clean"], rotation=45, ha="right")
    ax1.grid(axis="y", alpha=0.25)

    ax2 = ax1.twinx()
    ax2.plot(
        x,
        top["mean_core_variant_count"],
        marker="o",
        linewidth=2,
    )
    ax2.set_ylabel("Mean core variant count vs consensus")

    ax1.set_title(
        "Country-level Sampling Density and Mean Core-Genome Variant Burden",
        fontweight="bold",
    )

    caption = (
        "Bars show isolate counts per country; line shows mean core-genome "
        "variant burden among coordinate-resolved clinical isolates."
    )

    ax1.text(
        0.01,
        -0.27,
        textwrap.fill(caption, 135),
        transform=ax1.transAxes,
        fontsize=10,
        va="top",
    )

    fig.tight_layout()

    for ext in ["png", "svg", "pdf"]:
        path = figdir / f"Fig3_country_level_core_mutation_summary.{ext}"
        fig.savefig(path, bbox_inches="tight")
        print(f"[INFO] Saved: {path}")

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Generate publication-ready core mutation geographic figures."
    )

    parser.add_argument(
        "--results-dir",
        required=True,
        help="Path to clinical_core_mutation_bubble_results directory."
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="publication_core_mutation_figures",
        help="Output directory for publication figures."
    )

    parser.add_argument(
        "--top-label-n",
        type=int,
        default=10,
        help="Number of top locations to label on the world map."
    )

    parser.add_argument(
        "--top-panel-n",
        type=int,
        default=12,
        help="Number of locations/countries shown in summary panels."
    )

    parser.add_argument(
        "--score-column",
        default=None,
        help="Location-level score column for map color. Default: auto."
    )

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    tables_dir = results_dir / "tables"

    location_file = tables_dir / "core_mutation_bubble_location_summary.tsv"
    merged_file = tables_dir / "metadata_with_core_variant_burden.tsv"

    if not location_file.exists():
        raise FileNotFoundError(f"Missing file: {location_file}")

    if not merged_file.exists():
        raise FileNotFoundError(f"Missing file: {merged_file}")

    outdir = Path(args.outdir)

    for subdir in ["figures", "tables", "logs"]:
        (outdir / subdir).mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading location summary: {location_file}")
    location_df = safe_read_tsv(location_file)
    location_df = prepare_location_table(location_df)

    print(f"[INFO] Reading merged metadata: {merged_file}")
    merged_df = safe_read_tsv(merged_file)
    merged_df = prepare_merged_table(merged_df)

    if args.score_column is None:
        score_col = infer_score_column(location_df)
    else:
        score_col = args.score_column

    if score_col not in location_df.columns:
        raise ValueError(f"Score column not found in location table: {score_col}")

    print(f"[INFO] Using map color score column: {score_col}")

    location_df.to_csv(
        outdir / "tables" / "publication_location_table_used.tsv",
        sep="\t",
        index=False,
    )

    merged_df.to_csv(
        outdir / "tables" / "publication_merged_metadata_used.tsv",
        sep="\t",
        index=False,
    )

    create_publication_bubble_map(
        location_df=location_df,
        merged_df=merged_df,
        outdir=outdir,
        score_col=score_col,
        top_label_n=args.top_label_n,
    )

    create_summary_panels(
        location_df=location_df,
        merged_df=merged_df,
        outdir=outdir,
        top_n=args.top_panel_n,
    )

    create_country_level_figure(
        merged_df=merged_df,
        outdir=outdir,
        top_n=20,
    )

    print()
    print("[DONE] Publication-ready mutation figures generated.")
    print(f"[DONE] Output directory: {outdir}")


if __name__ == "__main__":
    main()
