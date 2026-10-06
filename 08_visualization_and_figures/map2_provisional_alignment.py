#!/usr/bin/env python3
"""Build an unmasked, alignment-coordinate Map 2 sensitivity analysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from core_utils import WorkflowError, fasta_records, load_config, write_tsv


BASE_CODE = {base: ord(base) for base in "ACGT"}
EARTH_RADIUS_KM = 6371.0088


def load_matrix(path: Path, selected: set[str]) -> tuple[list[str], np.ndarray]:
    records: list[tuple[str, bytes]] = []
    length = None
    for identifier, _description, sequence in fasta_records(path):
        if identifier not in selected:
            continue
        if length is None:
            length = len(sequence)
        elif len(sequence) != length:
            raise WorkflowError("Alignment sequence lengths differ")
        records.append((identifier, sequence))
    if {identifier for identifier, _seq in records} != selected or length is None:
        raise WorkflowError("Primary cohort does not match the Map 2 alignment")
    matrix = np.empty((len(records), length), dtype=np.uint8)
    for row, (_identifier, sequence) in enumerate(records):
        matrix[row] = np.frombuffer(sequence, dtype=np.uint8)
    return [identifier for identifier, _sequence in records], matrix


def max_haversine_km(latitudes: np.ndarray, longitudes: np.ndarray) -> float:
    if len(latitudes) < 2:
        return 0.0
    lat = np.radians(latitudes)
    lon = np.radians(longitudes)
    maximum = 0.0
    for start in range(0, len(lat), 256):
        lat1 = lat[start:start + 256, None]
        lon1 = lon[start:start + 256, None]
        dlat = lat1 - lat[None, :]
        dlon = lon1 - lon[None, :]
        a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
        distance = 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
        maximum = max(maximum, float(distance.max()))
    return maximum


@lru_cache(maxsize=None)
def continent(country: str) -> str:
    if not country:
        return ""
    try:
        import country_converter as coco
        value = coco.convert(names=[country], to="continent", not_found="__NOT_FOUND__")
        value = value[0] if isinstance(value, list) else str(value)
        return "" if value == "__NOT_FOUND__" else value
    except Exception:
        return ""


def make_small_multiples(top: pd.DataFrame, carriers: sparse.csc_matrix, isolate_ids: list[str],
                         geo: pd.DataFrame, output: Path) -> None:
    if top.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    id_to_row = {value: index for index, value in enumerate(isolate_ids)}
    geo = geo[geo["assembly_accession"].isin(id_to_row)].copy()
    geo["matrix_row"] = geo["assembly_accession"].map(id_to_row)
    location = (geo.groupby("location_cluster_id", as_index=False)
                .agg(longitude=("longitude", "mean"), latitude=("latitude", "mean"),
                     callable_isolates=("assembly_accession", "size")))
    n_panels = len(top)
    ncols = 2
    nrows = math.ceil(n_panels / ncols)
    fig = plt.figure(figsize=(13.4, 4.9 * nrows))
    max_carriers_at_location = 1
    panel_data: list[pd.DataFrame] = []
    for _, variant in top.iterrows():
        vector = carriers[:, int(variant["matrix_column_index"])].toarray().ravel().astype(bool)
        panel = geo.copy()
        panel["carrier"] = vector[panel["matrix_row"].to_numpy(int)]
        summary = (panel.groupby("location_cluster_id", as_index=False)
                   .agg(carriers=("carrier", "sum"), callable_isolates=("carrier", "size"),
                        longitude=("longitude", "mean"), latitude=("latitude", "mean")))
        summary = summary[summary["carriers"] > 0]
        summary["local_prevalence"] = summary["carriers"] / summary["callable_isolates"]
        if not summary.empty:
            max_carriers_at_location = max(max_carriers_at_location, int(summary["carriers"].max()))
        panel_data.append(summary)
    for panel_index, ((_, variant), summary) in enumerate(zip(top.iterrows(), panel_data), start=1):
        ax = fig.add_subplot(nrows, ncols, panel_index, projection=ccrs.Robinson())
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor="#f3f0e8", edgecolor="none")
        ax.add_feature(cfeature.OCEAN, facecolor="#f7fbfd", edgecolor="none")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.35, edgecolor="#777777")
        ax.add_feature(cfeature.BORDERS, linewidth=0.25, edgecolor="#aaaaaa")
        ax.scatter(location["longitude"], location["latitude"], s=6, color="#d0d0d0", alpha=0.65,
                   linewidths=0, transform=ccrs.PlateCarree())
        scatter = ax.scatter(summary["longitude"], summary["latitude"], s=16 * summary["carriers"],
                             c=summary["local_prevalence"], vmin=0, vmax=1, cmap="viridis",
                             edgecolor="#333333", linewidth=0.55, transform=ccrs.PlateCarree())
        ax.set_title(
            f"{variant['mutation_identifier']}\n"
            f"carriers={int(variant['carrier_count'])}; countries={int(variant['number_countries'])}; "
            f"locations={int(variant['number_locations'])}; gains=not estimated",
            fontsize=10,
        )
    for panel_index in range(n_panels + 1, nrows * ncols + 1):
        ax = fig.add_subplot(nrows, ncols, panel_index)
        ax.axis("off")
    colorbar = fig.colorbar(scatter, ax=fig.axes, orientation="horizontal", fraction=0.025, pad=0.025)
    colorbar.set_label("Local carrier prevalence among callable isolates")
    fig.suptitle("Provisional global distribution of recurrent alignment-coordinate core alleles\n"
                 "Unmasked sensitivity analysis; identifiers are not PAO1 chromosome coordinates", fontsize=14)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(output / f"PROVISIONAL_Map2_global_recurrent_alignment_variants.{suffix}",
                    dpi=600 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--eligible-sites", required=True)
    parser.add_argument("--cohort-ids", required=True)
    parser.add_argument("--geographic-metadata", required=True)
    parser.add_argument("--cluster-membership", required=True)
    parser.add_argument("--output-dir", default="results/map2/provisional_unfiltered_alignment_coordinates")
    args = parser.parse_args()
    config, project_root = load_config(args.config)
    output = (project_root / args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(args.cohort_ids, sep="\t")
    selected = set(cohort["assembly_accession"])
    isolate_ids, matrix = load_matrix(Path(args.alignment).resolve(), selected)
    sites = pd.read_csv(args.eligible_sites, sep="\t")
    if sites["recombination_status"].ne("unknown").any():
        raise WorkflowError("This script is restricted to the unmasked sensitivity analysis")
    positions = sites["alignment_position"].to_numpy(int) - 1
    major_codes = np.array([BASE_CODE[value] for value in sites["major_allele"]], dtype=np.uint8)
    minor_codes = np.array([BASE_CODE[value] for value in sites["minor_allele"]], dtype=np.uint8)
    selected_calls = matrix[:, positions]
    callable_matrix = np.isin(selected_calls, list(BASE_CODE.values()))
    if not callable_matrix.all():
        raise WorkflowError("Complete-site SNP input unexpectedly contains non-callable states")
    carrier_dense = selected_calls == minor_codes[np.newaxis, :]
    carrier_matrix = sparse.csc_matrix(carrier_dense)
    sparse.save_npz(output / "core_variant_isolate_matrix.npz", carrier_matrix.tocsr(), compressed=True)

    geo = pd.read_csv(args.geographic_metadata, sep="\t")
    membership = pd.read_csv(args.cluster_membership, sep="\t")
    geo = geo.merge(membership[["assembly_accession", "location_cluster_id"]], on="assembly_accession", how="inner", validate="one_to_one")
    geo_by_id = geo.set_index("assembly_accession")
    id_to_row = {identifier: row for row, identifier in enumerate(isolate_ids)}
    row_to_geo = {id_to_row[value]: geo_by_id.loc[value] for value in geo_by_id.index if value in id_to_row}
    metadata_rows: list[dict[str, Any]] = []
    recurrent_rows: list[dict[str, Any]] = []
    min_carriers = int(config["map2"]["minimum_carriers"])
    min_locations = int(config["map2"]["minimum_locations"])
    min_countries = int(config["map2"]["minimum_countries"])
    for column, site in enumerate(sites.itertuples(index=False)):
        carrier_indices = carrier_matrix[:, column].indices
        carrier_count = len(carrier_indices)
        mutation_id = f"ALN_{int(site.alignment_position):08d}_{site.major_allele}_{site.minor_allele}"
        base_row = {
            "matrix_column_index": column,
            "mutation_identifier": mutation_id,
            "alignment_position": int(site.alignment_position),
            "coordinate_system": "alignment",
            "baseline_allele": site.major_allele,
            "alternate_allele": site.minor_allele,
            "baseline_allele_definition": "deduplicated_cohort_majority_not_reference_or_ancestral",
            "carrier_count": carrier_count,
            "callable_isolate_count": len(isolate_ids),
            "global_carrier_frequency": carrier_count / len(isolate_ids),
            "recombination_status": "unknown_unfiltered",
            "analysis_status": "provisional_alignment_coordinate_sensitivity",
        }
        metadata_rows.append(base_row)
        if carrier_count < min_carriers:
            continue
        carrier_geo = [row_to_geo[index] for index in carrier_indices if index in row_to_geo]
        countries = sorted({str(row["country_standardized"]) for row in carrier_geo if pd.notna(row["country_standardized"]) and str(row["country_standardized"])})
        locations = sorted({str(row["location_cluster_id"]) for row in carrier_geo})
        if len(locations) < min_locations or len(countries) < min_countries:
            continue
        latitudes = np.array([float(row["latitude"]) for row in carrier_geo])
        longitudes = np.array([float(row["longitude"]) for row in carrier_geo])
        continents = {continent(country) for country in countries} - {""}
        recurrent_rows.append({
            **base_row,
            "coordinate_resolved_carrier_count": len(carrier_geo),
            "number_locations": len(locations),
            "number_countries": len(countries),
            "number_continents": len(continents),
            "maximum_pairwise_geographic_distance_km": max_haversine_km(latitudes, longitudes),
            "number_represented_lineages": "",
            "minimum_inferred_phylogenetic_gains": "",
            "phylogenetic_origin_status": "not_estimated_tree_requires_deduplicated_reconstruction",
        })
    metadata = pd.DataFrame(metadata_rows)
    recurrent = pd.DataFrame(recurrent_rows)
    metadata.to_csv(output / "core_variant_isolate_matrix_metadata.tsv", sep="\t", index=False)
    recurrent.to_csv(output / "globally_recurrent_core_variants.tsv", sep="\t", index=False)
    recurrent.to_csv(output / "variant_geographic_summary.tsv", sep="\t", index=False)

    recurrent_columns = set(recurrent["matrix_column_index"].astype(int)) if not recurrent.empty else set()
    if not recurrent.empty:
        long_rank = (np.log1p(recurrent["carrier_count"]) * np.log1p(recurrent["number_countries"]))
        callable_long_columns = set(recurrent.assign(_long_rank=long_rank)
                                    .nlargest(20, "_long_rank")["matrix_column_index"].astype(int))
    else:
        callable_long_columns = set()
    with gzip.open(output / "core_variant_carriers_long.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["mutation_identifier", "assembly_accession", "table_scope"])
        for column in sorted(recurrent_columns):
            mutation_id = metadata.iloc[column]["mutation_identifier"]
            for row in carrier_matrix[:, column].indices:
                writer.writerow([mutation_id, isolate_ids[row], "globally_recurrent_provisional_variants"])
    with gzip.open(output / "core_variant_callable_long.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["mutation_identifier", "assembly_accession", "callable", "table_scope"])
        for column in sorted(callable_long_columns):
            mutation_id = metadata.iloc[column]["mutation_identifier"]
            for isolate in isolate_ids:
                writer.writerow([mutation_id, isolate, 1, "top_20_geographically_ranked_provisional_variants"])

    location_rows: list[dict[str, Any]] = []
    country_rows: list[dict[str, Any]] = []
    prevalence_variants = recurrent[recurrent["matrix_column_index"].isin(callable_long_columns)]
    for variant in prevalence_variants.itertuples(index=False):
        vector = carrier_matrix[:, int(variant.matrix_column_index)].toarray().ravel().astype(bool)
        local = geo.copy()
        local["matrix_row"] = local["assembly_accession"].map(id_to_row)
        local = local.dropna(subset=["matrix_row"])
        local["carrier"] = vector[local["matrix_row"].astype(int).to_numpy()]
        for cluster_id, subset in local.groupby("location_cluster_id"):
            carriers_n = int(subset["carrier"].sum())
            if carriers_n:
                location_rows.append({"mutation_identifier": variant.mutation_identifier,
                                      "location_cluster_id": cluster_id,
                                      "carrier_count": carriers_n, "callable_isolate_count": len(subset),
                                      "local_carrier_prevalence": carriers_n / len(subset),
                                      "table_scope": "top_20_geographically_ranked_provisional_variants"})
        for country_name, subset in local.groupby("country_standardized"):
            carriers_n = int(subset["carrier"].sum())
            if carriers_n:
                country_rows.append({"mutation_identifier": variant.mutation_identifier,
                                     "country": country_name, "carrier_count": carriers_n,
                                     "callable_isolate_count": len(subset),
                                     "country_carrier_prevalence": carriers_n / len(subset),
                                     "table_scope": "top_20_geographically_ranked_provisional_variants"})
    write_tsv(output / "variant_location_prevalence.tsv", location_rows,
              ["mutation_identifier", "location_cluster_id", "carrier_count", "callable_isolate_count", "local_carrier_prevalence", "table_scope"])
    write_tsv(output / "variant_country_prevalence.tsv", country_rows,
              ["mutation_identifier", "country", "carrier_count", "callable_isolate_count", "country_carrier_prevalence", "table_scope"])

    if not recurrent.empty:
        recurrent["descriptive_geographic_ranking_score"] = np.log1p(recurrent["carrier_count"]) * np.log1p(recurrent["number_countries"])
        recurrent["ranking_score_definition"] = "log1p(carriers)*log1p(countries); origins omitted because not estimated"
        ranking = recurrent.sort_values(["descriptive_geographic_ranking_score", "mutation_identifier"], ascending=[False, True])
    else:
        ranking = recurrent.copy()
    ranking.to_csv(output / "map2_variant_ranking.tsv", sep="\t", index=False)
    ranking.head(20).to_csv(output / "map2_top_disseminated_variants.tsv", sep="\t", index=False)
    write_tsv(output / "map2_top_homoplasy_candidates.tsv", [], ["mutation_identifier", "status"])
    write_tsv(output / "variant_phylogenetic_origin_counts.tsv",
              ({"mutation_identifier": value, "status": "not_estimated_deduplicated_tree_required"}
               for value in recurrent.get("mutation_identifier", [])), ["mutation_identifier", "status"])
    write_tsv(output / "variant_homoplasy_qc.tsv",
              ({"mutation_identifier": value, "coordinate_mapping": "alignment_only", "recombination_status": "unknown",
                "artefact_review_status": "required"} for value in recurrent.get("mutation_identifier", [])),
              ["mutation_identifier", "coordinate_mapping", "recombination_status", "artefact_review_status"])
    write_tsv(output / "variant_functional_annotation.tsv",
              ({"mutation_identifier": value, "annotation_status": "blocked_no_validated_PAO1_mapping"}
               for value in recurrent.get("mutation_identifier", [])), ["mutation_identifier", "annotation_status"])
    write_tsv(output / "variant_lineage_distribution.tsv", [], ["mutation_identifier", "lineage", "carrier_count", "status"])
    write_tsv(output / "variant_study_platform_bias_checks.tsv",
              ({"mutation_identifier": value, "platform_status": "metadata_unavailable", "study_check_status": "pending"}
               for value in recurrent.get("mutation_identifier", [])),
              ["mutation_identifier", "platform_status", "study_check_status"])
    if not ranking.empty:
        top = ranking.head(6)
        top.to_csv(output / "PROVISIONAL_Map2_source_data.tsv", sep="\t", index=False)
        make_small_multiples(top, carrier_matrix, isolate_ids, geo, output)
    print(f"Eligible biallelic alignment variants: {len(sites)}")
    print(f"Provisional globally recurrent variants: {len(recurrent)}")
    print("Phylogenetic gain counts and PAO1 annotations remain blocked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
