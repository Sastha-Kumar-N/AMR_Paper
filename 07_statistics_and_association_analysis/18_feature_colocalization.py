#!/usr/bin/env python3
"""Evaluate physical proximity for statistically strong feature associations."""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

from analysis_common import ROOT, atomic_write_tsv


Coordinate = tuple[str, int, int, str]


def add_coordinate(store, category: str, feature: str, genome: str, contig: str, start, end, identifier: str):
    if not feature or not genome or not contig or str(start) == "" or str(end) == "":
        return
    key = (contig, int(float(start)), int(float(end)), identifier)
    store[category][feature][genome].add(key)


def minimum_same_contig_distance(left: set[Coordinate], right: set[Coordinate]):
    best = None
    best_pair = None
    for contig_1, start_1, end_1, id_1 in left:
        for contig_2, start_2, end_2, id_2 in right:
            if contig_1 != contig_2:
                continue
            if end_1 < start_2:
                distance = start_2 - end_1 - 1
            elif end_2 < start_1:
                distance = start_1 - end_2 - 1
            else:
                distance = 0
            if best is None or distance < best:
                best = distance
                best_pair = f"{contig_1}:{start_1}-{end_1}({id_1})|{start_2}-{end_2}({id_2})"
    return best, best_pair


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    coordinates = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

    amp_assign = pd.read_csv(ROOT / "03_matrices/AMP_sequence_clusters.tsv", sep="\t", dtype=str, keep_default_na=False)
    amp_assign = amp_assign[amp_assign["identity_threshold_percent"] == "50"]
    feature_to_cluster = dict(zip(amp_assign["feature_id"], amp_assign["cluster_id"]))
    amp = pd.read_csv(ROOT / "02_cleaned/AMP_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    for row in amp.itertuples(index=False):
        cluster = feature_to_cluster.get(row.feature_id, "")
        add_coordinate(coordinates, "AMP", cluster, row.genome_id, row.contig_id, row.CDS_start, row.CDS_end, row.feature_id)

    amr = pd.read_csv(ROOT / "02_cleaned/AMR_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    amr = amr[(amr["primary_analysis"] == "yes") & (amr["contig_id"] != "") & (amr["coordinate_start"] != "") & (amr["coordinate_end"] != "")]
    for row in amr.itertuples(index=False):
        add_coordinate(
            coordinates,
            "AMR",
            row.gene_family,
            row.genome_id,
            row.contig_id,
            row.coordinate_start,
            row.coordinate_end,
            row.feature_id,
        )

    bgc = pd.read_csv(ROOT / "02_cleaned/BGC_consensus_loci.tsv", sep="\t", dtype=str, keep_default_na=False)
    for row in bgc.itertuples(index=False):
        family = f"BGC_FAMILY::{row.BGC_similarity_family}"
        add_coordinate(coordinates, "BGC", family, row.genome_id, row.contig_id, row.BGC_start, row.BGC_end, row.consensus_BGC_id)
        for bgc_type in str(row.normalized_BGC_types).split(";"):
            if bgc_type:
                add_coordinate(
                    coordinates,
                    "BGC",
                    f"BGC_TYPE::{bgc_type}",
                    row.genome_id,
                    row.contig_id,
                    row.BGC_start,
                    row.BGC_end,
                    row.consensus_BGC_id,
                )

    inputs = [
        "AMP_ARG_associations.tsv",
        "AMP_BGC_associations.tsv",
        "ARG_BGC_associations.tsv",
        "AMR_gene_cooccurrence.tsv",
    ]
    output = []
    for filename in inputs:
        associations = pd.read_csv(ROOT / f"04_statistics/{filename}", sep="\t", dtype=str, keep_default_na=False)
        associations["abs_phi"] = pd.to_numeric(associations["Phi_coefficient"]).abs()
        selected = associations[associations["network_eligible"] == "yes"].sort_values("abs_phi", ascending=False)
        for rank, row in enumerate(selected.itertuples(index=False), start=1):
            left_by_genome = coordinates[row.category_1].get(row.feature_1, {})
            right_by_genome = coordinates[row.category_2].get(row.feature_2, {})
            shared_mappable = sorted(set(left_by_genome) & set(right_by_genome))
            distances = []
            examples = []
            for genome in shared_mappable:
                distance, pair = minimum_same_contig_distance(left_by_genome[genome], right_by_genome[genome])
                if distance is not None:
                    distances.append(distance)
                    examples.append((distance, f"{genome}:{distance}bp:{pair}"))
            examples.sort(key=lambda item: (item[0], item[1]))
            within_10kb = sum(distance <= 10000 for distance in distances)
            overlaps = sum(distance == 0 for distance in distances)
            output.append(
                {
                    "association_source": filename,
                    "source_rank_by_absolute_phi": rank,
                    "feature_1": row.feature_1,
                    "category_1": row.category_1,
                    "feature_2": row.feature_2,
                    "category_2": row.category_2,
                    "n_genomes_cooccurring": int(row.n_both),
                    "Phi_coefficient": row.Phi_coefficient,
                    "q_value": row.q_value,
                    "relationship": row.relationship,
                    "n_cooccurring_genomes_with_both_coordinate_sets": len(shared_mappable),
                    "n_genomes_same_contig": len(distances),
                    "percent_mappable_cooccurrences_same_contig": 100 * len(distances) / len(shared_mappable) if shared_mappable else "",
                    "n_genomes_distance_le_10kb": within_10kb,
                    "percent_mappable_cooccurrences_distance_le_10kb": 100 * within_10kb / len(shared_mappable) if shared_mappable else "",
                    "n_genomes_features_overlap": overlaps,
                    "minimum_same_contig_distance_bp": min(distances) if distances else "",
                    "median_same_contig_distance_bp": float(np.median(distances)) if distances else "",
                    "closest_examples": ";".join(text for _, text in examples[:5]),
                    "physical_linkage_interpretation": (
                        "frequent <=10 kb proximity among coordinate-mappable cooccurrences"
                        if shared_mappable and within_10kb / len(shared_mappable) >= 0.5
                        else "limited or no frequent <=10 kb proximity"
                        if shared_mappable
                        else "not evaluable from available coordinates"
                    ),
                    "scope_note": "BGC coordinates denote predicted regions; overlap does not prove biosynthetic function or direct interaction",
                }
            )
    output.sort(key=lambda row: (str(row["association_source"]), int(row["source_rank_by_absolute_phi"])))
    atomic_write_tsv(ROOT / "04_statistics/feature_colocalization.tsv", output, list(output[0]))
    print(
        f"Colocalization: {len(output)} q<0.05, |phi|>=0.20 associations evaluated; "
        f"{sum(row['physical_linkage_interpretation'].startswith('frequent') for row in output)} show frequent <=10 kb proximity"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
