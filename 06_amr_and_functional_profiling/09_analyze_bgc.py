#!/usr/bin/env python3
"""Construct coordinate consensus BGC loci, concordance, burdens, and domain families."""

from __future__ import annotations

import argparse
import hashlib
import sys
from collections import Counter, defaultdict

import pandas as pd

from analysis_common import ROOT, atomic_write_tsv, prevalence_band


N_EVALUABLE = 1350
RECIPROCAL_SHORT_THRESHOLD = 0.5


class UnionFind:
    def __init__(self, items: list[int]):
        self.parent = {x: x for x in items}

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def split_set(value: str) -> set[str]:
    return {x for x in str(value).split(";") if x}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    df = pd.read_csv(ROOT / "02_cleaned/BGC_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    df["BGC_start"] = pd.to_numeric(df["BGC_start"]).astype(int)
    df["BGC_end"] = pd.to_numeric(df["BGC_end"]).astype(int)
    records = df.to_dict("records")
    index_to_consensus: dict[int, str] = {}
    component_payloads: list[tuple[str, str, int, int, list[int]]] = []
    concordance_buffer: list[dict[str, object]] = []

    grouped_indices: dict[tuple[str, str], list[int]] = defaultdict(list)
    for idx, record in enumerate(records):
        grouped_indices[(record["genome_id"], record["contig_id"])].append(idx)
    for (genome, contig), indices in grouped_indices.items():
        indices.sort(key=lambda i: (int(records[i]["BGC_start"]), int(records[i]["BGC_end"])))
        uf = UnionFind(indices)
        active: list[int] = []
        for idx in indices:
            start = int(records[idx]["BGC_start"])
            active = [j for j in active if int(records[j]["BGC_end"]) >= start]
            for other in active:
                if records[other]["prediction_tool"] == records[idx]["prediction_tool"]:
                    continue
                left = max(start, int(records[other]["BGC_start"]))
                right = min(int(records[idx]["BGC_end"]), int(records[other]["BGC_end"]))
                intersection = max(0, right - left + 1)
                if not intersection:
                    continue
                len_a = int(records[idx]["BGC_end"]) - start + 1
                len_b = int(records[other]["BGC_end"]) - int(records[other]["BGC_start"]) + 1
                union = len_a + len_b - intersection
                reciprocal_short = intersection / min(len_a, len_b)
                if reciprocal_short >= RECIPROCAL_SHORT_THRESHOLD:
                    uf.union(idx, other)
                concordance_buffer.append(
                    {
                        "record_type": "overlapping_region_pair",
                        "genome_id": genome,
                        "contig_id": contig,
                        "region_A": records[other]["feature_id"],
                        "tool_A": records[other]["prediction_tool"],
                        "type_A": records[other]["normalized_BGC_types"],
                        "region_B": records[idx]["feature_id"],
                        "tool_B": records[idx]["prediction_tool"],
                        "type_B": records[idx]["normalized_BGC_types"],
                        "intersection_bp": intersection,
                        "union_bp": union,
                        "intersection_over_union": intersection / union,
                        "intersection_over_shorter": reciprocal_short,
                        "consensus_link_threshold": RECIPROCAL_SHORT_THRESHOLD,
                        "linked_into_same_consensus": "yes" if reciprocal_short >= RECIPROCAL_SHORT_THRESHOLD else "no",
                    }
                )
            active.append(idx)
        components: dict[int, list[int]] = defaultdict(list)
        for idx in indices:
            components[uf.find(idx)].append(idx)
        for members in components.values():
            component_payloads.append(
                (
                    genome,
                    contig,
                    min(int(records[i]["BGC_start"]) for i in members),
                    max(int(records[i]["BGC_end"]) for i in members),
                    members,
                )
            )

    component_payloads.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    consensus_rows: list[dict[str, object]] = []
    tool_priority = {"antiSMASH": 0, "GECCO": 1, "deepBGC": 2}
    for serial, (genome, contig, start, end, members) in enumerate(component_payloads, start=1):
        consensus_id = f"BGCC{serial:06d}"
        for idx in members:
            index_to_consensus[idx] = consensus_id
        tools = sorted({records[i]["prediction_tool"] for i in members})
        types = set().union(*(split_set(records[i]["normalized_BGC_types"]) for i in members))
        if len(types) > 1:
            types.discard("Unclassified")
        pfams = set().union(*(split_set(records[i]["PFAM_domain_set"]) for i in members))
        representative_idx = sorted(
            members,
            key=lambda i: (
                tool_priority.get(records[i]["prediction_tool"], 9),
                -(int(records[i]["BGC_end"]) - int(records[i]["BGC_start"])),
                records[i]["feature_id"],
            ),
        )[0]
        consensus_rows.append(
            {
                "consensus_BGC_id": consensus_id,
                "genome_id": genome,
                "contig_id": contig,
                "BGC_start": start,
                "BGC_end": end,
                "BGC_span": end - start + 1,
                "prediction_count": len(members),
                "prediction_tools": ";".join(tools),
                "cross_tool_support_count": len(tools),
                "confidence_class": "cross-tool consensus" if len(tools) >= 2 else "single-tool prediction",
                "normalized_BGC_types": ";".join(sorted(types)) if types else "Unclassified",
                "raw_product_classes": ";".join(sorted({records[i]["raw_product_class"] for i in members})),
                "PFAM_domain_set": ";".join(sorted(pfams)),
                "representative_region": records[representative_idx]["feature_id"],
                "representative_tool": records[representative_idx]["prediction_tool"],
                "member_region_ids": ";".join(records[i]["feature_id"] for i in members),
            }
        )

    feature_to_index = {record["feature_id"]: i for i, record in enumerate(records)}
    for row in concordance_buffer:
        a = feature_to_index[str(row["region_A"])]
        b = feature_to_index[str(row["region_B"])]
        row["consensus_BGC_A"] = index_to_consensus[a]
        row["consensus_BGC_B"] = index_to_consensus[b]
        row["same_consensus_BGC"] = "yes" if index_to_consensus[a] == index_to_consensus[b] else "no"
    atomic_write_tsv(
        ROOT / "04_statistics/BGC_tool_concordance.tsv",
        concordance_buffer,
        list(concordance_buffer[0]),
    )

    pair_summaries: list[dict[str, object]] = []
    for tools, group in pd.DataFrame(concordance_buffer).groupby(["tool_A", "tool_B"]):
        pair_summaries.append(
            {
                "tool_A": tools[0],
                "tool_B": tools[1],
                "overlapping_region_pairs": len(group),
                "pairs_meeting_consensus_threshold": (group["linked_into_same_consensus"] == "yes").sum(),
                "median_intersection_over_union": group["intersection_over_union"].median(),
                "median_intersection_over_shorter": group["intersection_over_shorter"].median(),
                "genomes_with_overlap": group["genome_id"].nunique(),
            }
        )
    atomic_write_tsv(
        ROOT / "04_statistics/BGC_tool_concordance_summary.tsv",
        pair_summaries,
        list(pair_summaries[0]),
    )

    # Conservative BGC similarity families: exact domain-content signatures plus
    # the normalized product type. These are not claimed as BiG-SCAPE GCFs.
    family_groups: dict[str, list[int]] = defaultdict(list)
    unresolved: list[int] = []
    for idx, row in enumerate(consensus_rows):
        domains = row["PFAM_domain_set"]
        if not domains:
            unresolved.append(idx)
            continue
        signature = f"{row['normalized_BGC_types']}|{domains}"
        family_groups[signature].append(idx)
    ordered_families = sorted(
        family_groups.items(),
        key=lambda item: (-len({consensus_rows[i]["genome_id"] for i in item[1]}), -len(item[1]), item[0]),
    )
    family_id_for_consensus: dict[str, str] = {}
    family_rows: list[dict[str, object]] = []
    representative_rows: list[dict[str, object]] = []
    similarity_rows: list[dict[str, object]] = []
    for number, (signature, indices) in enumerate(ordered_families, start=1):
        family_id = f"BGCF_DOMAIN_{number:05d}"
        rep_idx = sorted(
            indices,
            key=lambda i: (-int(consensus_rows[i]["cross_tool_support_count"]), -int(consensus_rows[i]["BGC_span"]), consensus_rows[i]["consensus_BGC_id"]),
        )[0]
        rep = consensus_rows[rep_idx]
        genomes = {consensus_rows[i]["genome_id"] for i in indices}
        for idx in indices:
            family_id_for_consensus[str(consensus_rows[idx]["consensus_BGC_id"])] = family_id
            similarity_rows.append(
                {
                    "family_id": family_id,
                    "representative_BGC": rep["consensus_BGC_id"],
                    "member_BGC": consensus_rows[idx]["consensus_BGC_id"],
                    "PFAM_Jaccard_similarity": 1.0,
                    "similarity_definition": "exact normalized type and exact PFAM-domain set",
                }
            )
        family_rows.append(
            {
                "family_id": family_id,
                "family_definition": "conservative exact product-type plus PFAM-domain signature",
                "number_BGCs": len(indices),
                "number_genomes": len(genomes),
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * len(genomes) / N_EVALUABLE,
                "predicted_product": rep["normalized_BGC_types"],
                "representative_cluster": rep["consensus_BGC_id"],
                "representative_region": rep["representative_region"],
                "PFAM_domain_signature": rep["PFAM_domain_set"],
                "within_family_similarity": 1.0,
                "cross_tool_supported_BGCs": sum(int(consensus_rows[i]["cross_tool_support_count"]) >= 2 for i in indices),
                "interpretation_scope": "domain-content similarity family; not a BiG-SCAPE gene-cluster family",
            }
        )
        representative_rows.append(
            {
                "family_id": family_id,
                "representative_cluster": rep["consensus_BGC_id"],
                "genome_id": rep["genome_id"],
                "coordinates": f"{rep['contig_id']}:{rep['BGC_start']}-{rep['BGC_end']}",
                "predicted_product": rep["normalized_BGC_types"],
                "PFAM_domain_signature": rep["PFAM_domain_set"],
            }
        )
    for serial, idx in enumerate(unresolved, start=1):
        cid = str(consensus_rows[idx]["consensus_BGC_id"])
        family_id_for_consensus[cid] = f"BGCF_UNRESOLVED_{serial:06d}"
    for row in consensus_rows:
        row["BGC_similarity_family"] = family_id_for_consensus[str(row["consensus_BGC_id"])]
    atomic_write_tsv(
        ROOT / "02_cleaned/BGC_consensus_loci.tsv",
        consensus_rows,
        list(consensus_rows[0]),
    )
    atomic_write_tsv(
        ROOT / "06_tables/BGC_gene_cluster_families.tsv",
        family_rows,
        list(family_rows[0]),
    )
    atomic_write_tsv(
        ROOT / "03_matrices/BGC_cluster_similarity.tsv",
        similarity_rows,
        list(similarity_rows[0]),
    )
    atomic_write_tsv(
        ROOT / "06_tables/BGC_cluster_representatives.tsv",
        representative_rows,
        list(representative_rows[0]),
    )

    type_rows: list[dict[str, object]] = []
    types = sorted(set().union(*(split_set(row["normalized_BGC_types"]) for row in consensus_rows)))
    for bgc_type in types:
        subset = [row for row in consensus_rows if bgc_type in split_set(row["normalized_BGC_types"])]
        genomes = {row["genome_id"] for row in subset}
        counts = Counter(row["genome_id"] for row in subset)
        type_rows.append(
            {
                "BGC_type": bgc_type,
                "n_genomes": len(genomes),
                "n_evaluable": N_EVALUABLE,
                "prevalence_percent": 100 * len(genomes) / N_EVALUABLE,
                "total_consensus_regions": len(subset),
                "mean_regions_per_positive_genome": sum(counts.values()) / len(counts),
                "cross_tool_supported_regions": sum(int(row["cross_tool_support_count"]) >= 2 for row in subset),
                "prevalence_band": prevalence_band(100 * len(genomes) / N_EVALUABLE),
            }
        )
    type_rows.sort(key=lambda r: (-int(r["n_genomes"]), -int(r["total_consensus_regions"]), str(r["BGC_type"])))
    atomic_write_tsv(ROOT / "06_tables/BGC_type_prevalence.tsv", type_rows, list(type_rows[0]))
    atomic_write_tsv(ROOT / "06_tables/BGC_top_prevalent_types.tsv", type_rows[:25], list(type_rows[0]))

    all_genomes = sorted(df["genome_id"].unique())
    burden_rows: list[dict[str, object]] = []
    for genome in all_genomes:
        loci = [row for row in consensus_rows if row["genome_id"] == genome]
        raw = df[df["genome_id"] == genome]
        burden_rows.append(
            {
                "genome_id": genome,
                "total_BGC_count": len(loci),
                "unique_BGC_types": len(set().union(*(split_set(row["normalized_BGC_types"]) for row in loci))),
                "total_BGC_span": sum(int(row["BGC_span"]) for row in loci),
                "biosynthetic_gene_count_if_available": "",
                "cross_tool_consensus_BGC_count": sum(int(row["cross_tool_support_count"]) >= 2 for row in loci),
                "raw_tool_region_predictions": len(raw),
                "antiSMASH_region_predictions": (raw["prediction_tool"] == "antiSMASH").sum(),
                "deepBGC_region_predictions": (raw["prediction_tool"] == "deepBGC").sum(),
                "GECCO_region_predictions": (raw["prediction_tool"] == "GECCO").sum(),
            }
        )
    atomic_write_tsv(ROOT / "03_matrices/BGC_burden_by_genome.tsv", burden_rows, list(burden_rows[0]))

    type_names = [row["BGC_type"] for row in type_rows]
    type_present = {(row["genome_id"], t) for row in consensus_rows for t in split_set(row["normalized_BGC_types"])}
    matrix = []
    for genome in all_genomes:
        item = {"genome_id": genome}
        item.update({name: int((genome, name) in type_present) for name in type_names})
        matrix.append(item)
    atomic_write_tsv(ROOT / "03_matrices/BGC_type_presence_absence.tsv", matrix, ["genome_id", *type_names])

    # A dense matrix of tens of thousands of singleton signatures would be
    # uninformative and wasteful. The association-ready family matrix follows
    # the predeclared >=1% prevalence rule (at least 14/1,350 genomes); the
    # complete family catalogue remains in BGC_gene_cluster_families.tsv.
    family_names = [row["family_id"] for row in family_rows if int(row["number_genomes"]) >= 14]
    family_present = {(row["genome_id"], row["BGC_similarity_family"]) for row in consensus_rows}
    family_matrix = []
    for genome in all_genomes:
        item = {"genome_id": genome}
        item.update({name: int((genome, name) in family_present) for name in family_names})
        family_matrix.append(item)
    atomic_write_tsv(ROOT / "03_matrices/BGC_family_presence_absence.tsv", family_matrix, ["genome_id", *family_names])
    atomic_write_tsv(
        ROOT / "04_statistics/BGC_family_matrix_filter.tsv",
        [
            {
                "minimum_genomes": 14,
                "minimum_prevalence_percent": 100 * 14 / N_EVALUABLE,
                "families_retained": len(family_names),
                "families_catalogued_total": len(family_rows),
                "rationale": "exclude extremely rare signatures from inferential association matrices",
            }
        ],
        ["minimum_genomes", "minimum_prevalence_percent", "families_retained", "families_catalogued_total", "rationale"],
    )
    print(
        f"BGC analysis: {len(records)} raw predictions -> {len(consensus_rows)} consensus loci; "
        f"{sum(int(r['cross_tool_support_count']) >= 2 for r in consensus_rows)} cross-tool-supported; "
        f"{len(family_rows)} conservative domain families"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
