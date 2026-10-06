#!/usr/bin/env python3
"""Create transparent, category-specific candidate rankings."""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from analysis_common import ROOT, atomic_write_tsv


def minmax(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").fillna(0.0)
    if numeric.max() == numeric.min():
        return pd.Series(np.zeros(len(numeric)), index=series.index)
    return (numeric - numeric.min()) / (numeric.max() - numeric.min())


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    nodes = pd.read_csv(ROOT / "06_tables/functional_network_nodes.tsv", sep="\t", dtype=str, keep_default_na=False)
    centrality = dict(zip(nodes["node"], pd.to_numeric(nodes["weighted_degree"], errors="coerce").fillna(0.0)))

    amp_summary = pd.read_csv(ROOT / "06_tables/AMP_cluster_summary.tsv", sep="\t")
    amp_summary = amp_summary[amp_summary["identity_threshold_percent"] == 50].copy()
    amp_assign = pd.read_csv(ROOT / "03_matrices/AMP_sequence_clusters.tsv", sep="\t", dtype=str, keep_default_na=False)
    amp_assign = amp_assign[amp_assign["identity_threshold_percent"] == "50"][["sequence_id", "cluster_id"]]
    amp_master = pd.read_csv(
        ROOT / "02_cleaned/AMP_master_long.tsv",
        sep="\t",
        usecols=["sequence_id", "sequence_length_aa", "hydrophobicity", "isoelectric_point", "molecular_weight"],
        low_memory=False,
    )
    amp_properties = amp_assign.merge(
        amp_master[["sequence_id", "sequence_length_aa", "hydrophobicity", "isoelectric_point", "molecular_weight"]],
        on="sequence_id",
        how="left",
    ).groupby("cluster_id").median(numeric_only=True)
    amp_summary = amp_summary.merge(amp_properties, left_on="cluster_id", right_index=True, how="left")
    amp_summary["prevalence_component"] = amp_summary["prevalence_percent"] / 100
    amp_summary["predictor_consensus_component"] = amp_summary["high_confidence_predictions"] / amp_summary["number_of_sequences"]
    amp_summary["sequence_conservation_component"] = pd.to_numeric(
        amp_summary["mean_identity_to_representative_percent"], errors="coerce"
    ).fillna(50) / 100
    amp_summary["database_support_component"] = amp_summary["database_supported_predictions"] / amp_summary["number_of_sequences"]
    amp_summary["mapping_component"] = amp_summary["Prokka_Panaroo_mapped_predictions"] / amp_summary["number_of_sequences"]
    amp_summary["priority_score"] = 100 * (
        0.25 * amp_summary["prevalence_component"]
        + 0.30 * amp_summary["predictor_consensus_component"]
        + 0.20 * amp_summary["sequence_conservation_component"]
        + 0.15 * amp_summary["database_support_component"]
        + 0.10 * amp_summary["mapping_component"]
    )
    amp_summary["physicochemical_scope"] = "reported for review but not scored because no universal favorable AMP profile was imposed"
    amp_summary["score_formula"] = "25% prevalence + 30% multi-predictor confidence + 20% sequence conservation + 15% database support + 10% defensible gene mapping"
    amp_summary["validation_scope"] = "computational AMP candidates; no experimental validation or curated database match in this dataset"
    amp_summary = amp_summary.sort_values(["priority_score", "prevalence_percent", "cluster_id"], ascending=[False, False, True])
    amp_summary.insert(0, "rank", range(1, len(amp_summary) + 1))
    atomic_write_tsv(ROOT / "06_tables/AMP_priority_candidates.tsv", amp_summary.to_dict("records"), list(amp_summary.columns))

    amr = pd.read_csv(ROOT / "06_tables/AMR_family_prevalence_all.tsv", sep="\t")
    amr["prevalence_informativeness_component"] = 1 - (amr["prevalence_percent"] - 50).abs() / 50
    amr["prevalence_informativeness_component"] = amr["prevalence_informativeness_component"].clip(lower=0)
    amr["network_weighted_degree"] = amr["gene_family"].map(centrality).fillna(0)
    amr["network_hub_component"] = minmax(amr["network_weighted_degree"])
    amr["copy_variation_component"] = minmax(amr["max_copies_per_genome"])
    amr["multi_tool_component"] = amr["prediction_tools"].fillna("").map(lambda text: min(1.0, len([x for x in text.split(";") if x]) / 4))
    amr["priority_score"] = 100 * (
        0.40 * amr["prevalence_informativeness_component"]
        + 0.30 * amr["network_hub_component"]
        + 0.15 * amr["copy_variation_component"]
        + 0.15 * amr["multi_tool_component"]
    )
    amr["score_formula"] = "40% prevalence informativeness + 30% stringent-network hub score + 15% copy variation + 15% multi-tool support"
    amr["intrinsic_acquired_scope"] = "ranking favors discriminatory distributions but does not classify a family as acquired or confirm resistance phenotype"
    amr = amr.sort_values(["priority_score", "n_genomes", "gene_family"], ascending=[False, False, True])
    amr.insert(0, "rank", range(1, len(amr) + 1))
    atomic_write_tsv(ROOT / "06_tables/AMR_priority_candidates.tsv", amr.to_dict("records"), list(amr.columns))

    bgc = pd.read_csv(ROOT / "06_tables/BGC_gene_cluster_families.tsv", sep="\t")
    bgc["network_node"] = "BGC_FAMILY::" + bgc["family_id"].astype(str)
    bgc["prevalence_component"] = bgc["prevalence_percent"] / 100
    bgc["cross_tool_component"] = (bgc["cross_tool_supported_BGCs"] / bgc["number_BGCs"]).clip(upper=1)
    bgc["domain_similarity_component"] = pd.to_numeric(bgc["within_family_similarity"], errors="coerce").fillna(0)
    bgc["network_weighted_degree"] = bgc["network_node"].map(centrality).fillna(0)
    bgc["network_hub_component"] = minmax(bgc["network_weighted_degree"])
    bgc["priority_score"] = 100 * (
        0.30 * bgc["prevalence_component"]
        + 0.30 * bgc["cross_tool_component"]
        + 0.20 * bgc["domain_similarity_component"]
        + 0.20 * bgc["network_hub_component"]
    )
    bgc["score_formula"] = "30% prevalence + 30% cross-tool support + 20% exact domain-signature similarity + 20% stringent-network hub score"
    bgc["family_scope"] = "conservative product/PFAM signature family, not a BiG-SCAPE gene-cluster family and not evidence of metabolite production"
    bgc = bgc.sort_values(["priority_score", "number_genomes", "family_id"], ascending=[False, False, True])
    bgc.insert(0, "rank", range(1, len(bgc) + 1))
    atomic_write_tsv(ROOT / "06_tables/BGC_priority_candidates.tsv", bgc.to_dict("records"), list(bgc.columns))

    criteria = [
        {
            "category": "AMP",
            "criteria": amp_summary.iloc[0]["score_formula"],
            "rationale": "prioritize distributed, conserved, multi-predictor candidates while retaining database and mapping evidence",
            "important_exclusion": "physicochemical values are reported but not assigned a universal desirability score",
        },
        {
            "category": "AMR",
            "criteria": amr.iloc[0]["score_formula"],
            "rationale": "prioritize variable, clinically informative genotypic patterns and association hubs",
            "important_exclusion": "no phenotypic MDR or acquired/intrinsic classification is inferred from score",
        },
        {
            "category": "BGC",
            "criteria": bgc.iloc[0]["score_formula"],
            "rationale": "prioritize recurrent, cross-tool-supported domain architectures with integrated association evidence",
            "important_exclusion": "product classes and domain signatures are not confirmed GCFs or produced metabolites",
        },
    ]
    atomic_write_tsv(ROOT / "06_tables/candidate_priority_scoring_criteria.tsv", criteria, list(criteria[0]))
    print(
        f"Candidate ranking: AMP={len(amp_summary)}, AMR={len(amr)}, BGC={len(bgc)}; "
        f"top AMP={amp_summary.iloc[0]['cluster_id']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
