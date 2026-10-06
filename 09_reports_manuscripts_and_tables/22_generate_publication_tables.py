#!/usr/bin/env python3
"""Assemble publication tables from traceable primary analysis outputs."""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from analysis_common import ROOT, atomic_write_tsv


def write_frame(path, frame):
    atomic_write_tsv(path, frame.to_dict("records"), list(frame.columns))


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    amp_master = pd.read_csv(ROOT / "02_cleaned/AMP_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    amr_master = pd.read_csv(ROOT / "02_cleaned/AMR_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    bgc_master = pd.read_csv(ROOT / "02_cleaned/BGC_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    bgc_consensus = pd.read_csv(ROOT / "02_cleaned/BGC_consensus_loci.tsv", sep="\t", dtype=str, keep_default_na=False)
    cazyme_master = pd.read_csv(ROOT / "02_cleaned/CAZyme_master_long.tsv", sep="\t", dtype=str, keep_default_na=False)
    cazyme_primary = cazyme_master[cazyme_master["analysis_tier"] == "primary_consensus"]
    amp_prev = pd.read_csv(ROOT / "06_tables/AMP_prevalence_all.tsv", sep="\t")
    amr_prev = pd.read_csv(ROOT / "06_tables/AMR_family_prevalence_all.tsv", sep="\t")
    bgc_types = pd.read_csv(ROOT / "06_tables/BGC_type_prevalence.tsv", sep="\t")
    bgc_families = pd.read_csv(ROOT / "06_tables/BGC_gene_cluster_families.tsv", sep="\t")
    cazyme_families = pd.read_csv(ROOT / "06_tables/CAZyme_family_prevalence.tsv", sep="\t")

    table1_rows = [
        ("Clinical genomes", 1350, "1350", "authoritative clinical metadata cohort", "02_cleaned/canonical_genome_metadata.tsv"),
        ("Raw files", sum(1 for _ in open(ROOT / "01_audit/funcscan_complete_file_manifest.tsv", encoding="utf-8")) - 1, "", "raw export files inventoried", "01_audit/funcscan_complete_file_manifest.tsv"),
        ("AMP evaluable genomes", 1350, "1350", "usable Ampcombi summary records", "01_audit/tool_output_completeness.tsv"),
        ("AMP predictions", len(amp_master), "", "all standardized prediction records", "02_cleaned/AMP_master_long.tsv"),
        ("High-confidence AMP predictions", int((amp_master["high_confidence"] == "yes").sum()), "", "Tier 1 or Tier 2", "02_cleaned/AMP_master_long.tsv"),
        ("Curated-database-supported AMP predictions", int((amp_master["database_match"] == "yes").sum()), "", "DRAMP/Swiss-Prot evidence in exported fields", "02_cleaned/AMP_master_long.tsv"),
        ("AMP 50%-identity families", len(amp_prev), "", "MMseqs2, 80% bidirectional coverage", "06_tables/AMP_prevalence_all.tsv"),
        ("AMR evaluable genomes", 1350, "1350", "all four AMR caller streams represented", "01_audit/tool_output_completeness.tsv"),
        ("Raw AMR tool calls", len(amr_master), "", "hamronization records", "02_cleaned/AMR_master_long.tsv"),
        ("Tier-A AMR tool calls", int((amr_master["primary_analysis"] == "yes").sum()), "", "curated caller or multi-tool-supported genome/family group", "02_cleaned/AMR_master_long.tsv"),
        ("Unique Tier-A AMR families", len(amr_prev), "", "conservative normalized family names", "06_tables/AMR_family_prevalence_all.tsv"),
        ("BGC evaluable genomes", 1350, "1350", "antiSMASH, DeepBGC, and GECCO collection-level summary", "01_audit/tool_output_completeness.tsv"),
        ("Raw BGC predictions", len(bgc_master), "", "tool-specific region calls", "02_cleaned/BGC_master_long.tsv"),
        ("Consensus BGC loci", len(bgc_consensus), "", "coordinate-overlap merged loci", "02_cleaned/BGC_consensus_loci.tsv"),
        ("Cross-tool-supported BGC loci", int((pd.to_numeric(bgc_consensus["cross_tool_support_count"]) >= 2).sum()), "", "at least two tools", "02_cleaned/BGC_consensus_loci.tsv"),
        ("Normalized BGC types", len(bgc_types), "", "only observed normalized categories", "06_tables/BGC_type_prevalence.tsv"),
        ("Conservative BGC domain-signature families", len(bgc_families), "", "not BiG-SCAPE GCFs", "06_tables/BGC_gene_cluster_families.tsv"),
        ("CAZyme evaluable genomes", 1350, "1350", "run_dbCAN overview available", "01_audit/tool_output_completeness.tsv"),
        ("Primary consensus CAZyme occurrences", len(cazyme_primary), "", "recommended/consensus CAZyme records", "02_cleaned/CAZyme_master_long.tsv"),
        ("Unique CAZyme families", len(cazyme_families), "", "family-level catalogue", "06_tables/CAZyme_family_prevalence.tsv"),
        ("CAZyme gene clusters", sum(1 for _ in open(ROOT / "02_cleaned/CAZyme_CGC_master_long.tsv", encoding="utf-8")) - 1, "", "run_dbCAN CGC records", "02_cleaned/CAZyme_CGC_master_long.tsv"),
        ("Clinical Panaroo clusters", 40811, "", "existing clinical Panaroo Rtab", "../../14.pangenome_panaroo/results/clinical/gene_presence_absence.Rtab"),
    ]
    table1 = pd.DataFrame(table1_rows, columns=["metric", "value", "denominator", "definition_or_scope", "primary_source"])
    write_frame(ROOT / "06_tables/Table_1_FuncScan_dataset_summary.tsv", table1)

    table2_cols = ["AMP", "annotation", "n_genomes", "n_evaluable", "prevalence_percent", "total_occurrences", "high_confidence_occurrences", "database_match", "prevalence_band"]
    write_frame(ROOT / "06_tables/Table_2_Most_prevalent_AMPs.tsv", amp_prev[table2_cols].head(20))
    table3_cols = ["gene_family", "annotation", "n_genomes", "n_evaluable", "prevalence_percent", "total_occurrences", "drug_class", "mechanism", "prediction_tools", "prevalence_band"]
    write_frame(ROOT / "06_tables/Table_3_Most_prevalent_AMR_genes.tsv", amr_prev[table3_cols].head(20))
    table4 = pd.read_csv(ROOT / "06_tables/AMR_drug_class_prevalence.tsv", sep="\t")
    write_frame(ROOT / "06_tables/Table_4_AMR_drug_classes.tsv", table4)

    type_rows = bgc_types.head(15).copy()
    type_rows.insert(0, "record_type", "normalized_BGC_type")
    type_rows = type_rows.rename(columns={"BGC_type": "feature", "total_consensus_regions": "total_predictions"})
    family_rows = bgc_families.head(15).copy()
    family_rows.insert(0, "record_type", "conservative_domain_signature_family")
    family_rows = family_rows.rename(columns={"family_id": "feature", "number_BGCs": "total_predictions"})
    common = ["record_type", "feature", "number_genomes", "n_evaluable", "prevalence_percent", "total_predictions"]
    type_rows = type_rows.rename(columns={"n_genomes": "number_genomes"})
    combined = pd.concat([type_rows.reindex(columns=common), family_rows.reindex(columns=common)], ignore_index=True)
    write_frame(ROOT / "06_tables/Table_5_Most_prevalent_BGCs.tsv", combined)

    burden_summary = pd.read_csv(ROOT / "06_tables/functional_burden_summary.tsv", sep="\t")
    write_frame(ROOT / "06_tables/Table_6_Functional_burden_summary.tsv", burden_summary)
    correlations = pd.read_csv(ROOT / "04_statistics/functional_burden_spearman.tsv", sep="\t")
    write_frame(ROOT / "06_tables/Table_7_Cross_functional_correlations.tsv", correlations)

    association_tables = []
    for label, filename in (
        ("AMP-AMR", "AMP_ARG_associations.tsv"),
        ("AMP-BGC", "AMP_BGC_associations.tsv"),
        ("AMR-BGC", "ARG_BGC_associations.tsv"),
        ("AMR-AMR", "AMR_gene_cooccurrence.tsv"),
    ):
        frame = pd.read_csv(ROOT / f"04_statistics/{filename}", sep="\t")
        frame = frame[frame["q_value"] < 0.05].copy()
        frame["abs_phi"] = frame["Phi_coefficient"].abs()
        frame = frame.sort_values(["abs_phi", "q_value"], ascending=[False, True]).head(10)
        frame.insert(0, "association_class", label)
        association_tables.append(frame)
    strongest = pd.concat(association_tables, ignore_index=True)
    keep = ["association_class", "feature_1", "category_1", "n_genomes_1", "feature_2", "category_2", "n_genomes_2", "n_both", "Phi_coefficient", "Fisher_odds_ratio", "OR_ci_95_low_Haldane_Anscombe", "OR_ci_95_high_Haldane_Anscombe", "q_value", "relationship"]
    write_frame(ROOT / "06_tables/Table_8_Strongest_feature_associations.tsv", strongest[keep])

    amp_clusters = pd.read_csv(ROOT / "06_tables/AMP_cluster_summary.tsv", sep="\t")
    amp_clusters = amp_clusters[amp_clusters["identity_threshold_percent"] == 50].head(30)
    write_frame(ROOT / "06_tables/Table_9_AMP_clusters.tsv", amp_clusters)
    amr_clusters = pd.read_csv(ROOT / "06_tables/AMR_homology_clusters.tsv", sep="\t").head(30)
    write_frame(ROOT / "06_tables/Table_10_AMR_clusters.tsv", amr_clusters)
    write_frame(ROOT / "06_tables/Table_11_BGC_cluster_families.tsv", bgc_families.head(30))

    pangenome_tables = []
    for category, filename in (
        ("AMP", "AMP_pangenome_distribution.tsv"),
        ("AMR", "AMR_pangenome_distribution.tsv"),
        ("BGC-associated gene", "BGC_gene_pangenome_distribution.tsv"),
        ("CAZyme", "CAZyme_pangenome_distribution.tsv"),
    ):
        frame = pd.read_csv(ROOT / f"06_tables/{filename}", sep="\t")
        frame.insert(0, "reported_category", category)
        pangenome_tables.append(frame)
    pangenome = pd.concat(pangenome_tables, ignore_index=True)
    write_frame(ROOT / "06_tables/Table_12_Pangenome_relationships.tsv", pangenome)

    manifest_rows = []
    for number, name in enumerate(
        [
            "Table_1_FuncScan_dataset_summary.tsv",
            "Table_2_Most_prevalent_AMPs.tsv",
            "Table_3_Most_prevalent_AMR_genes.tsv",
            "Table_4_AMR_drug_classes.tsv",
            "Table_5_Most_prevalent_BGCs.tsv",
            "Table_6_Functional_burden_summary.tsv",
            "Table_7_Cross_functional_correlations.tsv",
            "Table_8_Strongest_feature_associations.tsv",
            "Table_9_AMP_clusters.tsv",
            "Table_10_AMR_clusters.tsv",
            "Table_11_BGC_cluster_families.tsv",
            "Table_12_Pangenome_relationships.tsv",
        ],
        start=1,
    ):
        path = ROOT / f"06_tables/{name}"
        manifest_rows.append({"table_number": number, "filename": name, "rows": sum(1 for _ in open(path, encoding="utf-8")) - 1, "path": str(path.relative_to(ROOT))})
    atomic_write_tsv(ROOT / "06_tables/publication_table_manifest.tsv", manifest_rows, list(manifest_rows[0]))
    print("Publication tables: 12 main tables generated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
