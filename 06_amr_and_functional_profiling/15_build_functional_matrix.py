#!/usr/bin/env python3
"""Build the 1,350-genome functional burden matrix and phylogeny tip crosswalk."""

from __future__ import annotations

import argparse
import sys

import pandas as pd
from Bio import Phylo

from analysis_common import PHYLOGENY, ROOT, atomic_write_tsv, describe, load_clinical_metadata


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    _, raw_to_canonical, canonical_to_raw = load_clinical_metadata()
    genomes = sorted(canonical_to_raw)
    matrix = pd.DataFrame({"genome_id": genomes})

    inputs = [
        ROOT / "03_matrices/AMP_burden_by_genome.tsv",
        ROOT / "03_matrices/AMR_burden_by_genome.tsv",
        ROOT / "03_matrices/BGC_burden_by_genome.tsv",
        ROOT / "03_matrices/CAZyme_burden_by_genome.tsv",
        ROOT / "02_cleaned/genome_size_cds.tsv",
    ]
    for path in inputs:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if frame["genome_id"].duplicated().any():
            raise ValueError(f"Duplicate genome IDs in {path}")
        matrix = matrix.merge(frame, on="genome_id", how="left", validate="one_to_one")

    rename = {
        "total_ARG_count": "ARG_total",
        "unique_ARG_count": "ARG_unique",
        "drug_class_count": "AMR_drug_classes",
        "resistance_mechanism_count": "AMR_mechanisms",
        "total_BGC_count": "BGC_total",
        "unique_BGC_types": "BGC_types",
        "BGC_region_gene_count": "BGC_gene_count",
    }
    matrix = matrix.rename(columns=rename)
    required = [
        "AMP_total",
        "AMP_unique",
        "AMP_high_confidence",
        "ARG_total",
        "ARG_unique",
        "AMR_drug_classes",
        "AMR_mechanisms",
        "BGC_total",
        "BGC_types",
        "BGC_gene_count",
        "CAZyme_total",
        "CAZyme_families",
        "total_CDS",
        "genome_size_bp",
    ]
    missing = matrix[required].apply(lambda column: column.eq("") | column.isna()).sum()
    if int(missing.sum()) != 0:
        raise ValueError(f"Missing required burdens after merge: {missing[missing > 0].to_dict()}")
    for column in required:
        matrix[column] = pd.to_numeric(matrix[column], errors="raise")

    preferred = [
        "genome_id",
        "AMP_total",
        "AMP_unique",
        "AMP_high_confidence",
        "AMP_high_confidence_unique",
        "AMP_predictor_consensus_fraction",
        "ARG_total",
        "ARG_unique",
        "unique_specific_ARG_names",
        "AMR_drug_classes",
        "AMR_mechanisms",
        "coordinate_supported_ARG_loci",
        "BGC_total",
        "BGC_types",
        "total_BGC_span",
        "BGC_gene_count",
        "cross_tool_consensus_BGC_count",
        "CAZyme_total",
        "CAZyme_gene_count",
        "CAZyme_families",
        "CAZyme_subfamilies",
        "CAZyme_classes",
        "CAZyme_CGC_count",
        "genome_size_bp",
        "total_CDS",
        "contig_count",
        "raw_genome_id",
        "source_gff",
    ]
    ordered = [column for column in preferred if column in matrix.columns]
    extras = [column for column in matrix.columns if column not in ordered and column != "scope"]
    matrix = matrix[[*ordered, *extras]]
    atomic_write_tsv(
        ROOT / "03_matrices/genome_functional_burden_matrix.tsv",
        matrix.to_dict("records"),
        list(matrix.columns),
    )

    summary_rows = []
    for column in required:
        summary_rows.append({"metric": column, **describe(matrix[column].astype(float).tolist())})
    atomic_write_tsv(ROOT / "06_tables/functional_burden_summary.tsv", summary_rows, list(summary_rows[0]))

    tree = Phylo.read(PHYLOGENY, "newick")
    tip_rows = []
    seen = set()
    for index, tip in enumerate(tree.get_terminals(), start=1):
        original = tip.name or ""
        canonical = raw_to_canonical.get(original, original)
        status = "matched" if canonical in canonical_to_raw else "unmatched"
        if status == "matched":
            seen.add(canonical)
        tip_rows.append(
            {
                "tree_tip_order": index,
                "tree_tip_id": original,
                "canonical_genome_id": canonical if status == "matched" else "",
                "match_status": status,
                "source_tree": str(PHYLOGENY),
            }
        )
    for genome in sorted(set(genomes) - seen):
        tip_rows.append(
            {
                "tree_tip_order": "",
                "tree_tip_id": "",
                "canonical_genome_id": genome,
                "match_status": "genome_absent_from_tree",
                "source_tree": str(PHYLOGENY),
            }
        )
    atomic_write_tsv(ROOT / "02_cleaned/clinical_phylogeny_tip_crosswalk.tsv", tip_rows, list(tip_rows[0]))
    print(
        f"Functional matrix: {len(matrix)} genomes x {len(matrix.columns)-1} attributes; "
        f"phylogeny tips matched={sum(row['match_status']=='matched' for row in tip_rows)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
