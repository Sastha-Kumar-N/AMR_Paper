#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd

def read_tsv(path):
    path = Path(path)
    if path.exists() and path.stat().st_size > 0:
        try:
            return pd.read_csv(path, sep="\t", dtype=str).fillna("")
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    return pd.DataFrame()

def pick_col(df, candidates):
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    for c in df.columns:
        lc = c.lower()
        for cand in candidates:
            if cand.lower() in lc:
                return c
    return None

def normalize_source(df, source_name):
    if df.empty:
        return pd.DataFrame()

    asm_col = pick_col(df, [
        "input_assembly_accession",
        "assembly_accession",
        "ncbi_assembly_accession",
        "asm_acc",
        "bvbrc_genome.assembly_accession",
    ])
    bs_col = pick_col(df, [
        "input_biosample_accession",
        "biosample_accession",
        "query_biosample_accession",
        "ncbi_biosample_accession",
        "biosample_acc",
        "bvbrc_genome.biosample_accession",
    ])
    antibiotic_col = pick_col(df, [
        "antibiotic",
        "antibiotic_name",
        "bvbrc_genome_drug.antibiotic",
        "bvbrc_genome_drug.antibiotic_name",
        "drug",
    ])
    phenotype_col = pick_col(df, [
        "phenotype",
        "resistant_phenotype",
        "resistance_phenotype",
        "bvbrc_genome_drug.resistant_phenotype",
        "sir_label",
    ])
    mic_value_col = pick_col(df, [
        "mic",
        "mic_value",
        "measurement_value",
        "bvbrc_genome_drug.measurement_value",
    ])
    mic_sign_col = pick_col(df, [
        "measurement_sign",
        "mic_sign",
        "bvbrc_genome_drug.measurement_sign",
    ])
    mic_units_col = pick_col(df, [
        "measurement_units",
        "measurement_unit",
        "mic_units",
        "bvbrc_genome_drug.measurement_unit",
    ])
    method_col = pick_col(df, [
        "laboratory_typing_method",
        "lab_typing_method",
        "testing_method",
        "method",
        "bvbrc_genome_drug.laboratory_typing_method",
    ])
    standard_col = pick_col(df, [
        "testing_standard",
        "standard",
        "bvbrc_genome_drug.testing_standard",
    ])

    out = pd.DataFrame()
    out["source_database"] = source_name
    out["assembly_accession"] = df[asm_col] if asm_col else ""
    out["biosample_accession"] = df[bs_col] if bs_col else ""
    out["antibiotic"] = df[antibiotic_col] if antibiotic_col else ""
    out["phenotype_or_sir"] = df[phenotype_col] if phenotype_col else ""
    out["mic_value"] = df[mic_value_col] if mic_value_col else ""
    out["mic_sign"] = df[mic_sign_col] if mic_sign_col else ""
    out["mic_units"] = df[mic_units_col] if mic_units_col else ""
    out["testing_method"] = df[method_col] if method_col else ""
    out["testing_standard"] = df[standard_col] if standard_col else ""

    # Keep full source row for traceability.
    out["source_row_json"] = df.astype(str).apply(lambda r: r.to_json(), axis=1)

    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--results", required=True)
    args = ap.parse_args()

    results = Path(args.results)

    local = read_tsv(results / "local_metadata_ast_keyword_hits.tsv")
    ncbi_biosample = read_tsv(results / "ncbi_biosample_ast_keyword_hits.tsv")
    bvbrc = read_tsv(results / "bvbrc_genome_drug_amr.tsv")
    ncbi_ast1 = read_tsv(results / "ncbi_pathogen_isolates_ast_phenotypes.tsv")
    ncbi_ast2 = read_tsv(results / "ncbi_pathogen_ast_table_rows.tsv")

    parts = [
        normalize_source(local, "LOCAL_METADATA_KEYWORD"),
        normalize_source(ncbi_biosample, "NCBI_BIOSAMPLE"),
        normalize_source(bvbrc, "BV-BRC"),
        normalize_source(ncbi_ast1, "NCBI_PATHOGEN_ISOLATES_AST_PHENOTYPES"),
        normalize_source(ncbi_ast2, "NCBI_PATHOGEN_AST_TABLE"),
    ]

    combined = pd.concat([p for p in parts if not p.empty], ignore_index=True) if any(not p.empty for p in parts) else pd.DataFrame()

    combined.to_csv(results / "combined_ast_long_table.tsv", sep="\t", index=False)

    meta = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    meta["_row_key"] = range(len(meta))

    if combined.empty:
        meta["ast_any_source_count"] = 0
        meta["ast_sources"] = ""
        meta["ast_antibiotics"] = ""
        meta["ast_phenotypes"] = ""
    else:
        c = combined.copy()
        for col in ["assembly_accession", "biosample_accession", "antibiotic", "phenotype_or_sir", "source_database"]:
            if col not in c.columns:
                c[col] = ""

        summary_rows = []
        for i, row in meta.iterrows():
            asm = row.get("assembly_accession", "")
            bs = row.get("biosample_accession", "")

            hit = c[
                ((c["assembly_accession"] != "") & (c["assembly_accession"] == asm)) |
                ((c["biosample_accession"] != "") & (c["biosample_accession"] == bs))
            ]

            summary_rows.append({
                "_row_key": i,
                "ast_any_source_count": len(hit),
                "ast_sources": ";".join(sorted(set(hit["source_database"]))) if len(hit) else "",
                "ast_antibiotics": ";".join(sorted(set(x for x in hit["antibiotic"] if x))) if len(hit) else "",
                "ast_phenotypes": ";".join(sorted(set(x for x in hit["phenotype_or_sir"] if x))) if len(hit) else "",
            })

        summ = pd.DataFrame(summary_rows)
        meta = meta.merge(summ, on="_row_key", how="left")
        meta = meta.drop(columns=["_row_key"])

    meta.to_csv(results / "clinical_metadata_with_ast_summary.tsv", sep="\t", index=False)

    print(f"[MERGE] Wrote {results / 'combined_ast_long_table.tsv'}")
    print(f"[MERGE] Wrote {results / 'clinical_metadata_with_ast_summary.tsv'}")
    if not combined.empty:
        print(f"[MERGE] Combined AST-like rows: {len(combined)}")
    else:
        print("[MERGE] No AST-like rows found from completed sources.")

if __name__ == "__main__":
    main()
