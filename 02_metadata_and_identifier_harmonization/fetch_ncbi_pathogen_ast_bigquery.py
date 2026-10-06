#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd

def sql_quote_list(values):
    vals = []
    for v in values:
        v = str(v).strip()
        if not v:
            continue
        vals.append("'" + v.replace("'", "\\'") + "'")
    return ", ".join(vals) if vals else "''"

def write_sql(outdir, assemblies, biosamples):
    asm_sql = sql_quote_list(assemblies)
    bs_sql = sql_quote_list(biosamples)

    sql_isolates_ast = f"""
-- NCBI Pathogen Detection AST phenotypes from isolates.AST_phenotypes
-- Run in Google BigQuery.
-- Output target: ncbi_pathogen_isolates_ast_phenotypes.tsv

WITH matched_isolates AS (
  SELECT
    asm_acc,
    biosample_acc,
    target_acc,
    taxgroup_name,
    collection_date,
    geo_loc_name,
    isolation_source,
    host,
    AST_phenotypes
  FROM `ncbi-pathogen-detect.pdbrowser.isolates`
  WHERE asm_acc IN ({asm_sql})
     OR biosample_acc IN ({bs_sql})
)
SELECT
  asm_acc AS ncbi_assembly_accession,
  biosample_acc AS ncbi_biosample_accession,
  target_acc,
  taxgroup_name,
  collection_date,
  geo_loc_name,
  isolation_source,
  host,
  ast.*
FROM matched_isolates,
UNNEST(AST_phenotypes) AS ast
;
""".strip()

    sql_ast_table = f"""
-- NCBI Pathogen Detection standalone AST table joined to isolates by target_acc
-- Run in Google BigQuery.
-- Output target: ncbi_pathogen_ast_table_rows.tsv

WITH matched_targets AS (
  SELECT DISTINCT
    asm_acc,
    biosample_acc,
    target_acc,
    taxgroup_name
  FROM `ncbi-pathogen-detect.pdbrowser.isolates`
  WHERE asm_acc IN ({asm_sql})
     OR biosample_acc IN ({bs_sql})
)
SELECT
  mt.asm_acc AS ncbi_assembly_accession,
  mt.biosample_acc AS ncbi_biosample_accession,
  mt.taxgroup_name,
  a.*
FROM matched_targets mt
JOIN `ncbi-pathogen-detect.pdbrowser.ast` a
ON mt.target_acc = a.target_acc
;
""".strip()

    (outdir / "ncbi_pathogen_isolates_ast_phenotypes.sql").write_text(sql_isolates_ast + "\n")
    (outdir / "ncbi_pathogen_ast_table_rows.sql").write_text(sql_ast_table + "\n")

def run_bigquery(project, outdir):
    from google.cloud import bigquery

    client = bigquery.Client(project=project)

    queries = [
        (
            outdir / "ncbi_pathogen_isolates_ast_phenotypes.sql",
            outdir / "ncbi_pathogen_isolates_ast_phenotypes.tsv",
        ),
        (
            outdir / "ncbi_pathogen_ast_table_rows.sql",
            outdir / "ncbi_pathogen_ast_table_rows.tsv",
        ),
    ]

    for sql_file, out_file in queries:
        sql = sql_file.read_text()
        print(f"[BigQuery] Running {sql_file.name}")
        job = client.query(sql)
        df = job.result().to_dataframe()
        df.to_csv(out_file, sep="\t", index=False)
        print(f"[BigQuery] Wrote {out_file} rows={len(df)}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gcp-project", default="")
    ap.add_argument("--sql-only", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")

    assemblies = sorted(set(df["assembly_accession"].dropna().astype(str))) if "assembly_accession" in df.columns else []
    biosamples = sorted(set(df["biosample_accession"].dropna().astype(str))) if "biosample_accession" in df.columns else []

    assemblies = [x for x in assemblies if x.strip()]
    biosamples = [x for x in biosamples if x.strip()]

    write_sql(outdir, assemblies, biosamples)

    print("[NCBI Pathogen Detection] SQL files written:")
    print(f"  {outdir / 'ncbi_pathogen_isolates_ast_phenotypes.sql'}")
    print(f"  {outdir / 'ncbi_pathogen_ast_table_rows.sql'}")

    if args.sql_only or not args.gcp_project:
        print("[NCBI Pathogen Detection] BigQuery not run.")
        print("To run it, authenticate Google Cloud and rerun with:")
        print("  --gcp-project YOUR_GCP_PROJECT")
        return

    run_bigquery(args.gcp_project, outdir)

if __name__ == "__main__":
    main()
