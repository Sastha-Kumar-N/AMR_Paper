#!/usr/bin/env python3

import argparse
import time
import urllib.parse
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm
from tenacity import retry, stop_after_attempt, wait_exponential

API = "https://www.bv-brc.org/api"

def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def rql_value(value):
    value = str(value).strip()
    # Accession-like values are safe raw; quote anything more complex.
    if value.replace("_", "").replace(".", "").replace("-", "").isalnum():
        return value
    return '"' + value.replace('"', '\\"') + '"'

@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=1, min=1, max=20))
def api_query(session, dtype, rql):
    url = f"{API}/{dtype}/?{rql}"
    headers = {"Accept": "application/json"}
    r = session.get(url, headers=headers, timeout=90)
    r.raise_for_status()
    data = r.json()

    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        if "response" in data and isinstance(data["response"], dict):
            return data["response"].get("docs", [])
        if "results" in data:
            return data["results"]
        if "data" in data:
            return data["data"]
    return []

def first_nonempty(row, cols):
    for c in cols:
        if c in row and clean(row[c]):
            return clean(row[c])
    return ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sleep", type=float, default=0.25)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    session = requests.Session()

    genome_matches = []
    drug_rows = []
    status_rows = []

    seen_queries = {}

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="BV-BRC genome matching"):
        asm = first_nonempty(row, ["assembly_accession"])
        biosample = first_nonempty(row, ["biosample_accession"])
        bioproject = first_nonempty(row, ["bioproject_accession"])
        strain = first_nonempty(row, ["strain"])

        queries = []
        if asm:
            queries.append(("assembly_accession", asm))
        if biosample:
            queries.append(("biosample_accession", biosample))
        if bioproject:
            queries.append(("bioproject_accession", bioproject))

        matches = []
        match_basis = ""

        for field, value in queries:
            key = (field, value)
            if key in seen_queries:
                qrows = seen_queries[key]
            else:
                rql = f"eq({field},{rql_value(value)})&limit(100)"
                try:
                    qrows = api_query(session, "genome", rql)
                except Exception:
                    qrows = []
                seen_queries[key] = qrows
                time.sleep(args.sleep)

            if qrows:
                matches = qrows
                match_basis = f"{field}={value}"
                break

        if not matches:
            status_rows.append({
                "input_row_index": idx,
                "assembly_accession": asm,
                "biosample_accession": biosample,
                "bioproject_accession": bioproject,
                "strain": strain,
                "status": "no_bvbrc_genome_match",
                "match_basis": "",
                "n_genome_matches": 0,
                "n_drug_rows": 0,
                "error": "",
            })
            continue

        n_drug = 0

        for gm in matches:
            genome_id = clean(gm.get("genome_id", gm.get("genome.genome_id", "")))
            gm_rec = {
                "input_row_index": idx,
                "input_assembly_accession": asm,
                "input_biosample_accession": biosample,
                "input_bioproject_accession": bioproject,
                "input_strain": strain,
                "match_basis": match_basis,
                "source_database": "BV-BRC",
            }

            for k, v in gm.items():
                gm_rec[f"bvbrc_genome.{k}"] = v
            genome_matches.append(gm_rec)

            if not genome_id:
                continue

            try:
                drug_q = f"eq(genome_id,{rql_value(genome_id)})&limit(10000)"
                drugs = api_query(session, "genome_drug", drug_q)
                time.sleep(args.sleep)
            except Exception as e:
                drugs = []
                status_rows.append({
                    "input_row_index": idx,
                    "assembly_accession": asm,
                    "biosample_accession": biosample,
                    "bioproject_accession": bioproject,
                    "strain": strain,
                    "status": "bvbrc_genome_drug_error",
                    "match_basis": match_basis,
                    "n_genome_matches": len(matches),
                    "n_drug_rows": 0,
                    "error": repr(e),
                })

            for dr in drugs:
                rec = {
                    "input_row_index": idx,
                    "input_assembly_accession": asm,
                    "input_biosample_accession": biosample,
                    "input_bioproject_accession": bioproject,
                    "input_strain": strain,
                    "bvbrc_genome_id": genome_id,
                    "match_basis": match_basis,
                    "source_database": "BV-BRC",
                }
                for k, v in dr.items():
                    rec[f"bvbrc_genome_drug.{k}"] = v
                drug_rows.append(rec)
                n_drug += 1

        status_rows.append({
            "input_row_index": idx,
            "assembly_accession": asm,
            "biosample_accession": biosample,
            "bioproject_accession": bioproject,
            "strain": strain,
            "status": "matched",
            "match_basis": match_basis,
            "n_genome_matches": len(matches),
            "n_drug_rows": n_drug,
            "error": "",
        })

    pd.DataFrame(genome_matches).to_csv(outdir / "bvbrc_genome_matches.tsv", sep="\t", index=False)
    pd.DataFrame(drug_rows).to_csv(outdir / "bvbrc_genome_drug_amr.tsv", sep="\t", index=False)
    pd.DataFrame(status_rows).to_csv(outdir / "bvbrc_fetch_status.tsv", sep="\t", index=False)

    print(f"[BV-BRC] Genome matches: {len(genome_matches)}")
    print(f"[BV-BRC] Genome-drug AMR rows: {len(drug_rows)}")

if __name__ == "__main__":
    main()
