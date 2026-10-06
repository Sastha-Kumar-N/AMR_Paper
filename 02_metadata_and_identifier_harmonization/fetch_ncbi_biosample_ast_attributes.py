#!/usr/bin/env python3

import argparse
import time
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm
from tenacity import retry, stop_after_attempt, wait_exponential

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

AST_RE = re.compile(
    r"(antibiogram|antimicrobial|susceptib|resistan|resistance|phenotype|"
    r"\bMIC\b|minimum inhibitory|SIR|S/I/R|intermediate|"
    r"ciprofloxacin|levofloxacin|ofloxacin|moxifloxacin|"
    r"meropenem|imipenem|doripenem|ertapenem|carbapenem|"
    r"ceftazidime|cefepime|cefoperazone|ceftriaxone|cefotaxime|"
    r"piperacillin|tazobactam|aztreonam|ticarcillin|"
    r"amikacin|gentamicin|tobramycin|netilmicin|"
    r"colistin|polymyxin|fosfomycin|chloramphenicol|"
    r"beta-lactam|aminoglycoside|fluoroquinolone)",
    re.IGNORECASE,
)

@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=1, min=1, max=20))
def get_url(session, url, params):
    r = session.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r

def esearch_biosample(session, biosample, email=None, api_key=None):
    params = {
        "db": "biosample",
        "term": f"{biosample}[Accession]",
        "retmode": "json",
        "retmax": "5",
        "tool": "ast_metadata_fetcher",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    r = get_url(session, f"{EUTILS}/esearch.fcgi", params)
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])

def efetch_biosample_xml(session, uid, email=None, api_key=None):
    params = {
        "db": "biosample",
        "id": uid,
        "retmode": "xml",
        "tool": "ast_metadata_fetcher",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    r = get_url(session, f"{EUTILS}/efetch.fcgi", params)
    return r.text

def parse_biosample_xml(xml_text, input_biosample):
    root = ET.fromstring(xml_text)

    all_attrs = []
    ast_hits = []

    for bs in root.findall(".//BioSample"):
        biosample_accession = bs.attrib.get("accession", input_biosample)
        biosample_id = bs.attrib.get("id", "")

        title_el = bs.find(".//Description/Title")
        organism_el = bs.find(".//Organism/OrganismName")

        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        organism = "".join(organism_el.itertext()).strip() if organism_el is not None else ""

        for attr in bs.findall(".//Attribute"):
            attr_name = (
                attr.attrib.get("attribute_name")
                or attr.attrib.get("harmonized_name")
                or attr.attrib.get("display_name")
                or ""
            )
            harmonized_name = attr.attrib.get("harmonized_name", "")
            display_name = attr.attrib.get("display_name", "")
            unit = attr.attrib.get("unit", "")
            value = "".join(attr.itertext()).strip()

            rec = {
                "query_biosample_accession": input_biosample,
                "biosample_accession": biosample_accession,
                "biosample_internal_id": biosample_id,
                "title": title,
                "organism": organism,
                "attribute_name": attr_name,
                "harmonized_name": harmonized_name,
                "display_name": display_name,
                "unit": unit,
                "attribute_value": value,
                "source_database": "NCBI_BIOSAMPLE",
            }
            all_attrs.append(rec)

            test_text = f"{attr_name} {harmonized_name} {display_name} {value}"
            if AST_RE.search(test_text):
                ast_hits.append({
                    **rec,
                    "source_type": "biosample_ast_keyword_hit",
                })

    return all_attrs, ast_hits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--email", default="")
    ap.add_argument("--api-key", default="")
    ap.add_argument("--sleep", type=float, default=0.34)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    if "biosample_accession" not in df.columns:
        print("ERROR: biosample_accession column not found", file=sys.stderr)
        sys.exit(1)

    biosamples = sorted(set(x.strip() for x in df["biosample_accession"] if str(x).strip()))
    session = requests.Session()

    all_attrs = []
    ast_hits = []
    status_rows = []

    for biosample in tqdm(biosamples, desc="NCBI BioSample"):
        try:
            ids = esearch_biosample(session, biosample, args.email, args.api_key)
            time.sleep(args.sleep)

            if not ids:
                status_rows.append({
                    "biosample_accession": biosample,
                    "status": "not_found",
                    "uid": "",
                    "error": "",
                })
                continue

            for uid in ids:
                xml_text = efetch_biosample_xml(session, uid, args.email, args.api_key)
                attrs, hits = parse_biosample_xml(xml_text, biosample)
                all_attrs.extend(attrs)
                ast_hits.extend(hits)
                status_rows.append({
                    "biosample_accession": biosample,
                    "status": "fetched",
                    "uid": uid,
                    "error": "",
                })
                time.sleep(args.sleep)

        except Exception as e:
            status_rows.append({
                "biosample_accession": biosample,
                "status": "error",
                "uid": "",
                "error": repr(e),
            })

    pd.DataFrame(all_attrs).to_csv(
        outdir / "ncbi_biosample_all_attributes.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(ast_hits).to_csv(
        outdir / "ncbi_biosample_ast_keyword_hits.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(status_rows).to_csv(
        outdir / "ncbi_biosample_fetch_status.tsv",
        sep="\t",
        index=False,
    )

    print(f"[NCBI BioSample] BioSamples queried: {len(biosamples)}")
    print(f"[NCBI BioSample] AST-like attribute hits: {len(ast_hits)}")

if __name__ == "__main__":
    main()
