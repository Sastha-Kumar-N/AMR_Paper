#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from pathlib import Path

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    rows = []

    id_cols = [
        "genome_file",
        "assembly_accession",
        "biosample_accession",
        "bioproject_accession",
        "organism",
        "strain",
    ]

    for idx, row in df.iterrows():
        for col in df.columns:
            val = str(row.get(col, ""))
            if not val:
                continue
            if AST_RE.search(col) or AST_RE.search(val):
                rows.append({
                    "input_row_index": idx,
                    **{c: row.get(c, "") for c in id_cols if c in df.columns},
                    "matched_column": col,
                    "matched_text": val,
                    "source_database": "LOCAL_METADATA",
                    "source_type": "keyword_hit"
                })

    hits = pd.DataFrame(rows)
    hits.to_csv(outdir / "local_metadata_ast_keyword_hits.tsv", sep="\t", index=False)

    summary = {
        "metadata_rows": len(df),
        "rows_with_local_ast_keyword_hits": hits["input_row_index"].nunique() if not hits.empty else 0,
        "total_local_ast_keyword_hits": len(hits),
    }

    pd.DataFrame([summary]).to_csv(outdir / "local_metadata_ast_keyword_summary.tsv", sep="\t", index=False)

    print(f"[LOCAL] Wrote {outdir / 'local_metadata_ast_keyword_hits.tsv'}")
    print(f"[LOCAL] Rows with hits: {summary['rows_with_local_ast_keyword_hits']} / {summary['metadata_rows']}")

if __name__ == "__main__":
    main()
