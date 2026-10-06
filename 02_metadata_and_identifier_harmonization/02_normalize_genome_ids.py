#!/usr/bin/env python3
"""Construct an auditable crosswalk from every FuncScan identifier to metadata."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

from analysis_common import (
    METADATA,
    RAW,
    ROOT,
    atomic_write_json,
    atomic_write_tsv,
    canonicalize_genome_id,
    ensure_output_dirs,
    load_clinical_metadata,
)


def unique_values(path: Path, column: str) -> set[str]:
    values: set[str] = set()
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = (row.get(column) or "").strip()
            if value:
                values.add(value)
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    ensure_output_dirs()
    clinical, mapping, canonical_to_raw = load_clinical_metadata()
    rows: list[dict[str, str]] = []

    def add(original: str, source_file: str, source_tool: str, note: str = "") -> None:
        canonical, status = canonicalize_genome_id(original, mapping)
        rows.append(
            {
                "canonical_genome_id": canonical,
                "original_id": original,
                "source_file": source_file,
                "source_tool": source_tool,
                "match_status": status,
                "notes": note,
            }
        )

    for record in clinical:
        raw_id = Path(record["genome_file"]).stem
        add(raw_id, str(METADATA), "clinical metadata", "authoritative metadata mapping")
        add(record["assembly_accession"], str(METADATA), "clinical metadata", "canonical accession")

    sources = [
        (RAW / "AMP/Ampcombi_summary.tsv", "sample_id", "Ampcombi"),
        (RAW / "ARG/hamronization_combined_report.tsv", "input_file_name", "hamronization"),
        (RAW / "BGC/combgc_complete_summary.tsv", "sample_id", "comBGC"),
    ]
    for path, column, tool in sources:
        for value in sorted(unique_values(path, column)):
            add(value, str(path.relative_to(ROOT)), tool)

    for path in sorted((RAW / "ARG").glob("*.normalized.tsv")):
        add(path.name, str(path.relative_to(ROOT)), "argnorm", "identifier derived from filename")
    for path in sorted((RAW / "CAZYME/cazyme/dbcan/cazyme_annotation").glob("*/*_overview.tsv")):
        add(path.parent.name, str(path.relative_to(ROOT)), "run_dbCAN", "identifier derived from directory")

    # Deduplicate identical provenance records without losing any source surface.
    keyed = {
        (
            row["canonical_genome_id"],
            row["original_id"],
            row["source_file"],
            row["source_tool"],
            row["match_status"],
        ): row
        for row in rows
    }
    output = sorted(
        keyed.values(),
        key=lambda r: (r["canonical_genome_id"], r["source_tool"], r["original_id"], r["source_file"]),
    )
    atomic_write_tsv(
        ROOT / "02_cleaned/genome_id_crosswalk.tsv",
        output,
        ["canonical_genome_id", "original_id", "source_file", "source_tool", "match_status", "notes"],
    )

    canonical_rows = []
    for row in sorted(clinical, key=lambda r: r["assembly_accession"]):
        canonical = row["assembly_accession"]
        canonical_rows.append(
            {
                "canonical_genome_id": canonical,
                "raw_genome_id": canonical_to_raw[canonical],
                "genome_file": row["genome_file"],
                "organism": row["organism"],
                "strain": row["strain"],
                "biosample_accession": row["biosample_accession"],
                "assembly_level": row["assembly_level"],
                "isolation_source": row["isolation_source"],
                "host": row["host"],
                "geo_loc_name": row["geo_loc_name"],
                "collection_date": row["collection_date"],
                "clinical_source": row["clinical_source"],
            }
        )
    atomic_write_tsv(
        ROOT / "02_cleaned/canonical_genome_metadata.tsv",
        canonical_rows,
        list(canonical_rows[0]),
    )

    counts = Counter(row["match_status"] for row in output)
    summary = {
        "canonical_genomes": len(canonical_rows),
        "crosswalk_records": len(output),
        "number_exactly_matched": counts["exact"],
        "number_normalized": counts["normalized"],
        "number_ambiguous": counts["ambiguous"],
        "number_unmatched": counts["unmatched"],
    }
    atomic_write_json(ROOT / "01_audit/genome_id_normalization_summary.json", summary)
    atomic_write_tsv(
        ROOT / "01_audit/genome_id_normalization_summary.tsv",
        [summary],
        list(summary),
    )
    print(summary)
    if summary["number_ambiguous"] or summary["number_unmatched"]:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

