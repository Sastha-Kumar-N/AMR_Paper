#!/usr/bin/env python3
"""Validate paired GCA/GCF representations and prepare a deduplicated cohort manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path
from typing import Iterator

from core_utils import WorkflowError, fasta_records, write_tsv


def gff_fasta_records(path: Path) -> Iterator[tuple[str, bytes]]:
    identifier: str | None = None
    chunks: list[bytes] = []
    in_fasta = False
    with path.open("rb") as handle:
        for raw in handle:
            if raw.startswith(b"##FASTA"):
                in_fasta = True
                continue
            if not in_fasta:
                continue
            if raw.startswith(b">"):
                if identifier is not None:
                    yield identifier, b"".join(chunks).upper()
                identifier = raw[1:].decode("utf-8").strip().split()[0]
                chunks = []
            else:
                chunks.append(raw.strip())
    if identifier is not None:
        yield identifier, b"".join(chunks).upper()


def sequence_summary(path: Path) -> dict[str, object]:
    records = list(gff_fasta_records(path))
    if not records:
        raise WorkflowError(f"No embedded FASTA records in {path}")
    hashes = [hashlib.sha256(sequence).hexdigest() for _identifier, sequence in records]
    largest_index = max(range(len(records)), key=lambda i: len(records[i][1]))
    return {
        "record_count": len(records),
        "total_length": sum(len(sequence) for _identifier, sequence in records),
        "largest_length": len(records[largest_index][1]),
        "largest_sha256": hashes[largest_index],
        "all_sequence_hashes": sorted(hashes),
    }


def selected_fasta_hashes(path: Path, selected: set[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for identifier, _description, sequence in fasta_records(path):
        if identifier in selected:
            result[identifier] = hashlib.sha256(sequence).hexdigest()
    missing = selected - set(result)
    if missing:
        raise WorkflowError(f"Selected IDs missing from {path}: {sorted(missing)}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", default="results/qc/identifier_crosswalk.tsv")
    parser.add_argument("--metadata", default="clinical_metadata.tsv")
    parser.add_argument("--gff-list", required=True)
    parser.add_argument("--snp-alignment")
    parser.add_argument("--output-dir", default="results/qc")
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()

    with Path(args.crosswalk).open("r", encoding="utf-8", newline="") as handle:
        crosswalk = list(csv.DictReader(handle, delimiter="\t"))
    by_canonical: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in crosswalk:
        if row["deduplication_status"] in {"retained_canonical", "excluded_paired_representation"}:
            by_canonical[row["canonical_assembly_accession"]].append(row)
    if not by_canonical:
        raise WorkflowError("No canonical GCA/GCF pairs were identified in the crosswalk")

    gff_paths: dict[str, Path] = {}
    with Path(args.gff_list).open("r", encoding="utf-8") as handle:
        original_gff_lines = [line.strip() for line in handle if line.strip()]
    for value in original_gff_lines:
        gff_paths[Path(value).stem] = Path(value).resolve()

    selected_ids = {row["assembly_accession"] for rows in by_canonical.values() for row in rows}
    snp_hashes = selected_fasta_hashes(Path(args.snp_alignment).resolve(), selected_ids) if args.snp_alignment else {}
    pair_rows: list[dict[str, object]] = []
    excluded: set[str] = set()
    for canonical, rows in sorted(by_canonical.items()):
        if len(rows) != 2:
            raise WorkflowError(f"Expected two records for canonical pair {canonical}, observed {len(rows)}")
        retained = next((row for row in rows if row["deduplication_status"] == "retained_canonical"), None)
        removed = next((row for row in rows if row["deduplication_status"] == "excluded_paired_representation"), None)
        if retained is None or removed is None:
            raise WorkflowError(f"Incomplete retained/excluded designation for {canonical}")
        retained_id, removed_id = retained["assembly_accession"], removed["assembly_accession"]
        if retained_id not in gff_paths or removed_id not in gff_paths:
            raise WorkflowError(f"GFF input missing for pair {removed_id}, {retained_id}")
        retained_summary = sequence_summary(gff_paths[retained_id])
        removed_summary = sequence_summary(gff_paths[removed_id])
        chromosome_identical = retained_summary["largest_sha256"] == removed_summary["largest_sha256"]
        all_records_identical = retained_summary["all_sequence_hashes"] == removed_summary["all_sequence_hashes"]
        snp_identical = snp_hashes.get(retained_id) == snp_hashes.get(removed_id) if snp_hashes else None
        if not chromosome_identical or (snp_identical is False):
            status = "manual_review_sequence_disagreement"
        elif all_records_identical:
            status = "confirmed_exact_sequence_representation_pair"
        else:
            status = "confirmed_same_primary_chromosome_accessory_representation_differs"
        pair_rows.append({
            "biosample_accession": retained["biosample_accession"],
            "excluded_accession": removed_id,
            "retained_accession": retained_id,
            "excluded_GFF": str(gff_paths[removed_id]),
            "retained_GFF": str(gff_paths[retained_id]),
            "primary_chromosome_identical": chromosome_identical,
            "complete_site_SNP_sequence_identical": snp_identical,
            "all_assembly_records_identical": all_records_identical,
            "excluded_record_count": removed_summary["record_count"],
            "retained_record_count": retained_summary["record_count"],
            "excluded_total_length": removed_summary["total_length"],
            "retained_total_length": retained_summary["total_length"],
            "resolution_status": status,
        })
        if status.startswith("manual_review"):
            raise WorkflowError(f"Cannot safely deduplicate {removed_id}/{retained_id}: {status}")
        excluded.add(removed_id)

    write_tsv(output / "paired_gca_gcf_resolution.tsv", pair_rows,
              ["biosample_accession", "excluded_accession", "retained_accession", "excluded_GFF", "retained_GFF",
               "primary_chromosome_identical", "complete_site_SNP_sequence_identical", "all_assembly_records_identical",
               "excluded_record_count", "retained_record_count", "excluded_total_length", "retained_total_length",
               "resolution_status"])
    retained_crosswalk = [row for row in crosswalk if row["assembly_accession"] not in excluded]
    write_tsv(output / "primary_cohort_ids.tsv",
              ({"assembly_accession": row["assembly_accession"]} for row in retained_crosswalk),
              ["assembly_accession"])
    write_tsv(output / "excluded_duplicate_accessions.tsv",
              ({"assembly_accession": value, "reason": "paired_GCA_representation_prefer_GCF"} for value in sorted(excluded)),
              ["assembly_accession", "reason"])
    retained_gffs = [value for value in original_gff_lines if Path(value).stem not in excluded]
    output.mkdir(parents=True, exist_ok=True)
    (output / "deduplicated_clinical_gff_files.list").write_text("\n".join(retained_gffs) + "\n", encoding="utf-8")
    if len(retained_gffs) != len(original_gff_lines) - len(excluded):
        raise WorkflowError("GFF-list attrition does not equal the confirmed duplicate count")
    print(f"Confirmed {len(pair_rows)} GCA/GCF pairs; retained {len(retained_gffs)} primary GFF inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
