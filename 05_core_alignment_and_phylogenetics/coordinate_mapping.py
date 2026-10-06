#!/usr/bin/env python3
"""Validate Panaroo alignment blocks and optional one-to-one PAO1 maps."""

from __future__ import annotations

import argparse
import bisect
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from core_utils import WorkflowError, fasta_records, write_tsv


@dataclass(frozen=True)
class Block:
    start: int
    end: int
    label: str


def parse_panaroo_header(path: Path) -> list[Block]:
    blocks: list[Block] = []
    pending: tuple[int, int] | None = None
    location_re = re.compile(r"^FT\s+feature\s+(\d+)\.\.(\d+)\s*$")
    label_re = re.compile(r"^FT\s+/label=(.+?)\s*$")
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            location = location_re.match(line.rstrip("\n"))
            if location:
                pending = (int(location.group(1)), int(location.group(2)))
                continue
            label = label_re.match(line.rstrip("\n"))
            if label and pending:
                blocks.append(Block(pending[0], pending[1], label.group(1)))
                pending = None
    if not blocks:
        raise WorkflowError(f"No Panaroo blocks parsed from {path}")
    for index, block in enumerate(blocks):
        expected = 1 if index == 0 else blocks[index - 1].end + 1
        if block.start != expected or block.end < block.start:
            raise WorkflowError(
                f"Non-contiguous Panaroo block at index {index}: expected start {expected}, "
                f"observed {block.start}..{block.end} ({block.label})"
            )
    return blocks


def map_positions_to_blocks(positions: list[int], blocks: list[Block]) -> list[dict[str, Any]]:
    ends = [block.end for block in blocks]
    rows: list[dict[str, Any]] = []
    for position in positions:
        index = bisect.bisect_left(ends, position)
        if index >= len(blocks) or position < blocks[index].start:
            raise WorkflowError(f"Alignment position {position} does not map to one Panaroo block")
        block = blocks[index]
        rows.append({
            "alignment_position": position,
            "alignment_start": block.start,
            "alignment_end": block.end,
            "Panaroo_gene_cluster": block.label.removesuffix(".aln"),
            "within_gene_alignment_position": position - block.start + 1,
            "coordinate_system": "concatenated_core_alignment",
        })
    return rows


def parse_vcf_positions(path: Path) -> list[int]:
    positions: list[int] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            positions.append(int(fields[1]))
    return positions


def validate_reference_map(path: Path, fasta_path: Path) -> list[dict[str, Any]]:
    reference = {identifier: sequence.decode() for identifier, _desc, sequence in fasta_records(fasta_path)}
    seen_alignment: set[int] = set()
    seen_reference: set[tuple[str, int]] = set()
    validated: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"alignment_position", "reference_sequence", "reference_coordinate", "reference_allele"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise WorkflowError(f"Reference map lacks required columns: {sorted(required)}")
        for row_number, row in enumerate(reader, start=2):
            aln = int(row["alignment_position"])
            seqid = row["reference_sequence"]
            coord = int(row["reference_coordinate"])
            allele = row["reference_allele"].upper()
            if aln in seen_alignment:
                raise WorkflowError(f"Non-unique alignment position {aln} at mapping row {row_number}")
            if (seqid, coord) in seen_reference:
                raise WorkflowError(f"Non-unique reference coordinate {seqid}:{coord} at mapping row {row_number}")
            if seqid not in reference or not (1 <= coord <= len(reference[seqid])):
                raise WorkflowError(f"Invalid reference coordinate {seqid}:{coord} at mapping row {row_number}")
            observed = reference[seqid][coord - 1]
            if observed != allele:
                raise WorkflowError(
                    f"REF mismatch at mapping row {row_number}: table={allele}, FASTA={observed}, {seqid}:{coord}"
                )
            seen_alignment.add(aln)
            seen_reference.add((seqid, coord))
            validated.append(dict(row))
    return validated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panaroo-header", required=True)
    parser.add_argument("--alignment-length", type=int, required=True)
    parser.add_argument("--vcf")
    parser.add_argument("--reference-map")
    parser.add_argument("--reference-fasta")
    parser.add_argument("--output-dir", default="results/qc")
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    blocks = parse_panaroo_header(Path(args.panaroo_header).resolve())
    if blocks[-1].end != args.alignment_length:
        raise WorkflowError(f"Panaroo block extent {blocks[-1].end} != alignment length {args.alignment_length}")
    block_rows = ({"alignment_start": b.start, "alignment_end": b.end,
                   "Panaroo_gene_cluster": b.label.removesuffix(".aln"),
                   "block_length": b.end - b.start + 1} for b in blocks)
    write_tsv(output / "core_alignment_gene_blocks.tsv", block_rows,
              ["alignment_start", "alignment_end", "Panaroo_gene_cluster", "block_length"])
    status = [{"check": "Panaroo_block_contiguity", "status": "pass",
               "detail": f"{len(blocks)} contiguous blocks cover 1..{blocks[-1].end}"}]
    if args.vcf:
        mapped = map_positions_to_blocks(parse_vcf_positions(Path(args.vcf).resolve()), blocks)
        write_tsv(output / "snp_alignment_to_gene_blocks.tsv", mapped,
                  ["alignment_position", "alignment_start", "alignment_end", "Panaroo_gene_cluster",
                   "within_gene_alignment_position", "coordinate_system"])
        status.append({"check": "VCF_position_to_gene_block_mapping", "status": "pass",
                       "detail": f"{len(mapped)} positions mapped"})
    if args.reference_map or args.reference_fasta:
        if not (args.reference_map and args.reference_fasta):
            raise WorkflowError("--reference-map and --reference-fasta must be supplied together")
        validated = validate_reference_map(Path(args.reference_map).resolve(), Path(args.reference_fasta).resolve())
        status.append({"check": "one_to_one_reference_mapping_and_REF_validation", "status": "pass",
                       "detail": f"{len(validated)} rows validated"})
    else:
        status.append({"check": "one_to_one_reference_mapping_and_REF_validation", "status": "blocked",
                       "detail": "No PAO1 alignment-column mapping was supplied"})
    write_tsv(output / "reference_mapping_status.tsv", status, ["check", "status", "detail"])
    print(status[-1]["detail"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
