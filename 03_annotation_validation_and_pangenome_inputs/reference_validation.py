#!/usr/bin/env python3
"""Validate exact reference FASTA/annotation identity and coordinate bounds."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

from core_utils import WorkflowError, fasta_records, load_config, resolve_input, write_tsv


def parse_gff(path: Path) -> dict[str, Any]:
    sequence_regions: dict[str, int] = {}
    maximum_end: dict[str, int] = {}
    feature_counts: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##sequence-region"):
                fields = line.rstrip().split()
                if len(fields) != 4:
                    raise WorkflowError(f"Malformed ##sequence-region at {path}:{line_number}")
                sequence_regions[fields[1]] = int(fields[3])
            elif line.startswith("#") or not line.strip():
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise WorkflowError(f"Malformed GFF3 row at {path}:{line_number}")
                seqid, feature_type, end = fields[0], fields[2], int(fields[4])
                maximum_end[seqid] = max(maximum_end.get(seqid, 0), end)
                feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1
    return {"sequence_regions": sequence_regions, "maximum_end": maximum_end, "feature_counts": feature_counts}


def parse_genbank_summary(path: Path) -> dict[str, Any]:
    locus_length = None
    accession = None
    version = None
    feature_counts: dict[str, int] = {}
    in_features = False
    origin_bases = 0
    in_origin = False
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("LOCUS"):
                parts = line.split()
                locus_length = int(parts[2])
            elif line.startswith("ACCESSION"):
                accession = line.split()[1]
            elif line.startswith("VERSION"):
                version = line.split()[1]
            elif line.startswith("FEATURES"):
                in_features, in_origin = True, False
            elif line.startswith("ORIGIN"):
                in_features, in_origin = False, True
            elif line.startswith("CONTIG"):
                in_features = False
            elif line.startswith("//"):
                in_features = in_origin = False
            elif in_features and re.match(r"^     \S", line):
                key = line[5:21].strip()
                feature_counts[key] = feature_counts.get(key, 0) + 1
            elif in_origin:
                origin_bases += sum(base in "acgtACGT" for base in line)
    return {"locus_length": locus_length, "accession": accession, "version": version,
            "feature_counts": feature_counts, "origin_bases": origin_bases}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--fasta")
    parser.add_argument("--annotation")
    parser.add_argument("--output", default="results/qc/reference_validation.tsv")
    args = parser.parse_args()
    config, project_root = load_config(args.config)
    fasta_path = resolve_input(project_root, args.fasta or config["inputs"].get("reference_fasta"))
    annotation_path = resolve_input(project_root, args.annotation or config["inputs"].get("reference_annotation"))
    if fasta_path is None or not fasta_path.is_file():
        raise WorkflowError("Exact reference FASTA is missing")

    records = list(fasta_records(fasta_path))
    if len(records) != 1:
        raise WorkflowError(f"Expected one PAO1 chromosome in {fasta_path}, found {len(records)} records")
    seqid, _description, sequence = records[0]
    expected_accession = config["reference"]["accession"]
    expected_length = int(config["reference"]["expected_length"])
    rows = [
        {"check": "fasta_sequence_id", "observed": seqid, "expected": expected_accession,
         "status": "pass" if seqid == expected_accession else "fail"},
        {"check": "fasta_length", "observed": len(sequence), "expected": expected_length,
         "status": "pass" if len(sequence) == expected_length else "fail"},
        {"check": "fasta_alphabet", "observed": "".join(sorted(set(sequence.decode()))), "expected": "ACGT",
         "status": "pass" if set(sequence).issubset(set(b"ACGT")) else "fail"},
    ]
    if annotation_path is None or not annotation_path.is_file():
        rows.append({"check": "reference_annotation", "observed": "missing",
                     "expected": "matching GFF3 or full GenBank", "status": "fail"})
    elif annotation_path.suffix.lower() in {".gff", ".gff3"}:
        info = parse_gff(annotation_path)
        observed_len = info["sequence_regions"].get(seqid)
        max_end = info["maximum_end"].get(seqid, 0)
        rows.extend([
            {"check": "annotation_sequence_id", "observed": ",".join(info["sequence_regions"]),
             "expected": seqid, "status": "pass" if set(info["sequence_regions"]) == {seqid} else "fail"},
            {"check": "annotation_declared_length", "observed": observed_len, "expected": len(sequence),
             "status": "pass" if observed_len == len(sequence) else "fail"},
            {"check": "annotation_max_coordinate", "observed": max_end, "expected": f"<= {len(sequence)}",
             "status": "pass" if 0 < max_end <= len(sequence) else "fail"},
            {"check": "annotation_CDS_count", "observed": info["feature_counts"].get("CDS", 0), "expected": "> 0",
             "status": "pass" if info["feature_counts"].get("CDS", 0) > 0 else "fail"},
        ])
    else:
        info = parse_genbank_summary(annotation_path)
        rows.extend([
            {"check": "annotation_version", "observed": info["version"], "expected": expected_accession,
             "status": "pass" if info["version"] == expected_accession else "fail"},
            {"check": "annotation_locus_length", "observed": info["locus_length"], "expected": len(sequence),
             "status": "pass" if info["locus_length"] == len(sequence) else "fail"},
            {"check": "annotation_CDS_count", "observed": info["feature_counts"].get("CDS", 0), "expected": "> 0",
             "status": "pass" if info["feature_counts"].get("CDS", 0) > 0 else "fail"},
            {"check": "annotation_embedded_sequence_length", "observed": info["origin_bases"], "expected": len(sequence),
             "status": "pass" if info["origin_bases"] == len(sequence) else "fail"},
        ])

    write_tsv((project_root / args.output).resolve(), rows, ["check", "observed", "expected", "status"])
    failures = [row for row in rows if row["status"] == "fail"]
    if failures:
        raise WorkflowError("Reference validation failed: " + "; ".join(row["check"] for row in failures))
    print(f"Validated {seqid} ({len(sequence)} bp) against {annotation_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
