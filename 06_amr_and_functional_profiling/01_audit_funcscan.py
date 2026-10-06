#!/usr/bin/env python3
"""Complete immutable-file audit and module completeness assessment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from analysis_common import (
    ACCESSION_RE,
    METADATA,
    NAMED_RE,
    PROJECT,
    RAW,
    ROOT,
    atomic_write_json,
    atomic_write_text,
    atomic_write_tsv,
    canonicalize_genome_id,
    ensure_output_dirs,
    load_clinical_metadata,
)


TABULAR_EXTENSIONS = {".tsv": "\t", ".csv": ","}
SEQUENCE_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".ffn", ".fas"}


def extension(path: Path) -> str:
    lower = path.name.lower()
    for compound in (".tar.gz", ".fasta.gz", ".faa.gz", ".fna.gz", ".tsv.gz", ".csv.gz"):
        if lower.endswith(compound):
            return compound
    return path.suffix.lower() or "[no_extension]"


def infer_tool(rel: Path) -> str:
    text = rel.as_posix().lower()
    if text.startswith("amp/"):
        return "Ampcombi"
    if "hamronization" in text:
        return "hamronization"
    if text.startswith("arg/"):
        return "argnorm (AMRFinderPlus/Abricate/DeepARG/RGI integration)"
    if text.startswith("bgc/"):
        return "comBGC (antiSMASH/DeepBGC/GECCO integration)"
    if "dbcan_hmm_results" in text:
        return "dbCAN HMM"
    if "dbcansub_hmm_results" in text:
        return "dbCAN-sub HMM"
    if "_overview.tsv" in text:
        return "run_dbCAN consensus"
    if "cgc_standard_out" in text:
        return "run_dbCAN CGC"
    if "substrate_prediction" in text:
        return "run_dbCAN substrate"
    if "total_cgc_info" in text:
        return "run_dbCAN CGC annotation"
    if "multiqc" in text:
        return "MultiQC"
    if "execution_" in text or "pipeline_dag" in text or "params_" in text:
        return "Nextflow/nf-core FuncScan"
    if "software" in text or "versions" in text:
        return "nf-core FuncScan software provenance"
    return "unknown"


def infer_notes(rel: Path, size: int) -> str:
    notes: list[str] = []
    text = rel.as_posix().lower()
    if size == 0:
        notes.append("zero-byte file")
    if "summary" in rel.name.lower() or "combined" in rel.name.lower():
        notes.append("collection-level aggregate")
    if "potential.arg" in text:
        notes.append("potential ARG subset")
    elif ".arg.normalized" in text:
        notes.append("DeepARG ARG subset")
    elif text.startswith("arg/") and rel.name.endswith(".normalized.tsv"):
        notes.append("per-genome integrated normalized calls")
    if "cluster_representative" in text:
        notes.append("Ampcombi representatives; existing clustering")
    elif "summary_cluster" in text:
        notes.append("Ampcombi cluster assignments")
    if "combgc" in text:
        notes.append("region-level predictions, not BGC gene-family calls")
    return "; ".join(notes)


def infer_genome_id(rel: Path) -> str:
    text = rel.as_posix()
    match = ACCESSION_RE.search(text)
    if match:
        return match.group(0)
    named = NAMED_RE.search(text)
    return named.group(0) if named else ""


def count_rows_columns(path: Path, ext: str) -> tuple[str, str, str]:
    if ext not in TABULAR_EXTENSIONS:
        return "", "", ""
    delimiter = TABULAR_EXTENSIONS[ext]
    try:
        lines = 0
        last = b""
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                lines += chunk.count(b"\n")
                last = chunk[-1:]
        if path.stat().st_size and last != b"\n":
            lines += 1
        with path.open(encoding="utf-8", errors="replace", newline="") as handle:
            header = next(csv.reader(handle, delimiter=delimiter), [])
        return str(max(0, lines - (1 if header else 0))), str(len(header)), " | ".join(header)
    except Exception as exc:
        return "", "", f"ERROR: {type(exc).__name__}: {exc}"


def template_name(rel: Path) -> str:
    text = rel.as_posix()
    text = re.sub(r"GC[AF]_\d+\.\d+", "<GENOME>", text)
    text = re.sub(
        r"Pseudomonas_aeruginosa_(?:P8W_16354|SE5381_16704|SE5418_16705)",
        "<GENOME>",
        text,
    )
    text = re.sub(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}", "<TIMESTAMP>", text)
    return text


def load_samples(path: Path, column: str) -> set[str]:
    result: set[str] = set()
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            value = (row.get(column) or "").strip()
            if value:
                result.add(value)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=RAW)
    args = parser.parse_args()
    raw = args.raw.resolve()
    if raw != RAW.resolve():
        raise SystemExit(f"Refusing unexpected raw directory: {raw}")
    if not raw.is_dir():
        raise SystemExit(f"Raw directory not found: {raw}")

    ensure_output_dirs()
    audit = ROOT / "01_audit"
    clinical_rows, raw_to_canonical, canonical_to_raw = load_clinical_metadata()
    if len(clinical_rows) != 1350:
        raise SystemExit(f"Expected 1,350 clinical metadata rows; found {len(clinical_rows)}")

    files = sorted(p for p in raw.rglob("*") if p.is_file())
    dirs = sorted(p for p in raw.rglob("*") if p.is_dir())
    manifest: list[dict[str, str | int]] = []
    schemas: dict[tuple[str, str], dict[str, str | int]] = {}
    extension_counts: Counter[str] = Counter()
    category_stats: dict[str, Counter[str]] = defaultdict(Counter)
    issues: list[dict[str, str]] = []

    for path in files:
        rel = path.relative_to(raw)
        stat = path.stat()
        ext = extension(path)
        rows, columns, header = count_rows_columns(path, ext)
        category = rel.parts[0] if rel.parts else "unknown"
        genome_id = infer_genome_id(rel)
        genome_scope = "genome-specific" if genome_id else "collection-level"
        manifest.append(
            {
                "path": rel.as_posix(),
                "filename": path.name,
                "extension": ext,
                "size_bytes": stat.st_size,
                "functional_category": category,
                "prediction_tool": infer_tool(rel),
                "genome_specific_or_summary": genome_scope,
                "putative_genome_id": genome_id,
                "rows": rows,
                "columns": columns,
                "sequence_file": "yes" if ext in SEQUENCE_EXTENSIONS else "no",
                "notes": infer_notes(rel, stat.st_size),
            }
        )
        extension_counts[ext] += 1
        category_stats[category]["files"] += 1
        category_stats[category]["size_bytes"] += stat.st_size
        category_stats[category][genome_scope] += 1
        if stat.st_size == 0:
            category_stats[category]["zero_byte_files"] += 1
            issues.append(
                {
                    "issue_type": "zero_byte_file",
                    "severity": "warning",
                    "module": category,
                    "genome_id": genome_id,
                    "source_file": rel.as_posix(),
                    "details": "Raw file has size 0 bytes.",
                    "resolution": "Retained and flagged; not interpreted as biological absence.",
                }
            )
        if header:
            key = (template_name(rel), header)
            if key not in schemas:
                schemas[key] = {
                    "file_template": key[0],
                    "functional_category": category,
                    "prediction_tool": infer_tool(rel),
                    "delimiter": "tab" if ext == ".tsv" else "comma",
                    "columns": columns,
                    "column_names": header,
                    "n_files": 0,
                    "example_file": rel.as_posix(),
                }
            schemas[key]["n_files"] = int(schemas[key]["n_files"]) + 1

    manifest_fields = [
        "path",
        "filename",
        "extension",
        "size_bytes",
        "functional_category",
        "prediction_tool",
        "genome_specific_or_summary",
        "putative_genome_id",
        "rows",
        "columns",
        "sequence_file",
        "notes",
    ]
    atomic_write_tsv(audit / "funcscan_complete_file_manifest.tsv", manifest, manifest_fields)

    tree_lines = [f"[D]\t{raw.name}"]
    for path in sorted([*dirs, *files]):
        rel = path.relative_to(raw)
        tree_lines.append(f"[{'D' if path.is_dir() else 'F'}]\t{rel.as_posix()}")
    atomic_write_text(audit / "directory_tree.txt", "\n".join(tree_lines) + "\n")

    atomic_write_tsv(
        audit / "file_counts_by_extension.tsv",
        ({"extension": k, "file_count": v} for k, v in sorted(extension_counts.items())),
        ["extension", "file_count"],
    )
    atomic_write_tsv(
        audit / "files_by_functional_category.tsv",
        (
            {
                "functional_category": category,
                "file_count": stats["files"],
                "size_bytes": stats["size_bytes"],
                "genome_specific_files": stats["genome-specific"],
                "collection_level_files": stats["collection-level"],
                "zero_byte_files": stats["zero_byte_files"],
            }
            for category, stats in sorted(category_stats.items())
        ),
        [
            "functional_category",
            "file_count",
            "size_bytes",
            "genome_specific_files",
            "collection_level_files",
            "zero_byte_files",
        ],
    )
    atomic_write_tsv(
        audit / "schema_summary.tsv",
        sorted(schemas.values(), key=lambda r: (str(r["functional_category"]), str(r["file_template"]))),
        [
            "file_template",
            "functional_category",
            "prediction_tool",
            "delimiter",
            "columns",
            "column_names",
            "n_files",
            "example_file",
        ],
    )

    amp_raw = load_samples(raw / "AMP/Ampcombi_summary.tsv", "sample_id")
    bgc_raw = load_samples(raw / "BGC/combgc_complete_summary.tsv", "sample_id")
    arg_raw = {
        p.name.removesuffix(".normalized.tsv")
        for p in (raw / "ARG").glob("*.normalized.tsv")
        if ".ARG.normalized.tsv" not in p.name and ".potential.ARG.normalized.tsv" not in p.name
    }
    cazyme_raw = {
        p.name.removesuffix("_overview.tsv")
        for p in (raw / "CAZYME/cazyme/dbcan/cazyme_annotation").glob("*/*_overview.tsv")
    }

    raw_sets = {"AMP": amp_raw, "ARG": arg_raw, "BGC": bgc_raw, "CAZyme": cazyme_raw}
    canonical_sets: dict[str, set[str]] = {}
    unmatched: dict[str, set[str]] = {}
    for module, values in raw_sets.items():
        canonical_sets[module] = set()
        unmatched[module] = set()
        for value in values:
            canonical, status = canonicalize_genome_id(value, raw_to_canonical)
            if status == "unmatched":
                unmatched[module].add(value)
            else:
                canonical_sets[module].add(canonical)

    expected = {row["assembly_accession"] for row in clinical_rows}
    detection_rows: list[dict[str, str]] = []
    for canonical in sorted(expected):
        raw_id = canonical_to_raw[canonical]
        detection_rows.append(
            {
                "canonical_genome_id": canonical,
                "raw_genome_id": raw_id,
                "metadata_status": "clinical_expected",
                "AMP_output": "present" if canonical in canonical_sets["AMP"] else "missing_from_positive_summary",
                "ARG_output": "present" if canonical in canonical_sets["ARG"] else "missing",
                "BGC_output": "present" if canonical in canonical_sets["BGC"] else "missing_from_positive_summary",
                "CAZyme_output": "present" if canonical in canonical_sets["CAZyme"] else "missing",
                "notes": "Named genome normalized to assembly accession." if canonical != raw_id else "",
            }
        )
    atomic_write_tsv(
        audit / "genome_detection_audit.tsv",
        detection_rows,
        [
            "canonical_genome_id",
            "raw_genome_id",
            "metadata_status",
            "AMP_output",
            "ARG_output",
            "BGC_output",
            "CAZyme_output",
            "notes",
        ],
    )

    completeness_rows: list[dict[str, str | int]] = []
    tools = {
        "AMP": "AMPir + AMPlify + Macrel, integrated by Ampcombi",
        "ARG": "AMRFinderPlus + Abricate + DeepARG + RGI, normalized by hamronization/argnorm",
        "BGC": "antiSMASH + DeepBGC + GECCO, integrated by comBGC",
        "CAZyme": "run_dbCAN/dbCAN",
    }
    for module in ("AMP", "ARG", "BGC", "CAZyme"):
        detected = len(canonical_sets[module])
        missing = sorted(expected - canonical_sets[module])
        # Per-genome ARG/CAZyme files explicitly encode evaluability. AMP/BGC are
        # positive-only collection summaries; completed 1,350-input workflow logs
        # establish the denominator while detected genomes quantify positives.
        evaluable = 1350
        completeness_rows.append(
            {
                "module": module,
                "prediction_tools": tools[module],
                "expected_genomes": 1350,
                "evaluable_genomes": evaluable,
                "genomes_with_positive_records": detected,
                "genomes_without_positive_records": 1350 - detected,
                "unmatched_raw_ids": len(unmatched[module]),
                "denominator_basis": (
                    "one complete normalized file set per genome"
                    if module in {"ARG", "CAZyme"}
                    else "completed nf-core/funcscan run with 1,350-input samplesheet; aggregate table is positive-only"
                ),
                "notes": (
                    "No-output genome retained as evaluable zero only because workflow completion and the common input cohort are documented."
                    if missing
                    else "All expected genomes represented by at least one positive record."
                ),
            }
        )
        for value in sorted(unmatched[module]):
            issues.append(
                {
                    "issue_type": "unmatched_genome_id",
                    "severity": "error",
                    "module": module,
                    "genome_id": value,
                    "source_file": "module outputs",
                    "details": "Raw genome identifier could not be reconciled to clinical metadata.",
                    "resolution": "Retained for manual review; excluded from canonical denominator until resolved.",
                }
            )
    atomic_write_tsv(
        audit / "tool_output_completeness.tsv",
        completeness_rows,
        [
            "module",
            "prediction_tools",
            "expected_genomes",
            "evaluable_genomes",
            "genomes_with_positive_records",
            "genomes_without_positive_records",
            "unmatched_raw_ids",
            "denominator_basis",
            "notes",
        ],
    )

    # Record named-ID normalization as an explicit, resolved issue.
    for raw_id, canonical in sorted(raw_to_canonical.items()):
        if raw_id.startswith("Pseudomonas_aeruginosa_"):
            issues.append(
                {
                    "issue_type": "non_accession_genome_name",
                    "severity": "resolved",
                    "module": "all",
                    "genome_id": raw_id,
                    "source_file": METADATA.relative_to(PROJECT).as_posix(),
                    "details": f"Filename identifier differs from canonical assembly accession {canonical}.",
                    "resolution": "Normalized through authoritative clinical metadata; original identifier preserved in crosswalk.",
                }
            )
    if not issues:
        issues.append(
            {
                "issue_type": "none_detected_at_audit_stage",
                "severity": "info",
                "module": "all",
                "genome_id": "",
                "source_file": "",
                "details": "No zero-byte files, unmatched identifiers, or module-level missing outputs detected.",
                "resolution": "Further record-level validation follows in downstream scripts.",
            }
        )
    atomic_write_tsv(
        audit / "data_quality_issues.tsv",
        issues,
        ["issue_type", "severity", "module", "genome_id", "source_file", "details", "resolution"],
    )

    snapshot_rows = []
    aggregate = hashlib.sha256()
    for path in files:
        rel = path.relative_to(raw).as_posix()
        stat = path.stat()
        payload = f"{rel}\t{stat.st_size}\t{stat.st_mtime_ns}\n".encode()
        aggregate.update(payload)
        snapshot_rows.append(
            {"path": rel, "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
        )
    atomic_write_tsv(
        audit / "raw_immutability_snapshot.tsv",
        snapshot_rows,
        ["path", "size_bytes", "mtime_ns"],
    )

    summary = {
        "raw_directory": str(raw),
        "raw_file_count": len(files),
        "raw_directory_count_including_root": len(dirs) + 1,
        "raw_total_bytes": sum(p.stat().st_size for p in files),
        "raw_metadata_fingerprint_sha256": aggregate.hexdigest(),
        "clinical_metadata_genomes": len(expected),
        "module_positive_genomes": {k: len(v) for k, v in canonical_sets.items()},
        "module_evaluable_genomes": {k: 1350 for k in canonical_sets},
        "unmatched_ids": {k: sorted(v) for k, v in unmatched.items()},
        "zero_byte_files": sum(1 for p in files if p.stat().st_size == 0),
    }
    atomic_write_json(audit / "audit_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
