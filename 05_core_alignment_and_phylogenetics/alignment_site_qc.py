#!/usr/bin/env python3
"""Stream a FASTA twice to calculate site QC and isolate divergence metrics."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from core_utils import WorkflowError, fasta_records, load_config, write_tsv


BASE_CODES = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)
BASE_NAMES = np.array(["A", "C", "G", "T"])


def read_cohort_ids(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "assembly_accession" not in reader.fieldnames:
            raise WorkflowError(f"Cohort table must contain assembly_accession: {path}")
        return {row["assembly_accession"] for row in reader}


def count_alignment(path: Path, selected: set[str] | None) -> tuple[list[str], np.ndarray, int]:
    identifiers: list[str] = []
    counts: np.ndarray | None = None
    length = -1
    for identifier, _description, sequence in fasta_records(path):
        if selected is not None and identifier not in selected:
            continue
        if length < 0:
            length = len(sequence)
            counts = np.zeros((4, length), dtype=np.uint16)
        elif len(sequence) != length:
            raise WorkflowError(f"Non-uniform alignment: {identifier} has {len(sequence)} bases, expected {length}")
        assert counts is not None
        calls = np.frombuffer(sequence, dtype=np.uint8)
        for index, base in enumerate(BASE_CODES):
            counts[index] += calls == base
        identifiers.append(identifier)
    if counts is None or not identifiers:
        raise WorkflowError(f"No selected sequences were read from {path}")
    if selected is not None and set(identifiers) != selected:
        missing = sorted(selected - set(identifiers))
        extra = sorted(set(identifiers) - selected)
        raise WorkflowError(f"Cohort/alignment mismatch; missing={missing[:10]}, extra={extra[:10]}")
    return identifiers, counts, length


def site_rows(indices: Iterable[int], arrays: dict[str, np.ndarray], scope: str) -> Iterable[dict[str, Any]]:
    for j in indices:
        observed = [BASE_NAMES[k] for k in range(4) if arrays["counts"][k, j] > 0]
        yield {
            "alignment_position": j + 1,
            "valid_acgt_count": int(arrays["valid"][j]),
            "missing_or_ambiguous_count": int(arrays["missing"][j]),
            "gap_count": int(arrays["gap"][j]),
            "site_call_rate": float(arrays["call_rate"][j]),
            "observed_alleles": ",".join(observed),
            "major_allele": BASE_NAMES[arrays["major_index"][j]] if arrays["valid"][j] else "",
            "minor_allele": BASE_NAMES[arrays["minor_index"][j]] if arrays["allele_count"][j] > 1 else "",
            "minor_allele_count": int(arrays["minor_count"][j]),
            "biallelic": bool(arrays["allele_count"][j] == 2),
            "multiallelic": bool(arrays["allele_count"][j] > 2),
            "variable": bool(arrays["allele_count"][j] > 1),
            "invariant": bool(arrays["allele_count"][j] == 1),
            "recombination_status": "unknown",
            "primary_analysis_eligible": False,
            "unfiltered_sensitivity_eligible": bool(arrays["eligible_unfiltered"][j]),
            "coordinate_system": "alignment",
            "analysis_scope": scope,
        }


SITE_FIELDS = [
    "alignment_position", "valid_acgt_count", "missing_or_ambiguous_count", "gap_count", "site_call_rate",
    "observed_alleles", "major_allele", "minor_allele", "minor_allele_count", "biallelic",
    "multiallelic", "variable", "invariant", "recombination_status", "primary_analysis_eligible",
    "unfiltered_sensitivity_eligible", "coordinate_system", "analysis_scope",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--cohort-ids")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--analysis-label", default="unfiltered_sensitivity")
    parser.add_argument("--write-all-sites", action="store_true")
    args = parser.parse_args()
    config, _project_root = load_config(args.config)
    alignment_path = Path(args.alignment).resolve()
    cohort = read_cohort_ids(Path(args.cohort_ids).resolve()) if args.cohort_ids else None
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)

    identifiers, counts, alignment_length = count_alignment(alignment_path, cohort)
    n = len(identifiers)
    valid = counts.sum(axis=0, dtype=np.uint16)
    gap = np.zeros(alignment_length, dtype=np.uint16)
    # A second lightweight pass distinguishes gaps from other ambiguous calls.
    for identifier, _description, sequence in fasta_records(alignment_path):
        if cohort is not None and identifier not in cohort:
            continue
        calls = np.frombuffer(sequence, dtype=np.uint8)
        gap += (calls == ord("-")) | (calls == ord("."))
    missing = n - valid - gap
    call_rate = valid.astype(np.float64) / n
    allele_count = (counts > 0).sum(axis=0)
    major_index = counts.argmax(axis=0)
    order = np.argsort(counts, axis=0, kind="stable")
    minor_index = order[-2]
    minor_count = np.take_along_axis(counts, minor_index[np.newaxis, :], axis=0)[0]
    filters = config["site_filters"]
    high_quality = call_rate >= float(filters["site_call_rate"])
    if filters["exclude_gap_rich_sites"]:
        high_quality &= gap.astype(np.float64) / n <= float(filters["maximum_gap_fraction"])
    eligible_unfiltered = high_quality & (allele_count > 1) & (minor_count >= int(filters["minor_allele_count"]))
    if filters["exclude_multiallelic_primary"] or filters["prefer_biallelic"]:
        eligible_unfiltered &= allele_count == 2
    private = high_quality & (allele_count == 2) & (minor_count == 1)
    multiallelic = high_quality & (allele_count > 2)
    invariant_count = int((allele_count == 1).sum())
    scope = "variant_site_only" if invariant_count == 0 else "invariant_retaining_core_alignment"
    arrays = {
        "counts": counts, "valid": valid, "gap": gap, "missing": missing, "call_rate": call_rate,
        "allele_count": allele_count, "major_index": major_index, "minor_index": minor_index,
        "minor_count": minor_count, "eligible_unfiltered": eligible_unfiltered,
    }

    if args.write_all_sites:
        write_tsv(output / "all_core_site_qc.tsv", site_rows(range(alignment_length), arrays, scope), SITE_FIELDS)
    else:
        write_tsv(output / "all_core_site_qc.tsv",
                  [{"alignment_position": "not_written", "analysis_scope": scope}], SITE_FIELDS)
    write_tsv(output / "eligible_primary_core_sites.tsv", [], SITE_FIELDS)
    write_tsv(output / "eligible_unfiltered_sensitivity_core_sites.tsv",
              site_rows(np.flatnonzero(eligible_unfiltered), arrays, scope), SITE_FIELDS)
    write_tsv(output / "private_variant_sites.tsv", site_rows(np.flatnonzero(private), arrays, scope), SITE_FIELDS)
    write_tsv(output / "multiallelic_sites.tsv", site_rows(np.flatnonzero(multiallelic), arrays, scope), SITE_FIELDS)
    write_tsv(output / "recombination_excluded_sites.tsv", [], SITE_FIELDS)

    major_calls = BASE_CODES[major_index]
    metrics: list[dict[str, Any]] = []
    for identifier, _description, sequence in fasta_records(alignment_path):
        if cohort is not None and identifier not in cohort:
            continue
        calls = np.frombuffer(sequence, dtype=np.uint8)
        callable_mask = np.isin(calls, BASE_CODES)
        callable_all = callable_mask & high_quality
        callable_var = callable_mask & eligible_unfiltered
        nonconsensus = callable_mask & (calls != major_calls)
        raw_count = int((nonconsensus & eligible_unfiltered).sum())
        all_denom = int(callable_all.sum())
        var_denom = int(callable_var.sum())
        metrics.append({
            "assembly_accession": identifier,
            "alignment_sequence_count": n,
            "alignment_length": alignment_length,
            "alignment_scope": scope,
            "analysis_status": args.analysis_label,
            "recombination_status": "unknown_unfiltered",
            "isolate_site_call_rate": float(callable_mask.sum() / alignment_length),
            "raw_non_consensus_allele_count": raw_count,
            "callable_high_quality_core_positions": all_denom,
            "normalized_core_divergence_D_all": float((nonconsensus & high_quality).sum() / all_denom) if all_denom and scope != "variant_site_only" else None,
            "callable_eligible_variable_sites": var_denom,
            "normalized_variant_site_divergence_D_var": float(raw_count / var_denom) if var_denom else None,
            "private_non_consensus_allele_count": int((nonconsensus & private).sum()),
        })
    write_tsv(output / "isolate_core_divergence_metrics.tsv", metrics,
              ["assembly_accession", "alignment_sequence_count", "alignment_length", "alignment_scope",
               "analysis_status", "recombination_status", "isolate_site_call_rate",
               "raw_non_consensus_allele_count", "callable_high_quality_core_positions",
               "normalized_core_divergence_D_all", "callable_eligible_variable_sites",
               "normalized_variant_site_divergence_D_var", "private_non_consensus_allele_count"])
    summary = [
        {"category": "total_alignment_sites", "count": alignment_length},
        {"category": "invariant_sites", "count": invariant_count},
        {"category": "variable_sites", "count": int((allele_count > 1).sum())},
        {"category": "high_quality_sites", "count": int(high_quality.sum())},
        {"category": "eligible_unfiltered_sensitivity_sites", "count": int(eligible_unfiltered.sum())},
        {"category": "eligible_primary_recombination_filtered_sites", "count": 0},
        {"category": "private_variant_sites", "count": int(private.sum())},
        {"category": "multiallelic_sites", "count": int(multiallelic.sum())},
        {"category": "recombination_status", "count": "unknown"},
    ]
    write_tsv(output / "site_filtering_summary.tsv", summary, ["category", "count"])
    print(f"Processed {n} isolates x {alignment_length} columns ({scope})")
    print(f"Unfiltered sensitivity eligible sites: {int(eligible_unfiltered.sum())}")
    print("Primary eligible sites: 0 (recombination status unresolved)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
