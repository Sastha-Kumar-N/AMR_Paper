#!/usr/bin/env python3
"""Normalize run_dbCAN CAZyme calls and carbohydrate gene clusters."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from analysis_common import (
    RAW,
    ROOT,
    atomic_write_tsv,
    canonicalize_genome_id,
    clean_missing,
    ensure_output_dirs,
    extract_cazyme_families,
    load_clinical_metadata,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    ensure_output_dirs()
    _, mapping, _ = load_clinical_metadata()
    overview_root = RAW / "CAZYME/cazyme/dbcan/cazyme_annotation"
    master: list[dict[str, object]] = []
    feature_index = 0
    for path in sorted(overview_root.glob("*/*_overview.tsv")):
        raw_genome = path.parent.name
        canonical, status = canonicalize_genome_id(raw_genome, mapping)
        with path.open(encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                recommended = extract_cazyme_families(row.get("Recommend Results", ""))
                observed = sorted(
                    set(
                        extract_cazyme_families(row.get("dbCAN_hmm", ""))
                        + extract_cazyme_families(row.get("dbCAN_sub", ""))
                        + extract_cazyme_families(row.get("DIAMOND", ""))
                    )
                )
                selected = recommended or observed
                tier = "primary_consensus" if recommended else "single_tool_only"
                for token in selected:
                    feature_index += 1
                    base_match = re.match(r"((?:GH|GT|CE|PL|AA|CBM)\d+)", token, flags=re.I)
                    if not base_match:
                        continue
                    family = base_match.group(1).upper()
                    cazy_class = re.match(r"[A-Z]+", family).group(0)
                    master.append(
                        {
                            "genome_id": canonical,
                            "original_genome_id": raw_genome,
                            "genome_match_status": status,
                            "feature_id": f"CAZ:{feature_index:09d}",
                            "original_feature_name": row.get("Recommend Results", "") or ";".join(observed),
                            "normalized_feature_name": token,
                            "CAZyme_family": family,
                            "CAZyme_subfamily": token if token != family else "",
                            "CAZyme_class": cazy_class,
                            "prediction_tool": "run_dbCAN/dbCAN",
                            "coordinates_if_available": "",
                            "sequence_id": row.get("Gene ID", ""),
                            "cluster_id_if_available": "",
                            "functional_annotation": row.get("Substrate", ""),
                            "confidence_probability": "",
                            "analysis_tier": tier,
                            "tools_supporting": row.get("#ofTools", ""),
                            "EC_number": row.get("EC#", ""),
                            "dbCAN_hmm": row.get("dbCAN_hmm", ""),
                            "dbCAN_sub": row.get("dbCAN_sub", ""),
                            "DIAMOND": row.get("DIAMOND", ""),
                            "substrate": row.get("Substrate", ""),
                            "source_file": str(path.relative_to(ROOT)),
                        }
                    )
    atomic_write_tsv(ROOT / "02_cleaned/CAZyme_master_long.tsv", master, list(master[0]))

    # Aggregate the separate CGC output at cluster level; functional gene roles
    # remain the run_dbCAN labels and are not recast as BGC roles.
    cgc_root = RAW / "CAZYME/cazyme/dbcan/cgc"
    substrate_root = RAW / "CAZYME/cazyme/dbcan/substrate"
    cgc_rows: list[dict[str, object]] = []
    for path in sorted(cgc_root.glob("*/*_cgc_standard_out.tsv")):
        raw_genome = path.parent.name
        canonical, _ = canonicalize_genome_id(raw_genome, mapping)
        substrate_file = substrate_root / raw_genome / f"{raw_genome}_substrate_prediction.tsv"
        substrate_map: dict[str, list[str]] = defaultdict(list)
        if substrate_file.exists():
            with substrate_file.open(encoding="utf-8", errors="replace", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    cgcid = clean_missing(row.get("#cgcid"))
                    for value in (row.get("dbCAN-PUL substrate", ""), row.get("dbCAN-sub substrate", "")):
                        if clean_missing(value):
                            substrate_map[cgcid].append(clean_missing(value))
        groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        with path.open(encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                groups[(row["Contig ID"], row["CGC#"])].append(row)
        for (contig, cgc_id), members in sorted(groups.items()):
            full_id = f"{contig}|{cgc_id}"
            starts = [int(x["Gene Start"]) for x in members]
            stops = [int(x["Gene Stop"]) for x in members]
            cgc_rows.append(
                {
                    "genome_id": canonical,
                    "CGC_id": f"{canonical}:{full_id}",
                    "contig_id": contig,
                    "start": min(starts),
                    "end": max(stops),
                    "gene_count": len(members),
                    "gene_types": ";".join(sorted({x["Gene Type"] for x in members})),
                    "protein_ids": ";".join(x["Protein ID"] for x in members),
                    "gene_annotations": ";".join(sorted({x["Gene Annotation"] for x in members})),
                    "predicted_substrates": ";".join(sorted(set(substrate_map.get(full_id, [])))),
                    "source_file": str(path.relative_to(ROOT)),
                }
            )
    atomic_write_tsv(
        ROOT / "02_cleaned/CAZyme_CGC_master_long.tsv",
        cgc_rows,
        list(cgc_rows[0]) if cgc_rows else ["genome_id", "CGC_id"],
    )
    print(
        f"CAZyme master: {len(master)} gene-family records; "
        f"{sum(r['analysis_tier'] == 'primary_consensus' for r in master)} primary consensus records; "
        f"{len(cgc_rows)} CGCs"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
