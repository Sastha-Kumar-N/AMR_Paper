#!/usr/bin/env python3
"""Normalize comBGC region predictions and product classes without conflating tools."""

from __future__ import annotations

import argparse
import re
import sys

import pandas as pd

from analysis_common import (
    RAW,
    ROOT,
    atomic_write_tsv,
    canonicalize_genome_id,
    clean_missing,
    ensure_output_dirs,
    load_clinical_metadata,
)


def product_types(value: str) -> list[str]:
    text = clean_missing(value)
    if not text or text.lower() == "unknown":
        return ["Unclassified"]
    low = text.lower()
    types: set[str] = set()
    if "nrp" in low or "nrps" in low:
        types.add("NRPS/NRP")
    if "polyketide" in low or "pks" in low:
        types.add("Polyketide/PKS")
    if "ripp" in low or "ranthipeptide" in low or "bacteriocin" in low:
        types.add("RiPP")
    if "terpene" in low:
        types.add("Terpene")
    if "siderophore" in low or "metallophore" in low:
        types.add("Siderophore/metallophore")
    if "phenazine" in low:
        types.add("Phenazine")
    if "betalactone" in low:
        types.add("Betalactone")
    if "saccharide" in low or "naggn" in low:
        types.add("Saccharide")
    if "cdps" in low:
        types.add("CDPS")
    if "hserlactone" in low:
        types.add("Homoserine lactone")
    if "hydrogen-cyanide" in low:
        types.add("Hydrogen cyanide")
    if "redox-cofactor" in low:
        types.add("Redox cofactor")
    if not types:
        types.add(text if text.lower() != "other" else "Other")
    return sorted(types)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    ensure_output_dirs()
    _, mapping, _ = load_clinical_metadata()
    source = RAW / "BGC/combgc_complete_summary.tsv"
    df = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    master: list[dict[str, object]] = []
    for index, row in enumerate(df.to_dict("records"), start=1):
        canonical, status = canonicalize_genome_id(row["sample_id"], mapping)
        start = int(float(row["BGC_start"]))
        end = int(float(row["BGC_end"]))
        if end < start:
            start, end = end, start
        raw_product = clean_missing(row.get("Product_class")) or "Unclassified"
        normalized_types = product_types(raw_product)
        pfams = sorted(set(re.findall(r"PF\d+", row.get("PFAM_domains", ""), flags=re.I)))
        region_id = f"BGC:{canonical}:{row['contig_id']}:{row['Prediction_tool']}:{start}-{end}:{index}"
        master.append(
            {
                "genome_id": canonical,
                "original_genome_id": row["sample_id"],
                "genome_match_status": status,
                "feature_id": region_id,
                "original_feature_name": raw_product,
                "normalized_feature_name": ";".join(normalized_types),
                "prediction_tool": row["Prediction_tool"],
                "coordinates_if_available": f"{row['contig_id']}:{start}-{end}",
                "contig_id": row["contig_id"],
                "BGC_start": start,
                "BGC_end": end,
                "BGC_length": end - start + 1,
                "sequence_id": "",
                "cluster_id_if_available": "",
                "functional_annotation": raw_product,
                "raw_product_class": raw_product,
                "normalized_BGC_types": ";".join(normalized_types),
                "BGC_probability": row.get("BGC_probability", ""),
                "BGC_complete": row.get("BGC_complete", ""),
                "CDS_count": row.get("CDS_count", ""),
                "BGC_region_contig_ids": row.get("BGC_region_contig_ids", ""),
                "PFAM_domains": row.get("PFAM_domains", ""),
                "PFAM_domain_set": ";".join(pfams),
                "MIBiG_ID": row.get("MIBiG_ID", ""),
                "InterPro_ID": row.get("InterPro_ID", ""),
                "source_file": str(source.relative_to(ROOT)),
            }
        )
    atomic_write_tsv(ROOT / "02_cleaned/BGC_master_long.tsv", master, list(master[0]))
    print(
        f"BGC master: {len(master)} tool-specific regions, "
        f"{len({r['genome_id'] for r in master})} canonical genomes"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

