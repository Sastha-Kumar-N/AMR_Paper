#!/usr/bin/env python3
"""Normalize all integrated ARG calls while preserving evidence and ambiguity."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict

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


TOOL_PRIORITY = {"amrfinderplus": 0, "rgi": 1, "abricate": 2, "deeparg": 3}


def key_for_gene(symbol: str) -> str:
    text = clean_missing(symbol).upper()
    text = re.sub(r"^PSEUDOMONAS[_ -]+AERUGINOSA[_ -]+", "", text)
    return re.sub(r"[^A-Z0-9]+", "", text)


def family_key(gene_key: str) -> tuple[str, str]:
    key = gene_key
    if key.startswith("BLAPDC"):
        return "BLAPDC", "beta-lactamase allele collapsed to blaPDC family"
    if key.startswith("BLAOXA"):
        return "BLAOXA", "beta-lactamase allele collapsed to blaOXA family"
    for prefix in ("BLAVIM", "BLAIMP", "BLAGES", "BLAKPC", "BLANDM", "BLAPER", "BLACTXM", "BLATEM", "BLASHV"):
        if key.startswith(prefix):
            return prefix, "beta-lactamase allele collapsed to named family"
    if key.startswith("FOSA"):
        return "FOSA", "FosA database allele suffix collapsed"
    match = re.fullmatch(r"CATB\d+", key)
    if match:
        return "CATB", "catB allele collapsed"
    match = re.match(r"(APH\d+[A-Z]*)(?:\d+)?$", key)
    if match and key.startswith("APH"):
        return match.group(1), "aminoglycoside phosphotransferase family normalized conservatively"
    return key, "exact normalized symbol retained"


def broad_drug_classes(value: str) -> list[str]:
    text = clean_missing(value).lower()
    if not text:
        return []
    raw = [x.strip() for x in re.split(r"[;,]", text) if x.strip()]
    classes: set[str] = set()
    for item in raw:
        if any(x in item for x in ("beta-lactam", "cephalosporin", "cephamycin", "carbapenem", "monobactam", "penam", "penem")):
            classes.add("beta-lactam")
        elif "aminoglycoside" in item or item == "kanamycin":
            classes.add("aminoglycoside")
        elif "fluoroquinolone" in item or item == "quinolone":
            classes.add("fluoroquinolone")
        elif "phenicol" in item or "chloramphenicol" in item:
            classes.add("phenicol")
        elif "macrolide" in item or item == "mls":
            classes.add("macrolide/MLS")
        elif "tetracycline" in item:
            classes.add("tetracycline")
        elif "sulfonamide" in item:
            classes.add("sulfonamide")
        elif "diaminopyrimidine" in item or "trimethoprim" in item:
            classes.add("diaminopyrimidine")
        elif "fosfomycin" in item or "phosphonic acid" in item:
            classes.add("fosfomycin/phosphonic acid")
        elif "peptide" in item or "bacitracin" in item or "glycopeptide" in item:
            classes.add("peptide/glycopeptide")
        elif "rifamycin" in item:
            classes.add("rifamycin")
        elif "multidrug" in item:
            classes.add("multidrug")
        elif "disinfect" in item or "antiseptic" in item or item == "triclosan":
            classes.add("biocide/antiseptic")
        elif item not in {"unclassified", "na", "n/a"}:
            classes.add(item.replace(" antibiotic", ""))
    return sorted(classes)


def display_family(key: str, preferred_gene: str) -> str:
    special = {
        "BLAPDC": "blaPDC",
        "BLAOXA": "blaOXA",
        "BLAVIM": "blaVIM",
        "BLAIMP": "blaIMP",
        "BLAGES": "blaGES",
        "BLAKPC": "blaKPC",
        "BLANDM": "blaNDM",
        "BLAPER": "blaPER",
        "BLACTXM": "blaCTX-M",
        "BLATEM": "blaTEM",
        "BLASHV": "blaSHV",
        "FOSA": "fosA",
        "CATB": "catB",
    }
    return special.get(key, preferred_gene)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    ensure_output_dirs()
    _, mapping, _ = load_clinical_metadata()
    source = RAW / "ARG/hamronization_combined_report.tsv"
    df = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    records = df.to_dict("records")

    # Choose a stable display symbol for every punctuation/case-normalized key,
    # preferring curated AMRFinderPlus nomenclature.
    candidates: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in records:
        gene = clean_missing(row.get("gene_symbol")) or "unannotated_ARG"
        key = key_for_gene(gene)
        tool = clean_missing(row.get("analysis_software_name")).lower()
        candidates[key].append((TOOL_PRIORITY.get(tool, 9), len(gene), gene))
    preferred = {key: sorted(values)[0][2] for key, values in candidates.items()}

    provisional: list[dict[str, object]] = []
    group_tools: dict[tuple[str, str], set[str]] = defaultdict(set)
    group_has_amrfinder: set[tuple[str, str]] = set()
    group_has_nonpotential: set[tuple[str, str]] = set()
    for index, row in enumerate(records, start=1):
        canonical, match_status = canonicalize_genome_id(row["input_file_name"], mapping)
        tool = clean_missing(row.get("analysis_software_name")).lower()
        raw_gene = clean_missing(row.get("gene_symbol")) or "unannotated_ARG"
        gene_key = key_for_gene(raw_gene)
        normalized_gene = preferred[gene_key]
        fam_key, normalization_rule = family_key(gene_key)
        family = display_family(fam_key, normalized_gene)
        potential = ".potential." in row["input_file_name"].lower()
        call_tier = "potential_prediction" if potential else "nonpotential_prediction"
        drug_classes = broad_drug_classes(row.get("drug_class", ""))
        contig = clean_missing(row.get("input_sequence_id"))
        start = clean_missing(row.get("input_gene_start"))
        stop = clean_missing(row.get("input_gene_stop"))
        ordinal_match = re.fullmatch(r"(.+)_([0-9]+)", contig)
        if ordinal_match and tool == "deeparg":
            locus_key = f"ordinal:{ordinal_match.group(1)}:{ordinal_match.group(2)}"
        elif contig and start and stop:
            locus_key = f"coordinate:{contig}:{start}-{stop}"
        else:
            locus_key = f"unresolved:{tool}:{index}"
        key = (canonical, fam_key)
        if not potential:
            group_tools[key].add(tool)
            group_has_nonpotential.add(key)
            if tool == "amrfinderplus":
                group_has_amrfinder.add(key)
        provisional.append(
            {
                "genome_id": canonical,
                "original_genome_id": row["input_file_name"],
                "genome_match_status": match_status,
                "feature_id": f"AMR:{index:09d}",
                "original_feature_name": raw_gene,
                "normalized_feature_name": normalized_gene,
                "normalized_gene_key": gene_key,
                "gene_family": family,
                "gene_family_key": fam_key,
                "normalization_rule": normalization_rule,
                "prediction_tool": tool,
                "call_tier": call_tier,
                "genetic_variation_type": row.get("genetic_variation_type", ""),
                "functional_annotation": row.get("gene_name", ""),
                "reference_database_name": row.get("reference_database_name", ""),
                "reference_database_version": row.get("reference_database_version", ""),
                "reference_accession": row.get("reference_accession", ""),
                "sequence_identity": row.get("sequence_identity", ""),
                "coverage_percentage": row.get("coverage_percentage", ""),
                "raw_drug_class": row.get("drug_class", ""),
                "normalized_drug_classes": ";".join(drug_classes),
                "resistance_mechanism": row.get("resistance_mechanism", ""),
                "antimicrobial_agent": row.get("antimicrobial_agent", ""),
                "sequence_id": contig,
                "locus_key": locus_key,
                "contig_id": ordinal_match.group(1) if ordinal_match and tool == "deeparg" else contig,
                "coordinate_start": "" if ordinal_match and tool == "deeparg" else start,
                "coordinate_end": "" if ordinal_match and tool == "deeparg" else stop,
                "coordinates_if_available": f"{contig}:{start}-{stop}" if contig and start and stop and not (ordinal_match and tool == "deeparg") else "",
                "amino_acid_mutation": row.get("amino_acid_mutation", ""),
                "nucleotide_mutation": row.get("nucleotide_mutation", ""),
                "strand_orientation": row.get("strand_orientation", ""),
                "source_file": str(source.relative_to(ROOT)),
            }
        )

    for record in provisional:
        key = (str(record["genome_id"]), str(record["gene_family_key"]))
        tools = sorted(group_tools.get(key, set()))
        if record["call_tier"] == "potential_prediction" and key not in group_has_nonpotential:
            evidence = "Tier C - potential-only"
        elif key in group_has_amrfinder or len(tools) >= 2:
            evidence = "Tier A - curated or multi-tool"
        elif key in group_has_nonpotential:
            evidence = "Tier B - single-tool nonpotential"
        else:
            evidence = "Tier C - potential-only"
        record["gene_family_tool_support"] = ";".join(tools)
        record["gene_family_tool_support_count"] = len(tools)
        record["evidence_tier"] = evidence
        record["primary_analysis"] = "yes" if evidence.startswith("Tier A") else "no"

    atomic_write_tsv(
        ROOT / "02_cleaned/AMR_master_long.tsv",
        provisional,
        list(provisional[0]),
    )

    normalization_rows: list[dict[str, object]] = []
    grouped_raw: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    rules: dict[tuple[str, str], str] = {}
    for record in provisional:
        key = (str(record["original_feature_name"]), str(record["normalized_feature_name"]))
        grouped_raw[key][str(record["prediction_tool"])] += 1
        rules[key] = str(record["normalization_rule"])
    for (raw_gene, normalized_gene), tools in sorted(grouped_raw.items()):
        gene_key = key_for_gene(raw_gene)
        fam_key, _ = family_key(gene_key)
        normalization_rows.append(
            {
                "original_gene_name": raw_gene,
                "normalized_gene_name": normalized_gene,
                "normalized_gene_key": gene_key,
                "gene_family": display_family(fam_key, normalized_gene),
                "gene_family_key": fam_key,
                "normalization_rule": rules[(raw_gene, normalized_gene)],
                "source_tools": ";".join(sorted(tools)),
                "total_raw_calls": sum(tools.values()),
            }
        )
    atomic_write_tsv(
        ROOT / "02_cleaned/AMR_gene_normalization.tsv",
        normalization_rows,
        list(normalization_rows[0]),
    )
    print(
        f"AMR master: {len(provisional)} raw tool calls; "
        f"{sum(r['primary_analysis'] == 'yes' for r in provisional)} calls assigned to Tier-A genome/family groups"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

