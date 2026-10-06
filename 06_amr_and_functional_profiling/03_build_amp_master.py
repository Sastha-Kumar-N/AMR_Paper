#!/usr/bin/env python3
"""Normalize AMP predictions, evidence tiers, sequences, and predictor agreement."""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

import pandas as pd

from analysis_common import (
    RAW,
    ROOT,
    atomic_write_text,
    atomic_write_tsv,
    canonicalize_genome_id,
    clean_missing,
    ensure_output_dirs,
    load_clinical_metadata,
)


PREDICTORS = {"AMPir": "prob_ampir", "AMPlify": "prob_amplify", "Macrel": "prob_macrel"}
THRESHOLD = 0.6


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    ensure_output_dirs()
    _, mapping, _ = load_clinical_metadata()
    source = RAW / "AMP/Ampcombi_summary.tsv"
    clustered = RAW / "AMP/Ampcombi_summary_cluster.tsv"
    df = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    cdf = pd.read_csv(
        clustered,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=["sample_id", "CDS_id", "cluster_id"],
    )
    cluster_map = {
        (row.sample_id, row.CDS_id): row.cluster_id
        for row in cdf.itertuples(index=False)
        if clean_missing(row.cluster_id)
    }

    master: list[dict[str, object]] = []
    sequence_map: list[dict[str, str]] = []
    fasta_lines: list[str] = []
    for idx, row in enumerate(df.to_dict("records"), start=1):
        canonical, match_status = canonicalize_genome_id(row["sample_id"], mapping)
        sequence = clean_missing(row.get("aa_sequence")).upper().replace("*", "")
        if not sequence:
            raise ValueError(f"Missing AMP sequence at raw row {idx + 1}")
        seq_hash = hashlib.sha256(sequence.encode()).hexdigest()
        probabilities = {
            name: float(row.get(column) or 0)
            for name, column in PREDICTORS.items()
        }
        supported = [name for name, probability in probabilities.items() if probability >= THRESHOLD]
        db_evalue = None
        try:
            db_evalue = float(row["evalue"]) if clean_missing(row.get("evalue")) else None
        except ValueError:
            db_evalue = None
        strong_database = bool(clean_missing(row.get("DRAMP_ID")) or clean_missing(row.get("Swiss_Prot_Entry"))) and (
            db_evalue is not None and db_evalue <= 1e-5
        )
        if len(supported) == 3 or strong_database:
            tier = "Tier 1 - very high confidence"
        elif len(supported) >= 2:
            tier = "Tier 2 - high confidence"
        else:
            tier = "Tier 3 - putative single-predictor"
        annotation_candidates = [
            clean_missing(row.get("Name")),
            clean_missing(row.get("Gene")),
            clean_missing(row.get("DRAMP_ID")),
            clean_missing(row.get("Swiss_Prot_Entry")),
        ]
        annotation = next((x for x in annotation_candidates if x), "unannotated predicted AMP")
        normalized = annotation if annotation != "unannotated predicted AMP" else f"AMP_seq_{seq_hash[:12]}"
        feature_id = f"AMP:{canonical}:{row['CDS_id']}"
        mmseqs_id = f"AMPSEQ{idx:08d}"
        existing = cluster_map.get((row["sample_id"], row["CDS_id"]), "")
        record = {
            "genome_id": canonical,
            "original_genome_id": row["sample_id"],
            "genome_match_status": match_status,
            "feature_id": feature_id,
            "original_feature_name": annotation,
            "normalized_feature_name": normalized,
            "prediction_tool": ";".join(supported),
            "predictor_support_count": len(supported),
            "predictor_combination": " + ".join(supported) if supported else "none",
            "prob_ampir": row.get("prob_ampir", ""),
            "prob_amplify": row.get("prob_amplify", ""),
            "prob_macrel": row.get("prob_macrel", ""),
            "evidence_tier": tier,
            "high_confidence": "yes" if tier.startswith(("Tier 1", "Tier 2")) else "no",
            "database_match": "yes" if strong_database else "no",
            "DRAMP_ID": row.get("DRAMP_ID", ""),
            "Swiss_Prot_Entry": row.get("Swiss_Prot_Entry", ""),
            "database_target": row.get("target", ""),
            "database_evalue": row.get("evalue", ""),
            "database_percent_identity": row.get("pident", ""),
            "database_query_coverage": row.get("qcov", ""),
            "database_target_coverage": row.get("tcov", ""),
            "functional_annotation": annotation,
            "CDS_id": row["CDS_id"],
            "sequence_id": mmseqs_id,
            "aa_sequence": sequence,
            "sequence_sha256": seq_hash,
            "sequence_length_aa": len(sequence),
            "cluster_id_if_available": f"Ampcombi40_{existing}" if existing else "",
            "contig_id": row.get("contig_id", ""),
            "CDS_start": row.get("CDS_start", ""),
            "CDS_end": row.get("CDS_end", ""),
            "CDS_dir": row.get("CDS_dir", ""),
            "coordinates_if_available": (
                f"{row.get('contig_id', '')}:{row.get('CDS_start', '')}-{row.get('CDS_end', '')}"
            ),
            "molecular_weight": row.get("molecular_weight", ""),
            "helix_fraction": row.get("helix_fraction", ""),
            "turn_fraction": row.get("turn_fraction", ""),
            "sheet_fraction": row.get("sheet_fraction", ""),
            "isoelectric_point": row.get("isoelectric_point", ""),
            "hydrophobicity": row.get("hydrophobicity", ""),
            "transporter_protein": row.get("transporter_protein", ""),
            "CDS_stop_codon_found": row.get("CDS_stop_codon_found", ""),
            "source_file": str(source.relative_to(ROOT)),
        }
        master.append(record)
        sequence_map.append(
            {
                "sequence_id": mmseqs_id,
                "feature_id": feature_id,
                "genome_id": canonical,
                "CDS_id": row["CDS_id"],
                "sequence_sha256": seq_hash,
            }
        )
        fasta_lines.extend([f">{mmseqs_id}", sequence])

    fields = list(master[0])
    atomic_write_tsv(ROOT / "02_cleaned/AMP_master_long.tsv", master, fields)
    atomic_write_tsv(
        ROOT / "02_cleaned/AMP_sequence_id_map.tsv",
        sequence_map,
        list(sequence_map[0]),
    )
    atomic_write_text(ROOT / "02_cleaned/AMP_sequences.faa", "\n".join(fasta_lines) + "\n")

    concordance: list[dict[str, object]] = []
    by_tool = {name: {r["feature_id"] for r in master if name in str(r["prediction_tool"]).split(";")} for name in PREDICTORS}
    for name, ids in by_tool.items():
        genomes = {r["genome_id"] for r in master if r["feature_id"] in ids}
        concordance.append(
            {
                "record_type": "predictor_total",
                "predictor_1": name,
                "predictor_2": "",
                "predictor_combination": name,
                "n_predictions": len(ids),
                "n_genomes": len(genomes),
                "jaccard_similarity": "",
                "proportion_of_all_predictions": len(ids) / len(master),
                "threshold": THRESHOLD,
            }
        )
    names = list(PREDICTORS)
    for i, first in enumerate(names):
        for second in names[i + 1 :]:
            intersection = by_tool[first] & by_tool[second]
            union = by_tool[first] | by_tool[second]
            concordance.append(
                {
                    "record_type": "pairwise_overlap",
                    "predictor_1": first,
                    "predictor_2": second,
                    "predictor_combination": f"{first} + {second}",
                    "n_predictions": len(intersection),
                    "n_genomes": len({r["genome_id"] for r in master if r["feature_id"] in intersection}),
                    "jaccard_similarity": len(intersection) / len(union) if union else 0,
                    "proportion_of_all_predictions": len(intersection) / len(master),
                    "threshold": THRESHOLD,
                }
            )
    combination_counts = pd.Series([r["predictor_combination"] for r in master]).value_counts()
    for combination, count in combination_counts.items():
        members = {r["feature_id"] for r in master if r["predictor_combination"] == combination}
        concordance.append(
            {
                "record_type": "exclusive_combination",
                "predictor_1": "",
                "predictor_2": "",
                "predictor_combination": combination,
                "n_predictions": int(count),
                "n_genomes": len({r["genome_id"] for r in master if r["feature_id"] in members}),
                "jaccard_similarity": "",
                "proportion_of_all_predictions": int(count) / len(master),
                "threshold": THRESHOLD,
            }
        )
    concordance.append(
        {
            "record_type": "support_summary",
            "predictor_1": "",
            "predictor_2": "",
            "predictor_combination": ">=2 predictors",
            "n_predictions": sum(int(r["predictor_support_count"]) >= 2 for r in master),
            "n_genomes": len({r["genome_id"] for r in master if int(r["predictor_support_count"]) >= 2}),
            "jaccard_similarity": "",
            "proportion_of_all_predictions": sum(int(r["predictor_support_count"]) >= 2 for r in master) / len(master),
            "threshold": THRESHOLD,
        }
    )
    concordance.append(
        {
            "record_type": "support_summary",
            "predictor_1": "",
            "predictor_2": "",
            "predictor_combination": "all 3 predictors",
            "n_predictions": sum(int(r["predictor_support_count"]) == 3 for r in master),
            "n_genomes": len({r["genome_id"] for r in master if int(r["predictor_support_count"]) == 3}),
            "jaccard_similarity": "",
            "proportion_of_all_predictions": sum(int(r["predictor_support_count"]) == 3 for r in master) / len(master),
            "threshold": THRESHOLD,
        }
    )
    atomic_write_tsv(
        ROOT / "04_statistics/AMP_predictor_concordance.tsv",
        concordance,
        list(concordance[0]),
    )
    tier_rows = [
        {
            "tier": "Tier 1 - very high confidence",
            "criteria": "all three available predictors at probability >=0.6, or stringent curated-database support (E <= 1e-5)",
            "n_predictions": sum(r["evidence_tier"].startswith("Tier 1") for r in master),
        },
        {
            "tier": "Tier 2 - high confidence",
            "criteria": "two independent predictors at probability >=0.6 without Tier 1 database/all-predictor support",
            "n_predictions": sum(r["evidence_tier"].startswith("Tier 2") for r in master),
        },
        {
            "tier": "Tier 3 - putative single-predictor",
            "criteria": "one predictor at probability >=0.6 and no stringent curated-database support",
            "n_predictions": sum(r["evidence_tier"].startswith("Tier 3") for r in master),
        },
    ]
    atomic_write_tsv(
        ROOT / "04_statistics/AMP_confidence_tier_definition.tsv",
        tier_rows,
        ["tier", "criteria", "n_predictions"],
    )
    print(f"AMP master: {len(master)} predictions across {len({r['genome_id'] for r in master})} genomes")
    return 0


if __name__ == "__main__":
    sys.exit(main())

