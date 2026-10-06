#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import shutil
import re
import numpy as np

# ---------------- USER SETTINGS ----------------

QUAST_REPORT = Path("quast_qc_results/final_unique_genomes_quast/report.tsv")
INPUT_DIR = Path("4.final_unique_genomes")

PASS_DIR = Path("6.quast_pass_genomes")
FAIL_DIR = Path("6.quast_fail_genomes")
REVIEW_DIR = Path("6.quast_review_gc_peak_genomes")

SUMMARY_OUT = Path("quast_final_filter_report.csv")

# P. aeruginosa lower genome-size filter
# Do NOT remove large genomes at this stage.
MIN_GENOME_SIZE = 5_000_000

# Assembly quality thresholds
MAX_CONTIGS = 500
MIN_N50 = 50_000
MIN_LARGEST_CONTIG = 100_000
MAX_NS_PER_100KBP = 10

# GC peak detection
GC_WINDOW_MIN = 5000
GC_BIN_WIDTH = 2.0
MIN_SECOND_PEAK_FRACTION = 0.15
MIN_PEAK_DISTANCE = 6.0

# ------------------------------------------------

for d in [PASS_DIR, FAIL_DIR, REVIEW_DIR]:
    d.mkdir(exist_ok=True)

def clean_number(x):
    if pd.isna(x):
        return np.nan
    x = str(x).replace(",", "").replace("%", "").strip()
    if x in ["-", ""]:
        return np.nan
    return float(x)

def read_fasta_records(fasta):
    header = None
    seqs = []

    with open(fasta, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seqs).upper()
                header = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)

        if header is not None:
            yield header, "".join(seqs).upper()

def contig_gc_values(fasta):
    values = []

    for header, seq in read_fasta_records(fasta):
        length = len(seq)
        if length < GC_WINDOW_MIN:
            continue

        gc = seq.count("G") + seq.count("C")
        atgc = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")

        if atgc == 0:
            continue

        values.append((gc / atgc) * 100)

    return values

def detect_gc_peak_status(gc_values):
    if len(gc_values) < 5:
        return "NOT_ENOUGH_CONTIGS", np.nan, np.nan, 0

    bins = np.arange(30, 81, GC_BIN_WIDTH)
    hist, edges = np.histogram(gc_values, bins=bins)

    if hist.max() == 0:
        return "NO_GC_DATA", np.nan, np.nan, 0

    peak_indices = []

    for i in range(len(hist)):
        left = hist[i - 1] if i > 0 else 0
        right = hist[i + 1] if i < len(hist) - 1 else 0

        if hist[i] >= left and hist[i] >= right and hist[i] > 0:
            peak_indices.append(i)

    peak_indices = sorted(peak_indices, key=lambda i: hist[i], reverse=True)

    if len(peak_indices) < 2:
        main_peak = (edges[peak_indices[0]] + edges[peak_indices[0] + 1]) / 2
        return "SINGLE_PEAK", main_peak, np.nan, 1

    p1 = peak_indices[0]
    p2 = peak_indices[1]

    peak1_gc = (edges[p1] + edges[p1 + 1]) / 2
    peak2_gc = (edges[p2] + edges[p2 + 1]) / 2

    peak1_count = hist[p1]
    peak2_count = hist[p2]

    distance = abs(peak1_gc - peak2_gc)
    second_fraction = peak2_count / peak1_count

    if distance >= MIN_PEAK_DISTANCE and second_fraction >= MIN_SECOND_PEAK_FRACTION:
        return "DUAL_PEAK_REVIEW", peak1_gc, peak2_gc, 2

    return "SINGLE_PEAK", peak1_gc, peak2_gc, 1

def find_genome_file(genome_name):
    genome_name = str(genome_name)
    base = Path(genome_name).name

    candidates = [
        INPUT_DIR / genome_name,
        INPUT_DIR / base,
        INPUT_DIR / f"{genome_name}.fna",
        INPUT_DIR / f"{genome_name}.fa",
        INPUT_DIR / f"{genome_name}.fasta",
        INPUT_DIR / f"{base}.fna",
        INPUT_DIR / f"{base}.fa",
        INPUT_DIR / f"{base}.fasta",
    ]

    for c in candidates:
        if c.exists():
            return c

    matches = list(INPUT_DIR.glob(base + "*"))
    if matches:
        return matches[0]

    return None

def qc_decision(row):
    reasons = []

    if row["total_length"] < MIN_GENOME_SIZE:
        reasons.append("small_incomplete_genome")

    if row["contigs"] > MAX_CONTIGS:
        reasons.append("high_contig_count")

    if row["n50"] < MIN_N50:
        reasons.append("low_N50")

    if row["largest_contig"] < MIN_LARGEST_CONTIG:
        reasons.append("small_largest_contig")

    if not pd.isna(row["ns_per_100kbp"]) and row["ns_per_100kbp"] > MAX_NS_PER_100KBP:
        reasons.append("high_N_content")

    if reasons:
        return "FAIL", ";".join(reasons)

    return "PASS", "passed_quast_basic_filters"

# Read QUAST report
raw = pd.read_csv(QUAST_REPORT, sep="\t")

# QUAST report.tsv has metrics as rows and genomes as columns
df = raw.set_index(raw.columns[0]).T
df.index.name = "Genome"
df = df.reset_index()

required_metrics = {
    "# contigs": "contigs",
    "Total length": "total_length",
    "GC (%)": "gc",
    "N50": "n50",
    "Largest contig": "largest_contig",
    "# N's per 100 kbp": "ns_per_100kbp",
}

for old, new in required_metrics.items():
    if old in df.columns:
        df[new] = df[old].apply(clean_number)
    else:
        df[new] = np.nan

results = []

for _, row in df.iterrows():
    genome = row["Genome"]
    src = find_genome_file(genome)

    status, reason = qc_decision(row)

    gc_peak_status = "FILE_NOT_FOUND"
    gc_peak_1 = np.nan
    gc_peak_2 = np.nan
    gc_peak_count = 0

    if src is not None:
        gc_values = contig_gc_values(src)
        gc_peak_status, gc_peak_1, gc_peak_2, gc_peak_count = detect_gc_peak_status(gc_values)

    final_status = status

    if status == "PASS" and gc_peak_status == "DUAL_PEAK_REVIEW":
        final_status = "REVIEW"

    if src is None:
        final_status = "MISSING_FILE"
        reason = "file_not_found"

    out_row = row.to_dict()
    out_row["source_file"] = str(src) if src else "NOT_FOUND"
    out_row["basic_quast_status"] = status
    out_row["basic_quast_reason"] = reason
    out_row["gc_peak_status"] = gc_peak_status
    out_row["gc_peak_1"] = gc_peak_1
    out_row["gc_peak_2"] = gc_peak_2
    out_row["gc_peak_count"] = gc_peak_count
    out_row["final_status"] = final_status

    results.append(out_row)

    if src is not None:
        if final_status == "PASS":
            shutil.copy2(src, PASS_DIR / src.name)
        elif final_status == "FAIL":
            shutil.copy2(src, FAIL_DIR / src.name)
        elif final_status == "REVIEW":
            shutil.copy2(src, REVIEW_DIR / src.name)

final_df = pd.DataFrame(results)
final_df.to_csv(SUMMARY_OUT, index=False)

print("Filtering complete.")
print(final_df["final_status"].value_counts(dropna=False))
print()
print(f"PASS genomes copied to:   {PASS_DIR}")
print(f"FAIL genomes copied to:   {FAIL_DIR}")
print(f"REVIEW genomes copied to: {REVIEW_DIR}")
print(f"Full report saved as:     {SUMMARY_OUT}")
