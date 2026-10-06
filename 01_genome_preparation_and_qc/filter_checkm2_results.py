#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import shutil

BASE = Path("/mnt/phoenix/v4/1_raw_files")

REPORT = BASE / "7.checkm2/pass/quality_report.tsv"
INPUT_GENOMES = BASE / "6.quast_pass_genomes"

OUT_BASE = BASE / "8.checkm2_filtered_genomes"
PASS_DIR = OUT_BASE / "pass"
REVIEW_DIR = OUT_BASE / "review"
FAIL_DIR = OUT_BASE / "fail"
REPORT_DIR = OUT_BASE / "reports"

for d in [PASS_DIR, REVIEW_DIR, FAIL_DIR, REPORT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(REPORT, sep="\t")

# Detect genome/bin name column
possible_name_cols = ["Name", "Bin Id", "Bin_Id", "Genome", "genome"]
name_col = None
for col in possible_name_cols:
    if col in df.columns:
        name_col = col
        break

if name_col is None:
    raise ValueError(f"Could not identify genome name column. Columns found: {list(df.columns)}")

required_cols = ["Completeness", "Contamination"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing required column: {col}")

def classify(row):
    comp = float(row["Completeness"])
    cont = float(row["Contamination"])

    if comp >= 95 and cont <= 5:
        return "PASS"
    elif comp < 90 or cont > 10:
        return "FAIL"
    else:
        return "REVIEW"

def find_genome_file(name):
    name = str(name)

    candidates = [
        INPUT_GENOMES / name,
        INPUT_GENOMES / f"{name}.fna",
        INPUT_GENOMES / f"{name}.fa",
        INPUT_GENOMES / f"{name}.fasta",
    ]

    for c in candidates:
        if c.exists():
            return c

    # fallback: match by stem
    matches = list(INPUT_GENOMES.glob(f"{name}.*"))
    if matches:
        return matches[0]

    return None

df["checkm2_status"] = df.apply(classify, axis=1)

copied = {"PASS": 0, "REVIEW": 0, "FAIL": 0}
missing = []

for _, row in df.iterrows():
    name = str(row[name_col])
    status = row["checkm2_status"]

    src = find_genome_file(name)

    if src is None:
        missing.append(name)
        continue

    if status == "PASS":
        dest = PASS_DIR / src.name
    elif status == "REVIEW":
        dest = REVIEW_DIR / src.name
    else:
        dest = FAIL_DIR / src.name

    shutil.copy2(src, dest)
    copied[status] += 1

# Save full annotated report
full_report = REPORT_DIR / "checkm2_final_filter_report.csv"
df.to_csv(full_report, index=False)

# Save separate reports
df[df["checkm2_status"] == "PASS"].to_csv(REPORT_DIR / "checkm2_pass_genomes.csv", index=False)
df[df["checkm2_status"] == "REVIEW"].to_csv(REPORT_DIR / "checkm2_review_genomes.csv", index=False)
df[df["checkm2_status"] == "FAIL"].to_csv(REPORT_DIR / "checkm2_fail_genomes.csv", index=False)

# Save missing files if any
missing_file = REPORT_DIR / "checkm2_missing_genome_files.txt"
with open(missing_file, "w") as f:
    for m in missing:
        f.write(m + "\n")

print("CheckM2 filtering complete.")
print()
print("Input CheckM2 report:", REPORT)
print("Input genome folder:", INPUT_GENOMES)
print()
print("Copied genomes:")
print("PASS:", copied["PASS"])
print("REVIEW:", copied["REVIEW"])
print("FAIL:", copied["FAIL"])
print("Missing:", len(missing))
print()
print("Reports saved in:", REPORT_DIR)
print("Main report:", full_report)
