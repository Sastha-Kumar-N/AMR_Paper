#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import shutil

REPORT = Path("quast_quality_filter_report.csv")
INPUT_DIR = Path("4.final_unique_genomes")
GOOD_DIR = Path("5.good_quality_genomes")
BAD_DIR = Path("5.failed_quality_genomes")

GOOD_DIR.mkdir(exist_ok=True)
BAD_DIR.mkdir(exist_ok=True)

df = pd.read_csv(REPORT)

missing = []

for _, row in df.iterrows():
    genome = str(row["Genome"])
    status = row["QC_status"]

    possible_files = [
        INPUT_DIR / genome,
        INPUT_DIR / f"{genome}.fna",
        INPUT_DIR / f"{genome}.fa",
        INPUT_DIR / f"{genome}.fasta",
    ]

    src = None
    for p in possible_files:
        if p.exists():
            src = p
            break

    if src is None:
        matches = list(INPUT_DIR.glob(genome + "*"))
        if matches:
            src = matches[0]

    if src is None:
        missing.append(genome)
        continue

    dst_dir = GOOD_DIR if status == "PASS" else BAD_DIR
    shutil.copy2(src, dst_dir / src.name)

print("Copy correction complete.")
print("Good genomes:", len(list(GOOD_DIR.glob("*"))))
print("Failed genomes:", len(list(BAD_DIR.glob("*"))))
print("Total copied:", len(list(GOOD_DIR.glob("*"))) + len(list(BAD_DIR.glob("*"))))

if missing:
    print("\nStill missing:")
    for m in missing:
        print(m)
else:
    print("No missing files.")