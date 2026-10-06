#!/usr/bin/env python3

from pathlib import Path
import hashlib
import shutil
import csv
from collections import defaultdict

INPUT_DIR = Path("2.combined_input_files")
OUTPUT_DIR = Path("3.unique_fasta_files")

REPORT = "unique_fasta_selection_report.csv"

OUTPUT_DIR.mkdir(exist_ok=True)

def sequence_hash(fasta_path):
    h = hashlib.sha256()

    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            h.update(line.upper().encode())

    return h.hexdigest()

files = sorted(
    list(INPUT_DIR.glob("*.fna")) +
    list(INPUT_DIR.glob("*.fa")) +
    list(INPUT_DIR.glob("*.fasta"))
)

hash_to_files = defaultdict(list)

for file in files:
    h = sequence_hash(file)
    hash_to_files[h].append(file)

rows = []

for h, file_list in hash_to_files.items():
    file_list = sorted(file_list)
    representative = file_list[0]

    shutil.copy2(representative, OUTPUT_DIR / representative.name)

    for file in file_list:
        rows.append([
            file.name,
            h,
            representative.name,
            "KEPT" if file == representative else "REMOVED_DUPLICATE",
            len(file_list)
        ])

with open(REPORT, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "file",
        "sequence_hash",
        "representative_kept",
        "status",
        "duplicate_group_size"
    ])
    writer.writerows(rows)

print(f"Input files checked: {len(files)}")
print(f"Unique sequences found: {len(hash_to_files)}")
print(f"Files copied to clean folder: {len(list(OUTPUT_DIR.glob('*')))}")
print(f"Clean folder: {OUTPUT_DIR}")
print(f"Report saved: {REPORT}")
