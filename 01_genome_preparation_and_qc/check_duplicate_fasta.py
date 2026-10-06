#!/usr/bin/env python3

from pathlib import Path
import hashlib
import csv
from collections import defaultdict

INPUT_DIR = Path("2.combined_input_files")
OUT_REPORT = "duplicate_fasta_report.csv"
OUT_GROUPS = "duplicate_fasta_groups.txt"

def sequence_hash(fasta_path):
    """
    Hash only sequence content.
    Ignores FASTA headers, spaces, and line breaks.
    """
    h = hashlib.sha256()

    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            h.update(line.upper().encode())

    return h.hexdigest()

hash_to_files = defaultdict(list)

files = sorted(
    list(INPUT_DIR.glob("*.fna")) +
    list(INPUT_DIR.glob("*.fa")) +
    list(INPUT_DIR.glob("*.fasta"))
)

for file in files:
    h = sequence_hash(file)
    hash_to_files[h].append(file.name)

duplicate_groups = {
    h: names for h, names in hash_to_files.items()
    if len(names) > 1
}

with open(OUT_REPORT, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["file", "sequence_hash", "duplicate_status", "duplicate_group_size"])

    for h, names in hash_to_files.items():
        status = "DUPLICATE" if len(names) > 1 else "UNIQUE"
        for name in names:
            writer.writerow([name, h, status, len(names)])

with open(OUT_GROUPS, "w") as out:
    if duplicate_groups:
        for i, (h, names) in enumerate(duplicate_groups.items(), start=1):
            out.write(f"Duplicate group {i}\n")
            out.write(f"Hash: {h}\n")
            for name in names:
                out.write(f"  {name}\n")
            out.write("\n")
    else:
        out.write("No exact duplicate sequences found.\n")

total_files = len(files)
duplicate_file_count = sum(len(names) for names in duplicate_groups.values())
duplicate_group_count = len(duplicate_groups)
unique_count = total_files - duplicate_file_count

print(f"Checked files: {total_files}")
print(f"Unique files: {unique_count}")
print(f"Duplicate files: {duplicate_file_count}")
print(f"Duplicate groups: {duplicate_group_count}")
print(f"Report saved: {OUT_REPORT}")
print(f"Duplicate groups saved: {OUT_GROUPS}")

if duplicate_group_count == 0:
    print("\nNo exact duplicate sequences found.")
else:
    print("\nDuplicate sequences found. Check duplicate_fasta_groups.txt")
