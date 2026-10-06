#!/usr/bin/env python3

from pathlib import Path
import hashlib
import csv
from collections import defaultdict

INPUT_DIR = Path("3.unique_fasta_files")
REPORT = "unique_folder_duplicate_check.csv"
GROUPS = "unique_folder_duplicate_groups.txt"

def contig_set_hash(fasta_path):
    contig_hashes = []
    seq = []

    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq).upper()
                    contig_hashes.append(hashlib.sha256(sequence.encode()).hexdigest())
                    seq = []
            else:
                seq.append(line)

        if seq:
            sequence = "".join(seq).upper()
            contig_hashes.append(hashlib.sha256(sequence.encode()).hexdigest())

    contig_hashes.sort()

    final = hashlib.sha256()
    for h in contig_hashes:
        final.update(h.encode())

    return final.hexdigest(), len(contig_hashes)

files = sorted(
    list(INPUT_DIR.glob("*.fna")) +
    list(INPUT_DIR.glob("*.fa")) +
    list(INPUT_DIR.glob("*.fasta"))
)

hash_to_files = defaultdict(list)
contig_counts = {}

for file in files:
    h, n_contigs = contig_set_hash(file)
    hash_to_files[h].append(file.name)
    contig_counts[file.name] = n_contigs

duplicate_groups = {
    h: names for h, names in hash_to_files.items()
    if len(names) > 1
}

with open(REPORT, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["file", "contig_set_hash", "contig_count", "status", "group_size"])

    for h, names in hash_to_files.items():
        status = "DUPLICATE" if len(names) > 1 else "UNIQUE"
        for name in names:
            writer.writerow([name, h, contig_counts[name], status, len(names)])

with open(GROUPS, "w") as out:
    if duplicate_groups:
        for i, (h, names) in enumerate(duplicate_groups.items(), start=1):
            out.write(f"Duplicate group {i}\n")
            out.write(f"Hash: {h}\n")
            for name in names:
                out.write(f"  {name}\n")
            out.write("\n")
    else:
        out.write("No duplicate genomes found.\n")

print(f"Checked files: {len(files)}")
print(f"Unique genome sets: {len(hash_to_files)}")
print(f"Duplicate groups: {len(duplicate_groups)}")
print(f"Report saved: {REPORT}")
print(f"Groups saved: {GROUPS}")

if duplicate_groups:
    print("Duplicates still found. Check unique_folder_duplicate_groups.txt")
else:
    print("No duplicate genomes found in 3.unique_fasta_files.")
