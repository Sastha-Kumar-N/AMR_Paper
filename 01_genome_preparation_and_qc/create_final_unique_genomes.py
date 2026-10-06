#!/usr/bin/env python3

from pathlib import Path
import hashlib
import shutil
import csv
from collections import defaultdict

INPUT_DIR = Path("3.unique_fasta_files")
OUTPUT_DIR = Path("4.final_unique_genomes")

REPORT = "final_unique_genome_selection_report.csv"

OUTPUT_DIR.mkdir(exist_ok=True)

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

    return final.hexdigest()

files = sorted(
    list(INPUT_DIR.glob("*.fna")) +
    list(INPUT_DIR.glob("*.fa")) +
    list(INPUT_DIR.glob("*.fasta"))
)

hash_to_files = defaultdict(list)

for file in files:
    h = contig_set_hash(file)
    hash_to_files[h].append(file)

rows = []

for h, file_list in hash_to_files.items():
    file_list = sorted(file_list)

    # Prefer GCA/GCF accession filenames as representative
    accession_files = [
        f for f in file_list
        if f.name.startswith("GCA_") or f.name.startswith("GCF_")
    ]

    if accession_files:
        representative = sorted(accession_files)[0]
    else:
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
        "contig_set_hash",
        "representative_kept",
        "status",
        "duplicate_group_size"
    ])
    writer.writerows(rows)

print(f"Input files checked: {len(files)}")
print(f"Final unique genome sets: {len(hash_to_files)}")
print(f"Files copied to: {OUTPUT_DIR}")
print(f"Files copied count: {len(list(OUTPUT_DIR.glob('*')))}")
print(f"Report saved: {REPORT}")
