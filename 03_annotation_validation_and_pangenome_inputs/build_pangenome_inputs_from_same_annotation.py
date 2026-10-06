#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import os
import sys
from collections import Counter, defaultdict

ANNOT_DIR = Path("/mnt/phoenix_2/12.annotation/prokka_all")
META_FILE = Path("/mnt/phoenix_2/11.metadata/metadata/final_genome_metadata_classified.tsv")
OUT_DIR = Path("/mnt/phoenix_2/13.pangenome_inputs")

OUT_WHOLE = OUT_DIR / "whole"
OUT_CLINICAL = OUT_DIR / "clinical"
OUT_ENV = OUT_DIR / "environmental"
OUT_UNKNOWN = OUT_DIR / "unknown_other"

for d in [OUT_WHOLE, OUT_CLINICAL, OUT_ENV, OUT_UNKNOWN]:
    d.mkdir(parents=True, exist_ok=True)

manifest_file = OUT_DIR / "pangenome_input_manifest.tsv"
summary_file = OUT_DIR / "pangenome_input_summary.tsv"

acc_re = re.compile(r"(GC[AF]_\d+\.\d+)")

def normalize_class(x):
    x0 = str(x).strip()
    x = x0.lower().replace("-", "_").replace(" ", "_")

    if "clinical" in x:
        return "Clinical"
    if "environmental" in x:
        return "Environmental"
    if "unknown" in x or "other" in x:
        return "Unknown_Other"

    return x0

def detect_accession_column(headers, rows):
    preferred = [
        "assembly_accession",
        "accession",
        "genome_accession",
        "assembly_id",
        "genome_id",
        "file_stem",
    ]

    lower_map = {h.lower(): h for h in headers}

    for p in preferred:
        if p in lower_map:
            return lower_map[p]

    for h in headers:
        hits = 0
        for r in rows[:100]:
            if acc_re.search(str(r.get(h, ""))):
                hits += 1
        if hits >= 5:
            return h

    raise RuntimeError("Could not detect accession column in metadata file.")

def detect_class_column(headers, rows):
    preferred = [
        "final_class",
        "classification",
        "source_class",
        "isolate_class",
        "isolate_category",
        "category",
        "group",
        "metadata_class",
    ]

    lower_map = {h.lower(): h for h in headers}

    for p in preferred:
        if p in lower_map:
            return lower_map[p]

    best_col = None
    best_hits = 0

    for h in headers:
        hits = 0
        for r in rows:
            val = str(r.get(h, "")).strip().lower().replace("-", "_").replace(" ", "_")
            if (
                val in {"clinical", "environmental", "unknown_other", "unknown", "other"}
                or "clinical" in val
                or "environmental" in val
            ):
                hits += 1

        if hits > best_hits:
            best_hits = hits
            best_col = h

    if best_col and best_hits > 10:
        return best_col

    raise RuntimeError("Could not detect classification column in metadata file.")

def index_gff_files(annot_dir):
    gff_index = {}
    duplicate_sources = defaultdict(list)

    for gff in annot_dir.rglob("*.gff"):
        accession = None

        m = acc_re.search(gff.name)
        if m:
            accession = m.group(1)

        if accession is None:
            m = acc_re.search(str(gff))
            if m:
                accession = m.group(1)

        if accession:
            if accession in gff_index:
                duplicate_sources[accession].append(str(gff))
            else:
                gff_index[accession] = gff

    return gff_index, duplicate_sources

def force_symlink(src, dest):
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    os.symlink(src, dest)

if not ANNOT_DIR.exists():
    raise RuntimeError(f"Annotation directory not found: {ANNOT_DIR}")

if not META_FILE.exists():
    raise RuntimeError(f"Metadata file not found: {META_FILE}")

with META_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    headers = reader.fieldnames

if not rows:
    raise RuntimeError("Metadata file has no rows.")

acc_col = detect_accession_column(headers, rows)
class_col = detect_class_column(headers, rows)

print(f"Metadata file: {META_FILE}")
print(f"Annotation directory: {ANNOT_DIR}")
print(f"Detected accession column: {acc_col}")
print(f"Detected classification column: {class_col}")

gff_index, duplicate_sources = index_gff_files(ANNOT_DIR)

print(f"Indexed GFF files: {len(gff_index)}")

if duplicate_sources:
    print(f"Warning: duplicate GFF mappings found for {len(duplicate_sources)} accessions.")
    print("Using first discovered GFF for each accession.")

summary = Counter()
missing = []

with manifest_file.open("w", newline="") as mf:
    writer = csv.writer(mf, delimiter="\t")
    writer.writerow([
        "accession",
        "classification",
        "source_gff",
        "whole_link",
        "group_link",
        "status"
    ])

    for r in rows:
        raw_acc = str(r.get(acc_col, "")).strip()
        m = acc_re.search(raw_acc)

        if not m:
            continue

        acc = m.group(1)
        cls = normalize_class(r.get(class_col, ""))

        src = gff_index.get(acc)

        if not src:
            missing.append(acc)
            writer.writerow([acc, cls, "", "", "", "MISSING_GFF"])
            summary["missing_gff"] += 1
            continue

        whole_link = OUT_WHOLE / f"{acc}.gff"
        force_symlink(src, whole_link)
        summary["whole"] += 1

        if cls == "Clinical":
            group_link = OUT_CLINICAL / f"{acc}.gff"
            force_symlink(src, group_link)
            summary["clinical"] += 1

        elif cls == "Environmental":
            group_link = OUT_ENV / f"{acc}.gff"
            force_symlink(src, group_link)
            summary["environmental"] += 1

        else:
            group_link = OUT_UNKNOWN / f"{acc}.gff"
            force_symlink(src, group_link)
            summary["unknown_other"] += 1

        writer.writerow([acc, cls, src, whole_link, group_link, "OK"])

with summary_file.open("w") as sf:
    sf.write("dataset\tgff_count\n")
    for k in ["whole", "clinical", "environmental", "unknown_other", "missing_gff"]:
        sf.write(f"{k}\t{summary[k]}\n")

print("\nPangenome input generation complete.")
print(f"Manifest: {manifest_file}")
print(f"Summary:  {summary_file}")

print("\nCounts:")
for k in ["whole", "clinical", "environmental", "unknown_other", "missing_gff"]:
    print(f"{k}: {summary[k]}")

if missing:
    print("\nMissing GFF files:")
    for acc in missing[:50]:
        print(acc)
    if len(missing) > 50:
        print(f"... plus {len(missing) - 50} more")

if summary["missing_gff"] > 0:
    sys.exit(1)
