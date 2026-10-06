#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import defaultdict, Counter

ANNOT_DIR = Path("/mnt/phoenix_2/12.annotation/prokka_all")
META_FILE = Path("/mnt/phoenix_2/11.metadata/metadata/final_genome_metadata_classified.tsv")

acc_re = re.compile(r"(GC[AF]_\d+\.\d+)")

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

    raise RuntimeError("Could not detect accession column.")

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
            if "clinical" in val or "environmental" in val or "unknown" in val or "other" in val:
                hits += 1

        if hits > best_hits:
            best_hits = hits
            best_col = h

    if best_col and best_hits > 10:
        return best_col

    return None

def norm_class(x):
    x0 = str(x).strip()
    x = x0.lower().replace("-", "_").replace(" ", "_")
    if "clinical" in x:
        return "Clinical"
    if "environmental" in x:
        return "Environmental"
    if "unknown" in x or "other" in x:
        return "Unknown_Other"
    return x0

with META_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    headers = reader.fieldnames

acc_col = detect_accession_column(headers, rows)
class_col = detect_class_column(headers, rows)

metadata_accs = []
metadata_by_acc = {}

for r in rows:
    m = acc_re.search(str(r.get(acc_col, "")))
    if m:
        acc = m.group(1)
        metadata_accs.append(acc)
        metadata_by_acc[acc] = r

gff_by_acc = {}
gff_by_suffix = defaultdict(list)

for gff in ANNOT_DIR.rglob("*.gff"):
    m = acc_re.search(gff.name) or acc_re.search(str(gff))
    if not m:
        continue

    acc = m.group(1)
    suffix = acc[4:]

    gff_by_acc[acc] = gff
    gff_by_suffix[suffix].append((acc, gff))

used_gff_accs = set()
exact_matches = []
prefix_matches = []
missing = []

for acc in metadata_accs:
    if acc in gff_by_acc:
        exact_matches.append(acc)
        used_gff_accs.add(acc)
    else:
        suffix = acc[4:]
        candidates = gff_by_suffix.get(suffix, [])
        if len(candidates) == 1:
            source_acc, source_gff = candidates[0]
            prefix_matches.append((acc, source_acc, str(source_gff)))
            used_gff_accs.add(source_acc)
        else:
            missing.append(acc)

unused_gff_accs = sorted(set(gff_by_acc) - used_gff_accs)

print("=== BASIC COUNTS ===")
print(f"Metadata rows with accession: {len(metadata_accs)}")
print(f"Prokka GFF accessions:        {len(gff_by_acc)}")
print(f"Exact matches:                {len(exact_matches)}")
print(f"GCA/GCF prefix matches:       {len(prefix_matches)}")
print(f"Missing metadata accessions:  {len(missing)}")
print(f"Unused Prokka GFF accessions: {len(unused_gff_accs)}")

print("\n=== GCA/GCF PREFIX MATCHES ===")
for meta_acc, source_acc, source_gff in prefix_matches:
    print(f"{meta_acc}\t{source_acc}\t{source_gff}")

print("\n=== MISSING METADATA ACCESSIONS ===")
for acc in missing:
    r = metadata_by_acc.get(acc, {})
    cls = norm_class(r.get(class_col, "")) if class_col else ""
    print(f"\n{acc}")
    print(f"  classification: {cls}")
    for col in [
        "organism_name",
        "species_taxid",
        "biosample_accession",
        "BioSample",
        "isolation_source",
        "host",
        "geo_loc_name",
        "collection_date",
        "assembly_level",
        "refseq_category",
    ]:
        if col in r and str(r[col]).strip():
            print(f"  {col}: {r[col]}")

print("\n=== UNUSED PROKKA GFF ACCESSIONS ===")
for acc in unused_gff_accs:
    print(f"{acc}\t{gff_by_acc[acc]}")

print("\n=== NEARBY SEARCH FOR 021497 ===")
for acc, gff in sorted(gff_by_acc.items()):
    if "021497" in acc:
        print(f"{acc}\t{gff}")
