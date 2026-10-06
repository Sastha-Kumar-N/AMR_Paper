#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict

META_FILE = Path("/mnt/phoenix_2/11.metadata/metadata/final_genome_metadata_classified.tsv")
ANNOT_DIR = Path("/mnt/phoenix_2/12.annotation/prokka_all")

acc_re = re.compile(r"(GC[AF]_\d+\.\d+)")
root_re = re.compile(r"GC[AF]_(\d+)\.(\d+)")

def root_without_version(acc):
    m = root_re.search(acc)
    if not m:
        return None
    return m.group(1)

def suffix_with_version(acc):
    return acc[4:]

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
    return None

with META_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    headers = reader.fieldnames

acc_col = detect_accession_column(headers, rows)
class_col = detect_class_column(headers, rows)

meta_accs = []
meta_rows_by_acc = defaultdict(list)
meta_by_root = defaultdict(list)

for i, r in enumerate(rows, start=2):
    m = acc_re.search(str(r.get(acc_col, "")))
    if not m:
        continue
    acc = m.group(1)
    meta_accs.append(acc)
    meta_rows_by_acc[acc].append((i, r))
    root = root_without_version(acc)
    if root:
        meta_by_root[root].append((acc, i, r))

gff_files = list(ANNOT_DIR.rglob("*.gff"))
gff_by_acc = defaultdict(list)
gff_by_root = defaultdict(list)
unparsed_gff = []

for gff in gff_files:
    m = acc_re.search(gff.name) or acc_re.search(str(gff))
    if not m:
        unparsed_gff.append(str(gff))
        continue
    acc = m.group(1)
    gff_by_acc[acc].append(gff)
    root = root_without_version(acc)
    if root:
        gff_by_root[root].append((acc, gff))

print("=== TOTALS ===")
print(f"Metadata accession rows:       {len(meta_accs)}")
print(f"Metadata unique accessions:    {len(meta_rows_by_acc)}")
print(f"Prokka GFF files:              {len(gff_files)}")
print(f"Prokka parsed GFF files:       {sum(len(v) for v in gff_by_acc.values())}")
print(f"Prokka unique GFF accessions:  {len(gff_by_acc)}")
print(f"Unparsed GFF files:            {len(unparsed_gff)}")

print("\n=== DUPLICATE METADATA ACCESSIONS ===")
dup_meta = {acc: vals for acc, vals in meta_rows_by_acc.items() if len(vals) > 1}
if not dup_meta:
    print("None")
else:
    for acc, vals in dup_meta.items():
        print(f"{acc}\trows={','.join(str(v[0]) for v in vals)}")

print("\n=== DUPLICATE PROKKA GFF ACCESSIONS ===")
dup_gff = {acc: paths for acc, paths in gff_by_acc.items() if len(paths) > 1}
if not dup_gff:
    print("None")
else:
    for acc, paths in dup_gff.items():
        print(f"{acc}\tcount={len(paths)}")
        for p in paths:
            print(f"  {p}")

print("\n=== UNPARSED GFF FILES ===")
if not unparsed_gff:
    print("None")
else:
    for p in unparsed_gff:
        print(p)

targets = [
    "GCF_019738995.2",
    "GCA_019738995.2",
    "GCF_021497345.1",
    "GCA_021497345.2",
    "GCF_021497365.1",
    "GCA_021497365.2",
]

print("\n=== TARGET ACCESSION PRESENCE ===")
for acc in targets:
    print(f"\n{acc}")
    print(f"  in metadata: {'YES' if acc in meta_rows_by_acc else 'NO'}")
    if acc in meta_rows_by_acc:
        for line_no, r in meta_rows_by_acc[acc]:
            cls = r.get(class_col, "") if class_col else ""
            biosample = r.get("biosample_accession", r.get("BioSample", ""))
            host = r.get("host", "")
            source = r.get("isolation_source", "")
            print(f"    row={line_no}; class={cls}; biosample={biosample}; host={host}; source={source}")
    print(f"  in Prokka GFF: {'YES' if acc in gff_by_acc else 'NO'}")
    if acc in gff_by_acc:
        for p in gff_by_acc[acc]:
            print(f"    {p}")

print("\n=== ROOT-LEVEL CANDIDATES FOR MISSING ACCESSIONS ===")
for acc in ["GCF_019738995.2", "GCF_021497345.1", "GCF_021497365.1"]:
    root = root_without_version(acc)
    print(f"\n{acc}; root={root}")
    print("  metadata same root:")
    for meta_acc, line_no, r in meta_by_root.get(root, []):
        cls = r.get(class_col, "") if class_col else ""
        biosample = r.get("biosample_accession", r.get("BioSample", ""))
        host = r.get("host", "")
        source = r.get("isolation_source", "")
        print(f"    {meta_acc}; row={line_no}; class={cls}; biosample={biosample}; host={host}; source={source}")
    print("  GFF same root:")
    for gff_acc, gff in gff_by_root.get(root, []):
        print(f"    {gff_acc}; {gff}")
