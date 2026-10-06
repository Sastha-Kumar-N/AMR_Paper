#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import os
import sys
import shutil
from collections import Counter, defaultdict

ANNOT_DIR = Path("/mnt/phoenix_2/12.annotation/prokka_all")
META_FILE = Path("/mnt/phoenix_2/11.metadata/metadata/final_genome_metadata_classified.tsv")
OUT_DIR = Path("/mnt/phoenix_2/13.pangenome_inputs")

OUT_WHOLE = OUT_DIR / "whole"
OUT_CLINICAL = OUT_DIR / "clinical"
OUT_ENV = OUT_DIR / "environmental"
OUT_UNKNOWN = OUT_DIR / "unknown_other"

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

    raise RuntimeError("Could not detect classification column.")

def clean_output_dirs():
    for d in [OUT_WHOLE, OUT_CLINICAL, OUT_ENV, OUT_UNKNOWN]:
        d.mkdir(parents=True, exist_ok=True)
        for f in d.glob("*.gff"):
            f.unlink()

def extract_fna_stems_from_row(row):
    """
    Extract original genome file stems from any metadata field containing .fna/.fa/.fasta.
    This is needed for GFFs named like:
    Pseudomonas_aeruginosa_SE5381_16704.gff
    instead of accession-style names.
    """
    stems = []

    for value in row.values():
        value = str(value).strip()
        if not value:
            continue

        parts = re.split(r"\s+", value)

        for part in parts:
            if part.endswith((".fna", ".fa", ".fasta")):
                stems.append(Path(part).stem)

    # Preserve order but remove duplicates
    seen = set()
    unique = []
    for s in stems:
        if s not in seen:
            unique.append(s)
            seen.add(s)

    return unique

def index_gff_files(annot_dir):
    by_accession = defaultdict(list)
    by_suffix = defaultdict(list)
    by_stem = defaultdict(list)

    for gff in annot_dir.rglob("*.gff"):
        gff_resolved = gff.resolve()

        # Index by file stem and parent directory name
        by_stem[gff.stem].append(gff_resolved)
        by_stem[gff.parent.name].append(gff_resolved)

        m = acc_re.search(gff.name) or acc_re.search(str(gff))

        if m:
            acc = m.group(1)
            suffix = acc[4:]
            by_accession[acc].append(gff_resolved)
            by_suffix[suffix].append((acc, gff_resolved))

    return by_accession, by_suffix, by_stem

def choose_single(candidates, label):
    unique = []
    seen = set()

    for c in candidates:
        c = Path(c).resolve()
        if str(c) not in seen:
            unique.append(c)
            seen.add(str(c))

    if len(unique) == 1:
        return unique[0]

    if len(unique) > 1:
        raise RuntimeError(f"Ambiguous mapping for {label}: {unique}")

    return None

def resolve_gff(meta_acc, row, by_accession, by_suffix, by_stem):
    """
    Resolution order:
    1. Exact accession match.
    2. Original FASTA stem match from metadata.
    3. GCA/GCF same accession suffix match.

    The FASTA-stem step is essential because three valid Prokka outputs are named:
    Pseudomonas_aeruginosa_P8W_16354
    Pseudomonas_aeruginosa_SE5381_16704
    Pseudomonas_aeruginosa_SE5418_16705
    """
    exact = choose_single(by_accession.get(meta_acc, []), f"exact accession {meta_acc}")
    if exact:
        return meta_acc, exact, "exact_accession"

    for stem in extract_fna_stems_from_row(row):
        stem_match = choose_single(by_stem.get(stem, []), f"file stem {stem}")
        if stem_match:
            return stem, stem_match, "metadata_file_stem"

    suffix = meta_acc[4:]
    suffix_candidates = by_suffix.get(suffix, [])

    if len(suffix_candidates) == 1:
        source_acc, source_gff = suffix_candidates[0]
        return source_acc, source_gff, "gca_gcf_same_suffix"

    if len(suffix_candidates) > 1:
        raise RuntimeError(f"Ambiguous GCA/GCF suffix mapping for {meta_acc}: {suffix_candidates}")

    return "", None, "missing"

def link_file(src, dest):
    if dest.exists() or dest.is_symlink():
        dest.unlink()

    try:
        os.link(src, dest)
        return "hardlink"
    except OSError:
        shutil.copy2(src, dest)
        return "copy"

if not ANNOT_DIR.exists():
    raise RuntimeError(f"Annotation directory not found: {ANNOT_DIR}")

if not META_FILE.exists():
    raise RuntimeError(f"Metadata file not found: {META_FILE}")

clean_output_dirs()

with META_FILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    headers = reader.fieldnames

if not rows:
    raise RuntimeError("Metadata file has no rows.")

acc_col = detect_accession_column(headers, rows)
class_col = detect_class_column(headers, rows)

by_accession, by_suffix, by_stem = index_gff_files(ANNOT_DIR)

print(f"Metadata file: {META_FILE}")
print(f"Annotation directory: {ANNOT_DIR}")
print(f"Detected accession column: {acc_col}")
print(f"Detected classification column: {class_col}")
print(f"Metadata rows: {len(rows)}")
print(f"GFF files in annotation directory: {len(list(ANNOT_DIR.rglob('*.gff')))}")

summary = Counter()
missing = []
match_type_counts = Counter()
used_source_gffs = Counter()

with manifest_file.open("w", newline="") as mf:
    writer = csv.writer(mf, delimiter="\t")
    writer.writerow([
        "metadata_accession",
        "source_identifier",
        "classification",
        "source_gff",
        "whole_file",
        "group_file",
        "match_type",
        "file_mode",
        "status"
    ])

    for r in rows:
        raw_acc = str(r.get(acc_col, "")).strip()
        m = acc_re.search(raw_acc)

        if not m:
            continue

        meta_acc = m.group(1)
        cls = normalize_class(r.get(class_col, ""))

        source_id, src_gff, match_type = resolve_gff(meta_acc, r, by_accession, by_suffix, by_stem)

        if src_gff is None:
            missing.append(meta_acc)
            writer.writerow([meta_acc, "", cls, "", "", "", match_type, "", "MISSING_GFF"])
            summary["missing_gff"] += 1
            continue

        used_source_gffs[str(src_gff)] += 1
        match_type_counts[match_type] += 1

        whole_file = OUT_WHOLE / f"{meta_acc}.gff"
        file_mode = link_file(src_gff, whole_file)
        summary["whole"] += 1

        if cls == "Clinical":
            group_file = OUT_CLINICAL / f"{meta_acc}.gff"
            link_file(src_gff, group_file)
            summary["clinical"] += 1

        elif cls == "Environmental":
            group_file = OUT_ENV / f"{meta_acc}.gff"
            link_file(src_gff, group_file)
            summary["environmental"] += 1

        else:
            group_file = OUT_UNKNOWN / f"{meta_acc}.gff"
            link_file(src_gff, group_file)
            summary["unknown_other"] += 1

        writer.writerow([
            meta_acc,
            source_id,
            cls,
            src_gff,
            whole_file,
            group_file,
            match_type,
            file_mode,
            "OK"
        ])

duplicate_source_use = {k: v for k, v in used_source_gffs.items() if v > 1}

with summary_file.open("w") as sf:
    sf.write("dataset\tgff_count\n")
    for k in ["whole", "clinical", "environmental", "unknown_other", "missing_gff"]:
        sf.write(f"{k}\t{summary[k]}\n")

    sf.write("\nmatch_type\tcount\n")
    for k, v in sorted(match_type_counts.items()):
        sf.write(f"{k}\t{v}\n")

    sf.write(f"\nduplicate_source_gff_use\t{len(duplicate_source_use)}\n")

print("\nPangenome input generation complete.")
print(f"Manifest: {manifest_file}")
print(f"Summary:  {summary_file}")

print("\nCounts:")
for k in ["whole", "clinical", "environmental", "unknown_other", "missing_gff"]:
    print(f"{k}: {summary[k]}")

print("\nMatch types:")
for k, v in sorted(match_type_counts.items()):
    print(f"{k}: {v}")

if duplicate_source_use:
    print("\nERROR: Some source GFF files were used more than once:")
    for src, n in duplicate_source_use.items():
        print(f"{n}x {src}")
    sys.exit(1)

if missing:
    print("\nERROR: Missing GFF files:")
    for acc in missing:
        print(acc)
    sys.exit(1)

expected_total = summary["clinical"] + summary["environmental"] + summary["unknown_other"]

if expected_total != summary["whole"]:
    raise RuntimeError(
        f"Group totals do not equal whole total: groups={expected_total}, whole={summary['whole']}"
    )

print("\nValidation passed.")
print(f"Group total = {expected_total}")
print(f"Whole total = {summary['whole']}")
print("No duplicate source GFF usage detected.")
