#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import re
import sys

EXTENSIONS = [
    ".fasta.gz", ".fa.gz", ".fna.gz", ".faa.gz", ".gff.gz",
    ".fasta", ".fa", ".fna", ".faa", ".gff",
    ".gbk", ".gbff"
]

def clean_name(x):
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return ""

    name = Path(x).name.strip()

    for ext in EXTENSIONS:
        if name.endswith(ext):
            name = name[:-len(ext)]
            break

    return name.strip()

def strip_version(x):
    x = str(x).strip()
    return re.sub(r"\.\d+$", "", x)

def add_token(tokens, x):
    x = clean_name(x)
    if not x:
        return

    tokens.add(x)
    tokens.add(strip_version(x))

    # Capture GCF/GCA accession inside longer filenames
    for m in re.findall(r"(GC[AF]_\d+\.\d+)", x):
        tokens.add(m)
        tokens.add(strip_version(m))

def choose_class_col(cols):
    # Important: force isolate_class first
    for target in ["isolate_class", "isolate_category", "classification"]:
        for c in cols:
            if c.lower() == target:
                return c

    # Fallback only if isolate_class is absent
    for target in ["sample_type", "source_type"]:
        for c in cols:
            if c.lower() == target:
                return c

    return None

def load_clinical_tokens(metadata):
    clinical_tokens = set()
    clinical_rows = 0

    with open(metadata, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cols = reader.fieldnames or []

        class_col = choose_class_col(cols)
        if class_col is None:
            sys.exit(f"ERROR: Could not find classification column. Columns found: {cols}")

        print(f"Using classification column: {class_col}")

        for row in reader:
            cls = str(row.get(class_col, "")).strip().lower()

            if "clinical" not in cls:
                continue

            clinical_rows += 1

            # Preferred metadata columns
            preferred_cols = [
                "assembly_accession",
                "genome_file",
                "genome_path",
                "file_stem",
                "metadata_file_stem",
                "prokka_stem",
                "sample",
                "sample_id",
                "strain",
                "assembly_name",
                "biosample_accession",
            ]

            for c in preferred_cols:
                if c in row:
                    add_token(clinical_tokens, row.get(c, ""))

            # Also scan all row values for embedded GCF/GCA accessions
            for v in row.values():
                if v:
                    for m in re.findall(r"(GC[AF]_\d+\.\d+)", str(v)):
                        add_token(clinical_tokens, m)

    print(f"Clinical rows in metadata: {clinical_rows}")
    print(f"Clinical matching tokens collected: {len(clinical_tokens)}")

    return clinical_tokens, clinical_rows

def tokens_match(prokka_tokens, clinical_tokens):
    for p in prokka_tokens:
        if not p:
            continue

        p0 = p
        p1 = strip_version(p)

        if p0 in clinical_tokens or p1 in clinical_tokens:
            return True

        # Prefix match for names like GCF_000000000.1_ASM...
        for c in clinical_tokens:
            if not c:
                continue
            if p0.startswith(c) or c.startswith(p0):
                return True
            if p1.startswith(c) or c.startswith(p1):
                return True

    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prokka-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--missing-out", required=True)
    ap.add_argument("--all-prokka-out", default=None)
    args = ap.parse_args()

    prokka_dir = Path(args.prokka_dir)

    if not prokka_dir.exists():
        sys.exit(f"ERROR: Prokka directory does not exist: {prokka_dir}")

    clinical_tokens, clinical_rows = load_clinical_tokens(args.metadata)

    rows = []
    all_prokka = []
    incomplete = []

    faa_files = sorted(prokka_dir.rglob("*.faa"))

    print(f"Total .faa files found under Prokka directory: {len(faa_files)}")

    for faa in faa_files:
        folder = faa.parent
        sample = faa.stem

        same_fna = folder / f"{sample}.fna"
        same_gff = folder / f"{sample}.gff"

        fna_files = [same_fna] if same_fna.exists() else sorted(folder.glob("*.fna"))
        gff_files = [same_gff] if same_gff.exists() else sorted(folder.glob("*.gff"))

        if not fna_files or not gff_files:
            incomplete.append([sample, str(folder), "missing .fna or .gff"])
            continue

        fna = fna_files[0]
        gff = gff_files[0]

        prokka_tokens = set()
        for x in [
            sample,
            folder.name,
            faa.name,
            fna.name,
            gff.name,
            str(folder),
            str(faa),
            str(fna),
            str(gff),
        ]:
            add_token(prokka_tokens, x)

        is_clinical = tokens_match(prokka_tokens, clinical_tokens)

        all_prokka.append([
            sample,
            str(fna),
            str(faa),
            str(gff),
            "clinical" if is_clinical else "not_matched_to_clinical"
        ])

        if is_clinical:
            rows.append([sample, str(fna), str(faa), str(gff)])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample", "nucleotide_fasta", "protein_faa", "gff"])
        writer.writerows(rows)

    with open(args.missing_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample", "folder", "reason"])
        writer.writerows(incomplete)

    if args.all_prokka_out:
        with open(args.all_prokka_out, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["sample", "nucleotide_fasta", "protein_faa", "gff", "match_status"])
            writer.writerows(all_prokka)

    print()
    print(f"Manifest written: {args.out}")
    print(f"Samples in clinical manifest: {len(rows)}")
    print(f"Skipped incomplete Prokka folders: {len(incomplete)}")
    print(f"Expected clinical rows from metadata: {clinical_rows}")

    if len(rows) != clinical_rows:
        print()
        print("WARNING: Manifest count does not match clinical metadata count.")
        print("This means Prokka sample names and metadata names still do not fully match.")
        print("Inspect the diagnostic file if supplied with --all-prokka-out.")

if __name__ == "__main__":
    main()
