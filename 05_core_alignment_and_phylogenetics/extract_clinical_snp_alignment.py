from pathlib import Path
import csv
import re
import sys

meta = Path(sys.argv[1])
in_aln = Path(sys.argv[2])
out_aln = Path(sys.argv[3])

acc_re = re.compile(r'GC[AF]_\d+\.\d+')

# -----------------------------
# 1. Read metadata and detect clinical rows
# -----------------------------
clinical_keys = set()
clinical_rows = 0

with meta.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames or []

    for row in reader:
        values_lower = [str(v).strip().lower() for v in row.values() if v is not None]

        # A row is clinical if any metadata classification field says clinical.
        # This is intentionally broad because column names may vary.
        is_clinical = any(v == "clinical" or v.startswith("clinical") for v in values_lower)

        if not is_clinical:
            continue

        clinical_rows += 1

        for v in row.values():
            if not v:
                continue
            v = str(v).strip()
            if not v:
                continue

            clinical_keys.add(v)

            for acc in acc_re.findall(v):
                clinical_keys.add(acc)

            # Also store basename-like forms
            clinical_keys.add(Path(v).name)
            clinical_keys.add(Path(v).stem)

print(f"Clinical rows detected in metadata: {clinical_rows}")
print(f"Clinical matching keys collected: {len(clinical_keys)}")

if clinical_rows != 1350:
    print("WARNING: Metadata clinical row count is not 1350. Check classification column names.", file=sys.stderr)

# -----------------------------
# 2. Read FASTA alignment
# -----------------------------
records = []
header = None
seq_parts = []

with in_aln.open() as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts).upper()))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.strip())

if header is not None:
    records.append((header, "".join(seq_parts).upper()))

print(f"Input alignment sequences: {len(records)}")

if not records:
    raise SystemExit("No sequences found in input alignment.")

lengths = {len(seq) for _, seq in records}
if len(lengths) != 1:
    raise SystemExit(f"Alignment has inconsistent sequence lengths: {sorted(lengths)[:10]}")

aln_len = next(iter(lengths))
print(f"Input alignment length: {aln_len}")

# -----------------------------
# 3. Match clinical records
# -----------------------------
clinical_records = []
unmatched_headers = []

for h, seq in records:
    first_token = h.split()[0]
    h_stem = Path(first_token).stem
    accs = acc_re.findall(h)

    matched = False

    if first_token in clinical_keys or h in clinical_keys or h_stem in clinical_keys:
        matched = True

    if not matched and any(acc in clinical_keys for acc in accs):
        matched = True

    if matched:
        clinical_records.append((h, seq))
    else:
        unmatched_headers.append(h)

print(f"Clinical sequences extracted: {len(clinical_records)}")

if len(clinical_records) != 1350:
    print("WARNING: Extracted clinical sequence count is not 1350.", file=sys.stderr)
    debug = out_aln.with_suffix(".unmatched_headers.txt")
    with debug.open("w") as f:
        for h in unmatched_headers:
            f.write(h + "\n")
    print(f"Unmatched headers written to: {debug}", file=sys.stderr)
    raise SystemExit("Stopping because clinical alignment count is not exactly 1350.")

# -----------------------------
# 4. Remove columns invariant within clinical subset
#    This is required for GTR+ASC because ASC expects variable sites only.
# -----------------------------
valid_bases = set("ACGT")
keep_cols = []

for i in range(aln_len):
    observed = {seq[i] for _, seq in clinical_records if seq[i] in valid_bases}
    if len(observed) >= 2:
        keep_cols.append(i)

print(f"Clinical variable SNP columns retained: {len(keep_cols)}")

if len(keep_cols) == 0:
    raise SystemExit("No variable SNP columns retained.")

# -----------------------------
# 5. Write filtered clinical variable-site alignment
# -----------------------------
with out_aln.open("w") as f:
    for h, seq in clinical_records:
        new_seq = "".join(seq[i] for i in keep_cols)
        f.write(f">{h}\n")
        for j in range(0, len(new_seq), 80):
            f.write(new_seq[j:j+80] + "\n")

print(f"Written: {out_aln}")
