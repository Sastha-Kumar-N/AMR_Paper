#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import sys


def genome_from_filename(path: Path) -> str:
    name = path.name
    for suffix in [".amrfinder.tsv", ".tsv", ".txt"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def read_table(path: Path):
    """
    Reads one AMRFinderPlus TSV.
    Handles:
    - normal TSV header
    - header beginning with #
    - header-only files
    - empty files
    """
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    if not lines:
        return None, []

    header = lines[0]
    if header.startswith("#"):
        header = header[1:]

    columns = header.split("\t")
    rows = []

    for line in lines[1:]:
        parts = line.split("\t")

        if len(parts) < len(columns):
            parts = parts + [""] * (len(columns) - len(parts))
        elif len(parts) > len(columns):
            parts = parts[:len(columns)]

        rows.append(parts)

    return columns, rows


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-genome AMRFinderPlus TSV files into one combined table."
    )
    parser.add_argument("--raw-dir", required=True, help="Folder containing raw/*.amrfinder.tsv files")
    parser.add_argument("--out", required=True, help="Combined output TSV")
    parser.add_argument("--pattern", default="*.amrfinder.tsv", help="Input file pattern")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out_path = Path(args.out)
    summary_path = out_path.with_suffix(out_path.suffix + ".summary.txt")

    if not raw_dir.exists():
        sys.exit(f"ERROR: raw directory does not exist: {raw_dir}")

    files = sorted(raw_dir.glob(args.pattern))

    if not files:
        sys.exit(f"ERROR: no files matched {args.pattern} in {raw_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    master_columns = None
    total_rows = 0
    files_with_hits = 0
    files_without_hits = 0
    header_mismatch_files = []

    with open(out_path, "w", encoding="utf-8", newline="") as out_handle:
        writer = None

        for path in files:
            genome = genome_from_filename(path)
            columns, rows = read_table(path)

            if columns is None:
                files_without_hits += 1
                continue

            if master_columns is None:
                master_columns = columns
                writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["genome", "source_file"] + master_columns)

            if columns != master_columns:
                header_mismatch_files.append(str(path))

                # Align by column name when possible.
                col_to_index = {c: i for i, c in enumerate(columns)}
                aligned_rows = []

                for row in rows:
                    aligned = []
                    for col in master_columns:
                        if col in col_to_index:
                            idx = col_to_index[col]
                            aligned.append(row[idx] if idx < len(row) else "")
                        else:
                            aligned.append("")
                    aligned_rows.append(aligned)

                rows = aligned_rows

            if rows:
                files_with_hits += 1
            else:
                files_without_hits += 1

            for row in rows:
                writer.writerow([genome, str(path)] + row)
                total_rows += 1

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("Combined AMRFinderPlus results summary\n")
        summary.write("=====================================\n\n")
        summary.write(f"Raw directory: {raw_dir}\n")
        summary.write(f"Input pattern: {args.pattern}\n")
        summary.write(f"Input files detected: {len(files)}\n")
        summary.write(f"Files with AMR hits: {files_with_hits}\n")
        summary.write(f"Files without AMR hits or header-only: {files_without_hits}\n")
        summary.write(f"Combined AMR rows written: {total_rows}\n")
        summary.write(f"Combined output: {out_path}\n")
        summary.write(f"Header mismatch files: {len(header_mismatch_files)}\n\n")

        if header_mismatch_files:
            summary.write("Files with non-identical headers:\n")
            for f in header_mismatch_files:
                summary.write(f"{f}\n")

    print("Done.")
    print(f"Input files detected: {len(files)}")
    print(f"Combined AMR rows written: {total_rows}")
    print(f"Output: {out_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
