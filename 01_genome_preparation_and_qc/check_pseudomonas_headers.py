#!/usr/bin/env python3

from pathlib import Path
import csv

NCBI_DIR = Path(".")
TARGET = "pseudomonas aeruginosa"

report_file = "pseudomonas_header_check_report.csv"

def fasta_headers(file_path):
    with open(file_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                yield line.strip()

results = []

for fasta in sorted(NCBI_DIR.glob("*.fna")):
    headers = list(fasta_headers(fasta))

    if not headers:
        status = "NO_FASTA_HEADER_FOUND"
        first_header = ""
        matched_headers = 0
        total_headers = 0
    else:
        total_headers = len(headers)
        matched_headers = sum(TARGET in h.lower() for h in headers)
        first_header = headers[0]

        if matched_headers == total_headers:
            status = "PASS_ALL_HEADERS_MATCH"
        elif matched_headers > 0:
            status = "PARTIAL_MATCH"
        else:
            status = "FAIL_NO_MATCH"

    results.append([
        fasta.name,
        status,
        total_headers,
        matched_headers,
        first_header
    ])

with open(report_file, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "file",
        "status",
        "total_headers",
        "matched_headers",
        "first_header"
    ])
    writer.writerows(results)

total_files = len(results)
passed = sum(1 for r in results if r[1] == "PASS_ALL_HEADERS_MATCH")
failed = [r for r in results if r[1] != "PASS_ALL_HEADERS_MATCH"]

print(f"Checked files: {total_files}")
print(f"Passed: {passed}")
print(f"Failed/partial/no header: {len(failed)}")
print(f"Report saved: {report_file}")

if failed:
    print("\nProblem files:")
    for r in failed:
        print(f"{r[0]} --> {r[1]}")
else:
    print("\nAll files appear to be Pseudomonas aeruginosa based on FASTA headers.")
