from pathlib import Path
import re
import pandas as pd

GENOMES_DIR = Path("/mnt/phoenix_2/10.gtdbtk_filtered_genomes/pseudomonas_aeruginosa_filtered_final_list_after_gtdbk")
OUT = Path("accessions/genome_accession_table.tsv")

rows = []

for f in sorted(GENOMES_DIR.glob("*")):
    if f.suffix.lower() not in [".fna", ".fa", ".fasta"]:
        continue

    accession = None

    m = re.search(r"(GC[AF]_\d+\.\d+)", f.name)
    if m:
        accession = m.group(1)
    else:
        with open(f) as handle:
            for line in handle:
                if line.startswith(">"):
                    m = re.search(r"(GC[AF]_\d+\.\d+)", line)
                    if m:
                        accession = m.group(1)
                    break

    rows.append({
        "genome_file": f.name,
        "genome_path": str(f),
        "assembly_accession": accession
    })

df = pd.DataFrame(rows)
df.to_csv(OUT, sep="\t", index=False)

print(df["assembly_accession"].notna().sum(), "accessions found from", len(df), "genomes")
print("Saved:", OUT)
