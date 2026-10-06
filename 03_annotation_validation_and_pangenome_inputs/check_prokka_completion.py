import pandas as pd
from pathlib import Path

metadata = Path("/mnt/phoenix_2/11.metadata/metadata/final_genome_metadata_classified.tsv")
anno = Path("/mnt/phoenix_2/12.annotation/prokka_all")
out = Path("/mnt/phoenix_2/12.annotation/prokka_annotation_completion_summary.tsv")

df = pd.read_csv(metadata, sep="\t").fillna("")

rows = []

for _, row in df.iterrows():
    genome_file = row["genome_file"]
    stem = Path(genome_file).stem
    d = anno / stem

    rows.append({
        "genome_file": genome_file,
        "assembly_accession": row["assembly_accession"],
        "isolate_class": row["isolate_class"],
        "annotation_dir_exists": d.exists(),
        "gff_exists": (d / f"{stem}.gff").exists(),
        "gbk_exists": (d / f"{stem}.gbk").exists(),
        "faa_exists": (d / f"{stem}.faa").exists(),
        "ffn_exists": (d / f"{stem}.ffn").exists(),
        "fna_exists": (d / f"{stem}.fna").exists(),
        "tbl_exists": (d / f"{stem}.tbl").exists(),
        "sqn_exists": (d / f"{stem}.sqn").exists(),
        "err_exists": (d / f"{stem}.err").exists(),
        "log_exists": (d / f"{stem}.log").exists()
    })

res = pd.DataFrame(rows)
res.to_csv(out, sep="\t", index=False)

print("Saved:", out)
print("Expected genomes:", len(res))
print("GFF complete:", res["gff_exists"].sum())
print("GBK complete:", res["gbk_exists"].sum())
print("FAA complete:", res["faa_exists"].sum())
print("FFN complete:", res["ffn_exists"].sum())
print("Missing GFF:", (~res["gff_exists"]).sum())
print("Missing FAA:", (~res["faa_exists"]).sum())
print("Missing GBK:", (~res["gbk_exists"]).sum())

print("\nBy isolate class:")
print(res.groupby("isolate_class")[["gff_exists", "faa_exists", "gbk_exists"]].sum())
