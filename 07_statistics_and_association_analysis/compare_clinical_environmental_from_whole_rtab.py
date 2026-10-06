#!/usr/bin/env python3

import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


BASE = Path("/mnt/phoenix_2")

RTAB = BASE / "14.pangenome_panaroo/results/whole/gene_presence_absence.Rtab"
GPA_CSV = BASE / "14.pangenome_panaroo/results/whole/gene_presence_absence.csv"
META = BASE / "11.metadata/metadata/final_genome_metadata_classified.tsv"

OUTDIR = BASE / "15.pangenome_comparison/results"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_MAIN = OUTDIR / "clinical_environmental_gene_frequency_table.tsv"
OUT_CLIN_ENRICHED = OUTDIR / "clinical_enriched_genes.tsv"
OUT_ENV_ENRICHED = OUTDIR / "environmental_enriched_genes.tsv"
OUT_BOTH_CORE = OUTDIR / "conserved_in_both_core_genes.tsv"
OUT_CLIN_SPECIFIC = OUTDIR / "clinical_specific_common_genes.tsv"
OUT_ENV_SPECIFIC = OUTDIR / "environmental_specific_common_genes.tsv"
OUT_UNMATCHED = OUTDIR / "unmatched_samples.tsv"
OUT_SUMMARY = OUTDIR / "comparison_run_summary.txt"


def extract_accession(value):
    value = str(value)
    match = re.search(r"(GC[AF]_\d+\.\d+)", value)
    return match.group(1) if match else value


def bh_fdr(pvalues):
    """
    Benjamini-Hochberg FDR correction.
    """
    p = pd.Series(pvalues, dtype="float64")
    n = p.notna().sum()

    q = pd.Series([float("nan")] * len(p), index=p.index, dtype="float64")
    if n == 0:
        return q

    ranked = p.dropna().sort_values()
    ranks = pd.Series(range(1, len(ranked) + 1), index=ranked.index, dtype="float64")

    raw_q = ranked * n / ranks
    raw_q = raw_q.iloc[::-1].cummin().iloc[::-1]
    raw_q[raw_q > 1] = 1.0

    q.loc[raw_q.index] = raw_q
    return q


def assign_category(row):
    clinical_freq = row["clinical_frequency"]
    environmental_freq = row["environmental_frequency"]
    whole_freq = row["whole_frequency"]
    clinical_count = row["clinical_count"]
    environmental_count = row["environmental_count"]
    fdr = row["fisher_fdr"]
    delta = row["clinical_minus_environmental_frequency"]

    significant = pd.notna(fdr) and fdr <= 0.05

    if clinical_freq >= 0.99 and environmental_freq >= 0.99:
        return "conserved_in_both_core"

    if clinical_count > 0 and environmental_count == 0 and clinical_freq >= 0.15:
        return "clinical_specific_common"

    if environmental_count > 0 and clinical_count == 0 and environmental_freq >= 0.15:
        return "environmental_specific_common"

    if significant and delta >= 0.25 and clinical_freq >= 0.15:
        return "clinical_enriched"

    if significant and delta <= -0.25 and environmental_freq >= 0.15:
        return "environmental_enriched"

    if clinical_freq >= 0.95 and environmental_freq < 0.95:
        return "clinical_softcore_not_environmental"

    if environmental_freq >= 0.95 and clinical_freq < 0.95:
        return "environmental_softcore_not_clinical"

    if whole_freq < 0.15:
        return "rare_cloud"

    return "shared_accessory_or_variable"


print("Loading metadata...")
meta = pd.read_csv(META, sep="\t", dtype=str)

required_meta_cols = {"assembly_accession", "isolate_class"}
missing_meta_cols = required_meta_cols - set(meta.columns)
if missing_meta_cols:
    raise SystemExit(f"Missing metadata columns: {missing_meta_cols}")

meta["assembly_accession_clean"] = meta["assembly_accession"].map(extract_accession)
meta["isolate_class_clean"] = meta["isolate_class"].fillna("").str.lower()

clinical_accessions = set(
    meta.loc[meta["isolate_class_clean"].str.contains("clinical"), "assembly_accession_clean"]
)

environmental_accessions = set(
    meta.loc[meta["isolate_class_clean"].str.contains("environmental"), "assembly_accession_clean"]
)

unknown_accessions = set(
    meta.loc[meta["isolate_class_clean"].str.contains("unknown|other"), "assembly_accession_clean"]
)

print(f"Clinical metadata accessions: {len(clinical_accessions)}")
print(f"Environmental metadata accessions: {len(environmental_accessions)}")
print(f"Unknown/other metadata accessions: {len(unknown_accessions)}")

print("Loading whole pangenome Rtab...")
rtab = pd.read_csv(RTAB, sep="\t")

if "Gene" not in rtab.columns:
    raise SystemExit("Rtab file does not contain a 'Gene' column.")

matrix_samples = list(rtab.columns[1:])
matrix_sample_set = set(matrix_samples)

clinical_cols = [x for x in matrix_samples if x in clinical_accessions]
environmental_cols = [x for x in matrix_samples if x in environmental_accessions]
unknown_cols = [x for x in matrix_samples if x in unknown_accessions]

metadata_accessions = clinical_accessions | environmental_accessions | unknown_accessions
matrix_not_in_metadata = sorted(matrix_sample_set - metadata_accessions)
metadata_not_in_matrix = sorted(metadata_accessions - matrix_sample_set)

unmatched_rows = []
for sample in matrix_not_in_metadata:
    unmatched_rows.append({"sample": sample, "issue": "matrix_sample_not_in_metadata"})
for sample in metadata_not_in_matrix:
    unmatched_rows.append({"sample": sample, "issue": "metadata_sample_not_in_matrix"})

pd.DataFrame(unmatched_rows).to_csv(OUT_UNMATCHED, sep="\t", index=False)

print(f"Clinical columns matched in matrix: {len(clinical_cols)}")
print(f"Environmental columns matched in matrix: {len(environmental_cols)}")
print(f"Unknown/other columns matched in matrix: {len(unknown_cols)}")
print(f"Matrix samples not in metadata: {len(matrix_not_in_metadata)}")
print(f"Metadata samples not in matrix: {len(metadata_not_in_matrix)}")

if len(clinical_cols) != 1350:
    raise SystemExit(f"Expected 1350 clinical columns, found {len(clinical_cols)}")

if len(environmental_cols) != 99:
    raise SystemExit(f"Expected 99 environmental columns, found {len(environmental_cols)}")

all_cols = matrix_samples

print("Calculating gene frequencies...")
clinical_count = rtab[clinical_cols].sum(axis=1).astype(int)
environmental_count = rtab[environmental_cols].sum(axis=1).astype(int)
whole_count = rtab[all_cols].sum(axis=1).astype(int)

clinical_n = len(clinical_cols)
environmental_n = len(environmental_cols)
whole_n = len(all_cols)

result = pd.DataFrame({
    "gene_cluster": rtab["Gene"],
    "whole_count": whole_count,
    "whole_n": whole_n,
    "whole_frequency": whole_count / whole_n,
    "clinical_count": clinical_count,
    "clinical_n": clinical_n,
    "clinical_frequency": clinical_count / clinical_n,
    "environmental_count": environmental_count,
    "environmental_n": environmental_n,
    "environmental_frequency": environmental_count / environmental_n,
})

result["clinical_minus_environmental_frequency"] = (
    result["clinical_frequency"] - result["environmental_frequency"]
)

print("Running Fisher exact tests...")
pvalues = []
odds_ratios = []

for _, row in result.iterrows():
    c_present = int(row["clinical_count"])
    c_absent = clinical_n - c_present
    e_present = int(row["environmental_count"])
    e_absent = environmental_n - e_present

    table = [[c_present, c_absent], [e_present, e_absent]]
    odds_ratio, pvalue = fisher_exact(table, alternative="two-sided")

    odds_ratios.append(odds_ratio)
    pvalues.append(pvalue)

result["fisher_odds_ratio"] = odds_ratios
result["fisher_pvalue"] = pvalues
result["fisher_fdr"] = bh_fdr(result["fisher_pvalue"])

print("Adding gene annotations...")
try:
    ann = pd.read_csv(
        GPA_CSV,
        usecols=lambda c: c in {"Gene", "Annotation", "Non-unique Gene name"},
        dtype=str
    )
    result = result.merge(
        ann,
        left_on="gene_cluster",
        right_on="Gene",
        how="left"
    )
    result = result.drop(columns=["Gene"])
except Exception as exc:
    print(f"Warning: could not add annotations from gene_presence_absence.csv: {exc}")

print("Assigning categories...")
result["comparison_category"] = result.apply(assign_category, axis=1)

sort_cols = [
    "comparison_category",
    "fisher_fdr",
    "clinical_minus_environmental_frequency",
    "whole_frequency",
]

result = result.sort_values(
    by=sort_cols,
    ascending=[True, True, False, False],
    na_position="last"
)

result.to_csv(OUT_MAIN, sep="\t", index=False)

result.loc[result["comparison_category"] == "clinical_enriched"].to_csv(
    OUT_CLIN_ENRICHED, sep="\t", index=False
)

result.loc[result["comparison_category"] == "environmental_enriched"].to_csv(
    OUT_ENV_ENRICHED, sep="\t", index=False
)

result.loc[result["comparison_category"] == "conserved_in_both_core"].to_csv(
    OUT_BOTH_CORE, sep="\t", index=False
)

result.loc[result["comparison_category"] == "clinical_specific_common"].to_csv(
    OUT_CLIN_SPECIFIC, sep="\t", index=False
)

result.loc[result["comparison_category"] == "environmental_specific_common"].to_csv(
    OUT_ENV_SPECIFIC, sep="\t", index=False
)

summary = result["comparison_category"].value_counts().reset_index()
summary.columns = ["comparison_category", "gene_count"]

with open(OUT_SUMMARY, "w") as handle:
    handle.write("Clinical vs environmental pangenome comparison from whole Panaroo Rtab\n")
    handle.write("=====================================================================\n\n")
    handle.write(f"Whole Rtab: {RTAB}\n")
    handle.write(f"Metadata: {META}\n\n")
    handle.write(f"Total gene clusters analysed: {len(result)}\n")
    handle.write(f"Whole matrix genomes: {whole_n}\n")
    handle.write(f"Clinical genomes matched: {clinical_n}\n")
    handle.write(f"Environmental genomes matched: {environmental_n}\n")
    handle.write(f"Unknown/other genomes matched but not used for clinical-vs-environmental test: {len(unknown_cols)}\n")
    handle.write(f"Matrix samples not in metadata: {len(matrix_not_in_metadata)}\n")
    handle.write(f"Metadata samples not in matrix: {len(metadata_not_in_matrix)}\n\n")
    handle.write("Category counts:\n")
    for _, row in summary.iterrows():
        handle.write(f"{row['comparison_category']}\t{row['gene_count']}\n")

print(f"Wrote main table: {OUT_MAIN}")
print(f"Wrote summary: {OUT_SUMMARY}")
print("Done.")
