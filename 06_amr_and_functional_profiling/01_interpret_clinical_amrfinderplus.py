#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from io import StringIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def genome_from_filename(path: Path) -> str:
    name = path.name
    for suffix in [
        ".amrfinder.tsv",
        ".mutation_all.tsv",
        ".mutation.tsv",
        ".tsv",
        ".amr_proteins.faa",
        ".amr_regions.fna",
        ".faa",
        ".fna",
    ]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def read_amrfinder_table(path: Path) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = [line for line in handle if line.strip()]

    if not lines:
        return pd.DataFrame()

    header = lines[0].rstrip("\n")
    if header.startswith("#"):
        header = header[1:]

    columns = header.split("\t")
    body = "".join(lines[1:])

    if not body.strip():
        return pd.DataFrame(columns=columns)

    return pd.read_csv(
        StringIO(body),
        sep="\t",
        names=columns,
        dtype=str,
        keep_default_na=False,
    )


def find_col(df: pd.DataFrame, candidates):
    exact = {c.lower().strip(): c for c in df.columns}

    for cand in candidates:
        key = cand.lower().strip()
        if key in exact:
            return exact[key]

    for col in df.columns:
        low = col.lower().strip()
        for cand in candidates:
            words = cand.lower().strip().split()
            if all(w in low for w in words):
                return col

    return None


def series(df: pd.DataFrame, col, default=""):
    if col is None:
        return pd.Series([default] * len(df), index=df.index, dtype=str)
    return df[col].astype(str)


def standardise(df: pd.DataFrame, path: Path, source_kind: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    cols = {
        "name": find_col(df, ["Name", "Genome", "Sample"]),
        "protein_id": find_col(df, ["Protein id", "Protein identifier"]),
        "contig": find_col(df, ["Contig id", "Contig identifier", "Contig"]),
        "start": find_col(df, ["Start"]),
        "stop": find_col(df, ["Stop", "End"]),
        "strand": find_col(df, ["Strand"]),
        "element_symbol": find_col(df, ["Element symbol", "Gene symbol", "Gene"]),
        "element_name": find_col(df, ["Element name", "Sequence name", "Name of closest sequence"]),
        "scope": find_col(df, ["Scope"]),
        "type": find_col(df, ["Type", "Element type"]),
        "subtype": find_col(df, ["Subtype", "Element subtype"]),
        "drug_class": find_col(df, ["Class", "Drug class", "Antimicrobial class"]),
        "drug_subclass": find_col(df, ["Subclass", "Drug subclass"]),
        "method": find_col(df, ["Method"]),
        "coverage": find_col(df, ["Coverage", "% Coverage of reference sequence", "% Coverage"]),
        "identity": find_col(df, ["Identity", "% Identity to reference sequence", "% Identity"]),
        "alignment_length": find_col(df, ["Alignment length"]),
    }

    out = pd.DataFrame(index=df.index)

    if cols["name"] is None:
        out["genome"] = genome_from_filename(path)
    else:
        out["genome"] = series(df, cols["name"]).replace("", genome_from_filename(path))

    for c in [
        "protein_id", "contig", "start", "stop", "strand",
        "element_symbol", "element_name", "scope", "type", "subtype",
        "drug_class", "drug_subclass", "method", "coverage",
        "identity", "alignment_length",
    ]:
        out[c] = series(df, cols[c])

    out["source_file"] = str(path)
    out["source_kind"] = source_kind

    for c in ["start", "stop", "coverage", "identity", "alignment_length"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    return out


def load_folder(folder: Path, pattern: str, source_kind: str):
    frames = []
    genomes = []

    for path in sorted(folder.glob(pattern)):
        genomes.append(genome_from_filename(path))
        df = read_amrfinder_table(path)
        std = standardise(df, path, source_kind)
        if not std.empty:
            frames.append(std)

    if frames:
        return pd.concat(frames, ignore_index=True), genomes

    return pd.DataFrame(), genomes


def clean_text(x):
    x = str(x).strip()
    if x == "" or x.lower() in {"na", "nan", "none", "-"}:
        return "Unclassified"
    return x


def gene_family(symbol):
    s = str(symbol).strip()
    if not s:
        return "Unclassified"

    # Convert mutation names like gyrA_T83I, mexR_R83C, parC_S87L to gene family.
    s = re.sub(r"_[A-Z]\d+[A-Z*]+$", "", s)
    return s


def remove_non_resistance_rows(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    x = df.copy()

    if "type" in x.columns:
        type_upper = x["type"].astype(str).str.upper().str.strip()
        if type_upper.isin(["AMR", "STRESS", "VIRULENCE"]).any():
            x = x[type_upper.eq("AMR")].copy()

    joined = (
        x["element_symbol"].astype(str) + " " +
        x["element_name"].astype(str) + " " +
        x["scope"].astype(str) + " " +
        x["subtype"].astype(str) + " " +
        x["method"].astype(str)
    )

    non_resistant = joined.str.contains("WILDTYPE|UNKNOWN", case=False, regex=True, na=False)
    x = x[~non_resistant].copy()

    x["element_symbol"] = x["element_symbol"].replace("", np.nan)
    x["element_name"] = x["element_name"].replace("", np.nan)

    x["gene_exact"] = x["element_symbol"].fillna(x["element_name"]).fillna("Unclassified")
    x["gene_exact"] = x["gene_exact"].astype(str).str.strip().replace("", "Unclassified")

    x["gene_family"] = x["gene_exact"].apply(gene_family)
    x["drug_class_clean"] = x["drug_class"].apply(clean_text)
    x["drug_subclass_clean"] = x["drug_subclass"].apply(clean_text)
    x["method_clean"] = x["method"].apply(clean_text)

    return x


def save_barh(df, label_col, value_col, title, xlabel, outfile, top_n=30):
    if df.empty:
        return

    plot_df = df.head(top_n).copy()
    plot_df = plot_df.sort_values(value_col, ascending=True)

    plt.figure(figsize=(10, max(4, 0.35 * len(plot_df))))
    plt.barh(plot_df[label_col].astype(str), plot_df[value_col])
    plt.xlabel(xlabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.savefig(str(outfile).replace(".png", ".pdf"))
    plt.close()


def write_itol_binary(presence, top_genes, outfile: Path):
    if presence.empty or not top_genes:
        return

    labels = [re.sub(r"[^A-Za-z0-9_.()+-]", "_", str(g)) for g in top_genes]

    with open(outfile, "w") as out:
        out.write("DATASET_BINARY\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATASET_LABEL\tTop_AMR_determinants\n")
        out.write("COLOR\t#d73027\n")
        out.write("FIELD_SHAPES\t" + "\t".join(["1"] * len(labels)) + "\n")
        out.write("FIELD_LABELS\t" + "\t".join(labels) + "\n")
        out.write("FIELD_COLORS\t" + "\t".join(["#d73027"] * len(labels)) + "\n")
        out.write("DATA\n")

        for genome, row in presence[top_genes].iterrows():
            out.write(str(genome) + "\t" + "\t".join(str(int(row[g])) for g in top_genes) + "\n")


def write_itol_simplebar(summary, outfile: Path):
    if summary.empty:
        return

    with open(outfile, "w") as out:
        out.write("DATASET_SIMPLEBAR\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATASET_LABEL\tAMR_burden\n")
        out.write("COLOR\t#4575b4\n")
        out.write("DATA\n")

        for _, row in summary.iterrows():
            out.write(f"{row['genome']}\t{row['unique_amr_elements']}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="clinical_amrfinderplus result folder")
    parser.add_argument("--outdir", required=True, help="AMR interpretation output folder")
    parser.add_argument("--top-n", type=int, default=40)
    args = parser.parse_args()

    input_dir = Path(args.input)
    outdir = Path(args.outdir)

    tables = outdir / "tables"
    figures = outdir / "figures"
    itol = outdir / "itol"

    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    itol.mkdir(parents=True, exist_ok=True)

    raw_dir = input_dir / "raw"
    mut_dir = input_dir / "mutations_all"

    raw, raw_genomes = load_folder(raw_dir, "*.amrfinder.tsv", "raw_amrfinder")
    mutations_all, mutation_genomes = load_folder(mut_dir, "*.mutation_all.tsv", "mutation_all")

    all_genomes = sorted(set(raw_genomes) | set(mutation_genomes) | set(raw.get("genome", pd.Series(dtype=str)).astype(str)))

    raw.to_csv(tables / "combined_raw_amrfinder.tsv", sep="\t", index=False)
    mutations_all.to_csv(tables / "combined_mutation_all.tsv", sep="\t", index=False)

    amr = remove_non_resistance_rows(raw)
    amr.to_csv(tables / "combined_raw_amr_only_filtered.tsv", sep="\t", index=False)

    n_genomes = len(all_genomes)

    if amr.empty:
        raise SystemExit("No AMR-positive rows found after filtering. Check input raw/*.amrfinder.tsv files.")

    element_prev = (
        amr.groupby("gene_exact")
        .agg(
            total_hit_rows=("gene_exact", "size"),
            genomes_with_hit=("genome", "nunique"),
            gene_family=("gene_family", lambda x: ";".join(sorted(set(map(str, x))))),
            element_names=("element_name", lambda x: "; ".join(sorted(set(map(str, x)))[:8])),
            drug_classes=("drug_class_clean", lambda x: "; ".join(sorted(set(map(str, x))))),
            drug_subclasses=("drug_subclass_clean", lambda x: "; ".join(sorted(set(map(str, x))))),
            methods=("method_clean", lambda x: "; ".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )

    element_prev["prevalence_percent"] = 100 * element_prev["genomes_with_hit"] / max(n_genomes, 1)
    element_prev = element_prev.sort_values(["genomes_with_hit", "total_hit_rows"], ascending=False)
    element_prev.to_csv(tables / "amr_element_prevalence.tsv", sep="\t", index=False)

    family_prev = (
        amr.groupby("gene_family")
        .agg(
            total_hit_rows=("gene_family", "size"),
            genomes_with_hit=("genome", "nunique"),
            exact_elements=("gene_exact", lambda x: "; ".join(sorted(set(map(str, x)))[:25])),
            drug_classes=("drug_class_clean", lambda x: "; ".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )

    family_prev["prevalence_percent"] = 100 * family_prev["genomes_with_hit"] / max(n_genomes, 1)
    family_prev = family_prev.sort_values(["genomes_with_hit", "total_hit_rows"], ascending=False)
    family_prev.to_csv(tables / "amr_gene_family_prevalence.tsv", sep="\t", index=False)

    class_prev = (
        amr[amr["drug_class_clean"] != "Unclassified"]
        .drop_duplicates(["genome", "drug_class_clean"])
        .groupby("drug_class_clean")
        .agg(genomes_with_class=("genome", "nunique"))
        .reset_index()
        .sort_values("genomes_with_class", ascending=False)
    )

    class_prev["prevalence_percent"] = 100 * class_prev["genomes_with_class"] / max(n_genomes, 1)
    class_prev.to_csv(tables / "amr_drug_class_prevalence.tsv", sep="\t", index=False)

    summary = pd.DataFrame({"genome": all_genomes})
    grp = amr.groupby("genome")

    summary = summary.merge(grp.size().rename("amr_hit_rows").reset_index(), on="genome", how="left")
    summary = summary.merge(grp["gene_exact"].nunique().rename("unique_amr_elements").reset_index(), on="genome", how="left")
    summary = summary.merge(grp["gene_family"].nunique().rename("unique_amr_gene_families").reset_index(), on="genome", how="left")
    summary = summary.merge(grp["drug_class_clean"].nunique().rename("unique_drug_classes").reset_index(), on="genome", how="left")

    point_mask = (
        amr["method_clean"].str.contains("POINT", case=False, na=False)
        | amr["subtype"].str.contains("POINT", case=False, na=False)
    )

    point_counts = (
        amr[point_mask]
        .groupby("genome")["gene_exact"]
        .nunique()
        .rename("point_mutation_determinants")
        .reset_index()
    )

    gene_counts = (
        amr[~point_mask]
        .groupby("genome")["gene_exact"]
        .nunique()
        .rename("non_point_gene_determinants")
        .reset_index()
    )

    summary = summary.merge(point_counts, on="genome", how="left")
    summary = summary.merge(gene_counts, on="genome", how="left")

    for col in [
        "amr_hit_rows",
        "unique_amr_elements",
        "unique_amr_gene_families",
        "unique_drug_classes",
        "point_mutation_determinants",
        "non_point_gene_determinants",
    ]:
        summary[col] = summary[col].fillna(0).astype(int)

    drug_classes = (
        amr.groupby("genome")["drug_class_clean"]
        .apply(lambda x: "; ".join(sorted(set(v for v in x if v != "Unclassified"))))
        .rename("drug_classes_present")
        .reset_index()
    )

    summary = summary.merge(drug_classes, on="genome", how="left")
    summary["drug_classes_present"] = summary["drug_classes_present"].fillna("")

    summary = summary.sort_values(
        ["unique_amr_elements", "unique_drug_classes", "amr_hit_rows"],
        ascending=False,
    )

    summary.to_csv(tables / "per_genome_amr_burden.tsv", sep="\t", index=False)

    cutoff_90 = summary["unique_amr_elements"].quantile(0.90)
    high_burden = summary[summary["unique_amr_elements"] >= cutoff_90].copy()
    high_burden.to_csv(tables / "top_10_percent_high_amr_burden_genomes.tsv", sep="\t", index=False)

    presence = (
        amr.drop_duplicates(["genome", "gene_exact"])
        .assign(present=1)
        .pivot_table(index="genome", columns="gene_exact", values="present", fill_value=0, aggfunc="max")
        .astype(int)
    )

    for genome in all_genomes:
        if genome not in presence.index:
            presence.loc[genome] = 0

    presence = presence.loc[summary["genome"].tolist()]
    presence.to_csv(tables / "amr_gene_presence_absence_matrix_all.tsv", sep="\t")

    top_genes = element_prev.head(args.top_n)["gene_exact"].tolist()
    presence[top_genes].to_csv(tables / "amr_gene_presence_absence_matrix_top.tsv", sep="\t")

    edges = []
    X = presence[top_genes].values

    for i in range(len(top_genes)):
        for j in range(i + 1, len(top_genes)):
            a = X[:, i]
            b = X[:, j]
            both = int(np.sum((a == 1) & (b == 1)))
            union = int(np.sum((a == 1) | (b == 1)))
            if union == 0:
                continue

            jaccard = both / union

            if both >= max(5, int(0.01 * n_genomes)):
                edges.append(
                    {
                        "gene_a": top_genes[i],
                        "gene_b": top_genes[j],
                        "genomes_with_both": both,
                        "jaccard": jaccard,
                    }
                )

    edge_df = pd.DataFrame(edges)
    if not edge_df.empty:
        edge_df = edge_df.sort_values(["genomes_with_both", "jaccard"], ascending=False)

    edge_df.to_csv(tables / "top_amr_gene_cooccurrence_edges.tsv", sep="\t", index=False)

    save_barh(
        element_prev,
        "gene_exact",
        "prevalence_percent",
        "Top AMR determinants across clinical P. aeruginosa genomes",
        "Genome prevalence (%)",
        figures / "top_amr_determinants_prevalence.png",
        args.top_n,
    )

    save_barh(
        family_prev,
        "gene_family",
        "prevalence_percent",
        "Top AMR gene families across clinical P. aeruginosa genomes",
        "Genome prevalence (%)",
        figures / "top_amr_gene_families_prevalence.png",
        args.top_n,
    )

    save_barh(
        class_prev,
        "drug_class_clean",
        "prevalence_percent",
        "AMR drug-class prevalence",
        "Genome prevalence (%)",
        figures / "amr_drug_class_prevalence.png",
        30,
    )

    plt.figure(figsize=(8, 5))
    plt.hist(summary["unique_amr_elements"], bins=30)
    plt.xlabel("Unique AMR determinants per genome")
    plt.ylabel("Number of genomes")
    plt.title("Distribution of AMR burden across clinical genomes")
    plt.tight_layout()
    plt.savefig(figures / "per_genome_amr_burden_histogram.png", dpi=300)
    plt.savefig(figures / "per_genome_amr_burden_histogram.pdf")
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.hist(summary["unique_drug_classes"], bins=20)
    plt.xlabel("Drug classes represented per genome")
    plt.ylabel("Number of genomes")
    plt.title("Distribution of multidrug-class AMR burden")
    plt.tight_layout()
    plt.savefig(figures / "per_genome_drug_class_burden_histogram.png", dpi=300)
    plt.savefig(figures / "per_genome_drug_class_burden_histogram.pdf")
    plt.close()

    heat = presence[top_genes].T
    plt.figure(figsize=(16, max(6, 0.25 * len(top_genes))))
    plt.imshow(heat.values, aspect="auto", interpolation="nearest")
    plt.yticks(range(len(top_genes)), top_genes)
    plt.xticks([])
    plt.xlabel("Clinical genomes ordered by AMR burden")
    plt.ylabel("AMR determinants")
    plt.title("Presence/absence heatmap of top AMR determinants")
    plt.tight_layout()
    plt.savefig(figures / "top_amr_determinants_presence_absence_heatmap.png", dpi=300)
    plt.savefig(figures / "top_amr_determinants_presence_absence_heatmap.pdf")
    plt.close()

    write_itol_binary(presence, top_genes[:25], itol / "itol_top_amr_genes_binary.txt")
    write_itol_simplebar(summary, itol / "itol_amr_burden_simplebar.txt")

    report = outdir / "AMRFinderPlus_clinical_interpretation_summary.txt"

    with open(report, "w") as out:
        out.write("AMRFinderPlus clinical interpretation summary\n")
        out.write("============================================\n\n")
        out.write(f"Input folder: {input_dir}\n")
        out.write(f"Raw AMRFinder files detected: {len(raw_genomes)}\n")
        out.write(f"Mutation-all files detected: {len(mutation_genomes)}\n")
        out.write(f"Unique genomes detected: {n_genomes}\n")
        out.write(f"Raw rows loaded: {len(raw)}\n")
        out.write(f"Filtered AMR-positive rows used: {len(amr)}\n")
        out.write(f"90th percentile AMR burden cutoff: {cutoff_90}\n\n")

        out.write("Top 30 AMR determinants:\n")
        for _, row in element_prev.head(30).iterrows():
            out.write(
                f"{row['gene_exact']}\t"
                f"{row['genomes_with_hit']}/{n_genomes}\t"
                f"{row['prevalence_percent']:.2f}%\t"
                f"{row['drug_classes']}\n"
            )

        out.write("\nMain outputs:\n")
        out.write(f"Tables: {tables}\n")
        out.write(f"Figures: {figures}\n")
        out.write(f"iTOL: {itol}\n")

    print("Done.")
    print(f"Input: {input_dir}")
    print(f"Output: {outdir}")
    print(f"Tables: {tables}")
    print(f"Figures: {figures}")
    print(f"iTOL: {itol}")
    print(f"Report: {report}")


if __name__ == "__main__":
    main()
