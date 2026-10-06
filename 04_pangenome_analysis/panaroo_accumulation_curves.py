#!/usr/bin/env python3

import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def safe_label(label):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label)


def read_rtab(path):
    df = pd.read_csv(path, sep="\t")

    gene_col = df.columns[0]
    genome_cols = list(df.columns[1:])

    matrix = df[genome_cols].to_numpy()
    matrix = matrix > 0

    return gene_col, genome_cols, matrix


def compute_curves(matrix, permutations=100, seed=42):
    """
    matrix shape: genes x genomes
    values: True/False for presence/absence
    """
    rng = np.random.default_rng(seed)

    n_genes, n_genomes = matrix.shape

    pan_curves = np.zeros((permutations, n_genomes), dtype=np.int32)
    core_curves = np.zeros((permutations, n_genomes), dtype=np.int32)

    genome_counts = np.arange(1, n_genomes + 1, dtype=np.uint16)

    for i in range(permutations):
        order = rng.permutation(n_genomes)
        shuffled = matrix[:, order].astype(np.uint8)

        cumulative = np.cumsum(shuffled, axis=1, dtype=np.uint16)

        pan_curves[i, :] = np.sum(cumulative > 0, axis=0)
        core_curves[i, :] = np.sum(cumulative == genome_counts, axis=0)

        if (i + 1) % 10 == 0 or i == 0 or i + 1 == permutations:
            print(f"  Completed permutation {i + 1}/{permutations}")

    result = pd.DataFrame({
        "genomes_added": np.arange(1, n_genomes + 1),
        "pangenome_mean": pan_curves.mean(axis=0),
        "pangenome_sd": pan_curves.std(axis=0),
        "pangenome_q025": np.quantile(pan_curves, 0.025, axis=0),
        "pangenome_q975": np.quantile(pan_curves, 0.975, axis=0),
        "core_mean": core_curves.mean(axis=0),
        "core_sd": core_curves.std(axis=0),
        "core_q025": np.quantile(core_curves, 0.025, axis=0),
        "core_q975": np.quantile(core_curves, 0.975, axis=0),
    })

    return result


def plot_single_curve(result, label, outdir, prefix):
    label_safe = safe_label(label)

    x = result["genomes_added"]

    # Pangenome curve
    plt.figure(figsize=(8, 6))
    plt.plot(x, result["pangenome_mean"], linewidth=2)
    plt.fill_between(
        x,
        result["pangenome_q025"],
        result["pangenome_q975"],
        alpha=0.2
    )
    plt.xlabel("Number of genomes")
    plt.ylabel("Number of gene clusters")
    plt.title(f"Pangenome accumulation curve: {label}")
    plt.tight_layout()

    for ext in ["png", "pdf", "svg"]:
        plt.savefig(
            os.path.join(outdir, f"{prefix}_{label_safe}_pangenome_accumulation.{ext}"),
            dpi=600
        )
    plt.close()

    # Core genome curve
    plt.figure(figsize=(8, 6))
    plt.plot(x, result["core_mean"], linewidth=2)
    plt.fill_between(
        x,
        result["core_q025"],
        result["core_q975"],
        alpha=0.2
    )
    plt.xlabel("Number of genomes")
    plt.ylabel("Number of core gene clusters")
    plt.title(f"Core genome decay curve: {label}")
    plt.tight_layout()

    for ext in ["png", "pdf", "svg"]:
        plt.savefig(
            os.path.join(outdir, f"{prefix}_{label_safe}_core_genome_decay.{ext}"),
            dpi=600
        )
    plt.close()

    # Combined pan + core
    plt.figure(figsize=(8, 6))
    plt.plot(x, result["pangenome_mean"], linewidth=2, label="Pangenome")
    plt.plot(x, result["core_mean"], linewidth=2, label="Core genome")
    plt.xlabel("Number of genomes")
    plt.ylabel("Number of gene clusters")
    plt.title(f"Pangenome and core-genome curves: {label}")
    plt.legend()
    plt.tight_layout()

    for ext in ["png", "pdf", "svg"]:
        plt.savefig(
            os.path.join(outdir, f"{prefix}_{label_safe}_pan_core_combined.{ext}"),
            dpi=600
        )
    plt.close()


def plot_comparison(all_results, outdir, prefix):
    # Pangenome comparison
    plt.figure(figsize=(9, 6))
    for label, result in all_results.items():
        plt.plot(
            result["genomes_added"],
            result["pangenome_mean"],
            linewidth=2,
            label=label
        )
    plt.xlabel("Number of genomes")
    plt.ylabel("Number of gene clusters")
    plt.title("Pangenome accumulation curve comparison")
    plt.legend()
    plt.tight_layout()

    for ext in ["png", "pdf", "svg"]:
        plt.savefig(
            os.path.join(outdir, f"{prefix}_pangenome_accumulation_comparison.{ext}"),
            dpi=600
        )
    plt.close()

    # Core comparison
    plt.figure(figsize=(9, 6))
    for label, result in all_results.items():
        plt.plot(
            result["genomes_added"],
            result["core_mean"],
            linewidth=2,
            label=label
        )
    plt.xlabel("Number of genomes")
    plt.ylabel("Number of core gene clusters")
    plt.title("Core genome decay curve comparison")
    plt.legend()
    plt.tight_layout()

    for ext in ["png", "pdf", "svg"]:
        plt.savefig(
            os.path.join(outdir, f"{prefix}_core_genome_decay_comparison.{ext}"),
            dpi=600
        )
    plt.close()


def parse_rtab_argument(arg):
    if "=" not in arg:
        raise ValueError(
            "Each --rtab argument must be in the format Label=/path/to/gene_presence_absence.Rtab"
        )
    label, path = arg.split("=", 1)
    return label.strip(), path.strip()


def main():
    parser = argparse.ArgumentParser(
        description="Generate pangenome accumulation and core genome decay curves from Panaroo Rtab files."
    )

    parser.add_argument(
        "--rtab",
        action="append",
        required=True,
        help="Input in the format Label=/path/to/gene_presence_absence.Rtab. Can be used multiple times."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--prefix",
        default="panaroo",
        help="Output filename prefix."
    )

    parser.add_argument(
        "--permutations",
        type=int,
        default=100,
        help="Number of random genome-order permutations. Use 100 for draft, 500+ for final publication."
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed."
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    all_results = {}

    for rtab_arg in args.rtab:
        label, path = parse_rtab_argument(rtab_arg)

        print("=" * 70)
        print(f"Processing: {label}")
        print(f"Input Rtab: {path}")

        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

        gene_col, genome_cols, matrix = read_rtab(path)

        print(f"Gene clusters: {matrix.shape[0]}")
        print(f"Genomes: {matrix.shape[1]}")
        print(f"Permutations: {args.permutations}")

        result = compute_curves(
            matrix,
            permutations=args.permutations,
            seed=args.seed
        )

        label_safe = safe_label(label)
        table_path = os.path.join(
            args.outdir,
            f"{args.prefix}_{label_safe}_pan_core_accumulation_table.tsv"
        )

        result.to_csv(table_path, sep="\t", index=False)

        print(f"Saved table: {table_path}")

        plot_single_curve(result, label, args.outdir, args.prefix)

        all_results[label] = result

    if len(all_results) > 1:
        plot_comparison(all_results, args.outdir, args.prefix)

    print("=" * 70)
    print("Done.")
    print(f"Results saved in: {args.outdir}")


if __name__ == "__main__":
    main()
