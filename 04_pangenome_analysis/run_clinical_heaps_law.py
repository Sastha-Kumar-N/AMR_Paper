#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

rtab = Path("/mnt/phoenix_2/14.pangenome_panaroo/results/clinical/gene_presence_absence.Rtab")
outdir = Path("/mnt/phoenix_2/17.clinical_comparative_genomics/heaps_law/results")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(rtab, sep="\t")
gene_col = df.columns[0]
matrix = df.drop(columns=[gene_col])

X = matrix.values.astype(bool)
n_genes, n_genomes = X.shape

print(f"Clinical gene clusters: {n_genes}")
print(f"Clinical genomes: {n_genomes}")

n_permutations = 5000
rng = np.random.default_rng(42)

pan_curves = []
core_curves = []

for i in range(n_permutations):
    print(f"Permutation {i+1}/{n_permutations}", flush=True)

    order = rng.permutation(n_genomes)
    observed = np.zeros(n_genes, dtype=bool)
    core = np.ones(n_genes, dtype=bool)

    pan_counts = []
    core_counts = []

    for idx in order:
        observed |= X[:, idx]
        core &= X[:, idx]
        pan_counts.append(observed.sum())
        core_counts.append(core.sum())

    pan_curves.append(pan_counts)
    core_curves.append(core_counts)

pan_curves = np.array(pan_curves)
core_curves = np.array(core_curves)

n = np.arange(1, n_genomes + 1)

pan_mean = pan_curves.mean(axis=0)
pan_sd = pan_curves.std(axis=0)
core_mean = core_curves.mean(axis=0)
core_sd = core_curves.std(axis=0)

def heaps(n, kappa, gamma):
    return kappa * (n ** gamma)

params, _ = curve_fit(
    heaps,
    n,
    pan_mean,
    p0=[pan_mean[0], 0.2],
    maxfev=10000
)

kappa, gamma = params
pan_fit = heaps(n, kappa, gamma)

ss_res = np.sum((pan_mean - pan_fit) ** 2)
ss_tot = np.sum((pan_mean - np.mean(pan_mean)) ** 2)
r2 = 1 - (ss_res / ss_tot)

summary = pd.DataFrame({
    "metric": [
        "clinical_genomes",
        "clinical_gene_clusters",
        "heaps_kappa",
        "heaps_gamma",
        "heaps_r_squared",
        "n_permutations"
    ],
    "value": [
        n_genomes,
        n_genes,
        kappa,
        gamma,
        r2,
        n_permutations
    ]
})

curves = pd.DataFrame({
    "n_genomes": n,
    "pan_mean": pan_mean,
    "pan_sd": pan_sd,
    "core_mean": core_mean,
    "core_sd": core_sd,
    "heaps_fit": pan_fit
})

summary.to_csv(outdir / "clinical_heaps_law_summary.tsv", sep="\t", index=False)
curves.to_csv(outdir / "clinical_rarefaction_curves.tsv", sep="\t", index=False)

plt.figure(figsize=(8, 6))
plt.plot(n, pan_mean, label="Observed pan-genome")
plt.fill_between(n, pan_mean - pan_sd, pan_mean + pan_sd, alpha=0.2)
plt.plot(n, pan_fit, linestyle="--", label=f"Heaps fit: gamma={gamma:.4f}")
plt.xlabel("Number of clinical genomes sampled")
plt.ylabel("Pan-genome size")
plt.title("Clinical P. aeruginosa Pan-genome Openness")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "clinical_heaps_law_pan_genome.png", dpi=300)
plt.savefig(outdir / "clinical_heaps_law_pan_genome.pdf")
plt.close()

plt.figure(figsize=(8, 6))
plt.plot(n, core_mean, label="Core genome")
plt.fill_between(n, core_mean - core_sd, core_mean + core_sd, alpha=0.2)
plt.xlabel("Number of clinical genomes sampled")
plt.ylabel("Core genome size")
plt.title("Clinical P. aeruginosa Core-genome Decay")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "clinical_core_genome_decay.png", dpi=300)
plt.savefig(outdir / "clinical_core_genome_decay.pdf")
plt.close()

print("\n===== Heaps law summary =====")
print(summary.to_string(index=False))
print(f"\nResults saved to: {outdir}")
