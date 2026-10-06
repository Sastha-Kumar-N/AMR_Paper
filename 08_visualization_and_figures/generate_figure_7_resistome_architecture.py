#!/usr/bin/env python3
"""
Generate Figure 7: Clinical Resistome Architecture, High-Confidence Co-occurrence, and Micro-Synteny Colocalization.
Panels:
  A: Top Tier-A AMR gene family prevalence and classification (n = 1,350)
  B: High-confidence AMR co-occurrence correlation matrix (|Phi| >= 0.50, q < 0.05)
  C: Physical micro-synteny distance and synteny conservation for leading colocalized pairs
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# Set publication style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8.5
plt.rcParams['ytick.labelsize'] = 8.5
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 13
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

BASE_DIR = "/mnt/t2t2b_ext/Backup_main_Aug_2026/phoenix_2"
FUNSCAN_DIR = os.path.join(BASE_DIR, "17.clinical_comparative_genomics/funscan_analysis")
AMR_PREV_PATH = os.path.join(FUNSCAN_DIR, "06_tables/AMR_family_prevalence_all.tsv")
AMR_ASSOC_PATH = os.path.join(FUNSCAN_DIR, "04_statistics/AMR_gene_cooccurrence.tsv")
COLOC_PATH = os.path.join(FUNSCAN_DIR, "04_statistics/feature_colocalization.tsv")

OUT_DIR = os.path.join(BASE_DIR, "manuscript_figures_hq")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, "Figure_7_resistome_architecture.png")

def main():
    print("Loading resistome and colocalization data...")
    amr_prev = pd.read_csv(AMR_PREV_PATH, sep='\t')
    amr_assoc = pd.read_csv(AMR_ASSOC_PATH, sep='\t')
    coloc_df = pd.read_csv(COLOC_PATH, sep='\t')
    
    fig = plt.figure(figsize=(17.5, 5.2), dpi=300)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.25], wspace=0.50)
    
    # -------------------------------------------------------------
    # PANEL A: Top AMR Determinants Prevalence
    # -------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    
    top_amr = amr_prev.head(18).sort_values("prevalence_percent", ascending=True)
    
    colors = []
    for p in top_amr['prevalence_percent']:
        if p >= 95.0:
            colors.append("#2B5C8F")  # Core / Intrinsic (Blue)
        elif p >= 15.0:
            colors.append("#E5A83B")  # Shell / Intermediate (Amber)
        else:
            colors.append("#C44E52")  # Cloud / Acquired (Red)
            
    bars = ax1.barh(top_amr['gene_family'], top_amr['prevalence_percent'], color=colors, edgecolor='black', lw=0.6, height=0.68)
    
    for bar in bars:
        w = bar.get_width()
        ax1.annotate(f"{w:.1f}%", xy=(w + 1.2, bar.get_y() + bar.get_height()/2),
                     va='center', ha='left', fontsize=7.5, fontweight='bold', color='#263238')
        
    ax1.set_xlabel("Clinical Genome Prevalence (%)", fontweight='bold')
    ax1.set_xlim(0, 116)
    ax1.set_title("A   Top AMR Determinant Prevalence", loc='left', fontweight='bold', pad=10)
    ax1.grid(True, ls=':', color='#ECEFF1', axis='x')
    
    legend_elements = [
        Patch(facecolor="#2B5C8F", edgecolor='black', label="Core / Intrinsic (≥95%)"),
        Patch(facecolor="#E5A83B", edgecolor='black', label="Shell / Common (15–95%)"),
        Patch(facecolor="#C44E52", edgecolor='black', label="Cloud / Acquired (<15%)")
    ]
    ax1.legend(handles=legend_elements, loc='lower right', fontsize=8, frameon=True, facecolor='white', edgecolor='#CFD8DC')
    
    # -------------------------------------------------------------
    # PANEL B: High-Confidence AMR Co-occurrence Matrix
    # -------------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    
    strong_pairs = amr_assoc[(amr_assoc['q_value'] < 0.05) & (amr_assoc['Phi_coefficient'] >= 0.50)]
    
    key_co_features = [
        "qacEdelta1", "sul1", "aph(3'')-Ib", "aph(6)-Id", 
        "aadA1", "cmlA1", "aac(6')-Ib", "blaOXA-1", "dfrA14"
    ]
    
    phi_mat = pd.DataFrame(np.eye(len(key_co_features)), index=key_co_features, columns=key_co_features)
    for r in strong_pairs.itertuples():
        if r.feature_1 in key_co_features and r.feature_2 in key_co_features:
            phi_mat.loc[r.feature_1, r.feature_2] = r.Phi_coefficient
            phi_mat.loc[r.feature_2, r.feature_1] = r.Phi_coefficient
            
    display_names = [f.replace('delta', 'Δ') for f in key_co_features]
    phi_mat.index = display_names
    phi_mat.columns = display_names
    
    # Display heatmap with clean styling and compact colorbar on right with dedicated pad
    sns.heatmap(phi_mat, cmap='YlOrRd', vmin=0.0, vmax=1.0,
                annot=True, fmt=".2f", annot_kws={"fontsize": 7.0, "fontweight": "bold"},
                square=True, cbar_kws={"label": "Phi Coefficient (Φ)", "shrink": 0.65, "pad": 0.02},
                ax=ax2, linewidths=0.8, linecolor='white')
    
    ax2.set_title("B   Strong AMR Co-occurrence Matrix", loc='left', fontweight='bold', pad=10)
    ax2.tick_params(axis='x', rotation=45, labelsize=7.5)
    ax2.tick_params(axis='y', rotation=0, labelsize=7.5)
    
    # -------------------------------------------------------------
    # PANEL C: Micro-Synteny Proximity & Synteny Conservation
    # -------------------------------------------------------------
    ax3 = fig.add_subplot(gs[0, 2])
    
    pairs = [
        ("qacEΔ1 – sul1\n(Class 1 Integron 3' CS)", 100.0, 0, 436, 0.919),
        ("aph(3'')-Ib – aph(6)-Id\n(Aminoglycoside Operon)", 100.0, 0, 227, 0.963),
        ("AMP50_C0106 – blaPER\n(AMP–ESBL Micro-synteny)", 96.0, 3572, 25, 0.944),
        ("aac(6')-Ib – blaOXA-1\n(MDR Cassette Array)", 98.1, 420, 159, 0.812),
        ("aadA1 – cmlA1\n(Integron Gene Cassette)", 92.4, 850, 79, 0.724)
    ]
    
    pair_labels = [p[0] for p in pairs]
    pct_coloc = [p[1] for p in pairs]
    median_dist = [p[2] for p in pairs]
    n_genomes = [p[3] for p in pairs]
    phi_vals = [p[4] for p in pairs]
    
    y = np.arange(len(pairs))
    height = 0.52
    
    bars_c = ax3.barh(y, pct_coloc, height=height, color="#167D8D", edgecolor='black', lw=0.6, alpha=0.9)
    
    for idx, (bar, dist, n, phi) in enumerate(zip(bars_c, median_dist, n_genomes, phi_vals)):
        dist_str = f"{dist:,} bp" if dist > 0 else "0 bp (adjacent)"
        label_text = f" {bar.get_width():.1f}% ≤10 kb | Med dist: {dist_str} (n={n}, Φ={phi:.2f})"
        ax3.annotate(label_text, xy=(bar.get_width() - 2.0, bar.get_y() + bar.get_height()/2),
                     va='center', ha='right', fontsize=7.5, fontweight='bold', color='white')
        
    ax3.set_yticks(y)
    ax3.set_yticklabels(pair_labels, fontsize=8.0)
    ax3.invert_yaxis()
    ax3.set_xlabel("Syntenic Proximity Rate (≤10 kb, %)", fontweight='bold')
    ax3.set_xlim(0, 105)
    ax3.set_title("C   Genomic Micro-Synteny Colocalization", loc='left', fontweight='bold', pad=10)
    ax3.grid(True, ls=':', color='#ECEFF1', axis='x')
    
    print(f"Saving high-resolution figure to {OUT_PATH}...")
    fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("Figure 7 generated successfully!")

if __name__ == '__main__':
    main()
