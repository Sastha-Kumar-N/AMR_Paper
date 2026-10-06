#!/usr/bin/env python3
"""
Master script to reproduce and assemble all publication figures for:
'Global Pangenomic Architecture, Functional Diversification, and Genomic Fluidity of Clinical and Environmental Pseudomonas aeruginosa'
"""

import os, sys, subprocess, shutil
from PIL import Image

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
HQ_DIR = os.path.join(BASE_DIR, "manuscript_figures_hq")
SRC_IMG_DIR = "/home/igris/.gemini/antigravity-cli/brain/3a978a46-1db4-4ae1-9794-85febf1027bd/scratch/high_res_images"

os.makedirs(HQ_DIR, exist_ok=True)

def run_step(cmd, desc):
    print(f"--> {desc}...")
    res = subprocess.run([sys.executable, cmd], cwd=SCRIPT_DIR)
    if res.returncode != 0:
        raise RuntimeError(f"Step failed: {desc}")

def assemble_figures():
    print("--> Composing Figure 1 (Assembly Quality & Completeness)...")
    im1 = Image.open(os.path.join(SRC_IMG_DIR, 'image1.png')).convert('RGBA')
    im3 = Image.open(os.path.join(SRC_IMG_DIR, 'image3.png')).convert('RGBA')
    h1 = im1.height
    w3 = int(im3.width * (h1 / im3.height))
    fig1 = Image.new('RGBA', (im1.width + w3 + 40, h1), (255, 255, 255, 255))
    fig1.paste(im1, (0, 0))
    fig1.paste(im3.resize((w3, h1), Image.Resampling.LANCZOS), (im1.width + 40, 0))
    fig1.convert('RGB').save(os.path.join(HQ_DIR, 'Figure_1_assembly_quality.png'), dpi=(300, 300))

    print("--> Copying Figure 2 (Rare-cluster occupancy)...")
    shutil.copyfile(os.path.join(SRC_IMG_DIR, 'image7.png'), os.path.join(HQ_DIR, 'Figure_2_rare_clusters.png'))

    print("--> Composing Figure 3 (Accumulation curves and Heaps' law)...")
    im8 = Image.open(os.path.join(SRC_IMG_DIR, 'image8.png')).convert('RGBA')
    im9 = Image.open(os.path.join(SRC_IMG_DIR, 'image9.png')).convert('RGBA')
    im10 = Image.open(os.path.join(SRC_IMG_DIR, 'image10.png')).convert('RGBA')
    h3 = min(im8.height, im9.height, im10.height)
    w8 = int(im8.width * (h3 / im8.height))
    w9 = int(im9.width * (h3 / im9.height))
    w10 = int(im10.width * (h3 / im10.height))
    fig3 = Image.new('RGBA', (w8 + w9 + w10 + 60, h3), (255, 255, 255, 255))
    fig3.paste(im8.resize((w8, h3), Image.Resampling.LANCZOS), (0, 0))
    fig3.paste(im9.resize((w9, h3), Image.Resampling.LANCZOS), (w8 + 30, 0))
    fig3.paste(im10.resize((w10, h3), Image.Resampling.LANCZOS), (w8 + w9 + 60, 0))
    fig3.convert('RGB').save(os.path.join(HQ_DIR, 'Figure_3_accumulation_heaps.png'), dpi=(300, 300))

    print("--> Copying Figure 4 (Frequency classes across collections)...")
    shutil.copyfile(os.path.join(SRC_IMG_DIR, 'image11.png'), os.path.join(HQ_DIR, 'Figure_4_frequency_classes.png'))

    # Run scripts for Figures 5, 6, 7
    run_step(os.path.join(SCRIPT_DIR, "generate_figure_5_clinical_vs_environmental.py"), "Generating Figure 5 (Clinical vs Environmental)")
    run_step(os.path.join(SCRIPT_DIR, "generate_figure_6_functional_genomics.py"), "Generating Figure 6 (Multi-class Functional Genomics)")
    run_step(os.path.join(SCRIPT_DIR, "generate_figure_7_resistome_architecture.py"), "Generating Figure 7 (Resistome & Micro-synteny)")

    print("--> Copying Figure 8 (Genomic fluidity)...")
    shutil.copyfile(os.path.join(SRC_IMG_DIR, 'image18.png'), os.path.join(HQ_DIR, 'Figure_8_genomic_fluidity.png'))

    print("--> Copying Figure 9 (Clinical core-SNP phylogeny)...")
    shutil.copyfile(os.path.join(SRC_IMG_DIR, 'image16.png'), os.path.join(HQ_DIR, 'Figure_9_core_phylogeny.png'))

    print("--> Copying Figure 10 (Geographic divergence distribution)...")
    shutil.copyfile(os.path.join(SRC_IMG_DIR, 'image14.png'), os.path.join(HQ_DIR, 'Figure_10_geographic_divergence.png'))

    print("--> Copying Supplementary Figures S1 to S9...")
    supp_map = {
        'Supplementary_Figure_S1_gc_content.png': 'image2.png',
        'Supplementary_Figure_S2_gene_counts.png': 'image4.png',
        'Supplementary_Figure_S3_per_genome_distribution.png': 'image6.png',
        'Supplementary_Figure_S4_pangenome_expansion.png': 'image5.png',
        'Supplementary_Figure_S5_total_pangenome_size.png': 'image12.png',
        'Supplementary_Figure_S6_shared_framework.png': 'image13.png',
        'Supplementary_Figure_S7_fluidity_convergence.png': 'image17.png',
        'Supplementary_Figure_S8_country_sampling.png': 'image15.png',
        'Supplementary_Figure_S9_fluidity_formulation.png': 'image19.png',
    }
    for dst_name, src_name in supp_map.items():
        shutil.copyfile(os.path.join(SRC_IMG_DIR, src_name), os.path.join(HQ_DIR, dst_name))

    print("\nAll figures assembled and verified successfully!")

if __name__ == '__main__':
    assemble_figures()
