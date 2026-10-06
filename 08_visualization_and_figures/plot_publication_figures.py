#!/usr/bin/env python3
"""Render audited historical-cohort summaries; never modify project inputs.

All figures are descriptive. Sources and counting units are recorded per panel.
The retained historical source labels, caller definitions and cluster sets are
preserved. This script does not rerun Panaroo, predict functions or infer lineages.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import textwrap
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pa_publication_matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PUB = ROOT / "publication_manuscript"
TAB = PUB / "reconstructed_tables"
OUT = PUB / "reconstructed_figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5,
    "axes.titlesize": 12, "axes.labelsize": 10.5,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .8, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
CORE = json.loads((TAB / "verified_core_results.json").read_text())
FUNC = json.loads((TAB / "verified_functional_results.json").read_text())
EXACT = json.loads((TAB / "additional_validated_results.json").read_text())
ASSOC = json.loads((TAB / "current_association_integrity.json").read_text())
CAT_COLOURS = ["#1976A8", "#76AFC7", "#D8A04B", "#B4BDC7"]
GLOBAL, CLINICAL, ORANGE, GREY = "#1976A8", "#985592", "#D98235", "#677483"
manifest, legends, metadata = [], {}, {}
cache = {}


def source(figure, panel, path, fields, unit, n, status="RECONSTRUCTED_FROM_VERIFIED_OUTPUT"):
    p = Path(path)
    if not p.is_absolute():
        p = ROOT / p
    p = p.resolve()
    if p not in cache:
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for block in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(block)
        cache[p] = h.hexdigest()
    manifest.append({"figure": figure, "panel": panel,
                     "source_file": str(p.relative_to(ROOT)), "sha256": cache[p],
                     "plot_fields": fields, "unit": unit, "n": n,
                     "evidence_status": status,
                     "script": str(Path(__file__).relative_to(ROOT))})


def tsource(fig, panel, filename, fields, unit, n):
    source(fig, panel, TAB / filename, fields, unit, n)


def panel(ax, letter, title):
    width = max(24, int(ax.get_position().width * ax.figure.get_figwidth() * 7.6))
    ax.set_title(f"{letter}  " + textwrap.fill(title, width=width, break_on_hyphens=False),
                 loc="left", fontweight="bold", pad=12)


def tidy(ax, xgrid=False):
    ax.set_axisbelow(True)
    ax.grid(axis="x" if xgrid else "y", color="#e6e9ed", linewidth=.6)


def save(fig, number, suffix):
    stem = f"FIGURE_{number}_{suffix}"
    for extension in ("pdf", "png"):
        fig.savefig(OUT / f"{stem}.{extension}", dpi=220, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)
    metadata[str(number)] = {"stem": stem,
        "files": [f"{stem}.pdf", f"{stem}.png"],
        "canvas_inches": [float(v) for v in fig.get_size_inches()],
        "status": "AUDITED_HISTORICAL_COHORT_SUMMARY_NOT_SUBMISSION_FINAL"}


def box(ax, x, y, w, h, text, colour="#eef4f7", size=10.5):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle="round,pad=0.025,rounding_size=0.045",
        facecolor=colour, edgecolor="#607789", linewidth=.9))
    ax.text(x, y, text, ha="center", va="center", fontsize=size)


def arrow(ax, a, b):
    ax.annotate("", xy=b, xytext=a,
                arrowprops={"arrowstyle": "->", "color": "#607789", "lw": 1.1})


def figure_1():
    d, q = CORE["deduplication"], CORE["cohort"]
    fig, axes = plt.subplots(2, 1, figsize=(13.5, 7.2),
                             gridspec_kw={"height_ratios": [1.1, 1.9]})
    a, b = axes
    a.set(xlim=(0, 6), ylim=(0, 2)); a.axis("off")
    panel(a, "a", "Retained file manifests and revised curation branch")
    labels = [f"Input files\n{d['stage1_rows']:,}",
              f"Exact sequence\nhashes: {d['stage1_unique_hashes']:,}",
              f"Contig multiset\nhashes: {d['stage2_unique_hashes']:,}",
              f"Revised QUAST pass\n{q['quast_status']['PASS']:,}",
              f"CheckM2 pass\n{q['checkm_standard_count']:,}",
              f"P. aeruginosa\n{q['meta_n']:,}"]
    exclusions = ["2,113 exact duplicates", "22 exact contig-multiset\nduplicates",
                  "7 fail; 4 review excluded", "9 fail; 6 review excluded",
                  "9 P. paraeruginosa excluded"]
    xs = np.arange(6) + .5
    for i, (x, text) in enumerate(zip(xs, labels)):
        box(a, x, 1.25, .83, .72, text, "#dbeaf2" if i==5 else "#f3f5f7", 10)
        if i < 5:
            arrow(a, (x+.425, 1.25), (x+.575, 1.25))
            a.text(x+.5, .58, exclusions[i], ha="center", va="center", fontsize=8.6,
                   rotation=0, wrap=True)
    a.text(.5, .12, "Counts describe retained files; exact sequence deduplication does not establish isolate independence.",
           transform=a.transAxes, ha="center", va="center", fontsize=10, color=GREY)

    b.set(xlim=(0, 10.8), ylim=(-.5, 4.1)); b.axis("off")
    panel(b, "b", "Historical analysis cohorts and paired-assembly sensitivity")
    box(b, 1.3, 1.65, 2.15, 1.0, "Final collection\n1,595 assemblies", "#dbeaf2", 11)
    yv = [3.35, 2.2, 1.05, -.10]
    groups = ["Global\n1,595", "Clinical-labelled\n1,350",
              "Environmental-labelled\n99", "Unknown / other\n146"]
    analyses = ["Original global Panaroo matrix\n43,001 clusters",
                "Clinical Panaroo; AMR; AMP; CAZyme; core SNPs\nHistorical n = 1,350",
                "Separate environmental Panaroo matrix\nHistorical n = 99",
                "Retained in the global matrix\nNo separate functional cohort inferred"]
    for y, grp, ana in zip(yv, groups, analyses):
        box(b, 4.0, y, 2.0, .74, grp, "#f3f5f7", 10)
        arrow(b, (2.4, 1.65), (2.97, y))
        arrow(b, (5.03, y), (5.6, y))
        b.text(5.72, y, ana, va="center", fontsize=10.3)
    fig.text(.5, .007, "Three paired-assembly exclusions: fixed-cluster sensitivity n = 1,592 global / 1,347 clinical; later geographic subset n = 1,347, with 601 coordinate records.",
             ha="center", fontsize=9.5, color=GREY)
    fig.subplots_adjust(hspace=.44, bottom=.06, top=.94)
    tsource(1, "a,b", "verified_core_results.json", "deduplication; cohort; meta_groups", "assembly files", 1595)
    source(1, "a,b", PUB / "04_COHORT_FLOW.tsv", "all stages and branches", "assembly files", 1595)
    tsource(1, "b", "verified_phylogeny_geography.json", "paired_assembly_proof; geography", "assembly files / coordinate records", 1347)
    legends[1] = "Retained file manifests and analytical cohorts. (a) Two exact-sequence deduplication stages, revised QUAST selection, CheckM2 filtering and GTDB-Tk assignment retained 1,595 P. aeruginosa assembly files. The second hash compares a contig multiset, retaining repeated-contig multiplicity. (b) Historical source labels define 1,350 clinical-labelled, 99 environmental-labelled and 146 unknown/other assemblies. Clinical-labelled functional and core-SNP outputs retain historical n = 1,350. The footer reports fixed-cluster sensitivities excluding three paired representations, one from each of three BioSample pairs, and the later geographic branch. These are not primary Panaroo reruns. Exact file deduplication does not establish biological independence."
    save(fig, 1, "curation_and_cohorts")


def figure_2():
    w = pd.read_csv(TAB / "occupancy_weighted_architecture.tsv", sep="\t")
    w = w[w.population.eq("whole")].reset_index(drop=True)
    occ = pd.read_csv(TAB / "verified_occupancy_whole.tsv", sep="\t")
    per = pd.read_csv(TAB / "per_genome_frequency_classes.tsv", sep="\t")
    per = per[per.population.eq("whole")]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.0), gridspec_kw={"width_ratios":[1,1.1,1.1]})
    a, b, c = axes
    panel(a, "a", "Two counting scales")
    values = np.array([w.distinct_percent.to_numpy(), w.assignment_percent.to_numpy()]).T
    lower = np.zeros(2)
    for i, (vals, colour) in enumerate(zip(values, CAT_COLOURS)):
        a.bar([0,1], vals, bottom=lower, width=.62, color=colour,
              label=["Core ≥99%", "Soft core ≥95–<99%", "Shell ≥15–<95%", "Cloud >0–<15%"][i])
        for j in range(2):
            if vals[j] >= 7:
                a.text(j, lower[j]+vals[j]/2, f"{vals[j]:.2f}%", ha="center", va="center",
                       fontsize=10, color="white" if i==0 else "#142330", fontweight="bold")
        lower += vals
    a.set(ylim=(0,100), ylabel="Percentage of stated denominator",
          xticks=[0,1], xticklabels=["Distinct clusters\n43,001", "Pooled presences\n9,700,435"])
    a.legend(loc="upper left", bbox_to_anchor=(-.03,-.20), ncol=2, fontsize=9)
    tidy(a)
    panel(b, "b", "Global occupancy spectrum")
    spectrum = occ.prevalence.value_counts().sort_index()
    b.scatter(spectrum.index, spectrum.values, color=GLOBAL, s=14, alpha=.8)
    b.set(xscale="log", yscale="log", xlabel="Assemblies containing a cluster (log)",
          ylabel="Distinct clusters at that occupancy (log)")
    b.set_xticks([1,2,5,10,100,1595], labels=["1","2","5","10","100","1,595"])
    b.text(.05,.94, "12,293 singletons\n4,400 doubletons", transform=b.transAxes,
           va="top", fontsize=10, bbox={"facecolor":"white","alpha":.9,"edgecolor":"none"})
    tidy(b)
    panel(c, "c", "Class shares within assemblies")
    percent_cols = ["core99_percent", "soft95_percent", "shell15_percent", "cloud_percent"]
    if not set(percent_cols).issubset(per.columns):
        # The existing table preserves category labels; inspect their explicit names.
        category_cols = ["Core ≥99%", "Soft core ≥95–<99%", "Shell ≥15–<95%", "Cloud <15%"]
        total = per["cluster_count"].to_numpy()
        data = [100*per[k].to_numpy()/total for k in category_cols]
    else:
        data = [per[k].to_numpy() for k in percent_cols]
    bp = c.boxplot(data, tick_labels=["Core", "Soft core", "Shell", "Cloud"],
                   widths=.55, showfliers=False, patch_artist=True,
                   medianprops={"color":"#12202b","linewidth":1.4})
    for patch, colour in zip(bp["boxes"], CAT_COLOURS): patch.set_facecolor(colour)
    c.set(ylim=(0,100), ylabel="Percentage of per-assembly cluster set")
    c.text(.98,.93,"n = 1,595 assemblies", transform=c.transAxes, ha="right", fontsize=9.5)
    tidy(c)
    fig.subplots_adjust(wspace=.38, bottom=.24, top=.87)
    for pn, fn, fs, unit, n in [
        ("a","occupancy_weighted_architecture.tsv","distinct_percent; assignment_percent","distinct clusters / genome-cluster presences",1595),
        ("b","verified_occupancy_whole.tsv","prevalence","distinct clusters by assembly occupancy",43001),
        ("c","per_genome_frequency_classes.tsv","global category counts and total","assembly-level cluster shares",1595)]: tsource(2,pn,fn,fs,unit,n)
    legends[2] = "Global pangenome architecture at collective-richness and genome-content scales. The saved original binary matrix contains 43,001 clusters across 1,595 assemblies. (a) Cloud clusters represent 85.10% of distinct clusters but 7.26% of 9,700,435 pooled genome–cluster presences; ≥99% core clusters contribute 74.74% of those presences. (b) Distinct cluster counts by assembly occupancy are shown on logarithmic axes. (c) Per-assembly class shares use the same fixed global classes. Boxes show medians and IQRs, whiskers extend to the most extreme observations within 1.5 IQR of the quartiles, and outlier symbols are omitted for clarity. Presence weighting is not DNA length, expression or copy-number weighting."
    save(fig, 2, "cluster_richness_and_presence")


def figure_3():
    fig, axes = plt.subplots(1,3,figsize=(13.5,4.9),gridspec_kw={"width_ratios":[1.4,1.25,1]})
    a,b,c=axes
    panel(a,"a","Conditional pangenome accumulation")
    panel(b,"b","Expected marginal cluster discovery")
    for pop, colour, label in [("whole",GLOBAL,"Global"),("clinical",CLINICAL,"Clinical-labelled")]:
        t=pd.read_csv(TAB / f"{pop}_exact_rarefaction.tsv",sep="\t")
        a.plot(t.n_genomes,t.expected_pangenome,color=colour,lw=2.0,label=f"{label}: exact expectation")
        a.plot(t.n_genomes,t.fitted_pangenome,color=colour,lw=1.1,ls="--",label=f"{label}: cumulative fit")
        z=t.iloc[1:]
        b.plot(z.n_genomes,z.expected_new_clusters,color=colour,lw=1.8,label=label)
        b.scatter([z.n_genomes.iloc[-1]],[z.expected_new_clusters.iloc[-1]],color=colour,s=24)
        tsource(3,"a,b",f"{pop}_exact_rarefaction.tsv","n_genomes; expected_pangenome; fitted_pangenome; expected_new_clusters","conditional expected clusters",len(t))
    a.set(xlabel="Sampled assemblies",ylabel="Expected distinct clusters",ylim=(0,45000))
    a.yaxis.set_major_formatter(FuncFormatter(lambda y,_:f"{y/1000:.0f},000" if y else "0"))
    a.legend(fontsize=8.7,loc="lower right")
    a.text(.04,.95,"P(n) = κ nᵞ\nγ = 0.290002 global\nγ = 0.288462 clinical-labelled",
           transform=a.transAxes,va="top",fontsize=9.5)
    b.set(xscale="log",yscale="log",xlabel="Sampled assemblies (log)",
          ylabel="Expected new clusters (log)")
    b.legend(loc="upper right",fontsize=9.2)
    b.text(.04,.05,"Final increment: 7.71 global; 8.67 clinical\nConditional on the retained collection",
           transform=b.transAxes,fontsize=9.0,color=GREY)
    panel(c,"c","Paired-assembly sensitivity")
    y=np.array([1,0])
    for yy,pop,colour in zip(y,["whole","clinical"],[GLOBAL,CLINICAL]):
        full=EXACT[pop]["exact_fit"]["gamma"]
        sens=EXACT[pop]["paired_exclusion_sensitivity"]["exact_fit"]["gamma"]
        c.plot([full,sens],[yy,yy],color=colour,lw=2)
        c.scatter([full],[yy],s=52,color=colour,label="Historical cohort" if yy==1 else None)
        c.scatter([sens],[yy],s=62,facecolor="white",edgecolor=colour,lw=1.6,
                  label="Three representations\nexcluded (one per pair)" if yy==1 else None)
        c.text(full-.00012,yy+.15,f"{full:.6f}",ha="right",fontsize=9,color=colour)
        c.text(sens+.00012,yy-.22,f"{sens:.6f}",ha="left",fontsize=9,color=colour)
    c.set(yticks=y,yticklabels=["Global","Clinical-labelled"],ylim=(-.65,1.95),
          xlim=(.2875,.291),xlabel="Cumulative-fit exponent γ")
    c.legend(loc="upper center",fontsize=8.7)
    c.text(.5,-.22,"Fixed clusters; no inferential CI",transform=c.transAxes,
           ha="center",fontsize=9.2,color=GREY)
    for ax in axes:tidy(ax)
    fig.subplots_adjust(wspace=.44,bottom=.19,top=.86)
    tsource(3,"c","additional_validated_results.json","whole/clinical exact_fit and paired_exclusion_sensitivity","conditional cumulative-fit exponent",1595)
    legends[3] = "Exact finite-collection rarefaction and paired-representation sensitivity. (a) Curves integrate all equally weighted sampling orders conditional on the retained global (n = 1,595) and clinical-labelled (n = 1,350) occupancy spectra. Dashed lines fit P(n) = κn^γ across all sampled sizes; positive fitted exponents describe continued accumulation within that range. (b) Marginal discovery excludes the initial n = 1 genome from the log plot and ends at 7.71 and 8.67 expected clusters. (c) Points compare fitted cumulative exponents before and after excluding three paired representations, one from each BioSample pair; connecting lines are visual comparisons, not curves or uncertainty intervals. The sensitivity retains fixed original cluster definitions (n = 1,592 / 1,347), not a Panaroo rerun. No population confidence interval or extrapolation beyond the collection is implied."
    save(fig,3,"exact_accumulation_and_sensitivity")


def figure_4():
    h=pd.read_csv(TAB/"clinical_fluidity_histogram.tsv",sep="\t")
    s=pd.read_csv(TAB/"clinical_fluidity_sensitivity.tsv",sep="\t")
    raw=CORE["fluidity_recomputed"][0]
    fsummary=pd.read_csv(ROOT/"17.clinical_comparative_genomics/pangenome fluidity/fluidity_results/fluidity_summary.tsv",sep="\t").iloc[0]
    fig,axes=plt.subplots(1,3,figsize=(13.5,5.3),gridspec_kw={"width_ratios":[1.35,1.1,1]})
    a,b,c=axes
    panel(a,"a","Dependent pairwise dissimilarities")
    a.bar(h.lower,h.pairs,width=h.upper-h.lower,align="edge",color=GLOBAL,alpha=.9,edgecolor="none")
    a.axvline(raw["mean"],color="#202a33",ls="--",lw=1.5,label=f"Mean {raw['mean']:.5f}")
    a.axvline(raw["median"],color=ORANGE,ls=":",lw=1.7,label=f"Median {raw['median']:.5f}")
    a.set(xlabel="Pairwise genomic fluidity φ",ylabel="Unordered genome pairs",xlim=(0,.25))
    a.legend(fontsize=9,loc="upper left")
    a.text(.98,.91,"1,350 assemblies\n910,575 pairs\n(dependent)",transform=a.transAxes,
           ha="right",va="top",fontsize=9.2)
    tidy(a)
    panel(b,"b","Nominal assembly-resampling intervals")
    intervals=[(fsummary.jackknife_ci95_low,fsummary.jackknife_ci95_high),
               (fsummary.bootstrap_ci95_low,fsummary.bootstrap_ci95_high)]
    centers=[fsummary.mean_phi,fsummary.bootstrap_mean]
    for yy,(lo,hi),ctr,col in zip([1,0],intervals,centers,[GLOBAL,CLINICAL]):
        b.errorbar(ctr,yy,xerr=[[ctr-lo],[hi-ctr]],fmt="o",color=col,capsize=5,lw=1.8)
    b.set(yticks=[1,0],yticklabels=["Assembly jackknife", "Saved assembly\nbootstrap"],
          xlabel="Mean fluidity (nominal 95% interval)",xlim=(.1047,.1099),ylim=(-.6,1.65))
    b.text(.5,-.25,"Assembly exchangeability assumed\nSaved bootstrap: 1,000 resamples\nNot independently regenerated",
           transform=b.transAxes,ha="center",va="top",fontsize=8.5,color=GREY)
    tidy(b,True)
    panel(c,"c","Minimum-occupancy sensitivity")
    x=np.arange(len(s))
    c.plot(x,s["mean"],"o-",color=GLOBAL,lw=1.8)
    for xx,row in s.iterrows():
        c.text(xx,row["mean"]+.0035,f"{row['mean']:.5f}",
               ha="center",fontsize=8.8)
    c.set(xticks=x,xticklabels=["≥1\n40,811","≥2\n29,109","≥5\n20,625","≥14 (≥1%)\n15,280"],
          xlim=(-.5,3.5),ylim=(0,.125),xlabel="Minimum occupancy\nRetained cluster count below",ylabel="Mean fluidity")
    c.tick_params(axis="x",labelsize=8.7)
    tidy(c)
    fig.subplots_adjust(wspace=.62,bottom=.26,top=.85)
    for pn,fn,fields,unit,n in [
        ("a","clinical_fluidity_histogram.tsv","lower; upper; pairs","dependent unordered pairs",910575),
        ("c","clinical_fluidity_sensitivity.tsv","minimum_occupancy; retained_clusters; mean","clusters and pairwise mean",1350)]:tsource(4,pn,fn,fields,unit,n)
    source(4,"b",ROOT/"17.clinical_comparative_genomics/pangenome fluidity/fluidity_results/fluidity_summary.tsv",
           "jackknife_ci95; bootstrap_ci95; mean_phi; bootstrap_mean","nominal assembly-resampling estimate/interval under exchangeability",1350,"JACKKNIFE_INDEPENDENTLY_RECONCILED_BOOTSTRAP_SAVED_SUMMARY_NOT_REGENERATED")
    legends[4] = "Clinical-labelled genomic fluidity and sensitivity to rare clusters. (a) The saved 40,811-cluster matrix produces 910,575 dependent unordered pairs from 1,350 assemblies. Mean fluidity is 0.10736 and median 0.10907. (b) Delete-one-assembly jackknife and saved percentile intervals from 1,000 assembly bootstrap replicates are nominal 95% resampling intervals under assembly exchangeability. They do not account for relatedness or study-sampling imbalance. The jackknife was independently reconciled; the bootstrap marker is the saved bootstrap mean, and the bootstrap summary was not independently regenerated because its replicate values were not retained. (c) Minimum-occupancy filtering reduces mean fluidity to 0.09921 at ≥14 assemblies (≥1%). Removing rare clusters changes both numerator and denominator. Fluidity measures cluster-content dissimilarity, not gene transfer or evolutionary rate."
    save(fig,4,"clinical_fluidity_and_rare_sensitivity")


def figure_5():
    burden=pd.read_csv(TAB/"standalone_amr_burden.tsv",sep="\t")
    prev=pd.read_csv(TAB/"verified_amr_element_prevalence.tsv",sep="\t").set_index("element").genomes
    groups=FUNC["standalone_AMR"]["gene_groups"]
    labels=["fosA","aph(3′)-IIb","catB7","blaPDC family","blaOXA family",
            "blaVIM family","blaIMP family","blaKPC family","blaNDM family"]
    counts=[int(prev["fosA"]),int(prev["aph(3')-IIb"]),int(prev["catB7"])]+[
        int(groups[g]["genomes"]) for g in ["blaPDC","blaOXA","blaVIM","blaIMP","blaKPC","blaNDM"]]
    raw_amr_path=ROOT/"18.amrfinderplus/results/clinical_amrfinderplus/combined_amrfinderplus_results.tsv"
    raw_amr=pd.read_csv(raw_amr_path,sep="\t",keep_default_na=False)
    fos=raw_amr.loc[raw_amr["Element symbol"].eq("fosA")]
    fos_internal_stop=int(fos.loc[fos.Method.eq("INTERNAL_STOP"),"Name"].nunique())
    fos_partial=int(fos.loc[fos.Method.eq("PARTIALP"),"Name"].nunique())
    fig,(a,b)=plt.subplots(1,2,figsize=(12.8,5.2),gridspec_kw={"width_ratios":[1,1.3]})
    panel(a,"a","Distinct AMR-associated symbols")
    data=[burden[k].to_numpy() for k in ["unique_symbols","nonpoint_symbols","point_symbols"]]
    bp=a.boxplot(data,tick_labels=["All AMR\nsymbols\nMedian 10","Nonpoint\nsymbols\nMedian 5","Point / disruption\nsymbols\nMedian 4"],
                  showfliers=True,widths=.48,patch_artist=True,
                  flierprops={"marker":".","markersize":2.8,"alpha":.35},
                  medianprops={"color":"#17232e","linewidth":1.5})
    for patch,col in zip(bp["boxes"],[GLOBAL,"#8BB4CA",ORANGE]):patch.set_facecolor(col)
    a.set(ylabel="Distinct exact symbols per assembly",ylim=(0,38))
    a.text(.03,.96,"n = 1,350 assemblies",transform=a.transAxes,va="top",fontsize=10)
    tidy(a)
    panel(b,"b","Selected annotation-detection prevalence")
    yy=np.arange(len(labels))
    colours=[GLOBAL]*5+[ORANGE]*4
    b.barh(yy,np.array(counts)/len(burden)*100,color=colours,height=.63)
    for y,n in zip(yy,counts):
        percentage=n/len(burden)*100
        b.text(percentage-1.8 if percentage>70 else percentage+1.2,y,
               f"{n:,}/1,350 ({percentage:.1f}%)",va="center",fontsize=9.5,
               ha="right" if percentage>70 else "left",color="white" if percentage>70 else "#17232e")
    b.set(yticks=yy,yticklabels=labels,xlim=(0,102),xlabel="Assemblies with ≥1 annotation (%)")
    b.invert_yaxis();tidy(b,True)
    fig.text(.5,.007,"18,058 AMR records; 929 exact symbols = 413 nonpoint + 516 point/disruption. Prefix aggregates are not measured susceptibility.\n"
             f"fosA detections include {fos_internal_stop} INTERNAL_STOP and {fos_partial} PARTIALP calls; detection is not intact-locus prevalence.",
             ha="center",va="bottom",fontsize=8.9,color=GREY)
    fig.subplots_adjust(wspace=.40,bottom=.21,top=.88)
    tsource(5,"a","standalone_amr_burden.tsv","unique_symbols; nonpoint_symbols; point_symbols","distinct AMR annotation symbols per assembly",1350)
    tsource(5,"b","verified_amr_element_prevalence.tsv","element; genomes","positive assemblies per exact symbol",1350)
    tsource(5,"b","verified_functional_results.json","standalone_AMR.gene_groups","positive assemblies per prefix-defined family",1350)
    source(5,"b and qualification",raw_amr_path,"Name; Element symbol; Method for fosA","annotation-detection methods, including interrupted/partial matches",1350,"PRIMARY_CALLER_METHOD_COUNTS_INDEPENDENTLY_RECONCILED")
    legends[5] = "Standalone AMRFinderPlus resistance-associated annotations in 1,350 clinical-labelled assemblies. (a) Boxes show medians/IQRs of per-assembly distinct symbols, whiskers extend to the most extreme observations within 1.5 IQR of the quartiles, and dots show individual outliers. Nonpoint annotations and point/disruption symbols are separated; their medians are not additive. (b) Detection prevalence counts each positive assembly once, including partial/interrupted caller annotations. The 1,350 fosA detections include " + str(fos_internal_stop) + " INTERNAL_STOP and " + str(fos_partial) + " PARTIALP calls; 100% detection does not establish universal intact loci. bla-prefix aggregates can overlap within assemblies and contain multiple alleles with different biological contexts. The catalogue comprises 18,058 AMR-type rows and 929 distinct exact symbols (413 nonpoint, 516 point/disruption). STRESS/biocide rows and mutation_all wild-type/unknown baseline records are excluded from this burden. Genomic annotation does not establish measured antibiotic susceptibility or phenotypic MDR."
    save(fig,5,"standalone_AMR_architecture")


def figure_6():
    burden=pd.read_csv(TAB/"functional_burden.tsv",sep="\t")
    amp=FUNC["AMP"]
    fig,axes=plt.subplots(2,2,figsize=(12.8,7.6))
    a,b,c,d=axes.ravel()
    panel(a,"a","AMP prediction support")
    vals=[amp["tier_counts"]["Tier 3 - putative single-predictor"],
          amp["tier_counts"]["Tier 2 - high confidence"],
          amp["tier_counts"]["Tier 1 - very high confidence"]]
    a.bar(np.arange(3),vals,color=["#B4BDC7",GLOBAL,ORANGE],width=.57)
    for i,n in enumerate(vals):a.text(i,n*1.40,f"{n:,}",ha="center",fontsize=11,fontweight="bold")
    a.set(xticks=[0,1,2],xticklabels=["Single predictor","Multiple predictors","Database-supported"],
          yscale="log",ylim=(.65,110000),ylabel="Prediction occurrences (log)")
    a.text(.98,.94,"1,362 / 40,061 higher-support\n(3.40%; activity unvalidated)",transform=a.transAxes,
           ha="right",va="top",fontsize=10)
    tidy(a)
    panel(b,"b","AMP sequence-clustering sensitivity")
    idens=[50,70,90];families=[amp["families_per_identity"][str(i)] for i in idens]
    b.bar(np.arange(3),families,color=CLINICAL,width=.57)
    for i,n in enumerate(families):b.text(i,n+28,f"{n:,}",ha="center",fontsize=11)
    b.set(xticks=[0,1,2],xticklabels=["50%","70%","90%"],ylim=(0,1700),
          xlabel="Sequence identity threshold (80% coverage)",ylabel="Sequence families")
    b.text(.5,.98,"40,061 occurrence entries clustered\nRepresenting 2,109 unique sequences\nUniversal higher-support family: 1,350 / 1,350",
           transform=b.transAxes,ha="center",va="top",fontsize=9.2)
    tidy(b)
    panel(c,"c","CAZyme family/subfamily annotation classes")
    labels=list(FUNC["CAZyme"]["class_occurrences"])
    vals=list(FUNC["CAZyme"]["class_occurrences"].values())
    yy=np.arange(len(labels))
    c.barh(yy,vals,color=GLOBAL,height=.62)
    for y,n in zip(yy,vals):c.text(n+1200,y,f"{n:,}",va="center",fontsize=9.7)
    c.set(yticks=yy,yticklabels=labels,xlim=(0,80000),xlabel="Family/subfamily annotation tokens")
    c.invert_yaxis();tidy(c,True)
    c.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x/1000:.0f}k"))
    panel(d,"d","CAZyme annotation units")
    data=[burden.CAZyme_domain_annotations.to_numpy(),burden.CAZyme_unique_genes.to_numpy()]
    bp=d.boxplot(data,tick_labels=["Family/subfamily\nannotation tokens","Distinct annotated\ngenes"],
                  patch_artist=True,widths=.42,showfliers=False,
                  medianprops={"color":"#17232e","linewidth":1.4})
    for patch,col in zip(bp["boxes"],[GLOBAL,"#8BB4CA"]):patch.set_facecolor(col)
    for i,v in enumerate(data,1):d.text(i,np.percentile(v,75)+12,f"Median {np.median(v):.0f}",ha="center",fontsize=10)
    d.set(ylim=(0,140),ylabel="Annotation tokens / genes per assembly")
    d.text(.97,.06,"151,816 annotation tokens\n136,278 genes; 73 families\nn = 1,350 assemblies",
           transform=d.transAxes,ha="right",va="bottom",fontsize=9.7)
    tidy(d)
    fig.subplots_adjust(wspace=.30,hspace=.57,bottom=.08,top=.93)
    tsource(6,"a,b","verified_functional_results.json","AMP tier_counts; families_per_identity; unique_sequences","40,061 prediction-occurrence entries representing 2,109 unique sequences / operational families",1350)
    source(6,"b",ROOT/"17.clinical_comparative_genomics/funscan_analysis/06_tables/AMP_prevalence_all.tsv",
           "AMP50_C0001 n_genomes/high_confidence_occurrences","assemblies carrying predicted family",1350,"DERIVED_MATRIX_PREVALENCE_RECONCILED")
    tsource(6,"c","verified_functional_results.json","CAZyme class_occurrences","family/subfamily annotation tokens, not a physical-domain census",151816)
    tsource(6,"d","functional_burden.tsv","CAZyme_domain_annotations; CAZyme_unique_genes","per-assembly family/subfamily annotation tokens / distinct genes",1350)
    legends[6] = "Prediction support and functional counting units in 1,350 clinical-labelled assemblies. (a) Saved mutually exclusive AMP support tiers partition 40,061 prediction occurrences into 38,699 single-predictor, 1,361 multi-predictor and one database-supported occurrence. The 1,362 higher-support occurrences are computational support categories, not a validated true-positive fraction; separately reported predictor-specific counts overlap. (b) Clustering input comprised 40,061 occurrence entries representing 2,109 unique predicted sequences. Clustering at 50%, 70% and 90% identity with 80% coverage yields 1,051, 1,077 and 1,147 operational families; one higher-support family occurs in all 1,350 assemblies. (c) Consensus CAZyme family/subfamily annotation tokens span glycosyltransferases (GT), glycoside hydrolases (GH), carbohydrate-binding modules (CBM), auxiliary activities (AA), carbohydrate esterases (CE) and polysaccharide lyases (PL). (d) There are 151,816 family/subfamily annotation tokens assigned to 136,278 distinct genes. These tokens are not a sequence-coordinate-resolved physical-domain census. Boxes show medians/IQRs, whiskers extend to the most extreme observations within 1.5 IQR of the quartiles, and outlier symbols are omitted for clarity. Predictions and concordance do not establish antimicrobial or enzymatic activity."
    save(fig,6,"predicted_AMP_and_CAZyme_units")


def figure_7():
    burdenpath=ROOT/"17.clinical_comparative_genomics/funscan_analysis/03_matrices/genome_functional_burden_matrix.tsv"
    burden=pd.read_csv(burdenpath,sep="\t")
    gap=pd.read_csv(TAB/"standalone_qac_sul_raw_coordinate_sensitivity.tsv",sep="\t")
    x,y=burden.ARG_total.to_numpy(),burden.AMP_total.to_numpy()
    ranks=burden[["ARG_total","AMP_total","total_CDS"]].rank(method="average")
    design=np.column_stack([np.ones(len(ranks)),ranks.total_CDS.to_numpy()])
    residuals=np.column_stack([ranks.ARG_total.to_numpy(),ranks.AMP_total.to_numpy()])
    residuals-=design @ np.linalg.lstsq(design,residuals,rcond=None)[0]
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.9),gridspec_kw={"width_ratios":[1.05,1.05,1.0]})
    a,b,c=axes
    panel(a,"a","Unadjusted burden co-distribution")
    a.scatter(x,y,s=11,color=GLOBAL,alpha=.24,edgecolors="none")
    a.set(xlabel="Normalized resistance-associated\nannotation occurrence proxy",ylabel="Predicted AMP occurrences")
    a.text(.04,.95,f"Spearman ρ = {ASSOC['AMP_ARG_burden']['unadjusted_spearman']:.4f}\nn = 1,350 assemblies",
           transform=a.transAxes,va="top",fontsize=10.5,
           bbox={"facecolor":"white","alpha":.85,"edgecolor":"none"})
    tidy(a)
    panel(b,"b","After Prokka CDS-count rank adjustment")
    b.scatter(residuals[:,0],residuals[:,1],s=11,color=CLINICAL,alpha=.25,edgecolors="none")
    b.axhline(0,color="#aeb8c0",lw=.7);b.axvline(0,color="#aeb8c0",lw=.7)
    b.set(xlabel="Annotation-proxy rank residual",ylabel="AMP rank residual")
    b.text(.04,.95,f"Pearson r of rank residuals = {ASSOC['AMP_ARG_burden']['CDS_adjusted_rank_residual_correlation']:.5f}\nNo lineage / phylogeny adjustment",
           transform=b.transAxes,va="top",fontsize=9.1,
           bbox={"facecolor":"white","alpha":.85,"edgecolor":"none"})
    tidy(b)
    panel(c,"c","Raw qacEdelta1–sul1 proximity")
    counts=gap.minimum_intergenic_gap_bp.value_counts().sort_index()
    xx=np.arange(len(counts))
    c.bar(xx,counts.values,color=ORANGE,width=.58)
    for i,n in enumerate(counts.values):c.text(i,n*1.2,f"{n:,}",ha="center",fontsize=10.5)
    c.set(xticks=xx,xticklabels=[f"{int(v):,}" for v in counts.index],
          yscale="log",ylim=(.7,1000),xlabel="Minimum same-contig\nintergenic gap (bp)",
          ylabel="Co-positive assemblies (log)")
    c.text(.97,.94,"436 / 436 same-contig\nAt least one pair ≤1,950 bp",
           transform=c.transAxes,ha="right",va="top",fontsize=10,
           bbox={"facecolor":"white","alpha":.9,"edgecolor":"none"})
    tidy(c)
    fig.text(.5,.007,"Burden association, physical proximity and molecular interaction are different evidence types; none of these panels establishes interaction or cotransfer.",
             ha="center",fontsize=9.5,color=GREY)
    fig.subplots_adjust(wspace=.38,bottom=.19,top=.86)
    source(7,"a,b",burdenpath,"ARG_total; AMP_total; total_CDS (Prokka)","broad normalized resistance/stress/regulatory-associated annotation occurrence proxy / AMP predictions / Prokka CDS counts",1350,"HISTORICAL_BURDEN_INPUT_STATISTICS_INDEPENDENTLY_RECONCILED")
    tsource(7,"a,b","current_association_integrity.json","AMP_ARG_burden unadjusted and CDS-adjusted correlations","assembly-level Spearman correlation / Pearson correlation of rank residuals",1350)
    tsource(7,"c","standalone_qac_sul_raw_coordinate_sensitivity.tsv","assembly; same_contig; minimum_intergenic_gap_bp; starts/ends","minimum raw-caller intergenic gap per co-positive assembly",436)
    legends[7] = "Co-distribution, covariate sensitivity and independent raw-coordinate context. (a) Predicted AMP occurrence burden correlates with the saved custom FuncScan-normalized resistance/stress/regulatory-associated annotation occurrence proxy (Spearman ρ = 0.4816; n = 1,350). The proxy is broader than antibiotic-resistance genes, differs from standalone AMRFinderPlus exact symbols, and does not uniformly represent physical loci. (b) Pearson correlation of rank residuals is r = −0.06966 after residualizing each ranked burden on an intercept and ranked Prokka CDS count. Prokka CDS counts differ from the Pyrodigal predictions used by functional modules. This adjustment does not account for lineage, phylogeny or study sampling and does not establish causality or antagonism. (c) Separately, original standalone caller coordinates identify at least one qacEdelta1–sul1 same-contig pair within 1,950 bp in all 436 co-positive assemblies. qacEdelta1 is a STRESS/BIOCIDE annotation; sul1 is an AMR annotation. Minimum gaps are zero in 433 assemblies and 319, 1,756 and 1,950 bp in the remaining three; zero includes overlapping or abutting calls. The bars are gap categories rather than a distance density. This nearest-pair summary does not establish every copy's arrangement, mobile-element boundaries, cotransfer or molecular interaction."
    save(fig,7,"association_sensitivity_and_gene_proximity")


def main():
    for func in [figure_1,figure_2,figure_3,figure_4,figure_5,figure_6,figure_7]:
        func()
        print(f"Rendered {func.__name__}",flush=True)
    with (OUT/"FIGURE_SOURCE_MANIFEST.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(manifest[0]),delimiter="\t");w.writeheader();w.writerows(manifest)
    (OUT/"FIGURE_LEGENDS.md").write_text("# Authorable figure legends\n\nThese are audited historical-cohort summaries; submission-critical unresolved analyses are documented in the manuscript/audits.\n\n"+"\n\n".join(f"## Figure {k}\n\n{v}" for k,v in legends.items())+"\n")
    (OUT/"FIGURE_RENDER_METADATA.json").write_text(json.dumps({
        "software":{"Python":os.sys.version.split()[0],"matplotlib":matplotlib.__version__,"numpy":np.__version__,"pandas":pd.__version__},
        "figures":metadata,"source_count":len(manifest),
        "limitation":"Numerical rendering does not certify biological validity of unresolved metadata/provenance; visual inspection is documented separately."},indent=2)+"\n")


if __name__=="__main__":main()
