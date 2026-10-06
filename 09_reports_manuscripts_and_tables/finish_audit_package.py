#!/usr/bin/env python3
"""Assemble concise derivative tables and documentation inside the audit only."""
import sys
sys.dont_write_bytecode=True
import csv,json,re
from pathlib import Path
from audit_utils import ROOT,OUT,tsv,log

def records(n):return list(csv.DictReader((OUT/n).open(),delimiter='\t'))
allr=records('all_cluster_recalculation.tsv');by={r['cluster']:r for r in allr}
extra=['group_24031','group_11469','eta_1~~~eta_2~~~eta']
final=records('analysis_B_final_genes.tsv')
cases=final+[by[k] for k in extra if k not in {r['cluster'] for r in final}]
tsv('gene_case_studies_primary.tsv',cases)
def number(x):return 'NA' if x in ['','NA'] else f'{float(x):.6g}'
rows=['| Whole-run cluster (short label where indicated) | Clinical / 1,350 | Environmental / 99 | Delta (pp) | Fisher P | qA | qB |',
      '|---|---:|---:|---:|---:|---:|---:|']
for r in cases:
 name=r['cluster']
 if name.startswith('xerC_'):name='xerC/xerD cluster (full ID in A table)'
 elif name.startswith('sul1'):name='sul1~~~sul1_3'
 elif name.startswith('merT'):name='merT cluster (full ID in B table)'
 elif name.startswith('merP'):name='merP cluster (full ID in B table)'
 rows.append('| '+ ' | '.join([name,r['clinical_present'],r['environmental_present'],f"{float(r['delta'])*100:.5f}",number(r['fisher_p']),number(r['q_A_all']),number(r['q_B_minor_ge10'])])+' |')
p=OUT/'FINAL_AUDIT_REPORT.md';s=p.read_text()
if '<!-- GENE_TABLE -->' in s:p.write_text(s.replace('<!-- GENE_TABLE -->','\n'.join(rows)))

e=[]
def finding(name,num,interp,strength,caveat,section,sources,verification='New independent local calculation'):
 e.append(dict(finding=name,verified_numeric_evidence=num,biological_interpretation=interp,strength_of_evidence=strength,
   important_caveat=caveat,manuscript_section=section,source_paths=sources,verification_scope=verification))
whole='14.pangenome_panaroo/results/whole/gene_presence_absence.Rtab'
clinical='14.pangenome_panaroo/results/clinical/gene_presence_absence.Rtab'
func='17.clinical_comparative_genomics/funscan_analysis/'
amr='18.amrfinderplus/results/clinical_amrfinderplus/combined_amrfinderplus_results.tsv'
finding('Traceable assembly curation','3765 -> 1652 -> 1630 -> 1619 -> 1604 -> 1595; 9 P. paraeruginosa excluded',
 'Deduplication, assembly QC and taxonomy protect frequency and diversity estimates','HIGH',
 'Assembly representations are not guaranteed independent isolates; three BioSample pairs remain','Methods; Results; Discussion',
 'scripts /unique fasta files/unique_fasta_selection_report.csv;scripts /final unique file/final_unique_genome_selection_report.csv;quast_final_filter_report.csv;8.checkm2_filtered_genomes/reports/checkm2_final_filter_report.csv;9.gtdbtk_analysis/results/checkm2_pass_gtdbtk/classify/gtdbtk.bac120.summary.tsv')
finding('Validated common-space sample mapping','1595 = 1350 Clinical + 99 Environmental + 146 Unknown_Other; zero missing/extra/duplicate accession records',
 'Source comparison can be made within one clustering space','HIGH','Matching accuracy does not establish source-category accuracy or isolate independence','Methods; Results',
 whole+';11.metadata/metadata/final_genome_metadata_classified.tsv;13.pangenome_inputs/pangenome_input_manifest.tsv')
finding('Conserved backbone and expansive catalogue','43001 clusters; core>=99% 4551; soft-core 673; shell 1183; cloud 36594; strict-core 2352',
 'Shared backbone coexists with extensive accessory diversity','HIGH','Operational clustering and frequency definitions; no gain/transfer mechanism inferred','Abstract; Results; Discussion; Conclusion',whole)
finding('Per-genome cluster burden','median 6101; IQR 5773.5-6332; range 5304-7279; mean 6081.7774',
 'The species catalogue is distributed across genomes rather than contained within each genome','HIGH','Clusters are not necessarily one-to-one CDS or copy counts','Results',whole)
finding('Rare catalogue dominance','12293 singletons; 4400 doubletons; 22567 clusters <=5 genomes; 27754 <1% of genomes',
 'Rare clusters form a long tail of observed diversity','HIGH','Rare clusters are sensitive to annotation and clustering errors; rarity does not prove acquisition','Abstract; Results; Discussion',whole)
finding('Richness versus occupancy','cloud=85.1003% of clusters but 7.2609% of presences; core=10.5835% of clusters but 74.7380% of presences',
 'Rare diversity dominates catalogue richness while core content dominates individual repertoires','HIGH','Pooled presence assignments are not biological replicates','Results; Discussion; Conclusion',whole)
finding('Continued pangenome discovery','whole exact gamma 0.2900023644; terminal expected gain 7.707210; clinical gamma 0.2884624',
 'Gene discovery remains unsaturated within the observed collection','HIGH','Conditional finite-collection curve; not an acquisition rate or extrapolated infinite repertoire','Abstract; Results; Discussion; Conclusion',whole+';'+clinical)
finding('Heaps fit limitation','whole R2 0.9999772; first point underpredicted 16.6313%',
 'High aggregate fit quality does not validate all sampling scales','HIGH','Successive curve points are dependent; bounded fit is not a null test of openness','Methods; Discussion',whole+';publication_manuscript/reanalysis_scripts/audit_exact_extension.py')
finding('Clinical genomic fluidity','910575 pairs; mean 0.1073617094; median 0.1090704044; IQR 0.0887249-0.1275759',
 'Pairwise content variability complements global catalogue richness','HIGH','Dependent pairs; denominator includes all clusters; not Jaccard or transfer rate','Results; Discussion; Conclusion',clinical)
finding('Fluidity beyond the rarest clusters','mean after >=2 occupancy 0.1061258; >=5 0.1037009; >=14 0.0992106',
 'Variation is not solely driven by singletons','HIGH','Filtering changes both numerator and denominator; conditional on clinical sample','Results; Discussion',clinical)
finding('Shared clinical/environmental content','4174 core-level in both; 17872 detected in both; 2179 category shared_accessory_or_variable; env soft-core 29; clinical soft-core 25',
 'Extensive shared genomic framework with frequency variation','HIGH','Categories are ordered descriptive rules; sample-size imbalance affects observed sharing','Abstract; Results; Discussion',whole+';15.pangenome_comparison/scripts/compare_clinical_environmental_from_whole_rtab.py')
finding('Primary source associations','BH over 43001: 4 significant and effect-qualified; 2 clinical / 2 environmental',
 'Limited detectable frequency shifts under recorded source labels','MODERATE','No phylogenetic adjustment; metadata uncertainty; not clinical causation or specificity','Abstract; Results; Methods',whole+';audit_clinical_environmental_provenance/analysis_A_final_genes.tsv')
finding('Filtering sensitivity','minor count>=10: 12605 tests; 12 q<0.05; 9 effect-qualified; 7 clinical / 2 environmental',
 'Testing-family definition changes discovery threshold','MODERATE','Historical use verified; selection before examining results unverified','Results; Methods; Supplement',whole+';audit_clinical_environmental_provenance/analysis_B_final_genes.tsv')
finding('Historical numerical error resolved','178 wrongly assigned genomes; all 43001 historical source counts reproduced; 60 q-significant / 27 final; zero overlap with corrected sets',
 'Historical source interpretation must be retired','HIGH','Original generator missing; exact internal coding operation unverified','Audit; Manuscript corrections',
 '14.pangenome_panaroo/results/publication_clinical_vs_environmental/tables/all_gene_associations.tsv;audit_clinical_environmental_provenance/reconstruction_checks.json')
finding('Balanced output recovered','1000 draws of 99 without replacement; seed 20260721; 27 saved rows exactly reconstructed',
 'Equal numbers did not remedy wrong pools or adjust ancestry','HIGH','Invalid positional labels; balancing was a candidate display/stability exercise','Audit; Manuscript corrections',
 '14.pangenome_panaroo/results/publication_clinical_vs_environmental/tables/downsampling_sensitivity.tsv;audit_clinical_environmental_provenance/balanced_reconstruction_checks.json')
finding('Current toxA/eta occupancy','1321/1350 versus 95/99; delta 0.01892256; P 0.2782807; qA 1; qB 0.8014001',
 'No significant clinical enrichment in the current common matrix','HIGH','Existing local annotation only; no sequence, activity or mechanistic inference','Results correction; Supplement',whole)
finding('Clinical SNP-profile diversity','1350 samples; 130103 complete variable sites; 1195 distinct retained profiles',
 'Clinical collection is heterogeneous rather than one identical sampled clone','MODERATE',
 'Existing QC; clinical-only tree; bootstrap nonconvergence; no validated recombination masking or lineage assignment','Results; Discussion',
 'publication_manuscript/reconstructed_tables/verified_phylogeny_geography.json;17.clinical_comparative_genomics/core_genome_phylogeny/03_iqtree_clinical_only/clinical_only_complete_core_snps.log',
 'Verified existing local QC and tree log; no new sequence analysis')
finding('Accessory distance error','first two columns: true shared 837, uint8 shared 69; true Jaccard 0.442743, erroneous 0.969604; historical silhouette 0.0111432',
 'Saved distance/PCoA/silhouette cannot support ecological clustering conclusions','HIGH','Correct full ordination was not newly performed','Discussion limitation; Figure corrections',
 '14.pangenome_panaroo/results/whole/publication_deep_analysis/01_scripts/common.py;audit_clinical_environmental_provenance/integrated_verification.json')
finding('Gene-neighbourhood modularity','40811 nodes; 68454 edges; 260 components; 372 communities; weighted modularity 0.9707741',
 'Strong partition structure in a gene-neighbourhood graph','MODERATE',
 'Resolution/weight/component dependent; not pathways or ecological populations; modularity not rerun in continuation','Results; Discussion',
 '14.pangenome_panaroo/results/graph files pangenome/results_louvain_fdr/statistics_summary.json;14.pangenome_panaroo/results/graph files pangenome/results_louvain_fdr/statistics_validation.json',
 'Existing validated analysis plus independent node/community-table crosscheck')
finding('Standalone AMR units and burden','18058 AMR rows; 852 STRESS rows; 929 symbols=413 nonpoint+516 point/disruption; median symbols 10, IQR 8-16, range 4-35',
 'Resistance-associated repertoire has widespread and variable components','HIGH','Symbols, families, loci and raw calls are different units; not phenotype','Results; Discussion',amr)
finding('Layered resistance-associated distribution','fosA 1350; aph(3\')-IIb 1346; catB7 1317; blaPDC-prefix 1347; sul1 485; blaVIM-prefix 124',
 'Background components coexist with variably distributed annotations','MODERATE','fosA includes 90 INTERNAL_STOP and 1 PARTIALP; detection is not intact function or acquired origin','Results; Discussion; Conclusion',amr)
finding('AMP prediction support imbalance','40061 occurrences; 38699 single-predictor; 1362 higher-support; saved family counts 1051/1077/1147',
 'Predicted functional heterogeneity requires evidence-tier qualification','EXPLORATORY','No activity, secretion or therapeutic claims; no new sequence analysis','Results; Discussion',
 func+'02_cleaned/AMP_master_long.tsv;publication_manuscript/reconstructed_tables/verified_functional_results.json')
finding('CAZyme annotated repertoire','151816 recommended family/subfamily tokens; 136278 genome-gene identifiers; 73 families; 6 classes; median 112 tokens, 101 genes',
 'Common functional classes coexist with variable annotated family burdens','MODERATE','Tokens are not physical domain intervals or demonstrated metabolism','Results; Discussion',func+'02_cleaned/CAZyme_master_long.tsv;'+func+'03_matrices/genome_functional_burden_matrix.tsv')
finding('BGC integrity limitation','1355 of 1361 audited antiSMASH export rows have discontinuous ordinal sets; median 444 listed CDS / 22090 bp region',
 'BGC repertoire and dependent associations cannot yet support biological conclusions','HIGH','Saved integrity audit; this does not establish absence of BGCs','Limitations; withheld results',
 'publication_manuscript/reconstructed_tables/verified_antismash_export_audit.tsv;publication_manuscript/reanalysis_scripts/audit_additional.py')
finding('AMP-ARG co-distribution','138 AMP x 83 ARG =11454 tests; 1649 q<0.05; 494 also |phi|>=0.20',
 'Functional features show exploratory genome-wide co-occurrence','EXPLORATORY','No population adjustment; not physical linkage, biochemical interaction or causation','Results; Discussion',
 func+'03_matrices/AMP_presence_absence.tsv;'+func+'03_matrices/AMR_presence_absence.tsv;'+func+'04_statistics/AMP_ARG_associations.tsv')
finding('Burden confounding sensitivity','AMP-ARG Spearman 0.4815997; CDS-adjusted rank-residual correlation -0.0696633',
 'Aggregate relationship is highly sensitive to gene-content/annotation opportunity','MODERATE','Prokka CDS proxy differs from functional gene calling; no ancestry adjustment','Results; Discussion',func+'03_matrices/genome_functional_burden_matrix.tsv')
finding('qac/sul physical proximity','436 co-positive assemblies; 436 same-contig nearest pairs; gap 0-1950 bp; 433 zero-gap',
 'Local physical proximity is supported separately from genome-wide correlation','MODERATE','At least one pair per genome; no mobile-element identity, cotransfer or interaction established','Results; Discussion',amr)
finding('BGC-dependent associations withheld','saved AMP-BGC 70104 tests / 9470 q-significant; ARG-BGC 42164 / 3983; biological status WITHHELD',
 'Numerical significance cannot repair invalid feature provenance','EXPLORATORY','Unresolved BGC inputs; do not use counts as validated biological discoveries','Limitations',
 func+'04_statistics/AMP_BGC_associations.tsv;'+func+'04_statistics/ARG_BGC_associations.tsv',
 'Existing output counts inspected; no new inferential validation of BGC-dependent results')
finding('Ecological and surveillance model','Integrated common-backbone, openness, fluidity, source-overlap and functional evidence',
 'Consistent with ecological generalism and a shared genomic framework; surveillance should combine lineage and accessory information','MODERATE',
 'Interpretation, not causal proof; ecological substructure and host adaptation unresolved','Discussion; Conclusion',
 'audit_clinical_environmental_provenance/integrated_verification.json;audit_clinical_environmental_provenance/FINAL_AUDIT_REPORT.md',
 'Synthesis of verified local observations')
tsv('key_scientific_findings.tsv',e)
tsv('scientific_evidence_registry.tsv',[{'finding':x['finding'],'source_paths':x['source_paths'],'verification_scope':x['verification_scope']} for x in e])

readme='''# Clinical/environmental provenance and whole-genome scientific audit

Project root: `/mnt/phoenix_3/Backup_main_Aug_2026/phoenix_2`.
All generated output is in this directory; original matrices, metadata, scripts and manuscripts were read only. This is an audit, not a Panaroo rerun or manuscript edit.

## Read first

- [FINAL_AUDIT_REPORT.md](FINAL_AUDIT_REPORT.md): causal 27/9/4 explanation, scientific/statistical decisions, integrated interpretation and corrections.
- [MANUSCRIPT_CONCLUSION.md](MANUSCRIPT_CONCLUSION.md): manuscript-ready conclusion.
- [DISCUSSION_SCIENTIFIC_INTERPRETATION.md](DISCUSSION_SCIENTIFIC_INTERPRETATION.md): structured discussion with evidence labels.
- [key_scientific_findings.tsv](key_scientific_findings.tsv): evidence, interpretation, caveats, source paths and strength.

Decision: **4 primary; 9 filtering sensitivity; 27 numerically reproducible but invalid because 178 source labels were assigned to the wrong matrix columns.** Filtering was historically used; its selection before examination of results is not independently established. All source associations remain unadjusted for population structure.

## Tables and provenance

`matrix_inventory.tsv` covers 24 matrix exports/copies and five distinct binary matrices. `script_inventory.tsv` distinguishes configured analyses from executed outputs. `evidence_provenance.tsv` maps path-hash keys, original PDFs, tables and manuscript evidence. `manuscript_conflicts.tsv` contains passage-level document decisions; `manuscript_file_summary.tsv` groups them by file.

`all_cluster_recalculation.tsv` is the complete 43,001-row correct-label table. `analysis_A_final_genes.tsv` and `analysis_B_final_genes.tsv` contain four and nine rows. `statistical_variant_summary.tsv` includes evidenced historical variants. `historical_27_if_recovered.tsv` preserves historical values with explicit INVALID status and correct-label comparators. It is forensic evidence, not a discovery table. `gene_case_studies_primary.tsv` is a compact whole-matrix view; `gene_case_studies.tsv` retains extended run-qualified cases. `independent_run_id_audit.md` explains why arbitrary IDs cannot connect independent Panaroo runs.

`reconstruction_checks.json`, `balanced_reconstruction_checks.json`, `additional_integrity_checks.json`, and `integrated_verification.json` record numerical checks. `integrated_source_manifest.tsv` and `final_validation.json` record source integrity and delivery checks. `commands_and_checks.log` records searches, numerical commands and validations. Detailed document text and search intermediates are retained rather than deleted.

## Reproduce without modifying source files

Run from the project root. The existing graph environment supplies SciPy; statsmodels was installed only into this audit directory because it was absent. No system environment was modified. The selected interpreter inherits the system NumPy/pandas installation. Versions: Python 3.13.5, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.0, statsmodels 0.15.0. Package artifacts and caches are confined here.

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 '14.pangenome_panaroo/results/graph files pangenome/.venv/bin/python' -B audit_clinical_environmental_provenance/reproduce_source_association.py
PYTHONDONTWRITEBYTECODE=1 python3 -B audit_clinical_environmental_provenance/audit_balanced_resampling.py
OPENBLAS_NUM_THREADS=20 PYTHONDONTWRITEBYTECODE=1 '14.pangenome_panaroo/results/graph files pangenome/.venv/bin/python' -B audit_clinical_environmental_provenance/verify_integrated_results.py
```

The source-association script validates dimensions, exact accessions and group sizes; computes prevalence, raw two-sided Fisher OR/P and both BH families; exports deterministic TSVs; and reproduces the historical error as an explicitly invalid comparator. Identical contingency tables are cached, not removed from BH. The balanced script uses seed 20260721 and reconstructs only the invalid historical pools to establish provenance. It must not be presented as a corrected biological sensitivity.

The integrated script independently checks architecture, exact finite-collection rarefaction, all clinical fluidity pairs and rare-cluster sensitivities, curation tables, AMR aggregates, AMP-ARG tests, burden sensitivity and existing caller-coordinate proximity. Binary intersection arithmetic uses a safe dtype and validates integer diagonals. This continuation performs no sequence/protein/toxin/homology work and accesses no external biological source. Existing phylogeny QC and validated graph results are cited with their distinct verification scope.

## Boundaries and remaining work

The original historical comparison generator and complete original Panaroo invocation history are not recovered. Source metadata and isolate independence need curation; no lineage-adjusted source association has been demonstrated. Older overflowed accessory distances and BGC-dependent results are withheld. The 27-source error is resolved despite these separate gaps.

Earlier sequence-related scratch files and crosswalks are retained as prior audit artifacts. They are not inputs to the integrated biological interpretation or 27/9/4 decision; no further sequence analysis was performed after the scope change. Do not rerun `audit_cluster_identity.py` as part of this local computational continuation.
'''
(OUT/'README.md').write_text(readme)
log('WRITE final narrative package, compact gene table, scientific evidence and findings tables')
print('Assembled',len(e),'scientific findings and',len(cases),'primary gene cases')
