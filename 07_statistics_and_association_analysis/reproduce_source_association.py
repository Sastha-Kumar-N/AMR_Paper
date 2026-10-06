#!/usr/bin/env python3
"""Reconstruct source association from a common Panaroo matrix, with strict labels.

Outputs are confined to this script's audit directory. Historical positional
mislabeling is reconstructed only to identify the observed erroneous output.
No random resampling is needed for the deterministic Fisher/BH analyses.
"""
import sys, os
sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
# Locally installed optional dependencies; originals are never modified.
sys.path.insert(0,str(HERE/'python_packages'))
import argparse, csv, json, platform, importlib.metadata
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control
from statsmodels.stats.multitest import multipletests
from audit_utils import ROOT, OUT, sha, log, jsonout

def save(df, name):
    df.to_csv(OUT/name,sep='\t',index=False,float_format='%.17g',na_rep='NA',lineterminator='\n')
    log(f'WRITE {name}: {len(df)} rows')

def load_matrix(path):
    with path.open() as f: hdr=next(csv.reader(f,delimiter='\t'))
    assert hdr[0]=='Gene' and len(set(hdr))==len(hdr), 'duplicate matrix samples'
    m=pd.read_csv(path,sep='\t',dtype={c:('str' if c=='Gene' else 'uint8') for c in hdr},keep_default_na=False).set_index('Gene')
    assert m.index.is_unique and np.isin(m.to_numpy(),[0,1]).all()
    return m

def calculate(ids,a,b,nc=1350,ne=99):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    t=pd.DataFrame({'cluster':ids,'clinical_present':a,'clinical_absent':nc-a,'clinical_total':nc,'clinical_prevalence':a/nc,'environmental_present':b,'environmental_absent':ne-b,'environmental_total':ne,'environmental_prevalence':b/ne})
    t['delta']=a/nc-b/ne;t['combined_present']=a+b;t['combined_absent']=nc+ne-a-b
    t['combined_minor_state_count']=np.minimum(a+b,nc+ne-a-b)
    cache={}
    for ai,bi in sorted(set(zip(a,b))): cache[(ai,bi)]=fisher_exact([[int(ai),int(nc-ai)],[int(bi),int(ne-bi)]],alternative='two-sided')
    t['fisher_odds_ratio']=[cache[(ai,bi)].statistic for ai,bi in zip(a,b)]
    t['fisher_p']=[cache[(ai,bi)].pvalue for ai,bi in zip(a,b)]
    assert np.isfinite(t.fisher_p).all()
    t['q_A_all']=multipletests(t.fisher_p,method='fdr_bh')[1]
    eligible=t.combined_minor_state_count>=10
    t['eligible_B']=eligible;t['q_B_minor_ge10']=np.nan
    t.loc[eligible,'q_B_minor_ge10']=multipletests(t.loc[eligible,'fisher_p'],method='fdr_bh')[1]
    t['final_A']=(t.q_A_all<.05)&(t.delta.abs()>=.10)
    t['final_B']=(t.q_B_minor_ge10<.05)&(t.delta.abs()>=.10)
    t['direction']=np.where(t.delta>0,'Clinical-enriched',np.where(t.delta<0,'Environmental-enriched','Equal'))
    log(f'Fisher: {len(t)} cluster hypotheses, {len(cache)} distinct count pairs cached; duplicates retained in BH')
    return t

def summary(t,label,q,eligible,valid,evidence,threshold=.10):
    sig=t[q]<.05; final=sig&(t.delta.abs()>=threshold)
    return dict(analysis=label,total_clusters=len(t),total_tested=int(eligible.sum()),bh_family_size=int(eligible.sum()),q_lt_0_05=int(sig.sum()),q_and_effect=int(final.sum()),clinical_enriched=int((final&(t.delta>0)).sum()),environmental_enriched=int((final&(t.delta<0)).sum()),effect_threshold=threshold,scientific_status=valid,evidence=evidence)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matrix',type=Path,default=ROOT/'14.pangenome_panaroo/results/whole/gene_presence_absence.Rtab')
    ap.add_argument('--metadata',type=Path,default=ROOT/'11.metadata/metadata/final_genome_metadata_classified.tsv')
    args=ap.parse_args();log('COMMAND '+repr(sys.argv))
    meta=pd.read_csv(args.metadata,sep='\t',dtype=str,keep_default_na=False)
    assert meta.assembly_accession.is_unique
    assert meta.assembly_accession.str.fullmatch(r'GC[AF]_\d+\.\d+').all(), 'unexpected accession syntax; do not normalize silently'
    labels=meta.set_index('assembly_accession').isolate_class
    mat=load_matrix(args.matrix)
    assert mat.shape==(43001,1595)
    assert set(mat.columns)==set(labels.index),'missing matrix or metadata samples'
    labels=labels.loc[mat.columns]
    assert labels.value_counts().to_dict()=={'Clinical':1350,'Unknown_Other':146,'Environmental':99}
    manifest=pd.read_csv(ROOT/'13.pangenome_inputs/pangenome_input_manifest.tsv',sep='\t',dtype=str).set_index('metadata_accession')
    assert manifest.index.is_unique and set(manifest.index)==set(labels.index)
    assert (manifest.loc[labels.index,'classification']==labels).all()
    cl=labels.index[labels=='Clinical'].tolist();en=labels.index[labels=='Environmental'].tolist()
    for run,samples in [('clinical',cl),('environmental',en)]:
        with (ROOT/f'14.pangenome_panaroo/results/{run}/gene_presence_absence.Rtab').open() as f: cols=next(csv.reader(f,delimiter='\t'))[1:]
        assert len(cols)==len(set(cols)) and set(cols)==set(samples)
    save(pd.DataFrame({'matrix_column':mat.columns,'metadata_accession':mat.columns,'source':labels.values,'matrix_position_zero_based':np.arange(len(mat.columns))}),'sample_mapping.tsv')
    a=mat[cl].sum(axis=1).to_numpy();b=mat[en].sum(axis=1).to_numpy()
    t=calculate(mat.index,a,b)
    t['whole_present']=mat.sum(axis=1).to_numpy();t['whole_total']=len(mat.columns)
    t['unknown_other_present']=t.whole_present-t.combined_present
    ann=pd.read_csv(args.matrix.with_name('gene_presence_absence.csv'),usecols=['Gene','Non-unique Gene name','Annotation'],keep_default_na=False).rename(columns={'Gene':'cluster','Non-unique Gene name':'non_unique_gene_name','Annotation':'annotation'})
    assert ann.cluster.is_unique and set(ann.cluster)==set(t.cluster)
    t=t.merge(ann,on='cluster',validate='one_to_one',sort=False)
    save(t,'all_cluster_recalculation.tsv')
    for name,flag,q in [('A','final_A','q_A_all'),('B','final_B','q_B_minor_ge10')]:save(t[t[flag]].sort_values([q,'cluster'],kind='stable'),f'analysis_{name}_final_genes.tsv')
    # Exact reproduction of the recovered numerical fingerprint. This is an
    # explicitly invalid labeling scenario, not an alternative biological test.
    known=labels.index[labels.isin(['Clinical','Environmental'])].tolist()
    assumed=np.array(['Clinical']*1350+['Environmental']*99)
    actual=labels.loc[known].to_numpy()
    save(pd.DataFrame({'sample':known,'actual_source':actual,'positional_source':assumed,'mislabeled':actual!=assumed}),'historical_positional_label_audit.tsv')
    h=calculate(mat.index,mat[known[:1350]].sum(axis=1).to_numpy(),mat[known[1350:]].sum(axis=1).to_numpy())
    histpath=ROOT/'14.pangenome_panaroo/results/publication_clinical_vs_environmental/tables/all_gene_associations.tsv'
    old=pd.read_csv(histpath,sep='\t',keep_default_na=False).set_index('cluster').loc[mat.index]
    checks={}
    for col in ['clinical_present','environmental_present']:
        checks[f'historical_{col}_equals_positional_all_rows']=bool(np.array_equal(h[col],old[col]))
        checks[f'historical_{col}_differs_from_correct_rows']=int((t[col].to_numpy()!=old[col].to_numpy()).sum())
    eligible=h.eligible_B
    for newcol,oldcol in [('fisher_p','fisher_p_value'),('q_B_minor_ge10','fdr_bh')]:
        x=h.loc[eligible,newcol].to_numpy();y=old.loc[eligible.to_numpy(),oldcol].to_numpy(dtype=float)
        checks[f'historical_{oldcol}_max_abs_error']=float(np.max(np.abs(x-y)))
        checks[f'historical_{oldcol}_allclose']=bool(np.allclose(x,y,rtol=1e-10,atol=1e-15))
    saved27=pd.read_csv(histpath.with_name('publication_significant_gene_clusters.tsv'),sep='\t',keep_default_na=False)
    checks['historical_final_cluster_set_equal']=set(h.loc[h.final_B,'cluster'])==set(saved27.cluster)
    checks['historical_eligible_set_equal']=set(h.loc[h.eligible_B,'cluster'])==set(old.index[old.tested.astype(str).str.lower()=='true'])
    checks['historical_27_overlap_correct_A']=len(set(saved27.cluster)&set(t.loc[t.final_A,'cluster']))
    checks['historical_27_overlap_correct_B']=len(set(saved27.cluster)&set(t.loc[t.final_B,'cluster']))
    checks['historical_matrix_total_counts_match_all_rows']=bool(np.array_equal(h.combined_present,old.total_present))
    assert all(checks[k] for k in checks if k.endswith(('all_rows','allclose','set_equal')))
    saved27['audit_status']='INVALID_POSITIONAL_SOURCE_LABELS'
    correct=t.set_index('cluster')
    for c in ['clinical_present','environmental_present','delta','fisher_p','q_A_all','q_B_minor_ge10']:
        saved27['correct_'+c]=saved27.cluster.map(correct[c])
    save(saved27,'historical_27_if_recovered.tsv')
    save(h,'historical_positional_recalculation.tsv')
    variants=[summary(t,'A_correct_labels_BH_all','q_A_all',pd.Series(True,index=t.index),'Valid unadjusted exploratory association','15.pangenome_comparison/scripts/compare_clinical_environmental_from_whole_rtab.py'),summary(t,'B_correct_labels_minor_ge10_BH','q_B_minor_ge10',t.eligible_B,'Valid conditional low-information filtering; exploratory','14.pangenome_panaroo/results/publication_clinical_vs_environmental/manuscript_methods.md'),summary(h,'C_historical_positional_labels_minor_ge10_BH','q_B_minor_ge10',h.eligible_B,'INVALID source labels; forensic reconstruction only',str(histpath.relative_to(ROOT))),summary(t,'C_original_15_comparison_25pp_category_rule','q_A_all',pd.Series(True,index=t.index),'Different reporting effect threshold; same Fisher family','15.pangenome_comparison/scripts/compare_clinical_environmental_from_whole_rtab.py:53-89',.25)]
    save(pd.DataFrame(variants),'statistical_variant_summary.tsv')
    order=np.argsort(t.fisher_p.to_numpy(),kind='stable');p=t.fisher_p.to_numpy();manual=np.empty(len(p));manual[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    checks.update(dict(q_A_statsmodels_vs_scipy_max_abs_error=float(np.max(np.abs(t.q_A_all-false_discovery_control(p,method='bh')))),q_A_statsmodels_vs_manual_max_abs_error=float(np.max(np.abs(t.q_A_all-manual))),nonfinite_fisher_p=int((~np.isfinite(p)).sum()),duplicate_cluster_ids=int(t.cluster.duplicated().sum()),all_zero_compared=int((t.combined_present==0).sum()),all_one_compared=int((t.combined_absent==0).sum()),excluded_B=int((~t.eligible_B).sum()),eligible_B=int(t.eligible_B.sum()),mislabeled_samples=int((actual!=assumed).sum()),true_environmental_in_positional_clinical=int(((actual=='Environmental')&(assumed=='Clinical')).sum()),true_clinical_in_positional_environmental=int(((actual=='Clinical')&(assumed=='Environmental')).sum()),correct_labels_match_manifest=True,correct_labels_match_independent_run_headers=True,unknown_other_excluded=146,metadata_only_samples=0,matrix_only_samples=0,duplicate_matrix_samples=0,duplicate_metadata_samples=0,name_transformations=0))
    # All-zero/all-one rows receive p=1. A includes them; B excludes them.
    const=(t.combined_present==0)|(t.combined_absent==0)
    assert (t.loc[const,'fisher_p']==1).all()
    checks['constant_rows_final_A']=int(t.loc[const,'final_A'].sum())
    checks['constant_rows_final_B']=int(t.loc[const,'final_B'].sum())
    # Effect of removal is documented as a diagnostic, not another proposed analysis.
    q_nonconstant=multipletests(t.loc[~const,'fisher_p'],method='fdr_bh')[1]
    checks['diagnostic_only_exclude_constant_BH_family']=int((~const).sum())
    checks['diagnostic_only_exclude_constant_BH_significant']=int((q_nonconstant<.05).sum())
    checks['diagnostic_only_exclude_constant_final']=int(((q_nonconstant<.05)&(t.loc[~const,'delta'].abs()>=.1)).sum())
    legacy=pd.read_csv(ROOT/'15.pangenome_comparison/results/clinical_environmental_gene_frequency_table.tsv',sep='\t').set_index('gene_cluster').loc[t.cluster]
    checks['original_15_counts_equal_correct']=bool(np.array_equal(legacy.clinical_count,t.clinical_present) and np.array_equal(legacy.environmental_count,t.environmental_present))
    checks['original_15_q_max_abs_error']=float(np.max(np.abs(legacy.fisher_fdr.to_numpy()-t.q_A_all)))
    def original_category(r):
        cf=r.clinical_prevalence;ef=r.environmental_prevalence;wf=r.whole_present/r.whole_total
        if cf>=.99 and ef>=.99:return 'conserved_in_both_core'
        if r.clinical_present>0 and r.environmental_present==0 and cf>=.15:return 'clinical_specific_common'
        if r.environmental_present>0 and r.clinical_present==0 and ef>=.15:return 'environmental_specific_common'
        if r.q_A_all<=.05 and r.delta>=.25 and cf>=.15:return 'clinical_enriched'
        if r.q_A_all<=.05 and r.delta<=-.25 and ef>=.15:return 'environmental_enriched'
        if cf>=.95 and ef<.95:return 'clinical_softcore_not_environmental'
        if ef>=.95 and cf<.95:return 'environmental_softcore_not_clinical'
        if wf<.15:return 'rare_cloud'
        return 'shared_accessory_or_variable'
    categories=pd.Series([original_category(r) for r in t.itertuples()])
    checks['original_categories_equal_recomputed']=bool(np.array_equal(categories,legacy.comparison_category))
    assert checks['original_categories_equal_recomputed']
    counts=categories.value_counts().to_dict();jsonout('category_counts.json',counts)
    checks['software']={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','statsmodels']};checks['python']=platform.python_version()
    checks['inputs']={str(p.relative_to(ROOT)):sha(p) for p in [args.matrix,args.metadata,ROOT/'13.pangenome_inputs/pangenome_input_manifest.tsv',histpath,args.matrix.with_name('gene_presence_absence.csv')]}
    jsonout('reconstruction_checks.json',checks)
    print(pd.DataFrame(variants).to_string(index=False));print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
