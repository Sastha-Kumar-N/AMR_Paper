"""Read original results only. Exact finite-collection rarefaction and bounded sensitivity.
Run from project root with PYTHONDONTWRITEBYTECODE=1 python3 -B <this script>.
All outputs are below article_20260911/. No annotation/pangenome rerun occurs.
"""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np,pandas as pd
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript_reconstruction'
sys.path.append(str(ROOT/'17.clinical_comparative_genomics/map_genome/map_pan_3/map_pan_3/.venv/lib/python3.13/site-packages'))
from scipy import sparse
from scipy.optimize import curve_fit
R={};S={}
def read(p,**kw):
 p=ROOT/p;h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 S[str(p.relative_to(ROOT))]={'sha256':h.hexdigest(),'bytes':p.stat().st_size}
 return pd.read_csv(p,sep='\t',**kw)
def desc(x):
 x=np.asarray(x,float);return dict(zip(['n','min','q1','median','q3','max','mean'],[len(x),*np.quantile(x,[0,.25,.5,.75,1]),x.mean()]))
def exact_curve(f,n):
 hist=np.bincount(f,minlength=n+1).astype(float);hist[0]=0;freq=np.arange(n+1);absent=np.ones(n+1);core=np.ones(n+1);rows=[]
 for k in range(1,n+1):
  gain=float(hist@(absent*freq/(n-k+1)))
  absent*=np.maximum(n-freq-k+1,0)/(n-k+1)
  core*=np.maximum(freq-k+1,0)/(n-k+1)
  rows.append([k,float(hist@(1-absent)),float(hist@core),gain])
 t=pd.DataFrame(rows,columns=['n_genomes','expected_pangenome','expected_strict_core','expected_new_clusters'])
 assert np.isclose(t.expected_pangenome.iloc[-1],(f>0).sum())
 assert np.isclose(t.expected_new_clusters.iloc[-1],(f==1).sum()/n)
 pars,_=curve_fit(lambda x,k,g:k*x**g,t.n_genomes,t.expected_pangenome,p0=[5000,.3],bounds=(0,np.inf),maxfev=20000)
 fitted=pars[0]*t.n_genomes.to_numpy()**pars[1];res=t.expected_pangenome-fitted
 t['fitted_pangenome']=fitted
 nd=t.iloc[1:];pars2,_=curve_fit(lambda x,k,a:k*x**(-a),nd.n_genomes,nd.expected_new_clusters,p0=[1000,.5],bounds=(0,np.inf),maxfev=20000)
 pred2=pars2[0]*nd.n_genomes.to_numpy()**(-pars2[1]);r2=1-float((res**2).sum()/((t.expected_pangenome-t.expected_pangenome.mean())**2).sum())
 return t,{'kappa':float(pars[0]),'gamma':float(pars[1]),'R2':r2,'RMSE':float(np.sqrt(np.mean(res**2))),'new_kappa':float(pars2[0]),'new_alpha':float(pars2[1]),'new_R2':float(1-np.sum((nd.expected_new_clusters-pred2)**2)/np.sum((nd.expected_new_clusters-nd.expected_new_clusters.mean())**2)),'endpoint_expected_new':float(t.expected_new_clusters.iloc[-1]),'n100_expected_pan':float(t.loc[t.n_genomes==100,'expected_pangenome'].iloc[0]) if n>=100 else None,'n99_expected_pan':float(t.loc[t.n_genomes==99,'expected_pangenome'].iloc[0]),'CI':'Not estimated: deterministic finite-collection expectation; no population confidence interval implied.'}
def categories(f,n):
 return {'total':int((f>0).sum()),'core99':int((f>=.99*n).sum()),'strict100':int((f==n).sum()),'soft95':int(((f>=.95*n)&(f<.99*n)).sum()),'shell15':int(((f>=.15*n)&(f<.95*n)).sum()),'cloud':int(((f>0)&(f<.15*n)).sum()),'singletons':int((f==1).sum()),'doubletons':int((f==2).sum())}
def fluidity(x):
 xx=sparse.csr_matrix(x.astype(np.int32));n=x.shape[1];count=x.sum(axis=0,dtype=np.int64);shared=(xx.T@xx).toarray();den=count[:,None]+count[None,:];p=(den-2*shared)/den;v=p[np.triu_indices(n,1)];return {'mean':float(v.mean()),'pairs':len(v),'n_genomes':n},v
exclude={'GCA_019738995.2','GCA_021497345.2','GCA_021497365.2'}
weighted=[];pergenome=[]
for pop in ['whole','clinical','environmental']:
 print('Matrix and exact rarefaction:',pop,flush=True)
 t=read(f'14.pangenome_panaroo/results/{pop}/gene_presence_absence.Rtab',index_col=0);x=t.to_numpy(dtype=np.uint8);n=x.shape[1];f=x.sum(axis=1,dtype=np.int64);s=x.sum(axis=0,dtype=np.int64);total=int(f.sum());cats={'Core ≥99%':f>=.99*n,'Soft core ≥95–<99%':(f>=.95*n)&(f<.99*n),'Shell ≥15–<95%':(f>=.15*n)&(f<.95*n),'Cloud <15%':(f>0)&(f<.15*n)}
 pp=pd.DataFrame({'assembly_accession':t.columns,'population':pop,'cluster_count':s})
 for label,mask in cats.items():
  counts=x[mask].sum(axis=0,dtype=np.int64);pp[label]=counts
  weighted.append({'population':pop,'n_genomes':n,'category':label,'distinct_clusters':int(mask.sum()),'distinct_percent':float(mask.sum()/len(f)*100),'cluster_presence_assignments':int(f[mask].sum()),'assignment_percent':float(f[mask].sum()/total*100),'median_per_genome_count':float(np.median(counts)),'median_per_genome_percent':float(np.median(counts/s*100))})
 pergenome.append(pp)
 rare=(f<=5)&(f>0)
 R[pop]={'n_genomes':n,'categories':categories(f,n),'cluster_presence_assignments':total,'per_genome_clusters':desc(s),'rare_le5_distinct':int(rare.sum()),'rare_le5_distinct_percent':float(rare.sum()/len(f)*100),'rare_le5_assignment_percent':float(f[rare].sum()/total*100),'singleton_assignment_percent':float((f==1).sum()/total*100)}
 curve,fit=exact_curve(f,n);curve.to_csv(OUT/f'source_data/{pop}_exact_rarefaction.tsv',sep='\t',index=False);R[pop]['exact_fit']=fit
 if pop!='environmental':
  keep=np.array([v not in exclude for v in t.columns]);xx=x[:,keep];ff=xx.sum(axis=1,dtype=np.int64);cc,fit2=exact_curve(ff,int(keep.sum()));cc.to_csv(OUT/f'source_data/{pop}_paired_exclusion_exact_rarefaction.tsv',sep='\t',index=False)
  R[pop]['paired_exclusion_sensitivity']={'removed':sorted(set(t.columns)&exclude),'n_genomes':int(keep.sum()),'categories_fixed_clusters':categories(ff,int(keep.sum())),'exact_fit':fit2,'caveat':'Fixed original cluster definitions; this is not a Panaroo rerun.'}
 if pop=='clinical':
  v0,dist=fluidity(x);v1,_=fluidity(xx);R[pop]['fluidity_full']=v0;R[pop]['paired_exclusion_sensitivity']['fluidity']=v1
  h,e=np.histogram(dist,bins=np.linspace(0,.25,101));pd.DataFrame({'lower':e[:-1],'upper':e[1:],'pairs':h}).to_csv(OUT/'source_data/clinical_fluidity_histogram.tsv',sep='\t',index=False)
  rows=[]
  for z in [1,2,5,14]:
   v,_=fluidity(x[f>=z]);rows.append({'minimum_occupancy':z,'retained_clusters':int((f>=z).sum()),**v})
  pd.DataFrame(rows).to_csv(OUT/'source_data/clinical_fluidity_sensitivity.tsv',sep='\t',index=False)
 del t,x
pd.DataFrame(weighted).to_csv(OUT/'source_data/occupancy_weighted_architecture.tsv',sep='\t',index=False);pd.concat(pergenome).to_csv(OUT/'source_data/per_genome_frequency_classes.tsv',sep='\t',index=False)
print('AMR sensitivity and specimen-independent accounting',flush=True)
amr=read('18.amrfinderplus/results/clinical_amrfinderplus/combined_amrfinderplus_results.tsv',keep_default_na=False);amr=amr[amr.Type.eq('AMR')&~amr['Element name'].str.contains('WILDTYPE|UNKNOWN',case=False)]
point=amr.Subtype.str.startswith('POINT')|amr.Method.str.startswith('POINT');bur=amr.groupby('Name')['Element symbol'].nunique().rename('unique_symbols').to_frame();bur['nonpoint_symbols']=amr[~point].groupby('Name')['Element symbol'].nunique();bur['point_symbols']=amr[point].groupby('Name')['Element symbol'].nunique();bur=bur.fillna(0).astype(int);bur['paired_exclusion_retained']=~bur.index.isin(exclude);bur.to_csv(OUT/'source_data/standalone_amr_burden.tsv',sep='\t')
R['AMR_paired_exclusion']={'full':desc(bur.unique_symbols),'retained':desc(bur.loc[bur.paired_exclusion_retained,'unique_symbols']),'exact_symbols_retained':int(amr.loc[~amr.Name.isin(exclude),'Element symbol'].nunique())}
# Explicit source tables for predictors and CAZyme domains; distinct evidence units.
fa=ROOT/'17.clinical_comparative_genomics/funscan_analysis'
amp=read(str((fa/'02_cleaned/AMP_master_long.tsv').relative_to(ROOT)),keep_default_na=False)
a=amp.groupby('genome_id').size().rename('AMP_occurrences').to_frame();a['AMP_higher_support']=amp[amp.high_confidence.eq('yes')].groupby('genome_id').size();a=a.fillna(0).astype(int)
cazy=read(str((fa/'02_cleaned/CAZyme_master_long.tsv').relative_to(ROOT)),keep_default_na=False);cap=cazy[cazy.analysis_tier.eq('primary_consensus')]
a['CAZyme_domain_annotations']=cap.groupby('genome_id').size();a['CAZyme_unique_genes']=cap.groupby('genome_id').sequence_id.nunique();a.to_csv(OUT/'source_data/functional_burden.tsv',sep='\t')
R['functional_burden']={k:desc(a[k]) for k in a.columns}
R['amp_specificity_warning']={'higher_support_not_activity':True,'prediction_occurrences':len(amp),'unique_sequences':amp.aa_sequence.nunique()}
R['created_utc']=datetime.now(timezone.utc).isoformat()
(OUT/'source_data/additional_validated_results.json').write_text(json.dumps(R,indent=2));(OUT/'source_data/source_checksums.json').write_text(json.dumps(S,indent=2));print(json.dumps(R,indent=2),flush=True)
