"""Rebuild all new manuscript figures from verified numeric source tables.
Original figures/results remain read-only. All plotted data are saved separately.
"""
from pathlib import Path
import os,json,hashlib,shutil
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1];BASE=OUT/'source_data/audit_snapshot'
os.environ['MPLCONFIGDIR']=str(OUT/'source_data/matplotlib_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.4,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':400})
C=json.loads((BASE/'verified_core_results.json').read_text());F=json.loads((BASE/'verified_functional_results.json').read_text());E=json.loads((OUT/'source_data/additional_validated_results.json').read_text());D=OUT/'source_data';G=OUT/'figures';G.mkdir(exist_ok=True);REG=[]
TEAL='#007C83';BLUE='#3E6C9E';ORANGE='#CC7939';PURPLE='#8264A6';GREY='#71808D';COLORS=[BLUE,TEAL,ORANGE,GREY]
def read(p):return pd.read_csv(p,sep='\t')
def panel(ax,label,title):
 ax.set_title(title,loc='left',pad=12);ax.text(-.16,1.12,label,transform=ax.transAxes,fontweight='bold',fontsize=12,va='top')
def save(fig,name,sources):
 for ext in ['png','pdf','svg']:fig.savefig(G/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
 REG.append({'figure':name,'png':f'figures/{name}.png','pdf':f'figures/{name}.pdf','svg':f'figures/{name}.svg','source_tables':[str(Path(s).relative_to(ROOT)) for s in sources],'script':str(Path(__file__).relative_to(ROOT))})
 plt.close(fig)
q=read(BASE/'verified_cohort_qc.tsv')
q[['assembly_accession','Completeness','Contamination','total_length','contigs','gc','n50']].to_csv(D/'figure1_final_QC.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(7.1,5.05),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1.18,1]);a=fig.add_subplot(gs[:,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,1]);a.set_axis_off();a.set_xlim(0,1);a.set_ylim(0,1)
panel(a,'a','Traceable assembly collection')
labels=[('3,765 input files','3,130 + 635 preserved files'),('1,630 exact-hash representatives','2,113 + 22 duplicates removed'),('1,619 QUAST-passing assemblies','7 fail + 4 review excluded'),('1,604 CheckM2-passing assemblies','15 excluded; ≥95% / ≤5%'),('1,595 GTDB-defined P. aeruginosa','9 GTDB P. paraeruginosa excluded')]
for i,(t,sub) in enumerate(labels):
 y=.92-i*.176;a.add_patch(FancyBboxPatch((.02,y-.098),.95,.126,boxstyle='round,pad=0.012',facecolor='#EDF4F6',edgecolor='#CEDDE3',linewidth=.65));a.text(.495,y-.01,t,ha='center',va='center',fontsize=8.4,fontweight='bold',color='#17344B');a.text(.495,y-.057,sub,ha='center',va='center',fontsize=7.5,color='#536470')
 if i<4:a.annotate('',xy=(.495,y-.151),xytext=(.495,y-.115),arrowprops={'arrowstyle':'->','color':GREY,'lw':.8})
a.text(.495,.02,'Historical source labels\nClinical 1,350  |  Environmental 99  |  Other 146',ha='center',va='center',fontsize=8.3)
panel(b,'b','Assembly size in the final collection');b.hist(q.total_length/1e6,bins=26,color=TEAL,edgecolor='white',linewidth=.25);b.axvline(q.total_length.median()/1e6,color='#17344B',ls='--');b.set(xlabel='Assembly length (Mb)',ylabel='Assemblies');b.text(.97,.96,'n = 1,595\nMedian 6.758 Mb',transform=b.transAxes,ha='right',va='top',fontsize=8)
panel(c,'c','Predicted genome quality');c.scatter(q.Contamination,q.Completeness,s=9,c=TEAL,alpha=.35,edgecolors='none');c.axvline(2,color=ORANGE,ls='--',lw=.8);c.axhline(98,color=ORANGE,ls='--',lw=.8);c.set(xlabel='CheckM2 contamination (%)',ylabel='CheckM2 completeness (%)',xlim=(-.05,4.25),ylim=(95,100.3));c.text(.98,.09,'1,585 stringent\n10 standard-only',transform=c.transAxes,ha='right',va='bottom',fontsize=8)
save(fig,'Figure_1_collection',[D/'figure1_final_QC.tsv',BASE/'workflow_catalogue.json'])
weighted=read(D/'occupancy_weighted_architecture.tsv');pp=read(D/'per_genome_frequency_classes.tsv');occ=read(BASE/'verified_occupancy_whole.tsv');print('occupancy columns',occ.columns.to_list())
freq=occ.iloc[:,1].to_numpy(int);hist=np.bincount(freq,minlength=1596);pd.DataFrame({'occupancy':np.arange(1,1596),'clusters':hist[1:]}).to_csv(D/'figure2_occupancy_spectrum.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(7.1,5.1),layout='constrained');gs=fig.add_gridspec(2,2);a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,:]);panel(a,'a','Gene-cluster occupancy');xx=np.arange(1,1596);a.scatter(xx[hist[1:]>0],hist[1:][hist[1:]>0],s=6,color=BLUE,alpha=.7);a.set(xscale='log',yscale='log',xlabel='Assemblies containing cluster',ylabel='Distinct clusters');a.text(.98,.98,'43,001 clusters\n1,595 assemblies',transform=a.transAxes,ha='right',va='top',fontsize=8)
panel(b,'b','Rare clusters (occupancy ≤5)');rr=[E['whole']['rare_le5_distinct_percent'],E['whole']['rare_le5_assignment_percent']];b.bar([0,1],rr,color=[GREY,TEAL],width=.55);b.set(xticks=[0,1],xticklabels=['Distinct\nclusters','Presence\nassignments'],ylabel='Share of total (%)',ylim=(0,63));
for i,v in enumerate(rr):b.text(i,v+1.8,f'{v:.2f}%',ha='center',fontsize=9)
panel(c,'c','Repertoire composition and pooled genome content');dd=weighted[weighted.population.eq('whole')];left=np.zeros(2)
for (_,r),co in zip(dd.iterrows(),COLORS):
 vals=np.array([r.distinct_percent,r.assignment_percent]);c.barh([1,0],vals,left=left,color=co,height=.5,label=r.category)
 for j,v in enumerate(vals):
  if v>6:c.text(left[j]+v/2,[1,0][j],f'{v:.2f}%',ha='center',va='center',color='white',fontsize=8.5)
 left+=vals
c.set(yticks=[1,0],yticklabels=['Distinct clusters\n(n = 43,001)','Genome–cluster presences\n(n = 9,700,435)'],xlabel='Fraction (%)',xlim=(0,100));c.legend(ncol=4,loc='upper center',bbox_to_anchor=(.5,-.38),frameon=False,fontsize=7.6);save(fig,'Figure_2_pangenome_scales',[D/'figure2_occupancy_spectrum.tsv',D/'occupancy_weighted_architecture.tsv'])
curves={pop:read(D/f'{pop}_exact_rarefaction.tsv') for pop in ['whole','clinical']};sen=read(D/'clinical_fluidity_sensitivity.tsv')
fig,axes=plt.subplots(2,2,figsize=(7.1,5.2),layout='constrained')
for pop,co in [('whole',BLUE),('clinical',TEAL)]:
 t=curves[pop];lab=f'{pop.capitalize()} (n = {E[pop]["n_genomes"]:,})';axes[0,0].plot(t.n_genomes,t.expected_pangenome/1000,c=co,label=lab);axes[0,0].plot(t.n_genomes,t.fitted_pangenome/1000,c=co,ls='--',lw=.8);axes[0,1].plot(t.n_genomes[1:],t.expected_new_clusters[1:],c=co,label=lab);axes[1,0].plot(t.n_genomes,t.expected_strict_core/1000,c=co,label=lab)
for a,l,t in zip(axes.flat,'abcd',['Exact pangenome accumulation','Expected new clusters per addition','Strict-core expectation','Clinical fluidity sensitivity']):panel(a,l,t)
axes[0,0].set(xlabel='Assemblies sampled',ylabel='Expected clusters (×1,000)');axes[0,0].legend(frameon=False,loc='lower right');axes[0,0].text(.04,.94,'γ = 0.290002 / 0.288462\nDashed: power-law fits',transform=axes[0,0].transAxes,va='top',fontsize=7.5)
axes[0,1].set(xlabel='Assemblies sampled',ylabel='Expected new clusters',yscale='log');axes[0,1].text(.98,.96,'Endpoint: 7.707 / 8.668',transform=axes[0,1].transAxes,ha='right',va='top',fontsize=8)
axes[1,0].set(xlabel='Assemblies sampled',ylabel='Expected universal clusters (×1,000)');axes[1,0].text(.98,.95,'Endpoints: 2,352 / 2,468',transform=axes[1,0].transAxes,ha='right',va='top',fontsize=8)
a=axes[1,1];a.plot(np.arange(4),sen['mean'],marker='o',color=TEAL);mean=sen['mean'].iloc[0];a.errorbar([0],[mean],yerr=[[mean-.1053297665],[.1093551649-mean]],fmt='none',ecolor=TEAL,capsize=3);a.set(xticks=range(4),xticklabels=['≥1','≥2','≥5','≥14'],xlabel='Minimum cluster occupancy',ylabel='Mean genomic fluidity (φ)',ylim=(.095,.113));a.text(.03,.10,'n = 1,350 assemblies\n910,575 pairs at each threshold',transform=a.transAxes,fontsize=7.6);save(fig,'Figure_3_accumulation_fluidity',[D/'whole_exact_rarefaction.tsv',D/'clinical_exact_rarefaction.tsv',D/'clinical_fluidity_sensitivity.tsv'])
bur=read(D/'standalone_amr_burden.tsv');prev=read(BASE/'verified_amr_element_prevalence.tsv').set_index('element');fig,axes=plt.subplots(2,2,figsize=(7.1,5.45),layout='constrained')
for a,l,t in zip(axes.flat,'abcd',['Total resistance-associated symbols','Gene and point-variant contributions','Selected exact element symbols','Selected beta-lactamase groups']):panel(a,l,t)
a=axes[0,0];a.hist(bur.unique_symbols,bins=np.arange(3.5,36.5,1),color=BLUE,edgecolor='white',lw=.3);a.set(xlabel='Distinct symbols per assembly',ylabel='Assemblies');a.text(.98,.96,'n = 1,350\nMedian 10; IQR 8–16',transform=a.transAxes,ha='right',va='top',fontsize=8)
a=axes[0,1]
for col,lab,co in [('nonpoint_symbols','Non-point gene symbols',TEAL),('point_symbols','Point/disruption symbols',ORANGE)]:
 v=np.sort(bur[col]);a.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=lab,c=co)
a.set(xlabel='Distinct symbols per assembly',ylabel='Cumulative assembly fraction');a.legend(frameon=False,loc='lower right',fontsize=7.5)
a=axes[1,0];els=['fosA',"aph(3')-IIb",'catB7','nalC_G71E','gyrA_T83I','sul1'];els=[s for s in els if s in prev.index];vals=[int(prev.loc[s,'genomes']) for s in els];a.hlines(range(len(els)),0,np.array(vals)/1350*100,color='#D6E2EC');a.scatter(np.array(vals)/1350*100,range(len(els)),c=BLUE,s=22);a.set(yticks=range(len(els)),yticklabels=els,xlabel='Assemblies with detection (%)',xlim=(0,125));a.invert_yaxis()
for i,v in enumerate(vals):a.text(v/1350*100+3,i,str(v),va='center',fontsize=7.8)
a=axes[1,1];groups=['blaVIM','blaIMP','blaKPC','blaNDM'];vs=[F['standalone_AMR']['gene_groups'][s]['genomes'] for s in groups];a.hlines(range(4),0,np.array(vs)/1350*100,color='#D6E2EC');a.scatter(np.array(vs)/1350*100,range(4),c=ORANGE,s=24);a.set(yticks=range(4),yticklabels=groups,xlabel='Assemblies with prefix-group hit (%)',xlim=(0,12));a.invert_yaxis()
for i,v in enumerate(vs):a.text(v/1350*100+.35,i,str(v),va='center',fontsize=8)
pd.DataFrame({'group':groups,'genomes':vs,'denominator':1350}).to_csv(D/'figure4_selected_beta_lactamase_groups.tsv',sep='\t',index=False);save(fig,'Figure_4_resistance_detections',[D/'standalone_amr_burden.tsv',BASE/'verified_amr_element_prevalence.tsv',D/'figure4_selected_beta_lactamase_groups.tsv'])
fb=read(D/'functional_burden.tsv');fig,axes=plt.subplots(2,2,figsize=(7.1,5.4),layout='constrained')
for a,l,t in zip(axes.flat,'abcd',['AMP prediction evidence','Predicted peptide occurrences','Carbohydrate-active domain classes','Domains and distinct annotated genes']):panel(a,l,t)
a=axes[0,0];v=[38699,1361,1];a.bar(range(3),v,color=[GREY,TEAL,ORANGE]);a.set(yscale='log',ylim=(.7,200000),xticks=range(3),xticklabels=['Single\npredictor','≥2\npredictors','Database\nsupport'],ylabel='Peptide occurrences (log scale)')
for i,z in enumerate(v):a.text(i,z*1.45,f'{z:,}',ha='center',fontsize=8.5)
a=axes[0,1];a.hist(fb.AMP_occurrences,bins=np.arange(17.5,53.5,1),color=TEAL,edgecolor='white',lw=.2);a.set(xlabel='Predicted occurrences per assembly',ylabel='Assemblies');a.text(.98,.96,'n = 1,350\nMedian 29; IQR 26–33',transform=a.transAxes,ha='right',va='top',fontsize=8)
a=axes[1,0];cl=F['CAZyme']['class_occurrences'];keys=['GT','GH','CBM','AA','CE','PL'];vs=[cl[k] for k in keys];a.barh(keys,np.array(vs)/1000,color=BLUE);a.invert_yaxis();a.set(xlabel='Domain/family annotations (×1,000)',xlim=(0,83))
for i,z in enumerate(vs):a.text(z/1000+1,i,f'{z:,}',va='center',fontsize=8)
a=axes[1,1];arr=[fb.CAZyme_domain_annotations.to_numpy(),fb.CAZyme_unique_genes.to_numpy()];vp=a.violinplot(arr,showmedians=True,showextrema=True)
for pc in vp['bodies']:pc.set_facecolor(TEAL);pc.set_alpha(.35)
for k in ['cmedians','cmins','cmaxes','cbars']:vp[k].set_color(TEAL);vp[k].set_linewidth(.8)
a.set(xticks=[1,2],xticklabels=['Domain\nannotations','Distinct\ngenes'],ylabel='Count per assembly');a.text(.98,.96,'n = 1,350\nMedians 112 / 101',transform=a.transAxes,ha='right',va='top',fontsize=8);save(fig,'Figure_5_functional_predictions',[D/'functional_burden.tsv',BASE/'verified_functional_results.json'])
# Supplementary figure 1: byte-identical reuse of a verified original plot.
original=BASE/'original_corrected_fluidity.png';target=G/'Supplementary_Figure_1_existing_corrected_fluidity.png';shutil.copyfile(original,target);assert hashlib.sha256(original.read_bytes()).hexdigest()==hashlib.sha256(target.read_bytes()).hexdigest();REG.append({'figure':'Supplementary_Figure_1_existing_corrected_fluidity','png':str(target.relative_to(OUT)),'original_source':'17.clinical_comparative_genomics/pangenome fluidity/fluidity_results_corrected/figure_01_pairwise_fluidity_distribution.png','reproduction_input':str(original.relative_to(OUT)),'reuse':'Byte-identical; no image changes','sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
t=read(BASE/'verified_candidates_all43001.tsv');labels=['xerC/xerD cluster','sul1 cluster','group_5721','group_29837'];# use original row labels to avoid order assumptions
labels=[('xerC/xerD cluster' if 'xer' in x else 'sul1 cluster' if 'sul1' in x else x) for x in t.gene_cluster]
fig,a=plt.subplots(figsize=(7.1,3.4),layout='constrained');yy=np.arange(len(t));a.hlines(yy,t.environmental_frequency*100,t.clinical_frequency*100,color='#C4D2D9',lw=2);a.scatter(t.clinical_frequency*100,yy,c=BLUE,label='Clinical-labelled (n = 1,350)',s=35);a.scatter(t.environmental_frequency*100,yy,c=ORANGE,label='Environmental-labelled (n = 99)',s=35);a.set(yticks=yy,yticklabels=labels,xlabel='Assemblies containing the cluster (%)',xlim=(0,47));a.invert_yaxis()
for i,row in enumerate(t.itertuples()):a.text(33,i,f'q = {row.BH_all43001:.4f}',va='center',fontsize=8.5)
a.legend(loc='lower center',bbox_to_anchor=(.5,-.42),ncol=2,frameon=False);a.set_title('Exploratory source-label associations; metadata and lineage confounding remain',fontsize=9.5,pad=15);save(fig,'Supplementary_Figure_2_source_associations',[BASE/'verified_candidates_all43001.tsv'])
(G/'figure_manifest.json').write_text(json.dumps(REG,indent=2));print('Created',len(REG),'figure groups')
