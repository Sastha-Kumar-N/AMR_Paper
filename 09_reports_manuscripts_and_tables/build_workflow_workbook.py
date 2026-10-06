from pathlib import Path
import sys,json,math,hashlib,csv,collections,datetime
import pandas as pd
O=Path(__file__).resolve().parent;ROOT=O.parent
sys.path.append(str(ROOT/'18.amrfinderplus/amrfinderplus_interpretation/envs/amrinterp/lib/python3.13/site-packages'))
from openpyxl import Workbook,load_workbook
from openpyxl.styles import Font,PatternFill,Alignment,Border,Side
from openpyxl.worksheet.table import Table,TableStyleInfo
from openpyxl.utils import get_column_letter
from openpyxl.chart import BarChart,Reference
from study_workflow_data import C,F,P,A,S,W,FLOW,CORRECTIONS,FINDINGS,NEXT
NAVY='17344B';TEAL='007F86';LIGHT='E9F3F3';INK='203447';GREY='5B6D78'
STATUS={'Verified':('DDEFE6','176344'),'Qualified':('FFF0CC','886000'),'Preliminary':('FFE4CA','985000'),'Exploratory':('EDE4F5','67448C'),'Invalid':('FADDDD','A22828'),'Incomplete':('ECEFF2','626C75'),'Absent':('ECEFF2','626C75'),'Partial':('FFF0CC','886000')}
wb=Workbook();wb.remove(wb.active);sheet_rows={};metrics=[]
def metric(analysis,name,value,population,unit,source,status='Verified',note=''):
 metrics.append([f'N{len(metrics)+1:03}',analysis,name,value,population,unit,status,S.get(source,source),note])
for r in FLOW:
 metric('Dataset flow',r[0]+' output',r[4],r[0],'assemblies/files','acquisition' if r[0]=='Preserved input collection' else r[5],'Qualified' if 'branch' in r[0] or 'coordinate' in r[0] else 'Verified',r[6])
for label,val in [('Final stringent genomes',1585),('Final standard-only genomes',10),('Stringent before taxonomy',1592)]:metric('QC',label,val,'See metric label','assemblies','checkm')
for stage,key,source in [('QUAST','quast_final_stats','quast'),('CheckM2','checkm_final_stats','checkm')]:
 for k,d in C['cohort'][key].items():
  for stat in ['min','q1','median','q3','max','mean']:
   if stat in d:metric(stage,k+' '+stat,d[stat],'Final 1,595','bp' if k in ['total_length','n50','largest_contig','Genome_Size'] else ('%' if k in ['gc','Completeness','Contamination'] else 'metric units'),source)
for pop,d in C['pangenomes'].items():
 if pop.startswith('phylogeny'):continue
 for key,val in d.items():
  if isinstance(val,(int,float)) and not isinstance(val,bool):metric('Pangenome '+pop,key,val,f'{d["n_genomes"]} assemblies','clusters' if key!='n_genomes' else 'assemblies',pop)
for pop in ['whole','clinical5000']:
 d=C['heaps_refits'][pop];sk='heaps_whole' if pop=='whole' else 'heaps_clinical'
 for k in ['kappa','gamma','R2','last_step_gain']:metric('Heaps '+pop,k,d[k],'1,595 whole' if pop=='whole' else '1,350 clinical','parameter / clusters',sk,'Qualified','Refit to stored mean; fitted values may differ in final floating-point digits.')
for k,v in [('gamma_CI_low',0.28941053882335727),('gamma_CI_high',0.29104398871826487),('new_gene_alpha',0.6542100598013056),('alpha_CI_low',0.6501624114829013),('alpha_CI_high',0.6581479190844287),('permutations',500),('order_bootstraps',200)]:metric('Whole Heaps',k,v,'1,595 whole','parameter / replicates','heaps_whole','Qualified','CI resamples orders, not independent biological cohorts.')
for d in C['fluidity_recomputed']:
 metric('Clinical fluidity',f'Mean phi retaining occupancy ≥{d["threshold"]}',d['mean'],'1,350 clinical','fraction','fluidity');metric('Clinical fluidity',f'Clusters retaining occupancy ≥{d["threshold"]}',d['clusters'],'1,350 clinical','clusters','fluidity')
for k,v in [('pairs',910575),('bootstrap_CI_low',0.1053297665),('bootstrap_CI_high',0.1093551649),('jackknife_CI_low',0.1053113723),('jackknife_CI_high',0.1094120465)]:metric('Clinical fluidity',k,v,'1,350 clinical','pairs / fraction','fluidity')
for k in ['shared','clinical_only','environmental_only','absent_both','shared_core99','shared_below99_both','shared_not_core_in_both']:metric('Common-matrix sharing',k,C['common_matrix'][k],'1,350 clinical + 99 environmental','clusters','comparison','Qualified')
for k in ['raw_files','raw_rows','raw_genomes','filtered_rows','unique_exact_elements','unique_nonpoint_symbols','unique_point_symbols','point_positive_genomes','nonpoint_positive_genomes','distinct_class_labels','single_genome_elements','q90_burden','n_ge_q90']:metric('Standalone AMR',k,F['standalone_AMR'][k],'1,350 clinical','defined by metric','AMR','Verified','Detections are genomic, not AST phenotypes.')
for k in ['burden_unique_elements','burden_nonpoint','burden_point']:
 for stat in ['min','q1','median','q3','max','mean']:metric('Standalone AMR',k+' '+stat,F['standalone_AMR'][k][stat],'1,350 clinical','distinct element symbols/genome','AMR_burden','Qualified')
for family,d in F['standalone_AMR']['gene_groups'].items():metric('Standalone AMR',family+' detected genomes',d['genomes'],'1,350 clinical','genomes','AMR','Qualified','Prefix groups can overlap; no phenotype implied.')
for sec,keys,source in [
('AMP',['raw_records','unique_sequences','recomputed_high_confidence','multi_predictor','db_supported','three_predictor'],'AMP'),
('FuncScan_ARG',['raw_rows','primary_rows','primary_families','primary_exact_names','singleton_families'],'ARG'),
('BGC',['raw_rows','consensus_loci','cross_tool_loci','types','domain_signature_families','unclassified_only_loci'],'BGC'),
('CAZyme',['overview_files','overview_rows','recommended_raw_genes','recommended_raw_family_tokens','primary_records','primary_family_count','class_count','distinct_primary_genes','distinct_gene_family_pairs','CGC_rows'],'CAZyme')]:
 for k in keys:metric(sec,k,F[sec][k],'1,350 clinical','defined by metric',source,'Invalid' if sec=='BGC' else 'Qualified','BGC export counts only; biological catalogue withheld.' if sec=='BGC' else 'Prediction units; occurrence is not necessarily unique gene.')
for k,v in F['CAZyme']['class_occurrences'].items():metric('CAZyme',k+' occurrences',v,'1,350 clinical','domain/family annotations','CAZyme','Qualified')
for k,v in P['alignment'].items():
 if isinstance(v,(int,float)):metric('Clinical SNP alignment',k,v,'1,350 clinical','defined by metric','snp')
for k,v in [('UFBoot final correlation',0.930),('UFBoot replicates',1000),('SH-aLRT replicates',1000),('RF ML vs consensus',242)]:metric('Phylogeny',k,v,'1,350 clinical','defined by metric','phylogeny','Preliminary','Non-convergence explicitly warned in IQ-TREE log.')
for k,v in [('Old coordinate genomes',603),('Revised coordinate genomes',601),('Revised total',1347),('Eligible biallelic sites',84053),('High-divergence hybrid threshold genomes',10),('Significant locations',0)]:metric('Geography',k,v,'Population stated in metric','genomes/sites/locations','geography','Qualified','No AMR or epidemiological hotspot inference.')
for k,v in [('CDS sum',9833296),('CDS median',6184),('CDS min',5358),('CDS max',7526)]:metric('Prokka',k,v,'1,595 final','CDS','prokka')
for pop in ['whole','clinical','environmental']:metric('Representative sequences',pop+' reference FASTA records',A[pop+'_reference_fasta_records'],'Population-specific Panaroo run','FASTA records',pop,'Qualified','Representative FASTA counts do not equal cluster totals.')
# Machine-readable evidence references reused by Markdown documents.
(O/'workflow_catalogue.json').write_text(json.dumps(W,indent=2))
(O/'quantitative_result_catalogue.json').write_text(json.dumps(metrics,indent=2))

def sheet(name,title,subtitle,headers,rows,widths,status_col=None,height=70):
 ws=wb.create_sheet(name);n=len(headers);ws.sheet_view.showGridLines=False
 ws.merge_cells(start_row=1,start_column=1,end_row=1,end_column=n);c=ws.cell(1,1,title);c.font=Font(name='Aptos Display',size=20,bold=True,color='FFFFFF');c.fill=PatternFill('solid',fgColor=NAVY);c.alignment=Alignment(vertical='center');ws.row_dimensions[1].height=38
 for col in range(1,n+1):ws.cell(1,col).fill=PatternFill('solid',fgColor=NAVY)
 ws.merge_cells(start_row=2,start_column=1,end_row=2,end_column=n);ws.cell(2,1,subtitle).font=Font(name='Aptos',size=11,color=GREY);ws.cell(2,1).alignment=Alignment(wrap_text=True,vertical='center');ws.row_dimensions[2].height=36
 for j,h in enumerate(headers,1):c=ws.cell(4,j,h);c.font=Font(name='Aptos',bold=True,color='FFFFFF',size=11);c.fill=PatternFill('solid',fgColor=TEAL);c.alignment=Alignment(wrap_text=True,vertical='center')
 ws.row_dimensions[4].height=32
 for i,row in enumerate(rows,5):
  for j,v in enumerate(row,1):
   if isinstance(v,float) and not math.isfinite(v):v=None
   c=ws.cell(i,j,v);c.font=Font(name='Aptos',size=11,color=INK);c.alignment=Alignment(vertical='top',wrap_text=True);c.border=Border(bottom=Side(style='hair',color='DEE6EB'))
   if isinstance(v,float):c.number_format='0.000000' if abs(v)<2 else '#,##0.00'
   if isinstance(v,int):c.number_format='#,##0'
   if i%2==0:c.fill=PatternFill('solid',fgColor='F4F8FA')
  ws.row_dimensions[i].height=height
  if status_col:
   st=ws.cell(i,status_col).value
   if st in STATUS:
    bg,fg=STATUS[st];ws.cell(i,status_col).fill=PatternFill('solid',fgColor=bg);ws.cell(i,status_col).font=Font(name='Aptos',bold=True,color=fg,size=11)
 for j,w in enumerate(widths,1):ws.column_dimensions[get_column_letter(j)].width=w
 ws.freeze_panes='D5' if name=='Workflow' else 'B5';ws.sheet_properties.pageSetUpPr.fitToPage=True;ws.page_setup.orientation='landscape';ws.page_setup.paperSize=ws.PAPERSIZE_A3;ws.page_setup.fitToWidth=1;ws.page_setup.fitToHeight=0;ws.print_title_rows='1:4';ws.print_options.horizontalCentered=False;ws.oddFooter.center.text='P. aeruginosa | evidence reconstruction | 10 Sep 2026';ws.oddFooter.right.text='Page &P of &N'
 if rows:
  tab=Table(displayName='T'+''.join(ch for ch in name if ch.isalnum()),ref=f'A4:{get_column_letter(n)}{len(rows)+4}');tab.tableStyleInfo=TableStyleInfo(name='TableStyleMedium2',showFirstColumn=False,showLastColumn=False,showRowStripes=False,showColumnStripes=False);ws.add_table(tab)
 ws.sheet_properties.tabColor=TEAL;sheet_rows[name]=len(rows)
 return ws

readme=[
('Purpose','End-to-end study workflow with what was performed, verified results, interpreted outputs, distinctive findings, limitations and exact result sources.'),
('Start here','Open Workflow for all 39 analysis stages; Key Findings is the concise biological synthesis.'),
('Scope','This is a read-only reconstruction of completed project outputs, plus clearly identified numerical validation. No new expensive biological pipelines were run.'),
('Population rule','1,595 is the historical whole assembly collection; 1,350 clinical-labelled, 99 environmental-labelled and 146 other. The later geographic sensitivity uses 1,347, with 601 coordinates.'),
('Interpretation rule','Gene clusters, CDS, exact AMR symbols, gene families, peptide occurrences, BGC export groups and CAZyme domain annotations are different counting units.'),
('Unique finding','The Key Findings sheet describes distinctive observations in this dataset. Worldwide novelty has not been established by these files.'),
('Verified','Numbers independently reproduce underlying outputs or explicit recorded calculations. This status does not imply experimental validation.'),
('Qualified','Result is usable with the stated sample, method and biological limitations.'),
('Preliminary / Exploratory','A result or association exists, but unresolved confounding, source metadata or method issues prevent a strong conclusion.'),
('Invalid','Do not use for a biological claim. A demonstrated error or unresolved source-export defect affects the stated output; see Corrections.'),
('Incomplete / Absent / Partial','Required outputs are missing, retrieval failed, or only part of the analysis was completed. No results have been invented.'),
('Important corrections','Later habitat comparison: 178 mislabelled samples. Whole accessory distances: integer overflow. Whole cooccurrence: Boolean intersections. Phylogeny: non-converged UFBoot. BGC: region/CDS export inconsistency.'),
('Primary source','Project files take precedence over manuscripts. Source Registry records exact relative paths; Sources are relative to the project root.'),
('Preservation','Original project files were not altered. This workbook and audit scripts are under manuscript_reconstruction/.'),
('Preparation','10 September 2026. Prepared from the verified core, functional, phylogeny/geography and additional audit JSON/TSV files.'),
('Publication readiness','Core descriptive results are reproducible. Habitat inference, definitive phylogenetic support and BGC biology require the corrections listed in Next Analyses.'),
]
sheet('Read Me','P. aeruginosa | study workflow & evidence','An audit-ready workbook: observed results → interpretation → supported finding. Use the sheet tabs below.', ['Topic','Guide'],readme,[29,135],height=45)
wfheaders=['Step','Analysis stage','Experiment / analysis','What we have done','Input / population','Verified or audited result','Interpreted output','Distinctive finding','Evidence status','Limitations / publication boundary','Source key','Method / parameters']
sheet('Workflow','End-to-end experimental and computational workflow','39 stages. Computational analyses are not wet-lab experiments. Invalid and incomplete branches remain visible so they cannot enter the manuscript unnoticed.',wfheaders,[[w[k] for k in ['step','stage','experiment','what_done','population','result','interpreted_output','finding','status','limitation','source_key','method']] for w in W],[9,23,32,49,33,64,52,52,18,64,19,56],9,110)
sheet('Dataset Flow','Genome collection and analysis populations','Branches are separate populations, not sequential removals. “Outside branch” must not be read as a QC exclusion.', ['Stage','Input N','Filter / selection','Removed or outside branch','Output N','Source key','Population note'],FLOW,[36,14,55,23,14,20,68],height=52)
sheet('Key Findings','Distinctive observations supported by this collection','Dataset-specific findings, not claims of worldwide novelty. Each interpretation has a deliberate evidence boundary.', ['Finding ID','Finding','Quantitative support','Workflow steps','Evidence strength','Interpretation','Permitted wording','What the finding does not establish'],FINDINGS,[13,48,64,19,32,62,68,68],height=100)
sheet('Quantitative Results','Numerical reference catalogue','Numeric values preserve their underlying population and unit. BGC figures are retained solely as flagged export counts; do not publish them as a validated catalogue.', ['Metric ID','Analysis','Metric','Value','Population','Unit','Status','Exact source path(s)','Note'],metrics,[13,26,46,21,33,35,18,85,70],7,48)
assoc=pd.read_csv(O/'verified_candidates_all43001.tsv',sep='\t');print('association columns',assoc.columns.to_list())
# Preserve the audited association schema and add transparent population labels.
aheaders=list(assoc.columns);arows=assoc.astype(object).where(assoc.notna(),None).values.tolist()
aw=sheet('Clinical vs Env','Correctly matched source-label associations','Four BH-significant clusters from 43,001 tests. Historical labels: 1,350 clinical and 99 environmental. Metadata errors and phylogenetic confounding remain; this is preliminary.',aheaders,arows,[90]+[22]*(len(aheaders)-1),height=115)
for row in aw.iter_rows(min_row=5):
 for cell in row:
  if isinstance(cell.value,float) and abs(cell.value)<.001:cell.number_format='0.000E+00'
prev=pd.read_csv(O/'verified_amr_element_prevalence.tsv',sep='\t');amrows=[[x.element,int(x.genomes),int(x.genomes)/1350,'Standalone AMRFinderPlus 4.2.7; DB 2026-05-15.1','Exact element symbol includes point/disruption variants; not a phenotype'] for x in prev.itertuples()]
am=sheet('AMR Prevalence','Resistance-associated element prevalence','All 929 exact element symbols, including point/disruption labels. Denominator 1,350 assemblies; absence of AST prevents phenotypic claims.', ['Exact element symbol','Genomes detected','Prevalence','Analysis definition','Interpretation boundary'],amrows,[45,22,20,58,77],height=31)
for row in am.iter_rows(min_row=5):row[2].number_format='0.00%'
sheet('Corrections','Contradictions, errors and superseded results','Use the corrected/scoped interpretation; do not select numbers according to the most recent filename or old manuscript prose.', ['Issue ID','Analysis','Earlier / problematic claim','Verified or defensible interpretation','Why values differ','Required disposition','Source key'],CORRECTIONS,[13,25,70,77,70,73,20],height=105)
sheet('Next Analyses','Publication-critical work and future extensions','Priority is conditional on the claims retained. No expensive rerun has been launched automatically.', ['Priority','Module','Analysis / action needed','Why it is needed','Source key'],NEXT,[35,34,95,85,23],height=76)
source_rows=[]
for key,paths in S.items():
 for p in paths.split('; '):
  pp=ROOT/p;source_rows.append([key,p,'Directory/group' if pp.is_dir() else 'File','Present' if pp.exists() else 'Not found at exact path','; '.join(w['step'] for w in W if w['source_key']==key)])
sr=sheet('Source Registry','Exact sources and workflow cross-references','Paths are relative to the project root. File cells link to the preserved originals or new audit outputs. Directories are groups enumerated in source_file_manifest.tsv.', ['Source key','Exact relative path','Role / type','Path check','Workflow steps'],source_rows,[22,125,24,28,40],height=48)
for r in range(5,sr.max_row+1):
 c=sr.cell(r,2);c.hyperlink='../'+str(c.value);c.font=Font(name='Aptos',size=11,color='007F86',underline='single')
# Navigate source keys from Workflow to their first registry occurrence.
idx={}
for r in range(5,sr.max_row+1):idx.setdefault(sr.cell(r,1).value,r)
for r in range(5,wb['Workflow'].max_row+1):
 c=wb['Workflow'].cell(r,11);c.hyperlink=f"#'Source Registry'!A{idx[c.value]}";c.font=Font(name='Aptos',size=11,color=TEAL,underline='single')
# Compact QC and metadata review tabs support denominator and source-label assessment.
q=pd.read_csv(O/'verified_cohort_qc.tsv',sep='\t');keep=[x for x in ['assembly_accession','source_identifier','isolate_class','Completeness','Contamination','total_length','contigs','gc','n50'] if x in q.columns]
if len(keep)<5:keep=q.columns[:12].to_list()
sheet('Cohort QC','Assembly-level cohort and QC crosswalk','Derived audit table, not a new QC run. The complete verified_cohort_qc.tsv contains additional original metadata/QC columns.',keep,q[keep].astype(object).where(q[keep].notna(),None).values.tolist(),[28]*len(keep),height=32)
m=pd.read_csv(O/'metadata_manual_review_priority.tsv',sep='\t').fillna('')
sheet('Metadata Review','Priority examples for manual source recuration','This is a targeted review list, not an exhaustive list of misclassifications. Conflicting keyword hits alone do not prove a record is incorrect.',m.columns.to_list(),m.values.tolist(),[26,38,100,24,100],height=78)
# A directory is an analysis group, never an executable link. Scientific strings are plain values.
wb.active=0
name='P_aeruginosa_End_to_End_Workflow.xlsx';dest=O/name
if dest.exists():dest=O/f'P_aeruginosa_End_to_End_Workflow_{datetime.datetime.now():%Y%m%d_%H%M%S}.xlsx'
wb.save(dest)
check=load_workbook(dest,data_only=False)
errors=[]
for sn,rows in sheet_rows.items():
 if check[sn].max_row!=rows+4:errors.append(f'{sn}: unexpected row count')
 if check[sn].freeze_panes is None:errors.append(f'{sn}: missing freeze panes')
 if rows and len(check[sn].tables)!=1:errors.append(f'{sn}: missing filter table')
for ws in check:
 for row in ws:
  for c in row:
   if c.data_type=='f':errors.append(f'Unexpected formula: {ws.title}!{c.coordinate}')
assert not errors,errors
path_errors=[r for r in source_rows if r[3]!='Present']
assert not path_errors,path_errors
report={'workbook':str(dest.relative_to(ROOT)),'sheets':sheet_rows,'source_paths_checked':len(source_rows),'source_paths_missing':len(path_errors),'workflow_rows':len(W),'numeric_reference_rows':len(metrics),'amr_symbols':len(amrows),'validation':'PASS: saved and reopened; row counts, filters, freeze panes, numeric cells, no formulas and source existence checked','sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}
(O/'workflow_workbook_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
