#!/usr/bin/env python3
"""Write plain-text descriptions and legends for all manuscript figure images.
Inputs: registry, configs, existing legends, metrics, provenance and DOCX extracts.
Tools: Python standard library; date: 2026-10-05. No artwork or source data edited.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

S = Path(__file__).resolve().parents[1]
ROOT = S.parent
O = ROOT / 'Main_Pub_Figures'
sys.dont_write_bytecode = True
sys.path.insert(0, str(S / 'qc'))
from build_report import ADDENDA, CAPTIONS

# These descriptions explain the actual plotting functions. Scientific counts,
# thresholds and source labels remain in the source-derived legends/provenance.
PANELS = {
 'quast': {
  'a': 'Assembly-length distribution. The horizontal axis is total assembly length in megabases (Mb); the vertical axis counts assemblies. Blue bins contain the observed values, and the black dashed guide marks the median.',
  'b': 'N50 distribution. The horizontal axis is N50 in kilobases (kb); the vertical axis counts assemblies. N50 is a contiguity statistic: at least half the assembled length lies in contigs of this length or longer. The black dashed guide marks the median.',
  'c': 'Contig-count distribution. The horizontal axis counts contigs at or above the displayed minimum contig length; the vertical axis counts assemblies. The dashed guide marks the median.',
  'd': 'Assembly GC-content distribution. The horizontal axis is GC percentage; the vertical axis counts assemblies. The dashed guide marks the median.',
  'e': 'Largest-contig distribution. The horizontal axis gives each assembly’s longest contig in Mb; the vertical axis counts assemblies. The dashed guide marks the median.',
  'f': 'Relationship between contig count and N50. Each blue point represents an assembly at its observed values. Both axes use logarithmic scales. The printed Spearman coefficient summarizes rank association; it does not establish causation or a statistical significance test.'},
 'gc': {'a': 'Each translucent blue curve is one assembly’s retained GC-window distribution. The horizontal axis is GC percentage within a window; the vertical axis is the number of non-overlapping windows in that assembly. Overlapping curves show similar base-composition profiles; curve opacity is a display choice, not an additional quantitative measurement.'},
 'checkm': {
  'a': 'Distribution of model-predicted completeness among operationally passing assemblies. The horizontal axis is CheckM2 completeness (%), and the vertical axis counts assemblies.',
  'b': 'Distribution of model-predicted contamination among operationally passing assemblies. The horizontal axis is CheckM2 contamination (%), and the vertical axis counts assemblies.',
  'c': 'Model-predicted contamination versus completeness for the plotted operational-pass cohort. Filled blue points are stringent passes; open orange points are other operational passes. Dashed guides show the operational thresholds. Coincident points can overlap.',
  'd': 'Counts within the plotted operational-pass cohort. Blue denotes stringent passes, yellow-orange denotes other operational passes, and the outside-operational category is zero because excluded assemblies are not in this plotted subset. This zero must not be interpreted as zero exclusions from the evaluated cohort.'},
 'composition': {'a': 'Bars count Panaroo homologous gene clusters in the whole collection. Dark blue denotes core, light blue soft core, yellow-orange shell, and grey cloud. Bar labels are the exact category counts; the vertical axis is linear.'},
 'accumulation': {'a': 'The horizontal axis is the number of genomes added. The dark-grey curve is the mean cumulative number of distinct clusters; the blue curve is the mean number present in every genome of the current subset. Shaded envelopes show the saved percentile intervals across genome-order permutations. The curves summarize retained results rather than extrapolating beyond the analysed collection.'},
 'content': {'a': 'The vertical axis is the number of detected Panaroo clusters per assembly. Blue points are actual assembly observations at a common horizontal position. The box spans the first to third quartiles; its horizontal line marks the median. Whiskers reach the most extreme observed values within 1.5 interquartile ranges of the box. Points beyond the whiskers remain visible; horizontal position has no biological meaning.'},
 'rare': {
  'a': 'Horizontal bars count singleton clusters, doubletons, and clusters below the displayed occupancy thresholds. These definitions overlap: singleton and doubleton clusters also belong to the broader low-occupancy groups. Therefore these four bar counts must not be added as disjoint categories.',
  'b': 'Stacked bars show disjoint partitions using two different denominators: all clusters and cloud clusters only. Dark blue denotes singletons, yellow-orange doubletons, grey other cloud clusters, and green non-cloud clusters. Printed percentages use the denominator named beside each bar. Non-cloud clusters do not contribute to the cloud-only bar.'},
 'network': {
  'a': 'Bars count clinical Panaroo clusters in core, soft-core, shell and cloud categories. The vertical axis is logarithmic. Colours follow the shared frequency-category palette.',
  'b': 'The core/soft-core backbone is drawn using saved node coordinates. Nodes are homologous gene clusters and displayed edges are retained gene-adjacency relationships after the recorded support filter. Dark and light blue distinguish core and soft-core nodes. Spatial proximity in the layout is not a measured genomic distance.',
  'c': 'The full clinical node set is shown with edges passing the visualization support threshold. Dark blue, light blue, yellow-orange and grey denote core, soft-core, shell and cloud nodes. The caption under the graph reports displayed nodes and edges. Graph layout does not demonstrate transmission, mobile-element identity or functional interaction.',
  'd': 'Horizontal bars identify the five highest-degree clusters in the complete retained adjacency graph. Degree is the number of adjacent nodes. Italic names are annotation-inferred labels. This source-derived ranking differs from the manuscript prose and requires correction before figure substitution.',
  'e': 'The recorded radius-2 mexA neighbourhood shows nodes connected by gene-adjacency edges. Node outlines encode frequency category. Open nodes and italic gene labels distinguish annotation-inferred identities from experimentally demonstrated function or database-confirmed resistance. Label leader lines are guides, not additional graph edges.',
  'f': 'Connected components are ranked from largest to smallest. Each point represents one component; its vertical coordinate is the number of gene-cluster nodes on a logarithmic scale. A large component with many small components describes graph connectivity, not causal biological modules.'},
 'convergence': {'a': 'Pink points show mean pairwise genomic fluidity for each retained subset size. Error bars are the saved 2.5th–97.5th percentile bounds across replicate subset estimates. The dashed horizontal line marks the complete-cohort mean. These intervals describe the resampling distribution; they are not standard deviations or measured phenotype uncertainty.'},
 'fluidity': {'a': 'Pink histogram bars count genome pairs within the retained fluidity bins. The black dashed line marks the mean and the blue dotted line the median. A larger fluidity value denotes a greater fraction of summed gene-family occurrences unique to one member of a pair. It is not a nucleotide mutation rate.'},
 'map': {'a': 'Each bubble is a recorded coordinate-defined sampling location. Bubble area represents the number of assemblies, and its cividis fill represents mean core-genome differences from the cohort consensus. Blue outlines denote conserved-dominant locations, yellow-orange mixed locations, and orange candidate high-divergence clusters. The colour scale is clipped at the stated source-derived percentiles. Geographic gaps indicate missing sampling or coordinates, not biological absence.'},
 'country': {'a': 'Blue bars use the left axis to show mapped assembly counts for the most sampled countries. The orange line and points use the right axis for mean core-genome differences from the cohort consensus. Countries are ordered by sample count. Connections between countries are visual guides, not geographic trajectories or evolutionary paths; unequal sample sizes constrain comparisons.'},
 'tree': {'a': 'Grey branches show the native Newick topology, with lengths on the substitution-per-site scale. The outer ring provides one marker per tip: orange marks tips in the upper tail of terminal-branch length and blue denotes the others. Ring radius is a layout device, not a branch-length measurement. The labelled accession has the longest terminal branch; terminal length and total root-to-tip distance are distinct. The ring does not encode resistance, virulence or verified lineage identities.'},
 'comparison': {'a': 'Grouped bars compare frequency-category counts from the separate whole, clinical and environmental Panaroo runs. Dark blue, light blue, yellow-orange and grey denote core, soft core, shell and cloud. The vertical axis is logarithmic, and cohort genome counts appear beneath the groups. Frequency-category counts do not establish one-to-one cluster identity across independent runs.'},
 'totals': {'a': 'Bars show total Panaroo cluster counts for the whole, clinical and environmental collections. Blue, pink and green identify the respective cohorts. Direct labels report the source totals. Different cohort sizes prevent treating raw total differences as evidence of enrichment.'},
 'shared': {'a': 'Horizontal bars show retained comparison-category counts in the shared framework. Grey denotes rare cloud, blue conserved in both core, yellow-orange shared accessory/variable, green environmental soft-core but not clinical, and pink clinical soft-core but not environmental. Direct labels preserve the counts even where a small bar is difficult to see. These categories describe occupancy and are not significance-test results.'}
}

CAVEATS = {
 'gc': ['GC_SUBSET_UNSPECIFIED'], 'accumulation': ['ACCUMULATION_PERMUTATIONS'],
 'checkm': ['CHECKM_PLOTTED_DENOMINATOR'], 'network': ['NETWORK_RANKING','IMAGE_RELATIONSHIP_MISSING'],
}
FORMATS = {'.pdf', '.svg', '.tiff', '.png'}

def read(path):
 return json.loads(path.read_text(encoding='utf-8'))

def bullet_files(paths, base):
 return ['- '+str(p.relative_to(base)) for p in paths]

def main():
 registry = read(S/'figure_registry.json')
 paragraphs = {p['index']:p['text'] for p in read(O/'manuscript_extract/paragraphs.json')}
 issues = {x['id']:x for x in read(O/'NUMBER_QC_summary.json')['discrepancies']}
 records = []
 covered = []
 for row in registry:
  key = row['figure_id']; conf = read(S/key/'config.yaml'); folder = O/row['output_folder']; kind = conf['kind']
  provenance = read(folder/'provenance.json')
  lines = [key, '='*len(key), '', 'IMAGE DESCRIPTION AND LEGENDS', '',
    'Manuscript identifier: '+row['literal_number'], 'Analysis: '+row['analysis_type'],
    'Recovery status: '+('D — UNRECOVERABLE' if conf['classification']=='D' else 'A — regenerated from retained source data'), '',
    'This guide describes the files currently in this folder. Newly written explanatory legends are distinguished below from verbatim manuscript captions. Scientific corrections are flagged; the manuscript itself is unchanged.', '']
  final_files = sorted(p for p in (folder/'final').iterdir() if p.suffix in FORMATS)
  sources = ['SOURCE DATA AND REPRODUCIBILITY', '']
  for source in conf['inputs'].values():
   sources += ['- '+source['original_path'], '  SHA-256: '+source['sha256']]
  sources += ['', 'Canonical script: '+row['script_used'], 'Complete input and plotting provenance: provenance.json',
    'Element-to-source/column mapping: PLOTTED_ELEMENT_SOURCES.tsv',
    'Exact plotted observations and aggregates: panels/*.tsv', '']
  if conf['classification']=='D':
   sources[-3:-1] = ['Candidate files are evidence of the search only; no plotted data or final image is claimed.']
   lines += ['CURRENT IMAGE STATUS', '', ADDENDA[kind], '', 'The intended cohort panel breakdown is unspecified.' if key=='FigX_functional_cohort' else 'Intended panels a–c are described in the manuscript caption below; none was recovered.', '', 'MISSING INPUTS', '']
   lines += ['- '+x for x in conf['missing_inputs']]+['']
  else:
   assert set(PANELS[kind]) == set(conf['panels']), key
   legend = ['LEGEND FOR THE CURRENT REGENERATED FIGURE', '', ADDENDA[kind], '']
   for panel in conf['panels']:
    legend += ['Panel '+panel+'. '+PANELS[kind][panel], '']
   lines += ['FINAL IMAGE FILES', '', 'The formats below contain the same figure; this legend applies to all of them.', '']+bullet_files(final_files,folder)+['']+legend
   sidecar = folder/'final'/(key+'_legend.txt')
   sidecar.write_text('\n'.join([key,'','APPLIES TO:',*bullet_files(final_files,folder),'',*legend,*sources]),encoding='utf-8')
   records.append(sidecar)
   for p in final_files:covered.append({'image':str(p.relative_to(ROOT)), 'legend':str(sidecar.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
   # Standalone panel images have their own applicable legends beside them.
   for panel in conf['panels']:
    stem=key+'_panel_'+panel
    files=sorted(p for p in (folder/'panels').glob(stem+'*') if p.suffix in FORMATS)
    if not files:continue
    lines += ['STANDALONE PANEL '+panel, '']+bullet_files(files,folder)+['',PANELS[kind][panel],'']
    sidecar=folder/'panels'/(stem+'_legend.txt')
    ledger=[x for x in provenance['plotted_elements'] if x['panel']==panel]
    text=[key+' — panel '+panel,'','APPLIES TO:',*bullet_files(files,folder),'','PANEL LEGEND','',PANELS[kind][panel],'','FIGURE CONTEXT AND SCIENTIFIC CAVEATS','',ADDENDA[kind],'','PANEL DATA SOURCES','']
    for x in ledger:text += ['- '+x['source_file'], '  Columns: '+str(x['source_columns']), '  Exact plotted data: '+x['plot_data']]
    sidecar.write_text('\n'.join(text)+'\n',encoding='utf-8');records.append(sidecar)
    for p in files:covered.append({'image':str(p.relative_to(ROOT)),'legend':str(sidecar.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
  # Preserve original captions verbatim, rather than presenting corrections as quotes.
  caption_ids=CAPTIONS.get(kind,[])
  if key=='FigX_functional_cohort':caption_ids=[61]
  if key=='FigX_ARG_loci':caption_ids=[64,65]
  lines += ['ORIGINAL MANUSCRIPT CAPTION (VERBATIM; NOT AUTOMATICALLY CORRECTED)', '']
  if caption_ids:
   for i in caption_ids:
    caption=paragraphs[i]
    if i==64:caption=caption[caption.index('Figure X |'):]
    lines += ['Paragraph '+str(i)+':',caption,'']
  else:lines += ['No explicit, uniquely attributable caption was recovered beside this image. The current-artwork legend above is a source-grounded explanation, not a quotation from the manuscript.','']
  lines += ['SCIENTIFIC CAVEATS', '']
  for code in CAVEATS.get(kind,[]):
   issue=issues[code];lines += ['- '+issue['issue'], '  Action: '+issue['action']]
  lines += ['- Figure numbers repeat in the manuscript; the folder suffix identifies the analysis unambiguously.', '- Gene annotations, graph adjacency and genomic divergence do not by themselves establish measured resistance, virulence, transmission or causation.', '']
  original_files=sorted(p for p in (folder/'original_docx_image').iterdir() if p.suffix in FORMATS)
  reference=['ORIGINAL DOCX IMAGE — REFERENCE ONLY','']
  if not original_files:reference += ['The DOCX relationship references '+row['docx_image_file']+', but that archive member was absent. No substitute raster was invented.']
  else:
   reference += ['Files:',*bullet_files(original_files,folder),'']
   if key=='FigX_functional_cohort':reference += ['This embedded image actually shows fluidity convergence. It does not depict the cohort-composition analysis named in its surrounding caption. No cohort legend can accurately be assigned to this raster. See the regenerated FigS7_fluidity_convergence folder for the corresponding analysis.']
   elif key=='FigX_ARG_loci':reference += ['This embedded image actually shows a pairwise-fluidity distribution. It does not depict reconciled ARG loci. No ARG legend can accurately be assigned to this raster. See the regenerated FigUnnumbered_fluidity folder for the corresponding analysis.']
   else:reference += ['This is the unmodified low-resolution manuscript reference for '+row['analysis_type']+'. See the current-artwork legend and panel explanations in ../IMAGE_DESCRIPTIONS_AND_LEGENDS.txt. Exact agreement with unreadable raster labels is not asserted; source data and documented caveats govern the reconstruction.']
  reference += ['', 'Do not submit this reference raster as the regenerated final artwork.']
  sidecar=folder/'original_docx_image/REFERENCE_IMAGE_DESCRIPTION.txt';sidecar.write_text('\n'.join(reference)+'\n',encoding='utf-8');records.append(sidecar)
  lines += reference+['']+sources
  guide=folder/'IMAGE_DESCRIPTIONS_AND_LEGENDS.txt';guide.write_text('\n'.join(lines)+'\n',encoding='utf-8');records.append(guide)
  for p in original_files:covered.append({'image':str(p.relative_to(ROOT)),'legend':str(sidecar.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 with (O/'IMAGE_LEGEND_INDEX.tsv').open('w',encoding='utf-8',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=['image','legend','sha256'],delimiter='\t');writer.writeheader();writer.writerows(covered)
 # Coverage is checked against actual image files, not assumed format counts.
 actual={str(p.relative_to(ROOT)) for row in registry for part in ['final','panels','original_docx_image'] for p in (O/row['output_folder']/part).iterdir() if p.suffix in FORMATS}
 assert {r['image'] for r in covered}==actual
 for r in covered:assert (ROOT/r['legend']).is_file() and (ROOT/r['legend']).stat().st_size
 result={'figure_folders':len(registry),'text_files_written':len(records),'images_covered':len(covered),'artwork_modified':False,'files':[str(p.relative_to(ROOT)) for p in records]}
 (O/'IMAGE_GUIDES_QC.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(f"Image guides: {len(registry)} figure folders; {len(records)} text files; {len(covered)} images covered.")

if __name__=='__main__':main()
