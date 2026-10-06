"""Build a source-bounded bibliography registry and audit of the supplied drafts.

This script imports no original project code and changes no original files.
The curated corrections below were checked against the primary sources retained
in evidence_extracts/agent_reference_*.txt and reference_verification_*.txt.
"""
from pathlib import Path
import ast, csv, json, re

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'publication_manuscript'
TAB = OUT / 'reconstructed_tables'
DATE = '2026-10-01'
prior = ROOT / 'manuscript_draft/manuscript_reconstruction/article_20260911/scripts/literature.py'
tree = ast.parse(prior.read_text())
old = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
           and any(isinstance(t, ast.Name) and t.id == 'REFERENCES' for t in n.targets))
registry = {}
for key, (authors, title, details, doi, url) in old.items():
    registry[key] = dict(reference_id=key, authors=authors, title=title,
                         bibliographic_details=details, doi=doi, url=url,
                         status='VERIFIED', source_type='PRIMARY_RESEARCH',
                         checked_date=DATE, claim_scope='', limitations='',
                         recommended_use='METHODS_OR_CONTEXT', evidence_file='')

def update(key, **kw):
    registry[key].update(kw)

def add(key, authors, title, details, doi, url, scope, kind='PRIMARY_RESEARCH', use='CONTEXT', limitation=''):
    registry[key] = dict(reference_id=key, authors=authors, title=title,
        bibliographic_details=details, doi=doi, url=url, status='VERIFIED',
        source_type=kind, checked_date=DATE, claim_scope=scope,
        limitations=limitation, recommended_use=use, evidence_file='')

update('quast', title='QUAST: quality assessment tool for genome assemblies.')
update('funcscan', title='nf-core/funcscan v4.0.0: (Meta-)genome screening for functional and natural product gene sequences.',
       bibliographic_details='Official software release, 27 June 2026; commit aee3dc9. Accessed 1 October 2026.',
       url='https://github.com/nf-core/funcscan/releases/tag/4.0.0', source_type='OFFICIAL_SOFTWARE_DOCUMENTATION',
       limitations='The official 4.0.0 release and commit aee3dc9 match saved project manifests v4.0.0-gaee3dc9. A workflow citation does not establish that every optional module ran. Do not cite the unverified Fernandes preprint.')
update('bh', url='https://doi.org/10.1111/j.2517-6161.1995.tb02031.x', source_type='PRIMARY_METHODS')
update('ampir', bibliographic_details='Bioinformatics 36(21), 5262–5263 (2020).',
       limitations='Publisher issue and online-publication metadata give 2020; PubMed display has a January 2021 date. Use publisher year 2020, volume 36, pages 5262–5263, DOI btaa653.')
update('penaranda', bibliographic_details='mSystems 11(3), e01362-25 (2026).',
       url='https://journals.asm.org/doi/10.1128/msystems.01362-25')
add('who', 'World Health Organization.',
    'WHO bacterial priority pathogens list, 2024: bacterial pathogens of public health importance to guide research, development and strategies to prevent and control antimicrobial resistance.',
    'WHO, Geneva (2024). ISBN 9789240093461.', '',
    'https://www.who.int/publications/i/item/9789240093461',
    'The 2024 list places carbapenem-resistant P. aeruginosa in the HIGH-priority group.',
    'OFFICIAL_PUBLIC_HEALTH_REPORT', limitation='Priority applies to the specified resistance phenotype; it is not a susceptibility result for this genome collection.')
add('sati', 'Sati, H. et al.',
    'The WHO Bacterial Priority Pathogens List 2024: a prioritisation study to guide research, development, and public health strategies against antimicrobial resistance.',
    'Lancet Infect. Dis. 25(9), 1033–1043 (2025).', '10.1016/S1473-3099(25)00118-5',
    'https://pubmed.ncbi.nlm.nih.gov/40245910/',
    'Evidence and expert-input basis for the 2024 WHO prioritisation; high priority of CRPA.')
add('mosquera', 'Mosquera-Rendón, J. et al.',
    'Pangenome-wide and molecular evolution analyses of the Pseudomonas aeruginosa species.',
    'BMC Genomics 17, 45 (2016).', '10.1186/s12864-016-2364-4',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC4710005/',
    'Prior P. aeruginosa pangenome and molecular-evolution analysis.', limitation='Different collection, annotations and clustering; do not directly compare counts or exponents without aligning definitions.')
add('neu', 'Neu, H. C.', 'The role of Pseudomonas aeruginosa in infections.',
    'J. Antimicrob. Chemother. 11(Suppl B), 1–13 (1983).', '10.1093/jac/11.suppl_b.1',
    'https://pubmed.ncbi.nlm.nih.gov/6413479/', 'Historical broad clinical background.',
    'REVIEW', use='HISTORICAL_CONTEXT_ONLY', limitation='Does not establish contemporary global burden, ESKAPE nomenclature, current WHO ranking or present treatment recommendations.')
add('quinn', 'Quinn, J. P.',
    'Clinical problems posed by multiresistant nonfermenting gram-negative pathogens.',
    'Clin. Infect. Dis. 27(Suppl 1), S117–S124 (1998).', '10.1086/514912',
    'https://pubmed.ncbi.nlm.nih.gov/9710680/', 'Historical clinical background spanning several nonfermenting organisms.',
    'REVIEW', use='HISTORICAL_CONTEXT_ONLY', limitation='Not specific to P. aeruginosa; not evidence of susceptibility or prevalence in the current collection.')
add('fastani', 'Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. & Aluru, S.',
    'High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries.',
    'Nat. Commun. 9, 5114 (2018).', '10.1038/s41467-018-07641-9',
    'https://www.nature.com/articles/s41467-018-07641-9', 'FastANI methodology.',
    use='AUDIT_ONLY', limitation='No completed FastANI project output was located; a real paper cannot validate an undocumented run.')
add('modelfinder', 'Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A. & Jermiin, L. S.',
    'ModelFinder: fast model selection for accurate phylogenetic estimates.',
    'Nat. Methods 14, 587–589 (2017).', '10.1038/nmeth.4285',
    'https://www.nature.com/articles/nmeth.4285', 'Substitution-model selection implemented in the saved IQ-TREE run.', use='METHODS')
add('ufboot', 'Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q. & Vinh, L. S.',
    'UFBoot2: Improving the Ultrafast Bootstrap Approximation.',
    'Mol. Biol. Evol. 35, 518–522 (2018).', '10.1093/molbev/msx281',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC5850222/', 'Ultrafast bootstrap method.', use='METHODS',
    limitation='The publication does not override the project run-specific non-convergence warning.')
add('antismash', 'Blin, K. et al.',
    'antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation.',
    'Nucleic Acids Res. 51, W46–W50 (2023).', '10.1093/nar/gkad344',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC10320115/', 'antiSMASH 7 prediction method.',
    use='SUPPLEMENTARY_AUDIT_ONLY', limitation='The executed version must come from project logs, not this article title. Local region/CDS export remains unvalidated.')
add('deepbgc', 'Hannigan, G. D. et al.',
    'A deep learning genome-mining strategy for biosynthetic gene cluster prediction.',
    'Nucleic Acids Res. 47, e110 (2019).', '10.1093/nar/gkz654',
    'https://pubmed.ncbi.nlm.nih.gov/31400112/', 'DeepBGC method.', use='SUPPLEMENTARY_AUDIT_ONLY',
    limitation='Tool predictions are not experimental biosynthetic products; invalid downstream consensus cannot be repaired by citing the tool paper.')
add('gecco', 'Carroll, L. M. et al.',
    'Accurate de novo identification of biosynthetic gene clusters with GECCO.',
    'bioRxiv 2021.05.03.442509 (2021).', '10.1101/2021.05.03.442509',
    'https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1', 'GECCO method described in the cited preprint.',
    'PREPRINT', use='SUPPLEMENTARY_AUDIT_ONLY', limitation='No final journal replacement was verified; official GECCO documentation still cites this preprint. Identify it explicitly as a preprint.')
add('tettreview', 'Tettelin, H., Riley, D., Cattuto, C. & Medini, D.',
    'Comparative genomics: the bacterial pan-genome.', 'Curr. Opin. Microbiol. 11, 472–477 (2008).',
    '10.1016/j.mib.2008.09.006', 'https://pubmed.ncbi.nlm.nih.gov/19086349/',
    'Broad pangenome framework and accumulation formulations.', 'REVIEW', use='OPTIONAL_BACKGROUND',
    limitation='Cannot establish indefinite future expansion of the current finite sample. Prefer the 2005 primary paper for the origin of the pangenome concept.')
add('winsor', 'Winsor, G. L., Griffiths, E. J., Lo, R., Dhillon, B. K., Shay, J. A. & Brinkman, F. S. L.',
    'Enhanced annotations and features for comparing thousands of Pseudomonas genomes in the Pseudomonas genome database.',
    'Nucleic Acids Res. 44, D646–D653 (2016).', '10.1093/nar/gkv1227',
    'https://pubmed.ncbi.nlm.nih.gov/26578582/', 'Pseudomonas Genome Database resource.',
    'DATABASE_PUBLICATION', use='METHODS_RESOURCE', limitation='Does not prove the number or download date of files retrieved by this project.')
add('ncbi', 'Sayers, E. W. et al.', 'Database resources of the National Center for Biotechnology Information.',
    'Nucleic Acids Res. 50, D20–D26 (2022).', '10.1093/nar/gkab1112',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC8728269/', 'NCBI resources and database provenance.',
    'DATABASE_PUBLICATION', use='METHODS_RESOURCE', limitation='Does not establish the actual project download command, version or retrieval date.')

scopes = {
 'stover':'PAO1 genome as a foundation for understanding metabolic and regulatory versatility; one reference cannot represent all population variation.',
 'freschi':'Established population structure and pangenome diversity in P. aeruginosa; contextualises limitations of a first/unique-study claim.',
 'weimann':'Published global lineage and host-adaptation work; population structure is central to interpretation of source-associated genomic variation.',
 'tettelin':'Original population-genome analysis motivating the bacterial pangenome concept; not Pseudomonas-specific data.',
 'kislyuk':'Exact genomic-fluidity definition and gene-repertoire dissimilarity; not a horizontal-transfer or evolution-rate estimator.',
 'panaroo':'Panaroo pangenome reconstruction and the impact of annotation and assembly errors on accessory calls.',
 'penaranda':'Contemporary environmental/clinical comparison with paired phenotypic characterisation; challenges claims that ecological integration is wholly unexplored.',
 'amrfinder':'AMRFinderPlus identification of AMR, stress-response and virulence-associated genetic elements.',
 'ampir':'ampir prediction methodology; predictions require independent validation.',
 'amplify':'AMPlify method and experimental tests in that separate study; none transferred to the present project.',
 'macrel':'Macrel prediction methodology; no evidence that current bacterial predicted sequences have tested antimicrobial activity.',
 'dbcan':'dbCAN3 automated CAZyme/substrate annotation framework; executed run_dbCAN version and databases come from project records.',
 'pyseer':'An established structure-aware microbial association framework; cited only as a methodological option/limitation because it was not executed here.',
 'pizza':'The reverse-vaccinology concept; no completed vaccine-candidate pipeline exists in the current project.',
 'quast':'Assembly-quality measurement methodology.',
 'checkm':'CheckM2 completeness and contamination estimation methodology.',
 'gtdb':'GTDB-Tk v2 taxonomic classification methodology; actual project version/r232 and ANI pre-screen status derive from logs.',
 'prokka':'Prokka annotation methodology; exact command and database availability derive from local records.',
 'cdhit':'CD-HIT sequence-clustering algorithm; actual thresholds derive from project log.',
 'funcscan':'Workflow architecture and available modules; executed tool list must be demonstrated by local outputs.',
 'mmseqs':'MMseqs2 sequence search/clustering method; identity and coverage parameters derive from local scripts.',
 'deeparg':'DeepARG resistance-gene prediction methodology.',
 'card':'CARD resistance-ontology and curated reference resource used by RGI.',
 'bh':'Benjamini–Hochberg multiple-testing correction.',
 'iqtree':'IQ-TREE 2 maximum-likelihood phylogenetic methodology.'}
for key, scope in scopes.items():
    update(key, claim_scope=scope)
update('pyseer', recommended_use='DISCUSSION_FUTURE_ANALYSIS')
update('pizza', recommended_use='DISCUSSION_FUTURE_WORK')

evidence_assign = {
 'stover':'reference_verification_3.txt', 'freschi':'agent_reference_core_2.txt',
 'weimann':'agent_reference_core_2.txt', 'tettelin':'reference_verification_3.txt',
 'kislyuk':'reference_verification_0.txt', 'panaroo':'web_literature_batch_1.txt',
 'penaranda':'agent_reference_batch_0.txt', 'amrfinder':'reference_verification_0.txt',
 'ampir':'agent_reference_final_0.txt', 'amplify':'agent_reference_core_1.txt',
 'macrel':'reference_verification_0.txt', 'dbcan':'agent_reference_metadata.txt',
 'pyseer':'agent_reference_added.txt', 'pizza':'agent_reference_added.txt',
 'quast':'reference_verification_1.txt', 'checkm':'reference_verification_1.txt',
 'gtdb':'agent_reference_core_0.txt', 'prokka':'reference_verification_2.txt',
 'cdhit':'agent_reference_software.txt', 'funcscan':'reference_final_batch_1.txt',
 'mmseqs':'reference_verification_2.txt', 'deeparg':'agent_reference_amr_methods.txt',
 'card':'agent_reference_amr_methods.txt', 'bh':'reference_final_batch_3.txt',
 'iqtree':'agent_reference_software.txt', 'who':'agent_reference_metadata.txt',
 'sati':'agent_reference_core_2.txt', 'mosquera':'agent_reference_core_1.txt',
 'neu':'agent_reference_batch_0.txt', 'quinn':'agent_reference_verify_0.txt',
 'fastani':'agent_reference_core_1.txt', 'modelfinder':'agent_reference_core_0.txt',
 'ufboot':'agent_reference_core_0.txt', 'antismash':'agent_reference_final_0.txt',
 'deepbgc':'agent_reference_verify_0.txt', 'gecco':'agent_reference_batch_0.txt',
 'tettreview':'agent_reference_metadata.txt', 'winsor':'agent_reference_software.txt',
 'ncbi':'agent_reference_batch_1.txt'}
for key, f in evidence_assign.items():
    update(key, evidence_file='publication_manuscript/evidence_extracts/'+f)

def full(r):
    return f"{r['authors']} {r['title']} {r['bibliographic_details']}" + (f" doi:{r['doi']}." if r['doi'] else '')
for r in registry.values():
    r['full_citation'] = full(r)

(TAB/'reference_registry.json').write_text(json.dumps(registry, ensure_ascii=False, indent=2)+'\n')
with (TAB/'reference_registry.tsv').open('w',newline='') as fh:
    fields=list(next(iter(registry.values())).keys())
    w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(registry.values())

main=(OUT/'evidence_extracts/original_main_draft.txt').read_text()
comm=(OUT/'evidence_extracts/original_communication_draft.txt').read_text()
pre=main.split('\nReferences\n')[0]
counts={n:0 for n in range(1,30)}
for m in re.finditer(r'\[([\d,\s-]+)\]',pre):
    for part in m.group(1).split(','):
        part=part.strip()
        if '-' in part:
            a,b=map(int,part.split('-')); nums=range(a,b+1)
        else: nums=[int(part)]
        for n in nums:
            if n in counts:counts[n]+=1
entries={int(m.group(1)):m.group(0) for m in re.finditer(r'^(\d+)\. .+$',main,re.M) if 1<=int(m.group(1))<=29}
keymap={2:'neu',3:'quinn',4:'sati',5:'mosquera',6:'quast',7:'checkm',8:'gtdb',9:'prokka',10:'panaroo',11:'fastani',12:'iqtree',13:'modelfinder',14:'ufboot',15:'funcscan',16:'amrfinder',17:'deeparg',18:'antismash',19:'deepbgc',20:'gecco',21:'dbcan',22:'ampir',23:'amplify',24:'macrel',25:'kislyuk',26:'tettreview',27:'winsor',28:'ncbi',29:'bh'}
issues={
1:('UNTRACEABLE','Exact author/title/year/page combination could not be verified using searches; no DOI supplied.','Remove; use WHO and verified primary Pseudomonas literature for bounded background.'),
2:('CORRECTED','Journal, volume and pages are wrong; real article is J Antimicrob Chemother 11 Suppl B:1–13.','Correct; historical clinical review only.'),
3:('CORRECTED','Title narrowed incorrectly to P. aeruginosa; paper covers multiple nonfermenting pathogens.','Correct title and avoid attribution of all present resistance mechanisms to this review.'),
4:('CORRECTED_CLAIM_MISMATCH','Title and pages wrong; draft labels CRPA critical in updated WHO list although it is HIGH in 2024.','Correct reference and WHO wording.'),
5:('CORRECTED_CLAIM_MISMATCH','Title uses singular analysis; omnibus citation is overextended to essentiality, hypermutability, and indefinite future openness.','Correct plural analyses; cite only findings actually supported.'),
10:('VERIFIED_CLAIM_MISMATCH','Paper and metadata verified; >=99% core designation does not establish essentiality of every cluster or indefinite expansion.','Retain method citation; remove essentiality/unbounded-openness inference.'),
11:('CORRECTED_UNSUPPORTED_EXECUTION','Real paper has title High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries; no completed FastANI output located.','Correct in audit; omit FastANI as executed method.'),
15:('UNTRACEABLE_REPLACED','Fernandes/title/DOI preprint could not be verified; public official workflow repository exists.','Replace with official nf-core/funcscan software citation; describe executed version from local manifest.'),
16:('CORRECTED','Quoted title is wrong; paper addresses genomic links among antimicrobial resistance, stress response, and virulence.','Correct; tool records are not measured susceptibility phenotypes.'),
18:('VERIFIED_CONFLICTING_DOWNSTREAM','Tool publication is real; local biological BGC interpretation is blocked by region/CDS export inconsistency.','Keep audit citation; withhold downstream BGC biological results.'),
19:('CORRECTED_CONFLICTING_DOWNSTREAM','Wrong article title; exact title A deep learning genome-mining strategy for biosynthetic gene cluster prediction.','Correct in audit; BGC-derived findings withheld pending repair.'),
20:('CORRECTED_PREPRINT','Added phrase in bacterial and fungal genomes is not in verified title; source remains explicitly a preprint. No final journal replacement was verified.','Correct title and label preprint; BGC conclusions withheld.'),
22:('CORRECTED','Author initials, title, volume, pages and year wrong. Publisher gives Fingerhut L.C.H.W.; Bioinformatics36(21)5262–5263(2020).','Use publisher metadata and DOI btaa653.'),
23:('CORRECTED','Title ending application to loon genome is wrong; verified paper says effective against WHO priority pathogens.','Correct; its experiments do not validate current predictions.'),
24:('CORRECTED','Lead author wrong: Santos-Júnior, not Santos-Ferreira; title genomes and metagenomes.','Correct authors and title.'),
25:('CORRECTED','Title, year, volume and article number wrong; real fluidity paper BMCGenomics12:32(2011).','Replace with exact paper and equation; metric is repertoire dissimilarity.'),
26:('VERIFIED_SCOPE_LIMIT','Review verified; cumulative and new-gene exponent formulations must not be conflated.','Optional review; primary 2005 pangenome paper preferred for origin.'),
27:('CORRECTED','Title is enhanced annotations and features for comparing thousands of Pseudomonas genomes in the Pseudomonas genome database.','Correct; resource paper cannot validate collection counts.'),
28:('VERIFIED_SCOPE_LIMIT','Resource metadata verified; cannot support invented download version/date or counts.','Retain as resource citation only.')}
rows=[]
for n in range(1,30):
    key=keymap.get(n,'')
    r=registry.get(key,{})
    status, issue, action=issues.get(n,('VERIFIED','Bibliographic metadata verified against primary publication.','Retain when an executed method or bounded background claim is described.'))
    rows.append(dict(audit_id=f'R{n:03d}', original_document='manuscript_draft/Manuscript_draft.docx',
        original_reference_number=n, original_entry=entries[n], reference_id=key,
        verified_entry=r.get('full_citation','[NOT SUPPORTED BY PROJECT FILES]'),
        doi=r.get('doi',''), primary_verification_url=r.get('url',''), status=status,
        issue=issue, action=action, supported_claim_scope=r.get('claim_scope',''),
        unsupported_use=r.get('limitations',''), draft_in_text_citation_count=counts[n],
        recommended_master_use=r.get('recommended_use','EXCLUDE'), checked_date=DATE,
        evidence_file=r.get('evidence_file','publication_manuscript/evidence_extracts/agent_reference_verify_0.txt')))

extra=[
('R030','Communication Biology draft bibliography','TEMPLATE_RESIDUE','References section contains journal instructions and examples (Eigler/Schweizer, Iijima, Jones, GIDMaPS, quantum-scissors preprint), not an actual study bibliography.','Delete examples and reconstruct an actual bibliography.'),
('R031','Communication Biology draft author-year citations','CITED_BUT_BIBLIOGRAPHY_ABSENT','Asmare, Neu, Quinn, Sati, Mosquera-Rendón and other author-year mentions have no corresponding study bibliography in that draft.','Resolve each citation with verified registry; omit untraceable entries.'),
('R032','Original main numbered bibliography','NO_EXACT_DUPLICATES_DETECTED','29 distinct numbered entries; every number1–29 is cited at least once in extracted main body. This is an extraction check, not validation of attached claims.','Retain only references relevant to reconstructed scope; remove uncited references after final citation audit.'),
('R033','Unsupported rationale for novelty','MISSING_PRIOR_WORK','Original main list lacks central population-genomic studies Freschi2019, Weimann2024 and contemporary Peñaranda2026.','Add bounded context; do not claim first or largest study.'),
('R034','WHO terminology across drafts','OUTDATED_TERMINOLOGY','Updated 2024 WHO CRPA category is HIGH; draft critical/highest-priority phrases conflate species and resistance phenotype.','Use carbapenem-resistant P. aeruginosa is a high-priority pathogen in WHO2024.'),
('R035','Preprint replacement audit','NO_VERIFIED_FINAL_REPLACEMENT','GECCO preprint remains cited by official documentation; a final journal replacement was not verified. Fernandes funcscan preprint was untraceable.','Keep verified GECCO explicitly as a preprint in audit; replace funcscan with official software documentation.')]
for aid,entry,status,issue,action in extra:
    row={k:'' for k in rows[0]};row.update(audit_id=aid,original_document='manuscript_draft/pan/Pangenome_Pseudomonas Communication Biology.docx' if aid in ['R030','R031'] else 'Cross-draft audit',original_entry=entry,status=status,issue=issue,action=action,checked_date=DATE,evidence_file='publication_manuscript/evidence_extracts/original_communication_draft.txt' if aid in ['R030','R031'] else 'publication_manuscript/reconstructed_tables/reference_registry.json');rows.append(row)
for key in ['stover','freschi','weimann','penaranda','tettelin','who','pyseer','pizza','mmseqs','cdhit','card']:
    r=registry[key];row={k:'' for k in rows[0]};row.update(audit_id='NEW_'+key,original_document='Reconstructed master/context',reference_id=key,verified_entry=r['full_citation'],doi=r['doi'],primary_verification_url=r['url'],status='VERIFIED_ADDITION',issue='Verified primary or official source added for relevant bounded context/method.',action='Use only within claim_scope; no published result or experiment transferred to this collection.',supported_claim_scope=r['claim_scope'],unsupported_use=r['limitations'],recommended_master_use=r['recommended_use'],checked_date=DATE,evidence_file=r['evidence_file']);rows.append(row)
with (OUT/'07_REFERENCE_AUDIT.tsv').open('w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'registry_entries':len(registry),'audit_rows':len(rows),'main_refs':len(entries),'all_main_refs_cited':all(counts.values()),'outputs':['07_REFERENCE_AUDIT.tsv','reconstructed_tables/reference_registry.json','reconstructed_tables/reference_registry.tsv']}))
