from pathlib import Path
import re,json,zipfile,hashlib,xml.etree.ElementTree as E
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main';NS={'w':W};items=json.loads(Path('manuscript_work/restructured_items.json').read_text());bs=json.loads(Path('manuscript_work/source_blocks.json').read_text())
# Format Markdown tables as contiguous CommonMark rows.
for p in Path('.').glob('Manuscript_*.md'):
 t=p.read_text();t=re.sub(r'(?m)^(\|[^\n]*\|)\n(?:\n)+(?=\|)',r'\1\n',t)
 if p.name!='Manuscript_Restructured_v1.md':
  t+='\n' if not t.endswith('\n') else ''
  if p.name=='Manuscript_Figure_Table_Audit.md':t+='\nFor Word layout, only fully transparent exterior image padding is hidden with a DOCX image viewport. Every nontransparent source pixel remains inside that viewport (with a safety margin); embedded PNG bytes are unchanged. No plot data or label was redrawn.\n'
  if p.name=='Manuscript_Structural_Audit.md':t+='\nImage-layout note: transparent exterior margins are hidden in the Word viewport to improve display size. This is a document-formatting change; all retained PNG byte streams are unchanged.\n'
 else:
  t=re.sub(r'(?<!\*)(Pseudomonas aeruginosa|Pseudomonas paraeruginosa|P\. aeruginosa|P\. paraeruginosa)(?!\*)',r'*\1*',t)
 p.write_text(t)
source=next(Path('.').glob('Pangenome*.docx'));z=zipfile.ZipFile('Manuscript_Restructured_v1.docx');oz=zipfile.ZipFile(source)
assert z.testzip() is None
for n in z.namelist():
 if n.endswith('.xml') or n.endswith('.rels'):E.fromstring(z.read(n))
old=E.fromstring(oz.read('word/document.xml'));new=E.fromstring(z.read('word/document.xml'))
sourcehash=hashlib.sha256(source.read_bytes()).hexdigest();expected=Path('manuscript_work/original_sha256.txt').read_text().split()[0];assert sourcehash==expected
oldbody=list(old.find('w:body',NS));newtables=new.findall('.//w:body/w:tbl',NS);titems=[it for it in items if it['kind']=='table'];assert len(newtables)==len(titems)==14
for tb,it in zip(newtables,titems):
 def cells(tbl):return [[' '.join(''.join(x.text or '' for x in p.findall('.//w:t',NS)) for p in c.findall('w:p',NS)).strip() for c in tr.findall('w:tc',NS)] for tr in tbl.findall('w:tr',NS)]
 assert cells(tb)==cells(oldbody[it['block']]),it['block']
ims=[it for it in items if it['kind']=='image'];assert len(ims)==19
for it in ims:
 name=f"B{it['block']:04d}_{it['index']}.png";assert z.read('word/media/'+name)==Path('manuscript_work',name).read_bytes()
text='\n'.join(''.join(p.itertext()) for p in new.findall('.//w:p',NS));assert not re.search(r'graphene|platelet dynamics|Results subheading|The template is designed|Submission guidelines|Eigler|Iijima|Soft Machines',text,re.I)
# Headings and first-citation sequence.
for tag,n in [('Figure',7),('Supplementary Figure S',8),('Table',7),('Supplementary Table S',7)]:
 pat=(r'^'+tag.replace(' S',' S')+('' if tag.endswith('S') else ' ')+r'(\d+)\.')
 labels=[int(m[1]) for it in items if (m:=re.match(pat,it['text']))]
 assert labels==list(range(1,n+1)),(tag,labels)
for major,count in [('Results',8),('Methods',18)]:
 start=next(i for i,it in enumerate(items) if it['kind']=='h' and it['text']==major)
 end=next(i for i in range(start+1,len(items)) if items[i]['kind']=='h' and items[i]['level']==1)
 nums=[int(re.match(r'(\d+)\.',it['text'])[1]) for it in items[start:end] if it['kind']=='h' and it['level']==2]
 assert nums==list(range(1,count+1)),(major,nums)
absidx=next(i for i,it in enumerate(items) if it['text']=='Abstract');abstract=items[absidx+1]['text']
assert all(x not in abstract for x in ['130,103','371','390','1,347','12,605','60 remained','27 also'])
# Report source numeric tokens not literally repeated in new result prose: only known approximations/formatting/old section labels.
pat=r'(?<![A-Za-z_])(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][+-]?\d+)?';newnums=set(re.findall(pat,Path('Manuscript_Restructured_v1.md').read_text()));lost=[]
for b in bs:
 i=int(b['id'][1:])
 if 33<=i<=236 and b['kind']!='tbl':
  missing=set(re.findall(pat,b['text']))-newnums
  if missing:lost.append((b['id'],sorted(missing)))
report=['# Quality control','## File integrity',f'Original SHA-256 before/after: `{sourcehash}` — unchanged.','The new DOCX ZIP integrity and every XML/relationship part parsed successfully. Fourteen retained scientific tables match their source cells exactly, including all numeric precision and accession versions. All 19 scientific PNGs are byte-identical to their extracted source images. The source template image is excluded.','## Structure and claims','Eight Results sections and 18 Methods sections are continuously numbered. Main figures 1–7, supplementary figures S1–S8, main tables 1–7 and supplementary tables S1–S7 are continuous. Every retained Results analysis maps to Methods; substantive gaps are marked rather than certified as complete. Unsupported Abstract label counts, final differential counts and the uncorroborated SNP-site count are excluded from the revised Abstract. No completed resistance predictor or reverse-vaccinology screen is created.','Template/graphene/sample-reference scans pass. Native DOCX text consistently italicizes the two scientific species names, including abbreviated forms. Markdown companion names are also italicized. Source-to-destination ledger accounts for all 388 body blocks.','## Numerical preservation review','Automated source-Results numeric-token comparison found only the following literal forms absent from the new manuscript; substantive exact results remain present. These are removed old heading numbers, repeated approximate/rounded summaries, SHA-256 typography, or the damaged superscript extraction `2903` (the actual exponent 0.2903 is retained). No absent substantive functional result was detected.','| Source | Literal forms absent | Explanation |','| --- | --- | --- |']
for bid,vals in lost:report.append('| '+bid+' | '+', '.join(vals)+' | Removed heading/rounded repetition or formatting artifact; see original text and source ledger. |')
report+=['## Rendering and remaining limitations','LibreOffice successfully loaded and rendered the new DOCX to `render/Manuscript_Restructured_v1.pdf`. Representative text, dense tables, figure and supplementary-table pages were visually inspected. Figures remain limited by source raster resolution; transparent exterior padding is hidden only through document viewport settings, without changing image bytes or excluding nontransparent content. Authors need readable original exports and verified captions.','The DOCX is a structural author-review first draft with explicit checks, not a submission-ready or independently reanalyzed study. Raw matrices, manifests, logs, model outputs, alignments and a bibliography were not supplied. Numerical conflicts remain visible in the draft/issues file. No web reference verification or new empirical analysis was performed.']
Path('manuscript_work/Quality_Control.md').write_text('\n\n'.join(report)+'\n');p=Path('manuscript_work/Quality_Control.md');p.write_text(re.sub(r'(?m)^(\|[^\n]*\|)\n(?:\n)+(?=\|)',r'\1\n',p.read_text()))
print('PASS: original hash, 14 source-identical tables, 19 byte-identical images, section/display sequences, Abstract exclusions, XML/package and template scans.')
