from pathlib import Path
import zipfile,xml.etree.ElementTree as E,hashlib,json
p=next(Path('.').glob('*.docx')); b=p.read_bytes(); z=zipfile.ZipFile(p); ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
r=E.fromstring(z.read('word/document.xml')); out=[]
for i,e in enumerate(r.find('w:body',ns)):
 t=''.join(e.itertext()) if False else '\n'.join(''.join(p.itertext()) for p in [])
 paras=[''.join(n.text or '' for n in q.findall('.//w:t',ns)) for q in e.findall('.//w:p',ns)] if e.tag.endswith('}tbl') else [''.join(n.text or '' for n in e.findall('.//w:t',ns))]
 t='\n'.join(paras)
 out.append({'id':f'B{i:04d}','kind':e.tag.split('}')[-1],'text':t,'images':len(e.findall('.//{http://schemas.openxmlformats.org/drawingml/2006/main}blip'))})
Path('manuscript_work/source_blocks.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)); Path('manuscript_work/source_full.txt').write_text('\n\n'.join(f"[{x['id']}] ({x['kind']}, images={x['images']}) {x['text']}" for x in out));Path('manuscript_work/original_sha256.txt').write_text(hashlib.sha256(b).hexdigest()+'  '+str(p)+'\n')
print('blocks',len(out),'characters',sum(len(x['text']) for x in out));print('zip media', [n for n in z.namelist() if n.startswith('word/media/')]);print('other text parts',[n for n in z.namelist() if any(s in n for s in ['footnote','endnote','comment','header','footer'])])
