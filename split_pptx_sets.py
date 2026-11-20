"""Step 0: split one PowerPoint file that holds several specimens into the per-specimen input
folders that parse_traces.py expects (<input>/[Group][n]_raw/ppt/...).

Expected order of slides in the file
  * title slide of a specimen: text only (no picture, no line), e.g. 'D1_1'
  * then its section slides: one CT section picture each, with the traced red connectors
  * slides with traces but no picture (e.g. a merged overview of all sections) are skipped
  * slides before the first title slide are skipped
The specimen name is taken from the title text <code>_<n>; --map renames the code, e.g.
--map D1=Low D3=High turns 'D1_2' into the folder [Low][2]_raw (set name Low_2).

Each output folder holds presentation.xml with only the title slide and the section slides of
that specimen, those slides with their relationships, and the pictures they use (copied).

Usage:  python split_pptx_sets.py PPTX [--out raw/unpacked] [--map D1=Low D3=High]
"""
import os, re, shutil, zipfile, argparse
from xml.etree import ElementTree as ET

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('pptx')
ap.add_argument('--out', default=os.path.join('raw', 'unpacked'), help='input folder of parse_traces.py')
ap.add_argument('--map', nargs='*', default=[], help='rename title codes, e.g. D1=Low D3=High')
args = ap.parse_args()
NSP = 'http://schemas.openxmlformats.org/presentationml/2006/main'
NSA = 'http://schemas.openxmlformats.org/drawingml/2006/main'
NSR = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
for k, v in (('p', NSP), ('a', NSA), ('r', NSR)):
    ET.register_namespace(k, v)
RENAME = dict(m.split('=', 1) for m in args.map)

z = zipfile.ZipFile(args.pptx)
pres_xml = z.read('ppt/presentation.xml')
prels = ET.fromstring(z.read('ppt/_rels/presentation.xml.rels'))
rid2t = {r.get('Id'): r.get('Target') for r in prels}
pres = ET.fromstring(pres_xml)
order = [(s.get('{%s}id' % NSR), rid2t[s.get('{%s}id' % NSR)]) for s in pres.find('{%s}sldIdLst' % NSP)]

def slide_info(target):
    root = ET.fromstring(z.read('ppt/' + target))
    texts = [t.text.strip() for t in root.iter('{%s}t' % NSA) if t.text and t.text.strip()]
    return (len(list(root.iter('{%s}pic' % NSP))), len(list(root.iter('{%s}cxnSp' % NSP))), texts)

specimens, cur = [], None
for rid, t in order:
    npic, ncx, texts = slide_info(t)
    if npic == 0 and ncx == 0 and texts:
        m = re.match(r'^([A-Za-z0-9]+)_(\w+)$', texts[0])
        if m:
            cur = dict(title=texts[0], code=m.group(1), n=m.group(2), slides=[(rid, t)], skipped=[])
            specimens.append(cur); continue
    if cur is None:
        continue
    if npic == 1:
        cur['slides'].append((rid, t))
    else:
        cur['skipped'].append(os.path.basename(t))

for sp in specimens:
    group = RENAME.get(sp['code'], sp['code'])
    d = os.path.join(args.out, '[%s][%s]_raw' % (group, sp['n']))
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(os.path.join(d, 'ppt', '_rels'))
    os.makedirs(os.path.join(d, 'ppt', 'slides', '_rels'))
    with open(os.path.join(d, 'ppt', '_rels', 'presentation.xml.rels'), 'wb') as f:
        f.write(z.read('ppt/_rels/presentation.xml.rels'))
    p = ET.fromstring(pres_xml)
    lst = p.find('{%s}sldIdLst' % NSP)
    keep = {rid for rid, _ in sp['slides']}
    for s in list(lst):
        if s.get('{%s}id' % NSR) not in keep:
            lst.remove(s)
    ET.ElementTree(p).write(os.path.join(d, 'ppt', 'presentation.xml'), xml_declaration=True, encoding='UTF-8')
    n_media = 0
    for _, t in sp['slides']:
        name = os.path.basename(t)
        rel_name = 'slides/_rels/%s.rels' % name
        with open(os.path.join(d, 'ppt', 'slides', name), 'wb') as f:
            f.write(z.read('ppt/' + t))
        rels_xml = z.read('ppt/' + rel_name)
        with open(os.path.join(d, 'ppt', rel_name), 'wb') as f:
            f.write(rels_xml)
        for r in ET.fromstring(rels_xml):
            tgt = r.get('Target')
            if r.get('TargetMode') == 'External' or not tgt.startswith('../media/'):
                continue
            src = 'ppt/media/' + os.path.basename(tgt)
            dst = os.path.join(d, 'ppt', 'media', os.path.basename(tgt))
            if not os.path.exists(dst):
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                with open(dst, 'wb') as f:
                    f.write(z.read(src))
                n_media += 1
    print('%-18s <- title %-6s  section slides=%2d  pictures=%2d  skipped=%s' % (
        os.path.basename(d), sp['title'], len(sp['slides']) - 1, n_media, ', '.join(sp['skipped']) or '-'))
