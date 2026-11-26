#!/bin/bash
# Build the Collins checkpoint Word copy from the page source (pandoc HTML -> docx with the portrait reference doc, tables re-weighted).
# Self-contained: reads checkpoint_collins_v21.html and figures/ in this folder, and portrait.docx (shared with
# docs/jim_short_preread_src/build_short_preread.py); writes docs/checkpoint_collins_v21_editable.docx (v20 is never touched).
# Run: nice -n 19 bash docs/checkpoint_collins_src/build_v21.sh
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/../checkpoint_collins_v21_editable.docx"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
cd "$W"
# pandoc reads the HTML and the PNGs from one folder: stage them together
mkdir src; cp "$HERE/checkpoint_collins_v21.html" src/index.html; cp "$HERE"/figures/*.png src/
nice -n 19 pandoc src/index.html -f html+tex_math_dollars -t docx --resource-path=src --reference-doc="$HERE/portrait.docx" -o base.docx
mkdir d && cd d && unzip -q ../base.docx
python3 - <<'PY'
import re
from lxml import etree
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
ns={'w':W}
def q(t): return '{%s}%s'%(W,t)
TOTAL=10080   # text width in twips: 7.0 in (US Letter, 0.75 in left and right margins)
raw=open('word/document.xml','rb').read()
doc=etree.fromstring(raw)
def cells_of(tr, n):
    """(first column index, span, cell) for each cell of a row"""
    out=[]; c=0
    for tc in tr.findall('w:tc',ns):
        gs=tc.find('w:tcPr/w:gridSpan',ns); sp=int(gs.get(q('val'))) if gs is not None else 1
        out.append((c,sp,tc)); c+=sp
    return out
for tbl in doc.iter(q('tbl')):
    grid=tbl.find('w:tblGrid',ns); cols=grid.findall('w:gridCol',ns); n=len(cols)
    ptsz=7.0 if n>=12 else (8.0 if n>=8 else (8.5 if n>=6 else 10.5))
    marg=80 if n>=12 else 216                      # cell margin left + right, in twips (narrow for the 16-column tables)
    # column weights: wide enough for the longest word of every cell (no mid-word breaks) and roughly for the cell text
    wts=[4.0]*n; minw=[4.0]*n
    for tr in tbl.findall('w:tr',ns):
        for c,sp,tc in cells_of(tr,n):
            if sp!=1 or c>=n: continue             # a header cell that spans columns does not set the width of one column
            t=''.join(tc.itertext()).strip()
            lw=max([len(x) for x in t.split()] or [0])
            wts[c]=max(wts[c], lw*1.2+4.5, min(len(t),40)*(0.6 if n>=8 else 0.8)+3)
            minw[c]=max(minw[c], lw*0.6*ptsz*20+marg+110)
    if n>=9: wts[1]=max(wts[1],wts[0]*0.8)
    tot=sum(wts); gw=[TOTAL*wt/tot for wt in wts]
    for _ in range(8):    # raise narrow columns to their minimum, take the difference from the others in proportion
        low=[i for i in range(n) if gw[i]<minw[i]]
        if not low: break
        need=sum(minw[i]-gw[i] for i in low); rest=[i for i in range(n) if i not in low]
        if not rest: break
        sr=sum(gw[i] for i in rest)
        for i in rest: gw[i]-=need*gw[i]/sr
        for i in low: gw[i]=minw[i]
    gwi=[int(x) for x in gw]
    for gc,wd in zip(cols,gwi): gc.set(q('w'),str(wd))
    pr=tbl.find('w:tblPr',ns)
    for ch in list(pr):
        if ch.tag in (q('tblW'),q('tblLayout'),q('tblCellMar')): pr.remove(ch)
    etree.SubElement(pr,q('tblW')).attrib.update({q('type'):'dxa',q('w'):str(sum(gwi))})
    etree.SubElement(pr,q('tblLayout')).set(q('type'),'fixed')      # the column weights above are honoured (no autofit)
    if n>=12:
        cm=etree.SubElement(pr,q('tblCellMar'))
        for side in ('top','left','bottom','right'):
            e=etree.SubElement(cm,q(side)); e.set(q('w'),'40' if side in ('left','right') else '0'); e.set(q('type'),'dxa')
    order=['tblStyle','tblpPr','tblOverlap','bidiVisual','tblStyleRowBandSize','tblStyleColBandSize','tblW','jc','tblCellSpacing','tblInd','tblBorders','shd','tblLayout','tblCellMar','tblLook','tblCaption','tblDescription']
    kids=sorted(list(pr),key=lambda e: order.index(etree.QName(e).localname) if etree.QName(e).localname in order else 99)
    for k in kids: pr.remove(k)
    for k in kids: pr.append(k)
    for tr in tbl.findall('w:tr',ns):
        for c,sp,tc in cells_of(tr,n):
            tcpr=tc.find('w:tcPr',ns)
            if tcpr is None: tcpr=etree.Element(q('tcPr')); tc.insert(0,tcpr)
            for ch in list(tcpr):
                if ch.tag==q('tcW'): tcpr.remove(ch)
            tcw=etree.Element(q('tcW')); tcw.attrib.update({q('w'):str(sum(gwi[c:c+sp])),q('type'):'dxa'}); tcpr.insert(0,tcw)
    sz=14 if n>=12 else (16 if n>=8 else (17 if n>=6 else None))   # half-points: 7 pt for 12+ columns, 8 pt for 8-9, 8.5 pt for 6-7
    if sz:
        for r in tbl.iter(q('r')):
            rpr=r.find('w:rPr',ns)
            if rpr is None:
                rpr=etree.Element(q('rPr')); r.insert(0,rpr)
            for tag in ('sz','szCs'):
                e=rpr.find('w:'+tag,ns)
                if e is None: e=etree.SubElement(rpr,q(tag))
                e.set(q('val'),str(sz))
# every picture fits the text area (7.0 in wide, at most 9.7 in tall (US Letter, 0.6 in top and bottom margins) so a full-page figure stays on one page)
WP='http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'; A='http://schemas.openxmlformats.org/drawingml/2006/main'
EMU=914400; MAXW=int(7.0*EMU); MAXH=int(9.7*EMU)
for inl in doc.iter('{%s}inline'%WP):
    ext=inl.find('{%s}extent'%WP); cx=int(ext.get('cx')); cy=int(ext.get('cy'))
    f=min(1.0, MAXW/cx, MAXH/cy)
    if f<1.0:
        ncx,ncy=int(cx*f),int(cy*f); ext.set('cx',str(ncx)); ext.set('cy',str(ncy))
        for e in inl.iter('{%s}ext'%A):
            if e.get('cx') is not None: e.set('cx',str(ncx)); e.set('cy',str(ncy))
s=etree.tostring(doc,xml_declaration=True,encoding='UTF-8',standalone=True).decode('utf8')
s=re.sub(r'(<pic:cNvPr )descr="[^"]*" ?',r'\1',s)   # picture-level alt text is the source file name: drop it
open('word/document.xml','w',encoding='utf8').write(s)
st=open('word/styles.xml',encoding='utf8').read()
def sedit(sid,fn):
    global st
    m=re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>'%sid,st,re.S)
    if m: st=st.replace(m.group(0),fn(m.group(0)))
sedit('Caption',lambda x:x.replace('<w:pPr>','<w:pPr><w:keepLines/>',1))                          # a caption never splits across pages
sedit('VerbatimChar',lambda x:re.sub(r'<w:sz w:val="\d+"\s*/>','<w:sz w:val="16"/>',x))             # 8 pt code: pseudocode lines fit the 7.0 in text width
open('word/styles.xml','w',encoding='utf8').write(st)
PY
rm -f ../v21.docx; zip -q -r ../v21.docx .
cp ../v21.docx "$OUT"
ls -la "$OUT"
