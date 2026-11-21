#!/bin/bash
# Build the Collins checkpoint Word copy from the page source (pandoc HTML -> docx with the landscape reference doc, tables widened).
# Self-contained: reads checkpoint_collins_v21.html and figures/ in this folder, and landscape.docx (shared with
# docs/jim_short_preread_src/build_short_preread.py); writes docs/checkpoint_collins_v21_editable.docx (v20 is never touched).
# Run: nice -n 19 bash docs/checkpoint_collins_src/build_v21.sh
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/../checkpoint_collins_v21_editable.docx"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
cd "$W"
# pandoc reads the HTML and the PNGs from one folder: stage them together
mkdir src; cp "$HERE/checkpoint_collins_v21.html" src/index.html; cp "$HERE"/figures/*.png src/
nice -n 19 pandoc src/index.html -f html -t docx --resource-path=src --reference-doc="$HERE/landscape.docx" -o base.docx
mkdir d && cd d && unzip -q ../base.docx
python3 - <<'PY'
import re
s=open('word/document.xml',encoding='utf8').read()
TOTAL=14400
def fix(m):
    n=len(re.findall(r'<w:gridCol w:w="\d+" />',m.group(0)))
    w=[2.0]+[1.0]*(n-1) if n>2 else [1.0]*n
    if n>=9: w=[2.0,1.6]+[1.0]*(n-2)   # nine-column table: give the Truth column room
    t=sum(w)
    return '<w:tblGrid>'+''.join('<w:gridCol w:w="%d" />'%int(TOTAL*x/t) for x in w)+'</w:tblGrid>'
s=re.sub(r'<w:tblGrid>.*?</w:tblGrid>',fix,s,flags=re.S)
s=s.replace('<w:tblW w:type="auto" w:w="0" />','<w:tblW w:type="pct" w:w="5000" />')
s=re.sub(r'(<pic:cNvPr )descr="[^"]*" ?',r'\1',s)   # picture-level alt text is the source file name: drop it
open('word/document.xml','w',encoding='utf8').write(s)
PY
rm -f ../v21.docx; zip -q -r ../v21.docx .
cp ../v21.docx "$OUT"
ls -la "$OUT"
