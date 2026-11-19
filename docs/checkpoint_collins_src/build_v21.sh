#!/bin/bash
# Build the Collins checkpoint Word copy from the page source (same method as v20: pandoc HTML -> docx with the landscape reference doc, tables widened).
set -e
SP=/tmp/claude-1000/-home-boris-Projects-ClrAttention/4451c478-b0bc-4b2f-84d1-11ecf74c2d08/scratchpad
cd $SP/v21build; rm -rf d base.docx v21.docx
nice -n 19 pandoc ../v21/index.html -f html -t docx --resource-path=../v21 --reference-doc=../export/ref/landscape.docx -o base.docx
mkdir d && cd d && unzip -q ../base.docx
python3 - <<'PY'
import re
s=open('word/document.xml',encoding='utf8').read()
TOTAL=14400
def fix(m):
    n=len(re.findall(r'<w:gridCol w:w="\d+" />',m.group(0)))
    w=[2.0]+[1.0]*(n-1) if n>2 else [1.0]*n
    t=sum(w)
    return '<w:tblGrid>'+''.join('<w:gridCol w:w="%d" />'%int(TOTAL*x/t) for x in w)+'</w:tblGrid>'
s=re.sub(r'<w:tblGrid>.*?</w:tblGrid>',fix,s,flags=re.S)
s=s.replace('<w:tblW w:type="auto" w:w="0" />','<w:tblW w:type="pct" w:w="5000" />')
open('word/document.xml','w',encoding='utf8').write(s)
PY
zip -q -r ../v21.docx .
ls -la ../v21.docx
