#!/usr/bin/env python3
"""Build docs/jim_short_preread.docx from jim_short_preread.md (pandoc + post-processing).
Method names are temporary and live in the ONE block below (and once per figure script:
analysis/fig_graphical_abstract_short{,_p2,_p3}.R, variables ACLR_NAME and LEAD_NAME). A rename = edit this block,
edit those variables, re-run the three figure scripts, re-run this script.
Run: nice -n 19 python3 build_short_preread.py [source.md] [out.docx]
"""
import re, subprocess, sys, zipfile, shutil, os, tempfile
from lxml import etree

# ---- NAMES: DEFINED HERE ONLY (temporary) ------------------------------------------------------------------------------
ACLR_NAME = "aCLR"                                   # CLR + attention (earlier-adopted method)
LEAD_NAME = "tcaCLR"                                 # lead method ("twice cooked aCLR")
LEAD_DESC = f"{ACLR_NAME} followed by a second CLR on the log of the walk"
# -------------------------------------------------------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "jim_short_preread.md")
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/boris/Projects/ClrAttention/docs/jim_short_preread.docx"
REF = os.path.join(HERE, "..", "checkpoint_collins_src", "portrait.docx")   # Word reference template (US Letter portrait, 0.75 in side and 0.6 in top and bottom margins), shared with build_v21.sh
PAGEBREAK = "PAGEBREAKMARKER"

md = open(SRC, encoding="utf-8").read()
for k, v in {"ACLR_NAME": ACLR_NAME, "LEAD_NAME": LEAD_NAME, "LEAD_DESC": LEAD_DESC, "PAGEBREAK": PAGEBREAK}.items():
    md = md.replace("{{%s}}" % k, v)
assert "{{" not in md, re.findall(r"\{\{[A-Z_]+\}\}", md)
tmp = tempfile.mkdtemp()
from PIL import Image, ImageChops
def crop_figs(text):
    def sub(m):
        p = m.group(1)
        if "graphical_abstract_short" not in p: return m.group(0)
        im = Image.open(p).convert("RGB")
        bbox = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
        pad = 10
        bbox = (max(bbox[0]-pad, 0), max(bbox[1]-pad, 0), min(bbox[2]+pad, im.size[0]), min(bbox[3]+pad, im.size[1]))
        q_ = os.path.join(tmp, os.path.basename(p)); im.crop(bbox).save(q_)
        return "](" + q_ + ")"
    return re.sub(r"\]\((/[^)]*\.png)\)", sub, text)
md = crop_figs(md)
rendered = os.path.join(tmp, "doc.md"); open(rendered, "w", encoding="utf-8").write(md)
raw = os.path.join(tmp, "raw.docx")
subprocess.run(["pandoc", rendered, "-f", "markdown+pipe_tables+raw_attribute+footnotes", "-o", raw,
                "--reference-doc", REF, "--wrap=none"], check=True)

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
ns = {"w": W}
def q(t): return "{%s}%s" % (W, t)

zin = zipfile.ZipFile(raw)
files = {n: zin.read(n) for n in zin.namelist()}

# ---- styles ----
st = files["word/styles.xml"].decode("utf-8")
st = re.sub(r"<w:docDefaults>.*?</w:docDefaults>",
            '<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri" w:eastAsia="Calibri"/>'
            '<w:sz w:val="21"/><w:szCs w:val="21"/><w:lang w:val="en-US"/></w:rPr></w:rPrDefault>'
            '<w:pPrDefault><w:pPr><w:spacing w:after="80" w:line="259" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>', st, flags=re.S)
def style_edit(sid, fn):
    global st
    m = re.search(r'<w:style [^>]*w:styleId="%s".*?</w:style>' % sid, st, re.S)
    if m: st = st.replace(m.group(0), fn(m.group(0)))
def fonts(x): return re.sub(r'<w:rFonts [^>]*/>', '<w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri"/>', x)
def sp(before, after): return lambda x: re.sub(r'<w:spacing [^>]*/>', '<w:spacing w:before="%d" w:after="%d"/>' % (before, after), x)
for sid, b, a, sz in [("Heading1", 240, 80, 30), ("Heading2", 160, 60, 25), ("Heading3", 120, 40, 22)]:
    style_edit(sid, lambda x, b=b, a=a, sz=sz: re.sub(r'<w:sz w:val="\d+"\s*/>', '<w:sz w:val="%d"/>' % sz,
               sp(b, a)(fonts(x))) if '<w:sz ' in x else sp(b, a)(fonts(x)).replace('</w:rPr>', '<w:sz w:val="%d"/></w:rPr>' % sz))
style_edit("Title", lambda x: re.sub(r'<w:jc [^>]*/>', '<w:jc w:val="left"/>', re.sub(r'<w:sz w:val="\d+"\s*/>', '', sp(0, 60)(fonts(x))).replace("</w:rPr>", '<w:sz w:val="36"/></w:rPr>')))
style_edit("Subtitle", lambda x: sp(0, 40)(fonts(x)))
style_edit("BodyText", sp(60, 80)); style_edit("FirstParagraph", sp(60, 80)); style_edit("Compact", sp(20, 20))
style_edit("Caption", lambda x: x.replace('<w:rPr>', '<w:rPr><w:sz w:val="18"/><w:szCs w:val="18"/>', 1) if '<w:rPr>' in x else x)
style_edit("ImageCaption", lambda x: x)
style_edit("Caption", lambda x: x.replace('<w:pPr>', '<w:pPr><w:keepLines/>', 1))   # a caption never splits across pages
style_edit("VerbatimChar", lambda x: re.sub(r'<w:sz w:val="\d+"\s*/>', '<w:sz w:val="16"/>', x))   # 8 pt code: the pseudocode lines fit the 7.0 in text width
st = st.replace("Consolas", "Courier New")
files["word/styles.xml"] = st.encode("utf-8")

# ---- document: tables, images, code ----
doc = etree.fromstring(files["word/document.xml"])
body = doc.find("w:body", ns)
usable = 10080          # text width in twips: 7.0 in (US Letter, 0.75 in side margins)
for tbl in doc.iter(q("tbl")):
    rows = tbl.findall("w:tr", ns)
    ncol = len(tbl.find("w:tblGrid", ns).findall("w:gridCol", ns))
    wts = [4] * ncol
    for r in rows:
        for i, c in enumerate(r.findall("w:tc", ns)):
            t = "".join(c.itertext()).strip()
            lw = max([len(w_) for w_ in t.split()] or [0])       # longest word: the column must hold it without a mid-word break
            wts[i] = max(wts[i], lw * 1.2 + 4.5, min(len(t), 42) * (0.6 if ncol >= 8 else 1.0) + 3)
    tsz = 16 if ncol >= 8 else (17 if ncol >= 6 else 18)         # table font in half-points: 8 pt for 8 or 9 columns
    ptsz = tsz / 2
    minw = [0.0] * ncol      # absolute minimum width (twips): the longest word of the column plus the cell margins
    for r in rows:
        for i, c in enumerate(r.findall("w:tc", ns)):
            t = "".join(c.itertext()).strip()
            minw[i] = max(minw[i], max([len(w_) for w_ in t.split()] or [0]) * 0.6 * ptsz * 20 + 280)
    tot = sum(wts); gw = [usable * w / tot for w in wts]
    for _ in range(6):       # raise narrow columns to their minimum, take the difference from the others in proportion
        low = [i for i in range(ncol) if gw[i] < minw[i]]
        rest = [i for i in range(ncol) if i not in low]
        if not low or not rest: break
        need = sum(minw[i] - gw[i] for i in low); sr = sum(gw[i] for i in rest)
        for i in rest: gw[i] -= need * gw[i] / sr
        for i in low: gw[i] = minw[i]
    cols = [int(w) for w in gw]
    pr = tbl.find("w:tblPr", ns)
    for ch in list(pr):
        if ch.tag in (q("tblW"), q("tblLayout"), q("tblBorders"), q("tblCellMar")): pr.remove(ch)
    pr.append(etree.fromstring('<w:tblW xmlns:w="%s" w:type="dxa" w:w="%d"/>' % (W, usable)))
    pr.append(etree.fromstring('<w:tblBorders xmlns:w="%s">%s</w:tblBorders>' % (W, "".join(
        '<w:%s w:val="single" w:sz="4" w:space="0" w:color="BFBFBF"/>' % s for s in ("top", "left", "bottom", "right", "insideH", "insideV")))))
    pr.append(etree.fromstring('<w:tblLayout xmlns:w="%s" w:type="fixed"/>' % W))
    pr.append(etree.fromstring('<w:tblCellMar xmlns:w="%s"><w:top w:w="40" w:type="dxa"/><w:left w:w="100" w:type="dxa"/><w:bottom w:w="40" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tblCellMar>' % W))
    # order inside tblPr is schema-sensitive: tblStyle, tblW, jc, tblBorders, tblLayout, tblCellMar, tblLook
    order = ["tblStyle", "tblpPr", "tblOverlap", "bidiVisual", "tblStyleRowBandSize", "tblStyleColBandSize", "tblW", "jc", "tblCellSpacing",
             "tblInd", "tblBorders", "shd", "tblLayout", "tblCellMar", "tblLook"]
    kids = sorted(list(pr), key=lambda e: order.index(etree.QName(e).localname) if etree.QName(e).localname in order else 99)
    for k in kids: pr.remove(k)
    for k in kids: pr.append(k)
    for gc, w in zip(tbl.find("w:tblGrid", ns).findall("w:gridCol", ns), cols): gc.set(q("w"), str(w))
    for ri, r in enumerate(rows):
        trpr = r.find("w:trPr", ns)
        if trpr is None:
            trpr = etree.SubElement(r, q("trPr")); r.remove(trpr); r.insert(0, trpr)
        if ri == 0 and trpr.find("w:tblHeader", ns) is None: etree.SubElement(trpr, q("tblHeader"))
        etree.SubElement(trpr, q("cantSplit"))
        if ri < len(rows) - 1:
            for p in r.iter(q("p")):
                ppr = p.find("w:pPr", ns)
                if ppr is None: ppr = etree.Element(q("pPr")); p.insert(0, ppr)
                if ppr.find("w:keepNext", ns) is None: ppr.insert(0, etree.Element(q("keepNext")))
        for i, c in enumerate(r.findall("w:tc", ns)):
            tcpr = c.find("w:tcPr", ns)
            if tcpr is None:
                tcpr = etree.Element(q("tcPr")); c.insert(0, tcpr)
            for ch in list(tcpr):
                if ch.tag in (q("tcW"), q("shd")): tcpr.remove(ch)
            tcpr.insert(0, etree.fromstring('<w:tcW xmlns:w="%s" w:w="%d" w:type="dxa"/>' % (W, cols[i])))
            if ri == 0: tcpr.append(etree.fromstring('<w:shd xmlns:w="%s" w:val="clear" w:color="auto" w:fill="DCE6F2"/>' % W))
            for rr in c.iter(q("r")):
                rpr = rr.find("w:rPr", ns)
                if rpr is None: rpr = etree.Element(q("rPr")); rr.insert(0, rpr)
                if rpr.find("w:sz", ns) is None:
                    rpr.append(etree.fromstring('<w:sz xmlns:w="%s" w:val="%d"/>' % (W, tsz))); rpr.append(etree.fromstring('<w:szCs xmlns:w="%s" w:val="%d"/>' % (W, tsz)))
                if ri == 0 and rpr.find("w:b", ns) is None: rpr.insert(0, etree.Element(q("b")))
# page breaks: drop the marker paragraph, put pageBreakBefore on the next block's first paragraph
for p in list(body.findall("w:p", ns)):
    if "PAGEBREAKMARKER" in "".join(p.itertext()):
        nxt = p.getnext(); body.remove(p)
        while nxt.tag not in (q("p"), q("tbl")): nxt = nxt.getnext()
        tgt = nxt if nxt.tag == q("p") else nxt.find(".//" + q("p"))
        ppr = tgt.find("w:pPr", ns)
        if ppr is None: ppr = etree.Element(q("pPr")); tgt.insert(0, ppr)
        ppr.insert(0, etree.Element(q("pageBreakBefore")))
# bold-only short paragraphs (run-in headings such as "A.1 ...") stay with the next block
for p in body.findall("w:p", ns):
    runs = p.findall("w:r", ns)
    txt_ = "".join(p.itertext()).strip()
    if runs and txt_ and len(txt_) < 120 and all(r.find("w:rPr/w:b", ns) is not None for r in runs):
        ppr = p.find("w:pPr", ns)
        if ppr is None: ppr = etree.Element(q("pPr")); p.insert(0, ppr)
        if ppr.find("w:keepNext", ns) is None: ppr.insert(0, etree.Element(q("keepNext")))
# keep figure + caption together: drawing paragraphs keepNext
for p in doc.iter(q("p")):
    if p.find(".//" + q("drawing"), ) is not None:
        ppr = p.find("w:pPr", ns)
        if ppr is None: ppr = etree.Element(q("pPr")); p.insert(0, ppr)
        if ppr.find("w:keepNext", ns) is None: ppr.insert(0, etree.Element(q("keepNext")))
# drop the picture-level alt text pandoc fills with the source file path (the figure-level alt text and caption stay)
for e in doc.iter("{http://schemas.openxmlformats.org/drawingml/2006/picture}cNvPr"): e.attrib.pop("descr", None)
files["word/document.xml"] = etree.tostring(doc, xml_declaration=True, encoding="UTF-8", standalone=True)

# page size and margins come from the reference's final sectPr (portrait.docx); keep
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    for n, b in files.items(): z.writestr(n, b)
shutil.rmtree(tmp)
print("wrote", OUT)
