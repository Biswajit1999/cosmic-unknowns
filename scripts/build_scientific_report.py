"""Build the formal report from its Markdown source. Author: Biswajit Jana.

Run render_report_equations.py first in the numerical environment. This builder
needs reportlab, pypdf and Pillow, available in the bundled PDF runtime.
"""
import hashlib
import html
import json
from pathlib import Path
import re
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, KeepTogether, PageBreak
from PIL import Image as PILImage
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"reports/007-calibration-limited-identification.md"
OUTPUT = ROOT/"output/pdf/calibration-limited-identification.pdf"
EQUATIONS = ROOT/"tmp/pdfs/equations"


def rich(text):
    text = html.escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`([^`]+)`", r'<font name="Sans" size="8.8">\1</font>', text)
    text = re.sub(r"\b([pqfgY])_i\b", r"<i>\1</i><sub>i</sub>", text)
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r'<link href="\2" color="#175a73">\1</link>', text)
    # DOI URLs in the references stay clickable and visible, with break opportunities.
    text = re.sub(r"(?<![\"=])(https?://[^\s<]+)", lambda m: '<link href="'+m[1]+'" color="#175a73">'+m[1]+'</link>', text)
    return text


def main():
    fontdir = ROOT/".venv/Lib/site-packages/matplotlib/mpl-data/fonts/ttf"
    if not fontdir.exists():
        import matplotlib
        fontdir = Path(matplotlib.get_data_path())/"fonts/ttf"
    for name, filename in [("Body", "DejaVuSerif.ttf"), ("Body-Bold", "DejaVuSerif-Bold.ttf"),
                           ("Body-Italic", "DejaVuSerif-Italic.ttf"), ("Sans", "DejaVuSans.ttf")]:
        pdfmetrics.registerFont(TTFont(name, str(fontdir/filename)))
    pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold", italic="Body-Italic", boldItalic="Body-Bold")
    style = {
        "body": ParagraphStyle("body", fontName="Body", fontSize=10, leading=14.3, spaceAfter=8, alignment=TA_JUSTIFY),
        "title": ParagraphStyle("title", fontName="Body-Bold", fontSize=22, leading=28, spaceAfter=14, alignment=TA_CENTER),
        "front": ParagraphStyle("front", fontName="Body", fontSize=11, leading=15, spaceAfter=8, alignment=TA_CENTER),
        "h2": ParagraphStyle("h2", fontName="Body-Bold", fontSize=13, leading=17, spaceBefore=15, spaceAfter=8, keepWithNext=True),
        "h3": ParagraphStyle("h3", fontName="Body-Bold", fontSize=10.8, leading=15, spaceBefore=9, spaceAfter=6, keepWithNext=True),
        "caption": ParagraphStyle("caption", fontName="Body", fontSize=8.5, leading=11.4, spaceAfter=11),
        "cell": ParagraphStyle("cell", fontName="Body", fontSize=8.6, leading=11.5),
        "ref": ParagraphStyle("ref", fontName="Body", fontSize=8.4, leading=11.8, spaceAfter=8, wordWrap="CJK"),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    width = A4[0]-108
    doc = SimpleDocTemplate(str(OUTPUT), pagesize=A4, rightMargin=54, leftMargin=54,
        topMargin=57, bottomMargin=52, title="Calibration-Limited Model Identification in Scalar Transmission",
        author="Biswajit Jana", subject="Reproducible theoretical benchmark; not evidence of physical wormhole existence")
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    flow, i, eq, front, refs = [], 0, 0, True, False
    def image_flow(path, maxheight=365):
        with PILImage.open(path) as im:
            w, h = im.size
        factor = min(width/w, maxheight/h)
        return Image(str(path), width=w*factor, height=h*factor)
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith("# "):
            flow.append(Paragraph(rich(line[2:]), style["title"]))
            i += 1
        elif line.startswith("## "):
            front = False
            refs = line == "## References"
            if refs:
                flow.append(PageBreak())
            flow.append(Paragraph(rich(line[3:]), style["h2"]))
            i += 1
        elif line.startswith("### "):
            flow.append(Paragraph(rich(line[4:]), style["h3"]))
            i += 1
        elif line == "$$":
            eq += 1
            i += 1
            while i < len(lines) and lines[i].strip() != "$$":
                i += 1
            i += 1
            path = EQUATIONS/f"equation-{eq:02d}.png"
            with PILImage.open(path) as im:
                w, h = im.size
            factor = min(72/300, (width-38)/w)
            table = Table([[Image(str(path), width=w*factor, height=h*factor), Paragraph(f"({eq})", style["cell"])]], colWidths=[width-35, 35])
            table.setStyle(TableStyle([("VALIGN", (0,0), (-1,-1), "MIDDLE"), ("ALIGN", (0,0), (0,0), "CENTER"), ("ALIGN", (1,0), (1,0), "RIGHT"), ("TOPPADDING", (0,0), (-1,-1), 5), ("BOTTOMPADDING", (0,0), (-1,-1), 9)]))
            flow.append(table)
        elif line.startswith("!["):
            match = re.match(r"!\[(.*?)\]\((.*?)\)", line)
            path = (SOURCE.parent/match[2]).resolve()
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            caption = ""
            if i < len(lines) and lines[i].startswith("**Figure"):
                caption = lines[i]
                i += 1
            flow.append(KeepTogether([Spacer(1,7), image_flow(path), Spacer(1,5), Paragraph(rich(caption), style["caption"])]))
        elif line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = [v.strip() for v in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r"[-: ]+", v) for v in row):
                    rows.append([Paragraph(rich(v), style["cell"]) for v in row])
                i += 1
            tw = [width*.5]+[width*.5/(len(rows[0])-1)]*(len(rows[0])-1)
            table = Table(rows, colWidths=tw, repeatRows=1, hAlign="LEFT")
            table.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,0), colors.HexColor("#edf2f4")),
                ("LINEBELOW", (0,0), (-1,0), .8, colors.HexColor("#83939d")),
                ("LINEBELOW", (0,1), (-1,-1), .3, colors.HexColor("#d8dfe3")),
                ("VALIGN", (0,0), (-1,-1), "TOP"), ("TOPPADDING", (0,0), (-1,-1), 7),
                ("BOTTOMPADDING", (0,0), (-1,-1), 7)]))
            flow += [Spacer(1,5), table, Spacer(1,10)]
        else:
            words = [line]
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "$$", "![", "|")):
                words.append(lines[i].strip())
                i += 1
            flow.append(Paragraph(rich(" ".join(words)), style["front" if front else "ref" if refs else "body"]))
    def page(canvas, document):
        canvas.saveState()
        canvas.setFont("Sans", 8)
        canvas.setFillColor(colors.HexColor("#526370"))
        if document.page > 1:
            canvas.drawString(54, A4[1]-31, "Jana | Calibration-limited scalar-transmission identification")
        canvas.setStrokeColor(colors.HexColor("#c6d0d6"))
        canvas.line(54, 39, A4[0]-54, 39)
        canvas.drawString(54, 27, "Public working report | 10 November 2026 | Specialist review pending")
        canvas.drawRightString(A4[0]-54, 27, str(document.page))
        canvas.restoreState()
    doc.build(flow, onFirstPage=page, onLaterPages=page)
    reader = PdfReader(OUTPUT)
    text = "\n".join(p.extract_text() or "" for p in reader.pages)
    assert "Proposition 1" in text and "References" in text and len(reader.pages) >= 6
    assert reader.metadata.author == "Biswajit Jana"
    manifest = {"author": "Biswajit Jana", "source": SOURCE.relative_to(ROOT).as_posix(),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "pdf_sha256": hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        "pages": len(reader.pages), "equations": eq, "figures": 5,
        "build": "render_report_equations.py then build_scientific_report.py; numerical environment plus PDF runtime"}
    (OUTPUT.parent/"report-manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
