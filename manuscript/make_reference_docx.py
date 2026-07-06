#!/usr/bin/env python3
"""Generate reference.docx — the Word style template Quarto/Pandoc copies into
the manuscript's .docx output (wired via `reference-doc:` in manuscript.qmd).

Encodes a journal-standard manuscript spec:
  - Times New Roman, 12 pt body text
  - Double line spacing throughout (body, abstract, references, captions)
  - 1-inch margins on all four sides
  - Left-aligned paragraphs (never justified)
  - 0.5-inch first-line indent on body paragraphs
  - Continuous page numbers, top-right header
  - Continuous line numbering (for peer review)
  - Headings 13-20 pt via Word's built-in Heading styles
  - Tables with horizontal rules only (no vertical lines)

Regenerate with:  python3 make_reference_docx.py
It starts from Pandoc's own default reference doc (so every style id Pandoc
expects is present) and restyles it in place.
"""
import subprocess
import copy
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BODY_FONT = "Times New Roman"
BODY_PT = 12
INDENT = Inches(0.5)
DEFAULT_REF = "/private/tmp/claude-501/-Users-taavi-Projects-hpv-human-skin/c3a75b48-7370-49a1-8e1f-a4ec9cc83d67/scratchpad/ref_default.docx"
OUT = "reference.docx"

# Styles that get a first-line indent (running prose only).
INDENT_STYLES = {"Normal", "Body Text", "First Paragraph", "Block Text"}
# Prose-ish styles that stay double-spaced but with NO first-line indent.
NO_INDENT_STYLES = {
    "Compact", "Abstract", "Bibliography", "Caption", "Table Caption",
    "Image Caption", "Figure", "Figure with Caption", "Footnote Text",
    "Definition", "Definition Term", "Title", "Subtitle", "Author", "Date",
    "Heading 1", "Heading 2", "Heading 3", "Heading 4", "Heading 5", "Heading 6",
}
HEADING_PT = {
    "Title": 20, "Subtitle": 14,
    "Heading 1": 16, "Heading 2": 14, "Heading 3": 13,
    "Heading 4": 12, "Heading 5": 12, "Heading 6": 12,
}


def set_font_all_scripts(rpr_font_name, style):
    """Force Times New Roman for ascii/hAnsi/cs so nothing falls back."""
    rpr = style.element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.insert(0, rfonts)
    for attr in ("w:ascii", "w:hAnsi", "w:cs"):
        rfonts.set(qn(attr), rpr_font_name)


def style_paragraph(style, size_pt, indent, keep_alignment=False):
    f = style.font
    f.name = BODY_FONT
    f.size = Pt(size_pt)
    set_font_all_scripts(BODY_FONT, style)
    pf = style.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    if not keep_alignment:
        pf.alignment = WD_ALIGN_PARAGRAPH.LEFT
    pf.first_line_indent = indent if indent else Pt(0)


def ensure_body_sectPr(document):
    """Pandoc's default reference doc has an empty body with NO section
    properties, so python-docx reports zero sections and cannot attach a
    header/margins/line-numbering. Inject a US-Letter sectPr so there is a
    section to configure; Pandoc copies this sectPr (and the header part) into
    the rendered manuscript."""
    body = document.element.body
    if body.find(qn("w:sectPr")) is not None:
        return
    sectPr = OxmlElement("w:sectPr")
    pgSz = OxmlElement("w:pgSz")
    pgSz.set(qn("w:w"), "12240")   # 8.5" US Letter
    pgSz.set(qn("w:h"), "15840")   # 11"
    sectPr.append(pgSz)
    pgMar = OxmlElement("w:pgMar")
    for side, val in (("top", "1440"), ("right", "1440"), ("bottom", "1440"),
                      ("left", "1440"), ("header", "720"), ("footer", "720"),
                      ("gutter", "0")):
        pgMar.set(qn("w:" + side), val)
    sectPr.append(pgMar)
    body.append(sectPr)


def add_page_number_header(section):
    """Right-aligned { PAGE } field in the top header, Times New Roman 12 pt."""
    header = section.header
    header.is_linked_to_previous = False
    p = header.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE

    def run_with_font():
        r = OxmlElement("w:r")
        rpr = OxmlElement("w:rPr")
        rfonts = OxmlElement("w:rFonts")
        for a in ("w:ascii", "w:hAnsi", "w:cs"):
            rfonts.set(qn(a), BODY_FONT)
        sz = OxmlElement("w:sz")
        sz.set(qn("w:val"), str(BODY_PT * 2))
        rpr.append(rfonts)
        rpr.append(sz)
        r.append(rpr)
        return r, rpr

    # field: begin
    r_begin, _ = run_with_font()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    r_begin.append(fld_begin)
    # field: instruction " PAGE "
    r_instr, _ = run_with_font()
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    r_instr.append(instr)
    # field: end
    r_end, _ = run_with_font()
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    r_end.append(fld_end)

    for r in (r_begin, r_instr, r_end):
        p._p.append(r)


def add_line_numbering(section):
    """Continuous line numbering across the whole manuscript."""
    sectPr = section._sectPr
    existing = sectPr.find(qn("w:lnNumType"))
    if existing is not None:
        sectPr.remove(existing)
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")  # 0.25" gutter between numbers and text
    # w:lnNumType must precede certain trailing sectPr children; append is fine
    # for Word, which tolerates ordering here.
    sectPr.append(ln)


def define_table_style(document):
    """A 'Table' style with horizontal rules only — Pandoc applies it to every
    pipe table. Top/bottom of table + header underline; no vertical lines."""
    styles_el = document.styles.element
    # Remove any pre-existing Table style to avoid duplicates on re-run.
    for st in styles_el.findall(qn("w:style")):
        if st.get(qn("w:styleId")) == "Table":
            styles_el.remove(st)

    style = OxmlElement("w:style")
    style.set(qn("w:type"), "table")
    style.set(qn("w:styleId"), "Table")

    name = OxmlElement("w:name")
    name.set(qn("w:val"), "Table")
    style.append(name)
    based = OxmlElement("w:basedOn")
    based.set(qn("w:val"), "TableNormal")
    style.append(based)

    # Table-level run properties: Times New Roman 12 pt.
    rpr = OxmlElement("w:rPr")
    rfonts = OxmlElement("w:rFonts")
    for a in ("w:ascii", "w:hAnsi", "w:cs"):
        rfonts.set(qn(a), BODY_FONT)
    rpr.append(rfonts)
    sz = OxmlElement("w:sz")
    sz.set(qn("w:val"), str(BODY_PT * 2))
    rpr.append(sz)
    style.append(rpr)

    def border(kind, sz="4", val="single"):
        b = OxmlElement("w:" + kind)
        b.set(qn("w:val"), val)
        b.set(qn("w:sz"), sz)
        b.set(qn("w:space"), "0")
        b.set(qn("w:color"), "000000")
        return b

    tblPr = OxmlElement("w:tblPr")
    borders = OxmlElement("w:tblBorders")
    borders.append(border("top", sz="8"))
    borders.append(border("bottom", sz="8"))
    borders.append(border("insideH", sz="4"))
    borders.append(border("left", val="none", sz="0"))
    borders.append(border("right", val="none", sz="0"))
    borders.append(border("insideV", val="none", sz="0"))
    tblPr.append(borders)
    # Cell padding so text doesn't crowd the rules.
    margins = OxmlElement("w:tblCellMar")
    for side, w in (("top", "40"), ("bottom", "40"), ("left", "108"), ("right", "108")):
        m = OxmlElement("w:" + side)
        m.set(qn("w:w"), w)
        m.set(qn("w:type"), "dxa")
        margins.append(m)
    tblPr.append(margins)
    style.append(tblPr)

    # Header row: bold + bottom rule.
    tblStylePr = OxmlElement("w:tblStylePr")
    tblStylePr.set(qn("w:type"), "firstRow")
    h_rpr = OxmlElement("w:rPr")
    b = OxmlElement("w:b")
    h_rpr.append(b)
    tblStylePr.append(h_rpr)
    h_tcPr = OxmlElement("w:tcPr")
    h_borders = OxmlElement("w:tcBorders")
    h_borders.append(border("bottom", sz="8"))
    h_tcPr.append(h_borders)
    tblStylePr.append(h_tcPr)
    style.append(tblStylePr)

    styles_el.append(style)


def main():
    with open(DEFAULT_REF, "wb") as fh:
        subprocess.run(
            ["pandoc", "--print-default-data-file", "reference.docx"],
            check=True, stdout=fh,
        )

    doc = Document(DEFAULT_REF)

    by_name = {s.name: s for s in doc.styles if s.type is not None}
    for name, style in by_name.items():
        if name in INDENT_STYLES:
            style_paragraph(style, BODY_PT, INDENT)
        elif name in NO_INDENT_STYLES:
            size = HEADING_PT.get(name, BODY_PT)
            style_paragraph(style, size, None,
                            keep_alignment=(name in HEADING_PT and name != "Title"))

    # Title centered on the title page; headings stay left-aligned & bold.
    if "Title" in by_name:
        by_name["Title"].paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for h in ("Heading 1", "Heading 2", "Heading 3", "Heading 4", "Heading 5", "Heading 6"):
        if h in by_name:
            by_name[h].font.bold = True
            by_name[h].font.color.rgb = RGBColor(0, 0, 0)

    # Page geometry, header, line numbers on every section.
    ensure_body_sectPr(doc)
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
        section.header_distance = Inches(0.5)
        section.footer_distance = Inches(0.5)
        add_page_number_header(section)
        add_line_numbering(section)

    define_table_style(doc)

    doc.save(OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
