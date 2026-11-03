"""Post-render step for the manuscript docx: give the data tables horizontal rules.

Taken from the eluring manuscripts (adult_centenarian/reports/assets/docx_table_rules.py)
so this paper's Word copy matches theirs. Pandoc writes every docx table, figure wrappers
included, with the single "Table" style, whose only border is the line under the header row.
Each table with a header row (tblLook firstRow="1", which only data tables carry) gets a top
and a bottom rule and light grey lines between rows, no verticals. The paragraph right after
a data table gets the reference doc's 9 pt "Table Note" style.

Standalone, not a Quarto post-render hook: a _quarto.yml in manuscript/ would re-root
here::here() to this directory and break every here("results_...") path in _paths.R.

Usage, after `quarto render manuscript_short.qmd --to docx`:
    python3 docx_table_rules.py manuscript_short.docx
"""
import os
import re
import shutil
import tempfile
import zipfile
import xml.etree.ElementTree as ET

BORDERS = (
    '<w:tblBorders>'
    '<w:top w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
    '<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
    '<w:insideH w:val="single" w:sz="2" w:space="0" w:color="BFBFBF"/>'
    '</w:tblBorders>'
)
# tblBorders must precede these in tblPr (OOXML element order)
LATER = re.compile(r'<w:(shd|tblLayout|tblCellMar|tblLook)\b')
TBLPR = re.compile(r'<w:tblPr>.*?</w:tblPr>', re.S)


def add_rules(tblpr):
    if 'w:firstRow="1"' not in tblpr or '<w:tblBorders>' in tblpr:
        return tblpr
    m = LATER.search(tblpr)
    return tblpr[:m.start()] + BORDERS + tblpr[m.start():]


NOTE_PPR = '<w:pPr><w:pStyle w:val="TableNote" /></w:pPr>'
PPR = re.compile(r'<w:pPr>.*?</w:pPr>|<w:pPr\s*/>', re.S)


def style_notes(xml):
    """The paragraph right after a data table, inside its float cell, is the table note:
    give it the reference doc's 9 pt, left-aligned "Table Note" style (pandoc centres it)."""
    out, pos, n = [], 0, 0
    for m in TBLPR.finditer(xml):
        if 'w:firstRow="1"' not in m.group(0) or m.start() < pos:
            continue
        end = xml.find('</w:tbl>', m.end()) + len('</w:tbl>')
        if not xml.startswith('<w:p>', end):
            continue
        pend = xml.find('</w:p>', end)
        para = xml[end:pend]
        if PPR.match(para, len('<w:p>')):
            para = PPR.sub(NOTE_PPR, para, count=1)
        else:
            para = '<w:p>' + NOTE_PPR + para[len('<w:p>'):]
        out += [xml[pos:end], para]
        pos, n = pend, n + 1
    out.append(xml[pos:])
    return ''.join(out), n


CELL_SZ = 16   # half-points, so 8 pt
SZ = '<w:sz w:val="%d"/><w:szCs w:val="%d"/>' % (CELL_SZ, CELL_SZ)
TBL = re.compile(r'<w:tbl>.*?</w:tbl>', re.S)
# a run's opening tag together with its run-properties element, if it has one
RUN = re.compile(r'<w:r>(?:\s*(?:<w:rPr>(.*?)</w:rPr>|<w:rPr\s*/>))?', re.S)
# the rPr children OOXML requires to come after w:sz and w:szCs (CT_RPr order)
AFTER_SZ = re.compile(r'<w:(highlight|u|effect|bdr|shd|fitText|vertAlign|rtl|cs|em|lang'
                      r'|eastAsianLayout|specVanish|oMath)\b')


def set_size(inner):
    """Return a complete <w:rPr> carrying the 8 pt size, in schema order."""
    inner = re.sub(r'<w:sz\b[^>]*/>|<w:szCs\b[^>]*/>', '', inner)
    m = AFTER_SZ.search(inner)
    cut = m.start() if m else len(inner)
    return '<w:rPr>' + inner[:cut] + SZ + inner[cut:] + '</w:rPr>'


def shrink_cells(xml):
    """Set the text inside each data table to 8 pt.

    TP comment 101 on draft k: "use smaller font size in table cells to fit text,
    currently text is wrapped. also in other tables." The type catalogs are the widest,
    but every data table gains from it, so the size is set on each run inside a table
    rather than per table. It is a run property, not a style change, so the 9 pt table
    note and the 12 pt body text are untouched, and docx_jid_submit.py's 11-to-12 pt
    bump works on styles.xml and leaves these alone.

    Each run gets exactly one <w:rPr>: an existing one is rewritten in place (dropping
    any size it already carried), a run without one gets a new one. The first version of
    this function opened a <w:rPr> without closing it, which left document.xml unparseable
    and Word refusing to open draft l; the XML is checked after every run now.
    """
    def one_table(m):
        tbl = m.group(0)
        if 'w:firstRow="1"' not in tbl:
            return tbl
        return RUN.sub(lambda r: '<w:r>' + set_size(r.group(1) or ''), tbl)
    out = TBL.sub(one_table, xml)
    return out, sum(1 for a, b in zip(TBL.findall(xml), TBL.findall(out)) if a != b)


def process(path):
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.docx',
                                      dir=os.path.dirname(path))
    tmp.close()
    n = notes = shrunk = 0
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp.name, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == 'word/document.xml':
                xml = data.decode('utf-8')
                xml, notes = style_notes(xml)
                new = TBLPR.sub(lambda m: add_rules(m.group(0)), xml)
                n = sum(1 for a, b in zip(TBLPR.findall(xml), TBLPR.findall(new)) if a != b)
                new, shrunk = shrink_cells(new)
                # these edits are regex surgery on OOXML: refuse to write a
                # document.xml Word cannot parse (see shrink_cells' docstring)
                try:
                    ET.fromstring(new)
                except ET.ParseError as e:
                    raise SystemExit('docx_table_rules: edits broke document.xml '
                                     '(%s); %s left unchanged' % (e, path))
                data = new.encode('utf-8')
            zout.writestr(item, data)
    shutil.move(tmp.name, path)
    print(f'docx_table_rules: {n} tables ruled, {notes} notes styled, '
          f'{shrunk} tables set to 8 pt in {path}')


if __name__ == '__main__':
    import sys
    for f in sys.argv[1:]:
        process(f)
