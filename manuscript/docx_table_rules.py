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


def inner_table(tbl):
    """Split pandoc's float wrapper from the data table inside it.

    Pandoc puts each table and its caption in a one-cell outer <w:tbl>, so a
    non-greedy <w:tbl>...</w:tbl> match runs from the wrapper's opening tag to
    the inner table's closing one. Everything before the last <w:tbl> is
    wrapper and caption, and must not be restyled as table content.
    """
    i = tbl.rfind('<w:tbl>')
    return (tbl[:i], tbl[i:]) if i > 0 else ('', tbl)


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
        head, data = inner_table(tbl)
        return head + RUN.sub(lambda r: '<w:r>' + set_size(r.group(1) or ''), data)
    out = TBL.sub(one_table, xml)
    return out, sum(1 for a, b in zip(TBL.findall(xml), TBL.findall(out)) if a != b)


# Text width of the page, in twips: 8.5 in paper less two 1 in margins.
TEXT_TW = 9360
# Width of one character of 8 pt Times New Roman, in twips, and the default
# cell margin on each side. Deliberately generous: the cost of overestimating
# is a slightly wide column, the cost of underestimating is a wrapped word.
CHAR_TW = 95
CELL_MARGIN_TW = 108
GRIDCOL = re.compile(r'<w:gridCol w:w="(\d+)"\s*/>')
TC = re.compile(r'<w:tc>.*?</w:tc>', re.S)
TR = re.compile(r'<w:tr\b.*?</w:tr>', re.S)
TCW = re.compile(r'<w:tcW w:w="\d+" w:type="dxa"\s*/>')


def _longest_word_tw(cell_xml):
    """Twips needed by the longest unbreakable word in a cell."""
    # <w:t[^>]*> would also match <w:tc>, <w:tcPr/> and <w:tbl>, pulling markup
    # in as text; require a space or an immediate close after the element name
    text = ''.join(re.findall(r'<w:t(?: [^>]*)?>(.*?)</w:t>', cell_xml, re.S))
    text = re.sub(r'&[a-z]+;|&#\d+;', 'x', text)        # entities are one glyph
    words = text.split()
    return (max((len(w) for w in words), default=0) * CHAR_TW) + 2 * CELL_MARGIN_TW


def widen_columns(xml):
    """Give every column enough width for its longest word to fit on one line.

    TP comment 106 on draft o: "check column width, hpv type must fit into one
    line." Pandoc derives column widths from the width of the source markdown
    columns, so a table whose first column holds long group labels starves the
    rest: in Table S1 the first column took 51% of the table and the type
    columns were left with as little as 2%, which wraps HPV18 onto two lines.

    Each column is given max(its longest word, its pandoc share). The slack is
    handed back in proportion to the original widths, so a column that genuinely
    holds prose stays wide. If the minimums do not fit the page the whole set is
    scaled down, which is no worse than what pandoc wrote.
    """
    def one_table(whole):
        if 'w:firstRow="1"' not in whole:
            return whole, False
        head, tbl = inner_table(whole)
        grid = [int(g) for g in GRIDCOL.findall(tbl)]
        if not grid:
            return whole, False
        rows = TR.findall(tbl)
        if not rows or any(len(TC.findall(r)) != len(grid) for r in rows):
            return whole, False          # merged or spanned cells: leave alone
        need = [0] * len(grid)
        for r in rows:
            for i, c in enumerate(TC.findall(r)):
                need[i] = max(need[i], _longest_word_tw(c))
        total = sum(grid)
        rendered = [g * TEXT_TW / total for g in grid]
        if sum(need) >= TEXT_TW:                       # cannot satisfy: scale
            new = [n * TEXT_TW / sum(need) for n in need]
        else:
            base = [max(n, 0) for n in need]
            slack = TEXT_TW - sum(base)
            new = [b + slack * r / TEXT_TW for b, r in zip(base, rendered)]
        new = [max(1, round(v * total / TEXT_TW)) for v in new]   # back to grid units
        if new == grid:
            return whole, False
        it = iter(new)
        out = GRIDCOL.sub(lambda m: '<w:gridCol w:w="%d"/>' % next(it), tbl)
        # per-cell tcW must agree with the grid, or Word prefers the cell value
        def one_row(m):
            k = iter(range(len(new)))
            return TC.sub(lambda c: TCW.sub(
                '<w:tcW w:w="%d" w:type="dxa"/>' % new[next(k)], c.group(0), count=1),
                m.group(0))
        return head + TR.sub(one_row, out), True

    out, n = [], 0
    pos = 0
    for m in TBL.finditer(xml):
        tbl, changed = one_table(m.group(0))
        out += [xml[pos:m.start()], tbl]
        pos = m.end()
        n += changed
    out.append(xml[pos:])
    return ''.join(out), n


def process(path):
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.docx',
                                      dir=os.path.dirname(path))
    tmp.close()
    n = notes = shrunk = widened = 0
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp.name, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == 'word/document.xml':
                xml = data.decode('utf-8')
                xml, notes = style_notes(xml)
                new = TBLPR.sub(lambda m: add_rules(m.group(0)), xml)
                n = sum(1 for a, b in zip(TBLPR.findall(xml), TBLPR.findall(new)) if a != b)
                new, shrunk = shrink_cells(new)
                new, widened = widen_columns(new)
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
          f'{shrunk} tables set to 8 pt, {widened} rebalanced in {path}')


if __name__ == '__main__':
    import sys
    for f in sys.argv[1:]:
        process(f)
