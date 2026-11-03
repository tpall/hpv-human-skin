"""Turn the co-author (eluring-style) docx into the JID submission format.

The co-author drafts use the eluring reference.docx: Times New Roman 11 pt, double spacing,
continuous line numbering. JID's Information for Authors ("How to Format Your Manuscript",
read 2026-11-02) asks instead for Times New Roman 12 pt, double spacing, 1-inch margins,
page numbers starting with the title page as page 1, and NO line numbers. This script makes
exactly those changes and nothing else:

  - every 11 pt size in styles.xml becomes 12 pt (the 9 pt Table Note style is left alone);
  - every <w:lnNumType> is removed from the section properties;
  - a footer with a centred PAGE field is added to every section;
  - heading styles are set to JID's three levels ("Heading Styles", same page): first
    level ALL CAPS BOLD, second level Sentence Case Bold Underline, third level Sentence
    Case Bold Italics. All four pandoc heading styles are already bold and Heading4 is
    already bold italic, so this adds <w:caps/> to Heading1 and Heading2 and an underline
    to Heading3. Caps are applied as a run property rather than by uppercasing the text,
    so the heading text in the file stays as written and an editor can undo the look
    without retyping anything;
  - a page break is inserted before the References heading and before the supplementary
    material, as JID asks ("Begin a new page for References, Tables, and Figure Legends").
    Tables and figure legends sit inline in this manuscript rather than in sections of
    their own, so there is nothing for this script to break before; moving them to the
    end is a decision about the document, not about its formatting.

Margins are already 1 inch and spacing already double in the reference doc.

Which pandoc heading level maps to which JID level depends on this manuscript: Heading2
carries the main sections (Abstract, Introduction, ...), Heading3 the Results and
supplementary subsections, and Heading1 only the "Supplementary Material" divider. Both
Heading1 and Heading2 are therefore first-level headings in JID's sense.

Usage:
    python3 docx_jid_submit.py hpv_manuscript_DATE_jid.docx hpv_manuscript_DATE_jid_submit.docx
"""
import re
import sys
import zipfile

FOOTER = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
    '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate"/></w:r>'
    '<w:r><w:t>1</w:t></w:r>'
    '<w:r><w:fldChar w:fldCharType="end"/></w:r>'
    '</w:p></w:ftr>'
)
RID = 'rIdJidPageFooter'
REL = ('<Relationship Id="%s" '
       'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
       'Target="footer_jid.xml"/>' % RID)
CT = ('<Override PartName="/word/footer_jid.xml" '
      'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/>')


PAGE_BREAK = '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'
# Headings that must start on a new page. Matched on the heading's own text, because the
# References heading is emitted by citeproc and carries no identifier of its own.
BREAK_BEFORE = ('References', 'Tables', 'Figure Legends', 'Supplementary Material')


def add_page_breaks(xml):
    """Put a page break before each heading named in BREAK_BEFORE."""
    out, n = [], 0
    for para in re.split(r'(?=<w:p[ >])', xml):
        m = re.match(r'<w:p[ >]', para)
        if m and re.search(r'<w:pStyle w:val="Heading\d"', para):
            text = re.sub(r'<[^>]+>', '', re.sub(r'</w:p>.*', '', para, flags=re.S)).strip()
            if text in BREAK_BEFORE:
                out.append(PAGE_BREAK)
                n += 1
        out.append(para)
    return ''.join(out), n


def fix_document(xml):
    xml = re.sub(r'<w:lnNumType\b[^>]*/>', '', xml)
    # existing footer references would conflict; drop them, then add ours as the first
    # child of each sectPr (footerReference must precede pgSz in OOXML order)
    xml = re.sub(r'<w:footerReference\b[^>]*/>', '', xml)
    ref = '<w:footerReference w:type="default" r:id="%s"/>' % RID
    xml = re.sub(r'(<w:sectPr\b[^>]*>)', r'\1' + ref, xml)
    if 'xmlns:r=' not in xml[:3000]:
        xml = xml.replace('<w:document ', '<w:document xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" ', 1)
    return xml


def fix_styles(xml):
    # pandoc writes '<w:sz w:val="22" />' (space before the slash), Word writes it without
    xml = re.sub(r'<w:sz w:val="22"\s*/>', '<w:sz w:val="24"/>', xml)
    xml = re.sub(r'<w:szCs w:val="22"\s*/>', '<w:szCs w:val="24"/>', xml)
    return fix_heading_styles(xml)


# JID's three heading levels, as run properties to add to each style's <w:rPr>.
HEADING_RPR = {'Heading1': '<w:caps/>', 'Heading2': '<w:caps/>',
               'Heading3': '<w:u w:val="single"/>', 'Heading4': ''}


def fix_heading_styles(xml):
    for style, extra in HEADING_RPR.items():
        if not extra:
            continue          # Heading4 is already bold italic, which is JID's third level
        pat = re.compile(r'(<w:style [^>]*w:styleId="%s".*?<w:rPr>)' % style, re.S)
        if not pat.search(xml):
            raise SystemExit('docx_jid_submit: no rPr found for style %s' % style)
        if extra in pat.search(xml).group(1):
            continue          # already applied
        xml = pat.sub(lambda m: m.group(1) + extra, xml, count=1)
    return xml


def main(src, dst):
    breaks = 0
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == 'word/document.xml':
                xml, breaks = add_page_breaks(data.decode('utf-8'))
                data = fix_document(xml).encode('utf-8')
            elif item.filename == 'word/styles.xml':
                data = fix_styles(data.decode('utf-8')).encode('utf-8')
            elif item.filename == 'word/_rels/document.xml.rels':
                data = data.decode('utf-8').replace('</Relationships>', REL + '</Relationships>').encode('utf-8')
            elif item.filename == '[Content_Types].xml':
                data = data.decode('utf-8').replace('</Types>', CT + '</Types>').encode('utf-8')
            zout.writestr(item, data)
        zout.writestr('word/footer_jid.xml', FOOTER)
    print('docx_jid_submit: 12 pt, no line numbers, page-number footer, JID heading '
          'styles, %d page break(s) -> %s' % (breaks, dst))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
