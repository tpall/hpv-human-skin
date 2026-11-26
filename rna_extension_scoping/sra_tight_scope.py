#!/usr/bin/env python3
# Tightened, curation-aware SRA scoping for the two HPV-relevant RNA niches we
# lack: vaginal/cervical and anorectal. For each, count (a) the tissue-targeted
# RNA-seq set and (b) the fraction that trips cell-line / in-vitro indicator
# terms in metadata -- an upper-bound contamination burden BEFORE any download.
import urllib.parse, urllib.request, re, time, sys
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

def count(term):
    q = urllib.parse.urlencode({"db": "sra", "term": term, "retmax": 0})
    for _ in range(4):
        try:
            with urllib.request.urlopen(EUTILS + "?" + q, timeout=30) as r:
                xml = r.read().decode()
            m = re.search(r"<Count>(\d+)</Count>", xml)
            if m: return int(m.group(1))
            if "<ERROR>" in xml: return "ERR"
        except Exception:
            time.sleep(1.0)
    return "FAIL"

HUM = '"Homo sapiens"[Organism]'
RNA = '"RNA-Seq"[Strategy]'
PUB = '"public"[Access]'
BASE = f'{HUM} AND {RNA} AND {PUB}'

# tissue-tight positive terms (Text Word / MeSH)
CERVIX = ('("cervix"[MeSH Terms] OR cervical[Text Word] OR cervix[Text Word] OR '
          'ectocervix[Text Word] OR endocervix[Text Word] OR vaginal[Text Word] OR '
          'vagina[Text Word] OR vulvar[Text Word] OR vulva[Text Word] OR '
          '"cervical intraepithelial"[Text Word] OR CIN[Text Word] OR condyloma[Text Word])')
ANORECTAL = ('(anal[Text Word] OR anorectal[Text Word] OR perianal[Text Word] OR '
             '"anal canal"[Text Word] OR "anal squamous"[Text Word] OR '
             '"anal intraepithelial"[Text Word] OR AIN[Text Word] OR '
             '"anal cancer"[Text Word] OR "anal carcinoma"[Text Word] OR '
             '"anal condyloma"[Text Word])')
# broad gut (for contrast with the anorectal subset)
GUT = ('(gut[Text Word] OR intestin*[Text Word] OR colon[Text Word] OR colonic[Text Word] OR '
       'colorectal[Text Word] OR rectal[Text Word] OR rectum[Text Word] OR ileum[Text Word])')

# cell-line / in-vitro indicator terms (metadata upper bound; mirrors cell_line_patterns.py)
CL = ('(HeLa[Text Word] OR SiHa[Text Word] OR CaSki[Text Word] OR "ME-180"[Text Word] OR '
      '"MS751"[Text Word] OR "C-33A"[Text Word] OR "C33A"[Text Word] OR SW756[Text Word] OR '
      'W12[Text Word] OR "cell line"[Text Word] OR "cell lines"[Text Word] OR '
      'immortalized[Text Word] OR immortalised[Text Word] OR ATCC[Text Word] OR '
      'xenograft[Text Word] OR organoid[Text Word] OR organoids[Text Word] OR '
      'transfected[Text Word] OR transduced[Text Word] OR CRISPR[Text Word] OR '
      'shRNA[Text Word] OR siRNA[Text Word] OR "iPSC"[Text Word])')

# vaginal / gut metatranscriptome (metagenome-labelled)
VMETA = ('("vaginal metagenome"[Organism] OR "human vaginal metagenome"[Organism]) AND '
         '("metatranscriptomic"[Source] OR "RNA-Seq"[Strategy]) AND ' + PUB)
GMETA = ('("human gut metagenome"[Organism] OR "gut metagenome"[Organism]) AND '
         '("metatranscriptomic"[Source] OR "RNA-Seq"[Strategy]) AND ' + PUB)

rows = [
    ("VAGINAL/CERVICAL tissue RNA-seq",        f'{BASE} AND {CERVIX}'),
    ("  ...of which trip cell-line terms",     f'{BASE} AND {CERVIX} AND {CL}'),
    ("  ...cell-line-EXCLUDED (NOT)",          f'{BASE} AND {CERVIX} NOT {CL}'),
    ("VAGINAL metatranscriptome",              VMETA),
    ("ANORECTAL tissue RNA-seq (tight)",       f'{BASE} AND {ANORECTAL}'),
    ("  ...of which trip cell-line terms",     f'{BASE} AND {ANORECTAL} AND {CL}'),
    ("  ...cell-line-EXCLUDED (NOT)",          f'{BASE} AND {ANORECTAL} NOT {CL}'),
    ("ANORECTAL metatranscriptome",            f'{GMETA} AND {ANORECTAL}'),
    ("-- context: broad GUT RNA-seq",          f'{BASE} AND {GUT}'),
    ("  ...of which trip cell-line terms",     f'{BASE} AND {GUT} AND {CL}'),
]
print(f"{'query':44s} {'runs':>9s}")
print("-" * 56)
for label, term in rows:
    print(f"{label:44s} {str(count(term)):>9s}"); sys.stdout.flush(); time.sleep(0.4)
