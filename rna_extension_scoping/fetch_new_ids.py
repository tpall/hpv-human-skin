#!/usr/bin/env python3
# Fetch run accessions for the curation-ready queries and diff against the
# 12,200 SRRs already screened by hpv-human-skin -> how many are NEW.
import urllib.parse, urllib.request, re, time, sys
ESEARCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH ="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def get_runs(term):
    # esearch with history
    q=urllib.parse.urlencode({"db":"sra","term":term,"retmax":0,"usehistory":"y"})
    xml=urllib.request.urlopen(ESEARCH+"?"+q,timeout=60).read().decode()
    cnt=int(re.search(r"<Count>(\d+)</Count>",xml).group(1))
    webenv=re.search(r"<WebEnv>(\S+?)</WebEnv>",xml).group(1)
    qk=re.search(r"<QueryKey>(\d+)</QueryKey>",xml).group(1)
    runs=set()
    step=500
    for start in range(0,cnt,step):
        for _ in range(4):
            try:
                p=urllib.parse.urlencode({"db":"sra","query_key":qk,"WebEnv":webenv,
                    "rettype":"runinfo","retmode":"text","retstart":start,"retmax":step})
                txt=urllib.request.urlopen(EFETCH+"?"+p,timeout=120).read().decode(errors="ignore")
                for line in txt.splitlines():
                    acc=line.split(",",1)[0]
                    if acc[:3] in ("SRR","ERR","DRR"): runs.add(acc)
                break
            except Exception:
                time.sleep(2)
        time.sleep(0.4)
        print(f"    ...{min(start+step,cnt)}/{cnt}", end="\r", file=sys.stderr)
    print(" "*40, end="\r", file=sys.stderr)
    return cnt, runs

existing=set(l.strip() for l in open("/private/tmp/claude-501/-Users-taavi-Projects-eluring/b54054ee-044c-44ec-8fd7-c64fb7c316ce/scratchpad/existing_srrs.txt") if l.strip())
print("existing screened SRRs:", len(existing))

HUM='"Homo sapiens"[Organism]'; RNA='"RNA-Seq"[Strategy]'; PUB='"public"[Access]'
BASE=f'{HUM} AND {RNA} AND {PUB}'
CERVIX=('("cervix"[MeSH Terms] OR cervical[Text Word] OR cervix[Text Word] OR ectocervix[Text Word] OR '
        'endocervix[Text Word] OR vaginal[Text Word] OR vagina[Text Word] OR vulvar[Text Word] OR '
        'vulva[Text Word] OR "cervical intraepithelial"[Text Word] OR CIN[Text Word] OR condyloma[Text Word])')
ANORECTAL=('(anal[Text Word] OR anorectal[Text Word] OR perianal[Text Word] OR "anal canal"[Text Word] OR '
           '"anal squamous"[Text Word] OR "anal intraepithelial"[Text Word] OR AIN[Text Word] OR '
           '"anal cancer"[Text Word] OR "anal carcinoma"[Text Word] OR "anal condyloma"[Text Word])')
CL=('(HeLa[Text Word] OR SiHa[Text Word] OR CaSki[Text Word] OR "ME-180"[Text Word] OR "MS751"[Text Word] OR '
    '"C-33A"[Text Word] OR "C33A"[Text Word] OR SW756[Text Word] OR W12[Text Word] OR "cell line"[Text Word] OR '
    '"cell lines"[Text Word] OR immortalized[Text Word] OR immortalised[Text Word] OR ATCC[Text Word] OR '
    'xenograft[Text Word] OR organoid[Text Word] OR organoids[Text Word] OR transfected[Text Word] OR '
    'transduced[Text Word] OR CRISPR[Text Word] OR shRNA[Text Word] OR siRNA[Text Word] OR iPSC[Text Word])')

for label, term in [
    ("VAGINAL/CERVICAL (cell-line-excluded)", f'{BASE} AND {CERVIX} NOT {CL}'),
    ("ANORECTAL (cell-line-excluded)",        f'{BASE} AND {ANORECTAL} NOT {CL}'),
]:
    cnt, runs = get_runs(term)
    new = runs - existing
    have = runs & existing
    print(f"\n== {label} ==")
    print(f"  query count (esearch): {cnt}")
    print(f"  run accessions fetched: {len(runs)}")
    print(f"  already screened:       {len(have)}")
    print(f"  NEW (not yet screened): {len(new)}")
    fn="/private/tmp/claude-501/-Users-taavi-Projects-eluring/b54054ee-044c-44ec-8fd7-c64fb7c316ce/scratchpad/new_"+re.sub(r'\W+','_',label.split()[0]).lower()+".txt"
    open(fn,"w").write("\n".join(sorted(new)))
    print("  ->", fn)
