"""Resolve SRA study accessions to BioProject, linked publication and (where open
access) the verbatim ethics/consent sentences from PMC full text.

Writes ethics_lookup.tsv. Nothing is inferred: a field is either quoted from the
record or left as NOT_FOUND, so the gaps are visible.

Usage: python3 -I ethics_lookup.py ethics_studies.csv
"""
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "hpv-human-skin ethics audit (taavi.pall@ut.ee)"}


def get(url):
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if attempt == 3:
                return ""
            time.sleep(2 * (attempt + 1))
    return ""


def esearch(db, term):
    u = E + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmode": "json", "retmax": 5})
    try:
        return json.loads(get(u))["esearchresult"].get("idlist", [])
    except Exception:
        return []


def elink(dbfrom, db, uid):
    u = E + "elink.fcgi?" + urllib.parse.urlencode(
        {"dbfrom": dbfrom, "db": db, "id": uid, "retmode": "json"})
    try:
        d = json.loads(get(u))
        out = []
        for ls in d["linksets"]:
            for db_ in ls.get("linksetdbs", []):
                out += db_.get("links", [])
        return out
    except Exception:
        return []


def esummary(db, uid):
    u = E + "esummary.fcgi?" + urllib.parse.urlencode(
        {"db": db, "id": uid, "retmode": "json"})
    try:
        d = json.loads(get(u))["result"]
        return d[d["uids"][0]]
    except Exception:
        return {}


# Sentences that carry an approval or a consent statement.
APPROVAL = re.compile(
    r"[^.]*?\b("
    r"ethic\w*\s+(committee|board|approval|review)|"
    r"institutional\s+review\s+board|\bIRB\b|"
    r"research\s+ethics|ethical\s+(approval|committee|permission)|"
    r"approved\s+by|Helsinki|informed\s+consent|written\s+consent|"
    r"consent\s+was\s+obtained|waiver"
    r")\b[^.]*\.", re.I)


def pmc_ethics(pmcid):
    """Return verbatim sentences from the PMC full text, or []."""
    u = E + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "pmc", "id": pmcid, "retmode": "xml"})
    xml = get(u)
    if not xml or "<body" not in xml:
        return []
    txt = re.sub(r"<[^>]+>", " ", xml)
    txt = re.sub(r"&[a-z]+;", " ", txt)
    txt = re.sub(r"\s+", " ", txt)
    hits, seen = [], set()
    for m in APPROVAL.finditer(txt):
        s = m.group(0).strip()
        if 40 < len(s) < 600 and s.lower() not in seen:
            seen.add(s.lower())
            hits.append(s)
    return hits[:6]


def main(path):
    rows = list(csv.DictReader(open(path)))
    w = csv.writer(open("ethics_lookup.tsv", "w"), delimiter="\t")
    w.writerow(["study", "n_libs", "diagnoses", "bioproject", "bioproject_title",
                "pmid", "pmcid", "ethics_sentences"])
    for i, r in enumerate(rows, 1):
        srp = r["study"]
        bp = bptitle = pmid = pmcid = "NOT_FOUND"
        sents = []
        ids = esearch("sra", srp + "[All Fields]")
        bpids = elink("sra", "bioproject", ids[0]) if ids else []
        if bpids:
            s = esummary("bioproject", bpids[0])
            bp = s.get("project_acc") or "NOT_FOUND"
            bptitle = (s.get("project_title") or "")[:140] or "NOT_FOUND"
            pmids = elink("bioproject", "pubmed", bpids[0])
            if not pmids and bp != "NOT_FOUND":
                pmids = esearch("pubmed", bp + "[All Fields]")
            if pmids:
                pmid = pmids[0]
                pmcs = elink("pubmed", "pmc", pmid)
                if pmcs:
                    s2 = esummary("pmc", pmcs[0])
                    pmcid = "PMC" + str(s2.get("uid") or pmcs[0])
                    sents = pmc_ethics(pmcs[0])
        w.writerow([srp, r["n_libs"], r["diagnoses"], bp, bptitle, pmid, pmcid,
                    " || ".join(sents) if sents else "NOT_FOUND"])
        print(f"[{i}/{len(rows)}] {srp} bp={bp} pmid={pmid} pmc={pmcid} "
              f"sent={len(sents)}", flush=True)
        time.sleep(0.4)


if __name__ == "__main__":
    main(sys.argv[1])
