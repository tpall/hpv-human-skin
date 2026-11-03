"""Phase 2: for studies NCBI could not link to a paper, search Europe PMC, which
indexes data accessions in full text. Then pull ethics sentences from the open
access full text.

Nothing is inferred. Each row is either quoted from a record or NOT_FOUND.

Usage: python3 -I ethics_lookup2.py
Reads ethics_lookup.tsv, writes ethics_lookup2.tsv.
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
UA = {"User-Agent": "hpv-human-skin ethics audit (taavi.pall@ut.ee)"}

# Abbreviations that must not end a sentence, so "St. Louis" stays intact.
PROT = [(r"\bSt\.", "St<DOT>"), (r"\bDr\.", "Dr<DOT>"), (r"\bNo\.", "No<DOT>"),
        (r"\bRef\.", "Ref<DOT>"), (r"\bApprox\.", "Approx<DOT>"),
        (r"\bInc\.", "Inc<DOT>"), (r"\bLtd\.", "Ltd<DOT>"), (r"\bvs\.", "vs<DOT>"),
        (r"\be\.g\.", "e<DOT>g<DOT>"), (r"\bi\.e\.", "i<DOT>e<DOT>"),
        (r"\bet al\.", "et al<DOT>"), (r"\bU\.S\.", "U<DOT>S<DOT>")]

APPROVAL = re.compile(
    r"[^.]*?\b("
    r"ethic\w*\s+(committee|board|approval|review|statement)|"
    r"institutional\s+review\s+board|\bIRB\b|"
    r"research\s+ethics|ethical\s+(approval|committee|permission)|"
    r"approved\s+by\s+the|Helsinki|informed\s+consent|written\s+consent|"
    r"consent\s+was\s+obtained|consent\s+to\s+participate|waiver\s+of"
    r")\b[^.]*\.", re.I)


def get(url):
    for a in range(4):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if a == 3:
                return ""
            time.sleep(2 * (a + 1))
    return ""


def sentences(txt):
    for pat, rep in PROT:
        txt = re.sub(pat, rep, txt)
    txt = re.sub(r"<[^>]+>", " ", txt)
    txt = re.sub(r"&[a-z]+;", " ", txt)
    txt = re.sub(r"https?://\S+", " ", txt)
    txt = re.sub(r"\s+", " ", txt)
    out, seen = [], set()
    for m in APPROVAL.finditer(txt):
        s = m.group(0).replace("<DOT>", ".").strip()
        if 40 < len(s) < 600 and s.lower() not in seen:
            seen.add(s.lower())
            out.append(s)
    return out[:6]


def epmc_search(term):
    u = EPMC + "search?" + urllib.parse.urlencode(
        {"query": term, "format": "json", "pageSize": 10,
         "resultType": "core"})
    try:
        return json.loads(get(u)).get("resultList", {}).get("result", [])
    except Exception:
        return []


def epmc_fulltext(pmcid):
    return get(EPMC + f"{pmcid}/fullTextXML")


def main():
    rows = list(csv.DictReader(open("ethics_lookup.tsv"), delimiter="\t"))
    w = csv.writer(open("ethics_lookup2.tsv", "w"), delimiter="\t")
    w.writerow(["study", "n_libs", "diagnoses", "bioproject", "pmid", "pmcid",
                "source", "ethics_sentences"])
    for i, r in enumerate(rows, 1):
        pmid, pmcid = r["pmid"], r["pmcid"]
        sents = [s for s in r["ethics_sentences"].split(" || ")
                 if s and s != "NOT_FOUND"]
        source = "ncbi_elink" if pmid != "NOT_FOUND" else "NOT_FOUND"

        # Search Europe PMC when no paper, or when extraction produced nothing useful.
        need = pmid == "NOT_FOUND" or not sents
        if need and r["bioproject"] != "NOT_FOUND":
            for term in (f'"{r["bioproject"]}"', f'"{r["study"]}"'):
                res = epmc_search(term)
                if res:
                    hit = res[0]
                    pmid = hit.get("pmid") or pmid
                    pmcid = hit.get("pmcid") or pmcid
                    source = "europepmc:" + term.strip('"')
                    break
                time.sleep(0.3)

        # Pull full text if we have a PMC id and still no sentences.
        if pmcid != "NOT_FOUND" and not sents:
            xml = epmc_fulltext(pmcid)
            if xml:
                sents = sentences(xml)
                if sents and source == "NOT_FOUND":
                    source = "europepmc_fulltext"

        w.writerow([r["study"], r["n_libs"], r["diagnoses"], r["bioproject"],
                    pmid, pmcid, source,
                    " || ".join(sents) if sents else "NOT_FOUND"])
        print(f"[{i}/{len(rows)}] {r['study']} pmid={pmid} pmc={pmcid} "
              f"src={source} sent={len(sents)}", flush=True)
        time.sleep(0.3)


if __name__ == "__main__":
    main()
