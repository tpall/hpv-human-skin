"""Phase 3: verify each matched publication really belongs to the study, then
re-extract the ethics text with a sentence splitter that does not break on
abbreviations such as 'St.'.

The decisive check is whether the paper's own full text mentions the study or
BioProject accession. A paper that never names the accession is not evidence for
that dataset's approval, however well the search matched.

Usage: python3 -I ethics_verify.py
Reads ethics_lookup2.tsv, writes ethics_verified.tsv
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "hpv-human-skin ethics audit (taavi.pall@ut.ee)"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"

PROT = [(r"\bSt\.", "St<D>"), (r"\bDr\.", "Dr<D>"), (r"\bNo\.", "No<D>"),
        (r"\bRef\.", "Ref<D>"), (r"\bInc\.", "Inc<D>"), (r"\bLtd\.", "Ltd<D>"),
        (r"\bvs\.", "vs<D>"), (r"\be\.g\.", "e<D>g<D>"), (r"\bi\.e\.", "i<D>e<D>"),
        (r"\bet al\.", "et al<D>"), (r"\bU\.S\.", "U<D>S<D>"),
        (r"\bNo\s*\.", "No<D>"), (r"\bapprox\.", "approx<D>")]

# A committee statement, deliberately excluding the loose "approved by the",
# which matched drug-regulatory approvals in the first pass.
COMMITTEE = re.compile(
    r"\b(ethic\w*\s+(committee|board|review|approval|permission|statement)"
    r"|institutional\s+review\s+board|\bIRB\b|research\s+ethics"
    r"|ethical\s+(approval|committee|permission|review)"
    r"|institutional\s+animal\s+care)\b", re.I)
CONSENT = re.compile(r"\b(informed\s+consent|written\s+consent|consent\s+was"
                     r"|consent\s+to\s+participate|Helsinki|waiver)\b", re.I)
# Approvals that are not research-ethics approvals at all.
NOT_ETHICS = re.compile(
    r"Food\s+and\s+Drug\s+Administration|\bFDA\b|European\s+Medicines"
    r"|reimbursement|marketing\s+authoris", re.I)


def get(url):
    for a in range(3):
        try:
            with urllib.request.urlopen(
                    urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if a == 2:
                return ""
            time.sleep(2)
    return ""


def plain(xml):
    t = xml
    for pat, rep in PROT:
        t = re.sub(pat, rep, t)
    t = re.sub(r"<[^>]+>", " ", t)
    t = re.sub(r"&[a-z]+;", " ", t)
    t = re.sub(r"https?://\S+", " ", t)
    return re.sub(r"\s+", " ", t)


def pick(text, rx, limit=3):
    out, seen = [], set()
    for s in re.split(r"(?<=\.)\s+", text):
        s = s.strip()
        if not (40 < len(s) < 520):
            continue
        if rx.search(s) and not NOT_ETHICS.search(s):
            k = s.lower()
            if k not in seen:
                seen.add(k)
                out.append(s.replace("<D>", "."))
        if len(out) >= limit:
            break
    return out


def main():
    rows = list(csv.DictReader(open("ethics_lookup2.tsv"), delimiter="\t"))
    w = csv.writer(open("ethics_verified.tsv", "w"), delimiter="\t")
    w.writerow(["study", "bioproject", "n_libs", "diagnoses", "pmid", "pmcid",
                "accession_in_fulltext", "committee", "consent", "flag"])
    for i, r in enumerate(rows, 1):
        pmcid, flag = r["pmcid"], ""
        acc_ok, comm, cons = "no_fulltext", [], []
        if pmcid and pmcid != "NOT_FOUND":
            xml = get(EPMC + f"{pmcid}/fullTextXML")
            if xml and "<" in xml:
                txt = plain(xml)
                hits = [a for a in (r["study"], r["bioproject"])
                        if a != "NOT_FOUND" and a in txt]
                acc_ok = ",".join(hits) if hits else "ACCESSION ABSENT"
                comm = pick(txt, COMMITTEE)
                cons = pick(txt, CONSENT)
        if acc_ok == "ACCESSION ABSENT":
            flag = "VERIFY: paper never names this accession"
        elif acc_ok == "no_fulltext" and r["pmid"] != "NOT_FOUND":
            flag = "VERIFY: no open full text, check manually"
        elif r["pmid"] == "NOT_FOUND":
            flag = "no publication"
        if not comm and r["pmid"] != "NOT_FOUND" and not flag:
            flag = "VERIFY: no committee sentence found"
        w.writerow([r["study"], r["bioproject"], r["n_libs"], r["diagnoses"],
                    r["pmid"], pmcid, acc_ok,
                    " | ".join(comm) or "NONE", " | ".join(cons) or "NONE", flag])
        print(f"[{i}/{len(rows)}] {r['study']:12} acc={acc_ok[:28]:28} "
              f"comm={len(comm)} {flag}", flush=True)
        time.sleep(0.3)


if __name__ == "__main__":
    main()
