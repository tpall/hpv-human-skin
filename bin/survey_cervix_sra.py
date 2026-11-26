#!/usr/bin/env python3
"""Is there a PUBLIC, documented primary-tissue cervical RNA-seq cohort in SRA?

Motivation: the deck's Part III promises cervix as a "clean opposite to skin"
(presence without activity in skin vs presence with activity in cervix). The
vaginal baseline we had turned out 78% integrated cell-line material and a
random-1000 sample 39%, so that promise is not yet backed. TCGA-CESC cannot
back it either: all 2,512 read-level BAMs are dbGaP-controlled, and it holds
only 8 solid-tissue normals against 307 primary tumours.

So: survey what IS public. Pull every public human cervical RNA-seq run's full
metadata, then classify it with the SAME cell-line patterns and tissue keywords
the manuscript pipeline uses, so the answer is comparable to our own cohort.

Two known traps this must survive, both seen already in this project:
  - "cervical" also means the NECK/SPINE (cervical spine, cervical cord). Those
    are not cervix and must not be counted.
  - HeLa is a cervical adenocarcinoma line, so it saturates any text search for
    "cervical". That is the very contamination we are trying to escape.

Writes cervix_survey.tsv (one row per run, with the classification and the
evidence that drove it) and prints the composition breakdown.
"""
import collections
import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.expanduser("~/Projects/hpv-human-skin/bin"))
import cell_line_patterns as clp  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API_KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if API_KEY else 0.36
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cervix_survey.tsv")

QUERY = (
    '("cervix"[All Fields] OR "cervical"[All Fields] OR "ectocervix"[All Fields] '
    'OR "endocervix"[All Fields] OR "cervicovaginal"[All Fields] OR "CIN"[All Fields] '
    'OR "cervical intraepithelial neoplasia"[All Fields]) '
    'AND "rna seq"[Strategy] AND "public"[Access] AND txid9606[Organism]'
)

# The neck/spine sense of "cervical". If one of these fires and no true cervix
# term does, the run is anatomy-mismatched and gets dropped.
SPINE = re.compile(
    r"(?<![a-z])cervical[\s_-]*(spine|spinal|cord|vertebra\w*|disc|disk|"
    r"radiculopath\w*|myelopath\w*|lymph[\s_-]*node|esophag\w*|dorsal[\s_-]*root)"
    r"|(?<![a-z])(c[1-7][\s_-]*c[1-7]|spinal[\s_-]*cord)(?![a-z])",
    re.I,
)
CERVIX = re.compile(
    r"(?<![a-z])(cervix|ectocervi\w*|endocervi\w*|cervicovaginal|"
    r"cervical[\s_-]*(cancer|carcinoma|tumou?r|tissue|biops\w*|epitheli\w*|"
    r"cytolog\w*|smear|swab|screening|lesion|dysplasia|neoplasia|scrap\w*|"
    r"mucosa|stroma|squamous))(?![a-z])"
    r"|(?<![a-z])(cin[\s_-]?[123]|hsil|lsil|ascus)(?![a-z])",
    re.I,
)
# Primary material actually taken from a person, as opposed to a cultured line.
PRIMARY = re.compile(
    r"(?<![a-z])(biops\w*|tissue|patient\w*|surgical|resect\w*|specimen|"
    r"swab|smear|scrap\w*|brush\w*|autops\w*|donor|cohort|paraffin|ffpe|"
    r"fresh[\s_-]*frozen|hysterectom\w*|conization|cone[\s_-]*biops\w*|leep)(?![a-z])",
    re.I,
)


def get(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=180).read()
        except Exception as exc:  # noqa: BLE001
            if i == tries - 1:
                print(f"  !! giving up: {exc}", file=sys.stderr)
                return b""
            time.sleep(3 * (i + 1))
    return b""


def eurl(tool, **kw):
    if API_KEY:
        kw["api_key"] = API_KEY
    return EUTILS + tool + "?" + urllib.parse.urlencode(kw)


def main():
    raw = get(eurl("esearch.fcgi", db="sra", term=QUERY, usehistory="y",
                   retmode="json", retmax=0))
    if not raw:
        sys.exit("esearch failed")
    import json
    res = json.loads(raw)["esearchresult"]
    total, webenv, qkey = int(res["count"]), res["webenv"], res["querykey"]
    print(f"public human cervical-term RNA-seq experiments: {total}", flush=True)

    rows = []
    step = 200
    for start in range(0, total, step):
        xml = get(eurl("efetch.fcgi", db="sra", WebEnv=webenv, query_key=qkey,
                       retstart=start, retmax=step, rettype="full", retmode="xml"))
        time.sleep(PAUSE)
        if not xml:
            continue
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            continue
        for pkg in root.iter("EXPERIMENT_PACKAGE"):
            # every scrap of free text this run carries
            bits = []
            for tag in ("TITLE", "DESIGN_DESCRIPTION", "STUDY_TITLE",
                        "STUDY_ABSTRACT", "SCIENTIFIC_NAME", "VALUE", "TAG",
                        "LIBRARY_STRATEGY", "LIBRARY_SOURCE"):
                bits += [e.text or "" for e in pkg.iter(tag)]
            text = " ".join(bits)
            low = text.lower()

            runs = [r.get("accession") for r in pkg.iter("RUN") if r.get("accession")]
            stud = pkg.find(".//STUDY")
            proj = stud.get("alias") if stud is not None else ""
            for ext in pkg.iter("EXTERNAL_ID"):
                if ext.get("namespace") == "BioProject":
                    proj = ext.text or proj
                    break
            title = ""
            t = pkg.find(".//STUDY_TITLE")
            if t is not None and t.text:
                title = t.text.strip()

            hit_cx = bool(CERVIX.search(low))
            hit_sp = bool(SPINE.search(low))
            is_cl, cl_hit = clp.classify_cell_line(low)
            if not is_cl:
                # cultured-but-unnamed material counts as non-primary too
                iv, iv_hit = clp.classify_in_vitro(low)
                if iv:
                    is_cl, cl_hit = True, f"in_vitro:{iv_hit}"
            prim = bool(PRIMARY.search(low))

            if hit_sp and not hit_cx:
                verdict = "neck_or_spine"
            elif not hit_cx:
                verdict = "no_cervix_term"
            elif is_cl:
                verdict = "cell_line"
            elif prim:
                verdict = "PRIMARY_TISSUE"
            else:
                verdict = "cervix_undocumented"

            for rn in runs:
                rows.append(dict(run=rn, bioproject=proj, verdict=verdict,
                                 cell_line_evidence=(str(cl_hit)[:60] if cl_hit else ""),
                                 study_title=title[:150]))
        if start and start % 2000 == 0:
            print(f"  ...{start}/{total} experiments, {len(rows)} runs", flush=True)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["run", "bioproject", "verdict",
                                           "cell_line_evidence", "study_title"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\nruns surveyed: {len(rows)}   ->  {OUT}\n")
    c = collections.Counter(r["verdict"] for r in rows)
    for k, n in c.most_common():
        print(f"  {k:<24}{n:>6}  {100*n/max(len(rows),1):5.1f}%")

    prim_rows = [r for r in rows if r["verdict"] == "PRIMARY_TISSUE"]
    if prim_rows:
        print(f"\ntop BioProjects with documented primary cervical tissue "
              f"({len(prim_rows)} runs):")
        by = collections.Counter(r["bioproject"] for r in prim_rows)
        titles = {r["bioproject"]: r["study_title"] for r in prim_rows}
        print(f"  {'bioproject':<14}{'runs':>5}  study")
        for p, n in by.most_common(25):
            print(f"  {p:<14}{n:>5}  {titles.get(p,'')[:92]}")


if __name__ == "__main__":
    main()
