#!/usr/bin/env python3
"""Triage the documented-primary-tissue cervical studies into a usable shortlist.

The survey found 2,803 runs of documented primary cervical tissue in public SRA.
Not all of it is usable for HPV transcript typing, and the failure modes are
specific:

  - 10x/Smart-seq single-cell libraries are 3'-biased and shallow per cell; they
    are poor for viral read recovery and cannot support the E1:E2 ratio, which
    needs even coverage across the early region.
  - Some studies are PERIPHERAL BLOOD from cervical-cancer patients. The text says
    "cervical cancer" and "patient", so the survey's PRIMARY rule passes them, but
    the material is not cervix.
  - Some accessions are embargoed (GEO release date in the future) and will fail
    at download.
  - Amplicon/16S runs ride along in mixed submissions.

Emits cervix_shortlist.tsv with an arm assignment (normal / tumour) so the two
sides of the "presence vs activity" contrast can be pulled separately.
"""
import collections
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.join(HERE, "cervix_survey.tsv")
OUT = os.path.join(HERE, "cervix_shortlist.tsv")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API_KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if API_KEY else 0.36

SINGLE_CELL = re.compile(
    r"(?<![a-z])(10x|10[\s_-]?genomics|chromium|single[\s_-]?cell|scrna|snrna|"
    r"smart[\s_-]?seq|drop[\s_-]?seq|cel[\s_-]?seq|single[\s_-]?nucle\w*|"
    r"visium|spatial[\s_-]?transcriptom\w*|slide[\s_-]?seq)(?![a-z])", re.I)
BLOOD = re.compile(
    r"(?<![a-z])(peripheral[\s_-]*blood|pbmc|plasma|serum|whole[\s_-]*blood|"
    r"buffy[\s_-]*coat|leukocyte|circulating)(?![a-z])", re.I)
# material that is cervix itself, not blood/lymph node from a cervix patient
CERVICAL_MATERIAL = re.compile(
    r"(?<![a-z])(cervix|ectocervi\w*|endocervi\w*|cervical[\s_-]*"
    r"(tissue|biops\w*|tumou?r|cancer|carcinoma|mucosa|epitheli\w*|swab|smear|"
    r"scrap\w*|brush\w*|specimen|stroma|squamous|lesion|cytolog\w*))(?![a-z])",
    re.I)
TUMOUR = re.compile(
    r"(?<![a-z])(tumou?r|cancer|carcinoma|malignan\w*|squamous[\s_-]*cell|"
    r"adenocarcinoma|figo|stage[\s_-]*[i1-4]|neoplas\w*|hsil|cin[\s_-]?[23]|"
    r"chemoradi\w*|neoadjuvant|metasta\w*)(?![a-z])", re.I)
NORMAL = re.compile(
    r"(?<![a-z])(normal|healthy|non[\s_-]*malignant|benign|control|"
    r"disease[\s_-]*free|unaffected|proliferative[\s_-]*phase|secretory[\s_-]*phase)"
    r"(?![a-z])", re.I)


def get(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=180).read()
        except Exception:  # noqa: BLE001
            if i == tries - 1:
                return b""
            time.sleep(3 * (i + 1))
    return b""


def eurl(tool, **kw):
    if API_KEY:
        kw["api_key"] = API_KEY
    return EUTILS + tool + "?" + urllib.parse.urlencode(kw)


def main():
    rows = list(csv.DictReader(open(SURVEY), delimiter="\t"))
    prim = [r for r in rows if r["verdict"] == "PRIMARY_TISSUE"]
    by_proj = collections.Counter(r["bioproject"] for r in prim)
    # study-level triage: anything with >=10 documented primary runs is worth a look
    cands = [p for p, n in by_proj.most_common() if n >= 10 and p.startswith("PRJ")]
    print(f"{len(prim)} primary-tissue runs in {len(by_proj)} studies; "
          f"{len(cands)} studies carry >=10\n", flush=True)

    out = []
    for i, proj in enumerate(cands, 1):
        ids_raw = get(eurl("esearch.fcgi", db="sra", term=f"{proj}[BioProject]",
                           retmax=400, retmode="json"))
        time.sleep(PAUSE)
        if not ids_raw:
            continue
        ids = json.loads(ids_raw)["esearchresult"].get("idlist", [])
        if not ids:
            continue
        xml = get(eurl("efetch.fcgi", db="sra", id=",".join(ids[:400]),
                       rettype="full", retmode="xml"))
        time.sleep(PAUSE)
        if not xml:
            continue
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            continue

        title, n_rna, sc, blood, cxmat, tum, nor, n_runs = "", 0, 0, 0, 0, 0, 0, 0
        for pkg in root.iter("EXPERIMENT_PACKAGE"):
            bits = []
            for tag in ("TITLE", "DESIGN_DESCRIPTION", "STUDY_TITLE",
                        "STUDY_ABSTRACT", "VALUE", "TAG", "INSTRUMENT_MODEL",
                        "LIBRARY_CONSTRUCTION_PROTOCOL", "LIBRARY_SELECTION"):
                bits += [e.text or "" for e in pkg.iter(tag)]
            low = " ".join(bits).lower()
            st = pkg.find(".//STUDY_TITLE")
            if st is not None and st.text and not title:
                title = st.text.strip()
            strat = pkg.find(".//LIBRARY_STRATEGY")
            strat = (strat.text if strat is not None else "") or ""
            runs = [r for r in pkg.iter("RUN")]
            n_runs += len(runs)
            if strat.upper() == "RNA-SEQ":
                n_rna += len(runs)
            if SINGLE_CELL.search(low):
                sc += len(runs)
            if BLOOD.search(low) and not CERVICAL_MATERIAL.search(low):
                blood += len(runs)
            if CERVICAL_MATERIAL.search(low):
                cxmat += len(runs)
            if TUMOUR.search(low):
                tum += len(runs)
            if NORMAL.search(low):
                nor += len(runs)

        usable = n_rna - sc - blood
        if tum and not nor:
            arm = "tumour"
        elif nor and not tum:
            arm = "normal"
        elif nor and tum:
            arm = "both"
        else:
            arm = "unclear"
        verdict = ("DROP:single-cell" if sc > 0.5 * max(n_runs, 1) else
                   "DROP:blood" if blood > 0.5 * max(n_runs, 1) else
                   "DROP:not-cervical-material" if cxmat == 0 else
                   "DROP:too-few-bulk" if usable < 10 else "KEEP")
        out.append(dict(bioproject=proj, n_runs=n_runs, n_rnaseq=n_rna,
                        n_single_cell=sc, n_blood=blood, est_usable_bulk=max(usable, 0),
                        arm=arm, verdict=verdict, study_title=title[:120]))
        print(f"  [{i}/{len(cands)}] {proj:<14}{verdict:<28}{arm:<8}"
              f"bulk~{max(usable,0):<5}{title[:56]}", flush=True)

    out.sort(key=lambda r: (r["verdict"] != "KEEP", -r["est_usable_bulk"]))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    keep = [r for r in out if r["verdict"] == "KEEP"]
    print(f"\n=== SHORTLIST: {len(keep)} studies, "
          f"~{sum(r['est_usable_bulk'] for r in keep)} usable bulk RNA-seq runs ===")
    for arm in ("normal", "both", "tumour", "unclear"):
        sel = [r for r in keep if r["arm"] == arm]
        if not sel:
            continue
        print(f"\n  --- {arm} ({sum(r['est_usable_bulk'] for r in sel)} runs) ---")
        for r in sel:
            print(f"    {r['bioproject']:<14}{r['est_usable_bulk']:>5}  {r['study_title'][:88]}")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
