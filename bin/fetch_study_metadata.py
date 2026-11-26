#!/usr/bin/env python3
"""Build a curation worklist for libraries with no recorded tissue site.

`tissue_category == "muu"` means the SRA run fields carried no site keyword, not
that the site is something other than skin/anogenital/oral. Most of those runs
belong to studies whose *study-level* title or linked publication states the
tissue plainly -- the information exists, it just is not in the per-run fields
the pipeline parses.

This fetches, per study: the SRA study title, the BioProject title/description,
and the linked PubMed title/abstract. Output is sorted by how many unclassified
libraries the study would resolve, so curation effort goes where it pays.

The output is EVIDENCE, not a classification: fill in the `assign` column by
hand, then move the resolved rows into assets/tier_overrides.csv as
`study,<acc>,tissue_category,<nahk|anogenitaal|suuoos|muu>,<reason>` rows.

Usage:
    bin/fetch_study_metadata.py [--sheet results_full_v2/metadata/samplesheet_enriched.csv]
                                [--out study_curation_worklist.tsv] [--limit N]
"""
import argparse, csv, collections, html, re, sys, time, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PAUSE = 0.35          # NCBI allows 3 req/s without an API key


def get(url, tries=3):
    for _ in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
        except Exception:
            time.sleep(2)
    return ""


def flat(x, n=None):
    x = re.sub(r"\s+", " ", html.unescape(x or "")).strip()
    return x[:n] if n else x


def geo_series(gsm):
    """GEO sample accession -> parent series title + summary.

    GEO-submitted studies land in SRA with a per-sample title ("GSM8736547: 187F;
    Homo sapiens; RNA-Seq") that says nothing about tissue, while the parent GSE
    series title states it plainly. This recovers most of them.
    """
    ids = re.findall(r"<Id>(\d+)</Id>",
                     get(f"{EUTILS}esearch.fcgi?db=gds&term={gsm}[ACCN]&retmax=1&retmode=xml"))
    if not ids:
        return "", ""
    time.sleep(PAUSE)
    blob = html.unescape(get(f"{EUTILS}esummary.fcgi?db=gds&id={ids[0]}&retmode=xml"))
    t = re.search(r'Name="title" Type="String">(.*?)<', blob, re.S)
    s = re.search(r'Name="summary" Type="String">(.*?)<', blob, re.S)
    return flat(t.group(1) if t else "", 300), flat(s.group(1) if s else "", 400)


def study_info(acc):
    """SRA study accession -> title fields from SRA, GEO, BioProject and PubMed."""
    out = dict.fromkeys(
        ("sra_title", "bioproject", "bp_title", "bp_desc", "pmid", "pm_title",
         "geo_title", "geo_summary"), "")
    # Bare accession as the term; a [Study Accession] qualifier returns nothing.
    ids = re.findall(r"<Id>(\d+)</Id>", get(f"{EUTILS}esearch.fcgi?db=sra&term={acc}&retmax=1&retmode=xml"))
    if not ids:
        return out
    time.sleep(PAUSE)
    blob = html.unescape(get(f"{EUTILS}esummary.fcgi?db=sra&id={ids[0]}&retmode=xml"))
    m = re.search(r"<Title>(.*?)</Title>", blob, re.S)
    out["sra_title"] = flat(m.group(1) if m else "", 300)
    m = re.search(r"<Bioproject>(PRJ[A-Z0-9]+)</Bioproject>", blob)
    if m:
        out["bioproject"] = m.group(1)
    time.sleep(PAUSE)

    # GEO-submitted studies carry a per-sample GSM title that names no tissue;
    # the parent series does. Resolve it whenever the SRA title exposes a GSM.
    g = re.search(r"\b(GSM\d+)\b", out["sra_title"])
    if g:
        out["geo_title"], out["geo_summary"] = geo_series(g.group(1))

    if out["bioproject"]:
        bids = re.findall(r"<Id>(\d+)</Id>",
                          get(f"{EUTILS}esearch.fcgi?db=bioproject&term={out['bioproject']}[Project Accession]&retmode=xml"))
        time.sleep(PAUSE)
        if bids:
            bx = get(f"{EUTILS}efetch.fcgi?db=bioproject&id={bids[0]}&retmode=xml")
            m = re.search(r"<Title>(.*?)</Title>", bx, re.S)
            out["bp_title"] = flat(m.group(1) if m else "", 300)
            m = re.search(r"<Description>(.*?)</Description>", bx, re.S)
            out["bp_desc"] = flat(m.group(1) if m else "", 400)
            m = re.search(r"<Pubmed>(\d+)</Pubmed>", bx)
            if m:
                out["pmid"] = m.group(1)
            time.sleep(PAUSE)

    if out["pmid"]:
        px = get(f"{EUTILS}efetch.fcgi?db=pubmed&id={out['pmid']}&retmode=xml")
        m = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", px, re.S)
        out["pm_title"] = flat(re.sub(r"<[^>]+>", "", m.group(1)) if m else "", 300)
        time.sleep(PAUSE)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", default="results_full_v2/metadata/samplesheet_enriched.csv")
    ap.add_argument("--out", default="study_curation_worklist.tsv")
    ap.add_argument("--limit", type=int, default=0, help="only the top N studies by library count")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.sheet)))
    unc = [r for r in rows if r.get("tissue_category") == "muu"]
    by = collections.Counter(r["study"] for r in unc)
    studies = by.most_common(args.limit or None)
    total = len(unc)
    print(f"{total} unclassified libraries across {len(by)} studies; "
          f"querying {len(studies)}", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["study", "n_unclassified", "pct_of_unclassified", "assign", "note",
                    "best_title", "sra_title", "geo_title", "geo_summary",
                    "bioproject", "bp_title", "bp_desc", "pmid", "pm_title"])
        for i, (acc, n) in enumerate(studies, 1):
            info = study_info(acc)
            # The most informative title available, for one-glance curation.
            best = next((info[k] for k in ("geo_title", "bp_title", "pm_title", "sra_title")
                         if info[k] and not info[k].startswith("GSM")), info["sra_title"])
            w.writerow([acc, n, f"{100*n/total:.1f}", "", "", best,
                        info["sra_title"], info["geo_title"], info["geo_summary"],
                        info["bioproject"], info["bp_title"], info["bp_desc"],
                        info["pmid"], info["pm_title"]])
            fh.flush()
            print(f"  [{i}/{len(studies)}] {acc} ({n})", file=sys.stderr)

    print(f"wrote {args.out} -- fill in `assign`, then move rows into "
          f"assets/tier_overrides.csv", file=sys.stderr)


if __name__ == "__main__":
    main()
