#!/usr/bin/env python3
"""Find SRA specimens with BOTH a DNA and an RNA library, for a true paired
presence-vs-activity test.

Why this exists (PLAN.md section 5.5). Every presence-vs-activity comparison in
this project so far is a COMPOSITION comparison: the DNA side is metagenome
libraries and the RNA side is different libraries from the same body site. That
leaves a confound no amount of curation removes, because the two sets are
different specimens from different cohorts sequenced under different protocols.
A BioSample carrying both a DNA and an RNA library is the same specimen measured
both ways, which turns the comparison into a paired test.

Three tiers of pairing, weakest last, all reported separately:

  specimen  same BioSample carries both. The real thing.
  donor     same BioProject and a shared donor/subject token in the sample name,
            different BioSamples. Same person, probably not the same biopsy.
  study     same BioProject only. Weakest; useful for finding candidate studies
            to read by hand, not for a paired statistic.

Two modes:

  --niche anal|cervical|skin
      Site-anchored. Pulls every public human run for the niche across all
      library strategies, groups by BioSample, reports specimens with both.
      Finds pairs regardless of HPV status, which is what a paired test needs:
      HPV-negative DNA with HPV-negative RNA is an informative cell too.

  --shen PATH
      Catalogue-anchored. Takes the BioSamples of Shen/Logan HPV-POSITIVE DNA
      libraries and asks which of them also have an RNA-seq run. Higher yield
      per query because the DNA side is already known to carry HPV, but it can
      only ever find the HPV-positive half of the table.

Both modes write a TSV of candidate specimens with enough metadata to triage by
hand. Neither downloads sequence data.

Usage:
    presence_vs_activity/scope_paired_dna_rna.py --niche anal   --out anal_pairs.tsv
    presence_vs_activity/scope_paired_dna_rna.py --niche cervical --out cervix_pairs.tsv
    presence_vs_activity/scope_paired_dna_rna.py \
        --shen dna_rna_triangulation/cache/all_hits_info.list --out shen_pairs.tsv
"""
import argparse
import collections
import csv
import io
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# NCBI allows 3 requests/second without an API key, 10 with one. Stay under it.
DELAY = 0.40

HUMAN = '"Homo sapiens"[Organism]'
PUBLIC = '"public"[Access]'

# Site vocabularies. Deliberately WIDER than the RNA-side queries in
# rna_extension_scoping/, because a DNA library of the same biopsy is often
# annotated by the study's disease term rather than by its tissue.
NICHES = {
    "anal": ('(anal[Text Word] OR anorectal[Text Word] OR perianal[Text Word] '
             'OR "anal canal"[Text Word] OR "anal squamous"[Text Word] '
             'OR "anal intraepithelial"[Text Word] OR AIN[Text Word] '
             'OR "anal cancer"[Text Word] OR "anal carcinoma"[Text Word] '
             'OR "anal condyloma"[Text Word])'),
    "cervical": ('(cervical[Text Word] OR cervix[Text Word] OR ectocervix[Text Word] '
                 'OR endocervix[Text Word] OR vaginal[Text Word] OR vagina[Text Word] '
                 'OR vulvar[Text Word] OR vulva[Text Word] '
                 'OR "cervical intraepithelial"[Text Word] OR CIN[Text Word] '
                 'OR condyloma[Text Word])'),
    "skin": ('(skin[Text Word] OR cutaneous[Text Word] OR epidermis[Text Word] '
             'OR epidermal[Text Word] OR keratinocyte[Text Word] '
             'OR "actinic keratosis"[Text Word] OR cSCC[Text Word] '
             'OR verruca[Text Word] OR wart[Text Word] '
             'OR verruciformis[Text Word])'),
}

# LibrarySource is more reliable than LibraryStrategy for this split, because
# "OTHER" strategy is used by both metagenomic DNA and odd RNA protocols.
RNA_SOURCES = {"TRANSCRIPTOMIC", "METATRANSCRIPTOMIC",
               "TRANSCRIPTOMIC SINGLE CELL"}
DNA_SOURCES = {"GENOMIC", "METAGENOMIC", "GENOMIC SINGLE CELL",
               "SYNTHETIC", "VIRAL RNA"}
RNA_STRATEGIES = {"RNA-Seq", "miRNA-Seq", "ncRNA-Seq", "EST", "FL-cDNA",
                  "snRNA-seq", "scRNA-seq"}


def _get(url, tries=4, timeout=60):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # transient eutils failures are routine
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"eutils failed after {tries} tries: {last}")


def esearch_history(term):
    """Run an esearch and return (count, webenv, query_key)."""
    q = urllib.parse.urlencode({"db": "sra", "term": term, "retmax": 0,
                                "usehistory": "y"})
    xml = _get(f"{EUTILS}/esearch.fcgi?{q}")
    def grab(tag):
        m = re.search(rf"<{tag}>(.*?)</{tag}>", xml)
        return m.group(1) if m else None
    count = int(grab("Count") or 0)
    return count, grab("WebEnv"), grab("QueryKey")


# efetch rettype=runinfo under WebEnv paging returns CSV with NO header row, so
# the schema has to be supplied. This is the stable SRA RunInfo column order and
# it matches the header of Shen's all_hits_info.list exactly, which is the check
# that it has not drifted.
RUNINFO_COLS = [
    "Run", "ReleaseDate", "LoadDate", "spots", "bases", "spots_with_mates",
    "avgLength", "size_MB", "AssemblyName", "download_path", "Experiment",
    "LibraryName", "LibraryStrategy", "LibrarySelection", "LibrarySource",
    "LibraryLayout", "InsertSize", "InsertDev", "Platform", "Model",
    "SRAStudy", "BioProject", "Study_Pubmed_id", "ProjectID", "Sample",
    "BioSample", "SampleType", "TaxID", "ScientificName", "SampleName",
    "g1k_pop_code", "source", "g1k_analysis_group", "Subject_ID", "Sex",
    "Disease", "Tumor", "Affection_Status", "Analyte_Type",
    "Histological_Type", "Body_Site", "CenterName", "Submission",
    "dbgap_study_accession", "Consent", "RunHash", "ReadHash",
]


def fetch_runinfo(webenv, query_key, count, batch=500, label=""):
    """Page through efetch runinfo CSV, yielding dict rows."""
    header = None
    got = 0
    for start in range(0, count, batch):
        q = urllib.parse.urlencode({
            "db": "sra", "rettype": "runinfo", "retmode": "csv",
            "WebEnv": webenv, "query_key": query_key,
            "retstart": start, "retmax": batch})
        text = _get(f"{EUTILS}/efetch.fcgi?{q}")
        rdr = csv.reader(io.StringIO(text))
        for row in rdr:
            if not row or not row[0]:
                continue
            if row[0] == "Run":
                header = row          # some responses do carry one
                continue
            got += 1
            yield dict(zip(header or RUNINFO_COLS, row))
        print(f"    {label} runinfo {min(start + batch, count)}/{count}",
              file=sys.stderr)
        time.sleep(DELAY)


def classify(row):
    """Return 'rna', 'dna' or None for one runinfo row."""
    src = (row.get("LibrarySource") or "").strip().upper()
    strat = (row.get("LibraryStrategy") or "").strip()
    if src in RNA_SOURCES or strat in RNA_STRATEGIES:
        return "rna"
    if src in DNA_SOURCES:
        return "dna"
    return None


def donor_token(row):
    """A crude subject key: the sample name with trailing replicate/type
    suffixes stripped. Only used for the weaker 'donor' tier, never for the
    specimen tier, and always meant to be checked by hand afterwards."""
    name = (row.get("SampleName") or row.get("LibraryName") or "").strip()
    if not name:
        return None
    t = re.sub(r"[_\-. ]*(rna|dna|wgs|wxs|seq|rep\d*|r\d+|t\d+|n\d+)$", "",
               name, flags=re.I)
    return t.lower() or None


def summarise(rows, out_path, note, shen_runs=frozenset()):
    """Group runinfo rows by BioSample and report the three pairing tiers."""
    by_sample = collections.defaultdict(list)
    for r in rows:
        bs = (r.get("BioSample") or "").strip()
        if bs:
            by_sample[bs].append(r)

    specimen = {}
    for bs, rs in by_sample.items():
        kinds = {classify(r) for r in rs}
        if "rna" in kinds and "dna" in kinds:
            specimen[bs] = rs

    # donor tier: same BioProject, shared donor token, but DNA and RNA sit on
    # different BioSamples
    by_proj = collections.defaultdict(list)
    for r in rows:
        bp = (r.get("BioProject") or "").strip()
        if bp:
            by_proj[bp].append(r)

    donor_hits, study_hits = [], []
    for bp, rs in by_proj.items():
        kinds = {classify(r) for r in rs}
        if not ("rna" in kinds and "dna" in kinds):
            continue
        study_hits.append(bp)
        tok = collections.defaultdict(set)
        for r in rs:
            t = donor_token(r)
            k = classify(r)
            if t and k:
                tok[t].add(k)
        for t, ks in tok.items():
            if "rna" in ks and "dna" in ks:
                donor_hits.append((bp, t))

    print(f"\n=== {note} ===")
    print(f"  runs examined            {len(rows)}")
    print(f"  distinct BioSamples      {len(by_sample)}")
    print(f"  SPECIMEN pairs           {len(specimen)}   (same BioSample, DNA+RNA)")
    print(f"  DONOR pairs              {len(donor_hits)}   (same project + token)")
    print(f"  STUDY-level projects     {len(study_hits)}  (same project only)")

    cols = ["biosample", "bioproject", "n_dna", "n_rna", "dna_runs", "rna_runs",
            "dna_strategies", "rna_strategies", "organism", "sample_name",
            "dna_hpv_pos_shen", "rna_hpv_pos_shen", "priority"]
    n_priority = 0
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for bs, rs in sorted(specimen.items()):
            dna = [r for r in rs if classify(r) == "dna"]
            rna = [r for r in rs if classify(r) == "rna"]
            dna_pos = sum(1 for r in dna if r["Run"] in shen_runs)
            rna_pos = sum(1 for r in rna if r["Run"] in shen_runs)
            # The cell the paired test needs: HPV is already known to be present
            # in this specimen's DNA, so its matched RNA answers "present, but
            # is it active?" on the same material rather than across cohorts.
            priority = "1_dna_hpv_pos" if dna_pos else "2_unknown"
            if dna_pos:
                n_priority += 1
            w.writerow({
                "biosample": bs,
                "bioproject": rs[0].get("BioProject", ""),
                "n_dna": len(dna), "n_rna": len(rna),
                "dna_runs": ";".join(r["Run"] for r in dna[:8]),
                "rna_runs": ";".join(r["Run"] for r in rna[:8]),
                "dna_strategies": ";".join(sorted({r.get("LibraryStrategy", "") for r in dna})),
                "rna_strategies": ";".join(sorted({r.get("LibraryStrategy", "") for r in rna})),
                "organism": rs[0].get("ScientificName", ""),
                "sample_name": rs[0].get("SampleName", ""),
                "dna_hpv_pos_shen": dna_pos,
                "rna_hpv_pos_shen": rna_pos,
                "priority": priority,
            })
    print(f"  wrote {out_path}")
    if shen_runs:
        print(f"  of which DNA side is a KNOWN HPV+ Shen library: {n_priority}")

    if study_hits:
        print("  projects with both modalities (top 15 by run count):")
        ranked = sorted(study_hits, key=lambda b: -len(by_proj[b]))[:15]
        for bp in ranked:
            rs = by_proj[bp]
            nd = sum(1 for r in rs if classify(r) == "dna")
            nr = sum(1 for r in rs if classify(r) == "rna")
            print(f"    {bp}  dna={nd:5d} rna={nr:5d}  {rs[0].get('ScientificName','')}")
    return specimen


# SRA-side filters used to probe for a modality. These MUST match classify(),
# which keys on LibrarySource, or the query silently under-fetches: SAMD00005178
# carries 95 transcriptomic runs of which only 1 is LibraryStrategy "RNA-Seq"
# and 94 are "OTHER", so a strategy-only filter returned 1 of 95. Source is the
# authority; strategy is kept as a widener for records with an odd source.
DNA_STRAT_TERM = ('("genomic"[Source] OR "metagenomic"[Source] '
                  'OR "genomic single cell"[Source] '
                  'OR "WGS"[Strategy] OR "WXS"[Strategy] OR "Targeted-Capture"[Strategy] '
                  'OR "WGA"[Strategy] OR "AMPLICON"[Strategy] OR "ChIP-Seq"[Strategy] '
                  'OR "Bisulfite-Seq"[Strategy] OR "ATAC-seq"[Strategy] OR "Hi-C"[Strategy])')
RNA_STRAT_TERM = ('("transcriptomic"[Source] OR "metatranscriptomic"[Source] '
                  'OR "transcriptomic single cell"[Source] '
                  'OR "RNA-Seq"[Strategy] OR "miRNA-Seq"[Strategy] '
                  'OR "ncRNA-Seq"[Strategy] OR "FL-cDNA"[Strategy] OR "EST"[Strategy])')


def expand_projects(rows, batch=100):
    """Recover the modality a site query cannot see.

    A run only appears in the niche search if ITS OWN metadata carries the site
    term. Studies routinely label the RNA libraries "cervical tumour" and the
    matched DNA libraries something generic like "tumour DNA", so the project
    looks RNA-only and its pairs are invisible. For every project that currently
    shows one modality, ask SRA directly whether the project also holds the
    other, restricting the fetch to that modality so the cost stays bounded.
    """
    by_proj = collections.defaultdict(set)
    for r in rows:
        bp = (r.get("BioProject") or "").strip()
        k = classify(r)
        if bp and k:
            by_proj[bp].add(k)

    rna_only = [p for p, k in by_proj.items() if k == {"rna"}]
    dna_only = [p for p, k in by_proj.items() if k == {"dna"}]
    print(f"  projects: {len(by_proj)} total, {len(rna_only)} RNA-only, "
          f"{len(dna_only)} DNA-only", file=sys.stderr)

    extra = []
    for label, projects, want in (("DNA", rna_only, DNA_STRAT_TERM),
                                  ("RNA", dna_only, RNA_STRAT_TERM)):
        for i in range(0, len(projects), batch):
            chunk = projects[i:i + batch]
            term = ("(" + " OR ".join(f"{p}[BioProject]" for p in chunk) + ")"
                    f" AND {want} AND {PUBLIC}")
            try:
                count, we, qk = esearch_history(term)
            except RuntimeError as e:
                print(f"    probe failed, skipping chunk: {e}", file=sys.stderr)
                continue
            time.sleep(DELAY)
            if count:
                extra.extend(fetch_runinfo(we, qk, count,
                                           label=f"expand-{label}"))
            print(f"    expand {label}: {i + len(chunk)}/{len(projects)} "
                  f"projects probed, {len(extra)} extra runs", file=sys.stderr)
    print(f"  project expansion added {len(extra)} runs", file=sys.stderr)
    return extra


def load_shen_runs(path):
    """Run accessions of Shen/Logan HPV-POSITIVE libraries, for annotation.

    Kept in the tool rather than done ad hoc, because an unsaved post-processing
    step is exactly how the junction clustering statistic went missing.
    """
    if not path:
        return frozenset()
    with open(path) as fh:
        return frozenset(r["Run"] for r in csv.DictReader(fh) if r.get("Run"))


def mode_niche(niche, out_path, expand=False, shen_path=None):
    term = f"{HUMAN} AND {PUBLIC} AND {NICHES[niche]}"
    count, webenv, qk = esearch_history(term)
    print(f"[{niche}] {count} public human runs match the site vocabulary",
          file=sys.stderr)
    if not count:
        return
    rows = list(fetch_runinfo(webenv, qk, count, label=niche))
    shen_runs = load_shen_runs(shen_path)
    if expand:
        seen = {r["Run"] for r in rows}
        for r in expand_projects(rows):
            if r["Run"] not in seen:
                seen.add(r["Run"])
                rows.append(r)
    summarise(rows, out_path, f"niche = {niche}" + (" (project-expanded)" if expand else ""),
              shen_runs=shen_runs)


def mode_shen(path, out_path, batch=120):
    """BioSamples of Shen HPV-positive DNA libraries that also have RNA."""
    hpv_pos = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            bs = (r.get("BioSample") or "").strip()
            if bs and classify(r) == "dna":
                hpv_pos.setdefault(bs, r)
    samples = sorted(hpv_pos)
    print(f"[shen] {len(samples)} BioSamples with HPV-positive DNA libraries",
          file=sys.stderr)

    rows = []
    for i in range(0, len(samples), batch):
        chunk = samples[i:i + batch]
        term = " OR ".join(f"{s}[BioSample]" for s in chunk)
        count, webenv, qk = esearch_history(term)
        time.sleep(DELAY)
        if count:
            rows.extend(fetch_runinfo(webenv, qk, count,
                                      label=f"shen {i}-{i + len(chunk)}"))
        print(f"  [shen] {i + len(chunk)}/{len(samples)} biosamples, "
              f"{len(rows)} runs so far", file=sys.stderr)

    specimen = summarise(rows, out_path, "shen HPV-positive DNA BioSamples",
                         shen_runs=load_shen_runs(path))
    hits = [bs for bs in specimen if bs in hpv_pos]
    print(f"  of which the DNA side is a KNOWN HPV-POSITIVE library: {len(hits)}")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--niche", choices=sorted(NICHES))
    g.add_argument("--shen", metavar="PATH",
                   help="Shen all_hits_info.list (HPV-positive RunInfo)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--annotate-shen", metavar="PATH",
                    help="Shen all_hits_info.list; flags pairs whose DNA side is "
                         "already a known HPV-positive library (the priority set)")
    ap.add_argument("--expand-projects", action="store_true",
                    help="for every project showing one modality, ask SRA "
                         "whether it also holds the other (see expand_projects)")
    args = ap.parse_args()

    if args.niche:
        mode_niche(args.niche, args.out, expand=args.expand_projects,
                   shen_path=args.annotate_shen)
    else:
        mode_shen(args.shen, args.out)


if __name__ == "__main__":
    main()
