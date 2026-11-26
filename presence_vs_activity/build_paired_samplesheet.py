#!/usr/bin/env python3
"""Build the accession list for the paired presence-vs-activity RNA sets.

Takes the pair tables from scope_paired_dna_rna.py and emits the accession CSV
that bin/query_sra.py --accessions consumes, so the paired RNA libraries go
through exactly the same metadata and curation path as every other cohort in
this project rather than a bespoke one.

Two sets, both defined by their DNA side already being HPV-positive in Shen's
catalogue, which is what makes them a paired test rather than a survey:

  cervix  22 specimens whose BioSample carries HPV-positive DNA and matched,
          never-screened RNA. The mucosal pole, where presence is expected to
          predict activity.
  skin     3 specimens from PRJNA872116 (skin metagenome, WGS + RNA-seq). The
          cutaneous pole, where presence is expected NOT to predict activity.

Run accessions are re-queried from SRA rather than read out of the pair table,
because that table truncates its run lists to 8 per side for readability and
one specimen (SAMD00005178) carries 95 RNA runs.

Usage:
    presence_vs_activity/build_paired_samplesheet.py \
        --cervical presence_vs_activity/pairs_cervical.tsv \
        --shen     presence_vs_activity/pairs_shen.tsv \
        --out      presence_vs_activity/paired_accessions.csv
"""
import argparse
import collections
import csv
import importlib.util
import pathlib
import sys
import time

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "scope_paired", _HERE / "scope_paired_dna_rna.py")
_scope = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_scope)


# Protocol suitability, on the same principle as the single-cell exclusion in
# bin/cell_line_patterns.py: a library whose coverage is structurally uneven
# across the viral early region cannot carry an E1:E2 ratio, and truncating or
# downsampling it does not help. Only whole-transcript protocols qualify.
#
# The case that forced this: SAMD00005178 contributed 95 of 132 runs and is
# FANTOM5 CAGE on the TC-YIK cervical carcinoma line (RIKEN RCB0443). CAGE is a
# 5'-end tag protocol, so it reports promoter usage rather than transcript body
# coverage. It is also a cell line, so it fails the clinical-specimen
# requirement independently.
USABLE_STRATEGIES = {"RNA-Seq", "FL-cDNA"}

# Read length, which is a HARDER constraint than strategy and is not visible in
# it. The pipeline aligns with STAR after fastp, both short-read tools, so a
# long-read library is not slow or noisy here, it is fatal: STAR aborts with
# "quality string length is not equal to sequence length" on a 1,144 bp read.
#
# The case that forced this: job 520459 ran 11 hours and died at the STAR step
# on SRR23829220, an OXFORD_NANOPORE GridION run whose LibraryStrategy is the
# perfectly innocent "RNA-Seq". Strategy cannot be trusted for this; Platform
# can. Seven PacBio Sequel FL-cDNA runs in the same set had the same problem
# and would have failed next.
#
# If long-read libraries are ever wanted, they need their own aligner
# (minimap2) rather than a threshold here.
SHORT_READ_PLATFORMS = {"ILLUMINA", "BGISEQ", "DNBSEQ"}


def rna_runs_for(biosamples, batch=80):
    """Every public RNA run belonging to these BioSamples, from SRA."""
    out = []
    bs = sorted(biosamples)
    for i in range(0, len(bs), batch):
        chunk = bs[i:i + batch]
        term = ("(" + " OR ".join(f"{s}[BioSample]" for s in chunk) + ")"
                f" AND {_scope.RNA_STRAT_TERM} AND {_scope.PUBLIC}")
        count, we, qk = _scope.esearch_history(term)
        time.sleep(_scope.DELAY)
        if count:
            out.extend(_scope.fetch_runinfo(we, qk, count, label="rna-runs"))
        print(f"  {i + len(chunk)}/{len(bs)} biosamples, {len(out)} RNA runs",
              file=sys.stderr)
    return [r for r in out if _scope.classify(r) == "rna"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cervical", required=True)
    ap.add_argument("--shen", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cerv = {r["biosample"] for r in
            csv.DictReader(open(args.cervical), delimiter="\t")
            if r.get("priority") == "1_dna_hpv_pos"}
    skin = {r["biosample"] for r in
            csv.DictReader(open(args.shen), delimiter="\t")
            if "skin metagenome" in r.get("organism", "")}
    print(f"cervix: {len(cerv)} specimens, skin: {len(skin)} specimens",
          file=sys.stderr)

    label = {b: "cervix_paired" for b in cerv}
    label.update({b: "skin_paired" for b in skin})

    runs = rna_runs_for(set(label))

    dropped = collections.Counter()
    keep = []
    for r in runs:
        strat = (r.get("LibraryStrategy") or "").strip()
        plat = (r.get("Platform") or "").strip().upper()
        if strat not in USABLE_STRATEGIES:
            dropped[f"strategy={strat or '(blank)'}"] += 1
            continue
        if plat not in SHORT_READ_PLATFORMS:
            dropped[f"platform={plat or '(blank)'}"] += 1
            continue
        keep.append(r)
    if dropped:
        print("  dropped on protocol suitability (see USABLE_STRATEGIES):",
              file=sys.stderr)
        for k, n in dropped.most_common():
            print(f"    {k or '(blank)'}: {n}", file=sys.stderr)
    runs = keep

    by_set = collections.Counter()
    seen = set()
    rows = []
    for r in runs:
        srr = r["Run"]
        if srr in seen:
            continue
        seen.add(srr)
        bs = r.get("BioSample", "")
        s = label.get(bs, "unknown")
        by_set[s] += 1
        rows.append({
            "srr_id": srr,
            "srx_id": r.get("Experiment", ""),
            "study": r.get("SRAStudy", ""),
            "layout": r.get("LibraryLayout", ""),
            "biosample": bs,
            "bioproject": r.get("BioProject", ""),
            "paired_set": s,
            "strategy": r.get("LibraryStrategy", ""),
            "platform": r.get("Platform", ""),
            "avg_length": r.get("avgLength", ""),
            "sample_name": r.get("SampleName", ""),
            "n_bases": r.get("bases", ""),
        })

    rows.sort(key=lambda x: (x["paired_set"], x["bioproject"], x["srr_id"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {args.out}: {len(rows)} RNA runs")
    for s, n in by_set.most_common():
        print(f"  {s:16s} {n}")
    per_bs = collections.Counter(r["biosample"] for r in rows)
    big = [(b, n) for b, n in per_bs.most_common() if n > 8]
    if big:
        print("  specimens with many RNA runs (check these are not lines):")
        for b, n in big:
            print(f"    {b}  {n} runs")


if __name__ == "__main__":
    main()
