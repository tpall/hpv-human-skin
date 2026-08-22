#!/usr/bin/env python3
"""Classify each library's RNA-selection chemistry: poly-A vs rRNA-depleted.

Sequencing depth is modelled as a covariate throughout the manuscript;
library chemistry is not, and it is the obvious methodological counterpart.
HPV early transcripts are polyadenylated but late (L1) transcripts are
recovered less efficiently by poly-A capture, so the chemistry a library was
built with bounds what any viral-transcript survey can see in it.

Two evidence sources, in priority order:

  1. `LibrarySelection` from SRA RunInfo, when it carries a real value
     ("PolyA", "Oligo-dT", "Inverse rRNA", "RANDOM"). Authoritative but
     generic ("cDNA") for about three quarters of this cohort.
  2. The free-text `LIBRARY_CONSTRUCTION_PROTOCOL` / `DESIGN_DESCRIPTION`.
     This is what resolves the "cDNA" stratum: the graded cervical validation
     series is filed as "cDNA" while its protocol says "Ribosomal RNA was
     deleted".

Deliberately conservative. A library whose text supports neither call, or
supports both, stays `unspecified` rather than being guessed at: a
misclassified covariate is worse than an honest missing level. Amplification
chemistries that are neither (10x 3'-barcoded, SMART template-switching) get
their own levels so they are not silently pooled with poly-A capture.

Usage:
  bin/classify_library_chemistry.py \
      [--selection results_reframe_fixed/library_selection.tsv] \
      [--protocol  results_reframe_fixed/library_protocol.tsv] \
      [--out       results_reframe_fixed/library_chemistry.tsv]
"""
import argparse
import csv
import collections
import re
import sys

# --- Level 1: the LibrarySelection field ------------------------------------
SEL_MAP = {
    "PolyA": "polyA", "Oligo-dT": "polyA", "cDNA_oligo_dT": "polyA",
    "Inverse rRNA": "rRNA_depleted",
    "RANDOM": "random_primed", "cDNA_randomPriming": "random_primed",
    "RANDOM PCR": "random_primed",
}

# --- Level 2: free-text patterns --------------------------------------------
# Ordered within each class; the first match is reported as the evidence, so
# put the least ambiguous phrases first.
RIBO_PAT = [
    r"ribosomal rna was deleted", r"ribo-?zero", r"ribo-?minus",
    r"ribo-?erase", r"ribo-?off", r"rrna deplet", r"ribosomal rna deplet",
    r"deplet\w* (?:of|the) ribosom", r"remov\w* (?:of )?(?:the )?ribosom",
    r"ribosomal rna (?:was |were )?(?:remov|deplet)",
    r"rrna (?:was |were )?(?:remov|deplet)", r"depleted of ribosomal",
    r"stranded total rna prep with ribo", r"total rna.{0,40}ribo-?zero",
]
POLYA_PAT = [
    r"poly-?\(?a\)?[\s-]*(?:selec|enrich|captur|purif|tail|rna|containing)",
    r"oligo[\s-]*\(?dt\)?", r"poly-?\(?t\)?[\s-]*oligo",
    r"truseq stranded mrna", r"truseq rna (?:sample|library)",
    r"nebnext poly\(a\)", r"nebnext.{0,30}mrna magnetic",
    r"dynabeads.{0,20}mrna", r"magnetic beads with oligo",
    r"quantseq", r"lexogen",
    r"mrna (?:was |were )?(?:purif|isolat|enrich|select|captur)",
    r"purif\w+ .{0,25}mrna", r"enrich\w+ .{0,25}mrna",
    r"mrna[\s-]*seq(?:uencing)? librar",
]
# Neither chemistry, and not to be pooled with poly-A capture.
TENX_PAT = [r"10x genomics", r"10x chromium", r"chromium single cell",
            r"chromium next gem"]
SMART_PAT = [r"smart-?seq", r"smart-?er", r"template[\s-]switching"]


def _first(text, pats):
    for p in pats:
        m = re.search(p, text)
        if m:
            return m.group(0)[:60]
    return None


def classify(sel_row, prot_row):
    """Return (chemistry, evidence, source)."""
    sel = (sel_row.get("LibrarySelection") or "").strip()
    if sel in SEL_MAP:
        return SEL_MAP[sel], sel, "LibrarySelection"

    text = " || ".join([(prot_row.get("library_construction_protocol") or ""),
                        (prot_row.get("design_description") or ""),
                        (prot_row.get("experiment_title") or "")]).lower()
    if not text.strip(" |"):
        return "unspecified", "", "none"

    ribo = _first(text, RIBO_PAT)
    poly = _first(text, POLYA_PAT)
    if ribo and poly:
        # e.g. an Illumina Stranded Total RNA description that also names mRNA
        # in passing. Not resolvable from text, so do not guess.
        return "unspecified", f"conflict: {ribo} / {poly}", "protocol"
    if ribo:
        return "rRNA_depleted", ribo, "protocol"
    if poly:
        return "polyA", poly, "protocol"
    tenx = _first(text, TENX_PAT)
    if tenx:
        return "tenx_3prime", tenx, "protocol"
    smart = _first(text, SMART_PAT)
    if smart:
        return "template_switch", smart, "protocol"
    return "unspecified", "", "protocol"


def read_tsv(path, key="srr_id"):
    with open(path, newline="") as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", default="results_reframe_fixed/library_selection.tsv")
    ap.add_argument("--protocol", default="results_reframe_fixed/library_protocol.tsv")
    ap.add_argument("--out", default="results_reframe_fixed/library_chemistry.tsv")
    args = ap.parse_args()

    sel = read_tsv(args.selection)
    prot = read_tsv(args.protocol)

    rows = []
    for acc, srow in sel.items():
        chem, ev, src = classify(srow, prot.get(acc, {}))
        rows.append({"srr_id": acc, "chemistry": chem, "evidence": ev,
                     "evidence_source": src,
                     "library_selection": srow.get("LibrarySelection", ""),
                     "library_source": srow.get("LibrarySource", "")})

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    print(f"wrote {n} rows to {args.out}", file=sys.stderr)
    for k, v in collections.Counter(r["chemistry"] for r in rows).most_common():
        print(f"  {v:6d}  {v/n:5.1%}  {k}", file=sys.stderr)
    src = collections.Counter(r["evidence_source"] for r in rows
                              if r["chemistry"] != "unspecified")
    print(f"  resolved by: {dict(src)}", file=sys.stderr)


if __name__ == "__main__":
    main()
