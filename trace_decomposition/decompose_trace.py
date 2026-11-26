#!/usr/bin/env python3
"""Decompose the sub-threshold HPV trace into gene-less, early-only and productive.

Answers the question the Discussion used to concede: is the below-floor signal a
real low-level transcript of the commensal reservoir, or a k-mer artefact at the
sensitivity floor? Each HPV-positive library already has per-gene read counts, so
the trace can be split directly.

Emits `trace_rows.json` next to this script; `figure_trace_decomposition.R` turns
that into the figure, the composition table and the numbers the manuscript quotes.

Repo-local: paths are resolved from this file, the detection floor is read off the
sample-group table rather than pinned, and nothing is written to /tmp. All three
were hardcoded before, which survived only because the values happened to be right.
"""
import glob, csv, os, re, statistics as st, json
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
RPM = os.path.join(BASE, "results_reframe_fixed", "hpv_rpm_per_sample.tsv")
OUT = os.path.join(HERE, "trace_rows.json")


def quantile7(xs, p):
    """R's default quantile (type 7), so the floor here matches hpv_rpm_contrast.R."""
    xs = sorted(xs)
    if not xs:
        return float("nan")
    h = (len(xs) - 1) * p
    lo = int(h)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (h - lo) * (xs[hi] - xs[lo])


# ---- metadata + control-calibrated floor -----------------------------------
meta = {}
for r in csv.DictReader(open(RPM), delimiter="\t"):
    meta[r["srr_id"]] = dict(tier=r["tier"],
                             cl=r["is_cell_line"].upper() == "TRUE",
                             rpm=float(r["rpm"] or 0))

FLOOR = quantile7([m["rpm"] for m in meta.values()
                   if m["tier"] == "control_productive" and not m["cl"] and m["rpm"] > 0], 0.05)


# ---- per-library gene counts ------------------------------------------------
def index(suffix):
    d = {}
    for f in glob.glob(BASE + "/**/*" + suffix, recursive=True):
        s = os.path.basename(f).replace(suffix, "")
        d.setdefault(s, f)
    return d


TC = index("_transcript_classes.tsv")
TY = index("_hpv_types.tsv")


def genes(f):
    d = {r["gene"]: r["read_count"] for r in csv.DictReader(open(f), delimiter="\t")}
    gi = lambda k: (int(d[k]) if d.get(k, "").lstrip("-").isdigit() else 0)
    return dict(E6=gi("E6"), E7=gi("E7"), L1=gi("L1"), L2=gi("L2"),
                EARLY=gi("EARLY_TOTAL"), LATE=gi("LATE_TOTAL"))


def domtype(s):
    """Deepest assigned reference, used only to ask whether the early signal is HeLa."""
    f = TY.get(s)
    if not f:
        return None
    best, br = None, -1
    for r in csv.DictReader(open(f), delimiter="\t"):
        rc = float(r["read_count"]) if (r.get("read_count") or "").replace(".", "").isdigit() else 0
        if rc > br:
            br = rc
            m = re.search(r"HPV[- ]?(\d+[A-Za-z]?)", re.sub("REF$", "", r["hpv_reference"] or ""))
            best = m.group(1).upper() if m else "cand"
    return best


rows = []
for s, f in TC.items():
    m = meta.get(s)
    if not m or m["cl"]:
        continue
    g = genes(f)
    tot = g["EARLY"] + g["LATE"]
    rows.append(dict(tier=m["tier"], rpm=m["rpm"], **g, e67=g["E6"] + g["E7"],
                     latefrac=(g["LATE"] / tot if tot else 0.0), type=domtype(s)))


def grp(name, pred):
    g = [r for r in rows if pred(r)]
    if not g:
        print(f"\n{name}: n=0")
        return
    n = len(g)
    trans = [r for r in g if r["e67"] > 0]
    hela = sum(r["type"] in ("16", "18") for r in trans)
    print(f"\n{name}: n={n}")
    print(f"  no gene-mapped reads: {sum(r['EARLY'] == 0 and r['L1'] < 3 for r in g)}/{n}")
    print(f"  E6/E7 > 0:            {len(trans)}/{n} ({100*len(trans)/n:.0f}%)")
    print(f"  L1>=3 (productive):   {sum(r['L1'] >= 3 for r in g)}/{n} ({100*sum(r['L1']>=3 for r in g)/n:.0f}%)")
    print(f"  median late fraction: {st.median(r['latefrac'] for r in g):.3f}")
    if trans:
        print(f"  E6/E7+ dom types: HPV16/18={hela}/{len(trans)}; "
              f"{dict(Counter(r['type'] for r in trans).most_common(6))}")


print("=== EARLY-vs-LATE DECOMPOSITION | floor=%.2f RPM | non-cell-line ===" % FLOOR)
print("transcript_classes total:", len(TC), "| joined non-cell-line:", len(rows))
# The manuscript claim is about skin and neoplasia. `other` is the residual sample
# group whose site was never recorded, so it is reported for comparison only and
# is NOT part of the published denominator.
grp("SUB-THRESHOLD skin (unselected, rpm<floor)",
    lambda r: r["tier"] == "unselected_skin" and r["rpm"] < FLOOR)
grp("SUB-THRESHOLD skin+neoplasia (rpm<floor)  [PUBLISHED SET]",
    lambda r: r["tier"] in ("unselected_skin", "neoplasia") and r["rpm"] < FLOOR)
grp("SUB-THRESHOLD incl. site-not-recorded     [comparison only]",
    lambda r: r["tier"] in ("unselected_skin", "neoplasia", "other") and r["rpm"] < FLOOR)
grp("PRODUCTIVE CONTROLS (warts+EV)",
    lambda r: r["tier"] in ("control_productive", "EV"))

json.dump(dict(floor=FLOOR, rows=rows), open(OUT, "w"))
print("\nwrote %s (%d rows, floor %.2f)" % (OUT, len(rows), FLOOR))
