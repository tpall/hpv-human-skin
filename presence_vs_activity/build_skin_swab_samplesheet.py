#!/usr/bin/env python3
"""Select the skin-swab RNA arms to type, from the scoping table.

This is the compartment-matched test (PLAN.md 5.5b). The skin DNA reservoir is
100% surface swabs; the RNA axis it is compared against is 100% bulk biopsy.
These libraries are swab RNA, so they sit in the same compartment as the DNA,
and stratifying them by host immune status separates the two explanations for
the productive paired result:

  healthy productive        -> the reservoir was never silent, it was sampled in
                               the wrong compartment; R2 has to be rewritten
  healthy silent, DOCK8/RAG
  productive                -> the decoupling survives and gains a host axis

Emits the accession CSV for bin/query_sra.py --accessions, carrying the arm
label so the downstream grouping does not have to be reconstructed by hand.

The unstated arm (mostly SRP552405 diabetic foot ulcers) is deliberately left
out: a chronic wound is neither normal skin nor a clean disease arm, and it
would blur the contrast this run exists to make.

Usage:
    presence_vs_activity/build_skin_swab_samplesheet.py \
        --scope presence_vs_activity/skin_metatranscriptomes.tsv \
        --out   presence_vs_activity/skin_swab_accessions.csv
"""
import argparse
import collections
import csv

ARMS = {"immunocompromised", "healthy"}

# Runs that carry an arm label but are not the material this cohort is defined
# by. Excluded on cohort definition, not for convenience: the arm is the whole
# point of this run, so a mislabelled arm is worse than a missing sample.
#
# SRP687109's three "healthy" runs are HFF-1 dermal fibroblast CULTURES
# (isolation_source = "human dermal fibroblasts (HFF-1 cell line)"), which the
# host-status classifier called healthy because the study text says "Control".
# That is the control-arm conflation already documented for this project: in
# SRA, "control" means the untreated arm of an experiment at least as often as
# it means a healthy subject. Cultured fibroblasts have no stratum corneum, so
# they cannot answer a question about surface versus deep compartments.
NOT_SWAB = {
    "SRR37812028": "HFF-1 fibroblast culture, not a skin swab (SRP687109)",
    "SRR37812029": "HFF-1 fibroblast culture, not a skin swab (SRP687109)",
    "SRR37812030": "HFF-1 fibroblast culture, not a skin swab (SRP687109)",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--arms", default=",".join(sorted(ARMS)),
                    help="comma-separated host_status values to include")
    args = ap.parse_args()

    want = {a.strip() for a in args.arms.split(",") if a.strip()}
    rows, excluded = [], []
    for r in csv.DictReader(open(args.scope), delimiter="\t"):
        if r["host_status"] not in want:
            continue
        if r["srr_id"] in NOT_SWAB:
            excluded.append((r["srr_id"], NOT_SWAB[r["srr_id"]]))
            continue
        rows.append({
            "srr_id": r["srr_id"],
            "study": r["study"],
            "layout": r["layout"],
            "biosample": r["biosample"],
            "bioproject": r["bioproject"],
            "arm": r["host_status"],
            "arm_evidence": r["status_evidence"],
        })
    rows.sort(key=lambda x: (x["arm"], x["study"], x["srr_id"]))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    if excluded:
        print("excluded as not-this-material (see NOT_SWAB):")
        for srr, why in excluded:
            print(f"    {srr}  {why}")
    print(f"wrote {args.out}: {len(rows)} runs")
    for a, n in collections.Counter(r["arm"] for r in rows).most_common():
        print(f"  {a:20s} {n:4d}")
    print("  by study:")
    for (a, s), n in collections.Counter(
            (r["arm"], r["study"]) for r in rows).most_common():
        print(f"    {a:20s} {s:12s} {n:4d}")


if __name__ == "__main__":
    main()
