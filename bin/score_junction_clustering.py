#!/usr/bin/env python3
"""Score virus-host junctions by clustering, the statistic that separates
integrated from episomal libraries.

hpc/run_junction_pilot.sbatch emits per-library junction files and a summary
carrying n_junctions and n_distinct_sites. Neither is the discriminating
statistic. Raw junction COUNTS are uninformative and actively misleading: they
scale with viral depth, so the two epidermodysplasia-verruciformis libraries
produce the largest counts in the whole pilot (12,979 and 13,937) while being
episomal, because chimeric-read callers scatter noise in proportion to depth.

Real integration is a RECURRENT junction: many reads at one coordinate pair.
The statistic is therefore the fraction of a library's junctions falling at its
single most-supported site. Observed values are bimodal with an empty interval
between 0.03 and 0.33, and the manuscript places the threshold at 0.2 in it.

This step was originally run ad hoc and was not saved, so this script exists to
make the published junction_summary.tsv reproducible. It regenerates that file
byte for byte from the per-library junction files.

Junction file columns, written by the awk step in the sbatch:

    donor_contig  donor_pos  acceptor_contig  acceptor_pos  junction_type
    repeat_left   repeat_right

SITE DEFINITION. The published numbers group by the DONOR side alone
(donor_contig, donor_pos), which is what the sbatch's own n_distinct_sites
counts with `cut -f1,2`, so --site-key donor is the default and reproduces the
manuscript. Note this keys on one END of the junction rather than the pair, and
the donor is not always the viral side, since STAR orders by donor/acceptor
rather than by which reference a contig came from.

--site-key pair groups by the full (viral, human) coordinate pair, which is the
stricter reading of "one recurrent junction". The conclusion does not depend on
the choice: under pair the integrated libraries run 0.327 to 0.435 and the
episomal ones 0.000 to 0.034, so the gap stays empty and the 0.2 threshold
stays inside it. Pair is stricter on the noisy libraries in particular, pushing
the deepest EV library from 0.003 down to 0.000.

Usage:
    bin/score_junction_clustering.py \
        --junctions results_junction/junctions \
        --summary   results_junction/junction_summary.tsv \
        --out       results_junction/junction_summary.tsv
"""
import argparse
import csv
import collections
import pathlib
import sys


def score(path, site_key):
    """Return (n_junctions, n_sites, top_site_reads, frac_top) for one library."""
    sites = collections.Counter()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            key = (f[0], f[1]) if site_key == "donor" else (f[0], f[1], f[2], f[3])
            sites[key] += 1
    n = sum(sites.values())
    if n == 0:
        return 0, 0, 0, 0.0
    top = max(sites.values())
    return n, len(sites), top, top / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True,
                    help="directory of <SRR>.junction.tsv files")
    ap.add_argument("--summary", required=True,
                    help="summary emitted by run_junction_pilot.sbatch (supplies "
                         "srr_id, tier, hpv_reads, log2_e1e2, and the row order)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--site-key", choices=("donor", "pair"), default="donor",
                    help="donor (default, reproduces the published summary) or "
                         "pair (stricter; same conclusion, see module docstring)")
    args = ap.parse_args()

    jdir = pathlib.Path(args.junctions)
    rows = []
    with open(args.summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            srr = row["srr_id"]
            f = jdir / f"{srr}.junction.tsv"
            if not f.exists():
                print(f"  no junction file for {srr}, scoring as zero",
                      file=sys.stderr)
                n, ns, top, frac = 0, 0, 0, 0.0
            else:
                n, ns, top, frac = score(f, args.site_key)
            rows.append({
                "srr_id": srr,
                "tier": row["tier"],
                "hpv_reads": row["hpv_reads"],
                "log2_e1e2": row["log2_e1e2"],
                "n_junctions": n,
                "n_sites": ns,
                "top_site_reads": top,
                # 3 decimals: the bimodal gap is 0.03 to 0.33, so this is far
                # more precision than the call needs, and it keeps the file
                # byte-stable across reruns.
                "frac_top": f"{frac:.3f}",
            })

    # Row order follows the input summary, which follows the pilot samplesheet.
    cols = ["srr_id", "tier", "hpv_reads", "log2_e1e2",
            "n_junctions", "n_sites", "top_site_reads", "frac_top"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} libraries, site-key={args.site_key})")


if __name__ == "__main__":
    main()
