#!/usr/bin/env python3
"""
Post-hoc sensitivity sweep over per-sample HPV coverage tables.

Reads every ``*_hpv_coverage.tsv`` written by HPV_TYPING (the unfiltered
companion to ``*_hpv_types.tsv``) and produces:

  1. A combined long-format coverage table across all samples.
  2. A sensitivity sweep showing, for each (breadth, depth) cutoff pair,
     how many samples gain a type assignment and which references dominate.
  3. A per-sample "best hit" table at the most permissive cutoff, so
     low-titre candidates can be inspected manually.

Usage:
    sweep_hpv_coverage.py \\
        --input-dir results_full_v2/aggregated \\
        --output-dir results_full_v2/sweep \\
        [--breadth 0.10,0.05,0.02,0.01] \\
        [--depth   2,1,0.5,0.25]

Stdlib only — no pandas/numpy dependency.
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path


COVERAGE_COLS = [
    "sample_id",
    "hpv_reference",
    "ref_length",
    "read_count",
    "covered_bases",
    "coverage_breadth",
    "mean_depth",
]


def parse_float_list(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip()]


def load_coverage_files(input_dir: Path) -> list[dict]:
    rows = []
    files = sorted(input_dir.glob("*_hpv_coverage.tsv"))
    if not files:
        sys.exit(f"No *_hpv_coverage.tsv files under {input_dir}")
    for fp in files:
        with fp.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                r["ref_length"] = int(r["ref_length"])
                r["read_count"] = int(r["read_count"])
                r["covered_bases"] = int(r["covered_bases"])
                r["coverage_breadth"] = float(r["coverage_breadth"])
                r["mean_depth"] = float(r["mean_depth"])
                rows.append(r)
    print(f"Loaded {len(rows)} rows from {len(files)} sample(s)", file=sys.stderr)
    return rows


def write_combined(rows: list[dict], out: Path) -> None:
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COVERAGE_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def sweep(rows: list[dict], breadths: list[float], depths: list[float], out: Path) -> None:
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_sample[r["sample_id"]].append(r)

    summary_cols = [
        "min_breadth",
        "min_depth",
        "n_samples_typed",
        "n_assignments",
        "n_unique_refs",
        "top_refs",
    ]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_cols, delimiter="\t")
        w.writeheader()
        for b, d in product(breadths, depths):
            samples_typed = 0
            assignments = 0
            ref_counter: Counter[str] = Counter()
            for sample_id, sample_rows in by_sample.items():
                hits = [
                    r for r in sample_rows
                    if r["coverage_breadth"] >= b and r["mean_depth"] >= d
                ]
                if hits:
                    samples_typed += 1
                    assignments += len(hits)
                    for h in hits:
                        ref_counter[h["hpv_reference"]] += 1
            top = ", ".join(f"{ref}:{n}" for ref, n in ref_counter.most_common(5))
            w.writerow({
                "min_breadth": b,
                "min_depth": d,
                "n_samples_typed": samples_typed,
                "n_assignments": assignments,
                "n_unique_refs": len(ref_counter),
                "top_refs": top,
            })


def best_hits(rows: list[dict], min_breadth: float, min_depth: float, out: Path) -> None:
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        if r["coverage_breadth"] >= min_breadth and r["mean_depth"] >= min_depth:
            by_sample[r["sample_id"]].append(r)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COVERAGE_COLS, delimiter="\t")
        w.writeheader()
        for sample_id in sorted(by_sample):
            hits = sorted(
                by_sample[sample_id],
                key=lambda r: (r["coverage_breadth"], r["mean_depth"]),
                reverse=True,
            )
            w.writerow(hits[0])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", type=Path, required=True,
                   help="Directory containing *_hpv_coverage.tsv files")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Where to write combined.tsv, sweep.tsv, best_hits.tsv")
    p.add_argument("--breadth", default="0.10,0.05,0.02,0.01",
                   help="Comma-separated breadth cutoffs (default: 0.10,0.05,0.02,0.01)")
    p.add_argument("--depth", default="2,1,0.5,0.25",
                   help="Comma-separated mean-depth cutoffs (default: 2,1,0.5,0.25)")
    args = p.parse_args()

    breadths = parse_float_list(args.breadth)
    depths = parse_float_list(args.depth)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = load_coverage_files(args.input_dir)

    combined = args.output_dir / "combined_coverage.tsv"
    sweep_out = args.output_dir / "sweep_summary.tsv"
    best = args.output_dir / "best_hits_permissive.tsv"

    write_combined(rows, combined)
    sweep(rows, breadths, depths, sweep_out)
    best_hits(rows, min_breadth=min(breadths), min_depth=min(depths), out=best)

    print(f"Wrote {combined}", file=sys.stderr)
    print(f"Wrote {sweep_out}", file=sys.stderr)
    print(f"Wrote {best} (cutoff: breadth>={min(breadths)}, depth>={min(depths)})", file=sys.stderr)


if __name__ == "__main__":
    main()
