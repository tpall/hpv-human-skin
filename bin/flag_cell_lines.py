#!/usr/bin/env python3
"""
Post-hoc cell-line flagging over an existing raw samplesheet.

Scans the ``title`` and ``tissue_source`` fields of ``samplesheet_raw.csv``
written by ``query_sra.py`` and emits a TSV that can be joined back to HPV
results by ``srr_id``.

Usage:
    flag_cell_lines.py \\
        --input results_full_v2/metadata/samplesheet_raw.csv \\
        --output results_full_v2/metadata/cell_line_flags.tsv

Detection logic lives in ``cell_line_patterns.py``.
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from cell_line_patterns import classify_cell_line


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", type=Path, required=True,
                   help="Path to samplesheet_raw.csv from query_sra.py")
    p.add_argument("--output", "-o", type=Path, required=True,
                   help="Output TSV with cell-line flags")
    args = p.parse_args()

    n_total = 0
    n_flagged = 0
    label_counts: dict[str, int] = {}

    with args.input.open() as fin, args.output.open("w", newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["srr_id", "is_cell_line", "matched_pattern", "source_text"])

        for row in reader:
            n_total += 1
            title = row.get("title", "")
            tissue_source = row.get("tissue_source", "")
            combined = f"{title} {tissue_source}".strip()
            is_line, matched = classify_cell_line(combined)
            if is_line:
                n_flagged += 1
                label_counts[matched] = label_counts.get(matched, 0) + 1
            writer.writerow([
                row.get("srr_id", ""),
                "true" if is_line else "false",
                matched,
                combined[:200],
            ])

    print(f"Scanned {n_total} samples; flagged {n_flagged} as cell-line "
          f"({100 * n_flagged / n_total:.1f}%)" if n_total else "No rows.",
          file=sys.stderr)
    if label_counts:
        print("Top matches:", file=sys.stderr)
        for label, n in sorted(label_counts.items(), key=lambda kv: -kv[1])[:15]:
            print(f"  {label:40s} {n}", file=sys.stderr)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
