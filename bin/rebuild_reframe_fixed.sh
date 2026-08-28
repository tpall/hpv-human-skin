#!/usr/bin/env bash
# Rebuild results_reframe_fixed/ — the tier table the manuscript reads.
#
# The lesion cohort is NOT one results dir. The published table has 987 lesion
# libraries, which is the union of results_lesions_csc/aggregated (757) and
# results_lesions_wide/aggregated (230); results_lesions/aggregated (53) is a
# subset of those. hpv_rpm_contrast.R takes a single --lesion-agg, so the union
# has to be staged first. Running the script on any one dir silently produces a
# different cohort (csc alone loses every wart/EV control; the built-in default
# loses 413 neoplasia libraries), so always go through this script.
#
#   bin/rebuild_reframe_fixed.sh [outdir]
#
# Acceptance check: 10,834 rows = 972 lesion + 9,862 null. (11,282 before the
# two exclusions in tier_overrides.csv: 433 16S-amplicon libraries from
# SRP572392, and the 15 condyloma-adjacent CON_* libraries of
# SRP314958/SRP314961.)
set -euo pipefail

cd "$(dirname "$0")/.."
OUTDIR="${1:-results_reframe_fixed}"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

mkdir -p "$STAGE/aggregated"
for d in results_lesions_csc results_lesions_wide results_lesions; do
  [ -d "$d/aggregated" ] || continue
  for f in "$PWD/$d"/aggregated/*; do ln -sf "$f" "$STAGE/aggregated/"; done
done

# Merged samplesheet, csc winning ties (it is the deeper curation).
python3 - "$STAGE/lesions_samplesheet_full.csv" <<'PY'
import csv, sys
rows = {}
for p in ("results_lesions_csc/lesions_samplesheet_full.csv",
          "results_lesions_wide/lesions_samplesheet_full.csv"):
    for r in csv.DictReader(open(p)):
        rows.setdefault(r["srr_id"], r)
cols = list(next(iter(rows.values())).keys())
with open(sys.argv[1], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows.values():
        w.writerow({c: r.get(c, "") for c in cols})
print(f"merged lesion samplesheet: {len(rows)} rows", file=sys.stderr)
PY

Rscript bin/hpv_rpm_contrast.R \
  --lesion-agg   "$STAGE/aggregated" \
  --lesion-sheet "$STAGE/lesions_samplesheet_full.csv" \
  --outdir       "$OUTDIR"

n=$(( $(wc -l < "$OUTDIR/hpv_rpm_per_sample.tsv") - 1 ))
echo "rows: $n (expected 10834)"
[ "$n" -eq 10834 ] || echo "WARNING: row count differs from the published cohort" >&2

# Everything downstream of the tier table. These are NOT rebuilt by the
# manuscript render: figure2/figureS1 are copies of this script's plots,
# figure4 and figureS3 come from separate scripts, and bayes_summary.rds is read from disk by
# manuscript.qmd, so a stale one silently keeps the old posteriors.
if [ "$OUTDIR" = "results_reframe_fixed" ]; then
  cp "$OUTDIR/hpv_rpm_contrast.png"        manuscript/figure2.png
  cp "$OUTDIR/hpv_rpm_contrast.pdf"        manuscript/figure2.pdf
  cp "$OUTDIR/depth_matched_detection.png" manuscript/figureS1.png
  cp "$OUTDIR/depth_matched_detection.pdf" manuscript/figureS1.pdf
  echo "refreshed manuscript/figure2.{png,pdf}, figureS1.{png,pdf}"
  Rscript manuscript/figure_genus.R >/dev/null && echo "refreshed manuscript/figure4.{png,pdf}"
  Rscript manuscript/figure_type_grid.R >/dev/null && echo "refreshed manuscript/figureS3.{png,pdf}"
  Rscript manuscript/bayesian_models.R >/dev/null && echo "refit manuscript/bayes_summary.{rds,tsv}"
  Rscript manuscript/figure_bayes.R >/dev/null &&
    echo "refreshed manuscript/figureS5.{png,pdf}"
  # figureS2 reads the per-sample table and calibrates its own floor from it, and
  # figure5/figureS4 reach it through the type catalogue, so both go stale on a
  # cohort change exactly as silently as bayes_summary.rds does.
  python3 trace_decomposition/decompose_trace.py >/dev/null &&
    Rscript trace_decomposition/figure_trace_decomposition.R >/dev/null &&
    echo "refreshed manuscript/figureS2.{png,pdf} + trace_numbers.tsv"
  Rscript dna_rna_triangulation/shen_skin_dna_vs_rna.R >/dev/null 2>&1 &&
    Rscript dna_rna_triangulation/skin_dna_vs_rna_corrected.R >/dev/null 2>&1 &&
    Rscript manuscript/figure_reservoir.R >/dev/null &&
    echo "refreshed manuscript/figure5.{png,pdf}, figureS4.{png,pdf} + reservoir_numbers.tsv"
  echo "now re-render: quarto render manuscript/manuscript.qmd"
fi
