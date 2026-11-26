#!/bin/bash
#SBATCH --job-name=hpv-paired
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --output=hpc/logs/hpv-paired_%j.log
#
# Type the paired presence-vs-activity RNA sets.
#
# WHY THIS SET. Every presence-vs-activity comparison in this project so far is
# a composition comparison: the DNA side is metagenome libraries and the RNA
# side is different libraries from the same body site, so the two are different
# specimens from different cohorts under different protocols. These 20
# BioSamples each carry BOTH modalities, and their DNA is already HPV-positive
# in the Shen/Logan catalogue. Typing the matched RNA therefore asks "HPV is
# present here, is it active?" on the same material.
#
#   cervix_paired  17 specimens, 27 runs   mucosal pole, presence SHOULD predict activity
#   skin_paired     3 specimens,  3 runs   cutaneous pole, presence should NOT
#
# 5 cervical specimens were dropped after job 520459 because their only RNA is
# long-read (PacBio Sequel FL-cDNA, one Nanopore GridION). STAR is a short-read
# aligner and aborts outright on a 1,144 bp read, so they need minimap2 rather
# than a threshold. Recovering them is real work, not a filter change.
#
# Scoping and selection: presence_vs_activity/PLAN.md section 5.5,
# presence_vs_activity/scope_paired_dna_rna.py, pairs_cervical.tsv,
# pairs_shen.tsv. Samplesheet built by build_paired_samplesheet.py.
#
# FALSIFIABLE PREDICTION, written before the run so it cannot be fitted after:
#   - the cervical specimens should show productive alpha transcription, with
#     L1/L2 late reads present in a substantial fraction;
#   - the 3 skin specimens should be HPV-positive at trace level at most, with
#     no late transcription, because the cutaneous beta/gamma reservoir is
#     transcriptionally silent;
#   - if the skin specimens come back productively transcribing, R2 of the
#     paper is wrong and the composition result was an artefact of comparing
#     lesional RNA against commensal DNA.
#
# 30 runs is a single chunk, so this finishes in hours, not days. The point is
# to get the paired result BEFORE committing to the 6,396-library cervical
# cohort, both because it is publishable on its own and because it de-risks the
# big run.
#
# Usage (from the repo root on rocket, via login1 or login2):
#   sbatch hpc/run_paired_pilot.sbatch [outdir] [workdir]
set -euo pipefail

PROJ="${SLURM_SUBMIT_DIR:-$HOME/Projects/hpv-human-skin}"
cd "$PROJ"
mkdir -p hpc/logs

SHEET="$PROJ/presence_vs_activity/paired_samplesheet.csv"
# KEEP THE SAME OUTDIR ACROSS RERUNS. run_chunked.sh does `cd "$OUTDIR"` before
# launching, and Nextflow keeps its resume history in <launchDir>/.nextflow, so
# the cache belongs to the OUTDIR and not to the work dir. A new outdir means
# -resume finds nothing and every download runs again, even with the work dir
# handed over intact. That mistake cost a full re-download here.
#
# But run_chunked.sh also pins the cohort with chunks/.split_done and will
# silently rerun the OLD cohort if the samplesheet changed underneath it. So
# when the sheet changes, do NOT reach for a new outdir. Instead clear the pin
# and let the cache stand:
#
#     rm -f  <outdir>/chunks/.split_done
#     rm -rf <outdir>/chunks/chunk_*
#     # leave <outdir>/.nextflow alone
#
# That satisfies both constraints: the cohort is re-split from the new sheet,
# and unchanged samples still hit the resume cache.
OUTDIR="${1:-$PROJ/results_paired}"
WORKDIR="${2:-$PROJ/work_paired}"

[ -s "$SHEET" ] || { echo "ERROR: samplesheet missing: $SHEET" >&2; exit 1; }
echo "[$(date +%T)] samplesheet: $SHEET ($(( $(wc -l < "$SHEET") - 1 )) runs)"

# run_chunked.sh writes .split_done into the chunks dir and will then IGNORE a
# changed samplesheet in an existing outdir, silently rerunning the old cohort.
# A fresh outdir is the only safe way to start a new set.
if [ -e "$OUTDIR/chunks/.split_done" ]; then
    echo "ERROR: $OUTDIR already has chunks/.split_done." >&2
    echo "       Resume is fine, but if the samplesheet changed you must use a" >&2
    echo "       new outdir, because the existing split takes precedence." >&2
fi

# chunk size above the cohort size, so this is one chunk and one report
bash bin/run_chunked.sh "$SHEET" 50 "$OUTDIR" "$WORKDIR"

echo "[$(date +%T)] DONE"
echo "aggregated: $OUTDIR/aggregated"
