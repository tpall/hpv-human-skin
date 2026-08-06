#!/bin/bash
#SBATCH --job-name=hpv-skinswab
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --output=hpc/logs/hpv-skinswab_%j.log
#
# The compartment test: is the cutaneous HPV reservoir silent, or was it being
# sampled in the wrong compartment?
#
# WHY. Paper 2's skin arm compares a DNA reservoir that is 100% surface swabs
# (0 of 941 Shen skin BioSamples are biopsy or tissue) against an RNA axis that
# is 100% bulk tissue biopsy. Papillomavirus transcription is
# differentiation-coupled and completes in the upper epidermal layers, which is
# what a swab collects and what a biopsy dilutes with dermis and basal
# keratinocytes. The two axes therefore differ in the one variable that most
# directly controls the read-out, and the decoupling may be an artefact of that.
#
# The paired run (job 523244) made it concrete: both HPV-positive paired skin
# specimens were swab metatranscriptomes and both were PRODUCTIVE, with
# beta/gamma types at up to 99.5% genome coverage and 654x depth. But they are
# RAG-deficient patients, so host immunodeficiency was confounded with
# compartment at n=2. This run separates them.
#
#   immunocompromised  66 runs   SRP157584 DOCK8 (54), SRP403492 RAG (12)
#   healthy            40 runs   SRP046152 (24), SRP519513 (16)
#
# All are skin-swab metatranscriptomes, so compartment is held constant and host
# immune status is the variable. DOCK8 deficiency is the strongest available
# immunocompromised arm: it is the primary immunodeficiency defined by cutaneous
# viral disease including HPV.
#
# Scoping: presence_vs_activity/scope_skin_metatranscriptome.py, which discards
# 1,747 of 2,228 raw hits as DNA mislabelled RNA-Seq, 16S amplicon, or non-human
# hosts (horse, mouse, sea star, dog, deer, newt, lizard, frog) filed under the
# skin-metagenome taxid. Arms: build_skin_swab_samplesheet.py, which also drops
# three HFF-1 fibroblast cultures that the host classifier had called "healthy"
# because their study text says "Control".
#
# FALSIFIABLE PREDICTION, written before the run:
#   - if the HEALTHY arm is productive at rates near the immunocompromised arm,
#     the reservoir was never silent and R2 is wrong: HPV activity is stratified
#     by tissue depth and the biopsy axis simply could not see it;
#   - if the healthy arm is silent or trace while DOCK8/RAG are productive, the
#     decoupling survives and gains a host-status axis, which is a cleaner
#     result than the one currently planned;
#   - if BOTH arms are silent, the productive paired specimens were something
#     idiosyncratic and need re-examining before anything is concluded.
#
# Either of the first two is publishable. Running 106 libraries before the
# 6,396-library cervical cohort is what decides which paper is being written.
#
# Usage (from the worktree root on rocket, via login1 or login2):
#   sbatch hpc/run_skin_swab_pilot.sbatch [outdir] [workdir]
set -euo pipefail

PROJ="${SLURM_SUBMIT_DIR:-$HOME/Projects/hpv-human-skin}"
cd "$PROJ"
mkdir -p hpc/logs

SHEET="$PROJ/presence_vs_activity/skin_swab_samplesheet.csv"
# KEEP THE SAME OUTDIR ACROSS RERUNS. run_chunked.sh does `cd "$OUTDIR"` before
# launching and Nextflow keeps its resume history in <launchDir>/.nextflow, so
# the cache belongs to the outdir, not the work dir. If the samplesheet changes,
# do NOT switch outdir; clear the cohort pin in place instead:
#     rm -f  <outdir>/chunks/.split_done
#     rm -rf <outdir>/chunks/chunk_*
#     # leave <outdir>/.nextflow alone
OUTDIR="${1:-$PROJ/results_skin_swab}"
WORKDIR="${2:-$PROJ/work_skin_swab}"

[ -s "$SHEET" ] || { echo "ERROR: samplesheet missing: $SHEET" >&2; exit 1; }
echo "[$(date +%T)] samplesheet: $SHEET ($(( $(wc -l < "$SHEET") - 1 )) runs)"

if [ -e "$OUTDIR/chunks/.split_done" ]; then
    echo "NOTE: $OUTDIR/chunks/.split_done exists. Resuming the SAME cohort is" >&2
    echo "      fine; if the samplesheet changed, clear the pin as above." >&2
fi

bash bin/run_chunked.sh "$SHEET" 60 "$OUTDIR" "$WORKDIR"

echo "[$(date +%T)] DONE"
echo "aggregated: $OUTDIR/aggregated"
