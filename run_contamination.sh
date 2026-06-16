#!/usr/bin/env bash
#SBATCH --job-name=hpv-contam
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=contamination_%j.log
set -euo pipefail

# HPV contamination / cell-line flagging driver.
#
# Step 1 (cached): build the reference->genus lookup from PaVE taxonomy
#                  (bin/build_genus_lookup.py — needs outbound PaVE access).
# Step 2:          flag contamination / cell-line calls in base R
#                  (bin/flag_contamination.R — consumes is_cell_line /
#                  is_engineered from the enriched samplesheet + the genus
#                  lookup; flags, never excludes).
#
# Usage:
#   sbatch run_contamination.sh [diagnostic|flag] [breadth_threshold]
#   bash   run_contamination.sh diagnostic
#
# Env overrides (with defaults):
#   INDIR=results_full_v2/aggregated
#   REF_FAI=assets/hpv_references/hpv_all.fasta.fai
#   SAMPLESHEET=results_full_v2/metadata/samplesheet_enriched_v2.csv
#   GENUS_LOOKUP=results_full_v2/metadata/type_genus.csv
#   OUTDIR=results_full_v2/contamination
#   PAVE_JSON=...   # reuse a saved PaVE listing instead of querying (offline)
#   FORCE=1         # rebuild the genus lookup even if present

MODE="${1:-flag}"
THRESHOLD="${2:-0.15}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

INDIR="${INDIR:-${SCRIPT_DIR}/results_full_v2/aggregated}"
REF_FAI="${REF_FAI:-${SCRIPT_DIR}/assets/hpv_references/hpv_all.fasta.fai}"
SAMPLESHEET="${SAMPLESHEET:-${SCRIPT_DIR}/results_full_v2/metadata/samplesheet_enriched_v2.csv}"
GENUS_LOOKUP="${GENUS_LOOKUP:-${SCRIPT_DIR}/results_full_v2/metadata/type_genus.csv}"
OUTDIR="${OUTDIR:-${SCRIPT_DIR}/results_full_v2/contamination}"
ENV_YML="${SCRIPT_DIR}/envs/contamination.yml"
FORCE="${FORCE:-0}"

[[ -d "${INDIR}" ]]        || { echo "ERROR: --indir not found: ${INDIR}" >&2; exit 1; }
[[ -f "${SAMPLESHEET}" ]]  || { echo "ERROR: samplesheet not found: ${SAMPLESHEET}" >&2; exit 1; }
mkdir -p "${OUTDIR}" "$(dirname "${GENUS_LOOKUP}")"

echo "Mode:        ${MODE}   threshold: ${THRESHOLD}"
echo "Indir:       ${INDIR}"
echo "Samplesheet: ${SAMPLESHEET}"
echo "Genus lookup:${GENUS_LOOKUP}"
echo "Outdir:      ${OUTDIR}"
echo ""

# ── Conda env (idempotent, base R + stdlib python) ──────────────────────
CONDA_CMD="conda"
command -v mamba &>/dev/null && CONDA_CMD="mamba"
eval "$(conda shell.bash hook)"
ENV_NAME="contamination"
if conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Conda env '${ENV_NAME}' already exists — reusing"
else
    echo "Creating conda env '${ENV_NAME}' from ${ENV_YML}"
    ${CONDA_CMD} env create -n "${ENV_NAME}" -f "${ENV_YML}" -y
fi
conda activate "${ENV_NAME}"
# Beat the cluster's spack stack, which shadows the env's python / Rscript.
export PATH="${CONDA_PREFIX:+$CONDA_PREFIX/bin:}${PATH}"

# ── Step 1: genus lookup (cached) ───────────────────────────────────────
if [[ -s "${GENUS_LOOKUP}" && "${FORCE}" != "1" ]]; then
    echo "Genus lookup present — skipping build (set FORCE=1 to rebuild): ${GENUS_LOOKUP}"
else
    echo "=== Building genus lookup from PaVE ==="
    GENUS_ARGS=()
    [[ -f "${REF_FAI}" ]]            && GENUS_ARGS+=(--fai "${REF_FAI}")
    [[ -n "${PAVE_JSON:-}" ]]        && GENUS_ARGS+=(--pave-json "${PAVE_JSON}")
    python "${SCRIPT_DIR}/bin/build_genus_lookup.py" \
        --out "${GENUS_LOOKUP}" \
        ${GENUS_ARGS[@]+"${GENUS_ARGS[@]}"}
fi

# ── Step 2: contamination flagging ──────────────────────────────────────
echo ""
echo "=== Flagging (${MODE}) ==="
Rscript "${SCRIPT_DIR}/bin/flag_contamination.R" \
    --mode "${MODE}" \
    --indir "${INDIR}" \
    --genus-lookup "${GENUS_LOOKUP}" \
    --samplesheet "${SAMPLESHEET}" \
    --breadth-threshold "${THRESHOLD}" \
    --outdir "${OUTDIR}"

conda deactivate
echo ""
echo "==== Done. Outputs: ${OUTDIR}/ ===="
