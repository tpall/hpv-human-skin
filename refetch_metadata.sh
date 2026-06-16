#!/usr/bin/env bash
#SBATCH --job-name=hpv-refetch
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=refetch_%j.log
set -euo pipefail

# Re-fetch SRA/BioSample metadata for an EXISTING sample set and re-enrich it.
#
# The samplesheet shipped with results_full_v2 predates richer attribute
# capture: it lacks `genotype` and a full `characteristics` blob, so cell-line /
# engineered detection was blind to hints misfiled into off-spec fields. This
# driver re-fetches the exact same accessions (no fresh search → identical
# cohort) with all attributes, then re-runs tissue/cell-line enrichment.
#
# Needs outbound access to NCBI Entrez. Set NCBI_API_KEY for higher rate limits.
#
# Usage:
#   sbatch refetch_metadata.sh [input_samplesheet] [output_dir]
#   bash   refetch_metadata.sh [input_samplesheet] [output_dir]
#
# Env toggles:
#   NCBI_API_KEY=...   NCBI key (optional, faster)
#   NCBI_EMAIL=...     Entrez email (default tapa741@gmail.com)
#   FORCE=1            re-fetch even if the wide output already exists

# Under SLURM, BASH_SOURCE points to the spool copy — fall back to submit dir.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

INPUT="${1:-${SCRIPT_DIR}/results_full_v2/metadata/samplesheet_enriched.csv}"
OUTDIR="${2:-$(dirname "${INPUT}")}"
CATEGORIES="${SCRIPT_DIR}/assets/tissue_categories.csv"
ENV_YML="${SCRIPT_DIR}/modules/local/sra_discovery/environment.yml"

WIDE="${OUTDIR}/samplesheet_wide_v2.csv"
ENRICHED="${OUTDIR}/samplesheet_enriched_v2.csv"

NCBI_EMAIL="${NCBI_EMAIL:-tapa741@gmail.com}"
FORCE="${FORCE:-0}"

[[ -f "${INPUT}" ]]      || { echo "ERROR: input samplesheet not found: ${INPUT}" >&2; exit 1; }
[[ -f "${CATEGORIES}" ]] || { echo "ERROR: tissue categories not found: ${CATEGORIES}" >&2; exit 1; }
[[ -f "${ENV_YML}" ]]    || { echo "ERROR: env spec not found: ${ENV_YML}" >&2; exit 1; }
mkdir -p "${OUTDIR}"

echo "Project:     ${SCRIPT_DIR}"
echo "Input sheet: ${INPUT}"
echo "Wide out:    ${WIDE}"
echo "Enriched out:${ENRICHED}"
echo ""

# ── Conda env (idempotent, mirrors setup.sh) ────────────────────────────
CONDA_CMD="conda"
command -v mamba &>/dev/null && CONDA_CMD="mamba"
eval "$(conda shell.bash hook)"

ENV_NAME="sra_discovery"
if conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Conda env '${ENV_NAME}' already exists — reusing"
else
    echo "Creating conda env '${ENV_NAME}' from ${ENV_YML}"
    ${CONDA_CMD} env create -n "${ENV_NAME}" -f "${ENV_YML}" -y
fi
conda activate "${ENV_NAME}"

# This cluster's spack stack shadows the conda env's python on PATH; prepend
# the env's bin so `python` / `efetch` resolve to the activated env.
export PATH="${CONDA_PREFIX:+$CONDA_PREFIX/bin:}${PATH}"

# ── Step 1: re-fetch wide metadata for the exact accessions ─────────────
if [[ -s "${WIDE}" && "${FORCE}" != "1" ]]; then
    echo "Wide metadata already present — skipping fetch (set FORCE=1 to redo): ${WIDE}"
else
    echo "=== Re-fetching SRA metadata for accessions in ${INPUT} ==="
    API_ARGS=()
    [[ -n "${NCBI_API_KEY:-}" ]] && API_ARGS+=(--api-key "${NCBI_API_KEY}")
    python "${SCRIPT_DIR}/bin/query_sra.py" \
        --accessions "${INPUT}" \
        --output "${WIDE}" \
        --email "${NCBI_EMAIL}" \
        ${API_ARGS[@]+"${API_ARGS[@]}"}
fi

# ── Step 2: re-enrich (tissue category + cell-line/engineered flags) ─────
echo ""
echo "=== Enriching ${WIDE} → ${ENRICHED} ==="
python "${SCRIPT_DIR}/bin/parse_metadata.py" \
    --input "${WIDE}" \
    --categories "${CATEGORIES}" \
    --output "${ENRICHED}"

conda deactivate

echo ""
echo "==== Done ===="
echo "Wide (all attributes):    ${WIDE}"
echo "Enriched (slim + flags):  ${ENRICHED}"
