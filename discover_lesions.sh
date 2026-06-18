#!/usr/bin/env bash
#SBATCH --job-name=hpv-lesions
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=lesions_%j.log
set -euo pipefail

# Targeted discovery of cutaneous-HPV lesions, two axes the broad cohort
# under-samples: (1) productive lesions (epidermodysplasia verruciformis,
# warts, condyloma) carrying genuine L1+ productive HPV; (2) beta-HPV
# skin-oncogenesis (actinic keratosis, cutaneous SCC, Bowen disease) where
# beta-papillomavirus is hypothesised to act as a UV cofactor. Query terms
# were derived/refined by mining the metadata of an earlier broad pull.
#
# Steps (this job, lightweight):
#   1. query SRA for EV / wart RNA-seq accessions
#   2. enrich (tissue category + cell-line/in-vitro/diagnosis flags)
#   3. diff against the main cohort -> write a NEW-only samplesheet so we never
#      re-type samples already processed in results_full_v2
# Then (separate SLURM job, unless AUTO_TYPE=0): submit the typing pipeline on
# just the new samples via bin/run_chunked.sh.
#
# Usage:
#   sbatch discover_lesions.sh [outdir]
#   bash   discover_lesions.sh results_lesions
#
# Env overrides (with defaults):
#   NCBI_API_KEY=...     NCBI key (optional, faster)
#   NCBI_EMAIL=...       Entrez email (default tapa741@gmail.com)
#   QUERY=...            override the lesion search query
#   MAIN_COHORT=results_full_v2/metadata/samplesheet_enriched_v2.csv
#   EXCLUDE_COHORTS=...  extra already-typed samplesheets (comma/space list)
#                        whose srr_ids are ALSO excluded, on top of MAIN_COHORT
#                        — e.g. a prior lesions batch so it is not re-typed:
#                          EXCLUDE_COHORTS=results_lesions_wide/lesions_samplesheet_new.csv
#   CHUNK_SIZE=100
#   AUTO_TYPE=1          submit the typing run (set 0 to stop after discovery).
#                        run_chunked.sh self-activates the conda 'java' env, so
#                        the submitted job gets Java 17 regardless of this env.

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

OUTDIR="${1:-${SCRIPT_DIR}/results_lesions}"
CATEGORIES="${SCRIPT_DIR}/assets/tissue_categories.csv"
ENV_YML="${SCRIPT_DIR}/modules/local/sra_discovery/environment.yml"
MAIN_COHORT="${MAIN_COHORT:-${SCRIPT_DIR}/results_full_v2/metadata/samplesheet_enriched_v2.csv}"
NCBI_EMAIL="${NCBI_EMAIL:-tapa741@gmail.com}"

# Cohorts whose srr_ids (column 1) are excluded from the "new" samplesheet so
# they are not re-typed: MAIN_COHORT plus any files in EXCLUDE_COHORTS.
EXCLUDE_FILES=()
[[ -f "${MAIN_COHORT}" ]] && EXCLUDE_FILES+=("${MAIN_COHORT}")
if [[ -n "${EXCLUDE_COHORTS:-}" ]]; then
    IFS=', ' read -r -a _extra_cohorts <<< "${EXCLUDE_COHORTS}"
    for _f in "${_extra_cohorts[@]}"; do
        [[ -z "${_f}" ]] && continue
        if [[ -f "${_f}" ]]; then EXCLUDE_FILES+=("${_f}")
        else echo "WARN: EXCLUDE_COHORTS entry not found, skipping: ${_f}" >&2; fi
    done
fi
CHUNK_SIZE="${CHUNK_SIZE:-100}"
AUTO_TYPE="${AUTO_TYPE:-1}"
QUERY="${QUERY:-\"epidermodysplasia verruciformis\" OR verruca OR \"verruca plana\" OR \"flat wart\" OR wart OR condyloma OR papilloma OR \"cutaneous papillomavirus\" OR \"beta-papillomavirus\" OR betapapillomavirus OR \"skin papilloma\" OR \"Lewandowsky-Lutz\" OR \"cutaneous squamous cell carcinoma\" OR \"actinic keratosis\" OR \"Bowen disease\"}"

RAW="${OUTDIR}/lesions_raw.csv"
FULL="${OUTDIR}/lesions_samplesheet_full.csv"
NEW="${OUTDIR}/lesions_samplesheet_new.csv"

[[ -f "${CATEGORIES}" ]] || { echo "ERROR: tissue categories not found: ${CATEGORIES}" >&2; exit 1; }
mkdir -p "${OUTDIR}"

echo "Project:     ${SCRIPT_DIR}"
echo "Outdir:      ${OUTDIR}"
echo "Exclude:     ${EXCLUDE_FILES[*]:-(none found)}"
echo "Query:       ${QUERY}"
echo ""

# ── Conda env (idempotent) ──────────────────────────────────────────────
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
export PATH="${CONDA_PREFIX:+$CONDA_PREFIX/bin:}${PATH}"

# ── Step 1: discovery ───────────────────────────────────────────────────
echo "=== Discovering lesion RNA-seq accessions ==="
API_ARGS=()
[[ -n "${NCBI_API_KEY:-}" ]] && API_ARGS+=(--api-key "${NCBI_API_KEY}")
python "${SCRIPT_DIR}/bin/query_sra.py" \
    --query "${QUERY}" \
    --output "${RAW}" \
    --email "${NCBI_EMAIL}" \
    ${API_ARGS[@]+"${API_ARGS[@]}"}

# ── Step 2: enrich ──────────────────────────────────────────────────────
echo ""
echo "=== Enriching ==="
python "${SCRIPT_DIR}/bin/parse_metadata.py" \
    --input "${RAW}" --categories "${CATEGORIES}" --output "${FULL}"

conda deactivate

# ── Step 3: keep only samples NOT already in an excluded cohort ─────────
echo ""
echo "=== Filtering to samples new vs already-typed cohorts ==="
head -1 "${FULL}" > "${NEW}"
if [[ ${#EXCLUDE_FILES[@]} -gt 0 ]]; then
    echo "  excluding srr_ids already in: ${EXCLUDE_FILES[*]}"
    # Read every exclude file first (accumulate seen srr_ids from column 1),
    # then emit rows of FULL whose srr_id was not seen. FULL is passed last and
    # identified by FILENAME so the number of exclude files does not matter.
    awk -F, -v target="${FULL}" '
        FILENAME==target { if (FNR>1 && !($1 in seen)) print; next }
        FNR>1 { seen[$1]=1 }
    ' "${EXCLUDE_FILES[@]}" "${FULL}" >> "${NEW}"
else
    echo "  (no exclude cohorts found — treating all as new)" >&2
    tail -n +2 "${FULL}" >> "${NEW}"
fi
n_total=$(( $(wc -l < "${FULL}") - 1 ))
n_new=$(( $(wc -l < "${NEW}") - 1 ))
echo "  lesion samples found: ${n_total}   new (not in excluded cohorts): ${n_new}"

# ── Step 4: type the new samples ────────────────────────────────────────
if [[ "${n_new}" -le 0 ]]; then
    echo "No new lesion samples to type — done."
    exit 0
fi
if [[ "${AUTO_TYPE}" != "1" ]]; then
    echo ""
    echo "AUTO_TYPE=0 — skipping typing. To type the new samples, run:"
    echo "  sbatch ${SCRIPT_DIR}/bin/run_chunked.sh ${NEW} ${CHUNK_SIZE} ${OUTDIR} ${OUTDIR}/work"
    exit 0
fi
echo ""
echo "=== Submitting typing run on ${n_new} new lesion samples ==="
# run_chunked.sh activates the conda 'java' env itself, so the submitted job
# gets Java 17 regardless of this helper's environment.
sbatch "${SCRIPT_DIR}/bin/run_chunked.sh" "${NEW}" "${CHUNK_SIZE}" "${OUTDIR}" "${OUTDIR}/work"
echo "Submitted. After it finishes, flag with:"
echo "  INDIR=${OUTDIR}/aggregated SAMPLESHEET=${FULL} OUTDIR=${OUTDIR}/contamination \\"
echo "    sbatch ${SCRIPT_DIR}/run_contamination.sh flag 0.15"
