#!/usr/bin/env bash
#SBATCH --job-name=hpv-setup
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=setup_%j.log
set -euo pipefail

# HPV in Human Skin Pipeline — Setup
#
# Creates conda environments and builds all reference databases. Idempotent:
# on rerun, steps whose outputs already exist are skipped.
#
# Usage:
#   sbatch setup.sh               # submit to SLURM (uses #SBATCH defaults)
#   bash setup.sh [threads]       # run interactively
#   FORCE=1 sbatch setup.sh       # rebuild everything, even if outputs exist
#   FORCE_REFS=1 sbatch setup.sh  # rebuild only step 1 (HPV refs)
#   FORCE_KRAKEN=1 sbatch setup.sh  # rebuild only step 2 (Kraken2 DB)
#
# Prerequisites: conda or mamba must be in PATH

# Under SLURM, BASH_SOURCE points to the spool copy of the script, not the
# original location. Fall back to SLURM_SUBMIT_DIR in that case.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
THREADS="${SLURM_CPUS_PER_TASK:-${1:-8}}"

# Use mamba if available, fall back to conda
CONDA_CMD="conda"
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba for faster env creation"
fi

# Ensure we can activate envs in a script
eval "$(conda shell.bash hook)"

FORCE="${FORCE:-0}"
FORCE_REFS="${FORCE_REFS:-${FORCE}}"
FORCE_KRAKEN="${FORCE_KRAKEN:-${FORCE}}"

conda_env_exists() {
    conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "$1"
}

ensure_env() {
    # ensure_env <name> <yml> <force>
    local name="$1" yml="$2" force="$3"
    if conda_env_exists "${name}" && [[ "${force}" != "1" ]]; then
        echo "Conda env '${name}' already exists — reusing (set FORCE=1 to recreate)"
        return
    fi
    if conda_env_exists "${name}"; then
        ${CONDA_CMD} env remove -n "${name}" -y
    fi
    ${CONDA_CMD} env create -f "${yml}" -y
}

# ── Step 1: Build HPV references ────────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 1/2: Building HPV reference database"
echo "============================================"
echo ""

REFS_DIR="${SCRIPT_DIR}/assets/hpv_references"
# Sentinel files — if all present and non-empty, step 1 is complete.
REFS_SENTINELS=(
    "${REFS_DIR}/hpv_all.fasta"
    "${REFS_DIR}/hpv_all.fasta.fai"
    "${REFS_DIR}/hpv_gene_annotations.gff"
    "${REFS_DIR}/star_index/Genome"
    "${REFS_DIR}/hisat2_index/hpv_all.1.ht2"
)

refs_complete=1
for f in "${REFS_SENTINELS[@]}"; do
    [[ -s "$f" ]] || { refs_complete=0; break; }
done

if [[ "${refs_complete}" == "1" && "${FORCE_REFS}" != "1" ]]; then
    echo "HPV reference outputs already present — skipping (set FORCE_REFS=1 to rebuild):"
    for f in "${REFS_SENTINELS[@]}"; do echo "  ✓ $f"; done
else
    ensure_env build_refs "${SCRIPT_DIR}/envs/build_refs.yml" "${FORCE_REFS}"
    conda activate build_refs
    bash "${SCRIPT_DIR}/bin/build_hpv_refs.sh" \
        "${REFS_DIR}" \
        "${REFS_DIR}/custom" \
        "${THREADS}"
    conda deactivate
fi

# ── Step 2: Build Kraken2 database ──────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 2/2: Building Kraken2 database"
echo "============================================"
echo ""

K2_DIR="${SCRIPT_DIR}/assets/kraken2_db"
K2_SENTINELS=(
    "${K2_DIR}/hash.k2d"
    "${K2_DIR}/opts.k2d"
    "${K2_DIR}/taxo.k2d"
)

k2_complete=1
for f in "${K2_SENTINELS[@]}"; do
    [[ -s "$f" ]] || { k2_complete=0; break; }
done

if [[ "${k2_complete}" == "1" && "${FORCE_KRAKEN}" != "1" ]]; then
    echo "Kraken2 database already present — skipping (set FORCE_KRAKEN=1 to rebuild):"
    for f in "${K2_SENTINELS[@]}"; do echo "  ✓ $f"; done
else
    ensure_env build_kraken2_db "${SCRIPT_DIR}/envs/build_kraken2_db.yml" "${FORCE_KRAKEN}"
    conda activate build_kraken2_db
    bash "${SCRIPT_DIR}/bin/build_kraken2_db.sh" \
        "${K2_DIR}" \
        "${REFS_DIR}/hpv_all.fasta" \
        "${THREADS}"
    conda deactivate
fi

# ── Done ────────────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Setup complete!"
echo "============================================"
echo ""
echo "References: ${SCRIPT_DIR}/assets/hpv_references/"
echo "Kraken2 DB: ${SCRIPT_DIR}/assets/kraken2_db/"
echo ""
echo "Run the pipeline:"
echo "  nextflow run main.nf -profile conda,slurm --outdir results"
