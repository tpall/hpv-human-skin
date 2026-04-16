#!/usr/bin/env bash
#SBATCH --job-name=hpv-setup
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=setup_%j.log
set -euo pipefail

# HPV in Human Skin Pipeline — Setup
#
# Creates conda environments and builds all reference databases.
#
# Usage:
#   sbatch setup.sh           # submit to SLURM (uses #SBATCH defaults)
#   bash setup.sh [threads]   # run interactively
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

# ── Step 1: Build HPV references ────────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 1/2: Building HPV reference database"
echo "============================================"
echo ""

${CONDA_CMD} env remove -n build_refs -y 2>/dev/null || true
${CONDA_CMD} env create -f "${SCRIPT_DIR}/envs/build_refs.yml" -y
conda activate build_refs

bash "${SCRIPT_DIR}/bin/build_hpv_refs.sh" \
    "${SCRIPT_DIR}/assets/hpv_references" \
    "${SCRIPT_DIR}/assets/hpv_references/custom" \
    "${THREADS}"

conda deactivate

# ── Step 2: Build Kraken2 database ──────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 2/2: Building Kraken2 database"
echo "============================================"
echo ""

${CONDA_CMD} env remove -n build_kraken2_db -y 2>/dev/null || true
${CONDA_CMD} env create -f "${SCRIPT_DIR}/envs/build_kraken2_db.yml" -y
conda activate build_kraken2_db

bash "${SCRIPT_DIR}/bin/build_kraken2_db.sh" \
    "${SCRIPT_DIR}/assets/kraken2_db" \
    "${SCRIPT_DIR}/assets/hpv_references/hpv_all.fasta" \
    "${THREADS}"

conda deactivate

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
