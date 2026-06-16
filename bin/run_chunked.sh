#!/usr/bin/env bash
#SBATCH --job-name=hpv-chunked
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --output=chunked_%j.log
set -euo pipefail

# Nextflow creates conda envs in this driver process (not in SLURM task
# allocations), and builds several distinct envs concurrently. An 8G/2-CPU
# driver could OOM-kill a mamba extraction mid-flight, leaving a partial env
# directory that Nextflow then reuses — surfacing much later as a spurious
# "ModuleNotFoundError" in a task. 24G/4-CPU gives the concurrent builds room.

# Chunked driver for the HPV-in-skin pipeline.
#
# Splits a large samplesheet into fixed-size chunks, runs the pipeline
# once per chunk with its own work dir, aggregates per-sample summary
# artifacts into a shared directory, and wipes the chunk work dir on
# success. Completed chunks are marked with a `.done` sentinel so a
# rerun resumes at the first unfinished chunk. A final REPORT_ONLY run
# builds the report over all aggregated results.
#
# Usage:
#   sbatch bin/run_chunked.sh <samplesheet.csv> [chunk_size] [outdir] [workdir]
#   bash   bin/run_chunked.sh <samplesheet.csv> [chunk_size] [outdir] [workdir]
#
# Defaults: chunk_size=100, outdir=./results, workdir=./work

SAMPLESHEET="${1:?Usage: run_chunked.sh <samplesheet.csv> [chunk_size=100] [outdir=./results] [workdir=./work]}"
CHUNK_SIZE="${2:-100}"
OUTDIR="${3:-$(pwd)/results}"
WORKDIR="${4:-$(pwd)/work}"

# Resolve absolute paths so later pushd/popd or chunked paths stay valid
SAMPLESHEET="$(realpath "${SAMPLESHEET}")"
OUTDIR="$(realpath -m "${OUTDIR}")"
WORKDIR="$(realpath -m "${WORKDIR}")"

# Under SLURM, BASH_SOURCE points into /var/spool — fall back to SLURM_SUBMIT_DIR
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

# Nextflow needs Java 17+; this cluster's spack default is openjdk-11, which it
# rejects. Activate the conda 'java' env (openjdk 17) so the driver finds the
# right JVM regardless of the submit-time environment (e.g. when another batch
# job sbatch-submits this one), pinning JAVA_HOME/JAVA_CMD past any spack values.
eval "$(conda shell.bash hook)"
if conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx java; then
    conda activate java
    export JAVA_HOME="${CONDA_PREFIX}"
    export JAVA_CMD="${CONDA_PREFIX}/bin/java"
else
    echo "WARN: conda env 'java' not found — relying on ambient Java (need 17+)" >&2
fi

CHUNKS_DIR="${OUTDIR}/chunks"
AGG_DIR="${OUTDIR}/aggregated"
REPORT_DIR="${OUTDIR}/report"
mkdir -p "${CHUNKS_DIR}" "${AGG_DIR}" "${REPORT_DIR}" "${WORKDIR}"

[[ -f "${SAMPLESHEET}" ]] || { echo "ERROR: samplesheet not found: ${SAMPLESHEET}" >&2; exit 1; }

echo "Project:     ${PROJECT_DIR}"
echo "Samplesheet: ${SAMPLESHEET}"
echo "Chunk size:  ${CHUNK_SIZE}"
echo "Outdir:      ${OUTDIR}"
echo "Workdir:     ${WORKDIR}"
echo ""

# ── Split samplesheet (idempotent) ──────────────────────────────────────
# Slim the samplesheet to the 7 columns Nextflow's splitCsv consumes,
# dropping free-text fields (title, tissue_source, platform) whose embedded
# commas get column-shifted because splitCsv does not honour RFC 4180
# quoting. The original wide file stays intact for REPORT_ONLY, where
# R's readr parses it correctly.
SLIM_SHEET="${CHUNKS_DIR}/samplesheet_slim.csv"
python3 - "${SAMPLESHEET}" "${SLIM_SHEET}" <<'PY'
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
cols = ["srr_id", "srx_id", "study", "layout",
        "tissue_category", "diagnosis", "needs_curation"]
with open(src, newline="") as fi, open(dst, "w", newline="") as fo:
    r = csv.DictReader(fi)
    w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for row in r:
        w.writerow(row)
PY
HEADER="$(head -n 1 "${SLIM_SHEET}")"
SPLIT_MARKER="${CHUNKS_DIR}/.split_done"
if [[ ! -e "${SPLIT_MARKER}" ]]; then
    echo "Splitting samplesheet into chunks of ${CHUNK_SIZE} …"
    tail -n +2 "${SLIM_SHEET}" \
        | split -l "${CHUNK_SIZE}" --numeric-suffixes=1 --suffix-length=4 \
                --additional-suffix=.body - "${CHUNKS_DIR}/chunk_"
    for body in "${CHUNKS_DIR}"/chunk_*.body; do
        tag="$(basename "${body}" .body)"
        chunk_dir="${CHUNKS_DIR}/${tag}"
        mkdir -p "${chunk_dir}"
        { echo "${HEADER}"; cat "${body}"; } > "${chunk_dir}/samplesheet.csv"
        rm -f "${body}"
    done
    touch "${SPLIT_MARKER}"
fi

mapfile -t CHUNK_DIRS < <(find "${CHUNKS_DIR}" -mindepth 1 -maxdepth 1 -type d -name 'chunk_*' | sort)
N_CHUNKS="${#CHUNK_DIRS[@]}"
echo "Total chunks: ${N_CHUNKS}"
echo ""

# ── Process each chunk ──────────────────────────────────────────────────
aggregate_chunk() {
    local chunk_dir="$1"
    # Per-sample artifacts: just copy each file into AGG_DIR (unique SRR prefix).
    for sub in hpv_typing transcript_classification kraken2; do
        if [[ -d "${chunk_dir}/${sub}" ]]; then
            find "${chunk_dir}/${sub}" -maxdepth 1 -type f \
                -exec cp -a {} "${AGG_DIR}/" \;
        fi
    done
    # Concatenate hpv_status.tsv (keep header once).
    local status_src="${chunk_dir}/metadata/hpv_status.tsv"
    local status_dst="${AGG_DIR}/hpv_status.tsv"
    if [[ -f "${status_src}" ]]; then
        if [[ ! -s "${status_dst}" ]]; then
            head -n 1 "${status_src}" > "${status_dst}"
        fi
        tail -n +2 "${status_src}" >> "${status_dst}"
    fi
}

i=0
for chunk_dir in "${CHUNK_DIRS[@]}"; do
    i=$((i + 1))
    tag="$(basename "${chunk_dir}")"
    chunk_csv="${chunk_dir}/samplesheet.csv"
    chunk_work="${WORKDIR}/${tag}"
    done_sentinel="${chunk_dir}/.done"

    if [[ -e "${done_sentinel}" ]]; then
        echo "[${i}/${N_CHUNKS}] ${tag} already complete — skipping"
        continue
    fi

    echo ""
    echo "==== [${i}/${N_CHUNKS}] ${tag}: starting ===="
    n_samples="$(( $(wc -l < "${chunk_csv}") - 1 ))"
    echo "  samples: ${n_samples}   work: ${chunk_work}"

    if nextflow run "${PROJECT_DIR}/main.nf" \
           -profile conda,slurm \
           -w "${chunk_work}" \
           -resume \
           -ansi-log false \
           --samplesheet "${chunk_csv}" \
           --outdir "${chunk_dir}" \
           --skip_report true
    then
        echo "[${tag}] succeeded — aggregating and cleaning work dir"
        aggregate_chunk "${chunk_dir}"
        rm -rf "${chunk_work}"
        touch "${done_sentinel}"
    else
        echo "[${tag}] FAILED — leaving ${chunk_work} in place for inspection"
        exit 1
    fi
done

# ── Final aggregated report ─────────────────────────────────────────────
echo ""
echo "==== All chunks done — running REPORT_ONLY ===="
nextflow run "${PROJECT_DIR}/main.nf" \
    -profile conda,slurm \
    -w "${WORKDIR}/report" \
    -ansi-log false \
    --entry REPORT_ONLY \
    --samplesheet "${SAMPLESHEET}" \
    --agg_dir "${AGG_DIR}" \
    --outdir "${REPORT_DIR}"

rm -rf "${WORKDIR}/report"

echo ""
echo "==== Done. Aggregated data: ${AGG_DIR}/  Report: ${REPORT_DIR}/ ===="
