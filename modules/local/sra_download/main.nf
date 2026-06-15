/*
 * SRA FASTQ Download
 *
 * Download FASTQ files from NCBI SRA using fasterq-dump.
 */

process SRA_DOWNLOAD {
    tag "${meta.srr_id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/sra-tools:3.1.1--h4304569_0'

    input:
    val meta  // [srr_id, srx_id, study, tissue_category, diagnosis, layout, ...]

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads

    script:
    """
    set -euo pipefail

    # Prefer the conda env's binaries over the cluster's spack stack on PATH.
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

    # Two-step download: prefetch is resumable and handles retries itself;
    # fasterq-dump then extracts from the local .sra with no network in the
    # decompress path. Use the task workdir for temp — /tmp is too small on
    # many SLURM compute nodes for fasterq-dump scratch (2-3x output size).
    prefetch --max-size 100g --output-directory . ${meta.srr_id}

    fasterq-dump \\
        --split-3 \\
        --threads ${task.cpus} \\
        --temp . \\
        ./${meta.srr_id}

    pigz -p ${task.cpus} *.fastq

    # Drop the .sra cache dir to keep publish footprint small
    rm -rf ./${meta.srr_id}
    """
}
