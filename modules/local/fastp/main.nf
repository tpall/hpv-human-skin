/*
 * Read QC with fastp
 *
 * Adapter trimming, quality filtering, and QC reporting.
 */

process FASTP {
    tag "${meta.srr_id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_3'
    publishDir "${params.outdir}/qc/fastp", mode: params.publish_dir_mode, pattern: "*.json"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads
    path "${meta.srr_id}_fastp.json",             emit: json
    path "${meta.srr_id}_fastp.html",             emit: html

    script:
    if (meta.layout == "PAIRED") {
        """
        # Prefer the conda env's binaries over the cluster's spack stack on PATH.
        export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${meta.srr_id}_trimmed_R1.fastq.gz \\
            --out2 ${meta.srr_id}_trimmed_R2.fastq.gz \\
            --qualified_quality_phred ${params.fastp_min_quality} \\
            --length_required ${params.fastp_min_length} \\
            --detect_adapter_for_pe \\
            --thread ${task.cpus} \\
            --json ${meta.srr_id}_fastp.json \\
            --html ${meta.srr_id}_fastp.html
        """
    } else {
        """
        # Prefer the conda env's binaries over the cluster's spack stack on PATH.
        export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${meta.srr_id}_trimmed.fastq.gz \\
            --qualified_quality_phred ${params.fastp_min_quality} \\
            --length_required ${params.fastp_min_length} \\
            --thread ${task.cpus} \\
            --json ${meta.srr_id}_fastp.json \\
            --html ${meta.srr_id}_fastp.html
        """
    }
}
