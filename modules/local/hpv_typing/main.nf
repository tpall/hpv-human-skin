/*
 * HPV Type Assignment
 *
 * Determines HPV type(s) per sample based on alignment coverage
 * breadth and depth against the HPV reference panel.
 *
 * Coverage metrics come from `samtools coverage` (per-reference numreads,
 * covbases, coverage%, meandepth). This deliberately avoids pysam, whose
 * conda build proved unreliable on the cluster; samtools is the same package
 * the STAR step already depends on.
 */

process HPV_TYPING {
    tag "${meta.srr_id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "${params.outdir}/hpv_typing", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path hpv_ref_fasta
    path hpv_ref_fai

    output:
    tuple val(meta), path("${meta.srr_id}_hpv_types.tsv"), emit: types
    path "${meta.srr_id}_hpv_coverage.tsv",                emit: coverage

    script:
    def header = 'sample_id\\thpv_reference\\tref_length\\tread_count\\tcovered_bases\\tcoverage_breadth\\tmean_depth'
    """
    # This cluster prepends a large spack software stack (incl. its own python)
    # to PATH in every SLURM task, ahead of the activated conda env's bin. So a
    # bare 'samtools' would resolve to the spack build, not the pinned conda
    # one. CONDA_PREFIX is set reliably by activation, so call the env binary by
    # absolute path; under singularity/docker (CONDA_PREFIX unset) it falls back
    # to the container's samtools.
    SAMTOOLS="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin/}samtools"

    # Per-reference coverage; -q 10 mirrors the previous MAPQ>=10 filter.
    # samtools coverage lists every reference in the BAM header, including
    # zero-read ones, so the full panel is reported.
    "\$SAMTOOLS" coverage -q 10 ${bam} > coverage_raw.tsv

    # Full coverage report (all references).
    {
        printf '${header}\\n'
        awk -v s="${meta.srr_id}" 'BEGIN{FS="\\t";OFS="\\t"} !/^#/ {
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%.4f\\t%.2f\\n", s, \$1, \$3, \$4, \$5, \$6/100.0, \$7
        }' coverage_raw.tsv
    } > ${meta.srr_id}_hpv_coverage.tsv

    # Assigned types: breadth >= min_coverage AND depth >= min_depth,
    # sorted by coverage breadth (column 6) descending.
    {
        printf '${header}\\n'
        awk -v s="${meta.srr_id}" -v mc=${params.hpv_min_coverage} -v md=${params.hpv_min_depth} \\
            'BEGIN{FS="\\t";OFS="\\t"} !/^#/ {
                breadth = \$6/100.0
                if (breadth >= mc && \$7 >= md)
                    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%.4f\\t%.2f\\n", s, \$1, \$3, \$4, \$5, breadth, \$7
            }' coverage_raw.tsv | sort -t\$'\\t' -k6,6gr
    } > ${meta.srr_id}_hpv_types.tsv

    n_types=\$(( \$(wc -l < ${meta.srr_id}_hpv_types.tsv) - 1 ))
    echo "${meta.srr_id}: \${n_types} HPV type(s) assigned" >&2
    """
}
