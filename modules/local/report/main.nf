/*
 * Report Generation
 *
 * Aggregates all results and generates summary tables and HTML report.
 */

process REPORT {
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/r-tidyverse_r-rmarkdown_r-knitr_r-optparse:1.2.0--a6dcb8ce57f9c58f'
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    input:
    // stageAs keeps the in-workdir filenames stable and distinct even when
    // both inputs resolve to the same path (as in REPORT_ONLY's chunked
    // aggregation mode, where no separate raw samplesheet exists).
    path samplesheet,     stageAs: 'samplesheet_slim.csv'
    path raw_samplesheet, stageAs: 'samplesheet_raw.csv'
    path hpv_types_files
    path transcript_class_files
    path kraken_reports
    path hpv_status_csv

    output:
    path "hpv_skin_report.html",    emit: html
    path "summary_tables/",         emit: tables

    script:
    """
    # Prefer the conda env's binaries over the cluster's spack stack on PATH.
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

    shopt -s nullglob
    mkdir -p summary_tables

    # Merge all HPV type results (0 files when no sample was HPV+).
    types_files=( *_hpv_types.tsv )
    if (( \${#types_files[@]} > 0 )); then
        head -1 "\${types_files[0]}" > all_hpv_types.tsv
        for f in "\${types_files[@]}"; do
            tail -n +2 "\$f" >> all_hpv_types.tsv
        done
    else
        printf 'sample_id\\thpv_reference\\tref_length\\tread_count\\tcovered_bases\\tcoverage_breadth\\tmean_depth\\n' > all_hpv_types.tsv
    fi

    # Merge all transcript classification results.
    tx_files=( *_transcript_classes.tsv )
    if (( \${#tx_files[@]} > 0 )); then
        head -1 "\${tx_files[0]}" > all_transcript_classes.tsv
        for f in "\${tx_files[@]}"; do
            tail -n +2 "\$f" >> all_transcript_classes.tsv
        done
    else
        printf 'sample_id\\tgene\\tread_count\\tclass\\n' > all_transcript_classes.tsv
    fi

    # Generate report
    summarize_results.R \\
        --samplesheet samplesheet_slim.csv \\
        --raw-samplesheet samplesheet_raw.csv \\
        --hpv-types all_hpv_types.tsv \\
        --transcript-classes all_transcript_classes.tsv \\
        --hpv-status ${hpv_status_csv} \\
        --outdir summary_tables
    """
}
