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
    path samplesheet
    path raw_samplesheet
    path hpv_types_files
    path transcript_class_files
    path kraken_reports
    path hpv_status_csv

    output:
    path "hpv_skin_report.html",    emit: html
    path "summary_tables/",         emit: tables

    script:
    """
    mkdir -p summary_tables

    # Merge all HPV type results
    head -1 \$(ls *_hpv_types.tsv | head -1) > all_hpv_types.tsv 2>/dev/null || echo "sample_id\thpv_reference\tref_length\tread_count\tcovered_bases\tcoverage_breadth\tmean_depth" > all_hpv_types.tsv
    for f in *_hpv_types.tsv; do
        tail -n +2 "\$f" >> all_hpv_types.tsv
    done

    # Merge all transcript classification results
    head -1 \$(ls *_transcript_classes.tsv | head -1) > all_transcript_classes.tsv 2>/dev/null || echo "sample_id\tgene\tread_count\tclass" > all_transcript_classes.tsv
    for f in *_transcript_classes.tsv; do
        tail -n +2 "\$f" >> all_transcript_classes.tsv
    done

    # Generate report
    summarize_results.R \\
        --samplesheet ${samplesheet} \\
        --raw-samplesheet ${raw_samplesheet} \\
        --hpv-types all_hpv_types.tsv \\
        --transcript-classes all_transcript_classes.tsv \\
        --hpv-status ${hpv_status_csv} \\
        --outdir summary_tables
    """
}
