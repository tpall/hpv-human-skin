/*
 * SRA Dataset Discovery
 *
 * Query NCBI GEO/SRA for relevant RNA-seq datasets and produce
 * an enriched samplesheet with tissue categories.
 */

process SRA_DISCOVERY {
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/biopython_pysradb:2.2.1--e2decf4de5e5eb58'
    publishDir "${params.outdir}/metadata", mode: params.publish_dir_mode

    input:
    val query_terms

    output:
    path "samplesheet_enriched.csv", emit: samplesheet
    path "samplesheet_raw.csv",      emit: raw_samplesheet

    script:
    def max_samples = params.max_samples > 0 ? "--max-samples ${params.max_samples}" : ""
    """
    query_sra.py \\
        --query "${query_terms}" \\
        --output samplesheet_raw.csv \\
        ${max_samples}

    parse_metadata.py \\
        --input samplesheet_raw.csv \\
        --categories ${params.tissue_categories} \\
        --output samplesheet_enriched.csv
    """
}
