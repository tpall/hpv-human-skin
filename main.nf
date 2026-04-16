#!/usr/bin/env nextflow

/*
 * HPV in Human Skin Pipeline
 *
 * Identifies HPV types and their prevalence across human tissues
 * using public RNA-seq data from NCBI GEO/SRA.
 *
 * Research questions:
 *   1. Which HPV types prevail in healthy skin?
 *   2. Which HPV types prevail in different skin pathologies?
 *   3. In which non-traditional tissues can HPV transcripts be found?
 *   4. Which transcriptomes show productive viral infection (late transcripts)?
 */

nextflow.enable.dsl = 2

// Validate essential parameters
if (!params.samplesheet && !params.sra_query_terms) {
    error "Either --samplesheet or --sra_query_terms must be provided"
}

// Log pipeline info
log.info """
    ╔═══════════════════════════════════════════╗
    ║     HPV in Human Skin Pipeline v0.1.0     ║
    ╚═══════════════════════════════════════════╝

    Samplesheet     : ${params.samplesheet ?: 'Auto-discovery via SRA'}
    SRA query terms  : ${params.sra_query_terms}
    HPV references   : ${params.hpv_ref_dir}
    Kraken2 DB       : ${params.kraken2_db}
    HPV min reads    : ${params.hpv_min_reads}
    HPV min coverage : ${params.hpv_min_coverage}
    Output dir       : ${params.outdir}
    """.stripIndent()

// Include the main workflow
include { HPV_SKIN } from './workflows/hpv_skin'

workflow {
    HPV_SKIN()
}

workflow.onComplete {
    log.info(
        workflow.success
            ? "\nPipeline completed successfully!\nResults: ${params.outdir}\n"
            : "\nPipeline failed. Check .nextflow.log for details.\n"
    )
}
