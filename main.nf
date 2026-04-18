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
include { HPV_SKIN      } from './workflows/hpv_skin'
include { REPORT        } from './modules/local/report/main'
include { SRA_DISCOVERY } from './modules/local/sra_discovery/main'

workflow {
    HPV_SKIN()
}

/*
 * SRA_DISCOVERY_ONLY — run only the SRA query / samplesheet build step,
 * so the chunked driver can split the resulting CSV before running the
 * full pipeline.
 *
 *   nextflow run main.nf -entry SRA_DISCOVERY_ONLY \
 *     --sra_query_terms "..." --outdir results
 */
workflow SRA_DISCOVERY_ONLY {
    SRA_DISCOVERY(params.sra_query_terms)
}

/*
 * REPORT_ONLY — final aggregation step for chunked runs.
 *
 * The driver (bin/run_chunked.sh) collects per-chunk summary TSVs into
 * params.agg_dir and invokes this entry with:
 *   nextflow run main.nf -entry REPORT_ONLY \
 *     --samplesheet <path> --agg_dir <path> --outdir <path>
 */
workflow REPORT_ONLY {
    if (!params.agg_dir) {
        error "REPORT_ONLY requires --agg_dir pointing at aggregated per-chunk TSVs"
    }
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
    ch_hpv_types   = Channel.fromPath("${params.agg_dir}/*_hpv_types.tsv").collect()
    ch_tx_classes  = Channel.fromPath("${params.agg_dir}/*_transcript_classes.tsv").collect()
    ch_kraken      = Channel.fromPath("${params.agg_dir}/*_kraken2_report.txt").collect()
    ch_hpv_status  = Channel.fromPath("${params.agg_dir}/hpv_status.tsv", checkIfExists: true)

    // No separate raw samplesheet in chunked-aggregation mode; reuse the
    // slim one. The R report uses any_of() so missing free-text columns
    // (title/tissue_source) are tolerated — Table 3 just skips them.
    REPORT(ch_samplesheet, ch_samplesheet, ch_hpv_types, ch_tx_classes, ch_kraken, ch_hpv_status)
}

workflow.onComplete {
    log.info(
        workflow.success
            ? "\nPipeline completed successfully!\nResults: ${params.outdir}\n"
            : "\nPipeline failed. Check .nextflow.log for details.\n"
    )
}
