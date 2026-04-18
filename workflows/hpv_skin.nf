/*
 * HPV in Human Skin — Main Workflow
 *
 * Orchestrates the full pipeline from data discovery through
 * HPV typing and report generation.
 */

include { SRA_DISCOVERY        } from '../modules/local/sra_discovery/main'
include { SRA_DOWNLOAD         } from '../modules/local/sra_download/main'
include { FASTP                } from '../modules/local/fastp/main'
include { KRAKEN2_SCREEN       } from '../modules/local/kraken2_screen/main'
include { STAR_ALIGN           } from '../modules/local/star_align/main'
include { HPV_TYPING           } from '../modules/local/hpv_typing/main'
include { TRANSCRIPT_CLASSIFY  } from '../modules/local/transcript_classify/main'
include { REPORT               } from '../modules/local/report/main'

workflow HPV_SKIN {

    // ── Step 1: Get samplesheet ────────────────────────────────────────
    if (params.samplesheet) {
        // Use pre-built samplesheet — no raw/free-text companion, so reuse
        // the same file for the R report's Table 3 join.
        ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
        ch_raw_samplesheet = ch_samplesheet
    } else {
        // Discover datasets from SRA
        SRA_DISCOVERY(params.sra_query_terms)
        ch_samplesheet = SRA_DISCOVERY.out.samplesheet
        ch_raw_samplesheet = SRA_DISCOVERY.out.raw_samplesheet
    }
    ch_samples = ch_samplesheet
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                srr_id:          row.srr_id,
                srx_id:          row.srx_id ?: '',
                study:           row.study ?: '',
                tissue_category: row.tissue_category ?: 'muu',
                diagnosis:       row.diagnosis ?: 'unspecified',
                layout:          row.layout ?: 'SINGLE',
            ]
            meta
        }

    // ── Step 2: Download FASTQ ─────────────────────────────────────────
    SRA_DOWNLOAD(ch_samples)

    // ── Step 3: QC ─────────────────────────────────────────────────────
    FASTP(SRA_DOWNLOAD.out.reads)

    // ── Step 4: Kraken2 HPV screening ──────────────────────────────────
    ch_kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists: true)
    KRAKEN2_SCREEN(FASTP.out.reads, ch_kraken2_db.first())

    // ── Step 5: Filter to HPV+ samples ─────────────────────────────────
    // Collect HPV status for all samples (for reporting HPV- rates)
    ch_hpv_status = KRAKEN2_SCREEN.out.hpv_status
        .map { meta, status, count ->
            "${meta.srr_id}\t${meta.tissue_category}\t${meta.diagnosis}\t${status}\t${count}"
        }
        .collectFile(name: 'hpv_status.tsv', newLine: true,
                     storeDir: "${params.outdir}/metadata",
                     seed: "srr_id\ttissue_category\tdiagnosis\thpv_status\thpv_read_count")

    // Only HPV+ samples proceed to alignment
    ch_hpv_positive = KRAKEN2_SCREEN.out.hpv_reads
        .join(KRAKEN2_SCREEN.out.hpv_status)
        .filter { meta, reads, report, status, count ->
            status == "HPV+"
        }
        .map { meta, reads, report, status, count ->
            [meta, reads, report]
        }

    // ── Step 6: Align HPV reads ────────────────────────────────────────
    ch_star_index = Channel.fromPath("${params.hpv_ref_dir}/star_index", checkIfExists: true)
    STAR_ALIGN(ch_hpv_positive, ch_star_index.first())

    // ── Step 7: HPV typing ─────────────────────────────────────────────
    ch_hpv_ref = Channel.fromPath("${params.hpv_ref_dir}/hpv_all.fasta", checkIfExists: true)
    ch_hpv_fai = Channel.fromPath("${params.hpv_ref_dir}/hpv_all.fasta.fai", checkIfExists: true)
    HPV_TYPING(STAR_ALIGN.out.bam, ch_hpv_ref.first(), ch_hpv_fai.first())

    // ── Step 8: Early vs Late transcript classification ────────────────
    ch_gene_gff = Channel.fromPath(params.hpv_gene_gff, checkIfExists: true)
    TRANSCRIPT_CLASSIFY(STAR_ALIGN.out.bam, ch_gene_gff.first())

    // ── Step 9: Generate report ────────────────────────────────────────
    ch_all_types = HPV_TYPING.out.types
        .map { meta, tsv -> tsv }
        .collect()

    ch_all_classes = TRANSCRIPT_CLASSIFY.out.classes
        .map { meta, tsv -> tsv }
        .collect()

    ch_all_kraken = KRAKEN2_SCREEN.out.report.collect()

    if (!params.skip_report) {
        REPORT(
            ch_samplesheet,
            ch_raw_samplesheet,
            ch_all_types,
            ch_all_classes,
            ch_all_kraken,
            ch_hpv_status
        )
    }
}
