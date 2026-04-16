/*
 * STAR Alignment to HPV References
 *
 * Aligns HPV-classified reads against the full HPV reference panel
 * for precise type determination.
 */

process STAR_ALIGN {
    tag "${meta.srr_id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/star_samtools:2.7.11b--a43e00dda8a86e01'
    publishDir "${params.outdir}/alignments", mode: params.publish_dir_mode, pattern: "*.bam*"

    input:
    tuple val(meta), path(hpv_reads), path(kraken_report)
    path star_index

    output:
    tuple val(meta), path("${meta.srr_id}_Aligned.sortedByCoord.out.bam"),
                     path("${meta.srr_id}_Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(meta), path("${meta.srr_id}_Log.final.out"),                     emit: log

    script:
    def read_files = meta.layout == "PAIRED"
        ? "${hpv_reads[0]} ${hpv_reads[1]}"
        : "${hpv_reads[0]}"
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${read_files} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${meta.srr_id}_ \\
        --outSAMattributes NH HI AS NM MD \\
        --outFilterMultimapNmax 50 \\
        --outFilterMismatchNmax 3 \\
        --alignIntronMax 1 \\
        --alignEndsType EndToEnd

    # Index BAM
    samtools index ${meta.srr_id}_Aligned.sortedByCoord.out.bam
    """
}
