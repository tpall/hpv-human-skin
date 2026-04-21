/*
 * Kraken2 HPV Screening
 *
 * Fast taxonomic classification to identify samples containing HPV reads.
 * Extracts HPV-classified reads for downstream alignment.
 */

process KRAKEN2_SCREEN {
    tag "${meta.srr_id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/kraken2_seqtk:1.4--2093e3cfc1acab9f'
    publishDir "${params.outdir}/kraken2", mode: params.publish_dir_mode, pattern: "*_report.txt"

    input:
    tuple val(meta), path(reads)
    path kraken2_db

    output:
    tuple val(meta), path("*_hpv_reads*.fastq.gz"), path("${meta.srr_id}_kraken2_report.txt"), emit: hpv_reads
    tuple val(meta), env(HPV_STATUS), env(HPV_READ_COUNT),                                     emit: hpv_status
    path "${meta.srr_id}_kraken2_report.txt",                                                  emit: report

    script:
    def input_reads = meta.layout == "PAIRED"
        ? "--paired ${reads[0]} ${reads[1]}"
        : "${reads[0]}"
    """
    # Run Kraken2 classification
    kraken2 \\
        --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${meta.srr_id}_kraken2.out \\
        --report ${meta.srr_id}_kraken2_report.txt \\
        --classified-out ${meta.srr_id}_classified#.fastq \\
        ${input_reads}

    # Extract HPV-classified read IDs
    # Papillomaviridae taxid: 151340 (family), Alphapapillomavirus: 337043
    # Extract all reads classified under Papillomaviridae
    awk '\$3 ~ /^(151340|337043|337041|337042|337044|337045)/ || \$0 ~ /papilloma/ {print \$2}' \\
        ${meta.srr_id}_kraken2.out | sort -u > hpv_read_ids.txt

    HPV_READ_COUNT=\$(wc -l < hpv_read_ids.txt)

    if [ "\$HPV_READ_COUNT" -ge "${params.hpv_min_reads}" ]; then
        HPV_STATUS="HPV+"
    else
        HPV_STATUS="HPV-"
    fi

    # Extract HPV reads from FASTQ and drop malformed records (seq_len != qual_len).
    # SRA-downloaded FASTQs occasionally contain such records and STAR aborts on them.
    if [ "\$HPV_READ_COUNT" -gt 0 ]; then
        if [ "${meta.layout}" = "PAIRED" ]; then
            seqtk subseq ${reads[0]} hpv_read_ids.txt > r1.raw.fq
            seqtk subseq ${reads[1]} hpv_read_ids.txt > r2.raw.fq
            : > r1.clean.fq
            : > r2.clean.fq
            # Lockstep both mates through paste; drop the whole pair if either mate is malformed.
            paste r1.raw.fq r2.raw.fq | awk -F'\\t' '
                NR%4==1 { h1=\$1; h2=\$2 }
                NR%4==2 { s1=\$1; s2=\$2 }
                NR%4==3 { p1=\$1; p2=\$2 }
                NR%4==0 && length(s1)==length(\$1) && length(s2)==length(\$2) {
                    print h1"\\n"s1"\\n"p1"\\n"\$1 > "r1.clean.fq"
                    print h2"\\n"s2"\\n"p2"\\n"\$2 > "r2.clean.fq"
                }'
            gzip -c r1.clean.fq > ${meta.srr_id}_hpv_reads_R1.fastq.gz
            gzip -c r2.clean.fq > ${meta.srr_id}_hpv_reads_R2.fastq.gz
            rm -f r1.raw.fq r2.raw.fq r1.clean.fq r2.clean.fq
        else
            seqtk subseq ${reads[0]} hpv_read_ids.txt | awk '
                NR%4==1 { h=\$0 }
                NR%4==2 { s=\$0 }
                NR%4==3 { p=\$0 }
                NR%4==0 && length(s)==length(\$0) {
                    print h"\\n"s"\\n"p"\\n"\$0
                }' | gzip > ${meta.srr_id}_hpv_reads.fastq.gz
        fi
    else
        # Create empty files
        if [ "${meta.layout}" = "PAIRED" ]; then
            touch ${meta.srr_id}_hpv_reads_R1.fastq.gz ${meta.srr_id}_hpv_reads_R2.fastq.gz
        else
            touch ${meta.srr_id}_hpv_reads.fastq.gz
        fi
    fi

    # Cleanup large intermediate files
    rm -f ${meta.srr_id}_kraken2.out ${meta.srr_id}_classified*.fastq
    """
}
