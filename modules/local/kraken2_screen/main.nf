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
    tuple val(meta), env('HPV_STATUS'), env('HPV_READ_COUNT'),                                 emit: hpv_status
    path "${meta.srr_id}_kraken2_report.txt",                                                  emit: report

    script:
    def input_reads = meta.layout == "PAIRED"
        ? "--paired ${reads[0]} ${reads[1]}"
        : "${reads[0]}"
    """
    # Prefer the conda env's binaries over the cluster's spack stack on PATH.
    export PATH="\${CONDA_PREFIX:+\$CONDA_PREFIX/bin:}\$PATH"

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

    # Extract HPV reads and sanitise. STAR aborts on any malformed record
    # (seq_len != qual_len, missing/extra lines, truncated tail). Parse each
    # mate independently with a resyncing awk that only starts a record at
    # a line beginning with '@', then intersect read IDs across mates so
    # the emitted pair files stay in lockstep. A line-counted (NR%4) filter
    # is unsafe here: one missing line anywhere frame-shifts the rest of
    # the file and lets garbage through.
    if [ "\$HPV_READ_COUNT" -gt 0 ]; then
        cat > sanitise.awk <<'AWK'
state==0 && substr(\$0,1,1)=="@" { h=\$0; state=1; next }
state==1 { s=\$0; state=2; next }
state==2 { p=\$0; state=3; next }
state==3 {
    if (substr(p,1,1)=="+" && length(s)==length(\$0)) {
        print h; print s; print p; print \$0
    }
    state=0
    next
}
AWK

        if [ "${meta.layout}" = "PAIRED" ]; then
            seqtk subseq ${reads[0]} hpv_read_ids.txt | awk -f sanitise.awk > r1.valid.fq
            seqtk subseq ${reads[1]} hpv_read_ids.txt | awk -f sanitise.awk > r2.valid.fq

            # Intersect IDs (first whitespace token after '@') so mates stay paired.
            awk 'NR%4==1 { sub(/^@/,""); print \$1 }' r1.valid.fq | sort -u > r1.ids
            awk 'NR%4==1 { sub(/^@/,""); print \$1 }' r2.valid.fq | sort -u > r2.ids
            comm -12 r1.ids r2.ids > common.ids

            seqtk subseq r1.valid.fq common.ids | gzip -c > ${meta.srr_id}_hpv_reads_R1.fastq.gz
            seqtk subseq r2.valid.fq common.ids | gzip -c > ${meta.srr_id}_hpv_reads_R2.fastq.gz
            rm -f r1.valid.fq r2.valid.fq r1.ids r2.ids common.ids
        else
            seqtk subseq ${reads[0]} hpv_read_ids.txt \\
                | awk -f sanitise.awk \\
                | gzip -c > ${meta.srr_id}_hpv_reads.fastq.gz
        fi
        rm -f sanitise.awk
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
