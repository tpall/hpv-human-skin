/*
 * Early vs Late Transcript Classification
 *
 * Counts reads mapping to HPV early (E1-E7) and late (L1, L2) gene
 * regions to identify productive viral infection.
 */

process TRANSCRIPT_CLASSIFY {
    tag "${meta.srr_id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/subread_python:3.12--cb32a078e13625af'
    publishDir "${params.outdir}/transcript_classification", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path hpv_gene_gff

    output:
    tuple val(meta), path("${meta.srr_id}_transcript_classes.tsv"), emit: classes
    tuple val(meta), env(HAS_LATE_TRANSCRIPTS),                     emit: late_status

    script:
    """
    # Count reads per gene region using featureCounts
    featureCounts \\
        -a ${hpv_gene_gff} \\
        -o ${meta.srr_id}_featurecounts.txt \\
        -F GFF \\
        -t gene \\
        -g gene_name \\
        -T ${task.cpus} \\
        --minOverlap 20 \\
        -Q 10 \\
        ${bam}

    # Parse featureCounts output and classify early vs late
    python3 - <<'PYEOF'
import csv
import sys

sample_id = "${meta.srr_id}"
min_late_reads = ${params.late_transcript_min_reads}

early_genes = {"E1", "E2", "E4", "E5", "E6", "E7"}
late_genes = {"L1", "L2"}

gene_counts = {}
with open(f"{sample_id}_featurecounts.txt") as f:
    for line in f:
        if line.startswith("#") or line.startswith("Geneid"):
            continue
        parts = line.strip().split("\\t")
        gene_name = parts[0]
        count = int(parts[-1])
        gene_counts[gene_name] = count

# Aggregate early and late
early_total = sum(gene_counts.get(g, 0) for g in early_genes)
late_total = sum(gene_counts.get(g, 0) for g in late_genes)
l1_count = gene_counts.get("L1", 0)
l2_count = gene_counts.get("L2", 0)

has_late = late_total >= min_late_reads

# Write results
with open(f"{sample_id}_transcript_classes.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\\t")
    writer.writerow(["sample_id", "gene", "read_count", "class"])
    for gene, count in sorted(gene_counts.items()):
        gene_class = "late" if gene in late_genes else ("early" if gene in early_genes else "other")
        writer.writerow([sample_id, gene, count, gene_class])
    # Summary row
    writer.writerow([sample_id, "EARLY_TOTAL", early_total, "summary"])
    writer.writerow([sample_id, "LATE_TOTAL", late_total, "summary"])
    writer.writerow([sample_id, "PRODUCTIVE_INFECTION", "yes" if has_late else "no", "summary"])

print(f"{sample_id}: early={early_total}, late={late_total} (L1={l1_count}, L2={l2_count}), productive={'yes' if has_late else 'no'}", file=sys.stderr)
PYEOF

    # Set environment variable for Nextflow
    LATE_TOTAL=\$(grep "LATE_TOTAL" ${meta.srr_id}_transcript_classes.tsv | cut -f3)
    if [ "\$LATE_TOTAL" -ge "${params.late_transcript_min_reads}" ]; then
        HAS_LATE_TRANSCRIPTS="yes"
    else
        HAS_LATE_TRANSCRIPTS="no"
    fi
    """
}
