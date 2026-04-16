/*
 * HPV Type Assignment
 *
 * Determines HPV type(s) per sample based on alignment coverage
 * breadth and depth against the HPV reference panel.
 */

process HPV_TYPING {
    tag "${meta.srr_id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/pysam:0.22.1--py312hcfdcdd7_2'
    publishDir "${params.outdir}/hpv_typing", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path hpv_ref_fasta
    path hpv_ref_fai

    output:
    tuple val(meta), path("${meta.srr_id}_hpv_types.tsv"), emit: types
    path "${meta.srr_id}_hpv_coverage.tsv",                emit: coverage

    script:
    """
    #!/usr/bin/env python3
    import pysam
    import csv
    import sys

    bam_file = "${bam}"
    ref_fai = "${hpv_ref_fai}"
    sample_id = "${meta.srr_id}"
    min_coverage = ${params.hpv_min_coverage}
    min_depth = ${params.hpv_min_depth}

    # Load reference lengths
    ref_lengths = {}
    with open(ref_fai) as f:
        for line in f:
            parts = line.strip().split("\\t")
            ref_lengths[parts[0]] = int(parts[1])

    # Calculate coverage per reference
    results = []
    bam = pysam.AlignmentFile(bam_file, "rb")

    for ref_name, ref_len in ref_lengths.items():
        # Count aligned reads and coverage
        read_count = 0
        covered_positions = set()
        total_depth = 0

        for pileup_col in bam.pileup(ref_name, min_mapping_quality=10):
            covered_positions.add(pileup_col.pos)
            total_depth += pileup_col.n

        for read in bam.fetch(ref_name):
            if not read.is_unmapped and read.mapping_quality >= 10:
                read_count += 1

        coverage_breadth = len(covered_positions) / ref_len if ref_len > 0 else 0
        mean_depth = total_depth / ref_len if ref_len > 0 else 0

        results.append({
            "sample_id": sample_id,
            "hpv_reference": ref_name,
            "ref_length": ref_len,
            "read_count": read_count,
            "covered_bases": len(covered_positions),
            "coverage_breadth": round(coverage_breadth, 4),
            "mean_depth": round(mean_depth, 2),
        })

    bam.close()

    # Write full coverage report
    with open(f"{sample_id}_hpv_coverage.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter="\\t")
        writer.writeheader()
        writer.writerows(results)

    # Filter to assigned types
    assigned = [r for r in results
                if r["coverage_breadth"] >= min_coverage and r["mean_depth"] >= min_depth]

    # Sort by coverage breadth descending
    assigned.sort(key=lambda x: x["coverage_breadth"], reverse=True)

    with open(f"{sample_id}_hpv_types.tsv", "w", newline="") as f:
        if assigned:
            writer = csv.DictWriter(f, fieldnames=assigned[0].keys(), delimiter="\\t")
            writer.writeheader()
            writer.writerows(assigned)
        else:
            f.write("sample_id\\thpv_reference\\tref_length\\tread_count\\tcovered_bases\\tcoverage_breadth\\tmean_depth\\n")

    n_types = len(assigned)
    print(f"{sample_id}: {n_types} HPV type(s) assigned", file=sys.stderr)
    if assigned:
        for a in assigned[:5]:
            print(f"  {a['hpv_reference']}: breadth={a['coverage_breadth']}, depth={a['mean_depth']}", file=sys.stderr)
    """
}
