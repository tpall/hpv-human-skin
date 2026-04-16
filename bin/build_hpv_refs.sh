#!/usr/bin/env bash
set -euo pipefail

# Build HPV reference database from PaVE and NCBI RefSeq
#
# Usage: build_hpv_refs.sh <output_dir> [custom_fasta_dir]
#
# Downloads all HPV genomes from PaVE and RefSeq, merges, deduplicates,
# and creates indices for STAR alignment.

OUTPUT_DIR="${1:?Usage: build_hpv_refs.sh <output_dir> [custom_fasta_dir]}"
CUSTOM_DIR="${2:-}"
THREADS="${3:-8}"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "=== Downloading HPV genomes from PaVE ==="
# PaVE provides a bulk download of all papillomavirus reference genomes
PAVE_URL="https://pave.niaid.nih.gov/api/genomesequences/fasta"
curl -fsSL "${PAVE_URL}" -o pave_all_pv.fasta || {
    echo "Warning: PaVE bulk download failed, trying individual download..."
    # Fallback: download via the reference genomes page
    curl -fsSL "https://pave.niaid.nih.gov/api/hpv_ref_genomes/fasta" -o pave_all_pv.fasta
}

echo "=== Downloading HPV genomes from NCBI RefSeq ==="
# Search for Human papillomavirus complete genomes in RefSeq
esearch -db nucleotide -query '"Human papillomavirus"[Organism] AND refseq[filter] AND complete genome[Title]' \
    | efetch -format fasta > refseq_hpv.fasta || {
    echo "Warning: Entrez Direct not available. Skipping RefSeq download."
    echo "Install edirect: sh -c \"\$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)\""
    touch refseq_hpv.fasta
}

echo "=== Merging sequences ==="
cat pave_all_pv.fasta refseq_hpv.fasta > merged_raw.fasta

# Add custom sequences if provided
if [[ -n "${CUSTOM_DIR}" && -d "${CUSTOM_DIR}" ]]; then
    echo "=== Adding custom sequences from ${CUSTOM_DIR} ==="
    for f in "${CUSTOM_DIR}"/*.fasta "${CUSTOM_DIR}"/*.fa "${CUSTOM_DIR}"/*.fna; do
        [[ -f "$f" ]] && cat "$f" >> merged_raw.fasta
    done
fi

echo "=== Deduplicating by sequence ID ==="
# Use seqkit to remove duplicates by ID, keeping the first occurrence
if command -v seqkit &>/dev/null; then
    seqkit rmdup -n merged_raw.fasta -o hpv_all.fasta
else
    # Fallback: simple awk dedup by header
    awk '/^>/{id=$1; if(seen[id]++){skip=1}else{skip=0}} !skip' merged_raw.fasta > hpv_all.fasta
fi

echo "=== Creating samtools index ==="
samtools faidx hpv_all.fasta

echo "=== Building STAR genome index ==="
mkdir -p star_index
STAR --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles hpv_all.fasta \
     --genomeSAindexNbases 6 \
     --runThreadN "${THREADS}"

echo "=== Building HISAT2 index ==="
mkdir -p hisat2_index
hisat2-build -p "${THREADS}" hpv_all.fasta hisat2_index/hpv_all

# Count sequences
N_SEQ=$(grep -c '^>' hpv_all.fasta)
echo ""
echo "=== Done ==="
echo "Total HPV sequences: ${N_SEQ}"
echo "Reference FASTA:     ${OUTPUT_DIR}/hpv_all.fasta"
echo "STAR index:          ${OUTPUT_DIR}/star_index/"
echo "HISAT2 index:        ${OUTPUT_DIR}/hisat2_index/"

# Cleanup intermediate files
rm -f merged_raw.fasta
