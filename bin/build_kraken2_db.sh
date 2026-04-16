#!/usr/bin/env bash
set -euo pipefail

# Build custom Kraken2 database with human genome + all HPV genomes
#
# Usage: build_kraken2_db.sh <db_dir> <hpv_fasta> [threads]
#
# Prerequisites: kraken2, kraken2-build must be in PATH

DB_DIR="${1:?Usage: build_kraken2_db.sh <db_dir> <hpv_fasta> [threads]}"
HPV_FASTA="${2:?Must provide path to HPV reference FASTA}"
THREADS="${3:-8}"

echo "=== Building Kraken2 database ==="
echo "Database dir: ${DB_DIR}"
echo "HPV FASTA:    ${HPV_FASTA}"
echo "Threads:      ${THREADS}"

mkdir -p "${DB_DIR}"

echo "=== Downloading NCBI taxonomy ==="
kraken2-build --download-taxonomy --db "${DB_DIR}"

echo "=== Downloading human library (GRCh38) ==="
kraken2-build --download-library human --db "${DB_DIR}"

echo "=== Adding HPV sequences to library ==="
# Add each HPV sequence with proper taxonomy
# Kraken2 needs taxid in the FASTA header: >seq_id|kraken:taxid|TAXID
# Papillomaviridae taxid = 151340, but individual HPV types have specific taxids
# We add them under the Papillomaviridae family for classification

# Create a temporary FASTA with kraken-compatible headers
TMP_HPV="${DB_DIR}/hpv_kraken_format.fasta"
awk '
/^>/ {
    # Extract sequence ID
    split($0, a, " ")
    id = substr(a[1], 2)
    # Use Alphapapillomavirus taxid (337043) as default
    # This ensures HPV reads are classified within Papillomaviridae
    printf ">%s|kraken:taxid|337043 %s\n", id, substr($0, 2)
    next
}
{ print }
' "${HPV_FASTA}" > "${TMP_HPV}"

kraken2-build --add-to-library "${TMP_HPV}" --db "${DB_DIR}"
rm -f "${TMP_HPV}"

echo "=== Building database (this may take several hours) ==="
kraken2-build --build --db "${DB_DIR}" --threads "${THREADS}"

echo "=== Cleaning up intermediate files ==="
kraken2-build --clean --db "${DB_DIR}"

# Verify
echo ""
echo "=== Done ==="
echo "Database location: ${DB_DIR}"
echo "Database contents:"
ls -lh "${DB_DIR}"/*.k2d 2>/dev/null || echo "Warning: .k2d files not found"
