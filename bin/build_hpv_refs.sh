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
# PaVE 2.0 no longer exposes a bulk FASTA endpoint. Instead we enumerate all
# human-host papillomavirus genomes via /api/genome and fetch each one from
# /api/genome/{id}, assembling the FASTA ourselves.
PAVE_API="https://pave.niaid.nih.gov/api"
python3 - "${PAVE_API}" > pave_all_pv.fasta <<'PY'
import json
import sys
import time
import urllib.request

api = sys.argv[1].rstrip("/")

def fetch_json(url, retries=3):
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")

print(f"Listing genomes from {api}/genome ...", file=sys.stderr)
listing = fetch_json(f"{api}/genome?limit=9999&includeNonRef=true")
genomes = [g for g in listing.get("data", [])
           if (g.get("host_common_name") or "").lower() == "human"]
print(f"Found {len(genomes)} human papillomavirus genomes", file=sys.stderr)

out = sys.stdout
for i, meta in enumerate(genomes, 1):
    gid = meta["locus_id"]
    if i % 25 == 0 or i == len(genomes):
        print(f"  [{i}/{len(genomes)}] {gid}", file=sys.stderr)
    rec = fetch_json(f"{api}/genome/{gid}")
    seq = (rec.get("itemSequence") or {}).get("sequence")
    if not seq:
        print(f"  warning: no sequence for {gid}", file=sys.stderr)
        continue
    accession = (rec.get("genome") or {}).get("accession", "")
    desc = (rec.get("genome") or {}).get("description", meta.get("virus_name", ""))
    out.write(f">{gid} {accession} {desc}\n")
    for j in range(0, len(seq), 70):
        out.write(seq[j:j+70] + "\n")
PY

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
