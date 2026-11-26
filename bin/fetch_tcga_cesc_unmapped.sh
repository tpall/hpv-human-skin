#!/bin/bash
# Pull the UNMAPPED read fraction of TCGA-CESC RNA-seq BAMs and convert to FASTQ.
#
# Why this cohort: the vaginal SRA baseline we used as the cervical comparator is
# 78% integrated cell-line material (random-1000 subsample: 39%), because "cervix"
# in public SRA is saturated by HeLa. TCGA-CESC is the opposite situation --
# every library is documented primary tissue with a pathology diagnosis attached,
# so it can carry the "presence WITH activity in cervix" arm without a
# contamination caveat.
#
# Why unmapped only: GDC harmonized BAMs keep aligned AND unaligned reads, and the
# slicing API serves an arbitrary region -- including the keyword `unmapped`. HPV
# reads cannot align to a human-only reference, so they are all in that fraction.
# This turns a 2.27 TB download into a small one, and it is the only part we need.
#
# ACCESS: these are controlled (dbGaP phs000178). You need an approved Data Access
# Request and a GDC token. Without one the API returns 403 and this script stops.
# Everything upstream of the token -- the manifest, case IDs, sample types,
# diagnoses -- is open metadata and is already committed under assets/tcga/.
#
#   export GDC_TOKEN_FILE=~/.gdc/gdc-user-token.txt
#   bin/fetch_tcga_cesc_unmapped.sh [--limit N] [--outdir results_tcga_cesc]
set -euo pipefail

MANIFEST=${MANIFEST:-assets/tcga/tcga_cesc_genome_bams.tsv}
OUTDIR=results_tcga_cesc
LIMIT=0
while [[ $# -gt 0 ]]; do
  case $1 in
    --limit)  LIMIT=$2; shift 2 ;;
    --outdir) OUTDIR=$2; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

TOKFILE=${GDC_TOKEN_FILE:-$HOME/.gdc/gdc-user-token.txt}
[ -s "$TOKFILE" ] || { echo "No GDC token at $TOKFILE." >&2
  echo "Get one at portal.gdc.cancer.gov after the dbGaP phs000178 DAR is approved." >&2; exit 1; }
TOKEN=$(tr -d '\n' < "$TOKFILE")
[ -s "$MANIFEST" ] || { echo "missing manifest: $MANIFEST" >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not on PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"/{slices,fastq,logs}
SHEET=$OUTDIR/samplesheet.csv
[ -s "$SHEET" ] || echo "sample,fastq_1,fastq_2,case,sample_type,diagnosis" > "$SHEET"

n=0
# fields: file_id file_name file_size md5sum case sample sample_type diagnosis organ workflow
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r FID FNAME FSIZE MD5 CASE SAMPLE STYPE DX ORGAN WF; do
  [ -n "$FID" ] || continue
  [ "$LIMIT" -gt 0 ] && [ "$n" -ge "$LIMIT" ] && break
  n=$((n+1))
  ID="${SAMPLE:-$CASE}"
  R1=$OUTDIR/fastq/${ID}_1.fastq.gz
  [ -s "$R1" ] && { echo "  $ID already done"; continue; }

  BAM=$OUTDIR/slices/${ID}.unmapped.bam
  echo "[$(date +%T)] $ID  ($STYPE)  slicing unmapped from $FNAME"
  # 206 Partial Content is the success code for a slice; --fail catches 403/404.
  if ! curl -sS --fail --retry 3 --retry-delay 10 \
        --header "X-Auth-Token: $TOKEN" \
        "https://api.gdc.cancer.gov/slicing/view/${FID}?region=unmapped" \
        --output "$BAM" 2>>"$OUTDIR/logs/${ID}.err"; then
      echo "  !! slice failed for $ID (see logs/${ID}.err)"; rm -f "$BAM"; continue
  fi
  if ! samtools quickcheck "$BAM" 2>/dev/null; then
      echo "  !! truncated/invalid slice for $ID"; rm -f "$BAM"; continue
  fi

  # Unmapped mates are unordered in the slice; collate before fastq or pairing breaks.
  samtools collate -u -O "$BAM" "$OUTDIR/slices/${ID}.tmp" 2>/dev/null \
    | samtools fastq -N -0 /dev/null -s /dev/null \
        -1 "$OUTDIR/fastq/${ID}_1.fastq.gz" \
        -2 "$OUTDIR/fastq/${ID}_2.fastq.gz" - 2>>"$OUTDIR/logs/${ID}.err" || {
      echo "  !! fastq conversion failed for $ID"; continue; }

  SZ=$(du -sh "$OUTDIR/fastq/${ID}_1.fastq.gz" 2>/dev/null | cut -f1)
  echo "  unmapped R1: $SZ"
  printf '%s,%s,%s,%s,%s,"%s"\n' "$ID" \
    "$(cd "$OUTDIR" && pwd)/fastq/${ID}_1.fastq.gz" \
    "$(cd "$OUTDIR" && pwd)/fastq/${ID}_2.fastq.gz" \
    "$CASE" "$STYPE" "$DX" >> "$SHEET"
  rm -f "$BAM" "$OUTDIR/slices/${ID}.tmp"*
done

echo
echo "samplesheet: $SHEET  ($(( $(wc -l < "$SHEET") - 1 )) libraries)"
echo "next: run the typing pipeline on it, then"
echo "      bin/e1e2_integration_marker.R --agg $OUTDIR --no-tiers --outdir ${OUTDIR}_e1e2"
