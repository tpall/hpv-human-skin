#!/bin/bash
#SBATCH --job-name=hpv-junction
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.log
#
# Detect HPV integration from reads: virus-host chimeric junctions.
#
# Rationale: the assembly-based test failed (see
# dna_rna_triangulation/REPORT_e2_integration_negative.md) because an integrated
# genome does not assemble across the junction. Reads do span it. STAR's chimeric
# detection against a COMBINED human+HPV index calls exactly those reads, which is
# the standard route for virus-host fusion discovery in RNA-seq.
#
# Pre-registered prediction from the E1:E2 transcript marker:
#   cell lines (log2 E1/E2 ~ +10)  -> junctions expected
#   warts / EV  (log2 ~ -2 to -4)  -> none expected (episomal)
#   cSCC        (log2 ~ -3 to -6)  -> none expected, EXCEPT
#   SRR35133829 (log2 = +7.3)      -> junctions expected; this single library is
#                                     the marker's falsifiable call.
set -euo pipefail

PROJ=${SLURM_SUBMIT_DIR:-$HOME/Projects/hpv-human-skin}
WORK=$PROJ/results_junction_pilot
HUMAN=/gpfs/space/databases/HumanGenome_NCBI_RefSeq/GRCh38.15/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
HPV=$PROJ/assets/hpv_references/hpv_all.fasta
IDX=$WORK/star_index_human_hpv
SAMPLES=$WORK/junction_pilot.tsv

export PATH=$HOME/miniconda3/envs/build_refs/bin:$HOME/miniconda3/envs/durable_screen/bin:$PATH
mkdir -p "$WORK"/{fastq,align,junctions}
cd "$WORK"

# ── 1. combined reference (human + HPV panel as extra contigs) ───────────────
if [ ! -s "$IDX/SA" ]; then
  echo "[$(date +%T)] building combined human+HPV STAR index"
  mkdir -p "$IDX" ref
  [ -s ref/combined.fa ] || { zcat "$HUMAN" > ref/combined.fa; cat "$HPV" >> ref/combined.fa; }
  STAR --runMode genomeGenerate --runThreadN "${SLURM_CPUS_PER_TASK:-16}" \
       --genomeDir "$IDX" --genomeFastaFiles ref/combined.fa \
       --limitGenomeGenerateRAM 60000000000 --genomeSAsparseD 2
  echo "[$(date +%T)] index done"
fi

# HPV contig names, to tell viral from human sides of a junction
grep '^>' "$HPV" | sed 's/^>//; s/ .*//' | sort -u > hpv_contigs.txt
echo "HPV reference contigs: $(wc -l < hpv_contigs.txt)"

# ── 2. per library: fetch, align with chimeric detection, keep junctions ─────
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r SRR TIER READS LR; do
  [ -n "$SRR" ] || continue
  OUT=$WORK/junctions/${SRR}.junction.tsv
  [ -s "$OUT" ] && { echo "  $SRR already done"; continue; }
  echo "[$(date +%T)] $SRR  tier=$TIER  hpv_reads=$READS  log2(E1/E2)=$LR"

  fasterq-dump --split-3 --threads 8 --outdir fastq "$SRR" >/dev/null 2>&1 || {
      echo "  !! download failed: $SRR"; continue; }
  R1=$(ls fastq/${SRR}*_1.fastq 2>/dev/null || ls fastq/${SRR}.fastq 2>/dev/null || true)
  R2=$(ls fastq/${SRR}*_2.fastq 2>/dev/null || true)
  [ -n "$R1" ] || { echo "  !! no fastq: $SRR"; continue; }

  STAR --genomeDir "$IDX" --readFilesIn $R1 $R2 \
       --runThreadN "${SLURM_CPUS_PER_TASK:-16}" \
       --outFileNamePrefix align/${SRR}. \
       --outSAMtype None \
       --chimSegmentMin 20 --chimOutType Junctions \
       --chimJunctionOverhangMin 15 --chimScoreMin 1 \
       --outFilterMultimapNmax 20 >/dev/null 2>&1 || {
      echo "  !! STAR failed: $SRR"; rm -f fastq/${SRR}*; continue; }

  # keep only junctions with one side on an HPV contig and the other on human
  awk -v OFS='\t' 'NR>1 && $1!~"^#"' align/${SRR}.Chimeric.out.junction 2>/dev/null \
    | awk -v OFS='\t' 'NR==FNR{h[$1];next} {a=($1 in h); b=($4 in h);
        if (a!=b) print $1,$2,$4,$5,$7,$8,$9}' hpv_contigs.txt - \
    > "$OUT" || true
  echo "  virus-host junctions: $(wc -l < "$OUT")"
  rm -f fastq/${SRR}* align/${SRR}.Chimeric.out.junction
done

# ── 3. summary ──────────────────────────────────────────────────────────────
{
  printf "srr_id\ttier\thpv_reads\tlog2_e1e2\tn_junctions\tn_distinct_sites\n"
  tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r SRR TIER READS LR; do
    F=$WORK/junctions/${SRR}.junction.tsv
    if [ -s "$F" ]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SRR" "$TIER" "$READS" "$LR" \
        "$(wc -l < "$F")" "$(cut -f1,2 "$F" | sort -u | wc -l)"
    else
      printf "%s\t%s\t%s\t%s\t0\t0\n" "$SRR" "$TIER" "$READS" "$LR"
    fi
  done
} > "$WORK/junction_summary.tsv"
echo "[$(date +%T)] DONE"
cat "$WORK/junction_summary.tsv"
