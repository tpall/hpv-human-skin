#!/usr/bin/env bash
# Is the E2 ORF intact in the HPV DNA reservoir, and does it differ by body site?
#
# Integration of a papillomavirus genome characteristically disrupts E2, the
# transcriptional repressor of E6/E7. bin/e1e2_integration_marker.R shows that
# signature on the RNA side; this is the DNA-side complement, asking whether the
# contigs recovered by Shen et al.'s petabase survey carry an intact E2 -- and
# whether cutaneous beta/gamma (expected episomal) differs from vaginal alpha
# (expected to integrate in disease).
#
# Shen's tables carry no ORF annotation and no sequences, but the Zenodo deposit
# (10.5281/zenodo.19683048) holds the Logan contigs. Two caveats govern use:
#
#   RETRIEVAL IS NOT RANDOM. Only ~66% of contigs are in the deposit; the rest
#   carry Logan unitig-path identifiers ("Score=/Edges=") absent from BOTH S4
#   and S5. Retrieval varies by archive (DRR 41%, SRR 67%, ERR 75%), library
#   type (TRANSCRIPTOMIC 91% vs GENOMIC 62%), cluster status (blastn 46%) and
#   body site (skin 56%, vaginal 73%, gut 86%). If unitig-path reconstruction
#   was needed for fragmented assemblies, the missing set is enriched for
#   incomplete genomes and this analysis is biased toward completeness.
#
#   ASSEMBLY ASCERTAINMENT. An integrated genome need not assemble into a clean
#   papillomavirus contig, because the virus-host junction breaks the path. The
#   absolute rate of intact E2 is therefore an overestimate. The BETWEEN-TISSUE
#   contrast is the defensible readout: the same bias applies to both arms.
#
#   bin/e2_integrity_from_logan.sh <workdir>
set -euo pipefail
WORK="${1:?usage: e2_integrity_from_logan.sh <workdir>}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$WORK"; cd "$WORK"
BLAST="${BLAST_BIN:-$HOME/miniforge3/envs/hpvaln/bin}"

# ── 1. contig -> tissue map (ALL cluster members, not just centroids) ────────
REPO="$REPO" python3 "$REPO/bin/_e2_contig_map.py"

# ── 2. stream the Logan contigs, keep only ours (full file never stored) ─────
if [ ! -s contigs.fa ]; then
  curl -sL --max-time 1800 \
    "https://zenodo.org/records/19683048/files/S5.Logan_Outputs_R2.zst?download=1" \
    | zstd -d -c | python3 "$REPO/bin/_e2_extract.py"
fi

# ── 3. papillomavirus ORF proteins (RefSeq) as tblastn queries ──────────────
[ -s pv_orf_proteins.fa ] || python3 "$REPO/bin/_e2_fetch_orfs.py"

# ── 4. annotate and score E2 presence conditional on reaching E1 ────────────
"$BLAST/makeblastdb" -in contigs.fa -dbtype nucl -out contig_db >/dev/null
"$BLAST/tblastn" -query pv_orf_proteins.fa -db contig_db -evalue 1e-5 -num_threads 4 \
  -max_target_seqs 100000 \
  -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore" \
  -out orf_hits.tsv
REPO="$REPO" python3 "$REPO/bin/_e2_score.py"
