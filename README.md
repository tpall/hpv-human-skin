# HPV in Human Skin Pipeline

Nextflow DSL2 pipeline for identifying HPV types and their prevalence across human tissues using public RNA-seq data from NCBI GEO/SRA.

## Research Questions

1. Millised HPV tüübid prevalveerivad terves nahas?
2. Millised HPV tüübid prevalveerivad erinevates naha patoloogiate puhul?
3. Millistes mittetraditsioonilistest kudedes on võimalik leida HPV transkripte?
4. Millistes transkriptoomides esineb viiruse produktiivne nakkusfaas (hilised transkriptid)?

## Pipeline Overview

```
NCBI GEO/SRA → Metadata + FASTQ download → QC (fastp) → Kraken2 HPV screen
  → HPV+ samples: STAR alignment → HPV typing → Early/Late transcript classification
  → HPV- samples: Calculate HPV-negative rates per tissue
  → Report: Summary tables + figures
```

## Prerequisites

- Nextflow >= 23.04
- **One of:** Conda/Mamba, Singularity, or Docker
- SLURM cluster (or local execution)

Each pipeline process uses its own isolated environment (see `modules/local/*/environment.yml`). No monolithic environment needed — Nextflow resolves dependencies per-process. Standalone prep scripts have their own envs in `envs/`.

## Setup

```bash
# 1. Install Nextflow (if not already available)
curl -s https://get.nextflow.io | bash

# 2. Build reference databases (creates conda envs automatically)
sbatch setup.sh          # submit to SLURM (8 CPUs, 64GB, 24h)
# or: bash setup.sh 8    # run interactively
```

## Usage

Combine one **software profile** (`conda`, `singularity`, or `docker`) with an optional **executor profile** (`slurm` or `local`):

```bash
# --- SLURM cluster ---
# Conda environments (created per-process automatically)
nextflow run main.nf -profile conda,slurm --outdir results

# Singularity containers (recommended for HPC — no root needed)
nextflow run main.nf -profile singularity,slurm --outdir results

# --- Local machine ---
# Docker containers
nextflow run main.nf -profile docker --samplesheet my_samples.csv --outdir results

# Conda environments locally
nextflow run main.nf -profile conda --samplesheet my_samples.csv --outdir results

# --- Test run (small dataset) ---
nextflow run main.nf -profile test,conda,slurm

# --- Dry run (preview) ---
nextflow run main.nf -profile conda,slurm -preview
```

## Large runs (chunked execution)

For sample counts in the thousands (e.g. the full SRA query), a single Nextflow run is a poor fit — intermediate FASTQs in `work/` can reach tens of TB. Use the chunked driver instead, which processes the samplesheet in fixed-size batches, wipes each chunk's work dir once summaries are aggregated, and emits one combined report at the end.

```bash
# 1. Build the samplesheet once (no downloads yet)
nextflow run main.nf -profile conda,slurm -entry SRA_DISCOVERY_ONLY --outdir results_full

# 2. Run chunked (default: 100 samples per chunk)
sbatch bin/run_chunked.sh results_full/metadata/samplesheet_enriched.csv 100 results_full work_full
```

A failed chunk is not rerun automatically — its work dir is left in place for inspection, and a rerun of the driver will resume at that chunk. Per-sample failures within a chunk (withdrawn SRR, network blip) are tolerated by fail-soft `errorStrategy` on `SRA_DOWNLOAD` / `FASTP` / `KRAKEN2_SCREEN`.

## Samplesheet Format

CSV with columns: `srr_id, srx_id, study, tissue_category, diagnosis, layout`

```csv
srr_id,srx_id,study,tissue_category,diagnosis,layout
SRR1234567,SRX123456,SRP123456,nahk,normal,PAIRED
SRR2345678,SRX234567,SRP234567,nahk,wart,SINGLE
```

Tissue categories: `nahk` (skin), `anogenitaal`, `suuoos` (oral), `muu` (other)

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | null | Pre-built samplesheet (skip discovery) |
| `--sra_query_terms` | "skin RNA-seq Homo sapiens" | Search terms for SRA discovery |
| `--hpv_min_reads` | 10 | Min Kraken2 HPV reads for HPV+ call |
| `--hpv_min_coverage` | 0.10 | Min coverage breadth for type assignment |
| `--hpv_min_depth` | 2 | Min mean depth for type assignment |
| `--late_transcript_min_reads` | 3 | Min L1/L2 reads for productive infection |
| `--max_samples` | 0 | Limit number of samples (0 = all) |
| `--outdir` | results | Output directory |

## Output

```
results/
├── metadata/                  # Samplesheet and metadata
├── qc/fastp/                  # QC reports
├── kraken2/                   # Kraken2 classification reports
├── alignments/                # HPV-aligned BAMs
├── hpv_typing/                # Per-sample HPV type assignments
├── transcript_classification/ # Early/late transcript counts
├── report/
│   ├── hpv_skin_report.html   # Main HTML report
│   └── summary_tables/
│       ├── table1_hpv_healthy_skin.tsv
│       ├── table2_hpv_skin_pathology.tsv
│       ├── table3_hpv_nontraditional_tissues.tsv
│       ├── table4_productive_infection.tsv
│       ├── table5_hpv_negative_rates.tsv
│       ├── heatmap_hpv_tissue.png
│       └── barplot_hpv_status.png
└── pipeline_info/             # Nextflow execution reports
```

## Original Description

Meil oleks vaja teada millised HPV tüübid tegelikult millise sagedusega inimese nahas on. Probleem on vist selles, et enamustest transkriptoomidest visatakse mitte inimese transkriptid juba alguses välja ja tuleb minna tagasi raw data juurde.
