#!/usr/bin/env Rscript

# HPV in Human Skin — Summary Tables and Report
#
# Generates summary tables answering the four research questions and
# renders an HTML report.
#
# Usage:
#   summarize_results.R \
#     --samplesheet samplesheet_enriched.csv \
#     --hpv-types all_hpv_types.tsv \
#     --transcript-classes all_transcript_classes.tsv \
#     --hpv-status hpv_status.tsv \
#     --outdir summary_tables

suppressPackageStartupMessages({
  library(tidyverse)
  library(rmarkdown)
  library(knitr)
  library(optparse)
})

# ── Parse arguments ─────────────────────────────────────────────────────
option_list <- list(
  make_option("--samplesheet", type = "character", help = "Enriched samplesheet CSV"),
  make_option("--raw-samplesheet", type = "character", default = NULL,
              help = "Raw samplesheet CSV (with title/tissue_source free-text columns)"),
  make_option("--hpv-types", type = "character", help = "Merged HPV types TSV"),
  make_option("--transcript-classes", type = "character", help = "Merged transcript classes TSV"),
  make_option("--hpv-status", type = "character", help = "HPV status TSV (all samples)"),
  make_option("--late-min-reads", type = "integer", default = 3,
              help = "Min L1 (major capsid) reads to call productive infection"),
  make_option("--outdir", type = "character", default = "summary_tables", help = "Output directory")
)
opts <- parse_args(OptionParser(option_list = option_list))

outdir <- opts$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ───────────────────────────────────────────────────────────
samples <- read_csv(opts$samplesheet, show_col_types = FALSE)
# Slim samplesheets written before commit db18557 lack is_cell_line; default
# to "false" so older runs still render (every sample treated as non-line).
if (!"is_cell_line" %in% names(samples)) {
  samples$is_cell_line <- "false"
}
samples <- samples %>%
  mutate(is_cell_line = tolower(as.character(is_cell_line)) == "true")
# Raw samplesheet holds comma-prone free-text columns (title, tissue_source)
# kept out of the slim samplesheet so Nextflow's splitCsv can't mis-split.
raw_path <- opts$`raw-samplesheet`
if (!is.null(raw_path) && file.exists(raw_path)) {
  samples_raw <- read_csv(raw_path, show_col_types = FALSE)
} else {
  samples_raw <- tibble(srr_id = character(), title = character(),
                        tissue_source = character())
}
hpv_types <- read_tsv(opts$`hpv-types`, show_col_types = FALSE)
transcripts <- read_tsv(opts$`transcript-classes`, show_col_types = FALSE)
hpv_status <- read_tsv(opts$`hpv-status`, show_col_types = FALSE)

cat(sprintf("Loaded: %d samples, %d HPV type assignments, %d transcript records\n",
            nrow(samples), nrow(hpv_types), nrow(transcripts)))

# ── Merge metadata with results ─────────────────────────────────────────
# Join tissue category, diagnosis, and is_cell_line to HPV results.
# Tables 1-2 deliberately exclude cell lines so prevalence reflects
# real tissue; Tables 3-4 keep them visible since they're useful as
# positive controls.
hpv_full <- hpv_types %>%
  left_join(
    samples %>% select(srr_id, tissue_category, diagnosis, is_cell_line),
    by = c("sample_id" = "srr_id")
  )

# ── Table 1: HPV type prevalence in healthy skin (tissue only) ─────────
table1 <- hpv_full %>%
  filter(tissue_category == "nahk", !is_cell_line) %>%
  filter(str_detect(tolower(diagnosis), "normal|healthy|control") |
         diagnosis == "unspecified") %>%
  group_by(hpv_reference) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    mean_coverage = mean(coverage_breadth, na.rm = TRUE),
    mean_depth = mean(mean_depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

write_tsv(table1, file.path(outdir, "table1_hpv_healthy_skin.tsv"))
cat(sprintf("Table 1: %d HPV types in healthy skin\n", nrow(table1)))

# ── Table 2: HPV type prevalence by skin pathology (tissue only) ───────
table2 <- hpv_full %>%
  filter(tissue_category == "nahk", !is_cell_line) %>%
  group_by(diagnosis, hpv_reference) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    mean_coverage = mean(coverage_breadth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(diagnosis, desc(n_samples))

write_tsv(table2, file.path(outdir, "table2_hpv_skin_pathology.tsv"))
cat(sprintf("Table 2: %d rows (HPV type x pathology)\n", nrow(table2)))

# ── Table 3: HPV in non-traditional tissues ─────────────────────────────
table3 <- hpv_full %>%
  filter(tissue_category == "muu") %>%
  group_by(hpv_reference, sample_id) %>%
  slice_max(coverage_breadth, n = 1) %>%
  ungroup() %>%
  left_join(
    samples_raw %>% select(any_of(c("srr_id", "title", "tissue_source"))),
    by = c("sample_id" = "srr_id")
  ) %>%
  # Raw samplesheet may be slim in chunked-aggregation mode; materialise
  # the display columns as NA so the final select() always succeeds.
  mutate(title = if ("title" %in% names(.)) title else NA_character_,
         tissue_source = if ("tissue_source" %in% names(.)) tissue_source else NA_character_) %>%
  select(sample_id, is_cell_line, tissue_source, title,
         hpv_reference, coverage_breadth, mean_depth) %>%
  arrange(desc(coverage_breadth))

write_tsv(table3, file.path(outdir, "table3_hpv_nontraditional_tissues.tsv"))
cat(sprintf("Table 3: %d HPV+ samples in non-traditional tissues\n", nrow(table3)))

# ── Table 4: Productive infection (L1 major-capsid positive) ───────────
# Derive the productive call from the L1 read count rather than trusting the
# pipeline's PRODUCTIVE_INFECTION summary row: (a) older chunks wrote that row
# under an L1+L2 rule that over-calls (L2-only signal isn't productive — virion
# assembly needs the major capsid L1), so recomputing here keeps the whole
# dataset uniform regardless of which module version produced each chunk.
late_min <- opts$`late-min-reads`
productive <- transcripts %>%
  filter(gene == "L1") %>%
  mutate(l1_reads = suppressWarnings(as.numeric(read_count))) %>%
  filter(!is.na(l1_reads) & l1_reads >= late_min) %>%
  distinct(sample_id)
cat(sprintf("Productive (L1 >= %d reads): %d samples\n", late_min, nrow(productive)))

table4 <- productive %>%
  left_join(
    samples %>% select(srr_id, tissue_category, diagnosis, is_cell_line),
    by = c("sample_id" = "srr_id")
  ) %>%
  left_join(
    transcripts %>%
      filter(gene %in% c("L1", "L2")) %>%
      pivot_wider(id_cols = sample_id, names_from = gene,
                  values_from = read_count, names_prefix = "reads_"),
    by = "sample_id"
  ) %>%
  left_join(
    hpv_types %>%
      group_by(sample_id) %>%
      slice_max(coverage_breadth, n = 1) %>%
      ungroup() %>%
      select(sample_id, top_hpv_type = hpv_reference),
    by = "sample_id"
  ) %>%
  # When there are zero HPV+ samples the pivot_wider and hpv_types joins
  # yield no reads_L1 / reads_L2 / top_hpv_type columns; materialise them
  # as NA so the final select() always succeeds.
  mutate(
    reads_L1     = if ("reads_L1"     %in% names(.)) reads_L1     else NA_character_,
    reads_L2     = if ("reads_L2"     %in% names(.)) reads_L2     else NA_character_,
    top_hpv_type = if ("top_hpv_type" %in% names(.)) top_hpv_type else NA_character_
  ) %>%
  select(sample_id, tissue_category, is_cell_line, diagnosis,
         top_hpv_type, reads_L1, reads_L2) %>%
  arrange(tissue_category, is_cell_line,
          desc(suppressWarnings(as.numeric(reads_L1))))

write_tsv(table4, file.path(outdir, "table4_productive_infection.tsv"))
cat(sprintf("Table 4: %d samples with productive infection\n", nrow(table4)))

# ── Table 5: HPV-negative rates per tissue category ────────────────────
# hpv_status.tsv carries srr_id/tissue_category/diagnosis/hpv_status/count
# but no is_cell_line column — join from the slim samplesheet so the rates
# can be split below.
hpv_status_enriched <- hpv_status %>%
  left_join(samples %>% select(srr_id, is_cell_line),
            by = c("srr_id" = "srr_id")) %>%
  mutate(is_cell_line = ifelse(is.na(is_cell_line), FALSE, is_cell_line))

table5 <- hpv_status_enriched %>%
  group_by(tissue_category) %>%
  summarise(
    total_samples = n(),
    hpv_positive = sum(hpv_status == "HPV+"),
    hpv_negative = sum(hpv_status == "HPV-"),
    pct_positive = round(100 * hpv_positive / total_samples, 1),
    pct_negative = round(100 * hpv_negative / total_samples, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_positive))

write_tsv(table5, file.path(outdir, "table5_hpv_negative_rates.tsv"))
cat(sprintf("Table 5: HPV rates across %d tissue categories\n", nrow(table5)))

# ── Table 6: HPV rates stratified by cell-line vs primary tissue ───────
# Cell lines often dominate the HPV+ count; this split shows whether the
# rate in real tissue is meaningfully different from the cell-line subset
# (positive controls).
table6 <- hpv_status_enriched %>%
  mutate(sample_class = ifelse(is_cell_line, "cell_line", "tissue")) %>%
  group_by(tissue_category, sample_class) %>%
  summarise(
    total_samples = n(),
    hpv_positive = sum(hpv_status == "HPV+"),
    hpv_negative = sum(hpv_status == "HPV-"),
    pct_positive = round(100 * hpv_positive / total_samples, 1),
    .groups = "drop"
  ) %>%
  arrange(tissue_category, sample_class)

write_tsv(table6, file.path(outdir, "table6_hpv_rates_by_cellline.tsv"))
cat(sprintf("Table 6: HPV rates split by cell-line vs tissue (%d rows)\n",
            nrow(table6)))

# ── Figures ─────────────────────────────────────────────────────────────

# Heatmap: HPV type x tissue category
if (nrow(hpv_full) > 0) {
  heatmap_data <- hpv_full %>%
    group_by(tissue_category, hpv_reference) %>%
    summarise(n = n_distinct(sample_id), .groups = "drop")

  p_heatmap <- ggplot(heatmap_data, aes(x = tissue_category, y = hpv_reference, fill = n)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = "HPV Type Prevalence by Tissue Category",
         x = "Tissue Category", y = "HPV Type", fill = "# Samples") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))

  ggsave(file.path(outdir, "heatmap_hpv_tissue.pdf"), p_heatmap,
         width = 10, height = max(8, nrow(heatmap_data) * 0.15))
  ggsave(file.path(outdir, "heatmap_hpv_tissue.png"), p_heatmap,
         width = 10, height = max(8, nrow(heatmap_data) * 0.15), dpi = 150)
}

# Bar plot: HPV+/- per tissue category
if (nrow(table5) > 0) {
  p_bar <- table5 %>%
    pivot_longer(cols = c(hpv_positive, hpv_negative),
                 names_to = "status", values_to = "count") %>%
    mutate(status = ifelse(status == "hpv_positive", "HPV+", "HPV-")) %>%
    ggplot(aes(x = tissue_category, y = count, fill = status)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("HPV+" = "#e63946", "HPV-" = "#457b9d")) +
    labs(title = "HPV Status by Tissue Category",
         x = "Tissue Category", y = "Number of Samples", fill = "HPV Status") +
    theme_minimal()

  ggsave(file.path(outdir, "barplot_hpv_status.pdf"), p_bar, width = 8, height = 5)
  ggsave(file.path(outdir, "barplot_hpv_status.png"), p_bar, width = 8, height = 5, dpi = 150)
}

# Stratified bar plot: HPV+/- per tissue category, faceted by cell-line vs tissue
if (nrow(table6) > 0) {
  p_bar_strat <- table6 %>%
    pivot_longer(cols = c(hpv_positive, hpv_negative),
                 names_to = "status", values_to = "count") %>%
    mutate(status = ifelse(status == "hpv_positive", "HPV+", "HPV-")) %>%
    ggplot(aes(x = tissue_category, y = count, fill = status)) +
    geom_col(position = "stack") +
    facet_wrap(~ sample_class, scales = "free_y") +
    scale_fill_manual(values = c("HPV+" = "#e63946", "HPV-" = "#457b9d")) +
    labs(title = "HPV Status by Tissue Category — Cell Lines vs. Primary Tissue",
         x = "Tissue Category", y = "Number of Samples", fill = "HPV Status") +
    theme_minimal()

  ggsave(file.path(outdir, "barplot_hpv_status_by_cellline.pdf"),
         p_bar_strat, width = 10, height = 5)
  ggsave(file.path(outdir, "barplot_hpv_status_by_cellline.png"),
         p_bar_strat, width = 10, height = 5, dpi = 150)
}

# ── HTML Report ─────────────────────────────────────────────────────────
# Create a minimal Rmd for rendering
rmd_content <- sprintf('---
title: "HPV in Human Skin — Analysis Report"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(tidyverse)
library(knitr)
```

## Overview

This report summarizes HPV type prevalence across human tissues using public RNA-seq data.

- **Total samples analyzed:** `r nrow(read_tsv("%s", show_col_types = FALSE))`
- **HPV+ samples:** `r sum(read_tsv("%s", show_col_types = FALSE)$hpv_status == "HPV+")`

## 1. HPV Types in Healthy Skin

```{r}
kable(read_tsv("%s/table1_hpv_healthy_skin.tsv", show_col_types = FALSE))
```

## 2. HPV Types by Skin Pathology

```{r}
kable(read_tsv("%s/table2_hpv_skin_pathology.tsv", show_col_types = FALSE))
```

## 3. HPV in Non-Traditional Tissues

```{r}
kable(read_tsv("%s/table3_hpv_nontraditional_tissues.tsv", show_col_types = FALSE))
```

## 4. Productive Infection (Late Transcripts)

```{r}
kable(read_tsv("%s/table4_productive_infection.tsv", show_col_types = FALSE))
```

## 5. HPV-Negative Rates

```{r}
kable(read_tsv("%s/table5_hpv_negative_rates.tsv", show_col_types = FALSE))
```

## 6. HPV Rates: Cell Lines vs. Primary Tissue

Cell lines (HeLa, SiHa, CaSki, TIGK, HaCaT, etc.) are kept in the analysis as positive controls. This table separates the HPV+/- rate within cell-line samples from primary tissue so the latter can be interpreted on its own.

```{r}
kable(read_tsv("%s/table6_hpv_rates_by_cellline.tsv", show_col_types = FALSE))
```

## Figures

### HPV Type x Tissue Heatmap
![](summary_tables/heatmap_hpv_tissue.png)

### HPV Status by Tissue
![](summary_tables/barplot_hpv_status.png)

### HPV Status — Cell Lines vs. Primary Tissue
![](summary_tables/barplot_hpv_status_by_cellline.png)
',
  opts$`hpv-status`, opts$`hpv-status`,
  outdir, outdir, outdir, outdir, outdir, outdir
)

rmd_file <- "hpv_skin_report.Rmd"
writeLines(rmd_content, rmd_file)
rmarkdown::render(rmd_file, output_file = "hpv_skin_report.html", quiet = TRUE)

cat("\nReport generated: hpv_skin_report.html\n")
cat(sprintf("Summary tables: %s/\n", outdir))
