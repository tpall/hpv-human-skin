#!/usr/bin/env Rscript
# Quick-glimpse analysis of HPV typing sweep output.
#
# Reads the three TSVs written by bin/sweep_hpv_coverage.py:
#   results_full_v2/sweep/sweep_summary.tsv
#   results_full_v2/sweep/best_hits_permissive.tsv
#   results_full_v2/sweep/combined_coverage.tsv
#
# Optionally joins a slim samplesheet (results_full_v2/samplesheet.csv) to
# stratify by is_cell_line / tissue_category when present.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

sweep_dir <- here("results_full_v2", "sweep")
stopifnot(dir.exists(sweep_dir))

read_tsv_q <- function(p) read_tsv(p, show_col_types = FALSE)

sweep_summary <- read_tsv_q(file.path(sweep_dir, "sweep_summary.tsv"))
best_hits     <- read_tsv_q(file.path(sweep_dir, "best_hits_permissive.tsv"))
combined      <- read_tsv_q(file.path(sweep_dir, "combined_coverage.tsv"))

n_samples_total <- n_distinct(combined$sample_id)
message(sprintf("Coverage table: %d samples x %d refs (%d rows)",
                n_samples_total,
                n_distinct(combined$hpv_reference),
                nrow(combined)))

# Samplesheet join: enriched samplesheet from the pipeline + post-hoc
# cell-line flags. Older runs of parse_metadata.py didn't write is_cell_line
# inline, so we always pull it from cell_line_flags.tsv when available.
ss_path    <- here("results_full_v2", "metadata", "samplesheet_enriched.csv")
flags_path <- here("results_full_v2", "metadata", "cell_line_flags.tsv")

samples <- if (file.exists(ss_path)) {
  ss <- read_csv(ss_path, show_col_types = FALSE) |>
    select(any_of(c("srr_id", "study", "tissue_category",
                    "diagnosis", "is_cell_line")))
  if (file.exists(flags_path)) {
    flags <- read_tsv(
      flags_path, show_col_types = FALSE,
      col_types = cols(
        srr_id = "c", is_cell_line = "l", is_engineered = "l", .default = "c"
      )
    ) |> select(srr_id, is_cell_line, any_of("is_engineered"))
    ss <- ss |>
      select(-any_of(c("is_cell_line", "is_engineered"))) |>
      left_join(flags, by = "srr_id")
  }
  ss |> rename(sample_id = srr_id)
} else {
  message("No samplesheet at ", ss_path,
          " - skipping cell-line stratification.")
  NULL
}

# 1. Sensitivity sweep --------------------------------------------------------
sweep_summary |>
  arrange(desc(min_breadth), desc(min_depth)) |>
  print(n = Inf)

p_sweep <- sweep_summary |>
  mutate(
    min_breadth = factor(min_breadth, levels = sort(unique(min_breadth))),
    min_depth   = factor(min_depth,   levels = sort(unique(min_depth)))
  ) |>
  ggplot(aes(min_depth, min_breadth, fill = n_samples_typed)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n_samples_typed), size = 3) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  labs(
    title    = "HPV typing yield across breadth/depth thresholds",
    subtitle = sprintf("%d samples in coverage table", n_samples_total),
    x = "min mean depth", y = "min coverage breadth", fill = "samples typed"
  ) +
  theme_minimal()

# 2. Best per-sample hit at most permissive cutoff ----------------------------
best_annot <- if (!is.null(samples)) {
  best_hits |> left_join(samples, by = "sample_id")
} else {
  best_hits
}

# Per-reference counts
ref_counts <- best_annot |>
  count(hpv_reference, sort = TRUE)
message("Top references in best_hits:")
print(ref_counts, n = 20)

p_refs <- ref_counts |>
  mutate(hpv_reference = fct_reorder(hpv_reference, n)) |>
  ggplot(aes(n, hpv_reference)) +
  geom_col(fill = "steelblue") +
  labs(title = "Best-hit reference per sample (permissive cutoff)",
       x = "samples", y = NULL) +
  theme_minimal()

# Coverage vs depth, with default-threshold cross-hairs.
# Defaults from CLAUDE.md: hpv_min_coverage = 0.10, hpv_min_depth = 2.
if (!is.null(samples) && all(c("is_cell_line", "is_engineered") %in% names(best_annot))) {
  best_annot <- best_annot |>
    mutate(
      is_cell_line  = replace_na(is_cell_line,  FALSE),
      is_engineered = replace_na(is_engineered, FALSE),
      category = case_when(
        is_cell_line  ~ "cell_line",
        is_engineered ~ "engineered",
        TRUE          ~ "clinical"
      )
    )
  colour_by <- "category"
} else {
  colour_by <- "hpv_reference"
}

p_scatter <- best_annot |>
  ggplot(aes(coverage_breadth, mean_depth, colour = .data[[colour_by]])) +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, alpha = 0.85) +
  scale_y_log10() +
  labs(
    title    = "Best per-sample HPV hit",
    subtitle = "dashed lines = default thresholds (breadth 0.10, depth 2)",
    x = "coverage breadth", y = "mean depth (log10)", colour = colour_by
  ) +
  theme_minimal()

# 3. Headline numbers ---------------------------------------------------------
default_pass <- best_hits |>
  filter(coverage_breadth >= 0.10, mean_depth >= 2)

# Screening-stage HPV+ rate (Kraken2), stratified by cell-line if available.
status_path <- here("results_full_v2", "metadata", "hpv_status.tsv")
if (file.exists(status_path)) {
  status <- read_tsv(status_path, show_col_types = FALSE)
  if (!is.null(samples)) {
    status <- status |> left_join(samples, by = c("srr_id" = "sample_id"))
  }
  cat("\n--- screening (Kraken2) ---\n")
  cat(sprintf("Samples screened: %d\n", nrow(status)))
  status |> count(hpv_status) |> print()
  if (all(c("is_cell_line", "is_engineered") %in% names(status)) && !is.null(samples)) {
    status <- status |>
      mutate(
        is_cell_line  = replace_na(is_cell_line,  FALSE),
        is_engineered = replace_na(is_engineered, FALSE),
        category = case_when(
          is_cell_line  ~ "cell_line",
          is_engineered ~ "engineered",
          TRUE          ~ "clinical"
        )
      )
    cat_table <- status |>
      count(category, hpv_status) |>
      pivot_wider(names_from = hpv_status, values_from = n, values_fill = 0) |>
      mutate(
        n_screened = `HPV+` + `HPV-`,
        pct_pos    = round(100 * `HPV+` / n_screened, 2)
      )
    print(cat_table)

    # Post-stratified overall: weight per-category rates by full-cohort
    # proportions (not by screened-cohort proportions, which are biased by
    # uneven processing order). Per-category rates themselves are unbiased.
    cohort_totals <- samples |>
      mutate(
        is_cell_line  = replace_na(is_cell_line,  FALSE),
        is_engineered = replace_na(is_engineered, FALSE),
        category = case_when(
          is_cell_line  ~ "cell_line",
          is_engineered ~ "engineered",
          TRUE          ~ "clinical"
        )
      ) |>
      count(category, name = "n_cohort")

    weighted <- cat_table |>
      left_join(cohort_totals, by = "category") |>
      mutate(
        weight       = n_cohort / sum(n_cohort),
        contribution = (`HPV+` / n_screened) * weight
      )

    cat("\n--- post-stratified prevalence (weighted to full 12.2K cohort) ---\n")
    print(weighted |> select(category, n_cohort, weight, pct_pos, contribution))
    naive   <- sum(cat_table$`HPV+`) / sum(cat_table$n_screened)
    poststr <- sum(weighted$contribution)
    n_total <- sum(weighted$n_cohort)
    cat(sprintf("Naive overall:           %.3f%%  (%d/%d screened)\n",
                100 * naive, sum(cat_table$`HPV+`), sum(cat_table$n_screened)))
    cat(sprintf("Post-stratified overall: %.3f%%  (weighted to %d cohort)\n",
                100 * poststr, n_total))
    cat(sprintf("Projected total HPV+ at completion: %.0f\n",
                poststr * n_total))
  }
}

cat("\n--- typing (sweep) ---\n")
cat(sprintf("Samples in sweep:                   %d\n", n_samples_total))
cat(sprintf("Samples with any best-hit (perm.):  %d\n", nrow(best_hits)))
cat(sprintf("Samples passing default thresholds: %d\n", nrow(default_pass)))
if (!is.null(samples) && all(c("is_cell_line", "is_engineered") %in% names(best_annot))) {
  best_annot <- best_annot |>
    mutate(
      is_cell_line  = replace_na(is_cell_line,  FALSE),
      is_engineered = replace_na(is_engineered, FALSE),
      category = case_when(
        is_cell_line  ~ "cell_line",
        is_engineered ~ "engineered",
        TRUE          ~ "clinical"
      )
    )
  cat("Best-hit by category:\n")
  print(count(best_annot, category))
}

# 4. Save plots ---------------------------------------------------------------
ggsave(file.path(sweep_dir, "sweep_heatmap.png"),     p_sweep,
       width = 6, height = 4, dpi = 150)
ggsave(file.path(sweep_dir, "best_hits_refs.png"),    p_refs,
       width = 6, height = 4, dpi = 150)
ggsave(file.path(sweep_dir, "best_hits_scatter.png"), p_scatter,
       width = 6, height = 4, dpi = 150)

message("Wrote sweep_heatmap.png, best_hits_refs.png, best_hits_scatter.png to ",
        sweep_dir)
