#!/usr/bin/env Rscript
# Apply manual curation from the Excel workbook to produce a curated
# cell_line_flags table. Per-sample curation (if any) takes precedence
# over per-study curation, which takes precedence over the heuristic.
#
# Inputs:
#   reports/<date>_curation.xlsx
#   results_full_v2/metadata/cell_line_flags.tsv
#   results_full_v2/metadata/samplesheet_enriched.csv
#
# Output:
#   results_full_v2/metadata/cell_line_flags_curated.tsv
#     Same schema as cell_line_flags.tsv plus a `source` column tracking
#     which classifier produced the call (heuristic | study_curation |
#     sample_curation). Downstream (bin/sweep.R) prefers this file when
#     present.
#
# Usage:
#   Rscript bin/apply_curation.R [path/to/curation.xlsx]
#   Defaults to the most recent reports/*_curation.xlsx.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(openxlsx)
})

# Locate the curation workbook ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  curation_path <- args[1]
} else {
  candidates <- list.files(here("reports"), pattern = "_curation\\.xlsx$",
                           full.names = TRUE)
  if (!length(candidates))
    stop("No reports/*_curation.xlsx found and no path supplied.")
  curation_path <- tail(sort(candidates), 1)
}
message("Reading curation from: ", curation_path)

# Override-rule table: confirmed_* values that translate to definite booleans.
# Anything else (mixed / unknown / verify_individually / NA / "") leaves the
# heuristic call alone.
override_map <- tribble(
  ~confirmed,    ~ov_is_cell_line, ~ov_is_engineered,
  "cell_line",   TRUE,             FALSE,
  "engineered",  FALSE,            TRUE,
  "clinical",    FALSE,            FALSE,
)

normalize <- function(x) {
  x <- ifelse(is.na(x), "", str_trim(tolower(x)))
  x
}

# Load inputs -----------------------------------------------------------------
ss <- read_csv(here("results_full_v2/metadata/samplesheet_enriched.csv"),
               show_col_types = FALSE) |>
  select(srr_id, study)

heur <- read_tsv(here("results_full_v2/metadata/cell_line_flags.tsv"),
                 show_col_types = FALSE,
                 col_types = cols(srr_id = "c", is_cell_line = "l",
                                  is_engineered = "l", .default = "c"))

study_sheet  <- read.xlsx(curation_path, sheet = "per_study_curation")
sample_sheet <- read.xlsx(curation_path, sheet = "per_sample_curation")

# Build study- and sample-level override tables -------------------------------
study_ov <- study_sheet |>
  mutate(confirmed = normalize(confirmed_provenance)) |>
  inner_join(override_map, by = "confirmed") |>
  transmute(study,
            study_is_cell_line  = ov_is_cell_line,
            study_is_engineered = ov_is_engineered,
            study_note          = study_notes)

sample_ov <- sample_sheet |>
  mutate(confirmed = normalize(confirmed_category)) |>
  inner_join(override_map, by = "confirmed") |>
  transmute(srr_id,
            sample_is_cell_line  = ov_is_cell_line,
            sample_is_engineered = ov_is_engineered,
            sample_note          = notes)

message(sprintf("Study-level overrides: %d studies", nrow(study_ov)))
message(sprintf("Sample-level overrides: %d samples", nrow(sample_ov)))

# Apply --------------------------------------------------------------------
# Sample > study > heuristic. Track which source produced each row.
curated <- heur |>
  left_join(ss,        by = "srr_id") |>
  left_join(study_ov,  by = "study") |>
  left_join(sample_ov, by = "srr_id") |>
  mutate(
    new_is_cell_line  = coalesce(sample_is_cell_line,  study_is_cell_line,  is_cell_line),
    new_is_engineered = coalesce(sample_is_engineered, study_is_engineered, is_engineered),
    source = case_when(
      !is.na(sample_is_cell_line) ~ "sample_curation",
      !is.na(study_is_cell_line)  ~ "study_curation",
      TRUE                        ~ "heuristic"
    )
  )

# Stats: how many heuristic calls actually changed
changes <- curated |>
  mutate(
    cl_changed  = is_cell_line  != new_is_cell_line,
    eng_changed = is_engineered != new_is_engineered
  )

message(sprintf(
  "Heuristic vs curated: %d cell_line flips, %d engineered flips, %d unchanged",
  sum(changes$cl_changed, na.rm = TRUE),
  sum(changes$eng_changed, na.rm = TRUE),
  sum(!changes$cl_changed & !changes$eng_changed, na.rm = TRUE)
))
message("Calls by source:")
curated |> count(source) |> print()

# Build output table (preserve original schema + source column)
out <- curated |>
  transmute(
    srr_id,
    is_cell_line       = new_is_cell_line,
    cell_line_pattern,
    is_engineered      = new_is_engineered,
    engineered_pattern,
    source_text,
    source
  )

# Write
out_path <- here("results_full_v2/metadata/cell_line_flags_curated.tsv")
write_tsv(out, out_path, na = "")
message("Wrote: ", out_path)
