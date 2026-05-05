#!/usr/bin/env Rscript
# Build an Excel workbook for manual curation of typed/HPV+ samples.
#
# The workbook is structured for round-tripping: identifier and auto-classified
# columns are write-protected, decision columns are unlocked and have dropdown
# data validation. Once curated, feed it back through bin/apply_curation.R
# (TODO) to override is_cell_line / is_engineered for downstream analysis.
#
# Output: reports/<date>_curation.xlsx

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(openxlsx)
})

date_tag <- format(Sys.Date(), "%Y-%m-%d")
out_path <- here("reports", paste0(date_tag, "_curation.xlsx"))

# --- Load data ---------------------------------------------------------------
ss <- read_csv(here("results_full_v2/metadata/samplesheet_enriched.csv"),
               show_col_types = FALSE)
flags <- read_tsv(here("results_full_v2/metadata/cell_line_flags.tsv"),
                  show_col_types = FALSE,
                  col_types = cols(srr_id = "c", is_cell_line = "l",
                                   is_engineered = "l", .default = "c"))
status <- read_tsv(here("results_full_v2/metadata/hpv_status.tsv"),
                   show_col_types = FALSE)
best <- read_tsv(here("results_full_v2/sweep/best_hits_permissive.tsv"),
                 show_col_types = FALSE) |>
  rename(srr_id = sample_id, best_ref = hpv_reference,
         best_breadth = coverage_breadth, best_depth = mean_depth) |>
  select(-ref_length, -read_count, -covered_bases)

# Per-sample annotated table over all HPV+ samples (n = 90).
hpv_pos <- status |>
  filter(hpv_status == "HPV+") |>
  rename(hpv_read_count_kraken2 = hpv_read_count) |>
  left_join(ss |> select(srr_id, srx_id, study, title, tissue_source),
            by = "srr_id") |>
  left_join(flags |> select(srr_id, is_cell_line, is_engineered,
                            cell_line_pattern, engineered_pattern),
            by = "srr_id") |>
  left_join(best, by = "srr_id") |>
  mutate(
    is_cell_line  = replace_na(is_cell_line,  FALSE),
    is_engineered = replace_na(is_engineered, FALSE),
    auto_category = case_when(
      is_cell_line  ~ "cell_line",
      is_engineered ~ "engineered",
      TRUE          ~ "clinical"
    ),
    typed = !is.na(best_ref)
  )

# Per-study summary over the typed subset (n = 31 samples in 16 studies).
typed_only <- hpv_pos |> filter(typed)
per_study <- typed_only |>
  group_by(study) |>
  summarise(
    n_typed_samples = n(),
    n_total_samples_in_cohort = sum(ss$study == first(study)),
    refs_hit = paste(sort(unique(best_ref)), collapse = ", "),
    auto_categories = paste(sort(unique(auto_category)), collapse = ", "),
    example_titles = paste(unique(na.omit(title))[
      seq_len(min(3, length(unique(na.omit(title)))))
    ], collapse = " | "),
    .groups = "drop"
  ) |>
  mutate(
    is_alpha = str_detect(refs_hit, "HPV1[68]REF"),
    has_clinical = str_detect(auto_categories, "clinical"),
    priority = case_when(
      is_alpha & has_clinical               ~ "P1 (alpha-clinical, suspect)",
      !is_alpha & has_clinical              ~ "P2 (clinical β/γ, sanity check)",
      auto_categories == "cell_line"        ~ "P4a (cell_line, formality)",
      auto_categories == "engineered"       ~ "P4b (engineered, formality)",
      TRUE                                  ~ "P3 (mixed, review)"
    ),
    sra_url = case_when(
      str_starts(study, "ERP") ~ paste0("https://www.ebi.ac.uk/ena/browser/view/", study),
      TRUE                     ~ paste0("https://www.ncbi.nlm.nih.gov/Traces/study/?acc=", study)
    )
  ) |>
  select(priority, study, sra_url, n_typed_samples, n_total_samples_in_cohort,
         refs_hit, auto_categories, example_titles) |>
  arrange(priority, desc(n_typed_samples)) |>
  mutate(
    confirmed_provenance = "",
    study_notes          = "",
    curated_by           = "",
    curated_date         = ""
  )

# Per-sample curation table — order by priority then study.
per_sample <- hpv_pos |>
  arrange(desc(typed), auto_category, study, srr_id) |>
  select(srr_id, srx_id, study, title, tissue_source,
         tissue_category, diagnosis,
         hpv_read_count_kraken2, typed, best_ref, best_breadth, best_depth,
         is_cell_line, is_engineered, auto_category,
         cell_line_pattern, engineered_pattern) |>
  mutate(
    confirmed_category    = "",
    confirmed_subcategory = "",
    notes                 = "",
    curated_by            = "",
    curated_date          = ""
  )

# README sheet content
readme <- tibble(field = c(
  "Purpose",
  "Source snapshot",
  "",
  "Sheet: per_study_curation",
  "Sheet: per_sample_curation",
  "",
  "Editable columns (per_study)",
  "Editable columns (per_sample)",
  "",
  "Dropdown values: confirmed_provenance / confirmed_category",
  "",
  "Round-trip workflow",
  "",
  "Coverage thresholds in 'best_*' columns",
  "Auto-classification source"
), value = c(
  "Manual curation of typed and HPV+ samples for the full_v2 run.",
  paste0("Generated ", Sys.time(), " from results_full_v2/ snapshot."),
  "",
  "One row per study contributing typed HPV+ samples (16 studies). Use this when whole-study reclassification is appropriate.",
  "One row per HPV+ sample (90 rows: 31 typed + 59 untyped). Use this when per-sample decisions are needed.",
  "",
  "confirmed_provenance, study_notes, curated_by, curated_date",
  "confirmed_category, confirmed_subcategory, notes, curated_by, curated_date",
  "",
  "cell_line | engineered | clinical | mixed | unknown | verify_individually",
  "",
  "Edit decision columns and save. A future bin/apply_curation.R will read this back and override is_cell_line / is_engineered for re-analysis.",
  "",
  "Permissive cutoff: breadth >= 0.01, mean_depth >= 0.25 (all 31 typed samples). Default pipeline thresholds: 0.10 / 2.0 (only 14 samples).",
  "Heuristic patterns in bin/cell_line_patterns.py. cell_line_pattern / engineered_pattern columns show which regex matched (NA = not matched)."
))

# --- Build workbook ----------------------------------------------------------
wb <- createWorkbook()

# Styles
hdr_style <- createStyle(textDecoration = "bold", fgFill = "#1F4E79",
                        fontColour = "#FFFFFF", border = "Bottom",
                        halign = "center")
locked_style   <- createStyle(fgFill = "#F2F2F2", locked = TRUE)
editable_style <- createStyle(fgFill = "#FFF2CC", locked = FALSE,
                              border = "TopBottomLeftRight",
                              borderColour = "#BFBFBF")
priority_p1 <- createStyle(fgFill = "#F4B084", textDecoration = "bold")
priority_p2 <- createStyle(fgFill = "#FFE699")
priority_p3 <- createStyle(fgFill = "#C6E0B4")
priority_p4 <- createStyle(fgFill = "#D9D9D9")
url_style   <- createStyle(fontColour = "#0563C1", textDecoration = "underline")

# README
addWorksheet(wb, "README")
writeData(wb, "README", readme)
addStyle(wb, "README", hdr_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "README", cols = 1, widths = 50)
setColWidths(wb, "README", cols = 2, widths = 90)

# Per-study sheet
addWorksheet(wb, "per_study_curation", zoom = 110)
writeData(wb, "per_study_curation", per_study)

# Make sra_url clickable
class(per_study$sra_url) <- "hyperlink"
writeData(wb, "per_study_curation", per_study, withFilter = TRUE)

addStyle(wb, "per_study_curation", hdr_style,
         rows = 1, cols = seq_len(ncol(per_study)), gridExpand = TRUE)
freezePane(wb, "per_study_curation", firstActiveRow = 2, firstActiveCol = 3)

editable_cols_study <- which(names(per_study) %in%
  c("confirmed_provenance", "study_notes", "curated_by", "curated_date"))
locked_cols_study <- setdiff(seq_len(ncol(per_study)), editable_cols_study)
n_study_rows <- nrow(per_study)

addStyle(wb, "per_study_curation", locked_style,
         rows = 2:(n_study_rows + 1), cols = locked_cols_study,
         gridExpand = TRUE, stack = TRUE)
addStyle(wb, "per_study_curation", editable_style,
         rows = 2:(n_study_rows + 1), cols = editable_cols_study,
         gridExpand = TRUE, stack = TRUE)

# Priority colour bands (column 1)
for (i in seq_len(n_study_rows)) {
  pri <- per_study$priority[i]
  style <- switch(substr(pri, 1, 2),
                  P1 = priority_p1, P2 = priority_p2,
                  P3 = priority_p3, P4 = priority_p4, NULL)
  if (!is.null(style)) {
    addStyle(wb, "per_study_curation", style,
             rows = i + 1, cols = 1, stack = TRUE)
  }
}

# Dropdown for confirmed_provenance
prov_choices <- c("cell_line", "engineered", "clinical",
                  "mixed", "unknown", "verify_individually")
dataValidation(wb, "per_study_curation",
               cols = which(names(per_study) == "confirmed_provenance"),
               rows = 2:(n_study_rows + 1),
               type = "list",
               value = paste0('"', paste(prov_choices, collapse = ","), '"'))

setColWidths(wb, "per_study_curation",
             cols = seq_len(ncol(per_study)),
             widths = c(34, 12, 60, 8, 10, 28, 18, 60, 18, 30, 12, 12))

protectWorksheet(wb, "per_study_curation", protect = TRUE,
                 password = NULL,
                 lockSelectingLockedCells = FALSE,
                 lockSelectingUnlockedCells = FALSE,
                 lockFormattingCells = TRUE,
                 lockSorting = FALSE, lockAutoFilter = FALSE)

# Per-sample sheet
addWorksheet(wb, "per_sample_curation", zoom = 100)
writeData(wb, "per_sample_curation", per_sample, withFilter = TRUE)

addStyle(wb, "per_sample_curation", hdr_style,
         rows = 1, cols = seq_len(ncol(per_sample)), gridExpand = TRUE)
freezePane(wb, "per_sample_curation", firstActiveRow = 2, firstActiveCol = 4)

editable_cols_sample <- which(names(per_sample) %in%
  c("confirmed_category", "confirmed_subcategory", "notes",
    "curated_by", "curated_date"))
locked_cols_sample <- setdiff(seq_len(ncol(per_sample)), editable_cols_sample)
n_sample_rows <- nrow(per_sample)

addStyle(wb, "per_sample_curation", locked_style,
         rows = 2:(n_sample_rows + 1), cols = locked_cols_sample,
         gridExpand = TRUE, stack = TRUE)
addStyle(wb, "per_sample_curation", editable_style,
         rows = 2:(n_sample_rows + 1), cols = editable_cols_sample,
         gridExpand = TRUE, stack = TRUE)

# Dropdown
dataValidation(wb, "per_sample_curation",
               cols = which(names(per_sample) == "confirmed_category"),
               rows = 2:(n_sample_rows + 1),
               type = "list",
               value = paste0('"', paste(prov_choices, collapse = ","), '"'))

setColWidths(wb, "per_sample_curation",
             cols = seq_len(ncol(per_sample)),
             widths = c(13, 13, 12, 35, 30, 14, 14,
                        12, 8, 18, 12, 12, 12, 13, 13,
                        25, 25, 25, 25, 30, 12, 12))

protectWorksheet(wb, "per_sample_curation", protect = TRUE,
                 password = NULL,
                 lockSelectingLockedCells = FALSE,
                 lockSelectingUnlockedCells = FALSE,
                 lockFormattingCells = TRUE,
                 lockSorting = FALSE, lockAutoFilter = FALSE)

saveWorkbook(wb, out_path, overwrite = TRUE)
message(sprintf("Wrote %s\n  per_study_curation: %d rows\n  per_sample_curation: %d rows",
                out_path, n_study_rows, n_sample_rows))
