# Single source of truth for the depth-normalised burden analysis directory.
#
# `results_reframe_fixed` is the output of bin/hpv_rpm_contrast.R run with the
# curated tier corrections in assets/tier_overrides.csv. The earlier
# `results_reframe_csc` is retained on disk for provenance: it is the same run
# with the corrections disabled, and differs only in the tier assigned to 31
# libraries (EV biopsies, condyloma-adjacent mucosa, an HPV31+ cell line and a
# wart mis-flagged as a cell line).
REFRAME_DIR <- "results_reframe_fixed"
rpm_per_sample_tsv <- function() here::here(REFRAME_DIR, "hpv_rpm_per_sample.tsv")
type_catalogue_tsv <- function() here::here(REFRAME_DIR, "type_catalogue_by_tier.tsv")

# The same curated is_cell_line corrections that bin/hpv_rpm_contrast.R applies
# to the tier table must also reach Model 1, which classifies libraries by that
# flag. Study-level rules are expanded to run accessions via the null-sweep
# samplesheet; lesion-cohort accessions simply do not match and are ignored.
cellline_overrides <- function() {
  ov <- readr::read_csv(here::here("assets/tier_overrides.csv"),
                        comment = "#", show_col_types = FALSE) %>%
    dplyr::filter(field == "is_cell_line")
  if (!nrow(ov)) return(dplyr::tibble(srr_id = character(), ovr = logical()))
  ss <- readr::read_csv(here::here("results_full_v2/metadata/samplesheet_enriched.csv"),
                        show_col_types = FALSE) %>% dplyr::select(srr_id, study)
  dplyr::bind_rows(
    ov %>% dplyr::filter(scope == "study") %>%
      dplyr::inner_join(ss, by = c("id" = "study")) %>% dplyr::select(srr_id, value),
    ov %>% dplyr::filter(scope == "sample") %>% dplyr::transmute(srr_id = id, value)) %>%
    dplyr::distinct(srr_id, .keep_all = TRUE) %>%
    dplyr::transmute(srr_id, ovr = toupper(value) == "TRUE")
}
apply_cellline_overrides <- function(fl) {
  fl %>% dplyr::left_join(cellline_overrides(), by = "srr_id") %>%
    dplyr::mutate(is_cell_line = dplyr::coalesce(ovr, is_cell_line)) %>% dplyr::select(-ovr)
}

# Libraries whose Kraken2 report file is empty produced no screening result at
# all. Two causes, measured 2026-11-03: fastp discarded every read (410 in the
# sweep, all "too_short_reads", because params.fastp_min_length = 36 meets an R1
# of 26 to 33 bp, the barcode mate of 10x-style chemistry, and fastp in paired
# mode drops the pair when either mate fails), or the library was never
# downloaded (851, scattered 10 to 22 per 100-library chunk across 85 chunks,
# i.e. per-accession download failures rather than failed chunks).
#
# results_full_v2/metadata/hpv_status.tsv records all of them as "HPV-" with a
# read count of 0, so leaving them in counts a library that was never screened
# as a negative observation and deflates every positivity rate. They are
# therefore dropped from the positivity denominator. The depth-normalized burden
# analysis never had the problem: bin/hpv_rpm_contrast.R reads the report files
# themselves, and an empty one yields no RPM row.
unscreened_ids <- function(dir = "results_full_v2/aggregated") {
  f <- list.files(here::here(dir), pattern = "_kraken2_report\\.txt$", full.names = TRUE)
  sub("_kraken2_report\\.txt$", "", basename(f[file.size(f) == 0]))
}
