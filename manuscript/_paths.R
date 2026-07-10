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
