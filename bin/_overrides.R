# Shared curated-override layer (assets/tier_overrides.csv).
#
# The overrides are a curation layer applied at ANALYSIS time; they are
# deliberately not written back into the parsed samplesheets, so that what SRA
# said and what we judged stay separable (and so a re-run of parse_metadata.py
# cannot silently revert a curation). The consequence is that every consumer of
# tissue_category / is_cell_line has to apply them itself, which is what this
# file is for. bin/hpv_rpm_contrast.R applies them inline as part of its tier
# assignment; the standalone reporting scripts source this.
#
# Row format: scope(study|sample), id, field, value, reason
#   field = tissue_category | is_cell_line | tier | exclude
# Sample-scoped rows win over study-scoped rows for the same library.
#
# Usage:
#   source(file.path(dirname(this_script), "_overrides.R"))
#   d <- apply_overrides(d, c("tissue_category", "is_cell_line"))
#   d <- drop_excluded(d)

suppressPackageStartupMessages({library(dplyr); library(readr)})

OVERRIDES_PATH <- function() {
  p <- Sys.getenv("TIER_OVERRIDES", "")
  if (nzchar(p)) return(p)
  # repo root = parent of bin/
  here <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile %||% "."))),
                   error = function(e) ".")
  for (cand in c("assets/tier_overrides.csv",
                 file.path(here, "assets/tier_overrides.csv")))
    if (file.exists(cand)) return(cand)
  ""
}
`%||%` <- function(a, b) if (is.null(a)) b else a

read_overrides <- function(path = OVERRIDES_PATH()) {
  if (!nzchar(path) || !file.exists(path)) {
    message("  [overrides] none found; using automatic assignment only")
    return(tibble(scope = character(), id = character(), field = character(),
                  value = character(), reason = character()))
  }
  read_csv(path, comment = "#", show_col_types = FALSE)
}

# Apply one or more override fields to a data frame carrying srr_id (+ study
# for study-scoped rows). Columns absent from `d` are created.
apply_overrides <- function(d, fields, path = OVERRIDES_PATH(), quiet = FALSE) {
  ovr <- read_overrides(path)
  if (!nrow(ovr)) return(d)
  if (!"study" %in% names(d)) d$study <- NA_character_
  for (fld in fields) {
    o <- ovr %>% filter(field == fld)
    if (!nrow(o)) next
    by_sample <- o %>% filter(scope == "sample") %>% select(srr_id = id, .ovr = value)
    by_study  <- o %>% filter(scope == "study")  %>% select(study  = id, .ovr = value)
    d <- d %>%
      left_join(by_sample, by = "srr_id") %>%
      left_join(by_study,  by = "study", suffix = c("", ".study")) %>%
      mutate(.ovr = coalesce(.ovr, .ovr.study)) %>%
      select(-.ovr.study)
    if (!fld %in% names(d)) d[[fld]] <- NA
    n <- sum(!is.na(d$.ovr))
    # is_cell_line is logical in most callers; everything else is character.
    d[[fld]] <- if (is.logical(d[[fld]]))
      if_else(is.na(d$.ovr), d[[fld]], tolower(d$.ovr) %in% c("true", "1", "yes"))
    else
      coalesce(d$.ovr, as.character(d[[fld]]))
    d <- select(d, -.ovr)
    if (!quiet && n) message(sprintf("  [overrides] %s: %d libraries curated", fld, n))
  }
  d
}

# Drop libraries marked exclude=TRUE (assays that cannot detect HPV at all, e.g.
# 16S amplicon), which are structural zeros rather than negative observations.
drop_excluded <- function(d, path = OVERRIDES_PATH(), quiet = FALSE) {
  d <- apply_overrides(d, "exclude", path = path, quiet = TRUE)
  if (!"exclude" %in% names(d)) return(d)
  keep <- is.na(d$exclude) | !(tolower(as.character(d$exclude)) %in% c("true", "1", "yes"))
  if (!quiet && any(!keep))
    message(sprintf("  [overrides] excluded %d libraries", sum(!keep)))
  d[keep, ] %>% select(-exclude)
}
