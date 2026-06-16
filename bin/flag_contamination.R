#!/usr/bin/env Rscript

# HPV contamination / cell-line flagging — BASE R ONLY.
#
# The cluster's r-tidyverse module stack segfaults at library load, so this
# script deliberately uses only base R (read.csv/read.delim, aggregate, tapply).
# It does NOT re-implement cell-line detection: it consumes the is_cell_line /
# is_engineered columns already produced by parse_metadata.py / cell_line_patterns.py,
# and the genus assignment from build_genus_lookup.py. Flag, never exclude —
# every call is kept and labelled so contamination/cell-line positives stay as
# controls.
#
# Two findings drive the logic:
#   1. Coverage breadth is bimodal: an artifact "floor" cluster sits at the
#      0.10-0.12 typing threshold; genuine signal sits well above. The valley
#      (default 0.15) separates them.
#   2. Floor-cluster HPV16 calls cluster by study -> project-level contamination.
#
# Usage:
#   flag_contamination.R --mode diagnostic --indir results_full_v2/aggregated \
#       --genus-lookup type_genus.csv --samplesheet samplesheet_enriched_v2.csv
#   flag_contamination.R --mode flag --breadth-threshold 0.15 \
#       --indir results_full_v2/aggregated --genus-lookup type_genus.csv \
#       --samplesheet samplesheet_enriched_v2.csv --outdir contamination

# ── tiny base-R arg parser (no optparse dependency) ─────────────────────
parse_args <- function(args) {
  defaults <- list(
    mode = "flag", indir = "results_full_v2/aggregated", `hpv-types` = NA,
    `genus-lookup` = "type_genus.csv", samplesheet = NA,
    outdir = "contamination", `breadth-threshold` = "0.15"
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
      defaults[[key]] <- args[i + 1]; i <- i + 2
    } else { defaults[[key]] <- TRUE; i <- i + 1 }
  }
  defaults
}
opts <- parse_args(commandArgs(trailingOnly = TRUE))
breadth_thr <- as.numeric(opts$`breadth-threshold`)
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

# ── load HPV type calls ─────────────────────────────────────────────────
read_types <- function() {
  if (!is.na(opts$`hpv-types`) && file.exists(opts$`hpv-types`)) {
    return(read.delim(opts$`hpv-types`, stringsAsFactors = FALSE))
  }
  files <- list.files(opts$indir, pattern = "_hpv_types\\.tsv$", full.names = TRUE)
  if (length(files) == 0) stop("No *_hpv_types.tsv found in --indir and no --hpv-types given")
  parts <- lapply(files, function(f) {
    df <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) NULL else df
  })
  do.call(rbind, parts)
}
calls <- read_types()
cat(sprintf("Loaded %d HPV type calls from %d samples\n",
            nrow(calls), length(unique(calls$sample_id))))

# ── join genus ──────────────────────────────────────────────────────────
genus <- read.csv(opts$`genus-lookup`, stringsAsFactors = FALSE)
calls <- merge(calls, genus[, c("hpv_reference", "type", "genus", "genus_short", "species")],
               by = "hpv_reference", all.x = TRUE)
calls$genus_short[is.na(calls$genus_short)] <- "unknown"
# numeric HPV type where parseable (for alpha-9 / HPV16 / HPV18 selectors)
calls$type_num <- suppressWarnings(as.integer(calls$type))
calls$is_hpv16 <- calls$type_num == 16 & !is.na(calls$type_num)
calls$is_hpv18 <- calls$type_num == 18 & !is.na(calls$type_num)
calls$is_alpha <- calls$genus_short == "alpha"

# ── join sample metadata (tissue, cell-line, engineered, study) ─────────
if (!is.na(opts$samplesheet) && file.exists(opts$samplesheet)) {
  ss <- read.csv(opts$samplesheet, stringsAsFactors = FALSE)
} else {
  stop("--samplesheet is required")
}
norm_bool <- function(x) if (is.null(x)) FALSE else tolower(as.character(x)) %in% c("true", "1", "yes")
# Tolerate older samplesheets lacking the flag columns (treat as not-a-line).
for (col in c("is_cell_line", "is_engineered")) if (!col %in% names(ss)) ss[[col]] <- "false"
for (col in c("study", "tissue_category", "diagnosis")) if (!col %in% names(ss)) ss[[col]] <- NA
ss_keep <- ss[, c("srr_id", "study", "tissue_category", "diagnosis",
                  "is_cell_line", "is_engineered")]
calls <- merge(calls, ss_keep, by.x = "sample_id", by.y = "srr_id", all.x = TRUE)
calls$is_cell_line  <- vapply(calls$is_cell_line, norm_bool, logical(1))
calls$is_engineered <- vapply(calls$is_engineered, norm_bool, logical(1))

# ── ASCII histogram helper ──────────────────────────────────────────────
ascii_hist <- function(x, label, width = 50) {
  x <- x[!is.na(x)]
  cat(sprintf("\n  %s (n=%d)\n", label, length(x)))
  if (length(x) == 0) { cat("    (none)\n"); return(invisible()) }
  brks <- seq(0, 1, by = 0.05)
  h <- table(cut(x, breaks = brks, include.lowest = TRUE))
  mx <- max(h)
  for (k in seq_along(h)) {
    bar <- strrep("#", round(width * h[k] / mx))
    cat(sprintf("    %.2f-%.2f | %-*s %d\n", brks[k], brks[k + 1], width, bar, h[k]))
  }
}

# ════════════════════════════════════════════════════════════════════════
if (identical(opts$mode, "diagnostic")) {
  cat("\n=== DIAGNOSTIC: coverage breadth distributions ===\n")
  ascii_hist(calls$coverage_breadth, "ALL calls")
  ascii_hist(calls$coverage_breadth[calls$is_hpv16], "HPV16")
  ascii_hist(calls$coverage_breadth[calls$is_hpv18], "HPV18")
  ascii_hist(calls$coverage_breadth[calls$is_alpha], "alpha genus")
  q <- quantile(calls$coverage_breadth, probs = seq(0, 1, 0.1), na.rm = TRUE)
  qdf <- data.frame(quantile = names(q), coverage_breadth = as.numeric(q))
  write.table(qdf, file.path(opts$outdir, "breadth_quantiles.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\nWrote %s/breadth_quantiles.tsv\n", opts$outdir))
  cat(sprintf("Suggested threshold (current): %.2f\n", breadth_thr))
  quit(save = "no")
}

# ════════════════════════════════════════════════════════════════════════
# flag mode
cat(sprintf("\n=== FLAG mode (breadth threshold = %.2f) ===\n", breadth_thr))

calls$flag_low_breadth   <- calls$coverage_breadth < breadth_thr
calls$flag_hpv16_suspect <- calls$is_hpv16 & calls$flag_low_breadth
calls$flag_hpv18_suspect <- calls$is_hpv18 & calls$flag_low_breadth
calls$flag_alpha_suspect <- calls$is_alpha & calls$flag_low_breadth
# suspect_any: floor-cluster artifact OR known cell-line / engineered culture.
calls$suspect_any <- calls$flag_low_breadth | calls$is_cell_line | calls$is_engineered
# "clean" = survives breadth AND not a cultured/engineered sample.
calls$is_clean <- !calls$suspect_any

write.table(calls, file.path(opts$outdir, "hpv_calls_flagged.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# flag_summary
summ <- data.frame(
  metric = c("total_calls", "low_breadth", "cell_line", "engineered",
             "suspect_any", "clean", "hpv16_suspect", "hpv18_suspect", "alpha_suspect"),
  n = c(nrow(calls), sum(calls$flag_low_breadth), sum(calls$is_cell_line),
        sum(calls$is_engineered), sum(calls$suspect_any), sum(calls$is_clean),
        sum(calls$flag_hpv16_suspect), sum(calls$flag_hpv18_suspect),
        sum(calls$flag_alpha_suspect))
)
write.table(summ, file.path(opts$outdir, "flag_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# study clusters: where do low-breadth HPV16 calls concentrate? (project-level contamination)
sc16 <- calls[calls$flag_hpv16_suspect, ]
if (nrow(sc16) > 0) {
  study_clusters <- as.data.frame(table(study = sc16$study), stringsAsFactors = FALSE)
  names(study_clusters)[2] <- "n_hpv16_low"
  study_clusters <- study_clusters[order(-study_clusters$n_hpv16_low), ]
} else {
  study_clusters <- data.frame(study = character(), n_hpv16_low = integer())
}
write.table(study_clusters, file.path(opts$outdir, "study_clusters.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# genus breakdown among CLEAN calls
clean <- calls[calls$is_clean, ]
gb <- as.data.frame(table(genus = clean$genus_short), stringsAsFactors = FALSE)
names(gb)[2] <- "n_clean_calls"
gb <- gb[order(-gb$n_clean_calls), ]
write.table(gb, file.path(opts$outdir, "genus_breakdown_clean.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# genus x tissue (all calls)
gt <- as.data.frame(table(tissue = calls$tissue_category, genus = calls$genus_short),
                    stringsAsFactors = FALSE)
gt <- gt[gt$Freq > 0, ]
names(gt)[3] <- "n_calls"
gt <- gt[order(gt$tissue, -gt$n_calls), ]
write.table(gt, file.path(opts$outdir, "genus_breakdown_by_tissue.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── console summary ──────────────────────────────────────────────────────
cat("\n-- flag summary --\n")
print(summ, row.names = FALSE)
if (nrow(study_clusters) > 0) {
  cat("\n-- studies with low-breadth HPV16 (project-level contamination) --\n")
  print(head(study_clusters, 15), row.names = FALSE)
}

# SKIN-ONLY (nahk) headline — the biological question.
skin <- calls[!is.na(calls$tissue_category) & calls$tissue_category == "nahk", ]
skin_clean <- skin[skin$is_clean, ]
cat("\n-- SKIN (nahk) --\n")
cat(sprintf("  HPV calls in skin:            %d (across %d samples)\n",
            nrow(skin), length(unique(skin$sample_id))))
cat(sprintf("  cell-line / engineered:       %d\n",
            sum(skin$is_cell_line | skin$is_engineered)))
cat(sprintf("  low-breadth (floor artifact): %d\n", sum(skin$flag_low_breadth)))
cat(sprintf("  CLEAN skin calls remaining:   %d (across %d samples)\n",
            nrow(skin_clean), length(unique(skin_clean$sample_id))))
if (nrow(skin_clean) > 0) {
  st <- as.data.frame(table(genus = skin_clean$genus_short), stringsAsFactors = FALSE)
  st <- st[st$Freq > 0, ]; names(st)[2] <- "n"
  cat("  clean skin genus breakdown:\n")
  for (r in seq_len(nrow(st))) cat(sprintf("    %-8s %d\n", st$genus[r], st$n[r]))
}
cat(sprintf("\nOutputs written to %s/\n", opts$outdir))
