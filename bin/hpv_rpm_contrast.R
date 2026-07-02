#!/usr/bin/env Rscript
# hpv_rpm_contrast.R — depth-normalized HPV burden contrast + HPV-type catalogue.
#
# Reframe analysis ("where/when is cutaneous HPV transcriptionally active?"):
# treats the data as a case/control design with three tiers —
#   controls  : warts / condyloma (HPV-defined; validate the assay)
#   null      : unselected skin, with cell lines + contamination studies removed
#   open Q    : EV and AK->cSCC neoplasia (beta-HPV cofactor)
#
# The depth confound (lesions may just be sequenced deeper) is neutralised by
# scoring every screened sample as Papillomaviridae reads-per-million:
#   RPM = Papillomaviridae clade reads / (unclassified + root reads) * 1e6
# computed from Kraken2 reports so HPV-negative samples contribute too. The
# lesion-vs-null gap only matters if it survives this normalisation.
#
# It also emits a HPV-type x tier catalogue — the original descriptive question
# ("what HPV types are actually present") sharpened by population, which is what
# shows the contamination types (HPV16/18) sitting in cultures, not tissue.
#
# Every --*-agg and --*-sheet option accepts a COMMA-SEPARATED LIST, so extra
# lesion batches (e.g. the cSCC/AK run) fold in with no code change:
#   bin/hpv_rpm_contrast.R \
#     --lesion-agg   results_lesions_wide/aggregated,results_lesions_csc/aggregated \
#     --lesion-sheet results_lesions_wide/lesions_samplesheet_full.csv,results_lesions_csc/lesions_samplesheet_full.csv
#
# Usage (defaults point at the local result dirs):
#   bin/hpv_rpm_contrast.R [--outdir results_reframe]

suppressPackageStartupMessages({ library(tidyverse); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--lesion-agg",   default = "results_lesions_wide/aggregated"),
  make_option("--lesion-sheet", default = "results_lesions_wide/lesions_samplesheet_full.csv"),
  make_option("--null-agg",     default = "results_full_v2/aggregated"),
  make_option("--null-sheet",   default = "results_full_v2/metadata/samplesheet_enriched.csv"),
  make_option("--null-flags",   default = "results_full_v2/metadata/cell_line_flags_curated.tsv"),
  make_option("--contam-studies", default = "SRP563552,SRP516327"),
  make_option("--outdir",       default = "results_reframe")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
norm_bool <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes")
split_csv <- function(x) keep(str_trim(str_split(x, ",")[[1]]), nzchar)

# ── 1. RPM from Kraken2 reports (every screened sample contributes) ──────────
read_kraken_rpm <- function(f) {
  df <- tryCatch(
    read.delim(f, header = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = ""),
    error = function(e) NULL)
  if (is.null(df) || ncol(df) < 6) return(NULL)
  names(df)[1:6] <- c("pct", "clade", "taxon", "rank", "taxid", "name")
  nm    <- trimws(df$name)
  total <- sum(df$clade[df$rank == "U"]) + sum(df$clade[df$rank == "R"])
  hpv   <- sum(df$clade[nm == "Papillomaviridae"])
  tibble(total_reads = total, hpv_reads = hpv,
         rpm = if (total > 0) hpv / total * 1e6 else NA_real_)
}
gather_rpm <- function(dir, cohort) {
  files <- list.files(dir, pattern = "_kraken2_report\\.txt$", full.names = TRUE)
  message(sprintf("  %s [%s]: %d kraken2 reports", cohort, dir, length(files)))
  map_dfr(files, function(f) {
    r <- read_kraken_rpm(f); if (is.null(r)) return(NULL)
    mutate(r, srr_id = sub("_kraken2_report\\.txt$", "", basename(f)), cohort = cohort)
  })
}
gather_rpm_multi <- function(dirs, cohort)
  map_dfr(split_csv(dirs), ~ gather_rpm(.x, cohort))

message("Reading Kraken2 reports …")
rpm <- bind_rows(gather_rpm_multi(opt$`lesion-agg`, "lesion"),
                 gather_rpm_multi(opt$`null-agg`,   "null")) %>%
  distinct(srr_id, cohort, .keep_all = TRUE)

# ── 2. Sample labels (tissue / diagnosis / cell-line) ───────────────────────
les <- map_dfr(split_csv(opt$`lesion-sheet`), ~
    read_csv(.x, show_col_types = FALSE) %>%
      transmute(srr_id, cohort = "lesion", study, tissue_category, diagnosis,
                is_cell_line  = norm_bool(is_cell_line),
                is_engineered = norm_bool(is_engineered))) %>%
  distinct(srr_id, .keep_all = TRUE)

nullss <- map_dfr(split_csv(opt$`null-sheet`), ~
    read_csv(.x, show_col_types = FALSE) %>%
      select(srr_id, study, tissue_category, diagnosis)) %>%
  distinct(srr_id, .keep_all = TRUE)
flags <- map_dfr(split_csv(opt$`null-flags`), ~
    read_tsv(.x, show_col_types = FALSE) %>%
      transmute(srr_id, is_cell_line = norm_bool(is_cell_line),
                is_engineered = norm_bool(is_engineered))) %>%
  distinct(srr_id, .keep_all = TRUE)
nullss <- nullss %>% left_join(flags, by = "srr_id") %>%
  mutate(cohort = "null",
         is_cell_line  = coalesce(is_cell_line,  FALSE),
         is_engineered = coalesce(is_engineered, FALSE))

meta <- bind_rows(les, nullss)

# ── 3. Tier assignment (cultures & contamination subtracted from the null) ───
contam <- str_trim(str_split(opt$`contam-studies`, ",")[[1]])
dat <- rpm %>% left_join(meta, by = c("srr_id", "cohort"))

dat <- dat %>% mutate(
  dx = tolower(coalesce(diagnosis, "")),
  tc = tolower(coalesce(tissue_category, "")),
  tier = case_when(
    coalesce(is_cell_line, FALSE)                 ~ "cell_line",
    str_detect(dx, "verruciform")                 ~ "EV",
    str_detect(dx, "wart|condyloma")              ~ "control_productive",
    str_detect(dx, "squamous|carcinoma|scc|keratos|bowen") ~ "neoplasia",
    study %in% contam                             ~ "contamination_study",
    tc == "nahk"                                  ~ "unselected_skin",
    TRUE                                          ~ "other"
  ))

tier_levels <- c("unselected_skin", "contamination_study", "cell_line",
                 "control_productive", "EV", "neoplasia", "other")
dat <- dat %>% mutate(tier = factor(tier, levels = tier_levels))

# Detection floor calibrated on the productive controls: the 5th-percentile RPM
# of HPV-defined wart/condyloma tissue is the level we are demonstrably
# sensitive to. Anything below it in the null is "not detected at control level".
ctrl_rpm <- dat %>% filter(tier == "control_productive", !is_cell_line, rpm > 0) %>% pull(rpm)
floor_rpm <- if (length(ctrl_rpm)) as.numeric(quantile(ctrl_rpm, 0.05, na.rm = TRUE)) else NA_real_

# ── 4. Outputs ──────────────────────────────────────────────────────────────
write_tsv(dat %>% select(srr_id, cohort, study, tissue_category, diagnosis,
                         is_cell_line, is_engineered, total_reads, hpv_reads, rpm, tier),
          file.path(opt$outdir, "hpv_rpm_per_sample.tsv"))

tier_summary <- dat %>% group_by(tier) %>% summarise(
  n            = n(),
  median_total = median(total_reads, na.rm = TRUE),
  pct_any_hpv  = round(mean(rpm > 0, na.rm = TRUE) * 100, 1),
  pct_ge_floor = round(mean(rpm >= floor_rpm, na.rm = TRUE) * 100, 1),
  median_rpm   = round(median(rpm, na.rm = TRUE), 3),
  p95_rpm      = round(quantile(rpm, 0.95, na.rm = TRUE), 1),
  .groups = "drop") %>% arrange(tier)
write_tsv(tier_summary, file.path(opt$outdir, "rpm_by_tier.tsv"))

# Depth-normalised contrast plot
# Human-readable tier labels (match the manuscript wording); str_wrap keeps the
# longer ones to a few short lines so they sit horizontally under each box.
tier_display <- c(
  unselected_skin     = "Unselected skin",
  contamination_study = "Contamination studies",
  cell_line           = "Cell line",
  control_productive  = "Productive controls (warts/condyloma)",
  EV                  = "Epidermodysplasia verruciformis",
  neoplasia           = "Neoplasia (cSCC/AK)",
  other               = "Other tissues")
pdat <- dat %>% filter(!is.na(tier)) %>%
  mutate(rpm_plot = pmax(rpm, 0) + 0.05)
p <- ggplot(pdat, aes(tier, rpm_plot)) +
  geom_jitter(aes(colour = tier), width = 0.25, height = 0, alpha = 0.4, size = 0.7) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
  { if (!is.na(floor_rpm)) geom_hline(yintercept = floor_rpm + 0.05,
        linetype = "dashed", colour = "red") } +
  scale_y_log10() +
  scale_x_discrete(labels = function(lv) str_wrap(tier_display[lv], width = 16)) +
  labs(title = "Depth-normalised HPV transcript burden by population",
       subtitle = sprintf("Papillomaviridae reads/M; red = control-calibrated detection floor (%.1f RPM)", floor_rpm),
       x = NULL, y = "HPV reads per million (log10, +0.05)") +
  theme_minimal() + theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(file.path(opt$outdir, "hpv_rpm_contrast.pdf"), p, width = 9, height = 5.5)
ggsave(file.path(opt$outdir, "hpv_rpm_contrast.png"), p, width = 9, height = 5.5, dpi = 150)

# ── 4b. Depth-matched detection (pre-empt "it is still a depth effect") ─────
# Within each total-read bin, productive controls should still detect at ~100%
# while unselected skin stays near zero — i.e. the gap is not explained by depth.
depth_breaks <- c(0, 5e6, 1e7, 2e7, 4e7, Inf)
depth_labs   <- c("<5M", "5-10M", "10-20M", "20-40M", ">40M")
dmatch <- dat %>%
  filter(!is.na(tier), !is.na(total_reads), total_reads > 0) %>%
  mutate(depth_bin = cut(total_reads, breaks = depth_breaks, labels = depth_labs)) %>%
  group_by(tier, depth_bin) %>%
  summarise(n = n(),
            pct_ge_floor = round(mean(rpm >= floor_rpm, na.rm = TRUE) * 100, 1),
            median_rpm   = round(median(rpm, na.rm = TRUE), 2),
            .groups = "drop")
write_tsv(dmatch, file.path(opt$outdir, "depth_matched_detection.tsv"))

key_tiers <- c("unselected_skin", "control_productive", "EV", "neoplasia")
pd <- dmatch %>% filter(tier %in% key_tiers, n >= 3)
if (nrow(pd)) {
  pdepth <- ggplot(pd, aes(depth_bin, pct_ge_floor, colour = tier, group = tier)) +
    geom_line() + geom_point(aes(size = n)) +
    labs(title = "Control-level HPV detection within matched sequencing-depth bins",
         subtitle = sprintf("%% of samples >= %.1f RPM floor; flat-high controls vs flat-low skin = gap is not depth", floor_rpm),
         x = "Total reads (Kraken2-screened)", y = "% at/above detection floor", size = "n") +
    theme_minimal()
  ggsave(file.path(opt$outdir, "depth_matched_detection.pdf"), pdepth, width = 8, height = 5)
  ggsave(file.path(opt$outdir, "depth_matched_detection.png"), pdepth, width = 8, height = 5, dpi = 150)
}

# ── 5. HPV-type catalogue by tier (the original descriptive question) ───────
gather_types <- function(dir, cohort) {
  files <- list.files(dir, pattern = "_hpv_types\\.tsv$", full.names = TRUE)
  map_dfr(files, ~ suppressMessages(
            read_tsv(.x, show_col_types = FALSE,
                     col_types = cols(.default = col_character())))) %>%
    mutate(cohort = cohort)
}
gather_types_multi <- function(dirs, cohort)
  map_dfr(split_csv(dirs), ~ gather_types(.x, cohort))
types <- bind_rows(gather_types_multi(opt$`lesion-agg`, "lesion"),
                   gather_types_multi(opt$`null-agg`,   "null")) %>%
  left_join(dat %>% select(srr_id, cohort, tier), by = c("sample_id" = "srr_id", "cohort"))

type_catalogue <- types %>%
  count(hpv_reference, tier, name = "n_samples") %>%
  pivot_wider(names_from = tier, values_from = n_samples, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total))
write_tsv(type_catalogue, file.path(opt$outdir, "type_catalogue_by_tier.tsv"))

# ── 6. Console headline ─────────────────────────────────────────────────────
cat("\n================  RPM CONTRAST  ================\n")
print(tier_summary, n = Inf)
cat(sprintf("\nControl-calibrated detection floor: %.2f RPM\n", floor_rpm))
cat("\n================  DEPTH-MATCHED DETECTION (% >= floor)  ================\n")
print(dmatch %>%
        filter(tier %in% key_tiers) %>%
        select(tier, depth_bin, n, pct_ge_floor) %>%
        pivot_wider(names_from = depth_bin, values_from = c(n, pct_ge_floor)),
      n = Inf, width = Inf)
cat("\n================  TOP HPV TYPES x TIER  ================\n")
print(head(type_catalogue, 20), n = 20)
cat(sprintf("\nWrote: %s\n", normalizePath(opt$outdir)))
