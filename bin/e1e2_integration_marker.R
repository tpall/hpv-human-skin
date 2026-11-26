#!/usr/bin/env Rscript
# E1:E2 transcript ratio as a metadata-independent marker of HPV integration.
#
# Rationale (collaborator observation, 2026-08-03): E1 and E2 are encoded on the
# same early polycistronic transcript, so their read counts should track each
# other. They do not. Integration of the HPV genome characteristically disrupts
# the E2 ORF while leaving E1 intact, so integrated material loses E2 transcripts
# and the ratio inverts. That makes E1:E2 an orthogonal readout of integration --
# and therefore of cell-line contamination -- computed from the reads themselves,
# independent of the name- and metadata-based cell-line flagging the manuscript
# currently relies on.
#
# Two criteria, because integration presents in two ways:
#   (a) E2 DISRUPTED, E1 retained -> the E1:E2 ratio inverts. This is the
#       original marker and it works where both genes are still transcribed.
#   (b) EARLY REGION LOST except E6/E7 -> E1 and E2 are BOTH absent, so the
#       ratio is undefined and (a) is blind. The signature is instead E6/E7
#       transcription with no late (L1/L2) transcription at all, the canonical
#       integrated-carcinoma pattern. Found when the E1:E2 marker silently
#       dropped 52 HPV16 cervical-line libraries into "insufficient".
# A library is called integration-like if EITHER criterion fires.
#
# Outputs (into --outdir):
#   e1e2_by_library.tsv   per-library ratio, tier, flags, dominant type
#   e1e2_by_tier.tsv      tier-level summary + separation statistic
#   e1e2_discordant.tsv   libraries whose ratio contradicts their flags
#   figure_e1e2.png/pdf   ratio distribution by tier
#
#   bin/e1e2_integration_marker.R [--outdir results_e1e2]

suppressPackageStartupMessages({
  library(tidyverse); library(here); library(optparse); library(cowplot)
})
theme_set(theme_cowplot(font_size = 11))
source(here("manuscript/_paths.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--outdir", default = "results_e1e2"),
  make_option("--min-early", type = "integer", default = 5,
              help = "minimum E1+E2 reads for a library to yield a ratio [%default]"),
  make_option("--cut", type = "double", default = 2,
              help = "log2(E1/E2) above which a library is called integration-like [%default]"),
  make_option("--max-late-frac", type = "double", default = 0.05,
              help = "max share of viral reads that may be L1/L2 for the E6/E7-only call [%default]"),
  make_option("--agg", default = NULL,
              help = "comma-separated result dirs to read aggregated/ from; default is the main cohort"),
  make_option("--no-tiers", action = "store_true", default = FALSE,
              help = "cohort has no tier table (e.g. a standalone typing run)"),
  make_option("--groups", default = NULL,
              help = paste("TSV with columns sample_id/srr_id + group; used in place of",
                           "tiers so an external cohort can be summarised by its own design",
                           "arms (e.g. the cervical normal->HSIL->SCC series)")),
  make_option("--keep-single-cell", action = "store_true", default = FALSE,
              help = paste("keep single-cell/single-nucleus libraries (default: drop).",
                           "3'-biased barcoded reads cannot carry an E1:E2 ratio,",
                           "so including them adds noise, not signal"))
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

AGG <- if (is.null(opt$agg)) {
  c("results_full_v2", "results_lesions", "results_lesions_csc", "results_lesions_wide")
} else {
  str_trim(str_split(opt$agg, ",")[[1]])
}

# ── per-gene read counts ────────────────────────────────────────────────────
read_glob <- function(pattern) {
  fs <- unlist(lapply(AGG, function(d) list.files(here(d, "aggregated"), pattern, full.names = TRUE)))
  map_dfr(fs, ~ suppressMessages(read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c"))))
}
tc <- read_glob("_transcript_classes\\.tsv$") %>%
  filter(gene %in% c("E1", "E2", "E6", "E7", "L1", "L2")) %>%
  mutate(read_count = suppressWarnings(as.integer(read_count))) %>%
  filter(!is.na(read_count)) %>%
  distinct(sample_id, gene, .keep_all = TRUE) %>%
  pivot_wider(id_cols = sample_id, names_from = gene, values_from = read_count, values_fill = 0)

# Dominant assigned type, to test whether the effect is type-driven (HeLa/SiHa
# are HPV18/16) rather than integration-driven.
dom <- read_glob("_hpv_types\\.tsv$") %>%
  mutate(coverage_breadth = as.numeric(coverage_breadth), mean_depth = as.numeric(mean_depth)) %>%
  filter(coverage_breadth >= 0.10, mean_depth >= 2) %>%
  group_by(sample_id) %>% slice_max(mean_depth, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(sample_id, dom_type = str_remove(hpv_reference, "REF$"))

rp <- if (!is.null(opt$groups)) {
  # External cohort with its own design arms. Using them as `tier` means the
  # existing by-tier summary, separation statistic and figure all work unchanged,
  # and the arms come from the submitting study rather than from our tier rules --
  # so the marker is being tested against labels it had no part in assigning.
  g <- read_tsv(opt$groups, show_col_types = FALSE)
  id <- intersect(c("sample_id", "srr_id"), names(g))[1]
  if (is.na(id)) stop("--groups file needs a sample_id or srr_id column")
  if (!"group" %in% names(g)) stop("--groups file needs a `group` column")
  tibble(sample_id = g[[id]], study = NA_character_, tier = g$group,
         is_cell_line = FALSE, total_reads = NA_real_, rpm = NA_real_)
} else if (opt$`no-tiers`) {
  # standalone cohort: no tier assignment exists, so score every typed library
  tibble(sample_id = tc$sample_id, study = NA_character_, tier = "unassigned",
         is_cell_line = FALSE, total_reads = NA_real_, rpm = NA_real_)
} else {
  read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>%
    select(sample_id = srr_id, study, tier, is_cell_line, total_reads, rpm)
}

# Engineered flag is NOT consumed by the tier rules, so carry it explicitly.
flags <- tryCatch(
  read_tsv(here("results_full_v2/metadata/cell_line_flags_curated.tsv"),
           show_col_types = FALSE) %>%
    transmute(sample_id = srr_id,
              is_engineered = toupper(as.character(is_engineered)) == "TRUE"),
  error = function(e) tibble(sample_id = character(), is_engineered = logical()))

# Single-cell / single-nucleus libraries: 3'-biased and barcoded, so coverage
# is concentrated at the transcript end rather than spread across the early
# region. The E1:E2 ratio is then a statement about where the protocol put its
# reads, not about whether E2 is disrupted. Dropped by default; --keep-single-cell
# restores the old behaviour for comparison. Flag comes from parse_metadata.py.
sc <- tryCatch(
  read_csv(here("results_full_v2/metadata/samplesheet_enriched.csv"),
           show_col_types = FALSE, progress = FALSE) %>%
    transmute(sample_id = srr_id,
              is_single_cell = tolower(as.character(is_single_cell)) == "true") %>%
    distinct(sample_id, .keep_all = TRUE),
  error = function(e) tibble(sample_id = character(), is_single_cell = logical()))

lib <- tc %>%
  inner_join(rp, by = "sample_id") %>%
  left_join(dom, by = "sample_id") %>%
  left_join(flags, by = "sample_id") %>%
  left_join(sc, by = "sample_id") %>%
  mutate(is_single_cell = coalesce(is_single_cell, FALSE)) %>%
  { if (opt$`keep-single-cell`) . else {
      n_sc <- sum(.$is_single_cell)
      if (n_sc > 0) cat(sprintf(
        "Dropped %d single-cell/single-nucleus libraries (--keep-single-cell to retain)\n", n_sc))
      filter(., !is_single_cell)
    } } %>%
  mutate(is_engineered = coalesce(is_engineered, FALSE),
         early     = E1 + E2,
         oncogene  = E6 + E7,
         late      = L1 + L2,
         viral     = E1 + E2 + E6 + E7 + L1 + L2,
         late_frac = if_else(viral > 0, late / viral, NA_real_),
         log2_e1e2 = if_else(early >= opt$`min-early`,
                             log2((E1 + 0.5) / (E2 + 0.5)), NA_real_),
         # (a) E2 disrupted, E1 retained
         crit_ratio = !is.na(log2_e1e2) & log2_e1e2 >= opt$cut,
         # (b) early region lost but for E6/E7: oncogenes transcribed, no late
         #     transcription, and no usable E1/E2 signal to compute a ratio from
         crit_e67 = oncogene >= opt$`min-early` & early < opt$`min-early` &
                    !is.na(late_frac) & late_frac <= opt$`max-late-frac`,
         integration_like = crit_ratio | crit_e67,
         call = case_when(crit_ratio ~ "ratio (E2 disrupted)",
                          crit_e67   ~ "E6/E7-only (early region lost)",
                          !is.na(log2_e1e2) ~ "episomal-like",
                          TRUE ~ "uninformative")) %>%
  filter(viral >= opt$`min-early`)

write_tsv(lib %>% arrange(desc(log2_e1e2)) %>%
            select(sample_id, study, tier, is_cell_line, is_engineered, dom_type,
                   E1, E2, E6, E7, L1, L2, late_frac, log2_e1e2,
                   crit_ratio, crit_e67, integration_like, call, total_reads, rpm),
          file.path(opt$outdir, "e1e2_by_library.tsv"))

# ── tier summary ────────────────────────────────────────────────────────────
tiersum <- lib %>% group_by(tier) %>%
  summarise(n = n(), median_log2 = round(median(log2_e1e2, na.rm = TRUE), 2),
            fold = round(2^median(log2_e1e2, na.rm = TRUE), 2),
            n_ratio = sum(crit_ratio), n_e67 = sum(crit_e67),
            pct_integration_like = round(100 * mean(integration_like), 1),
            .groups = "drop") %>% arrange(desc(pct_integration_like))
write_tsv(tiersum, file.path(opt$outdir, "e1e2_by_tier.tsv"))
cat("\n=== E1:E2 by tier ===\n"); print(as.data.frame(tiersum))

# Separation: cell lines vs productive controls (AUC of the ranked ratio).
a <- lib$log2_e1e2[lib$tier == "cell_line"]
b <- lib$log2_e1e2[lib$tier == "control_productive"]
if (length(a) && length(b)) {
  auc <- mean(outer(a, b, ">") + 0.5 * outer(a, b, "=="))
  w <- suppressWarnings(wilcox.test(a, b))
  cat(sprintf("\ncell_line vs control_productive: AUC = %.3f, Wilcoxon P = %.3g\n", auc, w$p.value))
}

# Type confound: is the split explained by HPV16/18 rather than by integration?
cat("\n=== median log2(E1/E2) by dominant type (>=3 libraries) ===\n")
print(as.data.frame(lib %>% filter(!is.na(dom_type)) %>% group_by(dom_type) %>%
        filter(n() >= 3) %>%
        summarise(n = n(), median_log2 = round(median(log2_e1e2), 2),
                  n_cell_line = sum(is_cell_line), .groups = "drop") %>%
        arrange(desc(median_log2))))

# ── discordant libraries: the actionable output ─────────────────────────────
disc <- lib %>%
  filter(integration_like, !is_cell_line) %>%
  arrange(desc(coalesce(log2_e1e2, 99))) %>%
  select(sample_id, study, tier, is_engineered, dom_type, call,
         E1, E2, E6, E7, L1, L2, log2_e1e2)
write_tsv(disc, file.path(opt$outdir, "e1e2_discordant.tsv"))
cat(sprintf("\n=== %d integration-like libraries NOT flagged as cell line ===\n", nrow(disc)))
print(as.data.frame(disc %>% count(study, tier, is_engineered, sort = TRUE) %>% head(12)))

# ── figure ──────────────────────────────────────────────────────────────────
ord <- tiersum$tier
p <- lib %>% mutate(tier = factor(tier, ord)) %>%
  ggplot(aes(tier, log2_e1e2)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_hline(yintercept = opt$cut, linetype = 3, colour = "firebrick") +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey92") +
  geom_jitter(aes(colour = is_cell_line), width = 0.18, height = 0, size = 1.4, alpha = 0.8) +
  scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "firebrick"),
                      labels = c("primary / unflagged", "flagged cell line"), name = NULL) +
  labs(x = NULL, y = expression(log[2]~"(E1 / E2 reads)")) +
  coord_flip() + theme(legend.position = "bottom")
ggsave(file.path(opt$outdir, "figure_e1e2.png"), p, width = 8, height = 5, dpi = 200)
ggsave(file.path(opt$outdir, "figure_e1e2.pdf"), p, width = 8, height = 5)
cat(sprintf("\nWrote %s/{e1e2_by_library,e1e2_by_tier,e1e2_discordant}.tsv + figure_e1e2.{png,pdf}\n", opt$outdir))
