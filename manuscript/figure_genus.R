#!/usr/bin/env Rscript
# Builds manuscript/figure3.png: HPV genus prevalence by population tier.
#
#   Rscript manuscript/figure_genus.R
#
# Collapses the two type catalogues (Table S1, assigned types; Table S5,
# detection-tier trace types) from HPV type to HPV genus, so the type landscape
# can be read at a glance. Three rows:
#   Assigned  - library has a reference passing breadth >= 0.10 and depth >= 2
#   Trace     - HPV-positive, no reference passes, but a best-hit type resolves
#   Any       - union of the two (a library appears in at most one of them)
# A library can contribute to more than one genus if it carries more than one
# type, so the bars within a tier need not sum to the tier's positivity rate.
suppressPackageStartupMessages({library(tidyverse); library(cowplot); library(here)})
source(here("manuscript/_paths.R"))
old <- theme_set(theme_cowplot())

AGG   <- c("results_full_v2", "results_lesions", "results_lesions_csc", "results_lesions_wide")
TIERS <- c("unselected_skin", "control_productive", "EV", "peri_lesional",
           "non_permissive_tissue", "engineered",
           "neoplasia", "cell_line", "other")
GEN   <- c("alpha", "beta", "gamma", "mu", "nu")
TLAB  <- c(unselected_skin = "Unselected\nskin", control_productive = "Productive\ncontrols",
           EV = "EV", peri_lesional = "Peri-lesional\nmucosa", neoplasia = "Neoplasia",
           non_permissive_tissue = "Non-permissive\ntissue", engineered = "Engineered\ncells",
           cell_line = "Cell line",
           # Matches bin/hpv_rpm_contrast.R and the Results text. 99% of this
           # tier carries no site keyword in its SRA metadata rather than an
           # identified non-cutaneous site, so "Other" asserted something the
           # metadata never established.
           other = "Site not\nrecorded")

rp <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE)
tierN <- rp %>% count(tier) %>% { setNames(.$n, .$tier) }

gmap <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
  mutate(type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  select(type, genus_short) %>% distinct(type, .keep_all = TRUE)

# Assigned types (the Table S1 population): one row per library x genus.
assigned <- map_dfr(unlist(lapply(AGG, function(d)
    list.files(here(d, "aggregated"), "_hpv_types\\.tsv$", full.names = TRUE))),
  ~ suppressMessages(read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))) %>%
  distinct(sample_id, hpv_reference) %>%
  mutate(type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  left_join(gmap, by = "type") %>%
  left_join(rp %>% select(srr_id, tier), by = c("sample_id" = "srr_id")) %>%
  distinct(id = sample_id, tier, genus_short)

# Trace types (the Table S5 population): same quality filter as build_trace_table().
trace <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types = FALSE) %>%
  filter(!top_ref %in% c("NONE", "NO_BAM"), top_reads >= 5, max_mapq >= 20,
         top_reads / pmax(n_primary, 1) >= 0.5) %>%
  mutate(type = gsub("(REF|nr)$", "", top_ref)) %>%
  left_join(gmap, by = "type") %>%
  left_join(rp %>% select(srr_id, tier), by = "srr_id") %>%
  distinct(id = srr_id, tier, genus_short)

anyd <- bind_rows(assigned, trace) %>% distinct(id, tier, genus_short)

# Beta-Binomial posterior under the Jeffreys Beta(0.5, 0.5) prior, as in Tables
# S1 and S5, so figure and tables quote the same prevalence for the same cell.
summarise_level <- function(d, level) {
  d %>% filter(genus_short %in% GEN) %>% count(tier, genus_short, name = "k") %>%
    complete(tier = TIERS, genus_short = GEN, fill = list(k = 0)) %>%
    filter(!is.na(tier)) %>%
    mutate(N = tierN[tier], a = k + 0.5, b = N - k + 0.5,
           mean = 100 * a / (a + b), lo = 100 * qbeta(.025, a, b), hi = 100 * qbeta(.975, a, b),
           level = level)
}
pd <- bind_rows(
    summarise_level(assigned, "Assigned type (Table S1)"),
    summarise_level(trace,    "Trace type (Table S5)"),
    summarise_level(anyd,     "Any detection")) %>%
  mutate(tier  = factor(TLAB[tier], levels = TLAB[TIERS]),
         Genus = factor(str_to_title(genus_short), levels = str_to_title(GEN)),
         level = factor(level, levels = c("Assigned type (Table S1)",
                                          "Trace type (Table S5)", "Any detection")),
         # A zero count still has a small positive posterior mean. Drawing it
         # would read as a detection, so zeros are marked at the axis floor.
         mean_p = if_else(k > 0, mean, NA_real_),
         lo_p   = if_else(k > 0, lo,   NA_real_),
         hi_p   = if_else(k > 0, hi,   NA_real_))

FLOOR <- 0.01   # log axis needs a positive stand-in for "not detected"
p <- ggplot(pd, aes(Genus, mean_p, colour = Genus)) +
  geom_errorbar(aes(ymin = lo_p, ymax = hi_p), width = .18, linewidth = .4, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  geom_point(data = ~ filter(.x, k == 0), aes(y = FLOOR), shape = 4, size = 1.4,
             colour = "grey55", na.rm = TRUE) +
  geom_text(data = ~ filter(.x, k > 0), aes(y = hi_p, label = k), vjust = -.6,
            size = 2.3, show.legend = FALSE) +
  # Short row-strip labels so they fit the row height at the cowplot font size;
  # the data keys in pd$level stay long (the caption and _helpers.R use them).
  facet_grid(level ~ tier, labeller = labeller(level = c(
    "Assigned type (Table S1)" = "Assigned (S1)",
    "Trace type (Table S5)"    = "Trace (S5)",
    "Any detection"            = "Any (S1 ∪ S5)"))) +
  scale_y_log10(breaks = c(.01, .1, 1, 10, 100), labels = c("0", "0.1", "1", "10", "100"),
                limits = c(FLOOR, 400)) +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  labs(x = "HPV genus", y = "Prevalence within tier (%, log scale)") +
  # theme_cowplot() is set globally above; the only override is to angle the
  # genus labels, since five of them per facet collide across seven columns.
  # Column-strip labels are the width constraint: nine facets across 10.5in
  # leaves ~1.1in per strip, and at the default cowplot size "Non-permissive"
  # and "Unselected skin" were clipped at both ends rather than wrapped.
  # Every label is now at most two lines, and the strip text is sized to the
  # longest of them.
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 9, lineheight = 0.9),
        strip.text.y = element_text(size = 9))

ggsave(here("manuscript/figure3.png"), p, width = 10.5, height = 6.4, dpi = 200)
ggsave(here("manuscript/figure3.pdf"), p, width = 10.5, height = 6.4)

# The caption quotes prevalences from this figure. Persist them so _helpers.R
# reads the same numbers the figure plots, rather than having them hand-typed.
write_tsv(pd %>% mutate(tier = str_replace_all(as.character(tier), "\n", " ")) %>%
            select(level, tier, genus = Genus, k, N, mean, lo, hi),
          here("manuscript/genus_prevalence_by_tier.tsv"))

cat("\nAny-detection genus prevalence by tier (count, tier N, posterior % [95% CrI]):\n")
pd %>% filter(level == "Any detection", k > 0) %>%
  transmute(tier, Genus, k, N, pct = sprintf("%.2f [%.2f-%.2f]", mean, lo, hi)) %>%
  as.data.frame() %>% print(right = FALSE)
cat("\nGenera never detected in any tier:",
    paste(setdiff(str_to_title(GEN),
                  pd %>% filter(k > 0) %>% pull(Genus) %>% unique() %>% as.character()),
          collapse = ", "), "\n")
cat("Wrote manuscript/figure3.{png,pdf} and genus_prevalence_by_tier.tsv\n")
