#!/usr/bin/env Rscript
# Builds manuscript/figure4.png: HPV genus prevalence by population tier.
#
#   Rscript manuscript/figure_genus.R
#
# Collapses the two type catalogues (Table S1, assigned types; Table S6,
# detection-tier trace types) from HPV type to HPV genus, so the type landscape
# can be read at a glance. Three rows:
#   Assigned  - library has a reference passing breadth >= 0.10 and depth >= 2
#   Detection - HPV-positive, no reference passes, but a best-hit type resolves
#   Any       - union of the two (a library appears in at most one of them)
# A library can contribute to more than one genus if it carries more than one
# type, so the bars within a tier need not sum to the tier's positivity rate.
suppressPackageStartupMessages({library(tidyverse); library(cowplot); library(here)})
source(here("manuscript/_paths.R"))
old <- theme_set(theme_cowplot())

AGG   <- c("results_full_v2", "results_lesions", "results_lesions_csc", "results_lesions_wide")
TIERS <- c("unselected_skin", "control_productive", "EV",
           "non_permissive_tissue", "engineered",
           "neoplasia", "cell_line", "other")
GEN   <- c("alpha", "beta", "gamma", "mu", "nu")
TLAB  <- c(unselected_skin = "Unselected\nskin", control_productive = "Productive\ncontrols",
           EV = "EV", neoplasia = "Neoplasia",
           # Broken over three lines, not two: the strip is the width constraint
           # (see the theme note below), and "Non-permissive" on one line sets the
           # minimum strip width for the whole row.
           non_permissive_tissue = "Non-\npermissive\ntissue", engineered = "Engineered\ncells",
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

# Detection-tier types (the Table S6 population): same filter as build_trace_table().
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
    summarise_level(trace,    "Detection level (Table S6)"),
    summarise_level(anyd,     "Any detection")) %>%
  mutate(tier  = factor(TLAB[tier], levels = TLAB[TIERS]),
         Genus = factor(str_to_title(genus_short), levels = str_to_title(GEN)),
         level = factor(level, levels = c("Assigned type (Table S1)",
                                          "Detection level (Table S6)", "Any detection")),
         # A zero count still has a small positive posterior mean. Drawing it
         # would read as a detection, so zeros are marked at the axis floor.
         mean_p = if_else(k > 0, mean, NA_real_),
         lo_p   = if_else(k > 0, lo,   NA_real_),
         hi_p   = if_else(k > 0, hi,   NA_real_))

FLOOR <- 0.01   # log axis needs a positive stand-in for "not detected"
# Colour is deliberately not mapped. Genus is already the x axis, so colouring
# by it encodes nothing a reader cannot read off the position, which is why it
# carried no legend; a redundant scale without a key is noise rather than
# information. One ink colour throughout, with the grey crosses reserved for
# "not detected" so that distinction keeps a channel of its own.
INK <- "#2c6fbb"
p <- ggplot(pd, aes(Genus, mean_p)) +
  geom_errorbar(aes(ymin = lo_p, ymax = hi_p), width = .18, linewidth = .4,
                colour = INK, na.rm = TRUE) +
  geom_point(size = 2, colour = INK, na.rm = TRUE) +
  geom_point(data = ~ filter(.x, k == 0), aes(y = FLOOR), shape = 4, size = 1.4,
             colour = "grey55", na.rm = TRUE) +
  geom_text(data = ~ filter(.x, k > 0), aes(y = hi_p, label = k), vjust = -.6,
            colour = INK,
            size = 3, show.legend = FALSE) +
  # Strips name the level and nothing else. A supplementary table number burned
  # into the image cannot be corrected at proof and goes silently wrong whenever
  # an item is inserted ahead of it, which is exactly what happened when the
  # chemistry table took S5. The caption cross-references both tables instead,
  # and those renumber themselves.
  facet_grid(level ~ tier, labeller = labeller(level = c(
    "Assigned type (Table S1)"  = "Assigned",
    "Detection level (Table S6)" = "Detection level",
    "Any detection"             = "Any"))) +
  scale_y_log10(breaks = c(.01, .1, 1, 10, 100), labels = c("0", "0.1", "1", "10", "100"),
                limits = c(FLOOR, 400)) +
  labs(x = "HPV genus", y = "Prevalence within sample group (%, log scale)") +
  # theme_cowplot() is set globally above. Two overrides matter here.
  #
  # The genus labels are angled, since five of them per facet collide across
  # eight columns.
  #
  # Text size is set for the size this figure is READ at, not the size it is
  # drawn at (TP comment 25, 2026-11-02: "strip text too small, difficult to
  # read"). Word scales the image down to the ~6.5in text width, so on-page
  # point size is size * 6.5 / width: at the old 9 pt across 10.5in the strips
  # printed at about 5.6 pt. At 12 pt across 8in they print at about 9.8 pt.
  # Column-strip width is the constraint on how far this can go: 8in over eight
  # facets leaves about 1in per strip, which is why "Non-permissive" is broken
  # across lines above rather than set on one, and why 11 pt rather than 12 pt:
  # at 12 pt "Unselected" and "Engineered" were clipped at the strip edges.
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y  = element_text(size = 11),
        axis.title   = element_text(size = 13),
        strip.text.x = element_text(size = 11, lineheight = 0.9),
        strip.text.y = element_text(size = 11))

ggsave(here("manuscript/figure4.png"), p, width = 8, height = 6.0, dpi = 200)
ggsave(here("manuscript/figure4.pdf"), p, width = 8, height = 6.0)

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
cat("Wrote manuscript/figure4.{png,pdf} and genus_prevalence_by_tier.tsv\n")
