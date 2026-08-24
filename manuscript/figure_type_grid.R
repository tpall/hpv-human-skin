#!/usr/bin/env Rscript
# Builds manuscript/figureS2.{png,pdf}: the HPV type by tier grid, assigned
# against detection tier. Laid out landscape, types across and tiers down, so
# each tier reads as one line: the unselected-skin line of the detection panel
# is then the reservoir claim directly, a long run of pale cells across beta and
# gamma. Portrait it took 13 inches and the same claim had to be counted down a
# column.
#
#   Rscript manuscript/figure_type_grid.R
#
# WHY THIS FIGURE EXISTS. Tables S1 and S6 partition one cohort: a library is
# either typed above the assignment criterion or resolved below it, never both.
# Read as two tables that is hard to see, and the thing worth seeing is exactly
# the complementarity, because the healthy-skin beta/gamma reservoir lives
# almost entirely in the second one. Laid out as a single grid with genus as
# rows and the two catalogues as columns, that partition is the picture.
#
# COLOUR. Fill is log10, not linear: 37% of cells carry three or fewer
# libraries, and on a linear ramp those sit at the white end where the count
# labels become unreadable. The palette starts short of white for the same
# reason, and each label takes its colour from its own tile.
suppressPackageStartupMessages({library(tidyverse); library(here)})
old <- theme_set(theme_minimal(base_size = 10))
source(here("manuscript/_paths.R"))

TLAB <- c(unselected_skin = "Unselected skin", neoplasia = "Neoplasia", EV = "EV",
          control_productive = "Productive controls", peri_lesional = "Peri-lesional",
          cell_line = "Cell line", engineered = "Engineered",
          non_permissive_tissue = "Non-permissive", other = "Site not recorded")
TORD <- names(TLAB)
GLEV <- c("alpha", "beta", "gamma", "mu", "nu")

genus <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
  transmute(ref = hpv_reference, genus = genus_short)
tiers <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% select(srr_id, tier)
# Tier sizes are the denominators, and they span 10 to 7,800. Fill therefore
# encodes prevalence within tier rather than the raw count: a count of 3 is 30%
# of the EV tier and 0.04% of unselected skin, so a count-filled heatmap
# compares quantities that are not comparable, which is the very point Table S1's
# own caption makes. The count stays as the cell label.
tier_n <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>%
  count(tier, name = "tier_size")

assigned <- read_tsv(type_catalogue_tsv(), show_col_types = FALSE) %>%
  pivot_longer(-c(hpv_reference, total), names_to = "tier", values_to = "n") %>%
  filter(n > 0, tier %in% TORD) %>%
  left_join(genus, by = c("hpv_reference" = "ref")) %>%
  transmute(tier, type = sub("REF$", "", hpv_reference), genus, n)

detected <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types = FALSE) %>%
  filter(!top_ref %in% c("NONE", "NO_BAM"), top_reads >= 5, max_mapq >= 20,
         top_reads / pmax(n_primary, 1) >= 0.5) %>%
  left_join(tiers, by = "srr_id") %>% filter(tier %in% TORD) %>%
  count(tier, top_ref, name = "n") %>%
  left_join(genus, by = c("top_ref" = "ref")) %>%
  transmute(tier, type = sub("REF$", "", top_ref), genus, n)

both <- bind_rows(assigned = assigned, detected = detected, .id = "id") %>%
  mutate(genus = factor(replace_na(genus, "gamma"), levels = GLEV))
TYPE_ORDER <- both %>% group_by(type, genus) %>%
  summarise(tot = sum(n), .groups = "drop") %>%
  arrange(genus, desc(tot)) %>% pull(type)

d <- both %>%
  mutate(id = factor(id, levels = c("assigned", "detected"),
                     labels = c("Assigned", "Detection tier")),
         tier = factor(tier, levels = TORD)) %>%
  complete(id, tier, nesting(genus, type)) %>%
  left_join(tier_n %>% mutate(tier = factor(tier, levels = TORD)), by = "tier") %>%
  # Beta-Binomial posterior mean under Beta(0.5, 0.5), the same estimator the
  # tables use, so figure and table report one quantity rather than two.
  mutate(pct = 100 * (n + 0.5) / (tier_size + 1),
         type = factor(type, levels = TYPE_ORDER))
LIM <- range(d$pct, na.rm = TRUE)
d <- d %>%
  mutate(lab = ifelse(is.na(n), "", formatC(n, format = "d")),
         # Label colour from the tile's own position on the log ramp, so counts
         # stay legible at both ends instead of vanishing into a pale fill.
         lab_col = ifelse(is.na(pct) |
                          (log10(pct) - log10(LIM[1])) / (log10(LIM[2]) - log10(LIM[1])) < 0.55,
                          "grey15", "white"))

p <- ggplot(d, aes(type, fct_rev(tier), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = lab, colour = lab_col), size = 2.6, show.legend = FALSE) +
  facet_grid(id ~ genus, space = "free_x", scales = "free_x") +
  scale_fill_viridis_c("Prevalence\nwithin tier (%)", option = "mako", direction = -1,
                       trans = "log10", limits = LIM, begin = 0.12, end = 0.95,
                       breaks = c(0.03, 0.3, 3, 30),
                       labels = function(x) formatC(x, format = "g"),
                       na.value = "grey94") +
  scale_colour_identity() +
  scale_y_discrete(labels = function(x)
    sprintf("%s (n = %s)", TLAB[x],
            formatC(tier_n$tier_size[match(x, tier_n$tier)], format = "d", big.mark = ","))) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        axis.text.y = element_text(size = 8),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0), strip.text = element_text(size = 9),
        panel.spacing.x = unit(2, "pt"), legend.key.height = unit(1.1, "cm"))

ggsave(here("manuscript/figureS2.png"), p, width = 13.5, height = 6.2, dpi = 200, bg = "white")
ggsave(here("manuscript/figureS2.pdf"), p, width = 13.5, height = 6.2)
cat(sprintf("figureS2: %d types, %d assigned cells, %d detection cells; prevalence %.3f%% to %.1f%%\n",
            n_distinct(both$type), nrow(assigned), nrow(detected), LIM[1], LIM[2]))
