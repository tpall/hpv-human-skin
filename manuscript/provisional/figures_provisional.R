#!/usr/bin/env Rscript
# Provisional figures for Tables S1, S6 and S7, for inspection only.
#
#   Rscript manuscript/provisional/figures_provisional.R
#
# These are NOT manuscript figures. They exist to see whether the three
# supplementary tables are better read as pictures, and are deliberately built
# from the same data the tables are built from, so nothing here is a new result.
#
#   prov_S1_assigned_types.png   Table S1, assigned HPV type by tier
#   prov_S6_trace_types.png      Table S6, detection-tier type by tier
#   prov_S1S6_combined.png       the two together, since they partition one cohort
#   prov_S7_coinfection.png      Table S7, co-infected libraries
suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(here); library(patchwork); library(ggrepel)
})
old <- theme_set(theme_cowplot(font_size = 10))
source(here("manuscript/_paths.R"))
OUT <- here("manuscript/provisional")

GEN <- c(alpha = "#c0392b", beta = "#2c6fbb", gamma = "#27ae60",
         mu = "#8e44ad", nu = "#e67e22", other = "grey55")
TLAB <- c(unselected_skin = "Unselected\nskin", neoplasia = "Neoplasia", EV = "EV",
          control_productive = "Productive\ncontrols", peri_lesional = "Peri-\nlesional",
          cell_line = "Cell line", engineered = "Engineered",
          non_permissive_tissue = "Non-\npermissive", other = "Site not\nrecorded")
TORD <- names(TLAB)

genus <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
  transmute(ref = hpv_reference, genus = genus_short)
tiers <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% select(srr_id, tier, study)
tier_n <- tiers %>% count(tier, name = "tier_size")

## ── Table S1: assigned types ────────────────────────────────────────────────
s1 <- read_tsv(type_catalogue_tsv(), show_col_types = FALSE) %>%
  pivot_longer(-c(hpv_reference, total), names_to = "tier", values_to = "n") %>%
  filter(n > 0, tier %in% TORD) %>%
  left_join(genus, by = c("hpv_reference" = "ref")) %>%
  left_join(tier_n, by = "tier") %>%
  mutate(type = sub("REF$", "", hpv_reference),
         genus = replace_na(genus, "other"),
         pct = 100 * n / tier_size,
         tier = factor(tier, levels = TORD))

## ── Table S6: detection-tier (trace) types ──────────────────────────────────
s6 <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types = FALSE) %>%
  filter(!top_ref %in% c("NONE", "NO_BAM"), top_reads >= 5, max_mapq >= 20,
         top_reads / pmax(n_primary, 1) >= 0.5) %>%
  left_join(tiers, by = "srr_id") %>%
  filter(tier %in% TORD) %>%
  count(tier, top_ref, name = "n") %>%
  left_join(genus, by = c("top_ref" = "ref")) %>%
  left_join(tier_n, by = "tier") %>%
  mutate(type = sub("REF$", "", top_ref),
         genus = replace_na(genus, "other"),
         pct = 100 * n / tier_size,
         tier = factor(tier, levels = TORD))

# One genus factor and one type ordering for both panels: the two tables
# partition the same cohort (a library appears in exactly one), so they are only
# worth showing side by side if the reader can track a row across them.
# Only genera actually present, so the shared legend carries no empty swatches
# and both panels declare an identical scale, which is what lets patchwork
# collect them into one guide.
GLEV <- intersect(c("alpha", "beta", "gamma", "mu", "nu", "other"),
                  union(s1$genus, s6$genus))
TYPE_ORDER <- bind_rows(s1, s6) %>% group_by(type, genus) %>%
  summarise(tot = sum(n), .groups = "drop") %>%
  mutate(genus = factor(genus, levels = GLEV)) %>%
  arrange(genus, desc(tot)) %>% pull(type)

bubble <- function(d, ttl) {
  d <- d %>% mutate(type = factor(type, levels = rev(TYPE_ORDER)),
                    genus = factor(genus, levels = GLEV))
  # An invisible point per genus so both panels declare an identical colour
  # guide; without it the panel lacking a genus renders a different key and
  # patchwork cannot collect the two into one legend.
  ghost <- tibble(tier = factor(TORD[1], levels = TORD),
                  type = factor(TYPE_ORDER[1], levels = rev(TYPE_ORDER)),
                  n = 1, genus = factor(GLEV, levels = GLEV))
  ggplot(d, aes(tier, type, size = n, colour = genus)) +
    geom_point(data = ghost, alpha = 0) +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = n), size = 2.1, colour = "white", show.legend = FALSE) +
    scale_size_area(max_size = 8, guide = "none") +
    scale_colour_manual(values = GEN, name = "Genus", limits = GLEV, drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(labels = TLAB, drop = FALSE) +
    labs(x = NULL, y = NULL, title = ttl) +
    theme(axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(size = 8),
          panel.grid.major = element_line(colour = "grey92"),
          axis.line = element_blank(), axis.ticks = element_blank())
}

pS1 <- bubble(s1, "Table S1: assigned types (breadth ≥ 0.10, depth ≥ 2)")
pS6 <- bubble(s6, "Table S6: detection tier (below the assignment threshold)")
ggsave(file.path(OUT, "prov_S1_assigned_types.png"), pS1, width = 8, height = 8, dpi = 200, bg = "white")
ggsave(file.path(OUT, "prov_S6_trace_types.png"),   pS6, width = 8, height = 8, dpi = 200, bg = "white")
ggsave(file.path(OUT, "prov_S1S6_combined.png"),
       # patchwork will not collect the two colour guides (the panels' size
       # scales differ), so the right panel simply drops its legend.
       (pS1 + (pS6 + theme(axis.text.y = element_blank()) + guides(colour = "none"))) +
         plot_layout(guides = "collect") &
         theme(legend.position = "bottom") &
         guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))),
       width = 13, height = 9.5, dpi = 200, bg = "white")

## ── S1 + S6 as a faceted tile grid ──────────────────────────────────────────
## Structure follows the version drafted by hand: genus as facet rows, assigned
## against detected as facet columns, so the partition between the two tables is
## explicit rather than something the reader has to reconstruct.
##
## The colour problem in that draft is that fill was linear on a skewed count:
## 37% of cells carry n <= 3, and on a linear 1-35 ramp n = 3 sits at 6%, which
## is effectively white, so the white count labels vanished. Three changes:
##   1. log10 fill, which moves n = 3 to 31% of the ramp
##   2. a palette that does not start at white (begin = 0.12)
##   3. label colour chosen per tile from the fill, rather than always white
tile_labels <- function(n, lim) {
  # White on the dark half, near-black on the light half. Computed on the same
  # log scale as the fill so the switch happens where the palette actually
  # turns dark, not where the raw count happens to sit.
  pos <- (log10(n) - log10(lim[1])) / (log10(lim[2]) - log10(lim[1]))
  ifelse(is.na(pos) | pos < 0.55, "grey15", "white")
}

tiles <- function(d, opt = "mako", begin = 0.12, end = 0.95) {
  lim <- c(1, max(d$n, na.rm = TRUE))
  d <- d %>% mutate(type = factor(type, levels = rev(TYPE_ORDER)),
                    genus = factor(genus, levels = GLEV),
                    lab_col = tile_labels(n, lim))
  ggplot(d, aes(tier, type, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(n), "", n), colour = lab_col),
              size = 2.7, show.legend = FALSE) +
    facet_grid(genus ~ id, space = "free", scales = "free_y") +
    scale_fill_viridis_c("Libraries", option = opt, direction = -1,
                         trans = "log10", limits = lim,
                         begin = begin, end = end,
                         breaks = c(1, 3, 10, 30), na.value = "grey97") +
    scale_colour_identity() +
    scale_x_discrete(labels = function(x) gsub("\\n", " ", TLAB[x])) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7),
          axis.title = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text.y = element_text(angle = 0),
          panel.spacing.y = unit(2, "pt"),
          legend.key.height = unit(1.1, "cm"))
}

grid_d <- bind_rows(assigned = rename(s1, ref = hpv_reference),
                    detected = rename(s6, ref = top_ref), .id = "id") %>%
  select(id, tier, type, genus, n) %>%
  mutate(id = factor(id, levels = c("assigned", "detected"),
                     labels = c("Assigned (S1)", "Detection tier (S6)"))) %>%
  complete(id, nesting(tier, genus, type))

for (o in c("mako", "viridis", "rocket")) {
  ggsave(file.path(OUT, sprintf("prov_S1S6_tiles_%s.png", o)),
         tiles(grid_d, opt = o), width = 9, height = 12.5, dpi = 200, bg = "white")
}

## ── Table S7: co-infected libraries ─────────────────────────────────────────
## Every assigned type in every library carrying more than one, ordered by the
## deepest type. The three cautions the table's caption spells out are all
## visual here: how lopsided the pairs are, how close the secondary call sits
## to the assignment threshold, and how much of the set is one donor group.
dirs <- c("results_full_v2", "results_lesions", "results_lesions_csc", "results_lesions_wide")
ty <- map_dfr(unlist(lapply(dirs, function(d)
        list.files(here(d, "aggregated"), "_hpv_types\\.tsv$", full.names = TRUE))),
      ~ suppressMessages(read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))) %>%
  distinct(sample_id, hpv_reference, .keep_all = TRUE) %>%
  mutate(coverage_breadth = as.double(coverage_breadth), mean_depth = as.double(mean_depth),
         type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  filter(coverage_breadth >= 0.10, mean_depth >= 2)

co <- ty %>% add_count(sample_id, name = "n_types") %>% filter(n_types >= 2) %>%
  left_join(tiers, by = c("sample_id" = "srr_id")) %>%
  left_join(genus, by = c("hpv_reference" = "ref")) %>%
  mutate(genus = replace_na(genus, "other"))

ordr <- co %>% group_by(sample_id) %>%
  summarise(tier = first(tier), dom = max(mean_depth), .groups = "drop") %>%
  mutate(tier = factor(tier, levels = TORD)) %>% arrange(tier, dom)
co <- co %>% mutate(sample_id = factor(sample_id, levels = ordr$sample_id),
                    tier = factor(tier, levels = TORD),
                    pair_study = study %in% c("SRP314958", "SRP314961"))

pS7 <- ggplot(co, aes(mean_depth, sample_id)) +
  geom_line(aes(group = sample_id), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = genus, size = coverage_breadth), alpha = 0.9) +
  geom_text_repel(aes(label = type), size = 2, colour = "grey20",
                  direction = "x", nudge_x = 0.18, segment.size = 0.2,
                  segment.colour = "grey80", min.segment.length = 0.15,
                  box.padding = 0.12, max.overlaps = Inf, seed = 1) +
  geom_vline(xintercept = 2, linetype = 2, colour = "red") +
  scale_x_log10(labels = scales::comma, limits = c(1, 2e5),
                expand = expansion(mult = c(0.04, 0.10))) +
  scale_colour_manual(values = GEN, name = "Genus") +
  scale_size_area(max_size = 4, name = "Coverage\nbreadth") +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(tier = function(x) gsub("\n", " ", TLAB[x]))) +
  labs(x = "Mean depth of each assigned type (log scale)", y = NULL,
       title = "Table S7: libraries carrying more than one assigned HPV type",
       subtitle = "One row per library, a line joining the types called in it; red = the 2x assignment threshold") +
  theme(axis.text.y = element_text(size = 5.5), strip.text.y = element_text(angle = 0, size = 7),
        panel.grid.major.x = element_line(colour = "grey92"))
ggsave(file.path(OUT, "prov_S7_coinfection.png"), pS7, width = 9.5, height = 10, dpi = 200, bg = "white")

cat(sprintf("S1: %d type-tier cells, %d types\n", nrow(s1), n_distinct(s1$type)))
cat(sprintf("S6: %d type-tier cells, %d types\n", nrow(s6), n_distinct(s6$type)))
cat(sprintf("S7: %d libraries, %d type calls; %d libraries from the condyloma study pair\n",
            n_distinct(co$sample_id), nrow(co),
            n_distinct(co$sample_id[co$pair_study])))
cat("wrote 4 files to manuscript/provisional/\n")
