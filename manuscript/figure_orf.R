#!/usr/bin/env Rscript
# Builds manuscript/figure3.{png,pdf}: the viral transcriptional programme by tier.
#
#   Rscript manuscript/figure_orf.R
#
# WHY THIS FIGURE EXISTS. The paper's central criterion is the early-versus-late
# transcript programme: "productive" means late capsid transcripts are present.
# Every claim in the paper rests on that plane and the reader has never been
# shown it. Three things become visible at once:
#
#   A  programme breadth per ORF. Productive tiers run the complete genome,
#      capsid genes included; skin and neoplasia are uniformly low across all
#      of it (a trace, not a truncated programme); cell lines and engineered
#      material show the integrated signature directly, oncogenes on with E2
#      and the late genes collapsed.
#   B  programme composition. The same contrast in a different currency, and
#      the panel that speaks to cis-restriction of beta-HPV expression, which
#      is a claim about programme shape rather than about detection.
#   C  the criterion itself, early against late reads per library with the
#      productive threshold drawn.
#
# CONDITIONALITY. Every panel is conditional on HPV positivity and successful
# typing: these are the libraries for which a transcript classification exists.
# This is a figure about the programme of detected infections, never about
# prevalence. Tier n is printed on panel A for that reason.
suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(here); library(patchwork)
})
old <- theme_set(theme_cowplot(font_size = 11))
source(here("manuscript/_paths.R"))

DIRS <- c("results_full_v2", "results_lesions", "results_lesions_csc", "results_lesions_wide")
# E5 is dropped: beta and gamma types largely lack it and it is zero in every
# tier, so showing it adds an empty column and invites the wrong question.
ORFS <- c("E1", "E2", "E6", "E7", "L1", "L2")
LATE <- c("L1", "L2")
PRODUCTIVE_L1 <- 3          # the manuscript's productive criterion

TLAB <- c(control_productive = "Productive\ncontrols", EV = "EV",
          peri_lesional = "Peri-lesional", neoplasia = "Neoplasia",
          unselected_skin = "Unselected\nskin", cell_line = "Cell line",
          engineered = "Engineered")

fs <- unlist(lapply(DIRS, function(d)
  list.files(here(d, "aggregated"), "_transcript_classes[.]tsv$", full.names = TRUE)))
stopifnot(length(fs) > 0)
tc <- map_dfr(fs, ~ suppressMessages(
  read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c"))))

tiers <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% select(srr_id, tier)
# 53 libraries are processed in more than one results directory (a lesion cohort
# re-running an accession already in the null sweep). Every duplicated
# (library, gene) pair carries an identical count, so de-duplication is lossless;
# without it those libraries are weighted twice and the composition panel does
# not sum to one.
d <- tc %>%
  filter(class %in% c("early", "late"), gene %in% ORFS) %>%
  transmute(srr_id = sample_id, gene = factor(gene, levels = ORFS),
            n = as.numeric(read_count)) %>%
  distinct(srr_id, gene, .keep_all = TRUE) %>%
  inner_join(tiers, by = "srr_id") %>%
  filter(tier %in% names(TLAB)) %>%
  mutate(tier = factor(tier, levels = names(TLAB)))
stopifnot(nrow(d) == n_distinct(d$srr_id) * length(ORFS))

ntier <- d %>% distinct(srr_id, tier) %>% count(tier, name = "n_lib")
lab_n <- setNames(sprintf("%s\n(n = %d)", TLAB[as.character(ntier$tier)], ntier$n_lib),
                  as.character(ntier$tier))

## ── Panel A: programme breadth ──────────────────────────────────────────────
## Beta-Binomial posterior mean under Beta(0.5, 0.5), matching the type tables.
brd <- d %>% group_by(tier, gene) %>%
  summarise(k = sum(n > 0), n = n(), .groups = "drop") %>%
  mutate(pct = 100 * (k + 0.5) / (n + 1))

pA <- ggplot(brd, aes(gene, fct_rev(tier), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = round(pct)), size = 3,
            colour = ifelse(brd$pct > 55, "white", "grey20")) +
  scale_fill_gradient(low = "#f7f7f7", high = "#08519c", limits = c(0, 100),
                      name = "% of libraries\nwith ≥1 read") +
  scale_y_discrete(labels = lab_n) +
  labs(x = NULL, y = NULL, subtitle = "A  Programme breadth: which ORFs are transcribed at all") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.key.width = unit(0.35, "cm"))

## ── Panel B: programme composition ──────────────────────────────────────────
## Composition is only estimable where there are reads to apportion. Unselected
## skin has a median of ZERO viral reads across these six ORFs, so on a
## per-library mean-of-fractions every single read would swing a whole tier: a
## library with one L2 read would read as 100% late transcription. Restricting
## to libraries with at least MIN_COMP reads keeps the panel from appearing to
## show abundant late transcription in tiers the rest of the paper reports as
## quiet. The restriction is stated on the panel, not hidden.
MIN_COMP <- 20
comp_lib <- d %>% group_by(srr_id, tier) %>% mutate(tot = sum(n)) %>% ungroup() %>%
  filter(tot >= MIN_COMP)
comp <- comp_lib %>% mutate(frac = n / tot) %>%
  group_by(tier, gene) %>% summarise(share = mean(frac), .groups = "drop")
ncomp <- comp_lib %>% distinct(srr_id, tier) %>% count(tier, name = "n_lib")
lab_c <- setNames(sprintf("%s\n(n = %d)", TLAB[as.character(ncomp$tier)], ncomp$n_lib),
                  as.character(ncomp$tier))

pB <- ggplot(comp, aes(fct_rev(tier), share, fill = gene)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_discrete(labels = lab_c) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Mean share of a library's viral reads",
       subtitle = sprintf("B  Programme composition: where the viral reads fall (libraries with \u2265 %d viral reads)",
                          MIN_COMP)) +
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

## ── Panel C: the criterion ──────────────────────────────────────────────────
# The productive criterion is L1 >= 3 specifically, so the y axis is L1 rather
# than the whole late region: the dashed line is then exactly the rule the
# paper applies, not an approximation of it.
ec <- d %>%
  mutate(part = case_when(gene == "L1" ~ "L1",
                          gene %in% LATE ~ "L2", TRUE ~ "early")) %>%
  group_by(srr_id, tier, part) %>% summarise(n = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = part, values_from = n, values_fill = 0)

pC <- ggplot(ec, aes(early + 1, L1 + 1, colour = tier)) +
  geom_hline(yintercept = PRODUCTIVE_L1 + 1, linetype = 2, colour = "grey35") +
  annotate("text", x = Inf, y = PRODUCTIVE_L1 + 1, hjust = 1.05, vjust = -0.6,
           label = "productive threshold (L1 \u2265 3)", size = 3, colour = "grey35") +
  geom_point(alpha = 0.65, size = 1.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_brewer(palette = "Dark2", labels = TLAB, name = NULL) +
  labs(x = "Early-region reads + 1", y = "L1 major-capsid reads + 1",
       subtitle = "C  The criterion: early-region against L1 capsid transcription, per library") +
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

fig <- pA / pB / pC + plot_layout(heights = c(1.15, 0.85, 1.2))
ggsave(here("manuscript/figure3.png"), fig, width = 9, height = 11.5, dpi = 300, bg = "white")
ggsave(here("manuscript/figure3.pdf"), fig, width = 9, height = 11.5)

write_tsv(brd  %>% mutate(pct = round(pct, 1)), here("manuscript/orf_breadth_by_tier.tsv"))
write_tsv(comp %>% mutate(share = round(share, 4)) %>% left_join(ncomp, by = "tier"),
          here("manuscript/orf_composition_by_tier.tsv"))
cat(sprintf("figure3 written; %d libraries across %d tiers\n",
            n_distinct(d$srr_id), n_distinct(d$tier)))
print(as.data.frame(ntier))
