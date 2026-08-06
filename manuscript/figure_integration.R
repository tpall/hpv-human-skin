#!/usr/bin/env Rscript
# Builds manuscript/figure2.{png,pdf}: the integration read-out, three panels.
#
#   Rscript manuscript/figure_integration.R
#
# WHY THIS FIGURE EXISTS. The integration result was carried entirely by three
# dense paragraphs of numbers in the Results, which is the hardest passage in
# the manuscript to follow. All three claims are visual ones:
#
#   A  the E1:E2 transcript ratio separates integrated from episomal material
#   B  junctions confirm it, but only when scored as CLUSTERING, not as counts
#   C  the ratio moves with lesion grade in primary tissue, so it reads
#      integration as biology rather than culture as provenance
#
# Panel B is the one that most needs a picture: the claim "raw counts are
# actively misleading" is nearly impossible to state in prose without the reader
# having to hold two numbers per library in their head.
suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(here); library(patchwork)
})
old <- theme_set(theme_cowplot(font_size = 11))

E1E2   <- here("results_e1e2/e1e2_by_library.tsv")
JX     <- here("results_junction/junction_summary.tsv")
CXG    <- here("results_cervix_graded_e1e2/e1e2_by_library.tsv")
JXCUT  <- 0.2      # threshold inside the empty interval, see Methods

stopifnot(file.exists(E1E2), file.exists(JX), file.exists(CXG))

## ── Panel A: E1:E2 ratio by tier ────────────────────────────────────────────
TLAB <- c(cell_line = "Cell line", engineered = "Engineered", other = "Site not\nrecorded",
          unselected_skin = "Unselected\nskin", neoplasia = "Neoplasia",
          control_productive = "Productive\ncontrols", EV = "EV",
          peri_lesional = "Peri-lesional")
e <- read_tsv(E1E2, show_col_types = FALSE) %>%
  filter(tier %in% names(TLAB), !is.na(log2_e1e2)) %>%
  mutate(tier = factor(tier, levels = names(TLAB)),
         grp  = if_else(tier %in% c("cell_line", "engineered"),
                        "in vitro", "primary tissue"))

pA <- ggplot(e, aes(tier, log2_e1e2, colour = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey35", fill = NA) +
  geom_jitter(width = 0.18, height = 0, size = 1.1, alpha = 0.55) +
  scale_x_discrete(labels = TLAB) +
  scale_colour_manual(values = c(`in vitro` = "#c0392b",
                                 `primary tissue` = "#2c6fbb"), name = NULL) +
  labs(x = NULL, y = expression(log[2]~"(E1 / E2)"),
       title = "A  Integration inverts the E1:E2 transcript ratio") +
  theme(legend.position = c(0.72, 0.92),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        axis.text.x = element_text(size = 8.5))

## ── Panel B: junctions, counts vs clustering ────────────────────────────────
GLAB <- c(cell_line = "Cell line", neoplasia = "Neoplasia",
          control_productive = "Productive\ncontrols", EV = "EV")
j <- read_tsv(JX, show_col_types = FALSE) %>%
  mutate(state = if_else(frac_top >= JXCUT, "integrated", "episomal"),
         tier  = factor(tier, levels = names(GLAB)))

pB <- ggplot(j, aes(n_junctions, frac_top)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = JXCUT, ymax = Inf,
           fill = "#c0392b", alpha = 0.06) +
  geom_hline(yintercept = JXCUT, linetype = "dashed", colour = "#c0392b") +
  geom_point(aes(shape = tier, colour = state), size = 2.6) +
  scale_x_log10(labels = scales::comma) +
  scale_colour_manual(values = c(integrated = "#c0392b", episomal = "#2c6fbb"),
                      name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 3), labels = GLAB, name = NULL) +
  labs(x = "Virus-host junctions detected (n, log scale)",
       y = "Fraction at top site",
       title = "B  Clustering separates, counts do not") +
  theme(legend.position = "right", legend.key.height = unit(0.8, "lines"),
        legend.text = element_text(size = 8))

## ── Panel C: cervical progression ───────────────────────────────────────────
CLAB <- c(hpv16_normal_epithelium = "Normal\nepithelium",
          hpv16_hsil = "HSIL", hpv16_scc = "Carcinoma")
cx <- read_tsv(CXG, show_col_types = FALSE) %>%
  filter(tier %in% names(CLAB), !is.na(log2_e1e2)) %>%
  mutate(tier = factor(tier, levels = names(CLAB)))
cx_med <- cx %>% group_by(tier) %>% summarise(m = median(log2_e1e2), .groups = "drop")

pC <- ggplot(cx, aes(tier, log2_e1e2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(data = cx_med, aes(tier, m, group = 1), colour = "#c0392b", linewidth = 0.8) +
  geom_point(data = cx_med, aes(tier, m), colour = "#c0392b", size = 3) +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.6, colour = "grey30") +
  scale_x_discrete(labels = CLAB) +
  labs(x = NULL, y = expression(log[2]~"(E1 / E2)"),
       title = "C  Ratio tracks lesion grade") +
  theme(axis.text.x = element_text(size = 9))

fig <- (pA / (pB | pC)) + plot_layout(heights = c(1, 0.95))

ggsave(here("manuscript/figure2.png"), fig, width = 9.5, height = 7.6, dpi = 300, bg = "white")
ggsave(here("manuscript/figure2.pdf"), fig, width = 9.5, height = 7.6, bg = "white")

## ── numbers the caption quotes, printed so they are never hand-typed ────────
cat("\nPanel A, median fold E1:E2 by tier:\n")
e %>% group_by(tier) %>%
  summarise(n = n(), median_log2 = median(log2_e1e2), fold = 2^median(log2_e1e2)) %>%
  as.data.frame() %>% print(right = FALSE)
cat("\nPanel B, frac_top range by called state:\n")
j %>% group_by(state) %>%
  summarise(n = n(), min = min(frac_top), max = max(frac_top),
            junctions_min = min(n_junctions), junctions_max = max(n_junctions)) %>%
  as.data.frame() %>% print(right = FALSE)
cat("\n  the two EV libraries (episomal) carry the LARGEST junction counts:\n")
j %>% filter(tier == "EV") %>% select(srr_id, n_junctions, frac_top, state) %>%
  as.data.frame() %>% print(right = FALSE)
cat("\nPanel C, cervical progression:\n")
cx %>% group_by(tier) %>%
  summarise(n = n(), median_log2 = median(log2_e1e2), fold = 2^median(log2_e1e2)) %>%
  as.data.frame() %>% print(right = FALSE)
cat("\nWrote manuscript/figure2.{png,pdf}\n")
