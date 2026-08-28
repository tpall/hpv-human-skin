#!/usr/bin/env Rscript
# Builds manuscript/figureS5.{png,pdf}: the Bayesian re-analysis as a forest plot.
#
#   Rscript manuscript/figure_bayes.R
#
# The posterior summaries were reported only as a table (@supptbl-bayes) and as
# inline text. A table makes the reader compare nine intervals by reading digits;
# on a log axis the same nine intervals show their sizes, their overlap with the
# no-effect line, and the two that straddle it, at a glance. The table stays, for
# the exact numbers and the frequentist comparison.
#
# Rows and their order mirror build_bayes_table() in _helpers.R exactly, so the
# figure and the table cannot present the contrasts in different orders. The
# Student-t sensitivity fit, which the table carries only in its caption, is
# drawn here directly beneath the contrast it tests.
suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(here)
})
old <- theme_set(theme_cowplot(font_size = 11))

b <- read_tsv(here("manuscript/bayes_summary.tsv"), show_col_types = FALSE)

# Labels are the contrast strings build_bayes_table() prints, verbatim, so a row
# here and a row of @supptbl-bayes are recognisably the same contrast. The model
# tag rides in the label rather than in a facet strip: faceting by model looked
# tidier until the one-row M3 panel squeezed its own strip text to "producti".
rows <- tribble(
  ~model,                                ~label,
  "M1_cell_line",                        "Cell-line vs clinical positivity (M1)",
  "M1_engineered",                       "Engineered vs clinical positivity (M1)",
  "M2c_controls_vs_skin_floor",          "Productive controls vs unselected skin, floor (M2)",
  "M2c_controls_vs_skin_floor_studentt", "     the same, under a Student-t(3) prior",
  "M2b_neoplasia_vs_skin_floor",         "Neoplasia vs unselected skin, floor (M2)",
  "M2b_depth",                           "Sequencing depth, per 10× reads, floor (M2)",
  "M2d_depth_within_neoplasia",          "Sequencing depth within neoplasia, floor (M2)",
  "M2a_neoplasia_vs_skin_any",           "Neoplasia vs unselected skin, any-HPV / trace (M2)",
  "M3_controls_vs_neoplasia",            "Productive controls vs neoplasia, productive (M3)")

d <- rows %>% left_join(b, by = "model") %>%
  mutate(label  = fct_rev(fct_inorder(label)),
         # The two depth terms are the ones whose interval spans no effect, and
         # that is a result rather than a shortfall: deeper sequencing does not
         # buy productive-level detection. Draw them open so the eye finds them.
         decisive = if_else(lo > 1 | hi < 1, "credible interval excludes 1",
                                             "credible interval includes 1"))
stopifnot(!any(is.na(d$OR)))

# Same rounding as pdir() in _helpers.R: a posterior probability is never printed
# as 1.000 here when the text prints >0.999 for the same quantity.
d <- d %>% mutate(p_lab = if_else(p_gt1 >= 0.9995, ">0.999",
                                  formatC(p_gt1, format = "f", digits = 3)))

brk <- 10^(-1:5)
xmax <- max(d$hi)

# The P(OR > 1) column header, drawn once above the first row of the top facet.
hdr <- d %>% slice(1) %>% transmute(label, x = xmax * 12, txt = "P(OR > 1)")

p <- ggplot(d, aes(OR, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35") +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = decisive), linewidth = 0.6) +
  geom_point(aes(colour = decisive, shape = decisive), size = 2.4, fill = "white", stroke = 0.9) +
  geom_text(aes(x = xmax * 12, label = p_lab), size = 2.7, colour = "grey25", hjust = 0) +
  # One header, not one per facet. Annotating at a fixed y forced every panel's
  # range up to that value, which is what left the rows stranded at the bottom of
  # each block with dead space above them.
  geom_text(data = hdr, aes(x = x, y = label, label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = -1.4, size = 2.7, colour = "grey25", fontface = "bold") +
  scale_x_log10(breaks = brk,
                labels = function(x) formatC(x, format = "fg", big.mark = ",", drop0trailing = TRUE),
                expand = expansion(mult = c(0.04, 0.42))) +
  scale_colour_manual(values = c("credible interval excludes 1" = "#0072B2",
                                 "credible interval includes 1" = "grey55"), name = NULL) +
  scale_shape_manual(values = c("credible interval excludes 1" = 19,
                                "credible interval includes 1" = 21), name = NULL) +
  coord_cartesian(clip = "off") +
  labs(x = "Posterior odds ratio (log scale), with 95% credible interval", y = NULL) +
  theme(legend.position = "bottom",
        plot.margin = margin(t = 20, r = 8, b = 6, l = 6))

ggsave(here("manuscript/figureS5.png"), p, width = 9.0, height = 4.4, dpi = 200, bg = "white")
ggsave(here("manuscript/figureS5.pdf"), p, width = 9.0, height = 4.4, bg = "white")
cat(sprintf("Wrote manuscript/figureS5.{png,pdf}: %d contrasts, %d with a CrI excluding 1\n",
            nrow(d), sum(d$decisive == "credible interval excludes 1")))
