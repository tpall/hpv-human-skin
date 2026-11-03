#!/usr/bin/env Rscript
# Figure S6: how the two analysis populations are built, and where they differ.
#
# TP comment 82 on draft n: "maybe a box chart can visualize better how we came
# up with different numbers." Five counts are in play and the Methods carried all
# of them in prose. Drawn as a flow, the structure is simple and the thing a
# reader most needs to see becomes visible: the positivity denominator
# (@tbl-positivity) and the depth-normalized burden cohort are not nested and
# neither is a filtered version of the other. They share one stage, the screened
# libraries, and then diverge, because the burden analysis adds the targeted
# lesion cohorts and removes two sets of libraries that neither a null nor a
# control describes.
#
# Every count is computed here from the same sources the analysis scripts read,
# never typed in, so the figure cannot drift from the text.
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("manuscript/_paths.R"))

n_reports <- function(dir) {
  f <- list.files(here(dir, "aggregated"), "_kraken2_report\\.txt$", full.names = TRUE)
  c(submitted = length(f), screened = sum(file.size(f) > 0))
}
sweep  <- n_reports("results_full_v2")
les    <- n_reports("results_lesions_wide") + n_reports("results_lesions_csc")
burden <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% count(cohort)
b_null <- burden$n[burden$cohort == "null"]
b_les  <- burden$n[burden$cohort == "lesion"]

## the two causes of an empty report, told apart by whether fastp ran at all
fastp <- sub("_fastp\\.json$", "", basename(
  list.files(here("results_full_v2/chunks"), "_fastp\\.json$", recursive = TRUE)))
unscr <- unscreened_ids()
n_filt <- sum(unscr %in% fastp); n_nodl <- length(unscr) - n_filt

cm <- function(x) formatC(x, format = "d", big.mark = ",")
box <- function(x, y, label, kind) tibble(x, y, label, kind)

## main column stages, then the losses that sit between them
stages <- bind_rows(
  box(1, 5, sprintf("Untargeted prevalence sweep\n%s libraries submitted to the screen", cm(sweep["submitted"])), "stage"),
  box(4, 5, sprintf("Targeted lesion cohorts\n%s libraries submitted to the screen", cm(les["submitted"])), "stage"),
  box(1, 3, sprintf("%s screened\ndenominator of Table 1", cm(sweep["screened"])), "keep"),
  box(4, 3, sprintf("%s screened", cm(les["screened"])), "keep"),
  box(1, 1, cm(b_null), "keep"),
  box(4, 1, cm(b_les), "keep"),
  box(2.5, -0.42, sprintf("Depth-normalized burden cohort\n%s libraries", cm(b_null + b_les)), "stage"))

losses <- bind_rows(
  box(2.0, 4, sprintf("%s returned no screening\nresult: %s lost every read to\nthe %s nt length filter,\n%s never downloaded",
                      cm(length(unscr)), cm(n_filt), 36, cm(n_nodl)), "loss"),
  box(5.15, 4, sprintf("%s returned no\nscreening result", cm(les["submitted"] - les["screened"])), "loss"),
  box(2.0, 2, sprintf("%s 16S rRNA amplicon\nlibraries excluded:\nno detectable HPV at all", cm(sweep["screened"] - b_null)), "loss"),
  box(5.15, 2, sprintf("%s condyloma-adjacent\nmucosa excluded", cm(les["screened"] - b_les)), "loss"))

# the two columns converge on the burden cohort box rather than meeting in a
# corner, so the join reads as one step and not as a stray horizontal arrow
arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  1, 4.55, 1, 3.45,    4, 4.55, 4, 3.45,
  1, 2.55, 1, 1.45,    4, 2.55, 4, 1.45,
  1, 0.66, 2.32, 0.02, 4, 0.66, 2.68, 0.02)
side <- tribble(~x, ~y, ~xend, ~yend,
  1, 4.0, 1.55, 4.0,   4, 4.0, 4.55, 4.0,
  1, 2.0, 1.55, 2.0,   4, 2.0, 4.55, 2.0)

p <- ggplot() +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(5, "pt"), type = "closed"), linewidth = 0.4) +
  geom_segment(data = side, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(4, "pt"), type = "closed"),
               linewidth = 0.3, linetype = "22") +
  geom_label(data = stages, aes(x, y, label = label, fill = kind),
             size = 2.9, lineheight = 1.15, label.r = unit(3, "pt"),
             label.padding = unit(4.5, "pt"), colour = "black") +
  geom_label(data = losses, aes(x, y, label = label, fill = kind),
             size = 2.5, lineheight = 1.15, label.r = unit(3, "pt"),
             label.padding = unit(4, "pt"), colour = "grey25", hjust = 0) +
  scale_fill_manual(values = c(stage = "grey88", keep = "white", loss = "grey96"),
                    guide = "none") +
  # geom_label sizes are absolute, so the panel needs room the data range
  # does not imply: without this the leftmost boxes are clipped
  scale_x_continuous(expand = expansion(mult = c(0.36, 0.42))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.06))) +
  theme_void()

ggsave(here("manuscript/figureS6.png"), p, width = 6.5, height = 4.4, dpi = 300, bg = "white")
ggsave(here("manuscript/figureS6.pdf"), p, width = 6.5, height = 4.4)
cat(sprintf("figureS6: sweep %s->%s, lesion %s->%s, burden %s+%s=%s\n",
            cm(sweep["submitted"]), cm(sweep["screened"]), cm(les["submitted"]),
            cm(les["screened"]), cm(b_null), cm(b_les), cm(b_null + b_les)))
