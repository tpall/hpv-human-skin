## Supplementary figure: what the sub-threshold trace is made of.
##
## Splits the below-floor HPV-positive libraries of unselected skin and neoplasia
## into three exclusive classes and sets them against the productive controls.
## `other` (the residual sample group with no recorded site) is deliberately NOT
## in the denominator: the manuscript claim is about skin and neoplasia.
##
## Reads trace_decomposition/trace_rows.json (decompose_trace.py) and writes the
## figure, the composition table, and the numbers the manuscript quotes inline,
## so figure and text cannot drift apart.
suppressPackageStartupMessages({library(tidyverse); library(jsonlite); library(cowplot); library(here)})
theme_set(theme_cowplot(font_size = 12))

j     <- fromJSON(here("trace_decomposition/trace_rows.json"))
FLOOR <- j$floor
TRACE <- "Below-floor trace\n(skin and neoplasia)"
CTRL  <- "Productive controls\n(warts and EV)"

d <- as_tibble(j$rows) %>%
  mutate(group = case_when(
           tier %in% c("control_productive", "EV")                       ~ CTRL,
           tier %in% c("unselected_skin", "neoplasia") & rpm < FLOOR     ~ TRACE,
           TRUE                                                          ~ NA_character_)) %>%
  filter(!is.na(group)) %>%
  # Exclusive and ordered: the productive call wins, then anything that put a read
  # on a gene at all, and what is left carries no gene-mapped reads whatsoever.
  # The middle class is overwhelmingly early-gene-without-late (the latency
  # signature), but a handful carry only sub-threshold L2, so the label has to be
  # the weaker one; the split is reported in trace_numbers.tsv and in the text.
  mutate(cat = case_when(L1 >= 3                 ~ "Productive (L1 ≥ 3 reads)",
                         EARLY > 0 | LATE > 0    ~ "Gene-mapped, not productive",
                         TRUE                    ~ "No gene-mapped reads"),
         cat   = factor(cat, levels = c("No gene-mapped reads",
                                        "Gene-mapped, not productive",
                                        "Productive (L1 ≥ 3 reads)")),
         group = factor(group, levels = c(TRACE, CTRL)))

# The "no gene-mapped reads" label is only honest if those libraries really carry
# none, late region included. Fail loudly rather than mislabel the grey block.
stopifnot(all(d$LATE[d$cat == "No gene-mapped reads"] == 0),
          all(d$EARLY[d$cat == "No gene-mapped reads"] == 0))

comp <- d %>% count(group, cat) %>% group_by(group) %>%
  mutate(frac = n / sum(n), N = sum(n)) %>% ungroup()
write_csv(comp, here("trace_decomposition/trace_composition.csv"))

cols <- c("No gene-mapped reads"        = "grey70",
          "Gene-mapped, not productive" = "#D95F02",
          "Productive (L1 ≥ 3 reads)"   = "#1B9E77")

# No in-plot title or subtitle: what the panel shows belongs in the figure
# legend, matching every other figure in this manuscript.
p <- ggplot(comp, aes(group, frac, fill = cat)) +
  geom_col(width = .68, colour = "white", linewidth = .3) +
  geom_text(aes(label = ifelse(frac > .04, scales::percent(frac, 1), "")),
            position = position_stack(vjust = .5), size = 3.6) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .02))) +
  labs(x = NULL, y = "Percentage of HPV-positive libraries") +
  theme(legend.position = "right")
# Written twice on purpose: manuscript/figureS2.* is the supplementary figure,
# trace_decomposition/figure_trace_decomposition.* is what the slide deck points at.
# cairo_pdf, not the default device: the legend carries a U+2265 that the base
# pdf() device silently downgrades to ">=".
for (dest in c("trace_decomposition/figure_trace_decomposition",
               "manuscript/figureS2")) {
  ggsave(here(paste0(dest, ".png")), p, width = 9, height = 5, dpi = 200, bg = "white")
  ggsave(here(paste0(dest, ".pdf")), p, width = 9, height = 5, bg = "white", device = cairo_pdf)
}

## ── numbers quoted in the Results sentences and the caption ─────────────────
pc  <- function(g, k) round(100 * comp$frac[comp$group == g & comp$cat == k], 0)
nn  <- function(g, k) comp$n[comp$group == g & comp$cat == k]
med <- function(g) round(median(d$latefrac[d$group == g]), 3)
# Can the early-only signal be integrated-line contamination? That is mucosal
# HPV16/18, so ask what the E6/E7-transcribing trace libraries are typed as.
# Most are untypeable at this coverage, which bounds how far this can be pushed.
tr  <- d %>% filter(group == TRACE, e67 > 0)
# Within the middle class, how many are early-without-late (the interesting shape)
# versus late fragments below the productive bar.
mid <- d %>% filter(group == TRACE, cat == "Gene-mapped, not productive")
tibble(key = c("floor", "trace_n", "trace_nogene", "trace_mid", "trace_prod",
               "trace_pct_nogene", "trace_pct_mid", "trace_pct_prod", "trace_medlate",
               "trace_mid_early", "trace_mid_lateonly",
               "ctrl_n", "ctrl_nogene", "ctrl_mid", "ctrl_prod", "ctrl_pct_prod", "ctrl_medlate",
               "trace_e67_n", "trace_e67_typed", "trace_e67_hela"),
       value = c(round(FLOOR, 1),
                 comp$N[comp$group == TRACE][1], nn(TRACE, "No gene-mapped reads"),
                 nn(TRACE, "Gene-mapped, not productive"), nn(TRACE, "Productive (L1 ≥ 3 reads)"),
                 pc(TRACE, "No gene-mapped reads"), pc(TRACE, "Gene-mapped, not productive"),
                 pc(TRACE, "Productive (L1 ≥ 3 reads)"), med(TRACE),
                 sum(mid$EARLY > 0), sum(mid$EARLY == 0),
                 comp$N[comp$group == CTRL][1], nn(CTRL, "No gene-mapped reads"),
                 nn(CTRL, "Gene-mapped, not productive"), nn(CTRL, "Productive (L1 ≥ 3 reads)"),
                 pc(CTRL, "Productive (L1 ≥ 3 reads)"), med(CTRL),
                 nrow(tr), sum(!is.na(tr$type)), sum(tr$type %in% c("16", "18")))) %>%
  write_tsv(here("trace_decomposition/trace_numbers.tsv"))

cat("wrote figure + trace_composition.csv + trace_numbers.tsv\n")
print(comp)
cat(sprintf("\nE6/E7+ trace libraries: %d, of which typed %d, of which HPV16/18 %d\n",
            nrow(tr), sum(!is.na(tr$type)), sum(tr$type %in% c("16", "18"))))
