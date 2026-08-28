#!/usr/bin/env Rscript
# Builds manuscript/figure5.png: the skin DNA reservoir this survey is silent
# against, read from an independent petabase-scale DNA catalogue.
#
#   Rscript manuscript/figure_reservoir.R
#
# Inputs are written by dna_rna_triangulation/shen_reservoir_by_tissue.R (the
# per-body-site reservoir tabulation of Shen et al.'s Logan papillomavirus
# survey) and dna_rna_triangulation/skin_dna_vs_rna_corrected.R (the type-level
# join onto our own RNA tiers, assigned and trace).
#
# Figure 5 is skin only, on the principle that a skin paper's main figure should
# not lead with a vaginal panel:
#
#   A  the ten most frequent types in the skin DNA reservoir, each carrying the
#      RNA-side count from this study (assigned / trace), so the reservoir and
#      its transcriptional silence sit in one panel.
#   B  every type, DNA abundance against RNA productivity. The abundant-and-
#      productive corner is empty but for HPV2, which is the presence-is-not-
#      activity claim made as a picture rather than asserted.
#
# The body-site comparison (skin, vaginal, gut) moves to Figure S3, where it is
# the positive control for the read-out itself: if the composition were assembly
# noise it could not recover alpha tropism in the vaginal libraries, or the
# collection-route mixture in gut.
#
# Composition, not prevalence: the denominator is HPV+ libraries at that site,
# so nothing here depends on the unsettled DNA prevalence denominator.
suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(patchwork); library(here)
  library(ggrepel)
})
source(here("manuscript/_paths.R"))
old <- theme_set(theme_cowplot(font_size = 11))

SITES <- c("Skin", "Vaginal")
# Unresolved is not a genus but the residue, so it sits last in the legend and
# at the bottom of each stack, after the named genera and the novel clusters.
GEN_LEVELS <- c("Alpha", "Beta", "Gamma", "Mu", "Nu", "Novel", "Unresolved")
gcol <- c(Alpha = "#1B9E77", Beta = "#D95F02", Gamma = "#7570B3", Mu = "#E7298A",
          Nu = "#66A61E", Unresolved = "grey70", Novel = "grey40")

res_all <- read_csv(here("dna_rna_triangulation/reservoir_by_tissue_genus.csv"),
                    show_col_types = FALSE)
tt_all  <- read_csv(here("dna_rna_triangulation/reservoir_by_tissue_top_types.csv"),
                    show_col_types = FALSE)
cmp <- read_csv(here("dna_rna_triangulation/skin_dna_vs_rna_types_corrected.csv"),
                show_col_types = FALSE)

# RNA denominator for the annotation: our unselected-skin tier.
rpm <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE)
n_skin_rna <- sum(rpm$tier == "unselected_skin")
n_rna_all  <- nrow(rpm)

# Only the supplement is built from a site set now; Figure 5 has its own builder.
SITE_SETS <- list(figureS4 = c(SITES, "Gut"))

make_fig <- function(sites) {
  res <- res_all %>% filter(tissue %in% sites)
  tt  <- tt_all  %>% filter(tissue %in% sites) %>%
    group_by(tissue) %>% slice_max(libs, n = TOP_SUPP, with_ties = TRUE) %>% ungroup()
  gen_present <- intersect(GEN_LEVELS, union(res$genus, tt$genus))
  lab_map <- res %>% distinct(tissue, n_lib) %>%
    mutate(lab = sprintf("%s\n(n = %s HPV+ libraries)", tissue,
                         formatC(n_lib, format = "d", big.mark = ","))) %>%
    arrange(factor(tissue, sites)) %>% select(tissue, lab) %>% deframe()

  # A. Genus composition. A co-infected library contributes to more than one
  # genus, so library prevalences sum above 100%; normalising within site turns
  # them into the composition of type detections, which is what is compared here.
  pA <- res %>% filter(genus != "Novel_any") %>%
    group_by(tissue) %>% mutate(comp = libs / sum(libs)) %>% ungroup() %>%
    mutate(genus = factor(genus, gen_present),
           tissue = factor(lab_map[tissue], lab_map[sites])) %>%
    ggplot(aes(tissue, comp, fill = genus)) +
    geom_col(width = .6, colour = "white", linewidth = .2) +
    scale_fill_manual(values = gcol, name = "Genus", limits = gen_present, drop = FALSE) +
    scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .02))) +
    labs(x = NULL, y = "Genus composition\n(% of type detections)")

  # B. Top ten types per site. Skin rows carry our RNA counts, assigned / trace.
  d <- tt %>%
    mutate(type_n = str_remove(type, "^HPV")) %>%
    left_join(cmp %>% transmute(type_n = str_remove(type, "^HPV"),
                                rna_assigned_skin, rna_trace_skin),
              by = "type_n") %>%
    mutate(genus = factor(genus, gen_present),
           row = fct_reorder(paste0(type, "__", tissue), libs),
           rna_lab = if_else(tissue == "Skin",
                             sprintf("%d / %d", replace_na(rna_assigned_skin, 0L),
                                     replace_na(rna_trace_skin, 0L)), NA_character_),
           tissue = factor(lab_map[tissue], lab_map[sites]))

  pB <- ggplot(d, aes(prop, row, fill = genus)) +
    geom_col(width = .7, colour = "white", linewidth = .2) +
    geom_text(aes(label = rna_lab), hjust = 0, nudge_x = .004, size = 3,
              colour = "grey25", na.rm = TRUE) +
    facet_wrap(~ tissue, scales = "free_y", nrow = 1) +
    scale_y_discrete(labels = function(x) str_remove(x, "__.*")) +
    scale_fill_manual(values = gcol, name = "Genus", limits = gen_present, drop = FALSE) +
    scale_x_continuous(labels = scales::percent, expand = expansion(c(0, .18))) +
    guides(fill = "none") +
    labs(x = "% of that site's HPV+ DNA libraries", y = NULL)

  (pA / pB) + plot_layout(heights = c(1, 1.6), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"), legend.position = "right")
}

## ── Figure 5: skin only, reservoir against activity ─────────────────────────
# Panel B is lifted from dna_rna_triangulation/shen_skin_dna_vs_rna.R, which
# keeps building the two-panel triangulation figure for the companion paper.
# The shaded corner marks types that are both abundant in skin DNA and
# productive in our RNA; DNA_ABUNDANT / RNA_PRODUCTIVE are its two edges.
# How deep each panel goes into the ranked reservoir. The main figure runs to the
# tie group at rank 14 so that it reaches HPV2, the one type panel B identifies as
# both abundant in skin DNA and productive here; ranks 12 to 14 are tied on library
# count, so stopping short of 14 would split a tie. The supplement stays at ten per
# site, where its job is the cross-site comparison rather than depth.
TOP_MAIN <- 14
TOP_SUPP <- 10

DNA_ABUNDANT  <- 15
RNA_PRODUCTIVE <- 5
lg <- function(x) log10(x + 1)

make_fig5 <- function() {
  gen_present <- intersect(GEN_LEVELS, union(tt_all$genus[tt_all$tissue == "Skin"],
                                             as.character(cmp$genus)))
  n_skin <- res_all$n_lib[res_all$tissue == "Skin"][1]

  # A. Top ten skin types, with our RNA counts (assigned / trace) alongside.
  d <- tt_all %>% filter(tissue == "Skin") %>%
    slice_max(libs, n = TOP_MAIN, with_ties = TRUE) %>%
    mutate(type_n = str_remove(type, "^HPV")) %>%
    left_join(cmp %>% transmute(type_n = str_remove(type, "^HPV"),
                                rna_assigned_skin, rna_trace_skin),
              by = "type_n") %>%
    arrange(desc(libs)) %>%
    # Explicit levels, not fct_reorder: HPV107 and HPV173 are tied on library
    # count, and a reorder breaks that tie differently from the arrange() that
    # writes topn_types, which would put figure and Results text in different
    # orders for the same ten types.
    mutate(genus = factor(genus, gen_present),
           row = factor(type, levels = rev(type)))

  # The RNA side used to be a bare "0 / 1" beside each bar, which no reader can
  # decode without the caption. It is now a labelled two-column strip: a header
  # naming what the columns are, one cell per count, and empty cells left white
  # so that "the assigned column is blank all the way down" is the thing the eye
  # lands on. Tile fill is passed outside aes() because the genus fill scale is
  # already spoken for by the bars.
  xmax <- max(d$prop)
  colx <- c(assigned = xmax * 1.26, trace = xmax * 1.64)
  tw   <- xmax * 0.24
  nrow_d <- nrow(d)
  tiles <- d %>%
    transmute(row, assigned = replace_na(rna_assigned_skin, 0L),
              trace = replace_na(rna_trace_skin, 0L)) %>%
    pivot_longer(c(assigned, trace), names_to = "kind", values_to = "n") %>%
    mutate(x = unname(colx[kind]),
           fillcol = if_else(n > 0, "grey78", "white"),
           txtcol  = if_else(n > 0, "grey10", "grey60"))

  pA <- ggplot(d, aes(prop, row)) +
    geom_col(aes(fill = genus), width = .7, colour = "white", linewidth = .2) +
    geom_tile(data = tiles, aes(x = x, y = row), fill = tiles$fillcol,
              colour = "grey75", linewidth = .3, width = tw, height = .72) +
    geom_text(data = tiles, aes(x = x, y = row, label = n),
              colour = tiles$txtcol, size = 2.9) +
    # Headers sit ABOVE the panel, not inside it. Expanding the y scale to make
    # room for them stretched the axis line into a bare stub over the top bar;
    # with clip = "off" they draw in the top margin instead and the axis stops
    # where the data does.
    annotate("text", x = colx, y = nrow_d + 0.95, size = 2.5, colour = "grey25",
             fontface = "bold", label = c("assigned", "trace")) +
    annotate("text", x = mean(colx), y = nrow_d + 1.75, size = 2.6, colour = "grey25",
             label = "RNA detections in unselected skin") +
    scale_fill_manual(values = gcol, name = "Genus", limits = gen_present, drop = FALSE) +
    # Breaks stop where the bars stop. The tile strip lives further right on the
    # same axis, and ticks running under it would invite reading the cells off
    # the percentage scale. Expansion is small because the tiles already
    # extend the data range.
    scale_x_continuous(labels = scales::percent, breaks = c(0, .02, .04, .06),
                       expand = expansion(c(0, .06))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) +
    coord_cartesian(clip = "off") +
    guides(fill = "none") +
    labs(x = sprintf("%% of the %s HPV+ skin DNA libraries",
                     formatC(n_skin, format = "d", big.mark = ",")), y = NULL) +
    theme(plot.margin = margin(t = 30, r = 6, b = 6, l = 6))

  # B. Every type: how often it is present in skin DNA against how often it is
  # productive in our RNA. Labelled where either axis is at its extreme.
  lbrk <- c(0, 1, 5, 20, 100); lval <- lg(lbrk)
  lab_types <- cmp %>%
    filter(rank(-dna_skin_libs, ties.method = "first") <= 8 |
           rank(-rna_productive, ties.method = "first") <= 8) %>%
    mutate(genus = factor(as.character(genus), gen_present))

  pB <- cmp %>% mutate(genus = factor(as.character(genus), gen_present)) %>%
    ggplot(aes(lg(dna_skin_libs), lg(rna_productive), colour = genus)) +
    annotate("rect", xmin = lg(DNA_ABUNDANT), xmax = Inf,
             ymin = lg(RNA_PRODUCTIVE), ymax = Inf, fill = "grey85", alpha = .5) +
    # Seeded, for the same reason as the jitter in bin/hpv_rpm_contrast.R.
    geom_point(size = 1.9, alpha = .85,
               position = position_jitter(width = .02, height = .02, seed = 1)) +
    geom_text_repel(data = lab_types, aes(label = type), size = 2.9,
                    max.overlaps = 20, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = gcol, name = "Genus", limits = gen_present, drop = FALSE) +
    scale_x_continuous(breaks = lval, labels = lbrk) +
    scale_y_continuous(breaks = lval, labels = lbrk) +
    labs(x = "Skin DNA libraries carrying the type",
         y = "Productive RNA libraries (warts and EV)") +
    # Same top margin as panel A, whose margin holds the column headers, so the
    # two plotting regions start at the same height instead of stepping.
    theme(plot.margin = margin(t = 30, r = 6, b = 6, l = 6))

  (pA | pB) + plot_layout(widths = c(1, 1.35), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"), legend.position = "right")
}

ggsave(here("manuscript/figure5.png"), make_fig5(), width = 10.5, height = 4.6,
       dpi = 200, bg = "white")
ggsave(here("manuscript/figure5.pdf"), make_fig5(), width = 10.5, height = 4.6, bg = "white")

for (nm in names(SITE_SETS)) {
  p <- make_fig(SITE_SETS[[nm]])
  ggsave(here(paste0("manuscript/", nm, ".png")), p, width = 10.5, height = 7.6,
         dpi = 200, bg = "white")
  ggsave(here(paste0("manuscript/", nm, ".pdf")), p, width = 10.5, height = 7.6, bg = "white")
}

# Persist every number the caption and Results text quote, so they are read from
# the same object the figure plots rather than hand-typed.
gp <- function(tis, g) {
  x <- res_all$prop[res_all$tissue == tis & res_all$genus == g]
  if (length(x)) 100 * x else 0
}
top_skin <- tt_all %>% filter(tissue == "Skin") %>%
  slice_max(libs, n = TOP_MAIN, with_ties = TRUE) %>% arrange(desc(libs)) %>%
  mutate(type_n = str_remove(type, "^HPV")) %>%
  left_join(cmp %>% transmute(type_n = str_remove(type, "^HPV"),
                              rna_assigned_skin, rna_trace_skin), by = "type_n")
nums <- tibble::tibble(
  key = c("n_skin_lib", "n_vaginal_lib", "n_gut_lib", "n_skin_rna", "n_rna_all",
          "skin_beta", "skin_gamma", "skin_alpha", "skin_mu", "skin_unresolved",
          "skin_novel_any", "vaginal_alpha", "vaginal_beta", "vaginal_gamma",
          "vaginal_unresolved", "vaginal_novel_any", "gut_alpha", "gut_beta", "gut_gamma",
          "gut_novel_any", "topn_rna_assigned", "topn_rna_trace",
          "topn_types", "top_n", "top_type", "top_type_pct", "hpv8_dna", "hpv8_rna_all",
          "n_types_cmp", "n_dna_present", "n_dna_present_silent",
          "n_prod_types", "n_prod_dna_absent", "corner_types",
          "dna_abundant_cut", "rna_productive_cut"),
  value = c(res_all$n_lib[res_all$tissue == "Skin"][1],
            res_all$n_lib[res_all$tissue == "Vaginal"][1],
            res_all$n_lib[res_all$tissue == "Gut"][1],
            n_skin_rna, n_rna_all,
            gp("Skin", "Beta"), gp("Skin", "Gamma"), gp("Skin", "Alpha"),
            gp("Skin", "Mu"), gp("Skin", "Unresolved"), gp("Skin", "Novel_any"),
            gp("Vaginal", "Alpha"), gp("Vaginal", "Beta"), gp("Vaginal", "Gamma"),
            gp("Vaginal", "Unresolved"), gp("Vaginal", "Novel_any"),
            gp("Gut", "Alpha"), gp("Gut", "Beta"), gp("Gut", "Gamma"),
            gp("Gut", "Novel_any"),
            sum(top_skin$rna_assigned_skin, na.rm = TRUE),
            sum(top_skin$rna_trace_skin, na.rm = TRUE),
            NA, nrow(top_skin), NA, 100 * max(top_skin$prop),
            cmp$dna_skin_libs[cmp$type == "HPV8"],
            cmp$rna_any_all[cmp$type == "HPV8"],
            # Panel B: how empty the abundant-and-productive corner is.
            nrow(cmp),
            sum(cmp$dna_skin_libs > 0),
            sum(cmp$dna_skin_libs > 0 & cmp$rna_productive == 0),
            sum(cmp$rna_productive > 0),
            sum(cmp$rna_productive > 0 & cmp$dna_skin_libs == 0),
            NA, DNA_ABUNDANT, RNA_PRODUCTIVE) %>% as.character()) %>%
  mutate(value = case_when(key == "topn_types" ~ paste(top_skin$type, collapse = ", "),
                           key == "top_type" ~ top_skin$type[1],
                           key == "corner_types" ~ paste(sort(cmp$type[
                             cmp$dna_skin_libs >= DNA_ABUNDANT &
                             cmp$rna_productive >= RNA_PRODUCTIVE]), collapse = ", "),
                           TRUE ~ value))
write_tsv(nums, here("manuscript/reservoir_numbers.tsv"))

cat("Wrote manuscript/figure5.{png,pdf}, figureS4.{png,pdf} and reservoir_numbers.tsv\n")
print(as.data.frame(nums), right = FALSE)
