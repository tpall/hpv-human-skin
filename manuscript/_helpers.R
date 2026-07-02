# Formatting helpers + reproducible table builders for the manuscript inline R.
# Sourced by the setup chunk in manuscript.qmd after results_models.R (which
# provides `stats`). Keeps display conventions in one place so every quoted
# number is computed, not hand-typed.
suppressPackageStartupMessages({library(tidyverse); library(here)})

## ── Number formatting ───────────────────────────────────────────────────────
.superscript <- c(`0`="⁰",`1`="¹",`2`="²",`3`="³",`4`="⁴",`5`="⁵",`6`="⁶",`7`="⁷",`8`="⁸",`9`="⁹",`-`="⁻")
sup <- function(n) paste(.superscript[strsplit(as.character(n),"")[[1]]], collapse = "")

# P-value: decimal (2 sig figs) at or above 0.001, else a×10^b with unicode superscript.
fp <- function(p){
  if (is.na(p)) return("NA")
  if (p >= 0.001) return(formatC(signif(p, 2), format = "g"))
  e <- floor(log10(p)); m <- p / 10^e
  paste0(formatC(m, format = "f", digits = 1), "×10", sup(e))
}

# Odds ratio / OR-CI bound: <10 -> 1 decimal, <1000 -> integer, else comma integer.
fo <- function(x){
  if (is.na(x)) return("NA")
  if (x >= 1000) return(formatC(round(x), format = "d", big.mark = ","))
  if (x >= 10)   return(formatC(round(x), format = "d"))
  formatC(x, format = "f", digits = 1)
}
oci <- function(lo, hi) paste0(fo(lo), "–", fo(hi))          # OR confidence interval

# Percentages: f2 = 2 decimals (positivity rates + their Wilson CIs).
f2  <- function(x) formatC(x, format = "f", digits = 2)
pci <- function(lo, hi) paste0(f2(lo), "–", f2(hi))

# Tier percentages: 1 decimal, drop a trailing ".0" (e.g. 93.0 -> "93").
pp  <- function(x) sub("\\.0$", "", formatC(x, format = "f", digits = 1))
p1  <- function(x) formatC(x, format = "f", digits = 1)      # always 1 decimal (floor, depth)
comma <- function(x) formatC(round(x), format = "d", big.mark = ",")
mreads <- function(x) paste0(p1(x), "M")                     # x already in millions

## ── Bayesian posterior formatting ───────────────────────────────────────────
# Posterior OR with 95% credible interval, e.g. "1.9 [0.7, 4.3]".
bor <- function(b) paste0(fo(b$OR), " [", fo(b$lo), ", ", fo(b$hi), "]")
# Posterior probability of direction; report ">0.999" instead of a rounded 1.00.
pdir <- function(p) if (p >= 0.9995) ">0.999" else formatC(p, format = "f", digits = 3)

## ── Table S3: Bayesian re-analysis of the main contrasts (weakly-inf. priors) ─
# Frequentist estimate from `stats`, posterior summary from `bayes` (bget()).
build_bayes_table <- function(stats, bget){
  rows <- tibble::tribble(
    ~contrast, ~freqor, ~freqp, ~bmodel,
    "Cell-line vs clinical positivity (M1)",                 stats$cl_or,        stats$cl_or_p,      "M1_cell_line",
    "Engineered vs clinical positivity (M1)",                stats$eng_or,       stats$eng_or_p,     "M1_engineered",
    "Productive controls vs unselected skin, floor (M2)",    stats$fish_or,      stats$fish_p,       "M2c_controls_vs_skin_floor",
    "Neoplasia vs unselected skin, floor (M2)",              stats$m2b_or,       stats$m2b_p,        "M2b_neoplasia_vs_skin_floor",
    "Sequencing depth, per 10× reads, floor (M2)",           stats$m2b_l10_or,   stats$m2b_l10_p,    "M2b_depth",
    "Sequencing depth within neoplasia, floor (M2)",         stats$neo_depth_or, stats$neo_depth_p,  "M2d_depth_within_neoplasia",
    "Neoplasia vs unselected skin, any-HPV / trace (M2)",    stats$m2a_or,       stats$m2a_p,        "M2a_neoplasia_vs_skin_any",
    "Productive controls vs neoplasia, productive (M3)",     stats$m3_or,        stats$m3_p,         "M3_controls_vs_neoplasia")
  purrr::pmap_dfr(rows, function(contrast, freqor, freqp, bmodel){
    b <- bget(bmodel)
    tibble(`Contrast (model)` = contrast,
           `Frequentist OR (P)` = paste0(fo(freqor), " (", fp(freqp), ")"),
           `Posterior OR [95% CrI]` = bor(b),
           `P(OR > 1)` = pdir(b$p_gt1))
  })
}

## ── Table S1: HPV type-by-tier catalogue (types with >=2 assigned libraries) ─
build_type_table <- function(){
  read_tsv(here("results_reframe_csc/type_catalogue_by_tier.tsv"), show_col_types = FALSE) %>%
    filter(total >= 2) %>%
    transmute(`HPV type` = gsub("(REF|nr)$", "", hpv_reference),
              `Unselected skin` = unselected_skin, `Productive controls` = control_productive,
              `Neoplasia` = neoplasia, `Cell line` = cell_line, `Other` = other, `Total` = total) %>%
    arrange(desc(Total), desc(`Unselected skin`))
}

## ── Table S2: genuine non-cell-line productive infections ───────────────────
# Curated sample list + annotations; coverage breadth / L1 / L2 read from data.
.prod_find <- function(srr){
  for (d in c("results_lesions_csc", "results_lesions_wide", "results_lesions"))
    if (file.exists(here(d, "aggregated", paste0(srr, "_hpv_coverage.tsv")))) return(d)
  NA_character_
}
build_prod_table <- function(){
  curated <- tibble::tribble(
    ~sample,        ~tissue,                       ~type,
    "ERR2274877",   "Other (carcinoma; CIN)",      "HPV96 (+HPV182)",
    "ERR1971044",   "Other (actinic keratosis)",   "HPV96 (+HPV182)",
    "SRR19520005",  "Skin (cSCC; CIN)",            "HPV12",
    "SRR19520060",  "Skin (cSCC; CIN)",            "HPV16",
    "SRR19520029",  "Skin (cSCC; CIN)",            "HPV16",
    "SRR19520008",  "Skin (normal; CIN)",          "HPV105",
    "SRR24991592",  "Skin (unspecified)",          "HPV227")
  purrr::pmap_dfr(curated, function(sample, tissue, type){
    d   <- .prod_find(sample)
    cov <- read_tsv(here(d, "aggregated", paste0(sample, "_hpv_coverage.tsv")), show_col_types = FALSE)
    br  <- max(cov$coverage_breadth[cov$coverage_breadth >= 0.10 & cov$mean_depth >= 2])
    tcf <- read_tsv(here(d, "aggregated", paste0(sample, "_transcript_classes.tsv")),
                    show_col_types = FALSE, col_types = cols(.default = "c"))
    l1 <- as.integer(tcf$read_count[tcf$gene == "L1"]); l2 <- as.integer(tcf$read_count[tcf$gene == "L2"])
    tibble(Sample = sample, `Tissue (diagnosis)` = tissue, `HPV type(s)` = type,
           `Coverage breadth` = formatC(br, format = "f", digits = 2),
           `L1 reads` = formatC(l1, format = "d", big.mark = ","),
           `L2 reads` = formatC(l2, format = "d", big.mark = ","))
  })
}
