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
