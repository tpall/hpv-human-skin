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

## ── Per-tier prevalence for the type catalogues (Tables S1 & S5) ─────────────
# Raw counts across tiers are NOT comparable: tier sizes differ >100-fold, so 14
# skin vs 14 control libraries mean utterly different prevalences. Each catalogue
# cell therefore also carries the prevalence *within* its tier (count / tier N).
# Prevalence is a Beta-Binomial posterior under a weakly-informative Jeffreys
# Beta(0.5, 0.5) prior: conjugate, hence exact (no MCMC): posterior mean
# (k+0.5)/(N+1), 95% credible interval from qbeta(). N = libraries screened in
# the tier (from the per-sample RPM table); shown in each column header.
.tier_N <- local({
  read_tsv(here("results_reframe_csc/hpv_rpm_per_sample.tsv"), show_col_types = FALSE) %>%
    count(tier) %>% { setNames(.$n, .$tier) }
})
.type_cat <- read_tsv(here("results_reframe_csc/type_catalogue_by_tier.tsv"), show_col_types = FALSE)
.tier_lab <- c(unselected_skin = "Unselected skin", control_productive = "Productive controls",
               neoplasia = "Neoplasia", cell_line = "Cell line", other = "Other")
.hdr_N <- function(tc) paste0(.tier_lab[[tc]], " (N=", comma(.tier_N[[tc]]), ")")

# Prevalence % formatting: ~2 sig figs, decimals scaled by magnitude.
fpct <- function(p){
  if (is.na(p)) return("NA")
  if (p >= 10) formatC(p, format = "f", digits = 0)
  else if (p >= 1) formatC(p, format = "f", digits = 1)
  else formatC(p, format = "f", digits = 2)
}
.post <- function(k, N){ a <- k + 0.5; b <- N - k + 0.5
  list(mean = 100 * a / (a + b), lo = 100 * qbeta(0.025, a, b), hi = 100 * qbeta(0.975, a, b)) }
# Catalogue cell "n (p%)": count plus posterior prevalence within the tier. Empty
# cells stay a bare "0": the posterior mean of a zero count is a small non-zero
# number, so printing it in every empty cell is noise (and "0 (1.1%)" reads as a
# contradiction); the caption gives credible intervals where absence is the point.
.np <- function(counts, tc){ N <- .tier_N[[tc]]
  vapply(counts, function(k)
    if (k == 0) "0" else paste0(k, " (", fpct(.post(k, N)$mean), "%)"), character(1)) }
# Caption helper "p% [lo–hi]": posterior prevalence + 95% CrI for a count in a tier.
prev_ci <- function(k, tc){ p <- .post(k, .tier_N[[tc]])
  paste0(fpct(p$mean), "% [", fpct(p$lo), "–", fpct(p$hi), "]") }
# Assigned-count lookup (Table S1 source) for caption CrIs.
kA  <- function(type, tc){ v <- .type_cat[[tc]][gsub("(REF|nr)$", "", .type_cat$hpv_reference) == type]
  if (length(v)) as.integer(v[1]) else 0L }
ciA <- function(type, tc) prev_ci(kA(type, tc), tc)

## ── Table S1: HPV type-by-tier catalogue (types with >=2 assigned libraries) ─
build_type_table <- function(){
  raw <- .type_cat %>% filter(total >= 2) %>% arrange(desc(total), desc(unselected_skin))
  tibble(`HPV type` = gsub("(REF|nr)$", "", raw$hpv_reference)) %>%
    mutate(!!.hdr_N("unselected_skin")    := .np(raw$unselected_skin,    "unselected_skin"),
           !!.hdr_N("control_productive") := .np(raw$control_productive, "control_productive"),
           !!.hdr_N("neoplasia")          := .np(raw$neoplasia,          "neoplasia"),
           !!.hdr_N("cell_line")          := .np(raw$cell_line,          "cell_line"),
           !!.hdr_N("other")              := .np(raw$other,              "other"),
           Total = raw$total)
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

## ── Table S5: detection-tier (trace) HPV type catalogue ─────────────────────
# Best-hit type from STAR primary alignments (unique, MAPQ 255) for HPV-positive
# but sub-threshold samples, i.e. below the breadth>=0.10 & depth>=2 assignment
# cutoff. `trace_types.tsv` is a one-time BAM re-map; genus from the complete
# PaVE lookup (`type_genus.csv`). Detection-tier, low-coverage, not productive.
.trace_tiers <- c("unselected_skin","neoplasia","control_productive","other","cell_line")
# Raw trace count matrix (type × tier), filtered to Total >= 2. Shared by the
# display builder and the caption CrI lookups so both see identical counts.
.trace_counts <- function(){
  tt  <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types = FALSE)
  rp  <- read_tsv(here("results_reframe_csc/hpv_rpm_per_sample.tsv"), show_col_types = FALSE) %>%
    select(srr_id, tier)
  gen <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
    select(hpv_reference, genus_short)
  x <- tt %>%
    filter(!top_ref %in% c("NONE","NO_BAM"), top_reads >= 5, max_mapq >= 20,
           top_reads / pmax(n_primary, 1) >= 0.5) %>%
    left_join(rp, by = "srr_id") %>% left_join(gen, by = c("top_ref" = "hpv_reference")) %>%
    mutate(type  = gsub("(REF|nr)$", "", top_ref),
           Genus = tools::toTitleCase(coalesce(genus_short, "unclassified")),
           tier  = factor(tier, levels = .trace_tiers)) %>%
    count(type, Genus, tier, .drop = FALSE) %>%
    pivot_wider(names_from = tier, values_from = n, values_fill = 0)
  for (col in .trace_tiers) if (!col %in% names(x)) x[[col]] <- 0L
  x %>% mutate(Total = rowSums(across(all_of(.trace_tiers)))) %>%
    filter(Total >= 2) %>% arrange(desc(Total), desc(unselected_skin))
}
build_trace_table <- function(){
  x <- .trace_counts()
  tibble(`HPV type` = x$type, Genus = x$Genus) %>%
    mutate(!!.hdr_N("unselected_skin")    := .np(x$unselected_skin,    "unselected_skin"),
           !!.hdr_N("neoplasia")          := .np(x$neoplasia,          "neoplasia"),
           !!.hdr_N("control_productive") := .np(x$control_productive, "control_productive"),
           !!.hdr_N("other")              := .np(x$other,              "other"),
           !!.hdr_N("cell_line")          := .np(x$cell_line,          "cell_line"),
           Total = x$Total)
}
# Trace-count lookup for caption CrIs.
kT  <- function(type, tc){ x <- .trace_counts(); v <- x[[tc]][x$type == type]
  if (length(v)) as.integer(v[1]) else 0L }
ciT <- function(type, tc) prev_ci(kT(type, tc), tc)
