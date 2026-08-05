# Formatting helpers + reproducible table builders for the manuscript inline R.
# Sourced by the setup chunk in manuscript.qmd after results_models.R (which
# provides `stats`). Keeps display conventions in one place so every quoted
# number is computed, not hand-typed.
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("manuscript/_paths.R"))

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
# An unbounded value (complete separation) prints as ∞, never as a 16-digit number.
fo <- function(x){
  if (is.na(x)) return("NA")
  if (!is.finite(x)) return("∞")
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
  read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>%
    count(tier) %>% { setNames(.$n, .$tier) }
})
.type_cat <- read_tsv(type_catalogue_tsv(), show_col_types = FALSE)
.tier_lab <- c(unselected_skin = "Unselected skin", control_productive = "Productive controls",
               EV = "EV", peri_lesional = "Peri-lesional mucosa",
               non_permissive_tissue = "Non-permissive tissue",
               engineered = "Engineered cells",
               # "other" is the residual tier: 99% of it carries no site keyword
               # in its SRA metadata, so it is an absence of information rather
               # than an identified non-cutaneous site. Keep the label short for
               # table headers but do not call it a tissue.
               neoplasia = "Neoplasia", cell_line = "Cell line", other = "Site not recorded")
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
           !!.hdr_N("EV")                 := .np(raw$EV,                 "EV"),
           !!.hdr_N("peri_lesional")      := .np(raw$peri_lesional,      "peri_lesional"),
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
    # Tissue/diagnosis strings track the enriched samplesheets. The "; CIN"
    # tails present here previously were an artefact of bare-substring keyword
    # matching ("cin" inside "carCINoma"); none of these libraries is cervical.
    # ERR2274877/ERR1971044 are ERP106532/ERP022968, the two cutaneous AK/cSCC
    # studies whose SRA records carry no tissue text (see assets/tier_overrides.csv).
    "ERR2274877",   "Skin (carcinoma)",            "HPV96 (+HPV182)",
    "ERR1971044",   "Skin (actinic keratosis)",    "HPV96 (+HPV182)",
    "SRR19520005",  "Skin (cSCC)",                 "HPV12",
    "SRR19520060",  "Skin (cSCC)",                 "HPV16",
    "SRR19520029",  "Skin (cSCC)",                 "HPV16",
    "SRR19520008",  "Skin (normal)",               "HPV105",
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
# Peri-lesional mucosa contributes no trace calls (all 15 libraries clear the
# assignment threshold and therefore appear in Table S1), so it has no column here.
.trace_tiers <- c("unselected_skin","neoplasia","control_productive","EV","other","cell_line")
# Raw trace count matrix (type × tier), filtered to Total >= 2. Shared by the
# display builder and the caption CrI lookups so both see identical counts.
.trace_counts <- function(){
  tt  <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types = FALSE)
  rp  <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>%
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
           !!.hdr_N("EV")                 := .np(x$EV,                 "EV"),
           !!.hdr_N("other")              := .np(x$other,              "other"),
           !!.hdr_N("cell_line")          := .np(x$cell_line,          "cell_line"),
           Total = x$Total)
}
# Trace-eligible libraries per tier: HPV-positive but with no reference passing
# the breadth/depth assignment criterion. These are the only libraries that can
# appear in Table S5 at all, so the caption quotes them alongside the tier N
# (the denominator of the printed prevalences, which stays the tier N so that
# Tables S1 and S5 remain on one scale).
.trace_eligible_N <- local({
  assigned <- unique(unlist(lapply(
    c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide"),
    function(d) sub("_hpv_types\\.tsv$", "",
                    list.files(here(d, "aggregated"), "_hpv_types\\.tsv$")))))
  # a library has an assignment only if its types file is non-empty
  has_call <- function(srr) {
    for (d in c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide")) {
      f <- here(d, "aggregated", paste0(srr, "_hpv_types.tsv"))
      if (file.exists(f) && length(readLines(f, warn = FALSE)) > 1) return(TRUE)
    }
    FALSE
  }
  rp <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% filter(hpv_reads >= 10)
  rp %>% mutate(assigned = map_lgl(srr_id, ~ .x %in% assigned && has_call(.x))) %>%
    filter(!assigned) %>% count(tier) %>% { setNames(.$n, .$tier) }
})
trace_elig <- function(tc) unname(coalesce(.trace_eligible_N[tc], 0L))
# Human-readable "Unselected skin 196, Cell line 35, ..." for the caption, ordered
# as the table columns are, and using the display labels rather than the raw keys.
trace_elig_str <- function(tiers = .trace_tiers)
  paste(sprintf("%s %d", .tier_lab[tiers], vapply(tiers, trace_elig, integer(1))),
        collapse = ", ")
# Trace-count lookup for caption CrIs.
kT  <- function(type, tc){ x <- .trace_counts(); v <- x[[tc]][x$type == type]
  if (length(v)) as.integer(v[1]) else 0L }
ciT <- function(type, tc) prev_ci(kT(type, tc), tc)

## ── Table S6: multi-type (co-infection) libraries ───────────────────────────
# Every library with two or more assigned references, i.e. two or more types
# each clearing breadth >= 0.10 and depth >= 2. Coverage comes from
# `samtools coverage -q 10`, so only uniquely-mapped reads (STAR MAPQ 255)
# contribute: a read that maps equally well to two HPV references is a
# multimapper (MAPQ <= 3) and is discarded before these numbers are computed.
# Cross-mapping between types therefore cannot manufacture a co-infection call.
.coinf_dirs <- c("results_full_v2", "results_lesions", "results_lesions_csc", "results_lesions_wide")
.coinf <- local({
  ty <- map_dfr(unlist(lapply(.coinf_dirs, function(d)
      list.files(here(d, "aggregated"), "_hpv_types\\.tsv$", full.names = TRUE))),
    ~ suppressMessages(read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))) %>%
    distinct(sample_id, hpv_reference, .keep_all = TRUE) %>%
    mutate(coverage_breadth = as.double(coverage_breadth), mean_depth = as.double(mean_depth),
           type = gsub("(REF|nr)$", "", hpv_reference))
  rp <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% select(srr_id, tier, study)
  k  <- ty %>% count(sample_id, name = "n_types")
  list(types = ty %>% left_join(k, by = "sample_id") %>%
                      left_join(rp, by = c("sample_id" = "srr_id")),
       n_typed = nrow(k))
})
.coinf_tier_order <- c("control_productive", "EV", "peri_lesional", "neoplasia",
                       "non_permissive_tissue", "engineered",
                       "unselected_skin", "other", "cell_line")
build_coinf_table <- function(){
  .coinf$types %>% filter(n_types >= 2) %>%
    group_by(sample_id, tier, study) %>%
    arrange(desc(mean_depth), .by_group = TRUE) %>%
    summarise(dom_depth = max(mean_depth),
              types = paste(sprintf("%s (%.2f, %s×)", type, coverage_breadth,
                                    formatC(round(mean_depth), format = "d", big.mark = ",")),
                            collapse = " + "),
              .groups = "drop") %>%
    mutate(tier = factor(tier, levels = .coinf_tier_order)) %>%
    arrange(tier, desc(dom_depth)) %>%
    transmute(Library = sample_id, Tier = .tier_lab[as.character(tier)], Study = study,
              `Assigned types (coverage breadth, mean depth), ordered by depth` = types)
}
# Co-infection counts quoted in the caption and in the Results paragraph.
.CONDYLOMA_PAIR <- c("SRP314958", "SRP314961")
cf <- local({
  co  <- .coinf$types %>% filter(n_types >= 2)
  lib <- co %>% distinct(sample_id, study, n_types)
  sec <- co %>% group_by(sample_id) %>% arrange(desc(mean_depth), .by_group = TRUE) %>%
    slice(-1) %>% ungroup()                      # non-dominant calls
  # How lopsided is a typical co-infection? Depth of the deepest type over the
  # second deepest, per library.
  ratio <- co %>% group_by(sample_id) %>% arrange(desc(mean_depth), .by_group = TRUE) %>%
    summarise(r = first(mean_depth) / nth(mean_depth, 2), .groups = "drop") %>% pull(r)
  w <- .coinf$types %>% filter(type %in% c("HPV11", "HPV42")) %>%
    distinct(sample_id, type) %>% mutate(v = TRUE) %>%
    pivot_wider(names_from = type, values_from = v, values_fill = FALSE)
  # Donor identity for the condyloma study pair: SRP314958 is RNA-seq and
  # SRP314961 is MeRIP-seq of the same material, whose run titles are CA_n /
  # CON_n with an _input or _IP suffix on the MeRIP fractions.
  titles <- map_dfr(list.files(here(), "lesions_raw.csv", recursive = TRUE, full.names = TRUE),
      ~ suppressMessages(read_csv(.x, show_col_types = FALSE)) %>% mutate(across(everything(), as.character))) %>%
    distinct(srr_id, .keep_all = TRUE) %>%
    transmute(sample_id = srr_id, donor = str_remove(title, "_(IP|input)$"))
  pair <- co %>% filter(study %in% .CONDYLOMA_PAIR) %>% left_join(titles, by = "sample_id")
  # An HPV11 + HPV42 call is "unambiguous" only when *both* types clear a
  # coverage bar far above the assignment threshold (breadth 0.10, depth 2).
  strict <- .coinf$types %>% filter(type %in% c("HPV11", "HPV42"),
                                    coverage_breadth >= 0.70, mean_depth >= 5) %>%
    count(sample_id, name = "n") %>% filter(n == 2) %>%
    left_join(titles, by = "sample_id")
  list(libs = nrow(lib), typed = .coinf$n_typed,
       k2 = sum(lib$n_types == 2), k3 = sum(lib$n_types == 3),
       k4 = sum(lib$n_types == 4), k5 = sum(lib$n_types == 5),
       pair_libs = sum(lib$study %in% .CONDYLOMA_PAIR),
       pair_donors = n_distinct(pair$donor),
       sec_med_depth = median(sec$mean_depth),
       sec_shallow   = round(100 * mean(sec$mean_depth < 5)),
       dom_ratio_med = median(ratio),
       dom_ratio_lo  = unname(quantile(ratio, .25)),
       dom_ratio_hi  = unname(quantile(ratio, .75)),
       both1142 = sum(w$HPV11 & w$HPV42),
       only42   = sum(!w$HPV11 & w$HPV42),
       only11   = sum(w$HPV11 & !w$HPV42),
       strict_libs = nrow(strict), strict_donors = n_distinct(strict$donor))
})

## ── Figure S2: genus prevalence by tier ─────────────────────────────────────
# manuscript/figure_genus.R plots the figure and writes the counts it plotted.
# The caption reads that file, so figure and caption cannot drift apart.
.genus_prev <- read_tsv(here("manuscript/genus_prevalence_by_tier.tsv"), show_col_types = FALSE)
.gp_row <- function(genus, tier, level = "Any detection"){
  r <- .genus_prev[.genus_prev$genus == genus & .genus_prev$tier == tier &
                     .genus_prev$level == level, ]
  if (!nrow(r)) stop("no genus_prevalence row for ", genus, " / ", tier, " / ", level)
  r[1, ]
}
# Libraries detected: "10/10".
kG  <- function(genus, tier, level = "Any detection"){
  r <- .gp_row(genus, tier, level); paste0(r$k, "/", r$N) }
# Posterior prevalence with credible interval: "95% [78–100]".
ciG <- function(genus, tier, level = "Any detection"){
  r <- .gp_row(genus, tier, level)
  paste0(fpct(r$mean), "% [", fpct(r$lo), "–", fpct(r$hi), "]") }
