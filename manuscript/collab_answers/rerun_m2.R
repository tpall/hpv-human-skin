#!/usr/bin/env Rscript
# Refit Model 2 (depth-controlled detection across population tiers) on the
# tier-corrected RPM table and print it beside the published version.
#
#   Rscript manuscript/collab_answers/rerun_m2.R
#
# Published table : results_reframe_csc/hpv_rpm_per_sample.tsv
# Corrected table : results_reframe_fixed/hpv_rpm_per_sample.tsv
#   (bin/hpv_rpm_contrast.R + assets/tier_overrides.csv)
suppressPackageStartupMessages({library(tidyverse); library(here)})
options(width = 200)

ortab <- function(m){ s <- summary(m)$coefficients
  data.frame(term = rownames(s), OR = exp(s[,1]), lo = exp(s[,1]-1.96*s[,2]),
             hi = exp(s[,1]+1.96*s[,2]), p = s[,4], row.names = NULL) }
row_of <- function(o, tm) as.list(o[o$term == tm, ])
# Odds ratios here span 1 to 1e15 (near-complete separation), so switch to
# scientific notation once a bound stops being readable as a plain number.
num <- function(x) if (!is.finite(x)) "Inf" else if (abs(x) >= 1e5) formatC(x, format="e", digits=1) else
  if (abs(x) >= 100) formatC(x, format="d") else formatC(x, format="f", digits=2)
pval <- function(p) if (p < 1e-4) formatC(p, format="e", digits=1) else formatC(p, format="f", digits=4)
fmt <- function(r) if (is.null(r$OR) || !length(r$OR)) "n/a" else
  sprintf("OR %s [%s-%s], P = %s", num(r$OR), num(r$lo), num(r$hi), pval(r$p))
# Fisher exact, for contrasts where the glm Wald interval breaks down under
# separation (every library in the tier is above the floor).
fisher_vs_skin <- function(rp, other, floor_rpm){
  d <- rp %>% filter(tier %in% c("unselected_skin", other))
  if (n_distinct(d$tier) < 2) return(list())
  f <- fisher.test(table(tier = factor(as.character(d$tier),
                                       levels = c("unselected_skin", other)),
                         detect = d$detect))
  list(OR = unname(f$estimate), lo = f$conf.int[1], hi = f$conf.int[2], p = f$p.value)
}

# Dominant assigned type per library, for the CIN (mucosal high-risk) sensitivity.
.HR <- c(16,18,31,33,35,39,45,51,52,56,58,59,68,26,53,66,67,69,70,73,82)
domtype <- function(dirs){
  fs <- unlist(lapply(dirs, function(d) list.files(here(d,"aggregated"), "_hpv_types\\.tsv$", full.names = TRUE)))
  map_dfr(fs, ~ suppressMessages(read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))) %>%
    mutate(coverage_breadth = as.numeric(coverage_breadth), mean_depth = as.numeric(mean_depth)) %>%
    filter(coverage_breadth >= 0.10, mean_depth >= 2) %>%
    group_by(sample_id) %>% slice_max(mean_depth, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(srr_id = sample_id, dom_type = hpv_reference,
              mucosal_hr = suppressWarnings(as.integer(str_match(dom_type, "HPV[-]?(\\d+)")[,2])) %in% .HR)
}
DT <- domtype(c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide"))

fit <- function(path, label){
  rp <- read_tsv(here(path), show_col_types = FALSE) %>% filter(!is.na(total_reads), total_reads > 0)
  floor_rpm <- rp %>% filter(tier == "control_productive", !is_cell_line, rpm > 0) %>%
    pull(rpm) %>% quantile(.05, na.rm = TRUE) %>% as.numeric()
  rp <- rp %>% mutate(any_hpv = rpm > 0, detect = rpm >= floor_rpm, l10 = log10(total_reads)) %>%
    left_join(DT %>% select(srr_id, mucosal_hr), by = "srr_id") %>%
    mutate(mucosal_hr = coalesce(mucosal_hr, FALSE))

  tiersum <- rp %>% group_by(tier) %>%
    summarise(n = n(), med_reads_M = round(median(total_reads)/1e6, 1),
              pct_any = round(100*mean(any_hpv), 1), pct_floor = round(100*mean(detect), 1),
              med_rpm = round(median(rpm), 3), .groups = "drop")

  vs_skin <- function(other, outcome){
    d <- rp %>% filter(tier %in% c("unselected_skin", other)) %>%
      mutate(tier = relevel(factor(as.character(tier)), ref = "unselected_skin"))
    if (n_distinct(d$tier) < 2) return(list())
    f <- as.formula(paste(outcome, "~ tier + l10"))
    row_of(ortab(glm(f, binomial, d)), paste0("tier", other))
  }
  # Fisher for the near-separated controls-vs-skin contrast
  ct <- rp %>% filter(tier %in% c("control_productive","unselected_skin"))
  ft <- fisher.test(table(tier = factor(as.character(ct$tier),
                                        levels = c("unselected_skin","control_productive")),
                          detect = ct$detect))
  # depth only, within neoplasia
  mn <- row_of(ortab(glm(detect ~ l10, binomial, rp %>% filter(tier == "neoplasia"))), "l10")
  # CIN-cut: drop mucosal high-risk cases from the neoplasia tier
  rpc <- rp %>% filter(!(tier == "neoplasia" & mucosal_hr))
  subc <- rpc %>% filter(tier %in% c("unselected_skin","neoplasia")) %>%
    mutate(tier = relevel(factor(as.character(tier)), "unselected_skin"))
  m2b_cut <- row_of(ortab(glm(detect ~ tier + l10, binomial, subc)), "tierneoplasia")
  neo_cut <- rpc %>% filter(tier == "neoplasia") %>% summarise(n = n(), pct = 100*mean(detect))

  list(label = label, floor = floor_rpm, tiersum = tiersum,
       m2a_neo = vs_skin("neoplasia", "any_hpv"),
       m2b_neo = vs_skin("neoplasia", "detect"),
       # EV (8/10 above floor) and peri-lesional mucosa (15/15) are separated or
       # near-separated against skin, so the depth-adjusted glm is not usable.
       m2b_ev  = fisher_vs_skin(rp, "EV", floor_rpm),
       m2b_per = fisher_vs_skin(rp, "peri_lesional", floor_rpm),
       depth_any = row_of(ortab(glm(any_hpv ~ tier + l10, binomial,
                    rp %>% filter(tier %in% c("unselected_skin","neoplasia")) %>%
                    mutate(tier = relevel(factor(as.character(tier)), "unselected_skin")))), "l10"),
       fisher = ft, neo_depth = mn, m2b_cut = m2b_cut, neo_cut = neo_cut,
       n_cin = sum(rp$tier == "neoplasia" & rp$mucosal_hr))
}

A <- fit("results_reframe_csc/hpv_rpm_per_sample.tsv",   "A. published (results_reframe_csc)")
B <- fit("results_reframe_fixed/hpv_rpm_per_sample.tsv", "B. tier-corrected (results_reframe_fixed)")

hdr <- function(s) cat("\n", strrep("=", 78), "\n", s, "\n", strrep("=", 78), "\n", sep = "")
for (x in list(A, B)) {
  hdr(x$label)
  cat(sprintf("Control-calibrated detection floor: %.2f RPM\n\n", x$floor))
  print(as.data.frame(x$tiersum), right = FALSE)
}

hdr("MODEL 2 CONTRASTS, published vs tier-corrected")
cmp <- tribble(
  ~Contrast,                                              ~A,                    ~B,
  "M2a  neoplasia vs skin, ANY HPV read",                 fmt(A$m2a_neo),        fmt(B$m2a_neo),
  "M2a  sequencing depth (per 10x reads), any HPV",       fmt(A$depth_any),      fmt(B$depth_any),
  "M2b  neoplasia vs skin, at detection floor",           fmt(A$m2b_neo),        fmt(B$m2b_neo),
  "M2b  neoplasia vs skin, floor, CIN cases removed",     fmt(A$m2b_cut),        fmt(B$m2b_cut),
  "M2d  depth within neoplasia, at floor",                fmt(A$neo_depth),      fmt(B$neo_depth),
  "  EV vs skin, at floor (Fisher, separated)",           fmt(A$m2b_ev),         fmt(B$m2b_ev),
  "  peri-lesional vs skin, at floor (Fisher, separated)", fmt(A$m2b_per),       fmt(B$m2b_per))
names(cmp) <- c("Contrast", A$label, B$label)
print(as.data.frame(cmp), right = FALSE)

hdr("M2c  productive controls vs unselected skin, at floor (Fisher exact)")
for (x in list(A, B)) cat(sprintf("%-44s OR %s [%s - %s], P = %s\n", x$label,
  num(unname(x$fisher$estimate)), num(x$fisher$conf.int[1]), num(x$fisher$conf.int[2]),
  pval(x$fisher$p.value)))

hdr("Neoplasia tier after removing mucosal high-risk (CIN) cases")
for (x in list(A, B)) cat(sprintf("%-44s %d CIN cases; neoplasia at floor %.2f%% (n=%d)\n",
  x$label, x$n_cin, x$neo_cut$pct, x$neo_cut$n))

cat("\nNote: Model 1 (positivity by sample class) and Model 3 (L1-gated productivity)\n",
    "do not read the tier table, so they are unchanged by this correction.\n", sep = "")
