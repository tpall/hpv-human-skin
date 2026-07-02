#!/usr/bin/env Rscript
# Bayesian robustness check of the three frequentist premises (see results_models.R).
# Same data, same model structure, but fitted with rstanarm under EXPLICIT
# weakly-informative priors (NOT flat / not the improper default):
#   coefficients (log-odds scale):  Normal(0, 2.5)   [Gelman et al. 2008 weakly-informative]
#   intercept:                      Normal(0, 5)
# autoscale=FALSE so the priors act directly on the interpretable log-odds scale
# (this is what regularises the near-separated controls-vs-skin contrast).
# Continuous depth covariate is mean-centred (l10c) so the tier effect is read
# at mean sequencing depth and its prior is on a natural per-10x-reads scale.
suppressPackageStartupMessages({library(tidyverse); library(rstanarm); library(here)})
options(mc.cores = 4)
SEED <- 123
wk_coef <- normal(0, 2.5, autoscale = FALSE)   # weakly-informative slope prior
wk_int  <- normal(0, 5,   autoscale = FALSE)   # weakly-informative intercept prior

# Posterior OR summary for one coefficient: median, 95% credible interval, and
# posterior probability the effect is in the stated direction.
or_summ <- function(fit, term){
  b  <- as.matrix(fit)[, term]
  or <- exp(b)
  list(OR = median(or), lo = unname(quantile(or, .025)), hi = unname(quantile(or, .975)),
       p_gt1 = mean(b > 0), rhat = max(summary(fit)[, "Rhat"], na.rm = TRUE))
}
fo  <- function(x) if (x >= 1000) formatC(round(x), big.mark = ",", format = "d") else
                   if (x >= 10) formatC(round(x), format = "d") else formatC(x, format = "f", digits = 1)
line <- function(tag, freq, s)
  cat(sprintf("  %-26s freq OR %-7s | Bayes OR %s [%s, %s]  P(OR>1)=%.3f  (Rhat %.3f)\n",
              tag, freq, fo(s$OR), fo(s$lo), fo(s$hi), s$p_gt1, s$rhat))

## ── Data (mirrors results_models.R) ─────────────────────────────────────────
st <- read_tsv(here("results_full_v2/metadata/hpv_status.tsv"), show_col_types = FALSE)
fl <- read_tsv(here("results_full_v2/metadata/cell_line_flags_curated.tsv"), show_col_types = FALSE) %>%
  transmute(srr_id, is_cell_line = toupper(is_cell_line) == "TRUE", is_engineered = toupper(is_engineered) == "TRUE")
d1 <- st %>% left_join(fl, by = "srr_id") %>%
  mutate(is_cell_line = coalesce(is_cell_line, FALSE), is_engineered = coalesce(is_engineered, FALSE),
         hpv_pos = hpv_status == "HPV+",
         sample_class = factor(ifelse(is_cell_line, "cell_line", ifelse(is_engineered, "engineered", "clinical")),
                               levels = c("clinical", "engineered", "cell_line")))
rp <- read_tsv(here("results_reframe_csc/hpv_rpm_per_sample.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(total_reads), total_reads > 0)
floor_rpm <- rp %>% filter(tier == "control_productive", !is_cell_line, rpm > 0) %>% pull(rpm) %>% quantile(.05, na.rm = TRUE)
rp <- rp %>% mutate(any_hpv = rpm > 0, detect = rpm >= floor_rpm, l10 = log10(total_reads))

get_prod <- function(dir){
  fs <- list.files(here(dir, "aggregated"), "_transcript_classes\\.tsv$", full.names = TRUE)
  map_dfr(fs, ~ suppressMessages(read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))) %>%
    filter(gene == "PRODUCTIVE_INFECTION") %>% transmute(srr_id = sample_id, productive = read_count == "yes")
}
ss <- read_csv(here("results_lesions_csc/lesions_samplesheet_full.csv"), show_col_types = FALSE) %>%
  transmute(srr_id, diagnosis = tolower(coalesce(diagnosis, "")), is_cell_line = toupper(as.character(is_cell_line)) == "TRUE")
prod <- bind_rows(
  get_prod("results_lesions") %>% mutate(grp = "productive_controls"),
  get_prod("results_lesions_csc") %>% left_join(ss, by = "srr_id") %>%
    mutate(grp = case_when(is_cell_line ~ "cell_line",
                           str_detect(diagnosis, "squamous|carcinoma|scc|keratos|bowen") ~ "neoplasia",
                           TRUE ~ "other_lesion"))) %>%
  select(srr_id, productive, grp)

## ── Fits ────────────────────────────────────────────────────────────────────
message("Fitting Bayesian models (rstanarm, weakly-informative priors)…")
b1 <- stan_glm(hpv_pos ~ sample_class, binomial, d1, prior = wk_coef, prior_intercept = wk_int, seed = SEED, refresh = 0)

sub <- rp %>% filter(tier %in% c("unselected_skin", "neoplasia")) %>%
  mutate(tier = relevel(factor(as.character(tier)), ref = "unselected_skin"), l10c = l10 - mean(l10))
b2a <- stan_glm(any_hpv ~ tier + l10c, binomial, sub, prior = wk_coef, prior_intercept = wk_int, seed = SEED, refresh = 0)
b2b <- stan_glm(detect  ~ tier + l10c, binomial, sub, prior = wk_coef, prior_intercept = wk_int, seed = SEED, refresh = 0)

ctl <- rp %>% filter(tier %in% c("control_productive", "unselected_skin")) %>%
  mutate(tier = factor(as.character(tier), levels = c("unselected_skin", "control_productive")))
b2c <- stan_glm(detect ~ tier, binomial, ctl, prior = wk_coef, prior_intercept = wk_int, seed = SEED, refresh = 0)
# heavier-tailed sensitivity for the near-separated contrast
b2c_t <- stan_glm(detect ~ tier, binomial, ctl, prior = student_t(3, 0, 2.5, autoscale = FALSE),
                  prior_intercept = wk_int, seed = SEED, refresh = 0)

neo <- rp %>% filter(tier == "neoplasia") %>% mutate(l10c = l10 - mean(l10))
b2d <- stan_glm(detect ~ l10c, binomial, neo, prior = wk_coef, prior_intercept = wk_int, seed = SEED, refresh = 0)

pc <- prod %>% filter(grp %in% c("neoplasia", "productive_controls")) %>%
  mutate(grp = factor(grp, levels = c("neoplasia", "productive_controls")))
b3 <- stan_glm(productive ~ grp, binomial, pc, prior = wk_coef, prior_intercept = wk_int, seed = SEED, refresh = 0)

## ── Report ──────────────────────────────────────────────────────────────────
cat("\n================================================================\n")
cat("BAYESIAN PREMISE CHECKS  (weakly-informative Normal(0,2.5) coef priors)\n")
cat("================================================================\n")

cat("\nPREMISE 1 — Cell-line / engineered libraries are enriched for HPV positivity:\n")
line("engineered vs clinical", "5.4", or_summ(b1, "sample_classengineered"))
line("cell_line vs clinical",  "4.4", or_summ(b1, "sample_classcell_line"))

cat("\nPREMISE 2 — After depth adjustment, neoplasia is NOT elevated at productive level,\n")
cat("            but productive controls hugely exceed skin (the trace-not-productive story):\n")
line("neoplasia vs skin (floor)", "2.0", or_summ(b2b, "tierneoplasia"))
line("depth /10x reads (floor)",  "-",   or_summ(b2b, "l10c"))
line("controls vs skin (floor)",  "2,526", or_summ(b2c, "tiercontrol_productive"))
line("  ^ student-t(3) prior",    "2,526", or_summ(b2c_t, "tiercontrol_productive"))
line("depth within neoplasia",    "1.7", or_summ(b2d, "l10c"))

cat("\nPREMISE 2b — Neoplasia carries a modest TRACE (any-HPV) excess:\n")
line("neoplasia vs skin (any)",  "1.7", or_summ(b2a, "tierneoplasia"))
line("depth /10x reads (any)",   "4.6", or_summ(b2a, "l10c"))

cat("\nPREMISE 3 — Productive infection is more frequent in warts/EV than in neoplasia:\n")
line("controls vs neoplasia",    "3.8", or_summ(b3, "grpproductive_controls"))

## ── Persist posterior OR summaries (so a robustness paragraph can cite them) ──
grab <- function(fit, term, tag) tibble::tibble(model = tag, !!!or_summ(fit, term))
bayes_summary <- bind_rows(
  grab(b1, "sample_classengineered", "M1_engineered"),
  grab(b1, "sample_classcell_line",  "M1_cell_line"),
  grab(b2b, "tierneoplasia",          "M2b_neoplasia_vs_skin_floor"),
  grab(b2b, "l10c",                   "M2b_depth"),
  grab(b2c, "tiercontrol_productive", "M2c_controls_vs_skin_floor"),
  grab(b2c_t, "tiercontrol_productive","M2c_controls_vs_skin_floor_studentt"),
  grab(b2d, "l10c",                   "M2d_depth_within_neoplasia"),
  grab(b2a, "tierneoplasia",          "M2a_neoplasia_vs_skin_any"),
  grab(b3, "grpproductive_controls",  "M3_controls_vs_neoplasia"))
saveRDS(bayes_summary, here("manuscript/bayes_summary.rds"))
readr::write_tsv(bayes_summary %>% mutate(across(where(is.numeric), ~round(.x, 4))),
                 here("manuscript/bayes_summary.tsv"))
cat("\nDONE — posterior summaries saved to manuscript/bayes_summary.{rds,tsv}\n")
