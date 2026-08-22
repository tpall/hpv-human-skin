#!/usr/bin/env Rscript
# Library chemistry as a covariate alongside sequencing depth.
#
# The manuscript models depth (log10 total reads) throughout, because depth
# bounds what a viral-transcript survey can detect. RNA-selection chemistry is
# the obvious counterpart and was not modelled: poly-A capture recovers
# polyadenylated early transcripts but is inefficient for the rest of the viral
# pool, so a poly-A library and an rRNA-depleted library of equal depth do not
# have equal power to see HPV.
#
# Inputs (built by bin/fetch_library_selection.py + bin/classify_library_chemistry.py):
#   results_reframe_fixed/library_chemistry.tsv
#
# Outputs (results_reframe_fixed/chemistry_*.tsv), for the supplement:
#   chemistry_by_tier.tsv        composition of each tier
#   chemistry_within_skin.tsv    detection by chemistry inside one tier
#   chemistry_models.tsv         headline contrasts with and without the term
#
# Run: Rscript bin/library_chemistry_covariate.R
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("manuscript/_paths.R"))
options(width = 200)

CHEM_LEVELS <- c("unspecified", "polyA", "rRNA_depleted", "random_primed",
                 "template_switch", "tenx_3prime")

ortab <- function(m, label) {
  s <- summary(m)$coefficients
  ref <- tryCatch(levels(m$model$chem)[1], error = function(e) NA_character_)
  tibble(model = label, chem_reference = ref, term = rownames(s), OR = exp(s[, 1]),
         lo = exp(s[, 1] - 1.96 * s[, 2]), hi = exp(s[, 1] + 1.96 * s[, 2]),
         p = s[, 4])
}

rp <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>%
  filter(!is.na(total_reads), total_reads > 0)
ch <- read_tsv(here("results_reframe_fixed/library_chemistry.tsv"),
               show_col_types = FALSE)

# Same control-calibrated floor as manuscript/results_models.R, recomputed here
# so this script stands alone.
floor_rpm <- rp %>% filter(tier == "control_productive", !is_cell_line, rpm > 0) %>%
  pull(rpm) %>% quantile(.05, na.rm = TRUE)

d <- rp %>%
  left_join(ch %>% select(srr_id, chemistry, evidence_source), by = "srr_id") %>%
  mutate(any_hpv = rpm > 0, detect = rpm >= floor_rpm, l10 = log10(total_reads),
         chem = factor(coalesce(chemistry, "unspecified"), levels = CHEM_LEVELS),
         chem_known = chem != "unspecified")

cat(sprintf("%d libraries; chemistry resolved for %d (%.1f%%)\n",
            nrow(d), sum(d$chem_known), 100 * mean(d$chem_known)))
cat(sprintf("detection floor: %.2f RPM\n\n", floor_rpm))

## 1. Chemistry composition by tier -----------------------------------------
by_tier <- d %>% count(tier, chem) %>%
  group_by(tier) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
cat("Chemistry composition by tier (%):\n")
print(by_tier %>% select(-n) %>% pivot_wider(names_from = chem, values_from = pct,
                                             values_fill = 0), n = 30)
write_tsv(by_tier, here("results_reframe_fixed/chemistry_by_tier.tsv"))

## 2. Does chemistry move detection at all? ---------------------------------
## Asked inside unselected skin, the one tier large and heterogeneous enough to
## carry every chemistry across many studies, so the comparison is not a single
## study wearing a chemistry label.
sk <- d %>% filter(tier == "unselected_skin")
within_skin <- sk %>% group_by(chem) %>%
  summarise(n = n(), n_studies = n_distinct(study),
            med_reads_M = round(median(total_reads) / 1e6, 1),
            pct_any = round(100 * mean(any_hpv), 2),
            pct_floor = round(100 * mean(detect), 2), .groups = "drop")
cat("\nWithin unselected skin, by chemistry:\n"); print(within_skin)
write_tsv(within_skin, here("results_reframe_fixed/chemistry_within_skin.tsv"))

m_sk <- glm(any_hpv ~ chem + l10, binomial, sk)
o_sk <- ortab(m_sk, "within_skin: any_hpv ~ chem + l10")
cat("\nDepth-adjusted, within skin:\n"); print(as.data.frame(o_sk), digits = 3)

## 3. The headline contrast, with and without the chemistry term -------------
sub <- d %>% filter(tier %in% c("unselected_skin", "neoplasia")) %>%
  mutate(tier = relevel(factor(as.character(tier)), ref = "unselected_skin"))

fits <- list(
  ortab(glm(any_hpv ~ tier + l10, binomial, sub), "M2a  any ~ tier+depth"),
  ortab(glm(any_hpv ~ tier + l10 + chem, binomial, sub), "M2a+ any ~ tier+depth+chem"),
  ortab(glm(detect ~ tier + l10, binomial, sub), "M2b  floor ~ tier+depth"),
  ## Productive-floor detections occur only in the poly-A and unspecified
  ## strata, so the six-level factor is completely separated here and its
  ## coefficients are not estimable. Collapse to poly-A vs everything else,
  ## which is the contrast the chemistry objection is actually about.
  ortab(glm(detect ~ tier + l10 + polya, binomial,
            sub %>% mutate(polya = chem == "polyA")),
        "M2b+ floor ~ tier+depth+polyA"))

## Sensitivity 1: only libraries whose chemistry is actually known, so the
## adjustment is not carried by a 38% "unspecified" stratum.
subk <- sub %>% filter(chem_known) %>% mutate(chem = droplevels(chem))
fits <- c(fits, list(
  ortab(glm(any_hpv ~ tier + l10 + chem, binomial, subk),
        "S1   any ~ tier+depth+chem (chemistry-known only)")))

## Sensitivity 2: restrict to rRNA-depleted libraries, the chemistry most
## favourable to viral recovery. If the neoplasia deficit is a poly-A artefact
## it should disappear here.
subr <- sub %>% filter(chem == "rRNA_depleted")
cat(sprintf("\nrRNA-depleted subset: %d libraries (%s)\n", nrow(subr),
            paste(sprintf("%s=%d", names(table(subr$tier)), table(subr$tier)),
                  collapse = ", ")))
if (n_distinct(subr$tier) == 2) {
  fits <- c(fits, list(ortab(glm(any_hpv ~ tier + l10, binomial, subr),
                             "S2   any ~ tier+depth (rRNA-depleted only)")))
  tb <- table(tier = droplevels(subr$tier), detect = subr$detect)
  cat("\nrRNA-depleted subset, productive-floor detection:\n"); print(tb)
}

## Sensitivity 3: drop 3'-barcoded (10x) libraries. They are not excluded from
## the tier analysis upstream, and they sit at both the highest depth and the
## highest trace rate, so they are the obvious lever on the any-read contrast.
subn <- sub %>% filter(chem != "tenx_3prime")
fits <- c(fits, list(ortab(glm(any_hpv ~ tier + l10, binomial, subn),
                           "S3   any ~ tier+depth (no 10x)")))

models <- bind_rows(fits)
## Where the productive-floor detections actually sit
cat("\nProductive-floor detections by tier and chemistry:\n")
print(sub %>% group_by(tier, chem) %>%
        summarise(n = n(), n_floor = sum(detect), .groups = "drop") %>%
        filter(n_floor > 0))

cat("\nHeadline contrast across specifications (tier term only):\n")
print(as.data.frame(models %>% filter(str_detect(term, "tierneoplasia"))), digits = 3)
cat("\nFull coefficients:\n"); print(as.data.frame(models), digits = 3)
write_tsv(models, here("results_reframe_fixed/chemistry_models.tsv"))

cat("\nWrote results_reframe_fixed/chemistry_{by_tier,within_skin,models}.tsv\n")
