#!/usr/bin/env Rscript
# Fits the three inference models behind the manuscript and assembles a single
# `stats` list holding every number quoted in the text/tables, so the .qmd can
# render them inline (`r stats$...`). Run standalone to print the summary, or
# source() it (e.g. from the manuscript setup chunk) to obtain `stats`.
# All paths are resolved via here::here(), so it works from any CWD.
suppressPackageStartupMessages({library(tidyverse); library(here)})
options(width=200)
wilson <- function(x,n){ if(n==0) return(c(NA,NA)); p<-x/n; z<-1.96
  ctr<-(p+z^2/(2*n))/(1+z^2/n); hw<-z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n); c(ctr-hw,ctr+hw)}
ortab <- function(m){ s<-summary(m)$coefficients
  data.frame(term=rownames(s), OR=exp(s[,1]),
             lo=exp(s[,1]-1.96*s[,2]), hi=exp(s[,1]+1.96*s[,2]), p=s[,4], row.names=NULL)}

cat("######################################################################\n")
cat("# MODEL 1 — Cell-line / engineered enrichment of HPV positivity (full_v2 unbiased sweep)\n")
cat("######################################################################\n")
st <- read_tsv(here("results_full_v2/metadata/hpv_status.tsv"), show_col_types=FALSE)
fl <- read_tsv(here("results_full_v2/metadata/cell_line_flags_curated.tsv"), show_col_types=FALSE) %>%
  transmute(srr_id, is_cell_line=toupper(is_cell_line)=="TRUE", is_engineered=toupper(is_engineered)=="TRUE")
d1 <- st %>% left_join(fl, by="srr_id") %>%
  mutate(is_cell_line=coalesce(is_cell_line,FALSE), is_engineered=coalesce(is_engineered,FALSE),
         hpv_pos = hpv_status=="HPV+",
         sample_class = factor(ifelse(is_cell_line,"cell_line",
                                ifelse(is_engineered,"engineered","clinical")),
                               levels=c("clinical","engineered","cell_line")))
cat(sprintf("Screened libraries: %d\n", nrow(d1)))
rates <- d1 %>% group_by(sample_class) %>%
  summarise(n=n(), pos=sum(hpv_pos), pct=round(100*mean(hpv_pos),2), .groups="drop") %>%
  rowwise() %>% mutate(ci_lo=round(100*wilson(pos,n)[1],2), ci_hi=round(100*wilson(pos,n)[2],2)) %>% ungroup()
print(rates)
m1 <- glm(hpv_pos ~ sample_class, binomial, d1)
o1 <- ortab(m1)
cat("\nLogistic regression  hpv_pos ~ sample_class (ref=clinical):\n"); print(o1, digits=3)
lr_p <- anova(m1, test="Chisq")$`Pr(>Chi)`[2]
cat(sprintf("LR test vs null: p = %.3g\n", lr_p))

cat("\n######################################################################\n")
cat("# MODEL 2 — Depth-controlled HPV detection across population tiers (pooled reframe)\n")
cat("######################################################################\n")
rp <- read_tsv(here("results_reframe_csc/hpv_rpm_per_sample.tsv"), show_col_types=FALSE) %>%
  filter(!is.na(total_reads), total_reads>0)
floor_rpm <- rp %>% filter(tier=="control_productive", !is_cell_line, rpm>0) %>% pull(rpm) %>% quantile(.05, na.rm=TRUE)
cat(sprintf("Control-calibrated detection floor (5th pctile productive controls): %.2f RPM\n", floor_rpm))
rp <- rp %>% mutate(any_hpv = rpm>0, detect = rpm>=floor_rpm, l10=log10(total_reads))
tiersum <- rp %>% group_by(tier) %>% summarise(n=n(), med_reads=round(median(total_reads)/1e6,1),
      pct_any=round(100*mean(any_hpv),1), pct_floor=round(100*mean(detect),1),
      med_rpm=round(median(rpm),3), .groups="drop")
cat("\nTier detection summary:\n"); print(tiersum)

## Headline contrast: neoplasia vs unselected skin, depth-adjusted
sub <- rp %>% filter(tier %in% c("unselected_skin","neoplasia")) %>%
  mutate(tier=relevel(factor(as.character(tier)), ref="unselected_skin"))
cat(sprintf("\n[neoplasia vs unselected_skin] n=%d\n", nrow(sub)))
cat("\nM2a: ANY HPV (rpm>0) ~ tier + log10(reads):\n")
m2a <- glm(any_hpv ~ tier + l10, binomial, sub); o2a <- ortab(m2a); print(o2a, digits=3)
cat("\nM2b: PRODUCTIVE-LEVEL detection (rpm>=floor) ~ tier + log10(reads):\n")
m2b <- glm(detect ~ tier + l10, binomial, sub); o2b <- ortab(m2b); print(o2b, digits=3)

## Productive controls vs unselected skin (near-separation -> Fisher exact)
cat("\nControls vs unselected_skin, productive-level detection (Fisher exact):\n")
ct <- rp %>% filter(tier %in% c("control_productive","unselected_skin")) %>%
  mutate(tier=factor(as.character(tier), levels=c("unselected_skin","control_productive")))
tab <- table(tier=ct$tier, detect=ct$detect); print(tab)
ft <- fisher.test(tab); cat(sprintf("OR=%.1f  95%%CI %.1f-%s  p=%.3g\n",
    ft$estimate, ft$conf.int[1], ifelse(is.finite(ft$conf.int[2]),sprintf('%.1f',ft$conf.int[2]),'Inf'), ft$p.value))

## Depth-only check within neoplasia: does deeper sequencing buy productive-level detection?
cat("\nWithin neoplasia: detect ~ log10(reads):\n")
mn <- glm(detect ~ l10, binomial, rp %>% filter(tier=="neoplasia")); mn_o <- ortab(mn); print(mn_o, digits=3)

cat("\n######################################################################\n")
cat("# MODEL 3 — Productive infection (L1>=3) among HPV+ typed samples\n")
cat("######################################################################\n")
get_prod <- function(dir){
  fs <- list.files(here(dir,"aggregated"), "_transcript_classes\\.tsv$", full.names=TRUE)
  map_dfr(fs, ~ suppressMessages(read_tsv(.x, show_col_types=FALSE, col_types=cols(.default="c")))) %>%
    filter(gene=="PRODUCTIVE_INFECTION") %>% transmute(srr_id=sample_id, productive=read_count=="yes")
}
ctrl <- get_prod("results_lesions") %>% mutate(grp="productive_controls (targeted warts/EV)")
csc  <- get_prod("results_lesions_csc")
ss   <- read_csv(here("results_lesions_csc/lesions_samplesheet_full.csv"), show_col_types=FALSE) %>%
  transmute(srr_id, diagnosis=tolower(coalesce(diagnosis,"")),
            is_cell_line=toupper(as.character(is_cell_line))=="TRUE")
csc2 <- csc %>% left_join(ss, by="srr_id") %>%
  mutate(grp=case_when(is_cell_line ~ "cell_line",
                       str_detect(diagnosis,"squamous|carcinoma|scc|keratos|bowen") ~ "neoplasia (cSCC/AK)",
                       TRUE ~ "other_lesion"))
prod <- bind_rows(ctrl %>% select(srr_id,productive,grp), csc2 %>% select(srr_id,productive,grp))
cat("Productive (L1>=3) rate among HPV+ typed samples, by group:\n")
psum <- prod %>% group_by(grp) %>% summarise(n=n(), prod=sum(productive), pct=round(100*mean(productive),1), .groups="drop")
print(psum)
## Fisher: productive controls vs neoplasia
pc <- prod %>% filter(grp %in% c("productive_controls (targeted warts/EV)","neoplasia (cSCC/AK)"))
f3 <- NULL
if(n_distinct(pc$grp)==2){
  tb <- table(grp=pc$grp, productive=pc$productive); print(tb)
  f3 <- fisher.test(tb); cat(sprintf("Fisher OR=%.2f 95%%CI %.2f-%.2f p=%.3g\n",
       f3$estimate, f3$conf.int[1], f3$conf.int[2], f3$p.value))
}

## ── Reference panel size (built on the HPC; .fai not committed) ──────────────
.fai <- here("assets/hpv_references/hpv_all.fasta.fai")
panel_n <- if (file.exists(.fai)) length(readLines(.fai)) else 455L

## ── HPV type-by-tier catalogue (for HPV16/HPV18 contamination fractions) ─────
tc <- read_tsv(here("results_reframe_csc/type_catalogue_by_tier.tsv"), show_col_types=FALSE)
.h18 <- tc[tc$hpv_reference=="HPV18REF",]; .h16 <- tc[tc$hpv_reference=="HPV16REF",]

## ── Peripheral discovery-cohort counts ──────────────────────────────────────
.n_targeted <- nrow(read_tsv(here("results_lesions/aggregated/hpv_status.tsv"),     show_col_types=FALSE))
.n_csc       <- nrow(read_tsv(here("results_lesions_csc/aggregated/hpv_status.tsv"), show_col_types=FALSE))

## ── Assemble the flat `stats` list quoted throughout the manuscript ──────────
.rc <- as.list(rates[rates$sample_class=="clinical",]);   .re <- as.list(rates[rates$sample_class=="engineered",]); .rl <- as.list(rates[rates$sample_class=="cell_line",])
.oe <- as.list(o1[o1$term=="sample_classengineered",]);   .ol <- as.list(o1[o1$term=="sample_classcell_line",])
.tc_ctrl<-as.list(tiersum[tiersum$tier=="control_productive",]); .tc_neo<-as.list(tiersum[tiersum$tier=="neoplasia",])
.tc_uns <-as.list(tiersum[tiersum$tier=="unselected_skin",]);    .tc_cl <-as.list(tiersum[tiersum$tier=="cell_line",]); .tc_oth<-as.list(tiersum[tiersum$tier=="other",])
.a_neo<-as.list(o2a[o2a$term=="tierneoplasia",]); .a_l10<-as.list(o2a[o2a$term=="l10",])
.b_neo<-as.list(o2b[o2b$term=="tierneoplasia",]); .b_l10<-as.list(o2b[o2b$term=="l10",])
.mn_l10<-as.list(mn_o[mn_o$term=="l10",])
.m3c<-as.list(psum[psum$grp=="productive_controls (targeted warts/EV)",]); .m3n<-as.list(psum[psum$grp=="neoplasia (cSCC/AK)",])

stats <- lapply(list(
  n_screened=nrow(d1), n_targeted=.n_targeted, n_csc=.n_csc, panel_n=panel_n,
  clin_pct=.rc$pct, clin_lo=.rc$ci_lo, clin_hi=.rc$ci_hi, clin_n=.rc$n, clin_pos=.rc$pos,
  eng_pct=.re$pct,  eng_lo=.re$ci_lo,  eng_hi=.re$ci_hi,  eng_n=.re$n,  eng_pos=.re$pos,
  cl_pct=.rl$pct,   cl_lo=.rl$ci_lo,   cl_hi=.rl$ci_hi,   cl_n=.rl$n,   cl_pos=.rl$pos,
  eng_or=.oe$OR, eng_or_lo=.oe$lo, eng_or_hi=.oe$hi, eng_or_p=.oe$p,
  cl_or=.ol$OR,  cl_or_lo=.ol$lo,  cl_or_hi=.ol$hi,  cl_or_p=.ol$p,
  lr_p=lr_p,
  hpv18_cl=.h18$cell_line, hpv18_tot=.h18$total, hpv16_cl=.h16$cell_line, hpv16_tot=.h16$total,
  floor_rpm=floor_rpm,
  ctrl_n=.tc_ctrl$n, ctrl_pctfloor=.tc_ctrl$pct_floor, ctrl_medrpm=.tc_ctrl$med_rpm, ctrl_medreads=.tc_ctrl$med_reads,
  neo_n=.tc_neo$n,   neo_pctfloor=.tc_neo$pct_floor,   neo_medrpm=.tc_neo$med_rpm,   neo_medreads=.tc_neo$med_reads,
  uns_n=.tc_uns$n,   uns_pctfloor=.tc_uns$pct_floor,   cl2_n=.tc_cl$n, oth_n=.tc_oth$n,
  m2a_or=.a_neo$OR, m2a_lo=.a_neo$lo, m2a_hi=.a_neo$hi, m2a_p=.a_neo$p,
  m2a_l10_or=.a_l10$OR, m2a_l10_lo=.a_l10$lo, m2a_l10_hi=.a_l10$hi, m2a_l10_p=.a_l10$p,
  m2b_or=.b_neo$OR, m2b_lo=.b_neo$lo, m2b_hi=.b_neo$hi, m2b_p=.b_neo$p,
  m2b_l10_or=.b_l10$OR, m2b_l10_lo=.b_l10$lo, m2b_l10_hi=.b_l10$hi, m2b_l10_p=.b_l10$p,
  fish_or=ft$estimate, fish_lo=ft$conf.int[1], fish_hi=ft$conf.int[2], fish_p=ft$p.value,
  neo_depth_or=.mn_l10$OR, neo_depth_lo=.mn_l10$lo, neo_depth_hi=.mn_l10$hi, neo_depth_p=.mn_l10$p,
  m3_ctrl_pct=.m3c$pct, m3_ctrl_prod=.m3c$prod, m3_ctrl_n=.m3c$n,
  m3_neo_pct=.m3n$pct,  m3_neo_prod=.m3n$prod,  m3_neo_n=.m3n$n,
  m3_or=f3$estimate, m3_lo=f3$conf.int[1], m3_hi=f3$conf.int[2], m3_p=f3$p.value
), unname)

saveRDS(stats, here("manuscript/stats.rds"))
cat("\nDONE — stats saved to manuscript/stats.rds\n")
