#!/usr/bin/env Rscript
# Fits the three inference models behind the manuscript and assembles a single
# `stats` list holding every number quoted in the text/tables, so the .qmd can
# render them inline (`r stats$...`). Run standalone to print the summary, or
# source() it (e.g. from the manuscript setup chunk) to obtain `stats`.
# All paths are resolved via here::here(), so it works from any CWD.
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("manuscript/_paths.R"))
options(width=200)
wilson <- function(x,n){ if(n==0) return(c(NA,NA)); p<-x/n; z<-1.96
  ctr<-(p+z^2/(2*n))/(1+z^2/n); hw<-z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n); c(ctr-hw,ctr+hw)}
comma_ <- function(x) formatC(round(x), format="d", big.mark=",")
ortab <- function(m){ s<-summary(m)$coefficients
  data.frame(term=rownames(s), OR=exp(s[,1]),
             lo=exp(s[,1]-1.96*s[,2]), hi=exp(s[,1]+1.96*s[,2]), p=s[,4], row.names=NULL)}

cat("######################################################################\n")
cat("# MODEL 1: Cell-line / engineered enrichment of HPV positivity (full_v2 unbiased sweep)\n")
cat("######################################################################\n")
st <- read_tsv(here("results_full_v2/metadata/hpv_status.tsv"), show_col_types=FALSE)
fl <- read_tsv(here("results_full_v2/metadata/cell_line_flags_curated.tsv"), show_col_types=FALSE) %>%
  transmute(srr_id, is_cell_line=toupper(is_cell_line)=="TRUE", is_engineered=toupper(is_engineered)=="TRUE") %>%
  apply_cellline_overrides()
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
cat("# MODEL 2: Depth-controlled HPV detection across population tiers (pooled reframe)\n")
cat("######################################################################\n")
rp <- read_tsv(rpm_per_sample_tsv(), show_col_types=FALSE) %>%
  filter(!is.na(total_reads), total_reads>0)
floor_rpm <- rp %>% filter(tier=="control_productive", !is_cell_line, rpm>0) %>% pull(rpm) %>% quantile(.05, na.rm=TRUE)
cat(sprintf("Control-calibrated detection floor (5th pctile productive controls): %.2f RPM\n", floor_rpm))
rp <- rp %>% mutate(any_hpv = rpm>0, detect = rpm>=floor_rpm, l10=log10(total_reads))
tiersum <- rp %>% group_by(tier) %>% summarise(n=n(), med_reads=round(median(total_reads)/1e6,1),
      pct_any=round(100*mean(any_hpv),1), pct_floor=round(100*mean(detect),1),
      med_rpm=round(median(rpm),3), n_floor=sum(detect),
      max_reads=round(max(total_reads)/1e6,1), .groups="drop")
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
ft <- fisher.test(tab)
# Under near-complete separation `fisher.test` returns its internal search
# ceiling (~4.5e15) as the upper bound. That is an artefact of the optimiser,
# not an estimate, so report the interval as unbounded above.
fish_hi <- if (ft$conf.int[2] > 1e12) Inf else ft$conf.int[2]
cat(sprintf("OR=%.1f  95%%CI %.1f-%s  p=%.3g\n",
    ft$estimate, ft$conf.int[1], ifelse(is.finite(fish_hi),sprintf('%.1f',fish_hi),'Inf'), ft$p.value))

## EV and peri-lesional mucosa: separated against skin, so Fisher exact rather
## than a depth-adjusted glm (every peri-lesional library is above the floor).
fisher_vs_skin <- function(other){
  d <- rp %>% filter(tier %in% c("unselected_skin", other))
  if (n_distinct(d$tier) < 2) return(list(OR=NA, lo=NA, hi=NA, p=NA))
  f <- fisher.test(table(tier=factor(as.character(d$tier), levels=c("unselected_skin", other)),
                         detect=d$detect))
  list(OR=unname(f$estimate), lo=f$conf.int[1],
       hi=if (f$conf.int[2] > 1e12) Inf else f$conf.int[2], p=f$p.value)
}
ev_f  <- fisher_vs_skin("EV")
per_f <- fisher_vs_skin("peri_lesional")
cat(sprintf("\nEV vs unselected_skin (Fisher): OR=%.1f p=%.3g\n", ev_f$OR, ev_f$p))
cat(sprintf("Peri-lesional vs unselected_skin (Fisher): OR=%s p=%.3g\n",
            ifelse(is.finite(per_f$OR), sprintf('%.1f', per_f$OR), 'Inf'), per_f$p))

## Depth-only check within neoplasia: does deeper sequencing buy productive-level detection?
cat("\nWithin neoplasia: detect ~ log10(reads):\n")
mn <- glm(detect ~ l10, binomial, rp %>% filter(tier=="neoplasia")); mn_o <- ortab(mn); print(mn_o, digits=3)

cat("\n######################################################################\n")
cat("# MODEL 2b-split: Actinic keratosis vs cSCC within the neoplasia tier\n")
cat("######################################################################\n")
## The neoplasia tier pools AK with cSCC. That pooling has to be justified
## rather than assumed, because a published hit-and-run mechanism for HPV8 is
## built specifically on actinic keratosis. Split on the curated diagnosis text:
## libraries naming both AK and carcinoma are held out as ambiguous rather than
## assigned to either arm.
.dxsplit <- function(d) dplyr::case_when(
  d$tier != "neoplasia" ~ d$tier,
  stringr::str_detect(d$diagnosis, "actinic") &
    stringr::str_detect(d$diagnosis, "scc|carcinoma") ~ "neo_mixed",
  stringr::str_detect(d$diagnosis, "actinic") ~ "neo_AK",
  TRUE ~ "neo_cSCC")
rps <- rp %>% mutate(dx = .dxsplit(.))
.akc <- rps %>% filter(dx %in% c("neo_AK","neo_cSCC","neo_mixed")) %>% group_by(dx) %>%
  summarise(n=n(), n_studies=n_distinct(study), med_reads=median(total_reads)/1e6,
            pct_any=100*mean(any_hpv), n_floor=sum(detect), pct_floor=100*mean(detect),
            med_rpm=median(rpm), .groups="drop")
print(as.data.frame(.akc), digits=3)

## Each arm against the unselected-skin null, and then against each other.
.vs_skin <- function(arm, resp) {
  d <- rps %>% filter(dx %in% c("unselected_skin", arm)) %>%
    mutate(g = relevel(factor(dx), ref="unselected_skin"))
  as.list(ortab(glm(as.formula(paste(resp, "~ g + l10")), binomial, d)) %>%
          filter(term == paste0("g", arm)))
}
.ak_any <- .vs_skin("neo_AK","any_hpv");   .ak_fl <- .vs_skin("neo_AK","detect")
.sc_any <- .vs_skin("neo_cSCC","any_hpv"); .sc_fl <- .vs_skin("neo_cSCC","detect")
cat(sprintf("\nvs unselected skin: AK any OR %.2f (p=%.3g), floor OR %.1f (p=%.3g); cSCC any OR %.2f (p=%.3g), floor OR %.1f (p=%.3g)\n",
    .ak_any$OR, .ak_any$p, .ak_fl$OR, .ak_fl$p, .sc_any$OR, .sc_any$p, .sc_fl$OR, .sc_fl$p))

.ac <- rps %>% filter(dx %in% c("neo_AK","neo_cSCC")) %>%
  mutate(g = relevel(factor(dx), ref="neo_cSCC"))
.akvs <- as.list(ortab(glm(any_hpv ~ g + l10, binomial, .ac)) %>% filter(term=="gneo_AK"))
.akvf <- fisher.test(table(.ac$dx, .ac$detect))
cat(sprintf("AK vs cSCC directly: any OR %.2f (95%% CI %.2f-%.2f, p=%.3g); floor Fisher OR %.2f p=%.3g\n",
    .akvs$OR, .akvs$lo, .akvs$hi, .akvs$p, .akvf$estimate, .akvf$p.value))

## HPV8 specifically: the type the AK hit-and-run mechanism is built on. It is
## in the reference panel, so its absence is a result rather than a blind spot.
.tt_pre <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types=FALSE)
.h8_files <- unlist(lapply(c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide"),
  function(d) list.files(here(d,"aggregated"), "_hpv_types[.]tsv$", full.names=TRUE)))
.h8_all <- map_dfr(.h8_files, ~ suppressMessages(read_tsv(.x, show_col_types=FALSE, col_types=cols(.default="c"))))
.h8_typed <- sum(.h8_all$hpv_reference %in% c("HPV8REF","HPV8"))
.h8_trace <- sum(.tt_pre$top_ref %in% c("HPV8REF","HPV8"), na.rm=TRUE)
cat(sprintf("HPV8: %d assigned calls, %d trace best-hits across %d trace libraries\n",
            .h8_typed, .h8_trace, nrow(.tt_pre)))

cat("\n######################################################################\n")
cat("# MODEL 2c: Library chemistry alongside depth\n")
cat("######################################################################\n")
## Depth is modelled throughout because it bounds what a viral-transcript survey
## can see. RNA-selection chemistry is the counterpart: poly-A capture recovers
## polyadenylated early transcripts but is inefficient for the rest of the viral
## pool, so two libraries of equal depth need not have equal power to detect HPV.
## Chemistry is called from SRA `LibrarySelection` where it carries a real value
## and from the free-text construction protocol otherwise
## (bin/fetch_library_selection.py, bin/classify_library_chemistry.py).
.CHEMLV <- c("unspecified","polyA","rRNA_depleted","random_primed",
             "template_switch","tenx_3prime")
.chem <- tryCatch(read_tsv(here("results_reframe_fixed/library_chemistry.tsv"),
                           show_col_types=FALSE), error=function(e) NULL)
if (is.null(.chem)) {
  cat("library_chemistry.tsv absent; chemistry covariate skipped\n")
  .cq <- function(...) NA_real_
  chem_tbl <- tibble(Tier = character())   # @supptbl-chem renders only if filled
} else {
  rpc <- rp %>% left_join(.chem %>% select(srr_id, chemistry, library_selection), by="srr_id") %>%
    mutate(chem=factor(coalesce(chemistry,"unspecified"), levels=.CHEMLV),
           chem_known = chem != "unspecified")
  cat(sprintf("Chemistry resolved for %d/%d libraries (%.1f%%)\n",
      sum(rpc$chem_known), nrow(rpc), 100*mean(rpc$chem_known)))
  .chempct <- function(tr, cl) 100*mean(rpc$chem[rpc$tier==tr]==cl)
  cat("\nChemistry composition by tier (%):\n")
  print(round(100*prop.table(table(rpc$tier, rpc$chem),1),1))

  ## Chemistry effect on trace detection. Fitted among libraries whose chemistry
  ## is known, so the estimate is not carried by the unspecified stratum, and
  ## within the two tiers the headline contrast is drawn from, with a tier term:
  ## chemistry mix differs sharply across tiers (cell lines are half poly-A, EV
  ## is entirely random-primed), so an unstratified fit confounds the two.
  ## Reference is poly-A, since the objection is about poly-A vs rRNA depletion.
  .ck <- rpc %>% filter(chem_known, tier %in% c("unselected_skin","neoplasia")) %>%
    mutate(chem=droplevels(chem),
           tier=relevel(factor(as.character(tier)), ref="unselected_skin"))
  .ock <- ortab(glm(any_hpv ~ chem + l10 + tier, binomial, .ck))
  cat("\nany_hpv ~ chem + l10 + tier (skin+neoplasia, chemistry-known, ref = poly-A):\n")
  print(.ock, digits=3)
  ## Same term across the whole cohort, tier-adjusted, as a check that the null
  ## is not an artefact of restricting to two tiers. It is modestly above one
  ## there, carried by the non-cutaneous tiers.
  .ockall <- ortab(glm(any_hpv ~ chem + l10 + tier, binomial,
                       rpc %>% filter(chem_known) %>% mutate(chem=droplevels(chem))))
  cat("\nSame, all tiers:\n"); print(.ockall %>% filter(grepl("chem", term)), digits=3)

  ## Headline contrast re-fitted with the chemistry term.
  .subc2 <- rpc %>% filter(tier %in% c("unselected_skin","neoplasia")) %>%
    mutate(tier=relevel(factor(as.character(tier)), ref="unselected_skin"))
  .a2 <- as.list(ortab(glm(any_hpv ~ tier + l10 + chem, binomial, .subc2)) %>%
                 filter(term=="tierneoplasia"))
  ## Productive-floor detections occur only in the poly-A and unspecified strata,
  ## so the six-level factor is completely separated in the floor model. Collapse
  ## to poly-A vs everything else, which is the contrast at issue.
  .b2 <- as.list(ortab(glm(detect ~ tier + l10 + polya, binomial,
                           .subc2 %>% mutate(polya = chem=="polyA"))) %>%
                 filter(term=="tierneoplasia"))
  cat(sprintf("\nTier term with chemistry: any-read OR %.2f (p=%.3g); floor OR %.2f (p=%.3g)\n",
      .a2$OR, .a2$p, .b2$OR, .b2$p))

  ## The sensitivity that settles it: restrict to rRNA-depleted libraries, the
  ## chemistry most favourable to viral recovery. If the neoplasia deficit were a
  ## poly-A artefact it would disappear here.
  .rib <- .subc2 %>% filter(chem=="rRNA_depleted")
  .s2 <- as.list(ortab(glm(any_hpv ~ tier + l10, binomial, .rib)) %>%
                 filter(term=="tierneoplasia"))
  .ribneo <- .rib %>% filter(tier=="neoplasia")
  cat(sprintf("rRNA-depleted subset: %d libraries (%d skin, %d neoplasia from %d studies, median %.1fM, max %.1fM); %d above floor; tier OR %.2f (p=%.3g)\n",
      nrow(.rib), sum(.rib$tier=="unselected_skin"), nrow(.ribneo),
      n_distinct(.ribneo$study), median(.ribneo$total_reads)/1e6,
      max(.ribneo$total_reads)/1e6, sum(.rib$detect), .s2$OR, .s2$p))

  ## Supplementary table: chemistry composition by tier, in the display order and
  ## tier labels used elsewhere in the manuscript.
  .chemlab <- c(unspecified="Unspecified", polyA="Poly-A", rRNA_depleted="rRNA-depleted",
                random_primed="Random-primed", template_switch="Template-switching",
                tenx_3prime="3'-barcoded")
  .tierlab <- c(unselected_skin="Unselected skin", neoplasia="Neoplasia (cSCC/AK)",
                control_productive="Productive controls", EV="EV biopsies",
                peri_lesional="Peri-lesional mucosa", cell_line="Cell lines",
                engineered="Engineered", non_permissive_tissue="Non-permissive tissue",
                other="Other")
  chem_tbl <- rpc %>% count(tier, chem) %>% group_by(tier) %>%
    mutate(pct = 100*n/sum(n), N = sum(n)) %>% ungroup() %>%
    mutate(chem = .chemlab[as.character(chem)]) %>%
    select(-n) %>% pivot_wider(names_from=chem, values_from=pct, values_fill=0) %>%
    mutate(Tier = paste0(coalesce(.tierlab[tier], tier), " (", comma_(N), ")")) %>%
    arrange(match(tier, names(.tierlab))) %>%
    select(Tier, all_of(unname(.chemlab))) %>%
    mutate(across(-Tier, ~ sub("\\.0$", "", formatC(.x, format="f", digits=1))))

  .cqv <- list(
    chem_n=sum(rpc$chem_known), chem_pct=100*mean(rpc$chem_known),
    chem_total=nrow(rpc),
    # Share filed under the generic "cDNA": the reason the free-text protocol
    # pass is needed at all.
    chem_cdna_pct=100*mean(coalesce(rpc$library_selection,"")=="cDNA"),
    chem_neo_ribo=.chempct("neoplasia","rRNA_depleted"),
    chem_skin_ribo=.chempct("unselected_skin","rRNA_depleted"),
    chem_ctrl_ribo=.chempct("control_productive","rRNA_depleted"),
    chem_neo_polya=.chempct("neoplasia","polyA"),
    chem_skin_polya=.chempct("unselected_skin","polyA"),
    chem_ribo_or=.ock$OR[.ock$term=="chemrRNA_depleted"],
    chem_ribo_lo=.ock$lo[.ock$term=="chemrRNA_depleted"],
    chem_ribo_hi=.ock$hi[.ock$term=="chemrRNA_depleted"],
    chem_ribo_p =.ock$p [.ock$term=="chemrRNA_depleted"],
    chem_tenx_or=.ock$OR[.ock$term=="chemtenx_3prime"],
    chem_smart_or=.ock$OR[.ock$term=="chemtemplate_switch"],
    chem_rand_or=.ock$OR[.ock$term=="chemrandom_primed"],
    chem_ribo_or_all=.ockall$OR[.ockall$term=="chemrRNA_depleted"],
    chem_ribo_lo_all=.ockall$lo[.ockall$term=="chemrRNA_depleted"],
    chem_ribo_hi_all=.ockall$hi[.ockall$term=="chemrRNA_depleted"],
    m2a_or_chem=.a2$OR, m2a_p_chem=.a2$p,
    m2b_or_chem=.b2$OR, m2b_p_chem=.b2$p,
    ribo_n=nrow(.rib), ribo_skin_n=sum(.rib$tier=="unselected_skin"),
    ribo_neo_n=nrow(.ribneo), ribo_neo_studies=n_distinct(.ribneo$study),
    ribo_neo_medreads=median(.ribneo$total_reads)/1e6,
    ribo_neo_maxreads=max(.ribneo$total_reads)/1e6,
    ribo_floor=sum(.rib$detect),
    ribo_or=.s2$OR, ribo_lo=.s2$lo, ribo_hi=.s2$hi, ribo_p=.s2$p)
  .cq <- function(k) .cqv[[k]]
}

cat("\n######################################################################\n")
cat("# MODEL 3: Productive infection (L1>=3) among HPV+ typed samples\n")
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

## ── CIN data-hygiene sensitivity (Table S4) ─────────────────────────────────
## The broadened SCC/CIN query captured a few high-risk *mucosal* alpha (HPV16/18..)
## cases into the cutaneous neoplasia tier. Flag them by dominant assigned type and
## show that removing them only sharpens the trace-not-productive conclusion.
.HRTYPES <- c(16,18,31,33,35,39,45,51,52,56,58,59,68,26,53,66,67,69,70,73,82)
.domtype <- function(dirs){
  fs <- unlist(lapply(dirs, function(d) list.files(here(d,"aggregated"),"_hpv_types\\.tsv$",full.names=TRUE)))
  map_dfr(fs, ~ suppressMessages(read_tsv(.x, show_col_types=FALSE, col_types=cols(.default="c")))) %>%
    mutate(coverage_breadth=as.numeric(coverage_breadth), mean_depth=as.numeric(mean_depth)) %>%
    filter(coverage_breadth>=0.10, mean_depth>=2) %>%
    group_by(sample_id) %>% slice_max(mean_depth, n=1, with_ties=FALSE) %>% ungroup() %>%
    transmute(srr_id=sample_id, dom_type=hpv_reference,
              mucosal_hr = suppressWarnings(as.integer(str_match(dom_type,"HPV[-]?(\\d+)")[,2])) %in% .HRTYPES)
}
.dt  <- .domtype(c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide"))
.rp2 <- rp %>% left_join(.dt, by="srr_id") %>% mutate(mucosal_hr=coalesce(mucosal_hr,FALSE))
.cin <- .rp2 %>% filter(tier=="neoplasia", mucosal_hr) %>% arrange(desc(rpm))
# Estonian tissue categories -> English for display only (pipeline keeps Estonian).
# "muu" is the fall-through of the keyword classifier, i.e. no site keyword was
# found in the SRA text -- it means "site not recorded", not "a different site".
# Labelling it "Other" read as a positive finding about the tissue.
.tissue_en <- c(nahk="Skin", anogenitaal="Anogenital", suuoos="Oral", muu="Unclassified")
cin_tbl <- .cin %>% transmute(Accession=srr_id, Study=study, Type=sub("REF$","",dom_type),
                              Tissue=coalesce(.tissue_en[tissue_category], str_to_sentence(tissue_category)),
                              Diagnosis=diagnosis, RPM=round(rpm,1))
.rpc  <- .rp2 %>% filter(!(tier=="neoplasia" & mucosal_hr))
.neo_floor_cut <- .rpc %>% filter(tier=="neoplasia") %>% summarise(p=100*mean(detect)) %>% pull(p)
.subc <- .rpc %>% filter(tier %in% c("unselected_skin","neoplasia")) %>%
  mutate(tier=relevel(factor(as.character(tier)),"unselected_skin"))
.bcut <- as.list(ortab(glm(detect ~ tier + l10, binomial, .subc)) %>% filter(term=="tierneoplasia"))
.prodc <- prod %>% left_join(.dt %>% select(srr_id,mucosal_hr), by="srr_id") %>%
  mutate(mucosal_hr=coalesce(mucosal_hr,FALSE)) %>% filter(!(grp=="neoplasia (cSCC/AK)" & mucosal_hr))
.pcc  <- .prodc %>% filter(grp %in% c("productive_controls (targeted warts/EV)","neoplasia (cSCC/AK)"))
.m3nc <- .pcc %>% filter(grp=="neoplasia (cSCC/AK)") %>% summarise(n=n(), prod=sum(productive), pct=100*mean(productive))
.f3c  <- fisher.test(table(grp=.pcc$grp, productive=.pcc$productive))
cat(sprintf("\nCIN check: %d mucosal-HR neoplasia cases (%s); neo floor -> %.2f%%, M2b OR %.2f p=%.3g; M3 neo %.1f%% (%d/%d), OR %.2f p=%.4g\n",
    nrow(.cin), paste(.cin$srr_id,collapse=","), .neo_floor_cut, .bcut$OR, .bcut$p,
    .m3nc$pct, .m3nc$prod, .m3nc$n, .f3c$estimate, .f3c$p.value))

## ── Trace-type detection catalogue summary (Table S5) ───────────────────────
.tt <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types=FALSE) %>%
  mutate(trace_ok = !top_ref %in% c("NONE","NO_BAM") & top_reads>=5 & max_mapq>=20 &
                    top_reads/pmax(n_primary,1)>=0.5)
.tgen <- read_csv(here("manuscript/type_genus.csv"), show_col_types=FALSE) %>% select(hpv_reference, genus_short)
.tres <- .tt %>% filter(trace_ok) %>% left_join(rp %>% select(srr_id,tier), by="srr_id") %>%
  left_join(.tgen, by=c("top_ref"="hpv_reference"))
.tskin <- .tres %>% filter(tier=="unselected_skin")
cat(sprintf("\nTrace catalogue: %d/%d HPV+ untyped re-mapped to a confident type; %d distinct types; skin cutaneous(beta/gamma/nu) %d/%d\n",
    sum(.tt$trace_ok), nrow(.tt), n_distinct(.tres$top_ref),
    sum(.tskin$genus_short %in% c("beta","gamma","nu")), nrow(.tskin)))

## ── Reference panel size (built on the HPC; .fai not committed) ──────────────
.fai <- here("assets/hpv_references/hpv_all.fasta.fai")
panel_n <- if (file.exists(.fai)) length(readLines(.fai)) else 455L

## ── HPV type-by-tier catalogue (for HPV16/HPV18 contamination fractions) ─────
tc <- read_tsv(type_catalogue_tsv(), show_col_types=FALSE)
.h18 <- tc[tc$hpv_reference=="HPV18REF",]; .h16 <- tc[tc$hpv_reference=="HPV16REF",]

## ── Peripheral discovery-cohort counts ──────────────────────────────────────
.n_targeted <- nrow(read_tsv(here("results_lesions/aggregated/hpv_status.tsv"),     show_col_types=FALSE))
.n_csc       <- nrow(read_tsv(here("results_lesions_csc/aggregated/hpv_status.tsv"), show_col_types=FALSE))

## ── Library completeness: were whole libraries screened? ────────────────────
## Reading only the leading N reads of a library is a known way to miss viral
## signal, because a run submitted from a coordinate-sorted BAM emerges from the
## SRA in reference order and its unmapped (hence viral) reads sit at the tail.
## Nothing here subsamples, but the claim should be checkable rather than
## asserted, so compare reads actually screened against the SRA spot count.
.spots <- tryCatch(read_tsv(here("results_reframe_fixed/library_selection.tsv"),
                            show_col_types=FALSE) %>%
                     transmute(srr_id, spots=suppressWarnings(as.numeric(spots))),
                   error=function(e) NULL)
if (is.null(.spots)) {
  .sr <- function(k) NA_real_
} else {
  .srd <- rp %>% inner_join(.spots, by="srr_id") %>%
    filter(!is.na(spots), spots > 0) %>% mutate(ratio = total_reads/spots)
  .srtier <- .srd %>% group_by(tier) %>% summarise(m=median(ratio), .groups="drop")
  .srv <- list(n=nrow(.srd), med=100*median(.srd$ratio),
               p05=100*quantile(.srd$ratio, .05),
               tier_lo=100*min(.srtier$m), tier_hi=100*max(.srtier$m))
  .sr <- function(k) .srv[[k]]
  cat(sprintf("\nLibrary completeness: %d libraries, median %.1f%% of SRA spots screened (5th pctile %.1f%%); tier medians %.1f-%.1f%%\n",
      .srv$n, .srv$med, .srv$p05, .srv$tier_lo, .srv$tier_hi))
}

## ── ORF programme (figure 4; manuscript/figure_orf.R) ───────────────────────
## Read from disk so the caption quotes computed numbers. Absent on a machine
## where the figure has not been built, in which case the entries are NA.
.orfb <- tryCatch(read_tsv(here("manuscript/orf_breadth_by_tier.tsv"), show_col_types=FALSE),
                  error=function(e) NULL)
.orfc <- tryCatch(read_tsv(here("manuscript/orf_composition_by_tier.tsv"), show_col_types=FALSE),
                  error=function(e) NULL)
.ob <- function(tr, g) if (is.null(.orfb)) NA_real_ else
  .orfb$pct[.orfb$tier==tr & .orfb$gene==g]
.ocgrp <- function(tr, gs) if (is.null(.orfc)) NA_real_ else
  100*sum(.orfc$share[.orfc$tier==tr & .orfc$gene %in% gs])
.ocn <- function(tr) if (is.null(.orfc)) NA_integer_ else
  unique(.orfc$n_lib[.orfc$tier==tr])
.orfn <- if (is.null(.orfb)) NA_integer_ else sum(.orfb$n[.orfb$gene=="E1"])

## ── Assemble the flat `stats` list quoted throughout the manuscript ──────────
.rc <- as.list(rates[rates$sample_class=="clinical",]);   .re <- as.list(rates[rates$sample_class=="engineered",]); .rl <- as.list(rates[rates$sample_class=="cell_line",])
.oe <- as.list(o1[o1$term=="sample_classengineered",]);   .ol <- as.list(o1[o1$term=="sample_classcell_line",])
.tc_ctrl<-as.list(tiersum[tiersum$tier=="control_productive",]); .tc_neo<-as.list(tiersum[tiersum$tier=="neoplasia",])
.tc_uns <-as.list(tiersum[tiersum$tier=="unselected_skin",]);    .tc_cl <-as.list(tiersum[tiersum$tier=="cell_line",]); .tc_oth<-as.list(tiersum[tiersum$tier=="other",])
.tc_ev  <-as.list(tiersum[tiersum$tier=="EV",]);                 .tc_per<-as.list(tiersum[tiersum$tier=="peri_lesional",])
# Non-permissive tissue (cutaneous adipose): no stratified squamous compartment,
# so it is a negative control for the trace background rather than a candidate site.
.tc_npt<-as.list(tiersum[tiersum$tier=="non_permissive_tissue",])
.tc_eng<-as.list(tiersum[tiersum$tier=="engineered",])

## ── E1:E2 integration marker + E2 integrity in the DNA reservoir ─────────────
## E1 and E2 sit on the same early transcript, so their counts should track;
## integration disrupts E2 and inverts the ratio (bin/e1e2_integration_marker.R,
## bin/e2_integrity_from_logan.sh). Read from disk so the text stays live.
.e1e2 <- tryCatch(read_tsv(here("results_e1e2/e1e2_by_tier.tsv"), show_col_types=FALSE),
                  error=function(e) NULL)
.e1e2_disc <- tryCatch(nrow(read_tsv(here("results_e1e2/e1e2_discordant.tsv"), show_col_types=FALSE)),
                       error=function(e) NA_integer_)
.e2int <- tryCatch(read_tsv(here("dna_rna_triangulation/e2_integrity_by_site.tsv"), show_col_types=FALSE),
                   error=function(e) NULL)

## Virus-host junction pilot (rocket SLURM job 508305; results_junction/). Integration
## is called from chimeric reads against a combined GRCh38+HPV index. The
## discriminating statistic is CLUSTERING -- the fraction of a library's
## junctions falling at its single most-supported site -- not the raw count,
## which scales with viral depth and is highest in the episomal EV libraries.
.jx <- tryCatch(read_tsv(here("results_junction/junction_summary.tsv"), show_col_types=FALSE),
                error=function(e) NULL)
.JXCUT <- 0.2          # integrated vs episomal; observed gap is 0.034 to 0.327
.jxint <- function(tiers) if (is.null(.jx)) NA_integer_ else
  sum(.jx$tier %in% tiers & .jx$frac_top >= .JXCUT)
.jxn <- function(tiers) if (is.null(.jx)) NA_integer_ else sum(.jx$tier %in% tiers)
.jxrange <- function(tiers, integrated=TRUE) {
  if (is.null(.jx)) return("NA")
  v <- .jx$frac_top[.jx$tier %in% tiers & (if (integrated) .jx$frac_top >= .JXCUT else .jx$frac_top < .JXCUT)]
  if (!length(v)) return("NA")
  sprintf("%.2f\u2013%.2f", min(v), max(v))
}
## External validation of the ratio on a graded cervical series (PRJNA612785 /
## SRP252916, results_cervix_graded_e1e2/). Four arms of six, one rRNA-depleted
## protocol: HPV-negative normal epithelium, HPV16-positive normal epithelium,
## HPV16+ HSIL, HPV16+ squamous carcinoma. The arms come from the submitting
## study, so the marker is scored against labels it had no hand in assigning,
## in a tissue and a progression it was not developed on.
.cxg <- tryCatch(read_tsv(here("results_cervix_graded_e1e2/e1e2_by_tier.tsv"), show_col_types=FALSE),
                 error=function(e) NULL)
.cx <- function(grp, col="median_log2")
  if (is.null(.cxg) || !any(.cxg$tier==grp)) NA_real_ else .cxg[[col]][.cxg$tier==grp]

.ef <- function(tier) if (is.null(.e1e2)) NA_real_ else .e1e2$fold[.e1e2$tier==tier]
.epct <- function(tier) if (is.null(.e1e2)) NA_real_ else .e1e2$pct_integration_like[.e1e2$tier==tier]
.ecrit <- function(tier, col) if (is.null(.e1e2)) NA_integer_ else .e1e2[[col]][.e1e2$tier==tier]
.e2row <- function(site, set="ge7kb") if (is.null(.e2int)) NULL else
  as.list(.e2int[.e2int$site==site & .e2int$contig_set==set,])
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
  uns_n=.tc_uns$n,   uns_pctfloor=.tc_uns$pct_floor,   uns_pctany=.tc_uns$pct_any,
  cl2_n=.tc_cl$n, oth_n=.tc_oth$n,
  # Absolute above-floor counts: the odds ratios above are large only because the
  # skin baseline is near zero, so the text quotes the raw fractions too.
  ctrl_pos_floor=.tc_ctrl$n_floor, neo_pos_floor=.tc_neo$n_floor, uns_pos_floor=.tc_uns$n_floor,
  neo_maxreads=.tc_neo$max_reads,
  ev_n=.tc_ev$n,     ev_pctfloor=.tc_ev$pct_floor,     ev_medrpm=.tc_ev$med_rpm,   ev_medreads=.tc_ev$med_reads,
  ev_or=ev_f$OR,     ev_lo=ev_f$lo,   ev_hi=ev_f$hi,   ev_p=ev_f$p,
  npt_n=.tc_npt$n,   npt_pctfloor=.tc_npt$pct_floor,   npt_pctany=.tc_npt$pct_any,
  # NB: eng_n above is Model 1's engineered SAMPLE CLASS (whole sweep); this is
  # the engineered TIER in the RPM table. Distinct denominators, distinct names.
  engt_n=.tc_eng$n,  engt_pctfloor=.tc_eng$pct_floor,
  e1e2_cl=.ef("cell_line"), e1e2_ctrl=.ef("control_productive"),
  e1e2_neo=.ef("neoplasia"), e1e2_skin=.ef("unselected_skin"),
  e1e2_disc=.e1e2_disc,
  e1e2_cl_pct=.epct("cell_line"), e1e2_eng_pct=.epct("engineered"),
  e1e2_cl_ratio=.ecrit("cell_line","n_ratio"), e1e2_cl_e67=.ecrit("cell_line","n_e67"),
  e1e2_cl_n=.ecrit("cell_line","n"),
  e1e2_ev_pct=.epct("EV"), e1e2_per_pct=.epct("peri_lesional"),
  e1e2_ctrl_pct=.epct("control_productive"),
  # cervical graded series (external validation)
  cxg_norm=.cx("hpv16_normal_epithelium"), cxg_hsil=.cx("hpv16_hsil"),
  cxg_scc=.cx("hpv16_scc"),
  cxg_norm_n=.cx("hpv16_normal_epithelium","n"), cxg_hsil_n=.cx("hpv16_hsil","n"),
  cxg_scc_n=.cx("hpv16_scc","n"),
  cxg_norm_int=.cx("hpv16_normal_epithelium","n_ratio"),
  cxg_hsil_int=.cx("hpv16_hsil","n_ratio"), cxg_scc_int=.cx("hpv16_scc","n_ratio"),
  # Empirical detection limit: the shallowest library in which integration WAS
  # called bounds what a negative can mean. Every negative cSCC library carries
  # several-fold more viral reads than that, so the negatives are not depth-limited.
  jx_pos_minreads=if (is.null(.jx)) NA else min(.jx$hpv_reads[.jx$frac_top >= .JXCUT]),
  jx_neg_minfold=if (is.null(.jx)) NA else
    min(.jx$hpv_reads[.jx$frac_top < .JXCUT & .jx$tier=="neoplasia"]) /
    min(.jx$hpv_reads[.jx$frac_top >= .JXCUT]),
  jx_n=if (is.null(.jx)) NA_integer_ else nrow(.jx),
  jx_int_range=.jxrange(c("cell_line","neoplasia"), TRUE),
  jx_epi_range=.jxrange(c("control_productive","EV","neoplasia"), FALSE),
  jx_neo_n=.jxn("neoplasia"), jx_neo_int=.jxint("neoplasia"),
  jx_neo_epi=.jxn("neoplasia") - .jxint("neoplasia"),
  jx_cl_int=.jxint("cell_line"), jx_cl_n=.jxn("cell_line"),
  e2_full_n=.e2row("skin")$n_E1, e2_full_ok=.e2row("skin")$n_E2_given_E1,
  e2_vag_n=.e2row("vaginal")$n_E1, e2_vag_ok=.e2row("vaginal")$n_E2_given_E1,
  per_n=.tc_per$n,   per_pctfloor=.tc_per$pct_floor,   per_medrpm=.tc_per$med_rpm,
  per_or=per_f$OR,   per_lo=per_f$lo, per_hi=per_f$hi, per_p=per_f$p,
  m2a_or=.a_neo$OR, m2a_lo=.a_neo$lo, m2a_hi=.a_neo$hi, m2a_p=.a_neo$p,
  m2a_l10_or=.a_l10$OR, m2a_l10_lo=.a_l10$lo, m2a_l10_hi=.a_l10$hi, m2a_l10_p=.a_l10$p,
  m2b_or=.b_neo$OR, m2b_lo=.b_neo$lo, m2b_hi=.b_neo$hi, m2b_p=.b_neo$p,
  m2b_l10_or=.b_l10$OR, m2b_l10_lo=.b_l10$lo, m2b_l10_hi=.b_l10$hi, m2b_l10_p=.b_l10$p,
  fish_or=ft$estimate, fish_lo=ft$conf.int[1], fish_hi=fish_hi, fish_p=ft$p.value,
  neo_depth_or=.mn_l10$OR, neo_depth_lo=.mn_l10$lo, neo_depth_hi=.mn_l10$hi, neo_depth_p=.mn_l10$p,
  m3_ctrl_pct=.m3c$pct, m3_ctrl_prod=.m3c$prod, m3_ctrl_n=.m3c$n,
  m3_neo_pct=.m3n$pct,  m3_neo_prod=.m3n$prod,  m3_neo_n=.m3n$n,
  m3_or=f3$estimate, m3_lo=f3$conf.int[1], m3_hi=f3$conf.int[2], m3_p=f3$p.value,
  cin_n=nrow(.cin), neo_pctfloor_cut=round(.neo_floor_cut,2),
  m2b_or_cut=.bcut$OR, m2b_lo_cut=.bcut$lo, m2b_hi_cut=.bcut$hi, m2b_p_cut=.bcut$p,
  m3_neo_pct_cut=round(.m3nc$pct,1), m3_neo_n_cut=.m3nc$n, m3_neo_prod_cut=.m3nc$prod,
  m3_or_cut=.f3c$estimate, m3_lo_cut=.f3c$conf.int[1], m3_hi_cut=.f3c$conf.int[2], m3_p_cut=.f3c$p.value,
  trace_n=nrow(.tt), trace_resolved=sum(.tt$trace_ok), trace_ntypes=n_distinct(.tres$top_ref),
  trace_skin_n=nrow(.tskin), trace_skin_cut=sum(.tskin$genus_short %in% c("beta","gamma","nu")),
  trace_ctrl_resolved=sum(.tres$tier == "control_productive", na.rm=TRUE),
  # Library completeness
  spot_n=.sr("n"), spot_med=.sr("med"), spot_p05=.sr("p05"),
  spot_tier_lo=.sr("tier_lo"), spot_tier_hi=.sr("tier_hi"),
  # ORF programme (figure 4)
  orf_n=.orfn,
  orf_ctrl_l1=.ob("control_productive","L1"), orf_ctrl_e2=.ob("control_productive","E2"),
  orf_ev_l1=.ob("EV","L1"), orf_neo_l1=.ob("neoplasia","L1"), orf_neo_e7=.ob("neoplasia","E7"),
  orf_skin_l1=.ob("unselected_skin","L1"), orf_skin_e6=.ob("unselected_skin","E6"),
  orf_cl_e6=.ob("cell_line","E6"), orf_cl_e7=.ob("cell_line","E7"),
  orf_cl_e2=.ob("cell_line","E2"), orf_cl_l2=.ob("cell_line","L2"),
  orf_eng_e2=.ob("engineered","E2"), orf_eng_l2=.ob("engineered","L2"),
  orfc_ctrl_late=.ocgrp("control_productive", c("L1","L2")),
  orfc_ctrl_onco=.ocgrp("control_productive", c("E6","E7")),
  orfc_cl_late=.ocgrp("cell_line", c("L1","L2")),
  orfc_cl_onco=.ocgrp("cell_line", c("E6","E7")),
  orfc_neo_late=.ocgrp("neoplasia", c("L1","L2")),
  orfc_neo_onco=.ocgrp("neoplasia", c("E6","E7")),
  orfc_skin_late=.ocgrp("unselected_skin", c("L1","L2")),
  orfc_ctrl_n=.ocn("control_productive"), orfc_cl_n=.ocn("cell_line"),
  orfc_neo_n=.ocn("neoplasia"), orfc_skin_n=.ocn("unselected_skin"),
  # AK vs cSCC split (MODEL 2b-split)
  ak_n=.akc$n[.akc$dx=="neo_AK"], ak_studies=.akc$n_studies[.akc$dx=="neo_AK"],
  ak_medreads=.akc$med_reads[.akc$dx=="neo_AK"], ak_pctany=.akc$pct_any[.akc$dx=="neo_AK"],
  ak_floor=.akc$n_floor[.akc$dx=="neo_AK"], ak_pctfloor=.akc$pct_floor[.akc$dx=="neo_AK"],
  sc_n=.akc$n[.akc$dx=="neo_cSCC"], sc_studies=.akc$n_studies[.akc$dx=="neo_cSCC"],
  sc_pctany=.akc$pct_any[.akc$dx=="neo_cSCC"],
  sc_floor=.akc$n_floor[.akc$dx=="neo_cSCC"], sc_pctfloor=.akc$pct_floor[.akc$dx=="neo_cSCC"],
  ak_or=.ak_any$OR, ak_p=.ak_any$p, ak_fl_or=.ak_fl$OR, ak_fl_p=.ak_fl$p,
  sc_or=.sc_any$OR, sc_p=.sc_any$p, sc_fl_or=.sc_fl$OR, sc_fl_p=.sc_fl$p,
  akvs_or=.akvs$OR, akvs_lo=.akvs$lo, akvs_hi=.akvs$hi, akvs_p=.akvs$p,
  akvf_p=.akvf$p.value, mixed_n=.akc$n[.akc$dx=="neo_mixed"],
  hpv8_typed=.h8_typed, hpv8_trace=.h8_trace, trace_libs=nrow(.tt_pre),
  # Library chemistry (MODEL 2c)
  chem_n=.cq("chem_n"), chem_pct=.cq("chem_pct"),
  chem_total=.cq("chem_total"), chem_cdna_pct=.cq("chem_cdna_pct"),
  chem_neo_ribo=.cq("chem_neo_ribo"), chem_skin_ribo=.cq("chem_skin_ribo"),
  chem_ctrl_ribo=.cq("chem_ctrl_ribo"),
  chem_neo_polya=.cq("chem_neo_polya"), chem_skin_polya=.cq("chem_skin_polya"),
  chem_ribo_or=.cq("chem_ribo_or"), chem_ribo_lo=.cq("chem_ribo_lo"),
  chem_ribo_hi=.cq("chem_ribo_hi"), chem_ribo_p=.cq("chem_ribo_p"),
  chem_tenx_or=.cq("chem_tenx_or"), chem_smart_or=.cq("chem_smart_or"),
  chem_rand_or=.cq("chem_rand_or"),
  chem_ribo_or_all=.cq("chem_ribo_or_all"), chem_ribo_lo_all=.cq("chem_ribo_lo_all"),
  chem_ribo_hi_all=.cq("chem_ribo_hi_all"),
  m2a_or_chem=.cq("m2a_or_chem"), m2a_p_chem=.cq("m2a_p_chem"),
  m2b_or_chem=.cq("m2b_or_chem"), m2b_p_chem=.cq("m2b_p_chem"),
  ribo_n=.cq("ribo_n"), ribo_skin_n=.cq("ribo_skin_n"), ribo_neo_n=.cq("ribo_neo_n"),
  ribo_neo_studies=.cq("ribo_neo_studies"), ribo_neo_medreads=.cq("ribo_neo_medreads"),
  ribo_neo_maxreads=.cq("ribo_neo_maxreads"), ribo_floor=.cq("ribo_floor"),
  ribo_or=.cq("ribo_or"), ribo_lo=.cq("ribo_lo"), ribo_hi=.cq("ribo_hi"), ribo_p=.cq("ribo_p")
), unname)

saveRDS(stats, here("manuscript/stats.rds"))
cat("\nDONE: stats saved to manuscript/stats.rds\n")
