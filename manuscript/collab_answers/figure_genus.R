#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(tidyverse)})
setwd(here::here())
OUT <- here::here("manuscript/collab_answers")
AGG <- c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide")
rp <- read_tsv("results_reframe_csc/hpv_rpm_per_sample.tsv", show_col_types=FALSE)
types <- map_dfr(unlist(lapply(AGG,function(d) list.files(file.path(d,"aggregated"),"_hpv_types\\.tsv$",full.names=TRUE))),
  ~read_tsv(.x,show_col_types=FALSE,col_types=cols(.default="c"))) %>% distinct(sample_id,hpv_reference)
gmap <- read_csv("manuscript/type_genus.csv",show_col_types=FALSE) %>%
  mutate(type=gsub("(REF|nr)$","",hpv_reference)) %>% select(type,genus_short) %>% distinct(type,.keep_all=TRUE)
tierN <- rp %>% count(tier) %>% { setNames(.$n,.$tier) }
tiers <- c("unselected_skin","control_productive","neoplasia","cell_line","other")
GEN <- c("alpha","beta","gamma","mu","nu")

asg <- types %>% mutate(type=gsub("(REF|nr)$","",hpv_reference)) %>% left_join(gmap,by="type") %>%
  left_join(rp %>% select(srr_id,tier),by=c("sample_id"="srr_id")) %>% distinct(id=sample_id,tier,genus_short)
tt <- read_tsv("manuscript/trace_types.tsv",show_col_types=FALSE) %>%
  filter(!top_ref %in% c("NONE","NO_BAM"), top_reads>=5, max_mapq>=20, top_reads/pmax(n_primary,1)>=0.5) %>%
  mutate(type=gsub("(REF|nr)$","",top_ref)) %>% left_join(gmap,by="type") %>%
  left_join(rp %>% select(srr_id,tier),by="srr_id") %>% distinct(id=srr_id,tier,genus_short)
anyd <- bind_rows(asg,tt) %>% distinct(id,tier,genus_short)

mk <- function(d,level) d %>% filter(genus_short %in% GEN) %>% count(tier,genus_short,name="k") %>%
  complete(tier=tiers, genus_short=GEN, fill=list(k=0)) %>% filter(!is.na(tier)) %>%
  mutate(N=tierN[tier], a=k+.5, b=N-k+.5,
         mean=100*a/(a+b), lo=100*qbeta(.025,a,b), hi=100*qbeta(.975,a,b), level=level)

pd <- bind_rows(mk(asg,"Assigned type (Table S1)"), mk(tt,"Trace type (Table S5)"),
                mk(anyd,"Any detection (S1 ∪ S5)")) %>%
  mutate(tier=factor(recode(tier, unselected_skin="Unselected skin", control_productive="Productive\ncontrols",
                            neoplasia="Neoplasia", cell_line="Cell line", other="Other"),
                     levels=c("Unselected skin","Productive\ncontrols","Neoplasia","Cell line","Other")),
         Genus=factor(str_to_title(genus_short), levels=str_to_title(GEN)),
         level=factor(level, levels=c("Assigned type (Table S1)","Trace type (Table S5)","Any detection (S1 ∪ S5)")),
         # zero counts: no estimate, mark at the floor instead
         mean_p=ifelse(k>0, mean, NA_real_), lo_p=ifelse(k>0, lo, NA_real_), hi_p=ifelse(k>0, hi, NA_real_))
FLOOR <- 0.01
p <- ggplot(pd, aes(Genus, mean_p, colour=Genus)) +
  geom_errorbar(aes(ymin=lo_p, ymax=hi_p), width=.18, linewidth=.4, na.rm=TRUE) +
  geom_point(size=2.1, na.rm=TRUE) +
  geom_point(data=~filter(.x, k==0), aes(y=FLOOR), shape=4, size=1.5, colour="grey55", na.rm=TRUE) +
  geom_text(data=~filter(.x, k>0), aes(y=hi_p, label=k), vjust=-.6, size=2.4, show.legend=FALSE) +
  facet_grid(level ~ tier) +
  scale_y_log10(breaks=c(.01,.1,1,10,100), labels=c("0","0.1","1","10","100"),
                limits=c(FLOOR, 400)) +
  scale_colour_brewer(palette="Dark2", guide="none") +
  labs(x="HPV genus", y="Prevalence within tier (%, log scale)",
       title="HPV genus prevalence by population tier",
       subtitle=paste0("Beta-Binomial posterior mean with 95% CrI (Beta(0.5,0.5) prior); number above each point = libraries detected.\n",
                       "Grey × = no library detected in that tier (zero count). Tier N: skin ", tierN["unselected_skin"],
                       ", controls ", tierN["control_productive"], ", neoplasia ", tierN["neoplasia"],
                       ", cell line ", tierN["cell_line"], ", other ", tierN["other"], ".")) +
  theme_bw(base_size=9.5) +
  theme(strip.text.y=element_text(size=7.5), panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=7.3, colour="grey30"))
ggsave(file.path(OUT,"figure_genus_prevalence.png"), p, width=10.5, height=6.4, dpi=200)
cat("ok\n")
print(pd %>% filter(level=="Any detection (S1 ∪ S5)", k>0) %>%
  transmute(tier,Genus,k,N,pct=sprintf("%.2f [%.2f-%.2f]",mean,lo,hi)) %>% as.data.frame(), right=FALSE)
