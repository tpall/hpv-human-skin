#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Corrected type-level DNA-reservoir vs RNA-activity comparison.
#
# `skin_dna_vs_rna_types.csv` (from shen_skin_dna_vs_rna.R) builds its
# `rna_skin_present` column from the ASSIGNED type catalogue only, i.e. types
# clearing coverage breadth >= 0.10 and mean depth >= 2. Beta and gamma types
# in skin sit mostly below that bar, so the column reports zero for types the
# detection tier does in fact recover, and they must not be described as absent
# from skin RNA (PLAN.md 3.2).
#
# This rebuilds the comparison with the trace tier included, and separates the
# two questions that matter:
#   rna_*_skin    unselected-skin libraries only, assigned and trace
#   rna_any_all   any library in the whole RNA cohort, assigned or trace
# The second is what a claim of the form "type X is absent from our RNA" needs.
#
# Output: skin_dna_vs_rna_types_corrected.csv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("manuscript/_paths.R"))
source(here("dna_rna_triangulation/_type_norm.R"))

nm <- function(x) norm_type(str_remove(str_remove(x, "REF$"), "^HPV"))

## DNA side: keep the published per-type library counts as they stand.
dna <- read_csv(here("dna_rna_triangulation/skin_dna_vs_rna_types.csv"),
                show_col_types = FALSE) %>%
  transmute(type = nm(type), genus, dna_skin_libs, dna_rank,
            rna_assigned_skin_old = rna_skin_present, rna_productive)

## RNA side, assigned tier: unchanged definition, recomputed for the join.
assigned <- read_tsv(type_catalogue_tsv(), show_col_types = FALSE) %>%
  transmute(type = nm(hpv_reference), rna_assigned_skin = unselected_skin) %>%
  group_by(type) %>% summarise(rna_assigned_skin = sum(rna_assigned_skin), .groups = "drop")

## RNA side, assigned tier across every library, not only unselected skin.
assigned_all <- read_tsv(type_catalogue_tsv(), show_col_types = FALSE) %>%
  mutate(type = nm(hpv_reference)) %>%
  group_by(type) %>%
  summarise(rna_assigned_all = sum(across(where(is.numeric)), na.rm = TRUE), .groups = "drop")

## RNA side, trace tier. Same confidence rule as the manuscript's Table S6:
## at least five best-hit reads, MAPQ >= 20, and a majority of the library's
## primary alignments on that reference.
tiers <- read_tsv(rpm_per_sample_tsv(), show_col_types = FALSE) %>% select(srr_id, tier)
trace <- read_tsv(here("manuscript/trace_types.tsv"), show_col_types = FALSE) %>%
  filter(!top_ref %in% c("NONE", "NO_BAM"), top_reads >= 5, max_mapq >= 20,
         top_reads / pmax(n_primary, 1) >= 0.5) %>%
  mutate(type = nm(top_ref)) %>%
  left_join(tiers, by = "srr_id")
trace_skin <- trace %>% filter(tier == "unselected_skin") %>% count(type, name = "rna_trace_skin")
trace_all  <- trace %>% count(type, name = "rna_trace_all")

cmp <- dna %>%
  left_join(assigned, by = "type") %>% left_join(assigned_all, by = "type") %>%
  left_join(trace_skin, by = "type") %>% left_join(trace_all, by = "type") %>%
  mutate(across(starts_with("rna_"), ~ replace_na(., 0)),
         rna_any_skin = rna_assigned_skin + rna_trace_skin,
         rna_any_all  = rna_assigned_all + rna_trace_all) %>%
  arrange(dna_rank) %>%
  transmute(type = paste0("HPV", type), genus, dna_skin_libs, dna_rank,
            rna_assigned_skin, rna_trace_skin, rna_any_skin,
            rna_assigned_all, rna_trace_all, rna_any_all, rna_productive)

write_csv(cmp, here("dna_rna_triangulation/skin_dna_vs_rna_types_corrected.csv"))

cat("Top 10 skin DNA reservoir types against RNA recovery:\n")
print(as.data.frame(cmp %>% slice_head(n = 10) %>%
  select(type, genus, dna_skin_libs, rna_any_skin, rna_any_all, rna_productive)))

top <- cmp %>% slice_head(n = 10)
cat(sprintf("\nTop 10: DNA %d-%d libraries each; RNA %d-%d libraries each across the whole cohort.\n",
            min(top$dna_skin_libs), max(top$dna_skin_libs),
            min(top$rna_any_all), max(top$rna_any_all)))
cat(sprintf("Zero RNA evidence anywhere in the cohort: %d of the top 10 (%s)\n",
            sum(top$rna_any_all == 0), paste(top$type[top$rna_any_all == 0], collapse = ", ")))
cat(sprintf("Types the old column called absent from skin but which the trace tier recovers: %d (%s)\n",
            sum(cmp$rna_assigned_skin == 0 & cmp$rna_trace_skin > 0),
            paste(cmp$type[cmp$rna_assigned_skin == 0 & cmp$rna_trace_skin > 0], collapse = ", ")))
