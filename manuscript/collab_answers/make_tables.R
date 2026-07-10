#!/usr/bin/env Rscript
# Tables + figure answering the collaborator questions (July 2026).
# Run from the repo root:  Rscript manuscript/collab_answers/make_tables.R
suppressPackageStartupMessages({library(tidyverse); library(here)})
OUT <- here("manuscript/collab_answers")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

AGG <- c("results_full_v2","results_lesions","results_lesions_csc","results_lesions_wide")
rp  <- read_tsv(here("results_reframe_csc/hpv_rpm_per_sample.tsv"), show_col_types = FALSE)

gather <- function(pat)
  map_dfr(unlist(lapply(AGG, function(d) list.files(here(d,"aggregated"), pat, full.names = TRUE))),
          ~ read_tsv(.x, show_col_types = FALSE, col_types = cols(.default = "c")))

types <- gather("_hpv_types\\.tsv$") %>%
  mutate(across(c(ref_length, read_count, covered_bases), as.integer),
         across(c(coverage_breadth, mean_depth), as.double)) %>%
  distinct(sample_id, hpv_reference, .keep_all = TRUE)
tx <- gather("_transcript_classes\\.tsv$") %>% distinct(sample_id, gene, .keep_all = TRUE)

gmap <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
  mutate(type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  select(type, genus_short) %>% distinct(type, .keep_all = TRUE)

# featureCounts never emits E4 (PaVE names it "E1^E4") and E5 is always 0
# (PaVE names it E5_ALPHA/E5_BETA/E5_GAMMA/E5_DELTA), so only these 7 exist.
ORFS <- c("E1","E2","E5","E6","E7","L1","L2")
orf <- tx %>% filter(gene %in% ORFS) %>% mutate(read_count = as.integer(read_count)) %>%
  pivot_wider(names_from = gene, values_from = read_count, id_cols = sample_id)
prod <- tx %>% filter(gene == "PRODUCTIVE_INFECTION") %>% select(sample_id, productive = read_count)
nty  <- types %>% count(sample_id, name = "n_types_in_library")

## ── Table A: every above-threshold HPV call (type, reads, coverage) ──────────
tblA <- types %>%
  mutate(type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  left_join(gmap, by = "type") %>% left_join(nty, by = "sample_id") %>%
  left_join(rp %>% select(srr_id, tier, cohort, study, diagnosis, rpm, hpv_reads, total_reads),
            by = c("sample_id" = "srr_id")) %>%
  transmute(sample = sample_id, tier, cohort, study, diagnosis,
            hpv_type = type, genus = genus_short, n_types_in_library,
            aligned_reads = read_count, ref_length_bp = ref_length, covered_bases,
            coverage_breadth = round(coverage_breadth, 4), mean_depth = round(mean_depth, 2),
            kraken_hpv_reads = hpv_reads, total_reads, rpm = round(rpm, 3)) %>%
  arrange(tier, sample, desc(aligned_reads))
write_tsv(tblA, file.path(OUT, "tableA_assigned_type_calls.tsv"))

## ── Table B: type + reads + ORF profile, single-type libraries only ──────────
tblB <- tblA %>% filter(n_types_in_library == 1) %>%
  left_join(orf, by = c("sample" = "sample_id")) %>%
  left_join(prod, by = c("sample" = "sample_id")) %>%
  select(sample, tier, study, diagnosis, hpv_type, genus, aligned_reads,
         coverage_breadth, mean_depth, all_of(ORFS), productive) %>%
  arrange(tier, desc(aligned_reads))
write_tsv(tblB, file.path(OUT, "tableB_type_reads_orf.tsv"))

## ── Table C: co-infections ───────────────────────────────────────────────────
tblC <- types %>% mutate(type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  left_join(gmap, by = "type") %>%
  left_join(rp %>% select(srr_id, tier, study), by = c("sample_id" = "srr_id")) %>%
  group_by(sample_id, tier, study) %>% filter(n() >= 2) %>%
  arrange(desc(read_count), .by_group = TRUE) %>%
  summarise(n_types = n(),
            combination   = paste(type, collapse = " + "),
            genera        = paste(coalesce(genus_short, "?"), collapse = " + "),
            reads         = paste(read_count, collapse = " / "),
            breadth       = paste(sprintf("%.2f", coverage_breadth), collapse = " / "),
            mean_depth    = paste(sprintf("%.1f", mean_depth), collapse = " / "),
            .groups = "drop") %>%
  arrange(tier, desc(n_types))
write_tsv(tblC, file.path(OUT, "tableC_coinfections.tsv"))

## ── Table D: HPV-positive cell lines, named ─────────────────────────────────
fl <- read_tsv(here("results_full_v2/metadata/cell_line_flags_curated.tsv"), show_col_types = FALSE)
raw <- map_dfr(c("results_lesions_csc","results_lesions_wide","results_lesions_wide2"), function(d){
  f <- here(d, "lesions_raw.csv"); if (file.exists(f))
    read_csv(f, show_col_types = FALSE) %>% mutate(across(everything(), as.character)) else NULL }) %>%
  distinct(srr_id, .keep_all = TRUE)
meta <- read_csv(here("results_full_v2/metadata/samplesheet_enriched.csv"), show_col_types = FALSE) %>%
  mutate(across(everything(), as.character))
src <- bind_rows(
  fl   %>% transmute(srr_id, txt = source_text),
  meta %>% transmute(srr_id, txt = paste(coalesce(title,""), coalesce(tissue_source,""))),
  raw  %>% transmute(srr_id, txt = paste(coalesce(title,""), coalesce(tissue_source,""),
                                         coalesce(cell_line,""), coalesce(characteristics,"")))) %>%
  group_by(srr_id) %>% summarise(txt = paste(unique(na.omit(txt)), collapse = " | "), .groups = "drop")

LINES <- c("HeLa"="hela", "SiHa"="siha", "CaSki"="caski", "HaCaT"="hacat", "HaCaT31"="hacat31",
  "N/TERT"="n[/ -]?tert", "A431"="a431", "TIGK"="tigk", "MNT-1"="mnt[ -]?1", "A375"="a375",
  "HEK293"="hek[ -]?293|293t", "RWPE-1"="rwpe[ -]?1", "HK2"="\\bhk2\\b", "SCL-1"="scl[ -]?1",
  "HSC-1"="hsc[ -]?1", "W12"="w12", "iPSC"="\\bipsc\\b", "SK-MEL"="sk[ -]?mel")
name_of <- function(t){ if (is.na(t)) return(NA_character_)
  h <- names(LINES)[map_lgl(LINES, ~ str_detect(str_to_lower(t), .x))]
  if (!length(h)) NA_character_ else paste(h, collapse = ";") }

tblD <- types %>%
  semi_join(rp %>% filter(tier == "cell_line", hpv_reads >= 10), by = c("sample_id" = "srr_id")) %>%
  mutate(type = gsub("(REF|nr)$", "", hpv_reference)) %>%
  left_join(src, by = c("sample_id" = "srr_id")) %>%
  left_join(rp %>% select(srr_id, study, rpm), by = c("sample_id" = "srr_id")) %>%
  left_join(orf, by = "sample_id") %>%
  mutate(cell_line = coalesce(map_chr(txt, name_of), "unidentified")) %>%
  group_by(cell_line, study, hpv_type = type) %>%
  summarise(n_libraries = n(),
            coverage_breadth = sprintf("%.2f-%.2f", min(coverage_breadth), max(coverage_breadth)),
            mean_depth = sprintf("%.0f-%.0f", min(mean_depth), max(mean_depth)),
            median_rpm = round(median(rpm), 1),
            E1 = sum(E1), E2 = sum(E2), E6 = sum(E6), E7 = sum(E7), L1 = sum(L1), L2 = sum(L2),
            .groups = "drop") %>%
  arrange(desc(n_libraries))
write_tsv(tblD, file.path(OUT, "tableD_hpv_positive_cell_lines.tsv"))

cat("wrote:\n"); for (f in list.files(OUT, "\\.tsv$")) cat("  ", f, "\n")
cat("\nTable A rows:", nrow(tblA), " B:", nrow(tblB), " C:", nrow(tblC), " D:", nrow(tblD), "\n")
