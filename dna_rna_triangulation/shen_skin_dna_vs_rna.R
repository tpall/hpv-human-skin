#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Skin DNA reservoir (Shen et al., PVs-in-Logan) vs skin RNA activity (this study)
#
# Cross-modality triangulation of the "present but quiet" thesis at the TYPE level.
# The two datasets share no libraries (Logan = DNA metagenomes; ours = RNA-seq;
# disjoint SRRs), so this is a COMPOSITION comparison, not a per-library join:
#   - DNA presence : which HPV types dominate Shen's human skin-metagenome
#                    libraries (# libraries carrying each type/genus).
#   - RNA activity : per-type library counts in our tiers (unselected skin;
#                    productive = warts control_productive + EV).
#
# Inputs (Shen, auto-downloaded from github.com/syueqiao/pvs_in_logan):
#   files/all_hits_info.list                    SRA RunInfo for HPV+ libraries
#   outputs/2026.04.06.master_table_all_pvs.tsv library -> cluster
#   outputs/2026.04.13.cluster_taxonomy_demarcation.tsv  cluster -> best_hit_name
# Inputs (ours): manuscript/type_genus.csv, results_reframe_fixed tier catalogue.
# Outputs: skin_dna_vs_rna_{types,genus,summary}.csv + figure_dna_rna_triangulation.png
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(cowplot); library(patchwork); library(ggrepel)
})
theme_set(theme_cowplot(font_size = 11))
outdir <- here("dna_rna_triangulation"); cache <- file.path(outdir, "cache")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

# ---- fetch Shen data (cached) ----
base <- "https://raw.githubusercontent.com/syueqiao/pvs_in_logan/main"
shen <- c(hits = "files/all_hits_info.list",
          master = "outputs/2026.04.06.master_table_all_pvs.tsv",
          demarc = "outputs/2026.04.13.cluster_taxonomy_demarcation.tsv")
paths <- map_chr(shen, function(p) {
  dest <- file.path(cache, basename(p))
  if (!file.exists(dest)) download.file(file.path(base, p), dest, mode = "wb", quiet = TRUE)
  dest
})

# ---- type -> genus (ours), normalised to a short Greek label ----
GEN_LEVELS <- c("Alpha", "Beta", "Gamma", "Mu", "Nu", "Unresolved", "Novel")
tg <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
  transmute(type = toupper(as.character(type)),
            genus = str_to_title(genus_short))
type2genus <- deframe(tg %>% distinct(type, genus))

# ---- subtype-letter normalisation, applied to BOTH sides ----
# PaVE names variants of a type with a trailing letter: HPV38B is a subtype of
# HPV38, HPV14D of HPV14, HPV57B of HPV57. Our typing panel carries the subtype
# labels and Shen's best_hit_name carries the bare type, so without this the two
# sides never join and the subtype's libraries are silently lost to "Unresolved".
#
# It was not cosmetic. Before this fix HPV57 and HPV14 were reported as
# PRODUCTIVE TYPES WITH ZERO DNA DETECTIONS, one of the headline claims, while
# HPV57B (13 libraries) and HPV14D (58 libraries) sat unjoined in the same
# table. HPV38B (63 libraries) was likewise split from HPV38 (32).
#
# Only a single trailing letter directly after a digit is stripped, so novel
# candidate labels ending in letters (MKN1NR, MSK190NR) are left untouched.
norm_type <- function(x) str_remove(toupper(x), "(?<=[0-9])[A-Z]$")

# ---- parse Shen best_hit_name -> HPV type number + genus ----
parse_type <- function(x) {                     # only "human papillomavirus N" / "HPV N"
  m <- str_match(x, regex("human papillomavirus(?:\\s+type)?\\s+(\\d+[A-Za-z]?)|\\bHPV[- ]?(\\d+[A-Za-z]?)", ignore_case = TRUE))
  norm_type(coalesce(m[, 2], m[, 3]))
}
parse_genus_kw <- function(x) {                 # e.g. "Gammapapillomavirus 22"
  str_to_title(str_match(x, regex("(Alpha|Beta|Gamma|Mu|Nu)papillomavirus", ignore_case = TRUE))[, 2])
}
demarc <- read_tsv(paths["demarc"], show_col_types = FALSE) %>%
  mutate(type = parse_type(best_hit_name),
         genus_kw = parse_genus_kw(best_hit_name),
         is_novel = status == "novel" | taxonomy_level == "no_hit",
         genus = case_when(
           !is.na(type) & type %in% names(type2genus) ~ unname(type2genus[type]),
           !is.na(genus_kw) ~ genus_kw,
           is_novel ~ "Novel",
           TRUE ~ "Unresolved")) %>%
  select(cluster_number, type, genus, is_novel)

# ---- skin-metagenome HPV+ Runs (Shen RunInfo) ----
# The skin-metagenome taxid is also used for environmental transit surfaces and
# non-human animal skin; drop those (build_skin_run_exclusions.py writes the table).
excl_path <- file.path(outdir, "skin_run_exclusions.tsv")
skin_excl <- if (file.exists(excl_path)) {
  read_tsv(excl_path, show_col_types = FALSE)
} else {
  warning("skin_run_exclusions.tsv missing -- run build_skin_run_exclusions.py; ",
          "keeping non-human-skin libraries")
  tibble(run = character())
}
skin_sci <- c("human skin metagenome", "skin metagenome")
skin_runs <- read_csv(paths["hits"], show_col_types = FALSE) %>%
  filter(str_to_lower(ScientificName) %in% skin_sci) %>% pull(Run) %>% unique()
n_skin_runs_raw <- length(skin_runs)
skin_runs <- setdiff(skin_runs, skin_excl$run)
message(sprintf("skin runs: %d raw -> %d after dropping %d non-human-skin (%s)",
                n_skin_runs_raw, length(skin_runs),
                n_skin_runs_raw - length(skin_runs),
                paste(sort(unique(skin_excl$reason)), collapse = "; ")))

# ---- skin libraries -> clusters -> type/genus (master table) ----
lib_cluster <- read_tsv(paths["master"], show_col_types = FALSE) %>%
  filter(library %in% skin_runs) %>%
  distinct(library, cluster_number) %>%
  left_join(demarc, by = "cluster_number")
n_skin_lib <- n_distinct(lib_cluster$library)

# per-genus / per-type DNA frequency = # distinct skin libraries carrying it
dna_genus <- lib_cluster %>% distinct(library, genus) %>% count(genus, name = "dna_skin_libs")
dna_type  <- lib_cluster %>% filter(!is.na(type)) %>% distinct(library, type) %>%
  count(type, name = "dna_skin_libs")
novel_libs <- lib_cluster %>% filter(is_novel) %>% distinct(library) %>% nrow()

# ---- RNA tiers (ours): present in skin vs productive (warts + EV) ----
rna <- read_tsv(here("results_reframe_fixed/type_catalogue_by_tier.tsv"), show_col_types = FALSE) %>%
  transmute(type = norm_type(str_remove(str_remove(hpv_reference, "REF$"), "^HPV")),
            rna_skin_present = unselected_skin,
            rna_productive = control_productive + EV) %>%
  group_by(type) %>%                     # subtypes collapse onto their parent
  summarise(across(everything(), sum), .groups = "drop")

# ---- master per-type comparison table ----
cmp <- full_join(dna_type, rna, by = "type") %>%
  mutate(across(c(dna_skin_libs, rna_skin_present, rna_productive), ~ replace_na(., 0)),
         genus = factor(coalesce(unname(type2genus[type]), "Unresolved"), GEN_LEVELS)) %>%
  filter(dna_skin_libs > 0 | rna_skin_present > 0 | rna_productive > 0) %>%
  arrange(desc(dna_skin_libs)) %>%
  mutate(dna_rank = row_number())
write_csv(cmp %>% arrange(desc(dna_skin_libs)) %>%
            transmute(type = paste0("HPV", type), genus, dna_skin_libs,
                      rna_skin_present, rna_productive, dna_rank),
          file.path(outdir, "skin_dna_vs_rna_types.csv"))

# genus composition across the three views (proportions)
genus_view <- cmp %>%
  pivot_longer(c(dna_skin_libs, rna_skin_present, rna_productive),
               names_to = "view", values_to = "n") %>%
  group_by(view, genus) %>% summarise(n = sum(n), .groups = "drop") %>%
  group_by(view) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(view = recode(view, dna_skin_libs = "DNA skin\n(reservoir)",
                       rna_skin_present = "RNA skin\n(present)",
                       rna_productive = "RNA productive\n(warts + EV)"),
         view = factor(view, c("DNA skin\n(reservoir)", "RNA skin\n(present)", "RNA productive\n(warts + EV)")))
write_csv(genus_view, file.path(outdir, "skin_dna_vs_rna_genus.csv"))

# ---- overlap / summary scalars ----
topN <- 20
top_dna  <- cmp %>% slice_max(dna_skin_libs, n = topN) %>% pull(type)
top_prod <- cmp %>% filter(rna_productive > 0) %>% slice_max(rna_productive, n = topN) %>% pull(type)
deframe_p <- function(g) {
  x <- dna_genus %>% filter(genus == g) %>% pull(dna_skin_libs); if (length(x)) x else 0 }
summary_tbl <- tibble(
  n_skin_libs = n_skin_lib,
  n_types_dna = nrow(dna_type),
  pct_beta  = round(100 * deframe_p("Beta")  / n_skin_lib, 1),
  pct_gamma = round(100 * deframe_p("Gamma") / n_skin_lib, 1),
  pct_alpha = round(100 * deframe_p("Alpha") / n_skin_lib, 1),
  pct_novel = round(100 * novel_libs / n_skin_lib, 1),
  overlap_top20 = length(intersect(top_dna, top_prod)),
  overlap_types = paste(sort(paste0("HPV", intersect(top_dna, top_prod))), collapse = ", "),
  prod_zero_in_dna = paste(sort(paste0("HPV", setdiff(top_prod, cmp %>% filter(dna_skin_libs > 0) %>% pull(type)))), collapse = ", "))
write_csv(summary_tbl, file.path(outdir, "skin_dna_vs_rna_summary.csv"))

# ---- figure ----
gcol <- c(Alpha = "#1B9E77", Beta = "#D95F02", Gamma = "#7570B3", Mu = "#E7298A",
          Nu = "#66A61E", Unresolved = "grey70", Novel = "grey40")
lg <- function(x) log10(x + 1)
lbrk <- c(0, 1, 5, 20, 100); lval <- log10(lbrk + 1)

pA <- ggplot(genus_view, aes(view, prop, fill = genus)) +
  geom_col(width = .72, colour = "white", linewidth = .2) +
  scale_fill_manual(values = gcol, guide = "none", drop = TRUE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .02))) +
  labs(x = NULL, y = "Composition") +
  theme(axis.text.x = element_text(size = 9))

lab_types <- cmp %>% filter(type %in% union(top_dna[1:8], top_prod[1:8]))
pB <- ggplot(cmp, aes(lg(dna_skin_libs), lg(rna_productive), colour = genus)) +
  annotate("rect", xmin = lg(15), xmax = Inf, ymin = lg(5), ymax = Inf, fill = "grey85", alpha = .5) +
  annotate("text", x = lg(140), y = lg(24), label = "abundant AND productive:\nonly HPV2",
           hjust = 1, vjust = 1, size = 3, colour = "grey35", lineheight = .9) +
  geom_jitter(width = .02, height = .02, size = 1.9, alpha = .85) +
  geom_text_repel(data = lab_types, aes(label = paste0("HPV", type)),
                  size = 2.9, max.overlaps = 20, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = gcol, name = "Genus", drop = TRUE) +
  scale_x_continuous(breaks = lval, labels = lbrk) +
  scale_y_continuous(breaks = lval, labels = lbrk) +
  labs(x = "DNA skin-metagenome libraries (Shen)",
       y = "RNA productive libraries (warts + EV, this study)")

fig <- (pA + pB) + plot_layout(widths = c(1, 1.55)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "plain"))
ggsave(file.path(outdir, "figure_dna_rna_triangulation.png"), fig,
       width = 11, height = 5, dpi = 300, bg = "white")

cat(sprintf("skin libs=%d | Beta %.1f%% Gamma %.1f%% Alpha %.1f%% novel %.1f%% | top20 overlap=%d (%s)\n",
            summary_tbl$n_skin_libs, summary_tbl$pct_beta, summary_tbl$pct_gamma,
            summary_tbl$pct_alpha, summary_tbl$pct_novel, summary_tbl$overlap_top20, summary_tbl$overlap_types))
