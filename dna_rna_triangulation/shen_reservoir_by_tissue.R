#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# The HPV DNA reservoir across body sites (Shen et al., PVs-in-Logan)
#
# Companion to shen_skin_dna_vs_rna.R. That script contrasted the SKIN DNA
# reservoir with our skin RNA activity. Our RNA-seq survey is skin-only
# (tissue_category: nahk 9037; no vaginal/gut), so vaginal and gut cannot be
# given a matched RNA-activity axis. What they CAN show is the DNA reservoir
# itself: does the petabase-scale DNA reservoir recover the known genus tropism
# (skin = beta/gamma commensal; vaginal = alpha mucosal; gut = ?).
#
# This is a DNA-only, cross-body-site composition comparison built from the same
# Shen inputs (cached by shen_skin_dna_vs_rna.R). Composition, not prevalence:
# counts are # HPV+ libraries per body site carrying each type/genus (numerators).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(cowplot); library(patchwork)
})
theme_set(theme_cowplot(font_size = 11))
outdir <- here("dna_rna_triangulation"); cache <- file.path(outdir, "cache")

base <- "https://raw.githubusercontent.com/syueqiao/pvs_in_logan/main"
shen <- c(hits = "files/all_hits_info.list",
          master = "outputs/2026.04.06.master_table_all_pvs.tsv",
          demarc = "outputs/2026.04.13.cluster_taxonomy_demarcation.tsv")
paths <- map_chr(shen, function(p) {
  dest <- file.path(cache, basename(p))
  if (!file.exists(dest)) download.file(file.path(base, p), dest, mode = "wb", quiet = TRUE)
  dest
})

# ---- type -> genus (ours) ----
tg <- read_csv(here("manuscript/type_genus.csv"), show_col_types = FALSE) %>%
  transmute(type = toupper(as.character(type)), genus = str_to_title(genus_short))
type2genus <- deframe(tg %>% distinct(type, genus))

# ---- cluster -> type/genus from Shen demarcation ----
parse_type <- function(x) {
  m <- str_match(x, regex("human papillomavirus(?:\\s+type)?\\s+(\\d+[A-Za-z]?)|\\bHPV[- ]?(\\d+[A-Za-z]?)", ignore_case = TRUE))
  toupper(coalesce(m[, 2], m[, 3]))
}
parse_genus_kw <- function(x)
  str_to_title(str_match(x, regex("(Alpha|Beta|Gamma|Mu|Nu)papillomavirus", ignore_case = TRUE))[, 2])
demarc <- read_tsv(paths["demarc"], show_col_types = FALSE) %>%
  mutate(type = parse_type(best_hit_name), genus_kw = parse_genus_kw(best_hit_name),
         is_novel = status == "novel" | taxonomy_level == "no_hit",
         genus = case_when(
           !is.na(type) & type %in% names(type2genus) ~ unname(type2genus[type]),
           !is.na(genus_kw) ~ genus_kw, is_novel ~ "Novel", TRUE ~ "Unresolved")) %>%
  select(cluster_number, type, genus, is_novel)

runinfo <- read_csv(paths["hits"], show_col_types = FALSE) %>% transmute(Run, sci = str_to_lower(ScientificName))
master  <- read_tsv(paths["master"], show_col_types = FALSE) %>% distinct(library, cluster_number)

# ---- body-site definitions (ScientificName buckets) ----
TISSUES <- list(
  Skin    = c("human skin metagenome", "skin metagenome"),
  Vaginal = c("vaginal metagenome", "human vaginal metagenome"),
  Gut     = c("human gut metagenome", "gut metagenome"))
GEN_LEVELS <- c("Alpha", "Beta", "Gamma", "Mu", "Nu", "Unresolved", "Novel")

# The skin-metagenome taxid also covers environmental transit surfaces and
# non-human animal skin; drop those (build_skin_run_exclusions.py writes the table).
excl_path <- file.path(outdir, "skin_run_exclusions.tsv")
skin_excl <- if (file.exists(excl_path)) {
  read_tsv(excl_path, show_col_types = FALSE)$run
} else {
  warning("skin_run_exclusions.tsv missing -- run build_skin_run_exclusions.py; ",
          "keeping non-human-skin libraries")
  character()
}

reservoir <- imap_dfr(TISSUES, function(scis, tis) {
  runs <- runinfo %>% filter(sci %in% scis) %>% pull(Run) %>% unique()
  if (tis == "Skin") runs <- setdiff(runs, skin_excl)
  lc <- master %>% filter(library %in% runs) %>% left_join(demarc, by = "cluster_number")
  nlib <- n_distinct(lc$library)
  genus <- lc %>% distinct(library, genus) %>% count(genus, name = "libs") %>%
    mutate(tissue = tis, n_lib = nlib, prop = libs / nlib)
  novel <- tibble(genus = "Novel_any", libs = n_distinct(lc$library[lc$is_novel]),
                  tissue = tis, n_lib = nlib, prop = n_distinct(lc$library[lc$is_novel]) / nlib)
  bind_rows(genus, novel)
})

top_types <- imap_dfr(TISSUES, function(scis, tis) {
  runs <- runinfo %>% filter(sci %in% scis) %>% pull(Run) %>% unique()
  if (tis == "Skin") runs <- setdiff(runs, skin_excl)
  lc <- master %>% filter(library %in% runs) %>% left_join(demarc, by = "cluster_number")
  nlib <- n_distinct(lc$library)
  lc %>% filter(!is.na(type)) %>% distinct(library, type, genus) %>%
    count(type, genus, name = "libs") %>% mutate(tissue = tis, n_lib = nlib, prop = libs / nlib) %>%
    slice_max(libs, n = 10)
})

write_csv(reservoir %>% mutate(prop = round(prop, 4)), file.path(outdir, "reservoir_by_tissue_genus.csv"))
write_csv(top_types %>% mutate(prop = round(prop, 4)) %>%
            transmute(tissue, type = paste0("HPV", type), genus, libs, n_lib, prop),
          file.path(outdir, "reservoir_by_tissue_top_types.csv"))

# ---- figure ----
gcol <- c(Alpha = "#1B9E77", Beta = "#D95F02", Gamma = "#7570B3", Mu = "#E7298A",
          Nu = "#66A61E", Unresolved = "grey70", Novel = "grey40")
tis_lab <- reservoir %>% distinct(tissue, n_lib) %>%
  mutate(lab = sprintf("%s\n(n=%d)", tissue, n_lib)) %>%
  arrange(factor(tissue, names(TISSUES)))
lab_map <- deframe(tis_lab %>% select(tissue, lab))

# normalise each body site to 100%: composition of genus detections (a library
# with a co-infection contributes to more than one genus, so raw library-
# prevalences sum above 100%; here we show their relative composition).
pA <- reservoir %>% filter(genus != "Novel_any") %>%
  group_by(tissue) %>% mutate(comp = libs / sum(libs)) %>% ungroup() %>%
  mutate(genus = factor(genus, GEN_LEVELS),
         tissue = factor(lab_map[tissue], lab_map[names(TISSUES)])) %>%
  ggplot(aes(tissue, comp, fill = genus)) +
  geom_col(width = .72, colour = "white", linewidth = .2) +
  scale_fill_manual(values = gcol, name = "Genus", drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .02))) +
  labs(x = NULL, y = "Genus composition\n(% of type detections)",
       title = "DNA reservoir genus tropism by body site")

pB <- top_types %>%
  mutate(genus = factor(genus, GEN_LEVELS),
         row = fct_reorder(paste0("HPV", type, "__", tissue), libs),
         tissue = factor(lab_map[tissue], lab_map[names(TISSUES)])) %>%
  ggplot(aes(prop, row, fill = genus)) +
  geom_col(width = .78, colour = "white", linewidth = .2) +
  facet_wrap(~ tissue, scales = "free_y") +
  scale_y_discrete(labels = function(x) str_remove(x, "__.*")) +
  scale_fill_manual(values = gcol, name = "Genus", drop = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "% of body-site HPV+ libraries", y = NULL,
       title = "Top 10 types per body site")

fig <- (pA / pB) + plot_layout(heights = c(1, 1.35), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "plain"))
ggsave(file.path(outdir, "figure_reservoir_by_tissue.png"), fig,
       width = 10, height = 9, dpi = 300, bg = "white")

wide <- reservoir %>% filter(genus %in% c("Alpha", "Beta", "Gamma")) %>%
  select(tissue, genus, prop) %>% pivot_wider(names_from = genus, values_from = prop)
cat("=== DNA reservoir genus composition by body site ===\n"); print(as.data.frame(wide), digits = 3)
cat("\n=== novel fraction per body site ===\n")
print(reservoir %>% filter(genus == "Novel_any") %>% transmute(tissue, n_lib, novel_prop = round(prop, 3)) %>% as.data.frame())
