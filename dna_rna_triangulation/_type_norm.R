# Shared HPV type-label normalisation for the Shen/Logan comparisons.
#
# Lives in one file because it was duplicated before: shen_skin_dna_vs_rna.R was
# fixed on 2026-08-06 while shen_reservoir_by_tissue.R kept the broken copy, so
# the skin figure and the body-site figure disagreed about which genus HPV38B
# belongs to. Source this from both.
#
# PaVE names variants of a type with a trailing letter: HPV38B is a subtype of
# HPV38, HPV14D of HPV14, HPV57B of HPV57. Our typing panel carries the subtype
# label while Shen's `best_hit_name` carries the bare type, so without this the
# two never join, the subtype's libraries land in "Unresolved", and the type
# loses its genus.
#
# Strips a single trailing letter ONLY where it directly follows a digit, so
# novel candidate labels that legitimately end in letters (HPV-mKN1nr,
# HPV-mSK190nr) are left alone.
norm_type <- function(x) stringr::str_remove(toupper(x), "(?<=[0-9])[A-Z]$")
