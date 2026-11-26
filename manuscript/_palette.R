## Sample-group colours, defined once and used by every figure that colours by
## population group.
##
## Before this file each figure picked its own: Figure 2 and Figure S1 took
## ggplot's default hues, Figure 3C took a hand-written red/blue/grey set, so a
## reader who learned "productive controls are red" in 3C met them as gold in
## Figure 2. Sourced by bin/hpv_rpm_contrast.R and manuscript/figure_orf.R.
##
## The palette is Okabe-Ito, which is colour-blind safe by construction, mapped
## onto the semantics the figures argue: warm for productive disease, cool for
## the quiet tissues, neutral for material that is not tissue at all. The point
## of using Okabe-Ito rather than shades of one hue is separation within each of
## those families. The old set distinguished productive controls from EV, and
## neoplasia from unselected skin, by lightness alone, which is exactly what
## fails in a dense scatter with overplotting.
GROUP_COLS <- c(
  # productive disease: warm
  control_productive    = "#D55E00",  # vermillion
  EV                    = "#E69F00",  # orange
  # quiet tissue: cool
  neoplasia             = "#0072B2",  # blue
  unselected_skin       = "#56B4E9",  # sky blue
  # controls and residue
  non_permissive_tissue = "#009E73",  # bluish green
  other                 = "#CC79A7",  # reddish purple
  contamination_study   = "#F0E442",  # yellow
  # not tissue: neutral
  cell_line             = "#000000",
  engineered            = "#999999")

## Draw order for scatter plots: the two big masses (unselected skin, cell lines)
## go down first so the small groups that carry the argument are not buried under
## them. Any group missing from a given figure is simply skipped.
GROUP_DRAW_ORDER <- c("unselected_skin", "cell_line", "other",
                      "non_permissive_tissue", "engineered",
                      "neoplasia", "control_productive", "EV")
