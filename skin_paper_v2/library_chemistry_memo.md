# Library chemistry as a covariate: what the data say

Closes the "unflagged methodological gap" from `hpv_framing_review.md`: depth is
modelled throughout the manuscript, library chemistry is not, and poly-A
selection versus rRNA depletion is the obvious counterpart, since HPV late (L1)
transcripts are the ones poly-A capture recovers least well.

Result in one line: **chemistry is now measured, it does move trace detection,
and it does not touch any conclusion.** The productive-infection result is
strengthened, because the neoplasia tier turns out to have been built with the
*more* favourable chemistry.

## What was built

| file | what it does |
|---|---|
| `bin/fetch_library_selection.py` | pulls SRA RunInfo plus the full-XML free-text `LIBRARY_CONSTRUCTION_PROTOCOL` / `DESIGN_DESCRIPTION` for an existing accession list (no fresh search, so the cohort is unchanged) |
| `bin/classify_library_chemistry.py` | two-source, deliberately conservative chemistry call |
| `bin/library_chemistry_covariate.R` | composition tables, within-tier effect, headline contrast under five specifications |

Outputs: `results_reframe_fixed/library_{selection,protocol,chemistry}.tsv` and
`results_reframe_fixed/chemistry_{by_tier,within_skin,models}.tsv`.

## Why the SRA field alone was not enough

`LibrarySelection` is filed as the generic "cDNA" for 8,100 of 10,849 libraries
(75%). Resolving the rest needed the free-text construction protocol.

The proof that this matters is our own validation series: PRJNA612785 is filed
as `LibrarySelection = cDNA`, while its protocol text reads "Ribosomal RNA was
deleted". The methods statement that the cervical series is rRNA-depleted is
therefore correct, but it is not checkable from the field a reader would look
at. Running the classifier over those 15 libraries as a blind test recovers
rRNA-depleted for 15/15.

Coverage after both sources: **6,751/10,849 (62.2%)** resolved, against 18.8%
from the field alone. Libraries whose text supports both calls (103) or neither
are held at `unspecified` rather than guessed at.

| chemistry | n | % |
|---|---:|---:|
| polyA | 4,670 | 43.0 |
| unspecified | 4,098 | 37.8 |
| rRNA_depleted | 1,033 | 9.5 |
| template_switch (SMART) | 452 | 4.2 |
| random_primed | 304 | 2.8 |
| tenx_3prime | 292 | 2.7 |

## Finding 1: the tiers are not balanced, and the imbalance favours our claim

| tier | polyA | rRNA-depleted | unspecified |
|---|---:|---:|---:|
| unselected skin | 42.3% | **7.2%** | 40.2% |
| neoplasia (cSCC/AK) | 12.6% | **22.9%** | 60.6% |
| productive controls | 34.1% | **0%** | 65.9% |

The tier where we report no productive infection is three times more often
rRNA-depleted than unselected skin, and the tier where we report abundant
productive infection contains no rRNA-depleted libraries at all. The chemistry
imbalance runs against the negative result, not with it. That is the strongest
possible answer to the objection and it should be stated explicitly rather than
buried.

## Finding 2: poly-A versus rRNA depletion does not move HPV detection

Adjusted for depth and tier, restricted to libraries whose chemistry is known
and to the two tiers the headline contrast is drawn from (unselected skin and
neoplasia), with poly-A as reference. The tier restriction matters: chemistry
mix differs sharply between tiers (cell lines are half poly-A, EV is entirely
random-primed), so an unstratified fit confounds chemistry with tier.

| chemistry | OR for any HPV read | 95% CI | P |
|---|---:|---|---:|
| rRNA-depleted | **1.01** | 0.82 to 1.24 | 0.93 |
| random-primed | 1.97 | 1.45 to 2.68 | 1.6e-05 |
| SMART template-switching | 4.44 | 3.47 to 5.67 | 8.3e-33 |
| 10x 3'-barcoded | 3.63 | 2.64 to 4.98 | 1.5e-15 |

The axis the objection is about is flat. Across all tiers, tier-adjusted, the
rRNA-depleted term is modestly above one (OR 1.19, 1.00 to 1.40, P = 0.044),
carried by the non-cutaneous tiers; within the cutaneous population that the
analysis rests on it is flat (OR 0.94, 0.76 to 1.16 within unselected skin
alone). What does move trace detection, and moves it by three to four fold at
equal depth, is the amplification chemistries.
That is a contamination and artefact signal, not a recovery signal, and it sits
naturally with the paper's existing contamination argument.

Within unselected skin, trace positivity by chemistry (n studies in brackets):
poly-A 19.2% (122), unspecified 28.9% (139), rRNA-depleted 31.9% (18),
random-primed 47.1% (15), template-switch 50.9% (9), 10x 75.2% (12). Each level
spans many studies, so these are not single studies wearing a chemistry label.

## Finding 3: every conclusion survives, most get stronger

Neoplasia versus unselected skin, tier term only:

| specification | any-read OR | P |
|---|---:|---:|
| M2a depth only (published) | 1.95 | 5.5e-10 |
| + chemistry | 2.04 | 8.7e-11 |
| chemistry-known libraries only | 3.37 | 2.0e-11 |
| rRNA-depleted libraries only | 5.60 | 7.2e-12 |
| excluding 10x | 2.11 | 8.5e-12 |

Productive-floor detection: 14.36 (published) to 13.71 adding poly-A status.

**The decisive sensitivity analysis.** Restricting to the 657 rRNA-depleted
libraries (564 skin, 93 neoplasia from 3 studies, median 45.6M reads, deepest
96.8M): **not one library crossed the productive floor.** All 13 floor-crossing
libraries in the skin-versus-neoplasia contrast sit in the poly-A or unspecified
strata. Under the chemistry most favourable to viral recovery, at depth, the
productive signal is still absent from cutaneous neoplasia.

## One thing this exposed that is not about chemistry

3'-barcoded (10x) libraries are excluded from the E1:E2 marker
(`bin/e1e2_integration_marker.R`) but are **not** excluded from the RPM and tier
prevalence analysis. In unselected skin they are 258 libraries at a median 116M
reads carrying a 75.2% any-read rate, roughly four times the poly-A rate at
equal depth. They are almost certainly contributing barcode-swapping artefact to
the trace-signal denominator.

Dropping them does not change the conclusion (S3 above: OR 2.11 versus 1.95),
so this is a hygiene question, not a correction. But the exclusion is worth
making explicit and consistent, and the current asymmetry between the two
analyses will read as an oversight if a referee notices it.

## Wired into the manuscript

Done, and the manuscript renders clean.

- `manuscript/results_models.R`: new **MODEL 2c** block computes the chemistry
  composition, the chemistry effect, the re-fitted headline contrasts and the
  rRNA-depleted-only sensitivity, and exposes them as `stats$chem_*` /
  `stats$ribo_*`. It is guarded by a `tryCatch`, so the manuscript still renders
  where `library_chemistry.tsv` has not been built. It also builds `chem_tbl`
  for the supplement.
- `manuscript/_methods.qmd`: new subsection "Library-selection chemistry", with
  the two-source procedure and the 15/15 blind recovery of the validation series.
- `manuscript/_results.qmd`: new paragraph in the neoplasia section, placed
  directly after the depth argument it parallels and before the non-permissive
  tissue control.
- `manuscript/_discussion.qmd`: the one-line poly-A caveat in Limitations is
  replaced by the quantified version, which now states the residual limit
  honestly (37.8% unresolved; selection chemistry is not RNA integrity or input
  amount).
- `manuscript/_supplement.qmd`: new `@supptbl-chem`, chemistry composition by
  tier.

Still open, for the authors to decide: whether to exclude 3'-barcoded libraries
from the RPM analysis as they already are from the E1:E2 marker (see above).
