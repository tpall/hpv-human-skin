# Figure plan, 2026-08-27: trace decomposition in, Figure 5 rebuilt, vaginal demoted

Decided in discussion 2026-08-27. Two additions to paper 1, plus Marko's cut to
Figure 5. Follows the peri-lesional removal and the "tier" to "sample group"
rename made the same day.

## Context

Three candidate figures came off the collaboration slides
(`slides/hpv_collab_slaidid.qmd`):

| slide | figure | status |
|---|---|---|
| p15 | trace decomposition | genuinely new, draft text already written |
| p21 | DNA vs RNA triangulation | panel B claim already asserted in paper 1 prose, unillustrated |
| p23 | DNA reservoir tropism by body site | **already in paper 1 as Figure S3**, in a better version |

So p23 needs nothing added. What it triggered instead is Marko's point that the
vaginal arm does not belong in a skin paper's main figure.

### The paper 2 boundary, and how this settles it

`presence_vs_activity/PLAN.md` (2026-08-05) scopes paper 2 as "presence is not
activity": presence and productivity are decoupled, predictably, from genus
crossed with anatomical niche, with skin and cervix as the two poles. Its first
two results are exactly the figures under discussion:

- **R1**, DNA reservoir across body sites recovers textbook tropism (marked DONE)
- **R2**, skin: presence does not predict productivity

R1 has already crossed into paper 1 (Figure 5 exists, `@shenpv2026` cited), so
that boundary moved before today. This plan settles the rest as follows:

- Paper 1 takes **only the per-type scatter** (triangulation panel B), which is
  the visualization of a sentence paper 1 already carries.
- Paper 1 **gives R1 back** by demoting the body-site figure to the supplement,
  where it is a read-out control rather than a headline result.
- Paper 2 keeps the genus-shift panel (triangulation panel A), the whole mucosal
  pole, and the 6,396-library cervical run, which is where its novelty sits.
  Losing R2's skin half is survivable because paper 2 was never the skin result
  alone: it is skin contrasted with cervix, and a published, citable skin pole
  makes the methodological headline easier to land.

**Follow-up:** `presence_vs_activity/PLAN.md` R1/R2 need a note recording this,
and the memory note `project_paper2_presence_activity` still claims "no
Shen/Logan content is in `manuscript/`", which has been false since Figure 5
landed. Both should be corrected once this is applied.

---

## Part 1: trace decomposition to the supplement

Cheapest of the three placements: one supplementary figure, three sentences in
Results, one Discussion limitation rewritten. The argument lands in the text.

### Why it is worth doing at all

Paper 1's own Discussion currently *concedes* this ground. Limitation six
(`_discussion.qmd:63`) reads: "the same caveat attaches to distinguishing
low-level true infection from residual contaminant signal at the detection
floor." The decomposition answers it: of the below-floor HPV-positive libraries,
56% carry no reads on any HPV gene (a Kraken2 k-mer call that does not survive
alignment), 28% carry early-gene reads with the late programme absent, and 17%
reach the productive L1 threshold, against 91% productive in the warts/EV
controls. And the early-only libraries carry cutaneous types, not the HeLa
mucosal signature: 2 of 64 are HPV16/HPV18.

Note the draft calls this "limitation 4"; in the current text it is the sixth.

### Files

- `trace_decomposition/decompose_trace.py` produces the counts.
- `trace_decomposition/figure_trace_decomposition.R` draws the figure.
- `trace_decomposition/draft_manuscript_paragraph.md` already holds a drafted
  Results subsection and a drafted Discussion replacement. Use it as the source,
  cut to three sentences for Results.

### Work

1. **Fix the script's plumbing before trusting any number.**
   `decompose_trace.py` reads a hardcoded `/tmp/rpm.tsv` and writes
   `/tmp/trace_rows.json`, and hardcodes `FLOOR=13.2` at line 4. Point it at
   `results_reframe_fixed/hpv_rpm_per_sample.tsv`, read the floor from the tier
   table rather than pinning it, and write its outputs into
   `trace_decomposition/`. The floor value is currently correct, which is
   precisely why the hardcoding is dangerous.

2. **Settle which libraries the denominator covers.** The script's group at line
   48 is `unselected_skin, neoplasia, other`, but the draft paragraph describes
   "the 278 non-cell-line, below-floor HPV-positive libraries of skin/neoplasia".
   Those are different sets. Decide which one the paper claims and make script,
   figure and text agree. `other` is the residual sample group with no recorded
   site, so excluding it is the more defensible choice, but then the count must
   be recomputed rather than reused.

3. **Re-run and confirm the numbers are unchanged by today's edits.** They should
   be: the 15 excluded peri-lesional libraries were all above floor and were
   never in skin or neoplasia, and the control set is 44 + 10 = 54, which is
   `control_productive` plus `EV`, both untouched. Confirm 155 / 77 / 46 = 278
   and 2 / 3 / 49 = 54 rather than assuming it.

4. **Resolve the 16% versus 17% mismatch.** The figure rounds 0.1655 to 17%, the
   draft text says 16%. Pick one and use it in both.

5. **Rebuild the figure in English.** The current axis title is "Osakaal HPV+
   teekidest" and all three category labels are Estonian. Match the manuscript's
   figure style: no in-plot title or subtitle, since what the panel shows belongs
   in the legend (this is the convention already stated in
   `bin/hpv_rpm_contrast.R`).

6. **Place it as a supplementary figure directly after `@suppfig-depth`** in
   `manuscript/_supplement.qmd`. Both concern the detection floor, so they read
   as a pair. Quarto numbers float crossrefs by document order and everything
   references by label, so the existing `@suppfig-types` and
   `@suppfig-reservoir-gut` renumber themselves with no manual edits.

7. **Three sentences in `_results.qmd`.** Insert after the "Crucially, that
   marginal excess is not productive burden" paragraph and before the chemistry
   paragraph ("Nor is it explained by how the libraries were built"). That is
   where the trace has just been established and before the argument moves on.
   One sentence for the three-way split, one for the controls as mirror image,
   one for the cutaneous-not-HeLa point. Add a matching single sentence to
   `_results_short.qmd`.

8. **Rewrite limitation six in `_discussion.qmd:63`.** The draft's proposed
   replacement is usable nearly as written. Keep the cross-mapping half of the
   sentence, replace the concession half.

---

## Part 2: Figure 5 rebuilt, skin only

### Target shape

**Figure 5 (main), skin only, two panels:**

- **A**: the ten commonest skin DNA types, with their RNA assigned and trace
  counts annotated. This is the current Figure 5 panel B, skin facet only.
- **B**: per-type DNA abundance against RNA productivity, log-log, coloured by
  genus, with the empty "abundant AND productive" corner and HPV2 as the sole
  occupant. New, lifted from the triangulation.

**Figure S3 (supplement): unchanged.** It is already the three-site version
(skin, vaginal, gut) with both panels. Its *job* changes: it stops being "Figure
5 with gut added" and becomes the read-out control showing that the composition
measurement recovers textbook tropism from DNA alone. Only its caption is
rewritten.

That is the whole of Marko's cut: no vagina in the main figure, the validation
argument preserved one level down.

### Why this is a small edit

`manuscript/figure_reservoir.R` already reads
`dna_rna_triangulation/skin_dna_vs_rna_types_corrected.csv` as `cmp`, and that
file already carries both axes the new panel needs: `dna_skin_libs` and
`rna_productive`, plus `genus`. So panel B can be built inside
`figure_reservoir.R` from data it already loads. No cross-script surgery, and
`dna_rna_triangulation/shen_skin_dna_vs_rna.R` keeps producing the two-panel
triangulation figure for paper 2 unchanged.

Panel B's construction can be copied from that script (lines 163 to 176): the
`lg()` log1p transform, the shaded `annotate("rect", ...)` corner, `geom_jitter`,
`geom_text_repel` on the union of the top DNA and top productive types, and the
shared `gcol` genus palette, which `figure_reservoir.R` already defines
identically.

### Work

1. **`manuscript/figure_reservoir.R`**: keep `make_fig()` as it is for the
   supplement, and give `figure5` its own builder. `SITE_SETS` collapses to
   `figureS3` only; figure 5 no longer comes from a site set.
2. **`manuscript/reservoir_numbers.tsv`**: add the keys panel B's caption needs,
   written by the same script so figure and caption cannot drift. At minimum:
   how many types have DNA presence but zero productive libraries, and the HPV2
   exception. The existing `rv()`/`rvn()` helpers in `_helpers.R` pick them up
   with no change.
3. **`_results.qmd`, reservoir section (line 72 onward)**:
   - The vaginal validation sentence ("applying the identical procedure to
     vaginal metagenomes ... recovers the textbook mucosal niche from DNA
     alone") now points at `@suppfig-reservoir-gut` instead of standing inside
     `@fig-reservoir`. Merge it with the existing gut sentence, since both are
     now describing the same supplementary figure.
   - The `(@fig-reservoir B)` reference for the top-ten annotations becomes
     panel A.
   - Add a sentence for panel B. The paragraph already opens "The types are
     where the two modalities meet, and they do not agree", so the new sentence
     is its illustration, not a new claim.
   - Rewrite the `@fig-reservoir` caption for the new panels, and the
     `@suppfig-reservoir-gut` caption for its new job.
4. **`_results_short.qmd`**: same two adjustments, shorter.
5. **`_discussion.qmd:45`** cites `(@fig-reservoir)` for HPV8 being eighth
   commonest in the reservoir. Panel A retains the top ten skin types, so the
   reference still resolves. Confirm rather than assume.

---

## What is deliberately not being done

- **Triangulation panel A** (genus composition across DNA skin reservoir, RNA
  skin present, RNA productive) stays out of paper 1. It combines two things
  paper 1 already shows separately, in Figure 4 and in the reservoir section
  text, and it is the cleanest statement of paper 2's framing.
- **Promoting S3 to the main text.** The opposite is happening.
- **Any change to `tier` as a variable or column name.** Prose and display
  strings only, as settled earlier today.

---

## Rebuild sequence

Run from the repo root, in order.

```bash
# 1. refresh the DNA/RNA comparison against the current type catalogue
Rscript dna_rna_triangulation/shen_skin_dna_vs_rna.R
Rscript dna_rna_triangulation/skin_dna_vs_rna_corrected.R

# 2. trace decomposition, after the script plumbing is fixed
python3 trace_decomposition/decompose_trace.py
Rscript trace_decomposition/figure_trace_decomposition.R

# 3. figure 5 + figure S3 + reservoir_numbers.tsv
Rscript manuscript/figure_reservoir.R

# 4. text, then render
quarto render manuscript/manuscript.qmd
quarto render manuscript/manuscript_short.qmd
```

`bin/rebuild_reframe_fixed.sh` does **not** need re-running: the sample-group
table is unchanged by any of this.

Note that step 1 matters even though the numbers should not move.
`shen_skin_dna_vs_rna.R:101` reads
`results_reframe_fixed/type_catalogue_by_tier.tsv`, which was regenerated today
when the peri-lesional column was dropped and the stray `NA` column fixed. It
takes `unselected_skin` and `control_productive + EV`, none of which changed, so
the expectation is an identical figure. Verify that, do not assume it.

---

## Verification

1. **The scatter agrees with the prose it illustrates.** The types named in the
   Results sentence "the ten commonest skin types in the DNA reservoir are ..."
   must be the same ten labelled in panel A, and none of them may sit in panel
   B's shaded corner except HPV2.
2. **Figure 5 has no vaginal or gut data anywhere**, and Figure S3 still has all
   three sites with the skin RNA annotations intact.
3. **The trace decomposition percentages in the figure, the Results sentences
   and the Discussion rewrite are the same three numbers**, at the same
   rounding.
4. **Denominators are stated and consistent**: the below-floor set in the
   caption matches whatever was decided in Part 1 step 2, and the control set is
   named as productive controls plus EV, n = 54.
5. **Both manuscripts render with no `NA`, `NULL` or `Inf` in the prose**, and
   the supplementary figures number consecutively with every `@suppfig-`
   reference resolving.
6. **No em-dash characters** in any edited file, since the git hooks reject them.
