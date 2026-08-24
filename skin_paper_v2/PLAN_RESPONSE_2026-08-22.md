# Marko's structure, re-answered against the revised manuscript

Written 2026-08-22, revised the same day after checking the type-level numbers
against the data (the HPV8 claim in §4 was wrong; see there), and again on
2026-08-24, once the type grid stopped being a sketch and became a manuscript
figure. Supersedes the recommendation in
`PLAN.md` §5 only in its *reasons*; the recommendation itself is unchanged. The
proposal Marko has re-sent is verbatim the list answered in `PLAN.md` §2, but the
manuscript it would be folded into is not the same manuscript. Thirty-two commits
since 2026-08-07 rewrote the framing, added 30 references, added a
library-chemistry covariate, and added two figures. Two of those changes move the DNA question, and they move it in
opposite directions.

## 0. What we are adopting from your proposal

Short version, because the rest of this document is detail. Of the six analysis
blocks proposed, five are in or going in, and the one that is not is the one the
data cannot support at all.

| block | status |
|---|---|
| **R3, by ORF** | **Built. `manuscript/figure3.png`, three panels. See §3.** |
| R1, frequencies by sample type | In the manuscript, and stronger than in August: the contrast now survives adjustment for library chemistry |
| D1B, skin HPV types | Have the data; the palette and legend fix you asked for is still open |
| **R2A, types by sample type** | **Built as a figure. `manuscript/figureS2.png`, the type-by-tier grid. See §4** |
| D1A, genus composition | Doable now; the absolute-frequency half is blocked on the DNA denominator, not on the idea |
| R2B, co-infection | In the supplement as Table S7; promotable to main text if the hedge travels with it |
| D2, E2 ORF absence | Not doable, and the reason it fails became a section of its own (§5) |

The ORF analysis is the one worth talking about first. It was your idea, it was
the cheapest thing on the list, and it turned out to be the figure the manuscript
was most obviously missing: the paper's entire central criterion is the
early-versus-late transcript programme and there was no picture of it anywhere.
That figure now exists and it does more work than expected, so §3 is about it
rather than about whether to build it.

On structure, the proposal is for a different skeleton. What is being suggested
instead is that the blocks land inside the current one, which they do without much
strain, and that we talk through §7 in person.

## 1. What changed in the manuscript, and why it matters here

Three revisions bear directly on Marko's list.

**The activity result is now explicitly not novel.** The Discussion states, in
print, that transcriptional silence of HPV in cutaneous carcinoma has been
reported from four independent RNA-seq series (Arron 2011, Hultin 2019, Arroyo
Muhr 2021, Vandiver 2021) plus an orthogonal ISH cohort of 243 (Fujimoto 2019),
and that the depth argument "is no longer the discovery that depth was not the
explanation." What the paper now claims is the denominator, the controls and the
productive-replication criterion. That reframing is correct and it is stronger
than what it replaced, but it changes the calculus on the DNA arm: the paper's
novelty is now carried by *design*, and a second modality is a design argument.

**The paper now argues DNA-versus-RNA at length, using entirely other people's
DNA.** The new section "What a 2-4% transcriptome prevalence means against a 90%
DNA reservoir" rests on Weissenborn 2012 (viral copies per cell), Rollison 2019
(swab versus eyebrow-hair recovery), Ma 2014 (body-site metagenomic prevalence)
and Bzhalava 2014 (types with no reference genome). This is the single largest
change in the value of Marko's D1: the paper has committed to a quantitative
DNA/RNA reconciliation as a load-bearing section, and it currently sources every
number in it externally. Our own reservoir arm would convert that from a
literature argument into internal evidence.

**The paper has committed, twice and in print, to substrate as the explanation.**
The same section says "Our substrate is neither" of swabs and hair bulbs, and the
new HPV8/AK section says the Lrig1+ follicular junctional zone "is exactly the
compartment that our substrate under-samples." This is the inversion. In
`PLAN.md` §3.3 the position was that cohort size makes the swab-versus-biopsy
separation unlikely to matter. The manuscript now argues the opposite as one of
its explanations for its own headline number, and it cannot have this both ways.

## 2. The consequence for C1, stated plainly

`PLAN.md` §3.3 treated the cross-compartment problem as a stated limitation that
cohort size mitigates. That is no longer available. The DNA reservoir is 100%
surface swabs (0 of 941 skin BioSamples are biopsies); the RNA is bulk biopsy
tissue; and the manuscript now explains its own 2-4% by exactly that difference.
A "presence versus activity, per type" figure built across those two axes would
be read, correctly, as comparing two compartments and labelling the difference
activity.

The fix is not to drop C1 but to relabel it as what it is: **surface and
follicular DNA against intra-epidermal RNA**. Reframed that way it stops being a
confounded activity comparison and becomes a direct display of the paper's own
thesis, which is that the reservoir and the transcriptionally detectable
compartment are different things. That is a better figure than the one proposed,
and it is the only version consistent with the Discussion as it now stands.

**Decided (2026-08-22): substrate stays in the Discussion and C1 is relabelled as
compartments.** This reverses the annotation on `PLAN.md` §3.3, which held that
cohort size made the swab-versus-biopsy separation unlikely to matter. The two
positions are narrower than they look and both survive: that annotation is about
whether the compartments differ enough to *produce* the observed separation, which
the manuscript does not claim; what the manuscript claims is only that a
prevalence measured on swabs and one measured on biopsy tissue are not the same
quantity. The reconciliation section rests on three legs, abundance (Weissenborn),
substrate (Rollison) and the panel ceiling (Bzhalava), so it would survive dropping
substrate. It is kept because the abundance argument alone leaves open why a
near-ubiquitous reservoir yields 2-4%, and substrate is the concrete answer.

## 3. The ORF analysis, built

This was your item, it was the cheapest thing on the list, and it is now built:
`manuscript/figure_orf.R` writes `manuscript/figure3.{png,pdf}` plus
`orf_breadth_by_tier.tsv` and `orf_composition_by_tier.tsv`. It is worth leading
with, for a reason that was not obvious in August. The paper's central criterion
is the early-versus-late transcript programme, and until now there was no figure
of it anywhere: the whole argument rested on a threshold the reader had never been
shown. The Discussion has since also taken on a third model (Rehm 2022, in which
beta-HPV restricts its own expression through E8^E2), which is a claim about which
viral genes are transcribed, so an ORF-level panel now has two jobs rather than one.

*Since this was drafted* the panel titles have moved out of the image and into the
legend, and the figure renumbered from 4 to 3, so the version below is current.

No new computation was needed. The per-gene counts come from the 541
`*_transcript_classes.tsv` files the pipeline already writes; 488 of those join to
a tier: unselected skin 217, cell line 98, other 60, productive controls 44,
neoplasia 26, engineered 18, peri-lesional 15, EV 10.

### Figure 4: the viral transcriptional programme by population tier

![](../manuscript/figure3.png)

**Panel A, programme breadth.** Percentage of each tier's HPV-positive libraries
carrying at least one read on each ORF, as a Beta-Binomial posterior mean under
Beta(0.5, 0.5), the same treatment as the type tables:

| tier | E1 | E2 | E6 | E7 | L1 | L2 |
|------|---:|---:|---:|---:|---:|---:|
| peri-lesional (n = 15) | 97 | 97 | 97 | 97 | 97 | 97 |
| EV (n = 10) | 86 | 96 | 86 | 86 | 96 | 86 |
| productive controls (n = 44) | 90 | 94 | 86 | 90 | 90 | 90 |
| neoplasia (n = 26) | 57 | 65 | 50 | 68 | 50 | 54 |
| unselected skin (n = 217) | 22 | 23 | 10 | 14 | 19 | 21 |
| cell line (n = 98) | 53 | 34 | 90 | 93 | 36 | 11 |
| engineered (n = 18) | 66 | 18 | 71 | 76 | 40 | 13 |

Three shapes fall out of one panel. The productive tiers run the complete
programme, capsid genes included. Skin and neoplasia are uniformly low across the
whole genome, which is a trace signature rather than a truncated programme, and
that distinction is the paper's central claim made visible. The cell-line and
engineered rows are the canonical integrated pattern read straight off the
transcripts: oncogenes at 90% and 93%, E2 collapsed to 34% and 18%, L2 down to 11%
and 13%. That last row generalises the E1:E2 marker to the whole genome and makes
the contamination argument as a picture rather than as an odds ratio.

**Panel B, programme composition.** Mean share of a library's viral reads falling
on each ORF, stacked in genome order. Collapsed to the two groupings that matter:

| tier | late (L1+L2) | oncogenes (E6+E7) |
|---|---:|---:|
| productive controls (n = 40) | 34.9% | 6.4% |
| EV (n = 9) | 23.7% | 9.8% |
| peri-lesional (n = 15) | 18.0% | 3.2% |
| unselected skin (n = 36) | 27.8% | 9.3% |
| neoplasia (n = 18) | 11.7% | 21.6% |
| cell line (n = 87) | 2.2% | 66.6% |
| engineered (n = 12) | 3.6% | 63.2% |

The productive-versus-integrated inversion is complete and unambiguous: warts and
EV put a third of their viral reads into capsid genes and almost nothing into
oncogenes, and the cell lines do exactly the reverse. Neoplasia sits between,
shifted toward the oncogenes relative to skin (21.6% against 9.3%) and away from
late transcription (11.7% against 27.8%), which is the direction the paper argues
but on 18 libraries, so it is worth stating as suggestive and not more.

**Panel C, the criterion itself.** Early-region reads against L1 major-capsid
reads, per library, log-log, coloured by tier, with the L1 >= 3 productive
threshold drawn as a line. The y axis is L1 rather than the whole late region so
that the line is exactly the rule the paper applies rather than an approximation
of it. If only one panel survives review, this is the one.

Bayesian treatment, as you asked for, is what Panel A already uses and fits it
naturally, since each cell is a proportion.

### What building it turned up

Two things a proposal would have missed, both now fixed in the script.

The transcript-class tables **double-count 53 libraries** processed in more than
one results directory. Every duplicated (library, gene) pair carries an identical
count, so de-duplication is lossless, but without it those libraries take twice
their weight and the composition panel does not sum to 100%.

Composition is only estimable where there are reads to apportion, and **unselected
skin has a median of zero viral reads** across these six ORFs, against 60,963 in
the productive controls. On a per-library mean-of-fractions a single L2 read would
read as 100% late transcription, so an unrestricted Panel B would appear to show
abundant late transcription in exactly the tiers the paper reports as quiet. Panel
B is therefore restricted to libraries with at least 20 viral reads, with the
restriction and the surviving n printed on the panel. Only 36 of 217 skin
libraries qualify, which is itself worth saying out loud.

Two further things are built in rather than bolted on. E5 is dropped from the
display, since beta and gamma types largely lack it and it is zero in every tier.
And every panel is **conditional on HPV positivity and successful typing**: this
is a figure about the programme of detected infections and must never be read as a
prevalence statement. Panel A carries the prevalence story on its own.

## 4. The DNA reservoir against RNA activity, per type

*This section replaced an earlier version that made the same argument about HPV8
alone. Checked against the data, the HPV8-specific claim does not hold; the
general one is stronger. The corrected table is built by
`dna_rna_triangulation/skin_dna_vs_rna_corrected.R` and written to
`skin_dna_vs_rna_types_corrected.csv`. The RNA side of that comparison has since
been drawn as a manuscript figure, `manuscript/figureS2.png`, which is also the
answer to R2A and is shown at the end of this section.*

The manuscript's HPV8/AK section (added in response to Morgan 2024, which builds a
hit-and-run stem-cell mechanism on HPV8 in 70% of human AK) ends on an open
question: we detected no HPV8 anywhere in 10,849 RNA libraries, and cannot say
whether that is genuine absence, cohort composition, or DNA-versus-RNA
discordance. The DNA arm speaks to it, but not in the way first drafted here.

Two corrections had to be made before the numbers meant anything. The first is
the `rna_skin_present` threshold artefact already flagged in `PLAN.md` §3.2 and
still unfixed at the time of writing: that column reads the assigned catalogue
only, so it reports zero for types the detection tier does recover. Rebuilt with
the trace tier included, **15 types the old column called absent from skin RNA
are in fact recovered**, among them HPV12, HPV20, HPV38, HPV107 and HPV47. The
second is scope: a claim that a type is absent from our RNA has to be made
against the whole cohort, not against the unselected-skin tier alone.

With both fixed, the top of the skin DNA reservoir looks like this. DNA counts are
libraries out of 1,735 HPV-positive skin metagenomes; RNA counts are libraries out
of 10,849, assigned and trace tiers combined.

| rank | type | genus | DNA libs | RNA, skin | RNA, whole cohort |
|-----:|------|-------|---------:|----------:|------------------:|
| 1 | HPV24 | Beta | 128 | 0 | 1 |
| 2 | HPV12 | Beta | 96 | 3 | 7 |
| 3 | HPV38 | Beta | 94 | 1 | 2 |
| 4 | HPV20 | Beta | 81 | 3 | 7 |
| 5 | HPV107 | Beta | 80 | 1 | 4 |
| 6 | HPV173 | Gamma | 80 | 0 | 0 |
| 7 | HPV21 | Beta | 78 | 0 | 0 |
| 8 | **HPV8** | Beta | **63** | **0** | **0** |
| 9 | HPV5 | Beta | 62 | 0 | 3 |
| 10 | HPV75 | Beta | 61 | 0 | 0 |

**HPV8 is not specifically silent.** Four of the ten most prevalent skin
reservoir types have zero RNA evidence anywhere in the cohort: HPV173, HPV21,
HPV8 and HPV75. And HPV5, at essentially the same DNA prevalence as HPV8 (62
libraries against 63), appears in RNA three times out of 10,849, which is not
distinguishable from zero. A hero panel built on HPV8 invites the question "what
about HPV5?", and the answer deflates it.

The general result is the one to lead on, and it is a better one:

> Across the ten most prevalent types in the cutaneous DNA reservoir, each
> carried by 61 to 128 of 1,735 HPV-positive skin metagenomes, RNA recovery
> across our entire 10,849-library cohort ranges from 0 to 7 libraries per type.

That is the reservoir-versus-activity thesis stated at type resolution, on our own
data on both axes, and it does not depend on any single type behaving specially.
HPV8 then becomes one labelled point on that panel which happens to answer the
open question in the AK section: it is in the cutaneous DNA reservoir at a rate
comparable to the other beta-1 types and absent from the transcriptome, like three
of its neighbours. That is a caption, not the claim.

The compartment argument still corroborates, at the level it actually operates.
The Lrig1+ follicular junctional zone is sampled by swabs and hair bulbs and
diluted by biopsies, so DNA-positive-swab against RNA-negative-biopsy is what the
manuscript's own §2 argument predicts. But that prediction applies to every beta
type in the reservoir, not to HPV8 in particular, which is exactly why the general
framing and the mechanism fit each other and the HPV8-specific framing does not.

One thing this cannot settle, and the panel must not imply otherwise: the two arms
share no libraries. This is a composition comparison between two disjoint cohorts,
not a within-sample presence-versus-activity measurement. No paired skin specimens
exist in public data (`PLAN.md` §3.3), so the panel shows that the reservoir and
the transcriptome have different type compositions, and cannot show that any
individual genome is present-but-silent.

### The figure that came out of it

![](../manuscript/figureS2.png)

`manuscript/figureS2.png`, built by `manuscript/figure_type_grid.R`. It is the
answer to R2A, and it is in the supplement rather than the main text for a
reason worth stating, because it is the same reason the earlier draft of this
section was wrong. Of 66 types, 25 are seen exactly once. Thresholding the grid
to fit a main-text panel at three or more libraries would keep 86% of the
library calls but only 42% of the types, which removes the diversity evidence
and preserves the totals: precisely backwards for a claim about a reservoir that
is diverse and individually rare.

Three design points. Marko's palette complaint on D1B was that the legends do
not merge, which is the same bug this figure hit and fixed, so the fix transfers
straight to the DNA figure. Fill is prevalence *within* tier, not the count, because
tier sizes span 10 to 7,800 and three libraries is 23% of the EV tier against
0.04% of unselected skin; the raw count stays as the cell label so the figure
still ties back to the tables. The scale is log, because a third of non-empty
cells carry three or fewer libraries and a linear ramp puts them at the white
end where the labels vanish. And it is laid out landscape rather than portrait
so each tier is one line, which is what makes the unselected-skin line of the
lower panel legible as the reservoir claim rather than something to be counted
down a column.

The complementarity is the point. The upper panel, assigned types, is
concentrated in three rows, productive controls, peri-lesional mucosa and cell
lines, and almost entirely within alpha. The cutaneous beta and gamma types of
healthy skin appear only in the lower panel, because they sit below the coverage
threshold a type assignment requires, not because they are absent.

## 5. Block-by-block, updated

**D1A: HPV+ frequency and genus composition.** The denominator problem (§3.1) is
still unfixed and the instruction still stands: do not quote a percentage. But the
block is not blocked. Split it. The *genus composition of HPV-positive libraries*
(skin beta 57.9%, gamma 40.7%, alpha 7.4%, mu 1.8%, novel 14.5% of n = 1,735) is
conditional on positivity and needs no denominator at all; it can be modelled and
published now. The *absolute HPV+ frequency* needs the run-level rebuild with the
Logan Dec-2023 freeze on both numerator and denominator, using the
count-only-strategies-where-positives-occur approach. Marko's yes to a Bayesian
proportion model is right and is consistent with the manuscript, which already
carries a Bayesian re-analysis section.

**D1B: skin HPV types.** HAVE, and note this is the DNA side: `reservoir_by_tissue_top_types.csv` and `figure_reservoir_by_tissue.png`, not the RNA type grid of §4, which answers R2A. Agree on the palette and the colour/fill merge.
Agree on moving gut to the supplement, and the genus table supports the artefact
reading independently: gut returns beta at 34.9% and alpha at 30.5%, which is not
a credible gut virome and is what reagent or index contamination looks like.

**D2: E2 ORF absence.** Infeasible, as established, and already in the text. But
"we couldn't do it" undersells what happened, and Marko should be told the
stronger version: the reason it fails *is* the paper's most defensible novel
contribution. The dedicated prior-art sweep found the assembly-bias critique
(integrated genomes cannot assemble across the virus-host junction; ORF-truncation
frequency tracks distance from the L1 anchor rather than breakpoint position)
unclaimed in the literature. D2 did not fail, it turned into the section that
bounds what any L1-anchored petabase catalogue can claim.

**R1: frequencies by sample type.** HAVE, and now stronger than on 2026-08-07:
the headline contrast survives adjustment for library chemistry, and the decisive
sensitivity analysis (657 rRNA-depleted libraries, median 45.6M reads, none
crossing the productive floor) is in the paper.

**R2A: types by sample type.** HAVE, and now drawn: `figureS2` is exactly this block, the assigned catalogue against the detection tier as one grid (§4).

**R2B: co-infection.** Still not promoted; it remains supplement-only
(`@supptbl-coinf`, one reference in the Results). One caution before promoting it,
because the Results already hedge it hard and for good reasons: the co-infected
libraries are heavily pseudo-replicated, the secondary type is usually near the
detection threshold, and the strict-coverage set collapses to a handful of
libraries from a few donors. Promote it to main *text* with the hedge intact, or
to a supplementary figure. A headline co-infection panel would overstate what the
data support and would contradict the paper's own stated caveat.

**Protein analysis.** Closed, and stays closed. But it is worth sending Marko
`EORFEOME_SCREENS.md`, because the question has a partial answer from another
direction: mucosal high-risk E6 scores across p53 and NF-kB in a pooled effector
screen, HPV5 and HPV8 E6/E7 score in none of five pathways, and the single
beta-genus hit is HPV9 E6 on cGAS-STING. That is not protein evidence about our
samples and must not be presented as closing the block, but it does mean beta E6/E7
transcription should not be read as implying the mucosal oncogenic programme.

**C1: DNA/RNA comparison.** Keep, relabel as compartments (§2, now decided), and
lead on the per-type reservoir-versus-activity panel of §4 rather than on any single
type. The `rna_skin_present` threshold artefact from `PLAN.md` §3.2 is **fixed**:
`skin_dna_vs_rna_corrected.R` rebuilds the table with the trace tier included and
writes `skin_dna_vs_rna_types_corrected.csv`. Quote nothing from the old
`skin_dna_vs_rna_types.csv`, which called 15 types absent from skin RNA that the
detection tier recovers. The panel must be captioned as a composition comparison
between two disjoint cohorts, not a within-sample presence-versus-activity test.

## 6. Recommendation

Unchanged: **Option A**, fold the reservoir into the current manuscript as
supporting evidence rather than restructure around it. The reasons have shifted,
though. On 2026-08-07 the argument for A was that the RNA-only framing was tight
and adding DNA would import confounders. Half of that is now obsolete: the paper
has voluntarily taken on the DNA/RNA argument in its Discussion and sources it
entirely from the literature, so the confounders are already in the paper and the
only question is whether we answer them with our own data or someone else's.

What still argues against Option B is different from what argued against it
before. It is no longer that DNA imports the compartment problem, since the
manuscript has already conceded and used that problem. It is that a DNA-first
spine would make the reservoir the claim, and the reservoir claim needs the
denominator that §3.1 says is not yet defensible. Option A needs the genus
composition, which is conditional on positivity and needs no denominator. The
blocked number is load-bearing for B and optional for A.

Concretely, for Option A: one figure (D1A genus composition + D1B types as
panels), one Results subsection, the ORF panel from R3, HPV8 written into the
existing AK section, C1 relabelled as compartments, and the co-infection table
promoted to main text with its hedge. That lands five of Marko's six blocks, drops
only the infeasible one, and does not spend the manuscript's completeness.

## 7. Open, for the group

1. Option A versus B (§6). Unchanged from `PLAN.md`, still needs a decision.
2. The DNA denominator (§3.1). Blocking for B, optional for A.
3. ~~Substrate~~ **Decided (§2)**: substrate stays in the Discussion and C1 is
   relabelled as compartments. This reverses the annotation on `PLAN.md` §3.3 and
   the group should know that, though the two positions are compatible once the
   manuscript's narrower claim is stated (§2). Nothing further is blocked on it.
4. 3'-barcoded (10x) libraries are excluded from the E1:E2 marker but not from
   the RPM and prevalence analysis, where in unselected skin they are 258
   libraries at median 116M reads with a 75.2% any-read rate. Dropping them does
   not change any conclusion, but D1A and C1 would be built on that denominator,
   so the inconsistency is worth resolving before, not after.
