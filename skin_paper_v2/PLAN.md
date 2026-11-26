# Skin paper, Marko's structure: a figure-by-figure plan

Planning only. **The current manuscript in `../manuscript/` is not touched by
this.** Written 2026-08-07 in response to Marko's proposal to make one paper from
the skin data and leave mucosa for a second.

---

## 1. The decision this document is for

Marko's list is DNA-first: DNA reservoir, then integration, then RNA activity,
then a DNA/RNA comparison. The current manuscript is RNA-only and RNA-first. That
is not a rearrangement, it is a different paper with a different question:

| | current manuscript | Marko's structure |
|---|---|---|
| question | *Is HPV transcriptionally active in cSCC/AK, or only a footprint?* | *What cutaneous HPV is there, and is it doing anything?* |
| evidence | one modality (RNA), depth-controlled | two modalities (DNA + RNA) |
| the silence is | the **result** | the **punchline**, after presence is established |
| contamination | a co-equal result and methodological contribution | a filtering step |
| title now | "Cutaneous HPV is transcriptionally quiet in skin and carcinoma but productive in warts: a depth-controlled re-mining of public RNA-seq" | would need rewriting |

Both are publishable. The choice is not about which figures exist; it is about
whether the paper's claim is *about activity* (current) or *about the reservoir*
(Marko's). Section 5 gives my recommendation.

---

## 2. Figure-by-figure, with data status

Status key: **HAVE** = built and in the repo; **CHEAP** = one script, hours;
**NEEDS RUN** = new compute; **INFEASIBLE** = tested and cannot be done.

### DNA analysis 1

**Fig D1A. HPV+ frequency (%) across all skin DNA libraries + genus composition of the reservoir** (slide 19)

- Genus composition: **HAVE** (`dna_rna_triangulation/reservoir_by_tissue_genus.csv`;
  skin beta 57.9%, gamma 40.7%, alpha 7.4%, novel 14.5%, n = 1,735 HPV+ libraries).
- HPV+ *frequency*: **PROBLEM, see 3.1.** We have a numerator but the denominator
  is not settled, and my last three attempts at it were each wrong. Currently the
  most defensible internally consistent figure is 1,883/12,118 = **15.5%** of
  skin shotgun metagenomes, but this needs one careful pass before it is quoted.

**Fig D1B. Skin HPV types** (slide 23)

- **HAVE** (`reservoir_by_tissue_top_types.csv`, `figure_reservoir_by_tissue.png`).
- Corrected 2026-08-06 for the subtype-letter join, which moved HPV38 from 32 to
  94 libraries and pulled HPV14 and HPV57 out of the "no DNA reservoir" list.
  Any per-type ranking produced before that date is invalid.

### DNA analysis 2

**Fig D2. Which HPV types lack a detectable E2 ORF (integration from DNA)**

- **INFEASIBLE.** This was built, run, and is negative
  (`dna_rna_triangulation/REPORT_e2_integration_negative.md`, Table S7 in the
  current manuscript). Three independent framings fail for one reason: an
  integrated genome cannot assemble across the virus-host junction, so
  conditioning on near-complete contigs selects episomal genomes in every
  category. Cervical carcinoma, where integration is the rule, returns 99.3%
  intact E2 against 99.9% in skin.
- **The replacement is already done and is stronger:** integration read from RNA
  (E1:E2 ratio, AUC 0.968) and confirmed from reads (virus-host junctions,
  recovering HeLa's 8q24 site unprompted, every breakpoint in E1 or E2). That is
  Figure 2 of the current manuscript.
- Keep D2 in the paper as a **negative methods result** (Table S7), not as a
  planned analysis.

### RNA analysis 1

**Fig R1. HPV frequency by sample type** (skin, cell lines, papillomas, EV, skin cancer) (slide 10)

- **HAVE.** This is `figure1.png` / `@fig-rpm` in the current manuscript, 10,849
  libraries, depth-normalised RPM with a control-calibrated floor, plus the
  non-permissive-tissue and engineered tiers.

### RNA analysis 2

**Fig R2A. Types by sample type** and **Fig R2B. Co-infection**

- Types: **HAVE** (`@supptbl-types` assigned catalogue, `@supptbl-trace`
  detection tier, and `figure3.png` collapsing both to genus).
- Co-infection: **HAVE but unused in the text** (`build_coinf_table()`,
  `@supptbl-coinf`). Marko is right that it deserves promotion; it is currently
  a supplementary table nobody reads.

### RNA analysis 3

**Fig R3. By ORF**

- **HAVE** (`results_e1e2/e1e2_by_library.tsv` carries per-library E1, E2, E6,
  E7, L1, L2 counts and late fraction for every typed library).
- Currently used only for the E1:E2 integration marker. A per-tier ORF profile
  (early vs late programme by tier) is **CHEAP** and would be a genuinely new
  panel: it shows the productive tiers running the full early-to-late programme
  and neoplasia stopping at early.

### Protein analysis

- **DROPPED, and the question is closed.** Marko confirmed (2026-08-07) that he
  meant proteomics, and that he has already searched for suitable datasets and
  found none. That matches our position: RNA-seq gives transcripts, not protein,
  and no public cutaneous dataset carries mass spectrometry at the scale this
  would need.
- The nearest available thing is R3, which reports which ORFs are transcribed
  per library. It is not a substitute for protein evidence and should not be
  presented as one.
- Consequence for the paper: **two of Marko's blocks cannot be done**, this one
  and DNA analysis 2. Both should be stated in the Discussion as measurement
  limits rather than omitted silently, since both are the obvious next questions
  a reader will ask. The absence of protein evidence is worth one sentence:
  transcription is the closest available proxy for activity, and the paper's
  claims are therefore about transcription.

### DNA/RNA comparison

**Fig C1. Presence against activity, per type** (slide 21)

- **HAVE** (`figure_dna_rna_triangulation.png`).
- Carries two caveats that must be fixed or stated (3.2, 3.3).

---

## 3. What must be fixed before any of this is submitted

### 3.1 The DNA prevalence denominator
Three attempts, three different answers, because the denominator definition
moves it 14-fold and because `esearch` counts experiments while `efetch runinfo`
returns runs. The stratified approach (count only the strategies where positives
actually occur; `OTHER` yields 0 of 19,525 and is therefore 16S) is the right
method, but it must be redone entirely at run level, in one script, with the
Logan Dec-2023 freeze applied to both numerator and denominator. Until then, **do
not quote a percentage.**

### 3.2 "Absent from skin RNA" is a threshold artefact
`rna_skin_present` in the comparison table reads the **assigned** catalogue only,
so beta/gamma types show as zero when the detection tier in fact has HPV49 in 5
libraries, HPV20 in 4, HPV12 in 3. Either add the trace tier to that column or
relabel it *assigned-level presence*. **Do not describe those types as absent
from skin RNA.**

### 3.3 The DNA and RNA axes are different compartments
The DNA reservoir is 100% surface swabs (0 of 941 skin BioSamples are biopsies);
our RNA is bulk tissue. Papillomavirus late transcription completes in the upper
epidermis, so the axes differ in the variable that most controls the read-out.
Partly answered by the swab experiment (healthy skin swabs are 0/19 HPV+ while
DOCK8 swabs are 52/66), but it stays a stated limitation, not a solved problem.
No paired skin specimens exist in public data to fix it properly.

---

## 4. Framing: does it survive?

**Marko's structure changes the claim from an activity claim to a reservoir
claim.** That is defensible and arguably a bigger paper, but it inherits three
weaknesses the current manuscript does not have:

1. The DNA prevalence number is not yet solid (3.1).
2. The DNA/RNA comparison is cross-compartment and cross-pipeline (3.3).
3. DNA analysis 2, one of the six planned blocks, is infeasible (D2).

The current manuscript's framing is tight precisely because it stays in one
modality and controls depth explicitly. Adding DNA as the *spine* imports
confounders that the RNA-only design was built to avoid.

---

## 5. Two options, honestly

### Option A. Fold the DNA reservoir into the current manuscript
Add one figure (D1A + D1B as panels) and one Results subsection: *the reservoir
is real and independently confirmed*. It answers the sharpest objection to the
present paper, that the quiet majority might be a sampling artefact of RNA-seq.
The framing, title and argument are unchanged; the DNA becomes supporting
evidence rather than the spine.

- Cost: one figure, one subsection, plus 3.1 and 3.2 fixed.
- Risk: low. Nothing already written has to be rewritten.
- Marko's D1A, D1B and C1 all land; D2 becomes Table S7; R1, R2, R3 are already there.

### Option B. Build a new DNA-first skin paper
Restructure to Marko's order, with the reservoir established first and the
silence as the conclusion.

- Cost: new Introduction, new Results order, new title, plus 3.1, 3.2, 3.3.
- Risk: higher. The current manuscript is close to done; this restarts its
  argument. And it leaves the question of what the current manuscript then is.

**My recommendation is A**, with one addition: promote the co-infection table and
add the per-tier ORF panel, both of which Marko asks for and both of which are
cheap and already in hand. That delivers five of Marko's six blocks inside the
existing paper, drops only the one that is infeasible, and does not spend the
manuscript's completeness on a restructure.

If the group prefers B, the work is real but bounded, and this plan is the
figure list for it.

---

## 6. Open questions

None outstanding for Marko. The proteomics question was closed on 2026-08-07:
he meant proteomics, had already searched for datasets, and found none.

What remains is the group's decision between Option A and Option B in section 5,
and the three fixes in section 3, which are needed under either option.
