# Framing angles: how title/abstract/emphasis shift by target journal

The same dataset and results support three framings. The science is identical; what
changes is **what leads, what the "so what" is, and what gets foregrounded vs
backgrounded**. Pick one before finalising the abstract. The current draft is closest
to Angle A, though less exclusively so than when this note was first written.

*Revised 2026-08-23 against the current draft. The framing rewrite of 2026-08-22
(nine commits, 47 to 77 references) changed what two of these angles can honestly
claim. The stale parts are marked below where it matters.*

| Element | A. Tumor-virology | B. Dermatology / skin cancer | C. Data-reuse / methods |
|---|---|---|---|
| Target journals | J Virol, Tumour Virus Research, PLoS Pathogens | JID, Br J Dermatol, JAMA Dermatol | GigaScience, mBio, Genome Medicine |
| Core thesis | Beta-HPV is transcriptionally trace, not productive, in cSCC: the phenomenology every current account must accommodate | HPV is present but not productively replicating in cSCC and AK, so its clinical relevance is at the precursor and immune stage, not the tumour | Public RNA-seq can be re-mined for the cutaneous virome, but only with contamination, depth and chemistry controls |
| Hero display | Figure 2 (depth-normalised burden) | Figure 2, plus the AK/cSCC split and Table S2 exceptions | Table 1 + Figure 3 (programme), with @suppfig-depth promoted |
| Figure 3 (ORF programme) role | Supporting: the criterion made visible | Supporting | **Co-lead**: the assay's own read-out |
| Table 1 role | Supporting (contamination context) | Supporting | **Lead/centrepiece** |
| Backgrounded | pipeline mechanics | virology mechanism detail | the hit-and-run debate (becomes one application) |

**Figure numbers were reassigned on 2026-08-23** so filenames match display order:
1 integration, 2 depth-normalised burden, 3 transcriptional programme, 4 genus
prevalence. Earlier drafts of this note called the burden figure "Figure 1".

---

## Angle A: Tumor-virology (closest to the current draft)

**Thesis:** the productive/late-transcript phenotype settles where beta-HPV is and
isn't active across the skin disease spectrum, at a denominator and against controls
no prior series has had.

**Title (current):** *Cutaneous HPV is transcriptionally quiet in skin and carcinoma
but productive in warts: a depth-controlled re-mining of public RNA-seq*

**Abstract opens on:** HPV biology, the commensal beta/gamma reservoir, and the
unsettled role in cSCC.
**Abstract closes on:** "a trace, non-productive footprint that any account of
beta-HPV in skin must accommodate."
**Intro emphasis:** life cycle, productive vs latent, the beta-HPV/cSCC mechanism
literature, and the five prior transcriptomic and in-situ series.
**Discussion emphasis:** what the transcriptomic phenotype does and does not
adjudicate between the three live accounts.

**What changed, and what it costs this angle.** The draft no longer claims to test
the hit-and-run model, and the corresponding section heading was renamed. Son 2025
showed the protective-immunity account predicts the same declining viral signal, so
a decline discriminates nothing; Rehm 2022 supplies a third account in which
quiescence is the cis-encoded default. The old abstract close, "the phenomenology
both the hit-and-run and protective-immunity models must accommodate," is gone.

This is a real cost to Angle A: its "so what" used to be adjudication and can no
longer be. What replaces it is denominator, controls and criterion, which is a
methods-flavoured claim inside a virology framing. Reviewers at PLoS Pathogens in
particular are likely to want mechanism this angle cannot supply, so treat that
venue as a reach rather than a target.

---

## Angle B: Dermatology / skin cancer

**Thesis:** across the skin disease spectrum, HPV is productive in warts and EV and
present only as a non-productive trace in actinic keratosis and cSCC alike, which
locates the virus's clinical relevance at the precursor and immune stage rather than
in the established tumour.

**Title options:**
- *No productive HPV infection in cutaneous squamous cell carcinoma or actinic keratosis at population scale*
- *Human papillomavirus is transcriptionally silent in cutaneous squamous cell carcinoma: depth- and chemistry-controlled re-mining of 10,849 public transcriptomes*

**Abstract opens on:** the clinical problem, beta-HPV's contested role in
keratinocyte carcinoma and the immunosuppressed-patient risk.
**Abstract closes on:** where the virus can and cannot matter, given that it is not
productively replicating in the lesion.
**Intro emphasis:** cSCC epidemiology, organ-transplant risk, AK to cSCC
progression; compress the molecular life-cycle detail.
**Discussion emphasis:** the prior dermatology series (Arron, Vandiver, Fujimoto,
Purdie), the AK/HPV8 mechanism, and what VAXAK does and does not settle.
**Foreground:** the cSCC/AK result, the AK-versus-cSCC equivalence, and the
genuine-productive-lesion exceptions (Table S2). **Background:** pipeline internals
and contamination methodology.

**This angle got stronger, and its old thesis got falsified.**

Stronger, because nearly all the literature added on 2026-08-22 is dermatological:
Vandiver (*Exp Dermatol*), Fujimoto (*Am J Dermatopathol*), Morgan on actinic
keratosis (*JCI Insight*), Wenande/VAXAK and Oh 2026 (both *JAMA Dermatol*), the
Rollison cohorts. Both sides of the central literature sit in JID: Arron 2011 for the
negative, Purdie 2005 for a positive in-situ result. The conversation the paper now
joins is a dermatology conversation.

Falsified, because the previous thesis here was that HPV is "not a useful
therapeutic or screening target" in skin cancer. **Do not use that sentence.** The
VAXAK randomised trial reduced actinic-keratosis burden with nine-valent
vaccination, and the manuscript now cites it. The defensible version is narrower and
more interesting: an immunological intervention can reduce lesion burden without
productive viral replication being present in the lesion to target, which is
consistent with our finding rather than in tension with it. Note also that VAXAK's
vaccine is alpha-directed, so it cannot be neutralising the cutaneous beta
reservoir, and that it moved AK burden but not keratinocyte-carcinoma incidence over
12 months.

Two further additions serve this angle specifically: the AK arm and the cSCC arm are
statistically indistinguishable (trace OR 1.4, *P* = 0.28; floor detection Fisher
*P* = 1), which justifies pooling and answers the obvious objection; and the HPV8
actinic-keratosis mechanism is engaged directly rather than left for a referee.

---

## Angle C: Data-reuse / methods

**Thesis:** the petabase of public RNA-seq is a viable but booby-trapped substrate
for virome discovery; we provide a contamination-aware, depth- and
chemistry-normalised, reproducible framework, and HPV-in-skin is the demonstrating
case study.

**Title options:**
- *A contamination-aware, depth-normalised framework for mining the cutaneous virome from public RNA-seq*
- *Re-mining public transcriptomes for papillomaviruses: cell-line contamination, sequencing depth and library chemistry are the dominant confounders*

**Abstract opens on:** the opportunity and the trap. Non-human reads are discarded;
archives are huge but contaminated; naive mining overstates prevalence.
**Abstract closes on:** the generalisable framework and its reusability beyond HPV.
**Intro emphasis:** the discarded-reads opportunity, the prior archive-scale skin
minings (Leiendecker, Bhargava), the contamination literature, and the
methodological gap.
**Discussion emphasis:** the five-fold cell-line enrichment with Selitsky's traced
mechanism, the chemistry covariate, and the L1-anchored assembly-bias critique.
**Foreground:** Table 1, Figure 3, and Code availability; promote `@suppfig-depth`
into the main figures. **Background:** the hit-and-run debate as one application.

**What changed.** The "necessary minimum for credible virus mining" claim was
softened to an operationalisation claim, because Cantalupo and Pipas set out the
requirements in 2019 and Leiendecker and Bhargava had already mined skin archives.
That removes a priority claim this angle previously leaned on, and replaces it with a
quantified one: each safeguard is reported with the size of the artefact it guards
against.

Three things now strengthen this angle instead. Library chemistry is measured rather
than assumed, resolved for 62.2% of libraries from SRA metadata plus
library-construction protocol text, with the finding that poly-A versus rRNA
depletion does not move detection while amplification chemistries move it three- to
four-fold. Whole-library screening is stated and verified against SRA spot counts.
And the assembly-bias critique of L1-anchored contig catalogues appears to be
genuinely unclaimed, which is the most distinctive methodological contribution
available to this framing and is currently buried in a discussion subsection.

---

## Practical notes

- Switching angle = re-order the abstract, swap the lead Discussion paragraph, and
  move one or two display items between main text and supplement. The Methods,
  Results numbers, figures and `.bib` are unchanged.
- Both preconditions from the original note are now met: `full_v2` is complete, so
  the partial-denominator objection is gone, and the contamination result is stated
  as a contribution rather than a caveat.
- Quarto makes the swap cheap: edit `_abstract.qmd` ordering, retitle in
  `manuscript.qmd` YAML, and (for C) move the `@suppfig-depth` float from
  `_supplement.qmd` into `_results.qmd`.
- The reference list has grown from 47 to 77 entries and now covers all three
  audiences, so no angle requires new literature work before submission.
