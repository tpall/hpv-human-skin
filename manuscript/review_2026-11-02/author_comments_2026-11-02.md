# Author comments and tracked changes: hpv_manuscript_2026-11-02_TP.docx

Pulled 2026-11-02. All comments and tracked changes are by Taavi Päll and cover the Abstract,
Introduction and Results. Discussion and Methods were not yet reviewed. Applied in
`hpv_manuscript_2026-11-02_jid.docx` (sources: `_abstract_short.qmd`, `_introduction_short.qmd`,
`_results_short.qmd`, `_supplement.qmd`, `_methods_brief.qmd`, `_helpers.R`).

## Tracked changes

1. **Abstract:** delete "Because standard pipelines discard the non-human reads carrying viral
   signal," and change "re-mined" to "analyzed". Applied, as "reanalyzed" per comment 1 (JID closes
   the prefix: reanalysis, nonproductive).

## Comments

| # | Anchor | Comment | Action |
|---|---|---|---|
| 1, 2 | "Because standard pipelines discard ..." (Abstract) | Reader (a clinician) does not know what a standard pipeline is; use "re-analyzed" if apt | Clause deleted; "We reanalyzed ..." |
| 7 | "reached that floor" (Abstract) | What floor? | "a viral load threshold calibrated on [the wart controls]" |
| 8 | "a trace, non-productive footprint that any account ... must accommodate" | Rewrite plainly | "Together, these findings indicate that cutaneous HPV is productive in warts but present only at trace, transcriptionally inactive levels in healthy skin and carcinoma, which models of beta-HPV in skin cancer must account for." |
| 10 | "Dispensability is a claim about requirement rather than persistence ..." | Unwieldy start | Rewritten in 2 sentences and moved directly after the hit-and-run model it qualifies |
| 11 | "denominator" (Intro) | Readers may pause on what the denominator is | Word removed from Intro: "Most examined lesions alone, so they cannot say how often HPV is transcribed in carcinoma relative to skin without disease." |
| 12 | "Standard pipelines" (Intro) | Be clearer, even if more verbose | "RNA-seq studies of human tissue usually align reads to the human genome and discard the reads that do not match it. Viral reads are among those discarded, but they remain in the raw data ... and can be recovered by reanalysis." |
| 15 | "engineered libraries were positive in 9.28% and cell lines in 10.34%" | Explain engineered vs cell line at first use | New Results opening: the 3 searches (untargeted skin, targeted wart/condyloma/EV, broadened cSCC/AK) and the cell line / engineered / clinical classes in one sentence each |
| 17 | Table 1 caption | Short caption, rest as footnote, as in eluring | Table 1 and all 8 supplementary tables: bold one-sentence title above, note below in the 9 pt Table Note style (`docx_table_rules.py`) |
| 18 | "As an independent test" | Test of what? | "To test whether the ratio tracks integration in primary tissue, where integration increases with disease progression, ..." |
| 19 | "along the same axis" | Use "axis" consistently | Replaced by "across the same grades"; "axis" now only for plot axes |
| 20 | "The carcinoma median still sits well below the cell-line median." | Median of what? | "Even in carcinoma, however, the median E1:E2 ratio remained well below that of cell lines." |
| 21 | "The rise" | The rise of what? | "The rise of the ratio with grade in primary tissue ..." |
| 24 | "including negatives" | Implies negatives were not null | "including libraries below the HPV-positivity threshold of 10 reads" |
| 25 | "reached it" | Reached what? | "reached the floor"; floor defined in the preceding sentence (13.2 RPM, 5th percentile of primary-tissue productive controls) |
| 27 | "This gap" | What gap? | "Sequencing depth did not explain this difference between neoplasia and productive controls." |
| 28 | "the odds of any HPV read" | Finding any? | "the odds of finding any HPV read" |
| 29 | "either group" | Which groups? | "no neoplasia or unselected-skin library reached the productive floor" |
| 30 | "Pooling actinic keratosis with carcinoma is supported rather than assumed ..." | More direct | "AK and cSCC were indistinguishable, which justifies analyzing them together." |
| 31 | "Tissue that cannot host the virus fixes the meaning of that trace." | Rephrase | "To interpret this trace, we examined tissues that HPV cannot infect." |
| 32 | "puts a mechanism under that floor" | Mechanism must be assumed by the reader | "Splitting the below-floor libraries by viral gene shows what this trace consists of." No mechanism claimed |
| 33 | "All principal contrasts were refitted as Bayesian ..." | Out of place; final models should all use non-flat priors | Sentence removed from Results; Table S3 now cited from Methods. ⚠ See note below |
| 34 | "makes the distinction visible" | Distinction of what? | "The viral transcription program separates productive infection from this trace" |
| 37 | "the broadened query" | Methods come after Results | Defined in the new Results opening; here "the broadened cSCC/AK search" |
| 40 | "re-tabulated" | Check usage | "tabulated" (Results and Figure 5 legend); we tabulate the catalog by body site for the first time, nothing is re-done |
| 41 | "returns the mucosal niche" | Say alpha was detected in vaginal metagenomes | "Applied to vaginal metagenomes, the same procedure recovered mostly alpha types (84.5%), which are characteristic of the mucosal niche, ..." |
| 42 | "The 14 commonest reservoir types are beta and gamma but for HPV2 ..." | Hard to read | Split in 2: "Of the 14 most common types in the skin DNA reservoir, all except HPV2 are beta or gamma types. In our 7,800 unselected-skin RNA libraries, these types accounted for 0 assigned and 9 trace-level detections." |

⚠ **Comment 33 rests on an assumption that does not hold.** The principal contrasts (M1 to M3:
Table 1 odds ratios, the depth and chemistry models, the AK vs cSCC comparison) are frequentist
logistic regressions and Fisher tests, i.e. maximum likelihood with no prior. The Bayesian fits
with weakly informative priors are only the robustness check in Table S3. Removing the sentence is
fine, but making non-flat priors the default would mean refitting the main models and replacing
every OR/*P* in the text with posterior summaries. That is a decision, not an edit; not done.

## JID house style applied in the same sections

From the analysis of the 2026-11-02 JID issue (8 Original Articles):

- Numerals for numbers below 10 (5-fold, 4 times, 6 orders, 2 lesions, all 3 models, 4 of the 5).
- *P* values without a leading zero (`fp()` in `_helpers.R`, so everywhere: *P* = .28).
- "compared with" instead of "against"/"compared to" in prose; "vs" kept inside parentheses.
- Closed prefixes: nonprimary, nonproductive, nontumor, nonanogenital, posttransplant, coinfected,
  reanalyzed.
- HPV, RNA-seq, EV and AK defined at first use; the abstract now defines HPV.
- Register: "incredibly", "booby-trapped", "Critically", "Ultimately", "quiet" infections removed;
  fewer "claim: elaboration" colons.
- Abstract closes on "Together, these findings indicate ...", the dominant JID pattern.
- New citation: Hou et al. 2026 (WHIM syndrome, J Invest Dermatol 146:3142-54), Intro paragraph 1.

## Still open

- **Word count.** Main text (Intro + Results + Discussion, legends excluded) went from 3,878 to
  4,150 by the scratchpad counter, +272 words, so roughly 3,800 on the counter that gave 3,525.
  Most of it is the clarifications asked for in comments 12, 15, 18, 31. Cuts would have to come
  from the Discussion, which has not been reviewed yet. Abstract: 200 words.
- **Title:** changed to "transcriptionally inactive" (draft b). The Discussion still says "quiet".
- **Discussion and Methods:** not reviewed; one "compared to" and the house-style items remain there.
- **Citation format:** DONE in draft b (`hpv_manuscript_2026-11-02b_jid.docx`), via the CSL:
  in-text "(Hou et al, 2026; Doorbar, 2023)", same-author citations oldest first; bibliography
  "J Invest Dermatol 2026;146:3142-54." (no issue, closing period, "et al." kept). Journal
  abbreviations are NLM forms in a new `shortjournal` field on all 84 journal entries in
  `references.bib`, which the CSL reads via `container-title form="short"`.

## JID front matter and submission format (draft c)

From the JID Information for Authors, read 2026-11-02. Files: `hpv_manuscript_2026-11-02c_jid.docx`
(co-author copy, eluring style) and `hpv_manuscript_2026-11-02c_jid_submit.docx` (JID format).

- **Title** cut to 89 characters (limit 120): "Cutaneous HPV is transcriptionally inactive in skin and
  carcinoma but productive in warts". The internal subtitle is gone.
- **Title page** (`_titlepage.qmd`, then a page break): authors with numbered affiliations,
  place of work, corresponding author's address, telephone and email, short title "Cutaneous HPV
  transcription in skin and cSCC" (44 characters, limit 45), keywords not repeating the title
  (beta-papillomavirus; actinic keratosis; epidermodysplasia verruciformis; viral integration; cell
  line contamination), and the abbreviations list. Keywords removed from under the abstract.
- **Back matter** in JID order: Ethics, Data Availability, ORCiDs, Conflict of Interest,
  Acknowledgments, Author Contributions, AI declaration. Page break before References.
- **Preprints** may not be in the reference list. Shen et al. 2026 and Chikhi et al. 2024 (both
  still preprints on 2026-11-02; Shen is an eLife Reviewed Preprint, RP111778) are now a Word
  footnote in Results and "(…, bioRxiv preprint)" in the Figure 5 legend, Table S8 note and the
  supplementary Methods (`_methods.qmd`). Swap back to a citation once either is published.
- **SRA** defined at first use.
- **Submission format**: `python3 docx_jid_submit.py <draft>.docx <draft>_submit.docx` sets 12 pt,
  removes line numbers and adds page numbers. Co-author drafts keep the eluring 11 pt with line
  numbers.

Not done: the Ethics Statement does not name the ethics committees of the original studies or say
whether the original consent covered secondary use, which JID asks for in reanalyses of patient
data. Word count unchanged (about 3,800 on the old counter, against the suggested 3,500).

## Ethics Statement rebuilt for the JID secondary-analysis rule (draft d)

JID's Information for Authors (read 2026-11-02) requires, for "reanalysis, or secondary
analysis, of published patient data":

1. a list of all sources, with repository names, accession numbers and direct links, plus
   how the data were obtained and any preprocessing;
2. **the name of the ethics committee(s) that approved the original data collection**; and
3. a statement on whether informed consent was collected and whether that consent covered
   secondary use, or the exemption/waiver basis if not.

It also requires approvals and consent to be stated in the **Methods**, not only in the
Ethics Statement.

### Scope decision

Naming committees for every source study is not possible: 526 studies contribute, 426 of
them primary human material. Taavi chose (2026-11-02) to compile committee names for the
**36 studies behind the lesion and neoplasia claims** (559 libraries), and to state the
exemption basis for the remainder. The resolution of the 36:

| Status | Studies | Libraries |
|---|---|---|
| Committee named in the originating publication | 17 | 143 |
| Originating publication not identified | 12 | 284 |
| Committee named by reference to earlier trial papers only | 1 | 83 |
| Not patient material (cell line, xenograft, organotypic) | 3 | 29 |
| Publication names no committee | 2 | 13 |
| Committee not named specifically ("local Ethics Committees") | 1 | 7 |

So 284 of 559 libraries (51%) come from deposits with no recoverable committee name. The
statement says this in those terms rather than implying full coverage.

### How it was built

`manuscript/ethics_audit/` (four scripts, committed) resolves each SRA study to its
BioProject, then to a publication through NCBI elink and Europe PMC, then pulls the
approval and consent sentences verbatim from PMC full text. Output:
`manuscript/supplementary_data_ethics_provenance.tsv` (36 rows, 10 columns). The counts
quoted in the Ethics Statement are derived from that table by `eth()` in `_helpers.R`, so
prose and table cannot drift.

### Two errors the audit caught, worth keeping in mind

- **A reanalysis masquerading as the originating study.** Europe PMC matched ERP022968 (25
  AK/IEC/SCC libraries) to PMID 36225328, a keloid paper from Yonsei University. Its full
  text shows it *downloaded* ERP022968 and SRP227263 as control data. Citing its approval
  for our dataset would have attributed a keloid study's consent to skin carcinoma samples.
  Recorded as "originating publication not identified". Any accession-to-paper match needs
  this check.
- **A drug approval read as an ethics approval.** The pattern `approved by the` matched
  "approved by the United States Food and Drug Administration" in SRP448136. The verifier
  now excludes regulatory approvals.

### Data-quality finding, not an ethics matter

Three studies carrying a carcinoma or keratosis diagnosis are not patient material at all:
SRP520129 (HSC5 cSCC cell line plus mouse xenografts), SRP495320 (animal study) and
SRP528200 (organotypic 3D culture). Both publications with a human statement say
"Informed Consent Statement: Not applicable". **These look like `is_cell_line` /
`is_in_vitro` false negatives in `bin/cell_line_patterns.py` and should be checked against
the classifier**, since they carry 29 libraries into the lesion/neoplasia stratum.

### Still open

- **The UT waiver.** The statement says no approval was required and none was sought. An
  earlier draft asserted this as a matter of Estonian law, which we have not checked; the
  claim was removed and a VERIFY comment left in `_backmatter.qmd`. Some journals want the
  waiver in writing from the local committee.
- **ERP022968 and ERP106532** (49 libraries) have no identified originating publication.
  PRJEB24684 was submitted via EBI on 2018-07-10 by the Division of Genomics of Development
  and Disease, Institute for Molecular Bioscience, University of Queensland, which is a lead
  worth following; its description ("25 AK, IEC or SCC lesions from 17 individuals ... to
  detect and measure HPV E7 transcription") suggests a findable paper.

## Round 2 of TP comments, on draft c (13 comments + 6 insertions + 4 deletions) → draft e

Commented copy preserved as `hpv_manuscript_2026-11-02c_TP.docx`. Applied to the source and
rebuilt as draft e. Draft d was the ethics round and is left intact as the record.

| # | Anchor | Comment | Action |
|---|---|---|---|
| 0 | Author names on the title page | "this is pandoc extra, let's try to remove" | Removed `author:`/`affiliations:` from the YAML of `manuscript_short.qmd`, so the pandoc docx writer emits no byline. The names now appear once, from `_titlepage.qmd`, with the numbered affiliations JID wants. ORCiDs stay in the back matter. |
| 2 | "Cell-line and engineered libraries had 5-fold higher odds" (Abstract) | "this is framed here as our main finding" | Abstract results reordered. The warts/EV versus cSCC/AK contrast now leads; the cell-line enrichment follows as the confounder it is, introduced by "HPV signal instead concentrated in nonprimary material". No numbers changed. |
| 3 | "that threshold" (Abstract) | "productive infection threshold?" | Now "that productive threshold". |
| 10 | Table 1 note | "delete, restates results" | Deleted the sentence; the note keeps only the abbreviation expansions. |
| 16 | "with L1 represented by ≥ 3 reads" | below-floor and productive "collude mentally, needs to differentiate them more clearly" | Added two sentences to Results naming the two criteria explicitly: the detection floor is a burden threshold in RPM, productivity is a capsid criterion (L1 ≥ 3 reads), and a library can sit below the floor and still be scored for either. |
| 25 | Figure 4 | "strip text too small, difficult to read" | `figure_genus.R`: figure narrowed 10.5in → 8in and strip text 9 → 11 pt, so on-page size goes from about 5.6 pt to about 8.9 pt (Word scales the image to the 6.5in text width; on-page size is `size × 6.5 / width`, which the comment in the script now records). Axis text 11 pt, axis titles 13 pt, point labels 2.3 → 3. "Non-permissive tissue" broken over three lines, because strip width is the binding constraint; 12 pt clipped "Unselected" and "Engineered". |
| 27 | "public metagenomes (Methods)" | "there is also elife version now, still without official review, only public review" | Swapped to the eLife Reviewed Preprint. See the note below. |
| 28 | "in 0 of the 10,834 RNA" | "update inline R so that in appropriate cases 0 is replaced with none" | Added `n_or_none()` to `_helpers.R` (with a `cap` argument for sentence-initial use) and applied it to the HPV8 count, which now reads "in none of the 10,834 RNA libraries". |
| 31 | "silence in cutaneous carcinoma (Introduction)" | "if anything let's just add those references again" | "(Introduction)" replaced with the five citations themselves: @arron2011; @hultin2019; @arroyomuhr2021; @vandiver2021; @fujimoto2019. |
| 33 | The "broader objection" sentence | "too long sentence. use iso 24495 skill to split up such" | Split into three sentences and merged with the differentiation point during the trim. |
| 34 | "Our data measure this exact axis" | "rephrase without the use of axis pronoun" | Now "Distinguishing only transcription from its absence, our data cannot separate the three models." |
| 35 | "compared" (E2 integrity) | "why we want to compare them?" | Rewritten to state the purpose first: cervical carcinoma is where integration, and so E2 disruption, is expected, making it the natural test of whether a catalog can read genome state. It fails. Also fixes "compared to" → JID's "compared with" by removing the construction. |
| 37 | "per-group" (Methods) | "what groups?" | The groups are now named: productive controls, EV, cSCC and AK, unselected skin, and the cell-line and engineered comparators. |

Tracked changes, all accepted:

- Abstract: "in only 1.5%" → "in only 1.5% **of** libraries" (the insertion was the bare word "libraries"; "of" added for grammar).
- Abstract: deleted "Neither sequencing depth nor library chemistry explained this difference."
- Results: "reached the productive threshold" → "met the productive criterion of L1 represented by at least 3 reads".
- Figure 2 legend: deleted "Productive controls and EV biopsies lie far above it, whereas ...".
- Figure 3 legend: deleted "Productive groups transcribe the whole genome including the capsid genes, ...".
- Discussion: deleted "(Introduction)" (see comment 31).

### The eLife version (comment 27)

The two versions are the same paper. All seven figure titles match verbatim and every headline
number is identical (383 novel types, 105 host species, 65% re-identified = 646 of 992, 1,097
clusters, 21,319 L1 sequences, 34% expansion, 26.8M libraries, ~10-hour search). eLife's version
history: bioRxiv 22 April 2026, sent for review 1 May, Reviewed Preprint v1 on 16 July, not
revised. So v1 is the submitted manuscript with public reviews attached.

It is still a preprint, not a version of record, so JID's rule keeps it a **footnote**, not a
reference-list entry, either way. The footnote now cites eLife, with two deliberate choices:

- the **author order is corrected** to Shen, Korobeynikov, Chikhi, Babaian (the old footnote had
  Chikhi and Korobeynikov transposed);
- the DOI is given as `10.7554/eLife.111778` with "(Reviewed Preprint version 1, 16 July 2026;
  not yet a version of record)". eLife's own page says that DOI "represents all versions, and
  will always resolve to the latest one", so the version is stated in words rather than by
  inventing a versioned DOI we have not verified.

The Logan citation stays bioRxiv, which is still correct. The Zenodo record remains the data
source cited in the supplementary Methods, since that is what we actually tabulated.

### Discussion trim: target not reached

Cut 1,181 → 1,082 source words, and main text 4,145 → 4,073 rendered. That is 573 words over
JID's suggested 3,500, not the ~200 hoped for. Sentence-level and paragraph-level work yielded
only about 100 words because the Discussion was already dense: it carries no throat-clearing
openers, few nominalisations and almost no repeated sentences. Reaching 3,500 now means dropping
content, which is the author's call. See "Level 3 options" put to TP.

## Level 3 cuts and the classifier overrides (draft f)

TP chose, 2026-11-02: cut the petabase/E2 paragraph and the vaccine sentences; fix the
classifier by curation override only.

### Discussion, Level 3

| Cut | Words | Note |
|---|---|---|
| Petabase/E2 paragraph | 123 → 57 | Not deleted outright. Replaced by one sentence that keeps the claim and the `@supptbl-e2` crossref, so Table S8 stays cited and the argument survives in compressed form. |
| Nine-valent vaccine sentences | 45 → 0 | `@wenande2025` is now cited nowhere in the live manuscript (it survives only in the retired `_discussion.qmd`), so the VAXAK trial has **left the reference list**. Flagged to TP: it is a 2025 randomized trial on exactly our disease, and a reviewer may ask why it is absent. |

Discussion 1,181 → 980 source words. Main text 4,145 → **3,968** rendered, still 468 over the
suggested 3,500. Reaching 3,500 would now require cutting Results, which the condensing rules
put last for good reason.

### Classifier: curation overrides only

Four rows added to `assets/tier_overrides.csv` with a written reason each, following the
SRP560153 and SRP572392 precedents. No pattern edits, so no other library is reclassified.

| Study | Field | Why the classifier cannot reach it |
|---|---|---|
| SRP520129 | `is_cell_line` | `cell_line_patterns.py` has no HSC line name, and `ENGINEERED_PATTERNS` implements the lowercase-prefix convention for `sh` and `si` but not `sg`, so "HSC5 cells, sgWDR77#1" reads as naive tissue |
| SRP495320 | `exclude` | Every SRA field asserts human skin with actinic keratosis, organism *Homo sapiens*. Unreachable by any regex. Excluded rather than called a cell line, since it is neither patient tissue nor a human line |
| SRP528200 | `is_cell_line` | The existing `scc[\s-]?(4\|9\|13\|25)` pattern would match "RNA-Seq of 3D model of SCC13 cells", but that string is in the SRA `experiment_title`, which `parse_metadata.py` never reads; the per-run title is empty |
| SRP055713 | `is_cell_line` | Turku UT-SCC cSCC line panel. Nine other UT-SCC studies are already tiered `cell_line`; this one has no `cell_line` attribute and no siRNA label |

**These take effect only on a rerun**, and the expected direction is against us: the neoplasia
tier's any-HPV rate goes 54.7% → 57.0% and the cell-line tier's 26.4% → 25.9%. All 41 libraries
are HPV-negative or trace and none reaches the detection floor, so `pct_ge_floor` is unaffected.
State the direction explicitly when the numbers are refreshed.

### Deferred, not done

Two systemic defects found during the audit, both of which mean pattern edits alone would change
nothing downstream:

1. **`is_in_vitro` is computed but never consumed.** `bin/hpv_rpm_contrast.R` tiers on
   `is_cell_line` and `is_engineered` only, so 103 libraries flagged in-vitro are tiered as
   primary.
2. **`results_full_v2/metadata/cell_line_flags_curated.tsv` is stale** and is what the null
   cohort is tiered from: 158 libraries are `is_cell_line=true` in `samplesheet_enriched.csv`
   but `FALSE` there. Removing them from unselected skin moves its any-HPV rate 28.0% → 27.9%,
   so this is a defensibility fix rather than a results fix. Preserve the 10
   `source=study_curation` reversals (SRP545498, SRP514210) when regenerating.

Also deferred: feeding `experiment_title` to the classifier. It is the highest-value single
change (597 currently-unflagged libraries, including the 6 of SRP528200), but it is a pipeline
change, not a regex edit, and it reaches a 127-library study in the unselected-skin null, which
is the background every comparison is measured against.

**Do not add an `scc12` pattern.** ERP022968's `SCC_12` is patient 12 of a genuine clinical
AK/IEC/cSCC series (E-MTAB-5678), which is why the existing list skips 12. Note this is the same
ERP022968 whose originating publication the ethics audit could not identify.

## Vaccine clause restored (draft g)

TP: "use that single clause". The nine-valent vaccine finding is back as one clause
appended to the HPV8 sentence in the Son/reservoir paragraph, rather than the two
sentences cut at Level 3:

> The nine-valent vaccine, which reduced actinic keratosis burden in a randomized
> trial [@wenande2025], targets mucosal alpha types rather than this beta reservoir.

Kept: the trial finding, the design word ("randomized"), and the alpha/beta
mismatch, which is the reason the finding does not bear on our result.
Dropped from the original: the speculative second sentence attributing any
protective effect to cross-reactive or non-specific immunity. That was our
conjecture, not a reported result, and it is the part a reviewer would question.

@wenande2025 is back in the reference list. Main text 3,968 to 3,988 rendered words
(+20). Draft g supersedes f.

## Ethics Statement narrowed to open-access reuse (draft h)

TP, 2026-11-02, on JID's reanalysis rule: "this concerns reuse of patient data, we
are using public sequencing libs with public metadata, we dont need to cite original
works ethics permits."

Applied. The Ethics Statement no longer recites the originating studies' ethics
permits and no longer tallies committee coverage. It now says, in two paragraphs:
every library came from an open-access archive deposit with only archive-published
metadata; Supplementary Data lists all sources (repository, study and BioProject
accession, direct link, originating publication where identified, tissue and
diagnosis annotation) with retrieval and preprocessing in Methods; approval and
consent for the primary collection were obtained by the submitting investigators as
a condition of deposit; and the archive records do not state whether consent covered
secondary use, for which the deposit terms are the exemption basis.

Our own exemption basis (first paragraph) is unchanged and still carries the VERIFY
comment for the University of Tartu Research Ethics Committee.
[Superseded 2026-11-03: there is no such committee to ask. See Round 7.]

Also narrowed, for consistency: the Methods approvals paragraph no longer promises
"the per-study approvals of the primary studies", and the Data Availability
Statement describes the provenance table by what it contains rather than as a
committee list.

Kept deliberately: supplementary_data_ethics_provenance.tsv and ethics_audit/ stay in
the package, with the committee and consent columns intact, so the information can be
produced if an editor asks. The eth() helpers in _helpers.R are now unused by the
text; kept for the same reason.

Risk recorded, not hidden: the lesion libraries are human reads from AK and cSCC
biopsies with recorded diagnoses, so an editor may read them as published patient
data despite the deposits being open access, and the guide says manuscripts missing
required elements may be returned before review. TP's call, logged in the source
comment in _backmatter.qmd.

Counts: Ethics Statement 345 to 270 words, Methods 500 to 477. JID-countable
Introduction + Results + Discussion + Methods 3,870 to 3,847. Draft h supersedes g.

## Round 3 external review package, triaged and partly applied (draft i)

Source: ~/Downloads/hpv_review_r3_package.zip and ethics_statement_for_coauthors.md,
prepared elsewhere on draft g. 20 new points, 6 of them major. Verified against the
data and the rendered text rather than accepted as written.

### Ethics: GDPR basis adopted (TP's decision)

TP chose the package's proposal over keeping the earlier wording: the exemption now
rests on Regulation (EU) 2016/679 Article 2 with Recital 26, and "de-identified" is
replaced by "anonymous" in all three places (Ethics paragraphs 1 and 3, Methods).
Zero occurrences of "de-identified" remain in the manuscript.

Adapted rather than pasted, because the package predates TP's narrowing of the
Ethics Statement earlier the same day: its Methods replacement reinstated "the
per-study approvals of the primary studies", and its paragraph 3 ended on "no
per-study committee list was compiled". Both were dropped; the GDPR wording was
applied on top of draft h's narrowed text, not on top of draft g's.

Risk recorded in the source comment in _backmatter.qmd: "anonymous" is a stronger
characterisation than "de-identified" and raw human reads are widely treated as
personal data even stripped of identifiers, so this asserts a legal conclusion about
genomic data. Supporting facts now stated in the text: all deposits are open access
rather than controlled, and no individual-level information beyond archive annotation
was used. Confirming with the UT Research Ethics Committee is still recommended.
[Superseded 2026-11-03: there is no such committee to ask. See Round 7.]

### Majors verified and fixed

- **R3-6, three genus lists for the same libraries.** Confirmed from type_genus.csv:
  HPV12, HPV96 and HPV182 are beta, HPV16 is alpha, none is gamma. Table S2's legend
  said "beta- and gamma-HPV" and was simply wrong; fixed to beta. The Discussion
  silently excluded the two productive HPV16 cSCCs from its beta-only count; it now
  reads "mucosal HPV16 excepted (@supptbl-cin)" and "productive infection by a
  cutaneous type".
- **R3-4, Figure 3c sentence.** Confirmed: it asserted a clean partition that the
  paper's own 13% contradicts two paragraphs earlier. Rewritten to describe the
  overlap and to cite the 13% rather than contradict it.
- **R3-2, "sample-class definitions are identical across models".** Confirmed in
  _methods.qmd. Replaced: cell-line status is defined identically, but the two
  analyses are fitted on different populations with different class definitions and
  their counts are not interchangeable. Also fixed a dead crossref (@tbl-rates, no
  such label; the table is @tbl-positivity) and Figure 2's legend, which called the
  burden cohort "screened" where the paper uses that word for the sweep.

### R3-1: the package's diagnosis is wrong, and the real problem is larger

The package frames this as a filtering chain, 11,556 to 11,282 to 10,834 with an
unaccounted 274-library step. It is not a chain. Measured directly:

- results_full_v2/metadata/hpv_status.tsv (Table 1, positivity model): 11,556
- results_reframe_fixed/hpv_rpm_per_sample.tsv (burden cohort): 10,834
- in Table 1 but not in the burden cohort: 1,694
- in the burden cohort but not in Table 1: 972
- common to both: 9,862

Two overlapping populations, because Table 1 comes from the full_v2 sweep alone while
the burden cohort pools four result directories. Ruled out as explanations for the
1,694: the 433 amplicon libraries, the 15 adjacent-mucosa libraries, and single-cell
exclusion (46 of the 1,694 are single-cell, against 94 retained). They are 1,009
`nahk` and 685 `muu` with no library_selection recorded.

NOT FIXED, and not papered over: the reason 1,694 sweep libraries have no RPM record
needs to come from bin/hpv_rpm_contrast.R before submission. Writing a Methods
sentence that reconciles the two numbers without knowing why would be invention. The
Methods now states that the two counts are not interchangeable, which is true and is
the most that can be said honestly today.

### Still open from the package

- R3-3, "trace" carries three meanings about 250-fold apart. Real, not yet addressed;
  it reaches the title and the abstract, so it needs TP's decision on terminology.
- R3-7, the override rerun (29 lesion-stratum libraries, including the SRP495320
  animal study). Already known; needs the pipeline rerun.
- Publication DOIs or URLs for the 23 bare PMIDs in the provenance table.
- The data-minimisation sentence (age and sex deliberately not extracted).
- The sensitivity analysis on the 17 committee-named studies.

Draft i supersedes h. The review zip in manuscript/ is built on g and is now stale.

## Ethics material deduplicated across three sections (draft j)

TP: "now we have Approvals, consent and data provenance and Ethics Statement in
methods with largely overlapping content". Correct, and it was three-way: both of
those also re-described the Supplementary provenance table that the Data Availability
Statement already describes column by column. Each section now has one job.

- **Methods, retitled "Approvals and consent"** (was "Approvals, consent and data
  provenance"): two sentences. The approvals-and-consent fact that JID requires in
  Methods, pointing to the Ethics Statement for the basis, plus the retrieval and
  preprocessing sentence that JID's reanalysis rule asks for and that appears nowhere
  else in the main text. Dropped from here: the GDPR scope clause (Ethics
  Statement's job) and the Supplementary accession-list description (Data
  Availability's job).
- **Ethics Statement, paragraph 2:** no longer enumerates the provenance table's
  columns, which the Data Availability Statement lists in full. It now points there
  and keeps only the ethics content: the deposits are open access with archive
  metadata only, approval and consent for primary collection were the submitting
  investigators' and are documented in the originating studies, and where reported
  they are reproduced verbatim in the table.
- **Data Availability Statement:** unchanged. It is the one place the table is
  described in full, which is where a reader looking for the data will go.

Counts: Methods 477 to 430 words; Ethics Statement 343 to 303, the trim offsetting
most of the GDPR sentence added in draft i. Results and Discussion each grew slightly
in draft i (the R3-4 rewrite and the mucosal-HPV16 clause), so the JID-countable
Introduction + Results + Discussion + Methods is 3,843, against 3,847 in draft h.

No claim, accession, basis or disclosure was dropped; only restatements of material
that survives elsewhere in the same manuscript. Draft j supersedes i.

## JID heading styles, page breaks, and the ethics sections cut (draft k)

TP: "fix the heading styles and page breaks and approval and ethics sections are way
too verbose, cut them down, use iso-24495 on them too".

### docx_jid_submit.py

Two additions, both from "How to Format Your Manuscript" verbatim.

**Heading styles.** JID: "First level: ALL CAPS, BOLD; Second level: Sentence Case,
Bold Underline; Third level: Sentence Case, Bold Italics." All four pandoc heading
styles were already bold and Heading4 already bold italic, so the script adds
<w:caps/> to Heading1 and Heading2 and a single underline to Heading3. Caps are a run
property, not uppercased text, so the heading text in the file is unchanged and the
look can be undone without retyping. Level mapping for this manuscript is recorded in
the script's docstring: Heading2 carries the main sections, Heading3 the Results and
supplementary subsections, Heading1 only the "Supplementary Material" divider, so
Heading1 and Heading2 are both first-level in JID's sense.

**Page breaks.** JID: "Begin a new page for References, Tables, and Figure Legends."
The script now breaks before any heading named References, Tables, Figure Legends or
Supplementary Material; two fire in this manuscript (References, Supplementary
Material), on top of the two already in the document.

NOT done, because it is a decision about the document rather than its formatting:
tables and figure legends sit inline here, where Quarto puts them, not in sections of
their own. JID's own order of sections puts them at the end, after References. Moving
them is a restructuring choice for TP. JID also says formatting requirements "will be
imposed at the time of first revision", so nothing here blocks submission.

### Ethics Statement and Approvals paragraph, under iso-24495

Reader taken as the editor or reviewer checking compliance, so each JID requirement
has to be findable at a glance rather than embedded in a clause.

Ethics Statement 303 to 214 words (-29%); Methods paragraph 430 to 424 for the
section, the paragraph itself 70 to 62. Techniques that did the work: passive to
active with an explicit actor ("The authors had access to no identifiable private
information" to "We used no identifiable private information"), nominalisations to
verbs, one name per concept, and deletion of restatement ("this is the exemption basis
for our use of the data" after the sentence that is the basis; "Because the analysis
used only such data"; "including when processed for research purposes" to "even when
it is processed for research").

Every requirement survives, checked one by one against the rule: source list with
repository, accession and direct link (paragraph 2, pointing at Data Availability);
how obtained and preprocessed (paragraph 2, pointing at Methods); committee position
(paragraph 2); consent and secondary use with the exemption basis (paragraphs 1 and
3); anonymity, open-access tier, no controlled access, no re-identification attempt,
no individual-level data, GDPR Article 2 and Recital 26, no human subjects, none
sought (paragraph 1).

JID-countable Introduction + Results + Discussion + Methods: 3,837. Draft k
supersedes j.

## Round 3 of TP comments, on draft k (6 comments, no tracked changes) -> draft l

Commented copy preserved as hpv_manuscript_2026-11-02k_TP.docx.

| # | Anchor | Comment | What was done |
|---|---|---|---|
| 16 | Ethics Statement, para 2 | "can be removed, as all this is also in the next chapter. it's only function is to say that keep reading" | Paragraph deleted. Two clauses kept, as JID requirements appearing nowhere else: who approved the original collection, and whether consent covered secondary use. Both folded into para 1. |
| 17 | Ethics Statement, para 3 | "what's the function of this sentence, i think this can be remove as it's redundant with previous" | Paragraph deleted. The deposit-terms argument and "This basis applies to every source study" dropped as genuinely redundant. |
| 97 | Figure S3 | "axis text is too small to read" | Figure transposed to portrait. See below. |
| 101 | Table S1 | "use smaller font size in table cells to fit text, currently text is wrapped. also in other tables. footnote is way too extensive" | All 9 data tables set to 8 pt in docx_table_rules.py; note trimmed. |
| 106 | Table S5 note | "here and in other places remove results from the footnote. try to present elsewhere if not redundant." | Done for all six notes carrying results. |
| 109 | Table S7 note | "entire results section in the footnote!" | 484 words to 92. |

### Ethics Statement: now one paragraph, 303 to 136 words

Both deletions were right. The surviving paragraph still carries every JID clause:
anonymity, open-access tier, no controlled access, no re-identification attempt, no
individual-level data, GDPR Article 2 and Recital 26, no human subjects, none sought,
who approved the original collection, and the secondary-use position.

### Figure S3: the font was not the problem

Transposing it was the only fix available. Word scales an image to the ~6.5 in text
width, so 66 type columns get 0.0985 in each, about 7 pt. A label rotated 90 degrees
occupies roughly 1.2 x its font size in width, so nothing above about 5.9 pt can fit
whatever size the figure is drawn at, and widening the drawing makes it worse because
Word then shrinks it further. The 66 categories had to move to the axis with room for
them, which on a portrait page is the vertical one: 66 rows over about 8.5 in is
0.13 in per row. Types now run down, tiers across, 7 x 9.6 in. Checked by eye at print
scale. The reservoir claim still reads as one long run, now down the unselected-skin
column of the detection panel rather than across its row. The arithmetic is recorded in
the script header so this is not re-litigated at proof.

### Table notes: 1,881 to 714 words across the nine

Six notes carried results. Each note now describes its table (what the rows and columns
are, the estimator, the thresholds, the cross-references) and the results moved into the
supplementary narrative as six new subsections, where they are prose rather than legend:

| Table | Note words | New subsection |
|---|---|---|
| S1 types | 226 to 83 | What the assigned-type catalog shows |
| S4 mucosal HPV16 | 179 to 46 | Why the mucosal HPV16 cases are retained |
| S5 chemistry | 218 to 41 | Library chemistry does not manufacture the productivity contrast |
| S6 trace | 303 to 157 | What the detection level adds |
| S7 co-infection | 484 to 92 | How far the co-infection table can be read |
| S8 E2 | 227 to 51 | Why the E2 statistic cannot read genome state |

Nothing was deleted in the move except two sentences that restated their own table's
structure. The three untouched notes (37, 88, 119 words) are already descriptive.

Table cells are set to 8 pt as a run property, so the 9 pt table notes and the 12 pt
body are unaffected and docx_jid_submit.py's 11-to-12 pt bump, which works on
styles.xml, leaves them alone.

JID-countable Introduction + Results + Discussion + Methods: 3,837, unchanged, since
every cut was in material JID excludes. Draft l supersedes k.

## Positivity denominator corrected, draft l -> draft m (2026-11-03)

TP asked what caused the gap between Table 1's clinical count and the 10,834 of
the burden analysis. The answer turned out to be a defect in Table 1, not a
bookkeeping difference, so it is recorded here in full.

**What the two numbers are.** Not two ends of a filtering chain. Only 9,229
libraries are in both: 1,684 of Table 1's clinical set carry no burden row, and
1,605 burden rows are not clinical-in-the-sweep (mostly the targeted lesion
cohorts, which never entered the sweep). The near-match of 10,913 and 10,834 is
a coincidence.

**The defect.** 1,261 of the 11,556 sweep libraries produced an empty Kraken2
report and therefore no screening result, yet `hpv_status.tsv` records every one
of them as HPV-negative with a read count of 0. Table 1 was counting libraries
that were never screened as negative observations, and 1,251 of the 1,261 fall
in the clinical row, the reference group. Two causes:

- 410 were downloaded and lost every read to quality filtering, all of them
  `too_short_reads`: `fastp_min_length = 36` meets an R1 of 26 to 33 nt (median
  28), the barcode and UMI mate of droplet single-cell chemistry, and fastp in
  paired mode discards the pair, informative R2 included.
- 851 were never downloaded, spread 10 to 22 per 100-library batch across 85 of
  122 batches, i.e. per-accession download failures, not failed batches.

The depth-normalized burden analysis never carried the error: it reads the
report files themselves, and `read_kraken_rpm` returns NULL for an unparseable
one. The arithmetic closes exactly: 12,860 report files, minus 1,578 empty
(1,261 + 4 + 313 across the three aggregated dirs), gives 11,282; minus the 448
exclusions (433 amplicon, 15 adjacent mucosa) gives 10,834.

**What changed.** `unscreened_ids()` in `_paths.R` is the single source of truth;
`results_models.R` and `bayesian_models.R` both filter on it, so the frequentist
and Bayesian positivity models share a denominator.

| | draft l | draft m |
|---|---|---|
| Clinical | 10,913, 2.14% (1.88-2.42) | 9,662, **2.41% (2.12-2.74)** |
| Engineered | 194, 9.28% | 190, 9.47% |
| Cell line | 449, 9.58% | 443, 9.71% |
| OR engineered / cell line | 4.69 / 4.85 | 4.24 / 4.35 |
| Bayesian OR engineered / cell line | 4.57 / 5.25 | 4.09 / 4.30 |

Nothing reverses and no conclusion changes direction. Three prose claims had to
follow the numbers, because they no longer matched them:

- "roughly 5-fold higher odds" -> "roughly 4-fold" (Results, and the Abstract).
- Abstract "*P* < 10 to the -8" -> "*P* < 10 to the -7": the engineered *P* is
  now 1.8e-8, so the old bound had become false.
- "2-4% of primary skin libraries" -> "2-3%", in the Abstract, Results and
  Discussion. The clinical rate is 2.41% (CI 2.12-2.74) and skin specifically
  2.58%; nothing in the data reaches 4%, and it did not under the old
  denominator either (2.14%, CI to 2.42), so this corrects a pre-existing
  overstatement rather than one the fix introduced.

The attrition is now stated in both Methods: the brief one gives the submitted,
screened and unscreened counts, the supplementary one gives both causes and the
reasoning. `supplementary_data_unscreened_libraries.tsv`
(`bin/build_unscreened_table.py`) lists all 1,261 with their cause, batch, read
lengths and reads submitted, so a reviewer can reconstruct the denominator from
the accession lists; it is cited in the Data Availability Statement.

JID-countable Intro+Results+Discussion+Methods: 3,837 -> 3,889 (+52, all in the
brief Methods). Abstract 195 words, unchanged.

**Not done, and deliberately:** the 410 short-read libraries are not re-run with
a chemistry-aware minimum length, and the 851 downloads are not retried. Both
need a pipeline run and would change the denominator again; TP's call.

### HPV8 reservoir sentence clarified, draft m -> draft n (2026-11-03)

TP: the sentence was still unclear. It was, and the fault was a mismatched pair
of denominators that read as a rate comparison without being one. The DNA side
(63 of 1,735) is conditional on both HPV positivity and skin; the RNA side
(10,834) was every library of every sample group, and the paragraph had already
established a skin RNA denominator of 7,800 three sentences earlier, then
switched without saying so.

Now: "in 63 of the 1,735 HPV-positive skin DNA libraries (3.6%), yet absent at
any level, assigned or trace, from all 7,800 of our unselected-skin RNA
libraries and from all 10,834 we analyzed." The DNA proportion is stated so the
conditioning is visible, and both RNA denominators are given, because zero holds
on either: the site-matched comparison and the whole analyzed cohort. +2 words.

JID-countable total 3,889 -> 3,891.

## Round 4 of TP comments, on drafts l and n -> draft o (2026-11-03)

Three comments, applied in the order TP asked for.

**Comment 83 (draft n), on the new Methods attrition passage: "need to be
rephrased."** Correct, and it was my sentence. "That policy leaves a visible
trace, and it must be read before any rate is quoted" was vague about what the
trace is and told the reader what to do instead of what happened. Now the
paragraph simply states the fact: "1,261 of the 11,556 libraries that entered
the sweep returned an empty Kraken2 report, and so no screening result."

**Comment 82 (draft n): "maybe a box chart can visualize better how we came up
with different numbers."** Added as @suppfig-cohort-flow (Figure S6,
`manuscript/figure_cohort_flow.R`). Drawing it surfaced something the prose had
not made explicit: the positivity denominator and the burden cohort share a
stage. The sweep's 10,295 screened libraries are both the Table 1 denominator
and, minus the 433 amplicon libraries, the null half of the burden cohort; the
lesion cohorts contribute 987 screened, minus 15 adjacent mucosa, giving
9,862 + 972 = 10,834. So the two populations are not nested and neither is a
filtered version of the other, which is exactly the confusion the original
question came from. Every count in the figure is computed from the same files
the analysis scripts read, never typed in, so it cannot drift from the text.

**Comment 1 (draft l), on the Abstract: "can we also include disjoint between
DNA and transcriptionally active types? + which types are present, and whether
the infection is productive."** The Introduction promises three things and the
Abstract delivered two; the disjoint, a whole Results paragraph and a main
figure, was absent. Added: "The transcribed types were largely disjoint from the
skin DNA reservoir: of its 115 types, 102 were productive in no RNA library."

The Abstract was at 195 of 200 words, so this was a trade, not an addition.
Freed by tightening, with no claim, number or hedge dropped:

| Before | After | Saved |
|---|---|---|
| "typing HPV and scoring ... , controlling for cell-line origin, sequencing depth, and library chemistry" | "typing HPV, scoring ... , and controlling for cell-line origin, depth, and chemistry" | 3 |
| "exceeded that productive threshold in only" | "exceeded it in only" | 3 |
| "and did not differ from each other" | "and did not differ" | 4 |
| "A minority of carcinoma and keratosis specimens" | "A minority of lesions" | 3 |
| "had 4-fold higher odds of HPV positivity ... mostly from integrated" | "had 4-fold higher odds of positivity ... mostly integrated" | 3 |
| "Together, these findings indicate that cutaneous HPV is productive" | "Cutaneous HPV is therefore productive" | 4 |

Abstract now exactly 200 words. "Which types are present" is still carried by
the existing "beta- or alpha-HPV" and "integrated HPV16 and HPV18" clauses
rather than by a type list; there is no room for a list, and the main text and
@fig-reservoir carry it.

JID-countable total unchanged at 3,891: the Abstract is excluded from the count,
and the Methods rephrase was word-neutral.

## Round 5 of TP comments, on draft o (4 comments + 1 tracked change) -> draft p (2026-11-03)

**Tracked change, Figure S3 caption: "upper panel" -> "left", "lower panel" ->
"right".** Correct, and it exposed more than it fixed. The figure was transposed
to portrait for comment 97 but the caption still described the old orientation
throughout: it said "Columns are HPV types ... rows are sample groups" when the
figure now has types as rows, sample groups as columns, and the two catalogs
side by side. Four further sentences read the wrong way ("concentrated in two
rows", "the unselected-skin row is pale", "runs almost the full width", "the
length of that row"). All corrected against the rendered PNG.

**Comment 95 ("use more concise language") and comment 100 ("there is a simpler
and straighter way to say this, also other places in supp fig captions, use the
iso-24495 skill").** All six supplementary figure captions rewritten with the
skill, reader taken as a referee skimming the supplement. No number, condition
or hedge dropped.

| Caption | Before | After |
|---|---|---|
| suppfig-depth | 151 | 116 |
| suppfig-trace | 185 | 156 |
| suppfig-types | 357 | 312 |
| suppfig-reservoir-gut | 303 | 282 |
| suppfig-bayes | 185 | 164 |
| suppfig-cohort-flow | 146 | 140 |
| **total** | **1,327** | **1,170** |

One substantive fix fell out of the rewrite: in the reservoir caption "Its
composition is a mixture of mucosal alpha and cutaneous beta" followed "Only
skin is carried into the main figure", so "Its" pointed at skin when the
sentence is about the gut site. Now named outright.

**Comment 104: "good!"** No action. (Anchored on the cohort-flow figure.)

**Comment 106, on HPV18 in Table S1: "check column width, hpv type must fit into
one line."** Pandoc derives column widths from the width of the source markdown
columns, so Table S1's first column took 51% of the table and the type columns
were left with as little as 2%, which wrapped the type names. `widen_columns()`
in `docx_table_rules.py` now gives every column at least the width of its
longest unbreakable word, hands the slack back in proportion to the original
widths, and scales the whole set down if the minimums cannot fit. All 9 tables
rebalanced; every column now clears its longest word with margin (Table S1's
narrowest is 0.53 in for "Total", needing 0.43 in).

Two bugs of mine surfaced while doing it, both now fixed in
`docx_table_rules.py`:

- `<w:t[^>]*>` also matches `<w:tc>`, `<w:tcPr/>` and `<w:tbl>`, so the
  longest-word measurement was reading markup as text and every column came out
  needing 2,400 to 3,200 twips. The pattern now requires a space or an immediate
  close after the element name. The word counts reported earlier are unaffected:
  body paragraphs carry no other `w:t*` element, and the counts recompute
  identically.
- Pandoc wraps each table and its caption in a one-cell outer table, so a
  non-greedy `<w:tbl>...</w:tbl>` match runs from the wrapper to the inner
  table's close. `shrink_cells` was therefore setting the table captions to 8 pt
  as well as the cells, in every draft since l. Both passes now split the
  wrapper off with `inner_table()` and touch only the data table, so captions
  are back at 12 pt.

JID-countable total unchanged at 3,891 (supplementary captions are excluded).
Abstract 200 words.

### Draft p edit picked up, and the review pack built on draft q (2026-11-03)

While packaging p for review, two things turned up in the file itself.

First, TP had already made a tracked deletion in p at 14:59, cutting the last
sentence of the Figure S5 caption: "The axis spans five orders of magnitude
because that contrast genuinely does." Applied to `_supplement.qmd`. The pack is
therefore built on draft q, which is p plus that one cut and nothing else. The
commented copy of p is kept as `hpv_manuscript_2026-11-03p_TP.docx`.

Second, p's submission copy had been resaved by Word, which renames the
script-written page-number footer part from `footer_jid.xml` to `footer1.xml`
and reassigns its relationship id. Harmless in itself, but it means the file was
no longer the script's output. Worth knowing for future packs: build the archive
from a fresh render rather than from a docx that has been open in Word. The q
copies are clean (checked: `footer_jid.xml` present, no comments, no tracked
changes, line numbers in the editing copy and none in the submission copy).

`hpv_manuscript_2026-11-03q_review.zip`, 4.0 MB, 17 files: README, both docx,
this log, 11 figure PDFs, 2 Supplementary Data tables. The README gives the
denominator correction with its before/after table, the length position, and
six open points to comment on.

---

## Round 6: supplement cross-reference audit and legend copy-edit (draft r, 2026-11-03)

Source: `supplement_copyedit_q.md`, an external audit of draft q. Applied in full except
where noted. Verified every checkable claim against the source and the rendered docx before
acting; one claim was a false positive and was not acted on.

### Cross-references: two orphaned supplementary figures

Figure S5 (`@suppfig-bayes`) was cited in no file at all, and Figure S6
(`@suppfig-cohort-flow`) only in the supplementary Methods, so neither was reachable from
the main text. Both now have a main-text citation at the anchor each already had:

- Methods, "Cohort.": `... 10,295 returned a result (Figure S6).` The reader who meets the
  three cohort numbers for the first time is who that diagram is for.
- Methods, "Confounder control and statistics.": `... as a robustness check (Table S3, Figure S5).`

### Four errors corrected

1. **Figure S3's worked example was arithmetically wrong.** It read "three libraries is 23%
   of the EV group". EV is N=10, and under the caption's own estimator (Beta-Binomial
   posterior mean, Beta(0.5, 0.5) prior) three of ten is 31.8%. 23% is what *two* of ten
   gives, which is what Table S1's EV column prints in four rows (HPV14, HPV25, HPV47,
   HPV9). One count was carrying a percentage computed from a different count, in the
   sentence whose only job is to teach the reader how to read the fill scale. Now 32%.
2. **The Data Availability Statement miscounted its own tables**, introducing the third item
   named as "A second Supplementary Data table". Recast to "Two further Supplementary Data
   tables accompany them. The first ... The second ...".
3. **Two sentences in the supplementary Methods began with a numeral** (1,261 and 410),
   which JID style and the house rules both forbid. Recast as "Of the 11,556 libraries that
   entered the sweep, 1,261 returned ..." and "Of these, 410 were downloaded but ...".
4. **A third sentence-initial numeral, in the main text**, which the audit missed: Methods,
   "Cohort." opened "11,556 libraries entered the screen". Recast as "Of the 11,556
   libraries that entered the screen, 10,295 returned a result", which is also where the
   Figure S6 citation went.

### Not acted on: one false positive

The audit reported "space before a period, four instances in the main text, three in the
supplement", read as inline-R artifacts. All seven are `P = .49`, `P = .0055` and the like,
which is JID's own convention of dropping the leading zero from a p-value. The check caught
space-plus-period and misread it. Nothing to fix.

### Legends brought under JID's 125-word limit

Five of six supplementary legends were over. All six now comply:

| Legend | Before | After |
|---|---:|---:|
| Figure S1 | 106 | 106 |
| Figure S2 | 141 | 82 |
| Figure S3 | 302 | 125 |
| Figure S4 | 260 | 116 |
| Figure S5 | 145 | 116 |
| Figure S6 | 129 | 120 |
| **Total** | **1,083** | **665** |

418 words left the legends. 263 of them moved into prose rather than being cut, because
they were the only place their content appears:

- **Figure S2's late-read fractions** (57 gene-mapped non-productive libraries, 52 early-only,
  5 L2-only, median late-read fraction 0 against 0.367) into a new subsection, "What the
  gene-level split shows" (50 words). Put in the supplement rather than in the Results
  paragraph that cites Figure S2, since the supplement is outside JID's word limit and the
  main text is already over it.
- **Figure S4's whole interpretive half** (the vaginal positive control, the gut
  collection-route reading, the novel-cluster fractions) into a new subsection, "What the
  vaginal and gut sites establish" (144 words). Figure S4 was the one item with no prose
  section to defer to, so its caption had been carrying the entire argument.
- **Figure S3's commensal-reservoir reading** into the existing "What the detection level
  adds" (69 words), which already had the adjacent fact (40 of 51 resolved skin calls are
  cutaneous beta/gamma/nu) but not the interpretation.

The remaining 155 words were duplication:

- The Table S1 / Table S6 / Figure S3 relationship had been stated four times. The precise
  version (Table S6's note, with the 325-library count) stays; Figure S3's caption keeps the
  partition clause because a reader of the figure alone needs it; Table S1's note and Table
  S6's caption are now bare pointers.
- The estimator was defined three times. Defined once in Table S1's note; Table S6's note
  now says "the estimator of Table S1"; Figure S3's caption keeps the eight-word
  parenthetical, since the fill scale is unreadable without it.
- Figure S3's caption restated three things already in "What the assigned-type catalog
  shows" and in Figure 4's caption. Cut.
- Figure S5's Student-t(3) sentence was internally tangled: it called the item "the
  Student-t(3) row" while explaining that it is not a row, because "row" meant a row of the
  figure in its first use and a row of the table in its second. Replaced with the one thing
  a reader of the figure needs, which is why the figure has a line the table does not.
- Figure S6 lost "and neither is a filtered version of the other", which "not nested"
  already carries.

### Effect

Main text 3,891 to 3,897 words (+6, all in Methods, from the two citations and the
"Cohort." recast). Abstract unchanged at 200. Supplement 4,631 to 4,474 rendered words.
Legends and supplementary material are both outside JID's word limit, so the +6 is the only
change that touches the count.

### Still open from the audit

Two mechanical items in its §7 were out of scope for this round: inconsistent thousands
separators in running text (needs an eye, not a regex, since accessions are correctly left
bare), and seven supplementary sentences over 50 words, four of them in "Read download and
quality control".

### Packaged for review

`hpv_manuscript_2026-11-03r_review.zip`, 4.0 MB, 17 files: README, both docx,
this log, 11 figure PDFs, 2 Supplementary Data tables. The figure PDFs are
unchanged from the draft q pack, since no figure code changed this round, only
the captions. The README gives the two orphaned cross-references, the four
errors, the legend table before and after, and the open points.

---

## Round 7: regressions from round 6, a diction pass, and the ethics position (draft s, 2026-11-03)

Sources: `draft_r_fixes.md`, `diction_pass_r.md`, `ethics_readme_fix.md`, three external
reviews of draft r. Every checkable claim was verified against the source and the rendered
drafts before acting. One was taken as stated by TP, one was declined.

### Three defects introduced by round 6, all now fixed

These were mine, created an hour earlier by the legend cuts.

1. **A double colon in Figure S3's caption**: "as one grid: they partition one cohort: a
   library either clears ...". The conjunction was lost when three words were trimmed to
   reach 125. Restored ("one cohort, so a library"), paid for by "32% of the EV group,
   0.04% of unselected skin".
2. **The moved Figure S4 passage opened "This figure"**, which had an antecedent inside
   Figure S4's caption and none in a subsection that never names the figure.
3. **The moved Figure S2 passage opened "Of the 57 gene-mapped but non-productive
   libraries"**, a definite article on a set the subsection never introduces, and it named
   neither Figure S2 nor the 230-library denominator.

Both passages resolved only by sitting next to the caption they came from. Open point 6 is
whether to move figure legends to the end per JID's section order; if that happens, both
referents break. Each now names its figure, which costs three words and survives reordering.
"What the detection level adds" already did this and was the model.

### Diction pass

- **"licenses"** ("This figure licenses the skin composition ... to be read as biology").
  Replication licensing (origin licensing, MCM2-7, Cdt1, Cdc6) is standard vocabulary in
  this field, and this is a papillomavirus paper that reads genome replication state from
  E1:E2 with E1 the viral helicase. Using "license" for "warrants" attaches a second meaning
  to a word the field already owns. Now "@suppfig-reservoir-gut is the positive control for
  the composition read-out: ... so the skin composition in @fig-reservoir can be read as
  biology", which also fixes defect 2 above.
- **Two metaphorical uses of "axis"**, "the presence axis" (supplementary Methods) and "the
  DNA axis" ("Why the E2 statistic cannot read genome state"). Both breach the decision
  recorded at comment 19, "axis now only for plot axes", restated at comment 34. Now
  "presence data" and "DNA evidence". The four genuine plot-axis uses stay.
- **Five inanimate-agent constructions**: "This table exists because" (deleted, the
  substance starts at "Poly-A capture"), "Table S8 is reported because the inference is an
  inviting one to draw" (recast), "SRP565659 makes the point under internal control",
  "the capsid burst the controls are there to show", "which is the criterion working".
- Left alone: "manufacture" in "Library chemistry does not manufacture the productivity
  contrast", which is the precise verb for an artefactually created effect, and "buys" in
  the Figure S5 caption.

### Ethics: the position is settled, and the committee question is closed

TP, 2026-11-03: the University of Tartu research ethics committee no longer exists, and the
Estonian committee now sits under ETAG. Because this work reanalyzes published, anonymized
sequencing libraries, it needs no ethics approval and that committee has nothing to decide
on. The recommendation to "confirm the position with the University of Tartu Research Ethics
Committee" is therefore not pending, it is unavailable, and the comment block in
`_backmatter.qmd` now records that. The two earlier log lines carrying the old
recommendation are marked superseded rather than rewritten, since they are dated records.
The Ethics Statement itself needed no change: it names no committee and rests on GDPR
Article 2 with Recital 26.

### Only SRA is named as the source

TP asked to name only the SRA, and the code agrees. `modules/local/sra_download/main.nf`
fetches every library with `prefetch` followed by `fasterq-dump`; there is no ENA or EBI
endpoint anywhere in the pipeline. The cohorts do contain accessions issued by the other
INSDC members (1,194 ERR and 55 DRR in the prevalence sweep, 102 ERR and 19 DRR in the
lesion cohorts), but those runs are mirrored into the SRA and that is where we downloaded
them from. Naming the SRA alone is therefore both accurate and more complete than before,
since the old wording named the ENA but not the DDBJ. Changed in the Ethics Statement and in
the matching Methods sentence, which states the same fact and would otherwise contradict it.

### Declined

`ethics_readme_fix.md` proposed restoring the committee-coverage enumeration (36 source
studies, 17 naming a committee, 284 of 559 libraries with no recoverable committee) to the
Ethics Statement, reporting that the one-paragraph cut had removed it by accident. Verified:
the numbers are in draft g and absent from draft r. But the removal was deliberate and the
reasoning is in `_backmatter.qmd`, dated 2026-11-02: requirement (ii) is addressed to reuse
of patient data, so the statement "no longer tallies committee coverage ... without the
manuscript asserting coverage we do not have". The data still ship in
`supplementary_data_ethics_provenance.tsv` and the Data Availability Statement still points
at it. Not reversed, and now closed: TP confirmed on 2026-11-03 that the current wording
stands, namely "Approval and consent for the original collection were the submitting
investigators' responsibility; no archive record states whether that consent covered
secondary use, which anonymous data do not require." This is a decision, not an open point,
and should not be reopened by a later audit of the same section.

TP added, 2026-11-03: if this position does not fly with JID, the paper is published
elsewhere. The statement is not weakened and no approval is manufactured to satisfy a
reviewer. See the journal shortlist for the alternatives.

Two further notes on that report: its "51%" appears in neither draft (284 of 559 is 50.8%,
so the figure is right but was never in the text), and `cut_plan_r3.md`, cited as the
do-not-cut authority, is not in the repository. For the record, every number in the proposed
sentence does check out against supplementary_data_ethics_provenance.tsv (36 studies, 559
libraries, 17 with a named committee, 12 unrecoverable contributing 284 libraries), but the
proposed one-sentence form accounted for only 29 of the 36 studies, leaving 7 for a reader
to wonder about, which draft g's longer enumeration did not.
