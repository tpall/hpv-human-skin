# Author comments and tracked changes: hpv_manuscript_2026-11-11b_jid_claude-edits.docx

Pulled 2026-11-25. The file carries two rounds in one: the Nov 18 copy-edit pass over draft
2026-11-11b (177 insertions, 146 deletions, 10 queries) and TP's own edits and 8 comments of
2026-11-25. TP asked for all tracked edits to be kept and reattributed to his name, so every
revision in that file now carries the author "Taavi Päll" and the file is the record of one
accepted round. Applied in `hpv_manuscript_2026-11-25a_jid.docx` (sources:
`_introduction_short.qmd`, `_results_short.qmd`, `_discussion_short.qmd`, `_methods_brief.qmd`,
`_methods.qmd`).

## TP comments

| # | Anchor | Comment | Action |
|---|---|---|---|
| 402 | "HPV reference panel and Kraken2 database" | Last download date could be required, these databases are not static; a comment in the next section says the same and was deleted | The drifting clause ("exact composition depends on build date") is replaced by the date itself: "Because both sources are updated continuously, the panel is fixed at build time: PaVE and RefSeq were queried on 17 April 2026, and every analysis reported here uses that build." Date taken from the panel on the cluster (`assets/hpv_references/hpv_all.fasta`, 2026-04-17), which every run read |
| 414 | "Of the 11,556 libraries that entered the sweep ..." | Duplicates table in the main text | Both accounting paragraphs deleted (the empty-report and droplet split, and the never-screened denominator paragraph). The main Methods already give the same three facts: 11,556 entered, 10,295 returned a result, 1,261 counted in no denominator because a library never screened is not a negative observation, with Figure S6 cited there |
| 419 | "HPV screening and positivity calling" | Rather main text methods material | TP's own deletion of the closing "central branch point" sentence removes the only restatement; the remaining text is the Kraken2 version, the read-identifier rule, `hpv_min_reads` and the seqtk extraction, none of which the main Methods carry |
| 422 | "HPV typing" | Main Methods too, currently duplicates some part of the main methods | Dropped the guard-rationale sentence ("The dual breadth-and-depth criterion guards against ..."), since the main Methods state the criterion. Kept STAR version and flags, samtools, the per-reference metrics and the `hpv_min_coverage` / `hpv_min_depth` names |
| 424 | "Early/late transcript classification" | Main methods material too | The L1-gating justification cut from 54 words to 23: "Productivity was gated on L1 rather than on combined L1+L2 signal, since isolated L2 signal can arise from low-level transcription or cross-mapping." featureCounts settings and the `late_transcript_min_reads` threshold kept |
| 430 | "Statistical analysis" | Main methods material | Dropped the one clause the main Methods carry verbatim (Wilson intervals). The rest is detail that appears nowhere else (floor percentile, the five matched read bins, the model specifications with the log10-depth covariate, rstanarm chains, $\hat{R}$, the Normal(0, 2.5) / Normal(0, 5) priors and the Student-t check) and is kept |
| 433 | "Integration status: E1:E2 transcript ratio and virus-host junctions" | Compress to describe what was done, and minimally why | 404 to 331 words. Every parameter, threshold and script kept; the explanatory asides go ("which is the integrated-carcinoma pattern of oncogene transcription without a late program", "A genuine integration yields one recurrent junction, whereas chimeric-read artifacts are scattered", "because the first is blind to the material the second detects", the single-cell truncation aside). The premise sentence is left exactly as `a2d8f5c` set it |
| 435 | "Independent DNA reservoir from public metagenomes" | See previous comment (compress) | 394 to 337 words, same rule. The two caveats that are constraints rather than commentary are kept: composition conditional on positivity, and no prevalence from the catalog. Preprint citations stay inline, per the JID preprint rule |

## TP tracked changes

1. **Introduction, paragraphs 2 to 4 rewritten** in a more formal register ("Understanding this is
   critical", "account for this role", "acts synergistically", "actively cleared", "transcriptional
   silence is the viral default", "Standard RNA sequencing (RNA-seq) pipelines"). Applied verbatim
   with three slips corrected: "theories account this role" to "account for this role", "whether
   these infection are productive" to "these infections", and the misplaced comma in "Although,
   viral reads are excluded".
2. **Main Methods: "Whole-library screening" paragraph deleted** (87 words).
3. **Supplementary Methods deletions applied:** the per-process environment and default-parameter
   sentences in Overview; the PaVE REST-API enumeration sentence; the coordinate-sorted-order
   rationale and the spot-count medians; the two exclusion-accounting paragraphs (comment 414); the
   branch-point sentence; and the whole "Execution at scale" section (chunked driver, aggregation,
   R/tidyverse sentence).
4. **The Nov 18 copy-edits**, kept under TP's name, applied across Results and Discussion: "four-fold"
   for "five-fold" cell-line enrichment (now rendered from `stats`, as in the Abstract),
   "distorted by" for "dominated by" contamination, "We found no difference between AK and cSCC",
   the L1-criterion wording in the trace paragraph, "the types detected spanned the expected
   cutaneous repertoire", "a coinfected library counts under each genus", "the 130 types resolved in
   either modality, with one exception", "absent at both the assigned and detection levels", the
   split of the petabase-survey sentence, and the mechanics in the main Methods (positivity
   threshold in the pipeline sentence, "Group assignment used", "the L1 criterion", "Integration was
   inferred").

## One addition

The Results and the Figure 5 caption both cite "(Methods)" for the DNA reservoir, which the main
Methods did not describe (Nov 18 query 310). The paragraph freed by deleting "Whole-library
screening" is now an **Independent DNA reservoir** paragraph of the same length, so the main
Methods net +1 word.

## Effect

| | Nov 11 | draft a |
|---|---|---|
| Abstract | 200 | 200 |
| Introduction | 720 | 756 |
| Results | 1,600 | 1,623 |
| Discussion | 1,173 | 1,197 |
| Materials and Methods | 477 | 481 |
| **Main text** | **3,970** | **4,057** |
| Supplementary Methods | 3,541 | 2,901 |
| Supplementary Material | 5,222 | 4,582 |

No number, interval or result moves. Counts are rendered body paragraphs with tables, captions,
table notes and headings excluded, which reproduces the 3,970 this log quoted on 2026-11-11.

## Still open

- **The Introduction rewrite costs 83 words** in a main text already over JID's suggested 3,500.
  The standing decision is to submit over the suggested limit; if the register is to be kept at the
  old length, the rewrite can be trimmed by about 85 words without touching its claims.
- **The Results heading still reads "dominated by nonprimary material"** while the Abstract says
  "concentrated in". Table 1 puts 61 of 294 positives (21%) in cell-line or engineered libraries, a
  four-fold enrichment in rate, so "concentrated in" is the defensible word in both places.
- **Six Nov 18 queries are unanswered** and are factual, not editorial: the 398-versus-406 cSCC/AK
  denominator, which denominator each of the four L1 percentages (91%, 90%, 93.2%, 64%) uses,
  whether the two HPV96/HPV182 lesions are two or three, where the 10,834 burden set is defined in
  the main text, the missing citation for "most often proposed as a cSCC candidate", and whether the
  AK any-HPV rate should be stated against the 28.7% non-permissive background. The deleted
  whole-library passage also removes the only statement that libraries were screened in full to a
  median 98.6% of their spots; the one-sentence fact remains in the supplement, the evidence does not.

---

## Round 2 (draft b, 2026-11-25): the three follow-ups and the six open queries

Draft b supersedes draft a.

### TP's three follow-ups

1. **Introduction trimmed back.** 757 to 691 rendered words, 30 below the Nov 11 length. Every claim
   and citation is unchanged and the register is TP's; what went was the padding in his own rewrite,
   mostly stacked discourse markers and nominalisations ("Ultimately", "Similarly", "Furthermore",
   "However, the literature is not unanimous;", "Specifically", "meaning that viral activity is
   expected to decline only when", "demonstrating that", "sequenced to a median depth of", "found no
   viral transcripts despite"). The three limitation sentences of paragraph 3 became one list.
2. **"concentrated in" in the Results heading**, so the heading, the Abstract and Table 1 now agree
   (61 of 294 positives, four-fold in rate).
3. **The whole-library fact restored as one clause**, in the supplementary Methods where the
   subsampling sentence already sits: "; the reads reaching Kraken2 were a median 98.6% of the
   SRA-reported spot count, and 95.8% to 99.5% across the individual sample groups." The
   coordinate-sorted-order rationale TP cut stays cut.

### The six Nov 18 queries

| Query | Finding | Action |
|---|---|---|
| 398 vs 406 cSCC/AK | Neither number was wrong: `.akc` splits the group into `neo_AK` (72), `neo_cSCC` (326) and `neo_mixed` (8), the last naming both diagnoses in its metadata | "Both the 72 AK and the 326 cSCC libraries, the remaining 8 of the 406 naming both diagnoses, carried ..." |
| Four L1 percentages | Three different populations and one that is not L1 at all. 93.2% is the burden floor in the 44 primary-tissue productive controls, not a capsid figure (the Discussion sentence quoting it as "productivity" was corrected in round 1). 90% and 95.5% are Beta-Binomial posterior means over the ORF panel, 40 of 44 controls and 10 of 10 EV, so they are not k/n and 95.5% cannot be read as a count. 91% is 49 of the 54 control libraries with gene-level counts, a wider set than the 44. 64% is 34 of the 53 HPV-positive typed targeted libraries | The prose now quotes counts and leaves the posterior means to the figure: "L1 in 40 of 44 productive controls and 10 of 10 EV libraries", "49 of the 54 productive controls with gene-level counts", "34 of 53 ... compared with 7 of 22". `results_models.R` exports `orf_ctrl_l1_k/n` and `orf_ev_l1_k/n` for this, with a comment stating that `pct` is a posterior mean |
| Two lesions or three | Consistent as written. `build_prod_table()` in `_helpers.R` lists ERR2274877 (carcinoma) and ERR1971044 (actinic keratosis), each typed "HPV96 (+HPV182)", so both lesions carry both types, and with the HPV12 cSCC that is three lesions | One word for safety: "2 lesions **each** coinfected with the beta types HPV96 and HPV182" |
| 10,834 undefined in the main text | True; it was defined only in the supplement | Defined where RPM is introduced in the main Methods: "Every library that returned a Kraken2 report, including negatives, was scored as RPM ...: 10,834 libraries, the sweep pooled with the targeted lesion cohorts and minus two excluded sets (amplicon and lesion-adjacent mucosa), of which 9,862 are also among the 10,295 screened." |
| "most often proposed as a cSCC candidate" uncited | The bibliography has no serological or DNA-prevalence entry that singles out HPV8, so the claim as written is not citable from it | Narrowed to what the bibliography supports: "HPV8, the beta type most often implicated in cutaneous carcinogenesis [@meyers2017; @viarisio2018; @morgan2024]". Two of the three are new to the short version, so the main-text reference count goes from 41 to 43 |
| AK any-HPV against background | 62.5% in AK against the 28.7% any-HPV rate in non-permissive tissue is a real excess, and the AK-versus-skin model gives any-HPV OR 2.6 (*P* = 1.5×10⁻⁴) | "... recovering an any-HPV read from 62.5% of AK libraries, twice the 28.7% assay background (any-HPV OR 2.6 against unselected skin)." The sentence now reads as presence above background, not as presence tout court |

### Effect

| | Nov 11 | draft a | draft b |
|---|---|---|---|
| Abstract | 200 | 200 | 200 |
| Introduction | 720 | 756 | 690 |
| Results | 1,600 | 1,623 | 1,650 |
| Discussion | 1,173 | 1,197 | 1,208 |
| Materials and Methods | 477 | 481 | 514 |
| **Main text** | **3,970** | **4,057** | **4,062** |
| Supplementary Material | 5,222 | 4,582 | 4,604 |

No result changes. `results_models.R` gained four exported counts and no computation.

The Introduction trim gave back 66 words, but the six query fixes cost more than that and the
main text is 92 words above its Nov 11 length: counts and their denominators in the Results (+27),
the burden-set definition in Methods (+33), the background clause in the Discussion (+11), and
three author-date citations for HPV8, which JID counts in the limit. If that is not wanted, the
cheapest reversals are the parenthetical counts in the typed-sample sentence and the Methods
definition, which the supplement also carries.

### Left open

- Nov 18 query 233 ("dominated by") is now closed by the heading change, and 310 (DNA reservoir
  absent from the main Methods) by the paragraph added in round 1. Queries 262 (part two) and 290
  remain: the biopsy claim in the Discussion says deep biopsies "miss surface and hair follicle
  layers", which is doubtful for punch biopsies that include epidermis and follicles. The point
  TP may want there is that swabs and plucked hairs sample different material, not that biopsies
  miss layers.

---

## Round 3 (draft c, 2026-11-25): the sampling claim

Draft c supersedes draft b, and is draft b plus one sentence.

Nov 18 query 290 is now closed, as TP agreed with the reading. The Discussion said RNA-seq uses
"deep biopsies that miss surface and hair follicle layers", which is not true of a punch biopsy that
includes epidermis and follicles. The claim is now the one the comparison actually rests on, that
the two assays sample different material: "skin epidemiology studies test surface swabs and plucked
hairs (92% and 73% prevalence) [@rollison2019sites], which sample different material from the
biopsies RNA-seq is run on." "Hair bulbs" became "plucked hairs", which is what the cited study
collected.

Main text 4,062 to 4,061 rendered words. Abstract, Methods and supplement unchanged. No result moves.

## Round d: external review of draft c (hpv_review_2026-11-25c.zip)

Applied on TP's instruction; tracked copy `hpv_manuscript_2026-11-25d_jid_tracked.docx`
(draft c as the base, 10 paragraphs, all revisions under TP's name).

- **Abstract.** "We reanalyzed 10,834" is now "over 10,000": the burden rates run on 10,834 and
  the four-fold odds and 2-3% on Table 1's 10,295, and the two sets are not nested. The reviewer's
  own two-count sentence would have put the Abstract at 207. "genuine" dropped to hold 200.
- **Introduction and Discussion, the 87-word set**, taken whole: the TP53/NOTCH example and the
  Son mechanism sentence leave the Introduction (both stand in the Discussion), the limitation list
  stops making one point twice, and seven wording cuts.
- **Introduction P3 repair.** The reviewer's merged sentence read "four of the five examined
  lesions without a comparison group"; it is now "Across these studies, ... four of the five had
  no comparison group of unselected skin, so ...".
- **Not taken:** the 33-word risk-label cut in Discussion P2 (TP's own point from round 11).

| | draft c | draft d |
|---|---:|---:|
| Abstract | 200 | 200 |
| Introduction | 690 | 624 |
| Discussion | 1,207 | 1,187 |
| **Main text** | **4,061** | **3,975** |

Still open: "four of the five" sits in a paragraph that cites six studies.

## Round d, continued: two reviewer notes (pre_supplement_fixes.md, supplement_cut_chemistry_integration.md)

Same draft d, same tracked copy (19 paragraphs). Also in this draft, on TP's instruction: the
Methods approvals sentence cut to "This work required no ethics committee review or participant
consent (see Ethics Statement)." (-17 words).

- **Depth claim (Methods).** "Sequencing depth entered every model" was false: M1 has no depth
  term and M3 is a Fisher test. Now "was modeled ... in the detection models".
- **Abstract**, held at 200: "controlling for cell-line origin, depth, and chemistry" is now
  "modeling detection against depth and chemistry"; the four-fold odds briefly carried "n = 10,295", withdrawn by TP as unnecessary once the opening says "over 10,000";
  "and did not differ" is now "; the two did not differ"; "instead", "only"
  and "therefore" stay.
- **Supplementary Methods, chemistry and integration subsections**, reviewer's text taken with two
  departures: the sentence no longer opens on a numeral ("Neither 3'-barcoded droplet libraries nor
  SMART template-switching ..."), and the 3, 6 and 6 are described as libraries with enough viral
  gene reads to be scored, which is what `bin/e1e2_integration_marker.R` filters on (viral >= 5),
  not E1 and E2 reads as the reviewer wrote.
- **The 15/15 "blind test" sentence is deleted.** The reviewer is right that it is circular:
  `RIBO_PAT` opens with the literal protocol string of that series.
- **R-hat.** The clause was not broken in the source; it was a math object that text extraction
  drops. Written as plain text now so it survives conversion.
- "can fire" and "criterion fired" replaced.
- **Not done:** the depth-adjusted M1 sensitivity refit on the 9,862 intersection (offered by the
  reviewer, needs TP's call), and the near-separation clause (reviewer's own lowest priority).

## Round d, continued: Discussion cut for length (TP's decision, 2026-11-25)

Discussion 1,187 to 955 words; main text 3,961 to 3,729. Tracked in the same copy.

- **The contamination paragraph is deleted** ("Unadjusted analyses are distorted by
  contamination..."). The HeLa point stands in the Introduction and the enrichment in the Results.
- P2: the sentence repeating that the data cannot separate the models is cut, and the
  "absence does not imply clearance" caveat is reduced to one clause leading into the risk-label
  sentence, which stays (it needed that antecedent).
- P3: the nine-valent vaccine sentence is cut.
- P4: the species sentence and the novelty sentence are merged into one clause; the species
  numbers (Betapapillomavirus 1, 5, 2) leave the Discussion.
- P5: the eyebrow-hair concordance clause is cut.
- P9: the single-cell closing sentence is cut.

Seven references leave the list with these cuts: Wenande 2025, Oh 2026, Devitt 2019, Selitsky
2020, Asplund 2019, Gihawi 2023, Cantalupo 2019.

## Round d, continued: contamination is a control, not a finding (TP, 2026-11-25)

TP: cell-line contamination is known in the field and not novel, so it is not argued in depth.
The first Results heading is now "Cohort and HPV positivity by sample origin", and its two
closing commentary sentences ("This is the signature of ..." and "Much of the apparent HPV
signal ...") are cut. Results 1,650 to 1,614; main text 3,729 to 3,693.

## Draft e (2026-11-26): referee report on draft d (referee_report_25d.md, review_issues_25d.csv)

Twenty issues, all in text that predates draft d. Each was checked against the code, the result
tables or the cited paper before it was acted on. Tracked copy:
`hpv_manuscript_2026-11-26e_jid_tracked.docx` (draft c as the base).

| ID | Verdict | Action |
|---|---|---|
| R4-1 | Confirmed: `ft` in `results_models.R` is productive controls against unselected skin (41/44 vs 7/7,800) | The Fisher OR now sits in its own sentence naming that contrast. The later "almost clean null" sentence, which repeated the 7 of 7,800, is cut to pay for it |
| R4-2 | Confirmed from the PubMed abstract: Hultin 2019 found HPV RNA in 1 of 12 lesions | Introduction says so; Discussion P1 says "little or no transcription" |
| R4-3 | Confirmed from the paper's own title ("High-Risk Human Papillomavirus E6/E7 mRNA") | "high-risk HPV" added in the Introduction; Fujimoto dropped from the Discussion P1 list |
| R4-4 | The count cannot be settled from abstracts (Arron's transcriptome design is unverified) | "four of the five" is now "most", which holds under either reading; the last clause is narrowed to "none scored productive infection from capsid transcripts" |
| R4-5 | Confirmed against `bin/hpv_rpm_contrast.R` and the per-sample table | Group list rewritten from the code: eight groups with live sizes, engineered second, non-permissive tissue assigned through the override table, the contamination rule stated as matching no library |
| R4-6 | Confirmed, and the cause found: all 317 were length-filtered (read 1 of 26 to 35 nt), none a failed download; 280 belong to two studies without a diagnosis, 18 carried a lesion diagnosis | New paragraph in the supplementary Methods with live counts; `bin/build_unscreened_table.py` extended to the lesion cohorts and the Supplementary Data table rebuilt (1,261 sweep rows unchanged, 317 lesion rows added, new `cohort` column); Data Availability updated |
| R4-7 | Confirmed: SRR11310307 passes the viral-read filter but has no ratio | Legend n is now 2, 6 and 6; the supplement gives both counts; Results state n = 2 beside 0.36 |
| R4-8 | Confirmed: no integration-like call is reported anywhere | The two calling sentences are deleted; the ratio stays a continuous score |
| R4-9 | Confirmed | "49 of the 54 wart, condyloma and EV libraries"; the 34-of-53 sentence says its population is defined by search |
| R4-10 | Confirmed | "HPV positivity" in the Discussion; "Clinical primary tissue was HPV-positive in only 2-3%" in the Results; "primary tissue" in the Abstract |
| R4-11 | Confirmed | Trace defined with both bounds; any-HPV glossed on the reads; "To interpret the any-HPV rate" |
| R4-12 to R4-14 | Taken | E2 sentence rewritten; Discussion P5 split after the third reason; Introduction P2 ends at "to carcinoma." |
| R4-15, R4-16 | Taken | Register replacements; HPV8 wording; "one actinic keratosis and one cSCC" (ERR1971044, ERR2274877 checked); "residual viral genomes" |
| R4-17 | Mostly taken | Table S5 "Site not recorded"; Figure S3 key "within group" (figure rebuilt, the only figure that changed); 36 nt; seqkit citation; "Fisher exact tests run". Not done: P = .7 beside P = .49 |
| R4-18 | Not done | Renumbering supplementary items this late is the larger risk, as the referee says |
| R4-19 | Part | The three phrases are cut. The accession list the Data Availability Statement promised did not exist; it is now built by `bin/build_accession_table.py` (12,860 libraries, with screened status and sample group; its counts reproduce Figure S6) |
| R4-20 | Open | CRediT line, with the pending sign-off |

Main text 3,693 to 3,688 (Introduction 617, Results 1,616, Discussion 955, Methods 500). Abstract 200.

Follow-up on R4-6 (TP asked whether the 317 are single-cell libraries): mostly not. By read
geometry 280 have two mates of 27 to 28 nt and come from two spatial-profiling studies annotated by
area of illumination (SRP413119, SRP413116); 23 pair a 26 to 28 nt barcode read with a cDNA read,
the droplet single-cell layout; 14 are 35 nt single-end libraries. The supplement gives the three
classes. Tissue categories in both Supplementary Data tables are now in English (skin, anogenital,
oral, other) in place of the pipeline's Estonian strings.

TP, on draft e: "a library with any papillomavirus read is any-HPV" read as circular. any-HPV is a
defined label, so the sentence now introduces it as one: "Results counted at the level of a single
papillomavirus read, whether or not the library reaches the 10-read positivity threshold, are
termed any-HPV." Results 1,616 to 1,622; main text 3,694.
