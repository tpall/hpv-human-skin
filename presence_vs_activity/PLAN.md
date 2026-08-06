# Paper 2: presence is not activity

Planning document. Started 2026-08-05. Companion to the cSCC manuscript in
`../manuscript/`, which this paper does not repeat.

## 1. The claim

Presence and productivity are decoupled, and the decoupling is **predictable
from viral genus crossed with anatomical niche** rather than random. Skin and
cervix are the two poles of one axis: in skin the DNA reservoir tells you almost
nothing about transcriptional activity, in cervix it tells you almost everything.

The general form of the claim, which is what makes it a paper rather than a
survey: **petabase DNA catalogues measure presence, and presence is not
activity.** Logan, Serratus, Shen and the Babaian papillomavirus survey all
report counts of the form "virus X in N libraries". Every one of those is a
presence number. This paper is the first systematic RNA-activity layer over such
a catalogue, and it shows the presence-to-activity inference failing in a
specific and predictable way.

Framing decision (2026-08-05): lead methodological, with the biology as the
mechanism that explains the failure. Not a body map with methods attached.

## 2. Journal targets

Follows from the framing. The methodological lead widens the audience beyond
tumour virology, which the cSCC paper is aimed at.

- **PLoS Pathogens**: primary target. Hit-and-run work lives there, the
  catalogue critique is in scope, and it takes negative and reframing results.
- **Nature Communications**: reach, justified by the cross-catalogue scale and
  the field-wide methodological claim.
- **Genome Biology**: strong fit if the reusable activity-layer pipeline is
  foregrounded as a resource.
- **mBio / Journal of Virology**: safe specialty fallbacks.

Deliberately not the same tier list as paper 1 (`project_journal_targets`),
because this is a positive, larger-scale result rather than a depth-controlled
negative.

## 3. Results structure

**R1. The DNA reservoir across body sites recovers textbook tropism.**
Status: DONE (`../dna_rna_triangulation/`). Skin is beta/gamma commensal,
vaginal is 84% alpha mucosal, gut is a mixture. Recovered from DNA alone, which
validates the composition read-out before anything is asked of it.

**R2. Skin: presence does not predict productivity.**
Status: SURVIVES, with a host axis added. The objection was that the DNA side is
100% surface swabs while the RNA side is bulk tissue, so the two axes differ in
the variable that most directly controls late transcription (5.5b). It is
answered primarily by the scale argument in **5.2b**, which needs no swab data:
HPV presence is 9 to 38% in shotgun skin DNA against 2 to 4% in skin RNA, and
transcription amplifies, so an active infection would make RNA the more
sensitive detector rather than the less. The swab experiment (5.5d) supports
this on a small corner of the data: healthy skin swabs carry no HPV transcripts
at all, at greater depth than the immunocompromised arm that is 79% positive.

The composition finding stands: the DNA-abundant beta/gamma types show no
productive signal in our biopsy RNA, the productive types are alpha wart viruses
that are DNA-rare or absent, and only HPV2 occupies the abundant-and-productive
corner. It remains a composition comparison across cohorts, since no paired skin
specimens exist (5.5c), so the wording must stay at that strength.

**R3. Cervix and vagina: presence does predict productivity.**
Status: NEEDS THE RUN (section 6), but now has a paired core. 22 cervical
specimens carry HPV-positive DNA and unscreened matched RNA from the same
BioSample (section 5.5), which turns part of R3 into a genuine paired statistic
rather than a composition comparison. Run those first.

The full-cohort version is the within-tissue join that R2 could only do for
skin: are the DNA-frequent vaginal types (HPV91/51/66/67/16/90) the
productive ones? Prediction is yes, which makes the cervix the clean contrast
case to skin. Carries a second axis: high-risk alpha dominating productive
expression by lesion grade while low-risk and commensal mucosal types stay quiet.

**R4. Anorectal: both genera present in DNA, only one productive.**
Status: PILOT (681 typed, 230 HPV+ at 34%, productive signal alpha and HPV16-led
at 91 libraries). This removes the between-site confound: the cutaneous beta
types HPV24 and HPV124 that appear in anal and gut DNA are entirely silent in
the anal transcriptome. The strongest figure in the paper.

Honest limit: this is within-*niche*, not within-*sample*. The DNA is gut
metagenome and the RNA is anorectal RNA-seq, different libraries. Section 5.5
settles this: no paired anal specimens exist in SRA, so the limit is permanent
and should be stated in the paper rather than treated as future work.

**R5. Where the decoupling has a molecular correlate.**
Status: PARTIAL. The graded cervical series (PRJNA612785) shows E1:E2 rising
monotonically with grade with late transcription falling on the same axis. The
full cervical run extends this from 24 libraries to a real distribution. Care
needed: the marker itself belongs to paper 1, so R5 must use it rather than
re-derive it, and must cite paper 1 for its validation.

## 4. What is already in hand

| asset | location | state |
|---|---|---|
| DNA reservoir by body site | `../dna_rna_triangulation/` | done |
| Skin RNA activity axis, 10,849 libraries curated | `../results_reframe_fixed/` | done (paper 1) |
| Anal RNA pilot, 681 typed | `../dna_rna_triangulation/anorectal_pilot/` | untiered |
| Cervical graded series, 24 libraries | `../results_cervix_graded/` | done |
| E1:E2 marker and junction validation | `../bin/e1e2_integration_marker.R`, `../bin/score_junction_clustering.py` | done (paper 1) |
| Contamination machinery | `../bin/cell_line_patterns.py`, `../bin/_overrides.R` | done, reusable |
| Cervical/vaginal candidate accessions, 6,396 | `../rna_extension_scoping/candidate_new_vaginal_cervical.srr.txt` | scoped, not run |
| Anorectal candidate accessions, 651 | `../rna_extension_scoping/candidate_new_anorectal.srr.txt` | scoped, not run |

## 5. Fix before new compute

All cheap, all are the obvious referee objections, and none needs the cluster.

### 5.1 Restrict the skin RNA side to unselected skin
Our skin libraries skew clinical and lesional (warts, EV, cSCC), so "beta is
silent" is partly a statement about which skin got sequenced. Recompute R2 with
the RNA side restricted to the unselected-skin tier. The direction is expected to
hold; the magnitude of the silent column will shrink. Until this is done R2 is
weaker than it reads. (Caveat 1 of the triangulation report.)

### 5.2b The denominator: computed, and definition-dependent (2026-08-06)

**The freeze is Dec 2023.** Shen's hits run to 2023 (2,226 released in 2023, one
in 2024), which matches the Logan v1 freeze of 27 M accessions rather than the
v1.2 Dec 2025 freeze. Denominators must therefore be restricted to runs public
by 2023-12-31, or they are inflated by two years of subsequent deposition.

**The numerator** is 2,737 HPV-positive skin-metagenome libraries in Shen after
our non-human exclusions.

**The denominator depends on what counts as a shotgun library, and the choice
swings the answer 14-fold.** All restricted to <= 2023:

| denominator | n | HPV+ rate |
|---|---:|---:|
| A. `WGS` strategy only | 7,304 | **37.5%** |
| B. DNA source, any strategy except AMPLICON | 29,169 | **9.4%** |
| C. DNA source, any strategy (incl. amplicon) | 104,536 | 2.6% |

C is wrong: 75,866 of those runs are 16S amplicon and cannot contain HPV at any
infection rate. The real answer is between A and B, and the gap between them is
19,525 runs with `LibraryStrategy = OTHER`. That label is used both for genuine
shotgun metagenomics and for 16S submitted without the AMPLICON tag, and SRA
metadata cannot separate the two. Resolving it would need read-length or
library-name heuristics on 19.5k runs.

**So quote a range, not a point.** HPV DNA positivity in skin shotgun
metagenomes is **9 to 38%**, against **2 to 4% in unselected skin RNA**. The
collaborators' recollection of ~25% sits inside that range.

**One caveat that must travel with the comparison.** The DNA rate comes from
Shen's detection (profile-HMM on assembled contigs) and the RNA rate from ours
(Kraken2 screen on reads). These are different pipelines with different
sensitivities, so the ratio is not a calibrated fold-difference and should not
be reported as one. What survives is the direction and the order of magnitude.

**The argument the gap carries.** Transcription amplifies: one genome makes many
transcripts. If the cutaneous reservoir were transcriptionally active at any
normal rate, RNA-seq would be the *more* sensitive detector of its presence, not
the less. Seeing presence in 9 to 38% of DNA libraries and 2 to 4% of RNA
libraries points at low transcription directly. For compartment sampling to
explain it instead, biopsy dilution would have to overcome the transcriptional
amplification factor and then invert the ordering, across thousands of libraries
from hundreds of independent studies. That is not plausible at this scale, and
it is the primary argument for the decoupling. It needs no swab data.

### 5.2 Get a denominator
`files/all_hits_info.list` is RunInfo for HPV+ libraries only, so every DNA
number is a numerator among HPV+ libraries: "HPV24 is frequent" means frequent
among HPV-positive skin metagenomes. An SRA count of all skin-metagenome runs
converts composition into prevalence. Caveat to state explicitly: not every SRA
run is in Logan, so the denominator is the assembled subset and needs its own
coverage check. (Caveat 2.)

### 5.3 Fix the subtype-letter join
`HPV38B` against `HPV38`, `HPV14D` against `HPV14`. These fail to join and become
false zeros or land in `Unresolved`. Small effect, does not touch the
alpha-versus-beta/gamma contrast, but exactly the kind of thing a referee finds.
(Caveat 3.)

### 5.4 Tier the anorectal pilot
The 681 libraries carry automatic HPV-positive and productive calls without the
per-sample curation the skin cohort received, and 82% carry a generic tissue
label with an anal-carcinoma diagnosis rather than an explicit anal tissue field.
R4 cannot go in the paper on automatic calls. (Caveat 5.)

### 5.5 Paired DNA and RNA from the same specimen: DONE, and it exists

Run 2026-08-05 with `scope_paired_dna_rna.py`. Results in `pairs_*.tsv`.

Three tiers are reported separately and only the first is a paired test:
`specimen` (same BioSample carries both modalities), `donor` (same project plus
a shared subject token, different BioSamples), `study` (same project only).

| niche | runs | BioSamples | specimen | donor | study |
|---|---:|---:|---:|---:|---:|
| anal | 1,682 | 1,464 | **0** | 0 | 0 |
| anal, project-expanded | 2,541 | 2,320 | **0** | 0 | 1 |
| cervical | 54,932 | 44,723 | **64** | 58 | 62 |
| cervical, project-expanded | 64,393 | 48,163 | **64** | 59 | 94 |
| Shen HPV+ BioSamples | 10,657 | 6,614 | **61** | 53 | 21 |

**Anal is a clean negative.** Anal DNA studies and anal RNA studies are entirely
disjoint BioProjects, verified independently of the tool. R4 therefore cannot be
upgraded to a paired test and its within-niche limitation stands as written.

**Cervix is the win.** 64 same-BioSample pairs, of which **22 have a BioSample
already in Shen's HPV-positive catalogue** (29 of the 138 DNA runs ARE Shen HPV+
libraries), and **none of the 108 paired RNA runs has been screened by this
project**. Those 22 are the paired test in its strongest form: HPV is already
known to be present in the DNA, so typing the matched RNA asks "present, but is
it active?" on the same material instead of across cohorts. At ~108 runs this is
trivial compute and should run BEFORE the 6,396-library cohort, since it is a
publishable result on its own and it de-risks the big run.

Only 2 of the 64 are cell lines (one SiHa, one CaSki) against a 54% HeLa base
rate in the broad cervical space. Requiring a matched DNA library selects hard
for clinical studies, because nobody submits exome data for a cell line as a
patient series. The pairs concentrate in seven projects, so this is a
read-by-hand list rather than 64 separate curation problems.

**Paired skin specimens exist too**, which was not expected: 3 in PRJNA872116
(skin metagenome, WGS + RNA-seq, HPV-positive DNA). Small, but skin is the pole
where presence is NOT supposed to predict activity, so even three paired
specimens test R2 directly rather than by composition. Worth typing alongside
the cervical set since the marginal cost is nil.

The Shen-anchored mode also surfaced **13 paired vaginal-metagenome specimens**
(PRJNA797778, AMPLICON+WGS DNA with RNA-seq), which is a mucosal paired set in
metagenome rather than tissue material, and so a useful protocol contrast.

**Project expansion changes the weak tiers only.** Probing every single-modality
project for the modality the site vocabulary could not see added 10,123 cervical
runs and took study-level projects from 62 to 94, but specimen pairs stayed at
exactly 64. That is expected in hindsight: a BioSample carrying both modalities
shares its metadata across both runs, so the site query already sees both. Run
expansion for the study tier, not for pairs.

### 5.5b RESULT: the paired run, and the compartment problem it exposed

Job 523244, 2026-08-06, 29 libraries typed (`paired_results/`). The cervical arm
behaved as predicted. The skin arm did not, and the reason matters more than the
result.

**Cervix: prediction held.** 19 of 24 cervical libraries are productive, most
with substantial late transcription. Presence predicted activity in the mucosal
pole, which is R3 in miniature and on the same material rather than across
cohorts.

**Skin: prediction failed.** Of 3 paired skin specimens, 1 was HPV-negative and
2 were productive, not silent:

| library | types (breadth) | E1 | E2 | L1 | L2 | productive |
|---|---|---:|---:|---:|---:|---|
| SRR21981202 | HPV190 (99.5%), HPV105 (99.2%), HPV195 (99.2%), HPV47, HPV5, HPV20, HPV124 | 16,200 | 10,975 | 5,601 | 6,264 | yes |
| SRR21981198 | HPV107 (44%), HPV-mm090c10nr, HPV4 | 419 | 123 | 598 | 355 | yes |
| SRR21981204 | none | - | - | - | - | HPV-negative |

Beta and gamma types, at up to 654x depth and near-complete genome coverage,
with abundant capsid transcription. This is exactly what R2 says does not happen.

**Two confounders, and they are not separable at n=2.** PRJNA872116 is "The
microbiome of RAG-deficient patients": RAG1/RAG2 hypomorphic immunodeficiency.
HPV5, HPV20 and HPV47 are the classic EV beta types, and EV is itself a genetic
immunodeficiency, so productive beta-HPV is the expected phenotype in this host.
The libraries are also `env_medium = swab` metatranscriptomes, i.e. skin
SURFACE, whereas our RNA cohort is bulk tissue.

**The second confounder is systematic, and it is the real problem.** Checking
the cached BioSample attributes for Shen's skin reservoir
(`../dna_rna_triangulation/cache/skin_biosample_attrs.tsv`, n = 941):

| sampling | n | % |
|---|---:|---:|
| swab-like | 577 | 61.3 |
| biopsy or tissue | **0** | **0.0** |
| unstated | 364 | 38.7 |

Every named collection device is a swab, gauze or pore strip. The env_local_scale
values are the standard microbiome skin sites (volar forearm, retroauricular
crease, antecubital fossa). **The entire DNA reservoir is surface sampling, and
the entire RNA axis is bulk tissue biopsy.**

Papillomavirus transcription is differentiation-coupled: the late programme
completes in the upper, differentiated layers, which is what a swab collects and
what a biopsy dilutes with dermis and basal keratinocytes. So the two axes of
the central comparison differ in the one variable that most directly controls
the read-out. The decoupling in R2 may be a compartment artefact rather than
biology, and this paired result is the first direct evidence that it might be:
when the same compartment is sampled both ways, the beta/gamma reservoir is
loudly productive.

This does not overturn the cSCC manuscript, which compares RNA against RNA
within tissue and is unaffected. It goes at R2 of THIS paper.

**What it requires, and this now outranks section 6.** A compartment-matched RNA
axis for skin: metatranscriptomes from skin swabs, compared against the swab DNA
reservoir. Then the question becomes answerable rather than confounded:

- if swab RNA is productive and biopsy RNA is silent, the finding is that HPV
  activity is stratified by tissue depth, which is a better paper than the one
  planned and is a genuinely new statement about the reservoir;
- if swab RNA is silent in immunocompetent hosts and productive only in the
  RAG/EV patients, the decoupling survives and gains a host-status axis;
- either way the current R2, biopsy RNA against swab DNA, cannot be defended
  once a referee notices, and one will.

The 13 paired vaginal-metagenome specimens in PRJNA797778 (5.5) are the same
material class and give the mucosal counterpart for free.

### 5.5d RESULT: the compartment hypothesis is rejected. Host immunity is the variable.

Job 525959, completed 2026-08-06, 101 of 106 libraries typed. All are skin-swab
metatranscriptomes, so the sampling compartment is held constant and host immune
status is the only thing that varies. Results in `skin_swab_results/`.

| arm | typed | HPV+ | rate | productive | median depth |
|---|---:|---:|---:|---:|---:|
| healthy (SRP046152, pilosebaceous unit) | 19 | **0** | 0% | 0 | 1.63 Gb |
| immunocompromised (DOCK8 + RAG) | 66 | **52** | 79% | 37 | 0.95 / 17.48 Gb |
| ex-vivo skin models (SRP519513, excluded) | 16 | 0 | 0% | 0 | 18.40 Gb |

Complete separation, Fisher exact *P* = 1.9 × 10⁻¹⁰.

**It is not depth.** The two arms returning nothing are the more deeply
sequenced ones (healthy 1.63 Gb, skin models 18.40 Gb) while DOCK8, at 0.95 Gb,
is 85% positive. Depth does not drive positivity inside the immunocompromised
arm either: RAG at 17.48 Gb returns 50% against DOCK8's 85% at a twentieth of
the depth. An insensitivity artefact would have produced the opposite ordering
in both comparisons.

**It is the right virus.** Of 220 assigned type-calls in the immunocompromised
arm, 125 are Gammapapillomavirus and 87 Betapapillomavirus, i.e. 96% is the
cutaneous commensal repertoire (HPV208, HPV20, HPV12, HPV38, HPV17, HPV24,
HPV124, HPV107, plus novel mw/mSK/mKN candidate types). Only 4 alpha calls in 2
libraries, and **zero HPV16, HPV18, HPV31, HPV33 or HPV45**. This is not
cell-line or mucosal contamination; it is the beta/gamma reservoir itself,
transcribing.

**What it settles, and what it does not.** This is a control on a small and
atypical corner of the data, not the main argument. Transcriptomes are usually
run on tissue rather than on surface material, which is why the whole usable
healthy swab cohort came to 19 libraries; the evidence base for R2 is the
~10,849 skin RNA libraries, which are overwhelmingly tissue. What this shows is narrow and useful: surface
sampling per se does not reveal activity that tissue sampling misses, so the
swab-against-biopsy asymmetry is not the explanation for the silence.

The primary argument against the compartment explanation is the scale argument
in 5.2b, which needs no swab data: HPV presence is 17 to 25% in shotgun skin
DNA against 2 to 4% in skin RNA, and since transcription amplifies, an active
infection would make RNA the more sensitive detector rather than the less.

Together, R2 survives and gains a second axis:

> The cutaneous beta/gamma reservoir is transcriptionally silent in
> immunocompetent skin, and productivity is licensed by loss of host immune
> control rather than by anatomical site or sampling compartment.

DOCK8 deficiency is a positive control rather than an arbitrary cohort: it is
the primary immunodeficiency defined by cutaneous viral disease including HPV.

**How it must be worded.** These libraries are RNA only, with no DNA on the same
subjects (5.5c), so this is *absence of activity in immunocompetent skin*, not
*presence without activity*. The decoupling claim continues to rest on the
cross-cohort composition comparison. Stating it any more strongly would be
unsupported.

### 5.5c A paired skin test does not exist in public data. Stop looking.

Recorded so this is not re-attempted a fourth time. Every route to a paired
DNA/RNA test in skin was tried on 2026-08-06 and all of them fail, for one
structural reason.

| route | outcome |
|---|---|
| same-BioSample pairs, Shen-anchored | 3 specimens, all PRJNA872116 (RAG). Already typed. |
| same-BioSample pairs, skin niche | 604 found, but 327 with usable WGS are **all cultured cells**: 304 HipSci iPSC, 18 CCLE, 5 progeria fibroblasts. The non-WGS remainder adds CCLE, ENCODE and NCI PDX xenografts (PRJEB37099). Zero skin tissue. |
| donor-level linkage inside DOCK8 (SRP157584) | Impossible. All 879 DNA runs lack `host_subject_id`; all 54 RNA runs carry it. No join key, and the naming schemes differ (`VRNA_227` against `Dock8_001`). |

Two root causes, both structural rather than bad luck:

**Same-BioSample pairing selects for cell lines.** HipSci, CCLE, ENCODE and PDX
panels are precisely the projects that run both assays on identical material and
register it once. Clinical studies register DNA and RNA under separate
BioSamples. The cervical set, where only 2 of 64 pairs were cell lines, is the
exception rather than the rule; in skin the ratio is 100% the other way.

**Skin microbiome samples use MIMS environmental BioSample templates.** Those
carry `env_broad_scale`, `env_local_scale` and `env_medium` but have no slot for
subject linkage, and the `description` field is boilerplate ("MIMS
Environmental/Metagenome sample from human skin metagenome") for all 347 DOCK8
samples. The template is built to describe a site, not a patient. This is the
same wall that produced the non-human runs in the Shen reservoir, the horses and
sea stars in the metatranscriptome scoping, and the ex-vivo skin models in the
healthy arm.

Consequence for the paper: the skin DNA and RNA axes are cross-cohort and will
stay that way. R2 is a composition comparison, the swab-against-biopsy asymmetry
is a stated limitation rather than a solved problem, and the host-status contrast
(section 5.5b) is an RNA-only result that must be worded as absence of activity
rather than as presence without activity.

### 5.6 Triage still owed on the pair list

- **16 of the 64 cervical pairs have an AMPLICON DNA side.** If that amplicon is
  HPV genotyping it is ideal, a direct HPV DNA presence measure on the same
  specimen. If it is 16S it is useless here. Must be checked per project, not
  assumed.
- **The donor tier is unvalidated.** `donor_token` strips replicate suffixes off
  sample names, which is crude and will both over- and under-merge. The 58 to 59
  donor pairs are a lead list to read by hand, never a statistic.
- **Non-amniote entries in Shen's catalogue.** Among paired specimens alone the
  catalogue carries papillomavirus-positive libraries from rice, cacao,
  asparagus, sea bass, three butterflies and a springtail. Papillomaviruses
  infect amniotes only, so these are either false hits or human contamination in
  those libraries. Cow, pig and mouse also appear and are legitimate hosts, not
  errors. Do not call this k-mer noise without a protein-level check, which is
  the test that settled the deltavirus candidates in a sibling project. Either
  outcome is a small independent datapoint for the paper's central thesis about
  what a catalogue presence number means.

## 6. The cervical/vaginal run

6,396 curation-ready accessions, cell-line-term pre-filtered, none of them in the
current 12,200. This is the decisive piece and also the most dangerous cohort in
the archive for this question.

### Why it is dangerous
Cervical RNA-seq is **62% cell-line-tripped** (8,865 of 14,228), because cervical
carcinoma lines are among the most-sequenced human samples in SRA. Two specific
traps already hit in this project:

- The E1:E2 marker scored a naive public vaginal baseline at **78% integrated
  material**, and a random 1,000-library subsample of the same space at 39%.
  Re-stratifying that space does not fix it: HeLa is a cervical adenocarcinoma
  line, so any text search for cervix or cervical is saturated by it.
- **Control-arm conflation**: "normal" and "control" in SRA diagnosis fields can
  mean the untreated arm of a cell-culture experiment rather than healthy tissue.
  Caught by the E1:E2 marker, which scored such "normals" at 90% integrated.
  The cervical normal arm is where this will bite hardest.

Surviving this cohort cleanly is itself a contribution and should be written up
as one, not buried in methods.

### Operational notes
- Single-cell exclusion is mandatory and is now automatic
  (`classify_single_cell` in `bin/cell_line_patterns.py`). In a comparable
  cohort 301 of 1,000 runs were single-cell but carried over half the bases, and
  all 18 largest runs were single-cell. Exclude, do not cap reads.
- Use the chunked driver. Per-chunk work dirs have reached 0.7 to 1.3 TB and the
  home quota is the binding constraint; a failed chunk is fixed by removing its
  work dir and resuming, not by relaunching.
- `run_chunked.sh` silently ignores a new samplesheet if `.split_done` exists in
  the outdir. Use a fresh outdir.
- Submit through SLURM on rocket (login1 or login2), never on the login node.
  Scripts go in `hpc/`, logs in `hpc/logs/`.

### Order of work
The scoping recommendation is to start with the clean niches, and that still
holds even with the full run approved, because it validates the pipeline on the
mucosal material before committing to the hard cohort:

1. Anorectal 651 (6.6% cell-line-tripped) plus tiering the existing 681.
2. Vaginal metatranscriptome 2,837 (cell-line-free by construction).
3. Cervical/vaginal 6,396 with per-sample curation at full strength.

## 7. Open questions

- Does R5 have enough left in it once paper 1 has published the marker, or does
  it become a short section that mainly cites paper 1?
- Is the gut arm worth broadening beyond the tight anorectal terms? It would
  require the colorectal cell-line patterns in
  `../rna_extension_scoping/proposed_colorectal_patterns.py`, which are not
  applied. Probably no, since colorectal is not a papillomavirus niche and the
  cost is a large new contamination surface.
- Whether to report prevalence at all (5.2) or stay with composition and state
  the limitation. Depends on how well Logan's per-taxon coverage can be
  established.
