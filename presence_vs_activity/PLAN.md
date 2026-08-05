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
Status: DONE but needs the fix in section 5.1. The DNA-abundant beta/gamma types
are transcriptionally silent; the productive types are alpha wart viruses that
are DNA-rare or absent. Only HPV2 occupies the abundant-and-productive corner.

**R3. Cervix and vagina: presence does predict productivity.**
Status: NEEDS THE RUN (section 6). The within-tissue join that R2 could only do
for skin: are the DNA-frequent vaginal types (HPV91/51/66/67/16/90) the
productive ones? Prediction is yes, which makes the cervix the clean contrast
case to skin. Carries a second axis: high-risk alpha dominating productive
expression by lesion grade while low-risk and commensal mucosal types stay quiet.

**R4. Anorectal: both genera present in DNA, only one productive.**
Status: PILOT (681 typed, 230 HPV+ at 34%, productive signal alpha and HPV16-led
at 91 libraries). This removes the between-site confound: the cutaneous beta
types HPV24 and HPV124 that appear in anal and gut DNA are entirely silent in
the anal transcriptome. The strongest figure in the paper.

Honest limit: this is within-*niche*, not within-*sample*. The DNA is gut
metagenome and the RNA is anorectal RNA-seq, different libraries. See 5.5.

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

### 5.5 Scope paired DNA and RNA from the same specimen
Highest value per unit effort in this list, and new. Search SRA for BioSamples
carrying both a DNA library (WGS, metagenome, capture) and an RNA-seq library
from the same anal or cervical biopsy. If even a few dozen exist, R4 converts
from a composition comparison into a genuine paired test, which is a different
tier of evidence and would likely become the paper's centrepiece. If none exist,
say so in the limitations and move on.

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
