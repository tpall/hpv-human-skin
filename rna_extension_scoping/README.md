# Extending the RNA-activity survey to vaginal/cervical and anorectal niches

Scoping for a concrete next step: our RNA-seq survey is skin-only, so the DNA-vs-RNA
decoupling (`../dna_rna_triangulation/`) has a matched activity axis for skin but not
for the other two body sites where the Shen/Logan DNA reservoir is informative
(vaginal = α mucosal, gut/anorectal = mixed). This directory scopes the SRA RNA-seq
that would build those axes, with the cell-line contamination burden quantified up
front because new tissues drag in new cell lines.

## Headline

| curation-ready set (cell-line-excluded)          | runs  | already screened | **NEW** |
|--------------------------------------------------|-------|------------------|---------|
| Vaginal/cervical RNA-seq                          | 6,837 | 441              | **6,396** |
| Anorectal RNA-seq                                 |   651 |   0              | **651**   |
| Vaginal metatranscriptome (cell-line-free niche)  | 2,837 |   -              | ~2,837    |

The `anogenital` query already exists in `bin/query_sra.py` but was never run at scale
— only 441 cervical/vaginal libraries are in the current 12,200, and zero anorectal.
So ~7,000 new libraries turn up (before per-sample QC).

## Cell-line burden (the "be very careful" part)

Cervical/vaginal RNA-seq is **62% cell-line-tripped** in metadata (8,865 / 14,228):
cervical carcinoma lines are among the most-sequenced human samples on SRA. The good
news: the lines that actually appear are already covered by `bin/cell_line_patterns.py`
— HeLa (5,762 in the cervical context), SiHa (302), CaSki (85), C-33A (37), ME-180 (8),
W12 (14). Candidate gaps we probed (DoTc2, SNU-703/1005, GH354) return 0, so no cervical
additions are needed.

Anorectal is **6.6% cell-line-tripped** (37 / 560) — anal SCC lines are rare, so this
niche is clean and small: the better first target.

**The real new gap is colorectal**, which only matters if the gut query is ever
broadened beyond the tight anorectal terms. None of these are in the patterns yet, and
they are heavily sequenced (counts = libraries in a broad gut RNA-seq context):
HCT116 9,920 · SW480 2,032 · Caco-2 1,797 · HT-29 883 · DLD-1 881 · SW620 794 ·
RKO 750 · LoVo 379 · T84 273 · LS174T 175 · HCT-15 128 · COLO 205 41. A proposed
`NAMED_CELL_LINES` addition block is in `proposed_colorectal_patterns.py` — apply it
only if broadening to gut; the tight anorectal query below does not need it.

### Two caveats on the `NOT <cell-line-terms>` pre-filter

1. It **over-excludes**: it drops legitimate clinical studies that merely mention a
   line in their description ("compared against HeLa"), so the true clean clinical set
   is likely larger than the counts above.
2. It **under-excludes**: it only sees metadata text, so cell lines with uninformative
   per-run metadata survive it — which is exactly why the per-sample classifier
   (`flag_cell_lines.py`) and the downstream coverage/contamination checks
   (`flag_contamination.R`, the SNP cross-contamination step) must still run on every
   candidate. This pre-filter is a coarse first pass, not the curation.

## Curation-ready queries (drop into `bin/query_sra.py` TISSUE_QUERIES)

```python
# Tightened vaginal/cervical (mucosal alpha niche). Pair with BASE_FILTER
# (Homo sapiens + RNA-Seq + public) as query_sra.py already does.
"vaginal_cervical": [
    '("cervix"[MeSH Terms] OR cervical OR cervix OR ectocervix OR endocervix '
    'OR vaginal OR vagina OR vulvar OR vulva OR "cervical intraepithelial" '
    'OR CIN OR condyloma)',
],

# Anorectal subset of the gut (the HPV-relevant part; excludes colon/IBD/organoid).
"anorectal": [
    '(anal OR anorectal OR perianal OR "anal canal" OR "anal squamous" '
    'OR "anal intraepithelial" OR AIN OR "anal cancer" OR "anal carcinoma" '
    'OR "anal condyloma")',
],
```

Vaginal metatranscriptome is a separate, cell-line-free bonus source (microbiome RNA
of the vaginal community, organism = `vaginal metagenome`); worth a dedicated pass.

## Files here

- `candidate_new_vaginal_cervical.srr.txt` — 6,396 new run accessions (cell-line-term
  pre-filtered) not in the current 12,200. Feed to the pipeline, then run the full
  per-sample cell-line + contamination curation.
- `candidate_new_anorectal.srr.txt` — 651 new anorectal run accessions, none screened.
- `sra_tight_scope.py` — reproduces the counts (NCBI eutils, no download).
- `fetch_new_ids.py` — reproduces the run-accession fetch + diff against the existing set.
- `proposed_colorectal_patterns.py` — colorectal `NAMED_CELL_LINES` additions, apply
  only if broadening to gut.

## Recommendation

Start with **anorectal (651, clean, untouched)** as a small tractable pilot and the
**vaginal metatranscriptome (2,837, cell-line-free)** — both sidestep the 62% cervical
cell-line problem. Then take on the full cervical/vaginal 6,396 with the per-sample
curation running at full strength, since that is where the mucosal-α productive niche
lives and where the DNA-frequent vaginal types (HPV91/51/66/67/16) can be tested for
activity.
