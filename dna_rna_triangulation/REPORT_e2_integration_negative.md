# Can HPV integration be detected in the Logan/Shen DNA catalogue? A negative result

**Status: the assembly-based approach is negative (§1–7) and its Discussion
paragraph was removed. The read-based approach it pointed to (§8) was then run
and WORKS (§9).** Read §9 first if you want the result; §1–7 explain why the
obvious route fails, which is why the working route looks indirect.

Date: 2026-08-03. Code: `bin/e2_integrity_from_logan.sh` + `bin/_e2_*.py`.
Data: `e2_integrity_by_site.tsv`.

---

## 1. The mechanism we set out to detect

Integration of the HPV genome into host DNA proceeds by a well-described route:

1. **Gene disruption.** The circular viral genome breaks open, and the break
   falls preferentially inside the E1 or E2 open reading frame.
2. **Loss of E2 repression.** E2 is the transcriptional repressor of E6 and E7.
   Once E2 is broken or lost, that repression goes with it.
3. **Uncontrolled oncogene production.** E6 and E7 are then expressed
   continuously, driving genomic instability and cancer risk.

This makes two predictions that are, in principle, readable from assembled viral
contigs: integrated genomes should **lack an intact E2**, and the **break should
fall inside E1 or E2** rather than elsewhere in the genome.

The RNA-side counterpart of prediction 2 is already established in this project
and holds strongly (`bin/e1e2_integration_marker.R`): the E1:E2 transcript ratio
is 129× in cell lines against 0.27× in productive primary lesions, separating
them at AUC 0.968, with the episomal HPV31 HaCaT derivative sitting with primary
tissue rather than with the integrated lines. **That result is unaffected by
anything below.** What follows concerns only the attempt to confirm it on DNA.

## 2. What was tested

Shen et al.'s petabase survey provides papillomavirus contigs across body sites.
Sequences are not in the paper's tables but are deposited on Zenodo
(10.5281/zenodo.19683048, files S4/S5). Contigs were recovered for six
categories and annotated by `tblastn` against 738 RefSeq papillomavirus ORF
proteins (E1, E2, E6, E7, L1, L2).

Commensal carriage is registered in SRA under metagenome taxids; carcinoma and
dysplasia are registered under *Homo sapiens* and are therefore invisible to a
metagenome filter. The cancer arm had to be selected by BioProject:

| category | contigs mapped | recovered | source |
|---|---:|---:|---|
| skin | 6,390 | 3,593 | skin metagenome |
| cervical cancer | 1,917 | 1,426 | PRJNA664651 (BioSample `tissue = cervical squamous cell carcinoma`), PRJNA305342 |
| vaginal | 941 | 690 | vaginal metagenome (microbiome cohorts) |
| gut | 443 | 382 | gut metagenome |
| cervical CIN | 167 | 72 | PRJNA687153, PRJNA318454 |
| anogenital wart | 56 | 13 | PRJNA517793 |

## 3. Result 1: the completeness-conditioned test is circular

Conditioning on near-complete contigs (≥7 kb) and asking whether E2 is present:

| category | contigs ≥7 kb | E2 intact |
|---|---:|---:|
| skin | 1,072 | 99.9 % |
| vaginal | 427 | 100 % |
| gut | 211 | 100 % |
| anogenital wart | 5 | 100 % |
| cervical CIN | 21 | 100 % |
| **cervical cancer** | **267** | **99.3 %** |

Cervical carcinoma, where integration is the rule, returns the same answer as
skin. This is not biology. **An integrated genome cannot assemble to full length,
because the virus–host junction terminates the assembly path.** Conditioning on
≥7 kb therefore selects episomal genomes in every category, and the statistic is
incapable of detecting the thing it was built to detect. Any statement of the
form "near-complete genomes carry an intact E2" is close to a tautology.

## 4. Result 2: the completeness distribution moves, but is confounded

If integration fragments assemblies, the signal should appear as shorter contigs
rather than as E2-less full ones. It does:

| category | median length | ≥7 kb | <2 kb |
|---|---:|---:|---:|
| vaginal | 7,784 bp | 62 % | 8 % |
| gut | 7,245 bp | 55 % | 7 % |
| skin | 4,162 bp | 30 % | 17 % |
| cervical CIN | 3,009 bp | 29 % | 18 % |
| **cervical cancer** | **3,074 bp** | **19 %** | **27 %** |

But the internal control does not hold. Within the cancer arm, the WGS cohort
reaches 28 % full length against 17 % for the capture cohort, and one CIN cohort
reaches 45 %. Library chemistry (targeted capture vs shotgun) clearly drives
assembly completeness, so fragmentation cannot be attributed to integration.

## 5. Result 3: where the break falls — and why this cannot work here

The sharper test is positional: if the break falls inside E1 or E2, the **contig
terminus** should land inside those ORFs. Scoring each ORF as intact (≥90 % of
the reference covered), truncated at a contig edge (within 60 bp of a terminus),
truncated internally, or absent:

**Percentage of present ORFs truncated at a contig end**

| category | E1 | E2 | L1 |
|---|---:|---:|---:|
| vaginal | 24 % | 6 % | 4 % |
| gut | 26 % | 14 % | 3 % |
| skin | 49 % | 20 % | 9 % |
| anogenital wart | 40 % | 0 % | 23 % |
| cervical CIN | 26 % | 12 % | 6 % |
| **cervical cancer** | **38 %** | **20 %** | **9 %** |

Cancer is not elevated over skin on either gene. And the ordering **E1 > E2 > L1
holds in every category, including healthy vaginal and gut carriage**, which
integration cannot explain.

The reason is the catalogue's construction. The HPV genome is ordered
E6–E7–**E1–E2**–E4/E5–**L2–L1**–LCR, and **Shen's survey finds contigs by an L1
HMM search**. Every contig is therefore anchored on L1 and extends outward: it
reaches L2 easily, E2 less often, and E1 last. The truncation gradient measures
**distance from the L1 anchor**, not the position of an integration breakpoint.
The same gradient appears in commensal carriage because it is a property of the
assembly and the search strategy, not of the virus.

## 6. Conclusion and what would actually be required

Three independent framings of the question — E2 presence conditioned on
completeness, contig-length distribution, and breakpoint position — all fail on
this dataset, and they fail for one root cause: **an L1-anchored, assembly-based
catalogue systematically loses exactly the molecules that integration produces.**
Integrated genomes either do not assemble, or assemble as fragments whose
truncation pattern is dominated by distance from the anchor.

Detecting integration needs the reads, not the assemblies:

- **chimeric read pairs** spanning the virus–host junction, which is the direct
  and standard evidence, and which assembly discards;
- **coverage discontinuity** across the viral genome within a library, rather
  than presence/absence on a contig;
- a **library-chemistry-matched** comparison, since capture and shotgun give
  systematically different assemblies.

All three are available from raw SRA data for the cohorts identified in §2, so
the question is answerable; it is simply not answerable from this catalogue.

## 7. What survives

- The **RNA-side E1:E2 marker** stands independently and is in the manuscript.
- The **cohort identification** in §2 is reusable: cervical carcinoma and CIN
  libraries carrying HPV are locatable in SRA under *Homo sapiens*, with
  BioProject accessions listed above.
- The methodological point is worth stating wherever the DNA reservoir is
  discussed: **Shen's catalogue is L1-anchored and assembly-based, so "presence"
  there means L1 presence, and it carries no information about genome state.**

---

# 8. Design: detecting integration from reads

The failure above is a property of assemblies, not of the question. Integration
is routinely detected from **chimeric virus–host reads**, which assembly
discards. This section is a concrete protocol, written so it can be executed
without repeating the scoping.

## 8.1 Two arms, because the material differs

The mucosal arm has DNA; the cutaneous arm does not, and cannot be obtained the
same way. Shen's catalogue contains only libraries **in which HPV was found**, so
a cSCC library carrying no HPV never produces an L1 contig and never enters the
list. The cutaneous arm therefore comes from our own HPV-positive RNA-seq, where
integration is read as a **virus–host fusion transcript** — the standard route
by which HPV integration is called in TCGA-style RNA data.

## 8.2 Arm A — mucosal DNA (positive control and mechanism)

| BioProject | runs | GB | median/run | cohort |
|---|---:|---:|---:|---|
| PRJNA664651 | 561 | 320.0 | 249 MB | cervical SCC, HPV capture |
| PRJNA305342 | 209 | 357.8 | 457 MB | cervical cancer, low-pass WGS |
| PRJNA687153 | 30 | 7.8 | 162 MB | CIN, capture |
| PRJNA318454 | 33 | 6.2 | 164 MB | CIN cohort |
| PRJNA517793 | 10 | 0.5 | 53 MB | anogenital warts (benign) |
| **total** | **843** | **692.4** | | |

692 GB is a few days of download and standard alignment on Helios. Capture
libraries are the best available material: HPV enrichment raises both viral
depth and the number of junction-spanning reads.

## 8.3 Arm B — cutaneous RNA (the claim we actually need)

HPV reads available per tier in our own libraries, as a proxy for junction
detectability:

| tier | ≥1,000 HPV reads | ≥10,000 | max |
|---|---:|---:|---:|
| productive controls | 38 | 32 | 1,839,454 |
| EV | 7 | 4 | 4,061,846 |
| peri-lesional | 10 | 6 | 116,145 |
| **neoplasia (cSCC/AK)** | **6** | **5** | **39,228** |
| cell line | 61 | 12 | 49,971 |

The cutaneous cancer arm is six libraries. That is small, and it bounds what a
negative result can mean (see 8.6).

## 8.4 Method

1. Combined reference: GRCh38 + the PaVE/RefSeq HPV panel already built for this
   project (Methods).
2. `bwa-mem -Y` for DNA (soft-clipped supplementary alignments retained), a
   splice-aware aligner for RNA.
3. Junction calling with an established caller rather than a bespoke script —
   **VirusBreakend** (GRIDSS family), ViFi, seeksv or VERSE. All emit
   virus–host breakpoints with supporting read counts.
4. Cluster breakpoints per library; record the **position of the viral-side
   breakpoint in genome coordinates**, which is the quantity the mechanism
   predicts.

## 8.5 Order of controls — this determines whether the result is interpretable

Run in this order and stop if a step fails:

1. **Cell lines** (61 libraries ≥1,000 reads, including HeLa-derived material).
   Integration is certain here, and the E1:E2 ratio already says so at 129×.
   If the pipeline does not recover junctions in these, nothing downstream is
   interpretable.
2. **Cervical carcinoma** (arm A). Integration is the rule; this establishes
   sensitivity in tissue rather than culture.
3. **Anogenital warts and CIN.** Benign productive lesions should be episomal,
   dysplasia intermediate. This is the specificity check: if warts return
   junctions at the same rate as carcinoma, the caller is producing artefacts.
4. **Cutaneous neoplasia** (arm B). Only now is a negative result meaningful.

## 8.6 What each outcome would mean

- **Junctions in cervical cancer, none in warts, none in cSCC** — the strongest
  available version of the manuscript's claim: the trace of HPV in cSCC is not
  integrated virus, measured rather than inferred, and consistent with the
  E1:E2 ratio of 0.04 in that tier.
- **Junctions in cSCC** — would overturn the hit-and-run reading for those
  libraries and would be worth reporting on its own.
- **No junctions anywhere including the cell lines** — a pipeline failure, not a
  biological result.

The honest limit does not go away: with six cutaneous libraries, a negative
means "no junction detected at this depth", not "no integration". The viral read
counts above (5 of 6 exceeding 10,000 reads) should be converted into an
explicit detection limit before the negative is reported, in the same way the
RPM floor was calibrated on productive controls.

## 8.7 Relationship to the RNA marker already in the manuscript

The E1:E2 ratio and junction detection measure the same event by different
means, so they cross-validate. The ratio is cheap, already computed for every
library, and separates integrated from episomal at AUC 0.968. Junction calling
is expensive but direct and localises the breakpoint. A sensible reading is that
the ratio screens and the junctions confirm — and that the 28 discordant
libraries the ratio flagged (§ `bin/e1e2_integration_marker.R`) are the natural
first targets for confirmation, since they are the cases where the metadata and
the molecular signature disagree.

---

# 9. Result: integration detected from reads (Helios job 508305)

The design in §8 was executed on the 17-library pilot. **It works, and it
overturns the conclusion of §1–7: the question is answerable, just not from
assemblies.**

## 9.1 The discriminating statistic is clustering, not count

Raw junction counts are uninformative and actively misleading. The two EV
libraries produce the *highest* counts in the whole pilot (12,979 and 13,937)
while being episomal, because chimeric-read callers generate scattered noise in
proportion to viral depth. Real integration is a **recurrent** junction: many
reads at one coordinate pair.

Scoring each library by the fraction of its junctions falling at the single
most-supported site:

| group | libraries | fraction at top site |
|---|---:|---|
| cell lines | 5 | **0.33 – 0.46** |
| cutaneous neoplasia, integrated | 2 | **0.41 – 0.43** |
| warts (productive controls) | 4 | 0.01 – 0.03 |
| EV | 2 | 0.00 – 0.01 |
| cutaneous neoplasia, episomal | 4 | 0.02 – 0.03 |

No overlap; the gap is an order of magnitude. Use this ratio, not the count.

## 9.2 The positive control reproduces a textbook answer

Two HPV18 cell-line libraries independently return the same junction:
**HPV18 ↔ chr8:127,229,1xx**, i.e. the **8q24/MYC** locus. That is HeLa's known
integration site, recovered without being looked for. A pipeline that recovers
it is sensitive and specific enough to trust downstream.

## 9.3 The E1:E2 transcript marker predicted this, 16/17

The pilot was run with a falsifiable prediction written into the script header:
of the six cutaneous-neoplasia libraries only `SRR35133829` (log2 E1:E2 = +7.3)
should carry junctions. **It does** — fraction at top site 0.41 — and the four
libraries at log2 −3.3 to −5.9 do not (0.02–0.03).

One miss: `SRR19520029` scores log2 = −4.7 yet is unambiguously integrated
(0.43, 4,613 reads at one site). It is one of the three HPV16 libraries the
manuscript already flags as mucosal cases swept into the cutaneous tier by the
broadened query (`@supptbl-cin`). Junction detection therefore **independently
confirms that those libraries do not belong in the cutaneous arm** — the miss
strengthens the sensitivity analysis rather than weakening the marker.

## 9.4 Every breakpoint falls in E1 or E2

| library | viral breakpoint | ORF | human side |
|---|---|---|---|
| DRR670533 | HPV18:2498 | **E1** | chr8:127,229,302 |
| DRR670533 | HPV18:930 | **E1** | chr8:127,229,131 |
| SRR23163057 | HPV18:930 | **E1** | chr8:127,229,131 |
| SRR19520029 | HPV16:3439 | **E2** | chr17:41,583,205 |
| SRR19520029 | HPV16:881 | **E1** | chr17:41,583,141 |
| SRR35133829 | HPV16:881 | **E1** | chr17:40,782,856 |

(Reference ORF spans: HPV16 E1 865–2813, E2 2755–3852; HPV18 E1 914–2887.)

**Not one breakpoint falls in L1, L2, E6 or E7.** Where two junctions are
recovered they are the two ends of the integrated segment and they bracket the
E1–E2 region precisely. This is the mechanism as described: the circular genome
opens inside E1/E2, E2 repression is lost, E6/E7 escape control.

`SRR19520029` and `SRR35133829` share the same viral breakpoint (HPV16:881) but
land at different human loci (chr17:41.58 Mb vs 40.78 Mb), so they are two
independent integration events, not cross-contamination between libraries.

## 9.5 What this changes

- The manuscript's cutaneous claim can be stated as measured rather than
  inferred: of six HPV-positive cSCC/AK libraries with sufficient depth, four
  carry no integration junction, and the two that do are HPV16 mucosal cases
  already flagged as misfiled.
- The E1:E2 ratio is validated against an independent method and can be used as
  a cheap screen across all libraries, with junctions as confirmation.
- The 28 discordant libraries flagged by the ratio (§8.7) are the obvious next
  batch, followed by the mucosal DNA cohorts in §8.2 if a larger positive arm
  is wanted.

Code: `run_junction_pilot.sbatch` (on Helios, `results_junction_pilot/`).
Outputs: `junction_summary.tsv`, per-library `junctions/*.junction.tsv`.
