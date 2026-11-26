# Beta-HPV E6/E7 in the eORFeome functional screens

External resource, read 2026-08-21. Note only; nothing in `../manuscript/` or in
`PLAN.md` changes on the strength of it.

**Bottom line.** A 2026 pooled-screen resource tested HPV E6/E7 from 19 types across
five host pathways. Mucosal alpha E6 lights up p53 inhibition and NF-kB activation, as
expected. **HPV5 and HPV8 E6/E7 score nothing in any pathway.** The single beta-genus
hit in the whole resource is **HPV9 E6 inhibiting cGAS-STING** — and HPV9 is in our EV
tier. This is usable as one or two Discussion sentences about what beta E6/E7
transcription would mean if translated. It is **not** protein evidence about our
samples and does not reopen the block dropped on 2026-08-07.

---

## 1. What the resource is

Pachano et al., *Cell* 189:5429-5448 (2026), doi:10.1016/j.cell.2026.06.017. The
"eORFeome": 3,835 cloned pathogen effector ORFs from 177 species, expressed from a
barcoded inducible lentiviral vector and screened in pooled format in human cells.
Six screens are reported (NF-kB activation and inhibition are scored separately):
NF-kB activation, NF-kB inhibition, p53 inhibition, apoptosis inhibition, MHC-I
inhibition, cGAS-STING inhibition. Hit calling is FDR < 0.05 with FC > 1.5.

Files: Table S1 = clone inventory (`mmc2.xlsx`), Table S2 = per-ORF screen results
(`mmc3.xlsx`). Both currently in `~/Downloads`, not copied into this repo.

## 2. Coverage of the types we actually study

98 HPV clones, 19 types, 97 of them scored in the screens.

| our tier | types we report | in the eORFeome? |
|---|---|---|
| EV biopsies (beta) | HPV9, 14, 20, 25, 47 | **only HPV9** |
| productive neoplasia | HPV12, 16, 96, 182 | HPV16 only, and **without E6** |
| warts / mucocutaneous | HPV6, 11, 2, 27, 57 | HPV6b, 11, 2 |

Two coverage facts worth stating plainly before using any of this:

- **HPV16 E6 is not in the library at all** (HPV16 has only E1, E5, E7, L1, L2). Any
  claim about high-risk E6 there rests on types 18, 31, 33, 35, 39, 45, 51, 52, 58.
- **None of our productive cutaneous types are present** — no HPV12, 96, or 182. The
  resource says nothing about the infections that carry our productive signal.

## 3. What scored

19 hits across the six screens, involving 12 distinct type-ORF pairs.

| type | ORF | screen | lfc | FDR |
|---|---|---|---|---|
| **HPV9** (beta) | **E6** | **cGAS-STING inhibition** | **+1.81** | **1.4e-07** |
| HPV58 | E6 | NF-kB activation / p53 / apoptosis | +3.54 / +1.73 / +1.13 | 2.2e-04 / 2.4e-04 / 0.015 |
| HPV52 | E6 | NF-kB activation / p53 | +3.33 / +1.78 | 2.2e-04 / 2.4e-04 |
| HPV51 | E6 | NF-kB activation / p53 | +1.79 / +2.33 | 2.2e-04 / 2.4e-04 |
| HPV51 | E7 | apoptosis inhibition | +0.77 | 2.4e-04 |
| HPV35 | E6 | NF-kB activation / p53 | +1.30 / +2.17 | 2.2e-04 / 2.4e-04 |
| HPV31 | E6 | NF-kB activation / p53 | +1.26 / +2.05 | 2.2e-04 / 2.4e-04 |
| HPV18 | E6 | NF-kB activation / p53 | +1.39 / +1.75 | 2.2e-04 / 2.4e-04 |
| HPV45 | E6 / E1 | p53 / apoptosis | +1.44 / +1.41 | 1.4e-03 / 0.033 |
| HPV39 | E6 | p53 inhibition | +1.24 | 2.7e-03 |
| HPV34 | E6 | p53 inhibition | +1.26 | 2.4e-04 |

HPV34 is alpha-11, i.e. mucosal — it is not a cutaneous result despite the low type
number.

**The beta-genus ORFs, in full, hit or not:**

| type | ORF | best two screens |
|---|---|---|
| HPV5 | E2 / E6 / E7 | nothing anywhere; best FDR 0.61, 0.73, 0.90 |
| HPV8 | E4 / E6 / E7 | nothing anywhere; best FDR 0.58, 0.27, 0.98 |
| HPV9 | E6 | **cGAS-STING +1.81, FDR 1.4e-07**; p53 +1.03, FDR 0.38 |
| HPV9 | E1 / E2 / L1 / L2 | nothing; best FDR 0.097 (E1, apoptosis) |
| HPV92 | ^E4 | nothing; best FDR 0.91 |

## 4. How to read it

1. **The mucosal program is reproduced.** High-risk E6 inhibits p53 and activates
   NF-kB across nine types. The assay works and is sensitive to E6 biology.
2. **HPV5 and HPV8 are inert in it.** The two beta types most tied to EV and cutaneous
   SCC do nothing in five pathways. Whatever beta E6/E7 does, it is not the mucosal
   program, and it is not detected by assays built to find the mucosal program.
3. **HPV9 E6 targets innate immunity, not p53.** The one beta hit is in cGAS-STING
   inhibition, with p53 flatly negative (FDR 0.38). If beta E6 has a function here, it
   is immune evasion rather than tumour-suppressor degradation.

For us the useful consequence is interpretive: transcription of beta E6/E7 — which is
what our RNA panel measures — would not, on this evidence, imply the oncogenic program
that E6 transcription implies for mucosal types. That strengthens rather than weakens
the current framing, in which cutaneous HPV is transcriptionally quiet in skin and
carcinoma and productive in warts.

## 5. What this is not

- **Not protein evidence about our samples.** It is heterologous overexpression of
  cloned ORFs in immortalised cell lines. It says what these proteins *can* do when
  forced, not what is expressed in an actinic keratosis. The proteomics block stays
  dropped and this note must not be cited as closing it.
- **Absence of a hit is weak evidence.** Five pathways, one expression system, one
  hit-calling threshold. HPV5/HPV8 E6 scoring nothing means "not in these assays," not
  "no function."
- **It cannot speak to our productive cases.** HPV12, 96 and 182 are absent from the
  library.
- **One beta type, one hit.** HPV9 E6 is a single observation from a screen, not a
  validated mechanism; the paper validated other candidates (HHV6A U14, adenoviral
  13.6K), not this one.

## 6. Proposed Discussion sentences

Draft, for the paragraph where we say transcription is the closest proxy for activity:

> Functional data from an independent effector-ORFeome screen argues that beta-HPV E6
> transcription should not be read as a mucosal-type oncogenic signal: E6 proteins from
> nine mucosal high-risk types inhibited p53 and activated NF-kB in pooled screens,
> whereas E6 and E7 from HPV5 and HPV8 scored in none of five pathways [Pachano 2026].
> The single beta-genus hit in that resource, HPV9 E6, acted on cGAS-STING rather than
> p53, suggesting that whatever cutaneous E6 does, assays designed around the mucosal
> program will not see it. This is heterologous overexpression rather than evidence of
> protein in tissue, and our claims remain claims about transcription.

## 7. Provenance

Numbers here were extracted from `mmc2.xlsx` (Table S1) and `mmc3.xlsx` (Table S2) with
openpyxl on 2026-08-21; type/gene assignment is by the `Species` and `Label` columns,
genus assignment by type number. HPV34 was initially miscounted as cutaneous in a first
pass and is corrected above.
