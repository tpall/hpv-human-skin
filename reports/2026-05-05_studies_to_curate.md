# Studies to manually curate — 2026-05-05

Worklist for verifying sample provenance against the heuristic classification used in the intermediate report. Generated from the 16 studies that contribute typed HPV+ samples (n = 31, most permissive thresholds) in the in-progress full_v2 run.

For each study, decide whether the auto-assigned category (cell_line / engineered / clinical) is correct. Notes column flags the specific things to look for. ENA / SRA links open the study landing page.

## Priority 1 — alpha-type hits (HPV16 / HPV18) currently classified as **clinical**

These are the cases where unflagged in-vitro material would most distort the headline result. Top suspect: anything where the title or BioSample text reveals an experimental construct (sh/si/sg-RNA, lentivirus prefix, ASO, stable transfection) that the keyword list missed.

| ☐ | Study | n | Sample(s) | Ref | Breadth | Depth | Title / tissue text | What to confirm |
|---|---|---|---|---|---|---|---|---|
| ☐ | [ERP141290](https://www.ebi.ac.uk/ena/browser/view/PRJEB55636) | 2 | ERR10293422, ERR10293427 | HPV16 | 0.30 / 0.37 | 1.4 / 2.8 | titles `H04E`, `H05E`; tissue field empty | Is this a primary-tissue study or a cell-line / organoid study? Titles are uninformative — read the abstract. |
| ☐ | [SRP536293](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP536293) | 2 | SRR32387020, SRR32387021 | HPV18 | 0.25 / 0.18 | 1.0 / 0.5 | `DF H2AZ2 sh3`, `DF H2AZ2 sh2`; tissue `skin epidermis` | **Strong suspect.** `sh2` / `sh3` are numbered shRNA constructs against H2AZ2 — engineered. Heuristic deliberately doesn't match `sh\d+` (would false-positive on SH3 protein domains). Confirm and reclassify as engineered. |

## Priority 2 — beta / gamma clinical hits with in-vitro markers in the title

These hits are biologically plausible (β/γ types are cutaneous), but the title or tissue text contains experimental-construct hints that are worth a closer look.

| ☐ | Study | n | Sample(s) | Ref | Title / tissue text | What to confirm |
|---|---|---|---|---|---|---|
| ☐ | [SRP539392](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP539392) | 1 | SRR31039970 | HPV80 (depth 21!) | `FB_PLVX_80_Rep1`; tissue `Dermal-Connective` | **Strong suspect.** `PLVX` = pLVX, a Clontech lentiviral expression vector. Almost certainly engineered fibroblasts. Confirm and consider adding `pLV[A-Z]\\w+` (Clontech vector family) to the engineered patterns. |
| ☐ | [ERP149729](https://www.ebi.ac.uk/ena/browser/view/PRJEB67769) | 2 | ERR11758160, ERR11758161 | HPV208, HPV150 | `DRPLA-ASO2`; tissue `fibroblast` | DRPLA = dentatorubral-pallidoluysian atrophy patient line; ASO2 = antisense oligonucleotide #2. This is a primary patient fibroblast line under ASO treatment. Probably engineered (ASO is a manipulation), not "clinical" in the fresh-tissue sense. Decide whether ASO treatment counts as `is_engineered`. |

## Priority 3 — clinical hits that appear genuinely clinical

Sanity check that nothing in the study description undermines the clinical reading. These should be quick yes/no confirmations.

| ☐ | Study | n | Sample(s) | Ref | Title / tissue text | Notes |
|---|---|---|---|---|---|---|
| ☐ | [SRP498743](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP498743) | 1 | SRR28503328 | HPV190 (breadth 0.66, depth 2.4) | `SmartSticker, Forehead, Volunteer 5` | Looks like a wearable-sensor / volunteer skin study. Highest-confidence γ-type clinical hit in the cohort. |
| ☐ | [SRP548651](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP548651) | 1 | SRR31566149 | HPV182 (breadth 0.22, depth 7.0) | `NR_SK_19_012_Baseline_Uninvolved_..._S91`; tissue `Uninvolved` | Baseline / uninvolved-skin sample from a patient cohort (likely a dermatology study). |
| ☐ | [SRP565659](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP565659) | 3 | SRR32454238, SRR32454251, SRR32454283 | HPV65, 134, mw20c10anr | `WU1211_skin`, `WU1126_skin`, `WU1645_skin`; tissue `skin` | Three patient skin samples from a single study, three distinct cutaneous types — strong β/γ diversity signal in a clinical cohort. |
| ☐ | [SRP592803](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP592803) | 1 | SRR34022308 | HPV-mTVMBSFc09nr | title NA; tissue `Giant congenital melanocytic nevus tissue` | Real lesion tissue — high-confidence clinical. |
| ☐ | [SRP659354](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP659354) | 1 | SRR36658897 | HPV195 (breadth 0.06, depth 0.27) | `Pt 7 non-drainning tunnel`; tissue `Skin` | Likely hidradenitis suppurativa or similar wound / sinus-tract sample. Confirm. |

## Priority 4 — already-flagged samples (cell_line / engineered)

Confirm the heuristic call was right. Mostly a formality but quick.

### 4a. Flagged as `cell_line` (10 samples in 5 studies)

| ☐ | Study | n | Title pattern | Matched on |
|---|---|---|---|---|
| ☐ | [DRP013144](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DRP013144) | 1 | `MNT-1_WT1` (DRR670531) | MNT-1 (melanoma) — newly added pattern; HPV18 at depth 432 |
| ☐ | [SRP516327](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP516327) | 2 | `TIGK 24h-1`, `TIGK 6h-1` | TIGK (gingival keratinocyte line) |
| ☐ | [SRP554912](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP554912) | 3 | `SK-MEL-28-siLINC01291`, `A375-siLINC01291`, `A375-siNC` | "cell line" tissue field; SK-MEL-28 / A375 are melanoma lines (also engineered with siRNA) |
| ☐ | [SRP563552](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP563552) | 4 | `TIGK 60C-2`, `TIGK -25C-3`, `TIGK 37C-3`, `TIGK 37C-2` | TIGK (temperature-stress series; HPV16 depth 200–500) |

### 4b. Flagged as `engineered` (7 samples in 3 studies)

| ☐ | Study | n | Title pattern | Matched on |
|---|---|---|---|---|
| ☐ | [SRP514210](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP514210) | 2 | `LV-STK38-3`, `LV-NC3` | `LV-` lentiviral construct prefix |
| ☐ | [SRP545498](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP545498) | 2 | `IR_siSIRT7`, `IR_siControl` | `si` knockdown prefix |
| ☐ | [SRP552523](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP552523) | 3 | (no titles); tissue text `skin (shTFPI2-2-2)` etc. | `sh` knockdown prefix |

## Summary of expected reclassifications

If priorities 1–2 confirm as suspected:

| | Current | After curation (estimated) |
|---|---|---|
| cell_line   | 10 |  10 |
| engineered  |  7 |  9–11 (+SRP536293, +SRP539392, possibly +ERP149729) |
| clinical    | 14 | 10–12 |
| **Confirmed alpha-type clinical hits** | 4 | **0–2** |

Net effect: the "essentially zero confirmed alpha-type HPV in fresh tissue" claim in the report tightens from "essentially zero" to "zero", and the cohort's typed-clinical signal becomes purely β / γ (10–12 samples across 9 distinct cutaneous references).
