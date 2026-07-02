# Framing angles — how title/abstract/emphasis shift by target journal

The same dataset and results support three framings. The science is identical; what
changes is **what leads, what the "so what" is, and what gets foregrounded vs
backgrounded**. Pick one before finalising the abstract — the current draft is closest
to Angle A (tumor-virology).

| Element | A. Tumor-virology | B. Dermatology / skin cancer | C. Data-reuse / methods |
|---|---|---|---|
| Target journals | J Virol, Tumour Virus Research, PLoS Pathogens | JID, Br J Dermatol | GigaScience, mBio, Genome Medicine |
| Core thesis | Beta-HPV is transcriptionally trace, not productive, in cSCC — a population-scale test of hit-and-run | HPV does not drive established cSCC; it is a quiet bystander, not a therapeutic/screening target | Public RNA-seq can be re-mined for the cutaneous virome, but only with contamination + depth controls |
| Hero result | The depth-controlled neoplasia-vs-control contrast | The cSCC/AK trace finding + clinical interpretation | The contamination partition + reproducible pipeline |
| Table 1 role | Supporting (contamination context) | Supporting | **Lead/centrepiece** |
| Figure 1 role | **Lead** | **Lead** | Co-lead with Table 1 |
| Backgrounded | pipeline mechanics | virology mechanism detail | the hit-and-run debate (becomes one application) |

---

## Angle A — Tumor-virology (current draft)

**Thesis:** the productive/late-transcript phenotype settles where beta-HPV is and isn't active across the skin disease spectrum; the trace-in-neoplasia result is the population-scale transcriptomic test the hit-and-run model has lacked.

**Title (current):** *Cutaneous HPV is transcriptionally quiet in skin and carcinoma but productive in warts: a depth-controlled re-mining of public RNA-seq*

**Abstract opens on:** HPV biology — genera, commensal beta/gamma reservoir, the contested hit-and-run vs protective-immunity debate.
**Abstract closes on:** "the phenomenology both the hit-and-run and protective-immunity models must accommodate."
**Intro emphasis:** life cycle, productive vs latent, the beta-HPV/cSCC mechanism literature.
**Discussion emphasis:** Viarisio/Hasche/Arron; what the transcriptomic phenotype does and does not adjudicate.
→ Already written this way; minimal change needed.

---

## Angle B — Dermatology / skin cancer

**Thesis:** for the clinician/skin-cancer biologist, the actionable message is that HPV is **not a driver of established cSCC** and not a useful therapeutic or screening target there — while warts/EV are genuinely productive and a small subset of lesions carry real infection.

**Title options:**
- *Human papillomavirus is a transcriptionally silent bystander in cutaneous squamous cell carcinoma: evidence from depth-controlled re-mining of 1,000+ public transcriptomes*
- *No productive HPV infection in cutaneous squamous cell carcinoma at population scale*

**Abstract opens on:** the clinical problem — beta-HPV's contested role in keratinocyte carcinoma, immunosuppressed-patient risk, and why it matters whether HPV is a target.
**Abstract closes on:** the clinical implication — HPV-directed prevention/therapy is unlikely to benefit established cSCC; the virus's relevance, if any, is at the precursor/field stage.
**Intro emphasis:** cSCC epidemiology, organ-transplant risk, AK→cSCC progression, the screening/therapy stakes; compress the molecular life-cycle detail.
**Discussion emphasis:** lead with extending **Arron 2011 (JID)** at scale; clinical interpretation; what it means for HPV vaccination/antivirals in skin cancer; keep limitations about precursor lesions prominent.
**Foreground:** the cSCC/AK result and the genuine-productive-lesion exceptions (Table S2). **Background:** pipeline internals, contamination methodology (move more to supplement).

---

## Angle C — Data-reuse / methods

**Thesis:** the petabase of public RNA-seq is a viable but **booby-trapped** substrate for virome discovery; we provide a contamination-aware, depth-normalised, reproducible framework, and HPV-in-skin is the demonstrating case study.

**Title options:**
- *A contamination-aware, depth-normalised framework for mining the cutaneous virome from public RNA-seq*
- *Re-mining public transcriptomes for papillomaviruses: cell-line contamination and sequencing depth are the dominant confounders*

**Abstract opens on:** the opportunity and trap — non-human reads are discarded; archives are huge but contaminated (HeLa/HPV18); naïve mining overstates prevalence.
**Abstract closes on:** the generalisable framework (flagging + breadth/depth thresholds + RPM depth-normalisation) and reusability beyond HPV; HPV-in-skin as proof-of-principle.
**Intro emphasis:** the discarded-reads opportunity (Edgar/Serratus, VIRTUS), the contamination literature, and the methodological gap; compress the beta-HPV oncogenesis debate to one paragraph.
**Discussion emphasis:** lead with the ~7-fold cell-line enrichment and the necessary-minimum control set; position the trace-in-neoplasia and warts/EV findings as **what the method reveals once confounders are removed**.
**Foreground:** Table 1 (contamination) and the pipeline/Code-availability; promote the depth-matched detection (current Figure S1) into the main figures. **Background:** the hit-and-run mechanism becomes one downstream application.

---

## Practical notes
- Switching angle = re-order the abstract + swap the lead Discussion paragraph + move 1–2 display items between main/supplement. The Methods, Results numbers, Figure 1, and `.bib` are unchanged.
- All three benefit from completing full_v2 (removes the partial-denominator objection) and from stating the contamination result as a contribution.
- Quarto makes the swap cheap: edit `_abstract.qmd` ordering, retitle in `manuscript.qmd` YAML, and (for C) move the `@suppfig-depth` float from `_supplement.qmd` into `_results.qmd`.
