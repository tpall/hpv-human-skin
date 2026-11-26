# TCGA-CESC as the documented-tissue cervical comparator

## Why we went looking

The manuscript and the collaborator deck use cervix as the **clean opposite to
skin**: skin carries HPV DNA that is largely transcriptionally silent, cervix
carries HPV that is transcriptionally active. The cervical arm we had was a
public SRA "vaginal baseline", and the E1:E2 integration marker showed it is
**78 % integrated cell-line material**; a random 1 000-library subsample of the
same space ran at **39 %**. The cause is structural, not a filtering mistake:
HeLa is a cervical adenocarcinoma line, so any text search for "cervix" or
"cervical" in SRA is saturated by it. Re-stratifying that space does not fix it.

TCGA-CESC is the opposite situation. Every library is documented primary tissue
with an attached pathology diagnosis, submitted under a protocol, with the
sample type (`Primary Tumor` / `Solid Tissue Normal` / `Metastatic`) recorded as
structured metadata rather than inferred from free text.

## What is here (open metadata, no authorisation needed)

| file | contents |
|---|---|
| `tcga_cesc_rnaseq_manifest.tsv` | all 927 TCGA-CESC RNA-seq BAMs, 7.02 TB |
| `tcga_cesc_genome_bams.tsv` | the 309 we would actually pull, 2.27 TB |

927 is not 927 distinct libraries: GDC publishes each of 309 aliquots under
three workflows (`STAR 2-Pass Genome` / `Transcriptome` / `Chimeric`). Only the
**Genome** BAM retains unaligned reads, so that is the one we want, and the other
two are redundant for us.

Composition of the 309 (304 distinct cases):

| sample type | n |
|---|---|
| Primary Tumor | 304 |
| Solid Tissue Normal | 3 |
| Metastatic | 2 |

Diagnoses are squamous cell carcinoma (744 of the 927 file-rows) and endocervical
adenocarcinoma (138), i.e. the expected CESC mix.

## Feasibility — verified, not assumed

1. **The viral reads are still in the BAMs.** GDC's mRNA pipeline documentation:
   *"This workflow outputs a genomic BAM file, which contains both aligned and
   unaligned reads."* HPV reads cannot align to a human-only reference, so they
   are all in the unaligned fraction.
2. **We do not have to download 2.27 TB.** The GDC slicing API serves an
   arbitrary region, and accepts the keyword `unmapped`:
   `GET https://api.gdc.cancer.gov/slicing/view/{file_id}?region=unmapped`.
   That is exactly and only the fraction we need.
3. **The manifest IDs work.** An unauthenticated slice request against a file ID
   from `tcga_cesc_genome_bams.tsv` returns **HTTP 403**, not 404 — the endpoint
   and the ID are valid; the only missing piece is authorisation.

So the pipeline is verified end to end except for the access step.

## The one blocker: dbGaP

TCGA read-level data is controlled under **phs000178** (The Cancer Genome Atlas,
currently v11.p8). Metadata — everything in this directory — is open; reads are
not. Obtaining them requires:

- an **eRA Commons** account for the PI,
- a Data Access Request submitted through dbGaP Authorized Access,
- a **Data Use Certification** countersigned by the institution's **Signing
  Official** (University of Tartu has one; this is not something a single
  researcher can sign alone),
- approval by the Data Access Committee.

TCGA is broadly consented and approval for a legitimate research use is routine
rather than adversarial, but it is measured in weeks, not minutes. **This is a
human/administrative step — it cannot be automated from here.**

Once approved, download a token from portal.gdc.cancer.gov and run:

```bash
export GDC_TOKEN_FILE=~/.gdc/gdc-user-token.txt
bin/fetch_tcga_cesc_unmapped.sh --limit 5      # sanity check first
bin/fetch_tcga_cesc_unmapped.sh                # all 309
```

Then type the resulting samplesheet with the usual pipeline and score it with

```bash
bin/e1e2_integration_marker.R --agg results_tcga_cesc --no-tiers \
                              --outdir results_tcga_cesc_e1e2
```

running the marker **before** drawing composition conclusions, as with any new
cohort.

## What this fixes, and what it does not

**Fixes:** the *activity* arm. 304 documented primary cervical tumours is a solid
basis for "HPV present and transcriptionally active in cervix", with no
cell-line caveat, and the E1:E2 marker becomes a *finding* there rather than a
contamination check — integration is expected in cervical carcinoma, and the
marker should read positive. Recovering that is a positive control for the
method itself.

**Does not fix:** the *normal cervix* arm. TCGA-CESC has **3** solid-tissue
normals. If the deck's claim needs healthy cervix rather than tumour, TCGA
cannot supply it and neither can public SRA in its current state. GTEx has
ecto- and endocervix but is likewise dbGaP-controlled (phs000424). That arm
needs either a separate access request or a reframing of the claim to
tumour-vs-skin rather than normal-cervix-vs-skin.
