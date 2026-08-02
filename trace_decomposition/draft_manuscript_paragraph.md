# Draft manuscript text — early/late decomposition of the sub-threshold trace

*Status: DRAFT for the collab discussion; not yet wired into `_results.qmd`/`_discussion.qmd`
(manuscript rewrite deferred). Numbers are from `trace_composition.csv` /
`decompose_trace.py` and should be rendered inline (`r ...`) when integrated.*

---

## Proposed Results subsection

### The sub-threshold trace is part contamination, part latent early-gene transcription

The trace-level signal in unselected skin and neoplasia raises an obvious question: is
this residual HPV read a genuine, low-level transcript of the commensal reservoir, or a
contaminant/cross-mapping artefact at the sensitivity floor? Because our pipeline already
resolves each HPV-positive library into per-gene read counts, we could decompose the
sub-threshold signal directly, splitting each library's HPV reads into an early
(*E1*/*E2*/*E5*/*E6*/*E7*) and a late (*L1*/*L2*) total and comparing the below-floor
trace tissues against the productive controls. Among the 278 non-cell-line, below-floor
HPV-positive libraries of skin/neoplasia, the signal resolves into three parts. **56%
carry no reads mapping to any HPV gene at all** — a Kraken2 k-mer call that does not
survive alignment and gene assignment, i.e. cross-mapping or low-complexity noise rather
than transcription. A further **28% carry gene-mapped *early* reads with the late (L1/L2)
programme absent** (median late-read fraction 0.000), and only **16% reach the productive
L1 threshold** as low-burden infections. The productive controls (warts/EV, *n* = 54) are
the mirror image: **91% are productive** (L1 ≥ 3) with a median late fraction of 0.37, the
expected differentiation-coupled capsid burst.

Crucially, the early-biased trace is not cell-line contamination in disguise. Integrated
HPV in HeLa and related lines is itself early-biased (*E6*/*E7* expressed, *L1* silenced by
integration), so an early-only signature alone cannot exclude an index-hopped or
reagent-borne contaminant. But contamination of that origin is mucosal HPV16/HPV18,
whereas among the *E6*/*E7*-transcribing trace libraries only 2 of 64 were HPV16/HPV18;
the remainder were cutaneous or untypeable-at-low-coverage types. The early-gene trace
therefore carries the tissue-appropriate cutaneous reservoir, not the HeLa mucosal
signature.

This decomposition puts a mechanism under the "quiet reservoir": the below-floor signal is
a mixture of a dominant technical floor (k-mer-only calls with no transcript structure)
and a real biological minority showing exactly the transcriptional fingerprint of latent
or abortive infection — the early programme without the differentiation-dependent late
burst — on cutaneous types. It neither inflates the reservoir into productive infection nor
dismisses it as pure artefact.

## Proposed Discussion reframe

*Current limitation text (Discussion, limitation 4):* "…the same caveat attaches to
distinguishing low-level true infection from residual contaminant signal at the detection
floor."

*Proposed replacement / addition:* "We could partially resolve this at the transcript
level: decomposing the below-floor libraries into early and late viral reads shows that a
majority (56%) carry no gene-mapped reads and are best read as a k-mer/cross-mapping floor,
while a real minority carry early-gene (*E6*/*E7*) transcripts without the late programme —
the signature of latent or abortive infection — on cutaneous rather than HeLa mucosal
types. The sub-threshold trace is therefore neither wholly artefactual nor productive, but a
technical floor overlying genuine low-level early-gene transcription of the commensal
reservoir, consistent with the twilight-zone model of cutaneous HPV latency [@doorbar2023]."
