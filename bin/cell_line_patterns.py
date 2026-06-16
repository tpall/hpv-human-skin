"""
Cell-line detection patterns for HPV-relevant SRA metadata.

Single source of truth — imported by both the post-hoc flagging script
(``flag_cell_lines.py``) and the pipeline metadata enrichment
(``parse_metadata.py``) so detection logic stays consistent.

Usage:
    from cell_line_patterns import classify_cell_line
    is_line, matched = classify_cell_line("HeLa cells passage 12")
"""

import re

# Named cell lines relevant to HPV / skin / anogenital / head-and-neck research.
# Patterns are regex with word boundaries; matched case-insensitively.
# Ordering doesn't matter — first hit wins for the returned label.
NAMED_CELL_LINES = [
    # Cervical
    r"hela(?:[\s-]?(?:s3|229))?",
    r"siha",
    r"caski",
    r"me[\s-]?180",
    r"ms[\s-]?751",
    r"c[\s-]?33[\s-]?a",
    r"c[\s-]?4[\s-]?i{1,2}",
    r"ht[\s-]?3",
    r"w12",
    r"cin[\s-]?612",
    r"sw[\s-]?756",
    r"hcc[\s-]?94",
    # Head & neck SCC
    r"upci[\s-]?scc[\s-]?(?:090|152|154)",
    r"um[\s-]?scc[\s-]?(?:47|104|1|22a|22b)",
    r"ud[\s-]?scc[\s-]?2",
    r"93[\s-]?vu[\s-]?147t",
    r"fadu",
    r"detroit[\s-]?562",
    r"scc[\s-]?(?:4|9|13|25)\b",
    r"cal[\s-]?(?:27|33)",
    # Skin / keratinocyte / melanoma
    r"hacat",
    r"n[/\s-]?tert(?:[\s-]?1)?",
    r"niks",
    r"a431",
    r"tigk",        # telomerase-immortalized gingival keratinocyte
    r"hok[\s-]?16b",
    r"mnt[\s-]?1",  # melanoma line, used in melanocyte / pigmentation studies
    r"a375",        # melanoma line (often mislabeled "skin" in public RNA-seq)
    # Common producer / engineered backgrounds (often used as HPV VLP / model hosts)
    r"hek[\s-]?293(?:t|ft)?",
    r"293t",
    r"cos[\s-]?7",
]

# Generic indicators. Less specific than named lines but strong signal in
# combination. Matched case-insensitively; word boundaries used where needed.
GENERIC_PATTERNS = [
    r"\bcell[\s-]?lines?\b",
    r"\bimmortali[sz]ed\b",
    r"\bpassage[\s-]?\d+",
    r"\batcc\b",
    r"\bdsmz\b",
    r"\bjcrb\b",
    r"\becacc\b",
    r"\briken\b",
    r"\bipsc\b",
    r"\binduced pluripotent\b",
    r"\bxenograft\b",
    r"\bpdx\b",
    r"\bstably[\s-]?transfected\b",
    r"\bstable[\s-]?transfectant\b",
]


def _build_compiled(patterns: list[str]) -> list[tuple[re.Pattern, str]]:
    """Compile each pattern with alphanumeric-only boundaries.

    Plain ``\\b`` treats ``_`` as a word char, so ``\\bmnt-?1\\b`` won't match
    inside ``MNT-1_WT1``. We strip any author-supplied outer ``\\b`` and wrap
    with ``(?<![A-Za-z0-9])...(?![A-Za-z0-9])`` so the pattern terminates on
    underscores, hyphens, whitespace, punctuation, and string boundaries alike.
    """
    compiled = []
    for p in patterns:
        core = p
        if core.startswith(r"\b"):
            core = core[2:]
        if core.endswith(r"\b"):
            core = core[:-2]
        wrapped = rf"(?<![A-Za-z0-9]){core}(?![A-Za-z0-9])"
        compiled.append((re.compile(wrapped, re.IGNORECASE), p))
    return compiled


# Engineered / experimentally-manipulated cells: transduction, transfection,
# CRISPR / shRNA / siRNA knockdowns. Not strictly cell lines (the starting
# material may be primary), but also not naive clinical tissue — used to
# stratify in-vitro experiments from clinical signal in HPV reports.
#
# Gene-knockdown labels follow a lowercase-prefix + uppercase-target convention
# (shTFPI2, siSIRT7, shNC, siControl); (?-i:...) disables IGNORECASE locally so
# we don't match "shape", "side", "Skin", etc.
ENGINEERED_PATTERNS = [
    r"\bshRNA\b",
    r"\bsiRNA\b",
    r"\bsgRNA\b",
    r"\bknock[\s-]?(?:down|out)\b",
    r"\bcrispr\b",
    r"\bd?cas9\b",
    r"\btransduc(?:e|ed|ing|tion)\b",
    r"\btransfect(?:ed|ion|ing)\b",
    r"\blentivir(?:al|us)?\b",
    r"\bstably[\s-]?(?:expressing|transfected)\b",
    r"\b(?-i:sh[A-Z]\w*)\b",
    r"\b(?-i:si[A-Z]\w*)\b",
    r"\b(?-i:LV-[A-Z]\w*)\b",
    # Plasmid backbone vectors — strict lowercase p + 3+ uppercase letters,
    # optionally followed by digits or more uppercase chars. Matches pLVX /
    # pLKO / pBABE / pCDH / pCAGGS / pEGFP / pCMV3 etc.
    #
    # The case-sensitive p (via (?-i:)) is what makes this safe in skin /
    # immunology metadata: case-insensitive matching false-positives heavily
    # on uppercase abbreviations like PBMC, PATIENTS, POST, PURPL, and on
    # gene names PARP1 / POLE / PRDM1 / PCDH9 (all P + 3+ uppers). Mangled
    # all-uppercase plasmid names (e.g. PLVX in FB_PLVX_80_Rep1) are missed
    # by this pattern; rely on study-level curation for those.
    #
    # 3+ uppercase letters (not 2+) excludes pDC (plasmacytoid dendritic
    # cells, 25 samples in current cohort) and similar 2-letter cell-type
    # abbreviations. Trade-off: pET / pUC / pBR bacterial cloning vectors
    # would also be missed, but they don't appear in mammalian RNA-seq.
    r"\b(?-i:p[A-Z]{3,}[A-Z0-9]*)\b",
]

# In-vitro cultured material that is NOT an immortalized cell line and NOT
# (necessarily) genetically engineered: primary-cell cultures, organotypic /
# raft / air-liquid-interface skin models, organoids, differentiation
# time-courses. These are routinely mislabeled tissue="skin" yet are not
# clinical tissue — so a separate axis lets the skin analysis require genuine
# biopsies without violating "primary cells aren't cell lines".
#
# HIGH PRECISION by design: patterns must not fire on clinical biopsies.
# Notably we never match a bare "differentiated" (that collides with tumour
# grading: "well/poorly-differentiated carcinoma") — only "differentiated day N"
# and explicit culture vocabulary.
IN_VITRO_PATTERNS = [
    # cultured primary cells (the main gap): "primary keratinocyte",
    # "primary dermal fibroblast", "primary human melanocytes", "primary cells".
    r"primary[\s-]+(?:(?:human|epidermal|dermal|neonatal|adult|oral|gingival|foreskin)[\s-]+)*"
    r"(?:keratinocytes?|melanocytes?|fibroblasts?|cells?|cultures?)",
    # primary-keratinocyte product abbreviations (all primary cultures)
    r"\bn[/\s-]?heka?\b",        # NHEK / NHEKa normal human epidermal keratinocytes
    r"\bhek[\s-]?[an]\b",        # HEKn / HEKa
    # in-vitro skin / epithelial models
    r"\borganotypic\b",
    r"\braft[\s-]?culture",
    r"\bair[\s-]?liquid[\s-]?interface\b",
    r"\b3[\s-]?d[\s-]?culture",
    r"\bspheroids?\b",
    r"\borganoids?\b",
    # explicit culture vocabulary
    r"\bin[\s-]?vitro\b",
    r"\bcell[\s-]?culture\b",
    r"\bmonolayer\b",
    r"\bcultured\b",
    # differentiation/treatment time-course in culture (requires the "day N"
    # so it can't match "well-differentiated <tumour>")
    r"\bdifferentiat\w*[\s-]?day[\s-]?\d+",
    r"\bday[\s-]?\d+[\s-]?(?:post[\s-]?)?(?:differentiat|treatment|induction|culture)",
]

_NAMED = _build_compiled(NAMED_CELL_LINES)
_GENERIC = _build_compiled(GENERIC_PATTERNS)
_ENGINEERED = _build_compiled(ENGINEERED_PATTERNS)
_IN_VITRO = _build_compiled(IN_VITRO_PATTERNS)


def classify_cell_line(text: str) -> tuple[bool, str]:
    """Return (is_cell_line, matched_pattern_label).

    Named lines take precedence over generic indicators so the matched label
    is informative when both fire.
    """
    if not text:
        return (False, "")
    for rx, label in _NAMED:
        if rx.search(text):
            return (True, label)
    for rx, label in _GENERIC:
        if rx.search(text):
            return (True, label)
    return (False, "")


def classify_engineered(text: str) -> tuple[bool, str]:
    """Return (is_engineered, matched_pattern_label) for transduced/edited cells."""
    if not text:
        return (False, "")
    for rx, label in _ENGINEERED:
        if rx.search(text):
            return (True, label)
    return (False, "")


def classify_in_vitro(text: str) -> tuple[bool, str]:
    """Return (is_in_vitro, matched_pattern_label) for cultured/in-vitro material.

    Catches primary-cell cultures and in-vitro skin models that are not
    immortalized lines — distinct from classify_cell_line so the skin analysis
    can exclude cultured material without mislabeling primary cells as lines.
    """
    if not text:
        return (False, "")
    for rx, label in _IN_VITRO:
        if rx.search(text):
            return (True, label)
    return (False, "")


def has_explicit_cell_line_field(value: str) -> bool:
    """True if an SRA cell_line attribute carries a real immortalized-line value.

    SRA samples often have ``cell_line=not applicable`` / ``N/A`` / blank when
    the sample is primary tissue — treat those as negative.

    Researchers also occasionally mis-fill the field with ``primary
    keratinocytes`` or ``primary fibroblasts`` etc. Primary cultures are
    not cell lines for our purposes, so reject any value containing
    ``primary``.
    """
    if not value:
        return False
    v = value.strip().lower()
    if v in {"", "n/a", "na", "none", "not applicable", "not available", "missing"}:
        return False
    if "primary" in v:
        return False
    return True
