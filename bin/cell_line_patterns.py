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
    # Skin / keratinocyte
    r"hacat",
    r"n[/\s-]?tert(?:[\s-]?1)?",
    r"niks",
    r"a431",
    r"tigk",        # telomerase-immortalized gingival keratinocyte
    r"hok[\s-]?16b",
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
    """Compile each pattern with word boundaries; keep label = pattern source."""
    compiled = []
    for p in patterns:
        # Add word boundaries to named lines that don't already have them.
        # Generic patterns supply their own \b where needed.
        if not p.startswith(r"\b") and not p.endswith(r"\b"):
            wrapped = rf"\b{p}\b"
        else:
            wrapped = p
        compiled.append((re.compile(wrapped, re.IGNORECASE), p))
    return compiled


_NAMED = _build_compiled(NAMED_CELL_LINES)
_GENERIC = _build_compiled(GENERIC_PATTERNS)


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
