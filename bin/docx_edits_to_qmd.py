#!/usr/bin/env python3
"""Locate a collaborator's .docx edits in the Quarto source.

Do NOT convert an edited .docx back into .qmd. The rendered document has already
lost everything that makes the source reproducible: inline R expressions have
become literal numbers, `@fig-rpm` has become "Figure 2", and `[@arron2011]` has
become "Arron et al. 2011". Converting back freezes the manuscript.

The workable loop is instead:

    render .docx  ->  collaborator edits (track changes or not)  ->  this script
    ->  apply the reported changes to the .qmd by hand

Given the edited file and the version it was made from, this reports every
changed paragraph, and for each one the .qmd file and line where its source
lives, found by matching the paragraph back through the inline-R placeholders.

Usage:
  bin/docx_edits_to_qmd.py EDITED.docx [--from manuscript/manuscript.docx]
                           [--src manuscript]
"""
import argparse
import difflib
import pathlib
import re
import subprocess
import sys
import tempfile

QUARTO = "quarto"


def docx_text(path):
    """Rendered document as plain-text paragraphs."""
    with tempfile.TemporaryDirectory() as td:
        out = pathlib.Path(td) / "x.txt"
        subprocess.run([QUARTO, "pandoc", str(path), "-t", "plain",
                        "--wrap=none", "-o", str(out)],
                       check=True, capture_output=True)
        txt = out.read_text()
    paras = [re.sub(r"\s+", " ", p).strip() for p in txt.split("\n\n")]
    return [p for p in paras if len(p.split()) >= 6]


def qmd_paragraphs(srcdir):
    """Every source paragraph, as (file, line, raw text)."""
    out = []
    for f in sorted(pathlib.Path(srcdir).glob("_*.qmd")):
        line = 1
        for block in f.read_text().split("\n\n"):
            if block.strip() and not block.lstrip().startswith(("#", "|", ":", "```")):
                out.append((f.name, line, re.sub(r"\s+", " ", block).strip()))
            line += block.count("\n") + 2
    return out


def to_pattern(qmd_para):
    """A regex matching what this source paragraph renders to.

    Inline R, citations and cross-references all become unpredictable text in
    the rendered file, so each is replaced by a wildcard and the literal prose
    between them is what actually anchors the match.
    """
    parts = re.split(r"`r [^`]*`|\[@[^\]]*\]|@[a-z]+-[a-z0-9]+", qmd_para)
    lits = [p.strip() for p in parts if len(p.strip()) > 12]
    if not lits:
        return None
    # Two patterns, tried in order. Requiring every literal run to match makes
    # the search fail on exactly the paragraphs that were edited, since an edit
    # anywhere breaks one of the anchors. The opening run is usually unique on
    # its own, and the longest run is the fallback when the opening was the part
    # that got rewritten.
    head = re.compile(".*?".join(re.escape(x) for x in lits[:2]), re.S)
    longest = re.compile(re.escape(max(lits, key=len)), re.S)
    return (head, longest)


def locate(pats, text):
    """First source paragraph whose opening run, or failing that whose longest
    run, appears in the rendered text."""
    for idx in (0, 1):
        for f, ln, pat in pats:
            if pat and pat[idx].search(text):
                return f"{f}:{ln}"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("edited")
    ap.add_argument("--from", dest="orig", default="manuscript/manuscript.docx",
                    help="the rendered file the edits were made on")
    ap.add_argument("--src", default="manuscript", help="directory of _*.qmd sources")
    args = ap.parse_args()

    old, new = docx_text(args.orig), docx_text(args.edited)
    src = qmd_paragraphs(args.src)
    pats = [(f, ln, to_pattern(t)) for f, ln, t in src]

    changed = 0
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, old, new).get_opcodes():
        if op == "equal":
            continue
        for k in range(j1, j2) if op != "delete" else []:
            changed += 1
            para = new[k]
            # Locate against the ORIGINAL paragraph, not the edited one: the
            # source still says what the unedited render said, so matching the
            # edited text would fail on exactly the sentences that changed.
            anchor = old[i1 + (k - j1)] if op == "replace" and i1 + (k - j1) < i2 else para
            where = locate(pats, anchor) or (locate(pats, para) if anchor is not para else None)
            print(f"\n--- change {changed} [{op}] -> {where or 'SOURCE NOT LOCATED'}")
            if op == "replace" and i1 < i2:
                for d in difflib.unified_diff(old[i1].split(), para.split(),
                                              lineterm="", n=2):
                    if d.startswith(("+", "-")) and not d.startswith(("+++", "---")):
                        print("   ", d)
            else:
                print("   ", para[:300])
        if op == "delete":
            for k in range(i1, i2):
                changed += 1
                print(f"\n--- change {changed} [deleted]")
                print("   ", old[k][:300])

    print(f"\n{changed} changed paragraph(s); "
          f"{len(old)} -> {len(new)} paragraphs overall")
    if not changed:
        print("no textual differences found")


if __name__ == "__main__":
    main()
