# Proposed additions to bin/cell_line_patterns.py NAMED_CELL_LINES, needed ONLY
# if the gut RNA query is broadened beyond the tight anorectal terms (colon/IBD/
# organoid material drags in colorectal carcinoma lines). The tight anorectal
# query in this directory is 93% cell-line-free and does NOT require these.
#
# All are heavily-sequenced, unambiguous colorectal-line names (counts = libraries
# in a broad human gut RNA-seq context, 2026-07). None are currently in the
# patterns. Word-boundary wrapping is handled by _build_compiled(), same as the
# existing entries. Low false-positive risk (distinctive names), but review before
# merging per the project's cautious-curation stance.
COLORECTAL_CELL_LINES = [
    r"hct[\s-]?116",       # 9,920
    r"sw[\s-]?480",        # 2,032
    r"caco[\s-]?2",        # 1,797
    r"ht[\s-]?29",         #   883
    r"dld[\s-]?1",         #   881
    r"sw[\s-]?620",        #   794
    r"\brko\b",            #   750  (bare word; RKO is distinctive)
    r"lovo",               #   379
    r"\bt84\b",            #   273  (bare word to avoid partial hits)
    r"ls[\s-]?174t",       #   175
    r"ls[\s-]?180",
    r"hct[\s-]?15",        #   128
    r"colo[\s-]?205",      #    41
    r"colo[\s-]?320",
    r"hct[\s-]?8",
    r"sw[\s-]?48",
    r"lim[\s-]?1215",
]
# To apply: append COLORECTAL_CELL_LINES to NAMED_CELL_LINES in
# bin/cell_line_patterns.py (or import and extend), then re-run
# bin/flag_cell_lines.py over the gut samplesheet.
