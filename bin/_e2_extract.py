#!/usr/bin/env python3
"""Filter the streamed Logan FASTA down to the contigs we need."""
import csv, re, sys
strip = lambda c: re.sub(r':\d+-\d+$', '', c)
want = {strip(r['contig']) for r in csv.DictReader(open('contig_tissue.tsv'), delimiter='\t')}
out, n, keep = open('contigs.fa', 'w'), 0, False
for line in sys.stdin:
    if line.startswith('>'):
        keep = line.lstrip('>').strip() in want
        if keep:
            n += 1; out.write(line)
    elif keep:
        out.write(line)
out.close()
print(f"  extracted {n}/{len(want)} contigs", file=sys.stderr)
