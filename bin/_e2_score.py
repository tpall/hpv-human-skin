#!/usr/bin/env python3
"""Score E2 presence per contig, conditional on the contig reaching E1."""
import csv, collections, os
COV = 0.5                      # >=50% of the reference ORF covered = present
tis = {r['contig']: r['tissue']
       for r in csv.DictReader(open('contig_tissue.tsv'), delimiter='\t')}
import re
strip = lambda c: re.sub(r':\d+-\d+$', '', c)
tis = {strip(k): v for k, v in tis.items()}
seqs, h = {}, None
for l in open('contigs.fa'):
    if l.startswith('>'):
        h = l.strip().lstrip('>'); seqs[h] = 0
    elif h:
        seqs[h] += len(l.strip())
best = collections.defaultdict(dict)
for r in csv.reader(open('orf_hits.tsv'), delimiter='\t'):
    g = r[0].split('|')[1]; cov = int(r[3]) / int(r[4])
    if cov > best[r[1]].get(g, 0):
        best[r[1]][g] = cov
has = lambda c, g: best[c].get(g, 0) >= COV
rows = []
for site in ('skin', 'vaginal', 'gut', 'anogenital_wart', 'cervical_cin', 'cervical_cancer'):
    cs = [c for c in seqs if tis.get(c) == site]
    for lab, sel in (('all', lambda c: True), ('ge7kb', lambda c: seqs[c] >= 7000)):
        sub = [c for c in cs if sel(c)]
        e1 = [c for c in sub if has(c, 'E1')]
        both = sum(1 for c in e1 if has(c, 'E2'))
        rows.append(dict(site=site, contig_set=lab, n=len(sub), n_E1=len(e1),
                         n_E2=sum(1 for c in sub if has(c, 'E2')), n_E2_given_E1=both,
                         pct_E2_given_E1=round(100 * both / len(e1), 1) if e1 else ''))
out = os.path.join(os.environ["REPO"], "dna_rna_triangulation", "e2_integrity_by_site.tsv")
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(open(out).read())
