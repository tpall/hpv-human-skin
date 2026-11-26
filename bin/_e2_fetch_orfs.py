#!/usr/bin/env python3
"""Papillomavirus reference ORF proteins (RefSeq) used as tblastn queries."""
import urllib.request, urllib.parse, re, json, collections
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
g = lambda u, d=None: urllib.request.urlopen(
    urllib.request.Request(u, data=d), timeout=180).read().decode('utf8', 'replace')
ids = json.loads(g(E + "esearch.fcgi?db=protein&term=%s&retmax=900&retmode=json" %
    urllib.parse.quote('"Papillomaviridae"[Organism] AND refseq[filter]')
    ))['esearchresult']['idlist']
txt = "".join(g(E + "efetch.fcgi", urllib.parse.urlencode(
    {'db': 'protein', 'id': ','.join(ids[i:i+200]),
     'rettype': 'fasta', 'retmode': 'text'}).encode())
    for i in range(0, len(ids), 200))
GENE = [('E1', r'\bE1\b|replication protein'), ('E2', r'\bE2\b|regulatory protein'),
        ('E6', r'\bE6\b'), ('E7', r'\bE7\b'),
        ('L1', r'\bL1\b|major capsid'), ('L2', r'\bL2\b|minor capsid')]
out, c = [], collections.Counter()
for rec in txt.split('\n>'):
    if not rec.strip():
        continue
    h, _, seq = ('>' + rec.lstrip('>')).partition('\n')
    for gene, pat in GENE:
        if re.search(pat, h, re.I):
            out.append(f">{h.split()[0].lstrip('>')}|{gene}\n{seq.replace(chr(10), '').strip()}")
            c[gene] += 1
            break
open('pv_orf_proteins.fa', 'w').write('\n'.join(out) + '\n')
print("  ORF proteins:", dict(c))
