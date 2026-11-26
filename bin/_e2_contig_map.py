#!/usr/bin/env python3
"""contig -> body-site / disease-state map for the E2 integrity analysis.

Two selectors, because SRA labels the two arms differently:

  * commensal carriage is registered under a METAGENOME taxid
    ("human skin metagenome", "vaginal metagenome", ...), so it is selected by
    ScientificName;
  * carcinoma and dysplasia are registered under Homo sapiens, so they are
    invisible to that filter and must be selected by BioProject. Cohort
    identities were confirmed from BioSample attributes / study titles, quoted
    below; PRJNA664651 carries `tissue = cervical squamous cell carcinoma`.

Anogenital warts are kept separate: benign productive lesions, expected
episomal, so they act as a disease-state control rather than a cancer arm.
Cancer/dysplasia assignment wins over the metagenome one if a run has both.
"""
import csv, os, sys

d = os.path.join(os.environ["REPO"], "dna_rna_triangulation")

SITE = {'skin':    {'human skin metagenome', 'skin metagenome'},
        'vaginal': {'vaginal metagenome', 'human vaginal metagenome'},
        'gut':     {'human gut metagenome', 'gut metagenome'}}

PROJECT = {
    'cervical_cancer': {
        'PRJNA664651',   # "HPV"; BioSample tissue = cervical squamous cell carcinoma
        'PRJNA305342'},  # "Genome sequencing of cervical cancer in china"
    'cervical_cin': {
        'PRJNA687153',   # "CapNGS of human: cervical intraepithelial neoplasia"
        'PRJNA318454'},  # CIN cohort (sample names CIN1-16 / Control-NN)
    'anogenital_wart': {
        'PRJNA517793'},  # "virome of anogenital warts" -- benign, expected episomal
}

hits = list(csv.DictReader(open(f"{d}/cache/all_hits_info.list")))
runs = {t: {r['Run'] for r in hits if (r.get('ScientificName') or '').lower() in s}
        for t, s in SITE.items()}
# the non-human-skin exclusions (transit surfaces, bovine) apply to skin only
runs['skin'] -= {r['run'] for r in
                 csv.DictReader(open(f"{d}/skin_run_exclusions.tsv"), delimiter='\t')}
for cat, projs in PROJECT.items():
    runs[cat] = {r['Run'] for r in hits if r.get('BioProject') in projs}
# disease-state categories take precedence over the metagenome ones
claimed = set().union(*(runs[c] for c in PROJECT))
for t in SITE:
    runs[t] -= claimed

master = list(csv.DictReader(open(f"{d}/cache/2026.04.06.master_table_all_pvs.tsv"),
                             delimiter='\t'))
run2cat = {r: c for c, rs in runs.items() for r in rs}
seen, n = set(), {c: 0 for c in runs}
with open('contig_tissue.tsv', 'w') as fh:
    fh.write('contig\ttissue\n')
    for r in master:
        cat = run2cat.get(r['library'])
        if cat and (r['sequence'], cat) not in seen:
            seen.add((r['sequence'], cat)); n[cat] += 1
            fh.write(f"{r['sequence']}\t{cat}\n")
print("  contigs per category: " + ", ".join(f"{k} {v}" for k, v in n.items() if v),
      file=sys.stderr)
