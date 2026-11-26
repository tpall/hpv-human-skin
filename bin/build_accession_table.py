#!/usr/bin/env python3
"""List every library submitted to the screen, with its annotation and its fate.

Supports the Data Availability Statement: the per-cohort accession lists from
which the denominators of the paper can be reconstructed. One row per library
that has a Kraken2 report file, in either arm of the cohort flow (Figure S6):
the untargeted sweep (results_full_v2) and the targeted lesion cohorts
(results_lesions_wide, results_lesions_csc). The three sets are disjoint.

Columns
  cohort              sweep | lesion
  tissue_category     skin | anogenital | oral | other (the pipeline's Estonian
                      values nahk | anogenitaal | suuoos | muu, translated)
  is_cell_line,       as used in the burden analysis (after the override table)
  is_engineered       where the library is in it, otherwise as in the samplesheet
  screened            TRUE if the Kraken2 report is non-empty; FALSE rows are the
                      libraries of supplementary_data_unscreened_libraries.tsv
  in_burden_analysis  TRUE if the library is in the depth-normalized burden table
  sample_group        its group there (empty otherwise)

Usage: python3 bin/build_accession_table.py [outfile]
"""
import csv
import glob
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    ROOT, 'manuscript/supplementary_data_accessions.tsv')
BURDEN = os.path.join(ROOT, 'results_reframe_fixed/hpv_rpm_per_sample.tsv')

RUNS = [('sweep', 'results_full_v2', 'metadata/samplesheet_enriched.csv'),
        ('lesion', 'results_lesions_wide', 'lesions_samplesheet_full.csv'),
        ('lesion', 'results_lesions_csc', 'lesions_samplesheet_full.csv')]
TISSUE = {'nahk': 'skin', 'anogenitaal': 'anogenital', 'suuoos': 'oral', 'muu': 'other'}


def flag(v):
    return 'TRUE' if str(v).strip().lower() == 'true' else 'FALSE'


burden = {r['srr_id']: r for r in csv.DictReader(open(BURDEN), delimiter='\t')}

cols = ['srr_id', 'srx_id', 'cohort', 'study', 'layout', 'tissue_category', 'diagnosis',
        'is_cell_line', 'is_engineered', 'screened', 'in_burden_analysis', 'sample_group']
rows, seen = [], set()
for cohort, run, sheet in RUNS:
    base = os.path.join(ROOT, run)
    meta = {r['srr_id']: r for r in csv.DictReader(open(os.path.join(base, sheet)))}
    for f in sorted(glob.glob(os.path.join(base, 'aggregated/*_kraken2_report.txt'))):
        sid = os.path.basename(f)[:-len('_kraken2_report.txt')]
        assert sid not in seen, sid
        seen.add(sid)
        m, b = meta.get(sid, {}), burden.get(sid)
        src = b if b else m
        rows.append([sid, m.get('srx_id', ''), cohort, m.get('study', ''), m.get('layout', ''),
                     TISSUE.get(m.get('tissue_category', ''), m.get('tissue_category', '')),
                     m.get('diagnosis', ''), flag(src.get('is_cell_line')),
                     flag(src.get('is_engineered')),
                     'TRUE' if os.path.getsize(f) > 0 else 'FALSE',
                     'TRUE' if b else 'FALSE', b['tier'] if b else ''])

with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(cols)
    w.writerows(rows)

for cohort in ('sweep', 'lesion'):
    sub = [r for r in rows if r[2] == cohort]
    print('%s: %d submitted, %d screened, %d in the burden analysis'
          % (cohort, len(sub), sum(r[9] == 'TRUE' for r in sub), sum(r[10] == 'TRUE' for r in sub)))
print(OUT)
