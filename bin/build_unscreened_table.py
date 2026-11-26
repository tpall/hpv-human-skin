#!/usr/bin/env python3
"""List the libraries that returned no screening result, with the cause.

Supports the Materials and Methods statement that the positivity denominator is
the screened libraries, not every library submitted to the screen. A library
whose Kraken2 report file is empty produced no result; hpv_status.tsv records it
as HPV-negative with zero reads, which would count it as a negative observation.

Covers both arms of the cohort flow (Figure S6): the untargeted sweep
(results_full_v2) and the targeted lesion cohorts (results_lesions_wide,
results_lesions_csc). The `cohort` column says which.

Two causes, distinguished by whether the library has a fastp report:
  length_filtered  downloaded, but fastp kept 0 reads (R1 under --length_required 36)
  not_downloaded   no fastp output at all; the fail-soft download policy skipped it

The fastp reports of the lesion runs live on the cluster; those of the
unscreened libraries must be synced into results_lesions_*/chunks/ first, or
every lesion row is misread as not_downloaded.

Usage: python3 bin/build_unscreened_table.py [outfile]
"""
import csv
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    ROOT, 'manuscript/supplementary_data_unscreened_libraries.tsv')

# (cohort label, run directory, samplesheet relative to it)
RUNS = [('sweep', 'results_full_v2', 'metadata/samplesheet_enriched.csv'),
        ('lesion', 'results_lesions_wide', 'lesions_samplesheet_full.csv'),
        ('lesion', 'results_lesions_csc', 'lesions_samplesheet_full.csv')]

# The pipeline's tissue categories are Estonian strings; publish them in English.
TISSUE = {'nahk': 'skin', 'anogenitaal': 'anogenital', 'suuoos': 'oral', 'muu': 'other'}

cols = ['srr_id', 'cohort', 'study', 'tissue_category', 'diagnosis', 'layout', 'batch',
        'cause', 'reads_submitted', 'read1_mean_length', 'read2_mean_length',
        'reads_after_filtering']
rows = []
for cohort, run, sheet in RUNS:
    base = os.path.join(ROOT, run)
    empty = sorted(os.path.basename(f)[:-len('_kraken2_report.txt')]
                   for f in glob.glob(os.path.join(base, 'aggregated/*_kraken2_report.txt'))
                   if os.path.getsize(f) == 0)
    fastp, chunk = {}, {}
    for f in glob.glob(os.path.join(base, 'chunks/*/qc/fastp/*_fastp.json')):
        sid = os.path.basename(f)[:-len('_fastp.json')]
        chunk[sid] = f.split(os.sep)[-4]
        try:
            s = json.load(open(f))['summary']
            fastp[sid] = (s['before_filtering']['total_reads'],
                          s['before_filtering'].get('read1_mean_length'),
                          s['before_filtering'].get('read2_mean_length'),
                          s['after_filtering']['total_reads'])
        except Exception:
            fastp[sid] = None
    for f in glob.glob(os.path.join(base, 'chunks/*/samplesheet.csv')):
        c = os.path.basename(os.path.dirname(f))
        for row in csv.DictReader(open(f)):
            chunk[row['srr_id']] = c
    meta = {r['srr_id']: r for r in csv.DictReader(open(os.path.join(base, sheet)))}
    for sid in empty:
        m = meta.get(sid, {})
        fp = fastp.get(sid)
        cause = 'length_filtered' if sid in fastp else 'not_downloaded'
        rows.append([sid, cohort, m.get('study', ''),
                     TISSUE.get(m.get('tissue_category', ''), m.get('tissue_category', '')),
                     m.get('diagnosis', ''), m.get('layout', ''), chunk.get(sid, ''), cause]
                    + ([fp[0], fp[1], fp[2], fp[3]] if fp else ['', '', '', '']))

with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(cols)
    w.writerows(rows)

for cohort in ('sweep', 'lesion'):
    sub = [r for r in rows if r[1] == cohort]
    n_len = sum(1 for r in sub if r[7] == 'length_filtered')
    print('%s: %d libraries (%d length_filtered, %d not_downloaded)'
          % (cohort, len(sub), n_len, len(sub) - n_len))
print(OUT)
