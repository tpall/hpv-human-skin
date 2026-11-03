#!/usr/bin/env python3
"""List the sweep libraries that returned no screening result, with the cause.

Supports the Materials and Methods statement that the positivity denominator is
the screened libraries, not every library submitted to the screen. A library
whose Kraken2 report file is empty produced no result; hpv_status.tsv records it
as HPV-negative with zero reads, which would count it as a negative observation.

Two causes, distinguished by whether the library has a fastp report:
  length_filtered  downloaded, but fastp kept 0 reads (R1 under --length_required 36)
  not_downloaded   no fastp output at all; the fail-soft download policy skipped it

Usage: python3 bin/build_unscreened_table.py [outfile]
"""
import csv
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AGG = os.path.join(ROOT, 'results_full_v2/aggregated')
SHEET = os.path.join(ROOT, 'results_full_v2/metadata/samplesheet_enriched.csv')
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    ROOT, 'manuscript/supplementary_data_unscreened_libraries.tsv')

empty = sorted(os.path.basename(f)[:-len('_kraken2_report.txt')]
               for f in glob.glob(os.path.join(AGG, '*_kraken2_report.txt'))
               if os.path.getsize(f) == 0)

fastp = {}
for f in glob.glob(os.path.join(ROOT, 'results_full_v2/chunks/*/qc/fastp/*_fastp.json')):
    sid = os.path.basename(f)[:-len('_fastp.json')]
    try:
        s = json.load(open(f))['summary']
        fastp[sid] = (s['before_filtering']['total_reads'],
                      s['before_filtering'].get('read1_mean_length'),
                      s['before_filtering'].get('read2_mean_length'),
                      s['after_filtering']['total_reads'])
    except Exception:
        fastp[sid] = None

chunk = {}
for f in glob.glob(os.path.join(ROOT, 'results_full_v2/chunks/*/samplesheet.csv')):
    c = os.path.basename(os.path.dirname(f))
    for row in csv.DictReader(open(f)):
        chunk[row['srr_id']] = c

meta = {r['srr_id']: r for r in csv.DictReader(open(SHEET))}

cols = ['srr_id', 'study', 'tissue_category', 'diagnosis', 'layout', 'batch', 'cause',
        'reads_submitted', 'read1_mean_length', 'read2_mean_length', 'reads_after_filtering']
with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(cols)
    for sid in empty:
        m = meta.get(sid, {})
        fp = fastp.get(sid)
        cause = 'length_filtered' if sid in fastp else 'not_downloaded'
        w.writerow([sid, m.get('study', ''), m.get('tissue_category', ''),
                    m.get('diagnosis', ''), m.get('layout', ''), chunk.get(sid, ''), cause]
                   + ([fp[0], fp[1], fp[2], fp[3]] if fp else ['', '', '', '']))

n_len = sum(1 for s in empty if s in fastp)
print('%s: %d libraries (%d length_filtered, %d not_downloaded)'
      % (OUT, len(empty), n_len, len(empty) - n_len))
