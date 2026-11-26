#!/usr/bin/env python3
"""Summarise the read geometry that decides which libraries reach the screen.

Supports the Materials and Methods claim that droplet single-cell libraries are
excluded exactly rather than incidentally: fastp runs with --length_required 36
and discards a pair when either mate fails, so a library whose first mate is the
26-33 nt barcode-and-UMI read of droplet chemistry loses every read and returns
no screening result. The claim to support is that the cut is sharp, i.e. every
library in the sweep with a sub-36 nt first mate is unscreened and no screened
library has one.

Writes a small key/value table so the manuscript can quote these numbers without
re-reading 3,600 fastp reports at render time.

Usage: python3 bin/build_read_geometry_summary.py [outfile]
"""
import csv
import glob
import json
import os
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    ROOT, 'manuscript/read_geometry_summary.tsv')
UNSCR = os.path.join(ROOT, 'manuscript/supplementary_data_unscreened_libraries.tsv')
MINLEN = 36

geom = {}
for f in glob.glob(os.path.join(ROOT, 'results_full_v2/chunks/*/qc/fastp/*_fastp.json')):
    sid = os.path.basename(f)[:-len('_fastp.json')]
    try:
        s = json.load(open(f))['summary']['before_filtering']
        geom[sid] = (s.get('read1_mean_length'), s.get('read2_mean_length'))
    except Exception:
        pass

unscr = {r['srr_id']: r for r in csv.DictReader(open(UNSCR), delimiter='\t')}
filtered = [k for k, v in unscr.items() if v['cause'] == 'length_filtered']
screened = [k for k in geom if k not in unscr]

def med(vals):
    vals = [v for v in vals if v]
    return round(statistics.median(vals)) if vals else ''

r1_f = [geom[k][0] for k in filtered if k in geom]
r2_f = [geom[k][1] for k in filtered if k in geom]
r1_s = [geom[k][0] for k in screened if geom[k][0]]

rows = [
    ('droplet_n',             len(filtered)),
    ('droplet_paired',        sum(1 for k in filtered if unscr[k]['layout'] == 'PAIRED')),
    ('droplet_single',        sum(1 for k in filtered if unscr[k]['layout'] == 'SINGLE')),
    ('droplet_r1_med',        med(r1_f)),
    ('droplet_r2_med',        med(r2_f)),
    ('droplet_r2_ge_min',     sum(1 for v in r2_f if v and v >= MINLEN)),
    ('droplet_studies',       len({unscr[k]['study'] for k in filtered})),
    # the sharpness check: these two decide whether the exclusion is exact
    ('sweep_r1_below_min',    sum(1 for k in geom if geom[k][0] and geom[k][0] < MINLEN)),
    ('sweep_r1_below_min_unscreened',
                              sum(1 for k in geom if geom[k][0] and geom[k][0] < MINLEN and k in unscr)),
    ('screened_with_geometry', len(screened)),
    ('screened_r1_below_min', sum(1 for v in r1_s if v < MINLEN)),
    ('screened_r1_min',       round(min(r1_s)) if r1_s else ''),
]
with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['key', 'value'])
    for k, v in rows:
        w.writerow([k, v])
print('wrote %s (%d keys)' % (OUT, len(rows)))
for k, v in rows:
    print('  %-32s %s' % (k, v))
