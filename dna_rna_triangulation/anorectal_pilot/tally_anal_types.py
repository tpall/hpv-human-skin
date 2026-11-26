import glob, csv, collections, os, re
base="/gpfs/helios/home/taavi74/Projects/hpv-human-skin/results_anorectal_pilot"
# full type tally across all HPV+ samples (dominant type per sample by reads)
files=glob.glob(base+"/aggregated/*_hpv_types.tsv")
dom=collections.Counter()
for f in files:
    best=None; bestr=-1
    for r in csv.DictReader(open(f), delimiter="\t"):
        ref=re.sub("REF$","",(r.get("hpv_reference") or "")); 
        try: reads=float(r.get("read_count") or 0); br=float(r.get("coverage_breadth") or 0)
        except: reads=0; br=0
        if br>=0.10 and reads>=50 and reads>bestr: best=ref; bestr=reads
    if best: dom[best]+=1
print(f"HPV+ samples (real, dominant type): {sum(dom.values())}")
print("-- dominant type across HPV+ anal samples --")
for t,n in dom.most_common(15): print(f"  {n:4d}  {t}")
# productive infection (table4: L1/L2 late transcripts)
t4=base+"/report/report/summary_tables/table4_productive_infection.tsv"
rows=list(csv.DictReader(open(t4), delimiter="\t"))
prodc=collections.Counter(re.sub("REF$","",r["top_hpv_type"]) for r in rows if r.get("top_hpv_type"))
print(f"\nPRODUCTIVE-infection samples (L1/L2 late transcripts): {len(rows)}")
print("-- productive by top type --")
for t,n in prodc.most_common(15): print(f"  {n:4d}  {t}")
