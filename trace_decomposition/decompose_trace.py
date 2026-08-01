import glob, csv, os, re, statistics as st, json
from collections import Counter
base="/gpfs/helios/home/taavi74/Projects/hpv-human-skin"
FLOOR=13.2
# all transcript_classes + hpv_types across every result dir (first path wins)
def index(suffix):
    d={}
    for f in glob.glob(base+"/**/*"+suffix, recursive=True):
        s=os.path.basename(f).replace(suffix,""); d.setdefault(s,f)
    return d
TC=index("_transcript_classes.tsv"); TY=index("_hpv_types.tsv")
def genes(f):
    d={r["gene"]:r["read_count"] for r in csv.DictReader(open(f),delimiter="\t")}
    gi=lambda k:(int(d[k]) if d.get(k,"").lstrip("-").isdigit() else 0)
    return dict(E6=gi("E6"),E7=gi("E7"),L1=gi("L1"),L2=gi("L2"),EARLY=gi("EARLY_TOTAL"),LATE=gi("LATE_TOTAL"))
def domtype(s):
    f=TY.get(s)
    if not f: return None
    best=None;br=-1
    for r in csv.DictReader(open(f),delimiter="\t"):
        rc=float(r["read_count"]) if (r.get("read_count") or "").replace(".","").isdigit() else 0
        if rc>br:
            br=rc; m=re.search(r"HPV[- ]?(\d+[A-Za-z]?)",re.sub("REF$","",r["hpv_reference"] or "")); best=m.group(1).upper() if m else "cand"
    return best
meta={r["srr_id"]:dict(tier=r["tier"],cl=r["is_cell_line"].upper()=="TRUE",rpm=float(r["rpm"] or 0))
      for r in csv.DictReader(open("/tmp/rpm.tsv"),delimiter="\t")}
rows=[]
for s,f in TC.items():
    m=meta.get(s)
    if not m or m["cl"]: continue
    g=genes(f); tot=g["EARLY"]+g["LATE"]
    rows.append(dict(tier=m["tier"],rpm=m["rpm"],**g,e67=g["E6"]+g["E7"],
                     latefrac=(g["LATE"]/tot if tot else 0.0),type=domtype(s)))
def grp(name,pred):
    g=[r for r in rows if pred(r)]
    if not g: print(f"\n{name}: n=0"); return
    n=len(g); trans=[r for r in g if r["e67"]>0]
    hela=sum(r["type"] in("16","18") for r in trans)
    print(f"\n{name}: n={n}")
    print(f"  any early-gene reads: {sum(r['EARLY']>0 for r in g)}/{n} ({100*sum(r['EARLY']>0 for r in g)/n:.0f}%)")
    print(f"  E6/E7 > 0:            {len(trans)}/{n} ({100*len(trans)/n:.0f}%)")
    print(f"  L1>=3 (productive):   {sum(r['L1']>=3 for r in g)}/{n} ({100*sum(r['L1']>=3 for r in g)/n:.0f}%)")
    print(f"  median late fraction: {st.median(r['latefrac'] for r in g):.3f}  (E6/E7+ subset: {st.median(r['latefrac'] for r in trans) if trans else float('nan'):.3f})")
    if trans: print(f"  E6/E7+ dom types: HPV16/18={hela}/{len(trans)}; {dict(Counter(r['type'] for r in trans).most_common(6))}")
print("=== EARLY-vs-LATE DECOMPOSITION | floor=%.1f RPM | non-cell-line ==="%FLOOR)
print("transcript_classes total:",len(TC),"| joined non-cell-line:",len(rows))
grp("SUB-THRESHOLD skin (unselected, rpm<floor)",lambda r:r["tier"]=="unselected_skin" and r["rpm"]<FLOOR)
grp("SUB-THRESHOLD all trace (skin+neo+other)",lambda r:r["tier"] in("unselected_skin","neoplasia","other") and r["rpm"]<FLOOR)
grp("PRODUCTIVE CONTROLS (warts+EV)",lambda r:r["tier"] in("control_productive","EV"))
# export rows for local figure
json.dump(rows, open("/tmp/trace_rows.json","w"))
print("\nwrote /tmp/trace_rows.json ("+str(len(rows))+" rows)")
