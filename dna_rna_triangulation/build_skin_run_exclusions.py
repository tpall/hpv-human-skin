#!/usr/bin/env python3
"""Flag Shen skin-metagenome runs that are not human skin.

The Shen/Logan skin reservoir is selected by SRA ScientificName
("human skin metagenome" / "skin metagenome").  That taxid is also used by
studies that swab *environmental* surfaces (public-transit handrails) and
non-human animal skin (dairy-farm udder gauze).  Those libraries carry shed
human corneocytes, so they are real deposition -- but they are nobody's skin
and must not sit inside a "skin reservoir" prevalence.

Fetches BioSample attributes for every skin run, caches them, and writes an
auditable exclusion table consumed by shen_skin_dna_vs_rna.R and
shen_reservoir_by_tissue.R.

Outputs (both in this directory):
  cache/skin_biosample_attrs.tsv  raw attributes, one row per BioSample
  skin_run_exclusions.tsv         run, biosample, bioproject, reason
"""
import csv, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
HITS = os.path.join(CACHE, "all_hits_info.list")
ATTRS = os.path.join(CACHE, "skin_biosample_attrs.tsv")
OUT = os.path.join(HERE, "skin_run_exclusions.tsv")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SKIN_SCI = {"human skin metagenome", "skin metagenome"}

# fields that can reveal what was actually sampled
FIELDS = ["env_broad_scale", "env_local_scale", "env_medium", "isolation_source",
          "samp_collect_device", "host", "geo_loc_name", "body_site", "host_body_site"]

# non-human-skin signatures, checked against the joined FIELDS text
EXCLUDE = [
    ("environmental surface (transit)",
     r"transportation|subway|metro\b|\bcity\b|handrail|train\b|\bbus\b|public transport"),
    ("non-human animal",
     r"udder|dairy|\bfarm\b|cattle|bovine|\bcow\b|porcine|\bswine\b|canine|murine|\bmouse\b"),
]


def fetch_attrs(biosamples):
    """BioSample accession -> {attribute: value}, cached on disk."""
    if os.path.exists(ATTRS):
        out = {}
        with open(ATTRS) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                out[row["biosample"]] = {k: row.get(k, "") for k in FIELDS}
        missing = [b for b in biosamples if b not in out]
        if not missing:
            return out
        print(f"cache has {len(out)}, fetching {len(missing)} more", file=sys.stderr)
    else:
        out, missing = {}, list(biosamples)

    for i in range(0, len(missing), 150):
        chunk = missing[i:i + 150]
        data = urllib.parse.urlencode(
            {"db": "biosample", "id": ",".join(chunk),
             "rettype": "full", "retmode": "xml"}).encode()
        xml = ""
        for _ in range(4):
            try:
                xml = urllib.request.urlopen(
                    urllib.request.Request(EUTILS + "efetch.fcgi", data=data),
                    timeout=120).read().decode("utf-8", "replace")
                break
            except Exception:
                time.sleep(4)
        for m in re.finditer(r"<BioSample\b.*?</BioSample>", xml, re.S):
            try:
                node = ET.fromstring(m.group(0))
            except ET.ParseError:
                continue
            d = {}
            for a in node.findall(".//Attribute"):
                key = (a.get("harmonized_name") or a.get("attribute_name") or "").lower()
                if key in FIELDS and a.text:
                    d[key] = re.sub(r"\s+", " ", a.text.strip())
            out[node.get("accession")] = {k: d.get(k, "") for k in FIELDS}
        print(f"  fetched {min(i + 150, len(missing))}/{len(missing)}", file=sys.stderr)
        time.sleep(0.5)

    os.makedirs(CACHE, exist_ok=True)
    with open(ATTRS, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["biosample"] + FIELDS)
        for b in sorted(out):
            w.writerow([b] + [out[b].get(k, "") for k in FIELDS])
    return out


def main():
    if not os.path.exists(HITS):
        sys.exit(f"missing {HITS} -- run shen_skin_dna_vs_rna.R once to populate the cache")

    runs = defaultdict(list)          # biosample -> [(run, bioproject), ...]
    with open(HITS) as fh:
        for r in csv.DictReader(fh):
            if (r.get("ScientificName") or "").lower() in SKIN_SCI and r.get("BioSample"):
                runs[r["BioSample"]].append((r["Run"], r.get("BioProject", "")))

    attrs = fetch_attrs(sorted(runs))

    rows = []
    for bs, members in sorted(runs.items()):
        text = " ".join(attrs.get(bs, {}).get(k, "") for k in FIELDS).lower()
        for reason, pat in EXCLUDE:
            if re.search(pat, text):
                for run, proj in members:
                    rows.append((run, bs, proj, reason))
                break

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["run", "biosample", "bioproject", "reason"])
        w.writerows(sorted(set(rows)))

    n_run = len({r[0] for r in rows})
    print(f"skin BioSamples: {len(runs)}  attributes resolved: {len(attrs)}")
    print(f"excluded runs: {n_run}  BioSamples: {len({r[1] for r in rows})}")
    by = defaultdict(set)
    for run, bs, proj, reason in rows:
        by[reason].add(run)
    for reason, s in sorted(by.items()):
        print(f"  {reason:<32} {len(s):>4} runs")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
