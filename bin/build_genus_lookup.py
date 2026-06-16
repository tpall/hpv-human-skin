#!/usr/bin/env python3
"""
Build an HPV reference → genus lookup from authoritative PaVE taxonomy.

PaVE's /api/genome listing carries a `taxonomy` array per genome, e.g.
    ["HPV4REF", "Gammapapillomavirus 1", "Gammapapillomavirus", "Papillomaviridae"]
which gives the ICTV genus (and species) directly. This is the same source
`bin/build_hpv_refs.sh` uses to build the reference panel, so the genus
assignment stays consistent with — and as authoritative as — the references
themselves, rather than a hand-maintained type-number table that goes stale.

Genus is keyed on `locus_id`, which is exactly the FASTA header / .fai
reference name used downstream (HPV18REF, HPV-m7221nr, ...). The ~2 candidate
"HPV-m...nr" genomes PaVE marks "Unclassified", and any panel reference absent
from PaVE (e.g. RefSeq-merged NC_*) is emitted with genus "unknown" so the
lookup is complete with respect to the actual panel.

Usage:
    # Query PaVE live, emit a row for every reference in the panel .fai:
    build_genus_lookup.py --fai assets/hpv_references/hpv_all.fasta.fai \
                          --out type_genus.csv --save-json pave_genomes.json

    # Offline / firewalled compute: reuse a previously saved listing:
    build_genus_lookup.py --fai hpv_all.fasta.fai --pave-json pave_genomes.json \
                          --out type_genus.csv
"""

import argparse
import csv
import json
import re
import sys
import time
import urllib.request

PAVE_API = "https://pave.niaid.nih.gov/api"

# Full genus name (matches "...papillomavirus") → short label used by the
# downstream analysis (alpha-9 grouping, genus breakdown tables).
GENUS_SHORT = {
    "alphapapillomavirus": "alpha",
    "betapapillomavirus": "beta",
    "gammapapillomavirus": "gamma",
    "mupapillomavirus": "mu",
    "nupapillomavirus": "nu",
}

_GENUS_RE = re.compile(r"^[A-Za-z]+papillomavirus$")
_SPECIES_RE = re.compile(r"^[A-Za-z]+papillomavirus\s+\S+")


def fetch_json(url, retries=3):
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001 — network/parse, retry then raise
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")


def load_listing(args) -> list:
    """Return the PaVE genome listing, from a local JSON file or the live API."""
    if args.pave_json:
        with open(args.pave_json) as f:
            payload = json.load(f)
    else:
        url = f"{args.pave_api.rstrip('/')}/genome?limit=9999&includeNonRef=true"
        print(f"Querying PaVE: {url}", file=sys.stderr)
        payload = fetch_json(url)
    if args.save_json and not args.pave_json:
        with open(args.save_json, "w") as f:
            json.dump(payload, f)
        print(f"Saved listing to {args.save_json}", file=sys.stderr)
    return payload.get("data", [])


def genus_from_taxonomy(taxonomy) -> tuple[str, str, str]:
    """Extract (genus_full, genus_short, species) from a PaVE taxonomy array."""
    if not taxonomy:
        return ("unknown", "unknown", "")
    genus_full = ""
    species = ""
    for token in taxonomy:
        t = str(token).strip()
        if not genus_full and _GENUS_RE.match(t):
            genus_full = t
        elif not species and _SPECIES_RE.match(t):
            species = t
    if not genus_full:
        # PaVE marks candidate types "Unclassified".
        return ("Unclassified", "unclassified", species)
    return (genus_full, GENUS_SHORT.get(genus_full.lower(), genus_full.lower()), species)


def parse_type(ref: str) -> str:
    """Reference name → short type label (HPV18REF→18, HPV-m7221nr→m7221nr)."""
    t = re.sub(r"^HPV[-_]?", "", ref, flags=re.IGNORECASE)
    t = re.sub(r"REF$", "", t, flags=re.IGNORECASE)
    return t


def load_panel_refs(fai_path: str) -> list[str]:
    """First whitespace-delimited column of a .fai (or any ref list)."""
    refs = []
    with open(fai_path) as f:
        for line in f:
            line = line.strip()
            if line:
                refs.append(line.split()[0])
    return refs


def main():
    p = argparse.ArgumentParser(description="Build HPV reference → genus lookup from PaVE")
    p.add_argument("--fai", help="Reference .fai (or ref list); ensures a row per panel reference")
    p.add_argument("--out", "-o", default="type_genus.csv", help="Output CSV")
    p.add_argument("--pave-api", default=PAVE_API, help="PaVE API base URL")
    p.add_argument("--pave-json", help="Read PaVE listing from this JSON file instead of querying")
    p.add_argument("--save-json", help="Save the fetched PaVE listing here for offline reuse")
    args = p.parse_args()

    genomes = [g for g in load_listing(args)
               if (g.get("host_common_name") or "").lower() == "human"]
    print(f"PaVE human genomes: {len(genomes)}", file=sys.stderr)

    # locus_id → (genus_full, genus_short, species)
    by_ref = {}
    for g in genomes:
        ref = g.get("locus_id")
        if ref:
            by_ref[ref] = genus_from_taxonomy(g.get("taxonomy"))

    # Drive the row set from the panel when given, else from PaVE itself.
    if args.fai:
        refs = load_panel_refs(args.fai)
        print(f"Panel references: {len(refs)}", file=sys.stderr)
    else:
        refs = sorted(by_ref)

    n_unknown = 0
    rows = []
    for ref in refs:
        genus_full, genus_short, species = by_ref.get(ref, ("unknown", "unknown", ""))
        if genus_short in ("unknown", "unclassified"):
            n_unknown += 1
        rows.append({
            "hpv_reference": ref,
            "type": parse_type(ref),
            "genus": genus_full,
            "genus_short": genus_short,
            "species": species,
        })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["hpv_reference", "type", "genus",
                                          "genus_short", "species"])
        w.writeheader()
        w.writerows(rows)

    # Report genus tally (the headline split the analysis keys off).
    from collections import Counter
    tally = Counter(r["genus_short"] for r in rows)
    print(f"\nWrote {len(rows)} references to {args.out}", file=sys.stderr)
    for genus, n in tally.most_common():
        print(f"  {genus:13s}: {n}", file=sys.stderr)
    if n_unknown:
        print(f"  ({n_unknown} unknown/unclassified — candidate or non-PaVE refs)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
