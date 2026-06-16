#!/usr/bin/env python3
"""
Query NCBI GEO/SRA for RNA-seq datasets relevant to HPV in human tissues.

Searches for public RNA-seq datasets from human skin, anogenital, oral, and
other tissues, then outputs a samplesheet CSV for the Nextflow pipeline.

Usage:
    query_sra.py --query "skin RNA-seq Homo sapiens" --output samplesheet.csv [--max-samples 0]
"""

import argparse
import csv
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import Optional

try:
    from Bio import Entrez
except ImportError:
    sys.exit("Error: Biopython is required. Install with: pip install biopython")

try:
    import pysradb
    from pysradb.sraweb import SRAweb
    HAS_PYSRADB = True
except ImportError:
    HAS_PYSRADB = False

# NCBI requires an email
Entrez.email = "hpv-skin-pipeline@example.com"
Entrez.api_key = None  # Set via NCBI_API_KEY env var for higher rate limits

# Search terms for different tissue categories
TISSUE_QUERIES = {
    "skin": [
        '("skin"[MeSH Terms] OR "epidermis" OR "keratinocyte" OR "cutaneous")',
        '("wart" OR "verruca" OR "papilloma" OR "actinic keratosis")',
        '("basal cell carcinoma" OR "squamous cell carcinoma skin")',
    ],
    "anogenital": [
        '("cervix"[MeSH Terms] OR "cervical" OR "vaginal" OR "vulvar")',
        '("anal" OR "anogenital" OR "penile")',
    ],
    "oral": [
        '("oral mucosa" OR "oropharynx" OR "tonsil" OR "tongue")',
        '("head and neck" OR "larynx" OR "pharynx" OR "buccal")',
    ],
}

BASE_FILTER = (
    '"Homo sapiens"[Organism] AND '
    '"RNA-Seq"[Strategy] AND '
    '"public"[Access]'
)


def search_sra(query: str, max_results: int = 10000) -> list[str]:
    """Search SRA and return list of SRA run accessions."""
    full_query = f"({query}) AND {BASE_FILTER}"
    print(f"  Searching: {full_query[:120]}...", file=sys.stderr)

    handle = Entrez.esearch(db="sra", term=full_query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    uid_list = record.get("IdList", [])
    print(f"  Found {len(uid_list)} UIDs", file=sys.stderr)
    return uid_list


_WS_RE = re.compile(r"\s+")
_CTRL_RE = re.compile(r"[\x00-\x1f\x7f]")


def _clean(s: str) -> str:
    """Flatten whitespace and drop control chars so the value survives any CSV parser.

    SRA TITLE/attribute fields can contain embedded newlines, tabs, and stray
    quote characters; Nextflow's splitCsv tokenizes line-by-line and mis-aligns
    columns when a quoted field spans physical lines.
    """
    if not s:
        return ""
    s = _CTRL_RE.sub(" ", s)
    s = s.replace('"', "'")
    s = _WS_RE.sub(" ", s)
    return s.strip()


def _text(elem, path, default=""):
    """Safely extract text from an XML subelement, or return default."""
    if elem is None:
        return default
    child = elem.find(path)
    if child is None or child.text is None:
        return default
    return _clean(child.text)


def _parse_experiment_package(pkg) -> list[dict]:
    """Extract one sample row per RUN inside an <EXPERIMENT_PACKAGE>."""
    out = []
    exp = pkg.find("EXPERIMENT")
    sample = pkg.find("SAMPLE")
    study = pkg.find("STUDY")
    run_set = pkg.find("RUN_SET")
    if exp is None or run_set is None:
        return out

    srx = exp.get("accession", "")
    study_acc = study.get("accession", "") if study is not None else ""

    # Platform: the <PLATFORM> child has a single tagged child like <ILLUMINA>
    platform = ""
    plat_elem = exp.find("PLATFORM")
    if plat_elem is not None and len(plat_elem) > 0:
        platform = plat_elem[0].tag

    # Library layout: <LIBRARY_LAYOUT> contains either <PAIRED/> or <SINGLE/>
    layout = "SINGLE"
    lib_layout = exp.find("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT")
    if lib_layout is not None and lib_layout.find("PAIRED") is not None:
        layout = "PAIRED"

    # Sample title + tissue-ish attributes. Collect *all* relevant
    # SAMPLE_ATTRIBUTEs (don't break on first match) so a record with both
    # tissue=cervix and cell_type=HeLa surfaces both signals — needed for
    # cell-line vs. primary-tissue classification downstream.
    title = _text(sample, "TITLE")
    source_parts: list[str] = []
    all_attrs: list[str] = []
    cell_line = ""
    genotype = ""
    if sample is not None:
        for attr in sample.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
            tag = (_text(attr, "TAG") or "").lower().replace(" ", "_")
            value = _text(attr, "VALUE")
            if not value:
                continue
            # Keep EVERY attribute in `characteristics`: cell-line / culture
            # hints are routinely misfiled into off-spec tags (e.g.
            # genotype="TOP2A knockdown", note="HeLa-derived"), so downstream
            # detection scans the whole blob rather than a fixed whitelist.
            all_attrs.append(f"{tag}={value}")
            if tag in ("tissue", "cell_type", "source_name",
                       "sample_type", "isolation_source"):
                source_parts.append(f"{tag}={value}")
            if tag == "cell_line" and not cell_line:
                cell_line = value
            # genotype carries the strongest in-vitro-manipulation signals
            # ("TOP2A knockdown", "shNC", CRISPR edits) — kept as its own
            # column too so it stays visible in reports.
            if tag == "genotype" and not genotype:
                genotype = value
    source = " | ".join(source_parts)
    characteristics = " | ".join(all_attrs)

    for run in run_set.findall("RUN"):
        srr = run.get("accession", "")
        if not srr:
            continue
        out.append({
            "srr_id": srr,
            "srx_id": srx,
            "study": study_acc,
            "title": title,
            "tissue_source": source,
            "cell_line": cell_line,
            "genotype": genotype,
            "characteristics": characteristics,
            "platform": platform,
            "layout": layout,
        })
    return out


def fetch_sra_metadata(uid_list: list[str], batch_size: int = 200) -> list[dict]:
    """Fetch metadata for SRA UIDs in batches.

    Uses Entrez.efetch for XML retrieval but parses with ElementTree —
    recent Biopython versions refuse SRA XML because it lacks a DTD.
    """
    samples = []

    for i in range(0, len(uid_list), batch_size):
        batch = uid_list[i:i + batch_size]
        print(f"  Fetching metadata batch {i // batch_size + 1} "
              f"({len(batch)} records)...", file=sys.stderr)

        handle = Entrez.efetch(db="sra", id=",".join(batch), rettype="full", retmode="xml")
        try:
            root = ET.parse(handle).getroot()
        except ET.ParseError as e:
            print(f"  Warning: failed to parse batch XML: {e}", file=sys.stderr)
            handle.close()
            time.sleep(0.4)
            continue
        handle.close()

        for pkg in root.findall("EXPERIMENT_PACKAGE"):
            try:
                samples.extend(_parse_experiment_package(pkg))
            except Exception as e:
                print(f"  Warning: Failed to parse record: {e}", file=sys.stderr)
                continue

        # Rate limiting
        time.sleep(0.4)

    return samples


def search_geo_datasets(query: str, max_results: int = 500) -> list[str]:
    """Search GEO DataSets and return associated SRA UIDs."""
    geo_query = f"({query}) AND gse[Entry Type]"
    print(f"  Searching GEO: {geo_query[:100]}...", file=sys.stderr)

    handle = Entrez.esearch(db="gds", term=geo_query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    gds_ids = record.get("IdList", [])
    print(f"  Found {len(gds_ids)} GEO series", file=sys.stderr)

    # Link GEO → SRA
    if not gds_ids:
        return []

    sra_uids = []
    for i in range(0, len(gds_ids), 100):
        batch = gds_ids[i:i + 100]
        handle = Entrez.elink(dbfrom="gds", db="sra", id=batch)
        links = Entrez.read(handle)
        handle.close()

        for linkset in links:
            for link_db in linkset.get("LinkSetDb", []):
                for link in link_db.get("Link", []):
                    sra_uids.append(link["Id"])
        time.sleep(0.4)

    print(f"  Linked to {len(sra_uids)} SRA records", file=sys.stderr)
    return sra_uids


def load_accessions(path: str) -> tuple[list[str], set[str]]:
    """Read an existing samplesheet → (unique SRX ids to fetch, wanted SRR set).

    Rows lacking an srx_id fall back to fetching by their srr_id, so every
    requested run is reachable. efetch(db=sra) resolves both accession types.
    """
    srx_ids: list[str] = []
    seen_srx: set[str] = set()
    want_srr: set[str] = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            srr = (row.get("srr_id") or "").strip()
            srx = (row.get("srx_id") or "").strip()
            if srr:
                want_srr.add(srr)
            key = srx or srr
            if key and key not in seen_srx:
                seen_srx.add(key)
                srx_ids.append(key)
    return srx_ids, want_srr


def main():
    parser = argparse.ArgumentParser(description="Query NCBI SRA for HPV-relevant RNA-seq datasets")
    parser.add_argument("--query", type=str, default=None,
                        help="Custom search query (overrides built-in tissue queries)")
    parser.add_argument("--accessions", type=str, default=None,
                        help="CSV with srx_id (and srr_id) columns. Re-fetch metadata "
                             "for exactly these accessions instead of running a search "
                             "— reproduces an existing sample set with richer attributes.")
    parser.add_argument("--output", "-o", type=str, default="samplesheet.csv",
                        help="Output CSV file path")
    parser.add_argument("--max-samples", type=int, default=0,
                        help="Maximum samples to include (0 = all)")
    parser.add_argument("--email", type=str, default=None,
                        help="Email for NCBI Entrez (required by NCBI)")
    parser.add_argument("--api-key", type=str, default=None,
                        help="NCBI API key for faster queries")
    args = parser.parse_args()

    if args.email:
        Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # --- Accession re-fetch mode -------------------------------------------
    # Re-pull metadata for an existing sample set (e.g. samplesheet_enriched.csv)
    # so it gains the full attribute blob without a fresh search that would
    # return a different cohort. efetch resolves SRA accessions directly, so we
    # fetch whole experiments by SRX and keep only the runs we started with.
    if args.accessions:
        srx_ids, want_srr = load_accessions(args.accessions)
        print(f"Re-fetching {len(srx_ids)} experiments "
              f"({len(want_srr)} runs) from {args.accessions}", file=sys.stderr)
        samples = fetch_sra_metadata(srx_ids)
        samples = [s for s in samples if s["srr_id"] in want_srr]
        missing = want_srr - {s["srr_id"] for s in samples}
        if missing:
            print(f"  WARNING: {len(missing)} requested runs not returned "
                  f"(e.g. {', '.join(sorted(missing)[:5])})", file=sys.stderr)
    else:
        all_uids = set()
        if args.query:
            # Custom query mode
            print(f"Running custom query...", file=sys.stderr)
            uids = search_sra(args.query)
            all_uids.update(uids)
        else:
            # Search across all tissue categories
            for category, queries in TISSUE_QUERIES.items():
                print(f"\n=== Searching {category} datasets ===", file=sys.stderr)
                for q in queries:
                    uids = search_sra(q)
                    all_uids.update(uids)
                    time.sleep(0.5)

        print(f"\nTotal unique SRA UIDs: {len(all_uids)}", file=sys.stderr)

        if not all_uids:
            print("No datasets found.", file=sys.stderr)
            sys.exit(0)

        # Fetch metadata
        print("\n=== Fetching metadata ===", file=sys.stderr)
        samples = fetch_sra_metadata(list(all_uids))

    # Deduplicate by SRR
    seen = set()
    unique_samples = []
    for s in samples:
        if s["srr_id"] not in seen:
            seen.add(s["srr_id"])
            unique_samples.append(s)

    # Apply max samples limit
    if args.max_samples > 0:
        unique_samples = unique_samples[:args.max_samples]

    # Write samplesheet
    fieldnames = ["srr_id", "srx_id", "study", "title", "tissue_source",
                  "cell_line", "genotype", "characteristics", "platform", "layout"]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(unique_samples)

    print(f"\nWrote {len(unique_samples)} samples to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
