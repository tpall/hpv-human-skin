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
import sys
import time
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


def fetch_sra_metadata(uid_list: list[str], batch_size: int = 200) -> list[dict]:
    """Fetch metadata for SRA UIDs in batches."""
    samples = []

    for i in range(0, len(uid_list), batch_size):
        batch = uid_list[i:i + batch_size]
        print(f"  Fetching metadata batch {i // batch_size + 1} "
              f"({len(batch)} records)...", file=sys.stderr)

        handle = Entrez.efetch(db="sra", id=",".join(batch), rettype="full", retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for pkg in records:
            try:
                exp = pkg["EXPERIMENT_PACKAGE"]["EXPERIMENT"]
                run_set = pkg["EXPERIMENT_PACKAGE"].get("RUN_SET", {})
                sample = pkg["EXPERIMENT_PACKAGE"].get("SAMPLE", {})

                # Extract run accessions
                runs = run_set.get("RUN", []) if isinstance(run_set, dict) else []
                if not isinstance(runs, list):
                    runs = [runs]

                for run in runs:
                    srr = run.get("@accession", "")
                    if not srr:
                        continue

                    # Extract metadata
                    srx = exp.get("@accession", "")
                    platform = ""
                    try:
                        platform = exp["PLATFORM"]
                        if isinstance(platform, dict):
                            platform = list(platform.keys())[0]
                    except (KeyError, IndexError):
                        pass

                    # Library layout
                    layout = "SINGLE"
                    try:
                        lib_desc = exp["DESIGN"]["LIBRARY_DESCRIPTOR"]
                        layout_info = lib_desc.get("LIBRARY_LAYOUT", {})
                        if "PAIRED" in layout_info:
                            layout = "PAIRED"
                    except (KeyError, TypeError):
                        pass

                    # Sample attributes
                    title = ""
                    source = ""
                    try:
                        title = sample.get("TITLE", "")
                        attrs = sample.get("SAMPLE_ATTRIBUTES", {}).get("SAMPLE_ATTRIBUTE", [])
                        if not isinstance(attrs, list):
                            attrs = [attrs]
                        for attr in attrs:
                            tag = attr.get("TAG", "").lower()
                            val = attr.get("VALUE", "")
                            if tag in ("tissue", "cell_type", "source_name", "sample_type"):
                                source = val
                                break
                    except (KeyError, TypeError):
                        pass

                    # Study accession
                    study = ""
                    try:
                        study = pkg["EXPERIMENT_PACKAGE"]["STUDY"]["@accession"]
                    except (KeyError, TypeError):
                        pass

                    samples.append({
                        "srr_id": srr,
                        "srx_id": srx,
                        "study": study,
                        "title": title,
                        "tissue_source": source,
                        "platform": platform,
                        "layout": layout,
                    })

            except (KeyError, TypeError, IndexError) as e:
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


def main():
    parser = argparse.ArgumentParser(description="Query NCBI SRA for HPV-relevant RNA-seq datasets")
    parser.add_argument("--query", type=str, default=None,
                        help="Custom search query (overrides built-in tissue queries)")
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
    fieldnames = ["srr_id", "srx_id", "study", "title", "tissue_source", "platform", "layout"]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(unique_samples)

    print(f"\nWrote {len(unique_samples)} samples to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
