#!/usr/bin/env python3
"""
Parse and enrich samplesheet metadata with tissue category classification.

Reads the raw samplesheet from query_sra.py, classifies each sample into
tissue categories (nahk/anogenitaal/suuoos/muu), and extracts diagnosis info.

Usage:
    parse_metadata.py --input samplesheet.csv --categories tissue_categories.csv --output samplesheet_enriched.csv
"""

import argparse
import csv
import re
import sys


def load_tissue_categories(categories_file: str) -> dict[str, str]:
    """Load keyword → category mapping from CSV."""
    mapping = {}
    with open(categories_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            keyword = row["keyword"].strip().lower()
            category = row["category"].strip()
            mapping[keyword] = category
    return mapping


def classify_tissue(text: str, mapping: dict[str, str]) -> str:
    """Classify tissue text into category using keyword matching."""
    text_lower = text.lower()

    # Try exact substring matches, longest first for specificity
    sorted_keywords = sorted(mapping.keys(), key=len, reverse=True)
    for keyword in sorted_keywords:
        if keyword in text_lower:
            return mapping[keyword]

    return "muu"


def extract_diagnosis(title: str, tissue_source: str) -> str:
    """Extract disease/diagnosis information from sample metadata."""
    combined = f"{title} {tissue_source}".lower()

    # Common pathology patterns
    pathology_terms = [
        "carcinoma", "cancer", "tumor", "tumour", "malignant",
        "squamous cell carcinoma", "scc", "basal cell carcinoma", "bcc",
        "melanoma", "keratosis", "actinic keratosis",
        "wart", "verruca", "papilloma", "condyloma",
        "dysplasia", "cin", "neoplasia", "intraepithelial",
        "psoriasis", "eczema", "dermatitis", "lichen",
        "normal", "healthy", "control",
    ]

    found = []
    for term in pathology_terms:
        if term in combined:
            found.append(term)

    if not found:
        return "unspecified"

    # Prefer the most specific match
    return "; ".join(sorted(set(found), key=len, reverse=True))


def flag_for_curation(row: dict) -> bool:
    """Flag samples that may need manual review."""
    # Flag if tissue source is empty or ambiguous
    if not row.get("tissue_source", "").strip():
        return True
    # Flag if classified as 'muu' (other)
    if row.get("tissue_category") == "muu":
        return True
    # Flag if diagnosis is unspecified
    if row.get("diagnosis") == "unspecified":
        return True
    return False


def main():
    parser = argparse.ArgumentParser(description="Enrich samplesheet with tissue categories")
    parser.add_argument("--input", "-i", required=True, help="Input samplesheet CSV")
    parser.add_argument("--categories", "-c", required=True, help="Tissue categories CSV")
    parser.add_argument("--output", "-o", required=True, help="Output enriched samplesheet CSV")
    args = parser.parse_args()

    # Load tissue category mapping
    mapping = load_tissue_categories(args.categories)
    print(f"Loaded {len(mapping)} tissue keywords", file=sys.stderr)

    # Process samples
    enriched = []
    stats = {"nahk": 0, "anogenitaal": 0, "suuoos": 0, "muu": 0}
    flagged_count = 0

    with open(args.input) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Combine title and tissue_source for classification
            combined_text = f"{row.get('title', '')} {row.get('tissue_source', '')}"

            # Classify tissue
            category = classify_tissue(combined_text, mapping)
            row["tissue_category"] = category
            stats[category] += 1

            # Extract diagnosis
            row["diagnosis"] = extract_diagnosis(
                row.get("title", ""),
                row.get("tissue_source", "")
            )

            # Flag for manual curation
            row["needs_curation"] = flag_for_curation(row)
            if row["needs_curation"]:
                flagged_count += 1

            enriched.append(row)

    # Write enriched samplesheet
    if not enriched:
        print("No samples to process.", file=sys.stderr)
        sys.exit(0)

    fieldnames = list(enriched[0].keys())
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(enriched)

    # Report
    total = len(enriched)
    print(f"\nProcessed {total} samples:", file=sys.stderr)
    for cat, count in sorted(stats.items()):
        pct = 100 * count / total if total else 0
        print(f"  {cat:15s}: {count:5d} ({pct:.1f}%)", file=sys.stderr)
    print(f"  {'flagged':15s}: {flagged_count:5d} ({100 * flagged_count / total:.1f}%)", file=sys.stderr)
    print(f"\nOutput: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
