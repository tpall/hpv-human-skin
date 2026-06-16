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
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from cell_line_patterns import (
    classify_cell_line,
    classify_engineered,
    has_explicit_cell_line_field,
)


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
    cell_line_count = 0
    engineered_count = 0

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

            # Detection text for cell-line / engineered heuristics: title +
            # the FULL attribute blob (`characteristics`) plus cell_line /
            # genotype. Culture hints are routinely misfiled into arbitrary
            # tags (genotype="TOP2A knockdown", note="HeLa-derived"), so we
            # scan everything — not the tissue-classification whitelist. This
            # blob is deliberately kept out of tissue-category classification
            # above so off-spec values can't skew nahk/muu assignment.
            # `characteristics` is absent in pre-genotype samplesheets; the
            # cell_line/genotype/tissue_source fallbacks keep older inputs working.
            detect_text = " ".join(filter(None, (
                combined_text,
                row.get("characteristics", ""),
                row.get("cell_line", ""),
                row.get("genotype", ""),
            )))

            # Cell-line detection: trust the explicit SRA cell_line attribute
            # when present; otherwise fall back to keyword heuristic on
            # title + tissue_source + cell_line + genotype.
            explicit = has_explicit_cell_line_field(row.get("cell_line", ""))
            heuristic, _ = classify_cell_line(detect_text)
            row["is_cell_line"] = "true" if (explicit or heuristic) else "false"
            if row["is_cell_line"] == "true":
                cell_line_count += 1

            # Engineered cells (transduced / shRNA / siRNA / CRISPR). Flagged
            # separately so reports can split clinical signal from in-vitro
            # experiments without conflating with naming-based cell-line hits.
            is_eng, _ = classify_engineered(detect_text)
            row["is_engineered"] = "true" if is_eng else "false"
            if is_eng:
                engineered_count += 1

            # Flag for manual curation
            row["needs_curation"] = flag_for_curation(row)
            if row["needs_curation"]:
                flagged_count += 1

            enriched.append(row)

    # Write enriched samplesheet
    if not enriched:
        print("No samples to process.", file=sys.stderr)
        sys.exit(0)

    # Slim fieldset consumed by Nextflow's splitCsv — title/tissue_source/
    # platform are dropped because they contain commas that column-shift
    # the CSV when splitCsv doesn't honour RFC 4180 quoting. Free-text
    # columns live in the raw samplesheet, read separately by the R report.
    fieldnames = ["srr_id", "srx_id", "study", "layout",
                  "tissue_category", "diagnosis", "is_cell_line",
                  "is_engineered", "needs_curation"]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(enriched)

    # Report
    total = len(enriched)
    print(f"\nProcessed {total} samples:", file=sys.stderr)
    for cat, count in sorted(stats.items()):
        pct = 100 * count / total if total else 0
        print(f"  {cat:15s}: {count:5d} ({pct:.1f}%)", file=sys.stderr)
    print(f"  {'flagged':15s}: {flagged_count:5d} ({100 * flagged_count / total:.1f}%)", file=sys.stderr)
    print(f"  {'cell_line':15s}: {cell_line_count:5d} ({100 * cell_line_count / total:.1f}%)", file=sys.stderr)
    print(f"  {'engineered':15s}: {engineered_count:5d} ({100 * engineered_count / total:.1f}%)", file=sys.stderr)
    print(f"\nOutput: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
