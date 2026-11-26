#!/usr/bin/env python3
"""Fetch SRA library-construction metadata for an existing accession set.

Sequencing depth is modelled as a covariate throughout the manuscript, but
library chemistry is not, and poly-A selection versus rRNA depletion is
arguably the larger determinant of viral transcript recovery: HPV transcripts
are polyadenylated, but late (L1) transcripts and the overall recovery of a
low-abundance viral pool depend on how the library was built.

This pulls `LibrarySelection` / `LibrarySource` / `LibraryStrategy` from SRA
RunInfo for the exact accessions already analysed (no fresh search, so the
cohort is unchanged), so the field can be joined onto
`results_reframe_fixed/hpv_rpm_per_sample.tsv` and entered as a model term.

Caveat the caller must handle: `LibrarySelection` is submitter-supplied and
frequently uninformative ("cDNA", "RANDOM", "unspecified"). The values that
actually separate the two chemistries are "PolyA" and "Inverse rRNA". Expect a
large uninformative stratum and treat it as its own level rather than dropping
it.

Usage:
  bin/fetch_library_selection.py [--input TSV] [--out TSV] [--id-col srr_id]

Env:
  NCBI_API_KEY   raises the eutils rate limit from 3/s to 10/s (optional)
  NCBI_EMAIL     contact address sent with each request
"""
import argparse
import csv
import io
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "")
EMAIL = os.environ.get("NCBI_EMAIL", "tapa741@gmail.com")
DELAY = 0.15 if API_KEY else 0.40

# efetch rettype=runinfo returns CSV with no header row when queried by an
# explicit id list, so the schema has to be supplied. Stable SRA RunInfo order,
# shared with presence_vs_activity/scope_paired_dna_rna.py.
RUNINFO_COLS = [
    "Run", "ReleaseDate", "LoadDate", "spots", "bases", "spots_with_mates",
    "avgLength", "size_MB", "AssemblyName", "download_path", "Experiment",
    "LibraryName", "LibraryStrategy", "LibrarySelection", "LibrarySource",
    "LibraryLayout", "InsertSize", "InsertDev", "Platform", "Model",
    "SRAStudy", "BioProject", "Study_Pubmed_id", "ProjectID", "Sample",
    "BioSample", "SampleType", "TaxID", "ScientificName", "SampleName",
    "g1k_pop_code", "source", "g1k_analysis_group", "Subject_ID", "Sex",
    "Disease", "Tumor", "Affection_Status", "Analyte_Type",
    "Histological_Type", "Body_Site", "CenterName", "Submission",
    "dbgap_study_accession", "Consent", "RunHash", "ReadHash",
]

KEEP = ["Run", "LibraryStrategy", "LibrarySelection", "LibrarySource",
        "LibraryLayout", "LibraryName", "Platform", "Model", "avgLength",
        "spots", "bases", "Experiment", "SRAStudy"]


def _post(endpoint, params, tries=4, timeout=300):
    """eutils call with the parameters in the body.

    A 300-accession id list blows past the URI length limit as a GET (HTTP 414),
    and eutils accepts the same parameters by POST.
    """
    params = dict(params)
    if API_KEY:
        params["api_key"] = API_KEY
    params["email"] = EMAIL
    params["tool"] = "hpv-human-skin"
    data = urllib.parse.urlencode(params).encode()
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(endpoint, data=data)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:          # transient eutils failures are routine
            last = e
            time.sleep(2.0 * (i + 1))
    raise RuntimeError(f"eutils failed after {tries} tries: {last}")


def read_ids(path, id_col):
    delim = "\t" if path.endswith((".tsv", ".tab")) else ","
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        if id_col not in (rdr.fieldnames or []):
            sys.exit(f"no column '{id_col}' in {path}; have {rdr.fieldnames}")
        seen, out = set(), []
        for row in rdr:
            acc = (row[id_col] or "").strip()
            if acc and acc not in seen:
                seen.add(acc)
                out.append(acc)
    return out


def fetch_batch(accs):
    text = _post(f"{EUTILS}/efetch.fcgi",
                 {"db": "sra", "rettype": "runinfo", "retmode": "csv",
                  "id": ",".join(accs)})
    header = None
    rows = []
    for row in csv.reader(io.StringIO(text)):
        if not row or not row[0]:
            continue
        if row[0] == "Run":
            header = row               # some responses do carry one
            continue
        rows.append(dict(zip(header or RUNINFO_COLS, row)))
    return rows


# ---------------------------------------------------------------------------
# Free-text library-construction protocol.
#
# `LibrarySelection` is generic ("cDNA") for about three quarters of this
# cohort, including the graded cervical validation series, whose SRA record says
# "cDNA" while its construction protocol says "Ribosomal RNA was deleted". The
# protocol free text is therefore the more informative field, and the only way
# to resolve the poly-A / rRNA-depletion split for most libraries.
# ---------------------------------------------------------------------------

XML_FIELDS = ["design_description", "library_construction_protocol",
              "instrument_model", "experiment_title"]


def _tag(text, tag):
    m = re.search(rf"<{tag}>(.*?)</{tag}>", text, re.S)
    if not m:
        return ""
    v = re.sub(r"\s+", " ", m.group(1)).strip()
    for a, b in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"), ("&quot;", '"'),
                 ("&apos;", "'"), ("&#39;", "'")):
        v = v.replace(a, b)
    return v


def fetch_protocols(accs):
    """Return {run_accession: {field: text}} from full SRA XML.

    Parsed with regex over each EXPERIMENT_PACKAGE rather than a DOM parse:
    the packages are independent, the four fields are single-occurrence within
    one, and streaming keeps ~130 MB of XML out of memory.
    """
    text = _post(f"{EUTILS}/efetch.fcgi",
                 {"db": "sra", "rettype": "full", "retmode": "xml",
                  "id": ",".join(accs)})
    out = {}
    for pkg in text.split("<EXPERIMENT_PACKAGE>")[1:]:
        rec = {"design_description": _tag(pkg, "DESIGN_DESCRIPTION"),
               "library_construction_protocol":
                   _tag(pkg, "LIBRARY_CONSTRUCTION_PROTOCOL"),
               "instrument_model": _tag(pkg, "INSTRUMENT_MODEL"),
               "experiment_title": _tag(pkg, "TITLE")}
        runs = re.search(r"<RUN_SET.*?</RUN_SET>", pkg, re.S)
        for acc in re.findall(r'<RUN [^>]*accession="([DES]RR\d+)"',
                              runs.group(0) if runs else ""):
            out[acc] = dict(rec, srr_id=acc)
    return out


def run_protocol_pass(want, out_path, batch):
    have = {}
    if os.path.exists(out_path):
        with open(out_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                have[row["srr_id"]] = row
        print(f"{len(have)} protocols already cached in {out_path}",
              file=sys.stderr)

    todo = [a for a in want if a not in have]
    for start in range(0, len(todo), batch):
        chunk = todo[start:start + batch]
        have.update(fetch_protocols(chunk))
        print(f"  protocols {min(start + batch, len(todo))}/{len(todo)}"
              f" (have {len(have)})", file=sys.stderr)
        time.sleep(DELAY)

    cols = ["srr_id"] + XML_FIELDS
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for acc in want:
            if acc in have:
                w.writerow(have[acc])
    got = sum(1 for a in want if a in have)
    withtext = sum(1 for a in want if a in have and
                   (have[a].get("library_construction_protocol") or
                    have[a].get("design_description")))
    print(f"wrote {got}/{len(want)} to {out_path}; {withtext} carry protocol"
          f" or design free text", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="results_reframe_fixed/hpv_rpm_per_sample.tsv")
    ap.add_argument("--out", default="results_reframe_fixed/library_selection.tsv")
    ap.add_argument("--id-col", default="srr_id")
    ap.add_argument("--batch", type=int, default=300)
    ap.add_argument("--protocol-out",
                    default="results_reframe_fixed/library_protocol.tsv",
                    help="second pass: full XML free-text protocol fields")
    ap.add_argument("--skip-protocol", action="store_true")
    args = ap.parse_args()

    want = read_ids(args.input, args.id_col)
    print(f"{len(want)} accessions from {args.input}", file=sys.stderr)

    # Resume: keep whatever a previous (possibly interrupted) run already wrote.
    have = {}
    if os.path.exists(args.out):
        with open(args.out, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                have[row["srr_id"]] = row
        print(f"{len(have)} already cached in {args.out}", file=sys.stderr)

    todo = [a for a in want if a not in have]
    for start in range(0, len(todo), args.batch):
        chunk = todo[start:start + args.batch]
        for row in fetch_batch(chunk):
            rec = {"srr_id": row.get("Run", "")}
            rec.update({k: row.get(k, "") for k in KEEP if k != "Run"})
            if rec["srr_id"]:
                have[rec["srr_id"]] = rec
        print(f"  runinfo {min(start + args.batch, len(todo))}/{len(todo)}"
              f" (have {len(have)})", file=sys.stderr)
        time.sleep(DELAY)

    cols = ["srr_id"] + [k for k in KEEP if k != "Run"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for acc in want:
            if acc in have:
                w.writerow(have[acc])

    missing = [a for a in want if a not in have]
    print(f"wrote {len(want) - len(missing)}/{len(want)} to {args.out}",
          file=sys.stderr)
    if missing:
        print(f"MISSING {len(missing)}: {','.join(missing[:20])}"
              f"{' ...' if len(missing) > 20 else ''}", file=sys.stderr)

    if not args.skip_protocol:
        # Smaller batches: full XML is ~12 kB per run against ~0.4 kB for
        # runinfo, so 300 at a time is a large and slow response.
        run_protocol_pass(want, args.protocol_out, min(args.batch, 150))


if __name__ == "__main__":
    main()
