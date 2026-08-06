#!/usr/bin/env python3
"""Find skin metatranscriptomes: a compartment-matched RNA axis for the skin
DNA reservoir.

WHY (PLAN.md 5.5b). The skin arm of paper 2 compares a DNA reservoir that is
100% surface swabs (0 of 941 Shen skin BioSamples are biopsy or tissue) against
an RNA axis that is 100% bulk tissue biopsy. Papillomavirus transcription is
differentiation-coupled and completes in the upper epidermal layers, which is
exactly what a swab collects and what a biopsy dilutes with dermis and basal
keratinocytes. So the two axes differ in the one variable that most directly
controls the read-out, and the apparent decoupling may be a compartment
artefact rather than biology.

The paired run made that concrete: both HPV-positive paired skin specimens were
swab metatranscriptomes and both were PRODUCTIVE, with beta/gamma types at up to
99.5% genome coverage. But they come from RAG-deficient patients, so host
immunodeficiency is confounded with compartment at n=2.

This script builds the cohort that separates them: skin metatranscriptomes,
stratified by host immune status. The decision it serves:

  swab RNA productive in IMMUNOCOMPETENT hosts
      -> HPV activity is stratified by tissue depth, and the reservoir is not
         silent, it was being sampled in the wrong compartment. Better paper
         than the one planned, and R2 has to be rewritten.
  swab RNA silent in immunocompetent hosts, productive only in RAG/EV
      -> the decoupling survives and gains a host-status axis.

Host status is read from BioSample attributes AND the study title/abstract,
because immunodeficiency is usually a property of the STUDY (a cohort of
transplant recipients) rather than of a per-run field.

Outputs a per-run table with the host-status call and its evidence, plus a
per-study summary to pick cohorts from. Downloads no sequence data.

Usage:
    presence_vs_activity/scope_skin_metatranscriptome.py \
        --out presence_vs_activity/skin_metatranscriptomes.tsv
"""
import argparse
import collections
import csv
import importlib.util
import os
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "scope_paired", _HERE / "scope_paired_dna_rna.py")
_scope = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_scope)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE = _HERE / "cache"

SKIN_ORG = '("human skin metagenome"[Organism] OR "skin metagenome"[Organism])'
RNA_TERM = ('("metatranscriptomic"[Source] OR "transcriptomic"[Source] '
            'OR "RNA-Seq"[Strategy])')

ATTR_FIELDS = ["env_broad_scale", "env_local_scale", "env_medium",
               "isolation_source", "samp_collect_device", "host",
               "host_disease", "host_health_state", "disease", "phenotype",
               "body_site", "host_body_site", "host_subject_id", "treatment",
               "description"]

# Host immune status. Order matters: the first match wins, and the
# immunodeficiency patterns are checked before the healthy ones, because a
# study of transplant recipients will also say "healthy controls".
IMMUNOCOMPROMISED = [
    ("primary immunodeficiency (RAG/SCID)",
     r"\brag[-\s]?[12]?\b|\bscid\b|omenn|artemis[-\s]deficien|primary immunodeficien"),
    ("epidermodysplasia verruciformis",
     r"verruciform|\bev\b patient|tmc6|tmc8|\bever[12]\b|cib1"),
    ("transplant / iatrogenic",
     r"transplant|allograft|immunosuppress|azathioprine|ciclosporin|cyclosporin|tacrolimus"),
    ("HIV",
     r"\bhiv\b|aids\b|antiretroviral"),
    ("other inborn error / immune defect",
     r"immunodeficien|agammaglobulin|hyper[-\s]?ige|wiskott|dock8|stat[13]\b|gata2"),
    ("haematological malignancy / chemotherapy",
     r"leukaemia|leukemia|lymphoma|chemotherap|stem cell transplant|\bhsct\b"),
]
INFLAMMATORY = [
    ("atopic dermatitis / eczema", r"atopic dermatitis|\bad\b patient|eczema"),
    ("psoriasis", r"psorias"),
    ("acne / seborrhoeic", r"\bacne\b|seborrh"),
    ("other skin disease", r"rosacea|vitiligo|hidradenitis|ichthyos"),
]
HEALTHY = [
    ("healthy / control", r"healthy|\bcontrol\b|normal skin|unaffected|volunteer"),
]


# The "skin metagenome" taxid is a free-for-all, so the raw query is mostly
# unusable and every one of these filters was forced by something in the result:
#
#   non-transcriptomic source  1,126 runs carry LibrarySource METAGENOMIC or
#                              GENOMIC with LibraryStrategy "RNA-Seq". They are
#                              DNA. Strategy cannot be trusted here, Source can.
#   amplicon                     465 runs are 16S rRNA amplicon. A 16S library
#                              contains no HPV sequence whatever the infection
#                              status, so it is a structural zero, exactly like
#                              the SRP572392 urine cohort excluded from paper 1.
#   non-human host               150 runs are horse (equine pastern dermatitis),
#                              mouse, sea star, dog, fallow deer, newt,
#                              crocodile lizard and frog. All filed under
#                              "skin metagenome".
#   long-read platform         the pipeline is STAR after fastp; see
#                              build_paired_samplesheet.py, where a Nanopore run
#                              killed an 11-hour job.
#
# Without these the cohort is 2,228 runs of mostly animals and 16S. With them it
# is 481 human skin transcriptomes.
RNA_SOURCES_STRICT = {"TRANSCRIPTOMIC", "METATRANSCRIPTOMIC"}
SHORT_READ_PLATFORMS = {"ILLUMINA", "BGISEQ", "DNBSEQ"}
HUMAN_HOST = {"homo sapiens", "human"}


def usable(row, host):
    """(keep, reason_if_dropped) for one runinfo row plus its BioSample host."""
    if (row.get("LibrarySource") or "").strip().upper() not in RNA_SOURCES_STRICT:
        return False, "not transcriptomic source"
    if (row.get("LibraryStrategy") or "").strip() == "AMPLICON":
        return False, "amplicon (16S)"
    if (row.get("Platform") or "").strip().upper() not in SHORT_READ_PLATFORMS:
        return False, "long-read or other platform"
    h = (host or "").strip().lower()
    if not h:
        return False, "host unstated"
    if h not in HUMAN_HOST:
        return False, "non-human host"
    return True, ""


def _post(db, ids, rettype="full"):
    data = urllib.parse.urlencode(
        {"db": db, "id": ",".join(ids), "rettype": rettype,
         "retmode": "xml"}).encode()
    for _ in range(4):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(EUTILS + "efetch.fcgi", data=data),
                timeout=180).read().decode("utf-8", "replace")
        except Exception:
            time.sleep(4)
    return ""


def fetch_biosample_attrs(biosamples, cache_path):
    """BioSample -> {attr: value}, cached on disk like the triangulation tool."""
    out = {}
    if cache_path.exists():
        with open(cache_path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                out[row["biosample"]] = {k: row.get(k, "") for k in ATTR_FIELDS}
    missing = [b for b in biosamples if b not in out]
    print(f"  biosample attrs: cache has {len(out)}, fetching {len(missing)}",
          file=sys.stderr)
    for i in range(0, len(missing), 150):
        chunk = missing[i:i + 150]
        xml = _post("biosample", chunk)
        for mo in re.finditer(r"<BioSample\b.*?</BioSample>", xml, re.S):
            try:
                node = ET.fromstring(mo.group(0))
            except ET.ParseError:
                continue
            d = {}
            for a in node.findall(".//Attribute"):
                key = (a.get("harmonized_name") or a.get("attribute_name") or "").lower()
                if key in ATTR_FIELDS and a.text:
                    d[key] = re.sub(r"\s+", " ", a.text.strip())
            desc = node.find(".//Description/Title")
            if desc is not None and desc.text:
                d["description"] = re.sub(r"\s+", " ", desc.text.strip())
            out[node.get("accession")] = {k: d.get(k, "") for k in ATTR_FIELDS}
        print(f"    {min(i + 150, len(missing))}/{len(missing)}", file=sys.stderr)
        time.sleep(0.4)
    CACHE.mkdir(exist_ok=True)
    with open(cache_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["biosample"] + ATTR_FIELDS)
        for b in sorted(out):
            w.writerow([b] + [out[b].get(k, "") for k in ATTR_FIELDS])
    return out


def fetch_study_text(studies):
    """SRA study accession -> title + abstract, where immune status usually lives.

    ONE STUDY PER QUERY, deliberately. Batching study accessions into a single
    esearch returns every RUN in those studies (thousands), and an efetch with
    any sane retmax then truncates to the first few hundred experiment packages,
    so most studies silently come back with no text at all. That failure is
    invisible downstream: every unfetched study just reads as "unstated". With
    tens of studies the per-study query is cheap and cannot truncate.
    """
    out = {}
    accs = sorted(studies)
    for i, acc in enumerate(accs, 1):
        try:
            count, we, qk = _scope.esearch_history(f"{acc}[Accession]")
        except RuntimeError:
            continue
        time.sleep(_scope.DELAY)
        if not count:
            continue
        q = urllib.parse.urlencode({"db": "sra", "WebEnv": we, "query_key": qk,
                                    "retmax": 1, "rettype": "xml",
                                    "retmode": "xml"})
        try:
            doc = urllib.request.urlopen(f"{EUTILS}efetch.fcgi?{q}",
                                         timeout=180).read().decode("utf-8", "replace")
        except Exception:
            continue
        try:
            root = ET.fromstring(doc)
        except ET.ParseError:
            continue
        for pkg in root.iter("EXPERIMENT_PACKAGE"):
            title = pkg.find(".//STUDY_TITLE")
            abstract = pkg.find(".//STUDY_ABSTRACT")
            txt = " ".join(filter(None, [
                title.text if title is not None else "",
                abstract.text if abstract is not None else ""]))
            if txt:
                # keyed by the accession we ASKED for, not one parsed back out
                out[acc] = re.sub(r"\s+", " ", txt)[:1500]
            break
        if i % 10 == 0 or i == len(accs):
            print(f"  study text {i}/{len(accs)}", file=sys.stderr)
        time.sleep(_scope.DELAY)
    missing = [a for a in accs if a not in out]
    if missing:
        print(f"  WARNING: no study text for {len(missing)} of {len(accs)} "
              f"studies: {', '.join(missing[:8])}", file=sys.stderr)
    return out


def call_status(blob):
    for label, rx in IMMUNOCOMPROMISED:
        if re.search(rx, blob, re.I):
            return "immunocompromised", label
    for label, rx in INFLAMMATORY:
        if re.search(rx, blob, re.I):
            return "inflammatory_skin_disease", label
    for label, rx in HEALTHY:
        if re.search(rx, blob, re.I):
            return "healthy", label
    return "unstated", ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default=str(CACHE / "skin_rna_biosample_attrs.tsv"))
    args = ap.parse_args()

    term = f"{SKIN_ORG} AND {RNA_TERM} AND {_scope.PUBLIC}"
    count, we, qk = _scope.esearch_history(term)
    print(f"skin-metagenome RNA runs: {count}", file=sys.stderr)
    runs = list(_scope.fetch_runinfo(we, qk, count, label="skin-rna"))
    runs = [r for r in runs if _scope.classify(r) == "rna"]
    print(f"  RNA by LibrarySource/Strategy: {len(runs)}", file=sys.stderr)

    attrs = fetch_biosample_attrs(
        sorted({r.get("BioSample", "") for r in runs if r.get("BioSample")}),
        pathlib.Path(args.cache))
    study_text = fetch_study_text({r.get("SRAStudy", "") for r in runs
                                   if r.get("SRAStudy")})

    rows = []
    dropped = collections.Counter()
    for r in runs:
        a = attrs.get(r.get("BioSample", ""), {})
        ok, why = usable(r, a.get("host", ""))
        if not ok:
            dropped[why] += 1
            continue
        sample_blob = " ".join(a.get(k, "") for k in ATTR_FIELDS)
        study_blob = study_text.get(r.get("SRAStudy", ""), "")
        status, why = call_status(sample_blob + " " + study_blob)
        source = "sample" if call_status(sample_blob)[0] == status and sample_blob.strip() else "study"
        rows.append({
            "srr_id": r["Run"],
            "study": r.get("SRAStudy", ""),
            "bioproject": r.get("BioProject", ""),
            "biosample": r.get("BioSample", ""),
            "platform": r.get("Platform", ""),
            "strategy": r.get("LibraryStrategy", ""),
            "source": r.get("LibrarySource", ""),
            "layout": r.get("LibraryLayout", ""),
            "bases": r.get("bases", ""),
            "host_status": status,
            "status_evidence": why,
            "status_from": source if status != "unstated" else "",
            "env_medium": a.get("env_medium", ""),
            "body_site": a.get("env_local_scale", "") or a.get("body_site", ""),
            "host_disease": a.get("host_disease", "") or a.get("disease", ""),
        })

    rows.sort(key=lambda x: (x["host_status"], x["study"], x["srr_id"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print("\ndropped as unusable:", file=sys.stderr)
    for kk, nn in dropped.most_common():
        print(f"  {kk:30s} {nn}", file=sys.stderr)
    print(f"\nwrote {args.out}: {len(rows)} runs")
    print("\nhost status:")
    for k, n in collections.Counter(r["host_status"] for r in rows).most_common():
        print(f"  {k:28s} {n:5d}")
    print("\nplatform (short-read only is usable):")
    for k, n in collections.Counter(r["platform"] for r in rows).most_common():
        print(f"  {k:28s} {n:5d}")
    print("\nby study (top 20):")
    per = collections.defaultdict(lambda: collections.Counter())
    for r in rows:
        per[r["study"]][r["host_status"]] += 1
    ranked = sorted(per.items(), key=lambda kv: -sum(kv[1].values()))[:20]
    for st, c in ranked:
        tot = sum(c.values())
        detail = " ".join(f"{k}={v}" for k, v in c.most_common())
        print(f"  {st:12s} n={tot:5d}  {detail}")


if __name__ == "__main__":
    main()
