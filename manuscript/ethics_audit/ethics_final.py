"""Consolidate the ethics provenance table for the 36 studies contributing
lesion and neoplasia libraries.

Merges the NCBI pass (ethics_lookup.tsv, best coverage) with the verified pass
(ethics_verified.tsv, cleaner sentence boundaries), then applies the corrections
established by hand during the audit. Every correction is recorded inline with
its reason, so the table can be checked rather than trusted.

Usage: python3 -I ethics_final.py
Writes ethics_provenance_final.tsv
"""
import csv
import re

# Approvals that are not research-ethics approvals.
NOT_ETHICS = re.compile(
    r"Food\s+and\s+Drug\s+Administration|\bFDA\b|European\s+Medicines"
    r"|Creative\s+Commons|This\s+waiver\s+removes|RRID|Illumina|Agilent", re.I)
COMMITTEE = re.compile(
    r"\b(ethic\w*\s+(committee|board|review|approval|permission|statement)"
    r"|institutional\s+review\s+board|\bIRB\b|research\s+ethics"
    r"|ethical\s+(approval|committee|permission|review)"
    r"|institutional\s+animal\s+care|IRB-approved)\b", re.I)
CONSENT = re.compile(r"\b(informed\s+consent|written\s+consent|consent\s+was"
                     r"|consent\s+to\s+participate|Helsinki|waiver\s+of)\b", re.I)

# Corrections verified by hand during the audit. Each entry says why.
OVERRIDE = {
    "ERP022968": dict(
        status="ORIGINATING PUBLICATION NOT IDENTIFIED",
        committee=("Not recoverable. The only publication citing this accession "
                   "(PMID 36225328) is a reanalysis that downloaded it as control "
                   "data, not the originating study, so its approval does not "
                   "cover this dataset. Submitted via EBI 2017-05-16."),
        consent="Not stated in any public record", pub="Not identified"),
    "ERP106532": dict(
        status="ORIGINATING PUBLICATION NOT IDENTIFIED",
        committee=("Not recoverable from public records. Submitted via EBI "
                   "2018-07-10 by the Division of Genomics of Development and "
                   "Disease, Institute for Molecular Bioscience, University of "
                   "Queensland, Brisbane, Australia."),
        consent="Not stated in any public record", pub="Not identified"),
    "SRP448136": dict(
        status="PUBLICATION NAMES NO COMMITTEE",
        committee=("Not stated. The associated publication mentions only a drug "
                   "regulatory approval, which is not a research ethics approval."),
        consent="Not stated"),
    "SRP520129": dict(
        status="NOT PATIENT MATERIAL",
        committee=("Not applicable: HSC5 cSCC cell line and mouse xenografts. The "
                   "publication states 'Informed Consent Statement: Not applicable'."),
        consent="Not applicable"),
    "SRP495320": dict(
        status="NOT PATIENT MATERIAL",
        committee=("Not applicable: animal study. The publication states "
                   "'Informed Consent Statement: Not applicable'."),
        consent="Not applicable"),
    "SRP528200": dict(
        status="NOT PATIENT MATERIAL",
        committee="Not applicable: organotypic 3D culture, not patient tissue.",
        consent="Not applicable"),
    "SRP564228": dict(
        status="COMMITTEE NAMED BY REFERENCE ONLY",
        committee=("The publication states that institutional review board "
                   "approvals 'were previously reported' in its references 19 and "
                   "20, and names no committee itself. Trial conducted under the "
                   "Declaration of Helsinki and ICH Good Clinical Practice; all "
                   "patients provided written informed consent."),
        consent=("All patients provided written informed consent (clinical trial "
                 "cohorts).")),
    "SRP257684": dict(
        status="COMMITTEE NOT NAMED SPECIFICALLY",
        committee=("The publication states only 'Ethics approval was granted from "
                   "local Ethics Committees', without naming them."),
        consent=("Informed consent was obtained from all patients; the patients "
                 "provided written informed consent to participate.")),
}


def best(sentences, rx):
    out, seen = [], set()
    for s in sentences:
        s = s.strip()
        if not s or s == "NOT_FOUND" or NOT_ETHICS.search(s):
            continue
        if rx.search(s) and s.lower() not in seen:
            seen.add(s.lower())
            out.append(s)
    return out[:3]


def main():
    p1 = {r["study"]: r for r in
          csv.DictReader(open("ethics_lookup.tsv"), delimiter="\t")}
    p3 = {r["study"]: r for r in
          csv.DictReader(open("ethics_verified.tsv"), delimiter="\t")}
    base = list(csv.DictReader(open("ethics_studies.csv")))

    w = csv.writer(open("ethics_provenance_final.tsv", "w"), delimiter="\t")
    w.writerow(["Study accession", "BioProject", "Direct link", "Libraries used",
                "Diagnoses", "Originating publication",
                "Ethics committee, as reported in the original study",
                "Informed consent, as reported",
                "Consent covers secondary use?", "Status"])

    tally = {}
    libs = {}
    for b in base:
        st = b["study"]
        r1, r3 = p1.get(st, {}), p3.get(st, {})
        bp = r1.get("bioproject", "NOT_FOUND")
        pmid = r1.get("pmid", "NOT_FOUND")
        pub = f"PMID {pmid}" if pmid != "NOT_FOUND" else "Not identified"

        cands = []
        for src in (r3.get("committee", ""), r1.get("ethics_sentences", "")):
            cands += re.split(r"\s+\|\|\s+|\s+\|\s+", src or "")
        comm = best(cands, COMMITTEE)
        cons = best(cands, CONSENT)

        if pmid == "NOT_FOUND":
            status = "ORIGINATING PUBLICATION NOT IDENTIFIED"
            comm_t = ("Not recoverable: deposited without an associated "
                      "publication; the BioProject record states no approval.")
            cons_t = "Not stated in any public record"
        elif comm:
            status = "COMMITTEE NAMED"
            comm_t = " ".join(comm)[:460]
            cons_t = " ".join(cons)[:300] if cons else "Not stated"
        else:
            status = "PUBLICATION NAMES NO COMMITTEE"
            comm_t = "Not stated in the associated publication"
            cons_t = " ".join(cons)[:300] if cons else "Not stated"

        o = OVERRIDE.get(st)
        if o:
            status = o["status"]
            comm_t = o["committee"]
            cons_t = o["consent"]
            pub = o.get("pub", pub)

        sec = ("Not applicable" if status == "NOT PATIENT MATERIAL"
               else "Not stated in the original record")
        link = (f"https://www.ncbi.nlm.nih.gov/bioproject/{bp}"
                if bp != "NOT_FOUND" else
                "Not resolvable; no BioProject record")
        w.writerow([st, bp, link, b["n_libs"], b["diagnoses"], pub,
                    comm_t, cons_t, sec, status])
        tally[status] = tally.get(status, 0) + 1
        libs[status] = libs.get(status, 0) + int(b["n_libs"])

    print(f"{'status':42} studies  libraries")
    for k in sorted(tally, key=lambda x: -libs[x]):
        print(f"{k:42} {tally[k]:>7}  {libs[k]:>9}")
    print(f"{'TOTAL':42} {sum(tally.values()):>7}  {sum(libs.values()):>9}")


if __name__ == "__main__":
    main()
