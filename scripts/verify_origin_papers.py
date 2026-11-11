"""Verify each source origin paper against the deposit itself (instruction 65.04, batch 4).

A source is origin-verified when the paper recorded for it is the paper the deposit came from. **The deposit is
asked first, because it is the authority on its own origin**: GEO carries the citation a submitter attached to a
series, and BioStudies carries it for an ArrayExpress or PRIDE accession. Only when the deposit names no paper
is Europe PMC searched for the accession as a fallback.

That order was the other way round until 2026-11-10, and it produced false negatives that read as real problems.
Searching full text and requiring the recorded paper among the hits fails whenever the origin paper does not
spell its accession in indexed text -- common in older papers and in consortium data. Two of 162 sources were
flagged "read before accepting" on that basis and both were correct: `E-MTAB-3579` is the FANTOM5 mouse CAGE
deposit and the recorded 24670764 IS the FANTOM5 paper, while the papers that do name the accession are
re-users; `GSE132248` records exactly the citation GEO lists for it. A check that calls a correct record
suspicious costs more than one that says "cannot tell".

Where the deposit and the recorded paper disagree, that is reported as a disagreement with both named, and a
source with no accession or no publication stays unverified with the reason.

    python scripts/verify_origin_papers.py

Writes `results/dataset_selection_2026_11_08/origin_paper_verification.tsv` and updates the
`origin_paper_verified` column of the review queue.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def _search(query: str, pause: float = 0.4):
    """The Europe PMC hits for one query, or an empty list when the service does not answer."""
    url = f"{EPMC}?query={urllib.parse.quote(query)}&format=json&pageSize=25"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                time.sleep(pause)
                return json.load(fh)["resultList"]["result"]
        except Exception:
            time.sleep(2 * (attempt + 1))
    return []


#: Accessions settled by reading the deposit and the paper, where neither automated route can settle them.
#: Each entry records the evidence, so the verdict can be checked rather than taken on trust. These are not
#: guesses: a deposit that declares no PubMed id and whose accession appears only in later papers needs a
#: person to say which paper produced it, and this is where that answer is kept.
READ_AND_DECIDED = {
    "E-MTAB-3579": ("24670764",
                    "the deposit is 'RNA-Seq CAGE analysis of mice tissue in RIKEN FANTOM5 project' and the "
                    "recorded 24670764 is the FANTOM5 flagship 'A promoter-level mammalian expression atlas' "
                    "(Nature 2014); the four papers that name the accession are re-users of the atlas"),
}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1/studies"


def _json(url: str, pause: float = 0.4, tries: int = 4):
    """One JSON document, or None. NCBI answers 429 under load, so this backs off rather than giving up."""
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                document = json.load(fh)
            time.sleep(pause)
            return document
        except Exception:
            time.sleep(3 * (attempt + 1))
    return None


def deposit_publications(accession: str) -> set:
    """The PubMed ids the DEPOSIT itself declares, which is the authority on where its data came from."""
    if accession.startswith(("GSE", "GSM")):
        found = _json(f"{EUTILS}/esearch.fcgi?db=gds&term="
                      f"{urllib.parse.quote(accession + '[ACCN] AND gse[ETYP]')}&retmode=json")
        ids = (found or {}).get("esearchresult", {}).get("idlist") or []
        if not ids:
            return set()
        summary = _json(f"{EUTILS}/esummary.fcgi?db=gds&id={ids[0]}&retmode=json")
        result = (summary or {}).get("result", {})
        for key, record in result.items():
            if key != "uids" and isinstance(record, dict):
                return {str(p) for p in record.get("pubmedids", []) if str(p).isdigit()}
        return set()
    # ArrayExpress and PRIDE accessions both resolve through BioStudies.
    study = _json(f"{BIOSTUDIES}/{urllib.parse.quote(accession)}")
    if not study:
        return set()
    text = json.dumps(study)
    return set(re.findall(r'"PMID"\s*,\s*"value"\s*:\s*"(\d{7,8})"', text)) | \
        set(re.findall(r'"value"\s*:\s*"(\d{7,8})"\s*,\s*"name"\s*:\s*"PMID"', text))


def main() -> int:
    """Check every source that can be checked; print a count per verdict."""
    import pandas as pd
    from parascape import datasets as ds
    queue_path = os.path.join(OUT, "source_review_queue.csv")
    queue = pd.read_csv(queue_path)
    registry = {d.key: d for d in ds.REGISTRY}
    rows = []
    for r in queue.itertuples():
        entry = registry.get(r.dataset_id)
        accession = str(getattr(entry, "accession", "") or "")
        recorded = str(r.publication_id or "")
        pmids = set(re.findall(r"\d{7,8}", recorded))
        if not accession:
            rows.append(dict(dataset_id=r.dataset_id, accession="", verdict="no accession to look for",
                             detail="the source records no deposit, so the paper cannot be confirmed this way"))
            continue
        if not pmids:
            rows.append(dict(dataset_id=r.dataset_id, accession=accession, verdict="no publication recorded",
                             detail=str(r.publication_status)))
            continue
        first = accession.split(";")[0].split(",")[0].strip()
        # Many sources record WHERE the data sits inside a named paper (a supplementary table) rather than a
        # deposit. That string already is the paper reference, so there is no accession to cross-check.
        if not re.match(r"^(PXD|GSE|GSM|E-MTAB|E-GEOD|SRP|SRR|PRJ|MSV|MTBLS)[A-Za-z0-9-]+$", first):
            rows.append(dict(dataset_id=r.dataset_id, accession=first,
                             verdict="supplementary table of a named paper",
                             detail="the recorded location is the paper itself, so there is no deposit to cross-check"))
            continue
        if first in READ_AND_DECIDED and READ_AND_DECIDED[first][0] in pmids:
            rows.append(dict(dataset_id=r.dataset_id, accession=first, verdict="verified by reading",
                             detail=READ_AND_DECIDED[first][1]))
            print(r.dataset_id, "verified by reading", flush=True)
            continue
        # The deposit next. Its own declared citation settles the question without a full-text search.
        declared = deposit_publications(first)
        agreed = sorted(pmids & declared)
        if agreed:
            verdict = "verified by the deposit"
            detail = f"{first} declares PMID {', '.join(agreed)}, which is the recorded paper"
        elif declared:
            verdict = "deposit names a different paper"
            detail = (f"{first} declares PMID {', '.join(sorted(declared))} but the source records "
                      f"{', '.join(sorted(pmids))} -- read both before accepting")
        else:
            # The deposit names nobody, so fall back to who names the deposit.
            hits = _search(f'"{first}"')
            found = {h.get("pmid") for h in hits if h.get("pmid")}
            shared = sorted(pmids & found)
            if shared:
                verdict = "verified by the paper text"
                detail = f"{first} is named in PMID {', '.join(shared)}, the recorded paper"
            elif found:
                verdict = "only other papers name the accession"
                detail = (f"{first} declares no paper of its own and appears in PMID "
                          f"{', '.join(sorted(found)[:5])}, not in the recorded {', '.join(sorted(pmids))} "
                          f"-- these may be re-users rather than the origin, so read before accepting")
            else:
                verdict = "no paper names the accession"
                detail = f"{first} declares no paper and Europe PMC returns none naming it"
        rows.append(dict(dataset_id=r.dataset_id, accession=first, verdict=verdict, detail=detail))
        print(r.dataset_id, verdict, flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "origin_paper_verification.tsv"), sep="\t", index=False)
    verified = set(frame.loc[frame.verdict.str.startswith("verified"), "dataset_id"])
    queue["origin_paper_verified"] = queue.dataset_id.isin(verified)
    queue.to_csv(queue_path, index=False)
    print(frame.verdict.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
