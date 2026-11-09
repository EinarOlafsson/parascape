"""Verify each source origin paper by looking for its accession in the paper itself (instruction 65.04, batch 4).

A source is origin-verified when the paper recorded for it actually names the deposit the data came from. For
every registered source that has both an accession and a resolved publication, Europe PMC is searched for the
accession; if the recorded publication is among the papers that cite that accession, the origin paper is
confirmed. Where the accession is found only in other papers, those are listed instead of a verdict being
guessed, and a source with no accession or no publication stays unverified with the reason.

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
        hits = _search(f'"{first}"')
        found = {h.get("pmid") for h in hits if h.get("pmid")}
        shared = sorted(pmids & found)
        if shared:
            verdict, detail = "verified", f"{first} is named in PMID {', '.join(shared)}, the recorded paper"
        elif found:
            verdict = "accession found in other papers"
            detail = (f"{first} appears in PMID {', '.join(sorted(found)[:5])}, not in the recorded "
                      f"{', '.join(sorted(pmids))} -- read before accepting")
        else:
            verdict, detail = "accession not found in any paper text", f"Europe PMC returned no paper naming {first}"
        rows.append(dict(dataset_id=r.dataset_id, accession=first, verdict=verdict, detail=detail))
        print(r.dataset_id, verdict, flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "origin_paper_verification.tsv"), sep="\t", index=False)
    verified = set(frame.loc[frame.verdict == "verified", "dataset_id"])
    queue["origin_paper_verified"] = queue.dataset_id.isin(verified)
    queue.to_csv(queue_path, index=False)
    print(frame.verdict.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
