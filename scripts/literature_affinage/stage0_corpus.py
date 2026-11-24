"""Stage 0 of the free Toxoplasma Affinage run: each gene's papers and their abstracts, no model.

After Di Bernardo and Cheeseman 2026 (arXiv 2607.02217), adapted to ParaScape (Einar, 2026-11-24):
- the gene-to-paper links come from the release `mentions.parquet`: symbol and accession
  matches in the Toxoplasma literature, already disambiguated by `parascape.literature`;
- catalogue papers that name more than 50 genes are dropped, as in their precision layer;
- each gene's papers are ranked by how specific they are to it: named in the title first,
  then fewer genes per paper, then more mentions; at most 20 per gene are kept;
- abstracts are fetched once per paper from NCBI E-utilities (efetch, 200 ids per request,
  at most 3 requests a second, no key).

Each paper is read ONCE later, for every gene it names, rather than once per gene.

Outputs in --out:
- `abstracts.jsonl`: pmid, title, abstract, year, journal, publication types;
- `gene_papers.parquet`: gene_id, pmid, rank, in_title, genes_in_paper;
- `stage0_summary.json`.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MAX_GENES_PER_PAPER = 50
MAX_PAPERS_PER_GENE = 20


def gene_papers(mentions: pd.DataFrame) -> pd.DataFrame:
    """Rank each gene's papers by specificity and keep the top `MAX_PAPERS_PER_GENE`."""
    m = mentions.dropna(subset=["pmid"]).copy()
    m["pmid"] = m["pmid"].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
    m = m[m.pmid.str.fullmatch(r"\d+")]                     # full texts without a PubMed id
    fan = m.groupby("pmid").gene_id.nunique()
    m = m[m.pmid.map(fan) <= MAX_GENES_PER_PAPER]
    per = (m.assign(in_title=m.section.eq("title"))
           .groupby(["gene_id", "pmid"])
           .agg(in_title=("in_title", "max"), hits=("n_hits", "sum"))
           .reset_index())
    per["genes_in_paper"] = per.pmid.map(fan).astype(int)
    per = per.sort_values(["gene_id", "in_title", "genes_in_paper", "hits"],
                          ascending=[True, False, True, False])
    per["rank"] = per.groupby("gene_id").cumcount() + 1
    return per[per["rank"] <= MAX_PAPERS_PER_GENE].reset_index(drop=True)


def _text(node) -> str:
    return "".join(node.itertext()).strip() if node is not None else ""


def parse(xml: bytes) -> list[dict]:
    """Title, abstract, year, journal and publication types of each PubmedArticle."""
    out = []
    for art in ET.fromstring(xml).iter("PubmedArticle"):
        pmid = _text(art.find(".//PMID"))
        a = art.find(".//Article")
        parts = []
        for node in a.findall(".//Abstract/AbstractText") if a is not None else []:
            label = node.get("Label")
            parts.append((f"{label}: " if label else "") + _text(node))
        year = _text(art.find(".//PubDate/Year")) or _text(art.find(".//PubDate/MedlineDate"))[:4]
        out.append({"pmid": pmid, "title": _text(a.find("ArticleTitle")) if a is not None else "",
                    "abstract": "\n".join(parts), "year": year,
                    "journal": _text(art.find(".//Journal/Title")),
                    "types": [_text(t) for t in art.findall(".//PublicationType")]})
    return out


def fetch(pmids: list[str], raw: Path, log=print) -> list[dict]:
    """Fetch abstracts 200 at a time; each raw reply is kept so a rerun reads from disk."""
    raw.mkdir(parents=True, exist_ok=True)
    records = []
    for start in range(0, len(pmids), 200):
        chunk = pmids[start:start + 200]
        path = raw / f"efetch_{start:06d}.xml"
        if not path.exists():
            body = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(chunk), "retmode": "xml",
                                           "tool": "parascape"}).encode()
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(EFETCH, data=body, timeout=120) as r:
                        path.write_bytes(r.read())
                    break
                except Exception as e:                      # transient NCBI errors: back off
                    log(f"efetch {start}: {e}; retry {attempt + 1}")
                    time.sleep(5 * (attempt + 1))
            else:
                raise RuntimeError(f"efetch failed for papers {start}-{start + len(chunk)}")
            time.sleep(0.4)
        records += parse(path.read_bytes())
        log(f"{min(start + 200, len(pmids))}/{len(pmids)} papers")
    return records


def main(argv=None) -> int:
    """Build the per-gene corpus and fetch its abstracts."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mentions", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)                               # unknown arguments are an error
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    gp = gene_papers(pd.read_parquet(args.mentions))
    gp.to_parquet(out / "gene_papers.parquet", index=False)
    pmids = sorted(gp.pmid.unique(), key=int)
    recs = fetch(pmids, out / "pubmed_raw")
    with open(out / "abstracts.jsonl", "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    got = {r["pmid"] for r in recs}
    summary = {"genes": int(gp.gene_id.nunique()), "gene_paper_pairs": len(gp), "papers": len(pmids),
               "fetched": len(got), "with_abstract": sum(bool(r["abstract"]) for r in recs),
               "missing": sorted(set(pmids) - got)[:50], "mentions": args.mentions}
    (out / "stage0_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "missing"}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
