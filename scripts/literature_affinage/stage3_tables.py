"""Stage 3 of the free Toxoplasma Affinage run: findings into ParaScape tables, no model.

Inputs: the work folder, with batches/, findings/ and exclusions.tsv, and the audit rules in
audit.py. Outputs in <work>/tables/:

- `literature_findings.parquet`: one row per admitted finding, with pmid, year and journal.
  Kept out: findings failing A1, A3 or A4, and rows in exclusions.tsv. A partner id that is
  not on that paper's gene list (A2) is dropped from the ids and kept as a name. A finding
  whose protein the abstract does not name (A7: often a family member named only as a family)
  is kept with `named_in_abstract = False` and is not used for the labels or the pairs.
- `literature_gene_labels.parquet`: one row per gene:
  - `lit_compartment`: the compartment named by the most tier 1-2 established findings
    (ties joined with " | ");
  - `lit_host_facing`: true when an established host_facing finding exists;
  - counts per kind; `lit_n_papers`; `lit_first_year` and `lit_last_year`.
- `literature_partner_pairs.parquet`: gene pairs from established partner and complex
  findings, with the papers and the best method tier.
- `literature_doubts.parquet`: every reader doubt, a to-do list for the gene matching in
  `parascape.literature`.

Provenance for the circularity table: these are observed labels with source group
`literature_affinage`. A feature built from the same papers (literature co-mention,
attention depth) is excluded when scoring these labels.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import audit as A  # noqa: E402

SOURCE_GROUP = "literature_affinage"


def load(work: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Admitted findings and all doubts."""
    papers = {}
    for b in sorted((work / "batches").glob("batch_*.json")):
        for p in json.load(open(b)):
            papers[p["pmid"]] = p
    checks = A.audit(work)
    excl = pd.read_csv(work / "exclusions.tsv", sep="\t", dtype=str) if (work / "exclusions.tsv").exists() \
        else pd.DataFrame(columns=["pmid", "gene_id"])
    excluded = set(zip(excl.pmid, excl.gene_id))
    rows, doubts = [], []
    for f in sorted((work / "findings").glob("batch_[0-9][0-9][0-9].jsonl")):
        fails = checks[checks.batch == f.stem].set_index("line").fails.to_dict()
        for i, line in enumerate(open(f)):
            if not line.strip():
                continue
            r = json.loads(line)
            why = fails.get(i, "")
            if any(rule in why for rule in ("A1", "A3", "A4")) or (str(r["pmid"]), r["gene_id"]) in excluded:
                continue
            listed = {g["gene_id"]: g for g in papers[str(r["pmid"])]["genes"]}
            ids = [g for g in r.get("partner_gene_ids") or [] if g in listed]
            names = list(r.get("partner_names") or []) + [g for g in r.get("partner_gene_ids") or [] if g not in listed]
            p = papers[str(r["pmid"])]
            rows.append({**r, "pmid": str(r["pmid"]), "partner_gene_ids": ids, "partner_names": names,
                         "named_in_abstract": "A7" not in why,
                         "year": int(str(p["year"])[:4]) if str(p["year"])[:4].isdigit() else None,
                         "journal": p["journal"], "batch": f.stem, "source_group": SOURCE_GROUP})
        d = f.with_suffix(".doubts.jsonl")
        if d.exists():
            doubts += [{**json.loads(l), "batch": f.stem} for l in open(d) if l.strip()]
    return pd.DataFrame(rows), pd.DataFrame(doubts)


def gene_labels(F: pd.DataFrame) -> pd.DataFrame:
    """Observed labels per gene from established findings."""
    est = F[(F.status == "established") & F.named_in_abstract]
    out = []
    for gene, g in est.groupby("gene_id"):
        comp = Counter(c for c, t in zip(g.compartment, g.method_tier) if c and t <= 2)
        top = max(comp.values()) if comp else 0
        row = {"gene_id": gene,
               "lit_compartment": " | ".join(sorted(c for c, n in comp.items() if n == top)),
               "lit_host_facing": bool((g.kind == "host_facing").any()),
               "lit_n_papers": int(g.pmid.nunique()),
               "lit_first_year": g.year.min(), "lit_last_year": g.year.max()}
        for k in A.KINDS:
            row[f"lit_n_{k}"] = int((g.kind == k).sum())
        out.append(row)
    return pd.DataFrame(out)


def partner_pairs(F: pd.DataFrame) -> pd.DataFrame:
    """Unordered gene pairs from established partner and complex findings."""
    rows = []
    for r in F[(F.status == "established") & F.named_in_abstract & F.kind.isin(["partner", "complex"])].itertuples():
        for p in r.partner_gene_ids:
            if p != r.gene_id:
                a, b = sorted((r.gene_id, p))
                rows.append({"a": a, "b": b, "pmid": r.pmid, "tier": r.method_tier, "kind": r.kind})
    if not rows:
        return pd.DataFrame(columns=["a", "b", "papers", "best_tier", "kinds"])
    P = pd.DataFrame(rows)
    return (P.groupby(["a", "b"]).agg(papers=("pmid", lambda s: ";".join(sorted(set(s)))),
                                      best_tier=("tier", "min"), kinds=("kind", lambda s: ";".join(sorted(set(s)))))
            .reset_index())


def main(argv=None) -> int:
    """Write the four tables and print their sizes."""
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        raise SystemExit("usage: stage3_tables.py WORK_FOLDER")
    work = Path(args[0]); out = work / "tables"; out.mkdir(exist_ok=True)
    F, D = load(work)
    L, P = gene_labels(F), partner_pairs(F)
    F.to_parquet(out / "literature_findings.parquet", index=False)
    L.to_parquet(out / "literature_gene_labels.parquet", index=False)
    P.to_parquet(out / "literature_partner_pairs.parquet", index=False)
    D.to_parquet(out / "literature_doubts.parquet", index=False)
    print(json.dumps({"findings": len(F), "genes": int(F.gene_id.nunique()), "papers": int(F.pmid.nunique()),
                      "genes_with_compartment": int((L.lit_compartment != "").sum()),
                      "host_facing_genes": int(L.lit_host_facing.sum()), "partner_pairs": len(P),
                      "doubts": len(D), "by_kind": F.kind.value_counts().to_dict()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
