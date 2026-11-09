"""The verdict sheet for the biological questions (instruction 59), as a document Einar can answer in one pass.

27 of the 42 run questions carry no verdict. A verdict is a biological reading -- is this result true, useful,
or an artefact -- and it cannot be computed, so the useful thing is to put everything needed to make that call
in one place: the question in its own words, the anchor paper, how the run graded itself, what it measured
against chance where it tested itself, and the top genes it actually called, with their products.

The grade is the software's own statement about whether the method worked. The verdict is the biology, and the
two are deliberately separate columns: a reliable run can produce a biologically empty answer, and a weak run
can still point at the right gene.

    python scripts/build_verdict_sheet.py
"""
from __future__ import annotations

import datetime as _dt
import glob
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
BATCH = os.path.join(ROOT, "results", "questions_2026_11_06")
OUT = os.path.join(ROOT, "results", "verdict_sheet_2026_11_09")
TOP = 6


def _products():
    """gene id -> product, for both parasites, so a call can be read as biology."""
    names = {}
    for table in ("nodes.parquet", "pf_nodes.parquet"):
        path = os.path.join(ROOT, "parascape", "data", table)
        if not os.path.exists(path):
            continue
        frame = pd.read_parquet(path, columns=["gene_id", "product"])
        names.update(dict(zip(frame.gene_id.astype(str), frame["product"].astype(str))))
    return names


def _top_calls(qid: str, products: dict) -> str:
    """The first few genes this question called, named, from whichever table it wrote."""
    files = sorted(glob.glob(os.path.join(BATCH, f"{qid}_*.csv")))
    for path in files:
        try:
            frame = pd.read_csv(path, nrows=200)
        except Exception:
            continue
        column = next((c for c in frame.columns if c.lower() in ("gene_id", "gene", "id")), None)
        if column is None:
            continue
        genes = [str(g) for g in frame[column].head(TOP).tolist()]
        if not genes:
            continue
        named = [f"{g} ({products.get(g, 'product unknown')[:48]})" for g in genes]
        return os.path.basename(path) + ": " + "; ".join(named)
    return "no per-gene table written"


def main() -> int:
    """Write the sheet as a TSV and a readable Markdown document."""
    os.makedirs(OUT, exist_ok=True)
    summary = pd.read_csv(os.path.join(BATCH, "summary.csv"))
    pending = summary[summary.verdict.isna()].copy()
    products = _products()
    pending["top_calls"] = [_top_calls(q, products) for q in pending.qid]
    columns = ["qid", "organism", "strategy", "question", "anchor", "grade", "metric", "observed",
               "chance", "p", "skill", "summary", "top_calls"]
    sheet = pending[[c for c in columns if c in pending.columns]]
    sheet.to_csv(os.path.join(OUT, "verdict_sheet.tsv"), sep="\t", index=False)

    lines = ["# Verdict sheet: the biological questions awaiting Einar's reading (instruction 59)", "",
             f"Built {_dt.date.today().isoformat()} from `results/questions_2026_11_06/summary.csv`. "
             f"**{len(pending)} of {len(summary)} questions have no verdict.**", "",
             "`grade` is the software's own statement about whether the method worked. **The verdict is the "
             "biology** -- is the answer true, useful, or an artefact -- and only you can give it. Write PASS, "
             "FAIL or INCONCLUSIVE against each qid and I will record it.", ""]
    for row in sheet.itertuples(index=False):
        lines.append(f"## {row.qid} · {row.organism} · {row.strategy} · graded *{row.grade}*")
        lines.append("")
        lines.append(f"**{row.question}**")
        lines.append("")
        if getattr(row, "anchor", "") and str(row.anchor) != "nan":
            lines.append(f"* Anchor: {row.anchor}")
        measured = [(label, getattr(row, key, None)) for label, key in
                    (("measured", "observed"), ("chance", "chance"), ("p", "p"), ("skill", "skill"))]
        measured = [f"{label} {value}" for label, value in measured
                    if value is not None and str(value) != "nan"]
        if measured:
            metric = getattr(row, "metric", "")
            lines.append(f"* Self-test ({metric}): " + ", ".join(measured))
        lines.append(f"* What the run says: {str(row.summary)[:400]}")
        lines.append(f"* It called: {row.top_calls}")
        lines.append("")
        lines.append("**Verdict:** ______")
        lines.append("")
    with open(os.path.join(OUT, "README.md"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"{len(pending)} questions awaiting a verdict -> {OUT}/README.md")
    print(pending.grade.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
