"""The three layers of the gene card: what is observed, what is predicted, what is still open.

After Di Bernardo and Cheeseman 2026 (arXiv 2607.02217), as Einar chose on 2026-11-24. A sentence
about a gene is one of three different things, and the card must never let one pass for another:

- **observed**: a paper shows it by a direct experiment (the literature findings), or a dataset
  measures it (the existing "What is known" synthesis);
- **predicted**: a method infers it, said with how often such a call is right, against chance;
- **still open**: nothing measures it yet, or a paper names it as unknown.

Pure functions over tables, so they are tested without a window.
"""
from __future__ import annotations

from html import escape

import pandas as pd

#: The order kinds of finding are read in: where it is, what it does, what loses it, what it binds.
KIND_ORDER = ("localisation", "host_facing", "activity", "phenotype", "partner", "complex", "structure")
KIND_TITLE = {"localisation": "Where it is", "host_facing": "What it does to the host",
              "activity": "What it does", "phenotype": "What losing it does",
              "partner": "What it binds", "complex": "What complex it is in",
              "structure": "Its shape, solved"}
TIER_WORDS = {1: "genetic loss of function, an assay with mutants, or a structure",
              2: "imaging, co-immunoprecipitation, pulldown or proximity labelling",
              3: "another direct experiment"}


def _pubmed(pmid) -> str:
    p = escape(str(pmid))
    return f"<a href='https://pubmed.ncbi.nlm.nih.gov/{p}/'>PMID {p}</a>"


def _year(value) -> str:
    try:
        return str(int(value))
    except (TypeError, ValueError):
        return ""


def gene_findings(gene: str, findings: pd.DataFrame | None = None) -> pd.DataFrame:
    """The admitted literature findings of one gene: named in the abstract, newest first."""
    if findings is None:
        from . import literature_findings as LF
        findings = LF.findings()
    if findings is None or not len(findings):
        return pd.DataFrame()
    f = findings[findings["gene_id"].astype(str) == str(gene)]
    if "named_in_abstract" in f:
        f = f[f["named_in_abstract"].astype(bool)]
    return f.sort_values("year", ascending=False, na_position="last")


def observed_lines(gene: str, findings: pd.DataFrame | None = None, per_kind: int = 4) -> list:
    """One HTML block per kind of finding: each sentence with its method, year and paper.

    A refutation is said as one ("Shown not to"), because a paper ruling something out is
    knowledge too, and must not read like the claim it rules out.
    """
    f = gene_findings(gene, findings)
    if not len(f):
        return []
    out = []
    for kind in KIND_ORDER:
        part = f[f["kind"].astype(str) == kind]
        if not len(part):
            continue
        items = []
        for r in part.head(per_kind).itertuples(index=False):
            refuted = str(getattr(r, "status", "")) == "refuted"
            lead = "<b>Ruled out:</b> " if refuted else ""
            year = _year(getattr(r, "year", None))
            items.append(f"{lead}{escape(str(r.finding))} <span style='color:#999'>"
                         f"({escape(str(r.method))}; {year + ', ' if year else ''}{_pubmed(r.pmid)})</span>")
        more = len(part) - per_kind
        tail = f" <span style='color:#999'>and {more} more</span>" if more > 0 else ""
        out.append(f"<b>{KIND_TITLE[kind]}</b> ({len(part)} "
                   f"{'paper' if part['pmid'].nunique() == 1 else 'papers'}): "
                   + "; ".join(items) + tail)
    return out


def observed_summary(gene: str, findings: pd.DataFrame | None = None) -> str:
    """One line saying how much the papers show, and by what strength of experiment."""
    f = gene_findings(gene, findings)
    if not len(f):
        return ""
    tiers = f["method_tier"].dropna().astype(int)
    best = int(tiers.min()) if len(tiers) else 3
    return (f"{len(f)} findings shown by experiment in {f['pmid'].nunique()} papers; the strongest kind "
            f"of evidence is {TIER_WORDS.get(best, TIER_WORDS[3])}.")


def chance_rate(ledger: pd.DataFrame, target: str) -> tuple:
    """(share, answer): how often the commonest true answer of a label occurs among genes held out.

    The rate a guess of the commonest answer would reach -- what "how often right" is measured
    against.
    """
    g = ledger[ledger["target"].astype(str) == str(target)]
    if "truth" not in g or not len(g):
        return None, ""
    truth = g.drop_duplicates("gene_id")["truth"].astype(str)
    truth = truth[truth.ne("") & truth.ne("nan")]
    if not len(truth):
        return None, ""
    counts = truth.value_counts()
    return float(counts.iloc[0] / counts.sum()), str(counts.index[0])


def predicted_lines(reading, consensus: list | None = None, chance: dict | None = None,
                    limit: int = 8) -> list:
    """What is predicted for the gene, each with its own uncertainty, never phrased as knowledge.

    :param reading: a `gene_reading.GeneReading`.
    :param consensus: the held-out agreement sentences (`gene_reading.consensus_lines`), already
        carrying how often that much agreement was right.
    :param chance: target -> (share, answer) of the commonest answer, appended so each rate is read
        against what a guess would get.
    """
    from .gene_reading import value_text
    out = []
    predicted = getattr(reading, "predicted", None)
    if predicted is not None and len(predicted):
        for row in predicted.head(limit).itertuples(index=False):
            p = "" if row.probability != row.probability else f", {row.probability:.0%} certain"
            out.append(f"Predicted, not measured: {escape(str(row.description or row.column)[:70])} "
                       f"= <b>{escape(value_text(row.column, row.value))}</b>"
                       f" ({escape(str(row.source_tier))}{p}).")
    for line, target in consensus or []:
        share, answer = (chance or {}).get(target, (None, ""))
        if share is not None:
            line += (f" A guess of the commonest answer ({escape(answer)}) would be right "
                     f"{share:.0%} of the time.")
        out.append(line)
    return out


def open_lines(reading, gene: str, findings: pd.DataFrame | None = None, limit: int = 6) -> list:
    """What nothing measures yet, and what papers about the gene say is still unknown."""
    out = [f"Not measured: {escape(phrase)}." for _column, phrase in getattr(reading, "unknown", [])]
    f = gene_findings(gene, findings)
    if len(f) and "open_question" in f:
        q = f[f["open_question"].fillna("").astype(str).str.strip().ne("")]
        q = q.drop_duplicates("open_question")
        for r in q.head(limit).itertuples(index=False):
            year = _year(getattr(r, "year", None))
            out.append(f"Left open by {year + ' ' if year else ''}{_pubmed(r.pmid)}: "
                       f"{escape(str(r.open_question))}")
    return out
