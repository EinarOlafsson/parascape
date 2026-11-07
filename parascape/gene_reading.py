"""Reading one gene: what is known, from where, and what is not (the Gene module, instruction 93).

The Gene module leads with a SYNTHESIS -- a sentence per kind of knowledge (where the protein
is, what it does, what it binds, how much the parasite needs it, when it is made...) with the
sources of each claim -- then the evidence itself grouped the same way, then what is not known
and what can be inferred. The grouping is the information categories of
`parascape.information_categories` (place, function, interaction, fitness, transcript...), the
same ones that keep inference from being circular, so a sentence here and an exclusion there
mean the same thing by the same word.

Everything here is a pure function of a `PublishedMap` (or anything with its `gene_facts`,
`gene_evidence` and `gene_partners`), so every sentence is testable without a window.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from html import escape

import pandas as pd

#: The categories a synthesis reads, in the order a biologist asks about a gene, each with the
#: heading of its evidence section.
ORDER = (
    ("place", "Where it is"),
    ("function", "What it does"),
    ("interaction", "What it binds"),
    ("fitness", "How much the parasite needs it"),
    ("phenotype", "What happens when it is perturbed"),
    ("transcript", "When its RNA is made"),
    ("translation", "How it is translated"),
    ("protein_abundance", "How much protein there is"),
    ("modification", "How it is modified"),
    ("structure", "Its shape"),
    ("domain", "Its domains"),
    ("conservation", "Across species and strains"),
    ("host_response", "What it does to the host"),
    ("chromatin", "Its locus"),
    ("sequence", "Its sequence"),
)
TITLES = dict(ORDER)
#: Categories never shown as knowledge about the gene.
HIDDEN = ("identity", "attention")
#: Columns whose values are whole sequences or per-residue strings: shown by length, not text.
_LONG = re.compile(r"^(seq_|tdi_)(full|mature|nterm_flank|cterm_flank)$")
#: Source tiers that are predictions rather than observations.
PREDICTED = ("filled", "exploratory", "predicted")
#: The traits a synthesis reports as unknown when the gene has no observed value for them.
KEY_TRAITS = (("host_facing", "whether it faces the host"),
              ("compartment_hyperlopit", "its compartment by hyperLOPIT"),
              ("partner_best", "a recorded binding partner"))

_PMID = re.compile(r"PMID:?\s*(\d{5,9})")
_PMC = re.compile(r"\b(PMC\d{5,9})\b")


def category_of(column: str, organism: str | None = None) -> str:
    """The first information category of a column ('' when it has none)."""
    try:
        from . import information_categories as IC
        cats = IC.categories_of(column, organism) if organism else IC.categories_of(column)
    except Exception:
        cats = ()
    if cats:
        return cats[0]
    if column.startswith(("partner_", "host_partner")):
        return "interaction"
    if column == "host_facing" or column.startswith("compartment"):
        return "place"
    return ""


def source_links(text: str) -> str:
    """References in `text` as links (PubMed, PMC); other text escaped and kept short."""
    text = str(text or "")
    out = escape(text[:160])
    out = _PMID.sub(lambda m: f"<a href='https://pubmed.ncbi.nlm.nih.gov/{m.group(1)}/'>"
                              f"PMID:{m.group(1)}</a>", out)
    out = _PMC.sub(lambda m: f"<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/{m.group(1)}/'>"
                             f"{m.group(1)}</a>", out)
    return out


def value_text(column: str, value) -> str:
    """A cell value as it should be read: sequences by length, numbers rounded, booleans as words."""
    if value is None or (isinstance(value, float) and value != value):
        return "unknown"
    if _LONG.match(column):
        return f"{len(str(value))} residues"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.3g}"
    text = str(value)
    if text.lower() in ("true", "false"):
        return "yes" if text.lower() == "true" else "no"
    return text if len(text) <= 80 else text[:77] + "…"


@dataclass
class GeneReading:
    """One gene, read: facts and evidence by category, partners, and what is missing."""

    gene: str
    product: str = ""
    facts: dict = field(default_factory=dict)        # category -> DataFrame of observed facts
    evidence: pd.DataFrame | None = None            # every source row behind its facts
    partners: pd.DataFrame | None = None
    predicted: pd.DataFrame | None = None           # cells whose source is a prediction
    unknown: list = field(default_factory=list)     # (column, phrase) of key traits not observed

    def sources_for(self, column: str, limit: int = 4) -> list:
        """The references and datasets behind one trait, strongest first."""
        if self.evidence is None or not len(self.evidence):
            return []
        rows = self.evidence[self.evidence["trait"].astype(str) == column]
        rows = rows.sort_values("probability", ascending=False)
        out = []
        for row in rows.itertuples(index=False):
            ref = str(getattr(row, "reference", "") or "")
            if ref and ref not in out:
                out.append(ref)
            if len(out) >= limit:
                break
        return out


def _dataset_of(column: str) -> tuple:
    """(dataset name, "PMID:n" or "") for a gene-table column, from the dataset registry."""
    try:
        from . import datasets as D
        d = D.provenance(column, "Tg") or D.provenance(column)
    except Exception:
        d = None
    if d is None:
        return "", ""
    pmid = str(getattr(d, "pmid", "") or "")
    return str(getattr(d, "name", "") or ""), (f"PMID:{pmid}" if pmid else "")


def node_facts(row) -> pd.DataFrame:
    """The window's gene-table row (measurements the map is built from: fitness, expression,
    protein features...) in the shape of `PublishedMap.gene_facts`, observed values only."""
    if row is None:
        return pd.DataFrame(columns=["column", "tier", "table", "description", "observed",
                                     "value", "probability", "source_tier"])
    records = []
    for column, value in row.items():
        column = str(column)
        if column in ("gene_id", "product") or value is None:
            continue
        if isinstance(value, float) and value != value:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        name, reference = _dataset_of(column)
        records.append({"column": column, "tier": 0, "table": "gene table",
                        "description": name or column.replace("_", " "), "observed": True,
                        "value": value, "probability": float("nan"),
                        "source_tier": "measured (gene table)", "dataset": name,
                        "reference": reference})
    return pd.DataFrame(records)


def read(published, gene: str, product: str = "", organism: str | None = None,
         node_row=None) -> GeneReading:
    """Everything the Gene module shows about `gene`, grouped by information category.

    `node_row` adds the window's own measurements of the gene (its row of the node table),
    which the published map does not repeat: fitness screens, expression, protein features.
    """
    reading = GeneReading(gene=gene, product=product)
    facts = published.gene_facts(gene) if published is not None else None
    extra = node_facts(node_row)
    if facts is None or not len(facts):
        facts = extra
    elif len(extra):
        extra = extra[~extra["column"].isin(set(facts["column"].astype(str)))]
        facts = pd.concat([facts, extra], ignore_index=True)
    if facts is None or not len(facts):
        reading.unknown = list(KEY_TRAITS)
        return reading
    facts = facts.copy()
    facts["category"] = [category_of(str(c), organism) for c in facts["column"]]
    tiers = facts["source_tier"].astype(str).str.lower()
    is_predicted = tiers.str.startswith(PREDICTED)
    reading.predicted = facts[is_predicted]
    observed = facts[facts["observed"].astype(bool) & ~is_predicted]
    for category, rows in observed.groupby("category", sort=False):
        if category in HIDDEN:
            continue
        reading.facts[category or "other"] = rows
    try:
        reading.evidence = published.gene_evidence(gene) if published is not None else None
    except Exception:
        reading.evidence = None
    try:
        reading.partners = published.gene_partners(gene) if published is not None else None
    except Exception:
        reading.partners = None
    have = set(observed["column"].astype(str))
    reading.unknown = [(c, phrase) for c, phrase in KEY_TRAITS if c not in have]
    return reading


def _fact(reading: GeneReading, column: str):
    for rows in reading.facts.values():
        hit = rows[rows["column"].astype(str) == column]
        if len(hit):
            return hit.iloc[0]
    return None


def _cite(reading: GeneReading, column: str) -> str:
    refs = reading.sources_for(column)
    return f" ({'; '.join(source_links(r) for r in refs)})" if refs else ""


def informative(rows: pd.DataFrame) -> pd.DataFrame:
    """The facts worth a sentence: not whole-sequence strings, and not an absence recorded as
    "no" (a domain the protein does not have is not something known about it)."""
    if not len(rows):
        return rows
    keep = []
    for row in rows.itertuples(index=False):
        column = str(row.column)
        text = value_text(column, row.value)
        keep.append(not _LONG.match(column) and not column.startswith(("tdi_", "seq_"))
                    and text not in ("no", "unknown", "0", "False"))
    return rows[keep]


def _text(value) -> str:
    """A cell as text, with missing (None or NaN) as ""."""
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value)


def _ref_of(row) -> str:
    return _text(getattr(row, "reference", ""))


def _datasets(rows: pd.DataFrame, limit: int = 3) -> str:
    """The distinct datasets behind some facts, named with their sources."""
    seen, out = set(), []
    for row in rows.itertuples(index=False):
        name = _text(getattr(row, "dataset", "")) or _text(getattr(row, "description", ""))
        if not name or name in seen:
            continue
        seen.add(name)
        ref = _ref_of(row)
        out.append(escape(name) + (f" ({source_links(ref)})" if ref else ""))
    more = len(seen) - min(len(seen), limit)
    return "; ".join(out[:limit]) + (f"; and {more} more" if more > 0 else "")


#: Categories measured as many numbers (one per condition or replicate): summarised by which
#: datasets measured them, not by reciting values out of context.
MEASURED_MANY = ("transcript", "translation", "protein_abundance", "host_response", "chromatin",
                 "modification", "conservation")


#: A column must be measured on this many genes before a percentile on it means anything.
UNUSUAL_MIN_MEASURED = 500
#: Categories that describe how much a gene was studied or what it is called, not what it does.
UNUSUAL_SKIP = {"identity", "attention"}


def unusual(nodes: pd.DataFrame, gene: str, top: int = 10, organism: str | None = None,
            per_dataset: int = 2) -> pd.DataFrame:
    """The gene's most extreme measurements, as percentiles among every gene measured (instruction 60).

    "Of everything measured on this gene, what is unusual?" -- for each continuous column measured on at
    least `UNUSUAL_MIN_MEASURED` genes, the share of measured genes below this one's value (ties split),
    ranked by distance from the middle. Yes/no columns are left out (a percentile of a flag says
    nothing), as are identity and attention columns. At most `per_dataset` lines come from one
    dataset, so ten columns of one experiment do not crowd out everything else.
    """
    cols = ["column", "value", "percentile", "measured", "category", "dataset", "rank_high",
            "rank_low"]
    if nodes is None or "gene_id" not in nodes:
        return pd.DataFrame(columns=cols)
    hit = nodes.index[nodes["gene_id"].astype(str) == str(gene)]
    if not len(hit):
        return pd.DataFrame(columns=cols)
    i = hit[0]
    rows = []
    for column in nodes.columns:
        series = nodes[column]
        if not pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series):
            continue
        value = series.loc[i]
        if pd.isna(value):
            continue
        measured = series.dropna()
        if len(measured) < UNUSUAL_MIN_MEASURED or measured.nunique() <= 2:
            continue
        category = category_of(column, organism)
        if category in UNUSUAL_SKIP:
            continue
        below = float((measured < value).mean()) + 0.5 * float((measured == value).mean())
        rows.append((column, float(value), below, int(len(measured)), category,
                     _dataset_of(column)[0], int((measured > value).sum()) + 1,
                     int((measured < value).sum()) + 1))
    frame = pd.DataFrame(rows, columns=cols)
    if not len(frame):
        return frame
    frame["distance"] = (frame["percentile"] - 0.5).abs()
    frame = frame.sort_values("distance", ascending=False)
    frame = frame[frame.groupby(frame["dataset"].where(frame["dataset"] != "", frame["column"]))
                  .cumcount() < per_dataset]
    return frame.head(top).drop(columns="distance").reset_index(drop=True)


def _ordinal(n: int) -> str:
    return f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"


def unusual_sentences(frame: pd.DataFrame) -> list:
    """One line per unusual measurement: its rank among the genes measured, and where it is from.

    A rank rather than a rounded percentage: the top gene is "the highest of 8,131", never "higher
    than 100 %", which would be false.
    """
    out = []
    for r in frame.itertuples(index=False):
        high = r.percentile >= 0.5
        rank = r.rank_high if high else r.rank_low
        word = "highest" if high else "lowest"
        place = f"the {word}" if rank == 1 else f"the {_ordinal(rank)} {word}"
        share = rank / r.measured
        top = f", top {share:.1%}" if share >= 0.001 else ""
        where = f" ({escape(str(r.dataset))})" if r.dataset else ""
        out.append(f"<b>{escape(str(r.column))}</b> = {escape(value_text(r.column, r.value))}: "
                   f"{place} of the {r.measured:,} genes measured{top}{where}.")
    return out


def strategy_profile(gene: str, organism: str = "Tg", ledger: pd.DataFrame | None = None) -> pd.DataFrame:
    """Every catalogue strategy, and what it said about this gene when the gene was hidden (64.30).

    One row per strategy: its calibration grade, and for each label the hold-out record covers, the call it
    made with the gene's own answer hidden -- right, wrong or abstained -- at its recorded setting. A
    strategy that answers about structure rather than about one gene says so instead of being left out, so
    all of the catalogue is accounted for.
    """
    from . import calibration, strategies as S, track_record as T
    if ledger is None:
        ledger = T.shipped(organism)
    mine = ledger[(ledger["gene_id"].astype(str) == str(gene))
                  & (ledger["mode"].astype(str) == "together")] if len(ledger) else ledger
    recorded = set(T.supported(organism))
    rows = []
    for strategy in S.catalog():
        try:
            entry = calibration.entry(strategy.key, organism) or {}
            grade = calibration.grade(entry) if entry else "not calibrated"
        except Exception:
            grade = "not calibrated"
        part = mine[mine["strategy"].astype(str) == strategy.key] if len(mine) else mine
        if strategy.key not in recorded:
            calls, status = "", "answers about structure or a list, not about one gene"
        elif not len(part):
            calls, status = "", "no hold-out record for this gene (no label it was tested on)"
        else:
            bits = []
            for r in part.itertuples(index=False):
                if bool(r.abstained):
                    bits.append(f"{r.target}: abstained")
                else:
                    verdict = "right" if r.correct is True or r.correct == 1 else "wrong"
                    bits.append(f"{r.target}: {r.prediction} ({verdict}; truth {r.truth})")
            calls, status = "; ".join(dict.fromkeys(bits)), "held-out calls"
        rows.append(dict(strategy=strategy.name, key=strategy.key, grade=grade, status=status, calls=calls))
    return pd.DataFrame(rows)


def synthesis(reading: GeneReading) -> list:
    """(category, sentence as HTML) for each kind of knowledge the gene has, in `ORDER`.

    Plain sentences from the observed facts only -- a prediction is never phrased as knowledge
    here; it goes to `inferences`. Each sentence names its sources.
    """
    out = []
    facts = reading.facts
    if "place" in facts:
        parts = []
        hosted = _fact(reading, "host_facing")
        if hosted is not None and value_text("host_facing", hosted["value"]) == "yes":
            parts.append("it reaches beyond the parasite and faces the host"
                         + _cite(reading, "host_facing"))
        elif hosted is not None:
            parts.append("it is not recorded as facing the host" + _cite(reading, "host_facing"))
        rows = facts["place"]
        named = {}
        for row in rows.itertuples(index=False):
            text = value_text(row.column, row.value)
            if row.column == "host_facing" or not re.search(r"[A-Za-z]{3}", text) \
                    or text in ("yes", "no", "unknown"):
                continue
            ref = _ref_of(row) or "; ".join(reading.sources_for(row.column, 1))
            label = (_text(getattr(row, "dataset", "")) or _text(row.description) or row.column)[:50]
            if text.lower() in label.lower():
                continue                                  # a method name, not a place
            named.setdefault(text, []).append(escape(label) + (f", {source_links(ref)}" if ref else ""))
        for value, sources in list(named.items())[:3]:
            parts.append(f"<b>{escape(value)}</b> ({'; '.join(dict.fromkeys(sources))})")
        if parts:
            out.append(("place", "Where it is: " + "; ".join(parts) + "."))
    if "function" in facts:
        rows = informative(facts["function"]).head(3)
        if len(rows):
            items = [f"{escape(str(r.description or r.column)[:60])}: "
                     f"<b>{escape(value_text(r.column, r.value))}</b>"
                     for r in rows.itertuples(index=False)]
            out.append(("function", "What it does: " + "; ".join(items) + "."))
    partners = reading.partners
    if partners is not None and len(partners):
        host = partners[partners["side"].astype(str) == "host"]
        para = partners[partners["side"].astype(str) != "host"]
        strong = partners.sort_values("probability", ascending=False).head(3)
        named = "; ".join(f"<b>{escape(str(r.partner))}</b> ({escape(str(r.side))}, "
                          f"{r.probability:.0%}"
                          f"{', ' + source_links(r.reference) if r.reference else ''})"
                          for r in strong.itertuples(index=False))
        out.append(("interaction", f"What it binds: {len(partners)} recorded partners "
                                   f"({len(host)} in the host, {len(para)} in the parasite); "
                                   f"strongest {named}."))
    if "fitness" in facts:
        rows = informative(facts["fitness"])
        if len(rows):
            head = rows.head(1).iloc[0]
            name = _text(head.get("dataset", "")) or _text(head["description"])
            out.append(("fitness", f"{TITLES['fitness']}: <b>{escape(value_text(head['column'], head['value']))}</b> "
                                   f"in {escape(name)}"
                                   f"{' (' + source_links(_text(head.get('reference', ''))) + ')' if _text(head.get('reference', '')) else ''}; "
                                   f"{len(rows)} fitness measurements in all, from {_datasets(rows)}."))
    for category in ("phenotype", "structure", "domain", "function_extra"):
        if category not in facts:
            continue
        rows = informative(facts[category])
        if not len(rows):
            continue
        shown = rows.head(2)
        items = "; ".join(f"{escape((_text(getattr(r, 'dataset', '')) or _text(r.description) or r.column)[:50])}: "
                          f"{escape(str(r.column).replace('_', ' '))} "
                          f"<b>{escape(value_text(r.column, r.value))}</b>"
                          for r in shown.itertuples(index=False))
        more = f" (and {len(rows) - len(shown)} more)" if len(rows) > len(shown) else ""
        out.append((category, f"{TITLES[category]}: {items}{more}."))
    for category in MEASURED_MANY:
        if category not in facts:
            continue
        rows = informative(facts[category])
        if not len(rows):
            continue
        n = len(rows)
        out.append((category, f"{TITLES[category]}: {n} measurement{'s' if n != 1 else ''}, from "
                              f"{_datasets(rows)}."))
    order = [c for c, _ in ORDER]
    out.sort(key=lambda pair: order.index(pair[0]) if pair[0] in order else 99)
    return out


def inferences(reading: GeneReading) -> list:
    """What is not known, and what is predicted, each as an HTML sentence with its uncertainty."""
    out = []
    for _column, phrase in reading.unknown:
        out.append(f"Not known: {escape(phrase)}.")
    predicted = reading.predicted
    if predicted is not None and len(predicted):
        for row in predicted.head(8).itertuples(index=False):
            p = "" if row.probability != row.probability else f", {row.probability:.0%} certain"
            out.append(f"Predicted, not measured: {escape(str(row.description or row.column)[:70])} "
                       f"= <b>{escape(value_text(row.column, row.value))}</b>"
                       f" ({escape(str(row.source_tier))}{p}).")
    return out
