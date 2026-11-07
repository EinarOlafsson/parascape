"""A methods paragraph, a numbered reference list and BibTeX for any run (instruction 60, change 2).

Every citation comes from the dataset registry (`datasets.py`) by way of the columns a run touched:
`datasets.provenance` for each column, then the datasets of the columns it was declared to be derived
from, so a derived score cites the measurements behind it. Nothing is typed from memory. A dataset
whose publication is not recorded is listed as **not yet confirmed** and kept out of the BibTeX, the
same rule `datasets.unresolved()` exists for. The BibTeX carries only what the registry holds -- the
citation as recorded, the PMID, the accession and the URL -- and never a title, journal or DOI that
would have to be guessed.

    text = methods_paragraph(columns, "Tg", strategy="feature_knn", settings={"k": 15}, seed=0)
    refs = reference_list(sources(columns, "Tg"))
    bib = bibtex(sources(columns, "Tg"))

No window here; every rule is testable.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from . import datasets as D

#: Words for an organism code, for the sentence.
ORGANISM_WORDS = {"Tg": "Toxoplasma gondii (ME49)", "Pf": "Plasmodium falciparum (3D7)"}


@dataclass
class Source:
    """One dataset a run drew on, with the columns that brought it in."""

    dataset: object
    columns: list = field(default_factory=list)
    via: list = field(default_factory=list)     # derived columns that led here, when not used directly

    @property
    def confirmed(self) -> bool:
        """True when the registry records the originating publication."""
        return bool(getattr(self.dataset, "citation", None))

    @property
    def computed(self) -> bool:
        """True for a layer ParaScape computed from other columns: it cites its inputs, not a paper."""
        return bool(getattr(self.dataset, "derived_from", ())) and not self.confirmed


@dataclass
class Sources:
    """Every dataset behind a set of columns, in first-use order, and what could not be traced."""

    found: list = field(default_factory=list)
    unregistered: list = field(default_factory=list)

    def numbered(self) -> dict:
        """dataset key -> reference number, confirmed sources only, in order."""
        out, n = {}, 0
        for s in self.found:
            if s.confirmed:
                n += 1
                out[s.dataset.key] = n
        return out


def sources(columns, organism: str | None = None) -> Sources:
    """Resolve columns to the datasets that produced them, following declared derivations."""
    out = Sources()
    index: dict = {}

    def add(dataset, column, via=None):
        s = index.get(dataset.key)
        if s is None:
            s = index[dataset.key] = Source(dataset)
            out.found.append(s)
        if via is None:
            if column not in s.columns:
                s.columns.append(column)
        elif via not in s.via:
            s.via.append(via)

    for column in dict.fromkeys(str(c) for c in columns):
        dataset = D.provenance(column, organism)
        if dataset is None:
            out.unregistered.append(column)
            continue
        add(dataset, column)
        seen, queue = {column}, list(D.derived_sources(column, organism))
        while queue:
            parent = queue.pop(0)
            if parent in seen:
                continue
            seen.add(parent)
            origin = D.provenance(parent, organism)
            if origin is not None and origin.key != dataset.key:
                add(origin, parent, via=column)
            queue.extend(D.derived_sources(parent, organism))
    return out


def _settings_text(settings) -> str:
    if not settings:
        return "default settings"
    return ", ".join(f"{k} = {v}" for k, v in sorted(dict(settings).items()))


def methods_paragraph(columns, organism: str, *, strategy: str = "", settings=None, seed=None,
                      version: str = "", n_genes: int | None = None, held_out: str = "") -> str:
    """One paragraph a methods section can start from, with numbered citations."""
    src = sources(columns, organism)
    numbers = src.numbered()
    if not version:
        try:
            from . import __version__ as version
        except Exception:
            version = ""
    parts = []
    for s in src.found:
        n = numbers.get(s.dataset.key)
        mark = (f" [{n}]" if n else " (computed in ParaScape from the sources listed)" if s.computed
                else " (publication not yet confirmed)")
        parts.append(f"{s.dataset.name}{mark}")
    who = ORGANISM_WORDS.get(organism, organism)
    sentences = [f"Analyses were run in ParaScape{' ' + str(version) if version else ''} on "
                 f"{f'{n_genes:,} ' if n_genes else ''}{who} genes."]
    if parts:
        listing = parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + " and " + parts[-1]
        sentences.append(f"The evidence used was {listing}, {len(list(dict.fromkeys(map(str, columns))))} "
                         f"columns in all.")
    if src.unregistered:
        k = len(src.unregistered)
        sentences.append(f"{k} further column{'s have' if k != 1 else ' has'} no recorded source "
                         f"({', '.join(src.unregistered[:5])}{'...' if len(src.unregistered) > 5 else ''}).")
    if strategy:
        sentences.append(f"The strategy was {strategy} with {_settings_text(settings)}"
                         + (f", random seed {seed}" if seed is not None else "") + ".")
    if held_out:
        sentences.append(f"Evidence that restates the label {held_out} was held out before any "
                         "model was fitted.")
    return " ".join(sentences)


def reference_list(src: Sources) -> list:
    """Numbered references for the confirmed sources, then the unconfirmed ones, flagged."""
    numbers = src.numbered()
    lines = []
    for s in src.found:
        d = s.dataset
        bits = [d.citation] if d.citation else []
        if d.pmid:
            bits.append(f"PMID {d.pmid}")
        if d.accession:
            bits.append(d.accession)
        n = numbers.get(d.key)
        if n:
            lines.append(f"[{n}] " + "; ".join(bits))
    for s in src.found:
        if not s.confirmed and not s.computed:
            d = s.dataset
            extra = "; ".join(x for x in (d.accession, f"PMID {d.pmid}" if d.pmid else "") if x)
            lines.append(f"[--] {d.name}: publication not yet confirmed -- confirm before citing"
                         + (f" ({extra})" if extra else ""))
    return lines


def _bib_key(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", key).strip("_") or "dataset"


def _bib_escape(text: str) -> str:
    return str(text).replace("\\", "\\\\").replace("{", "\\{").replace("}", "\\}")


def bibtex(src: Sources) -> str:
    """BibTeX for the confirmed sources, holding only the fields the registry records."""
    entries = []
    for s in src.found:
        d = s.dataset
        if s.computed:
            continue
        if not s.confirmed:
            entries.append(f"% {d.key}: publication not yet confirmed -- not cited")
            continue
        fields = [f"  note = {{{_bib_escape(d.citation)}}}"]
        if d.pmid:
            fields.append(f"  pmid = {{{_bib_escape(d.pmid)}}}")
        if d.accession:
            fields.append(f"  howpublished = {{{_bib_escape(d.accession)}}}")
        if d.url:
            fields.append(f"  url = {{{_bib_escape(d.url)}}}")
        entries.append(f"@misc{{{_bib_key(d.key)},\n" + ",\n".join(fields) + "\n}")
    return "\n\n".join(entries) + ("\n" if entries else "")


def full_text(columns, organism: str, **kwargs) -> str:
    """The paragraph, the reference list and the BibTeX, ready to paste."""
    src = sources(columns, organism)
    refs = reference_list(src)
    return (methods_paragraph(columns, organism, **kwargs) + "\n\nReferences\n" + "\n".join(refs)
            + "\n\nBibTeX\n" + bibtex(src))
