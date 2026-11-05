"""A file of genes dropped on the window, and an honest account of what was in it.

Instruction 74: a gene list is one of the three things ParaScape can be asked about -- a block of
rows of the matrix -- and dropping a file is simply how one comes into being. This module is what a drop produces -- the file's
own facts, an identifier resolution report, and the extra per-gene columns -- with no Qt, so the
whole of it is testable without a window. `parascape.home_panel` draws it.

**Nothing here is a new reader.** `parascape.importer` already reads `.csv .tsv .txt .tab .xlsx
.xls .parquet`, sniffs the separator rather than trusting the extension (half the supplementary
tables in this field are `.txt` that are really TSV, and a `.csv` exported from Excel in a European
locale is semicolon-separated), lists a workbook's sheets, and describes each numeric column's
range. Its rule governs here too, and this module does not work around it:

    **THE RANGE DECIDES, NOT THE FILENAME.**

So an extra column arrives with the quantification its *numbers* suggest and the reason attached,
for the user to confirm -- never with a type inferred from what the column is called.

Two things this adds, because neither existed and both were being lost:

* **A resolution report in which nothing is silently dropped.** `tuning.import_table` keeps the rows
  that matched and says how many; a user who dropped 300 identifiers and got 250 genes has no way to
  see which 50 went and why. Here every input value comes back as resolved, unresolved or ambiguous,
  with its reason, and the unresolved ones are *listed*. A list that quietly shrinks is how a
  published screen lost every one of its rows once already.
* **A refusal that says what is accepted.** An unsupported extension raises `UnsupportedFormat`
  carrying `importer.READABLE`, rather than failing somewhere inside a parser.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import pandas as pd

from . import importer as IMP
from .tuning import IMPORT_RX, suggest_gene_column, suggest_repair

#: What may be dropped, from the reader itself rather than restated here.
READABLE = IMP.READABLE

#: The separators the reader considers, so the report and the read agree.
SEPARATORS = IMP.SEPARATORS


class UnsupportedFormat(ValueError):
    """A dropped file in a format no reader here accepts, carrying the list of those that are.

    A format that is simply not read must say so with the accepted list, because the alternative --
    an exception from inside a parser -- reads as a broken file rather than an unsupported one.
    """

    def __init__(self, path: str):
        self.path = path
        self.accepted = tuple(READABLE)
        super().__init__(f"{os.path.basename(path)}: ParaScape does not read "
                         f"{os.path.splitext(path)[1] or 'a file with no extension'}. "
                         f"Accepted: {', '.join(READABLE)}")


def accepts(path: str) -> bool:
    """Whether this file's extension is one of the readable ones."""
    return os.path.splitext(str(path))[1].lower() in READABLE


def separator_of(path: str) -> str:
    """The separator the reader will use, as a label: "" for Excel, parquet and a one-column file.

    Reported from `importer.sniff_separator`, the same call `importer.read_any` makes, so what the
    user is shown is what was actually used. A wrongly sniffed file gives one column with everything
    in it, which looks like a file with no data rather than a file read badly -- showing the
    separator is how a user can tell those two apart.
    """
    return IMP.sniff_separator(str(path))


@dataclass(frozen=True)
class Resolution:
    """One input identifier and what became of it. `state` is resolved, unresolved or ambiguous."""

    raw: str
    gene_id: str = ""
    state: str = "resolved"
    reason: str = ""
    candidates: tuple = ()


@dataclass
class Report:
    """What happened to every identifier in the dropped file, nothing dropped silently.

    `rows` holds one `Resolution` per input row in file order, so the user sees their own file back.
    """

    column: str = ""
    rows: list = field(default_factory=list)
    repair: str = ""
    checked_against_table: bool = True

    @property
    def resolved(self) -> list:
        """The resolutions that reached exactly one gene in the table."""
        return [r for r in self.rows if r.state == "resolved"]

    @property
    def unresolved(self) -> list:
        """The resolutions that reached no gene, each with the reason it did not."""
        return [r for r in self.rows if r.state == "unresolved"]

    @property
    def ambiguous(self) -> list:
        """The input values that named more than one gene, which are never guessed at."""
        return [r for r in self.rows if r.state == "ambiguous"]

    @property
    def genes(self) -> list:
        """The distinct gene ids the list resolved to, in file order."""
        out, seen = [], set()
        for r in self.resolved:
            if r.gene_id not in seen:
                seen.add(r.gene_id)
                out.append(r.gene_id)
        return out

    @property
    def duplicates(self) -> int:
        """How many resolved rows named a gene another row had already named."""
        return len(self.resolved) - len(self.genes)

    def summary(self) -> str:
        """One line: how many resolved, how many did not, and that the rest are listed not dropped."""
        parts = [f"{len(self.genes):,} genes from {len(self.resolved):,} of {len(self.rows):,} rows"]
        if self.duplicates:
            parts.append(f"{self.duplicates:,} repeated identifiers")
        if self.ambiguous:
            parts.append(f"{len(self.ambiguous):,} ambiguous")
        if self.unresolved:
            parts.append(f"{len(self.unresolved):,} unresolved and listed below, not discarded")
        if self.repair:
            parts.append(f"identifiers reshaped with {self.repair}")
        if not self.checked_against_table:
            parts.append("no node table to check against, so only the identifier form was verified")
        return " · ".join(parts)


def resolution_report(values, known_ids=None, resolve=None, column: str = "") -> Report:
    """Resolve every value of an identifier column, keeping the ones that did not.

    `known_ids` is the node table's own `gene_id` column: a well-formed accession that the table
    does not hold is unresolved **with that reason**, which is a different fact from a value that
    was never an accession at all. `resolve` is an identity-layer callable, used when one is
    available so that previous and strain accessions reach their current gene.
    """
    series = pd.Series(list(values)).astype(str)
    known = set(map(str, known_ids)) if known_ids is not None else None
    # The repair regex is suggested from the values that need one, and applied only to those: a file
    # mixing `TGME49_208830` with `TgME49.233460` is common, and judging the whole column by whether
    # ANY value is already canonical loses every malformed row in it.
    malformed = series[~series.str.contains(IMPORT_RX, regex=True, na=False)]
    repair = suggest_repair(malformed) if len(malformed) else None
    pattern, replacement = repair if repair else (None, None)
    report = Report(column=column, repair=pattern or "", checked_against_table=known is not None)
    for raw in series:
        text = raw
        if pattern and not re.search(IMPORT_RX.pattern, raw, flags=re.I):
            text = re.sub(pattern, replacement, raw, flags=re.I)
        found = []
        for match in re.finditer(IMPORT_RX.pattern, text, flags=re.I):
            candidate = match.group(0).upper()
            if resolve is not None:
                candidate = resolve(candidate) or candidate
            if candidate not in found:
                found.append(candidate)
        if not found:
            report.rows.append(Resolution(raw, state="unresolved",
                                          reason="no gene identifier in this value"))
            continue
        if len(found) > 1:
            report.rows.append(Resolution(raw, state="ambiguous", candidates=tuple(found),
                                          reason="names more than one gene; which one is not guessed"))
            continue
        gene = found[0]
        if known is not None and gene not in known:
            report.rows.append(Resolution(raw, gene_id=gene, state="unresolved",
                                          reason="a well-formed accession this table does not hold"))
            continue
        report.rows.append(Resolution(raw, gene_id=gene))
    return report


@dataclass
class Extra:
    """One per-gene column found beside the identifiers, and what it would be treated as.

    `treated_as` and `why` come from `importer.suggest_quantification`, which reads the numbers. The
    user confirms it; nothing is applied from this record alone.
    """

    column: str
    n: int
    minimum: float
    median: float
    maximum: float
    missing: float
    treated_as: str
    why: str

    def line(self) -> str:
        """The row a user reads: the range, then the treatment proposed for it, with the reason."""
        return (f"{self.column}: {self.n:,} values, {self.minimum:.3g} to {self.maximum:.3g} "
                f"(median {self.median:.3g}, {self.missing:.0%} missing) — would be treated as "
                f"{self.treated_as}, because {self.why}")


def extra_columns(table: pd.DataFrame, gene_column: str = "") -> list:
    """Every numeric column beside the identifier, with the quantification its numbers suggest."""
    described = IMP.describe_columns(table, IMP.numeric_columns(table, gene_column))
    out = []
    for row in described.itertuples(index=False):
        quant, why = IMP.suggest_quantification(table, [row.column])
        out.append(Extra(str(row.column), int(row.n), float(row.min), float(row.median),
                         float(row.max), float(row.missing), quant, why))
    return out


@dataclass
class Drop:
    """A dropped file, read: its own facts, its resolution report and its extra columns.

    Every one of the file's facts is changeable by the user, which is why they are fields and not
    only a message: the sheet, the identifier column and the treatment of each extra column are
    decisions, and a dropped file that silently decided them is the thing this project avoids.
    """

    path: str
    sheet: str = ""
    sheets: tuple = ()
    separator: str = ""
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    gene_column: str = ""
    report: Report = field(default_factory=Report)
    extras: list = field(default_factory=list)

    @property
    def name(self) -> str:
        """The file's name, which is what the user recognises it by."""
        return os.path.basename(self.path)

    @property
    def genes(self) -> list:
        """The gene ids this file resolved to."""
        return self.report.genes

    def headline(self) -> str:
        """The line naming the file, its sheet, its separator and the column taken as identifier."""
        bits = [self.name]
        if self.sheet:
            bits.append(f"sheet {self.sheet}")
        if self.separator:
            bits.append("separator " + ("tab" if self.separator == "\t" else repr(self.separator)))
        bits.append(f"identifier column {self.gene_column!r}")
        return " · ".join(bits)


def read_drop(path: str, sheet=None, known_ids=None, gene_column: str = "", resolve=None) -> Drop:
    """Read a dropped file into a `Drop`, or raise `UnsupportedFormat` saying what is accepted.

    The reading is `importer.read_any`, so the separator is sniffed and a workbook's sheet is
    honoured. What is added is the report: which identifiers resolved, which did not and why, and
    which extra columns came along.
    """
    path = str(path)
    if not accepts(path):
        raise UnsupportedFormat(path)
    sheets = tuple(IMP.sheets_in(path))
    chosen = sheet if sheet is not None else (sheets[0] if sheets else 0)
    table = IMP.read_any(path, sheet=chosen)
    if not isinstance(table, pd.DataFrame) or table.empty:
        raise ValueError(f"{os.path.basename(path)} holds no rows to read"
                         + (f" on sheet {chosen!r}" if sheets else ""))
    column = gene_column or suggest_gene_column(table)[0] or str(table.columns[0])
    report = resolution_report(table[column], known_ids=known_ids, resolve=resolve, column=column)
    return Drop(path=path, sheet=str(chosen) if sheets else "", sheets=sheets,
                separator=separator_of(path), table=table, gene_column=column, report=report,
                extras=extra_columns(table, column))


def annotated(drop: Drop, nodes: pd.DataFrame, columns=None) -> pd.DataFrame:
    """The resolved genes joined to the node table, unresolved rows kept and marked.

    The unresolved rows stay in the frame with an empty `gene_id` and their reason, because a table
    that silently holds only what worked is a table nobody can audit.
    """
    rows = [{"input": r.raw, "gene_id": r.gene_id, "state": r.state, "reason": r.reason}
            for r in drop.report.rows]
    frame = pd.DataFrame(rows, columns=["input", "gene_id", "state", "reason"])
    if nodes is None or "gene_id" not in getattr(nodes, "columns", ()):
        return frame
    keep = ["gene_id"] + [c for c in (columns or ()) if c in nodes.columns]
    return frame.merge(nodes[keep].astype({"gene_id": str}), on="gene_id", how="left")
