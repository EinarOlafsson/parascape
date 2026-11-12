"""A LIST of genes, read: what the genes of it share, what is measured, what is not.

The pure half of the gene list card (:mod:`parascape.list_module`), written apart from Qt so
every sentence the card shows can be tested on a made-up table. It is the list counterpart of
:mod:`parascape.gene_reading`, and it answers the four questions instruction 98 names for a
list, in this order:

1. **What do the genes of the list share that genes outside it do not?**
   :func:`shared_labels`. Every label a gene can carry is counted inside the list and outside
   it, and the ones the list holds far more of than the rest of the genome are reported with
   how surprising that is -- how often a list of the same size, picked at random, would hold
   that many or more. A share with nothing to compare it against is not a finding: 20 of 40
   genes in the nucleus means nothing until you know that 20 of every 40 genes are.
2. **Where are the genes of the list in the parasite?** :func:`spread`, which is
   :func:`parascape.gene_reading.placement` for many genes at once, so the schematic Einar
   asked the restructure to keep can be drawn for a list as well as for one gene.
3. **How much of the list is measured at all?** :func:`coverage`, per kind of knowledge, with
   how much of the rest of the genome is measured for the same kind -- so a list that is
   badly covered can be told from a kind of knowledge that nobody has measured for anything.
4. **What is unknown for most of the list?** :func:`unknown_for_most`, which is the list of
   places where inference is the only way forward, and is therefore what the card hands to
   Infer.

THE RULES THESE FUNCTIONS ARE WRITTEN TO
========================================
* **Against chance, with n.** No share is ever reported on its own. Every one of them carries
  the share outside the list and the number of genes both rest on.
* **Absence is unknown.** A gene with no value for something is counted as unknown and said
  so in words. It is never counted as a "no", and it is never quietly left out of a
  denominator, which would turn "we did not look" into "it is not there".
* **Whole numbers.** Percentages are whole numbers, as Einar asked.
* **Nothing is dropped silently.** An identifier that names no gene in the table is reported
  (:class:`Resolved`), not discarded, which is the rule `gene_drop` was written to.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from html import escape

from . import gene_reading as GR

#: What a value means when it means "nothing was measured", rather than naming an answer.
#: The same list the parasite drawing greys a gene out for, plus the spellings a table of
#: labels uses for an empty cell.
UNMEASURED = GR.UNPLACED + ("nd", "not measured", "not determined", "-", "?")

#: Kinds of knowledge never reported as something a list shares: what a gene is CALLED and how
#: much it was STUDIED are properties of the literature, not of the biology, and a list pulled
#: out of one paper trivially shares both. The same exclusion `gene_reading.unusual` makes.
SKIP_CATEGORIES = set(GR.UNUSUAL_SKIP) | {"sequence"}

#: A label column with more distinct answers than this is an identifier in disguise -- a
#: product name, an accession, a description -- and "these 40 genes share 40 different
#: products" is not a finding. Read against the number of genes as well, below.
MAX_DISTINCT = 80

#: A column must be measured on at least this many genes before a share outside the list means
#: anything. Below it the comparison is noise, and a noisy comparison reported as a finding is
#: worse than no comparison: it is a wrong one.
MIN_MEASURED = 30

#: How many of the list must carry a label before it is reported. Two, because one gene
#: carrying something is not a thing the LIST shares -- it is a thing one gene has, which the
#: gene card already says.
MIN_IN_LIST = 2

#: Below this share of the list measured, a kind of knowledge counts as "unknown for most of
#: it", which is what inference is for.
MOSTLY_UNKNOWN = 0.5

#: How much more often the list must carry a label than the other genes do before it is
#: reported: half as much again. Added after looking at the card, which offered "no: 28 of 28
#: (100 %) against 95 % of 8,112" as something the list shares -- true, and not a finding.
#: A bar on the SIZE of the difference rather than a list of columns to ignore, because the
#: columns that produced it (whether a transferred label was accepted, whether the scorers
#: agreed) are perfectly good evidence when a list really does differ in them.
MIN_TIMES = 1.5

#: And how rare it must be: a list of the same size picked at random must do it less than
#: once in this many tries. Twenty, so that a card full of things that chance would hand
#: anybody is not read as twenty findings.
MIN_ONE_IN = 20

#: What splits a typed list: commas, semicolons, newlines, tabs and spaces. All of them,
#: because people paste a column out of a spreadsheet as readily as they type three
#: accessions with commas, and asking which one they meant would be a dialog nobody wants.
_SPLIT = re.compile(r"[\s,;|]+")


def split_identifiers(text) -> list:
    """The identifiers in typed text, in the order typed, each once.

    :param text: what was typed or pasted: commas, semicolons, spaces, tabs or newlines
        between the identifiers, in any mixture.
    :returns: the identifiers, stripped, with repeats removed but the first of each kept
        where it was. Order is kept because a list has an order the person chose, and
        sorting it would quietly lose it.
    """
    out, seen = [], set()
    for part in _SPLIT.split(str(text or "")):
        token = part.strip().strip(",;").strip()
        if token and token.lower() not in seen:
            seen.add(token.lower())
            out.append(token)
    return out


def looks_like_a_list(text) -> bool:
    """Whether typed text is several identifiers rather than one thing to search for.

    Used by the Navigator to decide whether what was typed is a LIST, which is the whole
    difference between offering the gene card and offering the gene list card.

    :param text: what was typed.
    """
    return len(split_identifiers(text)) > 1


@dataclass(frozen=True)
class Resolved:
    """What a list of identifiers turned out to be, including what it did not.

    :param genes: the identifiers that name a gene of the table, in the order given.
    :param missing: the identifiers that name none, kept rather than dropped so somebody who
        pasted 300 rows and got 250 genes can see which 50 did not land.
    :param repeated: how many identifiers named a gene an earlier one had already named.
    """

    genes: tuple = ()
    missing: tuple = ()
    repeated: int = 0

    @property
    def count(self) -> int:
        """How many distinct genes of the table the list reached."""
        return len(self.genes)

    @property
    def asked(self) -> int:
        """How many identifiers were handed in, counting the ones that reached nothing."""
        return len(self.genes) + len(self.missing) + self.repeated

    def sentence(self) -> str:
        """One line saying what the list is, and what of it did not land. Plain words."""
        if not self.asked:
            return "No genes yet."
        parts = [f"{self.count:,} gene{'' if self.count == 1 else 's'} of "
                 f"{self.asked:,} identifier{'' if self.asked == 1 else 's'}"]
        if self.repeated:
            parts.append(f"{self.repeated:,} repeated")
        if self.missing:
            shown = ", ".join(self.missing[:4])
            more = f" and {len(self.missing) - 4:,} more" if len(self.missing) > 4 else ""
            parts.append(f"{len(self.missing):,} name no gene in this table ({shown}{more}) "
                         "-- listed, not dropped")
        return "; ".join(parts) + "."


def resolve(nodes, ids) -> Resolved:
    """Match identifiers against the gene table, reporting the ones that match nothing.

    Matched on the identifier exactly, ignoring case, because that is what a list of
    accessions is. A product word reaching several genes is the job of the search field of
    the gene card, not of a list: a list whose entries silently became three genes each
    would not be the list anybody handed in.

    :param nodes: the gene table, with a ``gene_id`` column.
    :param ids: the identifiers, from :func:`split_identifiers` or a file.
    :returns: a :class:`Resolved`.
    """
    wanted = [str(i).strip() for i in (ids or []) if str(i).strip()]
    if nodes is None or "gene_id" not in getattr(nodes, "columns", ()):
        return Resolved(missing=tuple(wanted))
    known = {str(g).lower(): str(g) for g in nodes["gene_id"].astype(str)}
    found, missing, seen, repeated = [], [], set(), 0
    for one in wanted:
        gene = known.get(one.lower())
        # One `seen` for both: an identifier given twice is counted once whether it names a
        # gene or not, so the count of what did not land is a count of DISTINCT things that
        # did not land rather than of how many times they were typed.
        token = gene or one.lower()
        if token in seen:
            repeated += 1
            continue
        seen.add(token)
        if gene is None:
            missing.append(one)
        else:
            found.append(gene)
    return Resolved(genes=tuple(found), missing=tuple(missing), repeated=repeated)


def _rows(nodes, genes):
    """The rows of the table for these genes, in table order; an empty frame when none."""
    ids = {str(g) for g in (genes or [])}
    if nodes is None or "gene_id" not in getattr(nodes, "columns", ()) or not ids:
        return None if nodes is None else nodes.iloc[:0]
    return nodes[nodes["gene_id"].astype(str).isin(ids)]


def _measured(series):
    """A boolean mask of the rows of `series` that hold a measured value.

    One definition, used by every reading here, so "measured" cannot come to mean two
    different things in two sentences on one card. A blank, a NaN and the words that spell
    an empty cell all count as unknown.
    """
    text = series.astype(str).str.strip().str.lower()
    return series.notna() & ~text.isin(UNMEASURED)


# --------------------------------------------------------------- how much is measured
@dataclass(frozen=True)
class Coverage:
    """How much of a list is measured for one kind of knowledge, against the rest of the genes.

    :param category: the information category, e.g. ``"place"``.
    :param title: what it is called in plain words, e.g. ``"Where it is"``.
    :param columns: how many measurements of that kind the table holds at all.
    :param in_list: how many genes of the list have at least one of them measured.
    :param listed: how many genes the list has in the table.
    :param outside_with: how many genes OUTSIDE the list have at least one measured.
    :param outside: how many genes are outside the list.
    """

    category: str
    title: str
    columns: int
    in_list: int
    listed: int
    outside_with: int
    outside: int

    @property
    def percent(self) -> int:
        """What share of the list is measured for this, as a whole number."""
        return int(round(100 * self.in_list / self.listed)) if self.listed else 0

    @property
    def outside_percent(self) -> int:
        """What share of the genes outside the list is measured for it, as a whole number."""
        return int(round(100 * self.outside_with / self.outside)) if self.outside else 0

    @property
    def unknown(self) -> int:
        """How many genes of the list have nothing of this kind measured."""
        return max(0, self.listed - self.in_list)

    @property
    def mostly_unknown(self) -> bool:
        """Whether this is unknown for most of the list, which is what inference is for."""
        return bool(self.listed) and self.in_list < MOSTLY_UNKNOWN * self.listed


def _columns_by_category(nodes, organism=None) -> dict:
    """Category -> the columns of the table that belong to it, skipping names and study effort."""
    out: dict = {}
    if nodes is None:
        return out
    for column in getattr(nodes, "columns", ()):
        name = str(column)
        if name == "gene_id":
            continue
        category = GR.category_of(name, organism)
        if not category or category in SKIP_CATEGORIES:
            continue
        out.setdefault(category, []).append(name)
    return out


def coverage(nodes, genes, organism=None) -> tuple:
    """How much of the list is measured for each kind of knowledge, against the other genes.

    Reported per KIND rather than per column, because "47 of the 130 columns hold something
    for this gene" is a fact about the table and "where they are is measured for 12 of your
    40 genes" is a fact about the list. The kinds and their wording are
    :data:`parascape.gene_reading.ORDER`, so the gene card and the list card ask about a gene
    in the same words.

    :param nodes: the gene table.
    :param genes: the genes of the list.
    :param organism: which organism the columns belong to, for the category lookup.
    :returns: one :class:`Coverage` per kind of knowledge the table holds, in the order of
        :data:`parascape.gene_reading.ORDER`; empty when there is no list or no table.
    """
    rows = _rows(nodes, genes)
    if rows is None or not len(rows):
        return ()
    inside = nodes["gene_id"].astype(str).isin({str(g) for g in genes})
    by_category = _columns_by_category(nodes, organism)
    out = []
    for category, title in GR.ORDER:
        columns = by_category.get(category) or []
        if not columns:
            continue
        any_measured = None
        for column in columns:
            mask = _measured(nodes[column])
            any_measured = mask if any_measured is None else (any_measured | mask)
        out.append(Coverage(category=category, title=title, columns=len(columns),
                            in_list=int((any_measured & inside).sum()),
                            listed=int(inside.sum()),
                            outside_with=int((any_measured & ~inside).sum()),
                            outside=int((~inside).sum())))
    return tuple(out)


def coverage_sentence(cov: Coverage) -> str:
    """One line: how much of the list this is measured for, against the other genes. HTML.

    :param cov: one reading from :func:`coverage`.
    """
    head = (f"<b>{escape(cov.title)}</b> is measured for <b>{cov.in_list:,}</b> of the "
            f"{cov.listed:,} genes of the list ({cov.percent} %), against "
            f"{cov.outside_percent} % of the {cov.outside:,} genes outside it")
    if cov.unknown:
        head += (f". It is <span style='color:#999'>unknown</span> for the other "
                 f"{cov.unknown:,}, which is not the same as there being nothing there")
    return head + "."


def unknown_for_most(coverages) -> tuple:
    """The kinds of knowledge unknown for most of the list, emptiest first.

    This is the list of places where inference is the only way forward, which is why the card
    puts it beside the button that opens Infer. A kind of knowledge nobody has measured for
    any gene is included and says so through `outside_percent`: the gap is then in the data
    rather than in the list, and that is a different problem with a different answer.

    :param coverages: the readings from :func:`coverage`.
    :returns: the ones measured for less than half the list, least measured first.
    """
    poor = [c for c in coverages if c.mostly_unknown]
    return tuple(sorted(poor, key=lambda c: (c.percent, -c.unknown)))


def unknown_sentence(cov: Coverage) -> str:
    """One line saying what is unknown for most of the list, and whether that is unusual. HTML.

    :param cov: one reading from :func:`unknown_for_most`.
    """
    head = (f"<b>{escape(cov.title)}</b> is <span style='color:#999'>unknown</span> for "
            f"<b>{cov.unknown:,}</b> of the {cov.listed:,} genes of the list "
            f"({100 - cov.percent} %)")
    if cov.outside_percent <= cov.percent:
        return (head + f", and for {100 - cov.outside_percent} % of the other genes as well "
                "-- the gap is in the data rather than in your list.")
    return (head + f", against {100 - cov.outside_percent} % of the other genes -- so this "
            "list is measured less than the genome is for it.")


# ------------------------------------------------------------ what the list shares
@dataclass(frozen=True)
class Shared:
    """One label the genes of a list hold more often than the genes outside it.

    Every number here is a count, and every share is reported beside the share outside the
    list, because a share on its own cannot be read: half the list in the nucleus is only a
    finding if half the genome is not.

    :param column: the measurement the label comes from.
    :param value: the label itself.
    :param in_list: how many genes of the list carry it.
    :param list_measured: how many genes of the list are measured for this column at all --
        the denominator, which is NOT the size of the list when part of it is unknown.
    :param outside_with: how many genes outside the list carry it.
    :param outside_measured: how many genes outside the list are measured for the column.
    :param surprise: the chance that a list of the same size, picked at random from the
        measured genes, would hold this many or more. Not shown as a number; turned into
        "once in so many" by :meth:`one_in`.
    :param also: the OTHER measurements that say the same thing about the list, named rather
        than repeated as rows of their own. Five measurements place a gene in this table and
        they mostly agree, so a list of dense-granule genes filled the top of the card with
        five rows all saying "dense granules". That is one finding, and five measurements
        agreeing on it is a fact about the finding -- which is what this field carries. See
        :func:`merge_by_value`.
    """

    column: str
    value: str
    in_list: int
    list_measured: int
    outside_with: int
    outside_measured: int
    surprise: float = 1.0
    also: tuple = ()

    @property
    def list_percent(self) -> int:
        """What share of the measured part of the list carries this label, whole number."""
        return int(round(100 * self.in_list / self.list_measured)) if self.list_measured else 0

    @property
    def outside_percent(self) -> int:
        """What share of the measured genes outside the list carries it, whole number."""
        return (int(round(100 * self.outside_with / self.outside_measured))
                if self.outside_measured else 0)

    @property
    def times(self) -> float:
        """How many times as often the list carries it as the other genes do."""
        if not self.outside_percent:
            return float("inf") if self.list_percent else 0.0
        return self.list_percent / self.outside_percent

    @property
    def one_in(self) -> int:
        """How rarely a list picked at random would hold this many or more: "once in n".

        :returns: n, or 0 when a random list would hold about this many -- which is the
            answer "this is not a surprise", and is said in those words rather than as a 1.
        """
        if self.surprise <= 0:
            return 1_000_000
        n = int(round(1 / self.surprise))
        return n if n >= 2 else 0


def _label_columns(nodes, organism=None) -> list:
    """The columns worth asking "what does the list share" of, in table order.

    A column qualifies when it holds LABELS: text, measured for enough genes to compare
    against (:data:`MIN_MEASURED`), with few enough distinct answers (:data:`MAX_DISTINCT`)
    that an answer is shared rather than unique, and belonging to a KIND OF KNOWLEDGE. A
    product name passes the first two and fails the third, which is the whole reason the
    third exists: "these 40 genes share 40 different products" is not a finding.

    The fourth was added after looking at the card. A column with no kind of knowledge is
    bookkeeping -- whether a transferred label was accepted, whether the scorers agreed --
    and the card reported "no: 28 of 28 (100 %) against 95 %" with "no" as the thing the
    list shares. A reader cannot tell what that is about, and it is not about the biology.

    :param nodes: the gene table.
    :param organism: which organism the columns belong to, for the category lookup.
    """
    if nodes is None:
        return []
    skip = SKIP_CATEGORIES
    out = []
    for column in getattr(nodes, "columns", ()):
        name = str(column)
        if name == "gene_id":
            continue
        category = GR.category_of(name, organism)
        if not category or category in skip:
            continue
        series = nodes[column]
        if series.dtype.kind not in ("O", "U", "S", "b"):
            continue
        try:
            values = series[_measured(series)].astype(str)
        except Exception:                      # a column of arrays holds no labels
            continue
        if len(values) < MIN_MEASURED:
            continue
        distinct = values.nunique()
        if distinct < 2 or distinct > min(MAX_DISTINCT, max(2, len(values) // 2)):
            continue
        out.append(name)
    return out


def _tail(in_list: int, with_label: int, measured: int, picked: int) -> float:
    """The chance a random list of `picked` measured genes holds `in_list` or more labelled ones.

    The hypergeometric tail -- drawing without replacement, which is what picking genes is.
    Named `_tail` and reported as "once in so many" rather than by the name of the
    distribution, because the rule for this application is to say what a number does rather
    than what it is called.

    :param in_list: how many of the list carry the label.
    :param with_label: how many measured genes carry it in all.
    :param measured: how many genes are measured for the column in all.
    :param picked: how many of the list are measured for it.
    :returns: a probability in 0..1; 1.0 when the question is not answerable.
    """
    if min(in_list, with_label, measured, picked) <= 0 or picked > measured:
        return 1.0
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(in_list - 1, measured, with_label, picked))
    except Exception:                          # a reading must not depend on scipy being there
        return 1.0


def shared_labels(nodes, genes, columns=(), min_in_list: int = MIN_IN_LIST,
                  limit: int = 12, organism=None) -> tuple:
    """What the genes of the list share that the genes outside it do not, most surprising first.

    This is the set-level reading the catalogue greyed the List module out for the absence of.
    For every label any gene can carry, it counts the genes of the list that carry it and the
    genes outside the list that carry it, and keeps the ones the list holds more of than
    chance would give it.

    Only labels the list holds MORE of are kept, and only by enough to be worth saying:
    :data:`MIN_TIMES` as often as the other genes, and rarely enough that a random list of
    the same size would do it less than once in :data:`MIN_ONE_IN` tries. A label the list
    holds LESS of is a real fact, but it is a fact about the other genes: "your 40 genes are
    not in the nucleus" is read off the places they ARE in, which the same table shows.

    :param nodes: the gene table.
    :param genes: the genes of the list.
    :param columns: which columns to ask of; empty asks :func:`_label_columns`.
    :param min_in_list: how many of the list must carry a label before it is reported.
    :param limit: how many to return.
    :param organism: which organism the columns belong to, for the category lookup.
    :returns: :class:`Shared` readings, most surprising first.
    """
    rows = _rows(nodes, genes)
    if rows is None or not len(rows):
        return ()
    inside = nodes["gene_id"].astype(str).isin({str(g) for g in genes})
    out = []
    for column in (list(columns) or _label_columns(nodes, organism)):
        if column not in nodes.columns:
            continue
        series = nodes[column]
        measured = _measured(series)
        values = series.astype(str)
        in_mask, out_mask = measured & inside, measured & ~inside
        list_measured, outside_measured = int(in_mask.sum()), int(out_mask.sum())
        if list_measured < min_in_list or outside_measured < MIN_MEASURED:
            continue
        counts = values[in_mask].value_counts()
        for value, n in counts.items():
            n = int(n)
            if n < min_in_list:
                continue
            outside_with = int((values[out_mask] == value).sum())
            # The shares, compared BEFORE the tail is computed: a label the list holds no
            # more of than the rest is not a thing the list shares, and asking how
            # surprising it is would be work done to throw away.
            if n / list_measured < MIN_TIMES * (outside_with / outside_measured or 1e-9):
                continue
            shared = Shared(
                column=str(column), value=str(value), in_list=n,
                list_measured=list_measured, outside_with=outside_with,
                outside_measured=outside_measured,
                surprise=_tail(n, outside_with + n, outside_measured + list_measured,
                               list_measured))
            if shared.one_in < MIN_ONE_IN:
                continue
            out.append(shared)
    out.sort(key=lambda s: (s.surprise, -s.in_list))
    return merge_by_value(out)[:limit]


def merge_by_value(shared) -> tuple:
    """One row per thing the list shares, naming the other measurements that agree.

    Written after looking at the card. A list of dense-granule genes filled the top of it
    with five rows -- "dense granules (lopit map)", "dense granules (compartment)", "dense
    granules (lopit mcmc)"... -- which is ONE finding said five times, and it pushed every
    other thing the list shares off the screen. The strongest of them is kept and the rest
    become :attr:`Shared.also`, so nothing is lost and the agreement is reported as the
    thing it is: corroboration.

    Only measurements that say the SAME label are merged. ``apical_secretory`` and ``dense
    granules`` are different answers from measurements of different coarseness, and merging
    those would be claiming an agreement nobody made.

    :param shared: readings from :func:`shared_labels`, most surprising first.
    :returns: the same readings, one per label, in the order given.
    """
    best: dict = {}
    order: list = []
    for s in shared:
        if s.value not in best:
            best[s.value] = s
            order.append(s.value)
        else:
            kept = best[s.value]
            best[s.value] = Shared(
                column=kept.column, value=kept.value, in_list=kept.in_list,
                list_measured=kept.list_measured, outside_with=kept.outside_with,
                outside_measured=kept.outside_measured, surprise=kept.surprise,
                also=kept.also + (s.column,))
    return tuple(best[value] for value in order)


def shared_sentence(s: Shared) -> str:
    """One line: what share of the list carries a label, against the other genes, and how
    surprising that is. HTML, plain words, whole-number percentages.

    :param s: one reading from :func:`shared_labels`.
    """
    what = escape(str(s.column).replace("_", " "))
    # Through `value_text`, so a boolean column reads "yes" rather than "True": the gene card
    # says yes and no, and one card saying True would be the same fact in two languages.
    head = (f"<b>{escape(GR.value_text(s.column, s.value))}</b> ({what}): "
            f"<b>{s.in_list:,}</b> of the "
            f"{s.list_measured:,} genes of the list measured for it ({s.list_percent} %), "
            f"against {s.outside_percent} % of the {s.outside_measured:,} other measured genes")
    one_in = s.one_in
    if not one_in:
        tail = " -- about what a list of this size picked at random would hold."
    elif one_in >= 1_000_000:
        tail = (" -- a list of this size picked at random would hold this many or more "
                "less than once in a million tries.")
    else:
        tail = (f" -- a list of this size picked at random would hold this many or more "
                f"about once in {one_in:,} tries.")
    if s.also:
        tail += (" The same is said by "
                 + escape(", ".join(c.replace("_", " ") for c in s.also)) + ".")
    return head + tail


def surprise_words(s: Shared) -> str:
    """How surprising a shared label is, in a few plain words fit for a table cell.

    The same judgement :func:`shared_sentence` makes, said shorter, so the cell and the
    sentence in the tooltip beside it cannot disagree.

    :param s: one reading from :func:`shared_labels`.
    """
    one_in = s.one_in
    if not one_in:
        return "about what chance would give"
    if one_in >= 1_000_000:
        return "less than once in a million random lists"
    return f"once in {one_in:,} random lists"


def shared_lead(shared, resolved: Resolved) -> str:
    """The one line above the shared labels, saying what was compared with what. HTML.

    :param shared: the readings from :func:`shared_labels`.
    :param resolved: the list itself, from :func:`resolve`.
    """
    if not resolved.count:
        return ("Type or choose a list of genes on the left. This reads what they share that "
                "other genes do not.")
    if not shared:
        return (f"Nothing is held by these {resolved.count:,} genes much more often than by "
                "the rest of them -- or too little of the list is measured to tell. That is "
                "an answer, not a failure: the list is not unusual in anything recorded here.")
    return (f"The {resolved.count:,} genes of the list, read against every other gene in this "
            f"table. {len(shared):,} label{'' if len(shared) == 1 else 's'} "
            f"{'is' if len(shared) == 1 else 'are'} held at least "
            f"{MIN_TIMES:g} times as often by the list as by the rest, and rarely enough "
            f"that a random list of {resolved.count:,} genes would do it less than once in "
            f"{MIN_ONE_IN:,} tries. Most surprising first.")


# ------------------------------------------------ where the list sits in the parasite
@dataclass(frozen=True)
class Spread:
    """Where the genes of a list sit in the parasite, by one placement measurement.

    The list counterpart of :class:`parascape.gene_reading.Placement`: the same measurement,
    the same idea of unknown, counted over many genes instead of one, so the parasite
    schematic Einar asked the restructure to keep can be drawn for a list.

    :param column: the placement measurement being read.
    :param counts: ``((compartment, how many of the list), ...)``, commonest first.
    :param measured: how many genes of the list have a measured compartment.
    :param listed: how many genes the list has in the table.
    :param outside: ``{compartment: how many genes outside the list}``.
    :param outside_measured: how many genes outside the list have one.
    """

    column: str = ""
    counts: tuple = ()
    measured: int = 0
    listed: int = 0
    outside: dict = None
    outside_measured: int = 0

    @property
    def commonest(self) -> str:
        """The compartment most of the measured part of the list is in, or ``""``."""
        return self.counts[0][0] if self.counts else ""

    @property
    def unknown(self) -> int:
        """How many genes of the list have no measured compartment."""
        return max(0, self.listed - self.measured)

    @property
    def compartments(self) -> tuple:
        """The compartments the list is found in, commonest first."""
        return tuple(name for name, _n in self.counts)

    def count_of(self, compartment: str) -> int:
        """How many genes of the list are in one compartment."""
        return next((n for name, n in self.counts if name == compartment), 0)

    def percent_of(self, compartment: str) -> int:
        """What share of the MEASURED part of the list is in one compartment, whole number."""
        return (int(round(100 * self.count_of(compartment) / self.measured))
                if self.measured else 0)

    def outside_percent_of(self, compartment: str) -> int:
        """What share of the measured genes outside the list is in it, whole number."""
        if not self.outside_measured:
            return 0
        return int(round(100 * (self.outside or {}).get(compartment, 0) / self.outside_measured))


def spread(nodes, genes, column: str) -> Spread:
    """Where the genes of a list are placed, and where the other genes are.

    :param nodes: the gene table.
    :param genes: the genes of the list.
    :param column: a placement column, from
        :func:`parascape.gene_reading.placement_columns`.
    :returns: a :class:`Spread`; ``measured`` is 0 when nothing of the list is placed.
    """
    if nodes is None or not column or column not in getattr(nodes, "columns", ()):
        return Spread(column=str(column or ""), outside={})
    inside = nodes["gene_id"].astype(str).isin({str(g) for g in (genes or [])})
    values = nodes[column].astype(str)
    placed = ~values.str.strip().str.lower().isin(GR.UNPLACED) & nodes[column].notna()
    mine = values[placed & inside].value_counts()
    theirs = values[placed & ~inside].value_counts()
    return Spread(column=str(column),
                  counts=tuple((str(k), int(v)) for k, v in mine.items()),
                  measured=int((placed & inside).sum()), listed=int(inside.sum()),
                  outside={str(k): int(v) for k, v in theirs.items()},
                  outside_measured=int((placed & ~inside).sum()))


def spread_sentence(sp: Spread, compartment: str = "") -> str:
    """One line for the drawing: how much of the list is in the organelle shown. HTML.

    Unknown is said in words. A grey parasite with no caption reads as "these genes are
    nowhere", which is a claim nobody made -- the same rule the gene card is written to.

    :param sp: the reading from :func:`spread`.
    :param compartment: the organelle coloured in; empty uses the commonest.
    """
    if not sp.listed:
        return "Type or choose a list on the left to see where in the cell the genes of it are."
    if not sp.measured:
        return ("Where these genes are in the cell is <b>unknown</b> &mdash; none of them has a "
                "measured compartment, which is not the same as them being nowhere. The "
                "parasite is drawn grey.")
    name = compartment or sp.commonest
    mine, theirs = sp.percent_of(name), sp.outside_percent_of(name)
    more = ("more often than" if mine > theirs else
            "less often than" if mine < theirs else "as often as")
    return (f"<b>{sp.count_of(name):,}</b> of the {sp.measured:,} placed genes of the list "
            f"({mine} %) are in the <b>{escape(name)}</b>, {more} the other genes "
            f"({theirs} %). The organelle is filled in the colour that compartment has on "
            f"the map ({escape(str(sp.column).replace('_', ' '))}).")


def spread_trust(sp: Spread) -> str:
    """How far to trust the drawing: how much of the list is placed at all, and by what. HTML.

    :param sp: the reading from :func:`spread`.
    """
    if not sp.listed:
        return ""
    head = (f"Where a gene is, is measured for <b>{sp.measured:,}</b> of the {sp.listed:,} "
            f"genes of the list")
    if sp.listed:
        head += f" ({int(round(100 * sp.measured / sp.listed))} %)"
    if sp.unknown:
        head += (f"; it is <span style='color:#999'>unknown</span> for the other "
                 f"{sp.unknown:,}, which are not drawn")
    if not sp.counts:
        return head + "."
    spread_over = len(sp.counts)
    return (head + f". The placed ones are spread over {spread_over:,} "
            f"compartment{'' if spread_over == 1 else 's'}, the commonest being "
            f"<b>{escape(sp.commonest)}</b> with {sp.percent_of(sp.commonest)} % of them.")
