"""ParaScape's modules, as the home screen's tiles.

This is the catalogue the home screen draws, and the one place that decides what
ParaScape *is* to somebody opening it. It replaces the Gene/List/Label node model
of instruction 74: that model was correct about the shape of the data -- a row, a
block of rows, a column -- and wrong as a front door, because nobody arrives
wanting "a column".

The bands are what the person came to do, not how the data is arranged:

* **Look up** -- I have a gene and I want to know about it.
* **Explore** -- I want to see the whole thing and find my way around.
* **Ask** -- I have a question and I want an answer with an uncertainty.
* **Bring your own** -- I have a screen, a list, or a table of my own.
* **Data** -- what is in here, where it came from, and what is running.

Each module declares what it needs. A module whose view or data is absent is
drawn greyed **with its reason** rather than hidden, so the home screen is an
honest inventory rather than a menu of promises. `availability` is computed at
open time, never hardcoded, which is why `Fill` greys itself out: it asks the
harness whether any question has cleared its bar, and today none has.

No Qt here, so the catalogue can be tested, listed and audited without a display.

ONE TABLE, AND WHY IT IS THIS ONE (instruction 98, WP3)
=======================================================
This file is the single declared catalogue: key, name, plain-words note,
section, icon, what opens it, how it is built, and -- through
:mod:`parascape.ports` -- what subject it takes and what it hands on.

It was extended rather than replaced by a new `module_catalogue.py`, and the
reason is the bug spaCR records as #16i: Home grouped by one table and the tabs
by another, so a section could gain a tab it had no tiles for. ParaScape had
the same shape of problem in four pieces -- this table, the icon map in
`navigator.py`, the shell routes in `landing_install.open_module`, and the tool
specs in `dock_module.SPECS`, which repeated the title and the description of
ten modules in different words. Adding a fifth table called "the catalogue"
would have made it five. So the fields the other tables held moved HERE, and
those tables are now computed from this one:

* ``dock_module.SPECS`` is built from the rows with ``builder == "tool"``;
* ``landing_install.open_module`` reads ``builder`` instead of a dict;
* the dock and the home bands already read :data:`MODULES`, and now read it
  through :func:`section_members` so both group the same way.

The ports are NOT repeated here. :mod:`parascape.ports` declares them and this
module reads them, never the other way round: ports imports nothing out of
ParaScape, so it stays testable and cheap, and `tests/test_ports.py` holds that
line.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import ports as PORTS


@dataclass(frozen=True)
class Module:
    """One tile on the home screen.

    :param key: stable identifier; what the home screen reports when clicked and
        what a saved layout or a search hit refers to.
    :param name: what the tile says. One or two words.
    :param note: one line under the name. Says what you get, not what it is.
    :param band: which group of tiles it belongs to.
    :param opens: where the module goes: a dock as `dock`, a tab in one as
        `dock:tab`, or a window method as `action:method` (the dataset browser and
        export are windows and dialogs, not docks). Empty means the view does not
        exist yet, and the module is unavailable for that reason alone.

        Every destination was checked against a live window, after the first
        version of this catalogue -- written from memory -- pointed six modules at
        tabs that did not exist. `tests/test_modules_in_window.py` keeps it checked.
    :param needs: what must be present for it to open, as free text used in the
        greyed-out reason. Empty means it only needs the app.
    :param icon: the piece of spaCR artwork this module is drawn with, named as
        spaCR names the file (``profiler``, ``umap``), without ``.png``.

        THE NAME OF THE ARTWORK, not a word for the module. It used to be the
        latter -- ``icon="gene"`` beside ``name="Gene"`` -- which said nothing
        the key did not already say, while the fact it was pretending to hold
        lived in a second table, ``iconset.MODULE_ICON``. Two tables holding one
        fact is the bug spaCR records as #16i, and instruction 98 named this
        pair as one of the two still open. :data:`parascape.iconset.MODULE_ICON`
        is now built from this column, so a module gets its symbol by having a
        row here and in no other way.
    :param instruction: what to DO first, in one line, shown under the name in
        the module shell. It lived in `dock_module.SPECS` for ten modules and
        nowhere for the rest; it belongs beside the name it goes under.
    :param builder: how the module is BUILT, which is the part that used to be
        a dict in `landing_install` and a second list in `dock_module`:

        * ``"shell:<module>"`` -- built by ``<module>.build(window)``, a module
          already through the mould. Several keys may name one builder:
          Partners is the Network module read by binding evidence, and
          Questions is the Infer module opened on its catalogue.
        * ``"tool"`` -- built by `dock_module` around the dock named in
          ``opens``: the same shell, borrowing a panel that already exists.
        * ``""`` -- not built at all. It opens whatever ``opens`` names, or
          nothing when ``opens`` is empty too.
    """

    key: str
    name: str
    note: str
    band: str
    opens: str = ""
    needs: tuple = ()
    icon: str = ""
    instruction: str = ""
    builder: str = ""

    @property
    def section(self) -> str:
        """The section this module is filed under.

        The same string as :attr:`band`. spaCR calls it a section and ParaScape
        has always called it a band; both names are used in the instructions, so
        both work here rather than one of them being a rename nobody finds.
        """
        return self.band

    @property
    def ports(self):
        """The subjects this module takes and hands on.

        :returns: the :class:`parascape.ports.ModulePorts` declared for this key.
        :raises parascape.ports.UnknownModule: when nothing is declared, which
            `tests/test_module_catalogue.py` makes impossible for a row that
            ships.
        """
        return PORTS.module_ports(self.key)

    @property
    def consumes(self) -> tuple:
        """The subjects this module takes, as kinds, in declared order."""
        return tuple(p.kind for p in self.ports.consumes)

    @property
    def produces(self) -> tuple:
        """The subjects this module hands on, as kinds, in declared order."""
        return tuple(p.kind for p in self.ports.produces)


#: The bands, in the order they are drawn, each with the question it answers.
BANDS = (
    ("Look up", "You have a gene."),
    ("Explore", "You want to see the whole thing."),
    ("Ask", "You have a question and want an answer with an uncertainty on it."),
    ("Bring your own", "You have a screen, a list or a table of your own."),
    ("Data", "What is in here, where it came from, and what is running."),
)

#: Every module. Order within a band is the order drawn.
MODULES = (
    # ---------------------------------------------------------------- Look up
    Module("gene", "Gene", "everything recorded about one gene, with its sources",
           "Look up", opens="evidence", icon="profiler",
           instruction="Name a gene, or step through a list one gene at a time.",
           builder="shell:gene_module"),
    Module("structure", "Structure", "the fold, its confidence, and what is known of it",
           "Look up", needs=("a structure view; AF3 summaries are on the evidence of a gene for now",),
           icon="layer_viewer",
           instruction="Name a gene to see its fold and how far the fold is to be trusted."),
    Module("partners", "Partners", "what it binds, and how strongly that is evidenced",
           "Look up", opens="evidence", icon="graph_builder",
           instruction="Name a gene; the links drawn are the physical ones, with their evidence.",
           builder="shell:network_module"),
    Module("expression", "Expression", "where and when its RNA and protein are measured",
           "Look up", opens="compare:Gene expression", icon="pca",
           instruction="Choose a gene; compare where and when it is made, across stages and strains.",
           builder="tool"),

    # ---------------------------------------------------------------- Explore
    Module("map", "Map", "the genome as a space you choose the axes of",
           "Explore", opens="maps", icon="umap",
           instruction="Put your genes on the map and read what surrounds them.",
           builder="shell:map_module"),
    Module("network", "Network", "a gene and its neighbourhood, edges you define",
           "Explore", opens="star map", icon="graph_builder",
           instruction="Name a gene to bring its region into focus; tick the links that count.",
           builder="shell:network_module"),
    Module("cell", "Cell", "the parasite drawn, with the genes placed on it",
           "Explore", opens="find", icon="image_scatter",
           instruction="Choose a compartment or a gene to see where it sits in the cell.",
           builder="tool"),
    Module("gallery", "Gallery", "every map built so far, to compare rather than rebuild",
           "Explore", opens="gallery", icon="trellis",
           instruction="Click a map to show it; compare maps rather than rebuild them.",
           builder="tool"),

    # -------------------------------------------------------------------- Ask
    Module("questions", "Questions", "the catalogue of what can be asked, and of what cannot",
           "Ask", opens="analysis:8 · Questions", icon="experiment_design",
           instruction="Read what can be asked, and what the answer would have to clear.",
           builder="shell:infer_module"),
    Module("infer", "Infer", "answer a question, choosing the evidence it may use",
           "Ask", opens="analysis", icon="ml_analyze",
           instruction="Choose a question, then choose what evidence it may use.",
           builder="shell:infer_module"),
    Module("strategies", "Strategies", "39 ways to turn the data into a claim, each graded",
           "Ask", opens="strategies", icon="model_zoo",
           instruction="Choose a strategy to see what it does, what it may use and how well it did.",
           builder="tool"),
    Module("scorecard", "Scorecard", "how well each answer did on genes it was not shown",
           "Ask", opens="discoveries", icon="model_compare",
           instruction="Choose a label to see which methods recover it, and how often they are right.",
           builder="tool"),
    Module("runs", "Runs", "every inference saved, with the data and settings it used",
           "Ask", needs=("a panel that reads results/harness_holdout_log.jsonl",),
           icon="run_history",
           instruction="Pick a saved run to see what it asked, what it used and what came back."),
    Module("fill", "Fill", "write answers into the map, where they clear the bar",
           "Ask", needs=("a question that has cleared its 90% bar",), icon="regression",
           instruction="Check the track record, then write in only the answers that clear the bar."),

    # --------------------------------------------------------- Bring your own
    Module("drop", "Drop", "a CSV, TSV, XLS or XLSX of genes, with or without values",
           "Bring your own", opens="gene list", icon="download",
           instruction="Drop a CSV, TSV or Excel file of genes here, or choose one.",
           builder="tool"),
    Module("screen", "Screen", "join your screen to everything already recorded",
           "Bring your own", opens="compare", icon="hit_list",
           instruction="Load your screen table and compare it with the recorded conditions.",
           builder="tool"),
    Module("list", "List", "what a set of genes has in common, and what it points to",
           "Bring your own", needs=("a set-level enrichment view",), icon="tabulate",
           instruction="Hand it a set of genes to see what they share, and how surprising that is."),
    Module("export", "Export", "take a list, a table or a figure out",
           "Bring your own", opens="action:export_visible", icon="methods_export",
           instruction="Choose the columns to take out, then where to write them."),

    # ------------------------------------------------------------------- Data
    Module("datasets", "Datasets", "all 130 sources, their coverage and their provenance",
           "Data", opens="action:open_dataset_browser", icon="data_manager",
           instruction="Browse the sources; each one says what it covers and where it came from."),
    Module("information", "Information", "the gene x trait map, and how sure each cell is",
           "Data", opens="information map", icon="feature_dict",
           instruction="Search a gene or choose a trait; every cell shows its source and certainty.",
           builder="tool"),
    Module("jobs", "Jobs", "what is running, what it is costing, and how to stop it",
           "Data", opens="jobs", icon="queue",
           instruction="Select a job to stop it; finished jobs keep their outcome here.",
           builder="tool"),
    Module("console", "Console", "what the program is doing, in its own words",
           "Data", opens="console", icon="report",
           instruction="Filter by level or text; copy what you need for a report.",
           builder="tool"),
)


#: Every section a module may be filed under, in the order they are drawn.
#: Computed from :data:`BANDS` so that the names and the order are declared
#: once -- spaCR keeps `SECTION_ORDER` beside the same table for the same
#: reason, and its bug #16i was two of these lists disagreeing.
SECTIONS = tuple(title for title, _note in BANDS)


def section_members(section: str) -> tuple:
    """The modules filed under one section, in catalogue order.

    The one place the grouping happens. The home bands, the dock and the tree
    all go through it, so a section cannot acquire a heading it has no modules
    for, or lose one it does.

    :param section: a section name, as in :data:`SECTIONS`.
    :returns: the modules in it, possibly empty.
    """
    return tuple(m for m in MODULES if m.band == section)


def sections() -> tuple:
    """`((section, note, (Module, ...)), ...)`, in drawing order, no empties.

    :returns: the same grouping :func:`bands` returns, under the spaCR name.
    """
    out = []
    for title, note in BANDS:
        members = section_members(title)
        if members:
            out.append((title, note, members))
    return tuple(out)


def bands() -> tuple:
    """The bands with their modules, in drawing order.

    The ParaScape name for :func:`sections`, kept because the home screen and
    its tests have used it since instruction 85. One implementation, two names:
    two implementations would be the thing this restructure exists to remove.

    :returns: `((title, note, (Module, ...)), ...)`, skipping an empty band.
    """
    return sections()


def get(key: str) -> Module:
    """One module by key.

    :raises KeyError: when no module has that key, rather than returning None and
        failing later at the point of use.
    """
    for m in MODULES:
        if m.key == key:
            return m
    raise KeyError(f"no module {key!r}; have {sorted(x.key for x in MODULES)}")


def tool_modules() -> tuple:
    """The modules built by `dock_module` around a dock that already exists.

    `dock_module.SPECS` was a second list of ten keys with their own titles and
    their own wording for the same ten modules -- "where and when a gene's RNA
    and protein are measured" here and "where and when its RNA and protein are
    measured" there. It is computed from this now.

    :returns: the rows whose builder is ``"tool"``, in catalogue order.
    """
    return tuple(m for m in MODULES if m.builder == "tool")


def shell_module(key: str) -> str:
    """The Python module that builds `key`, when one does, else "".

    Reads :attr:`Module.builder` so that `landing_install` no longer carries a
    dict of six keys that had to be kept in step with this table by hand.

    :param key: a module key.
    :returns: a name such as ``"map_module"``, relative to the parascape
        package, or ``""`` when this module is not built by one.
    :raises KeyError: when there is no such module.
    """
    builder = get(key).builder
    return builder.split(":", 1)[1] if builder.startswith("shell:") else ""


def shell_modules() -> tuple:
    """Every distinct Python module that builds a module shell, in catalogue order.

    Partners is the Network module read by binding evidence and Questions is
    the Infer module opened on its catalogue, so six keys name four builders --
    which is why this de-duplicates rather than returning one per key.
    """
    out = []
    for module in MODULES:
        name = shell_module(module.key)
        if name and name not in out:
            out.append(name)
    return tuple(out)


def unavailable_reason(module: Module, *, docks=(), fillable: int = 0) -> str:
    """Why this module cannot be opened right now, or empty when it can.

    Computed, never stored. `Fill` is the module this exists for: it is open only
    when some question has cleared its bar, and whether one has is a measurement
    that changes as evidence arrives -- so a hardcoded "coming soon" would be
    wrong in both directions.

    :param module: the module to judge.
    :param docks: the dock names the window actually built.
    :param fillable: how many questions have cleared the fill bar.
    :returns: a sentence naming what is missing, or "".
    """
    if module.key == "fill":
        if fillable:
            return ""
        return "no question has cleared its 90% bar yet, so there is nothing to write"
    if module.needs:
        return "needs " + "; ".join(module.needs)
    if not module.opens:
        return "this view has not been built yet"
    if module.opens.startswith("action:"):
        return ""
    if docks:
        dock = module.opens.split(":")[0]
        if dock not in docks:
            return f"the {dock} view is not available in this window"
    return ""


def _catalogue_order(keys) -> tuple:
    """Put module keys into the order the catalogue draws them, dropping any it has not got.

    The bands run Look up, Explore, Ask, Bring your own, Data, which is the order of
    somebody working: read the gene, see where it sits, ask what is missing. Sorting
    the offers alphabetically instead -- which is what :mod:`parascape.ports` returns,
    because a graph has no opinion about which answer to read first -- put Cell before
    Partners for no reason anybody could say out loud.

    :param keys: module keys in any order.
    :returns: the ones that are real modules, in catalogue order, each once.
    """
    wanted = set(keys)
    return tuple(m.key for m in MODULES if m.key in wanted)


def downstream_of_subject(kind: str, *, without: str = "") -> tuple:
    """The modules that can take one subject, ready ones first, in catalogue order.

    This is the answer Einar asked for by name: "all the modules that take 1 gene and
    do something should come downstream of the gene card module" is
    ``downstream_of_subject("gene", without="gene")``.

    Ready first and blocked last is spaCR's own order in the chaining strip
    (``spacr.chaining.next_steps`` sorts ``(not ok, module)``), and it matters here
    because the first offer under a gene card was otherwise Structure, which has no
    view built and is drawn greyed. An offer somebody cannot take should not be the
    first thing they are shown.

    :param kind: a subject from :mod:`parascape.ports`, such as ``"gene"``.
    :param without: a module key to leave out -- normally the one doing the offering,
        since the Gene card is itself a reader of genes and offering itself is noise.
    :returns: module keys, ready ones first.
    """
    keys = [k for k in PORTS.consumers_of(kind) if k != without]
    ordered = _catalogue_order(keys)
    return tuple(sorted(ordered, key=lambda k: (bool(unavailable_reason(get(k))),
                                                ordered.index(k))))


def downstream(key: str) -> tuple:
    """The modules that can work on something `key` hands on, best first.

    COMPUTED, not written down. Until instruction 98 this was ``DOWNSTREAM``: a
    hand-written dict of twenty-two lists saying which module follows which, with no
    reason recorded for any of it and nothing to check it against. A module added to
    :data:`MODULES` and forgotten there simply led nowhere, silently, and a module
    whose inputs changed went on being offered after the wrong thing.

    The answer now comes from the declared subjects (:mod:`parascape.ports`): for each
    subject this module hands on, in the order it declares them, the modules that read
    that subject. So "what follows the Gene card" is answered by what the Gene card
    produces -- a gene, the pairs it names, what could be inferred about it -- and by
    what every other module says it reads.

    Optional readers count as well as required ones, which is the difference between
    this and :func:`parascape.ports.next_modules`. Network reads a gene it can do
    without, and leaving Network out of what follows a gene would have been wrong in
    the one place Einar was looking.

    :param key: a module key. One with nothing declared leads nowhere rather than
        raising, because this is called while drawing.
    :returns: module keys, ready ones first, never including `key` itself.
    """
    try:
        spec = PORTS.module_ports(key)
    except PORTS.UnknownModule:
        return ()
    out: list = []
    for port in spec.produces:
        for candidate in downstream_of_subject(port.kind, without=key):
            if candidate not in out:
                out.append(candidate)
    return tuple(out)


def search(text: str) -> tuple:
    """Modules matching `text`, best first.

    Ranked rather than filtered, because the navigator shows modules alongside
    genes, questions and datasets, and a module that merely contains the word
    should not outrank one named by it.

    :returns: matching modules, most relevant first.
    """
    q = (text or "").strip().lower()
    if not q:
        return ()
    scored = []
    for m in MODULES:
        name, note = m.name.lower(), m.note.lower()
        if name == q:
            rank = 0
        elif name.startswith(q):
            rank = 1
        elif q in name:
            rank = 2
        elif q in note or q in m.band.lower():
            rank = 3
        else:
            continue
        scored.append((rank, m.name, m))
    return tuple(m for _, _, m in sorted(scored))
