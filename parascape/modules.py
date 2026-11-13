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
class View:
    """One named view that belongs to a module without being a module of its own.

    WHY THIS EXISTS (instruction 98, WP10). The inventory counted 29 tabs inside five
    docks and four windows that the navigator could not name at all. Each of them shows
    real measurements, so none could be dropped; and twenty-nine more tiles on Home
    would have answered "everything behind the navigator is unchanged" with a wall.

    A view is the middle answer: it is declared HERE, beside the module whose main
    section already contains it, so it has a name, a plain sentence about what it shows
    and one declared destination. The search field finds it by that name
    (`help_index.panel_entries` is computed from these rows) and
    `landing_install.open_view` opens it. It inherits the ports of the module that owns
    it, because a tab of a panel acts on whatever subject the panel is reading.

    :param name: what it is called, which is the tab label where there is one. What a
        person would type to look for it.
    :param shows: one plain line saying what data is on screen. Not what it is called
        again.
    :param opens: where it is, in the same spelling as :attr:`Module.opens`:
        ``dock``, ``dock:tab`` or ``action:method``.
    """

    name: str
    shows: str
    opens: str


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
        export are windows and dialogs, not docks). Empty means one of two things: the
        view does not exist yet, and the module is unavailable for that reason alone; or
        the module is built by a shell of its own (`builder` starts ``shell:``) and goes
        to no dock, which is Map, Network, Predict and Screen.

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
    :param views: the named views that belong to this module -- the tabs of the panel it
        shows, and the separate windows that read the same thing. See :class:`View`.
        They are declared here rather than in a second table so that a view cannot drift
        away from the module it belongs to, which is bug #16i again.
    :param words: the OTHER WORDS somebody would type for this module -- the ones a
        parasitologist uses rather than the name on the tile. "interactome" for Partners,
        "LOPIT" and "rhoptry" for Cell, "UMAP" for Map, "CRISPR" and "fitness" for Screen.

        WHY THEY ARE DECLARED HERE (instruction 98, WP4). The search field could already
        match the name and the one-line note, which is why typing "interactome" found
        nothing at all: no tile is called that and no note contains it. A synonym table of
        its own would be a second table about the same twenty-three modules, which is the
        failure this whole instruction exists to remove, so the words sit on the row they
        belong to. They are ranked BELOW every part of the name and above a hit in the
        note (:data:`OTHER_WORD`), because somebody who typed the name of a module means
        that module.

        Words that come free from the declared subjects are NOT repeated here: a module
        that reads one gene is already findable by "gene" through
        :func:`what_it_does`, which reads the ports.
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
    views: tuple = ()
    words: tuple = ()

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
    Module("gene", "Gene card", "everything recorded about one gene, with its sources",
           "Look up", opens="evidence", icon="profiler",
           instruction="Type a gene, or choose a list and step through it.",
           builder="shell:gene_module",
           words=("accession", "locus", "annotation", "product", "protein", "ToxoDB",
                  "orthologue", "gene card", "what is known"),
           views=(
               View("Every record behind one gene",
                    "one gene, column by column, including the empty ones, each with the "
                    "record it came from",
                    "action:open_gene_evidence"),
           )),
    Module("structure", "Structure", "the fold, its confidence, and what is known of it",
           "Look up", needs=("a structure view; AF3 summaries are on the evidence of a gene for now",),
           icon="layer_viewer",
           words=("fold", "AlphaFold", "AF3", "pLDDT", "domain", "PDB", "crystal",
                  "secondary structure"),
           instruction="Name a gene to see its fold and how far the fold is to be trusted."),
    Module("partners", "Partners", "what it binds, and how strongly that is evidenced",
           "Look up", opens="evidence", icon="graph_builder",
           words=("interactome", "interaction", "binding", "complex", "pulldown",
                  "immunoprecipitation", "IP-MS", "crosslink", "XL-MS", "two-hybrid",
                  "proximity labelling", "BioID"),
           instruction="Name a gene; the links drawn are the physical ones, with their evidence.",
           builder="shell:network_module"),
    Module("expression", "Expression", "where and when its RNA and protein are measured",
           "Look up", opens="compare:Gene expression", icon="pca",
           instruction="Choose a gene; compare where and when it is made, across stages and strains.",
           builder="tool",
           words=("RNA-seq", "transcriptome", "transcript", "mRNA", "proteome", "abundance",
                  "stage", "tachyzoite", "bradyzoite", "sporozoite", "life cycle",
                  "coexpression", "differential expression"),
           views=(
               View("Studies & samples",
                    "each quantified study, its samples and how they were grouped",
                    "compare:Studies & samples"),
               View("Coexpression",
                    "the genes whose expression rises and falls with the one chosen",
                    "compare:Coexpression"),
               View("Published expression",
                    "the values the original papers published, as published",
                    "compare:Published expression"),
               View("Deposited matrices",
                    "whole deposited measurement matrices, read from the source file",
                    "compare:Deposited matrices"),
               View("Other measurements",
                    "deposited values that are not expression, with their source and units",
                    "compare:Other measurements"),
               View("Sparse measurements",
                    "measurements recorded for only a few genes, kept sparse rather than "
                    "filled in with zeros",
                    "compare:Sparse measurements"),
               View("Original workbooks",
                    "the supplementary workbooks as the authors shipped them, sheet by sheet",
                    "compare:Original workbooks"),
               View("Reference identifiers",
                    "which identifier in a source maps to which gene here, and what did not map",
                    "compare:Reference identifiers"),
               View("Dataset coverage",
                    "how many genes each dataset actually covers, dataset by dataset",
                    "compare:Dataset coverage"),
               View("Condition comparison",
                    "two groups of samples from one study: which genes differ, by how much "
                    "and how sure",
                    "compare:Condition comparison"),
           )),

    # ---------------------------------------------------------------- Explore
    # `opens` is empty for Map and Network on purpose. It is where a module GOES, and
    # these two go to a shell of their own; it used to name the dock each was loosely
    # associated with (`maps`, `star map`), which was never where either opened and which
    # greyed Network out whenever the star map dock failed to build. Both of those docks
    # are declared views now -- of Gallery and of Network -- so each has exactly one name.
    Module("map", "Map", "every gene, placed by how alike its measurements are",
           "Explore", icon="umap",
           instruction="Put your genes on the map and read what surrounds them.",
           builder="shell:map_module",
           words=("UMAP", "embedding", "projection", "atlas", "latent space", "clusters",
                  "clustering", "neighbours", "dimensionality reduction", "PCA", "t-SNE"),
           views=(
               View("Choose individual measurements",
                    "every measurement block, how many genes it covers, and the scaling and "
                    "missing-value rules the map is built under",
                    "analysis:1 · Data"),
               View("Find the best map for a label",
                    "each combination of settings tried, scored by how well a label the map "
                    "never saw comes back",
                    "analysis:5 · Search"),
           )),
    Module("network", "Network", "every gene and the links between them, yours in focus",
           "Explore", icon="graph_builder",
           instruction="Name a gene to bring its region into focus; tick the links that count.",
           builder="shell:network_module",
           words=("graph", "edges", "neighbourhood", "star map", "interactome", "pathway",
                  "guilt by association", "connected"),
           views=(
               View("Links with their provenance",
                    "one gene at the centre and its links on rings around it, coloured by "
                    "which layer, which strategy or which run of yours says so -- the "
                    "inferred links the Network module does not draw",
                    "star map"),
           )),
    Module("cell", "Cell", "the parasite drawn, with the genes placed on it",
           "Explore", opens="find", icon="image_scatter",
           words=("compartment", "localisation", "localization", "organelle", "rhoptry",
                  "dense granule", "micronème", "microneme", "apicoplast", "inner membrane "
                  "complex", "parasitophorous vacuole", "LOPIT", "hyperLOPIT", "schematic",
                  "where it sits"),
           instruction="Choose a compartment or a gene to see where it sits in the cell.",
           builder="tool"),
    Module("gallery", "Gallery", "every map built so far, to compare rather than rebuild",
           "Explore", opens="gallery", icon="trellis",
           instruction="Click a map to show it; compare maps rather than rebuild them.",
           builder="tool",
           words=("thumbnails", "wall", "ready-made maps", "compare maps", "contact sheet"),
           views=(
               View("Maps shipped ready-made",
                    "the maps that ship with ParaScape, grouped by what each was built "
                    "from, with how well every label maps onto every one of them",
                    "maps"),
           )),

    # -------------------------------------------------------------------- Ask
    Module("questions", "Questions", "the catalogue of what can be asked, and of what cannot",
           "Ask", opens="analysis:8 · Questions", icon="experiment_design",
           instruction="Read what can be asked, and what the answer would have to clear.",
           builder="shell:infer_module",
           words=("hypothesis", "biological question", "catalogue", "what can be asked",
                  "question templates"),
           views=(
               View("Question scoreboard",
                    "every registered question scored on genes held out: its held-out AUC "
                    "against chance, and whether it may write into the map",
                    "action:open_question_scoreboard"),
           )),
    Module("infer", "Infer",
           "ask the data a question; get genes named, with how far to trust them",
           "Ask", opens="analysis", icon="ml_analyze",
           instruction="Choose a question, check what it may use, then answer it.",
           builder="shell:infer_module",
           words=("inference", "prediction", "annotate", "impute", "fill the gaps",
                  "machine learning", "held out", "guess what is missing"),
           views=(
               View("What the groups mean, held out",
                    "one row per measurement the map never saw, scored against the groups "
                    "the map found, with the best single group for each value",
                    "analysis:4 · Inference"),
               View("Hill-climb for claims",
                    "a search over maps and groupings for structure that predicts, read "
                    "back as ranked claims with their numbers and their caveats",
                    "analysis:7 · Discover"),
           )),
    Module("predict", "Predict", "predict any measured column, and measure how often that is right",
           "Ask", icon="ml_analyze",
           words=("classifier", "model", "cross-validation", "holdout", "folds", "AUC",
                  "random forest", "train", "abstention", "free-hand prediction"),
           instruction="Choose a column to predict, a method, and what to hold out.",
           builder="shell:predict_module"),
    Module("strategies", "Strategies", "39 ways to turn the data into a claim, each graded",
           "Ask", opens="strategies", icon="model_zoo",
           instruction="Choose a strategy to see what it does, what it may use and how well it did.",
           builder="tool",
           words=("methods", "approaches", "families", "mechanisms", "self-test",
                  "ways of inferring", "guide"),
           views=(
               View("Guide",
                    "what the chosen strategy does, in words, and what it may and may not use",
                    "strategies:Guide"),
               View("Settings",
                    "every setting of the chosen strategy, with the range each is allowed",
                    "strategies:Settings"),
               View("Results",
                    "what the chosen strategy returned, its verdict and its tables",
                    "strategies:Results"),
               View("Start here",
                    "one question at a time, from what you have to the strategies worth "
                    "running on it",
                    "start here"),
           )),
    Module("scorecard", "Scorecard", "how well each answer did on genes it was not shown",
           "Ask", opens="discoveries", icon="model_compare",
           instruction="Choose a label to see which methods recover it, and how often they are right.",
           builder="tool",
           words=("benchmark", "track record", "AUC", "precision", "recall", "validation",
                  "against chance", "how often right", "class balance", "labels"),
           views=(
               View("Labels and functions",
                    "every label, how many genes carry it and where that membership was "
                    "recorded",
                    "discoveries:Labels and functions"),
               View("Inferred claims",
                    "the claims already computed, each with its certainty and whether "
                    "anything independent agreed",
                    "discoveries:Inferred claims"),
               View("Functional tests",
                    "the frozen functional benchmarks and what each method scored on them",
                    "discoveries:Functional tests"),
               View("Strategy coverage",
                    "which strategies have been run on which labels, and which pairs have not",
                    "discoveries:Strategy coverage"),
               View("Host source status",
                    "each host-facing source, what it covers and whether it is loaded",
                    "discoveries:Host source status"),
               View("Hide labels and see if they come back",
                    "a share of the genes that already carry a label hidden, then scored "
                    "against the hidden ones, reported per category",
                    "analysis:6 · Validation"),
               View("Labels and classes",
                    "label membership, class balance, the source records and the held-out "
                    "scorecards recorded against each",
                    "action:open_label_browser"),
           )),
    Module("runs", "Runs", "every inference saved, with the data and settings it used",
           "Ask", needs=("a panel that reads results/harness_holdout_log.jsonl",),
           icon="run_history",
           words=("history", "log", "provenance", "what I ran", "reproduce", "saved runs"),
           instruction="Pick a saved run to see what it asked, what it used and what came back."),
    Module("fill", "Fill", "write answers into the map, where they clear the bar",
           "Ask", needs=("a question that has cleared its 90% bar",), icon="regression",
           words=("write in", "accept", "commit", "promote a claim", "fill the map"),
           instruction="Check the track record, then write in only the answers that clear the bar."),

    # --------------------------------------------------------- Bring your own
    Module("drop", "Drop", "a CSV, TSV, XLS or XLSX of genes, with or without values",
           "Bring your own", opens="gene list", icon="download",
           words=("upload", "import", "load my file", "CSV", "TSV", "Excel", "XLSX",
                  "spreadsheet", "paste genes"),
           instruction="Drop a CSV, TSV or Excel file of genes here, or choose one.",
           builder="tool"),
    Module("screen", "Screen", "join your screen to everything already recorded",
           "Bring your own", icon="hit_list",
           words=("CRISPR", "knockout", "fitness", "phenotype", "hits", "my screen",
                  "my data", "prioritise candidates", "join my table"),
           instruction="Open your measured table, say which columns it uses, then compare.",
           builder="shell:screen_module"),
    Module("list", "Gene list", "what a list of genes shares that other genes do not",
           "Bring your own", opens="gene list", icon="tabulate",
           words=("gene set", "set", "enrichment", "over-representation", "GO terms",
                  "what they share", "my candidates", "shared labels"),
           instruction="Paste or choose a list of genes to see what they share, and how "
                       "surprising that is.",
           builder="shell:list_module"),
    Module("export", "Export", "take a list, a table or a figure out",
           "Bring your own", opens="action:export_visible", icon="methods_export",
           words=("download", "save", "write out", "figure", "table", "CSV", "supplementary"),
           instruction="Choose the columns to take out, then where to write them."),

    # ------------------------------------------------------------------- Data
    Module("datasets", "Datasets", "all 130 sources, their coverage and their provenance",
           "Data", opens="action:open_dataset_browser", icon="data_manager",
           instruction="Browse the sources; each one says what it covers and where it came from.",
           words=("sources", "studies", "publications", "provenance", "coverage",
                  "where the data came from", "citations", "accessions"),
           views=(
               View("Slot tree",
                    "every information category, its address in each of the three "
                    "hierarchies, and what fills it -- the empty ones coloured",
                    "action:open_slot_tree"),
           )),
    Module("information", "Information", "the gene x trait map, and how sure each cell is",
           "Data", opens="information map", icon="feature_dict",
           instruction="Search a gene or choose a trait; every cell shows its source and certainty.",
           builder="tool",
           words=("what is measured", "gaps", "missing", "coverage", "certainty",
                  "information map", "trait", "matrix", "observed"),
           views=(
               View("One gene",
                    "one gene, trait by trait, with each cell marked measured, a measured "
                    "no, or not measured at all",
                    "information map:One gene"),
               View("One trait",
                    "one trait across every gene: how many carry an observed fact for it, "
                    "by source tier, and the values it takes",
                    "information map:One trait"),
           )),
    Module("jobs", "Jobs", "what is running, what it is costing, and how to stop it",
           "Data", opens="jobs", icon="queue",
           words=("running", "queue", "progress", "stop", "cancel", "memory", "busy"),
           instruction="Select a job to stop it; finished jobs keep their outcome here.",
           builder="tool"),
    Module("console", "Console", "what the program is doing, in its own words",
           "Data", opens="console", icon="report",
           words=("log", "output", "errors", "traceback", "messages", "what went wrong"),
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


def views() -> tuple:
    """Every declared view, as `((Module, View), ...)`, in catalogue order.

    The one place the views are gathered. The search index, the tree and the tests all
    read this rather than walking the rows themselves, so a view can only be reached by
    being declared on a module.
    """
    return tuple((m, v) for m in MODULES for v in m.views)


def module_views(key: str) -> tuple:
    """The views declared on one module.

    :param key: a module key.
    :raises KeyError: when there is no such module, rather than quietly returning none
        and leaving a mistyped key looking like a module with no views.
    """
    return tuple(get(key).views)


def find_view(name: str):
    """The module and the view with this name, matched without regard to case.

    :param name: a view name, as declared.
    :returns: `(Module, View)`.
    :raises KeyError: when nothing is called that.
    """
    wanted = (name or "").strip().lower()
    for module, view in views():
        if view.name.lower() == wanted:
            return module, view
    raise KeyError(f"no view {name!r}; have {sorted(v.name for _m, v in views())}")


def search_views(text: str) -> tuple:
    """The views matching `text`, best first, as `((Module, View), ...)`.

    Ranked by :func:`match_score`, the same bands :func:`search` ranks modules by, so a
    view named by what was typed comes before one that merely mentions it in the line
    about what it shows.
    """
    return tuple((m, v) for kind, m, v in search_all(text) if kind == "view")


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
    if module.builder.startswith("shell:"):
        # A module with its own shell is built by `landing_install`, not by a dock, so
        # `opens` says only where its controls came from. Predict and Screen came out of
        # the window instruction 98 retired and so came from nowhere; greying them for
        # want of a dock would grey the two modules WP10 exists to provide.
        return ""
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
    return tuple(m for kind, m, _v in search_all(text) if kind == "module")


# --------------------------------------------------------------------- the ranking
# WHAT THE SEARCH RANKS ON (instruction 98, WP4), and why in this order. The bands are
# spaCR's (`spacr.qt.help_index._term_score`), so one query ranks the same way in both
# applications, with ONE band added for the words declared on a row (`Module.words`).
#
# Einar asked that the navigator find a module by its name, by the name of a view behind
# it, by WHAT IT DOES, and by "the words a parasitologist would use rather than the
# internal name". Those are four different strengths of evidence that a person meant this
# module, and the order between them is the whole ranking: somebody who typed a name meant
# that thing, and must never have it buried under twenty rows whose description happens to
# mention the word.
#
#: Typed in full: "map" for Map. Nothing beats it.
NAME_EXACT = 100.0
#: A whole word of the name: "gene" for "Gene list".
NAME_WORD = 80.0
#: The name begins with it: "expr" for Expression. Plus a little for how much of the name
#: was typed, so a short name matched nearly in full beats a long one matched at the start.
NAME_STARTS = 60.0
#: A word of the name begins with it: "coex" for "Coexpression".
WORD_STARTS = 50.0
#: Somewhere inside the name: "ap" in Map.
NAME_INSIDE = 40.0
#: One of the OTHER WORDS declared for it (:attr:`Module.words`): "interactome" for
#: Partners. Below every part of the name, because a module NAMED by what was typed is a
#: better answer than one that merely answers to the word; above what it does, because
#: these words were declared for this module deliberately and a description hit was not.
OTHER_WORD = 30.0
#: The letters are all there in order and close enough together: "corexpr" for
#: Coexpression. spaCR measured the spread a fuzzy match is allowed.
FUZZY = 20.0
#: In what it DOES -- the one-line note, what to do first, the section, the plain names and
#: meanings of every subject it reads and hands on, and the line a view says about what is
#: on screen. Lowest on purpose.
WHAT_IT_DOES = 10.0


def what_it_does(module: "Module", view: "View" = None) -> str:
    """Everything about a module that is not its name, lowercased, for matching.

    Computed from the catalogue row and from the declared subjects, never written out a
    second time: the note, what to do first, the section it is filed under, the plain name
    and the meaning of every subject it reads and hands on, and -- for a view -- the line
    saying what is on screen. So "trait" finds every module that reads one without the
    word "trait" being typed anywhere on the row.

    :param module: the module.
    :param view: one of its views, when a view is being matched rather than the module.
    :returns: one lowercased string.
    """
    pieces = [module.note, module.instruction, module.band]
    try:
        spec = PORTS.module_ports(module.key)
    except PORTS.UnknownModule:
        spec = None
    if spec is not None:
        for port in spec.consumes + spec.produces:
            pieces.append(PORTS.describe_subject(port.kind))
            pieces.append(PORTS.explain_subject(port.kind))
            pieces.append(port.description)
    if view is not None:
        pieces.extend([view.shows, module.name])
    return " ".join(str(p) for p in pieces if p).lower()


def _term_score(term: str, name: str, words, does: str):
    """How well one typed word matches one candidate, or None when it does not match.

    :param term: one lowercased word of the query.
    :param name: the name of the module, or of the view.
    :param words: the other words declared for it.
    :param does: what it does, from :func:`what_it_does`, already lowercased.
    :returns: a score in the bands above, higher being better, or None.
    """
    from . import help_index as HI        # stdlib only, and it owns the name splitting
    whole, _tail, parts = HI._name_parts(name)
    if whole == term:
        return NAME_EXACT
    if term in parts:
        return NAME_WORD
    if whole.startswith(term):
        return NAME_STARTS + 10.0 * len(term) / max(len(whole), 1)
    if any(p.startswith(term) for p in parts):
        return WORD_STARTS + 10.0 * len(term) / max(len(whole), 1)
    if term in whole:
        return NAME_INSIDE
    best = None
    for word in words or ():
        low = str(word).lower()
        if low == term or term in low.split():
            best = OTHER_WORD + 5.0
            break
        if low.startswith(term) or term in low:
            best = OTHER_WORD
    if best is not None:
        return best
    span = HI._subsequence_span(term, whole) if len(term) >= 3 else None
    if span is not None:
        return FUZZY + 10.0 * len(term) / max(span, 1)
    if term in does:
        return WHAT_IT_DOES
    return None


def match_score(text: str, name: str, words=(), does: str = ""):
    """How well a whole query matches one candidate, or None when any word of it misses.

    The words of the query are ANDed, as spaCR ands them, so a second word narrows the
    answer rather than widening it: "gene list" is not everything about genes plus
    everything about lists.

    :param text: what was typed.
    :param name: the name of the module or of the view.
    :param words: the other words declared for it.
    :param does: what it does, from :func:`what_it_does`.
    :returns: the summed score, or None when a word of the query matches nothing.
    """
    terms = str(text or "").lower().split()
    if not terms:
        return None
    total = 0.0
    for term in terms:
        one = _term_score(term, name, words, does)
        if one is None:
            return None
        total += one
    return total


def search_all(text: str) -> tuple:
    """Modules AND views matching `text`, best first, from one ranking.

    One ranking over both rather than two lists stapled together, because the best answer
    to "coexpression" is a VIEW -- a tab of the panel the Expression module shows -- while
    the best answer to "expression" is the module itself. Two lists could not say which,
    and before instruction 98 the navigator could not name a view at all.

    A module and a view that score the same put the module first: a module is the coarser
    thing and is what somebody typing one word usually means.

    :param text: what was typed.
    :returns: ``(("module", Module, None) | ("view", Module, View), ...)``, best first.
    """
    if not str(text or "").strip():
        return ()
    scored = []
    for module in MODULES:
        value = match_score(text, module.name, module.words, what_it_does(module))
        if value is not None:
            scored.append((-value, 0, module.name.lower(), ("module", module, None)))
    for module, view in views():
        value = match_score(text, view.name, module.words, what_it_does(module, view))
        if value is not None:
            scored.append((-value, 1, view.name.lower(), ("view", module, view)))
    scored.sort(key=lambda row: row[:3])
    return tuple(row[3] for row in scored)


def takes_words(key: str) -> str:
    """What one module needs, in one plain phrase; "" when it needs nothing at all.

    Said in the words of :mod:`parascape.ports` rather than in words of its own, so the
    tree, the tooltips and the reason written on a greyed offer all name a subject the
    same way.

    :param key: a module key.
    :returns: e.g. ``"one gene"``, or ``"one gene, and a list of genes if you have one"``.
    """
    try:
        spec = PORTS.module_ports(key)
    except PORTS.UnknownModule:
        return ""
    need = [PORTS.describe_subject(k) for k in dict.fromkeys(spec.required_kinds())]
    spare = [PORTS.describe_subject(k) for k in dict.fromkeys(spec.optional_kinds())]
    if need and spare:
        return f"{PORTS.listing(need)}, and {PORTS.listing(spare)} if you have one"
    if need:
        return PORTS.listing(need)
    if spare:
        return f"{PORTS.listing(spare)} if you have one"
    return ""


def hands_on_words(key: str) -> str:
    """What one module hands on, in one plain phrase; "" when it hands on nothing.

    :param key: a module key.
    :returns: e.g. ``"a pair of genes or a list of genes"``.
    """
    try:
        spec = PORTS.module_ports(key)
    except PORTS.UnknownModule:
        return ""
    kinds = dict.fromkeys(p.kind for p in spec.produces)
    return PORTS.listing([PORTS.describe_subject(k) for k in kinds])
