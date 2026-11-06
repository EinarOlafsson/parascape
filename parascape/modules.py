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
"""
from __future__ import annotations

from dataclasses import dataclass, field


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
    :param icon: icon name, resolved by the home screen against the icon set.
    """

    key: str
    name: str
    note: str
    band: str
    opens: str = ""
    needs: tuple = ()
    icon: str = ""


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
           "Look up", opens="evidence", icon="gene"),
    Module("structure", "Structure", "the fold, its confidence, and what is known of it",
           "Look up", needs=("a structure view; AF3 summaries are on a gene's evidence for now",),
           icon="structure"),
    Module("partners", "Partners", "what it binds, and how strongly that is evidenced",
           "Look up", opens="evidence", icon="partners"),
    Module("expression", "Expression", "where and when its RNA and protein are measured",
           "Look up", opens="compare:Gene expression", icon="expression"),

    # ---------------------------------------------------------------- Explore
    Module("map", "Map", "the genome as a space you choose the axes of",
           "Explore", opens="maps", icon="map"),
    Module("network", "Network", "a gene and its neighbourhood, edges you define",
           "Explore", opens="star map", icon="network"),
    Module("cell", "Cell", "the parasite drawn, with the genes placed on it",
           "Explore", opens="find", icon="cell"),
    Module("gallery", "Gallery", "every map built so far, to compare rather than rebuild",
           "Explore", opens="gallery", icon="gallery"),

    # -------------------------------------------------------------------- Ask
    Module("questions", "Questions", "the catalogue of what can be asked, and of what cannot",
           "Ask", opens="analysis:8 · Questions", icon="questions"),
    Module("infer", "Infer", "answer a question, choosing the evidence it may use",
           "Ask", opens="analysis", icon="infer"),
    Module("strategies", "Strategies", "39 ways to turn the data into a claim, each graded",
           "Ask", opens="strategies", icon="strategies"),
    Module("scorecard", "Scorecard", "how well each answer did on genes it was not shown",
           "Ask", opens="discoveries", icon="scorecard"),
    Module("runs", "Runs", "every inference saved, with the data and settings it used",
           "Ask", needs=("a panel that reads results/harness_holdout_log.jsonl",),
           icon="runs"),
    Module("fill", "Fill", "write answers into the map, where they clear the bar",
           "Ask", needs=("a question that has cleared its 90% bar",), icon="fill"),

    # --------------------------------------------------------- Bring your own
    Module("drop", "Drop", "a CSV, TSV, XLS or XLSX of genes, with or without values",
           "Bring your own", opens="gene list", icon="drop"),
    Module("screen", "Screen", "join your screen to everything already recorded",
           "Bring your own", opens="compare", icon="screen"),
    Module("list", "List", "what a set of genes has in common, and what it points to",
           "Bring your own", needs=("a set-level enrichment view",), icon="list"),
    Module("export", "Export", "take a list, a table or a figure out",
           "Bring your own", opens="action:export_visible", icon="export"),

    # ------------------------------------------------------------------- Data
    Module("datasets", "Datasets", "all 130 sources, their coverage and their provenance",
           "Data", opens="action:open_dataset_browser", icon="datasets"),
    Module("information", "Information", "the gene x trait map, and how sure each cell is",
           "Data", opens="information map", icon="information"),
    Module("jobs", "Jobs", "what is running, what it is costing, and how to stop it",
           "Data", opens="jobs", icon="jobs"),
    Module("console", "Console", "what the program is doing, in its own words",
           "Data", opens="console", icon="console"),
)


def bands() -> tuple:
    """The bands with their modules, in drawing order.

    :returns: `((title, note, (Module, ...)), ...)`, skipping an empty band.
    """
    out = []
    for title, note in BANDS:
        members = tuple(m for m in MODULES if m.band == title)
        if members:
            out.append((title, note, members))
    return tuple(out)


def get(key: str) -> Module:
    """One module by key.

    :raises KeyError: when no module has that key, rather than returning None and
        failing later at the point of use.
    """
    for m in MODULES:
        if m.key == key:
            return m
    raise KeyError(f"no module {key!r}; have {sorted(x.key for x in MODULES)}")


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
