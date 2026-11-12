"""What ParaScape can be asked, as clickable nodes -- the model behind the home screen.

Instruction 74 phase 1. The complaint it answers is Einar's own, from instruction 70: "as it is now,
I get lost right away." The application already does nearly everything on this list; what was
missing is a statement of what can be asked at all, in the order a person asks it.

**Level 1 is the query unit, because the map is a matrix.** A query is one row (`Gene`), a block of
rows (`List`), or one column (`Label`) -- and there is no fourth, because provenance is an attribute
of every value rather than a place to visit, and a UMAP is a view of the rows rather than a kind of
query. That is the data model rather than a guess at the user's mood, and it is why three categories
replaced the six intention categories this instruction started with.

Nesting is three deep and no deeper: **category -> node -> destination view**. A level-3 node (one
question of the label catalogue) is reached from its parent and is never a top-level card.

Five rules are structural here rather than left to whoever draws the screen:

* **A node says what clicking it does.** `summary` is required, and a title is one or two words --
  the card carries both, and a title that explains itself is a summary in the wrong place.
* **Trust is not a destination.** The "how much can this be trusted" category was deliberately
  dropped: a node that shows a claim must show that claim's probability, its uncertainty and its
  source tier *inline*, and `Node.__post_init__` refuses a claim-bearing node that does not declare
  all three. If reliability lives in a separate room, the default experience is untrustworthy.
* **A question carries its own scoreboard, coverage and state.** `Scoreboard` holds the
  cross-validated estimate *and* the frozen-test result side by side, the genes answerable at the
  bar, and the coverage sentence that says an unmeasured cell is not a negative. "Cannot be answered
  yet" is a state on a question, not a page.
* **Unavailable is a state, not a crash.** A node whose data is absent, or whose view has not been
  built yet, is shown greyed with the reason ("needs a background table"), so the absence is
  visible on the screen instead of being discovered by clicking -- or, worse, by the card being
  left out and the gap never being noticed.
* **Provenance is baked in.** There is no `Sources` category: `Gene ▸ Card` carries the datasets a
  gene appears in and the papers that mention it, and every value behind these cards carries its own
  source tier, which is what those nodes were for.

No Qt is imported here, so the catalogue, the resolver and the scoreboard reader are testable
without an application. `parascape.home_panel` draws these nodes; `parascape.help_index` indexes
them as the ``node`` kind so the search field ranks them beside everything else.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, replace

#: The kinds of thing a node can be. `help_index` gains one entry kind, ``node``, for all of them;
#: this is what the node itself is *about*.
KINDS = ("gene", "question", "strategy", "dataset", "panel", "guide")

#: How a node opens. `resolve` turns a node into an `Action` with one of these kinds, and
#: `home_panel.open_node` is the only place that acts on them.
#:
#: ``dock``       a dock of the main window, by its title, optionally on a named tab
#: ``window``     a method of the main window that opens something (`open_gene_evidence`)
#: ``dialog``     a dialog built on demand (the data downloader)
#: ``catalogue``  the label catalogue, drawn by the home screen itself
#: ``question``   one question of that catalogue, by its registry id
#: ``drop``       the home screen's own drop target, on one of its steps
#: ``unavailable`` nothing happens except saying why, which is the point of it
#: ``gene``       one gene, by its identifier: the gene is selected and its card opened
#: ``entry``      a `help_index` row the search returned -- a strategy, a dataset, a panel, a guide
#:                page -- opened by `help_search.open_entry`, so the navigator reaches everything
#:                the search field already reached
DESTINATION_KINDS = ("dock", "window", "dialog", "catalogue", "question", "drop",
                     "gene", "entry", "unavailable")

#: What must be shown beside a claim, in the same view as the claim itself. Instruction 74, after
#: Einar dropped "how much can this be trusted?" as a category: trust is integrated, not visited.
CLAIM_FIELDS = ("probability", "uncertainty", "source_tier")

#: The bar a frozen-test lower bound has to clear before an answer may be written into the matrix.
#: Instruction 73; repeated here rather than imported so this module stays free of pandas.
FILL_BAR = 0.90

#: Genes in the ME49 information map, for the coverage sentence. The number a coverage statement is
#: a fraction OF matters as much as the numerator.
GENOME_GENES = 8140

#: Longest a node title may be, in words. One is preferred; the card's summary does the explaining.
MAX_TITLE_WORDS = 2


@dataclass(frozen=True)
class Category:
    """One level-1 section: the unit of the matrix being asked about."""

    id: str
    title: str
    purpose: str


#: `Map` is one node name under all three categories: the same embedding view entered with
#: different context -- a gene marked, a set projected, a column coloured. It is a lens on the rows
#: rather than a category, which is why it is not a fourth one, and the embedding *strategies* stay
#: under `Label` (families 1 and 2 are ways to infer, not ways to browse).
MAP_LENS = ("gene_map", "list_map", "label_map")

#: The three query units, in the order they are drawn: exactly the three shapes of the matrix. There
#: is no `Sources` and no `Space`, deliberately -- provenance is an attribute of every value rather
#: than a place, and a UMAP is a view of the rows rather than a fourth kind of query.
CATEGORIES = (
    Category("gene", "Gene", "one protein: one row of the matrix"),
    Category("list", "List", "many proteins: a block of rows, yours or the map's"),
    Category("label", "Label", "one question across all proteins: a column"),
)

#: The category ids, in order, for a caller that only wants the order.
CATEGORY_IDS = tuple(c.id for c in CATEGORIES)


@dataclass(frozen=True)
class Destination:
    """Where a node goes, named by what already exists rather than by a new screen.

    :ivar kind: one of `DESTINATION_KINDS`.
    :ivar target: a dock title, a window method, a question id, a step name.
    :ivar tab: the tab inside the destination dock, when it has tabs.
    :ivar section: which part of the destination to point at, for the status line.
    """

    kind: str
    target: str = ""
    tab: str = ""
    section: str = ""

    def __post_init__(self):
        if self.kind not in DESTINATION_KINDS:
            raise ValueError(f"destination kind must be one of {list(DESTINATION_KINDS)}")


@dataclass(frozen=True)
class Scoreboard:
    """What is known about how well one question is answered, as both numbers, never one.

    Instruction 73 ships two numbers per question or none: the nested cross-validated estimate (how
    a strategy compares, over many folds) and the frozen-test result (what it achieved, once). A cell
    may be filled only on the frozen-test lower bound, so `may_fill` reads that one alone.
    """

    question_id: str
    n_labels: int = 0
    bar: float = FILL_BAR
    cross_validated_lower: float | None = None
    frozen_test_lower: float | None = None
    coverage_at_bar: float = 0.0
    genes_at_bar: int = 0
    may_fill: bool = False
    flags: tuple = ()
    strategy: str = ""
    source: str = ""

    @property
    def scored(self) -> bool:
        """Whether this question has been through the harness at all."""
        return self.cross_validated_lower is not None

    def inline(self) -> str:
        """The scoreboard as one line, both numbers side by side, for the card itself."""
        if not self.scored:
            return "not yet scored: neither a cross-validated estimate nor a frozen test exists"
        cv = f"right ≥ {self.cross_validated_lower:.1%} of the time when the hidden genes rotate"
        if self.frozen_test_lower is None:
            from parascape.harness import plain
            why = '; '.join(plain(f) for f in self.flags) or 'not enough known answers to grade'
            frozen = f"not graded on hidden genes ({why})"
        else:
            frozen = f"right ≥ {self.frozen_test_lower:.1%} of the time on hidden genes"
        verdict = (f"{self.genes_at_bar:,} genes answerable at the {self.bar:.0%} bar"
                   if self.may_fill
                   else f"below the {self.bar:.0%} bar, so it ranks genes and writes nothing")
        return f"{cv} · {frozen} · {verdict}"

    def coverage(self, genome: int = GENOME_GENES) -> str:
        """The coverage sentence, which has to say what an empty cell is not."""
        if not self.n_labels:
            return ("no label count recorded for this question yet; an unmeasured cell is not a "
                    "negative")
        return (f"{self.n_labels:,} of {genome:,} genes have an observed fact for this question; "
                f"the rest are not measured, which is not a negative")


@dataclass(frozen=True)
class Node:
    """One clickable thing, and everything a card needs to draw it honestly.

    :ivar id: stable identifier, used by the resolver and by the search index.
    :ivar title: what a biologist calls it, in one or two words.
    :ivar kind: one of `KINDS` -- what the node is about.
    :ivar category: a `CATEGORY_IDS` member.
    :ivar summary: one line stating exactly what clicking it does. Required.
    :ivar destination: how it opens.
    :ivar available: whether its data is present. False means shown and greyed, never hidden.
    :ivar unavailable_reason: why not, in words, required when `available` is False.
    :ivar advanced: rendered last in its category.
    :ivar level: 2 for a card in a category, 3 for a node reached from one.
    :ivar parent: the level-2 node a level-3 node hangs under.
    :ivar shows_claim: whether the destination displays a claim about a gene. When it does,
        `claim_fields` must name all of `CLAIM_FIELDS`: a claim without its probability, its
        uncertainty and its source tier in the same view is the thing instruction 74 forbids.
    :ivar produces_inference: whether what this node shows is an inference rather than a
        measurement. The hook for the user's own layer: instruction 74 gives every such node a
        **Keep** action that writes the claim -- with its probability, uncertainty, strategy and
        version, the question, the score it was accepted at, whether it cleared the frozen-test bar
        and the date -- to the user's own layer and never to the shipped map. Keep itself is a later
        phase, since the harness is what produces the inference; the flag is here so the layer can
        be added without reshaping the model, and so nothing can quietly route an inference into
        `Facts`, which shows measurements only.
    :ivar question_id: the registry id, for a question node.
    :ivar scoreboard: both numbers, for a question node.
    :ivar keywords: extra words the search should match, so a question is found by its subject.
    :ivar from_search: whether the search made this node for one result rather than it being a card
        of the catalogue. Such a node carries the name of the thing it found -- a gene accession, a
        strategy's own name -- so the one-or-two-word title rule does not apply to it: the rule is
        about what the catalogue CALLS its nodes, not about what a search result is named.
    :ivar entry: the `help_index.HelpEntry` behind an ``entry`` destination.
    """

    id: str
    title: str
    kind: str
    category: str
    summary: str
    destination: Destination
    available: bool = True
    unavailable_reason: str = ""
    advanced: bool = False
    level: int = 2
    parent: str = ""
    shows_claim: bool = False
    claim_fields: tuple = ()
    produces_inference: bool = False
    question_id: str = ""
    scoreboard: Scoreboard | None = None
    keywords: tuple = ()
    from_search: bool = False
    entry: object = None

    def __post_init__(self):
        if not self.id or not self.title:
            raise ValueError("a node needs an id and a title")
        if len(self.title.split()) > MAX_TITLE_WORDS and not self.from_search:
            raise ValueError(f"{self.id}: a node title is at most {MAX_TITLE_WORDS} words, "
                             f"not {self.title!r} -- the summary does the explaining")
        if not self.summary.strip():
            raise ValueError(f"{self.id}: a node must say what clicking it does")
        if self.kind not in KINDS:
            raise ValueError(f"{self.id}: kind must be one of {list(KINDS)}")
        if self.category not in CATEGORY_IDS:
            raise ValueError(f"{self.id}: category must be one of {list(CATEGORY_IDS)}")
        if not self.available and not self.unavailable_reason.strip():
            raise ValueError(f"{self.id}: an unavailable node must say why")
        if self.produces_inference and not self.shows_claim:
            raise ValueError(f"{self.id}: a node that shows an inference shows its probability, "
                             "uncertainty and source tier with it")
        if self.id == "gene_facts" and self.produces_inference:
            raise ValueError("Facts shows measurements only; an inference belongs in Predicted")
        if self.shows_claim and tuple(self.claim_fields) != CLAIM_FIELDS:
            raise ValueError(f"{self.id}: a node that shows a claim shows {list(CLAIM_FIELDS)} "
                             f"inline; it declares {list(self.claim_fields)}")
        if self.kind == "question" and not self.question_id:
            raise ValueError(f"{self.id}: a question node must name its registry question")
        if self.level not in (2, 3):
            raise ValueError(f"{self.id}: nesting is category -> node -> view, so level is 2 or 3")
        if self.level == 3 and not self.parent:
            raise ValueError(f"{self.id}: a level-3 node must name the node it hangs under")

    @property
    def haystack(self) -> str:
        """Title, summary and keywords, for the search index."""
        return " ".join([self.title, self.summary, " ".join(self.keywords)])

    def claim_line(self) -> str:
        """The line a claim-bearing card prints, naming the three things shown beside the claim."""
        if not self.shows_claim:
            return ""
        return "shown with " + ", ".join(CLAIM_FIELDS).replace("source_tier", "source tier")

    def unavailable(self, reason: str) -> "Node":
        """The same node, marked unavailable for `reason` -- how a probe reports absent data."""
        return replace(self, available=False, unavailable_reason=reason)


@dataclass(frozen=True)
class Action:
    """What opening a node actually does, named so it can be checked without a window."""

    node_id: str
    kind: str
    target: str = ""
    tab: str = ""
    section: str = ""
    describe: str = ""


def resolve(node: Node) -> Action:
    """The action that opens `node`, or an ``unavailable`` action carrying the reason.

    The resolver is deliberately total: every node resolves to something, and a node whose data is
    missing resolves to the statement of why rather than to a destination that would fail.
    """
    if not node.available:
        return Action(node.id, "unavailable", describe=node.unavailable_reason)
    d = node.destination
    where = {"dock": lambda: f"the {d.target} panel" + (f", {d.tab}" if d.tab else ""),
             "window": lambda: f"{d.target}()",
             "dialog": lambda: f"the {d.target} dialog",
             "catalogue": lambda: "the label catalogue",
             "question": lambda: f"question {d.target}",
             "drop": lambda: f"the drop target ({d.target or 'file'})",
             "gene": lambda: f"gene {d.target}, its card",
             "entry": lambda: f"{d.section or d.target} ({d.target})",
             "unavailable": lambda: (node.unavailable_reason
                                     or "no destination has been built for this node yet")}[d.kind]
    return Action(node.id, d.kind, d.target, d.tab, d.section, f"{node.title}: {where()}")


# --------------------------------------------------------------------------- scoreboards
#: Where a harness run writes its results, relative to the repository root.
RESULTS_DIR = "results"
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def results_root(root=None) -> str:
    """The directory harness results are written under."""
    return os.path.join(str(root) if root else _REPO, RESULTS_DIR)


def _scoreboard_from(result: dict, source: str) -> Scoreboard:
    """One `Scoreboard` from one `harness.HarnessResult.to_json` record."""
    cv = result.get("cross_validated") or {}
    frozen = result.get("frozen_test") or None
    bar = float(result.get("bar", FILL_BAR))
    may_fill = bool(frozen and frozen.get("reaches_bar"))
    coverage = float((frozen or cv).get("coverage_at_bar", 0.0) or 0.0)
    n_labels = int(result.get("n_labels", 0) or 0)
    return Scoreboard(
        question_id=str(result.get("question_id", "")),
        n_labels=n_labels, bar=bar,
        cross_validated_lower=(None if "precision_lower90" not in cv
                               else float(cv["precision_lower90"])),
        frozen_test_lower=(None if not frozen or "precision_lower90" not in frozen
                           else float(frozen["precision_lower90"])),
        coverage_at_bar=coverage,
        genes_at_bar=int(round(coverage * n_labels)) if may_fill else 0,
        may_fill=may_fill,
        flags=tuple(str(f) for f in result.get("flags", ())),
        strategy=str(result.get("strategy", "")),
        source=source)


def scoreboards(root=None) -> dict:
    """Every question's last harness scoreboard, read from ``results/*/results.json``.

    Newest run wins, because a question re-scored with a corrected evidence list has two files and
    only the later one is true. A results file that is not a harness run is skipped rather than
    guessed at.
    """
    out: dict = {}
    directory = results_root(root)
    if not os.path.isdir(directory):
        return out
    files = []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name, "results.json")
        if os.path.isfile(path):
            try:
                files.append((os.path.getmtime(path), path))
            except OSError:
                continue
    for _when, path in sorted(files):
        try:
            with open(path, encoding="utf8") as handle:
                payload = json.load(handle)
        except (OSError, ValueError):
            continue
        if not isinstance(payload, dict):
            continue
        for key, result in payload.items():
            if isinstance(result, dict) and "cross_validated" in result:
                out[str(key)] = _scoreboard_from(result, os.path.basename(os.path.dirname(path)))
    return out


# --------------------------------------------------------------------------- the catalogue
@dataclass
class Presence:
    """What is actually on this machine, so a node can be greyed with a reason instead of failing.

    Each field is a probe answered once per catalogue build. They are fields rather than calls
    inside `catalogue` so a test can build the unavailable rendering without hiding a release.
    """

    information_map: bool = True
    information_map_reason: str = ""
    gene_list_loaded: bool = False
    pulldown_background: bool = False
    fillable_questions: int = 0

    @classmethod
    def probe(cls, scores=None, gene_list_loaded: bool = False, root=None) -> "Presence":
        """Read the machine: is a map release published, is a background table here, can any
        question fill a cell."""
        from . import information_map as IM
        try:
            root_path = IM.published_root(root)
        except Exception:                                   # a malformed release must not stop the home screen
            root_path = None
        scores = scoreboards(root) if scores is None else scores
        return cls(information_map=root_path is not None,
                   information_map_reason=("needs the information map release; none is published in "
                                           "this checkout"),
                   gene_list_loaded=bool(gene_list_loaded),
                   pulldown_background=False,
                   fillable_questions=sum(1 for s in scores.values() if s.may_fill))


#: Why the pull-down specificity node cannot be opened yet. Named here because the same sentence is
#: the one `information_map_panel` prints beside a pull-down partner.
BACKGROUND_REASON = ("needs a background table, which has not been acquired; until it exists a "
                  "pull-down member is provisional and frequency alone must not be used as a filter")

#: Why no answer may be written into the matrix. Recomputed every build from the frozen-test bounds.
NOTHING_TO_FILL = (f"no question has cleared the {FILL_BAR:.0%} frozen-test bound, so there is "
                   "nothing to write; filling is gated on that bound, never on the cross-validated one")

#: Why a `List` node is greyed before anything has been dropped. Public, because the drop view has
#: to recognise the one obstacle a loaded list removes, as distinct from the obstacles it does not
#: (no background table, no mechanism view).
NEEDS_LIST = "needs a gene list: drop or paste one first"


#: Why a node that needs a view nobody has built yet is greyed. Phase 1 is a navigation layer: a
#: card for work the application cannot do must say that, rather than being left off the screen and
#: rediscovered as a gap.
_NEEDS_SET_ENRICHMENT = ("needs the set-level enrichment view; what exists today is Answer, which "
                         "hides the labels of your list and asks the strategies for them back")
_NEEDS_MECHANISM = ("needs the mechanism view: a process or complex fitted to a set, with the "
                    "explanations it cannot rule out stated. No panel does this yet")
_NEEDS_SET_RANKING = ("needs genome-wide ranking from a set; the strategies rank from labels "
                      "rather than from a set, and are reachable from Label \u25b8 Strategies")
_NEEDS_NEIGHBOUR_SPACE = ("needs the gene-centred neighbour space of instruction 74 phase 3; the "
                          "neighbour strategies themselves are under Label \u25b8 Strategies")


def _gene_nodes() -> list:
    """`Gene`: one protein, one row of the matrix.

    Provenance is baked in rather than given a category of its own: `Card` carries the datasets the
    gene appears in and the papers that mention it, which is what the retired `Sources` nodes were
    for, and every value shown behind these cards carries its own source tier.
    """
    claim = dict(shows_claim=True, claim_fields=CLAIM_FIELDS)
    return [
        Node("gene_card", "Card", "gene", "gene",
             "its structure, the datasets it appears in, the papers that mention it, and what makes "
             "this one unusual",
             Destination("dock", "evidence", section="gene card"),
             keywords=("summary", "overview", "datasets", "references", "provenance", "papers",
                       "unusual", "outlier")),
        Node("gene_facts", "Facts", "gene", "gene",
             "every trait ever recorded for it, with how certain each is and who measured it",
             Destination("dock", "information map", section="gene view"),
             keywords=("information map", "what is known", "traits", "observed", "not measured",
                       "source tier"), **claim),
        Node("gene_predicted", "Predicted", "gene", "gene",
             "what the strategies think is true but nobody has measured, and whether that claim is "
             "strong enough to keep",
             Destination("dock", "discoveries", tab="Inferred claims"),
             produces_inference=True,
             keywords=("inference", "claims", "discoveries", "exploratory", "filled", "keep",
                       "calibrated"), **claim),
        Node("gene_similar", "Similar", "gene", "gene",
             "the genes that most resemble it, by whichever evidence you choose to judge resemblance",
             Destination("dock", "star map", section="neighbour space"),
             available=False, unavailable_reason=_NEEDS_NEIGHBOUR_SPACE,
             keywords=("neighbours", "resemblance", "graphspace", "nearest", "similarity")),
        Node("gene_map", "Map", "gene", "gene",
             "the space of all genes with this one marked, and its neighbourhood around it",
             Destination("dock", "maps", section="one gene marked"),
             keywords=("umap", "embedding", "projection", "gallery", "space", "neighbourhood")),
        Node("gene_network", "Network", "gene", "gene",
             "its neighbours as a map you can walk, with you deciding what counts as a connection",
             Destination("dock", "star map"),
             keywords=("star map", "links", "graph", "crosslink", "co-expression", "edges")),
        Node("gene_partners", "Partners", "gene", "gene",
             "the proteins it has been shown to bind, and how good that evidence is",
             Destination("dock", "information map", section="partners"),
             keywords=("binding", "pairs", "crosslink", "pull-down", "interactome"), **claim),
        Node("gene_expression", "Expression", "gene", "gene",
             "where and when it is switched on, preferring measurements made against a control",
             Destination("dock", "compare", tab="Gene expression"),
             keywords=("rna", "transcript", "profile", "samples", "normalization", "stage")),
        Node("gene_structure", "Structure", "gene", "gene",
             "the numbers behind the fold: confidence, topology, signal peptide and the regions "
             "either side of the membrane",
             Destination("window", "open_gene_evidence", section="structure columns"),
             keywords=("alphafold", "plddt", "3di", "foldseek", "topology", "signal peptide")),
    ]


def _list_nodes(presence: Presence) -> list:
    """`List`: many proteins. Dropping a file is simply how a list comes into being."""
    needs = {} if presence.gene_list_loaded else {"available": False,
                                                  "unavailable_reason": NEEDS_LIST}
    return [
        Node("list_drop", "Drop", "panel", "list",
             "bring your own genes from a spreadsheet or a paste, with any data you have for them",
             Destination("drop", "file"),
             keywords=("import", "paste", "upload", "excel", "csv", "tsv", "parquet", "spreadsheet",
                       "gene list", "identifiers", "drag", "drop")),
        Node("list_profile", "Profile", "panel", "list",
             "what this set has in common that a random set would not: places, domains, conditions, "
             "phenotypes",
             Destination("dock", "analysis", tab="7 \u00b7 Discover"),
             available=False, unavailable_reason=_NEEDS_SET_ENRICHMENT,
             keywords=("enrichment", "over-represented", "shared", "domains", "compartment")),
        Node("list_mechanism", "Mechanism", "panel", "list",
             "the process or complex that best explains the set, and the explanations it cannot "
             "rule out",
             Destination("dock", "analysis", tab="4 \u00b7 Inference"),
             available=False, unavailable_reason=_NEEDS_MECHANISM,
             keywords=("process", "complex", "pathway", "explanation", "hypothesis")),
        Node("list_rank", "Rank", "panel", "list",
             "the rest of the genome ordered by how well each gene would fit this set",
             Destination("dock", "strategies", section="ranking"),
             available=False, unavailable_reason=_NEEDS_SET_RANKING,
             keywords=("prioritize", "ranking", "candidates", "fit", "nearest")),
        Node("list_answer", "Answer", "guide", "list",
             "one question answered for every gene in the set, each with its own certainty",
             Destination("window", "run_my_list"),
             shows_claim=True, claim_fields=CLAIM_FIELDS, produces_inference=True,
             keywords=("predict", "label", "held out", "track record", "keep",
                       "would they have found"), **needs),
        Node("list_pulldown", "Pulldown", "dataset", "list",
             "which of your hits are specific and which turn up in everybody's pull-downs",
             Destination("unavailable", "pulldown_background"),
             available=presence.pulldown_background,
             unavailable_reason="" if presence.pulldown_background else BACKGROUND_REASON,
             keywords=("pulldown_background", "specificity", "ip-ms", "contaminant", "background")),
        Node("list_map", "Map", "panel", "list",
             "the set projected onto the space of all genes, and what surrounds it",
             Destination("dock", "maps", section="the set projected"),
             keywords=("umap", "embedding", "project", "projection", "gallery", "space",
                       "surroundings"), **needs),
        Node("list_compare", "Compare", "panel", "list",
             "two sets against each other: what they share, and what separates them",
             Destination("dock", "compare", tab="Condition comparison"),
             keywords=("overlap", "difference", "two sets", "screen", "contrast")),
        Node("list_export", "Export", "panel", "list",
             "everything attached, with the origin of every value",
             Destination("drop", "export"),
             keywords=("csv", "save", "download", "provenance", "record", "resolution"), **needs),
    ]


#: Why the run ledger cannot be opened yet. The harness already writes one -- every frozen-partition
#: use is logged with its strategy and evidence hashes -- but nothing reads it back into the window.
_NEEDS_RUN_LEDGER = ("needs the run ledger view; the harness records every frozen-partition use in "
                     "results/harness_holdout_log.jsonl with its strategy and evidence hash, and no "
                     "panel reads it back yet")


def _label_nodes(presence: Presence) -> list:
    """`Label`: one question across all proteins, in the order Einar settled on.

    `Strategies` is third rather than last: it stopped being the "advanced" afterthought the first
    draft made of it when he required that nothing of the previous ParaScape be lost -- all 39
    strategies in their nine mechanism families, each keeping its self-test and its scorecard.
    """
    return [
        Node("label_catalogue", "Catalogue", "panel", "label",
             "every question the data can be asked, how well each is currently answered, and for "
             "how many genes",
             Destination("catalogue", "questions"),
             keywords=("questions", "matrix", "registry", "scoreboard", "coverage",
                       "effective independent axes", "cannot be answered")),
        Node("label_infer", "Infer", "panel", "label",
             "run the strategies on one question and predict it for the genes nobody has measured",
             Destination("dock", "discoveries", tab="Inferred claims"),
             shows_claim=True, claim_fields=CLAIM_FIELDS, produces_inference=True,
             keywords=("inference", "predict", "claims", "unmeasured", "discoveries", "keep")),
        Node("label_strategies", "Strategies", "strategy", "label",
             "every way to infer, grouped by mechanism, each stating what it may use and how it was "
             "tested",
             Destination("dock", "strategies"),
             keywords=("strategies", "all strategies", "39", "random forest", "propagation",
                       "baselines", "method", "catalogue", "families", "self-test")),
        Node("label_scorecard", "Scorecard", "panel", "label",
             "how a strategy did on answers it was not allowed to see, against the simple baselines",
             Destination("dock", "discoveries", tab="Strategy coverage"),
             keywords=("scorecard", "baselines", "held out", "precision", "recall", "track record",
                       "frozen test", "cross-validated")),
        Node("label_map", "Map", "panel", "label",
             "the space coloured by the answer, showing where it concentrates",
             Destination("dock", "maps", section="coloured by the answer"),
             keywords=("clusters", "hdbscan", "regions", "map", "umap", "embedding", "lights up",
                       "concentrates")),
        Node("label_runs", "Runs", "panel", "label",
             "every inference you have run, what data it used, and the result -- re-runnable",
             Destination("dock", "discoveries", tab="Inferred claims"),
             available=False, unavailable_reason=_NEEDS_RUN_LEDGER,
             keywords=("runs", "history", "ledger", "metadata", "reproducible", "holdout log",
                       "seed", "partition")),
        Node("label_fill", "Fill", "panel", "label",
             "write only the answers strong enough to keep; everything weaker stays marked "
             "exploratory",
             Destination("catalogue", "fillable"),
             available=presence.fillable_questions > 0,
             unavailable_reason="" if presence.fillable_questions else NOTHING_TO_FILL,
             shows_claim=True, claim_fields=CLAIM_FIELDS, produces_inference=True,
             keywords=("fill", "deploy", "inference", "frozen test", "bound", "calibrated",
                       "keep", "my own layer")),
        Node("label_compare", "Compare", "panel", "label",
             "two questions side by side, including whether they are secretly the same question",
             Destination("dock", "compare", tab="Condition comparison"),
             keywords=("two labels", "redundancy", "agreement", "contrast", "condition",
                       "correlation", "axes")),

    ]


def _question_nodes(questions=None, scores=None) -> list:
    """One level-3 node per registered question, under `label_catalogue`.

    Derived from `question_registry`, so a question that is registered appears here and a question
    that is only talked about does not. A declared gap renders as unanswerable with its reason,
    which is instruction 74's replacement for a "what cannot be answered yet" page.
    """
    from . import question_registry as QR
    scores = scoreboards() if scores is None else scores
    registry = QR.QUESTIONS.values() if questions is None else questions
    out = []
    for q in sorted(registry, key=lambda q: (q.question_class, q.id)):
        if q.status == "retired":
            continue                       # a rejected question is a worked example, not a destination
        board = scores.get(q.id) or Scoreboard(question_id=q.id)
        state = {}
        if not q.answerable or q.status == "gap":
            state = {"available": False,
                     "unavailable_reason": f"cannot be answered yet: {q.label_source}. {q.note}"}
        out.append(Node(
            f"question_{q.id}", _short_title(q.id), "question", "label",
            f"{q.text} — opens its scoreboard, coverage and the genes it reaches",
            Destination("question", q.id),
            level=3, parent="label_catalogue", question_id=q.id, scoreboard=board,
            shows_claim=True, claim_fields=CLAIM_FIELDS,
            keywords=(q.text, q.id, q.answer_type, q.label_source, q.label_family, q.note),
            **state))
    return out


#: One- or two-word names for the registered questions. Derived from the id rather than hand
#: written, so a question registered tomorrow gets a name without anybody editing this file; the
#: exceptions are the ones whose id does not shorten to anything a biologist would say.
_TITLE_WORDS = {"required_without_glucose": "Glucose", "screen_egress_phenotype": "Egress",
                "host_facing": "Host-facing", "function_class": "Function",
                "binding_partners": "Partners",
                "favourable_immune_response_human": "Human immunity",
                "autophagy_involvement": "Autophagy"}


def _short_title(question_id: str) -> str:
    """A one- or two-word name for a question id, which is all a card title may be."""
    if question_id in _TITLE_WORDS:
        return _TITLE_WORDS[question_id]
    words = [w for w in question_id.replace("-", "_").split("_") if w not in
             ("is", "the", "a", "of", "for", "in", "to", "and", "without", "involvement")]
    return " ".join(w.capitalize() for w in words[:MAX_TITLE_WORDS]) or question_id


def catalogue(questions=None, scores=None, presence=None, gene_list_loaded: bool = False,
              root=None) -> tuple:
    """Every node, level 2 then level 3, in the order the home screen draws them.

    `presence` is injectable so the unavailable rendering can be built in a test without removing a
    release from the machine; `scores` and `questions` likewise, so the question nodes can be built
    from a fixture registry.
    """
    scores = scoreboards(root) if scores is None else scores
    if presence is None:
        presence = Presence.probe(scores=scores, gene_list_loaded=gene_list_loaded, root=root)
    nodes = list(_gene_nodes()) + _list_nodes(presence) + _label_nodes(presence)
    if not presence.information_map:
        nodes = [n.unavailable(presence.information_map_reason)
                 if n.destination.target == "information map" else n for n in nodes]
    nodes += _question_nodes(questions, scores)
    order = {c: i for i, c in enumerate(CATEGORY_IDS)}
    nodes.sort(key=lambda n: (order[n.category], n.level, n.advanced))
    return tuple(nodes)


def by_category(nodes=None) -> dict:
    """The level-2 cards of each category, in drawing order, advanced ones last.

    A dict in `CATEGORY_IDS` order, so iterating it draws the home screen top to bottom.
    """
    nodes = catalogue() if nodes is None else nodes
    out = {c: [] for c in CATEGORY_IDS}
    for node in nodes:
        if node.level == 2:
            out[node.category].append(node)
    for cards in out.values():
        cards.sort(key=lambda n: n.advanced)
    return out


def children(node_id: str, nodes=None) -> tuple:
    """The level-3 nodes reached from one level-2 node."""
    nodes = catalogue() if nodes is None else nodes
    return tuple(n for n in nodes if n.level == 3 and n.parent == node_id)


def get(node_id: str, nodes=None) -> Node:
    """One node by id."""
    nodes = catalogue() if nodes is None else nodes
    for node in nodes:
        if node.id == node_id:
            return node
    raise KeyError(f"no such node: {node_id!r}")


def category(category_id: str) -> Category:
    """One category by id."""
    for c in CATEGORIES:
        if c.id == category_id:
            return c
    raise KeyError(f"no such category: {category_id!r}")


# --------------------------------------------------------------------------- catalogue header
@dataclass(frozen=True)
class Axes:
    """How many independent things a question list really holds, and whether that was measured."""

    questions: int
    answerable: int
    effective: int
    measured: bool
    explanation: str


def effective_independent_axes(questions=None) -> Axes:
    """The count that qualifies the whole catalogue, and an honest statement of its status.

    Instruction 73: a thousand questions is not a thousand pieces of knowledge, so the catalogue
    publishes effective independent axes rather than a question count. The redundancy report that
    would measure it -- questions clustered by the correlation of their *labels* -- has not been
    computed, so what is reported until then is an **upper bound**, and it says so. Declared
    redundancy (`redundant_with`) is subtracted, which is the only part that is known today.
    """
    from . import question_registry as QR
    registry = list(QR.QUESTIONS.values() if questions is None else questions)
    live = [q for q in registry if q.status != "retired"]
    answerable = [q for q in live if q.answerable and q.status != "gap"]
    declared = {q.id for q in answerable if q.redundant_with}
    effective = len(answerable) - len(declared)
    measured = bool(declared)
    explanation = (
        "measured only where a question declares what it restates; the redundancy report of "
        "instruction 73, which clusters questions by the correlation of their labels, has not been "
        "computed, so this is an upper bound"
        if not measured else
        "declared redundancy subtracted; the label-correlation clustering is still outstanding")
    return Axes(len(live), len(answerable), effective, measured, explanation)


def catalogue_header(questions=None) -> str:
    """The line above the label catalogue: how many questions, and how many independent axes.

    It sits in the header rather than in a footnote because it qualifies every row beneath it: a
    list of questions read as a list of facts is exactly the failure mode instruction 73 names.
    """
    axes = effective_independent_axes(questions)
    return (f"{axes.questions} questions registered, {axes.answerable} with labels in hand. "
            f"Effective independent axes: {axes.effective} ({axes.explanation}). "
            f"An unmeasured cell is not a negative.")


# --------------------------------------------------------------------------- search
def search_rows(nodes=None) -> list:
    """Every node as a `help_index` entry of kind ``node``, ranked by the same code as everything else.

    The node model extends the index rather than duplicating it: one ranking answers "glucose",
    "star map" and "pulldown_background" the same way, whichever kind the best answer happens to be.
    """
    from . import help_index as HI
    nodes = catalogue() if nodes is None else nodes
    out = []
    for node in nodes:
        where = category(node.category).title
        out.append(HI.entry("node", node.title, f"Home ▸ {where}",
                            f"{node.summary} {' '.join(node.keywords)}"
                            + ("" if node.available else f" unavailable: {node.unavailable_reason}"),
                            node=node.id))
    return out


def rank(query: str, nodes=None, limit: int = 60) -> list:
    """The nodes matching `query`, best first; an empty query is every node in drawing order.

    Used by the home screen's search box, so typing filters and ranks the cards with the ranking
    the search field beside Help already uses.
    """
    from . import help_index as HI
    nodes = catalogue() if nodes is None else nodes
    text = str(query or "").strip()
    if not text:
        return list(nodes)
    by_id = {n.id: n for n in nodes}
    found = HI.search(search_rows(nodes), text, limit=limit, per_kind=None)
    return [by_id[e.data["node"]] for e in found if e.data.get("node") in by_id]


# --------------------------------------------------------------------------- the navigator
#: How many genes one query may return. A query like "GRA" matches dozens; the list is a way in,
#: not a browse.
GENE_LIMIT = 8

#: How many rows of the help index one query may return beside the cards.
ENTRY_LIMIT = 12

#: What a `help_index` kind is about, so a found row carries an honest node kind and lands in a
#: category. A panel, a command or a setting is a way INTO the program, which is `Label`'s business
#: only by accident -- they are filed where they are used from, which is the gene the user is on.
_ENTRY_KIND = {"strategy": ("strategy", "label"), "dataset": ("dataset", "gene"),
               "slot": ("dataset", "gene"), "doc": ("guide", "gene"),
               "tutorial": ("guide", "gene"), "panel": ("panel", "gene"),
               "action": ("panel", "gene"), "setting": ("panel", "gene")}


def gene_node(gene_id: str, product: str = "") -> Node:
    """A found gene, as the node that opens its card.

    The navigator's promise: searching a gene and clicking the result lands on that gene's card with
    the gene selected, not on a generic page. The title is the accession, which is what was typed.
    """
    gene = str(gene_id).strip()
    if not gene:
        raise ValueError("a gene node needs an identifier")
    text = " ".join(str(product or "").split())
    return Node(f"gene:{gene}", gene, "gene", "gene",
                text or "everything known about this gene, with the source of each fact",
                Destination("gene", gene, section="card"),
                from_search=True, keywords=(text, gene))


def entry_node(entry) -> Node:
    """A `help_index` row -- a strategy, a dataset, a panel, a guide page -- as a node card.

    The navigator extends the existing index rather than replacing it, so one field answers with
    genes, questions, strategies, datasets, panels and guide pages together, ranked by the same
    code. Opening one is `help_search.open_entry`, which already knows where each kind lives.
    """
    kind, category = _ENTRY_KIND.get(entry.kind, ("panel", "gene"))
    where = entry.subtitle or HI_KIND_LABEL.get(entry.kind, entry.kind)
    summary = f"{where} — {entry.description}" if entry.description else where
    return Node(f"entry:{entry.kind}:{entry.title}", entry.title, kind, category,
                summary[:400] or f"open {entry.title}",
                Destination("entry", entry.kind, section=entry.title),
                from_search=True, entry=entry, keywords=(entry.subtitle, entry.description))


#: `help_index.KIND_LABEL`, imported lazily so this module stays free of its import at load time.
HI_KIND_LABEL: dict = {}


def genes_matching(query: str, genes, limit: int = GENE_LIMIT) -> list:
    """The genes whose identifier, product or alias matches `query`, as (gene_id, product) pairs.

    `genes` is whatever the window has: a callable taking the query, or an iterable of pairs. It is
    passed in rather than read from a table here, because the node model holds no data.
    """
    text = str(query or "").strip().lower()
    if not text or genes is None:
        return []
    if callable(genes):
        try:
            rows = list(genes(text))
        except Exception:                    # a bad lookup must not break the search box
            return []
    else:
        rows = [(g, p) for g, p in genes
                if text in str(g).lower() or text in str(p).lower()]
    return [(str(g), str(p)) for g, p in rows][:limit]


def search(query: str, nodes=None, genes=None, index=None, limit: int = 40) -> list:
    """Everything one search bar returns, ranked and interleaved: cards, genes and index rows.

    This is the navigator of instruction 70 and 74: one field for genes, questions, strategies,
    datasets, panels and guide pages, every result a node. The ranking is `help_index.search` over
    all of them at once rather than three lists stapled together, so the best answer wins whatever
    kind it happens to be.
    """
    from . import help_index as HI
    HI_KIND_LABEL.update(HI.KIND_LABEL)
    nodes = catalogue() if nodes is None else nodes
    text = str(query or "").strip()
    if not text:
        return list(nodes)
    found = [gene_node(gene, product) for gene, product in genes_matching(text, genes)]
    if index is not None:
        rows = index() if callable(index) else list(index)
        for entry in HI.search([e for e in rows if e.kind != "node"], text, limit=ENTRY_LIMIT):
            found.append(entry_node(entry))
    pool = list(nodes) + found
    by_id = {n.id: n for n in pool}
    ranked = HI.search(search_rows(pool), text, limit=limit, per_kind=None)
    out, seen = [], set()
    for entry in ranked:
        node_id = entry.data.get("node")
        if node_id in by_id and node_id not in seen:
            seen.add(node_id)
            out.append(by_id[node_id])
    return out
