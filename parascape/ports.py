"""Typed subjects: what each ParaScape module takes, and what it hands on.

WHY THIS EXISTS
---------------
spaCR gives every pipeline module a declared contract (``spacr/ports.py``): the
*kind* of data it consumes and the kind it produces, so that "can this run?" and
"what comes next?" are answered from a table rather than discovered when a run
fails. ParaScape had nothing equivalent. What it had instead was
``modules.DOWNSTREAM`` -- twenty-two hand-written lists saying which module
follows which -- which is the answer without the reason, and which nobody can
check. A module added to ``modules.MODULES`` and forgotten in ``DOWNSTREAM``
simply led nowhere, silently.

That list is **gone** (instruction 98, WP8, 2026-11-12).
:func:`parascape.modules.downstream` computes the same answer from the
declarations below, and ``tests/test_module_catalogue.py`` fails if anything
writes one down again.

THE TRANSLATION FROM spaCR
--------------------------
A spaCR port kind is a FOLDER OF DATA ON DISK: merged arrays, a measurements
database, crops. ParaScape has no pipeline folders. Its equivalent is the
**subject** a module acts on -- one gene, a list of genes, a pair, a map, a
screen of your own, a trait, a question, a claim. Everything the restructure
asks for falls out of that one substitution:

* "the modules downstream of one gene" is :func:`consumers_of` on :data:`GENE`;
* "the modules downstream of a list" is :func:`consumers_of` on
  :data:`GENE_LIST`;
* the tree of what ParaScape can do is :func:`next_modules` walked from a
  subject the user has;
* "why is this greyed out?" is :func:`check_ready`, answered before anything is
  opened.

WHAT A SUBJECT IS, AND WHAT IT IS NOT
-------------------------------------
A subject is a thing somebody can HOLD and hand to a module. It is not a file
format, not a panel, and not a measurement table. The seven Einar named -- a
gene, a gene list, a gene pair, a space, a screen, a trait, a compound -- are
the spine. Six more are declared because the inventory forced them: a module
that takes a question, a claim, a track record, a saved run, a dataset or a
running job has a real subject, and pretending otherwise would leave six
modules with empty contracts.

Two kinds are deliberately kept apart where it would be tempting to merge them,
for the reason spaCR keeps object counts apart from measurements: merging them
lets one thing pass itself off as another. A :data:`CLAIM` is an answer about a
gene; a :data:`SCORECARD` is how often that kind of answer is right on genes it
was not shown. Fill needs BOTH, and a Fill that accepted a claim alone would
write unmeasured guesses into the map -- which is the one thing instruction 98
says must not happen.

The unit words ``gene`` and ``pair`` match :data:`parascape.capabilities.UNITS`,
which is the vocabulary the strategy catalogue already uses for what an output
is about. They are the same words on purpose.

PURITY
------
Standard library only. No Qt, no numpy, no pandas, and nothing else out of
ParaScape -- not even :mod:`parascape.modules`, which imports THIS module rather
than the other way round. A readiness answer has to cost less than the thing it
is protecting, and a question like "what follows a gene?" has to be answerable
in a test, in a script and on a machine with no display.

Public API
----------
``Subject`` / ``Port`` / ``ModulePorts``
    The declarations.
``PORTS``, ``module_ports``, ``register_module_ports``, ``known_modules``
    The registry and the seam a later module joins through.
``producers_of``, ``consumers_of``, ``next_modules``, ``upstream_modules``
    The module graph: the navigator, the tree and the chaining strip read this.
``check_ready``, ``format_readiness``, ``describe_ports``, ``route``
    "Can this module work on what I have, and if not, what gets me there?"
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple, Union

__all__ = [
    "ALL_SUBJECTS",
    "CLAIM",
    "COMPOUND",
    "DATASET",
    "ERROR",
    "GENE",
    "GENE_LIST",
    "GENE_PAIR",
    "JOB",
    "ModulePorts",
    "PORTS",
    "Port",
    "Problem",
    "QUESTION",
    "RUN",
    "Readiness",
    "SCORECARD",
    "SCREEN",
    "SPACE",
    "SUBJECT_MEANINGS",
    "SUBJECT_NOTES",
    "SUBJECT_SOURCES",
    "Subject",
    "TRAIT",
    "UnknownModule",
    "WARNING",
    "check_ready",
    "consumers_of",
    "describe_ports",
    "describe_subject",
    "explain_subject",
    "format_readiness",
    "known_modules",
    "module_ports",
    "next_modules",
    "producers_of",
    "register_module_ports",
    "route",
    "subject_kind",
    "upstream_modules",
]


# --------------------------------------------------------------- the subjects

#: One gene: an accession, with everything recorded against it.
GENE = "gene"
#: Several genes held together: a dropped file, a gated selection, the members
#: of a cluster, the answer to a question.
GENE_LIST = "gene_list"
#: Two genes and the link between them, with the evidence for the link.
GENE_PAIR = "gene_pair"
#: A label or measured column read across genes: a compartment, a functional
#: class, host-facing, a fitness score. The thing a map is coloured by and the
#: thing an inference is about.
TRAIT = "trait"
#: A table of your own: genes with values, from your screen or your collaborator.
#: Distinct from :data:`GENE_LIST` because it carries the values, and a module
#: that needs the values cannot work on a bare list of identifiers.
SCREEN = "screen"
#: A map: every gene placed, so that genes with similar measurements sit close
#: together. Built by Map, kept in the gallery, reusable as a layout.
SPACE = "space"
#: A drug or small molecule. Named by Einar as one of the subjects.
#: NOTHING CONSUMES IT TODAY and :func:`consumers_of` returns an empty tuple for
#: it. That is the inventory finding, not an oversight: drug sensitivity is in
#: ParaScape only as per-gene columns (``datasets.curated_drug_sensitivity``),
#: so there is no module you can hand a compound to. The kind is declared so
#: that the module which does take one joins the graph by declaring a port
#: rather than by editing anything here.
COMPOUND = "compound"

#: A registered question with the bar it has to clear, from the question
#: catalogue. What Infer is pointed at.
QUESTION = "question"
#: An answer about a gene -- a predicted trait with its uncertainty. NEVER
#: phrased as knowledge, and never the same thing as a :data:`SCORECARD`.
CLAIM = "claim"
#: How often a kind of call is right on genes it was not shown, against chance,
#: with the number it was measured on. The thing that turns a claim into
#: something you may act on.
SCORECARD = "scorecard"
#: One saved inference: what was asked, the data and settings it used, and what
#: came back.
RUN = "run"
#: A registered source, its coverage and its provenance.
DATASET = "dataset"
#: Work running now, or finished this session with its outcome kept. Produced
#: by the three modules that hand work to the job runner (`app.Window.jobs`):
#: Map, Infer and Strategies. Declaring it with no producer left Jobs and
#: Console as a dead end whose only advice was "this cannot be reached from
#: here", while every long build in the application was making one.
JOB = "job"

#: Every subject ParaScape declares, in the order they are introduced above:
#: the things you bring first, then the things the application makes.
#: Not a closed set -- :func:`register_module_ports` accepts any kind, so a
#: module written later can bring its own.
ALL_SUBJECTS: Tuple[str, ...] = (
    GENE, GENE_LIST, GENE_PAIR, TRAIT, SCREEN, SPACE, COMPOUND,
    QUESTION, CLAIM, SCORECARD, RUN, DATASET, JOB,
)

#: What to CALL each subject in a sentence: a short noun phrase with its
#: article, and no comma in it. The first version of this table held the full
#: explanation instead, and every message built from it came out as "no how
#: often that kind of answer is right, against chance in hand". A name and an
#: explanation are two different things and are kept in two tables.
SUBJECT_NOTES: Dict[str, str] = {
    GENE: "one gene",
    GENE_LIST: "a list of genes",
    GENE_PAIR: "a pair of genes",
    TRAIT: "a trait",
    SCREEN: "a screen of your own",
    SPACE: "a map",
    COMPOUND: "a compound",
    QUESTION: "a question",
    CLAIM: "a claim",
    SCORECARD: "a track record",
    RUN: "a saved run",
    DATASET: "a dataset",
    JOB: "a job",
}

#: What each subject MEANS, in one line, for a tooltip, a tree node or a
#: report. Longer than the name and allowed its commas.
SUBJECT_MEANINGS: Dict[str, str] = {
    GENE: "one gene, with everything recorded against it",
    GENE_LIST: "several genes held together",
    GENE_PAIR: "two genes and the link between them, with its evidence",
    TRAIT: "a label or measured column read across genes",
    SCREEN: "a table of your own: genes with values",
    SPACE: "every gene placed by how alike the measurements are",
    COMPOUND: "a drug or small molecule",
    QUESTION: "a registered question and the bar it has to clear",
    CLAIM: "an answer about a gene, with its uncertainty",
    SCORECARD: "how often that kind of answer is right, against chance",
    RUN: "one saved inference, with the data and settings it used",
    DATASET: "a registered source, its coverage and its provenance",
    JOB: "work running now, or finished this session",
}

#: How to get a subject WITHOUT going through a module, where there is a way.
#:
#: Written after the first readiness check told somebody with no gene chosen
#: that "nothing in ParaScape hands you one gene yet, so this cannot be reached
#: from here" -- while the navigator field at the top of the page does exactly
#: that. The module graph cannot know about typing, dropping or clicking, so
#: the ways in are named here and put first in the fix, ahead of the modules.
SUBJECT_SOURCES: Dict[str, str] = {
    GENE: "Type a gene in the navigator, or click one on the map.",
    GENE_LIST: "Drop a file of genes, or draw a gate on the map.",
    SCREEN: "Drop your own table of genes and values.",
    TRAIT: "Choose one in the colour-by control, or on the information map.",
    QUESTION: "Choose one from the question catalogue.",
    JOB: "Start something that runs: build a map, answer a question, test a strategy.",
}


def describe_subject(kind: str) -> str:
    """Return what to call ``kind`` in a sentence, or the kind when none is known.

    Falls back rather than raising because this is called to BUILD a message --
    usually a message about something already going wrong -- and a lookup
    failure inside an error report hides the error it was reporting.

    :param kind: a subject name such as :data:`GENE_LIST`.
    """
    return SUBJECT_NOTES.get(kind, kind)


def explain_subject(kind: str) -> str:
    """Return what ``kind`` means, in one line, falling back to the short name.

    :param kind: a subject name such as :data:`GENE_LIST`.
    """
    return SUBJECT_MEANINGS.get(kind, describe_subject(kind))


def _listing(words: Sequence[str]) -> str:
    """Join names the way a sentence does: "a", "a or b", "a, b or c".

    Written because ``" or ".join`` over seven module names produced "cell or
    expression or gene or information or network or partners or structure",
    which nobody reads to the end of.
    """
    words = list(words)
    if not words:
        return ""
    if len(words) == 1:
        return words[0]
    return ", ".join(words[:-1]) + " or " + words[-1]


# ------------------------------------------------------------------- problems

#: A problem that stops the module working.
ERROR = "error"
#: A problem worth saying out loud that does not stop anything.
WARNING = "warning"


@dataclass(frozen=True)
class Problem:
    """One thing wrong, in the shape spaCR uses for a settings pre-flight.

    Deliberately the same four fields as ``spacr.validate.Problem`` -- level,
    what it is about, what is wrong, what to do -- so that a ParaScape
    pre-flight written later can hand its problems and these to one formatter.
    ParaScape has no ``validate.Problem`` of its own to import (``validate.py``
    here is the hold-out check, a different thing with a confusingly similar
    name), so the type is declared rather than borrowed.

    :param level: :data:`ERROR` or :data:`WARNING`.
    :param subject: the subject or role the problem is about.
    :param message: what is wrong, in one line.
    :param fix: what to do about it. Never empty: a problem a user cannot act
        on is a problem reported to the wrong person.
    """

    level: str
    subject: str
    message: str
    fix: str

    @property
    def is_error(self) -> bool:
        """Whether this stops the module working, rather than merely warning."""
        return self.level == ERROR


class UnknownModule(KeyError):
    """No subjects are declared for the requested module key."""


# --------------------------------------------------------------- declarations

@dataclass(frozen=True)
class Subject:
    """A subject with something actually in it.

    The difference between "Infer takes a question" (a :class:`Port`, a
    declaration) and "I have the host-facing question in my hand" (this, a
    value). :func:`check_ready` accepts either, because both questions get
    asked: the navigator asks the first to decide what to offer, and a module
    about to run asks the second.

    :param kind: one of :data:`ALL_SUBJECTS`, or any kind a module declared.
    :param value: what is in it -- an accession, a list of them, a question key.
        An empty value is the useful case: it means "a gene is what I would
        hand you, but I have not named one yet", which is a different answer
        from "I have no gene at all" and is reported differently.
    :param label: how to name it to somebody, when the value does not say it.
    """

    kind: str
    value: object = None
    label: str = ""

    def __bool__(self) -> bool:
        """True when something is actually in it.

        An empty list, an empty string and ``None`` are all emptiness. A zero
        is not: a trait whose value is 0 is a trait that was chosen.
        """
        if self.value is None:
            return False
        if isinstance(self.value, (str, bytes, list, tuple, set, dict, frozenset)):
            return len(self.value) > 0
        return True

    def describe(self) -> str:
        """Return one line naming what this is, for a message."""
        if self.label:
            return self.label
        if not self:
            return f"{describe_subject(self.kind)} (nothing chosen)"
        if isinstance(self.value, (list, tuple, set, frozenset)):
            return f"{len(self.value)} {describe_subject(self.kind)}"
        return f"{describe_subject(self.kind)} {self.value}"


@dataclass(frozen=True)
class Port:
    """One subject a module consumes or produces.

    :param kind: the subject, e.g. :data:`GENE_LIST`. Ports match across modules
        by kind, which is the whole reason the graph can be computed.
    :param role: a short name, unique within the module, naming what the module
        uses it FOR. Two ports of the same kind in one module -- Network takes a
        space to lay genes out and a list to put in focus -- are told apart by
        this and by nothing else.
    :param required: False for a subject the module can work without. Map works
        with no genes chosen at all; it just shows every gene.
    :param one_of: the name of a set of ports of which ANY ONE will do, so that
        a module can say "one gene, OR a list of them" instead of having to
        pick which of the two it needs.

        Added for instruction 98, WP9, and the reason is a wrong answer it gave
        Einar. The Gene card declared one gene as required and a list as
        optional, which reads as "it needs a gene and may also be given a
        list". What it actually does is step THROUGH a list -- its own port
        says so, "a list to step through, one gene at a time" -- and it has a
        button for taking the gene list whole. So somebody holding a list of
        forty genes was told the Gene card could not work on it, in the one
        place a list is the thing in hand. The same was true of Partners,
        whose port says "several genes, whose partners are shown together",
        and of Expression, whose field reads "one or more exact gene
        accessions".

        Only :func:`check_ready` reads it: a required port whose set has any
        member filled counts as satisfied. :func:`required_kinds` is unchanged,
        because what the module NEEDS has not changed -- what changed is the
        admission that two different things can meet the need. With nothing at
        all in hand nothing in the set is filled, so such a module is still
        reported as not ready, with its own reason.
    :param description: one line, in plain words, for :func:`describe_ports` and
        for the reason text on a greyed offer.
    """

    kind: str
    role: str
    required: bool = True
    description: str = ""
    one_of: str = ""


@dataclass(frozen=True)
class ModulePorts:
    """Everything one module takes and everything it hands on.

    :param key: the module key, matching ``parascape.modules.Module.key``. The
        two tables are keyed the same on purpose and
        ``tests/test_module_catalogue.py`` checks that neither grows a key the
        other does not have.
    :param consumes: what you hand it.
    :param produces: what you can carry out of it into the next module. A
        module that only displays what it was given produces nothing, and
        saying so is what keeps :func:`next_modules` from answering
        "everything".
    :param summary: one line describing the module, for :func:`describe_ports`.
    """

    key: str
    consumes: Tuple[Port, ...] = ()
    produces: Tuple[Port, ...] = ()
    summary: str = ""

    def port(self, role: str) -> Port:
        """Return the consumed or produced port named ``role``.

        :param role: the role name.
        :raises KeyError: when this module declares no such role.
        """
        for candidate in self.consumes + self.produces:
            if candidate.role == role:
                return candidate
        raise KeyError(f"{self.key} declares no port with role {role!r}")

    def required_kinds(self) -> Tuple[str, ...]:
        """The subjects this module cannot work without, in declared order."""
        return tuple(p.kind for p in self.consumes if p.required)

    def optional_kinds(self) -> Tuple[str, ...]:
        """The subjects this module can use but does not need."""
        return tuple(p.kind for p in self.consumes if not p.required)


#: Declared ports, keyed by module. Mutated only through
#: :func:`register_module_ports`, so that two declarations claiming one key is
#: an error rather than a silent last-one-wins.
PORTS: Dict[str, ModulePorts] = {}


def register_module_ports(ports: ModulePorts, *,
                          overwrite: bool = False) -> ModulePorts:
    """Add or replace the port declaration of one module.

    The seam a module written after this one joins the graph through, so
    :data:`PORTS` never has to be edited by hand to add a module.

    :param ports: the declaration. Its key is lowercased before storage.
    :param overwrite: allow replacing an existing declaration. Off by default.
    :returns: the stored declaration, for chaining.
    :raises ValueError: when the key is empty, when it is already declared and
        ``overwrite`` is False, or when two of its ports share a role.
    """
    key = str(ports.key).strip().lower()
    if not key:
        raise ValueError("a module port declaration needs a non-empty key")
    if key in PORTS and not overwrite:
        raise ValueError(
            f"ports for {key!r} are already declared; pass overwrite=True to "
            f"replace them")
    roles = [p.role for p in ports.consumes + ports.produces]
    duplicates = sorted({r for r in roles if roles.count(r) > 1})
    if duplicates:
        raise ValueError(
            f"{key}: port roles must be unique, repeated: "
            f"{', '.join(duplicates)}")
    stored = ports if ports.key == key else ModulePorts(
        key=key, consumes=ports.consumes, produces=ports.produces,
        summary=ports.summary)
    PORTS[key] = stored
    return stored


# The declarations themselves. One per module in `parascape.modules.MODULES`,
# in the same band order, so the two tables can be read side by side.
#
# Every one of these was written by opening the module and reading what it
# actually does with what it is given -- not from the name of the module and not
# from `DOWNSTREAM`. The first version of the module catalogue was written from
# memory and pointed six modules at views that did not exist
# (tests/test_modules_in_window.py records it); guessing at an input would be
# the same mistake one layer down, where no render check would catch it.
_DECLARATIONS = (
    # ---------------------------------------------------------------- Look up
    ModulePorts(
        key="gene",
        summary="one gene, read as sentences, with its sources",
        consumes=(
            # One gene OR a list: the card steps through a list and has a button for
            # taking the gene list whole, so being handed a list is not a reason it
            # cannot work. See `Port.one_of`.
            Port(GENE, "gene", one_of="the genes to read",
                 description="the gene to read"),
            Port(GENE_LIST, "list", required=False, one_of="the genes to read",
                 description="a list to step through, one gene at a time"),
        ),
        produces=(
            Port(GENE, "stepped", required=False,
                 description="the gene you stopped at, stepping through a list"),
            Port(GENE_PAIR, "partners", required=False,
                 description="the partners it names, each with its evidence"),
            Port(CLAIM, "inferable", required=False,
                 description="what is not known about it and could be inferred"),
        )),
    ModulePorts(
        key="structure",
        summary="the fold of one gene, its confidence, and what is known of it",
        consumes=(Port(GENE, "gene", description="the gene whose fold to show"),)),
    ModulePorts(
        key="partners",
        summary="what one gene binds, and how strongly that is evidenced",
        consumes=(
            Port(GENE, "gene", one_of="the genes to show partners of",
                 description="the gene whose partners to show"),
            Port(GENE_LIST, "list", required=False, one_of="the genes to show partners of",
                 description="several genes, whose partners are shown together"),
        ),
        produces=(
            Port(GENE_PAIR, "pairs", description="each partner with its evidence"),
            Port(GENE_LIST, "partner_list", required=False,
                 description="the partners as a list, to carry on"),
        )),
    ModulePorts(
        key="expression",
        summary="where and when the RNA and protein of a gene are measured",
        consumes=(
            # "Enter one or more exact reference gene IDs" is what the field of the
            # expression panel says, so a list is as good an input as a gene.
            Port(GENE, "gene", one_of="the genes to profile",
                 description="the gene to profile"),
            Port(GENE_LIST, "list", required=False, one_of="the genes to profile",
                 description="several genes, profiled together"),
        ),
        produces=(Port(TRAIT, "profile", required=False,
                       description="expression as a column you can read genes by"),)),

    # ---------------------------------------------------------------- Explore
    ModulePorts(
        key="map",
        summary="every gene placed by how alike the measurements are",
        consumes=(
            Port(GENE_LIST, "yours", required=False,
                 description="the genes to put on the map and read around"),
            Port(TRAIT, "lens", required=False,
                 description="the label the map is coloured and read by"),
        ),
        produces=(
            Port(SPACE, "map", description="the map itself, reusable as a layout"),
            Port(GENE_LIST, "gated", required=False,
                 description="the genes inside a gate you drew"),
            # Clicking a gene is how most people get one, and leaving it
            # undeclared made `producers_of(GENE)` empty -- so a reader with no
            # gene chosen was told one could not be reached from here, while
            # the map in front of them was handing out genes all day.
            Port(GENE, "picked", required=False,
                 description="the gene you clicked on the map"),
            Port(JOB, "work", required=False,
                 description="the build, while it is running"),
        )),
    ModulePorts(
        key="network",
        summary="the whole network with your gene in focus, on the links you choose",
        consumes=(
            Port(GENE, "focus", required=False,
                 description="the gene the network is sharp around"),
            Port(GENE_LIST, "yours", required=False,
                 description="several genes to bring into focus"),
            Port(SPACE, "layout", required=False,
                 description="a map to place the genes by, instead of the links"),
        ),
        produces=(
            Port(GENE_PAIR, "links", description="the links drawn, with their layer"),
            Port(GENE_LIST, "neighbourhood", required=False,
                 description="the genes one and two steps from yours"),
            Port(GENE, "picked", required=False,
                 description="the gene you travelled to"),
        )),
    ModulePorts(
        key="cell",
        summary="the parasite drawn, with genes placed where they are found",
        consumes=(
            Port(GENE, "gene", required=False,
                 description="a gene, to show where in the cell it sits"),
            Port(TRAIT, "compartment", required=False,
                 description="a compartment, to show what is in it"),
        ),
        produces=(
            Port(GENE_LIST, "in_compartment", required=False,
                 description="the genes placed in the chosen compartment"),
            # The field beside the drawing, not the drawing itself: clicking an
            # organelle ticks a COMPARTMENT (`app.Window.select_compartment`),
            # and only the search field selects a gene. Written after opening
            # the file, because "clicked in the drawing" was the obvious guess
            # and was wrong.
            Port(GENE, "found", required=False,
                 description="the gene found with the field beside the drawing"),
        )),
    ModulePorts(
        key="gallery",
        summary="every map built so far, to compare rather than rebuild",
        consumes=(Port(SPACE, "maps",
                       description="the maps already built, to choose among"),),
        produces=(Port(SPACE, "chosen",
                       description="the map you picked, shown in the main view"),)),

    # -------------------------------------------------------------------- Ask
    ModulePorts(
        key="questions",
        summary="the catalogue of what can be asked, and of what cannot",
        produces=(Port(QUESTION, "catalogue",
                       description="every registered question and the bar it has to clear"),)),
    ModulePorts(
        key="infer",
        summary="answer a question, choosing the evidence it may use",
        consumes=(
            Port(QUESTION, "question", description="the question to answer"),
            Port(GENE_LIST, "yours", required=False,
                 description="your own genes, read against the answer"),
            Port(TRAIT, "target", required=False,
                 description="the label being inferred, when the question names one"),
        ),
        produces=(
            Port(CLAIM, "answers", description="the genes named, each with its uncertainty"),
            Port(GENE_LIST, "candidates", required=False,
                 description="the named genes as a list, to carry on"),
            Port(RUN, "run", description="the inference saved, with its data and settings"),
            Port(JOB, "work", required=False,
                 description="the inference, while it is running"),
        )),
    ModulePorts(
        key="predict",
        # The free-hand counterpart of `infer`: there the question arrives with its
        # holdout already chosen and checked, here the user chooses the column, the
        # method and what to hold out. It was the "Predict a trait" tab of the window
        # instruction 98 retired, and it is declared as its own module because nothing
        # else in ParaScape runs `prediction.run`.
        summary="predict any measured column, with how often such a call is right",
        consumes=(
            Port(TRAIT, "target", description="the measured column to predict"),
            Port(GENE_LIST, "yours", required=False,
                 description="your own genes, to read the hypotheses against"),
        ),
        produces=(
            Port(CLAIM, "hypotheses",
                 description="a call for each gene that carries no measurement, with its support"),
            Port(SCORECARD, "held_out",
                 description="how the method did on the groups it was not trained on"),
            Port(RUN, "run", description="the whole evaluated run, with its splits and provenance"),
            Port(JOB, "work", required=False,
                 description="the evaluation, while it is running"),
        )),
    ModulePorts(
        key="strategies",
        summary="the ways of turning the data into a claim, each graded",
        consumes=(
            Port(TRAIT, "target", required=False,
                 description="the label a strategy is asked to recover"),
            Port(GENE_LIST, "yours", required=False,
                 description="your own genes, to test a strategy on"),
        ),
        produces=(
            Port(CLAIM, "claims", required=False,
                 description="what the chosen strategy claims"),
            Port(SCORECARD, "grades",
                 description="how well each strategy did on genes it was not shown"),
            Port(JOB, "work", required=False,
                 description="the test, while it is running"),
        )),
    ModulePorts(
        key="scorecard",
        summary="how well each answer did on genes it was not shown",
        consumes=(
            Port(TRAIT, "label", required=False,
                 description="the label to score the methods against"),
            Port(CLAIM, "claims", required=False,
                 description="the claims being scored"),
        ),
        produces=(Port(SCORECARD, "record",
                       description="how often that kind of call is right, against chance"),)),
    ModulePorts(
        key="runs",
        summary="every inference saved, with the data and settings it used",
        consumes=(Port(RUN, "runs", description="the saved inferences to look through"),),
        produces=(Port(SCORECARD, "history", required=False,
                       description="how the saved runs did, over time"),)),
    ModulePorts(
        key="fill",
        # Fill needs BOTH, and this is the declaration that says so. A claim on
        # its own is a guess; the track record is what says how often a guess of
        # that kind is right. Writing a claim into the map without one is the
        # failure instruction 98 names: "no analysis, column, calibration or
        # frozen result may be dropped", and an unmeasured claim written in
        # beside measured ones drops the distinction between them.
        summary="write answers into the map, where they clear the bar",
        consumes=(
            Port(CLAIM, "claims", description="the answers to write in"),
            Port(SCORECARD, "record",
                 description="how often that kind of answer is right, which is the bar"),
        ),
        produces=(Port(TRAIT, "filled",
                       description="the trait with the cleared answers written in, marked as inferred"),)),

    # --------------------------------------------------------- Bring your own
    ModulePorts(
        key="drop",
        # Consumes nothing: a file on disk is not a ParaScape subject, it is
        # where one comes from. Declaring a `file` kind would put a subject in
        # the vocabulary that no second module could ever take.
        summary="a file of genes read and resolved, with or without values",
        produces=(
            Port(GENE_LIST, "genes", description="the genes in the file, resolved to accessions"),
            Port(SCREEN, "values", required=False,
                 description="the values beside them, when the file carries any"),
        )),
    ModulePorts(
        key="screen",
        summary="your screen joined to everything already recorded",
        consumes=(Port(SCREEN, "screen",
                       description="your table of genes and their values"),),
        produces=(
            Port(GENE_LIST, "genes", description="the genes in your screen"),
            Port(TRAIT, "measure", required=False,
                 description="your measurement, as a column to read genes by"),
        )),
    ModulePorts(
        key="list",
        summary="what a set of genes has in common, and what it points to",
        consumes=(Port(GENE_LIST, "genes", description="the set to characterise"),),
        produces=(
            # THE LIST ITSELF FIRST, and the order matters: `modules.downstream` walks the
            # produced ports in declared order, so the first column of the module network
            # under this card is the modules that read a list -- which is the answer Einar
            # asked for, "all the moduals that come downstream of gene lists". Declaring
            # the trait first put the readers of a trait there instead.
            Port(GENE_LIST, "same",
                 description="the list itself, to carry into the next module"),
            Port(TRAIT, "shared", required=False,
                 description="what the set shares, as labels with how surprising each is"),
            Port(GENE, "picked", required=False,
                 description="the gene of the list you opened, read on its own"),
        )),
    ModulePorts(
        key="export",
        # Produces nothing: what leaves the application is a file, and a file is
        # not a subject (see `drop`). Export is a sink on purpose, which is why
        # `next_modules("export")` is empty and the tree stops there.
        summary="take a list, a table or a figure out",
        consumes=(
            Port(GENE_LIST, "genes", required=False, description="a list to write out"),
            Port(SPACE, "map", required=False, description="a map to write out"),
            Port(TRAIT, "column", required=False, description="a column to write out"),
            Port(CLAIM, "claims", required=False, description="claims to write out"),
        )),

    # ------------------------------------------------------------------- Data
    ModulePorts(
        key="datasets",
        summary="every source, its coverage and its provenance",
        produces=(Port(DATASET, "sources",
                       description="the registered sources and the computed layers"),)),
    ModulePorts(
        key="information",
        summary="the gene by trait map, and how sure each cell is",
        consumes=(
            Port(GENE, "gene", required=False, description="a gene, to read its row"),
            Port(TRAIT, "trait", required=False, description="a trait, to read its column"),
            Port(DATASET, "sources", required=False,
                 description="the sources behind each cell"),
        ),
        produces=(
            Port(TRAIT, "measured",
                 description="what is measured, with the source and certainty of each cell"),
            Port(CLAIM, "gaps", required=False,
                 description="the cells marked not measured, which are what inference is for"),
            Port(GENE, "picked", required=False,
                 description="the gene whose row you opened"),
        )),
    ModulePorts(
        key="jobs",
        summary="what is running, what it is costing, and how to stop it",
        consumes=(Port(JOB, "jobs", description="the work started this session"),)),
    ModulePorts(
        key="console",
        summary="what the program is doing, in its own words",
        consumes=(Port(JOB, "log", description="what the running work reports"),)),
)

for _declaration in _DECLARATIONS:
    register_module_ports(_declaration)


def known_modules() -> Tuple[str, ...]:
    """Return every module key with declared subjects, sorted."""
    return tuple(sorted(PORTS))


def module_ports(module: str) -> ModulePorts:
    """Return the port declaration for ``module``.

    :param module: the module key.
    :returns: the :class:`ModulePorts` declared for it.
    :raises UnknownModule: when nothing is declared for it. Raising rather than
        returning an empty declaration, because an empty one reads as "this
        module takes nothing", which is a real and different answer that
        Questions, Datasets and Drop give.
    """
    key = str(module).strip().lower()
    if key not in PORTS:
        raise UnknownModule(
            f"no subjects declared for {module!r}; known: "
            f"{', '.join(known_modules())}")
    return PORTS[key]


# ----------------------------------------------------------------- the graph

def producers_of(kind: str) -> Tuple[str, ...]:
    """Return the modules that hand you ``kind``, sorted.

    The answer to "I need a list of genes -- where do I get one?", and the fix
    line on every missing-subject problem :func:`check_ready` reports.

    :param kind: a subject such as :data:`GENE_LIST`.
    """
    return tuple(sorted(
        key for key, spec in PORTS.items()
        if any(p.kind == kind for p in spec.produces)))


def consumers_of(kind: str, *, required_only: bool = False) -> Tuple[str, ...]:
    """Return the modules that take ``kind``, sorted.

    This is the function instruction 98 is written around: "the modules
    downstream of one gene" is ``consumers_of(GENE)`` and "downstream of a
    list" is ``consumers_of(GENE_LIST)``.

    :param kind: a subject such as :data:`GENE`.
    :param required_only: only the modules that NEED it. The navigator wants
        both answers: everything that can use a gene when it is offering, and
        only what needs one when it is explaining why something is greyed.
    """
    return tuple(sorted(
        key for key, spec in PORTS.items()
        if any(p.kind == kind and (p.required or not required_only)
               for p in spec.consumes)))


def next_modules(module: str) -> Tuple[str, ...]:
    """Return the modules that can work on what ``module`` hands you.

    The answer to "you have this -- what now?", which is the strip spaCR shows
    above Run and the module network the navigator draws.

    Only REQUIRED consumption counts. Including optional consumption made
    almost every module follow almost every other, because Map, Network, Cell,
    Infer and Export all take a gene list they can do without -- an offer of
    twenty next steps is the same as no offer at all.

    ONE required subject matching is enough to be offered, which is deliberate
    and is why this is not the last word. Gene hands on a claim, so Fill is
    offered after Gene -- and Fill also needs a track record, so
    :func:`check_ready` greys that offer and says which one is missing. That is
    the division instruction 98 asks for: successors offered, each checked, the
    ones that cannot work shown greyed with the reason. An offer list that
    silently dropped Fill would leave somebody with no idea it existed.

    :param module: the module key.
    :raises UnknownModule: when nothing is declared for it.
    """
    spec = module_ports(module)
    produced = {p.kind for p in spec.produces}
    return tuple(sorted(
        key for key, candidate in PORTS.items()
        if key != spec.key
        and any(p.kind in produced and p.required for p in candidate.consumes)))


def upstream_modules(module: str) -> Tuple[str, ...]:
    """Return the modules that hand you what ``module`` needs.

    The answer to "Fill needs a track record -- who makes one?".

    :param module: the module key.
    :raises UnknownModule: when nothing is declared for it.
    """
    spec = module_ports(module)
    needed = set(spec.required_kinds())
    return tuple(sorted(
        key for key, candidate in PORTS.items()
        if key != spec.key and any(p.kind in needed for p in candidate.produces)))


def route(kind: str, module: str) -> Tuple[str, ...]:
    """Return a shortest chain of modules that turns ``kind`` into what
    ``module`` needs, or ``()`` when there is none.

    The reason this exists rather than a bare "that will not work": somebody
    holding a gene and wanting Fill should be told HOW to get there, not only
    that they cannot. The chain ends with ``module`` itself, so
    ``route(GENE, "fill")`` reads as a set of steps to take.

    Shortest by number of steps, and ties broken alphabetically so the answer
    is the same on every run -- an offer that reorders itself between two
    openings of the same page looks like a bug to whoever is reading it.

    The search is over what you HOLD, not over single subjects, and that is the
    part a first version got wrong. Searching one subject at a time found a way
    to the first thing Fill needs -- a claim -- and called the job done, while
    Fill also needs the track record saying how often such a claim is right.
    A route that stops at half of what a module needs is worse than no route,
    because it reads as an instruction that will work.

    :param kind: the subject in hand. Empty means you hold nothing.
    :param module: the module key you want to reach.
    :returns: module keys in the order they would be used, ending with
        ``module``, or ``()`` when no chain of modules gets there.
    :raises UnknownModule: when nothing is declared for ``module``.
    """
    spec = module_ports(module)
    wanted = frozenset(spec.required_kinds())
    held = frozenset([kind]) if kind else frozenset()
    if wanted <= held:
        return (spec.key,)
    seen = {held}
    frontier: List[Tuple[frozenset, Tuple[str, ...]]] = [(held, ())]
    while frontier:
        nxt: List[Tuple[frozenset, Tuple[str, ...]]] = []
        for have, steps in frontier:
            for key in sorted(PORTS):
                if key == spec.key or key in steps:
                    continue
                candidate = PORTS[key]
                if not frozenset(candidate.required_kinds()) <= have:
                    continue
                gained = have | {p.kind for p in candidate.produces}
                if gained == have:
                    continue                 # a step that adds nothing is not a step
                if wanted <= gained:
                    return steps + (key, spec.key)
                if gained not in seen:
                    seen.add(gained)
                    nxt.append((gained, steps + (key,)))
        frontier = nxt
    return ()


# -------------------------------------------------------------- is it ready?

@dataclass(frozen=True)
class Readiness:
    """Whether a module can work on what you have, and why not when it cannot.

    :param module: the module key that was checked.
    :param ok: True when nothing blocking was found. ``bool(readiness)`` is the
        same answer, so a caller can write ``if check_ready(...)``.
    :param problems: :class:`Problem` instances, errors and warnings mixed.
    :param satisfied: the roles that were filled.
    :param missing: the subjects that were needed and not supplied.
    """

    module: str
    ok: bool
    problems: Tuple[Problem, ...] = ()
    satisfied: Tuple[str, ...] = ()
    missing: Tuple[str, ...] = ()

    def __bool__(self) -> bool:
        """True when the module can work on what was offered."""
        return self.ok

    @property
    def errors(self) -> Tuple[Problem, ...]:
        """The problems that stop it."""
        return tuple(p for p in self.problems if p.is_error)

    @property
    def warnings(self) -> Tuple[Problem, ...]:
        """The problems that do not stop it."""
        return tuple(p for p in self.problems if not p.is_error)

    @property
    def reason(self) -> str:
        """One line saying why the answer is what it is.

        This is the sentence a greyed offer carries. Instruction 98 asks for
        "an offer that cannot work shown greyed WITH THE REASON", so the reason
        is part of the answer rather than something a caller assembles.
        """
        if self.ok:
            found = ", ".join(self.satisfied) if self.satisfied else "nothing needed"
            return f"{self.module} can work on this: {found}"
        first = self.errors[0]
        more = len(self.errors) - 1
        tail = f" (and {more} more)" if more else ""
        return f"{self.module} cannot: {first.message}{tail}"

    def __str__(self) -> str:
        """The full report; see :func:`format_readiness`."""
        return format_readiness(self)


def subject_kind(subject: Union[str, Subject]) -> str:
    """Return the kind of ``subject``, whether it is a bare kind or a value.

    :param subject: a kind name such as :data:`GENE`, or a :class:`Subject`.
    """
    return subject.kind if isinstance(subject, Subject) else str(subject)


def _as_subjects(subject, have) -> Tuple[Subject, ...]:
    """Normalise what a caller offered into subjects, keeping the order given."""
    offered: List[Subject] = []
    items: List[object] = []
    if subject is not None:
        items.append(subject)
    if isinstance(have, (str, Subject)):
        items.append(have)
    elif have:
        items.extend(have)
    for item in items:
        offered.append(item if isinstance(item, Subject)
                       else Subject(str(item), value=True))
    return tuple(offered)


def _fix_for(kind: str, module: str) -> str:
    """Return what to do about a missing ``kind``, naming the modules that make one."""
    direct = SUBJECT_SOURCES.get(kind, "")
    makers = [key for key in producers_of(kind) if key != module]
    if not makers:
        return direct or (f"Nothing in ParaScape hands you {describe_subject(kind)} "
                          f"yet, so this cannot be reached from here.")
    opening = f"Open {_listing(makers)} first to get {describe_subject(kind)}."
    return f"{direct} {opening}".strip() if direct else opening


def check_ready(module: str,
                subject: Union[str, Subject, None] = None,
                *,
                have: Union[str, Subject, Iterable, None] = ()) -> Readiness:
    """Answer "can ``module`` work on what I have?" before anything is opened.

    Every required subject is looked for among what was offered. A subject that
    was offered but is empty -- "a gene is what I would hand you, but I have not
    named one" -- is reported differently from one not offered at all, because
    the two have different fixes: type a gene, versus open the module that
    produces genes.

    A subject offered that the module does not take is a WARNING, not an error.
    Handing Map a question does not stop Map working; it only means the question
    is not what Map will use, and saying so out loud is what stops somebody
    believing their question was taken into account.

    :param module: the module key.
    :param subject: the main thing in hand -- a kind name or a :class:`Subject`.
    :param have: anything else in hand, for a module that needs more than one.
        Accepts a single subject or an iterable of them.
    :returns: a :class:`Readiness`; ``bool(result)`` is the yes or no and
        :attr:`Readiness.reason` is the sentence to show somebody.
    :raises UnknownModule: when nothing is declared for ``module``.
    """
    spec = module_ports(module)
    offered = _as_subjects(subject, have)
    by_kind: Dict[str, List[Subject]] = {}
    for item in offered:
        by_kind.setdefault(item.kind, []).append(item)

    # The sets of ports of which any one will do (`Port.one_of`), resolved before the
    # ports are judged: a required port is satisfied either by its own kind being in
    # hand or by a sibling of its set being in hand. Resolved first because the sibling
    # may be declared after it -- the Gene card names the gene before the list.
    filled_sets = {
        port.one_of for port in spec.consumes
        if port.one_of and any(c for c in by_kind.get(port.kind, []))}

    problems: List[Problem] = []
    satisfied: List[str] = []
    missing: List[str] = []
    for port in spec.consumes:
        candidates = by_kind.get(port.kind, [])
        filled = [c for c in candidates if c]
        if filled or (port.one_of and port.one_of in filled_sets):
            satisfied.append(port.role)
            continue
        if candidates:
            # Offered but empty. A different answer from "not offered", and a
            # different fix: name one, rather than go and get one.
            message = (f"{spec.key} reads {describe_subject(port.kind)}, "
                       f"and none is chosen")
            fix = (f"Choose {describe_subject(port.kind)} first: "
                   f"{port.description or port.role}.")
        else:
            # Nothing of this kind offered at all. When something else would have done
            # instead (`Port.one_of`), say so: "one gene is needed" is a problem reported
            # in a way that hides half of its own fix, when a list of them would do.
            instead = [p for p in spec.consumes
                       if port.one_of and p.one_of == port.one_of and p.kind != port.kind]
            kinds = [describe_subject(port.kind)] + [describe_subject(p.kind) for p in instead]
            message = f"{_listing(kinds)} is needed and none is in hand"
            fix = " ".join([_fix_for(port.kind, spec.key)]
                           + [_fix_for(p.kind, spec.key) for p in instead])
        if port.required:
            missing.append(port.kind)
            problems.append(Problem(ERROR, port.role, message, fix))
        else:
            problems.append(Problem(
                WARNING, port.role,
                message + f", so {spec.key} works without it",
                fix))

    takes = {p.kind for p in spec.consumes}
    for item in offered:
        if item.kind not in takes:
            problems.append(Problem(
                WARNING, item.kind,
                f"{spec.key} does not read {describe_subject(item.kind)}",
                _unused_fix(item.kind, spec.key)))

    return Readiness(module=spec.key, ok=not any(p.is_error for p in problems),
                     problems=tuple(problems), satisfied=tuple(satisfied),
                     missing=tuple(missing))


def _unused_fix(kind: str, module: str) -> str:
    """Return what to do with a subject the module will not read."""
    takers = [key for key in consumers_of(kind) if key != module]
    if not takers:
        return f"Nothing reads {describe_subject(kind)} yet, so it is ignored here."
    verb = "reads" if len(takers) == 1 else "read"
    return (f"{_listing(takers)} {verb} {describe_subject(kind)}. "
            f"Open one of those to use it.")


def format_readiness(readiness: Readiness) -> str:
    """Render a :class:`Readiness` as a block of text.

    Errors first and each with its fix on the line under it, because a report
    that lists five problems and then five fixes makes the reader do the
    pairing.

    :param readiness: the result of :func:`check_ready`.
    """
    verdict = "READY" if readiness.ok else "NOT READY"
    lines = [f"{verdict}: {readiness.module}"]
    if readiness.satisfied:
        lines.append(f"  has: {', '.join(readiness.satisfied)}")
    for problem in readiness.errors:
        lines.append(f"  error   [{problem.subject}] {problem.message}")
        lines.append(f"          fix: {problem.fix}")
    for problem in readiness.warnings:
        lines.append(f"  warning [{problem.subject}] {problem.message}")
        lines.append(f"          fix: {problem.fix}")
    return "\n".join(lines)


def describe_ports(module: str) -> str:
    """Render the declared contract of one module as text.

    :param module: the module key.
    :returns: a multi-line description of what it takes and what it hands on.
    :raises UnknownModule: when nothing is declared for ``module``.
    """
    spec = module_ports(module)
    lines = [f"{spec.key} — {spec.summary}" if spec.summary else spec.key]
    for label, ports in (("takes", spec.consumes), ("hands on", spec.produces)):
        lines.append(f"  {label}:")
        if not ports:
            lines.append("    (nothing)")
            continue
        for port in ports:
            flag = "" if port.required else "  [optional]"
            lines.append(f"    {port.role}: {explain_subject(port.kind)}{flag}")
            if port.description:
                lines.append(f"      {port.description}")
    return "\n".join(lines)
