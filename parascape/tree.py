"""The tree of what ParaScape can do, computed from the catalogue and the ports.

Einar asked for this by name: *"there should be a tree of things that parascape can do,
and this progression should guide the user to the modual they want. just like in spacr.
the one in spacr works really well, so study it and implement conceptually the same
thing in parascape."*

WHAT SPACR DOES, AND WHAT IS KEPT
=================================
spaCR files every app under a declared section and then **computes** its home bands, its
tabs and its dock from that one table (``spacr.qt.app.section_members``,
``home_categories``, ``home_bands``). The comments beside them record the bug this
avoids, #16i: Home grouped by one table and the tabs by another, so a section could
acquire a tab it had no apps for. Everything here is computed for the same reason -- a
tree written down a second time would say a module leads somewhere it does not the first
day a port declaration changes, and nobody would find out.

WHAT IS DIFFERENT, AND WHY
==========================
**The headings are ParaScape's own** -- Start, Look up, Explore, Ask, Bring your own,
Data -- which Einar settled on 2026-11-12. They say what the user is trying to do;
spaCR's say which step of an imaging pipeline the app belongs to, and there are no
plates or segmentation models here to have steps of.

**The progression is the subject, not the pipeline stage.** spaCR progresses because its
sections are stages of one run: masks, then measurements, then a plot. ParaScape has no
single run. What it has instead is the subject in hand -- one gene, a list, your own
screen, a trait, a question -- and what becomes possible once you hold one. So the
**Start** heading is a progression over subjects: pick what you have, and the tree names
the modules that read it, ready ones first, each with the reason when one cannot be
taken. That is :func:`parascape.modules.downstream_of_subject`, the same answer the gene
card and the gene list card offer, from the same function.

THE SHAPE
=========
Four levels, and nothing deeper::

    Start                                      a heading
      Home                                     back to the Navigator and the tiles
      One gene                                 a subject you can bring
        Gene          ─┐
        Partners       ├ the modules that read one gene, ready ones first
        Structure ·   ─┘ (greyed, with its reason)
    Look up                                    a heading
      Gene                                     a module, filed where the catalogue files it
        Every record behind one gene           a view behind it
        Then: Partners                         what follows it
    ...

Each module is filed under **exactly one heading**, the one its catalogue row declares.
The rows under *Start* are REFERENCES to those same modules, marked as such
(:attr:`Node.reference`), not a second filing: `tests/test_tree.py` holds that every
module has exactly one node of kind ``module`` and that every reference points at one.
`tests/test_module_catalogue.py::test_no_module_is_in_two_sections` is untouched and
still passes, because it is about the catalogue and the catalogue has not changed.

No Qt here. The tree can be printed, counted and checked without a display;
:mod:`parascape.module_tree` is the widget that draws it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import modules as MOD
from . import ports as PORTS

#: The heading above the progression: not a band of the catalogue, and deliberately so.
#:
#: It holds no module of its own -- Home is not a module and a subject is not a module --
#: and `modules.BANDS` is the declaration that a section HAS modules
#: (`test_no_section_is_empty` says so by name). spaCR kept `SECTION_HELP` out of its own
#: `SECTION_ORDER` for exactly this reason, and its comment records that putting it in was
#: tried first and reverted: "Every Help module is TILELESS by construction, so Help can
#: never satisfy that and should not pretend to."
START = "Start"

#: What the Start heading is for, in one line.
START_NOTE = "What you have, and what it lets you do."

#: Every heading of the tree, in the order they are drawn: the progression first, then the
#: five bands of the catalogue, read from it rather than listed again.
SECTIONS = (START,) + tuple(MOD.SECTIONS)

#: What each heading answers, by name. The five bands bring their own note with them.
SECTION_NOTES = dict([(START, START_NOTE)] + [(title, note) for title, note in MOD.BANDS])


@dataclass(frozen=True)
class Node:
    """One row of the tree.

    :param kind: ``section``, ``home``, ``subject``, ``module``, ``view`` or ``step``.

        * ``section`` -- one of :data:`SECTIONS`. Folds.
        * ``home`` -- back to the Navigator and the module tiles.
        * ``subject`` -- something you can bring: one gene, a list, your own screen.
          Folds, and is not a destination: you cannot open a gene, you type one.
        * ``module`` -- a module, filed under the heading its catalogue row declares.
          Opens, and folds to show what is behind it and what follows it.
        * ``view`` -- a named view declared on that module: a tab of the panel it shows,
          or a window that reads the same thing. Opens.
        * ``step`` -- a module offered as what comes next, or as a reader of a subject.
          Opens the module it names. A reference, not a filing -- see :attr:`reference`.
    :param key: what opening it needs. A module key for ``module`` and ``step``, a subject
        kind for ``subject``, the heading for ``section``, the view name for ``view``.
    :param name: what the row says.
    :param note: one plain line about it, for the tooltip.
    :param owner: for a ``view``, the key of the module it belongs to -- which is what
        opening it needs, since a view is addressed by two names. For a ``step``, the row
        it hangs under, spelled ``module:<key>`` or ``subject:<kind>``.

        Spelled that way rather than as a bare key because one bare key means two things:
        the module ``gene`` and the subject ``gene`` are both called "gene", so "the steps
        under the Gene module" and "the modules that read one gene" could not be told
        apart -- and the second of those includes the first, which made the mistake look
        like an answer.
    :param blocked: why this cannot be opened right now, or "" when it can. A blocked row
        is kept, greyed, with the reason -- never dropped, which is spaCR's rule for its
        chaining strip and Einar's for the module tiles.
    :param reference: True when this row points at a module filed elsewhere. The one thing
        that keeps "a module is in exactly one section" true while still letting the
        progression name it.
    :param children: the rows under it.
    """

    kind: str
    key: str
    name: str
    note: str = ""
    owner: str = ""
    blocked: str = ""
    reference: bool = False
    children: tuple = field(default_factory=tuple)

    @property
    def opens(self) -> bool:
        """Whether clicking this row opens something, rather than only folding it."""
        return self.kind in ("home", "module", "view", "step")

    @property
    def folds(self) -> bool:
        """Whether this row has anything under it to show or hide."""
        return bool(self.children)


def subjects_you_can_bring() -> tuple:
    """The subjects a person can arrive holding, in the order the ports declare them.

    Read from :data:`parascape.ports.SUBJECT_SOURCES`, which is the table that already
    says how to GET a subject without going through a module -- type a gene, drop a file,
    draw a gate. A subject with no way in and no module to hand it to would be a dead row,
    and a compound is exactly that today, which is why it is absent here and recorded in
    `tests/test_ports.py` rather than quietly left out.

    :returns: subject kinds, e.g. ``("gene", "gene_list", "trait", ...)``.
    """
    return tuple(kind for kind in PORTS.ALL_SUBJECTS
                 if kind in PORTS.SUBJECT_SOURCES and MOD.downstream_of_subject(kind))


def _reason(key: str, reasons=None) -> str:
    """Why a module cannot be opened, from the window when it said, else the catalogue.

    :param key: a module key.
    :param reasons: ``{key: reason}`` from `landing_install._reasons`, which knows which
        docks this window really built. Omitted, the catalogue answers on its own.
    """
    if reasons and key in reasons:
        return str(reasons[key] or "")
    try:
        return MOD.unavailable_reason(MOD.get(key))
    except KeyError:
        return ""


def _module_note(key: str) -> str:
    """One line saying what a module is, what it takes and what it hands on.

    Three facts in one sentence because the rail has room for a name and a tooltip, and
    the tooltip is the only place "what each module takes and what follows it" can be read
    without opening anything. Every part of it is computed: the note from the catalogue
    row, the subjects from the declared ports.

    :param key: a module key.
    """
    module = MOD.get(key)
    parts = [module.note.rstrip(".") + "."]
    takes = MOD.takes_words(key)
    parts.append(f"Takes {takes}." if takes else "Needs nothing to open.")
    hands = MOD.hands_on_words(key)
    if hands:
        parts.append(f"Hands on {hands}.")
    return " ".join(parts)


def _steps_after(key: str, reasons=None) -> tuple:
    """The modules that can work on something this module hands on, as ``step`` rows.

    :func:`parascape.modules.downstream` computes them from the declared subjects, which
    is what makes the tree agree with the row of next steps at the foot of the gene card
    and with the module network of the navigator: all three read one function.

    :param key: the module the steps follow.
    :param reasons: ``{key: reason}`` for the modules this window cannot open.
    """
    out = []
    for nxt in MOD.downstream(key):
        module = MOD.get(nxt)
        out.append(Node("step", nxt, module.name, _module_note(nxt),
                        owner=f"module:{key}", blocked=_reason(nxt, reasons),
                        reference=True))
    return tuple(out)


def module_node(key: str, reasons=None, *, with_children: bool = True) -> Node:
    """One module as a tree row, with its views and its next steps under it.

    :param key: a module key.
    :param reasons: ``{key: reason}`` for the modules this window cannot open.
    :param with_children: False for a bare row, which is what a reference wants.
    :returns: the row.
    """
    module = MOD.get(key)
    children = ()
    if with_children:
        views = tuple(Node("view", view.name, view.name, view.shows, owner=key,
                           blocked=_reason(key, reasons))
                      for view in module.views)
        children = views + _steps_after(key, reasons)
    return Node("module", key, module.name, _module_note(key),
                blocked=_reason(key, reasons), children=children)


def subject_node(kind: str, reasons=None) -> Node:
    """One subject you can bring, with the modules that read it under it.

    The rows under it are REFERENCES: the same modules, filed elsewhere under the heading
    their catalogue row declares. Ready ones come first and a blocked one keeps its place
    with its reason, which is the order and the rule
    :func:`parascape.modules.downstream_of_subject` already applies everywhere else.

    :param kind: a subject from :mod:`parascape.ports`.
    :param reasons: ``{key: reason}`` for the modules this window cannot open.
    """
    name = PORTS.describe_subject(kind)
    note = PORTS.explain_subject(kind).rstrip(".") + "."
    how = PORTS.SUBJECT_SOURCES.get(kind, "")
    children = []
    for key in MOD.downstream_of_subject(kind):
        module = MOD.get(key)
        children.append(Node("step", key, module.name, _module_note(key),
                             owner=f"subject:{kind}", blocked=_reason(key, reasons),
                             reference=True))
    return Node("subject", kind, name[:1].upper() + name[1:],
                f"{note} {how}".strip(), children=tuple(children))


def tree(reasons=None) -> tuple:
    """The whole tree, as the six heading rows with everything under them.

    :param reasons: ``{key: reason}`` for the modules this window cannot open, from
        `landing_install._reasons`. Omitted, the catalogue answers on its own, which is
        what a test with no window wants.
    :returns: one :class:`Node` of kind ``section`` per heading, in drawing order.
    """
    out = [Node("section", START, START, START_NOTE, children=(
        Node("home", "home", "Home", "Back to the Navigator and the module tiles. "
                                     "Everything underneath stays as you left it."),
    ) + tuple(subject_node(kind, reasons) for kind in subjects_you_can_bring()))]
    for section in MOD.SECTIONS:
        members = MOD.section_members(section)
        if not members:
            continue
        out.append(Node("section", section, section, SECTION_NOTES.get(section, ""),
                        children=tuple(module_node(m.key, reasons) for m in members)))
    return tuple(out)


def walk(nodes=None) -> tuple:
    """Every row of the tree, depth first, as ``((depth, Node), ...)``.

    What the widget draws and what the tests count. Depth 0 is a heading.

    :param nodes: the tree, from :func:`tree`; built with no reasons when omitted.
    """
    nodes = tree() if nodes is None else nodes
    out = []

    def descend(items, depth):
        """Add one level and then everything under each of its rows."""
        for node in items:
            out.append((depth, node))
            descend(node.children, depth + 1)
    descend(nodes, 0)
    return tuple(out)


def modules_in(nodes=None) -> tuple:
    """The module keys the tree files, in drawing order, references excluded.

    Separate from :func:`walk` because "is every module in the tree exactly once" is the
    question bug #16i was, and it is asked of the FILINGS rather than of the references.

    :param nodes: the tree, from :func:`tree`.
    """
    return tuple(node.key for _depth, node in walk(nodes)
                 if node.kind == "module" and not node.reference)


def views_in(nodes=None) -> tuple:
    """Every ``(module key, view name)`` the tree reaches, in drawing order.

    The 32 declared views, which is the point of including them: a view not in here is a
    view behind the navigator that nothing can name, which is what instruction 98 was
    written about.

    :param nodes: the tree, from :func:`tree`.
    """
    return tuple((node.owner, node.key) for _depth, node in walk(nodes)
                 if node.kind == "view")


def find(key: str, nodes=None):
    """The row that FILES `key`, or None.

    A reference is never returned, so "where does this module live" has one answer.

    :param key: a module key.
    :param nodes: the tree, from :func:`tree`.
    """
    for _depth, node in walk(nodes):
        if node.kind == "module" and not node.reference and node.key == key:
            return node
    return None


def section_of(key: str, nodes=None) -> str:
    """The heading a module is filed under in the tree, or "".

    :param key: a module key.
    :param nodes: the tree, from :func:`tree`.
    """
    nodes = tree() if nodes is None else nodes
    for section in nodes:
        for child in section.children:
            if child.kind == "module" and not child.reference and child.key == key:
                return section.key
    return ""
