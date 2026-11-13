"""The left rail, drawn as the tree of what ParaScape can do (instruction 98, WP7).

Einar: *"there should be a tree of things that parascape can do, and this progression
should guide the user to the modual they want. just like in spacr. the one in spacr works
really well, so study it and implement conceptually the same thing in parascape."*

WHAT THIS REPLACES. The rail was a flat list: six headings, and under each the modules of
that band, which is `spacr.qt.app`'s dock. It answered "what is there" and not "what do I
do next" -- and it could not name the 29 tabs and four windows behind the navigator at
all, which is the complaint instruction 98 was written about. The rows are the same rows;
what is new is that they have a level, that a module unfolds to the views behind it and
the modules that follow it, and that the first heading is a progression over the SUBJECT
in hand rather than a list of modules.

WHAT IS COMPUTED, AND WHY NONE OF IT IS WRITTEN DOWN HERE. Every row comes from
:mod:`parascape.tree`, which computes the whole tree from the one catalogue
(:mod:`parascape.modules`) and the declared subjects (:mod:`parascape.ports`). This file
only turns rows into widgets. The reason is the bug spaCR records as #16i -- Home grouped
by one table and the tabs by another -- and the shape it would take here: a rail with its
own idea of which module follows which would start lying the first day a port declaration
changed, silently, in the one place a new user looks.

BUILT ON THE DOCK RATHER THAN BESIDE IT. :class:`ModuleTree` subclasses
:class:`parascape.panes.dock.Dock`, the vendored spaCR rail, so the rounded panel, the
width the user drags, the fold marks, the hover ink, the icon sizes and the accessibility
of the edge are the ones spaCR already has, unchanged. What is overridden is which rows
exist, how deep each one sits and when it is shown. The vendored file is not edited: its
own docstring records that a second indented level was REMOVED from spaCR on request, and
ParaScape asking for one is not a reason to put it back into spaCR's copy.

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

NO SMALL CROSS, anywhere: the only marks are the disclosure glyphs ▸ and ▾ that the
folding headings of the Home side panel already use, and an arrow on a next step.
`tests/test_symbols_and_no_cross.py` walks the real window and would fail otherwise.
"""
from __future__ import annotations

from PySide6 import QtCore

from . import tree as TREE
from .panes.dock import Dock, DockRow, SectionHeader
from .panes.foldable import OPEN_MARK, SHUT_MARK

#: How far each level is indented, in pixels before the interface scale. A heading sits at
#: the edge; a module is indented once; a view or a next step twice, so the eye can see at
#: a glance which rows belong to the module above them.
INDENT = (10, 22, 38)

#: What goes before a next step, so a row that LEADS somewhere is told apart from a row
#: that IS something, without colour being the only difference.
STEP_MARK = "→ "


class TreeHeader(SectionHeader):
    """A folding heading of the tree: one of the six, or a subject under Start.

    spaCR's :class:`~parascape.panes.dock.SectionHeader` paints the disclosure glyph and
    the dock toggles it; both are kept. What this adds is a KEY separate from the text, so
    that a subject heading can say "One gene" while folding under the key ``subject:gene``
    -- the heading text is what a person reads and the key is what the fold state is
    remembered by, and making one serve as the other is how a renamed heading loses its
    fold.

    :param key: what the fold state is remembered under.
    :param text: what the heading says.
    :param level: 0 for one of the six headings, 1 for a subject under Start.
    :param parent: parent widget; ownership only.
    """

    def __init__(self, key: str, text: str, level: int = 0, parent=None):
        super().__init__(text, parent)
        self.section = key
        self.setProperty("sectionName", key)
        self.setProperty("treeLevel", str(level))
        self.setAccessibleName(text)


class ModuleTree(Dock):
    """The rail: the six headings, and a progression under each, down to a module.

    ``nav_selected(key)`` fires for Home and for every module, exactly as the flat rail
    fired it, so everything already wired to the rail keeps working.
    ``view_selected(module key, view name)`` is new and fires for one of the 32 declared
    views.

    :param icon_for: ``key -> QIcon | None``, the one symbol provider.
    :param reasons: ``{module key: why it cannot open}`` for this window, from
        `landing_install._reasons`. A module that cannot open is drawn greyed with the
        reason in its tooltip rather than hidden, which is Einar's rule for the tiles and
        spaCR's for its chaining strip.
    :param parent: parent widget; ownership only.
    """

    #: ``(module key, view name)`` when one of the declared views is chosen.
    view_selected = QtCore.Signal(str, str)

    #: Class defaults, because the vendored dock calls :meth:`refresh_visibility`,
    #: :meth:`refresh_icons` and :meth:`apply_theme` from inside its own ``__init__`` --
    #: before anything of this class exists. Overriding them means they run then, and a
    #: Qt widget cannot be given attributes before its base class is constructed, so the
    #: defaults live here instead of being guarded at every reader.
    _entries: tuple = ()
    _extra_rows: tuple = ()
    _tree_built = False

    def __init__(self, icon_for=None, reasons=None, parent=None):
        super().__init__([], icon_for=icon_for, parent=parent)
        self._reasons = dict(reasons or {})
        #: Every row widget, in drawing order, as ``(node, depth, folds above it, widget)``.
        self._entries = []
        #: The view and next-step rows. Kept apart from ``_rows`` on purpose: ``_rows`` is
        #: "one row per module, and Home", which is what the window tests count and what
        #: "every module is reachable" is asked of. A step row is the SAME module named a
        #: second time, and counting it would make that question unanswerable.
        self._extra_rows = []
        self._build()
        self._tree_built = True
        self.refresh_icons()
        self.apply_theme()
        self.refresh_visibility()
        self.setToolTip("What ParaScape can do, from what you have to the module that "
                        "does it. Click a heading to open or close it; click a module to "
                        "open it and see what is behind it and what follows it.")

    # ------------------------------------------------------------------- building
    def _column(self):
        """The layout the rows go into: the one the vendored dock already made.

        Found rather than kept, because the dock builds it without naming it and finding
        it is one line, while a copy of the dock made to name it would be the whole file
        again.
        """
        return self._scroll.widget().layout()

    def _build(self) -> None:
        """Put one widget in the rail for every row of :func:`parascape.tree.tree`.

        Walks the tree rather than the catalogue, so the rail cannot show a different
        shape from the one the tests count.
        """
        column = self._column()
        self._nodes = TREE.tree(self._reasons)
        # The fold key of the row at each depth, so a row knows every heading above it.
        # Kept as a stack while walking rather than searched for afterwards: the widget
        # list is flat and the question asked of it at refresh time is flat too -- "is
        # every heading above me open" -- so the answer is recorded once, here.
        stack: list = []
        for depth, node in TREE.walk(self._nodes):
            del stack[depth:]
            above = tuple(k for k in stack if k)
            stack.append(self.fold_key(node))
            widget = self._widget_for(node, depth)
            if widget is None:
                continue
            column.insertWidget(column.count() - 1, widget)
            self._entries.append((node, depth, above, widget))
        # The first heading open, the rest shut: Start is the progression, so it is the
        # one a person arriving should see without clicking anything.
        self._open = {TREE.START}

    @staticmethod
    def fold_key(node) -> str:
        """What this row is folded under, or "" when it folds nothing.

        One function so the fold state, the visibility pass and the accordion all name a
        fold the same way.

        :param node: a :class:`parascape.tree.Node`.
        """
        if not node.folds:
            return ""
        if node.kind == "section":
            return node.key
        return f"{node.kind}:{node.key}"

    def _widget_for(self, node, depth: int):
        """One row of the tree as a widget, or None when the tree holds nothing for it.

        :param node: a :class:`parascape.tree.Node`.
        :param depth: how deep it sits.
        """
        if node.kind in ("section", "subject"):
            header = TreeHeader(self.fold_key(node), node.name, depth)
            self._watch_header(header)
            header.setToolTip(self._header_tip(node))
            self._headers[self.fold_key(node)] = header
            return header
        row = DockRow(node.key, self._row_text(node), node.note)
        row.node = node
        row.icon_key = node.owner if node.kind == "view" else node.key
        row.setProperty("treeLevel", str(min(depth, len(INDENT) - 1)))
        row.setProperty("treeKind", node.kind)
        row.setProperty("treeBlocked", "true" if node.blocked else "false")
        row.setToolTip(self._row_tip(node))
        row.setAccessibleDescription(self._row_tip(node))
        row.hovered.connect(self._on_row_hovered)
        row.clicked.connect(lambda _checked=False, n=node: self._row_clicked(n))
        if node.kind in ("module", "home"):
            self._rows.append(row)
            self._section_of[node.key] = ""
        else:
            self._extra_rows.append(row)
        return row

    def _watch_header(self, header) -> None:
        """Let the dock light and toggle this heading, as it does its own.

        The dock watches Enter, Leave, focus and a left-button release on every heading it
        makes; a heading made here has to be watched the same way or it would draw its
        glyph and never fold.
        """
        from .panes._compat import _watch_object_events
        _watch_object_events(header, self, (
            QtCore.QEvent.Type.Enter, QtCore.QEvent.Type.Leave,
            QtCore.QEvent.Type.FocusIn, QtCore.QEvent.Type.FocusOut,
            QtCore.QEvent.Type.KeyPress, QtCore.QEvent.Type.MouseButtonRelease))

    def _row_text(self, node) -> str:
        """What a row says: its name, with a mark for what clicking it will do.

        A module that unfolds carries the same disclosure glyph the headings carry, so
        "there is more under this" is said once, in one way, everywhere in the rail. A
        next step carries an arrow, because it leads out of the branch it is drawn in
        rather than deeper into it.

        :param node: a :class:`parascape.tree.Node`.
        """
        if node.kind == "step":
            return STEP_MARK + node.name
        if node.folds:
            mark = OPEN_MARK if self.fold_key(node) in getattr(self, "_open", ()) else SHUT_MARK
            return f"{node.name}  {mark}"
        return node.name

    @staticmethod
    def _header_tip(node) -> str:
        """What a heading of the tree says on hover.

        :param node: the heading row.
        """
        n = len(node.children)
        if node.kind == "subject":
            return (f"{node.note} {n} module{'' if n == 1 else 's'} read it; click to "
                    "open or close them.")
        return f"{node.note} {n} row{'' if n == 1 else 's'}; click to open or close them."

    @staticmethod
    def _row_tip(node) -> str:
        """What one row of the tree says on hover, including why it cannot be opened.

        Every control in ParaScape carries a tooltip, and in the rail the tooltip is the
        only place what a module TAKES and what follows it can be read without opening
        anything -- which is what instruction 98 asks the tree to show.

        IT OPENS WITH THE NAME, which is not the repetition it looks like. The rail is
        220 to 320 px wide and a name longer than that is elided: "Hide labels and see if
        they come back" is drawn as "Hide labels and see if they com…", and the tooltip
        is then the only place the whole of it can be read.

        :param node: the row.
        """
        return f"{node.name} — {ModuleTree._row_body(node)}"

    @staticmethod
    def _row_body(node) -> str:
        """What a row tooltip says after the name of the row.

        :param node: the row.
        """
        from . import modules as MOD
        if node.kind == "home":
            return node.note
        if node.blocked:
            return f"{node.note}\n\nNot available: {node.blocked}."
        if node.kind == "view":
            owner = MOD.get(node.owner).name
            return (f"A view of {owner}: {node.note} Opens it, inside the module where "
                    "that module shows it, and as a panel where it lives somewhere else.")
        if node.kind == "step":
            return f"{node.note} Opens it; it is filed under {MOD.get(node.key).band}."
        n = len(node.children)
        if not n:
            return node.note
        return (f"{node.note} Opens it, and shows the {n} row{'' if n == 1 else 's'} "
                "under it: what is behind it, and what follows it.")

    # -------------------------------------------------------------------- clicking
    def _row_clicked(self, node) -> None:
        """Open what a row names, and unfold a module so what follows it is on screen.

        A module row does BOTH: it opens the module and leaves the rail showing the views
        behind it and the modules that follow it. One click rather than two because the
        rail is the way to a module, not a thing to browse -- and because the alternative,
        a second small control on the row to fold it with, is another mark on a rail Einar
        asked to keep clean.

        Only ONE module per heading stays unfolded. Without that, opening the Expression
        module put its ten views and its next steps between every other module of Look up
        and the heading below, and the rail became the wall of rows the module shell
        exists to avoid.

        :param node: the row that was clicked.
        """
        if node.kind == "module":
            self._unfold_only(node)
            self.nav_selected.emit(node.key)
            return
        if node.kind == "view":
            self.view_selected.emit(node.owner, node.key)
            return
        self.nav_selected.emit(node.key)          # home, and a next step

    def _unfold_only(self, node) -> None:
        """Unfold one module and fold away the others under the same heading.

        :param node: the module row to leave open.
        """
        key = self.fold_key(node)
        if not key:
            return
        section = TREE.section_of(node.key, self._nodes)
        for other, _depth, above, _widget in self._entries:
            if (other.kind == "module" and self.fold_key(other)
                    and section and section in above):
                self._open.discard(self.fold_key(other))
        self._open.add(key)
        self.refresh_visibility()

    def toggle_section(self, section: str) -> bool:
        """Open a shut heading or shut an open one; returns the new state.

        Shutting a heading also FORGETS which module under it was unfolded. Without that,
        shutting Look up and opening it again brought back a module unfolded from minutes
        ago, under a heading the person had just asked to see the whole of.

        :param section: a fold key, from :meth:`fold_key`.
        """
        if section in self._open:
            self._open.discard(section)
            for node, _depth, above, _widget in self._entries:
                if section in above and self.fold_key(node):
                    self._open.discard(self.fold_key(node))
        else:
            self._open.add(section)
        self.refresh_visibility()
        return section in self._open

    def open_to(self, key: str) -> None:
        """Open every heading above one module, and unfold the module itself.

        Used when a module is opened from somewhere else -- a tile, the navigator, a menu
        -- so the rail shows where in the tree the person now is instead of staying shut.

        :param key: a module key.
        """
        for node, _depth, above, _widget in self._entries:
            if node.kind == "module" and node.key == key:
                self._open.update(k for k in above if k)
                self._unfold_only(node)
                # Again, and not only inside `_unfold_only`: a module that hands nothing
                # on and has no view -- Export is one -- folds nothing, so `_unfold_only`
                # has nothing to do and returns, and the heading opened just above would
                # never have been drawn.
                self.refresh_visibility()
                return

    # ------------------------------------------------------------------ appearance
    def _on_row_hovered(self, key: str, entered: bool) -> None:
        """Light the row under the pointer, and name the module it is about.

        The dock lights by KEY, which is right while every row has its own; in the tree a
        next step carries the key of the module it leads to, so lighting by key would
        light that module row as well, somewhere else in the rail. This lights the row
        that was actually entered.

        The explanation stays in the tooltip rather than going to a strip along the foot
        of the window: spaCR has such a strip and ParaScape does not.

        :param key: the key of the row.
        :param entered: whether the pointer arrived or left.
        """
        self._light_row(self.sender() if entered else None)
        if entered:
            self.module_hovered.emit(key)

    def _light_row(self, wanted) -> None:
        """Ink one row and no other; None clears every one.

        :param wanted: the row widget to light, or None.
        """
        for row in self.all_rows():
            want = row is wanted
            if bool(row.property("hovered")) == want:
                continue
            row.setProperty("hovered", want)
            style = row.style()
            if style is not None:
                style.unpolish(row)
                style.polish(row)

    def leaveEvent(self, event):                             # noqa: N802
        """No row is lit once the pointer leaves the rail.

        :param event: the leave event, passed on after every row is unlit.
        """
        self._light_row(None)
        super().leaveEvent(event)

    def refresh_icons(self) -> None:
        """Re-ask the one provider for the symbol of every row, views and steps included.

        A view is drawn with the symbol of the module it belongs to, which is the honest
        answer: a view is a part of that module and has no identity of its own to have a
        symbol for.
        """
        if not getattr(self, "_tree_built", False):
            super().refresh_icons()
            return
        from .panes._compat import _set_scaled_icon_size
        from .panes.dock import ICON_PX
        for row in self.all_rows():
            node = getattr(row, "node", None)
            px = ICON_PX if node is None or node.kind in ("module", "home") else ICON_PX - 4
            _set_scaled_icon_size(row, px)
            if self._icon_for is None:
                continue
            icon = self._icon_for(getattr(row, "icon_key", row.key))
            if icon is not None:
                row.setIcon(icon)

    def refresh_visibility(self) -> None:
        """Show a row when every heading above it is open, and redraw the fold marks.

        The one place the tree decides what is on screen. A blocked module is still shown:
        Einar asked for the tiles to be an honest inventory rather than a menu of
        promises, and a rail that hid what is not built yet would be the opposite.
        """
        if not getattr(self, "_tree_built", False):
            super().refresh_visibility()
            return
        for node, _depth, above, widget in self._entries:
            widget.setVisible(all(k in self._open for k in above))
            if isinstance(widget, SectionHeader):
                widget.set_open(self.fold_key(node) in self._open)
                widget.style().unpolish(widget)
                widget.style().polish(widget)
            elif node.folds:
                widget.setText(self._row_text(node))
        self.setFixedWidth(self.column_width())

    def apply_theme(self) -> None:
        """Paint the rail: the dock panel of spaCR, plus one indent rule per level.

        The vendored sheet is applied first and added to rather than replaced, so a change
        to how spaCR paints its rail reaches this one too.
        """
        super().apply_theme()
        if not getattr(self, "_tree_built", False):
            return
        from .panes._compat import active_palette, scaled_px
        palette = active_palette()
        rules = [
            f'QLabel#SidebarSection[treeLevel="1"] {{ padding-left: {scaled_px(INDENT[1])}px; '
            "font-weight: 400; }",
        ]
        for level, pad in enumerate(INDENT):
            rules.append(f'QPushButton#SidebarItem[treeLevel="{level}"] '
                         f"{{ padding-left: {scaled_px(pad)}px; }}")
        rules.append('QPushButton#SidebarItem[treeKind="view"], '
                     'QPushButton#SidebarItem[treeKind="step"] '
                     f'{{ color: {palette["fg_muted"]}; }}')
        rules.append('QPushButton#SidebarItem[treeBlocked="true"] '
                     f'{{ color: {palette["fg_muted"]}; }}')
        self.setStyleSheet(self.styleSheet() + "\n" + "\n".join(rules))

    # ---------------------------------------------------------------------- reading
    def all_rows(self) -> list:
        """Every clickable row of the tree: the modules, Home, the views and the steps.

        Separate from :meth:`rows`, which is the modules and Home alone. Both exist
        because two different questions are asked of the rail -- "what modules can be
        opened from here", which the window tests ask, and "paint, measure and unlight
        everything", which the widget asks -- and answering both with one list made the
        first one unanswerable.
        """
        return list(self._rows) + list(self._extra_rows)

    def rows(self) -> list:
        """Every module row and the Home row, in drawing order.

        Deliberately NOT the view and step rows: this answers "what modules can be opened
        from the rail", which is the question the window tests ask of it.
        """
        return list(self._rows)

    def view_rows(self) -> list:
        """Every declared view as a row, in drawing order."""
        return [r for r in self._extra_rows if getattr(r, "node", None) is not None
                and r.node.kind == "view"]

    def step_rows(self) -> list:
        """Every next-step row, in drawing order."""
        return [r for r in self._extra_rows if getattr(r, "node", None) is not None
                and r.node.kind == "step"]

    def visible_rows(self) -> list:
        """Every row the rail is showing right now, headings included."""
        return [w for _n, _d, _a, w in self._entries if not w.isHidden()]

    def fitting_width(self) -> int:
        """Width that shows the longest visible name in full, within the dock bounds.

        Two things the dock does not have to allow for, both found by rendering the rail
        and reading it rather than by a test:

        * the VIEW and STEP rows, which the dock has none of. A rail measured on the
          module names alone is sized for "Gallery" and elides "Condition comparison" the
          moment a module unfolds.
        * the INDENT. A row three levels in starts 38 px from the edge, and a width that
          ignores that is 38 px short of the text it was measuring -- which is what drew
          "Published expressi…", "Other measureme…" and "Sparse measureme…" in a rail
          wide enough for every one of them.
        """
        from PySide6.QtGui import QFontMetrics
        from .panes._compat import scaled_px
        from .panes.dock import ICON_PX
        widest = 0
        for row in self.all_rows():
            if row.isHidden():
                continue
            metrics = QFontMetrics(row.font())
            text = getattr(row, "full_text", lambda: row.text())()
            indent = INDENT[min(int(row.property("treeLevel") or 0), len(INDENT) - 1)]
            widest = max(widest, metrics.horizontalAdvance(str(text)) + scaled_px(indent))
        room = widest + scaled_px(ICON_PX) + scaled_px(44)
        return max(scaled_px(self.WIDTH_MIN), min(room, scaled_px(self.WIDTH_MAX)))

    def clipped_items(self) -> list:
        """Rows whose name had to be shortened to fit."""
        return [r for r in self.all_rows() if r.is_elided()]
