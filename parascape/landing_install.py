"""Lay the landing view over the whole main window, as spaCR's home is.

Einar's direction: the home screen should be spaCR's, and the UMAP and the parasite
drawing -- "beautiful, and I don't want to discard them" -- should be part of a module
rather than the thing every launch opens on.

The landing is an OVERLAY: a child of the main window covering everything between the
menu bar and the status bar. Nothing under it changes. Two earlier designs did change
things, and both failed in the same way:

1. A central stack (landing in front, map behind) left every dock in place, so the
   landing got a 470 px strip between dock columns, beside the old node home dock.
2. Hiding the docks while the landing showed fixed the look and broke the meaning: a
   hidden dock leaves `tabifiedDockWidgets` and its children stop being `isVisibleTo`
   the window. The full suite caught three layout invariants failing, and the same
   change would have left any panel that skips work while hidden stale on return.

Covering is not hiding. The docks keep their tabbing, visibility and state; the map
stays the central widget; the landing simply sits in front until something is opened.
Anything that brings a panel forward -- a tile, a Navigator result, a menu command, the
search field, a dock raised by any code -- takes the overlay away. View ▸ Home (Alt+Home)
puts it back.

Genes and index entries go through `home_panel.open_node`, because they genuinely are
nodes and it already does the subtle part -- selecting a gene before showing its card,
since the other order flashes the previous gene's card. Modules are not nodes (a node's
kind is gene/question/strategy/dataset/panel/guide and its category gene/list/label, and
a module is none of those), so they get `raise_dock`.
"""
from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from . import modules as MOD
from . import nodes as ND


def _counts(window) -> dict:
    """What the Measured panel shows, read from the window rather than hard-coded."""
    out = {}
    table = getattr(window, "nodes", None)
    if table is not None and "gene_id" in getattr(table, "columns", ()):
        try:
            species = getattr(window, "species", "") or ""
            label = "P. falciparum genes" if "falc" in str(species).lower() else "T. gondii genes"
        except Exception:
            label = "Genes"
        out[label] = (f"{len(table):,}", "genes in the map on screen")
    try:
        from . import datasets as DS
        n = len(DS.registry()) if hasattr(DS, "registry") else None
        if n:
            out["Datasets"] = (f"{n:,}", "registered sources and computed layers")
    except Exception:
        pass
    return out


def _reasons(window) -> dict:
    """Which modules the window cannot actually open, by the docks it really built."""
    docks = tuple(d.windowTitle() for d in window.findChildren(QtWidgets.QDockWidget))
    out = {}
    for module in MOD.MODULES:
        reason = MOD.unavailable_reason(module, docks=docks)
        if reason:
            out[module.key] = reason
    return out


def _subject_in_hand(landing) -> tuple:
    """What the Navigator is holding, as ``(kind, value)``; ``("", "")`` when nothing.

    Asked of the landing at the moment a module is opened rather than carried on the
    signal, because the module tiles open modules too and they hold no subject -- a second
    signal for the same event would be a second thing to keep in step.
    """
    navigator = getattr(landing, "navigator", None)
    return getattr(navigator, "subject", ("", "")) if navigator is not None else ("", "")


def _takes_a_subject(shell_module) -> bool:
    """Whether a module's `open_in` accepts a subject as well as the window.

    Read from the signature rather than listed here: a list of which shells take one would
    be a table of the same facts as the code, kept in step by hand, which is the one thing
    instruction 98 is about.
    """
    import inspect
    try:
        return len(inspect.signature(shell_module.open_in).parameters) > 1
    except (TypeError, ValueError):
        return False


def carry_subject(window, subject: str, value: str) -> bool:
    """Put a subject into the window before a module is opened on it. Returns whether it took.

    This is the pre-filling half of an offer (instruction 98, WP8): a next step that opened
    an empty module would leave the gene to be typed a second time, which is the whole thing
    spaCR's chaining strip exists to stop.

    One gene is carried by SELECTING it in the window -- the same route the search field
    uses -- because that is what every panel already reads: the evidence dock, the
    information map, the compare tabs and the find panel all follow the window selection.
    A module built on its own shell is additionally handed the gene by name, below.

    A LIST is carried by GATING it (instruction 98, WP9). The gate is what the application
    already means by "this set of genes is the one in hand": the map marks it, `describe_gate`
    reports what is in it, the strategy panel has a button that fills a gene box from it, and
    File ▸ Export gated selection writes it. Inventing a second window-wide notion of a
    chosen list beside that one would be two tables of the same fact, which is what this
    instruction exists to remove.

    :param window: the main window.
    :param subject: a subject from :mod:`parascape.ports`; ``gene`` and ``gene_list`` are
        carried, and anything else is reported as not taken rather than silently ignored.
    :param value: what is in it. A list is identifiers separated by commas.
    :returns: True when the window was moved onto the subject.
    """
    from . import ports as PORTS
    if not value:
        return False
    if subject == PORTS.GENE_LIST:
        return _carry_list(window, value)
    if subject != PORTS.GENE:
        return False
    select = getattr(window, "_workflow_gene", None)
    if not callable(select):
        return False
    try:
        select(str(value))
    except Exception:                   # a selection that fails must not stop the module
        return False
    return True


def _carry_list(window, value: str) -> bool:
    """Gate the genes of a list in the window. Returns whether the gate took.

    :param window: the main window.
    :param value: the identifiers, separated by commas.
    """
    nodes = getattr(window, "nodes", None)
    if nodes is None or "gene_id" not in getattr(nodes, "columns", ()):
        return False
    from . import list_reading as LR
    wanted = {g.lower() for g in LR.split_identifiers(value)}
    if not wanted:
        return False
    try:
        import numpy as np
        ids = nodes["gene_id"].astype(str).str.lower().to_numpy()
        index = np.flatnonzero(np.isin(ids, np.array(sorted(wanted))))
        if not len(index):
            return False
        window.on_gated(index)
    except Exception:                   # a gate that fails must not stop the module
        return False
    return True


def open_module(window, key: str, subject: str = "", value: str = "") -> str:
    """Open one module in this window. Returns what the status bar was told.

    Which Python module builds which key is read from the catalogue
    (`modules.Module.builder`) rather than from a dict kept here. The dict was
    a sixth table of the same facts, and a module added to `modules.py` without
    a line added here simply opened nothing.

    :param window: the main window.
    :param key: the module key.
    :param subject: the kind of thing being carried in, e.g. ``"gene"``. Empty opens the
        module on whatever it was already showing, which is what every caller before
        instruction 98 wanted.
    :param value: what is in the subject -- the gene identifier.
    """
    module = MOD.get(key)
    built_by = MOD.shell_module(key)
    carried = carry_subject(window, subject, value)
    message = None
    if built_by and getattr(window, built_by, None) is not None:
        shell_module = __import__(f"{__package__}.{built_by}", fromlist=["open_in"])
        # Only the shells that READ a subject take one; `map_module.open_in` and
        # `infer_module.open_in` take the window alone, and passing an argument they do
        # not accept would turn an offer into a TypeError in front of Einar.
        #
        # On the VALUE rather than on `carried`, which is what WP8 asked: a window that
        # cannot hold the subject itself is no reason to open the module empty. The gene
        # card handed a list has the list even where `_workflow_gene` is missing, and a
        # module opened empty leaves the thing to be typed a second time, which is what
        # the whole offer exists to stop.
        if value and _takes_a_subject(shell_module):
            message = shell_module.open_in(window, str(value))
        else:
            message = shell_module.open_in(window)
        if key == "partners":
            # Partners is the network read by binding evidence: the physical links ticked.
            window.network_module.module.show_partners()
            message = "Partners"
    elif key in getattr(window, "tool_modules", {}):
        from .dock_module import open_in as open_tool
        message = open_tool(window, key)
    if message is not None:
        try:
            window.statusBar().showMessage(message, 6000)
        except Exception:
            pass
        return message
    if module.opens and not MOD.unavailable_reason(module):
        leave_landing(window)
    if key == "map":
        message = "Map"
    elif not module.opens:
        message = f"{module.name}: {MOD.unavailable_reason(module)}"
    elif module.opens.startswith("action:"):
        name = module.opens.split(":", 1)[1]
        # Export handed a LIST writes the list, not the whole visible map (instruction 98,
        # WP9). Export declares that it reads a list of genes, and `export_visible` -- the
        # action on its catalogue row -- writes every gene passing the class filter, so the
        # offer "Export, with your 40 genes already in it" wrote eight thousand of them.
        # `carry_subject` has just gated the list, and `export_gated` is the action that
        # writes a gated set: a separate action on purpose, because the two answer different
        # questions. So the subject in hand decides which of the two is meant.
        from . import ports as PORTS
        if (name == "export_visible" and subject == PORTS.GENE_LIST and carried
                and callable(getattr(window, "export_gated", None))):
            name = "export_gated"
        method = getattr(window, name, None)
        if not callable(method):
            message = f"{module.name} cannot be opened in this window"
        else:
            try:
                method()
                message = module.name
            except Exception as error:      # a broken action must not take the landing
                message = f"could not open {module.name}: {type(error).__name__}: {error}"
    else:
        dock, _, tab = module.opens.partition(":")
        message = raise_dock(window, dock, tab)
    try:
        window.statusBar().showMessage(message, 6000)
    except Exception:
        pass
    return message


def _say(window, message: str) -> str:
    """Put one line in the status bar, and never fail for want of a status bar.

    :param window: the main window.
    :param message: what to say.
    :returns: the message, so a caller can report it.
    """
    try:
        window.statusBar().showMessage(message, 6000)
    except Exception:                       # a window without a bar must still work
        pass
    return message


def open_view(window, key: str, name: str) -> str:
    """Open one declared view of a module. Returns what the status bar was told.

    WHY A VIEW IS OPENED THROUGH ITS MODULE (instruction 98, WP10). The inventory found
    29 tabs and four windows the navigator could not name. Each is declared on the module
    whose panel holds it (`modules.View`), and this is the one route to them, so a view
    cannot be reachable without being declared and cannot be declared without being
    reachable -- which is what `tests/test_views.py` holds in both directions.

    A view that lives in the dock its module borrows is opened IN the module, with its tab
    brought forward, because that is the shell Einar asked every module to look like. A
    view in some other dock, or in a window of its own, leaves the landing and opens there:
    honest about being a panel rather than pretending to be part of the module.

    :param window: the main window.
    :param key: the key of the module the view belongs to.
    :param name: the name of the view, as declared.
    :returns: a line for the status bar; never raises, because this is called from a click.
    """
    try:
        module = MOD.get(key)
    except KeyError:
        return f"there is no {key} module"
    view = next((v for v in module.views if v.name.lower() == name.strip().lower()), None)
    if view is None:
        return f"{module.name} has no view called {name}"
    if view.opens.startswith("action:"):
        method = getattr(window, view.opens.split(":", 1)[1], None)
        if not callable(method):
            return f"{view.name} cannot be opened in this window"
        leave_landing(window)
        try:
            method()
        except Exception as error:               # a broken view must not take the window
            return f"could not open {view.name}: {type(error).__name__}: {error}"
        return view.name
    dock, _, tab = view.opens.partition(":")
    inside = module.opens.partition(":")[0]
    if dock == inside and key in getattr(window, "tool_modules", {}):
        # The module already shows this dock: open the module and turn to the tab, so the
        # view arrives inside the shell rather than as a bare panel beside the map.
        from .dock_module import open_in as open_tool
        message = open_tool(window, key)
        if _turn_to_tab(window, key, dock, tab):
            return f"{message} ▸ {view.name}"
        return message
    return raise_dock(window, dock, tab)


def _turn_to_tab(window, key: str, dock: str, tab: str) -> bool:
    """Bring a tab forward in the panel a tool module borrows. Did it work?

    The panel is LENT: while the module is on screen it lives under the shell, and the rest
    of the time under its dock. It is the SAME object either way, so this asks the module
    for it and falls back to the dock -- which matters because a module takes its panel on
    the show event, and a module opened in the same breath as a tab being chosen has not
    had one yet.
    """
    if not tab:
        return False
    shell = getattr(window, "tool_modules", {}).get(key)
    borrowed = getattr(getattr(shell, "module", None), "main", None)
    panel = getattr(borrowed, "panel", None)
    if panel is None:
        holder = next((d for d in window.findChildren(QtWidgets.QDockWidget)
                       if d.windowTitle() == dock), None)
        panel = holder.widget() if holder is not None else None
    if panel is None:
        return False
    # The Strategies tabs sit BEHIND the card of the chosen strategy, in a stack. Asking the
    # panel to show its details first is what the search field has always done; without it
    # the tab is selected under a card nobody has moved out of the way.
    if hasattr(panel, "show_details"):
        panel.show_details()
    for tabs in [panel] + panel.findChildren(QtWidgets.QTabWidget):
        if not isinstance(tabs, QtWidgets.QTabWidget):
            continue
        for i in range(tabs.count()):
            if tabs.tabText(i) == tab:
                tabs.setCurrentIndex(i)
                return True
    return False


def raise_dock(window, title: str, tab: str = "") -> str:
    """Show a dock, bring it forward, and select a tab in it. Reports, never raises:
    a dock that failed to install must leave the landing usable."""
    dock = next((d for d in window.findChildren(QtWidgets.QDockWidget)
                 if d.windowTitle() == title), None)
    if dock is None:
        return f"the {title} panel is not available in this window"
    dock.show()
    dock.raise_()
    if not tab:
        return title
    tabs = dock.findChild(QtWidgets.QTabWidget)
    if tabs is not None:
        for i in range(tabs.count()):
            if tabs.tabText(i) == tab:
                tabs.setCurrentIndex(i)
                return f"{title} ▸ {tab}"
    return f"{title}: {tab} is not a tab here"


def open_item(window, kind: str, key: str) -> str:
    """Open a gene, question, method, dataset, panel or guide the Navigator found.

    A GENE OPENS THE GENE CARD. Einar, 2026-11-12: "now when i put in a gene and pressed
    the modual that came up for the gene, i was brought to the old version of parascape.
    this also needs to be overhauled." It went to `home_panel.open_node`, which took the
    landing away and raised the evidence dock -- the pre-navigator window. It now selects
    the gene in the window, exactly as before, so every panel behind the card is reading
    the same gene, and then shows the gene card over it.

    A LIST OF GENES OPENS THE GENE LIST CARD. Einar, 2026-11-12: "if i type in a list of
    genes, the same logic but this time not gene card instead a gene list module". The kind
    is ``gene_list`` and the key is the identifiers separated by commas, which is what the
    Navigator resolved them to, so nothing is typed a second time.
    """
    from .home_panel import _index, open_node
    if kind == "gene":
        from . import ports as PORTS
        return open_module(window, "gene", subject=PORTS.GENE, value=key)
    if kind == "gene_list":
        from . import ports as PORTS
        return open_module(window, "list", subject=PORTS.GENE_LIST, value=key)
    if kind == "view":
        # One of the 32 declared views, named in the Navigator (instruction 98, WP4). The
        # key is the view name, which the catalogue keeps unique, and `find_view` answers
        # which module it belongs to -- so the Navigator does not have to carry two names
        # through a signal that takes one.
        try:
            module, view = MOD.find_view(key)
        except KeyError:
            return _say(window, f"{key}: there is no view by that name")
        return _say(window, open_view(window, module.key, view.name))
    # Everything else opens in a dock, so the overlay has to go: a panel raised behind the
    # landing is one nobody sees.
    leave_landing(window)
    entry = next((e for e in _index(window)
                  if getattr(e, "title", None) == key and getattr(e, "kind", None) == kind), None)
    if entry is None:
        entry = next((e for e in _index(window) if getattr(e, "title", None) == key), None)
    if entry is None:
        message = f"{key}: not found in this window's index"
        window.statusBar().showMessage(message, 6000)
        return message
    return open_node(window, ND.entry_node(entry))


def _docks(window):
    return [d for d in window.findChildren(QtWidgets.QDockWidget)]


class _Cover(QtCore.QObject):
    """Keeps the overlay covering the window's content area as the window resizes."""

    def __init__(self, window):
        super().__init__(window)
        self._window = window

    def eventFilter(self, obj, event):                    # noqa: N802
        """Refit on a resize or show of the window; never consume the event."""
        if obj is self._window and event.type() in (QtCore.QEvent.Type.Resize,
                                                    QtCore.QEvent.Type.Show):
            fit_landing(self._window)
            _arm_later(self._window)
        return False


ARM_AFTER_MS = 400


def _arm_later(window) -> None:
    """Disarm the dock watcher, and arm it once the window has settled.

    When the window first appears every dock reports that it became visible. Read as "the
    person opened a panel", that would take the landing away the instant the application
    opened. No test catches it -- the suite never shows a window -- so it is guarded here:
    visibility changes while the window is still being put on screen are ignored.

    A DEADLINE READ WHEN THE SIGNAL ARRIVES, NOT A TIMER THAT ARMS. This was
    ``QTimer.singleShot(400, ...)``, and the two were in a race: the platform delivers
    those dock visibility changes asynchronously, so whether the landing survived a show
    depended on whether the timer fired before they arrived. Instruction 98 records the
    failure (`test_the_dock_watcher_ignores_the_window_appearing`) being reported as a
    possible defect and passing on committed code, which is what a race looks like. WP7
    settled it: the rail became the tree, the window got slower to put on screen by a
    tenth of a second, and the race turned into a reliable failure.

    Read as a deadline there is nothing to race. Every show AND every resize pushes it
    forward, because a window being put on screen resizes several times while it settles,
    and the question being asked is "has anything stopped happening yet".

    :param window: the main window.
    """
    import time
    window._landing_quiet_until = time.monotonic() + ARM_AFTER_MS / 1000.0
    # Kept because it is what the attribute was called, and because reading it is how
    # anything outside asks whether the watcher is awake.
    window._landing_armed = False
    QtCore.QTimer.singleShot(ARM_AFTER_MS, window,
                             lambda: setattr(window, "_landing_armed", True))


def _watcher_awake(window) -> bool:
    """Whether a dock becoming visible should now be read as the person opening a panel.

    :param window: the main window.
    """
    import time
    return time.monotonic() >= getattr(window, "_landing_quiet_until", 0.0)


def fit_landing(window) -> None:
    """Cover everything between the menu bar and the status bar -- docks included."""
    landing = getattr(window, "landing", None)
    if landing is None:
        return
    # isHidden, not isVisible: isVisible is False for every child until the window itself
    # is on screen, so the first fit -- made at install, before the window shows -- put the
    # overlay over the menu bar. isHidden asks whether something deliberately hid the bar.
    menu = window.menuBar()
    top = menu.geometry().bottom() + 1 if not menu.isHidden() else 0
    bar = window.statusBar()
    bottom = bar.geometry().top() if bar is not None and not bar.isHidden() else window.height()
    landing.setGeometry(0, top, window.width(), max(0, bottom - top))
    if landing.isVisible():
        landing.raise_()


def landing_up(window) -> bool:
    """Whether the landing is covering the window right now."""
    landing = getattr(window, "landing", None)
    return landing is not None and not landing.isHidden()


def enter_landing(window) -> None:
    """Cover the window with the landing. Nothing underneath changes."""
    landing = window.landing
    fit_landing(window)
    # Einar: "navigator should be reset every time I go back to it".
    landing.navigator.reset()
    landing.show()
    landing.raise_()
    if landing.page == "Navigator":
        landing.navigator.field.setFocus()


def leave_landing(window, show=None, keep=None) -> None:
    """Take the landing away, revealing the window as it was.

    `keep` is a dock that was just brought forward; it is raised again after the overlay
    goes, so it is the tab on top. `show` is accepted for callers that name a view; the
    map is the central widget and is already there underneath.
    """
    landing = getattr(window, "landing", None)
    if landing is None or landing.isHidden():
        return
    landing.hide()
    if keep is not None:
        keep.raise_()


def _watch(window) -> None:
    """Take the landing away whenever anything brings a panel forward.

    Only the landing's own tiles knew to leave it at first, so a panel opened from a menu
    or the search field appeared UNDER the overlay. Three signals cover every route: a
    dock becoming visible (raising a tab that was not on top), any menu command other than
    Home and the animation choices, and the search field opening a result.
    """
    for dock in _docks(window):
        dock.visibilityChanged.connect(
            lambda visible, d=dock: leave_landing(window, keep=d)
            if visible and landing_up(window) and _watcher_awake(window)
            else None)

    stay = set()

    def on_menu(action):
        if action in stay or action.parent() in (getattr(window, "landing_animation_menu", None),
                                                 getattr(window, "landing_opacity_menu", None)):
            return
        leave_landing(window)

    for top in window.menuBar().actions():
        if top.menu() is not None:
            top.menu().triggered.connect(on_menu)
    window._landing_stay_actions = stay

    try:
        from . import help_search as HS
        field = HS.field_of(window)
        if field is not None:
            field.opened.connect(lambda *_: leave_landing(window))
    except Exception:                                     # the overlay works without it
        pass


def go_home(window) -> None:
    """Back to the landing, from anywhere -- and off a module shell onto the home page."""
    landing = getattr(window, "landing", None)
    if landing is not None and landing.page == "Module":
        landing.show_page(landing.home_page())
    enter_landing(window)


def install(window):
    """Lay the landing over the window and wire it up. Returns the landing."""
    from .home_panel import _index, gene_lookup
    from .landing import Landing

    def resolve_list(ids):
        """The identifiers that are genes of this table, as the table spells them.

        Given to the Navigator so a typed LIST can be recognised in one vectorised match
        (instruction 98, WP9). The gene lookup beside it is a SEARCH taking the whole query,
        which asked for "TGME49_1, TGME49_2" looks for one gene whose product contains that
        whole string, and finds nothing.
        """
        from . import list_reading as LR
        return LR.resolve(getattr(window, "nodes", None), ids).genes

    landing = Landing(counts=_counts(window), reasons=_reasons(window),
                      genes=gene_lookup(window), index=lambda: _index(window),
                      resolve_list=resolve_list, parent=window)
    window.landing = landing
    window.landing_stack = None                 # the stack design is gone; see the docstring

    # The subject the Navigator is holding travels with the module it opens, so a module
    # reached by typing a gene opens already reading that gene rather than empty.
    landing.module_chosen.connect(
        lambda key: open_module(window, key, *_subject_in_hand(landing)))
    landing.item_chosen.connect(lambda kind, key: open_item(window, kind, key))
    # One of the 32 declared views, chosen in the tree down the left or named in the
    # Navigator. `open_view` is the one route to a view, so both go through it.
    landing.view_chosen.connect(lambda key, name: _say(window, open_view(window, key, name)))

    home = QtGui.QAction("Home", window)
    # Alt+Home, the browser convention. Ctrl+Shift+H was the first choice and was already
    # Help ▸ Search ParaScape: taking it silently stopped that shortcut focusing the search
    # field. Ctrl+Home was the next obvious one and would steal start-of-line from every
    # text field in the window.
    home.setShortcut(QtGui.QKeySequence("Alt+Home"))
    home.setToolTip("Back to the Navigator and the module tiles, from anywhere. Everything "
                    "underneath stays exactly as you left it.")
    home.setStatusTip("Return to the Navigator and the module tiles")
    home.triggered.connect(lambda: go_home(window))
    window.addAction(home)
    window.landing_home_action = home

    navigator = QtGui.QAction("Navigator", window, checkable=True)
    navigator.setChecked(landing.navigator_enabled())
    navigator.setToolTip("Show the Navigator -- one field for a gene, a question, a module or a "
                         "dataset -- as the first page of Home. Off, Home is the module tiles alone.")
    navigator.setStatusTip("Turn the Navigator on or off")

    def toggle(on):
        landing.set_navigator_enabled(on)
        enter_landing(window)
    navigator.toggled.connect(toggle)
    window.landing_navigator_action = navigator

    menu = _view_menu(window)
    if menu is not None:
        menu.addSeparator()
        menu.addAction(home)
        menu.addAction(navigator)
        menu.addMenu(_animation_menu(window, landing))
        menu.addMenu(_opacity_menu(window))
    # Also in Help, which is where spaCR keeps it; Einar looked for it in View.
    for top in window.menuBar().actions():
        if top.text().replace("&", "") == "Help" and top.menu() is not None:
            top.menu().addAction(navigator)

    # The modules that have moved onto spaCR's module shell (instruction 86). Built here,
    # after every dock, because each one borrows the console and the assistant and moves
    # its controls out of the panel they lived in. A failure leaves the old route.
    #
    # WHICH ones is read from the catalogue (`modules.shell_modules`), not listed here:
    # the list and the table had to agree, and Map was built by its own copy of this block
    # purely because it was written first. The window attribute carries the same name as
    # the Python module, which is what `open_module` looks the shell up by.
    #
    # Every attribute is set to None BEFORE any of them is built, because each
    # shell reads the others through `getattr(window, "map_module", None)` and
    # a shell built first would otherwise see an attribute that does not exist
    # rather than one that is not ready. The default makes both read the same,
    # which is exactly the kind of agreement that stops holding the day somebody
    # drops the default.
    for name in MOD.shell_modules():
        setattr(window, name, None)
    for name in MOD.shell_modules():
        try:
            module = __import__(f"{__package__}.{name}", fromlist=["build"])
            shell = module.build(window)
            setattr(window, name, shell)
            if shell is not None:
                # Into the landing at once: a shell left unparented until first opened is a
                # top-level window of its own, deleted out of order with the window it borrows from.
                landing.stack.addWidget(shell)
        except Exception as error:              # the window must open without it
            window.statusBar().showMessage(f"{name.split('_')[0].title()} module unavailable: "
                                           f"{error}", 8000)
    # Every other tile with a panel: the same shell around the panel, lent while on screen.
    window.tool_modules = {}
    try:
        from . import dock_module
        window.tool_modules = dock_module.build_all(window, landing)
    except Exception as error:                  # the window must open without them
        window.statusBar().showMessage(f"tool modules unavailable: {error}", 8000)
    window._landing_cover = _Cover(window)
    window.installEventFilter(window._landing_cover)
    _watch(window)
    _arm_later(window)
    window._landing_stay_actions.add(home)
    window._landing_stay_actions.add(navigator)
    enter_landing(window)
    return landing


#: The container opacities offered, as spaCR's preference slider's common stops.
OPACITY_CHOICES = (30, 45, 60, 75, 90, 100)


def _opacity_menu(window) -> QtWidgets.QMenu:
    """View ▸ Container opacity: how solid every module's containers are, as spaCR's
    pane-opacity preference. Lower lets more of the animated field through."""
    from .backdrops import _compat as AC
    from .shell import ModuleShell
    menu = QtWidgets.QMenu("Container opacity", window)
    menu.setToolTip("How solid the containers in every module are. Lower lets more of the "
                    "animated background show through, as in spaCR.")
    menu.setToolTipsVisible(True)
    group = QtGui.QActionGroup(menu)
    group.setExclusive(True)
    current = int(round(AC.get_pane_opacity() * 100))

    def choose(pct):
        AC.set_pane_opacity(pct / 100)
        for shell in window.findChildren(ModuleShell):
            shell.refresh_style()
    for pct in OPACITY_CHOICES:
        action = QtGui.QAction(f"{pct} %", menu, checkable=True)
        action.setChecked(pct == current)
        action.setToolTip(f"Containers {pct} % solid.")
        action.setStatusTip(f"Set every module's containers to {pct} % opacity")
        action.triggered.connect(lambda _=False, p=pct: choose(p))
        group.addAction(action)
        menu.addAction(action)
    window.landing_opacity_menu = menu
    return menu


def _view_menu(window):
    for action in window.menuBar().actions():
        if action.text().replace("&", "") == "View" and action.menu() is not None:
            return action.menu()
    return None


def _animation_menu(window, landing) -> QtWidgets.QMenu:
    """View ▸ Animation: spaCR's themes, one checked, changed live and remembered."""
    from .backdrops import _compat as AC
    from .backdrops import engines as AE
    menu = QtWidgets.QMenu("Animation", window)
    group = QtGui.QActionGroup(menu)
    group.setExclusive(True)
    current = AC.get_ambient_theme() if AC.get_ambient_enabled() else AE.NO_ANIMATION
    for theme in AE.ANIMATION_CHOICES:
        try:
            label = AE.animation_label(theme)
        except Exception:
            label = theme
        action = QtGui.QAction(label, menu, checkable=True)
        action.setChecked(theme == current)
        try:
            note = str(AE.animation_note(theme) or "").strip()
        except Exception:
            note = ""
        if theme == AE.NO_ANIMATION:
            note = note or ("No animated backdrop. Nothing is drawn and no timer runs, so it "
                            "costs nothing at all.")
        action.setToolTip(note or f"The {label} animated backdrop, from spaCR, behind the "
                                  "Navigator and the module tiles.")
        action.setStatusTip(f"Animated backdrop: {label}")
        action.triggered.connect(lambda _=False, t=theme: landing.set_theme(t))
        group.addAction(action)
        menu.addAction(action)
    window.landing_animation_menu = menu
    return menu
