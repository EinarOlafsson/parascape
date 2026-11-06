"""What the vendored spaCR pane widgets ask of the rest of spaCR, answered by ParaScape.

The four files beside this one -- `collapsible_splitter`, `collapsible_section`,
`foldable` and `collapse_arrow` -- are spaCR's own (`spacr/qt/widgets/`), copied with
only their imports pointed here. They reached into five spaCR modules:

``gil_priority``  an event hub that keeps Python filters off paint events and tells a
                  watcher about shows and hides afterwards. Here a small relay per
                  watched object does the second part.
``i18n``          translation. ParaScape is English only; `tr` returns its text.
``preferences``   the QSettings store and the "which panels are folded" record. The
                  same keys, in ParaScape's own settings file.
``gui_scale``     spaCR's whole-interface zoom. ParaScape has none, so the scale is 1.
``theme`` / ``widgets.ambient``  the painting palette and the ripple a pane sends
                  through the animated field when it folds; both answered from
                  `parascape.backdrops`, which is the same spaCR animation.

Licensed under spaCR's terms -- LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt.
"""
from __future__ import annotations

import json

from PySide6.QtCore import QEvent, QObject, QSettings, Qt, QTimer
from shiboken6 import isValid

_KEY_FOLDED = "ui/folded_panels"


# ------------------------------------------------------------------ gil_priority
#: Events a watcher is told about on the next turn of the event loop rather than at once,
#: as spaCR's hub tells them ("a move, resize, show or hide ... is told afterwards"). A
#: watcher is then never run in the middle of a pane being deleted, and is not run at all
#: once either side is gone.
_LATER = frozenset((QEvent.Type.Show, QEvent.Type.Hide, QEvent.Type.ShowToParent,
                    QEvent.Type.HideToParent, QEvent.Type.LayoutRequest))


def _deliver(watcher, watched, kind) -> None:
    if isValid(watcher) and isValid(watched):
        watcher.eventFilter(watched, QEvent(kind))


class _Relay(QObject):
    """Stands between a watched object and its watcher, as spaCR's event hub does."""

    def __init__(self, watched, watcher, kinds):
        super().__init__(watched)                     # goes when the watched object goes
        self.watcher = watcher
        self.kinds = frozenset(kinds) if kinds else None

    def eventFilter(self, watched, event):            # noqa: N802 - Qt name
        kind = event.type()
        if self.kinds is not None and kind not in self.kinds:
            return False
        watcher = self.watcher
        if not isValid(watcher):
            return False
        if kind in _LATER:
            QTimer.singleShot(0, watcher, lambda: _deliver(watcher, watched, kind))
            return False
        return bool(watcher.eventFilter(watched, event))


def _watch_object_events(watched, watcher, kinds=(), replay=True) -> bool:
    """Tell `watcher.eventFilter` about `watched`'s events of `kinds` (spaCR's API)."""
    if watched is None or watcher is None:
        return False
    relay = _Relay(watched, watcher, kinds)
    watched.installEventFilter(relay)
    return True


def _stop_watching_object_events(watched, watcher) -> bool:
    """Undo `_watch_object_events`. True if it was watching."""
    if watched is None or watcher is None or not isValid(watched):
        return False
    found = False
    for relay in watched.findChildren(_Relay, options=Qt.FindChildOption.FindDirectChildrenOnly):
        if relay.watcher is watcher:
            watched.removeEventFilter(relay)
            relay.deleteLater()
            found = True
    return found


# ------------------------------------------------------------------ i18n
def tr(text, /, language=None, **_kwargs) -> str:
    """ParaScape is English only: the text itself."""
    return str(text)


# ------------------------------------------------------------------ preferences
def _settings() -> QSettings:
    """ParaScape's settings store; the backdrop shim owns the organisation name."""
    from ..backdrops import _compat as AC
    return AC._settings()


def get_folded_panels() -> dict:
    """``{"<module>/<pane>": folded}`` for every pane left out of its default state."""
    raw = _settings().value(_KEY_FOLDED, "")
    if not raw:
        return {}
    try:
        value = json.loads(raw)
    except (TypeError, ValueError):
        return {}
    return {str(k): bool(v) for k, v in value.items()} if isinstance(value, dict) else {}


def set_folded_panel(key: str, shut: bool, *, default_shut: bool = False) -> None:
    """Remember a pane folded or opened; a pane back in its default state is forgotten."""
    key = str(key or "").strip()
    if not key:
        return
    state = get_folded_panels()
    if bool(shut) == bool(default_shut):
        state.pop(key, None)
    else:
        state[key] = bool(shut)
    _settings().setValue(_KEY_FOLDED, json.dumps(state))


# ------------------------------------------------------------------ gui_scale
QWIDGETSIZE_MAX = (1 << 24) - 1


def add_listener(_callback) -> None:
    """No interface zoom in ParaScape, so nothing ever changes scale."""


def scale_int(value, factor=None) -> int:
    return int(value)


def _factor_for(_widget) -> float:
    return 1.0


def _actual_sizes(splitter) -> list:
    return [int(v) for v in splitter.sizes()]


def _set_actual_sizes(splitter, sizes) -> None:
    splitter.setSizes([max(0, int(v)) for v in sizes])


def _set_actual_dim(widget, key: str, value: int) -> None:
    setter = {"minh": "setMinimumHeight", "minw": "setMinimumWidth",
              "maxh": "setMaximumHeight", "maxw": "setMaximumWidth"}[key]
    getattr(widget, setter)(max(0, min(QWIDGETSIZE_MAX, int(value))))


# ------------------------------------------------------------------ theme, ambient
def _painting_palette() -> dict:
    """The colours of the ParaScape theme on screen."""
    from ..backdrops import _compat as AC
    return AC.palette_for(AC._painting_theme())


def field_ripple_for_widget(widget, edge=None, rect=None, strength=1.0) -> None:
    """Send a ripple through the animated field, when there is one. Never raises."""
    try:
        from ..backdrops.engines import field_ripple_for_widget as ripple
        ripple(widget, edge=edge, rect=rect, strength=strength)
    except Exception:                                 # decoration must never break a fold
        pass
