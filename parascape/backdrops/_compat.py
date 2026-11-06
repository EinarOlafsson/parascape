"""What spaCR's ambient engines expect from the rest of spaCR, supplied by ParaScape.

PORTED FROM spaCR, MODIFIED
===========================
The functions here stand in for parts of ``spacr.qt.preferences``,
``spacr.qt.theme``, ``spacr.qt.gil_priority``, ``spacr.qt.hidpi`` and
``spacr.qt._after``, and several are close transcriptions of spaCR's own bodies
(``field_fade_alpha``, ``field_chrome``, ``relative_luminance``, ``_after``).

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

Why a shim rather than a rewrite. ``engines.py`` is spaCR's ``ambient.py`` -- fifteen
thousand lines of animation -- and it reaches into about a dozen spaCR modules for
roughly sixty names. Nearly all of those are preferences and theme lookups, not
animation. Supplying them here leaves the animation code byte-for-byte spaCR's,
which is the reason ParaScape moved to PySide6: so this could transfer.

Three rules govern what is here:

1. **Preferences are ParaScape's own.** They are read from ParaScape's QSettings
   under ``ambient/...``, so a choice made in ParaScape persists in ParaScape and
   never touches spaCR's settings. Every default is spaCR's shipped default, read
   from the engine module itself where the engine defines it, so the first run
   looks exactly like spaCR.
2. **Performance hooks degrade to the plain Qt call.** spaCR's ``gil_priority``
   routes event filters through an application hub so a paint never re-enters
   Python; here they are ordinary ``installEventFilter`` calls. Correct, and
   slower only on the very largest screens.
3. **Anything spaceout- or fractal-only raises if reached.** ParaScape does not
   offer those themes, so reaching one is a bug worth a loud error rather than a
   silent blank backdrop.
"""
from __future__ import annotations

import contextlib
import math
from functools import lru_cache
from types import SimpleNamespace

from PySide6.QtCore import QObject, QSettings, QTimer
from PySide6.QtGui import QColor

from .. import theme as _TH

# ------------------------------------------------------------------- settings
_ORG, _APP = "ParaScape", "ParaScape"


def _settings() -> QSettings:
    """ParaScape's own settings store. Never spaCR's."""
    return QSettings(_ORG, _APP)


def _as_bool(raw, default: bool) -> bool:
    """A QSettings value as a bool. QSettings returns strings from INI files."""
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, str):
        return raw.strip().lower() in ("1", "true", "yes", "on")
    return bool(raw)


def _float(key: str, default: float, lo: float = -math.inf, hi: float = math.inf) -> float:
    try:
        value = float(_settings().value(key, default))
    except (TypeError, ValueError):
        return default
    return max(lo, min(hi, value)) if math.isfinite(value) else default


def _engines():
    """The engine module, imported late: it imports this one at load."""
    from . import engines
    return engines


#: The preference keys, all under ``ambient/`` in ParaScape's own store.
_KEY_SPACEOUT_ANIMATION = "ambient/spaceout_animation"
_KEY_THEME = "ambient/theme"
_KEY_PALETTE = "ambient/palette"
_KEY_ENABLED = "ambient/enabled"


# --------------------------------------------------------------- preferences
def get_ambient_theme() -> str:
    """The chosen theme, or spaCR's default (the spaCR field) when unset or unknown."""
    E = _engines()
    raw = str(_settings().value(_KEY_THEME, E.DEFAULT_THEME) or E.DEFAULT_THEME)
    return raw if raw in E.ANIMATION_CHOICES else E.DEFAULT_THEME


def get_ambient_palette() -> str:
    """The chosen palette, falling back to the theme's default palette."""
    E = _engines()
    theme = get_ambient_theme()
    fallback = E.default_palette_for(theme) if theme in E.AMBIENT_THEMES else E.DEFAULT_PALETTE
    raw = str(_settings().value(_KEY_PALETTE, fallback) or fallback)
    try:
        return raw if raw in E.palettes_for(theme) else fallback
    except Exception:
        return fallback


def get_ambient_enabled() -> bool:
    """Off when the user chose No animation, or switched it off."""
    if get_ambient_theme() == _engines().NO_ANIMATION:
        return False
    return _as_bool(_settings().value(_KEY_ENABLED, True), True)


def get_ambient_speed() -> float:
    return _float("ambient/speed", _engines().DEFAULT_SPEED, 0.0, 10.0)


def get_ambient_size() -> float:
    return _float("ambient/size", _engines().DEFAULT_SIZE, 0.0, 10.0)


def get_ambient_resolution() -> float:
    return _float("ambient/resolution", _engines().DEFAULT_RESOLUTION, 0.05, 4.0)


def get_ambient_density(theme=None) -> float:
    E = _engines()
    default = E._default_density_for(get_ambient_theme() if theme is None else theme)
    return _float("ambient/density", default, 0.0, 10.0)


def get_ambient_drift_direction() -> str:
    E = _engines()
    raw = str(_settings().value("ambient/drift_direction", E.DEFAULT_DRIFT_DIRECTION))
    return raw if raw in getattr(E, "DRIFT_DIRECTIONS", (raw,)) else E.DEFAULT_DRIFT_DIRECTION


def get_reduce_motion() -> bool:
    """Off by default, as in spaCR. On, the backdrop holds still."""
    return _as_bool(_settings().value("ambient/reduce_motion", False), False)


def _hold_still_when_motion_is_reduced() -> bool:
    return False


def get_field_fade_enabled() -> bool:
    return _as_bool(_settings().value("ambient/field_fade", True), True)


def _ambient_blink_percent() -> float:
    return _float("ambient/blink_percent", _engines().DEFAULT_BLINK_PERCENT, 0.0, 100.0)


def _field_ripple_intensity() -> float:
    """Zero is off, one is spaCR's default."""
    return _float("ambient/field_ripples", 1.0, 0.0, 2.0)


def _field_ripples_enabled() -> bool:
    return _field_ripple_intensity() > 0.0


def _field_click_ripples() -> bool:
    return _as_bool(_settings().value("ambient/click_ripples", True), True)


def _field_popup_wave_frequency() -> float:
    return _float("ambient/popup_wave_frequency", 0.0, 0.0, 60.0)


def _ambient_gravity_radius() -> float:
    return 0.0


def _dynamic_animation_enabled() -> bool:
    """Whether running analysis may pause the animation. On, as in spaCR."""
    return _as_bool(_settings().value("ambient/dynamic", True), True)


#: spaCR requests the GPU path by default and falls back when moderngl or a
#: headless context is missing. ParaScape asks for "auto" so a machine without
#: moderngl gets the CPU path without a failed attempt logged on every launch.
GRAPHICS_BACKENDS = ("cpu", "auto", "gpu")


def _flow_graphics_backend() -> str:
    raw = str(_settings().value("ambient/graphics_backend", "auto"))
    return raw if raw in GRAPHICS_BACKENDS else "auto"


_SPACEOUT_FIELD_EFFECT_KEYS = (
    "attractors", "relaxation", "elastic_release", "vortex",
    "density_pulses", "density_waves", "color_waves", "spirals",
)


def _spaceout_field_effects() -> dict:
    return {k: True for k in _SPACEOUT_FIELD_EFFECT_KEYS}


def _ambient_background_choice():
    """None means "follow the page colour", which is what ParaScape always wants."""
    return None


def _effective_ambient_background():
    return QColor(active_page_colour())


def _ambient_custom_colors():
    return ()


def get_fractal_settings() -> dict:
    raise RuntimeError("get_fractal_settings: not part of ParaScape (the fractal backdrop is spaCR-only)")


# --------------------------------------------------------------------- theme
#: ParaScape's palette for the animation. spaCR's theme names do not exist here,
#: so anything spaCR would have asked about a named theme is answered from the
#: one ParaScape is showing.
_PAINT_THEME = ["dark"]


def set_painting_theme(name: str) -> None:
    """Tell the animation which ParaScape theme is on screen."""
    _PAINT_THEME[0] = name if name in _TH.THEMES else "dark"


def _painting_theme() -> str:
    return _PAINT_THEME[0]


def _forget_painting_theme() -> None:
    pass


def palette_for(theme: str = "dark") -> dict:
    return _TH.palette_for(theme if theme in _TH.THEMES else _painting_theme())


def page_colour(theme: str = "dark") -> str:
    """The flat colour behind the panels. spaCR's rule: `page`, else `bg`."""
    p = palette_for(theme)
    return p.get("page") or p["bg"]


def active_page_colour() -> str:
    return page_colour(_painting_theme())


@lru_cache(maxsize=8192)
def _luminance_of(color: str) -> float:
    c = QColor(color)
    def lin(v):
        v = v / 255.0
        return v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4
    return 0.2126 * lin(c.red()) + 0.7152 * lin(c.green()) + 0.0722 * lin(c.blue())


def relative_luminance(color: str) -> float:
    """WCAG relative luminance in [0, 1]. Decides additive versus multiply blending:
    additive glows on a dark page and clips to white on a light one."""
    return _luminance_of(str(color))


def spaceout_enabled() -> bool:
    return False


def advance_spaceout_drift(dt: float) -> float:
    return 0.0


FIELD_FADE_EXPONENT = 3.0
FIELD_FADE_STOPS = 17


def field_fade_alpha(t: float) -> float:
    """Alpha of a field container at fraction t across its width: 1 at the start,
    0 at the end. spaCR's curve, so fields dissolve into the animation the same way."""
    t = max(0.0, min(1.0, float(t)))
    return 1.0 - t ** FIELD_FADE_EXPONENT


def field_chrome(theme: str = "dark") -> dict:
    base = palette_for(theme)
    return {
        "radius": float(_TH.RADIUS["sm"]),
        "fill": (base["surface_alt"], 1.0),
        "fill_disabled": (base["surface"], 1.0),
        "border": (base["border"], 1.0),
        "border_focus": (base["accent"], 1.0),
        "border_disabled": (base["border_soft"], 1.0),
    }


_REGISTERED_QSS: dict = {}


def register_widget_qss(name: str, fn, *, replace: bool = False):
    """Recorded, not applied: ParaScape styles the home screen from `tiles.tile_qss`.
    Kept so engine code that registers a block still runs."""
    if replace or name not in _REGISTERED_QSS:
        _REGISTERED_QSS[name] = fn
    return fn


# --------------------------------------------------------- event watchers etc.
def _watch_object_events(watched, watcher, kinds=(), replay=True) -> bool:
    """spaCR's hub-routed filter, as the plain Qt call it optimises."""
    if watched is None or watcher is None:
        return False
    watched.installEventFilter(watcher)
    return True


def _stop_watching_object_events(watched, watcher) -> bool:
    if watched is None or watcher is None:
        return False
    try:
        watched.removeEventFilter(watcher)
    except RuntimeError:            # either side already destroyed
        return False
    return True


def _watch_application_events(app, watcher, kinds=()) -> bool:
    return _watch_object_events(app, watcher, kinds)


def _stop_watching_application_events(app, watcher) -> bool:
    return _stop_watching_object_events(app, watcher)


def active() -> bool:
    """Whether a heavy job holds the interpreter. ParaScape does not track that."""
    return False


def screen_for_widget(widget):
    """The screen a widget is on, falling back to the primary screen."""
    from PySide6.QtGui import QGuiApplication
    try:
        screen = widget.screen() if widget is not None else None
    except RuntimeError:
        screen = None
    return screen or QGuiApplication.primaryScreen()


def _after(msec: int, owner: QObject, callback):
    """Call `callback` once after `msec`, unless `owner` is deleted first.

    spaCR's helper, re-derived: a single-shot timer parented to the owner, so
    deleting the owner deletes the timer and the call never runs against a dead
    widget. The callback is dropped after it returns, not inside the emission --
    a callback holding the last reference to its owner would otherwise delete
    the emitting timer mid-signal.
    """
    import shiboken6
    if owner is None or not shiboken6.isValid(owner):
        return None
    timer = QTimer(owner)
    timer.setSingleShot(True)
    holder = [callback]

    def fire():
        fn = holder[0]
        holder[0] = None
        if fn is not None:
            fn()

    timer.timeout.connect(fire)
    timer.start(int(msec))
    return timer


# ------------------------------------------------------ timing, misc, stubs
class _Timing:
    """spaCR's frame profiler, as a no-op. `span` wraps every ambient paint, so it
    must exist and must cost nothing."""

    @staticmethod
    @contextlib.contextmanager
    def span(*_a, **_k):
        yield


timing = _timing = _Timing()


def _metal_gpu_name():
    return None


def _numba_cache_folder():
    """Where precompiled animation kernels are cached. A `Path`, as spaCR's is --
    the engine joins onto it with `/`. ParaScape's own folder, so the cache is
    neither shared with spaCR nor scattered across the temp directory root."""
    from pathlib import Path
    from .. import paths
    try:
        root = Path(paths.state_dir()) if hasattr(paths, "state_dir") else None
    except Exception:
        root = None
    if root is None:
        import tempfile
        root = Path(tempfile.gettempdir()) / "parascape"
    folder = root / "ambient-kernels"
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def registry():
    """spaCR's job registry, empty: ParaScape's jobs live in its own runner."""
    return SimpleNamespace(running=lambda: (), jobs={}, active=lambda: ())


class _Unavailable:
    """Stands in for a spaCR subsystem ParaScape deliberately does not ship. Any use
    raises, naming what was reached, because a silent blank backdrop is worse."""

    def __init__(self, what: str):
        self._what = what

    def __getattr__(self, name):
        raise RuntimeError(f"{self._what}.{name}: not part of ParaScape "
                           f"(spaceout and fractal backdrops are spaCR-only)")

    def __call__(self, *a, **k):
        raise RuntimeError(f"{self._what}: not part of ParaScape")


resonance = rs = _Unavailable("resonance")
BANDS = ()


def silence(*a, **k):
    return None


def playing_moment(*a, **k):
    return None


CascadeEngine = SpaceEngine = OrbitEngine = _Unavailable("fractal engine")
RuntimeControls = Settings = _Unavailable("fractal travel")
create_fractal_widget = render_space_frame = _Unavailable("fractal")
_add_galaxies = _galaxies_at = _galaxy_atlas = _Unavailable("fractal space")


class _HeavyImportInProgress(Exception):
    pass


_heavy_import_lock = contextlib.nullcontext()
