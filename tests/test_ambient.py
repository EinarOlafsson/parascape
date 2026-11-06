"""spaCR's animated backdrops, vendored, behind ParaScape's home screen.

The engines are spaCR's code transferred nearly verbatim; what these tests guard
is the seam -- the shim standing in for the rest of spaCR -- because that is
where the port can be wrong without any engine code changing. Two of the three
bugs found while porting were in exactly that seam: a `str` returned where the
engine joins a `Path` onto it, and the module naming itself by spaCR's file name.
"""
import pytest
from PySide6 import QtCore, QtGui


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """ParaScape's ambient preferences in a throwaway INI file, never the user's."""
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "ambient.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))
    return AC


def test_the_menu_offers_spacrs_five_themes_with_the_field_as_default():
    from parascape.backdrops import engines as E
    assert E.AMBIENT_THEMES == ("data_art_impulse_lens", "data_art_genetic_advection",
                                "data_art_point_atlas", "blobs", "drift")
    assert E.DEFAULT_THEME == "data_art_impulse_lens", "the spaCR field is the default"
    assert E.NO_ANIMATION in E.ANIMATION_CHOICES


@pytest.mark.parametrize("theme", ["data_art_impulse_lens", "data_art_genetic_advection",
                                   "data_art_point_atlas", "blobs", "drift"])
def test_every_menu_theme_constructs_and_paints(qapp, theme):
    from parascape import theme as TH
    from parascape.backdrops import engines as E
    bg = QtGui.QColor(TH.palette_for("dark")["page"])
    engine = E.make_engine(theme, E.default_palette_for(theme), bg, seed=3)
    for _ in range(4):
        engine.advance(1 / 30)
    image = QtGui.QImage(320, 200, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(bg)
    painter = QtGui.QPainter(image)
    engine.paint(painter, 320, 200)        # (painter, width, height) -- not a QRect
    painter.end()


def test_the_data_art_themes_actually_draw_something(qapp):
    """A theme that constructs and paints a blank frame would pass the test above."""
    from parascape import theme as TH
    from parascape.backdrops import engines as E
    bg = QtGui.QColor(TH.palette_for("dark")["page"])
    for theme in ("data_art_impulse_lens", "data_art_genetic_advection", "data_art_point_atlas"):
        engine = E.make_engine(theme, E.default_palette_for(theme), bg, seed=3)
        engine.advance(0.2)
        image = QtGui.QImage(320, 200, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        image.fill(bg)
        painter = QtGui.QPainter(image)
        engine.paint(painter, 320, 200)
        painter.end()
        drawn = sum(image.pixel(x, y) != bg.rgba() for x in range(0, 320, 2) for y in range(0, 200, 2))
        assert drawn > 0, f"{theme} painted nothing"


# ------------------------------------------------------------------- the shim


def test_preferences_default_to_spacrs_shipped_animation(settings):
    from parascape.backdrops import engines as E
    assert settings.get_ambient_theme() == E.DEFAULT_THEME
    assert settings.get_ambient_enabled() is True
    assert settings.get_ambient_speed() == E.DEFAULT_SPEED


def test_an_unknown_saved_theme_falls_back_rather_than_failing(settings):
    from parascape.backdrops import engines as E
    settings._settings().setValue("ambient/theme", "a_theme_from_a_newer_build")
    assert settings.get_ambient_theme() == E.DEFAULT_THEME


def test_choosing_no_animation_turns_the_backdrop_off(settings):
    from parascape.backdrops import engines as E
    settings._settings().setValue("ambient/theme", E.NO_ANIMATION)
    assert settings.get_ambient_enabled() is False


def test_the_kernel_cache_folder_is_a_path():
    """The engine does `folder / name`. A `str` here broke the background kernel
    compiler with a TypeError in a thread nobody was watching."""
    from pathlib import Path
    from parascape.backdrops import _compat as AC
    assert isinstance(AC._numba_cache_folder(), Path)


def test_the_engine_stamps_its_own_file_not_spacrs():
    """The precompiled-kernel cache fingerprints its own source by file name. Left as
    spaCR's `ambient.py`, the stat failed in the compile thread."""
    import os
    from parascape.backdrops import engines as E
    here = os.path.dirname(E.__file__)
    for name in E._PRIMED_SOURCES:
        assert os.path.isfile(os.path.join(here, name + ".py")), name
    for kernel in E._PRECOMPILED_KERNELS:
        assert kernel.startswith("parascape.backdrops.engines."), kernel


def test_a_spaceout_only_subsystem_fails_loudly(settings):
    """A blank backdrop is worse than an error naming what was reached."""
    with pytest.raises(RuntimeError, match="not part of ParaScape"):
        settings.get_fractal_settings()
    with pytest.raises(RuntimeError, match="not part of ParaScape"):
        settings.resonance.mode_blend


def test_after_never_fires_against_a_deleted_owner(qapp):
    import shiboken6
    from PySide6.QtCore import QObject
    from parascape.backdrops import _compat as AC
    owner = QObject()
    fired = []
    AC._after(5, owner, lambda: fired.append(1))
    shiboken6.delete(owner)
    for _ in range(20):
        qapp.processEvents()
        QtCore.QThread.msleep(2)
    assert fired == []


# --------------------------------------------------------- behind the home page


def test_the_home_screen_animates_with_transparent_containers(qtbot, settings):
    from parascape.home import HomeScreen
    home = HomeScreen()
    qtbot.addWidget(home)
    assert home.animated
    for name in ("HomeScreen", "HomeScroll", "HomeScrollViewport", "HomeBands"):
        widget = home if name == "HomeScreen" else home.findChild(QtCore.QObject, name)
        assert widget is not None and widget.property("ambient") == "on", \
            f"{name} stays opaque, so the backdrop runs invisibly underneath it"


def test_with_animation_off_the_page_keeps_its_own_colour(qtbot, settings):
    """The other half: transparent containers with NO backdrop paint the Base colour,
    which is near-white -- the bug the very first render of this page caught."""
    from parascape.backdrops import engines as E
    from parascape.home import HomeScreen
    settings._settings().setValue("ambient/theme", E.NO_ANIMATION)
    home = HomeScreen()
    qtbot.addWidget(home)
    assert not home.animated
    assert home.property("ambient") != "on"


def test_a_backdrop_failure_never_stops_the_page_opening(qtbot, settings, monkeypatch):
    from parascape.backdrops import engines as E
    from parascape.home import HomeScreen

    def boom(*a, **k):
        raise RuntimeError("no GL today")

    monkeypatch.setattr(E, "install_ambient", boom)
    home = HomeScreen()
    qtbot.addWidget(home)
    assert not home.animated and "no GL today" in home._ambient_error
    assert len(home._tiles) > 0
