"""Recipes of module settings (instruction 98, work package 6).

A recipe is a named bundle of the settings of one module, kept as one file, that knows its module
and the version of ParaScape it was saved with. These tests mirror what spaCR guarantees for its
own recipes: a round trip, a refusal across modules, a version gap that names the settings that
have gone, a damaged file that is reported rather than ignored, an overridable folder, and a name
with awkward characters that gets a safe file without taking another recipe file.

They are not about :mod:`parascape.recipes`, which holds the Infer recipes of a biological question.
"""
import json
import os

import pytest

from parascape import module_recipes as R


@pytest.fixture
def folder(tmp_path, monkeypatch):
    """Keep every recipe in a scratch folder so no test touches the home folder of the user."""
    root = tmp_path / "recipes"
    monkeypatch.setenv(R.RECIPE_DIR_ENV, str(root))
    return root


def make(name="Tight map", module="map", **settings):
    """A recipe with a few settings, captured the way a module captures its own."""
    return R.capture_recipe(module, name, settings or {"n_neighbors": 15, "min_dist": 0.05})


# ------------------------------------------------------------------ round trip and format
def test_a_recipe_survives_a_save_and_a_load(folder):
    recipe = R.capture_recipe("map", "Tight map", {"n_neighbors": 15, "min_dist": 0.05,
                                                  "measurements": ["a", "b"], "live": True},
                              notes="for the screen")
    path = R.save_recipe(recipe)
    assert path == str(folder / "map" / "tight-map.json")
    loaded = R.load_recipe(path)
    assert loaded.name == "Tight map" and loaded.module == "map"
    assert loaded.settings == recipe.settings
    assert loaded.notes == "for the screen"
    assert loaded.parascape_version == R.parascape_version()
    assert loaded.created and loaded.path == path
    assert [r.name for r in R.list_recipes("map")] == ["Tight map"]


def test_the_file_is_one_json_with_the_documented_keys(folder):
    path = R.save_recipe(make())
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    assert set(data) == {"parascape_recipe", "name", "module", "parascape_version", "created",
                         "notes", "settings"}
    assert data["parascape_recipe"] == R.FORMAT_VERSION == 1
    assert data["module"] == "map" and data["parascape_version"] == R.parascape_version()
    assert data["settings"] == {"n_neighbors": 15, "min_dist": 0.05}
    assert not any(f.startswith(".recipe-") for f in os.listdir(folder / "map")), "temp file left"


def test_the_format_version_is_not_the_application_version(folder):
    """A file written by any application version is read when its format is one this reads."""
    path = R.save_recipe(make())
    data = json.load(open(path))
    data["parascape_version"] = "9.9.9"
    json.dump(data, open(path, "w"))
    assert R.load_recipe(path).parascape_version == "9.9.9"


def test_numpy_values_and_sets_are_kept_as_plain_data(folder):
    np = pytest.importorskip("numpy")
    recipe = R.capture_recipe("map", "np", {"a": np.float64(0.5), "b": np.int64(3), "c": {"x", "y"}})
    loaded = R.load_recipe(R.save_recipe(recipe))
    assert loaded.settings == {"a": 0.5, "b": 3, "c": ["x", "y"]}


# ------------------------------------------------------------------ the refusals
def test_a_recipe_from_another_module_is_refused_and_nothing_is_written(folder):
    recipe = make(module="map")
    written = []
    with pytest.raises(R.WrongModuleError) as error:
        R.apply_recipe(recipe, "network", {"n_neighbors", "min_dist"},
                       lambda values: written.append(values) or len(values))
    assert "'map'" in str(error.value) and "'network'" in str(error.value)
    assert written == [], "a refused recipe must not be partly applied"
    with pytest.raises(R.WrongModuleError):
        R.plan_apply(recipe, "network", {"n_neighbors"})


def test_importing_a_recipe_made_in_another_module_is_refused(folder, tmp_path):
    shared = tmp_path / "from_a_colleague.json"
    R.export_recipe(make(module="map"), str(shared))
    with pytest.raises(R.WrongModuleError):
        R.import_recipe(str(shared), "network")
    assert R.list_recipes("network") == []


def test_a_version_gap_names_the_settings_the_module_no_longer_has(folder):
    recipe = R.capture_recipe(
        "map", "old", {"n_neighbors": 15, "min_dist": 0.05, "retired_b": 1, "retired_a": 2})
    recipe.parascape_version = "0.0.9.0"
    plan = R.plan_apply(recipe, "map", {"n_neighbors", "min_dist", "seed"}, running="0.0.11.0")
    assert plan.version_differs and plan.needs_confirmation
    assert plan.missing == ("retired_a", "retired_b")
    assert plan.settings == {"n_neighbors": 15, "min_dist": 0.05}
    note = plan.note()
    assert "0.0.9.0" in note and "0.0.11.0" in note
    assert "retired_a" in note and "retired_b" in note, "the settings are named, not counted"
    seen = []
    report = R.apply_recipe(recipe, "map", {"n_neighbors", "min_dist", "seed"},
                            lambda values: seen.append(values) or len(values), running="0.0.11.0")
    assert seen == [{"n_neighbors": 15, "min_dist": 0.05}], "the missing settings are not written"
    assert report.applied == 2


def test_a_missing_setting_is_reported_even_when_the_version_is_the_same(folder):
    recipe = R.capture_recipe("map", "r", {"n_neighbors": 15, "gone": 1})
    recipe.parascape_version = "1.0"
    plan = R.plan_apply(recipe, "map", {"n_neighbors"}, running="1.0")
    assert not plan.version_differs and plan.missing == ("gone",) and plan.needs_confirmation


def test_a_matching_recipe_needs_no_question(folder):
    recipe = R.capture_recipe("map", "r", {"n_neighbors": 15})
    recipe.parascape_version = "1.0"
    plan = R.plan_apply(recipe, "map", {"n_neighbors", "seed"}, running="1.0")
    assert not plan.needs_confirmation and plan.note() == ""


def test_a_damaged_file_is_reported_not_ignored(folder):
    R.save_recipe(make("Good"))
    mod = folder / "map"
    (mod / "broken.json").write_text("{not json at all", encoding="utf-8")
    (mod / "package.json").write_text(json.dumps({"name": "x", "version": "1"}), encoding="utf-8")
    (mod / "newer.json").write_text(json.dumps(
        {"parascape_recipe": 2, "module": "map", "settings": {}}), encoding="utf-8")
    (mod / "no_module.json").write_text(json.dumps(
        {"parascape_recipe": 1, "settings": {}}), encoding="utf-8")
    (mod / "no_settings.json").write_text(json.dumps(
        {"parascape_recipe": 1, "module": "map"}), encoding="utf-8")
    (mod / "binary.json").write_bytes(b"\xff\xfe\x00\x80")
    other = make("Misfiled", module="network")
    other.path = ""
    R.save_recipe(other, directory=str(mod))
    listing = R.scan_recipes("map")
    assert [r.name for r in listing.recipes] == ["Good"], "good recipes are still listed"
    reasons = {os.path.basename(p): why for p, why in listing.problems}
    assert set(reasons) == {"broken.json", "package.json", "newer.json", "no_module.json",
                            "no_settings.json", "binary.json", "misfiled.json"}
    assert "not valid JSON" in reasons["broken.json"]
    assert "not a ParaScape recipe" in reasons["package.json"]
    assert "newer" in reasons["newer.json"] and "2" in reasons["newer.json"]
    assert "module" in reasons["no_module.json"]
    assert "settings" in reasons["no_settings.json"]
    assert "'network'" in reasons["misfiled.json"]
    for why in reasons.values():
        assert why, "every problem says why"
    with pytest.raises(R.NotARecipeError) as error:
        R.load_recipe(str(mod / "broken.json"))
    assert "broken.json" in str(error.value)
    with pytest.raises(R.NotARecipeError):
        R.load_recipe(str(mod / "missing.json"))


def test_saving_never_overwrites_a_damaged_file(folder):
    mod = folder / "map"
    mod.mkdir(parents=True)
    (mod / "tight-map.json").write_text("garbage", encoding="utf-8")
    path = R.save_recipe(make("Tight map"))
    assert os.path.basename(path) == "tight-map-2.json"
    assert (mod / "tight-map.json").read_text(encoding="utf-8") == "garbage"


# ------------------------------------------------------------------ storage
def test_the_folder_override_is_honoured(tmp_path, monkeypatch):
    elsewhere = tmp_path / "shared" / "lab"
    monkeypatch.setenv(R.RECIPE_DIR_ENV, str(elsewhere))
    path = R.save_recipe(make())
    assert path.startswith(str(elsewhere / "map"))
    assert R.recipes_dir() == str(elsewhere)
    assert R.recipes_dir("map") == str(elsewhere / "map")


def test_without_an_override_the_folder_is_under_the_home_folder(tmp_path, monkeypatch):
    monkeypatch.delenv(R.RECIPE_DIR_ENV, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    assert R.recipes_dir("map", create=False) == str(tmp_path / ".parascape" / "recipes" / "map")


def test_looking_creates_no_folder(folder):
    assert R.scan_recipes("map").recipes == []
    assert R.find_recipe("map", "nothing") is None
    assert not folder.exists()


def test_each_module_has_its_own_folder_and_list(folder):
    R.save_recipe(make("A", module="map"))
    R.save_recipe(make("B", module="network"))
    assert [r.name for r in R.list_recipes("map")] == ["A"]
    assert [r.name for r in R.list_recipes("network")] == ["B"]
    assert (folder / "map").is_dir() and (folder / "network").is_dir()


# ------------------------------------------------------------------ names and file names
def test_awkward_names_give_safe_file_names_that_stay_inside_the_folder(folder):
    names = ["../../etc/passwd", "a/b\\c:d*e?f", "Café ☕ crème", "日本語", "???", "CON",
             "x" * 500, ".hidden", "tab\tname\nnewline"]
    for name in names:
        path = R.save_recipe(R.capture_recipe("map", name, {"k": 1}))
        assert os.path.dirname(path) == str(folder / "map"), name
        stem = os.path.basename(path)[:-5]
        assert stem and all(c.isalnum() or c == "-" for c in stem), (name, stem)
        assert len(stem) <= 70
        assert R.load_recipe(path).name == name.strip(), "the display name is kept verbatim"
    assert R.slug("Café ☕ crème") == "cafe-creme"
    assert R.slug("CON") == "con-recipe"
    assert R.slug("日本語").startswith("recipe-") and R.slug("日本語") != R.slug("???")
    assert not (folder / "etc").exists() and not (folder.parent / "etc").exists()


def test_names_that_slug_alike_do_not_collide(folder):
    first = R.save_recipe(R.capture_recipe("map", "Tight map", {"k": 1}))
    second = R.save_recipe(R.capture_recipe("map", "tight-map", {"k": 2}))
    third = R.save_recipe(R.capture_recipe("map", "Tight  Map!", {"k": 3}))
    assert len({first, second, third}) == 3
    by_name = {r.name: r.settings["k"] for r in R.list_recipes("map")}
    assert by_name == {"Tight map": 1, "tight-map": 2, "Tight  Map!": 3}, "nothing was overwritten"


def test_two_names_with_no_plain_letters_do_not_collide(folder):
    a = R.save_recipe(R.capture_recipe("map", "日本語", {"k": 1}))
    b = R.save_recipe(R.capture_recipe("map", "中文", {"k": 2}))
    assert a != b and len(R.list_recipes("map")) == 2


def test_saving_the_same_name_again_replaces_it(folder):
    R.save_recipe(R.capture_recipe("map", "Tight map", {"k": 1}))
    R.save_recipe(R.capture_recipe("map", "tight MAP ", {"k": 2}))
    listing = R.list_recipes("map")
    assert len(listing) == 1 and listing[0].settings == {"k": 2}


def test_a_recipe_needs_a_name_and_a_module(folder):
    with pytest.raises(R.RecipeError):
        R.capture_recipe("map", "   ", {})
    with pytest.raises(R.RecipeError):
        R.capture_recipe("", "x", {})


# ------------------------------------------------------------------ rename, delete, share
def test_rename_moves_the_file_and_refuses_a_taken_name(folder):
    a = make("A")
    R.save_recipe(a)
    b = make("B")
    R.save_recipe(b)
    with pytest.raises(R.RecipeExistsError):
        R.rename_recipe(b, "a")
    assert b.name == "B" and os.path.isfile(b.path), "a refused rename changes nothing"
    old = a.path
    new = R.rename_recipe(a, "Wider map")
    assert os.path.basename(new) == "wider-map.json" and not os.path.exists(old)
    assert sorted(r.name for r in R.list_recipes("map")) == ["B", "Wider map"]
    with pytest.raises(R.RecipeError):
        R.rename_recipe(a, "  ")


def test_rename_changing_only_the_case_keeps_one_file(folder):
    recipe = make("tight map")
    R.save_recipe(recipe)
    R.rename_recipe(recipe, "Tight Map")
    assert [r.name for r in R.list_recipes("map")] == ["Tight Map"]
    assert len(os.listdir(folder / "map")) == 1


def test_delete_removes_the_file(folder):
    recipe = make()
    R.save_recipe(recipe)
    assert R.delete_recipe(recipe) is True
    assert R.list_recipes("map") == [] and R.delete_recipe(recipe) is False


def test_sharing_is_one_file_that_another_person_can_import(folder, tmp_path, monkeypatch):
    recipe = R.capture_recipe("map", "Tight map", {"n_neighbors": 15}, notes="n")
    sent = tmp_path / "outbox" / "tight.json"
    R.export_recipe(recipe, str(sent))
    monkeypatch.setenv(R.RECIPE_DIR_ENV, str(tmp_path / "their_folder"))
    got = R.import_recipe(str(sent), "map")
    assert got.name == "Tight map" and got.settings == {"n_neighbors": 15} and got.notes == "n"
    assert [r.name for r in R.list_recipes("map")] == ["Tight map"]
    with pytest.raises(R.RecipeExistsError):
        R.import_recipe(str(sent), "map")
    again = R.import_recipe(str(sent), "map", replace=True)
    assert len(R.list_recipes("map")) == 1 and again.path


def test_importing_a_damaged_file_says_why(folder, tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("[1, 2", encoding="utf-8")
    with pytest.raises(R.NotARecipeError) as error:
        R.import_recipe(str(bad), "map")
    assert "bad.json" in str(error.value) and "JSON" in str(error.value)


# ------------------------------------------------------------------ the panel
from PySide6 import QtCore, QtWidgets  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def screen(qapp, folder):
    """A small module with one of each kind of control, and a panel on it."""
    from parascape.recipes_panel import RecipePanel, widget_setting
    holder = QtWidgets.QWidget()
    spin = QtWidgets.QSpinBox()
    spin.setRange(2, 100)
    spin.setValue(25)
    dspin = QtWidgets.QDoubleSpinBox()
    dspin.setRange(0, 1)
    dspin.setValue(0.1)
    combo = QtWidgets.QComboBox()
    combo.addItems(["hdbscan", "dbscan"])
    tick = QtWidgets.QCheckBox()
    text = QtWidgets.QLineEdit("5, 15")
    for w in (spin, dspin, combo, tick, text):
        w.setParent(holder)
    applied = []
    panel = RecipePanel("map", [widget_setting("n", spin), widget_setting("d", dspin),
                                widget_setting("method", combo), widget_setting("tick", tick),
                                widget_setting("grid", text)],
                        after_apply=lambda: applied.append(1))
    panel.answers = []
    panel.warnings = []
    panel.ask_yes_no = lambda title, text: (panel.answers.append((title, text)) or True)
    panel.warn = lambda title, text: panel.warnings.append((title, text))
    panel.holder = holder
    panel.controls = dict(spin=spin, dspin=dspin, combo=combo, tick=tick, text=text)
    panel.applied_calls = applied
    yield panel
    holder.deleteLater()
    panel.deleteLater()


def test_the_panel_saves_and_applies_a_recipe_round_trip(screen):
    c = screen.controls
    c["spin"].setValue(15)
    c["dspin"].setValue(0.05)
    c["combo"].setCurrentText("dbscan")
    c["tick"].setChecked(True)
    c["text"].setText("1, 2")
    saved = screen.save_as("Tight")
    assert saved.settings == {"n": 15, "d": 0.05, "method": "dbscan", "tick": True, "grid": "1, 2"}
    assert [r.name for r in screen.recipes()] == ["Tight"] and screen.selected().name == "Tight"
    c["spin"].setValue(40)
    c["combo"].setCurrentText("hdbscan")
    c["tick"].setChecked(False)
    c["text"].setText("x")
    report = screen.apply_selected()
    assert report.applied == 5 and screen.answers == [], "same version, nothing missing: no question"
    assert (c["spin"].value(), c["dspin"].value(), c["combo"].currentText(),
            c["tick"].isChecked(), c["text"].text()) == (15, 0.05, "dbscan", True, "1, 2")
    assert screen.applied_calls == [1]
    assert "Applied" in screen.detail_text() and "5 settings" in screen.detail_text()


def test_the_panel_asks_before_applying_an_old_recipe_and_names_what_is_gone(screen, folder):
    old = R.capture_recipe("map", "Old", {"n": 10, "retired_one": 1, "retired_two": 2})
    old.parascape_version = "0.0.1"
    R.save_recipe(old)
    screen.reload()
    assert "retired_one" in screen.detail_text() and "retired_two" in screen.detail_text()
    assert "0.0.1" in screen.detail_text(), "told before applying, not after"
    screen.apply_selected()
    assert len(screen.answers) == 1
    assert "retired_one" in screen.answers[0][1] and "0.0.1" in screen.answers[0][1]
    assert screen.controls["spin"].value() == 10
    screen.answers.clear()
    screen.controls["spin"].setValue(50)
    screen.ask_yes_no = lambda t, x: False
    assert screen.apply_selected() is None and screen.controls["spin"].value() == 50, "declined"
    assert screen.apply_selected(confirmed=True).applied == 1


def test_the_panel_refuses_a_recipe_of_another_module(screen, folder):
    wrong = R.capture_recipe("network", "Net", {"n": 5})
    path = R.save_recipe(wrong)
    # put it where the map recipes are, as a person copying files by hand would
    target = folder / "map" / "net.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(path, target)
    screen.reload()
    assert screen.recipes() == [], "listed nowhere as usable"
    assert "'network'" in screen.problems_text()
    assert not screen.problems.isHidden()
    # and applying one that was selected from elsewhere is refused with nothing written
    screen._listed = [wrong]
    screen.choice.clear()
    screen.choice.addItem("Net")
    screen.controls["spin"].setValue(33)
    assert "network" in screen.detail_text() and not screen.buttons["Apply"].isEnabled()
    assert screen.apply_selected() is None
    assert screen.controls["spin"].value() == 33
    assert screen.warnings and "network" in screen.warnings[-1][1]


def test_a_choice_that_is_not_offered_any_more_is_reported_and_the_rest_apply(screen, folder):
    R.save_recipe(R.capture_recipe("map", "Gone", {"n": 11, "method": "optics", "d": 7.0}))
    screen.reload()
    report = screen.apply_selected()
    assert report.applied == 1 and screen.controls["spin"].value() == 11
    text = screen.detail_text()
    assert "method" in text and "optics" in text and "d" in text
    assert screen.controls["combo"].currentText() == "hdbscan"


def test_the_panel_asks_before_replacing_and_before_deleting(screen):
    screen.save_as("Same")
    screen.controls["spin"].setValue(77)
    screen.ask_yes_no = lambda t, x: False
    assert screen.save_as("same") is None
    assert screen.recipes()[0].settings["n"] == 25, "declined, so not replaced"
    screen.ask_yes_no = lambda t, x: True
    assert screen.save_as("same").settings["n"] == 77 and len(screen.recipes()) == 1
    screen.ask_yes_no = lambda t, x: False
    assert screen.delete() is False and len(screen.recipes()) == 1
    screen.ask_yes_no = lambda t, x: True
    assert screen.delete() is True and screen.recipes() == []
    assert not screen.buttons["Apply"].isEnabled(), "nothing to apply"
    assert "No recipes yet" in screen.detail_text()


def test_the_panel_renames_shares_and_imports(screen, tmp_path):
    screen.save_as("One")
    assert screen.rename("Two") and [r.name for r in screen.recipes()] == ["Two"]
    out = tmp_path / "share.json"
    screen.choose_save = lambda suggested: str(out)
    screen.share()
    assert json.load(open(out))["name"] == "Two"
    screen.delete()
    screen.choose_open = lambda: str(out)
    assert screen.import_file().name == "Two"
    bad = tmp_path / "bad.json"
    bad.write_text("nope", encoding="utf-8")
    screen.choose_open = lambda: str(bad)
    assert screen.import_file() is None and "bad.json" in screen.warnings[-1][1]


def test_every_control_of_the_panel_has_a_tooltip(screen):
    controls = [screen.choice, screen.problems] + list(screen.buttons.values())
    assert [w for w in controls if not w.toolTip()] == []
    assert set(screen.buttons) == {"Apply", "Save current settings…", "Share…", "Import…",
                                   "Rename…", "Delete"}


def test_a_control_a_recipe_cannot_keep_is_a_mistake_at_build_time(qapp):
    from parascape.recipes_panel import widget_setting
    with pytest.raises(TypeError):
        widget_setting("x", QtWidgets.QLabel("no"))


def test_a_number_outside_the_range_of_the_control_is_refused(qapp):
    from parascape.recipes_panel import widget_setting
    spin = QtWidgets.QSpinBox()
    spin.setRange(2, 10)
    setting = widget_setting("n", spin)
    with pytest.raises(ValueError):
        setting.set(500)
    with pytest.raises(ValueError):
        setting.set("a lot")
    with pytest.raises(ValueError):
        setting.set(True)


def test_the_shell_gives_a_module_a_recipes_group_that_starts_folded(qapp, folder, tmp_path,
                                                                    monkeypatch):
    from parascape.backdrops import _compat as AC
    ini = str(tmp_path / "shell.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(ini, QtCore.QSettings.Format.IniFormat))
    from parascape.recipes_panel import widget_setting
    from parascape.shell import ModuleShell
    shell = ModuleShell("t_recipes", "Test")
    box = QtWidgets.QSpinBox()
    shell.add_setting("A group", box)
    panel = shell.add_recipes([widget_setting("n", box)])
    assert shell.recipes is panel and "Recipes" in shell.sections
    assert not shell.sections["Recipes"].is_expanded()
    names = list(shell.sections)
    assert names.index("Recipes") > names.index("A group")
    shell.close()


# ------------------------------------------------------------------ in the real window
@pytest.fixture(scope="module")
def window(qapp):
    from parascape import app as A
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                      # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    from parascape.map_module import wait_for_projections
    wait_for_projections()
    win.close()
    qapp.processEvents()


def test_the_map_module_keeps_and_restores_its_settings_as_a_recipe(window, folder):
    shell = window.map_module
    assert shell is not None and shell.recipes is not None
    assert not shell.sections["Recipes"].is_expanded()
    names = list(shell.sections)
    assert names.index("Recipes") > names.index("Advanced")
    panel, module = shell.recipes, shell.module
    module.live.setChecked(False)
    window.panel.nn.setValue(17)
    window.panel.md.setValue(0.3)
    window.panel.seed.setValue(7)
    before = sorted(window.panel.blocks.checked())
    panel.save_as("Map recipe test")
    window.panel.nn.setValue(60)
    window.panel.md.setValue(0.0)
    window.panel.seed.setValue(1)
    window.panel.blocks.set_checked([])
    panel.apply_selected(confirmed=True)
    assert (window.panel.nn.value(), window.panel.md.value(), window.panel.seed.value()) == (17, 0.3, 7)
    assert sorted(window.panel.blocks.checked()) == before, "the exact measurements come back"
    saved = R.find_recipe("map", "Map recipe test")
    assert saved.module == "map" and "measurements" in saved.settings
    module.live.setChecked(True)
    window.panel.nn.setValue(25)
    window.panel.md.setValue(0.1)
    window.panel.seed.setValue(42)


def test_a_map_recipe_is_refused_by_the_network_module(window, folder):
    network = window.network_module
    assert network is not None and network.recipes is not None
    R.save_recipe(R.capture_recipe("map", "From map", {"n_neighbors": 9}))
    recipe = R.find_recipe("map", "From map")
    with pytest.raises(R.WrongModuleError):
        R.apply_recipe(recipe, "network", network.recipes.settings, lambda v: len(v))
    assert R.list_recipes("network") == []


def test_the_network_module_keeps_its_links_and_strength_as_a_recipe(window, folder):
    network = window.network_module
    module, panel = network.module, network.recipes
    assert not network.sections["Recipes"].is_expanded()
    keys = set(panel.settings)
    assert {"colour_by", "link_strength", "arrange_genes", "sharp_link_steps",
            "name_linked_genes"} <= keys
    assert all(f"links_{layer}" in keys for layer in module.layer_boxes)
    layer = next(iter(module.layer_boxes))
    before = {name: box.isChecked() for name, box in module.layer_boxes.items()}
    module.strength.setValue(70)
    panel.save_as("Network recipe test")
    module.strength.setValue(20)
    module.layer_boxes[layer].setChecked(not before[layer])
    panel.apply_selected(confirmed=True)
    assert module.strength.value() == 70
    assert {name: box.isChecked() for name, box in module.layer_boxes.items()} == before
    module.strength.setValue(50)
