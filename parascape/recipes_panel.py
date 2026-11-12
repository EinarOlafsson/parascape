"""The recipe controls of a module: save its settings under a name, apply one, share one.

The logic lives in :mod:`parascape.module_recipes` and has no Qt. This file only shows it. A module
declares which of its controls are settings (:class:`Setting`, or :func:`widget_setting` for an
ordinary spin box, drop-down, tick box or text field) and hands them to
``ModuleShell.add_recipes``, which puts this panel in the Settings column as a section that
starts folded, like the Save section of every module.

Every question that would open a dialog goes through a small set of replaceable functions
(:attr:`RecipePanel.ask_name` and the others), so the tests and any other host can answer them
without a modal window.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Optional

from PySide6 import QtCore, QtWidgets

from . import module_recipes as R
from . import theme as TH


@dataclass
class Setting:
    """One setting of a module that a recipe can keep.

    :param key: the name stored in the recipe. It must stay the same between versions for old
        recipes to keep applying, which is why a setting that goes away is reported by name.
    :param get: returns the value now. It must be plain data (number, text, list, true/false).
    :param set: puts a value back. It raises :class:`ValueError` with a plain reason when the
        value cannot be used any more (a choice the module no longer offers). It may return a
        note to show the user when it could only use part of the value.
    :param label: plain words for the setting, shown instead of the key where one is needed.
    """

    key: str
    get: Callable[[], Any]
    set: Callable[[Any], Optional[str]]
    label: str = ""


def widget_setting(key: str, widget: QtWidgets.QWidget, label: str = "") -> Setting:
    """A :class:`Setting` that reads and writes an ordinary Qt control.

    This is what most settings are, so a module lists its controls instead of writing a getter
    and a setter for each. A drop-down refuses a value it does not offer, rather than showing a
    blank, so the user is told that the choice has gone.

    :param key: the name stored in the recipe.
    :param widget: a spin box, double spin box, check box, drop-down, text field or slider.
    :param label: plain words for the setting.
    :raises TypeError: for a control that is not one of those, so the mistake shows at build
        time rather than as a recipe that quietly keeps nothing.
    """
    if isinstance(widget, QtWidgets.QComboBox):
        def get() -> Any:
            return widget.currentText()

        def put(value: Any) -> Optional[str]:
            index = widget.findText(str(value))
            if index < 0:
                raise ValueError(f"{value!r} is not offered any more")
            widget.setCurrentIndex(index)
            return None
    elif isinstance(widget, QtWidgets.QCheckBox):
        def get() -> Any:
            return bool(widget.isChecked())

        def put(value: Any) -> Optional[str]:
            widget.setChecked(bool(value))
            return None
    elif isinstance(widget, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox, QtWidgets.QSlider)):
        def get() -> Any:
            return widget.value()

        def put(value: Any) -> Optional[str]:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{value!r} is not a number")
            if not (widget.minimum() <= value <= widget.maximum()):
                raise ValueError(
                    f"{value!r} is outside {widget.minimum()} to {widget.maximum()}")
            widget.setValue(type(widget.value())(value))
            return None
    elif isinstance(widget, QtWidgets.QLineEdit):
        def get() -> Any:
            return widget.text()

        def put(value: Any) -> Optional[str]:
            widget.setText(str(value))
            return None
    else:
        raise TypeError(f"{type(widget).__name__} is not a control that a recipe can keep")
    return Setting(key=key, get=get, set=put, label=label)


class RecipePanel(QtWidgets.QWidget):
    """Save the settings of a module as a named recipe, and apply, share or delete one.

    :param module: key of the module. It scopes which recipes are listed and is what a recipe
        is checked against when applied.
    :param settings: the settings of the module that a recipe keeps.
    :param after_apply: called after a recipe has been applied, so the module can redraw.
    :param parent: parent widget.
    """

    applied = QtCore.Signal(str)

    def __init__(self, module: str, settings: Iterable[Setting],
                 after_apply: Optional[Callable[[], None]] = None, parent=None):
        """Build the panel and read the recipe folder.

        :param module: key of the module.
        :param settings: the settings a recipe of this module keeps; the keys are also the list
            of what the module has, against which an old recipe is checked.
        :param after_apply: called once after a successful apply.
        :param parent: parent widget.
        """
        super().__init__(parent)
        self.module = module
        self.settings: Dict[str, Setting] = {s.key: s for s in settings}
        self.after_apply = after_apply
        self._listed: List[R.Recipe] = []
        self._skipped: Dict[str, str] = {}
        self._notes: List[str] = []
        self.setObjectName("RecipePanel")

        # The questions that would open a dialog. A host or a test replaces any of them.
        self.ask_name: Callable[[str], Optional[str]] = self._ask_name
        self.ask_yes_no: Callable[[str, str], bool] = self._ask_yes_no
        self.choose_open: Callable[[], str] = self._choose_open
        self.choose_save: Callable[[str], str] = self._choose_save
        self.warn: Callable[[str, str], None] = self._warn

        column = QtWidgets.QVBoxLayout(self)
        column.setContentsMargins(6, 2, 6, 6)
        column.setSpacing(6)

        intro = QtWidgets.QLabel(
            "A recipe is the settings of this module under a name you chose. Save one when the "
            "module is set up the way you want, then apply it next time instead of setting it "
            "up again.")
        intro.setObjectName("ShellMuted")
        intro.setWordWrap(True)
        column.addWidget(intro)

        self.choice = QtWidgets.QComboBox()
        self.choice.setObjectName("RecipeChoice")
        self.choice.setToolTip(TH.tip(
            "Your saved recipes for this module, newest first, each with the version of "
            "ParaScape it was saved with."))
        from .shell import solid_dropdown
        solid_dropdown(self.choice)
        self.choice.currentIndexChanged.connect(lambda *_: self._describe())
        column.addWidget(self.choice)

        self.detail = QtWidgets.QLabel("")
        self.detail.setObjectName("ShellMuted")
        self.detail.setWordWrap(True)
        self.detail.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.detail.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        column.addWidget(self.detail)

        self.problems = QtWidgets.QLabel("")
        self.problems.setObjectName("ShellMuted")
        self.problems.setWordWrap(True)
        self.problems.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.problems.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.problems.setToolTip(TH.tip(
            "Files in the recipe folder that could not be read. They are left where they are; "
            "fix or remove them."))
        column.addWidget(self.problems)

        grid = QtWidgets.QGridLayout()
        grid.setSpacing(6)
        self.buttons: Dict[str, QtWidgets.QPushButton] = {}
        specs = [
            ("Apply", self.apply_selected,
             "Put the settings of the chosen recipe into this module. If it was saved with "
             "another version, or names settings this module no longer has, you are told which "
             "before anything changes."),
            ("Save current settings…", self.save_dialog,
             "Keep the settings this module has now under a name you choose. A recipe with the "
             "same name is replaced after you say yes."),
            ("Share…", self.share,
             "Write the chosen recipe to one file you can send to someone."),
            ("Import…", self.import_file,
             "Add a recipe somebody sent you. A recipe made in another module is refused."),
            ("Rename…", self.rename,
             "Give the chosen recipe a different name."),
            ("Delete", self.delete,
             "Remove the chosen recipe. This cannot be undone."),
        ]
        for i, (text, slot, tip) in enumerate(specs):
            button = QtWidgets.QPushButton(text)
            button.setToolTip(TH.tip(tip))
            button.clicked.connect(lambda _=False, f=slot: f())
            if text == "Apply":
                button.setProperty("primary", True)
            grid.addWidget(button, i // 2, i % 2)
            self.buttons[text] = button
        column.addLayout(grid)
        self.reload()

    # ------------------------------------------------------------------ what is listed
    def reload(self, select: str = "") -> None:
        """Read the recipe folder again and show what is in it.

        Files that could not be read are shown under the list instead of being left out.

        :param select: path of a recipe to select afterwards; the newest is selected otherwise.
        """
        listing = R.scan_recipes(self.module)
        self._listed = listing.recipes
        self.choice.blockSignals(True)
        self.choice.clear()
        for recipe in self._listed:
            self.choice.addItem(f"{recipe.name}   ·   {recipe.parascape_version}", recipe.path)
        row = 0
        if select:
            for i, recipe in enumerate(self._listed):
                if os.path.abspath(recipe.path) == os.path.abspath(select):
                    row = i
        if self._listed:
            self.choice.setCurrentIndex(row)
        self.choice.blockSignals(False)
        if listing.problems:
            lines = [f"{len(listing.problems)} file(s) in the recipe folder could not be used:"]
            lines += [f"  {reason}" for _path, reason in listing.problems]
            self.problems.setText("\n".join(lines))
        else:
            self.problems.setText("")
        self.problems.setVisible(bool(listing.problems))
        self._describe()

    def recipes(self) -> List[R.Recipe]:
        """The recipes currently listed, newest first."""
        return list(self._listed)

    def selected(self) -> Optional[R.Recipe]:
        """The chosen recipe, or ``None`` when there is none."""
        row = self.choice.currentIndex()
        return self._listed[row] if 0 <= row < len(self._listed) else None

    def detail_text(self) -> str:
        """The text under the list. Public so that tests read what the user reads."""
        return self.detail.text()

    def problems_text(self) -> str:
        """The note about files that could not be read. Public for the same reason."""
        return self.problems.text()

    def _describe(self) -> None:
        """Say what the chosen recipe holds and what would not carry over, before it is applied."""
        recipe = self.selected()
        for name, button in self.buttons.items():
            if name not in ("Save current settings…", "Import…"):
                button.setEnabled(recipe is not None)
        if recipe is None:
            self.detail.setText(
                "No recipes yet for this module. Set it up the way you want, then choose "
                "“Save current settings…”.")
            return
        try:
            plan = R.plan_apply(recipe, self.module, self.settings)
        except R.WrongModuleError as exc:
            self.detail.setText(str(exc))
            self.buttons["Apply"].setEnabled(False)
            return
        parts = [f"{len(recipe.settings)} settings, saved {recipe.created}."]
        if plan.note():
            parts.append(plan.note())
        if recipe.notes:
            parts.append(recipe.notes)
        text = " ".join(parts)
        if recipe.path:
            text += "\nSaved as: " + os.path.abspath(recipe.path)
        self.detail.setText(text)

    # ------------------------------------------------------------------ what the buttons do
    def current_settings(self) -> Dict[str, Any]:
        """The settings of the module right now, by key."""
        return {key: setting.get() for key, setting in self.settings.items()}

    def save_as(self, name: str, replace: bool = False) -> Optional[R.Recipe]:
        """Save the settings of the module under a name.

        A recipe of the same name is replaced only when ``replace`` is true or the user agrees
        when asked, so a recipe is never lost to a slip of the finger.

        :param name: the name to save under.
        :param replace: replace a recipe of that name without asking.
        :returns: the saved recipe, or ``None`` when the user declined to replace.
        """
        existing = R.find_recipe(self.module, name)
        if existing is not None and not replace:
            if not self.ask_yes_no(
                    "Replace this recipe?",
                    f"A recipe called {existing.name!r} already exists for this module. "
                    f"Replace it with the settings as they are now?"):
                return None
        recipe = R.capture_recipe(self.module, name, self.current_settings())
        if existing is not None:
            recipe.path = existing.path
        R.save_recipe(recipe)
        self.reload(select=recipe.path)
        return recipe

    def save_dialog(self) -> None:
        """Ask for a name and save the current settings under it; report a failure in a message."""
        name = self.ask_name("")
        if not name or not str(name).strip():
            return
        try:
            self.save_as(str(name).strip())
        except (R.RecipeError, OSError) as exc:
            self.warn("Could not save the recipe", str(exc))

    def apply_selected(self, confirmed: Optional[bool] = None) -> Optional[R.ApplyReport]:
        """Apply the chosen recipe to the module.

        A recipe from another version of ParaScape, or one that names settings the module no
        longer has, is applied only after the user agrees, with the reason shown so that the
        choice is informed. A recipe for another module is refused outright.

        :param confirmed: answer the question in advance (``True`` or ``False``) instead of
            asking; ``None`` asks when a question is needed.
        :returns: what was done, or ``None`` when nothing was applied.
        """
        recipe = self.selected()
        if recipe is None:
            return None
        try:
            plan = R.plan_apply(recipe, self.module, self.settings)
        except R.WrongModuleError as exc:
            self.warn("This recipe is for another module", str(exc))
            self.detail.setText(str(exc))
            return None
        if plan.needs_confirmation:
            ok = confirmed if confirmed is not None else self.ask_yes_no(
                "Apply this recipe?", plan.note() + "\n\nApply it anyway?")
            if not ok:
                return None
        self._skipped = {}
        self._notes = []
        try:
            report = R.apply_recipe(recipe, self.module, self.settings, self._write)
        except R.RecipeError as exc:
            self.warn("Could not apply the recipe", str(exc))
            return None
        lines = [f"Applied {recipe.name!r}: {report.applied} settings written."]
        if plan.missing:
            lines.append("Left out, not in this module any more: " + ", ".join(plan.missing) + ".")
        if self._skipped:
            lines.append("Could not be set: " + "; ".join(
                f"{self.settings[k].label or k} ({why})" for k, why in sorted(self._skipped.items()))
                + ".")
        lines.extend(self._notes)
        self.detail.setText(" ".join(lines))
        if self.after_apply is not None:
            self.after_apply()
        self.applied.emit(recipe.name)
        return report

    def _write(self, values: Dict[str, Any]) -> int:
        """Put values into the settings, remembering the ones that could not be used.

        One value that cannot be set (a choice that has gone) must not stop the others, so each is
        tried alone and the reasons are kept for the message to the user.

        :param values: setting key mapped to the value to put in.
        :returns: how many were set.
        """
        count = 0
        for key, value in values.items():
            try:
                note = self.settings[key].set(value)
            except (ValueError, TypeError) as exc:
                self._skipped[key] = str(exc)
                continue
            if note:
                self._notes.append(str(note))
            count += 1
        return count

    def share(self) -> None:
        """Write the chosen recipe to a file chosen by the user, for sending."""
        recipe = self.selected()
        if recipe is None:
            return
        path = self.choose_save(f"{R.slug(recipe.name) or 'recipe'}.json")
        if not path:
            return
        try:
            R.export_recipe(recipe, path)
        except OSError as exc:
            self.warn("Could not write the file", str(exc))

    def import_file(self, replace: bool = False) -> Optional[R.Recipe]:
        """Add a recipe from a file somebody sent, saying why when it cannot be used.

        :param replace: replace a recipe of the same name without asking.
        :returns: the imported recipe, or ``None``.
        """
        path = self.choose_open()
        if not path:
            return None
        try:
            try:
                recipe = R.import_recipe(path, self.module, replace=replace)
            except R.RecipeExistsError as exc:
                if not self.ask_yes_no("Replace this recipe?", f"{exc}. Replace it?"):
                    return None
                recipe = R.import_recipe(path, self.module, replace=True)
        except R.RecipeError as exc:
            self.warn("Could not import the recipe", str(exc))
            return None
        self.reload(select=recipe.path)
        return recipe

    def rename(self, new_name: str = "") -> Optional[str]:
        """Give the chosen recipe a new name.

        :param new_name: the name; the user is asked when it is empty.
        :returns: the new path, or ``None``.
        """
        recipe = self.selected()
        if recipe is None:
            return None
        name = new_name or self.ask_name(recipe.name)
        if not name or not str(name).strip() or str(name).strip() == recipe.name:
            return None
        try:
            path = R.rename_recipe(recipe, str(name))
        except (R.RecipeError, OSError) as exc:
            self.warn("Could not rename the recipe", str(exc))
            return None
        self.reload(select=path)
        return path

    def delete(self) -> bool:
        """Remove the chosen recipe after the user agrees. ``True`` when it was removed."""
        recipe = self.selected()
        if recipe is None:
            return False
        if not self.ask_yes_no("Delete this recipe?",
                               f"Delete the recipe {recipe.name!r}? This cannot be undone."):
            return False
        try:
            removed = R.delete_recipe(recipe)
        except OSError as exc:
            self.warn("Could not delete the recipe", str(exc))
            return False
        self.reload()
        return removed

    # ------------------------------------------------------------------ the real dialogs
    def _ask_name(self, current: str) -> Optional[str]:
        """Ask the user for a recipe name with a text box.

        :param current: the name to start from; empty for a new recipe.
        """
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Name the recipe", "A name for these settings, for example “Tight map”:",
            text=current)
        return name if ok else None

    def _ask_yes_no(self, title: str, text: str) -> bool:
        """Ask a yes or no question in a message box.

        :param title: heading of the box.
        :param text: the question.
        """
        answer = QtWidgets.QMessageBox.question(
            self, title, text,
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No)
        return answer == QtWidgets.QMessageBox.StandardButton.Yes

    def _choose_open(self) -> str:
        """Ask the user for a recipe file to read."""
        path, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self, "Import a recipe", "", "ParaScape recipe (*.json);;All files (*)")
        return path

    def _choose_save(self, suggested: str) -> str:
        """Ask the user where to write a recipe file.

        :param suggested: the file name to offer.
        """
        path, _filter = QtWidgets.QFileDialog.getSaveFileName(
            self, "Share a recipe", suggested, "ParaScape recipe (*.json)")
        return path

    def _warn(self, title: str, text: str) -> None:
        """Tell the user something went wrong, in a message box.

        :param title: heading of the box.
        :param text: what happened and why.
        """
        QtWidgets.QMessageBox.warning(self, title, text)
