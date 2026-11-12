"""Recipes: named, reusable and shareable bundles of module settings.

This is the ParaScape counterpart of ``spacr/qt/recipes.py``, and it exists for the same reason.
A lab or a person often returns to one set of settings, such as "Map by expression and
chromatin, tight spacing". A table of settings saved in a folder has no name, no record of the
module it belongs to and no record of the version of ParaScape that wrote it.

A recipe records all four:

* it has a **name** that the user chose;
* it knows its **module**, so applying a Map recipe to Network is REFUSED rather than writing
  keys that do not apply there (the keys that two modules happen to share are exactly the
  generic ones, so a cross-module apply would write only those and still report success);
* it records the **version** of ParaScape it was captured with, and applying it under a
  different one says so and lists, by name, the settings that the module no longer has;
* it is **one file**, so sharing a recipe is sending a file.

Storage is ``~/.parascape/recipes/<module>/<slug>.json``, or the folder named by the
environment variable :data:`RECIPE_DIR_ENV`, with one subfolder for each module so a listing is
already scoped to that module and a whole module of recipes can be handed over by copying one
folder.

The file is plain JSON::

    {
      "parascape_recipe": 1,
      "name": "Tight map by expression",
      "module": "map",
      "parascape_version": "0.0.11.0",
      "created": "2026-11-12T09:41:07",
      "notes": "",
      "settings": {"n_neighbors": 15, "min_dist": 0.05, ...}
    }

``parascape_recipe`` is the version of the file format and ``parascape_version`` is the version
of the application. They change for unrelated reasons, so a reader accepts any file whose format
version it can read, whatever application version wrote it.

A recipe is always a bundle of settings, here and in spaCR. Until 2026-11-12 this word also named
the way of asking one biological question in the Infer module, with the measurements allowed to
answer it. That is now a **question template** and lives in :mod:`parascape.question_templates`;
the old names still answer, with a message, from the bottom of this file.

Nothing here imports Qt, so the logic is tested without a display and a script can read, apply
or share recipes. The panel that shows them is :mod:`parascape.recipes_panel`.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
import tempfile
import unicodedata
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

#: Environment variable that moves the recipe folder, so a lab can keep its recipes on shared
#: storage and a test can keep them out of the home folder.
RECIPE_DIR_ENV = "PARASCAPE_RECIPE_DIR"

#: Version of the file format. Raised only when a reader would need to behave differently.
FORMAT_VERSION = 1

#: Longest stem of a file name taken from a recipe name, so a pasted paragraph cannot make a path
#: the file system refuses.
_SLUG_LIMIT = 60

_NOT_ALNUM = re.compile(r"[^a-z0-9]+")

#: File names that some systems refuse whatever their extension, which matters because a recipe
#: file is meant to be sent to other people.
_RESERVED = frozenset(
    {"con", "prn", "aux", "nul"}
    | {f"com{n}" for n in range(1, 10)}
    | {f"lpt{n}" for n in range(1, 10)})


class RecipeError(ValueError):
    """Anything wrong with a recipe: unreadable, not a recipe, or not for this module.

    A subclass of :class:`ValueError` so a caller that only wants to show the message can catch
    one type, while a caller that must tell the cases apart can catch the narrower ones.
    """


class NotARecipeError(RecipeError):
    """The file is unreadable, is not valid JSON or is not a ParaScape recipe."""


class WrongModuleError(RecipeError):
    """The recipe was captured in a different module and is refused here."""


class RecipeExistsError(RecipeError):
    """A recipe of that name is already saved for this module."""


def parascape_version() -> str:
    """The running version of ParaScape, or ``"unknown"`` when it cannot be read.

    Read here rather than imported at the top so that a recipe can still be loaded and shown in
    an environment where the package metadata is broken, which is when a record of the version
    matters most.
    """
    try:
        from . import __version__
        return str(__version__)
    except Exception:
        return "unknown"


def recipes_dir(module: Optional[str] = None, create: bool = True) -> str:
    """The folder that holds recipes, optionally the one for a single module.

    ``~/.parascape/recipes``, or the folder named by :data:`RECIPE_DIR_ENV`. There is one
    subfolder for each module, so a listing is scoped without reading every file and a user can
    hand over the recipes of one module by copying one folder.

    :param module: key of the module whose subfolder is wanted; omit for the top folder.
    :param create: make the folder when it is missing. A reader that only wants to look passes
        ``False`` so that looking never creates anything.
    """
    override = os.environ.get(RECIPE_DIR_ENV, "").strip()
    root = (os.path.abspath(os.path.expanduser(override)) if override
            else os.path.join(os.path.expanduser("~"), ".parascape", "recipes"))
    if module:
        root = os.path.join(root, slug(module) or "unknown")
    if create:
        os.makedirs(root, exist_ok=True)
    return root


def slug(text: str) -> str:
    """A file name stem that is safe for any file system, made from a recipe name.

    The name the user typed is kept inside the file and is what they read, so the stem only has
    to be safe and typeable. Accented letters are folded to plain ones ("Café" becomes "cafe"),
    every other run of characters becomes one hyphen, path separators and dots can never survive
    (so a name cannot reach outside the folder), and the length is limited. A name with no
    plain letters at all ("日本語", "???") gets a short fingerprint of itself instead of an
    empty stem, so two such names do not share a file. A reserved system name gets a suffix.

    Two different names can still give the same stem ("A b" and "a-b"); :func:`save_recipe`
    resolves that, so this function alone is not a unique key.

    :param text: the name to turn into a stem.
    """
    raw = str(text or "")
    folded = unicodedata.normalize("NFKD", raw).encode("ascii", "ignore").decode("ascii")
    stem = _NOT_ALNUM.sub("-", folded.lower()).strip("-")[:_SLUG_LIMIT].strip("-")
    if not stem:
        if not raw.strip():
            return ""
        stem = "recipe-" + hashlib.sha1(raw.strip().encode("utf-8")).hexdigest()[:8]
    if stem in _RESERVED:
        stem += "-recipe"
    return stem


def _same_name(a: str, b: str) -> bool:
    """Whether two recipe names are the same name, ignoring case and surrounding spaces.

    Used so that "Tight map" and "tight map " are one recipe rather than two files that differ
    only in how they were typed.
    """
    return str(a).strip().casefold() == str(b).strip().casefold()


def _jsonable(value: Any) -> Any:
    """Turn the odd values a settings widget can hand over into plain JSON values.

    Sets become sorted lists (tuples already become lists in JSON) and numpy scalars become Python numbers. Anything else becomes
    its text, which is better than refusing to save a recipe over one stray value.
    """
    if isinstance(value, (set, frozenset)):
        return sorted(value, key=str)
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return item()
        except Exception:
            pass
    return str(value)


@dataclass
class Recipe:
    """One named bundle of settings for one module.

    :param name: the name the user chose; its slug becomes the stem of the file.
    :param module: key of the module the settings were captured in; it picks the subfolder and
        is checked when the recipe is applied.
    :param settings: setting key mapped to the value it had.
    :param parascape_version: version of ParaScape that captured it; filled in on save when empty.
    :param created: ISO timestamp, to the second; filled in on save when empty and used to list
        the newest first.
    :param notes: free text kept with the recipe.
    :param path: the file it was read from or written to. Not written into the file, because a
        file that moves would then lie about where it is.
    """

    name: str
    module: str
    settings: Dict[str, Any] = field(default_factory=dict)
    parascape_version: str = ""
    created: str = ""
    notes: str = ""
    path: str = ""

    def to_json(self) -> Dict[str, Any]:
        """The mapping that is written to disk.

        The version and the time are stamped here when missing, so a recipe built by hand and a
        recipe captured from a module are written the same way.
        """
        return {
            "parascape_recipe": FORMAT_VERSION,
            "name": self.name,
            "module": self.module,
            "parascape_version": self.parascape_version or parascape_version(),
            "created": self.created or _dt.datetime.now().isoformat(timespec="seconds"),
            "notes": self.notes,
            "settings": json.loads(json.dumps(self.settings, default=_jsonable)),
        }

    @classmethod
    def from_json(cls, data: Any, path: str = "") -> "Recipe":
        """Build a recipe from parsed JSON, refusing anything that is not one.

        An arbitrary JSON file, such as a ``package.json``, must never be applied as settings, so
        the format marker is required. A format version newer than this build reads is refused
        with a message that says so, rather than guessed at. A file that does not say which
        module it is for is refused too, because a recipe that cannot be checked against a
        module cannot be applied safely.

        :param data: the parsed JSON.
        :param path: the file it came from, kept on the recipe.
        :raises NotARecipeError: when the data is not a usable recipe.
        """
        if not isinstance(data, dict) or "parascape_recipe" not in data:
            raise NotARecipeError("this is not a ParaScape recipe (it has no recipe format marker)")
        version = data["parascape_recipe"]
        if isinstance(version, bool) or not isinstance(version, int) or version < 1:
            raise NotARecipeError(f"the recipe format marker {version!r} is not a format number")
        if version > FORMAT_VERSION:
            raise NotARecipeError(
                f"this recipe uses file format {version}, which is newer than this ParaScape "
                f"reads (it reads formats up to {FORMAT_VERSION}); update ParaScape to use it")
        module = data.get("module")
        if not isinstance(module, str) or not module.strip():
            raise NotARecipeError("this recipe does not say which module it is for")
        settings = data.get("settings")
        if not isinstance(settings, dict):
            raise NotARecipeError("this recipe has no settings")
        return cls(
            name=str(data.get("name") or "Untitled"),
            module=module.strip(),
            settings=dict(settings),
            parascape_version=str(data.get("parascape_version") or "unknown"),
            created=str(data.get("created") or ""),
            notes=str(data.get("notes") or ""),
            path=path,
        )


# ---------------------------------------------------------------------------- files

def _write_text(path: str, text: str) -> None:
    """Write a file so that it is either the old file or the new one, never half of one.

    The text goes to a temporary file in the same folder and is moved into place, because a
    recipe cut short by a crash would be the corrupt file that the next listing has to report.

    :param path: destination file.
    :param text: what to write.
    """
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".recipe-", suffix=".tmp", dir=folder)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _read_name(path: str) -> Optional[str]:
    """The name stored in a recipe file, or ``None`` when the file cannot be read as a recipe.

    :param path: the file to look in.
    """
    try:
        return load_recipe(path).name
    except RecipeError:
        return None


def _free_path(folder: str, name: str, own_path: str = "") -> str:
    """The file that a recipe of this name should be written to, without taking anyone else.

    The stem of the name is tried first. A file already there is replaced only when it is the
    same recipe: the file this recipe was loaded from, or one that holds the same name. A file
    that holds a different name (two names that slug alike, such as "Tight map" and
    "tight-map"), or one that cannot be read at all, is left alone and the next free stem
    ("-2", "-3", ...) is used. Replacing a corrupt file would destroy the evidence that the user
    needs to see.

    :param folder: the module folder.
    :param name: the display name of the recipe.
    :param own_path: file the recipe already lives in, if any.
    """
    base = slug(name) or "recipe"
    own = os.path.abspath(own_path) if own_path else ""
    n = 1
    while True:
        stem = base if n == 1 else f"{base}-{n}"
        path = os.path.join(folder, f"{stem}.json")
        if not os.path.exists(path):
            return path
        if own and os.path.abspath(path) == own:
            return path
        existing = _read_name(path)
        if existing is not None and _same_name(existing, name):
            return path
        n += 1


def save_recipe(recipe: Recipe, directory: Optional[str] = None) -> str:
    """Write a recipe to its module folder and return the path of the file.

    A recipe with the same name as one already saved replaces it (saving again is how a recipe is
    updated). A different name that happens to give the same file stem is given its own file, so
    saving never silently overwrites another recipe. The version and the time are stamped on the
    recipe when empty.

    :param recipe: the bundle to write; its ``path`` is set to the file written.
    :param directory: write here instead of the module folder, for sharing.
    :raises RecipeError: when the recipe has no name or no module.
    """
    if not str(recipe.name or "").strip():
        raise RecipeError("a recipe needs a name")
    if not str(recipe.module or "").strip():
        raise RecipeError("a recipe needs the module it is for")
    recipe.name = str(recipe.name).strip()
    recipe.parascape_version = recipe.parascape_version or parascape_version()
    recipe.created = recipe.created or _dt.datetime.now().isoformat(timespec="seconds")
    folder = directory or recipes_dir(recipe.module)
    os.makedirs(folder, exist_ok=True)
    path = _free_path(folder, recipe.name, recipe.path if not directory else "")
    _write_text(path, json.dumps(recipe.to_json(), indent=2, sort_keys=True) + "\n")
    recipe.path = path
    return path


def load_recipe(path: str) -> Recipe:
    """Read one recipe file, saying plainly what is wrong with it if it cannot be read.

    :param path: the file.
    :raises NotARecipeError: when the file cannot be opened, is not valid JSON, or is not a
        usable recipe. The message names the file and the reason, so a listing can show it.
    """
    name = os.path.basename(str(path))
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except json.JSONDecodeError as exc:
        raise NotARecipeError(f"{name} is not valid JSON ({exc.msg}, line {exc.lineno})") from exc
    except UnicodeDecodeError as exc:
        raise NotARecipeError(f"{name} is not a text file") from exc
    except OSError as exc:
        raise NotARecipeError(f"{name} could not be opened ({exc.strerror or exc})") from exc
    try:
        return Recipe.from_json(data, path=str(path))
    except NotARecipeError as exc:
        raise NotARecipeError(f"{name}: {exc}") from exc


@dataclass
class Listing:
    """The result of reading a recipe folder: what was usable and what was not.

    :param recipes: the readable recipes of the module, newest first.
    :param problems: ``(path, reason)`` for every file that was not usable. Reported rather than
        dropped, because a recipe that has silently vanished from the list looks like a recipe
        that was never saved.
    """

    recipes: List[Recipe] = field(default_factory=list)
    problems: List[Tuple[str, str]] = field(default_factory=list)


def scan_recipes(module: str) -> Listing:
    """Read the recipe folder of a module and report every file that could not be used.

    A file is a problem when it is not valid JSON, is not a recipe, or is a recipe that says it
    is for a different module than the folder it sits in (typically copied into the wrong
    folder by hand). Looking never creates the folder.

    :param module: key of the module.
    """
    out = Listing()
    folder = recipes_dir(module, create=False)
    if not os.path.isdir(folder):
        return out
    for entry in sorted(os.listdir(folder)):
        if not entry.endswith(".json") or entry.startswith("."):
            continue
        path = os.path.join(folder, entry)
        if not os.path.isfile(path):
            continue
        try:
            recipe = load_recipe(path)
        except RecipeError as exc:
            out.problems.append((path, str(exc)))
            continue
        if recipe.module != module:
            out.problems.append((
                path, f"{entry} is for the {recipe.module!r} module but is filed under "
                      f"{module!r}; it was not listed"))
            continue
        out.recipes.append(recipe)
    out.recipes.sort(key=lambda r: (r.created, r.name), reverse=True)
    return out


def list_recipes(module: str) -> List[Recipe]:
    """The readable recipes of a module, newest first.

    Use :func:`scan_recipes` when the files that could not be read must be reported too.

    :param module: key of the module.
    """
    return scan_recipes(module).recipes


def find_recipe(module: str, name: str) -> Optional[Recipe]:
    """The saved recipe of a module with this name, or ``None``.

    Used to ask "replace the existing one?" before a save, instead of replacing it unasked.

    :param module: key of the module.
    :param name: the name to look for; case and surrounding spaces are ignored.
    """
    for recipe in list_recipes(module):
        if _same_name(recipe.name, name):
            return recipe
    return None


def delete_recipe(recipe: Recipe) -> bool:
    """Remove the file of a recipe. ``True`` when a file was removed.

    :param recipe: the recipe whose ``path`` is deleted; with no path nothing is removed.
    """
    if not recipe.path or not os.path.isfile(recipe.path):
        return False
    os.remove(recipe.path)
    return True


def rename_recipe(recipe: Recipe, new_name: str) -> str:
    """Give a saved recipe a new name and move its file to match.

    The new file is written before the old one is removed, so a failure part way leaves the
    recipe under its old name rather than under neither.

    :param recipe: the saved recipe; its ``name`` and ``path`` are updated.
    :param new_name: the new name; surrounding spaces are dropped.
    :returns: the new path.
    :raises RecipeError: when the name is empty.
    :raises RecipeExistsError: when another recipe of the module already has that name.
    """
    name = str(new_name or "").strip()
    if not name:
        raise RecipeError("a recipe needs a name")
    if not _same_name(name, recipe.name):
        other = find_recipe(recipe.module, name)
        if other is not None:
            raise RecipeExistsError(f"a recipe called {name!r} already exists for this module")
    old_name, old_path = recipe.name, recipe.path
    recipe.name = name
    folder = os.path.dirname(old_path) if old_path else recipes_dir(recipe.module)
    try:
        # The old file is not offered as "own" so that a rename never lands on it by the stem
        # alone: it is replaced only if the stem is genuinely free or holds this same name.
        recipe.path = ""
        path = save_recipe(recipe, directory=folder)
    except Exception:
        recipe.name, recipe.path = old_name, old_path
        raise
    if old_path and os.path.abspath(old_path) != os.path.abspath(path) and os.path.isfile(old_path):
        os.remove(old_path)
    return path


def capture_recipe(module: str, name: str, settings: Dict[str, Any], notes: str = "") -> Recipe:
    """Build a recipe from the settings a module has right now.

    Every value is kept, including the ones the user never touched, so the recipe reproduces the
    same result after a default changes.

    :param module: key of the module the settings come from.
    :param name: the name the user chose.
    :param settings: setting key mapped to its current value.
    :param notes: free text to keep with it.
    :raises RecipeError: when the name or the module is empty.
    """
    if not str(name or "").strip():
        raise RecipeError("a recipe needs a name")
    if not str(module or "").strip():
        raise RecipeError("a recipe needs the module it is for")
    return Recipe(name=str(name).strip(), module=str(module).strip(), settings=dict(settings),
                  notes=str(notes or ""))


def export_recipe(recipe: Recipe, path: str) -> str:
    """Write one recipe as a single file anywhere, for sending to someone.

    :param recipe: the recipe to share.
    :param path: destination file.
    :returns: the path written.
    """
    _write_text(path, json.dumps(recipe.to_json(), indent=2, sort_keys=True) + "\n")
    return path


def import_recipe(path: str, module: str, replace: bool = False) -> Recipe:
    """Read a shared recipe file and save it for a module.

    A recipe for another module is refused, because it would be listed here and then be refused
    again on every apply. An unreadable file is reported with the reason.

    :param path: the file that was sent.
    :param module: key of the module importing it.
    :param replace: replace a recipe of the same name instead of refusing.
    :raises NotARecipeError: when the file is not a usable recipe.
    :raises WrongModuleError: when it is for a different module.
    :raises RecipeExistsError: when the name is taken and ``replace`` is false.
    """
    recipe = load_recipe(path)
    if recipe.module != module:
        raise WrongModuleError(
            f"{recipe.name!r} is a recipe for the {recipe.module!r} module, not {module!r}")
    existing = find_recipe(module, recipe.name)
    if existing is not None and not replace:
        raise RecipeExistsError(f"a recipe called {recipe.name!r} already exists for this module")
    recipe.path = existing.path if existing is not None else ""
    save_recipe(recipe)
    return recipe


# ---------------------------------------------------------------------------- applying

@dataclass(frozen=True)
class ApplyPlan:
    """What applying a recipe to a module would do, worked out before anything is written.

    :param name: name of the recipe.
    :param made_with: version of ParaScape that captured it.
    :param running: version of ParaScape that is running.
    :param settings: the settings the module has, which are the ones that will be written.
    :param missing: names of settings in the recipe that the module no longer has, sorted. They
        are skipped, and listed in full so the user can judge whether the recipe still means
        what it did.
    """

    name: str
    made_with: str
    running: str
    settings: Dict[str, Any]
    missing: Tuple[str, ...] = ()

    @property
    def version_differs(self) -> bool:
        """Whether the recipe was captured with another version than the one running."""
        return self.made_with != self.running

    @property
    def needs_confirmation(self) -> bool:
        """Whether the user should be told something before it is applied.

        True for a version difference or for any setting that no longer exists. A recipe from the
        same version with every setting present is applied without a question.
        """
        return self.version_differs or bool(self.missing)

    def note(self) -> str:
        """A plain sentence for the user: the version gap and the settings that are gone.

        Empty when there is nothing to say. The settings are named in full rather than counted,
        because "3 settings are ignored" gives no way to judge whether they matter.
        """
        parts = []
        if self.version_differs:
            parts.append(
                f"This recipe was saved with ParaScape {self.made_with}; you are running "
                f"{self.running}.")
        if self.missing:
            names = ", ".join(self.missing)
            many = len(self.missing) != 1
            parts.append(
                f"{'These settings are' if many else 'This setting is'} not in this module any "
                f"more and will be left out: {names}.")
        return " ".join(parts)


def plan_apply(recipe: Recipe, module: str, known_keys: Iterable[str],
               running: Optional[str] = None) -> ApplyPlan:
    """Work out what applying a recipe to a module would do, or refuse it.

    The refusal comes first and is absolute: a recipe from another module is never partly
    applied.

    :param recipe: the recipe to apply.
    :param module: key of the module it would be applied to.
    :param known_keys: the names of every setting the module has now.
    :param running: override the running version (tests).
    :raises WrongModuleError: when the recipe was captured in a different module.
    """
    if recipe.module != module:
        raise WrongModuleError(
            f"This recipe is for the {recipe.module!r} module, not {module!r}, so it was not "
            f"applied.")
    known = set(known_keys)
    keep = {k: v for k, v in recipe.settings.items() if k in known}
    missing = tuple(sorted(k for k in recipe.settings if k not in known))
    return ApplyPlan(
        name=recipe.name,
        made_with=recipe.parascape_version or "unknown",
        running=running if running is not None else parascape_version(),
        settings=keep,
        missing=missing)


@dataclass(frozen=True)
class ApplyReport:
    """What applying a recipe did.

    :param plan: the plan that was carried out.
    :param applied: how many settings the module reported as written.
    """

    plan: ApplyPlan
    applied: int


def apply_recipe(recipe: Recipe, module: str, known_keys: Iterable[str],
                 write: Callable[[Dict[str, Any]], int],
                 running: Optional[str] = None) -> ApplyReport:
    """Apply a recipe to a module through the function the module gives for writing settings.

    Planning and writing are one call so that no caller can write the settings and skip the
    module check. The module is passed in as a key and a writer rather than as a window, which
    keeps this function free of Qt.

    :param recipe: the recipe to apply.
    :param module: key of the module it is applied to.
    :param known_keys: the names of every setting the module has now.
    :param write: called once with the settings to write; returns how many it wrote.
    :param running: override the running version (tests).
    :raises WrongModuleError: when the recipe belongs to a different module; nothing is written.
    """
    plan = plan_apply(recipe, module, known_keys, running=running)
    applied = int(write(dict(plan.settings)) or 0)
    return ApplyReport(plan=plan, applied=applied)


# --------------------------------------------------------------------------- names that moved
#: What this module used to export, before it held settings recipes. Each name maps to where the same
#: thing is now. ``Recipe`` is missing on purpose: it still exists here with a different meaning, so
#: it cannot be redirected, and a call written for the old one fails on its own arguments.
_MOVED = {
    "RecipeResult": "TemplateResult", "RecipeStore": "TemplateStore",
    **{name: name for name in (
        "Closure", "CONTROL_ENRICHMENT", "GOOD_REACH", "IMPORT_GROUPS", "MIN_LABEL",
        "STRONG_ENRICHMENT", "address_blocks", "address_label", "build_map", "build_model", "close",
        "dominant_by_cluster", "infer", "input_categories", "label_series", "recovery_by_cluster",
        "register_import", "relevance", "run", "same_category", "same_import")},
}


def __getattr__(name: str):
    """Say where a name that used to live here has gone, instead of reporting it as unknown.

    Code written before 2026-11-12 imports ``run`` or ``close`` from here. It must fail, because
    the settings recipes have no such thing, and the failure should say what to import instead.
    """
    if name in _MOVED:
        raise ImportError(
            f"parascape.recipes.{name} moved: the recipes of a biological question are now question "
            f"templates, so use parascape.question_templates.{_MOVED[name]}. A recipe in "
            f"parascape.recipes is a bundle of module settings.")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
