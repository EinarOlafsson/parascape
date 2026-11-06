# 83 · Move ParaScape from PyQt6 to PySide6

**Status: done (2026-11-05).** The precondition for the real overhaul of instruction 74.

## Why

Einar asked for the home screen and navigator to be spaCR's, with ParaScape's modules mapped onto them.
Scoping that turned up the thing that decides how it can be done: **spaCR is PySide6 and ParaScape is PyQt6**
-- 154 PySide6 imports in the five files in question alone. So spaCR's visual code cannot be copied across;
on PyQt6 it has to be reimplemented from the structure.

Two further facts settled it:

1. **Licence.** PyQt6 is **GPL-3.0-only** (read from the installed package, 6.11.0) or a paid commercial
   licence from Riverbank. PySide6 is LGPL-3.0 / GPL-2.0 / GPL-3.0. ParaScape declares itself MIT, so an
   MIT package that cannot run without a GPLv3 library already put anyone redistributing the combination
   under GPL terms. Moving to PySide6 removes that, and removes the awkwardness of shipping ported spaCR
   files -- spaCR's own licence is BSD-3 with noncommercial additional terms.
2. **Nothing visual is lost.** Verified rather than assumed: with `PYQTGRAPH_QT_LIB=PySide6`, pyqtgraph
   genuinely binds PySide6 and `GLViewWidget`, `GLScatterPlotItem`, `pyqtgraph.opengl.shaders` and
   PyOpenGL all come up. ParaScape's effects live below the binding -- pyqtgraph.opengl, PyOpenGL, and
   38 lines of its own GLSL in `sprite.py`. The strongest evidence is on the same machine: spaCR does
   more visually (ambient backdrops, fractals, live video, 3D) entirely on PySide6.

## Measured scope

| | count |
|---|---|
| files importing PyQt6 | 113 (parascape 40, tests 51, scripts 22) |
| `from PyQt6` | 496 |
| `pyqtSignal` | 423 |
| `pyqtSlot` / `pyqtProperty` | 4 / 3 |
| `PYQTGRAPH_QT_LIB` references | 59 |
| **`sip.` call sites** | **4** |

What makes this tractable: **no `exec_()`, no `toPyObject`, no `QVariant`, no `QStringList`, and zero
unscoped enums.** PyQt6 forbids the unscoped form, so everything is already fully scoped -- which is
exactly what PySide6 accepts. The port is a rename plus four hand edits.

## The hand edits

`PyQt6.sip` has no PySide6 equivalent by that name; `shiboken6` does the same job with inverted sense:

| PyQt6 | PySide6 |
|---|---|
| `sip.unwrapinstance(obj)` | `shiboken6.getCppPointer(obj)[0]` |
| `sip.isdeleted(obj)` | `not shiboken6.isValid(obj)` |

Sites: `parascape/glass.py:250,251`, `tests/test_functional_agreement_ui.py:24`,
`tests/test_app_gene_evidence.py:126`.

## Stages

1. Install PySide6 **alongside** PyQt6, so nothing breaks mid-flight. Done: 6.12.0.
2. Codemod `parascape/`, `tests/`, `scripts/`: imports and the `pyqt*` factory names.
3. The four `sip` sites, by hand, with the inverted sense checked.
4. `PYQTGRAPH_QT_LIB` to `PySide6` everywhere, including CI workflows and docs.
5. `pyproject.toml`: `PyQt6` out, `PySide6` in; `pytest-qt` works with either.
6. Full suite. Only then remove PyQt6 from the environment.

## Known risks, to watch rather than assume away

- **Ownership and garbage collection.** PySide6 can collect a parentless widget that PyQt6 kept alive.
  This is the likeliest source of new failures and it will not show up as an import error.
- **Flag construction.** `QtCore.Qt.Edge(0)` at `glass.py:352` builds a flag from an int; PySide6's true
  enums should accept it, but it is on the list to confirm.
- **Enum/int comparisons.** PySide6 6.4+ uses real Python enums, so anything comparing a flag to a bare
  integer can change behaviour silently.
- **The display painting test** (`test_display.py::test_the_background_never_climbs_on_top_of_a_panel`)
  is already failing in CI before this work. It must not be used to judge the migration, and the
  migration must not be used to excuse it.

## Gate

The full suite passes at its pre-migration count (6,505 passed / 41 skipped as of 2026-11-05, minus the
one already-red display test), `python -m parascape` opens with its 3D map rendering, and PyQt6 appears
nowhere in `pyproject.toml` or the imports.

## Then, and only then

Instruction 74's overhaul proper: port spaCR's `HomePage`, `AppTile`, theme tokens, navigator aside and
the `AppScreen` module shell -- which on PySide6 transfer rather than need rewriting -- map ParaScape's
modules onto the tiles, move the UMAP and the parasite schematic off the landing page into a module, and
rebuild the network as a navigable universe.


## Done 2026-11-05: 6,550 passed, 41 skipped, 0 failed, single process, 30m24s

Against the pre-migration baseline of 6,505 passed / 41 skipped. The extra 45 are the tests added
today elsewhere (contrast ladder, the derived label floor, the download fixes, the import resolver).

### Five real defects, every one a PyQt6 courtesy that PySide6 does not extend

1. **Completion handlers ran on the worker thread.** `_run` documents GUI-thread delivery and
   silently stopped providing it: `done.connect(relay)` bound a plain closure, which carries no
   QObject receiver, so PySide6 has no receiving thread to queue to and ran it on the emitter. PyQt6
   wrapped connected callables in a proxy belonging to the connecting thread. Qt was reporting it as
   `QThread::wait: Thread tried to wait on itself` -- `relay` joining the thread it was running on.
   **This would have shipped**: touching widgets from a worker is the crash that happens on a slow
   machine and not a fast one. Fixed by making the receiver a bound `@Slot` on the panel.
2. **`_Chrome.eventFilter` received an event after its Python attributes were gone**, raising inside
   the Qt event loop where no caller can catch it. Guarded.
3. **Modal menus that tests believed they had stubbed**, at three call sites. `QMenu.exec` has a
   static overload, so PySide6 exposes it as a mixed static/instance method that a plain assignment
   cannot replace; the real modal loop then ran with no user to dismiss it, so the suite HUNG rather
   than failed. Verified by probe that `QDialog.exec` and static patches
   (`QFileDialog.getSaveFileName`, `QApplication.clipboard`, `QDesktopServices.openUrl`,
   `QGuiApplication.primaryScreen`) are all unaffected, so the blast radius was three sites and not
   the eleven first feared. Fixed with a `_show_menu` seam on the panel and the window.
4. **`QWidget.render(painter)`** -- PySide6 has no single-argument overload; the targetOffset is
   required.
5. **`shiboken6` imported unguarded and undeclared.** A packaging bug introduced by this migration
   and caught by `tests/test_publishing.py`, which exists for exactly that. Declared in pyproject.

### The premise that was wrong

Single-process whole-suite runs were believed impossible on this machine -- a 30-minute wedge, zero
CPU, main thread in `do_poll`. That was defect 3: a modal menu loop. The slicing built around it was
scaffolding for a problem that did not exist, and it is why the same wall was hit three times before.

**`do_poll` with no CPU in a Qt suite should read as "a modal loop is open", not "the machine cannot
do this".**

### Two environment lessons, both already written down and both applied too late

- Long jobs need their own `systemd-run --user --scope` with `MemoryMax`. Inheriting the VS Code
  scope means `systemd-oomd` kills them under pressure from other sessions, and two kills at matching
  elapsed times are indistinguishable from a reproducible hang. Two investigations were lost to this.
- A command piped through `tail` buffers everything, so its output file reads as empty until it
  finishes. Progress read from such a file is not progress.
