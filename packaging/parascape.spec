# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec shared by the Windows and macOS builds.

    pyinstaller packaging/parascape.spec

The committed cache is the whole point of bundling: data/ is 11 MB and makes the installed application
work with no network and no dataset, so it is collected rather than fetched on first run.
"""
from __future__ import annotations

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

ROOT = Path(SPECPATH).resolve().parent
ENTRY = str(ROOT / "packaging" / "parascape_launcher.py")

# PyInstaller's static analysis misses these: pyqtgraph and umap import lazily, and sklearn's
# HDBSCAN and the numba-compiled umap internals are only reachable at run time.
hiddenimports = []
for pkg in ("parascape", "pyqtgraph", "pyqtgraph.opengl", "OpenGL", "OpenGL.platform",
            "sklearn.cluster", "sklearn.manifold", "sklearn.metrics",
            "umap", "numba", "pynndescent", "scipy.stats", "pyarrow"):
    try:
        hiddenimports += collect_submodules(pkg)
    except Exception:
        hiddenimports.append(pkg)

# The committed cache, which is what makes the installed application standalone. It lives INSIDE the
# package (parascape/data), and it is bundled to the same relative place so paths.data_dir() finds it
# by its normal first candidate rather than by a special case for frozen builds.
#
# This pointed at ROOT/data until 2026-08-12. After the cache moved into the package that directory
# stopped existing, so PyInstaller would have bundled nothing and produced an installer whose
# application cannot find its own data -- the failure with no error message that the whole resolver
# exists to prevent.
CACHE = ROOT / "parascape" / "data"
if not CACHE.is_dir():
    raise SystemExit(f"cache not found at {CACHE}; build it with python -m parascape.build_graph "
                     "before packaging, or the installer will ship an application with no data")
datas = [(str(CACHE), "parascape/data")]
for pkg in ("pyqtgraph", "umap", "pynndescent"):
    try:
        datas += collect_data_files(pkg)
    except Exception:
        pass

a = Analysis([ENTRY], pathex=[str(ROOT)], binaries=[], datas=datas,
             hiddenimports=hiddenimports, hookspath=[], runtime_hooks=[],
             excludes=["tkinter", "tensorflow", "torch", "PySide6", "PyQt5"],
             noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="parascape",
          console=False, icon=None)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="parascape")

app = BUNDLE(coll, name="parascape.app", icon=None, bundle_identifier="io.github.einarolafsson.parascape",
             info_plist={"NSHighResolutionCapable": True, "LSMinimumSystemVersion": "11.0"})
