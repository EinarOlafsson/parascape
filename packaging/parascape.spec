# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec shared by the Windows and macOS builds.

    pyinstaller packaging/parascape.spec

Application icons and the data descriptor ship with the executable. Scientific
data is downloaded from GitHub on first launch and then works offline.
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

# Dataset files are downloaded from a versioned GitHub release on first launch.
ASSETS = ROOT / "parascape" / "data" / "icons"
DESCRIPTOR = ROOT / "parascape" / "data_manifest.json"
if not ASSETS.is_dir() or not DESCRIPTOR.is_file():
    raise SystemExit("Missing icons or data release descriptor; restore the source assets before packaging")
datas = [(str(ASSETS), "parascape/data/icons"), (str(DESCRIPTOR), "parascape")]
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
