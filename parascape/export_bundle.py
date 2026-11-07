"""One export bundle for a query: tables, figures, methods text and a hashed manifest (64.39).

Everything a reader needs to check a result outside ParaScape, in one folder:

* each table as CSV, marked **measured**, **inferred** or **computed** in the manifest -- an inference
  never leaves the program looking like a measurement;
* each figure as a vector SVG (`figure_export`);
* `methods.md`: the methods paragraph, numbered references and BibTeX, every citation from the dataset
  registry, unconfirmed sources flagged rather than cited (`methods_text`);
* `manifest.json`: the query, organism, ParaScape version, the table fingerprint, and the SHA-256 of
  every file written, so `load` can tell a bundle that was edited afterwards from one that was not.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import pandas as pd

#: measured = observed values; inferred = a claim about a gene that was not observed; computed = derived
#: from measurements (map coordinates, groups, ranks) without being a claim about any gene.
KINDS = ("measured", "inferred", "computed")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(name)).strip("_") or "table"


def write(folder, *, organism: str, query: dict, tables: dict, kinds: dict, columns=(),
          figures: dict | None = None, nodes: pd.DataFrame | None = None, **methods) -> Path:
    """Write the bundle into `folder` (created) and return the manifest path.

    `tables` maps a name to a DataFrame; `kinds` says for each whether it is measured or inferred
    (required, so nothing leaves unlabelled); `columns` are the evidence columns the result rests
    on, which the methods text cites; `figures` maps a name to a widget to draw as SVG.
    """
    from . import __version__, methods_text
    from .searches import fingerprint
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    missing = sorted(set(tables) - set(kinds))
    bad = sorted(n for n, k in kinds.items() if k not in KINDS)
    if missing or bad:
        raise ValueError(f"every table must be marked measured, inferred or computed: missing {missing}, "
                         f"unknown {bad}")
    files = {}
    for name, frame in tables.items():
        path = folder / f"{_safe(name)}.csv"
        pd.DataFrame(frame).to_csv(path, index=False)
        files[path.name] = {"table": name, "kind": kinds[name], "rows": int(len(frame))}
    for name, widget in (figures or {}).items():
        from .figure_export import save
        path = folder / f"{_safe(name)}.svg"
        save(widget, str(path), title=str(name))
        files[path.name] = {"figure": name}
    src = methods_text.sources(columns, organism)
    text = methods_text.full_text(columns, organism, **methods)
    (folder / "methods.md").write_text(text)
    files["methods.md"] = {"methods": True, "cited": len(src.numbered()),
                           "unconfirmed": [s.dataset.key for s in src.found
                                           if not s.confirmed and not s.computed],
                           "unregistered_columns": list(src.unregistered)}
    for name in files:
        files[name]["sha256"] = _sha256(folder / name)
    manifest = {"query": query, "organism": organism, "parascape_version": __version__,
                "table": fingerprint(nodes) if nodes is not None else None,
                "evidence_columns": list(map(str, columns)), "written_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "files": files}
    path = folder / "manifest.json"
    path.write_text(json.dumps(manifest, indent=1, default=str))
    return path


def load(folder) -> dict:
    """Read a bundle back: the manifest, after checking every file still matches its hash."""
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    changed = [name for name, meta in manifest["files"].items()
               if not (folder / name).exists() or _sha256(folder / name) != meta["sha256"]]
    if changed:
        raise ValueError(f"changed or missing since the bundle was written: {changed}")
    return manifest


def ask_folder(parent) -> str:
    """A folder dialog for a bundle; "" when cancelled."""
    from PySide6 import QtWidgets
    return QtWidgets.QFileDialog.getExistingDirectory(parent, "Folder for the bundle") or ""


def reopen(module, folder: str = "") -> dict:
    """Read a bundle back into a Map, Gene or Network module: verify it, then restore its query.

    Map: the selected genes and the colour-by label; Gene: the gene; Network: the focus gene. What was
    computed is not re-run; the bundle files are the record, and the restored view is where to look again.
    """
    if not folder:
        folder = ask_folder(module.shell)
        if not folder:
            return {}
    try:
        manifest = load(folder)
    except (OSError, ValueError) as error:
        module.shell.say(f"bundle refused: {error}")
        return {}
    query = manifest.get("query", {})
    kind = query.get("module")
    if kind == "gene" and query.get("gene"):
        module.show(str(query["gene"]))
    elif kind == "network" and query.get("gene"):
        module.field.setText(str(query["gene"]))
        module._typed()
    elif kind == "map":
        ids = module.window.nodes["gene_id"].astype(str).tolist()
        wanted = [ids.index(g) for g in query.get("selected", []) if g in ids]
        if query.get("colour_by"):
            module._colour_chosen(str(query["colour_by"]))
        if wanted:
            module.select(wanted, "from the bundle")
    else:
        module.shell.say(f"this bundle is for the {kind or 'unknown'} module")
        return manifest
    module.shell.say(f"bundle opened: {folder}")
    return manifest
