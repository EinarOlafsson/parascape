#!/usr/bin/env python3
"""Build the wheel from the committed tree, install it apart from the checkout, and walk it (instruction 64.40).

`git archive HEAD` -> `pip wheel` (no dependencies) -> `pip install --target` into a scratch folder, then a fresh
interpreter in another directory (isolated mode, the checkout not on the path) imports the INSTALLED package and
walks the goals: data resolves, both parasite tables and the held-out record load, the Gene module reads a gene,
the gene and class records open, and a strategy card renders. The data cache comes from `PARASCAPE_CACHE` (default:
this checkout's `parascape/data`), standing in for `parascape-data --download`, whose download path has its own
tests. Uses the current Python environment for dependencies; creates no environment.

    python scripts/check_installed_wheel.py [--work DIR] [--cache DIR]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

WALK = r'''
import os, sys, json
from PySide6 import QtWidgets
app = QtWidgets.QApplication([])
import parascape
from parascape import paths, organisms as O, track_record as T, claims as C, strategies as S
out = {"parascape_from": os.path.dirname(parascape.__file__), "data_dir": paths.data_dir()}
assert os.path.realpath(out["parascape_from"]).startswith(os.path.realpath(sys.argv[1])), "imported the checkout"
import pandas as pd
out["tg_genes"] = len(pd.read_parquet(O.nodes_path("Tg")))
out["pf_genes"] = len(pd.read_parquet(O.nodes_path("Pf")))
L = T.shipped("Tg"); out["track_record_rows_tg"] = len(L)
out["claims_tg"] = int((C.shipped()["organism"].astype(str) == "Tg").sum())
from parascape.app import Window
from parascape.landing_install import open_module
w = Window(); w.show(); app.processEvents()
open_module(w, "gene"); m = w.gene_module.module
gene = str(L["gene_id"].iloc[0]); m.show(gene)
out["gene_card_has_gene"] = gene in m.known.toHtml()
out["gene_methods_section"] = "methods were tested" in m.methods.toHtml()
out["held_out_record_for_gene"] = bool(T.gene_html(gene, "Tg"))
label = str(L[L["target"].astype(str) == "compartment"]["truth"].iloc[0])
out["class_page"] = bool(T.class_html("compartment", label, "Tg"))
out["strategy_card"] = "Scorecard" in w.strategy_panel.scorecard_html(S.catalog()[0])
w.close()
print(json.dumps(out))
'''


def main(argv=None) -> int:
    """Build, install apart, walk; print the walk record and exit non-zero on any failed step."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", default=None)
    ap.add_argument("--cache", default=str(ROOT / "parascape" / "data"))
    args = ap.parse_args(argv)
    work = Path(args.work or tempfile.mkdtemp(prefix="parascape-wheel-"))
    src, wheel, site, run = (work / n for n in ("src", "wheel", "site", "run"))
    for d in (src, wheel, site, run):
        d.mkdir(parents=True, exist_ok=True)
    archive = subprocess.run(["git", "archive", "HEAD"], cwd=ROOT, check=True, capture_output=True).stdout
    subprocess.run(["tar", "-x", "-C", str(src)], input=archive, check=True)
    subprocess.run([sys.executable, "-m", "pip", "wheel", ".", "--no-deps", "--no-build-isolation", "-w", str(wheel),
                    "-q"], cwd=src, check=True)
    (whl,) = list(wheel.glob("*.whl"))
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-compile", "--target", str(site),
                    str(whl), "-q"], check=True)
    env = dict(os.environ, PARASCAPE_CACHE=args.cache, PARASCAPE_STATE=str(run / "state"),
               XDG_CONFIG_HOME=str(run / "config"), XDG_CACHE_HOME=str(run / "cache"), PARASCAPE_AMBIENT_DEFAULT="0",
               QT_QPA_PLATFORM=os.environ.get("QT_QPA_PLATFORM", "offscreen"))
    code = f"import sys; sys.path.insert(0, {str(site)!r}); sys.argv=['walk', {str(site)!r}]; exec({WALK!r})"
    done = subprocess.run([sys.executable, "-I", "-c", code], cwd=run, env=env, capture_output=True, text=True)
    line = next((l for l in done.stdout.splitlines()[::-1] if l.startswith("{")), "")
    record = json.loads(line) if line else {"error": done.stderr[-2000:]}
    record.update(wheel=whl.name, wheel_bytes=whl.stat().st_size)
    print(json.dumps(record, indent=1))
    ok = done.returncode == 0 and all(v for k, v in record.items() if isinstance(v, bool))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
