"""Every paper panel, exported with what made it: PNG, vector PDF, the plotted data, the statistics,
and the source data (or where it lives, when too large to copy).

Einar, 2026-11-20: the parascape-paper repository holds, per figure, a folder per panel with
`stats.csv` (every statistical analysis done for the panel), `data.csv` (the data plotted),
`panel_x.png`, `panel_x.pdf`, and `source_data/` -- or a `SOURCE_DATA_LOCATION.txt` saying where
the source is on disk and why it could not be included.

Figure builders call `panel(...)` once per panel and `manifest(...)` once per figure. The
repository is then laid out by parascape-paper/tools/build_repo.py from the manifests written to
results/paper_export/<figure>/.
"""
from __future__ import annotations

import json
import os
import shutil

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPORT = os.path.join(ROOT, "results", "paper_export")
#: A source file larger than this is not copied; its location is written down instead. GitHub
#: refuses files over 100 MB and a paper repository should stay cloneable.
MAX_SOURCE_MB = 10.0

_REGISTRY: dict = {}


def _stats_frame(stats) -> pd.DataFrame:
    if stats is None:
        return pd.DataFrame(columns=["statistic", "value", "note"])
    if isinstance(stats, pd.DataFrame):
        return stats
    rows = []
    for k, v in dict(stats).items():
        if isinstance(v, tuple):
            rows.append({"statistic": k, "value": v[0], "note": v[1]})
        else:
            rows.append({"statistic": k, "value": v, "note": ""})
    return pd.DataFrame(rows)


def panel(figure: str, letter: str, fig=None, png: str | None = None, data=None, stats=None,
          sources=(), scale: float = 1.0, what: str = "") -> str:
    """Export one panel; returns the PNG path, for composing the figure.

    `fig` is a matplotlib figure drawn here, saved as PNG (600 dpi) and PDF (vector, with dense
    point clouds left rasterised by the drawing code). `png` instead names a panel drawn
    elsewhere (the bradyzoite paper), whose sibling PDF is copied when it exists. `scale` carries
    the enlargement compose() applies to a panel drawn at another paper type size.
    """
    out = os.path.join(EXPORT, figure, f"panel_{letter}")
    os.makedirs(out, exist_ok=True)
    stem = os.path.join(out, f"panel_{letter}")
    if fig is not None:
        import matplotlib.pyplot as plt
        fig.savefig(stem + ".png", dpi=600, bbox_inches="tight", pad_inches=0.02)
        fig.savefig(stem + ".pdf", bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
    elif png is not None:
        shutil.copyfile(png, stem + ".png")
        pdf = os.path.splitext(png)[0] + ".pdf"
        if os.path.exists(pdf):
            shutil.copyfile(pdf, stem + ".pdf")
    else:
        raise ValueError("panel() needs a figure or a png")

    if data is not None:
        (data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)).to_csv(
            os.path.join(out, "data.csv"), index=isinstance(data, pd.DataFrame)
            and not isinstance(data.index, pd.RangeIndex))
    _stats_frame(stats).to_csv(os.path.join(out, "stats.csv"), index=False)

    src_dir = os.path.join(out, "source_data")
    if os.path.isdir(src_dir):
        shutil.rmtree(src_dir)
    os.makedirs(src_dir)
    too_big = []
    for s in sources:
        if isinstance(s, tuple):
            # (table, columns): the columns this panel used, extracted from a large shared table
            # (the 20 MB gene table is read by most panels; one copy per panel would not do).
            path, cols = s
            frame = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(
                path, sep="\t")
            keep = [c for c in (["gene_id"] + list(cols)) if c in frame.columns]
            base = os.path.splitext(os.path.basename(path))[0]
            frame[list(dict.fromkeys(keep))].to_csv(
                os.path.join(src_dir, f"{base}__columns_used.csv"), index=False)
            too_big.append((path, os.path.getsize(path) / 1e6))
            continue
        if not s or not os.path.exists(s):
            continue
        size_mb = os.path.getsize(s) / 1e6 if os.path.isfile(s) else sum(
            os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(s) for f in fs) / 1e6
        if os.path.isfile(s) and size_mb <= MAX_SOURCE_MB:
            shutil.copyfile(s, os.path.join(src_dir, os.path.basename(s)))
        else:
            too_big.append((s, size_mb))
    if too_big:
        with open(os.path.join(src_dir, "SOURCE_DATA_LOCATION.txt"), "w") as fh:
            fh.write("Source data not copied whole into the repository. Where a file named "
                     "<table>__columns_used.csv sits beside this note, it holds every column of "
                     "that table this panel used.\n\n")
            for s, mb in too_big:
                why = (f"{mb:,.1f} MB, above the {MAX_SOURCE_MB:.0f} MB limit for this repository"
                       if os.path.isfile(s) else f"a directory of {mb:,.1f} MB")
                fh.write(f"{s}\n  on disk at the path above (workstation 'carruthers');\n"
                         f"  not included: {why}.\n\n")
    _REGISTRY.setdefault(figure, []).append({"letter": letter, "png": stem + ".png",
                                             "pdf": stem + ".pdf" if os.path.exists(stem + ".pdf")
                                             else None, "scale": scale, "what": what})
    return stem + ".png"


def manifest(figure: str, kind: str, title: str, composed_png: str, layout: str | None = None):
    """Record a finished figure: its kind (main / supplementary / development), its one-line
    claim, its panels, and the composed figure with its panel positions."""
    out = os.path.join(EXPORT, figure)
    os.makedirs(out, exist_ok=True)
    m = {"figure": figure, "kind": kind, "title": title, "panels": _REGISTRY.get(figure, []),
         "composed_png": composed_png, "layout": layout}
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump(m, fh, indent=2)
    return m


def compose(figure: str, kind: str, title: str, rows=None, per_row: int = 4) -> str:
    """Compose the exported panels of `figure` (A4, true size) and write its manifest."""
    import assemble_paper_figures as A
    A.PER_ROW = per_row
    panels = [{"png": p["png"], "letter": p["letter"], "scale": p["scale"]}
              for p in _REGISTRY.get(figure, [])]
    out = os.path.join(EXPORT, figure)
    made = A.compose(figure, panels, out, rows=rows)          # writes figure<figure>.*
    stem = os.path.splitext(made)[0]
    final = os.path.join(out, figure)
    for ext in (".png", ".pdf", "_layout.json"):
        os.replace(stem + ext, final + ext)
    manifest(figure, kind, title, final + ".png", layout=final + "_layout.json")
    return final + ".png"
