"""Build all twelve ParaScape figures from `paper_figure_spec.FIGURES`.

For each panel in the spec, in this order of preference:
  1. a REAL panel already drawn by an analysis script (reused as its PNG);
  2. a REAL panel drawn here from a result table that has no drawing yet;
  3. PLACEHOLDER DATA in the panel's declared graph form, seeded, and marked as such in the
     panel's own legend -- Einar, 2026-11-19: "wherever data is missing make placeholder data so
     i can evaluate the figure structure".

Then every figure is assembled as a grid. Figure 5 comes from `build_figure05_bradyzoite.py`.
"""
from __future__ import annotations

import glob
import os
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from parascape import figstyle as FS        # noqa: E402
import paper_figure_spec as SPEC            # noqa: E402

OUT = os.path.join(ROOT, "results", "paper_figures_v3")
PANELS = os.path.join(OUT, "panels")
P = os.path.join(ROOT, "results")

#: Real panels already drawn, by (figure, letter). The flag says whether the PNG lacks a letter,
#: in which case the assembler adds it.
REUSE = {
    (1, "a"): ("figure1_schematic_*/panel_schematic.png", True),
    (1, "e"): ("corpus_inventory_*/panel_coverage.png", True),
    (1, "f"): ("paper_panels_*/figure12b_categories_per_gene.png", False),
    (2, "b"): ("paper_panels_*/figure02_comparability.png", True),
    (3, "a"): ("paper_panels_*/figure03_modules.png", True),
    (3, "b"): ("paper_panels_*/figure03b_vs_baseline.png", False),
    (3, "c"): ("paper_panels_*/figure03c_sizes.png", False),
    (3, "d"): ("paper_panels_*/figure03d_detection.png", False),
    (4, "a"): ("paper_panels_*/figure04_what_structures.png", False),
    (6, "a"): ("coexpression_components_localization_*/panel_which_evidence_places_a_protein.png",
               False),
    (6, "b"): ("paper_panels_*/figure06b_protein_abundance.png", False),
    (6, "c"): ("paper_panels_*/figure06c_explains_which.png", False),
    (7, "b"): ("paper_panels_*/figure07_network_layers.png", False),
    (8, "a"): ("paper_panels_*/figure08_strain_variation.png", False),
    (8, "b"): ("paper_panels_*/figure08b_foldability.png", False),
    (9, "a"): ("paper_panels_*/figure09_species_transfer.png", True),
    (10, "a"): ("paper_panels_*/figure10_host_shared.png", False),
    (10, "b"): ("paper_panels_*/figure10b_single_cell.png", False),
    (11, "a"): ("paper_panels_*/figure11_filling_the_map.png", False),
}


def newest(pattern: str):
    found = sorted(glob.glob(os.path.join(P, pattern)))
    return found[-1] if found else None


def save(fig, figno: int, letter: str) -> str:
    os.makedirs(PANELS, exist_ok=True)
    path = os.path.join(PANELS, f"fig{figno:02d}{letter}.png")
    fig.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return path


# ------------------------------------------------------------------ real panels drawn here
def bars_from(values: dict, xlabel: str, ylabel: str, letter: str, accent=()):
    names = list(values)
    FS.use()
    fig, ax = FS.figure(len(names), vertical=True)
    FS.bars(ax, [values[n] for n in names], names, highlight=set(accent))
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, len(names) - 0.4)
    FS.pad_limits(ax, list(values.values()) + [0])
    FS.letter(fig, letter)
    return fig


def real_1b():
    d = pd.read_csv(newest("corpus_inventory_*/datasets.tsv"), sep="\t")
    v = d.groupby("organism")["datasets"].sum().sort_values(ascending=False).to_dict()
    return bars_from(v, "Organism", "Datasets", "b", accent={next(iter(v))})


def real_1c():
    d = pd.read_csv(newest("corpus_inventory_*/datasets.tsv"), sep="\t")
    v = d.groupby("level")["datasets"].sum().sort_values(ascending=False)
    v.index = [i.replace("_", " ") for i in v.index]
    return bars_from(v.to_dict(), "Level of biology", "Datasets", "c")


def real_7(col, ylabel, letter):
    import draw_remaining_panels as D
    L = pd.read_csv(newest("network_layers_*/layers.tsv"), sep="\t")
    L = L[L["enough"].astype(bool)].sort_values(col, ascending=False)
    keys = [s.replace(" AT ONCE", "").lower() for s in L["layer"]]
    # Display names, not record keys: 'xlms' and 'struct' are keys.
    names = [D.LAYER_NAMES.get(k, k.replace("_", " ")) for k in keys]
    vals = L[col].to_numpy()
    if col == "accuracy":
        vals = vals * 100          # per cent, to share an axis unit with panel B
    v = dict(zip(names, vals))
    return bars_from(v, "Edge layer", ylabel, letter,
                     accent={n for n in names if n.startswith("all layers")})


def real_10c():
    d = pd.read_csv(newest("host_response_shared_*/effectors_per_host_gene.tsv"), sep="\t")
    col = [c for c in d.columns if "effector" in c.lower() and d[c].dtype != object]
    s = d[col[0]] if col else d.iloc[:, 1]
    v = s.value_counts().sort_index()
    v = {str(k): int(n) for k, n in v.items() if k <= 12}
    return bars_from(v, "Effectors moving the gene", "Host genes", "c")


def real_10f():
    files = sorted(glob.glob(os.path.join(P, "dual_perturb_single_cell_*/*_targets.tsv")))
    t = max((pd.read_csv(f, sep="\t") for f in files), key=len)
    FS.use()
    fig, ax = FS.figure(2, vertical=True)
    FS.violin(ax, [t["cells"].to_numpy()], ["Knockouts"], control=["Knockouts"], vertical=True)
    ax.set_xticks([0]); ax.set_xticklabels(["Knockouts"])
    ax.set_ylabel("Cells per knockout"); ax.set_xlim(-0.75, 1.25)
    ax.tick_params(axis="x", length=0)
    FS.letter(fig, "f")
    return fig


def real_11b():
    a = pd.read_csv(newest("filling_the_map_*/by_target_and_strategy.tsv"), sep="\t")
    order = a.groupby("target")["lift_over_baseline"].max().sort_values(ascending=False).index
    groups = [a.loc[a["target"] == t, "lift_over_baseline"].to_numpy() for t in order]
    FS.use()
    fig, ax = FS.figure(len(groups), vertical=True)
    bp = ax.boxplot(groups, positions=range(len(groups)), widths=FS.bar_width(),
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color=FS.INK, lw=0.9),
                    whiskerprops=dict(color=FS.INK, lw=0.6), capprops=dict(color=FS.INK, lw=0.6))
    for b in bp["boxes"]:
        b.set_facecolor(FS.CONTROL); b.set_alpha(FS.ALPHA); b.set_linewidth(0)
    ax.axhline(0, color=FS.INK, lw=0.6)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([t.replace("_", " ")[:16] for t in order], rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel("Accuracy above baseline"); ax.tick_params(axis="x", length=0)
    FS.letter(fig, "b")
    return fig


def real_11c():
    b = pd.read_csv(newest("filling_the_map_*/best_per_target.tsv"), sep="\t")
    FS.use()
    fig, ax = FS.square_figure()
    ax.plot(b["answer_rate"], b["right_when_answered"], "o", color=FS.ACCENT, ms=4,
            alpha=FS.ALPHA, lw=0)
    ax.set_xlabel("Answer rate"); ax.set_ylabel("Accuracy when answered")
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
    FS.letter(fig, "c")
    return fig


REAL = {(1, "b"): real_1b, (1, "c"): real_1c,
        (7, "c"): lambda: real_7("accuracy", "Accuracy (%)", "c"),
        (7, "d"): lambda: real_7("answered", "Genes answered", "d"),
        (10, "c"): real_10c, (10, "f"): real_10f,
        (11, "b"): real_11b, (11, "c"): real_11c}


# ------------------------------------------------------------------ placeholder data
def placeholder(panel, figno: int):
    """Seeded placeholder data in the panel's REAL graph form, marked in its own legend."""
    rng = np.random.default_rng(figno * 100 + ord(panel.letter))
    FS.use()
    form = panel.form
    mark = dict(color="none", label="Placeholder data")
    if form in ("bar", "stacked bar"):
        n = 6
        fig, ax = FS.figure(n, vertical=True)
        vals = np.sort(rng.uniform(0.1, 1.0, n))[::-1]
        if form == "bar":
            FS.bars(ax, list(vals), [f"cond. {i + 1}" for i in range(n)], highlight={"cond. 1"})
        else:
            low = vals * rng.uniform(0.2, 0.6, n)
            ax.bar(range(n), low, width=FS.bar_width(), color=FS.ACCENT, alpha=FS.ALPHA, lw=0)
            ax.bar(range(n), vals - low, bottom=low, width=FS.bar_width(), color=FS.CONTROL,
                   alpha=FS.ALPHA, lw=0)
            ax.set_xticks(range(n)); ax.set_xticklabels([f"cond. {i + 1}" for i in range(n)],
                                                        rotation=30, ha="right",
                                                        rotation_mode="anchor")
        ax.set_xlim(-0.6, n - 0.4); ax.set_ylim(0, 1.1)
    elif form in ("box", "violin"):
        n = 4
        fig, ax = FS.figure(n, vertical=True)
        groups = [rng.normal(0.5 + 0.1 * i, 0.15, 60) for i in range(n)]
        names = [f"group {i + 1}" for i in range(n)]
        FS.violin(ax, groups, names, control=names[1:], vertical=True, draw_points=False)
        ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=30, ha="right",
                                                    rotation_mode="anchor")
        ax.set_xlim(-0.6, n - 0.4); ax.tick_params(axis="x", length=0)
    elif form == "scatter":
        fig, ax = FS.square_figure()
        x = rng.uniform(0, 1, 80); y = 0.7 * x + rng.normal(0, 0.12, 80)
        ax.plot(x, y, "o", color=FS.CONTROL, ms=3, alpha=FS.ALPHA, lw=0)
        ax.plot([0, 1], [0, 1], "-", color=FS.BACKGROUND, lw=0.7)
    elif form == "line":
        fig, ax = FS.square_figure()
        x = np.arange(8); y = np.cumsum(rng.uniform(0.2, 1.0, 8))
        ax.plot(x, y, "-o", color=FS.ACCENT, ms=3, lw=1.0, alpha=FS.ALPHA)
        ax.set_ylim(0, y.max() * 1.1)
    elif form == "heatmap":
        fig = plt.figure(figsize=(2.4, 2.2))
        ax = fig.add_axes([0.18, 0.18, 0.62, 0.70])
        m = rng.normal(0, 1, (6, 6))
        im = ax.imshow(m, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="equal")
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.outline.set_linewidth(0.6)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])
    elif form == "dot plot":
        fig, ax = FS.square_figure()
        xs, ys = np.meshgrid(range(5), range(5))
        ax.scatter(xs.ravel(), ys.ravel(), s=rng.uniform(5, 60, 25), color=FS.ACCENT,
                   alpha=FS.ALPHA, lw=0)
        ax.set_xlim(-0.6, 4.6); ax.set_ylim(-0.6, 4.6)
    else:   # schematic: a node-link sketch, so the slot reads as a schematic and not as nothing
        fig, ax = FS.square_figure()
        ax.axis("off")
        pts = rng.uniform(0.1, 0.9, (8, 2))
        for i in range(8):
            for j in range(i + 1, 8):
                if rng.random() < 0.3:
                    ax.plot(*zip(pts[i], pts[j]), color=FS.BACKGROUND, lw=0.8, zorder=1)
        ax.scatter(pts[:, 0], pts[:, 1], s=120, color="#E4E4E4", edgecolor="#9A9A9A",
                   lw=0.6, zorder=2)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.plot([], [], " ", **mark)
    # A heatmap has no empty corner, so its key goes outside; elsewhere it stays with the panel.
    if form == "heatmap":
        FS.legend(ax, loc="upper left", bbox_to_anchor=(0.0, 1.18))
    else:
        FS.legend(ax, loc="upper left")
    if form not in ("heatmap", "schematic"):
        ax.set_ylabel("Value")
    FS.letter(fig, panel.letter)
    return fig


# ------------------------------------------------------------------ build
def build_panel(figno: int, panel) -> tuple:
    key = (figno, panel.letter)
    if key in REUSE:
        src = newest(REUSE[key][0])
        if src:
            dst = os.path.join(PANELS, f"fig{figno:02d}{panel.letter}.png")
            os.makedirs(PANELS, exist_ok=True)
            shutil.copyfile(src, dst)
            return dst, "real", REUSE[key][1]
    if key in REAL:
        try:
            return save(REAL[key](), figno, panel.letter), "real", False
        except Exception as e:
            print(f"    real {figno}{panel.letter} failed ({type(e).__name__}: {e}); placeholder")
    return save(placeholder(panel, figno), figno, panel.letter), "placeholder", False


def main() -> None:
    import assemble_paper_figures as A
    A.PER_ROW = 3
    os.makedirs(OUT, exist_ok=True)
    totals = {"real": 0, "placeholder": 0}
    for f in SPEC.FIGURES:
        if f.number == 5:
            pngs = sorted(glob.glob(os.path.join(P, "figure05_bradyzoite", "*.png")))
            panels = [{"png": p, "letter": os.path.basename(p)[0]} for p in pngs]
            totals["real"] += len(pngs)
        else:
            panels = []
            for panel in f.panels:
                path, kind, _needs = build_panel(f.number, panel)
                totals[kind] += 1
                panels.append({"png": path, "letter": panel.letter})
        A.PER_ROW = 4 if len(panels) > 6 else 3
        out = A.compose(f"{f.number:02d}", panels, OUT)
        print(f"  figure {f.number:2d}: {len(panels):2d} panels -> {os.path.basename(out)}")
    print(f"\n{totals['real']} real panels, {totals['placeholder']} placeholder")


if __name__ == "__main__":
    main()
