"""Figure 1, one claim: ParaScape gathers the published Toxoplasma data into one gene space, and
that space shows what has never been measured: what a gene does to the host, and what happens
when it is lost, are measured for a few percent of genes.

Einar, 2026-11-19: 1e "phenotype" means nothing (say what was measured), 1f (categories per gene,
one violin) is pointless, the placeholder line in 1d says nothing. Every panel here is real data:

  A  what ParaScape does                                                 (schematic, reused)
  B  datasets by organism and level of biology                           (what was gathered)
  C  datasets by year of publication: all vs those measuring host or mutant phenotype
  D  coverage map: kinds of knowledge x genes                          (where the space is empty)
  E  the thin kinds, opened up: genes measured by each specific assay    (what "phenotype" is)
  F  the same coverage for secreted proteins vs all genes. Result (2026-11-19): secreted proteins
     are better measured in nearly every kind -- host response 64.0 % against 5.3 % genome-wide,
     because the effector screens target them -- but mutant phenotype (4.6 %) and enzyme
     activity (14.5 %) stay thin even for them. Localisation is left out: "secreted" is itself
     a localisation call.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import datasets as DS                    # noqa: E402
from parascape import figstyle as FS                    # noqa: E402
from parascape import information_categories as IC      # noqa: E402
import build_corpus_inventory as CI                     # noqa: E402

OUT = os.path.join(ROOT, "results", "figure01_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
NOT_KNOWLEDGE = CI.NOT_KNOWLEDGE

#: What each kind of knowledge is, in words a reader does not have to look up.
CATEGORY_NAMES = {
    "transcript": "mRNA level", "chromatin": "chromatin state", "translation": "translation",
    "structure_computed": "predicted structure", "fitness": "fitness (CRISPR)",
    "sequence": "sequence variation", "protein_abundance": "protein level",
    "modification": "protein modification", "place": "localisation",
    "structure_measured": "measured structure", "interaction": "interaction partners",
    "function": "enzyme activity (EC)", "phenotype": "mutant phenotype",
    "host_response": "host response",
}
THIN = ("function", "phenotype", "host_response")
LEVEL_NAMES = {"post_translation": "protein", "transcription": "RNA", "translation": "translation",
               "DNA": "DNA", "reference": "reference"}

#: The specific assays inside the thin categories, by column prefix. Einar: "phenotype" alone
#: says nothing -- the reader needs to know it means an imaging screen, a drug, a reporter.
ASSAYS = [
    ("hosttx_", "host mRNA after effector loss"),
    ("escrt_recruitment", "host ESCRT recruitment"),
    ("iedb_epitope", "immune epitopes"), ("n_bcell_epitopes", "immune epitopes"),
    ("screen_", "imaging screen"),
    ("drug_", "drug sensitivity of mutants"),
    ("resistance_", "drug-resistance alleles"),
    ("diff_reporter", "differentiation reporter"),
    ("ec_number", "EC number"), ("has_ec", "EC number"),
]
SECRETED = ("dense granules", "rhoptries 1", "rhoptries 2", "micronemes")


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def knowledge_columns(nodes):
    """Category -> columns, the same rule as the corpus inventory: a column with a value for every
    gene is a derivation or a zero-filled count, not a measurement, and is set aside."""
    derived = set(getattr(IC, "DERIVED_FROM", {}))
    out = {}
    for c in nodes.columns:
        if c in derived or not nodes[c].isna().any():
            continue
        for cat in IC.categories_of(c):
            if cat not in NOT_KNOWLEDGE:
                out.setdefault(cat, []).append(c)
    return out


#: What ParaScape answers, each tied to the figure of this paper that shows it answered.
QUESTIONS = [
    ("Can these studies be compared?", "removes cell cycle and laboratory first", "Figs 2, 4"),
    ("Which genes work together?", "co-expression modules", "Fig. 3"),
    ("Which cyst state does a drug reach?", "single-cell bradyzoite map", "Fig. 5"),
    ("Where is this protein?", "placed from RNA and from the graph", "Figs 6, 7"),
    ("What does it do to the host?", "host response, strain variation, co-folds", "Figs 8, 10"),
    ("Does it hold in Plasmodium?", "transfer between species", "Fig. 9"),
    ("What is unknown, and what next?", "filled-in map, missing measurements", "Figs 11, 12"),
]


def panel_a(nodes, cols):
    """What ParaScape does, drawn from its own numbers: what it gathered, the one gene space
    those land in (the real coverage of 8,140 genes, holes and all), and the questions it
    answers, each pointed at the figure that answers it.

    Einar, 2026-11-19: the Gene/List/Map/Ask/Export flow chart was ugly and said nothing about
    what ParaScape can do. Every number below is read from the registry or the gene table.
    """
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    FS.use()
    W, H = 180.0, 62.0
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    bold = dict(fontsize=FS.FONT_PT, fontweight="bold", color=FS.INK)
    body = dict(fontsize=FS.FONT_PT, color=FS.INK)

    # 1. gathered
    t = CI.by_organism_and_level()
    lev = t.groupby("level")["datasets"].sum().sort_values(ascending=False)
    org = t.groupby("organism")["datasets"].sum().sort_values(ascending=False)
    ax.text(2, H - 3, f"{int(lev.sum())} published datasets", va="top", **bold)
    y = H - 10
    top = float(lev.max())
    for level, n in lev.items():
        ax.add_patch(FancyBboxPatch((2, y - 3.2), 30 * n / top, 3.6,
                                    boxstyle="round,pad=0,rounding_size=0.8", lw=0,
                                    fc=FS.CONTROL, alpha=FS.ALPHA))
        ax.text(2 + 30 * n / top + 1.2, y - 1.4, f"{LEVEL_NAMES.get(level, level)}  {int(n)}",
                va="center", **body)
        y -= 5.2
    y -= 1.5
    ax.text(2, y, "from", va="top", color=FS.DARK_GREY, fontsize=FS.FONT_PT)
    y -= 4.2
    for name, n in org.items():
        ax.text(2, y, f"{name}  {int(n)}", va="top", fontstyle="italic"
                if name[0].isupper() else "normal", **body)
        y -= 3.9

    # 2. one gene space: the real coverage, 14 kinds x every gene, holes and all
    x0, x1 = 58.0, 104.0
    ax.add_patch(FancyArrowPatch((47, H / 2), (x0 - 2, H / 2), arrowstyle="-|>",
                                 mutation_scale=7, lw=0.8, color=FS.INK))
    share = {c: nodes[v].notna().any(axis=1).mean() for c, v in cols.items()}
    order = sorted(share, key=share.get, reverse=True)
    m = np.column_stack([nodes[cols[c]].notna().any(axis=1).to_numpy() for c in order]).T
    m = m[:, np.argsort(-m.sum(axis=0), kind="stable")]
    rgb = np.ones(m.shape + (3,))
    thin = np.array([c in THIN for c in order])
    for k, colour in ((True, FS.ACCENT), (False, FS.DARK_GREY)):
        rgb[m & (thin[:, None] == k)] = matplotlib.colors.to_rgb(colour)
    ax.text(x0, H - 3, "One gene space", va="top", **bold)
    ax.text(x0, H - 7.2, f"{len(nodes):,} genes \u00d7 {len(order)} kinds of knowledge", va="top",
            **body)
    ax.imshow(rgb, extent=(x0, x1, 14, H - 13), aspect="auto", interpolation="nearest")
    ax.add_patch(plt.Rectangle((x0, 14), x1 - x0, H - 27, fill=False, lw=0.6, ec=FS.INK))
    ax.text(x0, 11.5, "white: never measured", va="top", color=FS.DARK_GREY,
            fontsize=FS.FONT_PT)
    ax.text(x0, 7.5, "blue: host, mutant phenotype, EC", va="top",
            color=FS.ACCENT, fontsize=FS.FONT_PT)
    ax.text(x0, 3.6, "bring a gene, a list or your own data", va="top",
            color=FS.INK, fontsize=FS.FONT_PT)

    # 3. the questions it answers
    q0 = 116.0
    ax.add_patch(FancyArrowPatch((x1 + 2, H / 2), (q0 - 3, H / 2), arrowstyle="-|>",
                                 mutation_scale=7, lw=0.8, color=FS.INK))
    ax.text(q0, H - 3, "Questions it answers", va="top", **bold)
    y = H - 10
    for q, how, figs in QUESTIONS:
        ax.plot([q0 + 0.8], [y - 1.6], "o", color=FS.ACCENT, ms=3)
        ax.text(q0 + 2.6, y, q, va="top", **body)
        ax.text(q0 + 2.6, y - 3.4, f"{how}  ({figs})", va="top", color=FS.DARK_GREY,
                fontsize=FS.FONT_PT)
        y -= 7.4
    return fig


def panel_b():
    t = CI.by_organism_and_level()
    m = t.pivot_table(index="organism", columns="level", values="datasets", aggfunc="sum",
                      fill_value=0)
    m = m.loc[m.sum(axis=1).sort_values(ascending=False).index,
              m.sum(axis=0).sort_values(ascending=False).index]
    m.columns = [LEVEL_NAMES.get(c, c) for c in m.columns]
    FS.use()
    nr, nc = m.shape
    cell = FS.SLOT_IN
    w, h = 1.0 + nc * cell + 0.9, 0.3 + nr * cell + 0.75
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([1.0 / w, 0.75 / h, nc * cell / w, nr * cell / h])
    im = ax.imshow(m.to_numpy(), cmap="Blues", vmin=0, aspect="equal")
    ax.set_xticks(range(nc))
    ax.set_xticklabels(m.columns, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(nr)); ax.set_yticklabels(m.index)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.05)
    cb.set_label("Datasets"); cb.outline.set_linewidth(0.6)
    return fig


def dataset_years():
    rows = []
    for d in DS.REGISTRY:
        y = re.findall(r"\b(19[5-9]\d|20[0-2]\d)\b", str(d.citation or ""))
        cats = set()
        for c in d.columns:
            try:
                cats |= set(IC.categories_of(c, d.organism))
            except Exception:
                pass
        rows.append((int(y[0]) if y else np.nan, bool(cats & {"phenotype", "host_response"})))
    return pd.DataFrame(rows, columns=["year", "thin"])


def panel_c(years):
    y = years.dropna(subset=["year"])
    span = np.arange(int(y["year"].min()), int(y["year"].max()) + 1)
    allc = [int((y["year"] <= s).sum()) for s in span]
    thin = [int(((y["year"] <= s) & y["thin"]).sum()) for s in span]
    FS.use()
    fig, ax = FS.square_figure()
    ax.step(span, allc, where="post", color=FS.CONTROL, lw=1.2, alpha=1.0,
            label="all datasets (ref.)")
    ax.step(span, thin, where="post", color=FS.ACCENT, lw=1.2, alpha=1.0,
            label="host response or\nmutant phenotype")
    ax.set_xlabel("Year of publication"); ax.set_ylabel("Datasets, cumulative")
    FS.pad_limits(ax, allc + [0], axis="y", frac=0.05)
    FS.legend(ax, loc="upper left")
    return fig


def panel_d(nodes, cols):
    """Coverage map: one row per kind of knowledge, one column per gene, filled where measured.

    Einar, 2026-11-19: wide graphs need another form. Fourteen bars in a row made a strip twice
    as wide as tall; the map is square, reads its row names flat, and shows per GENE what the
    old violin of categories-per-gene tried to (and what was pointless as a lone violin).
    """
    share = {c: float(nodes[v].notna().any(axis=1).mean()) * 100 for c, v in cols.items()}
    order = sorted(share, key=share.get, reverse=True)
    m = np.column_stack([nodes[cols[c]].notna().any(axis=1).to_numpy() for c in order]).T
    m = m[:, np.lexsort(m[::-1])[::-1]]                 # genes ordered by what is known
    m = m[:, np.argsort(-m.sum(axis=0), kind="stable")]
    FS.use()
    fig, ax = FS.square_figure()
    from matplotlib.colors import ListedColormap
    rgb = np.ones(m.shape + (3,))
    thin = np.array([c in THIN for c in order])
    for k, colour in ((True, FS.ACCENT), (False, FS.DARK_GREY)):
        rows = m & (thin[:, None] == k)
        rgb[rows] = matplotlib.colors.to_rgb(colour)
    ax.imshow(rgb, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([CATEGORY_NAMES.get(c, c) for c in order])
    ax.set_xticks([0, m.shape[1] - 1]); ax.set_xticklabels(["1", f"{m.shape[1]:,}"])
    ax.set_xlabel("Genes, most to least measured")
    ax.tick_params(axis="y", length=0)
    ax.plot([], [], "s", color=FS.DARK_GREY, ms=5, label="measured")
    ax.plot([], [], "s", color=FS.ACCENT, ms=5, label="measured, thin kind")
    ax.plot([], [], "s", color="white", mec=FS.INK, mew=0.4, ms=5, label="never measured")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig, share


def dot_rows(n, value_mm=46.0):
    """Rows of named dots: names read flat down the left axis, values along x.

    Einar, 2026-11-19: wide graphs need another form. A row of 8-14 bars with tilted 30-character
    names was twice as wide as tall and put the names at an angle; named rows are compact and
    keep every word on the axis.
    """
    fig, ax = FS.figure(n, vertical=False, value_mm=value_mm)
    ax.set_ylim(n - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", color=FS.BACKGROUND, lw=0.4)
    ax.set_axisbelow(True)
    return fig, ax


HOST_ASSAYS = {"host mRNA after effector loss", "host ESCRT recruitment", "immune epitopes"}


def panel_e(nodes, cols):
    """The thin kinds opened up: how many genes each specific assay has measured."""
    thin_cols = sorted({c for cat in THIN for c in cols.get(cat, [])})
    groups = {}
    for c in thin_cols:
        name = next((lab for pre, lab in ASSAYS if c.startswith(pre)), c)
        groups.setdefault(name, []).append(c)
    count = {k: int(nodes[v].notna().any(axis=1).sum()) for k, v in groups.items()}
    order = sorted(count, key=count.get, reverse=True)
    FS.use()
    fig, ax = dot_rows(len(order))
    y = np.arange(len(order))
    host = np.array([k in HOST_ASSAYS for k in order])
    vals = np.array([count[k] for k in order])
    ax.plot(vals[~host], y[~host], "o", color=FS.CONTROL, ms=4, alpha=FS.ALPHA,
            label="mutant phenotype or EC")
    ax.plot(vals[host], y[host], "o", color=FS.ACCENT, ms=4, alpha=FS.ALPHA,
            label="host response")
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.set_xscale("log"); ax.set_xlim(0.7, len(nodes) * 1.3)
    ax.set_xticks([1, 10, 100, 1000]); ax.set_xticklabels(["1", "10", "100", "1,000"])
    ax.minorticks_off()
    ax.axvline(len(nodes), color=FS.INK, lw=0.6, ls="--", label="all genes")
    ax.set_xlabel("Genes measured")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig, {k: 100 * v / len(nodes) for k, v in count.items()}


def panel_f(nodes, cols):
    """Per kind of knowledge, genes measured among all genes vs among secreted proteins."""
    sec = nodes["lopit_map"].astype(str).isin(SECRETED).to_numpy()
    rows = []
    for c, v in cols.items():
        if c == "place":       # "secreted" IS a localisation call: 100 % by definition
            continue
        m = nodes[v].notna().any(axis=1).to_numpy()
        rows.append((c, m.mean() * 100, m[sec].mean() * 100))
    d = pd.DataFrame(rows, columns=["cat", "all", "sec"]).sort_values("all", ascending=False)
    FS.use()
    fig, ax = dot_rows(len(d))
    y = np.arange(len(d))
    ax.hlines(y, d["all"], d["sec"], color=FS.CONTROL, lw=0.6)
    ax.plot(d["all"], y, "o", color=FS.CONTROL, ms=4, alpha=FS.ALPHA, label="all genes (ref.)")
    ax.plot(d["sec"], y, "o", color=FS.ACCENT, ms=4, alpha=FS.ALPHA, label="secreted proteins")
    ax.set_yticks(y); ax.set_yticklabels([CATEGORY_NAMES.get(c, c) for c in d["cat"]])
    ax.set_xlim(-4, 104); ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Genes measured (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig, d, int(sec.sum())


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    cols = knowledge_columns(nodes)
    years = dataset_years()
    figd, share = panel_d(nodes, cols)
    fige, assay = panel_e(nodes, cols)
    figf, sec, n_sec = panel_f(nodes, cols)
    panels = [("a", save(panel_a(nodes, cols), "a_schematic")), ("b", save(panel_b(), "b_datasets")),
              ("c", save(panel_c(years), "c_years")), ("d", save(figd, "d_coverage")),
              ("e", save(fige, "e_assays")), ("f", save(figf, "f_secreted"))]
    print(f"{len(nodes):,} genes; {len(DS.REGISTRY)} datasets, "
          f"{years['year'].notna().sum()} with a year in the citation")
    print("coverage %:", {k: round(v, 1) for k, v in share.items()})
    print("assays %:", {k: round(v, 1) for k, v in assay.items()})
    print(f"secreted proteins: {n_sec}")
    print(sec.round(1).to_string(index=False))
    import assemble_paper_figures as A
    A.PER_ROW = 3
    out = A.compose("01", [{"png": p, "letter": l} for l, p in panels], OUT, rows=[1, 3, 2])
    print("figure 1 ->", out)


if __name__ == "__main__":
    main()
