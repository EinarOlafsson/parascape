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
    ("Can these studies be compared?", "removes cell cycle and laboratory first", "Figs 2, 3"),
    ("Which genes work together?", "co-expression modules", "Fig. 4"),
    ("Where is this protein?", "placed from RNA and from the graph", "Figs 5, 6"),
    ("Does it hold in Plasmodium?", "transfer between species", "Fig. 7"),
    ("What can be filled in?", "validated prediction", "Fig. 8"),
    ("What does it do to the host?", "strain variation, host response", "Figs 9, 10"),
    ("Which cyst state does a drug reach?", "single-cell bradyzoite map", "Fig. 11"),
    ("What should be measured next?", "missing measurements", "Fig. 12"),
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
    W, H = 180.0, 66.0
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    bold = dict(fontsize=FS.FONT_PT, fontweight="bold", color=FS.INK)
    body = dict(fontsize=FS.FONT_PT, color=FS.INK)

    # 1. gathered
    t = CI.by_organism_and_level()
    lev = t.groupby("level")["datasets"].sum().sort_values(ascending=False)
    org = t.groupby("organism")["datasets"].sum().sort_values(ascending=False)
    ax.text(2, H - 3, f"{int(lev.sum())} published datasets", va="top", **bold)
    # Levels as labelled tiles, not a miniature bar chart (Einar 2026-11-20: no graphs inside
    # the schematic).
    y = H - 10
    for level, n in lev.items():
        ax.add_patch(FancyBboxPatch((2, y - 3.4), 36, 3.8,
                                    boxstyle="round,pad=0,rounding_size=0.8", lw=0.5,
                                    ec=FS.DARK_GREY, fc=FS.BACKGROUND, alpha=0.5))
        ax.text(3.2, y - 1.5, LEVEL_NAMES.get(level, level), va="center", **body)
        ax.text(37, y - 1.5, f"{int(n)}", va="center", ha="right", **body)
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
    # The gene space as an icon: genes as rows, each with a slot per kind of knowledge, filled
    # where measured. Illustrative, not data (the real coverage is panel D).
    share = {c: nodes[v].notna().any(axis=1).mean() for c, v in cols.items()}
    order = sorted(share, key=share.get, reverse=True)
    ax.text(x0, H - 3, "One gene space", va="top", **bold)
    ax.text(x0, H - 7.2, f"{len(nodes):,} genes \u00d7 {len(order)} kinds of knowledge",
            va="top", **body)
    rng = np.random.default_rng(3)
    nk, ng = len(order), 7
    cw, ch = (x1 - x0) / nk, 3.4
    for gi in range(ng):
        yy = H - 14 - gi * (ch + 0.9)
        for ki, c in enumerate(order):
            # thin kinds drawn as measured in two of the seven icon genes so the key reads
            filled = (gi in (1, 4)) if c in THIN else rng.random() < share[c]
            fc = (FS.ACCENT if c in THIN else FS.DARK_GREY) if filled else "white"
            ax.add_patch(plt.Rectangle((x0 + ki * cw + 0.25, yy - ch), cw - 0.5, ch, lw=0.4,
                                       ec=FS.DARK_GREY, fc=fc, alpha=0.8 if filled else 1))
    yb = H - 14 - ng * (ch + 0.9) - 1.5
    for k, (fc, lab) in enumerate(((FS.DARK_GREY, "measured"), (FS.ACCENT, "host, phenotype"),
                                  ("white", "never"))):
        yy = yb - k * 3.6
        ax.add_patch(plt.Rectangle((x0, yy - 2.6), 2.4, 2.4, lw=0.4, ec=FS.DARK_GREY, fc=fc))
        ax.text(x0 + 3.4, yy - 1.4, lab, va="center", color=FS.DARK_GREY, fontsize=FS.FONT_PT)
    ax.text(x1, yb - 1.4, "bring a gene,\na list or\nyour own data", va="top", ha="right", **body)

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
        y -= 7.0
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


SHORT_KIND = {"mRNA level": "mRNA", "chromatin state": "chromatin", "translation": "translation",
              "predicted structure": "pred. struct.", "fitness (CRISPR)": "fitness",
              "sequence variation": "variation", "protein level": "protein",
              "protein modification": "PTM", "localisation": "location",
              "measured structure": "meas. struct.", "interaction partners": "interaction",
              "enzyme activity (EC)": "enzyme (EC)", "mutant phenotype": "phenotype",
              "host response": "host response"}


def kind_name(c):
    return SHORT_KIND.get(CATEGORY_NAMES.get(c, c), CATEGORY_NAMES.get(c, c))


def panel_d(nodes, cols):
    """Genes measured per kind of knowledge: bars, the thin kinds in the accent."""
    share = {c: float(nodes[v].notna().any(axis=1).mean()) * 100 for c, v in cols.items()}
    order = sorted(share, key=share.get, reverse=True)
    FS.use()
    fig, ax = FS.ref_bars([share[c] for c in order], [kind_name(c) for c in order],
                          highlight=[kind_name(c) for c in THIN], ylabel="Genes measured (%)")
    ax.set_ylim(0, 100)
    return fig, share


HOST_ASSAYS = {"host mRNA after effector loss", "host ESCRT recruitment", "immune epitopes"}
ASSAY_SHORT = {"host mRNA after effector loss": "host mRNA", "host ESCRT recruitment":
               "host ESCRT", "immune epitopes": "epitopes", "imaging screen": "imaging",
               "drug sensitivity of mutants": "drug sens.", "drug-resistance alleles":
               "resistance", "differentiation reporter": "diff. reporter", "EC number": "EC number"}


def panel_e(nodes, cols):
    """The thin kinds opened up: genes measured by each assay."""
    thin_cols = sorted({c for cat in THIN for c in cols.get(cat, [])})
    groups = {}
    for c in thin_cols:
        name = next((lab for pre, lab in ASSAYS if c.startswith(pre)), c)
        groups.setdefault(name, []).append(c)
    count = {k: int(nodes[v].notna().any(axis=1).sum()) for k, v in groups.items()}
    order = sorted(count, key=count.get, reverse=True)
    FS.use()
    fig, ax = FS.ref_bars([max(count[k], 0.8) for k in order],
                          [ASSAY_SHORT.get(k, k) for k in order],
                          highlight=[ASSAY_SHORT[k] for k in HOST_ASSAYS],
                          ylabel="Genes measured", log=True)
    ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["1", "10", "100", "1,000"])
    ax.minorticks_off(); ax.set_ylim(0.7, len(nodes))
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=FS.ACCENT, alpha=FS.ALPHA, label="host response"),
                       Patch(color=FS.CONTROL, alpha=FS.ALPHA, label="phenotype or enzyme")],
              frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2,
              fontsize=FS.NOTE_PT, handlelength=1.0, borderaxespad=0.0)
    return fig, {k: 100 * v / len(nodes) for k, v in count.items()}


def panel_f(nodes, cols):
    """Per kind, genes measured among secreted proteins (bars) and among all genes (tick)."""
    sec = nodes["lopit_map"].astype(str).isin(SECRETED).to_numpy()
    rows = []
    for c, v in cols.items():
        if c == "place":       # "secreted" IS a localisation call: 100 % by definition
            continue
        m = nodes[v].notna().any(axis=1).to_numpy()
        rows.append((c, m.mean() * 100, m[sec].mean() * 100))
    d = pd.DataFrame(rows, columns=["cat", "all", "sec"]).sort_values("all", ascending=False)
    FS.use()
    names = [kind_name(c) for c in d["cat"]]
    fig, ax = FS.ref_bars(d["sec"].to_numpy(), names, refs=d["all"].to_numpy(),
                          highlight=names, ref_label="all genes",
                          ylabel="Secreted proteins measured (%)")
    ax.set_ylim(0, 105)
    return fig, d, int(sec.sum())


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    cols = knowledge_columns(nodes)
    years = dataset_years()
    figd, share = panel_d(nodes, cols)
    fige, assay = panel_e(nodes, cols)
    figf, sec, n_sec = panel_f(nodes, cols)
    import paper_panels as PP
    from scipy.stats import fisher_exact
    F = "figure_1"
    NODES = os.path.join(CACHE, "nodes.parquet")
    REG = os.path.join(ROOT, "parascape", "datasets.py")
    inv = CI.by_organism_and_level()
    lev = inv.groupby("level")["datasets"].sum(); org = inv.groupby("organism")["datasets"].sum()
    PP.panel(F, "a", fig=panel_a(nodes, cols),
             data=pd.concat([lev.rename("datasets").rename_axis("group").reset_index().assign(by="level"),
                             org.rename("datasets").rename_axis("group").reset_index().assign(by="organism")]),
             stats={"datasets": int(lev.sum()), "genes": len(nodes),
                    "kinds of knowledge": len(cols)}, sources=[REG],
             what="what ParaScape gathers, the gene space, and the questions it answers")
    PP.panel(F, "b", fig=panel_b(), data=inv, stats={"datasets": int(inv.datasets.sum()),
             "organisms": inv.organism.nunique(), "levels": inv.level.nunique()},
             sources=[REG], what="datasets by organism and level of biology")
    PP.panel(F, "c", fig=panel_c(years), data=years, stats={
        "datasets with a citation year": int(years.year.notna().sum()),
        "datasets without": int(years.year.isna().sum()),
        "dated datasets measuring host response or mutant phenotype": int(
            (years.thin & years.year.notna()).sum())}, sources=[REG],
             what="datasets by year of publication")
    presence = pd.DataFrame({CATEGORY_NAMES.get(c, c): nodes[v].notna().any(axis=1)
                             for c, v in cols.items()})
    presence.insert(0, "gene_id", nodes["gene_id"].astype(str).to_numpy())
    PP.panel(F, "d", fig=figd, data=presence, stats=pd.DataFrame(
        [{"kind": CATEGORY_NAMES.get(c, c), "columns": len(v), "genes_measured_pct": share[c]}
         for c, v in cols.items()]), sources=[NODES],
             what="coverage map: kinds of knowledge x genes")
    PP.panel(F, "e", fig=fige, data=pd.Series(assay, name="genes_measured_pct")
             .rename_axis("assay").reset_index(),
             stats={"genes": len(nodes)}, sources=[NODES],
             what="the thin kinds opened up: genes measured per assay")
    secm = nodes["lopit_map"].astype(str).isin(SECRETED).to_numpy()
    rows = []
    for c, v in cols.items():
        if c == "place":
            continue
        m = nodes[v].notna().any(axis=1).to_numpy()
        t = [[int((m & secm).sum()), int((~m & secm).sum())],
             [int((m & ~secm).sum()), int((~m & ~secm).sum())]]
        odds, pv = fisher_exact(t)
        rows.append({"kind": CATEGORY_NAMES.get(c, c), "secreted_measured": t[0][0],
                     "secreted_not": t[0][1], "other_measured": t[1][0], "other_not": t[1][1],
                     "fisher_odds_ratio": odds, "fisher_p": pv})
    PP.panel(F, "f", fig=figf, data=sec, stats=pd.DataFrame(rows), sources=[NODES],
             what="coverage among secreted proteins vs all genes")
    print(f"{len(nodes):,} genes; {len(DS.REGISTRY)} datasets, "
          f"{years['year'].notna().sum()} with a year in the citation")
    print("coverage %:", {k: round(v, 1) for k, v in share.items()})
    print("assays %:", {k: round(v, 1) for k, v in assay.items()})
    print(f"secreted proteins: {n_sec}")
    print(sec.round(1).to_string(index=False))
    out = PP.compose(F, "main", "ParaScape gathers the published Toxoplasma data into one gene "
                     "space; what a gene does to the host, and what happens when it is lost, are "
                     "measured for a few percent of genes.", rows=[1, 3, 2])
    print("figure 1 ->", out)


if __name__ == "__main__":
    main()
