"""Development figures: the observations of the main figures at the grain of single units.

Einar, 2026-11-20: development figures show the same as the main and supplementary figures "on a
more granular level where i can see the data more granularily". Each is a small multiple, one
panel per unit, every point drawn:

  development_figure_2   per study: its control profile vs the median of the others, before
                         and after removing the cell-cycle axis
  development_figure_3   per module: the cell-cycle trajectory of every gene
  development_figure_6   per compartment: every gene trajectory
  development_figure_8   per class: N/S against pLDDT, every gene
  development_figure_10  per effector: the distribution of host changes
  development_figure_11  per target: every strategy, answer rate against accuracy when answered
The unit of each panel is named in its legend (text only on axes and legends).
"""
from __future__ import annotations

import os
import string
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import paper_panels as PP                               # noqa: E402

CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
SIDE = float(os.environ.get("DEV_SIDE_MM", "30"))   # smaller for 20+ panels, to fit A4
LETTERS = list(string.ascii_lowercase) + [f"a{c}" for c in string.ascii_lowercase]


def small():
    FS.use()
    return FS.square_figure(side_mm=SIDE)


def name_legend(ax, text):
    ax.plot([], [], " ", label=text)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))


def dev2():
    import build_control_comparability as CC
    import build_figure02_v2 as F2
    F = "development_figure_2"
    raw, _ = CC.build_profiles(F2.OMICS)
    frame = CC.shared_genes(raw)
    ph = CC.phase_table(F2.CACHE); groups = {p: list(g.index) for p, g in ph.groupby(ph)}
    after = CC.correct_for(frame, groups)
    srcs = [os.path.join(F2.VALID, "agreement_before.tsv"), os.path.join(F2.OMICS, "matrices")]
    for i, st in enumerate(frame.columns):
        others_b = frame.drop(columns=[st]).median(axis=1)
        others_a = after.drop(columns=[st]).median(axis=1)
        fig, ax = small()
        ax.plot(others_b, frame[st], "o", color=FS.CONTROL, ms=0.8, alpha=0.3, lw=0,
                rasterized=True, label="before (ref.)")
        ax.plot(others_a, after[st], "o", color=FS.ACCENT, ms=0.8, alpha=0.3, lw=0,
                rasterized=True, label="cell cycle removed")
        ax.set_xlabel("Others, median (log$_2$)"); ax.set_ylabel("Study (log$_2$)")
        name_legend(ax, st)
        from scipy.stats import spearmanr
        PP.panel(F, LETTERS[i], fig=fig, data=pd.DataFrame({"study_before": frame[st],
                 "others_before": others_b, "study_after": after[st], "others_after": others_a}),
                 stats={"Spearman before": spearmanr(frame[st], others_b, nan_policy="omit")
                        .statistic, "Spearman after": spearmanr(after[st], others_a,
                                                                nan_policy="omit").statistic},
                 sources=srcs, what=f"{st}: control profile vs the other studies")
    PP.compose(F, "development", "Figure 2 per study: each control profile against the others, "
               "before and after removing the cell-cycle axis.")


def dev3():
    import build_figure03_v2 as F3
    F = "development_figure_3"
    gm = pd.read_csv(os.path.join(F3.MODS, "gene_modules.tsv"), sep="\t", index_col=0)["module"]
    mods = pd.read_csv(os.path.join(F3.MODS, "modules.tsv"), sep="\t")
    mods = mods[mods.genes >= 20].sort_values("over_detection_limit", ascending=False)
    n = pd.read_parquet(os.path.join(CACHE, "nodes.parquet")); n = n.set_index(n.gene_id.astype(str))
    H = list(range(1, 13))
    prof = pd.DataFrame({h: n[[f"cellcycle19092_{h}h_r1", f"cellcycle19092_{h}h_r2"]].mean(axis=1)
                         for h in H}).dropna()
    z = prof.sub(prof.mean(axis=1), axis=0).div(prof.std(axis=1).replace(0, np.nan), axis=0).dropna()
    # Two pages, so each fits A4: 3a the named modules, 3b those that do not clear chance.
    count = {"development_figure_3a": 0, "development_figure_3b": 0}
    for r in mods.itertuples():
        F = "development_figure_3a" if str(r.named) == "True" else "development_figure_3b"
        i = count[F]; count[F] += 1
        g = z.loc[z.index.intersection(gm.index[gm == r.module])]
        fig, ax = small()
        named = str(r.named) == "True"
        col = FS.ACCENT if named else FS.CONTROL
        for _, row in g.head(150).iterrows():
            ax.plot(H, row.to_numpy(), color=col, lw=0.3, alpha=0.15)
        ax.plot(H, g.mean(), color=FS.INK, lw=1.0, label="mean")
        ax.set_xticks([3, 6, 9, 12]); ax.set_xlabel("Time (h)"); ax.set_ylabel("Expression (z)")
        name_legend(ax, f"M{r.module} {r.compartment if named else '(not named)'}")
        PP.panel(F, LETTERS[i], fig=fig, data=g.rename_axis("gene_id"),
                 stats={"genes": len(g), "named": named, "compartment": r.compartment,
                        "share": r.share, "q": r.q}, sources=[os.path.join(F3.MODS, "modules.tsv")],
                 what=f"module M{r.module}: every gene across the cell cycle")
    PP.compose("development_figure_3a", "development", "Figure 3 per named module: the "
               "cell-cycle trajectory of every gene.")
    PP.compose("development_figure_3b", "development", "Figure 3 per module that does not clear "
               "chance: the cell-cycle trajectory of every gene.")


def dev6():
    import build_figure06_v2 as F6
    F = "development_figure_6"
    z = F6.profiles()
    for i, k in enumerate(F6.WAVE + F6.CONTINUOUS):
        g = F6.group(z, k)
        fig, ax = small()
        col = FS.ACCENT if k in F6.WAVE else FS.CONTROL
        for _, row in g.iterrows():
            ax.plot(F6.HOURS, row.to_numpy(), color=col, lw=0.25, alpha=0.12)
        ax.plot(F6.HOURS, g.mean(), color=FS.INK, lw=1.0, label="mean")
        ax.set_xticks([3, 6, 9, 12]); ax.set_xlabel("Time (h)"); ax.set_ylabel("Expression (z)")
        name_legend(ax, f"{k} ({len(g)} genes)")
        PP.panel(F, LETTERS[i], fig=fig, data=g.rename_axis("gene_id"), stats={"genes": len(g)},
                 sources=[(os.path.join(CACHE, "nodes.parquet"), ["lopit_map"])],
                 what=f"{k}: every gene across the cell cycle")
    PP.compose(F, "development", "Figure 6 per compartment: every gene trajectory.")


def dev8():
    F = "development_figure_8"
    SV = os.path.join(ROOT, "results", "strain_variation_by_class_2026_11_19")
    pg = pd.read_csv(os.path.join(SV, "per_gene.tsv"), sep="\t")
    n = pd.read_parquet(os.path.join(CACHE, "nodes.parquet")); n = n.set_index(n.gene_id.astype(str))
    from scipy.stats import spearmanr
    for i, c in enumerate(pg["class"].drop_duplicates()):
        s = pg[pg["class"] == c].set_index("gene_id")
        s = s.assign(plddt=n["mean_plddt"].reindex(s.index)).dropna(subset=["plddt"])
        s = s[s.ratio > 0]
        fig, ax = small()
        ax.plot(s.plddt, s.ratio, "o", color=FS.ACCENT, ms=1.5, alpha=0.5, lw=0, rasterized=True)
        ax.set_yscale("log"); ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(["0.1", "1", "10"])
        ax.minorticks_off(); ax.set_xlim(20, 100)
        ax.set_xlabel("pLDDT"); ax.set_ylabel("N/S")
        r = spearmanr(s.plddt, s.ratio) if len(s) > 3 else None
        name_legend(ax, c.split(" (")[0])
        PP.panel(F, LETTERS[i], fig=fig, data=s[["ratio", "plddt"]],
                 stats={"genes": len(s), "Spearman rho": r.statistic if r else np.nan,
                        "p": r.pvalue if r else np.nan}, sources=[os.path.join(SV, "per_gene.tsv")],
                 what=f"{c}: N/S against pLDDT")
    PP.compose(F, "development", "Figure 8 per class: strain variation against order, every gene.")


def dev10():
    F = "development_figure_10"
    ht = pd.read_parquet(os.path.join(ROOT, "results", "hosttx_per_host_gene_2026_11_18",
                                      "hosttx_host_gene_effects.parquet"))
    order = ht.groupby("effector").size().sort_values(ascending=False).index
    for i, e in enumerate(order):
        s = ht[ht.effector == e]
        fig, ax = small()
        ax.hist(s.hosttx_log2fc.clip(-4, 4), bins=40, color=FS.ACCENT, alpha=FS.ALPHA, lw=0)
        ax.axvline(0, color=FS.INK, lw=0.6)
        ax.set_xlabel("Host change (log$_2$ FC)"); ax.set_ylabel("Host genes")
        name_legend(ax, f"{e} ({len(s)} genes)")
        PP.panel(F, LETTERS[i], fig=fig, data=s[["host_id", "host_symbol", "hosttx_log2fc",
                                                 "hosttx_q"]],
                 stats={"host genes moved": len(s), "up (%)": float((s.hosttx_log2fc > 0).mean()
                                                                    * 100)},
                 sources=[os.path.join(ROOT, "results", "hosttx_per_host_gene_2026_11_18")],
                 what=f"{e}: host genes moved and their change")
    PP.compose(F, "development", "Figure 10 per effector: every host gene moved.")


def dev11():
    import build_figure11_v2 as F11
    F = "development_figure_11"
    t = pd.read_csv(os.path.join(F11.FM, "by_target_and_strategy.tsv"), sep="\t")
    for i, ((org, tgt), (name, kind)) in enumerate(F11.TARGETS.items()):
        s = t[(t.organism == org) & (t.target == tgt)]
        fig, ax = small()
        ax.plot(s.answer_rate * 100, s.right_when_answered * 100, "o", color=FS.ACCENT, ms=3,
                alpha=FS.ALPHA, lw=0)
        ax.axhline(float(s.commonest_answer_baseline.iloc[0]) * 100, color=FS.INK, lw=0.6,
                   ls="--")
        ax.set_xlim(0, 103); ax.set_ylim(0, 103)
        ax.set_xlabel("Answered (%)"); ax.set_ylabel("Correct (%)")
        name_legend(ax, name)
        PP.panel(F, LETTERS[i], fig=fig, data=s, stats={"strategies": len(s),
                 "commonest answer": float(s.commonest_answer_baseline.iloc[0])},
                 sources=[os.path.join(F11.FM, "by_target_and_strategy.tsv")],
                 what=f"{name}: every strategy")
    PP.compose(F, "development", "Figure 11 per target: every strategy, answer rate against "
               "accuracy; dashed line, commonest answer.")


if __name__ == "__main__":
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    which = sys.argv[1:] or ["2", "3", "6", "8", "10", "11"]
    for w in which:
        {"2": dev2, "3": dev3, "6": dev6, "8": dev8, "10": dev10, "11": dev11}[w]()
        print("done", w)
