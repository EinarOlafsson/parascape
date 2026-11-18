"""The corpus, assembled: what ParaScape has imported and how much of the genome it reaches.

Figure 1 of the ParaScape paper (instruction 101). The schematic around it is a drawing; this is
the part that has to be measured, because a figure claiming a corpus has to say how big it is and
how much of it is empty.

WHAT IS COUNTED. Datasets by organism and by level of biology, the slots of the catalogue that
carry something, and -- the number that matters -- how much of the genome each kind of knowledge
actually reaches. A dataset list is a claim about effort; coverage is a claim about knowledge.

A COLUMN COMPUTED FOR EVERY GENE IS NOT COVERAGE. The same rule as figure 12: a count that is zero
where nobody looked, or a flag that is False for an unassayed gene, would report a category as
fully covered. Those columns are set aside and counted.

Writes `results/corpus_inventory_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from datetime import date

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import datasets as DS                   # noqa: E402
from parascape import figstyle as FS                   # noqa: E402
from parascape import information_categories as IC     # noqa: E402
from parascape import organisms as ORGANISMS           # noqa: E402

#: Categories that are not knowledge about a gene.
NOT_KNOWLEDGE = ("identity", "attention")

def organism_name(code: str) -> str:
    """How an organism code reads on a figure, from `parascape.organisms`.

    Read from the organism table rather than written out here. A second copy of the codes in a
    figure script is the bug the project calls #16i, and `test_organisms` refuses new literals.
    """
    space = ORGANISMS.SPACES.get(code)
    if space is not None:
        return space.species.split()[0]
    # The host entries are table filenames, not species, so the reading word is supplied here --
    # keyed on the constants `organisms` already exports, so no organism code is written out.
    hosts = {ORGANISMS.HUMAN: "human", ORGANISMS.MOUSE: "mouse",
             ORGANISMS.ANOPHELES: "Anopheles"}
    return hosts.get(code, code)


def by_organism_and_level() -> pd.DataFrame:
    """Datasets counted by organism and by level of biology."""
    rows = Counter((organism_name(d.organism), d.level) for d in DS.REGISTRY)
    frame = pd.DataFrame([{"organism": o, "level": lv, "datasets": n}
                          for (o, lv), n in rows.items()])
    return frame.sort_values("datasets", ascending=False)


def coverage(cache: str) -> pd.DataFrame:
    """Per category of knowledge, the share of genes anything has ever measured."""
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    derived = set(getattr(IC, "DERIVED_FROM", {}))
    columns, set_aside = {}, 0
    for column in nodes.columns:
        if column in derived or not nodes[column].isna().any():
            set_aside += 1
            continue
        for category in IC.categories_of(column):
            if category not in NOT_KNOWLEDGE:
                columns.setdefault(category, []).append(column)
    rows = [{"category": c, "columns": len(cols),
             "share_measured": float(nodes[cols].notna().any(axis=1).mean())}
            for c, cols in columns.items()]
    table = pd.DataFrame(rows).sort_values("share_measured", ascending=False)
    table.attrs["genes"] = len(nodes)
    table.attrs["set_aside"] = set_aside
    return table


def draw(table: pd.DataFrame, letter: str = "b"):
    """How much of the genome each kind of knowledge reaches.

    Bars, in the house geometry: one slot per category, uniform width and gap. The categories that
    reach nearly everything are the reference against which the empty ones are the point, so they
    are grey and the three least covered carry the first series colour.
    """
    FS.use()
    order = table.sort_values("share_measured")
    worst = set(order["category"].head(3))
    fig, ax = FS.figure(len(order), height_mm=56.0)
    y = list(range(len(order)))
    ax.barh(y, order["share_measured"].to_numpy(), height=FS.bar_width(), linewidth=0,
            color=[FS.SERIES[0] if c in worst else FS.CONTROL for c in order["category"]],
            alpha=FS.ALPHA)
    for i, row in enumerate(order.itertuples()):
        ax.annotate(f"{row.share_measured:.1%}", (row.share_measured, i), xytext=(3, 0),
                    textcoords="offset points", va="center", ha="left", fontsize=5.6,
                    color=FS.SERIES[0] if row.category in worst else FS.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ") for c in order["category"]], fontsize=6.0)
    ax.set_xlabel(f"share of {table.attrs['genes']:,} genes anything has measured", fontsize=6.0,
                  labelpad=2)
    ax.set_xlim(0, 1.16)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "50%", "100%"], fontsize=5.5)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.spines["left"].set_visible(False)
    FS.letter(fig, letter)
    return fig


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    datasets = by_organism_and_level()
    print(f"{len(DS.REGISTRY)} datasets")
    print(datasets.groupby("organism")["datasets"].sum().sort_values(ascending=False).to_string())
    print()
    print(datasets.groupby("level")["datasets"].sum().sort_values(ascending=False).to_string())

    table = coverage(cache)
    print(f"\n{table.attrs['genes']:,} genes; {table.attrs['set_aside']} columns set aside as "
          f"derivations, flags or counts of zero\n")
    for row in table.itertuples():
        print(f"  {row.category:18s} {row.share_measured:6.1%} of genes   "
              f"{row.columns:3d} columns")

    out = os.path.join(ROOT, "results", f"corpus_inventory_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    datasets.to_csv(os.path.join(out, "datasets.tsv"), sep="\t", index=False)
    table.to_csv(os.path.join(out, "coverage.tsv"), sep="\t", index=False)
    fig = draw(table)
    for extension in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"panel_coverage.{extension}"),
                    dpi=600 if extension == "png" else None, bbox_inches="tight", pad_inches=0.02)
    json.dump({"datasets": len(DS.REGISTRY), "genes": int(table.attrs["genes"]),
               "columns_set_aside": int(table.attrs["set_aside"]),
               "by_organism": datasets.groupby("organism")["datasets"].sum().to_dict(),
               "by_level": datasets.groupby("level")["datasets"].sum().to_dict(),
               "coverage": table.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
