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
    order = table.sort_values("share_measured", ascending=False)
    worst = set(order["category"].tail(3))
    fig, ax = FS.figure(len(order), vertical=True)
    names = [c.replace("_", " ") for c in order["category"]]
    FS.bars(ax, order["share_measured"].to_numpy(), names, highlight=worst and
            {n for n, c in zip(names, order["category"]) if c in worst})
    ax.set_ylabel("share measured")
    ax.set_ylim(0, 1.06)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "50%", "100%"])
    ax.set_xlim(-0.6, len(order) - 0.4)
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
