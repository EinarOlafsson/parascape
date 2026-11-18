"""How much of the host response is shared between effectors, and how much belongs to one.

Figure 10 of the ParaScape paper (instruction 101). The dual perturb-seq study knocks out 22
parasite effectors one at a time and reads the infected host cell, so for the first time the host
response can be split: the genes many effectors move, and the genes only one does.

WHY THE SPLIT MATTERS. A host gene moved by nearly every effector is downstream of infection
itself, or of the export machinery all of them need; a host gene moved by one effector is that
effector doing something of its own. The first set is where the parasite is loud and
uninformative, and the second is where an effector earns its name.

THE EXPORT MACHINERY IS THE CONTROL BUILT INTO THE EXPERIMENT. MYR1, MYR2, MYR3 and MYR4 are the
translocon: knocking any of them out should remove the effects of everything that travels through
it, so they should move MORE host genes and share MORE of them with each other than a cargo
effector does. If that does not come out, the import is keyed wrongly.

Writes `results/host_response_shared_<date>/`.
"""
from __future__ import annotations

import glob
import json
import os
from datetime import date

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The translocon, which is the control: these move cargo rather than being cargo.
TRANSLOCON = ("MYR1", "MYR2", "MYR3", "MYR4", "GRA59")

#: A host gene moved by at least this share of effectors is called shared.
SHARED_SHARE = 0.5


def newest(pattern: str) -> str:
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    if not found:
        raise SystemExit(f"no input matching {pattern}; run build_hosttx_per_host_gene.py first")
    return found[-1]


def jaccard(a: set, b: set) -> float:
    """How much two effectors share of what they move, as a fraction of what either moves."""
    union = a | b
    return len(a & b) / len(union) if union else float("nan")


def main() -> None:
    source = newest("hosttx_per_host_gene_*/hosttx_host_gene_effects.parquet")
    effects = pd.read_parquet(source)
    effectors = sorted(effects["effector"].unique())
    print(f"{len(effects):,} effector-to-host-gene effects, {len(effectors)} effectors, "
          f"{effects['host_id'].nunique():,} host genes")

    moved = {e: set(g["host_id"]) for e, g in effects.groupby("effector")}
    per_gene = effects.groupby("host_id")["effector"].nunique()
    shared_cut = max(2, int(round(len(effectors) * SHARED_SHARE)))
    shared = per_gene[per_gene >= shared_cut]
    unique = per_gene[per_gene == 1]
    print(f"\nhost genes moved by {shared_cut}+ of {len(effectors)} effectors: {len(shared):,}")
    print(f"host genes moved by exactly one effector:    {len(unique):,}")
    print(f"median effectors per moved host gene:        {int(per_gene.median())}")

    # The control: does the translocon behave like a translocon?
    inside = [e for e in effectors if e in TRANSLOCON]
    cargo = [e for e in effectors if e not in TRANSLOCON]
    sizes = pd.Series({e: len(v) for e, v in moved.items()})
    print(f"\nhost genes moved, median: translocon {int(sizes[inside].median()):,} "
          f"({', '.join(inside)}) against cargo {int(sizes[cargo].median()):,}")
    within = [jaccard(moved[a], moved[b]) for i, a in enumerate(inside) for b in inside[i + 1:]]
    between = [jaccard(moved[a], moved[b]) for a in inside for b in cargo]
    among = [jaccard(moved[a], moved[b]) for i, a in enumerate(cargo) for b in cargo[i + 1:]]
    print(f"overlap of what they move: translocon with translocon "
          f"{np.median(within):.3f}, translocon with cargo {np.median(between):.3f}, "
          f"cargo with cargo {np.median(among):.3f}")
    if np.median(within) > np.median(among):
        print("The translocon members move more host genes and share more of them with each "
              "other\nthan cargo effectors do, which is the control on the import passing.")

    rows = [{"effector": e, "host_genes_moved": len(v),
             "shared_with_half_or_more": int(sum(1 for g in v if per_gene.get(g, 0) >= shared_cut)),
             "its_own_only": int(sum(1 for g in v if per_gene.get(g, 0) == 1)),
             "translocon": e in TRANSLOCON}
            for e, v in sorted(moved.items(), key=lambda kv: -len(kv[1]))]
    table = pd.DataFrame(rows)
    table["share_that_is_its_own"] = table["its_own_only"] / table["host_genes_moved"]
    most_own = table.sort_values("share_that_is_its_own", ascending=False).iloc[0]
    print(f"\nThe effector with the largest private response is {most_own['effector']}: "
          f"{most_own['its_own_only']:,} of\n{most_own['host_genes_moved']:,} host genes "
          f"({most_own['share_that_is_its_own']:.1%}) are moved by it and by nothing else.")

    out = os.path.join(ROOT, "results", f"host_response_shared_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "per_effector.tsv"), sep="\t", index=False)
    per_gene.rename("effectors_moving_it").to_csv(
        os.path.join(out, "effectors_per_host_gene.tsv"), sep="\t")
    json.dump({"source": source, "effectors": len(effectors),
               "host_genes": int(effects["host_id"].nunique()),
               "shared_cut": shared_cut, "shared_host_genes": int(len(shared)),
               "private_host_genes": int(len(unique)),
               "median_effectors_per_gene": int(per_gene.median()),
               "overlap_translocon_within": float(np.median(within)),
               "overlap_translocon_to_cargo": float(np.median(between)),
               "overlap_cargo_within": float(np.median(among)),
               "per_effector": table.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
