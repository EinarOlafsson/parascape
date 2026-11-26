"""Does the host-facing cluster hold up among little-studied genes? (Einar 2026-11-26: is study attention a problem?)

The E8 all-views clustering (Louvain 1.0) is built twice: with every view, and without the views that count study
effort (the interaction view: numbers of BioID / IP-MS studies and crosslink partners). For each, per attention
stratum (papers on the gene: none, 1-2, 3-10, more than 10): the measured host-facing genes ("effect on host"),
how many land in the most host-facing cluster, and the enrichment over the stratum's own rate of landing there.

    python scripts/infer_bench/host_attention_check.py OUT_DIR
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e8_selected_clustering as E8  # noqa: E402
import heldout_sets as HS  # noqa: E402


def main(argv=None) -> int:
    out = (argv if argv is not None else sys.argv[1:])[0]
    os.makedirs(out, exist_ok=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    eff = np.asarray(L["effect on host"][0], dtype=object) == "host-facing"
    papers = (nodes["n_papers_focal"].fillna(0) + nodes["n_papers_substantive"].fillna(0)
              + nodes["n_papers_incidental"].fillna(0))
    stratum = HS.attention_stratum(papers).to_numpy()
    rows = []
    for variant, keep in (("all views", lambda v: True), ("without study-effort counts", lambda v: v != "interaction")):
        feats = {v: c for v, c in views.items() if len(c) >= 2 and keep(v)}
        E, idx = E8.embed(nodes, views, extra, feats, list(extra))
        c = np.full(len(nodes), -1)
        c[idx] = E5.louvain(E, 1.0, 0)
        top = int(pd.Series(c[eff & (c >= 0)]).value_counts().index[0])
        inside = c == top
        for s in ("none", "1-2", "3-10", "more than 10", "all"):
            m = (stratum == s) if s != "all" else np.ones(len(nodes), bool)
            pos = m & eff
            base = inside[m & (c >= 0)].mean() if (m & (c >= 0)).any() else np.nan
            share = inside[pos].mean() if pos.any() else np.nan
            rows.append({"variant": variant, "papers": s, "genes": int(m.sum()), "host_facing": int(pos.sum()),
                         "in_host_cluster": int((pos & inside).sum()), "share": share, "stratum_base": base,
                         "enrichment": share / base if base else np.nan, "cluster_size": int(inside.sum())})
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(out, "attention_check.csv"), index=False)
    print(T.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
