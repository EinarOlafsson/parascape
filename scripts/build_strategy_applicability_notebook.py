"""Which strategies run in which space (instruction 53, WP8), as an executed notebook.

Every strategy declares what it needs in words. `parascape.strategy_applicability` checks those
words against each space own data, so a space that is not Toxoplasma arrives with an honest list of
what it can and cannot do instead of a wall of failures. Nothing is run here: the question is only
whether it could, and when not, which requirement is missing.

    python scripts/build_strategy_applicability_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"strategy_applicability_{DATE}")
NB = os.path.join(REPO, "notebooks", f"strategy_applicability_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Which strategies run in which space", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this answers (instruction 53, WP8).** ParaScape has 39 strategies, written against the "
          "Toxoplasma space. Which of them can run in another space at all?",
          "**How.** Each strategy declares what it needs -- \"a categorical column\", \"two or more edge "
          "layers\", \"the other organism's table\" -- and those words are checked against what each space "
          "actually holds. A strategy that cannot run names the requirement it is missing.",
          "**What this is not.** A measurement of how well a strategy works in a space. That is the track "
          "record, and a space with no track record is never called calibrated.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, numpy as np, pandas as pd",
            "import parascape.strategy_graph, parascape.strategy_catalog, parascape.strategy_learning  # register them",
            "from parascape import strategies, strategy_applicability as SA",
            "SPACES = (('Tg', 'parascape/data/nodes.parquet', 'parascape/data/graph.npz'),",
            "          ('Pf', 'parascape/data/pf_nodes.parquet', 'parascape/data/pf_graph.npz'),",
            "          ('Cp', 'results/cryptosporidium_space_2026_11_08/cp_gene_nodes.parquet', None))",
            "def layer_names(relative):",
            "    if not relative:",
            "        return ()",
            "    path = os.path.join(REPO, relative)",
            "    if not os.path.exists(path):",
            "        return ()",
            "    return tuple(np.load(path, allow_pickle=True).files)",
            "print(len(strategies.REGISTRY), 'strategies registered')",
            "rows = []",
            "for code, table, graph in SPACES:",
            "    nodes = pd.read_parquet(os.path.join(REPO, table))",
            "    partner = None if code != 'Cp' else False   # no Cryptosporidium counterpart table is built",
            "    rows.append(SA.evidence_of(code, nodes, layer_names(graph), partner=partner))",
            "E = pd.DataFrame([vars(e) | {'edge_layers': len(e.edge_layers)} for e in rows])",
            "E.to_csv(os.path.join(OUT, 'space_evidence.tsv'), sep='\\t', index=False)",
            "E")
    nb.md("## 1 · The matrix",
          "One row per strategy and space. A space with no graph loses the network strategies; a space with "
          "no counterpart table loses the transfer strategy; a space with no literature column loses the two "
          "that read attention.")
    nb.code("frames = [SA.applicability(e) for e in rows]",
            "A = pd.concat(frames, ignore_index=True)",
            "A.to_csv(os.path.join(OUT, 'applicability.tsv'), sep='\\t', index=False)",
            "counts = A.groupby('space').applicable.agg(['sum', 'count'])",
            "counts.columns = ['applicable', 'strategies']",
            "counts['share'] = (100 * counts.applicable / counts.strategies).round(0)",
            "json.dump({r.Index: int(r.applicable) for r in counts.itertuples()},",
            "          open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "counts")
    nb.md("## 2 · What each space cannot do, and why",
          "The useful half of the answer: the missing requirement names what to build next for that space.")
    nb.code("blocked = A[~A.applicable][['space', 'number', 'strategy', 'family', 'missing']]",
            "blocked.to_csv(os.path.join(OUT, 'blocked.tsv'), sep='\\t', index=False)",
            "print(blocked.groupby(['space', 'missing']).size().to_string())",
            "blocked")
    nb.md("## 3 · Reading it",
          "For Cryptosporidium every missing requirement is an edge layer or a counterpart table, not a "
          "strategy that fails to apply to the organism: build its graph and its Toxoplasma counterpart "
          "table and most of the blocked strategies open. That is a build order, which is what WP8 was for.")
    nb.code("cp = blocked[blocked.space == 'Cp']",
            "needs_graph = cp[cp.missing.str.contains('layer')]",
            "print(len(needs_graph), 'of', len(cp), 'blocked Cryptosporidium strategies wait on a graph')",
            "print(sorted(set(cp.missing) - set(needs_graph.missing)))",
            "cp.groupby('family').size().to_frame('blocked')")
    nb.write(NB)


if __name__ == "__main__":
    main()
