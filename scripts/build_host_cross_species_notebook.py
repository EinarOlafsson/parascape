"""Crossing human and mouse in the host response atlas (instruction 71, product 3), as an executed notebook.

The host atlas keeps the species apart because human and mouse signatures share no gene identifiers. This maps
both onto the gene symbol they share (`perturbation_atlas.host_symbol_map`) and asks the only question that
matters: does the host response to the same parasite perturbation look the same in a mouse cell as in a human
cell? One perturbation is stated in both species -- the BFD1 knockout -- so it is the test, and its answer is
reported against what chance would give.

The mapping is by symbol, not tested orthology: a gene family where one species has several members is dropped
rather than guessed. What that costs is counted in the notebook.

    python scripts/build_host_cross_species_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"host_cross_species_{DATE}")
NB = os.path.join(REPO, "notebooks", f"host_cross_species_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Host response: does a mouse cell answer like a human cell?",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**The open point this closes (instruction 71, product 3).** The host response atlas held human and "
          "mouse signatures apart, because they share no gene identifiers. Mapped onto the gene symbol they "
          "share, they can be compared, and the question becomes answerable: does the host response to the same "
          "parasite perturbation look the same in a mouse cell as in a human cell?",
          "**The mapping, and its limit.** Symbols only -- mouse `Bfd1` meets human `BFD1`. A symbol that names "
          "more than one gene in a species is dropped rather than guessed, so a gene family where the species "
          "differ in member count is simply absent. This is not tested orthology and is not offered as such.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, numpy as np, pandas as pd",
            "from parascape import perturbation_atlas as PA",
            "store = PA.host_store()",
            "crossed = PA.cross_species(store)",
            "reg = crossed['registry']",
            "human = [c for c in crossed['lfc'].columns if reg.loc[c, 'species'] == 'human']",
            "mouse = [c for c in crossed['lfc'].columns if reg.loc[c, 'species'] == 'mouse']",
            "both = int((crossed['lfc'][human].notna().any(axis=1) & crossed['lfc'][mouse].notna().any(axis=1)).sum())",
            "print(len(human), 'human and', len(mouse), 'mouse signatures over', crossed['lfc'].shape[0], 'symbols')",
            "print(both, 'symbols measured in both species;', crossed['dropped']['rows_without_symbol'],",
            "      'gene-signature values dropped for want of a one-to-one symbol')",
            "pd.DataFrame([crossed['dropped']])")
    nb.md("## 1 · The one perturbation both species state",
          "The treated label is normalised the way `replicate_groups` normalises it, then grouped across "
          "species. A perturbation named by only one species cannot be crossed and is listed.")
    nb.code("norm = reg.treated.astype(str).str.lower().str.replace(r'[^a-z0-9]+', '', regex=True)",
            "reg = reg.assign(norm=norm)",
            "def family(text):",
            "    for key in ('bfd1', 'bfd2', 'gra18', 'myr1', 'lps', 'ifng'):",
            "        if key in text:",
            "            return key",
            "    return ''",
            "reg['family'] = reg.norm.map(family)",
            "shared = sorted({f for f in reg.loc[mouse, 'family'] if f} & {f for f in reg.loc[human, 'family'] if f})",
            "print('perturbations stated in both species:', shared)",
            "print('mouse only:', sorted({f for f in reg.loc[mouse, 'family'] if f} - set(shared)))",
            "reg.loc[reg.family != '', ['study', 'species', 'treated', 'family']]")
    nb.md("## 2 · Does the mouse signature find the human one?",
          "Each mouse signature is correlated (Spearman, on the symbols both measure) with every human "
          "signature. The rank of the best human signature **of the same perturbation** is the answer; chance is "
          "the rank an unrelated signature would get, which is half the number of human signatures.")
    nb.code("MIN_SHARED = 500",
            "lfc = crossed['lfc']",
            "rows = []",
            "for m in mouse:",
            "    fam = reg.loc[m, 'family']",
            "    scores = {}",
            "    for h in human:",
            "        pair = lfc[[m, h]].dropna()",
            "        if len(pair) < MIN_SHARED:",
            "            continue",
            "        scores[h] = pair.corr(method='spearman').iloc[0, 1]",
            "    if not scores:",
            "        rows.append(dict(mouse_signature=m, perturbation=fam, comparable=0, verdict='no human signature shares enough symbols'))",
            "        continue",
            "    ranked = pd.Series(scores).sort_values(ascending=False)",
            "    same = [h for h in ranked.index if reg.loc[h, 'family'] == fam and fam]",
            "    best = int(min(ranked.index.get_loc(h) for h in same) + 1) if same else None",
            "    rows.append(dict(mouse_signature=m, perturbation=fam or '(not one of the shared kinds)',",
            "                     comparable=len(ranked), same_perturbation_signatures=len(same),",
            "                     best_same_rank=best, chance_rank=round((len(ranked) + 1) / 2, 1),",
            "                     best_rho=round(float(ranked.iloc[0]), 3),",
            "                     best_same_rho=round(float(ranked[same].max()), 3) if same else None,",
            "                     top_match=ranked.index[0], top_match_label=reg.loc[ranked.index[0], 'contrast_label']))",
            "X = pd.DataFrame(rows)",
            "X.to_csv(os.path.join(OUT, 'cross_species_recovery.tsv'), sep='\\t', index=False)",
            "X")
    nb.md("## 3 · What it says",
          "A mouse signature whose same-perturbation human partner ranks near the top means the host response "
          "crosses the species; one that ranks at chance means the two cells answer differently, or that the "
          "symbol mapping loses what matters. Both readings are stated, never only the hopeful one.")
    nb.code("tested = X[X.best_same_rank.notna()]",
            "summary = dict(mouse_signatures=len(X), testable=len(tested),",
            "               median_best_same_rank=(float(tested.best_same_rank.median()) if len(tested) else None),",
            "               median_chance_rank=(float(tested.chance_rank.median()) if len(tested) else None),",
            "               symbols_in_both=both)",
            "json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "print(summary)",
            "tested[['mouse_signature', 'perturbation', 'best_same_rank', 'chance_rank', 'best_same_rho', 'top_match_label']]")
    nb.write(NB)


if __name__ == "__main__":
    main()
