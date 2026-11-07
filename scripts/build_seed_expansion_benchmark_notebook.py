"""Held-out member recovery for seed expansion (instruction 64.12, the first ranking adapter), as a notebook.

For each localization niche, its curation-grade hyperLOPIT markers are split by orthogroup into seeds and hidden
members (seed 0, half each). Seed expansion runs from the seeds with the hyperLOPIT compartment excluded, and the
hidden members are looked for in the ranked genome. Depths (50, 100, 200) are fixed here, before any run.

    python scripts/build_seed_expansion_benchmark_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"seed_expansion_benchmark_{DATE}")
NB = os.path.join(REPO, "notebooks", f"seed_expansion_benchmark_{DATE}.ipynb")


def one_niche(niche: str, out_path: str, strategy: str = "seed_expansion", split_seed: int = 0,
              top: int = 2000) -> None:
    """Run seed expansion for one niche and write its recovery row (run in its own process)."""
    import json
    import numpy as np
    import pandas as pd
    sys.path.insert(0, REPO)
    from parascape import strategies as S
    ctx = S.Context.shipped("Tg")
    nodes = ctx.nodes.set_index("gene_id")
    m = pd.read_csv(os.path.join(REPO, "parascape/data/lopit_markers_coarse.tsv"), sep="\t")
    members = [g for g in m[m.coarse_niche == niche].gene_id if g in nodes.index]
    rng = np.random.default_rng(split_seed)
    og = nodes.loc[members, "orthogroup"].fillna(pd.Series(members, index=members))
    groups = sorted(pd.unique(og), key=str)
    rng.shuffle(groups)
    seed_groups = set(groups[:len(groups) // 2])
    seeds = [g for g in members if og[g] in seed_groups]
    hidden = [g for g in members if og[g] not in seed_groups]
    extra = {"top": top} if any(p.name == "top" for p in S.get(strategy).params) else {}
    r = S.get(strategy).run(ctx, genes="\n".join(seeds), exclude="compartment", **extra)
    table = next(t for t in r.tables.values() if "gene_id" in getattr(t, "columns", []))
    ranked = [g for g in table.gene_id.astype(str) if g not in set(seeds)]
    hits = np.isin(ranked, hidden)
    prevalence = len(hidden) / (len(nodes) - len(seeds))
    row = dict(strategy=strategy, niche=niche, split_seed=split_seed, seeds=len(seeds), hidden=len(hidden),
               prevalence=round(prevalence, 5))
    # Average precision over the ranking; hidden members the ranking never reaches count as misses at the end.
    found = np.flatnonzero(hits)
    row["average_precision"] = round(float(sum((i + 1) / (pos + 1) for i, pos in enumerate(found)) / len(hidden)), 4)
    for k in (50, 100, 200):
        row[f"hidden_in_top{k}"] = int(hits[:k].sum())
        row[f"recall_at_{k}"] = round(float(hits[:k].sum() / len(hidden)), 3)
        row[f"lift_at_{k}"] = round(float(hits[:k].mean() / prevalence), 1)
    with open(out_path, "w") as fh:
        json.dump(row, fh)


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Does seed expansion find the members it was not given?",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**The test (instruction 64.12).** Half of a niche's hyperLOPIT markers, whole orthogroups at a time, "
          "are given as seeds; the other half are hidden. Seed expansion (random walk with restart over the "
          "measured networks and the measurement space) ranks the genome with the hyperLOPIT compartment "
          "excluded, and the hidden members are counted at fixed depths. The candidate universe is every gene "
          "except the seeds. **What this is not:** the genes that are not markers are not verified negatives, "
          "so precision here is a lower bound on how often a ranked gene is truly in the niche.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, subprocess, sys",
            "from concurrent.futures import ThreadPoolExecutor",
            "import pandas as pd",
            "NICHES = ['rhoptries', 'micronemes', 'dense granules', 'IMC', 'apicoplast', 'mitochondrion']",
            "STRATEGIES = ['seed_expansion', 'positive_unlabeled', 'set_enrichment', 'geneset_hunt']",
            "SCRIPT = os.path.join(REPO, 'scripts', 'build_seed_expansion_benchmark_notebook.py')",
            "def run(job):",
            "    strategy, niche = job",
            "    tag = niche.replace(' ', '_')",
            "    out = os.path.join(OUT, f'row_{tag}.json' if strategy == 'seed_expansion' else f'row_{strategy}_{tag}.json')",
            "    if not os.path.exists(out):      # a finished niche is reused, never silently recomputed",
            "        subprocess.run([sys.executable, SCRIPT, '--niche', niche, '--out', out, '--strategy', strategy],",
            "                       check=True, env={**os.environ, 'OMP_NUM_THREADS': '2'})",
            "    row = json.load(open(out)); row.setdefault('strategy', strategy)",
            "    return row",
            "jobs = [(s, n) for s in STRATEGIES for n in NICHES]",
            "with ThreadPoolExecutor(max_workers=len(NICHES)) as pool:",
            "    rows = list(pool.map(run, jobs))",
            "table = pd.DataFrame(rows)",
            "table.to_csv(os.path.join(OUT, 'recovery.tsv'), sep='\\t', index=False)",
            "table.pivot(index='niche', columns='strategy', values='recall_at_100')")
    nb.md("### Five splits, full ranking: positive-unlabelled learning",
          "The fast strategy repeated over five orthogroup splits (seeds 0-4), ranking every gene, with average "
          "precision over the whole ranking (hidden members it never reaches count as misses). Mean and range over "
          "splits; prevalence is the average precision of a random ranking.")
    nb.code("REP = os.path.join(OUT, 'repeats'); os.makedirs(REP, exist_ok=True)",
            "def rep(job):",
            "    seed, niche = job",
            "    out = os.path.join(REP, f\"pu_{seed}_{niche.replace(' ', '_')}.json\")",
            "    if not os.path.exists(out):",
            "        subprocess.run([sys.executable, SCRIPT, '--niche', niche, '--out', out, '--strategy', 'positive_unlabeled',",
            "                        '--split-seed', str(seed), '--top', '8140'], check=True, env={**os.environ, 'OMP_NUM_THREADS': '2'})",
            "    return json.load(open(out))",
            "with ThreadPoolExecutor(max_workers=len(NICHES)) as pool:",
            "    reps = pd.DataFrame(list(pool.map(rep, [(s, n) for s in range(5) for n in NICHES])))",
            "summary = reps.groupby('niche').agg(splits=('split_seed', 'size'), ap_mean=('average_precision', 'mean'),",
            "    ap_min=('average_precision', 'min'), ap_max=('average_precision', 'max'),",
            "    recall100_mean=('recall_at_100', 'mean'), prevalence=('prevalence', 'mean')).round(3)",
            "summary['ap_lift'] = (summary.ap_mean / summary.prevalence).round(1)",
            "summary.to_csv(os.path.join(OUT, 'pu_repeats_summary.tsv'), sep='\\t')",
            "summary")
    nb.write(NB)


if __name__ == "__main__":
    if "--niche" in sys.argv:
        arg = lambda name, default: sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default
        one_niche(arg("--niche", ""), arg("--out", ""), arg("--strategy", "seed_expansion"),
                  int(arg("--split-seed", 0)), int(arg("--top", 2000)))
    else:
        main()
