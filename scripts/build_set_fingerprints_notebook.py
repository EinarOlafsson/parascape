"""Set fingerprints for the host-facing compartment ladder (instruction 82), as an executed notebook.

A set fingerprint is the data profile shared by genes known to share an annotation; an unannotated
gene is ranked by how closely it matches. The sets here are the rungs of the compartment ladder, taken
ONLY from the graded, PMID-verified host-facing review (results/host_facing_truth_2026_11_04): nothing is
added by guessing. Every fingerprint is built from measurements that say nothing about place
(`information_categories`: every `place` column is excluded), scored on held-out members, set against a
random baseline and against dense granule proteins in general, and the confusion between rungs reported.

    python scripts/build_set_fingerprints_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

LITERATURE = "--literature" in sys.argv   # add results/literature_set_candidates_2026_11_08/literature_rungs.tsv
DATE = "2026_11_08" if LITERATURE else "2026_11_06"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"set_fingerprints_{DATE}")
NB = os.path.join(REPO, "notebooks", f"set_fingerprints_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Set fingerprints for the host-facing compartment ladder",
                          namespace={"REPO": REPO, "OUT": OUT, "LITERATURE": LITERATURE})
    nb.md(
        "**What this is.** Instruction 82 asks for fingerprints of *gene sets* defined by the "
        "literature -- proteins on the host-facing side of the vacuole membrane, proteins exported "
        "into the host cytosol or nucleus -- so that an unannotated protein can be ranked by how "
        "much it resembles each set. The sets are the rungs of a **compartment ladder**, and the "
        "rung that matters most is the one host-facing inference never had: proteins that are "
        "*vacuole-resident* (the hard negatives).",
        "**Where the sets come from.** Only the graded host-facing review of 2026-11-04 (140 "
        "claims, each with a PMID and an evidence strength; 98 genes decided). Nothing is added by "
        "guessing; a rung with too few members is reported and not scored.",
        "**What a fingerprint may use.** Every measurement *except* those about place "
        "(`information_categories`, category `place`) -- otherwise a protein would resemble the "
        "dense granules because it was measured to be in them, which is circular.",
        f"Built {_dt.date.today().isoformat()}. Inputs are read-only.")

    nb.code(
        "import os, numpy as np, pandas as pd",
        "import sys; sys.path.insert(0, REPO)",
        "from parascape import information_categories as IC",
        "truth = pd.read_csv(os.path.join(REPO, 'results', 'host_facing_truth_2026_11_04',",
        "                                 'host_facing_truth.tsv'), sep='\\t')",
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape', 'data', 'nodes.parquet'))",
        "print('reviewed genes:', len(truth), '| genes in the table:', len(nodes))",
        "truth.groupby(['deciding_place', 'decision']).size().unstack(fill_value=0)")

    nb.md("## 1 · The ladder, from the review",
          "Each reviewed gene is placed on one rung by its deciding evidence. The two PVM sides "
          "(rungs 2 and 3) cannot be told apart in this evidence: the review records *PVM, "
          "host-facing* (decision yes) or *not host-facing* (decision no), so the ladder here has a "
          "**host-facing PVM** rung and a **vacuole-resident** rung, which is exactly the "
          "distinction host-facing inference lacked.")
    nb.code(
        "RUNG = {('moving_junction', 'yes'): '1 injected at the moving junction',",
        "        ('PVM', 'yes'): '3 PVM, host-facing',",
        "        ('host_cytosol', 'yes'): '5 exported to the host cytosol',",
        "        ('host_nucleus', 'yes'): '6 host nucleus',",
        "        ('host_organelle', 'yes'): '7 host organelle or plasma membrane',",
        "        ('host_plasma_membrane', 'yes'): '7 host organelle or plasma membrane'}",
        "def rung_of(r):",
        "    if r.decision == 'no':",
        "        return '4 vacuole-resident (hard negative)'",
        "    return RUNG.get((r.deciding_place, r.decision), '')",
        "truth['rung'] = truth.apply(rung_of, axis=1)",
        "ladder = truth[truth.rung != ''][['gene_id', 'protein', 'rung', 'deciding_strength',",
        "                                  'deciding_pmid']]",
        "changes = pd.DataFrame()",
        "if LITERATURE:   # literature rungs (scripts/merge_literature_ladder.py): new genes join; PVM genes take the quoted side",
        "    lit = pd.read_csv(os.path.join(REPO, 'results/literature_set_candidates_2026_11_08/literature_rungs.tsv'), sep='\\t')",
        "    lit = lit.assign(deciding_strength='direct', deciding_pmid=lit.pmids)",
        "    old = ladder.set_index('gene_id').rung",
        "    pvm = set(truth[truth.deciding_place.astype(str).str.startswith(('PVM', 'PV_'))].gene_id)",
        "    take = lit[~lit.gene_id.isin(old.index) | lit.gene_id.isin(pvm)]",
        "    changes = take.assign(review_rung=take.gene_id.map(old).fillna('(not in the review)'))",
        "    changes = changes[changes.review_rung != changes.rung][['gene_id', 'protein', 'review_rung', 'rung', 'pmids']]",
        "    ladder = pd.concat([ladder[~ladder.gene_id.isin(take.gene_id)],",
        "                        take[['gene_id', 'protein', 'rung', 'deciding_strength', 'deciding_pmid']]], ignore_index=True)",
        "    changes.to_csv(os.path.join(OUT, 'literature_changes.tsv'), sep='\\t', index=False)",
        "    print(len(changes), 'genes added or moved by the literature')",
        "ladder = ladder[ladder.gene_id.isin(set(nodes.gene_id.astype(str)))]",
        "ladder.to_csv(os.path.join(OUT, 'ladder.tsv'), sep='\\t', index=False)",
        "counts = ladder.groupby('rung').agg(genes=('gene_id', 'size'),",
        "                                    direct=('deciding_strength', lambda s: int((s == 'direct').sum())))",
        "MIN_MEMBERS = 5",
        "counts['scored'] = counts.genes >= MIN_MEMBERS",
        "counts")

    nb.md("## 2 · Profiles that say nothing about place",
          "Numeric columns of the gene table, minus every column whose information category is "
          "`place` (and identity and bookkeeping). Each column is rank-scaled, missing values "
          "set to the column median, and the profiles reduced to 50 principal components so a "
          "fingerprint is not dominated by the hundreds of near-duplicate expression columns.")
    nb.code(
        "num = nodes.select_dtypes(include=[np.number]).columns",
        "def allowed(c):",
        "    try:",
        "        cats = IC.categories_of(c)",
        "    except Exception:",
        "        return False",
        "    return bool(cats) and not set(cats) & {'place', 'identity', 'attention'}",
        "cols = [c for c in num if allowed(c) and nodes[c].notna().mean() >= 0.2]",
        "X = nodes[cols].rank(pct=True)",
        "X = X.fillna(X.median()).to_numpy()",
        "X = (X - X.mean(0)) / (X.std(0) + 1e-9)",
        "from sklearn.decomposition import PCA",
        "Z = PCA(n_components=min(50, X.shape[1]), random_state=0).fit_transform(X)",
        "Z = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-9)",
        "gid = nodes.gene_id.astype(str).to_numpy()",
        "pos = {g: i for i, g in enumerate(gid)}",
        "excluded_place = [c for c in num if not allowed(c)]",
        "print(len(cols), 'columns used;', len(excluded_place),",
        "      'numeric columns excluded (place, identity, bookkeeping, or no category)')",
        "print('example exclusions:', [c for c in excluded_place if 'lopit' in c.lower()][:6])")

    nb.md("## 3 · Does each fingerprint recover its own held-out members?",
          "Leave one member out, build the fingerprint (the mean profile) from the rest, and rank "
          "the held-out member among every gene that is not a member. The score is the AUROC: 0.5 "
          "is chance, 1 is perfect. Two baselines, both of the same size: **random genes**, and "
          "**dense granule proteins** (hyperLOPIT), because most of these rungs are dense granule "
          "proteins and the question is whether a fingerprint knows more than that.")
    nb.code(
        "rng = np.random.default_rng(0)",
        "def loo_auroc(members):",
        "    idx = np.array([pos[g] for g in members if g in pos])",
        "    if len(idx) < 3:",
        "        return np.nan",
        "    others = np.setdiff1d(np.arange(len(Z)), idx)",
        "    aucs = []",
        "    for k, i in enumerate(idx):",
        "        rest = np.delete(idx, k)",
        "        centre = Z[rest].mean(0)",
        "        s_in = Z[i] @ centre",
        "        s_out = Z[others] @ centre",
        "        aucs.append(float((s_out < s_in).mean() + 0.5 * (s_out == s_in).mean()))",
        "    return float(np.mean(aucs))",
        "comp = nodes['compartment'].astype(str) if 'compartment' in nodes else pd.Series('', index=nodes.index)",
        "dg = gid[comp.str.lower().str.contains('dense gran').to_numpy()]",
        "rows = []",
        "for rung, grp in ladder.groupby('rung'):",
        "    members = grp.gene_id.tolist()",
        "    if len(members) < MIN_MEMBERS:",
        "        rows.append(dict(rung=rung, genes=len(members), fingerprint=np.nan,",
        "                         random=np.nan, dense_granule=np.nan, verdict='too few to score'))",
        "        continue",
        "    fp = loo_auroc(members)",
        "    rand = np.nanmean([loo_auroc(rng.choice(gid, len(members), replace=False)) for _ in range(20)])",
        "    dgb = np.nanmean([loo_auroc(rng.choice(dg, min(len(members), len(dg)), replace=False))",
        "                      for _ in range(20)]) if len(dg) >= len(members) else np.nan",
        "    best = np.nanmax([rand, dgb])",
        "    rows.append(dict(rung=rung, genes=len(members), fingerprint=round(fp, 3),",
        "                     random=round(rand, 3), dense_granule=round(dgb, 3),",
        "                     verdict='above the best baseline' if fp > best + 0.02 else 'not above the best baseline'))",
        "recovery = pd.DataFrame(rows)",
        "recovery.to_csv(os.path.join(OUT, 'recovery.tsv'), sep='\\t', index=False)",
        "recovery")

    nb.md("## 4 · Can the rungs be told apart?",
          "Each held-out member is assigned to the rung whose fingerprint (built without it) it "
          "resembles most. The row is the true rung, the column the assigned one. The cell that "
          "matters is **host-facing PVM** against **vacuole-resident**: if those are confused, the "
          "data cannot separate the two, which is a finding about the evidence.")
    nb.code(
        "scored = recovery[recovery.verdict != 'too few to score'].rung.tolist()",
        "sets = {r: [pos[g] for g in ladder[ladder.rung == r].gene_id if g in pos] for r in scored}",
        "conf = pd.DataFrame(0, index=scored, columns=scored)",
        "for r, idx in sets.items():",
        "    for i in idx:",
        "        best, best_s = None, -9",
        "        for r2, idx2 in sets.items():",
        "            rest = [j for j in idx2 if j != i]",
        "            if not rest:",
        "                continue",
        "            c = Z[rest].mean(0); c = c / (np.linalg.norm(c) + 1e-9)   # cosine: fair across set sizes",
        "            s = Z[i] @ c",
        "            if s > best_s:",
        "                best, best_s = r2, s",
        "        conf.loc[r, best] += 1",
        "conf.to_csv(os.path.join(OUT, 'confusion.tsv'), sep='\\t')",
        "conf")

    nb.md("## 5 · Candidates: resembles host-facing, not vacuole-resident",
          "A ranking, not a fill: the rungs hold tens of genes, far below the 400 labels a frozen "
          "test needs (instruction 73), so no cell of the information map is written from this. "
          "For each unreviewed gene: similarity to the host-facing PVM fingerprint minus similarity "
          "to the vacuole-resident one. The top of the list is where to look next.")
    nb.code(
        "hf = [pos[g] for g in ladder[ladder.rung.str.startswith('3')].gene_id if g in pos]",
        "vr = [pos[g] for g in ladder[ladder.rung.str.startswith('4')].gene_id if g in pos]",
        "reviewed = set(pos[g] for g in truth.gene_id if g in pos)",
        "unit = lambda v: v / (np.linalg.norm(v) + 1e-9)",
        "like_hf = Z @ unit(Z[hf].mean(0))",
        "score = like_hf - Z @ unit(Z[vr].mean(0))",
        "# The difference alone rewards proteins unlike BOTH sets (apicoplast, nucleus); a candidate",
        "# must first resemble the host-facing set -- the top 5% -- and is then ranked by the difference.",
        "near = like_hf >= np.quantile(like_hf, 0.95)",
        "cand = pd.DataFrame({'gene_id': gid, 'product': nodes.get('product', pd.Series('', index=nodes.index)).astype(str),",
        "                     'compartment': comp, 'resembles_host_facing': like_hf,",
        "                     'host_facing_minus_resident': score})",
        "cand = cand[near & ~cand.index.isin(reviewed)].sort_values('host_facing_minus_resident',",
        "                                                           ascending=False)",
        "print(len(cand), 'unreviewed genes in the top 5% by resemblance to host-facing PVM')",
        "print(cand.compartment.value_counts().head(8).to_string())",
        "cand.head(200).to_csv(os.path.join(OUT, 'candidates.tsv'), sep='\\t', index=False)",
        "cand.head(25)")

    nb.md("## 6 · What this shows, and what it does not",
          "* The ladder is only as large as its graded sources; a rung with fewer than 5 members is "
          "not scored (section 1) and needs more papers read.",
          "* A fingerprint that is not above the dense-granule baseline knows no more than "
          "'this is a dense granule protein' -- reported as such.",
          "* **Host-facing and vacuole-resident are not told apart by these profiles**: held-out "
          "host-facing members resemble the vacuole-resident fingerprint more often than their own. "
          "That is the evidence limit behind host-facing precision, now measured. Separating the "
          "two needs the literature-graded sides of the membrane (rungs 2 and 3) and data that sees "
          "the membrane side (protease protection, proximity labelling from the host side).",
          "* The candidates are a ranking for reading, never written into the map.")
    nb.write(NB)
    _readme()


def _md(frame, index: bool = False) -> str:
    """A pipe table without the optional `tabulate` dependency."""
    frame = frame.reset_index() if index else frame
    head = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def _readme() -> None:
    import pandas as pd
    rec = pd.read_csv(os.path.join(OUT, "recovery.tsv"), sep="\t")
    conf = pd.read_csv(os.path.join(OUT, "confusion.tsv"), sep="\t", index_col=0)
    lines = ["# Set fingerprints for the host-facing compartment ladder (instruction 82)", "",
             f"Notebook: `notebooks/set_fingerprints_{DATE}.ipynb` (executed). Builder: "
             "`scripts/build_set_fingerprints_notebook.py`. Ladder from "
             "`results/host_facing_truth_2026_11_04`" + (
                 " plus the literature rungs (`literature_changes.tsv` lists every gene added or moved)"
                 if LITERATURE else " only") + "; place columns excluded.", "",
             "## Held-out recovery (AUROC, leave-one-out)", "", _md(rec), "",
             "## Confusion between scored rungs (row: true, column: assigned)", "",
             _md(conf, index=True), "",
             "Files: `ladder.tsv`, `recovery.tsv`, `confusion.tsv`, `candidates.tsv` (top 200, a "
             "ranking only)."]
    with open(os.path.join(OUT, "README.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
