"""Pilot of the resumable production track record (instruction 64.25), as an executed notebook.

Three methods on one Toxoplasma label, run as fingerprinted partitions: stopped after one partition, resumed,
replayed with nothing to build, then one method setting changed so only that partition rebuilds. Runtime, peak
memory and bytes are reported per partition, and the calls are compared with the shipped track record.

    python scripts/build_precompute_track_record_pilot_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"precompute_track_record_pilot_{DATE}")
NB = os.path.join(REPO, "notebooks", f"precompute_track_record_pilot_{DATE}.ipynb")
TARGET = "screen_egress_phenotype"
METHODS = ["feature_knn", "map_neighbours", "layer_vote"]


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Resumable track record: pilot", namespace={"REPO": REPO, "OUT": OUT, "TARGET": TARGET,
                                                                     "METHODS": METHODS})
    nb.md(f"**The pilot (instruction 64.25).** The production hold-out computation "
          f"(`track_record.evaluate` + set hold-outs, exactly as `build_track_record.py`) for {len(METHODS)} methods "
          f"on the Toxoplasma label `{TARGET}`, cut into one partition per method. Each partition carries a "
          "fingerprint of the node table bytes, the method settings for this label and the code that computes "
          "the calls.",
          "**Pre-declared steps:** (1) run with a budget of one partition, so it stops after the first; (2) resume: "
          "only the other two are built; (3) replay: nothing is built; (4) change one setting of one method: only "
          "that partition rebuilds; (5) compare with the shipped track record.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, json, shutil, dataclasses, pandas as pd",
            "from parascape import precompute_track_record as P, track_record as T",
            "work = os.path.join(OUT, 'run'); shutil.rmtree(work, ignore_errors=True)",
            "parts = P.plan('Tg', [TARGET], METHODS)",
            "pd.DataFrame([dict(method=p.strategy, settings=p.settings, fingerprint=p.fingerprint[:16]) for p in parts])")
    nb.md("## 1-3 · Stop, resume, replay")
    nb.code("first = P.run(parts, work, P.Budget(max_jobs=1, max_seconds=3600, memory_limit_bytes=12 << 30))",
            "second = P.run(parts, work, P.Budget(max_seconds=3600, memory_limit_bytes=12 << 30))",
            "third = P.run(parts, work)",
            "steps = pd.DataFrame({'stop after one': first, 'resume': second, 'replay': third})",
            "assert list(first.values()) == ['built', 'pending', 'pending']",
            "assert list(second.values()) == ['reused', 'built', 'built'] and set(third.values()) == {'reused'}",
            "steps")
    nb.md("## 4 · One changed setting rebuilds one partition",
          "The `layer_vote` partition is given a different settings record (as if its default changed); the other "
          "two keep their fingerprints and are reused.")
    nb.code("changed = [dataclasses.replace(p, settings=p.settings[:-1] + ', \"pilot_change\": 1}') if p.strategy == 'layer_vote' else p for p in parts]",
            "built = []",
            "def recording(p):",
            "    built.append(p.strategy)",
            "    return P.build_partition(dataclasses.replace(p, settings=[q for q in parts if q.strategy == p.strategy][0].settings))",
            "fourth = P.run(changed, work, builder=recording)",
            "assert built == ['layer_vote'], built",
            "fourth")
    nb.md("## 5 · Cost per partition, and agreement with the shipped record",
          "Peak memory is the process high-water mark (VmHWM) after each partition, so it includes earlier work in "
          "the same process; it is not an isolated per-partition peak.")
    nb.code("cost = P.report(work)",
            "cost['peak_rss_MiB'] = (cost.peak_rss_bytes / 2**20).round(1); cost['KiB'] = (cost.bytes / 1024).round(1)",
            "table = P.assemble(work, parts)",
            "shipped = T.shipped('Tg')",
            "rows = []",
            "for m in METHODS:",
            "    a = table[(table.strategy.astype(str) == m) & (table['mode'].astype(str) == 'together')]",
            "    s = shipped[(shipped.strategy.astype(str) == m) & (shipped.target.astype(str) == TARGET) & (shipped['mode'].astype(str) == 'together')]",
            "    s = s[s.setting_key.astype(str) == s.setting_key.astype(str).iloc[0]] if len(s) else s",
            "    j = a.merge(s, on='gene_id', suffixes=('', '_shipped'))",
            "    both_missing = j.prediction.isna() & j.prediction_shipped.isna()  # an abstention is a missing call on both sides",
            "    same = ((j.prediction.astype(str) == j.prediction_shipped.astype(str)) | both_missing).mean() if len(j) else float('nan')",
            "    rows.append(dict(method=m, genes_now=a.gene_id.nunique(), genes_shipped=s.gene_id.nunique(), both=len(j),",
            "                     same_call=round(float(same), 3), accuracy_now=round(float(a.correct.fillna(False).mean()), 3),",
            "                     accuracy_shipped=round(float(s.correct.fillna(False).mean()), 3) if len(s) else None))",
            "agree = pd.DataFrame(rows)",
            "cost.to_csv(os.path.join(OUT, 'cost.tsv'), sep='\\t', index=False); agree.to_csv(os.path.join(OUT, 'agreement.tsv'), sep='\\t', index=False)",
            "storage = sum(os.path.getsize(os.path.join(dp, f)) for dp, _d, fs in os.walk(work) for f in fs)",
            "json.dump(dict(target=TARGET, methods=METHODS, steps={k: dict(v) for k, v in dict(first=first, second=second, third=third, fourth=fourth).items()},",
            "               storage_bytes=storage, total_seconds=float(cost.seconds.sum())), open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "print('storage', storage, 'bytes; total build seconds', round(float(cost.seconds.sum()), 1))",
            "print(cost[['partition', 'status', 'rows', 'seconds', 'peak_rss_MiB', 'KiB']].to_string(index=False))",
            "agree")
    nb.write(NB)


if __name__ == "__main__":
    main()
