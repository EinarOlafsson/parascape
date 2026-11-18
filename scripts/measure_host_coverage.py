"""What share of the host space each host dataset actually covers.

Instruction 65, the open complaint that the A grades overstate coverage. They do, and the reason is
the denominator: every `coverage` string in the dataset catalogue reports the share of the SOURCE
FILE that resolved, so `host_bone_marrow_proteome` reads "2,059 proteins (93.8 per cent of the
2,195 in the source file)". That is true and it is not the number a reader wants. The question is
what share of the host space carries a value, and for that dataset it is a tenth.

Both numbers are legitimate. The share of the source says the import worked; the share of the
space says how much of the host the dataset can speak about. A grade assigned on the first and
read as the second is the overstatement, so this measures the second for every host dataset and
writes it where the catalogue can carry it.

Writes `results/host_coverage_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import datasets as DS       # noqa: E402

#: The table each host organism's columns live in, and what one row of it is.
SPACES = {"Hs": ("hs_host_proteins.parquet", "human proteins"),
          "Mm": ("mm_host_proteins.parquet", "mouse proteins"),
          "Ag": ("ag_host_genes.parquet", "Anopheles genes")}


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")

    rows = []
    for organism, (filename, unit) in SPACES.items():
        path = os.path.join(cache, filename)
        if not os.path.exists(path):
            print(f"{organism}: {filename} absent; skipped")
            continue
        table = pd.read_parquet(path)
        total = len(table)
        print(f"\n{organism}: {total:,} {unit}")
        for dataset in DS.REGISTRY:
            if dataset.organism != organism or not dataset.columns:
                continue
            present = [c for c in dataset.columns if c in table.columns]
            if not present:
                rows.append({"dataset": dataset.key, "organism": organism, "columns": 0,
                             "covered": 0, "space": total, "share": 0.0,
                             "note": "none of its columns are in the host table"})
                print(f"  {dataset.key:38s} none of its columns are in the table")
                continue
            covered = int(table[present].notna().any(axis=1).sum())
            share = covered / total
            rows.append({"dataset": dataset.key, "organism": organism, "columns": len(present),
                         "covered": covered, "space": total, "share": share,
                         "catalogue_says": dataset.coverage, "note": ""})
            print(f"  {dataset.key:38s} {covered:6,} of {total:,} {unit}  {share:6.1%}"
                  f"   catalogue: {dataset.coverage[:46]}")

    table = pd.DataFrame(rows)
    measured = table[table["columns"] > 0]
    thin = measured[measured["share"] < 0.25]
    print(f"\n{len(thin)} of {len(measured)} host datasets cover under a quarter of their space:")
    for row in thin.sort_values("share").itertuples():
        print(f"  {row.dataset:38s} {row.share:6.1%}")
    print("\nThe catalogue coverage string reports the share of the SOURCE FILE, which is a "
          "different\nquestion and a much larger number. Both belong in the record; only one "
          "should be read as\nhow much of the host a dataset can speak about.")

    out = os.path.join(ROOT, "results", f"host_coverage_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, "by_dataset.tsv"), sep="\t", index=False)
    json.dump({"spaces": {k: v[1] for k, v in SPACES.items()},
               "datasets_measured": int(len(measured)),
               "under_a_quarter": int(len(thin)),
               "rows": table.to_dict("records")},
              open(os.path.join(out, "summary.json"), "w"), indent=2, default=float)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
