# Set fingerprints for the host-facing compartment ladder (instruction 82)

Notebook: `notebooks/set_fingerprints_2026_11_08.ipynb` (executed). Builder: `scripts/build_set_fingerprints_notebook.py`. Ladder from `results/host_facing_truth_2026_11_04` plus the literature rungs (`literature_changes.tsv` lists every gene added or moved); place columns excluded.

## Held-out recovery (AUROC, leave-one-out)

| rung | genes | fingerprint | random | dense_granule | verdict |
|---|---|---|---|---|---|
| 1 injected at the moving junction | 5 | 0.954 | 0.523 | 0.778 | above the best baseline |
| 2 PVM, vacuole-facing | 14 | 0.975 | 0.513 | 0.826 | above the best baseline |
| 3 PVM, host-facing | 39 | 0.941 | 0.487 | 0.826 | above the best baseline |
| 4 vacuole-resident (hard negative) | 26 | 0.939 | 0.482 | 0.834 | above the best baseline |
| 5 exported to the host cytosol | 4 | nan | nan | nan | too few to score |
| 6 host nucleus | 9 | 0.976 | 0.518 | 0.816 | above the best baseline |
| 7 host organelle or plasma membrane | 5 | 0.98 | 0.483 | 0.818 | above the best baseline |

## Confusion between scored rungs (row: true, column: assigned)

| index | 1 injected at the moving junction | 2 PVM, vacuole-facing | 3 PVM, host-facing | 4 vacuole-resident (hard negative) | 6 host nucleus | 7 host organelle or plasma membrane |
|---|---|---|---|---|---|---|
| 1 injected at the moving junction | 3 | 0 | 0 | 1 | 0 | 1 |
| 2 PVM, vacuole-facing | 0 | 3 | 1 | 3 | 5 | 2 |
| 3 PVM, host-facing | 2 | 8 | 9 | 8 | 6 | 6 |
| 4 vacuole-resident (hard negative) | 1 | 6 | 2 | 12 | 1 | 4 |
| 6 host nucleus | 0 | 2 | 2 | 0 | 5 | 0 |
| 7 host organelle or plasma membrane | 1 | 0 | 1 | 0 | 0 | 3 |

Files: `ladder.tsv`, `recovery.tsv`, `confusion.tsv`, `candidates.tsv` (top 200, a ranking only).
