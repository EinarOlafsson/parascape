# Set fingerprints for the host-facing compartment ladder (instruction 82)

Notebook: `notebooks/set_fingerprints_2026_11_06.ipynb` (executed). Builder: `scripts/build_set_fingerprints_notebook.py`. Ladder from `results/host_facing_truth_2026_11_04` only; place columns excluded.

## Held-out recovery (AUROC, leave-one-out)

| rung | genes | fingerprint | random | dense_granule | verdict |
|---|---|---|---|---|---|
| 1 injected at the moving junction | 3 | nan | nan | nan | too few to score |
| 3 PVM, host-facing | 45 | 0.944 | 0.525 | 0.828 | above the best baseline |
| 4 vacuole-resident (hard negative) | 32 | 0.947 | 0.5 | 0.831 | above the best baseline |
| 5 exported to the host cytosol | 4 | nan | nan | nan | too few to score |
| 6 host nucleus | 9 | 0.976 | 0.518 | 0.808 | above the best baseline |
| 7 host organelle or plasma membrane | 5 | 0.98 | 0.559 | 0.807 | above the best baseline |

## Confusion between scored rungs (row: true, column: assigned)

| index | 3 PVM, host-facing | 4 vacuole-resident (hard negative) | 6 host nucleus | 7 host organelle or plasma membrane |
|---|---|---|---|---|
| 3 PVM, host-facing | 10 | 18 | 9 | 8 |
| 4 vacuole-resident (hard negative) | 5 | 18 | 4 | 5 |
| 6 host nucleus | 2 | 2 | 5 | 0 |
| 7 host organelle or plasma membrane | 1 | 0 | 0 | 4 |

Files: `ladder.tsv`, `recovery.tsv`, `confusion.tsv`, `candidates.tsv` (top 200, a ranking only).
