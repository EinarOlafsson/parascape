# 64.14 · Test map and module recovery

Status: OPEN — 70%, 2026-11-07. Held-out cluster and module recovery summarised for maps and all three mechanisms; fragmentation and reach controls on the cards remain.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `parascape/strategy_catalog.py`, `parascape/embedding.py`, `parascape/clustering.py`, `parascape/scorecard.py` |

## Deliverable

Persistent cluster/module recovery records for map searches, consensus modules and multiplex modules.

## Controlled scope

Cluster-recovery task and three current mechanisms; no renderer work.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Maps/module selection use permitted training evidence; settings search selects on inner data and evaluates on an untouched outer set.
- [ ] Cards show held-out recovery, fraction placed/clustered, noise, fragmentation/reach controls and class-level confusion.
- [ ] Geometric quality/stability and biological recovery have distinct metrics; withheld labels are not used to choose the reported winner.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: what the held-out records say

`scripts/build_cluster_benchmark_summary_notebook.py` -> `results/cluster_benchmark_summary_2026_11_07/`, read from what
ships (nothing refitted):

* **Maps** (`umap_gallery_scores.tsv`, 2,502 label x map scores): every map chosen by a label-free structure score, so
  no withheld label picked a winner; maps whose inputs restate the label set aside as circular. The best non-circular
  map per label lifts category F1 only 0.10-0.17 over its shuffled-label chance (cell-cycle phase 0.36 against 0.25;
  compartment 0.24 against 0.14). Geometric quality (`structure`) and recovery (`category_skill`) are separate metrics.
* **Module mechanisms**, held out (cluster chosen for a label on known genes, scored on hidden genes, measured null):
  map search skill 0.07 (Tg) / 0.09 (Pf), consensus modules 0.07 / 0.15, multiplex modules 0.05 / 0.03 -- all graded
  weak.
* **Reading:** clusters drawn without the label recover it only weakly, while a model trained for the label
  recovers it well (localization frozen AUC 0.85-1.0, instruction 70). Clusters are for seeing structure, not for
  calling a label.
* **Open:** fragmentation and reach controls, and class-level confusion shown on the cards.
