# 64.14 · Test map and module recovery

Status: DONE 2026-11-07. Held-out recovery, noise and fragmentation measured and shown on the Map module cards.

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

**Fragmentation and reach (same day, `fragmentation.tsv`):** on each label's best non-circular map, the best category
is 81-97 % clustered (noise is not the problem: maps leave 6-20 % unclustered, one 44 %) but covering 80 % of its
clustered members takes **3-9 clusters** -- the classes are split across several clusters, which is why their F1
stays low. Open: class-level confusion and these controls shown on the result cards (today they are tables).

**On the cards (same day):** the Map module **Labels** table now shows, per label, "groups for 80 %" (how many groups
hold 80 % of its grouped genes) and "in no group" (its unclustered share), beside the best group's precision, recall
and F1 (`map_reading.LabelRow.groups_for_80`, `.unclustered`; test in `tests/test_map_reading.py`). The **Groups**
table shows each group's majority label with its precision -- the confusion seen from the group side. All done-when
points are met: label-free selection, held-out recovery against chance, noise and fragmentation, geometric and
biological metrics kept apart.
