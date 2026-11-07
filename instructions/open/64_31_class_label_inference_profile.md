# 64.31 · Compare inferences for a label or protein class

Status: OPEN — 35%, 2026-11-07. Per-class and per-label recovery pages exist (track record, item 62); inferred-member lists and numeric-context cards remain.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.17, 64.20, 64.26, 64.27, 64.30 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `parascape/track_record.py`, `parascape/strategy_panel.py`, `parascape/app.py`, `parascape/guided.py` |

## Deliverable

Class/label inference profiles listing measured members, predicted candidates, all applicable methods and per-class/target performance.

## Controlled scope

Class and label profiles over declared annotation schemas; no silent new ontology mappings.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Known membership and inferred membership are distinguished; class-specific precision/recall, priors/lift and missed/confused classes remain accessible.
- [ ] Numeric contexts use their own quantities and value-task cards rather than being collapsed into categorical calls.
- [ ] Query-dependent set results are tied to the exact input members/context; a gene-set scorecard cannot be mistaken for global method accuracy.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: what already exists, read against this contract

* **Per-class performance and confusions -- exist** (instruction 62): `track_record.class_html` shows, for one class
  of one label, how often each strategy recovers it with members hidden, with Wilson intervals, abstentions, and
  the classes it is mistaken for; `track_record.target_html` ranks every class of a label by recoverability and
  names classes no strategy recovers. Reached from the Gene module and strategy cards.
* **Known against inferred membership -- partly:** the Map module "Groups" section lists, per group, the known
  label share and the unknown genes that become candidates (`map_reading.learn_from_groups`, candidates only in
  groups at least 60 % pure); kept claims are shown as kept, never as facts (instruction 78). A class page that
  lists measured members beside predicted candidates with their per-class precision is still to build.
* **Open:** numeric-context cards (value tasks), and tying query-dependent set results to their exact input
  members on the class page.
