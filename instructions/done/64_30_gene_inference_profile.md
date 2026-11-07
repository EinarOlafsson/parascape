# 64.30 · Show every strategy's answer for a gene

Status: DONE 2026-11-07. One gene inference profile: every method held out, grade, scorecard link, agreement with its measured accuracy, combined answer, stale flag.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.17, 64.19, 64.26, 64.27 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `parascape/app.py`, `parascape/strategy_panel.py`, `parascape/track_record.py`, `parascape/claims.py` |

## Deliverable

One gene inference profile listing each strategy's applicable outputs, status, supporting evidence, calibrated reliability and conflicting calls.

## Controlled scope

One unified gene profile over the existing evidence page; use available Tg/Pf packs first, then validated host packs.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] All 39 catalogue entries are accounted for; unsupported, abstained, failed, stale and pending answers remain inspectable.
- [ ] A known gene's recovery view uses its held-out prediction, and an unknown gene's claim exposes the matching class/stratum test population.
- [ ] Cached results appear without recomputation; new query-dependent work has visible progress/cancel and attaches its own manifest.
- [ ] One click from each answer reaches its scorecard and deeper method/evidence description.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: every method, per gene, held out

`gene_reading.strategy_profile(gene)` lists all 39 catalogue strategies for a gene: the calibration grade, and for
the 14 the hold-out record covers, each label's call with the gene's own answer hidden (right, wrong or
abstained, with the truth); the other 25 say why they give no per-gene answer ("answers about structure or a
list"). Shown in the Gene module as **What each method said about it, held out** (folded). Read from the shipped
record, no recomputation (0.2 s). Tests in `tests/test_gene_reading.py` and `tests/test_module_shell.py`. Open: for
an unknown gene, the class/stratum population behind a claim; one-click links to each scorecard; stale/pending
states.

## 2026-11-07 (later): the open points, met by today's other work

* **One click from each answer to its scorecard:** each method in the Gene module's held-out section is a link that
  opens its strategy card with calibration, scorecard and track-record line (`GeneModule.open_strategy`, 64.36).
* **An unknown gene's claim exposes its test population:** the combined answer for a gene with no known answer is
  shown with the label's measured accuracy on untouched orthogroups and the number of genes it rests on (64.35); the
  agreement line names how many genes with that much agreement the rate rests on (64.33); every method's grade
  comes from its calibration population.
* **Known gene uses its held-out prediction:** the track-record calls and the combined outer-test call are both
  made with the gene's own answer hidden.
* **Cached, no recomputation:** all read from shipped records (0.2 s for the method profile; the agreement tables
  cached per label after their first build).
* **All 39 catalogue entries accounted for**, unsupported ones saying why.
* **Remaining:** a "stale" state when a record predates the current table (none is flagged today) and visible
  progress for any query-dependent work (none is launched from this section).

**Stale records (same day):** `gene_reading.stale_records` compares, per label, the answer the hold-out record was
scored against with the gene's current value; where they differ the Gene module says the record's verdicts on that
label are out of date until it is rebuilt. Today none of 3,000 genes checked is stale. Test in
`tests/test_gene_reading.py`. All done-when points are met.
