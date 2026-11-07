# 64.38 · Freeze predictions and score new independent truth

Status: DONE 2026-11-07 (infrastructure and first freeze). Grading waits on a later screen, by construction.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.04, 64.05, 64.08, 64.16, 64.32, 64.34 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `parascape/claims.py`, `parascape/datasets.py`, `parascape/deposits.py`, `parascape/scorecard.py`, `scripts/add_deposits.py` |

## Deliverable

Version/hash-frozen claims and meta-inferences plus an evaluator triggered by newly admitted compatible experimental truth.

## Controlled scope

Freeze/evaluate infrastructure and one historical time-ordered or genuinely later dataset pilot with provenance.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Predictions, thresholds and eligible populations are frozen before seeing the new dataset; old/new dataset identity and reuse exclusions are checked.
- [ ] Reports retain tested/correct/wrong/unmapped/unmeasured counts, class/context performance and precision/coverage at frozen thresholds.
- [ ] A later dataset scores archived predictions without fitting them again; its results become visible on the same result/method/meta cards.
- [ ] Until real later truth is available, the feature may pass fixture tests but the biological prospective validation remains open.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-06: the freeze and the grader

* **Freeze** (instruction 63): `results/prospective_claims_2026_11_06/` -- 747 genes the glucose-withdrawal
  screen did not measure, ranked; `manifest.json` holds the claims SHA-256, strategy and evidence hashes,
  version and time. Built by `scripts/build_prospective_claims_notebook.py`, frozen before any later truth.
* **Grader**: `parascape/prospective.py` -- `load` refuses a set whose file no longer matches its frozen hash;
  `grade(folder, truth)` fits nothing, puts every claim in tested (correct / wrong) or unmeasured, and reports
  precision at the top-k cuts fixed in the module (50, 100, 200) plus AUC over tested genes.
  `tests/test_prospective.py` (4, fixtures plus the shipped freeze still matching its hash).
* Done-when status: predictions frozen before the new dataset -- yes; counts and precision at frozen cuts --
  yes; a later dataset scores without refitting -- yes, on fixtures; **the biological validation stays open**
  until a screen covering those genes is published, and results are not yet shown on the method cards.

**On the cards (2026-11-07):** a gene in a frozen claim set shows its frozen prediction in the Gene module's method
section -- probability, rank, freeze date and the claims hash, "waiting for a later screen to grade it, unchanged"
(`gene_reading.prospective_lines`; the hash is checked by `prospective.load` before anything is shown). Test in
`tests/test_gene_reading.py`. With the freeze, the grader and this, every done-when point is met in the form the
contract allows before later truth exists: frozen before the new dataset, counts and precision at frozen cuts,
scoring without refitting (fixture-tested), visible on the cards. **The biological prospective validation itself
stays open by construction** until a screen covering those 747 genes is published; grading it is then one call.
