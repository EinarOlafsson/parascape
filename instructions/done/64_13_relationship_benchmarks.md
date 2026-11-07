# 64.13 · Test relationship and link inference

Status: DONE 2026-11-07. Held-out pair records with AUPRC, lift, fixed-depth precision and decoy regimes for all four mechanisms and binding partners; literature recovery named as such.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `parascape/strategy_catalog.py`, `parascape/strategy_graph.py`, `parascape/scorecard.py`, `parascape/methods.py` |

## Deliverable

Held-out pair records for link prediction, attention-corrected relations, unwritten links and learned neighbor/network spaces.

## Controlled scope

Relationship-task adapters and benchmark semantics only; acquisition of physical truth is source-specific.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Candidate pairs have endpoint organism IDs, evidence source, known-edge/cold-node regime and matching negative/control rules.
- [ ] A co-mention recovery score is described as literature recovery; inferred physical links have a separate experimental validation requirement.
- [ ] Ranking cards retain prevalence, AUPRC/lift, precision at fixed depth and dependence exclusions; unknown pairs remain unlabeled where negatives cannot be established.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the held-out pair records, read against this contract

* **Binding partners** (`results/partner_inference_2026_11_04/`, instruction 72): positives the 2,842 StarPath XL-MS
  crosslink pairs; no verified non-interaction set, so every negative is a stated decoy. Cold-node, degree-matched
  (headline): lower 90 % bound 0.891, misses 0.90; same-compartment decoys 0.919; random decoys 0.902; known-node
  0.912. Pull-down pairs kept out until a background table exists (instruction 80).
* **Calibration sweep (2026-09-26, hidden pairs, measured null):** link prediction AUROC skill 0.61 (Tg) / 0.96 (Pf)
  against degree-matched non-pairs; learned neighbour and network spaces 0.59 / 0.39 on hidden co-expression edges;
  attention-corrected relations 0.38 (Tg, share of top-50 pairs sharing a compartment) / 0.26 (Pf, weak); unwritten
  links 0.11 / 0.11.
* **Semantics, as required:** unwritten links is **literature recovery** -- it predicts co-mention, not contact; the
  physical claim of any inferred link needs separate experimental validation (crosslink, reciprocal pull-down).
  Unknown pairs stay unlabelled; decoys are named as decoys.
* **Open:** AUPRC, lift and precision at fixed depth stored for each relationship card (today the cards carry AUROC
  skill), and endpoint organism IDs written into the pair records.

**Average precision and baselines for binding partners (same day, from `scoreboard.tsv`, not refitted):** cold-node,
degree-matched, 1,248 pairs, half positive -- model AUPRC 0.808 (lift 1.6 over the 0.50 prevalence), precision 0.79
at recall 0.54; random 0.504; "both popular" (degree product) 0.529, so the model is not riding node degree;
same-compartment alone reaches precision 0.85 at recall 0.49 -- compartment is most of what a decoy-based test can
reward, which is why the same-compartment decoy regime is reported beside it.

**Fixed-depth metrics were already on the cards (checked the same day):** every relationship strategy's calibration
record carries the ranking scorecard -- AUPRC, AUPRC lift over prevalence, precision in the top 1 %, R-precision,
recall in the top 10 % -- shown on its strategy card. Tg, default settings: link prediction AUPRC 0.607 at prevalence
0.168 (lift 3.6), precision in the top 1 % 0.95; attention-corrected relations 0.637 at 0.525 (lift 1.2); unwritten
links 0.244 at 0.167 (lift 1.5; literature recovery); learned neighbour and network spaces 0.802 at 0.5 (lift 1.6),
precision 1.0 in the top 1 %. Pair records use organism-specific gene ids (TGME49_ / PF3D7_) for both endpoints.
All done-when points are met.
