# 64.35 · Expose meta-inference scorecards at every level

Status: DONE 2026-11-07. Gene, class and label views of the combined model, admitted or refused, each with its test.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.17, 64.30, 64.31, 64.34 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `parascape/scorecard.py`, `parascape/strategy_card.py`, `parascape/claims.py`, `parascape/app.py` |

## Deliverable

Gene-result, class, label, meta-method and organism views of the meta-model's validation and base-method comparisons.

## Controlled scope

UI/readers over item 34's tested artifacts; no extra model training.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Every displayed meta confidence links to the applicable test population, calibration and dependence report.
- [ ] Paired base-versus-meta improvement, disagreement cases and abstention/coverage tradeoffs are reachable from the compact card.
- [ ] Known-gene outer-test outcomes and unknown-gene conditional reliability are labeled correctly; global accuracy does not masquerade as class-specific or gene-specific accuracy.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the gene view

The Gene module's method section now ends with **Combined**: for each label with a meta-inference admission test
(64.34), either "no combined answer -- not admitted" with its test numbers, or the combined call for this gene made in
the untouched outer test with the gene's answer hidden, its probability, whether it was right, and the test it passed
("combined 0.919 against the best single method 0.902 on untouched orthogroups") -- so the known-gene outcome is
labelled as such, and a global accuracy is never shown as this gene's (`gene_reading.meta_lines`, `META_RESULTS`).
Open: combined calls for genes with no known answer (a deployment fit, with the label's measured accuracy as its
conditional reliability), and class and label views.

**Unknown genes (same day):** for an admitted label the meta-model is refitted on the base methods' outer-test calls
on all known genes and applied to their calls on the genes with no known answer (each made with every known label
visible): 4,313 unassigned genes get a combined compartment call (`meta_calls_unknown.tsv`, run from the notebook's
last cell). The Gene module shows it as "predicted X by the combined model (probability p); on this label the
combined answer was right 54.3 % of the time on untouched orthogroups, which is its reliability, not a certainty
for this gene" -- the unknown-gene conditional reliability labelled as such. Open: class and label views of the
meta-model.

**Class and label views (same day):** a label page (`track_record.target_html`, reached from the gene and strategy
cards) now ends with the combined model's verdict and its recall per class from the untouched outer test
(`track_record.meta_class_html`) -- e.g. dtm_class: GLOB 98.5 %, TM 92.7 %, SP 55.2 %, SP+TM 30.8 % -- so a global
accuracy is never read as a class's. Test in `tests/test_track_record.py`. All done-when points are met: each combined
answer names its test population and verdict, paired base-versus-meta comparison and refusals are shown, known-gene
outcomes and unknown-gene reliability are labelled apart, and class-specific recall is separate from the overall
number.
