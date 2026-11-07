# 64.12 · Test gene rankings and set retrieval

Status: DONE 2026-11-07. All seven ranking tasks have held-out records; full-ranking average precision over five splits for the fast strategy.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `parascape/strategy_catalog.py`, `parascape/scorecard.py`, `parascape/track_record.py` |

## Deliverable

Persistent held-out member/ranking records for gene-set search, recoverability, outlier ranking, PU learning, enrichment, seed expansion and paralog profiles.

## Controlled scope

Gene/list ranking and retrieval tasks only; one native adapter per commit.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Seeds and validation members are disjoint at the declared group scope; scores include candidate-universe and query-set identities.
- [ ] Cards expose precision/recall at fixed depths, AUPRC and prevalence lift where defined; thresholds are fixed on tuning data.
- [ ] Corrupted-label detection is marked as a synthetic task; PU/background and paralog-label proxies never claim verified biological negatives or adjudicated misannotation.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the first adapter, seed expansion

`results/seed_expansion_benchmark_2026_11_07/README.md`. Seeds and hidden members disjoint by orthogroup; depths
fixed before the run; hyperLOPIT compartment excluded. Recall at 100: rhoptries 0.81, micronemes 0.86, dense
granules 0.57, IMC 0.53, mitochondrion 0.34, apicoplast 0.29 (lift 24-70 over the prevalence). Non-markers are not
verified negatives (stated in the notebook). Next adapters: positive-unlabelled learning, gene-set search,
recoverability, outlier ranking, enrichment, paralog profiles.

**Second adapter (same day): positive-unlabelled learning** on the identical splits -- better than seed expansion for
dense granules (recall at 100: 0.76 against 0.57) and the apicoplast (0.46 against 0.29), worse for the invasion
organelles. Table in the results README.

**Third adapter (same day): set enrichment** -- never the best of the three (recall at 100 0.23-0.69); its profile
includes publication counts, so it partly ranks by study effort.

**Fourth adapter (same day): outlier ranking, marked SYNTHETIC** (`scripts/build_label_outlier_benchmark_notebook.py`,
`results/label_outlier_benchmark_2026_11_07/`): 5 % of a label's values swapped to another class in a copy of the
table, label-outlier ranking run on the corrupted label. Precision among its top 50: dtm_class 0.74, compartment
0.34, cell-cycle phase 0.32 (base rate 0.05); cell-cycle phase recall 0.51 at 100. It finds planted mistakes;
it adjudicates no real annotation. Placeholder values (`unassigned`, `unknown`) are not labels and are never
planted (`compartment`: 3,827 real labels of 8,140).

**Fifth adapter (same day): gene-set search** -- a set-retrieval task (one cluster of about 16 candidates), scored as
hidden members inside it: 2-3 for four niches, none for micronemes and apicoplast. Remaining: recoverability atlas
(category-level, not per gene) and paralog profiles (pairs).

**Recoverability and paralog profiles (same day): their held-out records already exist** in the shipped calibration
(`strategy_calibration.json`, sweep of 2026-09-26: settings chosen on seeds 1-3, reported on seeds 4-5, against a
measured null). Recoverability atlas: hidden-gene AUROC of the categories it ranks in its top half, skill 0.50 (Tg,
25 conclusive runs) and 0.69 (Pf), graded reliable. Paralog divergence: AUROC of profile divergence for paralogs
of different compartment, skill 0.16 (Tg, reliable) and 0.11 (Pf, weak) -- a paralog-label proxy, never a
verified biological difference, as this contract requires it to be called.

**Where this leaves the contract:** all seven named ranking tasks have held-out records (five built today on
orthogroup-disjoint seed / hidden splits at depths fixed in advance, outlier ranking marked synthetic, two from the
calibration sweep). Not yet: AUPRC over the full ranking (the rankings stop at 2,000 genes) and repeats over several
splits for intervals.

**Full ranking and repeats (same day):** positive-unlabelled learning over five orthogroup splits, ranking the whole
genome: average precision rhoptries 0.49 (0.41-0.55) down to mitochondrion 0.08 (0.04-0.12), lift 20-160 over a
random ranking (`pu_repeats_summary.tsv`). The slow network walk (seed expansion, about 2.5 h per niche set) keeps its
single split. With this the done-when points are met: disjoint seed / hidden groups, depths fixed beforehand,
recall, lift and average precision, the synthetic and proxy tasks named as such.
