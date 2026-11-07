# 64.12 · Test gene rankings and set retrieval

Status: OPEN — 15%, 2026-11-07. First adapter: seed expansion, held-out member recovery on six localization niches; the other ranking strategies remain.

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
