# 78 · A user's own layer: keeping an inference without contaminating the map
**Status: done 2026-11-06 (moved to done).** `parascape/kept_claims.py`: a separate store in the user
cache folder (`kept_claims.jsonl`), never the release; every claim records value, probability, uncertainty,
method and version, question, accepted score, whether it cleared the 90 % bar, the data release and the date;
refused without its method. Shown in the Gene module under "What can be inferred" as **Kept by you** (never
among the facts) with a *forget* link, and **Now measured -- agrees / DISAGREES** when a release measures the
cell (`check_against`). Kept from the Map-neighbour suggestion (*keep this*, with its hit rate) and from the
Infer module (*Keep the chosen answers*); exported with provenance (Infer ▸ Save the answer). Gate tests in
`tests/test_kept_claims.py` (never facts, no training code imports the store, survives a release upgrade,
superseded loudly, revert per claim and per method version).

Requested 2026-11-05. The user-facing half of inference; specified in instruction 74 and
lifted here because it is infrastructure, not layout.

## The request (Einar)
"There should be mechanisms that connect the one protein and labels, where if inference is done the user can
choose to add it to their own ParaScape dataset and then it should turn up in `Predicted`."

## The plan
1. **A separate store**, never the shipped map: the application writes nothing to the release.
2. **A kept claim records**: the value, its probability and uncertainty, the strategy and version, the question,
   the score it was accepted at, whether it cleared the 90 % frozen-test bar, the data release it was computed
   against, and the date. Reproducible from its record or not kept.
3. **Where it appears**: `Gene -> Predicted`, visibly distinct from measured facts and from predictions nobody
   accepted. **Never in `Facts`.** `Label -> Fill` is the bulk form of the same action.
4. **Rules.** A kept inference is never training data -- now enforceable because the layer is separate. It can
   be reverted individually or by strategy version. It exports with its provenance so it can be shared or cited.
   A later release that *measures* a kept claim supersedes it **and says so**, rather than silently agreeing:
   without that, a kept guess that happens to be right becomes indistinguishable from evidence.

## Gate
A test that a kept claim never reaches a training set, never appears in `Facts`, survives a release upgrade with
its provenance intact, and is superseded loudly when the same cell becomes measured.
