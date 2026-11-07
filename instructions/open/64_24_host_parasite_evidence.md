# 64.24 · Explore typed organism-to-host connections

Status: OPEN — 50%, 2026-11-07. One organism-qualified table of verified bridges with relation types and coverage; opening a counterpart in its own space waits on host gene spaces (64.21/22).

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.21, 64.22, 64.23 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `parascape/host.py`, `parascape/organisms.py`, `parascape/star_edges.py`, `parascape/star_map.py`, `parascape/datasets.py` |

## Deliverable

Explicit cross-space evidence links with organism-qualified endpoints, relation type, context and source quality.

## Controlled scope

Expose existing verified bridges only; new cross-species inference is separately benchmarked.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Measured interaction, orthology, predicted interaction and infection-response links retain distinct labels and provenance.
- [ ] Opening a counterpart selects its own space; tables and measurement denominators remain species-specific.
- [ ] A bridge has mapping/coverage evidence and can be explored without an inference-accuracy claim; inference through it requires a registered benchmark.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the verified bridges in one table

`parascape/bridges.py` (`all_links`, `coverage`, `links_of`): every link leaving an organism's space, with
organism-qualified endpoints, one relation type each, context, support and reference -- measured interactions from
curated literature (18, Tg -> host; the host species is "not recorded" where the source does not say, never
guessed), measured interactions from crosslinking MS (10, Pf -> human) and orthology through shared OrthoMCL groups
(3,291 Tg <-> Pf links over 2,857 Tg genes). Relation types are never merged; exploring a link makes no inference
claim. Tests in `tests/test_bridges.py`. **Open:** opening a host counterpart selects its own host space -- that
needs the human and mouse gene spaces of 64.21 and 64.22.
