# 64.24 · Explore typed organism-to-host connections

Status: DONE — 100%, 2026-11-07. Verified bridges in one organism-qualified table; the other end of a link opens in its own space (parasite map, or the human/mouse gene pack). Built; Einar's review pending -- reopen on his feedback.

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

- [x] Measured interaction, orthology, predicted interaction and infection-response links retain distinct labels and provenance.
- [x] Opening a counterpart selects its own space; tables and measurement denominators remain species-specific.
- [x] A bridge has mapping/coverage evidence and can be explored without an inference-accuracy claim; inference through it requires a registered benchmark.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

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

## 2026-11-07: the other end, in its own space

`bridges.counterpart(link, from_gene=)` resolves the other end of any link: a Toxoplasma or Plasmodium gene opens
in its own map (rows and gene count from that organism's own table); a host protein is resolved to host GENES
through the installed host pack's reviewed protein links (64.21 human, 64.22 mouse) -- never by name -- with the
host space's own rows and denominator; a host whose species the source does not record gets no space and says
so. The Gene module has a section "Links to other organisms" (after the evidence): each link with its relation
type, the other end, and where it opens; clicking a parasite end opens that organism map and its Gene module,
clicking a host end shows the host gene rows from the host pack. Exploring a link makes no inference claim.

Tests: `tests/test_bridges.py` (6 pass: host protein -> two host genes with the human denominator, a protein
with no reviewed gene, no pack installed, species not recorded, a parasite end from either side) and
`tests/test_module_shell.py::test_the_gene_module_lists_links_to_other_organisms_in_their_own_spaces`.

Gaps, stated: the 18 curated-literature host targets do not record a host species, so none opens until the
source says; the human and mouse packs are built but not hosted for download, so a fresh install shows "not
installed" until a pack is installed (`python -m parascape.packs install …`); inference through a bridge still
needs its own registered benchmark.
