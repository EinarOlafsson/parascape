# 64.21 · Build the human gene information space

Status: DONE — 100%, 2026-11-07. Foundation, documentary terms, HPA annotations, reviewed protein links and a distributable pack that installs and opens on its own. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.02, 64.03, 64.04, 64.05, 64.08 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `parascape/organisms.py`, `parascape/host.py`, `parascape/spaces/`, `parascape/packs.py`, `parascape/datasets.py` |

## Deliverable

One reproducible Hs gene-space builder using already acquired, verified reference and selected measurement sources.

## Controlled scope

One host species and a named pilot source set selected from the inventory; use instruction 53's builders/pack work.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Declare canonical gene/protein mappings and per-column species, units, licenses, provenance and targets; quantify one-to-many/lost joins.
- [x] A notebook reproduces the build and per-source sanity checks; unique IDs, non-loss and graph-order checks pass.
- [x] A validated distributable human pack opens independently; protein reference rows are not relabeled as genes.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## Frozen first partition: canonical source genes and explicit protein links

2026-10-08, before full source parsing. Reuse the verified GTEx v10 processed
median-TPM source and verified UniProt 2026_03 reviewed human entry/idmapping
snapshots already acquired in the host-expression review. Retain all 59,033
source records and both selected numeric fields. Build canonical unversioned
Ensembl gene rows only from ordinary source IDs; retain 45 `_PAR_Y` records as
separate qualified source evidence, never collapse them into X-gene values.
The resulting 58,988-gene universe is this source's annotation scope, not a
claim of genome-wide completeness. Preserve versioned source IDs, names and
source row positions; zero median TPM is an observed summary, not missingness.

Use round-trip float parsing and verify each retained measurement against its
original source text. Map explicit reviewed UniProt Ensembl cross-reference
pairs without resolving one-to-many relations by row order, symbol guesses or
averaging. Gene/protein mappings do not authorize automatic protein-measurement
projection. Quantify unmapped, multi-protein genes, multi-gene proteins and
qualified records. Declare source hashes, contexts, median-TPM units, curation
versus measured-evidence grades and redistribution status/gaps per column.
Cultured adult GTEx fibroblasts are not independently validated HFF measurements.

Execute a reproducible builder/replay notebook in a new
`results/human_gene_space_foundation_2026_10_08/` packet. This first partition
does not register/promote an installed host gene space or claim a distributable
pack before reference, license, non-loss/graph-order/opening and benchmark
gates pass. Keep installed protein references and parasite data unchanged;
complete pack/opening and additional selected sources remain controlled follow-up
partitions of this item. GTEx v11 stays a pending source-comparison candidate
under 65.04; version releases do not reset publication citation age.

Foundation validated in the executed build/replay packet. All 59,033 original
source records and 118,066 selected median-TPM cells replay exactly from source
text. 58,988 canonical gene rows retain 117,976 ordinary-gene cells; 45 qualified
Y records retain their 90 cells separately. No source record is lost, no zero is
turned missing, and no protein reference row is renamed as a gene. Reviewed
cross-references retain 22,320 candidate pairs: 39,630 genes unmapped, 19,208
with one reviewed protein, 150 with multiple proteins, and 1,713 source genes
with a protein shared with another reference gene. All reverse mappings remain
explicit; no protein measurements are projected. Source positions/versions and
mapping counts stay outside the biological-feature table.

Per-column context, species, units, source hashes, evidence grade, redistribution
gaps and typed measurement provenance are retained. 618 relevant checks pass;
two wrapper provenance-role failures were recorded and resolved before output
writes. All installed host/parasite input hashes remain unchanged. Build/opening,
graph-order, distributable-pack licensing and reference/biological admission
still require their gates; no additional complete action earns a tick.

## Frozen documentary partition: reference and public redistribution terms

2026-10-08, before audit jobs. Verify only the existing human foundation and its
pinned GTEx/UniProt source hashes; capture authoritative GTEx public-data terms
(actual rendered/frontend text, not a JavaScript shell), UniProt license JSON,
and GENCODE 39/47 release-to-assembly statements. Keep URL, retrieval date, bytes,
SHA-256 and bounded-request outcomes in a new executed
`results/human_reference_terms_2026_10_08/` packet. Requests are limited to 8 MiB
per response and the worker to 400 MiB. Unreadable pages remain explicit gaps
while independent pages continue. This documentary review neither changes the
foundation, admits a benchmark, promotes v11 nor declares a distributable pack
or genome-wide completeness; downstream licensing/reference/opening gates remain
separate controlled work.

Documentary partition accepted in canonical
`results/human_reference_terms_2026_10_08_v5/`: 25 foundation/source receipts
reverified twice, five authoritative captures and five GTEx V10→GENCODE39 checks.
Official frontend license component/route text confirms conditional public-data
redistribution: source/date acknowledgment, no implied endorsement, current data
or a visible stale-version notice. Future pinned v10 packs must implement those
conditions. UniProt copyrightable database parts use CC BY 4.0 with attribution
and a separate other-rights disclaimer. GTEx's own release table/query setting
links V10 to GENCODE39/GRCh38; the official GENCODE39 page identifies GRCh38.p13.
GENCODE47 is GRCh38.p14; the captured frontend does not establish V11's reference.
Sources: [GTEx terms](https://gtexportal.org/home/license),
[UniProt license JSON](https://rest.uniprot.org/help/license),
[GENCODE39](https://www.gencodegenes.org/human/release_39.html).

19 final output receipts verify. Executed final review took 0.44 s, 111.2 MiB
peak under 400 MiB, with zero new requests (verified captures reused). Earlier
wrapper-import failure, request-ceiling and text-review gaps remain preserved;
V2 inherited-launcher high-water measurement is explicitly superseded by actual
executed-process VmHWM. No original foundation or admission fields were mutated.
Pack notices/attribution, reference-gene coverage/mapping, graph/order/opening,
release-specific publication lineage and biological benchmarks remain open.
64.21 stays 30%; this documentary gate does not admit a distributable host pack.

## Second partition: HPA annotations and shipped-protein links (2026-11-07)

`scripts/build_human_gene_space_notebook.py` → `notebooks/human_gene_space_2026_11_07.ipynb`,
`results/human_gene_space_2026_11_07/`. Built on the foundation without changing it (checked: same genes, same
order, same values).

- Human Protein Atlas annotations (RNA tissue specificity/distribution, single-cell specificity, subcellular
  location and reliability, secretome location, protein class) attached by Ensembl id: 19,706 of 58,988 genes.
  1,035 HPA rows have no Ensembl id and 456 HPA genes are outside the GTEx v10 universe; both are listed, none
  added as genes. HPA licence CC BY-SA 4.0 recorded per column (attribution and share-alike on a pack).
- The shipped protein table (`hs_host_proteins.parquet`, 20,989 UniProt rows) linked through reviewed
  cross-references only: 17,714 proteins to one gene, 1,572 to several genes, 1,703 to none
  (`proteins_without_gene.tsv`); 145 genes have more than one shipped protein. The links are their own table
  (`hs_protein_gene_links.parquet`); no protein measurement is copied onto a gene.
- `column_metadata.json` describes every column (species, unit, evidence grade, source hash, licence,
  redistribution, missingness, allowed use). Checks pass: unique ids, every foundation gene kept in order,
  foundation values unchanged, no protein id used as a gene, every shipped protein in the links, every column
  described.

Remaining: the distributable pack (pack format, GTEx/UniProt/HPA attribution notices, opening independently in
ParaScape) and its benchmark gate.

## Last partition: the distributable pack (2026-11-07)

`scripts/build_human_pack_notebook.py` → `notebooks/human_pack_2026_11_07.ipynb`, `results/human_pack_2026_11_07/`
(`parascape_Hs_2026.11.07.zip`, 2.5 MB, SHA256 in `SHA256`). Members: `hs_gene_nodes.parquet` (58,988 genes),
`hs_gene_graph.npz` (3-D display layout: UMAP of log GTEx expression + HPA category indicators, cosine, seed 0;
display only, never an input), `hs_protein_gene_links.parquet`, `column_metadata.json`, `NOTICE.txt`.

- Licence per column: GTEx and UniProt-derived columns O, HPA columns SA. `NOTICE.txt` carries what the terms
  require: GTEx source/date acknowledgment, no endorsement, a stale-version pointer; UniProt CC BY 4.0; HPA
  CC BY-SA 4.0 with the share-alike condition.
- `spaces.build_pack` gained `extra=` so notices, column descriptions and link tables travel hash-checked in the
  pack and cannot replace the node table or graph (`tests/test_packs.py`, 16 pass).
- Opened on its own: installed into an empty cache with its checksum, read back through `organisms.nodes_path` /
  `graph_path`: identical table, valid graph order, notice present, 0 protein rows used as genes.

Gaps, stated: Hs is not yet registered in `organisms.SPACES` (the app does not list it until 64.24 wires host
spaces in), the pack is not hosted for download, and biological benchmarks on human labels belong to 64.28 (the
human inference pack). No PyPI or pack release before the licence review in item 90.
