# 64.22 · Build the mouse gene information space

Status: DONE — 100%, 2026-11-07. Mouse gene table, protein links and a distributable pack that opens on its own; gaps stated below. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.02, 64.03, 64.04, 64.05, 64.08 |
| Estimated remaining engineering time | 8–16 h |
| Existing code to inspect | `parascape/organisms.py`, `parascape/host.py`, `parascape/spaces/`, `parascape/packs.py`, `parascape/datasets.py` |

## Deliverable

One reproducible Mm gene-space builder using already acquired, verified reference and selected measurement sources.

## Controlled scope

One host species and a named pilot source set; share builder contracts with item 21 while keeping species data separate.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Apply explicit gene/protein mappings and source-specific mouse tissue/condition semantics with measured mapping loss.
- [x] A notebook reproduces the build and source sanity checks; unique IDs, non-loss and graph-order checks pass.
- [x] A validated distributable mouse pack opens independently with its own available targets and benchmark gaps.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

2026-11-07. `scripts/build_mouse_gene_space_notebook.py` → `notebooks/mouse_gene_space_2026_11_07.ipynb`,
`results/mouse_gene_space_2026_11_07/` (`parascape_Mm_2026.11.07.zip`, 0.8 MB, SHA256 in `SHA256`); the space is
declared in `parascape/spaces/mouse.py` (`MOUSE_SPACE`, Ensembl mouse gene ids, pack distribution).

- 18,428 mouse genes (Ensembl) from the 22,132 source rows; 3,704 rows without an Ensembl id are listed
  (`rows_without_gene_id.tsv`), not guessed.
- Columns: Bgee best score, best tissue and tissues present (Bgee scores are per tissue, 0-100); human ortholog,
  its type and confidence (Ensembl Compara); and subcellular location TRANSFERRED from the human ortholog (HPA),
  named `location_from_human_ortholog*`, graded inferred and marked "never an observed label".
- Shipped mouse proteins (`mm_host_proteins.parquet`, 15,590): 14,439 link to one gene, 1,151 to none
  (`proteins_without_gene.tsv`); links are their own table, no protein row is used as a gene.
- Licence per column: Bgee CC0 and Ensembl open (O); transferred HPA columns CC BY-SA 4.0 (SA); `NOTICE.txt`
  in the pack says so.
- Checks pass: unique ids, every source gene kept, sorted order, graph order, no protein ids as genes, every
  shipped protein in the links, pack installed into an empty cache and read back identical with its notice.
  `tests/test_packs.py` (17 pass, including the mouse space declaration).

Gaps, stated: the source files (assembled 2026-07-17) do not record Bgee/Ensembl release versions and their
assembling script is not in the tree; gene->protein comes from one source column, not reviewed many-to-many
cross-references; there are no observed mouse labels yet, so benchmarks belong to 64.29; Mm is not registered
in `organisms.SPACES` until 64.24 wires host spaces into the app.
