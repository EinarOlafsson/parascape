# 64.11a · Take the GSE108740 expression log only once

Status: DONE — 100%, 2026-11-08. Table, map, track record, claims, Tg frozen functional recoveries and Tg calibration rebuilt on the corrected expression columns. Built; Einar's review pending -- reopen on his feedback. Registered from 64.11 (units curation).

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

## Finding

The GSE108740 columns in `nodes.parquet` (`rna108740_*_FPKM`) are labelled FPKM but range 0 to 8.9 (median 2.7):
they are already on a log scale. `build_graph.py` then computes `expr_tachy`, `expr_cyst` and `expr_max` as
log2(mean + 1), so these three columns are doubly logged. Rankings and correlations are unaffected (monotone);
values, MAE/RMSE, value-strategy calibration and any comparison with other expression columns are compressed.

## Done when

- [x] The scale of the GSE108740 series is confirmed from GEO / the paper (log2 FPKM, or log2(FPKM + 1)).
- [x] `build_graph.py` takes the log once (or not at all) for these summaries; `units.CURATED` updated to match.
- [x] The shipped table, value-strategy calibration and the track record entries that use these columns are
      rebuilt (`precompute_track_record` reuses everything not depending on the node table only when the table is
      unchanged, so this is a full Tg rebuild), and the change is recorded in the CHANGELOG.
- [x] The result is committed and pushed to nightly.

## Evidence (2026-11-07)

GEO series GSE108740 ("Transcriptomics of Toxoplasma gondii enteroepithelial stages"; SOFT file via the GEO FTP
mirror): every sample's data processing says "Supplementary_files_format_and_content: Excel table with FPKM
values", FPKM from featureCounts with ToxoDB 24 transcript lengths; the supplementary file is
`GSE108740_FPKM.xlsx` (8,322 genes, 12 samples, plain FPKM, e.g. 6,955 for TGME49_270250). The shipped
`rna108740_*_FPKM` columns equal ln(FPKM + 1) of that table exactly (largest difference 2e-15 over 7,739 shared
genes; log2(x+1) and log10(x+1) do not fit). So the node table carries natural-log values under an FPKM name, and
`expr_tachy`/`expr_cyst`/`expr_max` = log2(mean(ln(FPKM + 1)) + 1). `units.CURATED` now states this exactly.
The ln is applied upstream of `build_graph.py` (where the columns are read in); the fix belongs there plus a
rename (or a unit) for the `rna108740_*` columns, and the rebuilds listed above.

## Fixed (2026-11-07, Einar: fix and rebuild now)

`build_graph.load_nodes` now undoes the stored natural log (expm1) of the GSE108740 columns and takes log2 once:
`expr_tachy`/`expr_cyst` = log2(mean FPKM + 1), and `expr_max` is the highest GSE108740 sample only (it used to
take the maximum over every `rna*` family, most of which are centred log values with median 0). The raw tree
`build_graph` reads is not on this machine, so `scripts/patch_expression_summaries.py` applies the identical
formula to the shipped `nodes.parquet` and recomputes the packaged layout with `build_graph.embed`, replacing only
`xyz` and `layout_metadata` (with the new node-file hash) in `graph.npz`; edge layers are untouched. Record in
`results/expression_fix_2026_11_07/record.json`. New ranges: expr_tachy 0-12.8 (median 4.0), expr_cyst 0-13.7,
expr_max 0-14.0; gene order nearly unchanged (Spearman with the old values 0.99999 / 0.997 / 0.985).
`units.CURATED` states the corrected units, and that the `rna108740_*` columns themselves are ln(FPKM + 1).
`expr_sporulated` (GSE206344, 0.23-2.01) is left as is: its series scale is not yet checked.
529 data-dependent tests pass, including the packaged-layout hash replay.

Still to do: rebuild the Toxoplasma track record and claims (their features include these columns), refresh the
value-strategy calibration for Toxoplasma, CHANGELOG entry.

## Rebuilt (2026-11-07/08)

- **Track record (Tg):** rebuilt with `build_track_record.py --organism Tg --workers 4` and merged with the
  unchanged Pf rows (1,038,372 rows). Accuracy per label moves by at most 0.002 -- expected, since the fix keeps
  each gene's rank.
- **Claims:** rebuilt from the new record: 19,798 (was 19,641).
- **Frozen functional recoveries (Tg):** bound to the old node bytes, so refused; re-run on fresh splits (Einar:
  "Re-run on a fresh split") -- EC FN-EC-01 seed 24 (184 test genes, native replay exact) and Pfam seed 20261107
  (686 test genes, 317 with training-unsupported profiles); Pf entries unchanged
  (`scripts/replace_tg_functional_results.py`).
- **Calibration:** a fresh Tg sweep (4,275 runs, `results/calibration_2026-11-07/`) published together with the
  unchanged Pf runs of the shipped 2026-09-26b sweep (Pf calibration identical to before). Tg changes: conjunctions
  default skill 0.49 -> 0.41 (still reliable); condition_shift and map_neighbours graded "works when tuned" (were
  "weak"); everything else within 0.03. The sweep cannot separate the expression fix from strategy code changes
  since 2026-09-26.
- **Slot table:** expr_max now covers the 7,739 GSE108740 genes (it used to reach 8,131 by mixing other series).
- The published data asset (0.0.9.0) still carries the old table; the next release publishes the corrected one
  (64.40).
