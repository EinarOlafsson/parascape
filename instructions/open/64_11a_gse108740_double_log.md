# 64.11a · Take the GSE108740 expression log only once

Status: OPEN — 60%, 2026-11-07. Scale confirmed from GEO; Einar chose fix and rebuild; table and map fixed, track record/claims/calibration rebuild next. Registered from 64.11 (units curation).

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

## Finding

The GSE108740 columns in `nodes.parquet` (`rna108740_*_FPKM`) are labelled FPKM but range 0 to 8.9 (median 2.7):
they are already on a log scale. `build_graph.py` then computes `expr_tachy`, `expr_cyst` and `expr_max` as
log2(mean + 1), so these three columns are doubly logged. Rankings and correlations are unaffected (monotone);
values, MAE/RMSE, value-strategy calibration and any comparison with other expression columns are compressed.

## Done when

- [x] The scale of the GSE108740 series is confirmed from GEO / the paper (log2 FPKM, or log2(FPKM + 1)).
- [x] `build_graph.py` takes the log once (or not at all) for these summaries; `units.CURATED` updated to match.
- [ ] The shipped table, value-strategy calibration and the track record entries that use these columns are
      rebuilt (`precompute_track_record` reuses everything not depending on the node table only when the table is
      unchanged, so this is a full Tg rebuild), and the change is recorded in the CHANGELOG.
- [ ] The result is committed and pushed to nightly.

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
