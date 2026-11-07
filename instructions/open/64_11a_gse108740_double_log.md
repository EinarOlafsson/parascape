# 64.11a · Take the GSE108740 expression log only once

Status: OPEN — 25%, 2026-11-07. Scale confirmed from GEO. Registered from 64.11 (units curation).

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

## Finding

The GSE108740 columns in `nodes.parquet` (`rna108740_*_FPKM`) are labelled FPKM but range 0 to 8.9 (median 2.7):
they are already on a log scale. `build_graph.py` then computes `expr_tachy`, `expr_cyst` and `expr_max` as
log2(mean + 1), so these three columns are doubly logged. Rankings and correlations are unaffected (monotone);
values, MAE/RMSE, value-strategy calibration and any comparison with other expression columns are compressed.

## Done when

- [x] The scale of the GSE108740 series is confirmed from GEO / the paper (log2 FPKM, or log2(FPKM + 1)).
- [ ] `build_graph.py` takes the log once (or not at all) for these summaries; `units.CURATED` updated to match.
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
