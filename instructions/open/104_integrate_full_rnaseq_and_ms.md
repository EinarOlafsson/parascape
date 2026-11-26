# 104 · Integrate all RNA-seq and the uniform MS re-search into the app

**Einar, 2026-11-26** (asked "Are all the rnaseq and ms datasets downloaded/ redownloaded and integrated?", then chose
"Both, RNA first").

## State found (2026-11-26 ~08:30)

- RNA-seq: 2,535 runs planned (priority plan), 2,498 quantified, 2,129 verified by the queue. The app's
  snapshot (`app_snapshots/current`) was 20261121T120206Z-2002: 2,002 runs, 135 studies. The snapshot service
  was killed by the 2026-11-21 09:49 reboot mid-build and its timer never fired again.
- MS: 207 PRIDE projects selected; 156 searched (Sage, uniform), 13 unsupported formats, 38 left (3 held:
  PXD000297, PXD078724, PXD052656). Nothing from the uniform search reaches the gene table.
- PRIDE deposited tables: 3,474 of 3,519 files acquired; measurements not admitted.

## Done

- [x] Snapshot service restarted (timer + one build of 2,129 verified runs); later batches follow automatically.
- [x] `parascape/ms_genes.py` + `scripts/ms/build_ms_gene_tables.py` + `tests/test_ms_genes.py`: per project a
      Toxoplasma gene x raw-file table (q <= 0.01, peptides unique to one parasite gene, Top3, log2,
      median-centred) and `prepared_ms_genes/summary.csv`.

## To do

- [ ] Verify the new snapshot validated and replaced `current`; release it to the app (release build).
- [ ] MS admission: per-sample design (PRIDE designs pipeline) for the raw files, an admission rule (e.g. >= 100
      genes quantified, parasite share of peptides), then a node/expression-style pack for the app.
- [ ] Retry PXD000297 at 48 GB; re-download the corrupt PXD052656 file.
