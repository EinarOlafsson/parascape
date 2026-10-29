# RNA expression and GitHub data release audit

Prepared ParaScape0.0.2.0 at the user’s request. All296 completed study/run entries
across23 studies are retained, with8,164 genes in the common ME49_host reference.
Both quantifiers and all already-verified reference variants are included.
Sample design covers all1,011 planned entries, source annotations for256 BioSamples
and explicit pending coverage; Trx1 has four author-declared biological replicates
per group and was prioritized in the live worker after its current completed run.

Scientific application data is outside the PyPI distribution: a versioned GitHub
asset is pinned by parascape/data_manifest.json with archive and per-file SHA256.
First-run setup is consent-based, cancellable, resumable and staged; damaged or
incomplete archives never replace working data. Icons remain application assets.

Evidence: notebooks/rna_snapshot_release_2026_10_29.ipynb; source/output hashes
in the snapshot; scripts/build_rna_snapshot.py, extend_rna_references.py,
enrich_rna_samples.py and build_data_release.py. No feature/model/calibration
values were changed. Raw reads and quantifier index binaries are not application
browsing data and stay in the external processing archive.

Validation before publication:655 focused release gates; clean wheel imports
without datasets; missing-data popup and cancellation; verified local asset
installation; offline app load and one gene across all296 runs. Full regression
and publication verification results will be appended after completion.
