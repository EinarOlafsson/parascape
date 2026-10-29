# RNA expression and GitHub data release audit

Prepared ParaScape 0.0.2.0 at the user’s request. All 296 completed study/run entries
across 23 studies are retained, with 8,164 genes in the common ME49_host reference.
Both quantifiers and all already-verified reference variants are included.
Sample design covers all 1,011 planned entries, source annotations for 256 BioSamples
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

Validation before publication: 655 focused release gates; clean wheel imports
without datasets; missing-data popup and cancellation; verified local asset
installation; offline app load and one gene across all 296 runs. Full regression
and publication verification results will be appended after completion.


Final local full regression run:5,901 passed,41 skipped,7 deselected. Two documentation
checks failed: short module explanations (including pre-existing modules) and a stale
cache-size statement. The explanations and cache/download documentation were fixed;
495 focused documentation/number/docstring checks then passed. The cloud full-suite
Tests and Documentation checks remain mandatory publication gates.

GitHub data draft release ID428109947 contains the221,482,867-byte asset. GitHub
reported SHA256891b93e69bff81e9f2d0985bb35792cec4e7a63c3d653e5d7fb22be58ae24f74,
matching the package descriptor. Workflow push is blocked by the CLI OAuth login
lacking workflow scope. Native service parascape-publish-prepared.service waits
for `gh auth refresh -h github.com -s workflow`, then pushes the frozen commit,
requires green CI, publishes the verified asset, fast-forwards main without force
and verifies the resulting PyPI wheel. It makes no model calls.
