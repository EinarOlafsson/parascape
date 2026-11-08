# 90 · Release 0.0.10.0 under the new licence

**Status: DONE 2026-11-08.** Einar approved the release under the declared licence; 0.0.10.0 is on PyPI and GitHub (instruction 64.40). Previously: waits on Einar. The spaCR port put noncommercial-licensed files into ParaScape (Einar chose to
ship them under spaCR's terms). `pyproject` declares `MIT AND LicenseRef-spaCR-BSD-3-Clause-Additional-Terms`
and LICENSE states that commercial use is no longer granted by MIT alone. The standing permission to publish
once CI is green was given when ParaScape was wholly MIT, so publishing waits for Einar's review.

Also outstanding before a release: `main` is far behind `nightly`, and the data manifest still points at the
v0.0.9.0 data release while the package is 0.0.10.0.
