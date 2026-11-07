# Host response atlas (instruction 71, product 3)

Notebook `notebooks/host_response_atlas_2026_11_07.ipynb`; builder `scripts/build_host_response_atlas_notebook.py`.
The contrast-ladder contrasts read on the **host** side of the competitive salmon quantification: host
transcripts summed to genes, CPM over the run's own host species, log2(CPM + 1), treated minus baseline means,
with a standard error where both arms are replicated. **106 host signatures across 18 studies** (human or mouse,
never pooled): `host_signatures.parquet`, `host_signature_se.parquet`, `registry.tsv`.

**The check (GSE125120, infected against bystander cells):** wild type changes **344** host genes at least
2-fold; the MYR1 knockout **173** -- about half, as expected if MYR1-exported effectors drive much of the
response -- and the two responses correlate at 0.536 (`myr1_check.tsv`). MYC moves in neither (log2 change 0.05
and 0.08); its MYR1-dependent induction is described for tachyzoites, and these samples are labelled as
bradyzoite infections (not checked further).

Not done: human-mouse orthologue mapping for cross-species comparison, a query function like the perturbation
atlas, and contrasts against mock-infected cells (the zero-parasite runs are outside the ladder).
