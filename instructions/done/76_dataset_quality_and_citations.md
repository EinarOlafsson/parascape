# 76 · Dataset quality, including citations normalised by age
**Status: done 2026-11-06.** Built and measured: sample quality predicts replicate agreement (0.47 over 41 studies); citations do not (−0.19), so they describe a dataset and do not weight it.
(`results/rna_sample_quality_2026_11_04/`) upward to the dataset, and adds a bibliometric component.

## The request (Einar)
A quality score for each dataset, and "please try to keep track of how many citations the study that the data
came from has, this should be normalized also by time the study has been out, this can be used as an additional
quality metric."

## The plan
1. **Citations per year since publication.** Fetch the citation count for each dataset's publication and divide
   by years since it appeared. Source: Europe PMC's `citedByCount` keyed by PMID (it covers PubMed and gives a
   count without scraping), with NCBI E-utilities as the cross-check on the publication date. Store the raw
   count, the retrieval date and the derived rate -- never the rate alone, because it changes under you.
2. **Treat it as weak evidence and say so.** Citations measure attention, which this project has already caught
   predicting things it should not: an attention-only model reached AUC 0.902 on host-facing. So the citation
   rate is **a dataset-level weight, never a per-gene feature**, and it is in the `attention` information
   category, excluded by default from any inference about a gene.
3. **Combine with the measurable components.** Per dataset: the aggregate of its samples' technical quality,
   replicate agreement, the number of independent studies that reproduce its direction where testable, how
   completely it was processed by our own uniform pipeline, and the citation rate. Report the components
   separately as well as combined, because a well-cited dataset with poor replicate agreement should look
   different from a quiet dataset that reproduces cleanly.
4. **Where it is used.** Weighting samples when a column is built (already done per sample), ranking which
   datasets to acquire next, and as a stated caveat on a claim that rests on one dataset.

## Gate
The score is published only with its components visible, and with a measured statement of whether it predicts
replicate agreement on held-out datasets -- the same check the per-sample score had to pass (it reached 0.63,
falling to 0.40 within study).

## Done 2026-11-06

`scripts/build_dataset_quality_notebook.py` -> `notebooks/dataset_quality_2026_11_06.ipynb`,
`results/dataset_quality_2026_11_06/` (`dataset_quality.tsv`, `citations.tsv`, `gate.tsv`, README).

* **Citations**: Europe PMC `citedByCount` for 279 of the 281 catalogue datasets with a PubMed id (195
  papers), stored with the raw count, the retrieval date and the rate. Publication dates are checked
  against NCBI E-utilities; where the two years differ the earlier date is used and the row is flagged.
  Median 4.4 citations a year.
* **Components, shown separately**: sample quality (mean per-sample score without its replicate term, from
  the rebuilt `rna_sample_quality_2026_11_06` over 717 runs), split-half replicate agreement, per-run
  replicate correlation, reproduced elsewhere (share of a study contrasts that another study reproduces at
  Spearman at least 0.3), processing completeness, citations per year. Combined score = mean of ranks.
* **Gate, measured.** The score has no fitted weights, so every study is held out by construction.

| predictor | against per-run replicate correlation (41 studies) | against split-half agreement (7) |
|---|---|---|
| sample quality | **0.47** (p 0.002) | 0.86 (p 0.014) |
| combined score without replicate terms | 0.32 (p 0.045) | 0.75 |
| reproduced elsewhere | 0.24 (25 studies) | 0.47 |
| processing completeness | −0.11 | constant |
| citations per year | −0.19 (16 studies) | 2 studies only |

**Decision**: where a dataset weight is needed, use sample quality. The combined score is published as a
description, with its components; the citation rate stays in the `attention` category and weights nothing.

Found on the way: the Europe PMC query needs brackets (`(a OR b) AND SRC:MED`), or the source filter binds to
the last id only and 6 of 195 papers come back; the per-sample quality builder now takes its date as one
constant and was rebuilt on the 735-run snapshot (`results/rna_sample_quality_2026_11_06/`).
