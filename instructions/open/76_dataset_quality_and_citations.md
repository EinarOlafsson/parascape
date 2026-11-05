# 76 · Dataset quality, including citations normalised by age
**Status: open (requested 2026-11-05).** Extends the per-sample RNA quality score
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
