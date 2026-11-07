# Dataset quality, with citations normalised by age (instruction 76)

Notebook: `notebooks/dataset_quality_2026_11_06.ipynb` (executed). Builder: `scripts/build_dataset_quality_notebook.py`.

692 datasets; 279 with a citation count (Europe PMC, retrieved 2026-11-06); median 4.4 citations a year.

## Gate: predicting replicate agreement (Spearman across studies)

| target | predictor | studies | spearman | p |
|---|---|---|---|---|
| replicate_agreement | score_without_replicate | 7 | 0.75 | 0.0522 |
| replicate_agreement | sample_quality | 7 | 0.857 | 0.0137 |
| replicate_agreement | reproduced_elsewhere | 6 | 0.471 | 0.346 |
| replicate_agreement | processing_completeness | 7 | nan | nan |
| replicate_agreement | citations_per_year | 2 | nan | nan |
| replicate_rho | score_without_replicate | 41 | 0.315 | 0.0451 |
| replicate_rho | sample_quality | 41 | 0.472 | 0.0018 |
| replicate_rho | reproduced_elsewhere | 25 | 0.239 | 0.2501 |
| replicate_rho | processing_completeness | 41 | -0.107 | 0.5059 |
| replicate_rho | citations_per_year | 16 | -0.186 | 0.4904 |

**Reading.** Over the 41 studies with replicate groups, the component that predicts replicate agreement best is `sample_quality` (Spearman 0.47, p 0.002). The combined score without replicate terms reaches 0.32: averaging in the other components dilutes it. Citations per year reach -0.19 over 16 studies -- attention does not predict reproducibility here, so the citation rate is a description of a dataset, not a weight on it. Use `sample_quality` where a dataset weight is needed.

Files: `dataset_quality.tsv` (every component and the combined score), `citations.tsv` (raw counts, dates, rates), `gate.tsv`. The citation rate is a dataset weight in the `attention` category, never a gene feature.
