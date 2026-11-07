# Perturbation atlas: recovery benchmark (instruction 71, product 2)

`parascape/perturbation_atlas.py` (store, query, recovery); notebook
`notebooks/perturbation_atlas_2026_11_07.ipynb`; builder `scripts/build_perturbation_atlas_notebook.py`.
Signatures: the within-study contrast columns of `expression_traits_2026_11_06` (shrunken log2 fold change against
each study's own control). Query: Spearman over shared genes; the query's own study is always held out.

| perturbation (curated, independent studies) | queries | top 1 | top 5 | median best partner rank | chance rank |
|---|---|---|---|---|---|
| BFD1 stabilised by Shield-1 (GSE139537, GSE223877) | 2 | 100 % | 100 % | 1 | 98.8 |
| BFD1 lost under stress (GSE132238, GSE223819, GSE218921) | 5 | 40 % | 60 % | 3 | 59.3 |
| Alkaline stress vs normal, wild type (GSE132238, GSE218921, GSE223869) | 5 | 0 % | 80 % | 3 | 59.3 |
| AP2XII-1 knockdown (GSE224934, GSE228580; one lab, two deposits) | 5 | 60 % | 100 % | 1 | 66.3 |

Every query finds another lab's experiment of the same perturbation within its top 10, against chance ranks of
37-99. Alkaline stress is never first; why has not been looked at. `example_query.tsv`: the top ten
matches for the GSE218921 BFD1 knockout under stress.
