# Co-expression places proteins: localization from RNA alone (instruction 71 product 1; 70.08)

Notebook `notebooks/coexpression_components_localization_2026_11_06.ipynb`; builder
`scripts/build_coexpression_components_localization_notebook.py`. Complements the first product 1 build
(`coexpression_localization_2026_11_04`, nearest-neighbour votes, coarse accuracy 0.466 against 0.132).

**Built:** within-study Spearman co-expression for 57 studies with at least 6 bulk runs (735-run snapshot),
ranked within each study and averaged; the consensus reduced to 20 leading components -- label-free, so no
frozen-test label can leak in. **Scored:** the 13 coarse hyperLOPIT niches (718 curation-grade markers), the
20 components as the only evidence, boosted trees, nested cross-validation, a frozen partition per niche
(seed 2).

| niche | co-expression only (cross-validated AUC) | all allowed evidence |
|---|---|---|
| ribosome | 0.999 | 0.998 |
| rhoptries | **0.975** | 0.974 |
| apical | **0.958** | 0.943 |
| proteasome | 0.958 | 0.978 |
| micronemes | **0.936** | 0.932 |
| IMC | **0.913** | 0.879 |
| nucleus | 0.897 | 0.940 |
| dense granules | 0.872 | 0.957 |
| mitochondrion | 0.813 | 0.844 |
| apicoplast | 0.803 | 0.927 |
| cytosol | 0.799 | 0.913 |
| ER | 0.780 | 0.831 |
| plasma membrane | 0.669 | 0.934 |

**Median 0.897 from co-expression alone against 0.934 from everything allowed.** For the organelles made in
waves during division -- rhoptries, micronemes, the apical complex, the IMC -- co-expression alone matches or
beats all other evidence together; it is weakest for the plasma membrane and the ER, which are not made on a
shared schedule. Nothing may fill (4-13 positives per frozen partition). The cross-validated column is the
steadier one; frozen-test numbers are in `scores.tsv`.
