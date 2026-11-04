# Host-facing inference, stages 0 and 1 (2026-11-04)

Built by `scripts/build_host_facing_inference_notebook.py`; the executed notebook is
`notebooks/host_facing_inference_2026_11_04.ipynb`. Instruction 72 inference plan, stages 0
(baselines) and 1 (one gradient-boosted-tree model), scored by the instruction 70 rules.

## The question

Can the `host_facing` column of the information map be filled by prediction at the 90 % bar
Einar set -- that is, is there a calibrated probability threshold at which the **lower 90 %
confidence bound** of held-out precision reaches 0.90? Only at or above such a threshold may a
predicted cell be written as a fact.

## The answer

YES. The lowest calibrated threshold whose lower 90 % bound of precision reaches 0.90 is 0.90: 55 held-out calls, 53 correct, precision 0.964, lower bound 0.906. At that threshold 93 currently-empty genes would be filled.


Restricted to proteins with a predicted signal peptide -- the question that matters, because the overall negatives are mostly proteins that live inside the parasite and are separable by signal peptide alone -- **never reaches the bar; best lower bound 0.868**.


**And the overall answer is conditional on two feature families that sit close to the conclusion.** The most important features are a screen for effectors that alter host transcription (`hosttx_padj`, `hosttx_T2`) and a publication count (`n_papers_incidental`). Removing them, and the other proximity-style assays:


```
                                              feature_set  n_features  roc_auc  average_precision  best_lower90  at_threshold  calls  correct  reaches_bar
                                      all tier-2 features         409   0.9946             0.9625        0.9061          0.90     55       53         True
                                          minus attention         402   0.9945             0.9472        0.8771          0.85     63       59        False
                               minus host-effector screen         385   0.9950             0.9635        0.8872          0.97     33       32        False
                 minus host-effector screen and attention         378   0.9940             0.9360        0.9027          0.95     53       51         True
minus host-effector, attention and proximity-style assays         371   0.9938             0.9319        0.9027          0.95     53       51         True
```


So: the bar is still reached with the host-effector screen, the attention columns and the proximity-style assays all removed.


**But it is met by a hair, and not robustly.** Across those five feature sets the best lower bound ranges 0.877 to 0.906, and 2 of 5 fall below 0.90 -- including the set that removes only the host-effector screen. The bar being cleared at 0.903-0.906 while a neighbouring feature set gives 0.887 means the result sits inside the noise of a 73-positive label set. The right reading is "this is at the edge of what the existing labels can demonstrate", not "the fill bar is cleared".

## Method

* Labels: `host_facing` from `results/information_map_v0_2026_11_04/values.parquet`, observed
  facts only: 73 True and 486 False after excluding the
  8 genes in `conflicts.tsv` (TGME49_227620, TGME49_251540, TGME49_263090, TGME49_269690, TGME49_277270, TGME49_289380, TGME49_312420, TGME49_316250) from both training and
  scoring. 7573 genes are untested and stay empty.
* Features: 409 numeric/boolean columns of `nodes.parquet` plus nine sequence and
  topology traits. Excluded as leakage: 8 columns matching
  lopit|compartment|ortholopit|tagm, plus `pvm_proximity_positive` and `n_host_targets`, plus
  all identity/free-text columns. Full lists are printed in the notebook and in
  `importance.tsv` (family per feature).
* Folds: 5, grouped by OrthoMCL orthogroup so paralogs never split across folds, seed 0.
* Model: `HistGradientBoostingClassifier`, missing values passed through with no imputation,
  `nterm_side`/`cterm_side` as declared categoricals, 300 rounds at learning rate 0.05, the same
  for every fold, early stopping OFF so the fit is deterministic and no training gene is held
  back. No hyperparameter search; nothing was tuned on any held-out fold.
* Calibration: isotonic (`parascape.claims._isotonic`), cross-fitted across folds; expected
  calibration error 0.0069 over 10 quantile bins (uncalibrated
  0.0133).
* Bound: exact one-sided Clopper-Pearson, `scipy.stats.beta.ppf(0.10, k, n-k+1)`.

## Headline numbers

```
metric                                       accuracy  precision  recall  roc_auc  average_precision  calls  coverage
strategy                                                                                                             
GBT: all tier-2 features                        0.977      0.926   0.899    0.995              0.956   14.0       1.0
GBT: attention features only                    0.916      0.699   0.628    0.902              0.700   12.8       1.0
GBT: no attention features                      0.977      0.953   0.886    0.996              0.966   13.6       1.0
GBT: topology features only                     0.904      0.643   0.549    0.903              0.674   13.0       1.0
baseline: ortholog copy (grouped, abstains)     0.870        NaN   0.000      NaN                NaN    0.0       0.0
baseline: prevalence-weighted random            0.782      0.101   0.079    0.500              0.130   12.2       1.0
baseline: signal peptide alone                  0.934      0.753   0.751    0.856              0.601   14.8       1.0
baseline: training majority                     0.870        NaN   0.000      NaN                NaN    0.0       1.0
```

Secretory-pathway proteins only (74 genes, 55 positives):

```
metric                                            accuracy  average_precision  calls  precision  recall  roc_auc
strategy                                                                                                        
secretory · GBT                                      0.887              0.970   11.2      0.912   0.960    0.938
secretory · prevalence-weighted random               0.673              0.753   11.0      0.773   0.786    0.500
secretory · signal peptide alone (constant True)     0.753              0.753   14.8      0.753   1.000    0.500
secretory · training majority                        0.753                NaN   14.8      0.753   1.000      NaN
```

Publication-bias ablation:

```
                   roc_auc  average_precision  best_lower90
all features        0.9946             0.9625        0.9061
without attention   0.9945             0.9472        0.8771
drop                0.0000             0.0153           NaN
```

Top 10 features by permutation importance (held-out folds, average precision):

```
                                                 family  importance_mean
feature                                                                 
hosttx_padj                            hosttx_effectors           0.0232
n_papers_incidental                              pubmed           0.0201
ees_vs_tachyzoite_log2          toxodb_enteroepithelial           0.0139
has_signal_peptide                         unregistered           0.0029
hosttx_T2                              hosttx_effectors           0.0022
fit_no_glucose                crispr_glucose_limitation           0.0019
cellcycle19092_9h_r2                           gse19092           0.0008
surface_mitochondrion_log2fc  organelle_surface_turboid           0.0006
rna206344_Unsporulated_R1                     gse206344           0.0006
cellcycle19092_5h_r1                           gse19092           0.0004
```

## Files

| File | Contents |
|---|---|
| `scoreboard.tsv` | strategy x metric, mean plus every fold |
| `scoreboard_secretory.tsv` | the same restricted to predicted signal-peptide proteins |
| `thresholds.tsv` | threshold, calls, precision, Clopper-Pearson lower 90 % bound, genes fillable |
| `thresholds_secretory.tsv` | the same, secretory proteins only |
| `calibration.tsv` | reliability table: predicted range, n, observed positive rate |
| `by_label_source.tsv` | performance per positive class x negative class |
| `ablations.tsv` | the lower 90 % bound with the host-effector screen, the attention columns and the proximity-style assays removed |
| `importance.tsv` | permutation importance per feature with its provenance family |
| `predictions.tsv` | all 8,140 genes: calibrated probability, out-of-fold probability, fold, label, per-gene occlusion attributions |


## Limitations (none of these are softened)

1. **73 positives.** ~14 per fold. The Clopper-Pearson bound cannot
   reach 0.90 on fewer than ~22 calls even with perfect precision, so over most of the threshold
   range the bar is unreachable as arithmetic, before any question of model quality.
2. **The negatives are easy.** 486 negatives, almost all hyperLOPIT markers of
   compartments inside the parasite. "Has a signal peptide" nearly separates the classes on its
   own, so the overall numbers overstate usefulness; the secretory-only numbers are the honest
   ones.
3. **Most positives are one assay.** 52 of
   73 positives are PVM proximity-labelling hits and nothing else (tier probability
   0.70). Proximity to the PVM is not proof of export, and the 8 excluded conflicts are exactly
   where curation disagrees with that screen.
4. **Publication bias is only partly removable.** The ablation drops the counting columns, but a
   well-studied gene is also the gene with a melting temperature, a palmitome entry and a
   phosphoproteome row. Missingness itself carries fame.
5. **DeepTMHMM misses signal peptides** (GRA16 / TGME49_208830 is called GLOBULAR). The
   secretory-restricted analysis is therefore run on a set missing some true secreted proteins,
   and `mature_len` is wrong wherever the prediction misses. SignalP would fix it; open gap.
6. **The negatives are an inference, not an export test** (residence in an internal compartment,
   probability 0.80). Precision against them is precision against a curation decision.
7. **The ortholog-copy baseline has coverage 0** under orthogroup-grouped folds, by construction.
8. **The all-gene probability column is not independently validated**: model fitted on all labels,
   isotonic curve fitted in-sample on pooled out-of-fold scores. Only the cross-fitted
   out-of-fold numbers are honest.
9. **One model, fixed settings.** Stages 2-5 of the instruction 72 plan are not in this notebook.
10. **The signal is concentrated in two families that sit close to the conclusion.** The top
    permutation-importance features are `hosttx_padj` / `hosttx_T2`, from a screen for effectors
    that target host-cell transcription, and `n_papers_incidental`, a publication count. Neither
    defines the label, so neither is leakage under instruction 72, but a model built on them is
    partly re-reading one host-effector screen and partly reading how famous a gene is.
    `ablations.tsv` is the honest version of the headline result.
11. **Permutation importance is noisy here.** With ~112 genes and ~15 positives per held-out
    fold the individual importances are small and their signs flip between folds
    (`hosttx_padj` is -0.08 in one fold and +0.05 in another). Read the top of that table as a
    direction, not as a per-feature ranking to be trusted.
