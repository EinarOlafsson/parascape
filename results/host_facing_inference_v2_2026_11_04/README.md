# Host-facing inference v2 on Einar's reviewed labels (2026-11-04)

Built by `scripts/build_host_facing_inference_v2_notebook.py`; the executed notebook is
`notebooks/host_facing_inference_v2_2026_11_04.ipynb`. Rerun of
`results/host_facing_inference_2026_11_04/` after Einar reviewed all 140 published claims.
Run 1's outputs are untouched.

## Why this rerun

Run 1 reported ROC-AUC 0.995 and a lower bound of 0.906, and it was not trustworthy: its
negatives were 477 markers of compartments inside the parasite against 9 curated vacuole
residents, so it separated "inside the parasite" from "secreted" rather than "stays in the
vacuole" from "exported into the host"; against the 9 hard negatives its AUC fell to 0.855 with
a 33 % false-positive rate; and its top two features were a host-transcription effector screen
and a publication count. Einar then rejected 22 dense granule proteins the PVM proximity screen
had called positive, so the map now holds 43 vacuole-resident or curated negatives.

## The three questions and the answer to each

```
          reaches_bar       model  threshold  calls  correct  precision  lower90  best_lower90  fillable_all  fillable_secretory
question                                                                                                                        
hard            False  all-labels       0.65     45       37     0.8222   0.7273        0.7273            66                  46
reviewed        False      within       0.50     98       66     0.6735   0.6054        0.6054          7555                 875
easy            False      within       0.60     49       38     0.7755   0.6809        0.6809            66                  46
```


**1 · HARD (the headline).** NO for the HARD question. No calibrated threshold from 0.50 to 0.99, under either training variant, reaches a lower 90 % confidence bound of 0.90. The best bound anywhere on the curve is 0.7273 (all-labels model, threshold 0.65: 45 calls, 37 correct, point-estimate precision 0.822). Among secretory-pathway proteins we still cannot say at 90 % which are exported, so nothing may be written into the information map as a fact for this question.


**2 · HARDEST (Einar-reviewed genes only).** NO for the REVIEWED question. Best lower bound anywhere on the curve 0.6054 (within, threshold 0.50: 98 calls, 66 correct, precision 0.673).


**3 · easy (comparison only, the number run 1 reported).** NO even for the EASY question. Best lower bound 0.6809 at threshold 0.60 (49 calls, 38 correct).


Robustness of the hard question across feature sets: hard question: 0 of 14 feature sets x training variants reach the bar; best lower bound ranges 0.625 to 0.727.

With the two contested genes flipped (GRA9, NTPase II): with GRA9 and NTPase II as positives the hard question still does not reach the bar (best lower bound 0.7540).


**More labels cannot fix this one.** The best precision the model reaches anywhere on the hard
question — at any calibrated threshold, under either training variant — is **0.822** (45 calls, 37
correct), which is below the 0.90 bar itself. The Clopper-Pearson lower bound converges upward to
the true precision, so at 0.822 no number of additional labels can ever bring the bound to 0.90;
more labels would only measure 0.822 more precisely. The binding constraint for the hard question
is therefore **the evidence, not the label count**: against the 43 hard negatives the model has
ROC-AUC 0.667, against the hyperLOPIT internal markers 0.984, and that gap is the whole result.
Only once a model clears 0.90 precision does the label count become binding — at a true precision
of 0.95 the bound needs about 65 held-out calls, i.e. roughly 150 labelled secretory genes at the
current call rate against the 105 that exist today.

## Definitions

* **Secretory-pathway set (the hard question):** a predicted DeepTMHMM signal peptide OR a
  hyperLOPIT marker of dense granules / rhoptries 1 / rhoptries 2 / micronemes. 980 genes
  genome-wide, 105 of them labelled (58 host-facing / 47 not), 875 untested. 7 of the 65
  host-facing positives are outside it (no predicted signal peptide, no secretory marker).

* **Reviewed set (the hardest question):** the 98 genes Einar ruled on, with his decisions as
  labels (66 yes / 32 no).

* **Easy set:** all 585 labelled genes (65 / 520), hyperLOPIT internal markers included.

* **Label provenance** travels with every gene: hyperlopit_internal (477), einar_direct (48), einar_curated (48), gold_negative (11), host_interaction (1).

* **Conflicts:** all 33 genes in `conflicts.tsv` have Einar's review as a source, so all 33
  are resolved and kept. None dropped (run 1 dropped 8).

* **Two genes where the map and his review disagree:** GRA9 (TGME49_251540) and NTPase II
  (TGME49_277270). The map says not host-facing (curated gold negative, 0.90); his review
  accepted the PVM proximity claim (0.90). Questions 1 and 3 use the map column as the
  instruction names it; question 2 uses his decisions; the flip test is in the notebook.

## Method

* Features: **398** — numeric/boolean `nodes.parquet` columns, the DeepTMHMM sequence and
  topology traits, the 3Di quality columns (`tdi_coverage`, `tdi_align_identity`,
  `tdi_align_coverage`, `tdi_n_fragments`, `tdi_source`) and a new `foldseek_3di` family of
  21 3Di composition features. Full list with families is printed in the notebook and in
  `importance.tsv`.
* Held out as declared dependencies of `host_facing` (**47 columns**): every column of
  `hosttx_effectors` (24), `organelle_surface_turboid` (6), `pubmed` (4) and
  `pvm_proximity` (1); `n_host_targets`; the attention columns (`n_bioid_studies`, `n_fulltext`, `n_ipms_studies`, `n_papers_focal`, `n_papers_incidental`, `n_papers_substantive`, `n_publications`); and everything matching
  `lopit|compartment|ortholopit|tagm` (8).
* Folds: 5, grouped by OrthoMCL orthogroup so paralogs never split across folds, seed 0,
  assigned once over all labelled genes and inherited by every question.
* Model: `HistGradientBoostingClassifier`, missing values passed through with no imputation,
  300 rounds at learning rate 0.05, 15 leaves, early stopping OFF, identical for every fold
  and question. No hyperparameter search; nothing tuned on a held-out fold.
* Two training variants per question: `within` (fitted on that question's genes only) and
  `all-labels` (fitted on all 585 map labels, scored only on that question's genes). Both are
  out-of-fold; both are reported.
* Calibration: isotonic (`parascape.claims._isotonic`), cross-fitted across folds.
* Bound: exact one-sided Clopper-Pearson, `scipy.stats.beta.ppf(0.10, k, n-k+1)`.

## Scoreboards (mean over 5 folds)

### HARD (headline): among secretory-pathway proteins, which are exported?

```
metric                                       accuracy  precision  recall  roc_auc  average_precision  calls  coverage
strategy                                                                                                             
GBT (trained on all labels, scored here)        0.694      0.791   0.644    0.804              0.835    9.8       1.0
GBT (trained within the question set)           0.751      0.780   0.768    0.774              0.802   11.4       1.0
baseline: ortholog copy (grouped, abstains)     0.459        NaN   0.000      NaN                NaN    0.0       0.0
baseline: prevalence-weighted random            0.510      0.549   0.485    0.500              0.541   10.0       1.0
baseline: signal peptide alone                  0.502      0.516   0.835    0.459              0.524   18.6       1.0
baseline: training majority                     0.541      0.541   1.000      NaN                NaN   21.0       1.0
```

### HARDEST: the 98 genes Einar reviewed, his rulings as labels

```
metric                                       accuracy  precision  recall  roc_auc  average_precision  calls  coverage
strategy                                                                                                             
GBT (trained on all labels, scored here)        0.561      0.772   0.560    0.691              0.824   10.4       1.0
GBT (trained within the question set)           0.548      0.665   0.720    0.558              0.810   15.0       1.0
baseline: ortholog copy (grouped, abstains)     0.301        NaN   0.000      NaN                NaN    0.0       0.0
baseline: prevalence-weighted random            0.504      0.690   0.602    0.500              0.699   12.0       1.0
baseline: signal peptide alone                  0.644      0.712   0.766    0.515              0.712   14.8       1.0
baseline: training majority                     0.699      0.699   1.000      NaN                NaN   19.6       1.0
```

### easy (comparison only): all labels, internal markers included

```
metric                                       accuracy  precision  recall  roc_auc  average_precision  calls  coverage
strategy                                                                                                             
GBT (trained within the question set)           0.928      0.757   0.569    0.959              0.777   10.6       1.0
baseline: ortholog copy (grouped, abstains)     0.889        NaN   0.000      NaN                NaN    0.0       0.0
baseline: prevalence-weighted random            0.817      0.120   0.106    0.500              0.111   11.2       1.0
baseline: signal peptide alone                  0.900      0.516   0.744    0.830              0.428   18.6       1.0
baseline: training majority                     0.889        NaN   0.000      NaN                NaN    0.0       1.0
```

## Calibration

```
                        ece  ece_uncalibrated
question model                               
easy     within      0.0178            0.0379
hard     all-labels  0.1492            0.1937
         within      0.1937            0.1478
reviewed all-labels  0.1813            0.2759
         within      0.1377            0.3243
```

## Top 10 features, hard question, permutation importance on held-out folds

```
                                                    family  importance_mean
feature                                                                    
crispr_gra12s1_disco                                 gra12           0.0206
fit_naive_bmdm                           crispr_macrophage           0.0172
stress_5-1_unstress_red_24                gse132248_stress           0.0149
te245775_eif12ko_tachy_r2                        gse245775           0.0109
fit_density_dependence_log10padj   crispr_parasite_density           0.0104
fit_serum_differential_p8         crispr_serum_restriction           0.0101
codon_cai_ribosomal                     toxodb_codon_usage           0.0088
brain168465_1d_base_mean                         gse168465           0.0079
protein_ibaq_log2                             proteome_pru           0.0077
cellcycle19092_7h_r1                              gse19092           0.0068
```

## Performance by label provenance (does it do worse on the hard negatives?)

```
    positives                       negatives  n_pos  n_neg  roc_auc  average_precision  recall_at_0_5  false_positive_rate_at_0_5  median_score_pos  median_score_neg
 einar_direct                   einar_curated     46     30    0.606              0.701          0.609                       0.467             0.665             0.373
 einar_direct                   gold_negative     46     11    0.852              0.963          0.609                       0.091             0.665             0.013
 einar_direct             hyperlopit_internal     46    477    0.987              0.950          0.609                       0.000             0.665             0.000
 einar_direct einar+gold (the hard negatives)     46     43    0.685              0.684          0.609                       0.349             0.665             0.217
 einar_direct                   all negatives     46    520    0.962              0.665          0.609                       0.029             0.665             0.000
einar_curated                   einar_curated     18     30    0.554              0.446          0.556                       0.467             0.652             0.373
einar_curated                   gold_negative     18     11    0.813              0.896          0.556                       0.091             0.652             0.013
einar_curated             hyperlopit_internal     18    477    0.974              0.884          0.556                       0.000             0.652             0.000
einar_curated einar+gold (the hard negatives)     18     43    0.636              0.426          0.556                       0.349             0.652             0.217
einar_curated                   all negatives     18    520    0.946              0.403          0.556                       0.029             0.652             0.000
all positives                   einar_curated     65     30    0.586              0.751          0.585                       0.467             0.657             0.373
all positives                   gold_negative     65     11    0.838              0.970          0.585                       0.091             0.657             0.013
all positives             hyperlopit_internal     65    477    0.984              0.951          0.585                       0.000             0.657             0.000
all positives einar+gold (the hard negatives)     65     43    0.667              0.734          0.585                       0.349             0.657             0.217
all positives                   all negatives     65    520    0.958              0.711          0.585                       0.029             0.657             0.000
```

## Robustness of the hard question

```
                                                                    feature_set      model  n_features  roc_auc  average_precision  best_lower90  at_threshold  calls  correct  reaches_bar
                                                                   all features     within         398   0.7542             0.7330        0.6851          0.60     57       44        False
                                                                   all features all-labels         398   0.7755             0.7641        0.7273          0.65     45       37        False
                                                             minus foldseek_3di     within         377   0.7619             0.7491        0.6715          0.60     51       39        False
                                                             minus foldseek_3di all-labels         377   0.7678             0.7643        0.6784          0.65     45       35        False
                                                 minus sequence/topology traits     within         384   0.7509             0.7315        0.6851          0.60     57       44        False
                                                 minus sequence/topology traits all-labels         384   0.7674             0.7553        0.6246          0.70     35       26        False
minus attention-adjacent (BioID/IP-MS corpus membership, PMC, proximity counts)     within         396   0.7542             0.7330        0.6851          0.60     57       44        False
minus attention-adjacent (BioID/IP-MS corpus membership, PMC, proximity counts) all-labels         396   0.7755             0.7641        0.7273          0.65     45       37        False
                                                   minus structure availability     within         392   0.7542             0.7330        0.6851          0.60     57       44        False
                                                   minus structure availability all-labels         392   0.7751             0.7637        0.7273          0.65     45       37        False
                                             minus host-facing-adjacent screens     within         381   0.7450             0.7454        0.6717          0.65     58       44        False
                                             minus host-facing-adjacent screens all-labels         381   0.7223             0.7096        0.6770          0.70     59       45        False
                                                         minus all of the above     within         373   0.7450             0.7454        0.6717          0.65     58       44        False
                                                         minus all of the above all-labels         373   0.7216             0.7109        0.6770          0.70     59       45        False
```

## Files

| File | Contents |
|---|---|
| `scoreboard.tsv` | question x strategy x metric, mean plus every fold |
| `thresholds_hard.tsv` | the headline: threshold, calls, correct, precision, lower 90 % bound, genes fillable, reaches 0.90 |
| `thresholds_reviewed.tsv` | the same on the 98 Einar-reviewed genes |
| `thresholds_easy.tsv` | the same on all labels (comparison only) |
| `calibration.tsv` | reliability table and ECE per question and training variant |
| `importance.tsv` | permutation importance per feature with its provenance family, hard and easy questions |
| `by_label_provenance.tsv` | AUC / FPR for each positive provenance against each negative provenance |
| `by_label_provenance_errors.tsv` | per-provenance median score and error rate, every question |
| `robustness_hard.tsv` | the hard-question bound with each worrying feature block removed |
| `predictions.tsv` | all 8,140 genes: calibrated probability per question model, out-of-fold probability, label, provenance, secretory flag, per-gene occlusion attributions |


## Limitations (not softened)

1. **The hard set is 105 genes, 58 positives.** ~11 positives per fold. The Clopper-Pearson
   bound cannot reach 0.90 on fewer than ~22 calls even with perfect precision, so over most of
   the threshold range the bar is unreachable as arithmetic, before any question of model
   quality. This is a data-quantity result before it is a modelling result.

2. **Two of the hard negatives are contested** (GRA9, NTPase II), and the contested direction
   is the one that decides the question. Flip test in the notebook.
3. **The 43 hard negatives are a curation decision about published papers, not an export
   assay.** Usually immunogold EM or IFA in one strain at one time point. A protein vacuolar in
   a tachyzoite and exported in a bradyzoite would be invisible to these labels.
4. **DeepTMHMM misses signal peptides** (GRA16 / TGME49_208830 is called GLOBULAR). The
   secretory-marker arm recovers some, but a protein that is neither predicted nor a hyperLOPIT
   marker is absent from the hard question entirely. SignalP would fix it; open gap.
5. **The positives are still mostly one assay** read by one human: the PVM proximity screen
   (PMID 34749525) plus the literature.
6. **Publication bias is reduced, not removed.** The declared dependencies are held out, but a
   famous gene is also the gene with a melting temperature, a palmitome entry and a deposited
   structure. Missingness itself carries fame.
7. **`foldseek_3di` is a crude shape summary**, 637 genes have no 3Di and 520 have one copied
   by alignment; structure availability is fame-correlated.
8. **The all-gene probability columns are not independently validated** (isotonic curve fitted
   in-sample on pooled out-of-fold scores). Only the cross-fitted out-of-fold numbers are honest.
9. **The ortholog-copy baseline has coverage 0** under orthogroup-grouped folds, by construction.
10. **One model, fixed settings.** Stages 2-5 of the instruction 72 plan are not in this notebook.
11. **Permutation importance is noisy at this size** (~20 held-out genes per fold in the hard
    question); signs flip between folds. Read families, not features.
