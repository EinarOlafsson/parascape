# Function inference on the 242-gene truth set (2026-11-04)

Built by `scripts/build_function_inference_notebook.py`; the executed notebook is
`notebooks/function_inference_2026_11_04.ipynb`. Instruction 72, target trait 2 of 3.
Labels: `results/function_truth_2026_11_04/`.

## The answer

**NO for the POOLED multi-label question. No calibrated threshold from 0.20 to 0.99 reaches a lower 90 % confidence bound of 0.90. The best bound anywhere on the curve is 0.5254 (threshold 0.60: 27 calls, 18 correct, point-estimate precision 0.667, IC-weighted 0.734). Nothing may be written into the information map as a fact for the function trait.**


Non-screen-tier variant: does not reach the bar (best bound 0.6167).


Not a calibration artefact: ranking the pooled pairs by raw out-of-fold score and taking
the top k directly, for every k from 10 up - which is what any threshold rule amounts to -
the bar is cleared by no k. all-labels: top 25, 17 correct, precision 0.680, lower bound 0.5327; non-screen: top 25, 20 correct, precision 0.800, lower bound 0.6603. (`precision_at_k.tsv`)


Robustness across feature sets: 0 of 10 feature sets reach the bar; best bound ranges 0.395 to 0.550; pooled AUC ranges 0.663 to 0.792. (`robustness.tsv`)


## Labels or evidence?

WHICH CASE IS THIS? BOTH, at different places on the curve, and the two have to be separated.

(a) AT THE POOLED OPERATING POINT IT IS THE EVIDENCE. The best precision this model reaches anywhere - at any calibrated threshold and at any point on the pooled ranking - is 0.680, which is BELOW the 0.90 bar itself. The lower bound converges upward to the TRUE precision, so if this model is what we have, no number of labels brings the bound to 0.90: more labels would measure 0.680 more precisely.

(b) BUT THE MODEL IS NOT FIXED, AND THE LEARNING CURVE SAYS LABELS STILL BUY ACCURACY: from 80 % to 100 % of the labels pooled average precision moves +0.0242 and precision of the top 25 calls moves +0.1000; over the whole curve average precision goes 0.278 -> 0.365. So (a) is a statement about the model trained on 242 genes, not a ceiling on the problem. This is where the truth-set README is right: it argued the constraint for function is the label count, and the curve has not flattened.

(c) ON THE HARDEST HONEST CUT THE BOUND IS CAPPED BY ARITHMETIC, NOT BY ERRORS. Restricted to the 145 non-screen genes and the 11 non-screen-supported classes, the best threshold makes 7 calls and 7 are correct - precision 1.000 - and the lower 90 % bound is 0.720 ONLY because 7 calls cannot support more. That cut is label-limited in the plainest possible way: the bound needs about 22 calls before it can reach 0.90 at all, and 65 at a true precision of 0.95.

ACTIONABLE NUMBER. At the current call rate (0.7 % of held-out pairs called at the best threshold) a 65-call demonstration needs about 583 labelled genes against the 242 there are now - roughly 341 more. Einar's 254-gene review queue is larger than that, and its top 100 alone would take the labelled set from 242 to about 342 genes and the non-screen set from 145 to 245. That is worth doing. What it will NOT do on its own is lift a model whose precision is 0.680 to 0.90: the queue buys the measurement and some accuracy, and the rest has to come from evidence that distinguishes the classes better than expression, fitness and domain flags do.


### Learning curve (mean over draws; training labels dropped by orthogroup)

```
          train_genes  roc_auc  average_precision  precision_top25  precision_top50  best_lower90
fraction                                                                                         
0.4              94.0   0.7152             0.2776             0.46             0.44        0.3422
0.6             148.5   0.7454             0.3037             0.46             0.49        0.3483
0.8             192.5   0.7732             0.3412             0.58             0.57        0.4438
1.0             242.0   0.7900             0.3654             0.68             0.54        0.5254
```

## Why the headline is pooled

The truth-set README does the arithmetic: a Clopper-Pearson lower 90 % bound needs ~65
held-out calls even at a true precision of 0.95, and the largest non-screen class has 37
genes. So no class can be scored alone, and the headline is the pooled multi-label
question over all 3630 (gene, class) pairs (469 positive). Per-class numbers are in
`per_class.tsv` with every call count, so the under-powering is visible rather than
averaged away: 1 of 15 classes have >= 22 held-out calls at the 0.5 cut, which is the
minimum at which the bound could reach 0.90 even with perfect precision.

## Unweighted or information-content weighted?

Both, at every threshold. The **unweighted** Clopper-Pearson bound is the number the fill
decision turns on: the fill rule is a statement about cells, and only the unweighted bound
is exact for counting cells. The **IC-weighted** precision is the headline for what the
filled cells are worth - it says whether a model passes by calling specific things or by
calling `apicoplast_mitochondrion` (IC 1.60) on everything. The weighted bound uses a Kish
effective sample size and is an approximation, labelled as such in every table.

## Pooled scoreboard (pooled over held-out pairs, 0.5 cut)

```
metric                                                    precision  precision_lower90  precision_ic_weighted  recall  roc_auc  average_precision  calls  coverage
variant    strategy                                                                                                                                               
all-labels GBT one-vs-rest (all labels)                       0.479              0.427                  0.504   0.171    0.790              0.365  167.0       1.0
           baseline: best single domain column per class      0.179              0.158                  0.172   0.194    0.528              0.136  507.0       1.0
           baseline: class prior (most common class)          0.331              0.291                  0.331   0.171    0.560              0.164  242.0       1.0
           baseline: ortholog copy (grouped, abstains)          NaN                NaN                    NaN   0.000    0.500              0.129    0.0       0.0
           baseline: prevalence-weighted random               0.178              0.155                  0.153   0.177    0.528              0.138  467.0       1.0
non-screen GBT one-vs-rest (non-screen genes only)            0.519              0.452                  0.539   0.188    0.728              0.363  104.0       1.0
```

Best baseline: **baseline: class prior (most common class)** (pooled average precision 0.164). Model 0.365, improvement +0.202.


Per-fold values for every strategy and metric are in `scoreboard.tsv`.

## The pooled threshold table (the headline)

```
 threshold  calls  correct  precision  precision_lower90  precision_ic_weighted  recall  fillable_cells  fillable_genes  reaches_bar
      0.20    895      295     0.3296             0.3092                 0.3188  0.6290           23918            7769        False
      0.30    604      229     0.3791             0.3533                 0.3732  0.4883           12214            6281        False
      0.40    129       69     0.5349             0.4746                 0.5591  0.1471            2256            1879        False
      0.50     94       49     0.5213             0.4501                 0.5521  0.1045            2240            1867        False
      0.55     41       21     0.5122             0.4015                 0.5476  0.0448             294             287        False
      0.60     27       18     0.6667             0.5254                 0.7336  0.0384             294             287        False
      0.65     21       12     0.5714             0.4095                 0.6355  0.0256             292             285        False
      0.70      8        3     0.3750             0.1469                 0.3888  0.0064              61              60        False
      0.75      6        2     0.3333             0.0926                 0.3213  0.0043              46              45        False
      0.80      1        0     0.0000             0.0000                 0.0000  0.0000              36              36        False
      0.85      0        0        NaN                NaN                    NaN  0.0000              30              30        False
      0.90      0        0        NaN                NaN                    NaN  0.0000              22              22        False
      0.95      0        0        NaN                NaN                    NaN  0.0000              13              13        False
      0.97      0        0        NaN                NaN                    NaN  0.0000               8               8        False
      0.99      0        0        NaN                NaN                    NaN  0.0000               5               5        False
```

## Per class

```
                          cls  class_genes  ic_corpus  screen_supported  calls  correct  precision  precision_lower90  recall  roc_auc  average_precision  reaches_bar
     apicoplast_mitochondrion           80       1.60              True     33       13      0.394              0.278   0.162    0.654              0.427        False
          cell_cycle_division           62       1.96              True     13        6      0.462              0.264   0.097    0.678              0.370        False
                 cytoskeleton           58       2.06              True     17       11      0.647              0.463   0.190    0.704              0.452        False
                       egress           42       2.53              True      7        4      0.571              0.279   0.095    0.695              0.302        False
            invasion_motility           37       2.71             False      4        1      0.250              0.026   0.027    0.700              0.313        False
protein_trafficking_secretion           32       2.92             False      2        2      1.000              0.316   0.062    0.745              0.376        False
        transport_homeostasis           28       3.11             False      9        6      0.667              0.401   0.214    0.835              0.545        False
                   metabolism           26       3.22             False      1        1      1.000              0.100   0.038    0.759              0.281        False
                 proteostasis           23       3.40             False      0        0        NaN                NaN   0.000    0.656              0.232        False
                   signalling           20       3.60             False      4        2      0.500              0.143   0.100    0.717              0.336        False
       lipid_membrane_biology           19       3.67             False      0        0        NaN                NaN   0.000    0.456              0.077        False
         host_cell_modulation           16       3.92             False      0        0        NaN                NaN   0.000    0.780              0.335        False
       stress_differentiation            9       4.75             False      0        0        NaN                NaN   0.000    0.668              0.147        False
                  translation            9       4.75             False      4        3      0.750              0.320   0.333    0.809              0.544        False
      transcription_chromatin            8       4.92             False      0        0        NaN                NaN   0.000    0.634              0.309        False
```

## CAFA no-knowledge split

**Definition.** A gene is *no-knowledge* when it has no prior functional annotation of any
kind in the feature set: (1) `domain_novel` is True - no InterPro and no Pfam entry, so all
313 domain columns are False - and (2) its ToxoDB `product` is "hypothetical protein" or an
unspecified product. A looser definition (domain-free only, any name) is reported beside it.
Strict: 20 of the 242 labelled genes. Loose: 40. Genome-wide the strict definition covers
2664 of 8140 genes, which is the population the map most needs filled.


no-knowledge (strict) 20 genes, 300 pooled pairs, 38 positives; pooled AUC 0.898 against 0.771 on the 222 annotated genes. Calls at 0.5: 10, correct 5.

```
                                          split  genes  threshold  pairs  positives  calls  correct  precision  precision_lower90  precision_ic_weighted  recall  roc_auc  average_precision  reaches_bar
       no-knowledge strict (no domain, no name)     20        0.5    300         38     10        5      0.500              0.267                  0.525   0.132    0.898              0.433        False
       no-knowledge strict (no domain, no name)     20        0.7    300         38      1        0      0.000              0.000                  0.000   0.000    0.898              0.433        False
       no-knowledge strict (no domain, no name)     20        0.9    300         38      0        0        NaN                NaN                    NaN   0.000    0.898              0.433        False
       no-knowledge strict (no domain, no name)     20        NaN    300         38     61       29        NaN              0.387                    NaN     NaN    0.898              0.433        False
                 no-knowledge loose (no domain)     40        0.5    600         66     17        8      0.471              0.297                  0.470   0.121    0.850              0.349        False
                 no-knowledge loose (no domain)     40        0.7    600         66      1        0      0.000              0.000                  0.000   0.000    0.850              0.349        False
                 no-knowledge loose (no domain)     40        0.9    600         66      0        0        NaN                NaN                    NaN   0.000    0.850              0.349        False
                 no-knowledge loose (no domain)     40        NaN    600         66    105       43        NaN              0.345                    NaN     NaN    0.850              0.349        False
has prior annotation (the complement of strict)    222        0.5   3330        431     84       44      0.524              0.448                  0.555   0.102    0.771              0.344        False
has prior annotation (the complement of strict)    222        0.7   3330        431      7        3      0.429              0.170                  0.435   0.007    0.771              0.344        False
has prior annotation (the complement of strict)    222        0.9   3330        431      0        0        NaN                NaN                    NaN   0.000    0.771              0.344        False
has prior annotation (the complement of strict)    222        NaN   3330        431     23       17        NaN              0.587                    NaN     NaN    0.771              0.344        False
                             all labelled genes    242        0.5   3630        469     94       49      0.521              0.450                  0.552   0.104    0.784              0.350        False
                             all labelled genes    242        0.7   3630        469      8        3      0.375              0.147                  0.389   0.006    0.784              0.350        False
                             all labelled genes    242        0.9   3630        469      0        0        NaN                NaN                    NaN   0.000    0.784              0.350        False
                             all labelled genes    242        NaN   3630        469     27       18        NaN              0.525                    NaN     NaN    0.784              0.350        False
```

## Tier split, and the function analogue of the easy contrast

screen-only genes: pooled AUC 0.893, best bound 0.4668 over 65 calls. Non-screen genes: pooled AUC 0.701, best bound 0.6402. The four screen-supported classes: AUC 0.697; the other eleven: AUC 0.755. Non-screen genes restricted to the other eleven classes (the hardest honest cut): AUC 0.720, best bound 0.7197.

```
   variant                                                   split  genes  classes  pairs  positives  calls  correct  precision  precision_lower90  roc_auc  average_precision  best_lower90  reaches_bar
all-labels non-screen genes (experimental GO + curated Swiss-Prot)    145       15   2175        288     47       23      0.489              0.387    0.701              0.303         0.640        False
all-labels                                       screen-only genes     97       15   1455        181     47       26      0.553              0.449    0.893              0.434         0.467        False
all-labels               all genes, the 4 screen-supported classes    242        4    968        242     70       34      0.486              0.403    0.697              0.387         0.424        False
all-labels                         all genes, the other 11 classes    242       11   2662        227     24       15      0.625              0.474    0.755              0.313         0.720        False
all-labels                  non-screen genes, the other 11 classes    145       11   1595        227     20       15      0.750              0.585    0.720              0.376         0.720        False
non-screen non-screen genes (experimental GO + curated Swiss-Prot)    145       15   2175        288     53       34      0.642              0.545    0.713              0.336         0.617        False
non-screen               all genes, the 4 screen-supported classes    145        4    580         61      9        7      0.778              0.510    0.692              0.300         0.562        False
non-screen                         all genes, the other 11 classes    145       11   1595        227     44       27      0.614              0.506    0.715              0.349         0.552        False
non-screen                  non-screen genes, the other 11 classes    145       11   1595        227     44       27      0.614              0.506    0.715              0.349         0.552        False
```

## Method

* **Labels**: 242 genes x 15 classes = 3630 (gene, class) pairs, 469 positive. Unlabelled genes
  are never negatives. A negative pair means "not recorded", not "not true".
* **Features**: 732, in 86 families - numeric/boolean `nodes.parquet` columns, the DeepTMHMM
  sequence/topology traits, the 313 domain columns, and 35 co-expression neighbour-label
  features. Every feature with its family is in `features.tsv` and printed in the notebook.
* **Held out**: the whole `splitcas9_imaging_screen` family (the 6 `screen_*` columns),
  dropped THROUGHOUT rather than conditionally, because those four phenotype columns are the
  source of the 97 screen-tier labels; every `lopit|compartment|ortholopit|tagm` column
  (localization is a tier-1 conclusion and restates some classes); `hosttx_effectors`,
  `pvm_proximity_positive`, `n_host_targets`; and every free-text or identity column
  (`sequence`, `product`, `interpro_id`, `interpro_desc`, `pfam_id`, `ec_number`,
  `orthogroup`, ...). 73 columns held out in total.
* **Folds**: 5, grouped by OrthoMCL orthogroup so paralogs never split, seed 0, assigned once
  over the 242 genes and shared by all 15 class models.
* **Model**: one `HistGradientBoostingClassifier` per class per fold (one-vs-rest), missing
  values passed through with no imputation, 300 rounds at learning rate 0.05, 15 leaves,
  early stopping off, no hyperparameter search. Two variants: `all-labels` (242 genes) and
  `non-screen` (145 genes, scored on 145).
* **Co-expression features are built inside the fold loop**: partner labels come only from
  the training folds, so no test gene sees a same-fold partner label.
* **Calibration**: isotonic (`parascape.claims._isotonic`), cross-fitted over the pooled
  pairs. ECE in `calibration.tsv`.
* **Bound**: exact one-sided Clopper-Pearson, `scipy.stats.beta.ppf(0.10, k, n-k+1)`.

## Calibration

```
               ece  ece_uncalibrated
variant                             
all-labels  0.0100            0.0847
non-screen  0.0134            0.0829
```

## Top 10 features by permutation importance (held-out folds, pooled average precision)

```
                  feature                    family  importance_mean  importance_sd
 fit_complete_medium_2025 crispr_glucose_limitation           0.0143         0.0015
                   has_ec         toxodb_ec_numbers           0.0124         0.0023
                     n_tm           sequence_traits           0.0119         0.0015
            dom_IPR000719                    domain           0.0100         0.0007
h4_acetylation_chip_score     toxodb_h4_acetylation           0.0080         0.0015
               tm_kd_mean              unregistered           0.0062         0.0010
   melting_temperature_tm       melting_temperature           0.0062         0.0008
               n_fulltext                    pmc_oa           0.0059         0.0029
                cterm_len           sequence_traits           0.0052         0.0018
          expr_sporulated                 gse206344           0.0044         0.0007
```


By family, which is the unit of evidence:

```
                     family  n_features  importance_mean  importance_sd
            sequence_traits          14           0.0239         0.0039
                   gse19092          28           0.0164         0.0039
  crispr_glucose_limitation           5           0.0138         0.0045
               unregistered          12           0.0127         0.0026
          toxodb_ec_numbers           1           0.0113         0.0014
                  gse206344           7           0.0111         0.0016
   crispr_serum_restriction           6           0.0093         0.0066
                     domain         313           0.0087         0.0017
                     pmc_oa           1           0.0073         0.0020
        melting_temperature           1           0.0071         0.0010
               coexpression          35           0.0059         0.0035
             crispr_invitro           1           0.0055         0.0054
                  gse108740          15           0.0051         0.0028
                   gse51780           6           0.0050         0.0020
                  gse168465          10           0.0046         0.0036
      toxodb_h4_acetylation           1           0.0044         0.0017
           gse277553_cuttag           1           0.0038         0.0018
          gra17_synthlethal           3           0.0038         0.0035
 invivo_brain_transcriptome          14           0.0027         0.0014
             gse313048_atac           1           0.0026         0.0010
  organelle_surface_turboid           6           0.0026         0.0017
                   gse99395          12           0.0022         0.0061
             proteome_total          15           0.0021         0.0018
                   interpro           1           0.0020         0.0024
    crispr_parasite_density           5           0.0019         0.0014
        thermal_shift_cetsa           1           0.0017         0.0006
                     pubmed           4           0.0013         0.0023
          crispr_macrophage           2           0.0012         0.0013
               oocyst_itraq           8           0.0011         0.0021
         toxodb_strain_snps           5           0.0011         0.0016
            antisense_level           1           0.0010         0.0004
                  alphafold           1           0.0010         0.0023
                  gse200962          16           0.0008         0.0021
             morc_depletion          18           0.0008         0.0008
    crispr_invivo_composite           6           0.0007         0.0015
            foldseek_struct           1           0.0005         0.0002
   toxodb_nanopore_isoforms           1           0.0005         0.0003
          pride_acetylation           1           0.0002         0.0009
         toxodb_codon_usage           3           0.0002         0.0014
         gse223620_bfd2_rip           1           0.0001         0.0008
             xue_singlecell           1           0.0001         0.0004
    toxodb_enteroepithelial           1           0.0000         0.0004
              starpath_xlms           2           0.0000         0.0002
     ipms_corpus_membership           1           0.0000         0.0001
                  m6a_peaks           1           0.0000         0.0001
            toxodb_epitopes           1           0.0000         0.0001
     differentiation_screen           1           0.0000         0.0000
            pride_proximity           1           0.0000         0.0000
toxodb_arginine_methylation           1           0.0000         0.0000
  second_background_fitness           1           0.0000         0.0000
               myristoylome           1           0.0000         0.0000
           crispr_young2019           1           0.0000         0.0000
      cyst_wall_interactome           2           0.0000         0.0000
 curated_resistance_alleles           2           0.0000         0.0000
                 iedb_bcell           1           0.0000         0.0000
             mrna_stability           1           0.0000         0.0000
                 ipms_baits           1           0.0000         0.0000
        secretome_partition           1          -0.0000         0.0000
          pride_lactylation           1          -0.0000         0.0001
            invivo_platform           1          -0.0000         0.0000
        pride_glycosylation           1          -0.0000         0.0000
       pride_ubiquitination           1          -0.0000         0.0000
      crosslink_interactome           1          -0.0000         0.0003
             stage_enriched           1          -0.0001         0.0001
               proteome_pru           1          -0.0001         0.0002
          toxodb_macrophage           1          -0.0001         0.0000
          gse302107_riboseq           2          -0.0002         0.0003
                  local_af3           9          -0.0002         0.0002
           toxodb_palmitome           2          -0.0003         0.0002
           cdpk1_substrates           1          -0.0003         0.0001
    bioid_corpus_membership           1          -0.0004         0.0005
               phosphosites           2          -0.0004         0.0006
        pride_nitrosylation           1          -0.0004         0.0004
              sexual_stages           1          -0.0004         0.0004
                      gra12           6          -0.0004         0.0003
       phospho_quantitative           5          -0.0006         0.0005
                  gse168155           8          -0.0008         0.0008
       mrna_decay_gse329845           1          -0.0009         0.0010
    iron_depletion_proteome           3          -0.0014         0.0046
                   gse22258           2          -0.0020         0.0016
    oxidative_stress_screen           1          -0.0023         0.0010
                  gse245775          36          -0.0024         0.0024
             brady_subtypes           5          -0.0024         0.0021
           gse132248_stress          10          -0.0032         0.0012
                  gse129869          18          -0.0035         0.0014
             gse302108_utr5           6          -0.0049         0.0020
```

## Files

| File | Contents |
|---|---|
| `scoreboard.tsv` | strategy x variant x metric, pooled and per fold |
| `thresholds_pooled.tsv` | the headline: threshold, calls, correct, precision, lower 90 % bound, IC-weighted precision, cells fillable, reaches 0.90 |
| `per_class.tsv` | the same per class with its call count, at 0.5 and at the pooled operating threshold |
| `per_class_label_need.tsv` | labelled genes each class would need before its own bound could be demonstrated |
| `no_knowledge.tsv` | the CAFA no-knowledge split, both definitions, against the annotated complement |
| `by_tier.tsv` | non-screen vs screen-only genes, and screen-supported vs other classes |
| `calibration.tsv` | reliability table and ECE over the pooled pairs |
| `importance.tsv` | permutation importance per feature and per family |
| `features.tsv` | every feature with its family |
| `robustness.tsv` | the bound with each feature family removed, and with each one alone |
| `learning_curve.tsv` | pooled metrics when trained on 40/60/80/100 % of the labels, dropped by orthogroup |
| `precision_at_k.tsv` | precision and bound of the top k pooled pairs, every k from 10 up |
| `predictions.tsv` | all 8,140 genes x 15 classes: calibrated probability, out-of-fold probability for labelled genes, label, tier, fold, no-knowledge flags |


## Limitations (not softened)

1. **Negative pairs are "not recorded", not "not true".** Open-world annotation: a correct
   prediction of an unrecorded function scores as a false positive, so precision is biased
   downward. A passing bound is trustworthy; a failing bound may be pessimistic.
2. **The classes are coarse and of unequal specificity** (IC 1.60 to 4.92); the pooled number
   is dominated by the big vague classes, which is why the IC-weighted figure sits beside it.
3. **The labels come from regular expressions over curator text and GO descendant sets.** The
   domain columns read the same annotation act. Not label leakage, but not independent
   evidence either; the no-knowledge split and the annotation-free feature set measure it.
4. **The screen tier is a different experiment** with four process classes and a non-random
   319-gene list. Its four source columns are dropped throughout; the gene selection cannot be.
5. **242 genes, 145 at the non-screen tiers.** Per class this is too small to score at all.
6. **One-vs-rest ignores the structure between classes**; no joint model here.
7. **Co-expression features are fold-safe but not study-held-out** - the same expression
   studies feed the expression family.
8. **The all-gene probability columns are not independently validated** (isotonic curve fitted
   in-sample on pooled out-of-fold pairs). Only the out-of-fold numbers are honest.
9. **The ortholog-copy baseline has coverage 0.00** under orthogroup-grouped folds, by
   construction.
10. **One model, fixed settings, no search.** Stages 2-5 of the instruction 72 plan are absent.
11. **Permutation importance is noisy** at ~48 held-out genes per fold; read families.
