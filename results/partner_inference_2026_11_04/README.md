# Binding-partner inference, first scoring run (instruction 72, trait 2)
Notebook: `notebooks/partner_inference_2026_11_04.ipynb`. Builder: `scripts/build_partner_inference_notebook.py`. Seed 0.
## Question
Positives: the **2842 StarPath XL-MS crosslink pairs** over 1630 proteins, probability 1 - local FDR (median 0.9996). Direct physical evidence.
The 33596 bait-member pull-down pairs in the information map are **not** headline positives: they sit at a provisional 0.5 because no Toxoplasma background table exists, and without a background frequency a sticky contaminant is indistinguishable from a real partner. They are scored separately (§12).
**There is no verified non-interaction set.** Negatives are decoys — pairs with no recorded evidence — built three ways, and the headline is the hardest: for each positive, one endpoint kept and the other replaced by a protein of similar crosslink degree from the same held-out protein block, which removes the popular-protein shortcut and the "detectable by XL-MS at all" shortcut.
## Scoreboard (lower 90 % Clopper-Pearson bound on held-out precision)
```
                          question  bound  reaches    n
HEADLINE cold-node, degree-matched 0.8907    False 1248
cold-node, same-compartment decoys 0.9192     True 1199
    easy: cold-node, random decoys 0.9022     True 1185
  easy: known-node, degree-matched 0.9120     True 1959
```

**Verdict: partner prediction does NOT reach the 0.90 bar on the cold-node, degree-matched question.** Best lower bound **0.891**.
Pairs that may be written as facts: **0**.
## Negative-set and split gap (the pair analogue of the easy contrast)
```
regime                                       cold_node                                                      known_node                         
negatives                               degree_matched pulldown_degree_matched  random same_compartment degree_matched  random same_compartment
strategy                                                                                                                                       
baseline: both popular (degree product)         0.5122                     NaN  0.9954           0.9775         0.5433  0.9842           0.9653
baseline: random                                0.4979                     NaN  0.5117           0.4961         0.4978  0.5088           0.5005
baseline: same compartment alone                0.7050                     NaN  0.7366           0.2461         0.6916  0.7246           0.2331
baseline: shared orthogroup alone               0.5173                     NaN  0.5173           0.5078         0.5036  0.5030           0.4997
model: all families                             0.7879                   0.643  0.9328           0.9381         0.8002  0.9546           0.9526
```

## Feature families, hard question only
```
               family  n_features  n_scored  auc_alone  ap_alone  auc_without  coverage  auc_drop_when_removed
compartment_predicted           3      1248     0.7047    0.6954       0.7477    1.0000                 0.0210
             sequence           8      1248     0.6596    0.6509       0.7537    1.0000                 0.0150
 compartment_observed           4      1248     0.6589    0.6544       0.7686    1.0000                 0.0001
               domain           5      1248     0.5813    0.6146       0.7655    1.0000                 0.0032
              fitness          60      1248     0.5481    0.5479       0.7688    1.0000                -0.0001
         coexpression           4      1248     0.5395    0.5715       0.7741    1.0000                -0.0054
          host_facing           3      1248     0.5317    0.5255       0.7746    0.6363                -0.0059
           orthogroup           2      1248     0.5164    0.5142       0.7725    1.0000                -0.0038
        structure_3di           4      1248     0.5021    0.5220       0.7699    1.0000                -0.0012
```

## Top 10 features by permutation importance (held-out folds, hard question)
```
                        feature                family  importance
          same_comp_pred_coarse compartment_predicted      0.0975
                    dom_jaccard                domain      0.0216
                  same_comp_obs  compartment_observed      0.0113
                 same_comp_pred compartment_predicted      0.0073
           signal_peptide_agree              sequence      0.0065
           same_comp_obs_coarse  compartment_observed      0.0056
            both_signal_peptide              sequence      0.0043
                tdi_kmer_cosine         structure_3di      0.0029
       fitdiff_fit_density_high               fitness      0.0021
fitmin_fit_complete_medium_2025               fitness      0.0016
```

## Known complexes
```
        complex        pair  crosslink_positive  in_pulldown  probability  percentile_among_candidates
MyoA glideosome   MyoA-MLC1                True        False          NaN                          NaN
MyoA glideosome   MyoA-ELC1                True        False       0.8361                       0.9951
MyoA glideosome  MyoA-GAP45                True        False          NaN                          NaN
MyoA glideosome  MyoA-GAP50               False        False       0.8197                       0.8654
MyoA glideosome  MyoA-GAP40               False        False       0.8197                       0.8654
MyoA glideosome   MLC1-ELC1                True        False       0.9464                       0.9963
MyoA glideosome  MLC1-GAP45                True        False       0.7576                       0.7449
MyoA glideosome  MLC1-GAP50                True        False       0.9275                       0.9963
MyoA glideosome  MLC1-GAP40                True        False          NaN                          NaN
MyoA glideosome  ELC1-GAP45                True        False       0.9275                       0.9963
MyoA glideosome  ELC1-GAP50               False        False       0.8197                       0.8654
MyoA glideosome  ELC1-GAP40               False        False       0.8197                       0.8654
MyoA glideosome GAP45-GAP50                True        False       0.9275                       0.9963
MyoA glideosome GAP45-GAP40               False        False       0.8197                       0.8654
MyoA glideosome GAP50-GAP40               False        False       0.9737                       0.9963
      AMA1-RON2   AMA1-RON2               False        False       0.2667                       0.0053
      MIC2-M2AP   MIC2-M2AP                True        False          NaN                          NaN
```

## Calibration
Expected calibration error, cross-fitted isotonic:

```
       negatives     regime    n    ECE
  degree_matched  cold_node 1248 0.0524
  degree_matched known_node 1959 0.0254
          random  cold_node 1185 0.0256
          random known_node 1019 0.0112
same_compartment  cold_node 1199 0.0178
same_compartment known_node 1353 0.0141
```

## What the signal is, and what to check
* **Co-compartment carries most of it.** The strongest family on the hard question is the predicted compartment (AUC alone 0.705) with the observed marker compartment close behind (0.659); permutation importance is dominated by `same_comp_pred_coarse`. Crosslink positives share a compartment in 51 % of pairs against 10 % of degree-matched decoys. Two things follow: the model is largely a co-localisation predictor, and **whether the pair FDR used by StarPath relied on localisation agreement must be checked** — if it did, the strongest feature is partly circular. The crosslink table does carry `loc_a`/`loc_b` columns, which is why they are held out here.
* **The CRISPR-fitness pair features do not reproduce the StarPath finding.** The fitness family — the absolute difference and the minimum of each of the 28 `fit_*` columns, plus the profile correlation — reaches AUC 0.548 alone and removing it changes the full model by +0.0001. Instruction 72 records the difference as highly predictive; on this question, with degree-matched decoys and cold-node folds, it is not. The mouse-immune-dependence note in 72 found the same thing for single genes: subtracting two relative scores amplifies noise.
* **The known-complex check is mixed and one miss is flat.** 6 of the 12 glideosome pairs that could be scored sit above the 90th percentile of candidates, but **AMA1-RON2 is placed at the 0.5th percentile** — the model puts one of the best-established complexes in Toxoplasma near the bottom, because AMA1 (microneme) and RON2 (rhoptry neck) are in different compartments and the model is mostly a co-localisation predictor. That is the clearest single demonstration of what it has and has not learned.
## What constrains the answer
* Best precision at any threshold with at least 20 held-out calls: **0.957**. A 0.90 lower bound at that precision needs ~52 held-out calls, against 47 scored pairs in the headline question.
* Cold-node scoring keeps only 22 % of the labelled pairs (1248 of 5684): a pair is testable only when both endpoints land in the same held-out protein block. That is inherent to the regime, not a bug.
* The decoys are assumed non-interacting. Some fraction of them interact and were never reported, which depresses measured precision by an unknown amount. Nothing in the data bounds that fraction.
## What would change it
* **A Toxoplasma background table** (instruction 71, product 4) is the obvious candidate and the largest single gain available: it would turn 33596 provisional bait-member pairs into graded co-purification evidence. The positives would grow from 2842 to the specific subset of those 33596 — and, more important for the scoring, background frequency-filtered co-purification is the one *independent* assay that could serve as held-out truth for crosslink-trained predictions, which nothing in this notebook has. How much it would add cannot be stated as a number before the background exists; what can be stated is the ceiling: 33596 pairs are currently unusable as labels for want of it.
* Reciprocal pull-downs or solved complexes would give the first **verified negatives**; every number here rests on decoys instead.
* AF3 pair confidence on a shortlist is the per-pair evidence the feature table lacks entirely: every feature here is a property of the two proteins separately, combined. The crosslink-derived Chai models exist but are held out as label-derived.
## Files

| file | contents |
|---|---|
| `scoreboard.tsv` | every strategy x negative set x split regime, per-fold spread |
| `thresholds.tsv` | threshold sweep with the Clopper-Pearson lower 90 % bound |
| `verdicts.tsv` | the lowest threshold reaching the bar, per question |
| `calibration.tsv` | reliability bins and expected calibration error |
| `importance.tsv` | permutation importance per feature with its family |
| `family_ablation.tsv` | each family alone and held out, hard question only |
| `known_complexes.tsv` | the MyoA glideosome, AMA1-RON2, MIC2-M2AP pair by pair |
| `predictions.tsv` | top 500 candidate pairs with features and fillable/exploratory status |
| `fill_counts.tsv` | candidate pairs above each threshold |
| `features.tsv` | every feature with its family |
