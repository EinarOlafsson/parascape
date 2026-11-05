# Expression tier-2 columns of the information map (2026-11-05)

Built by `scripts/build_expression_traits_notebook.py`; notebook `notebooks/expression_traits_2026_11_05.ipynb`.
Snapshot `20261105T124050Z-498` (498 verified runs of 2227 planned).
Values come from the uniform reprocessing only: salmon, ME49_host reference, parasite genes, TPM reclosed
to a million within the parasite genes, `log2(TPM + 1)`. Author-deposited matrices are not mixed in.

## The rules, verbatim

* **Normalization preference (instruction 72, Einar):** `relative_to_control` > `internal_reference` >
  `library_size` > `raw`. Where both a relative and an absolute form of one quantity exist, learners take
  the relative one; absolute values are used only where no control exists.
* **A control group** is a condition group in which *every* run is marked `control (source label)` or
  `candidate baseline; review required` in `sample_design.tsv`. Every other condition group of that study
  becomes one contrast against the pooled control runs of the same study.
* **`control_link` carries no pairing.** For all 450 kept runs its value is the sentence "No verified
  pairing; choose the matching control within this study", and the snapshot records `control_pairing:
  unverified; requires within-study review`. The pooled baseline is therefore the control-marked *part of a
  study*, not the matched arm, and in at least one study (GSE158231) that pool mixes drug-treated and
  untreated baselines. Read `control_groups` and `contrast_label` in `column_registry.tsv` before using a column.
* **Internal reference rule:** group weighted mean `log2(TPM + 1)` minus the median of that gene across the
  condition groups of the same study. It removes gene length, composition and study batch; it does not
  remove a shift shared by every condition of the study. No housekeeping reference was used: TPM is already
  reclosed within the parasite genes and no invariant gene set is validated for these conditions.
* **Absolute rule:** weighted mean `log2(TPM + 1)` over all kept runs of the study. Last resort.
* **Weighting:** run weight = its quality score (`results/rna_sample_quality_2026_11_04/`). A group mean is
  weighted; `effective_n = (sum w)^2 / sum w^2` (Kish); `mean_quality` is the weighted mean of those weights.
  A run with no score takes its study median quality, or the global median (0.898) if its study has
  none scored. Both medians are taken over scored runs that are not themselves excluded, so a study whose
  only scored runs are the zero-parasite ones (GSE245775) does not inherit their weight.
  `n_quality_imputed` counts imputed-weight runs per column.
* **Standard error:** within-group weighted variances pooled over the groups of the study with >= 2 runs,
  moderated towards the study median variance with a prior worth 4 degrees of freedom, then
  `se = sqrt(s2_mod * (1/n_eff_treatment + 1/n_eff_control))`. Studies with no replicated group fall back to
  the variance across all their runs and are flagged `se_approximate`.
* **Shrinkage:** `lfc = lfc_raw * tau2 / (tau2 + se^2)` with `tau2 = var(lfc_raw) - mean(se^2)` over the genes
  of that contrast, floored just above zero (normal-normal posterior mean).
* **Probability:** `base x mean_quality x replication`, clipped to [0.05, 0.95]. `base` = 0.90
  (`relative_to_control`), 0.85 (`internal_reference`), 0.80 (`library_size`); `replication` = 1.0 when the
  treatment arm has `effective_n >= 2`, else 0.85. A relative column does not get the 0.95 of a direct
  experiment because its control pairing is unverified.
* **Traits:** every column is registered with `information_map.register`, tier 2, kind `measure`, role
  `measurement`, `family` = the study accession, and the normalization above. `parascape/information_map.py`
  was not edited; the lines are in `register_lines.txt`.
* **Missing is empty, never zero.** A gene with no value in a study was not measured there.

## Exclusions

| exclusion | runs | detail |
|---|---|---|
| design rows outside the studies that have a matrix | 1737 | runs planned or verified but with no salmon ME49_host matrix in this snapshot |
| single-cell collections and the GSE229505 quarantine | 0 | instruction 69: per-cell 10x libraries; 4 quarantined runs on the NAS |
| zero parasite reads (untrimmed short-insert libraries) | 40 | studies GSE243206, GSE243207, GSE245775 |

Single-cell collections are out by instruction 69 (per-cell 10x libraries; deposited single-cell matrices are
used instead), including the 4 GSE229505 runs quarantined under
`parascape_omics/rnaseq_excluded_single_cell/`. Runs with zero parasite reads are the untrimmed
short-insert libraries of GSE243206, GSE243207 and GSE245775; they return after the planned fastp
re-quantification. GSE243206 and GSE243207 lose every run and so contribute no column at all.

## What was built

| normalization | quantity | columns |
|---|---|---|
| internal_reference | log2_vs_study_median | 106 |
| library_size | mean_log2_tpm | 36 |
| relative_to_control | log2_fold_change | 62 |
| relative_to_control | standard_error | 62 |

17 of 36 usable studies have a control group; 19 do not.
Coverage: 8140 of 8140 genes carry at least one value;
the gene x column grid is 97.7 % filled.

## How comparable are these columns, honestly

**Columns with no signal.** 34 of 62 fold-change columns move fewer than 50 genes by at least
half a log2 unit (`genes_abs_lfc_ge_0_5` in the registry). Eight of them are GSE132238, which deposits one run
per condition group: with no replicated group its standard error falls back to about 3.9 log2 units and every
fold change is shrunk to zero. That is the correct answer for an unreplicated design, and those columns must
not be read as "nothing changed in the experiment".

**Two contrasts the control marking hides, with the numbers.** In GSE218921 the keyword rule marked the
wild-type alkaline-switch arm as a candidate baseline, so the best alkaline contrast in the snapshot never
became a column. Computed directly from the group means it is exactly as expected: BFD1 +1.36, BAG1 +4.75,
ENO1 +4.91, SAG1 -1.34. The same thing happens to WT_stress in GSE132238 (BAG1 +6.92, ENO1 +6.99, SAG1 -5.25).
The data hold the biology; the control marking loses it. These rows are in `known_biology.tsv` under
`hidden by the control marking` and are the clearest argument for a within-study control review.

**Known biology** (`known_biology.tsv`): 30 pass,
7 weak pass,
5 inconclusive,
3 miss. Every miss is listed in the table with its value.

**Replicate agreement** (`replicate_agreement.tsv`): of 62
fold-change columns, 26 have at least two runs in both arms and can be split in half; median
split-half Spearman 0.615, and 8 fall
below 0.3. A low split-half correlation usually means the perturbation changed little, so the column is mostly
noise -- correct arithmetic, no evidence. Those columns should not be used on their own.

**Cross-study agreement** (`cross_study_agreement.tsv`, `cross_study_named_pairs.tsv`): over
932 same-normalization pairs inside a keyword theme the median Spearman is only
0.036. **That number must not be read as "the columns are incomparable"**: most
pairs inside a theme perturb different genes and have no reason to agree. The pairs that should agree do:
BFD1 stabilised in two independent studies (GSE139537 against GSE223877) reaches
0.856, BFD1 lost in two
independent studies (GSE132238 against GSE223819) reaches
0.452, and the
opposite perturbations correlate negatively
(-0.545), while
Shield-1 given to a line with no BFD1 to stabilise -- the control on the test itself -- correlates
-0.089, i.e. not at all.
Where it does **not** work is alkaline stress: two independent alkaline-stress studies agree only at
0.236. Both of
those are `internal_reference` columns whose centre is a different mixture of conditions in each study, so
part of that disagreement is the normalization and not the biology -- which is itself the argument for
getting controls reviewed so those studies can carry `relative_to_control` columns instead.

**Absolute profiles** agree well: median pairwise Spearman 0.882 over 630 study
pairs. Studies do agree about which genes are highly expressed; the hard part is the differences.

**Summary judgement.** A perturbation column is trustworthy where the perturbation is strong and the control
arm is the right one -- BFD1 gain and loss reproduce across labs on these values, which is what the uniform
reprocessing was for. Columns from weak perturbations, from studies with no replicated group, and from the
studies whose baseline mixes conditions are not yet evidence on their own. The binding constraint is not the
quantification; it is the **control assignment**, which no code here can fix.
