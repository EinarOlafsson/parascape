% ParaScape: one gene space for the apicomplexan literature% 18 November 2026## Figures

**Figure 1. What ParaScape is, and how much of the genome it reaches**

![](results/figure1_schematic_2026_11_18/panel_schematic.pdf){width=6.5in}

**(a)** 177 datasets imported into one gene space of 8,140 genes and 362 slots, made comparable across studies, and read through 24 modules. Counts are read from the dataset catalogue at draw time. **(b)** The share of genes anything has ever measured, per category of knowledge. Transcript reaches 99.9% and chromatin 98.5%, while function reaches 16.1%, phenotype 6.9% and host response 5.3%. 51 columns are set aside as derivations, flags or counts in which absence was written as zero.


**Figure 2. Which studies can be compared at all**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

Agreement between the untreated arms of 12 independent studies, over 8,157 genes. Removing the cell-cycle axis raises the median agreement from 0.832 to 0.887. A random axis of the same size and shape gains +0.018 (sd 0.005, 200 draws), so the cell-cycle excess is +0.038 at p = 0.0050. A study whose control profile carries one value is refused rather than scored.


**Figure 3. The co-expression modules, and what each one is**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

8,018 genes cut into 24 modules on co-expression components from 57 independent studies. 14 modules hold a compartment above their own matched baseline; the rest hold none, which is reported rather than forced. The baseline is a random gene set matched to each module in size and in how much is already known about its members (1000 draws), never the prior size of the compartment.


**Figure 4. The laboratory structures the corpus, not the perturbation**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

20,100 pairs of stored contrasts. Two from the same study agree at a median absolute Spearman of 0.284, two from different studies at 0.085 -- 3.3 times. Of 178 parsed contrasts over 19 studies, 1 perturbation is stated by two independent laboratories, so cross-laboratory replication cannot be asked of this corpus.


**Figure 5. The bradyzoite subtype nothing published reaches**

![](results/bradyzoite_reach_2026_11_18/panel_reach_by_subtype.pdf){width=6.5in}

Five subtypes of in vivo bradyzoite tested against all 201 stored within-study contrasts, 1,000 expression-matched random gene sets behind each. 135 contrasts reach at least one subtype. Subtype C is reached by 4, against 86 for the most reached. Each subtype is marked by 150 genes.


**Figure 6. Which measurement explains which other measurement**

![](results/coexpression_components_localization_2026_11_06/panel_which_evidence_places_a_protein.pdf){width=6.5in}

**(a)** The 13 coarse hyperLOPIT compartments recovered from RNA co-expression alone against all allowed evidence, as cross-validated AUC. Five are placed as well or better by transcription alone, and each is built in a wave during division. **(b)** What sets protein abundance in *Plasmodium*: transcript abundance explains 48.6% of the variance on the 439 genes carrying every layer, and every layer together 50.5%. On the wider 695-gene set transcript explains 24.3%. Comparable with the human literature, one averaged mRNA measure against protein gives r2 = 0.396 over 2,611 genes. **(c)** Every pair of categories: 11 of 14 recover more than five points above their own shuffled baseline. Pairs in which the predictors determine the target by arithmetic are excluded by name (translation->transcript).


**Figure 7. The whole network against its best single layer**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

Compartment recovered from a gene's neighbours, layer by layer, over 3,855 labelled genes against a commonest-answer baseline of 32.7%. The best single layer is comention at 71.3% on 157 genes; every layer at once is 49.2% on 3,837 genes. Weighting each layer by what it is worth, fitted on half the labelled genes and scored on the other half, gives 54.8% against 49.1% unweighted. The `compartment` layer is excluded because it joins genes that share the label being recovered.


**Figure 8. What differs between strains is what goes into the host**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

Nonsynonymous over synonymous SNPs across the sequenced strains, by class of gene, over 7,939 genes carrying a ratio. Dense granules reach 2.06 and apical secretory 1.73, against a length-matched baseline near 1.36; nucleus 0.89 and cytosol 0.72 fall below it. The null is matched gene by gene on coding length (2,000 draws), because a longer gene carries more SNPs of both kinds. The ratio is an index and not dN/dS: the counts are not normalised by the number of sites at which each change could occur.


**Figure 9. Which layers transfer between Toxoplasma and Plasmodium**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

2,221 one-to-one orthologue pairs; 305 orthogroups are set aside as many-to-many rather than averaged. Where the protein is transfers at 0.630 against a 0.374 baseline; whether losing it costs fitness at 0.314. Protein length transfers at 0.907, which is the control on the pairing: a sequence property must transfer between true orthologues, so a layer that does not is reporting biology rather than a mapping fault.


**Figure 10. The host response is mostly private to one effector**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

**(a)** 7,625 host genes moved by at least one of 21 effectors: 3,435 are moved by exactly one and 68 by 10 or more. The translocon is the shared part, overlapping itself at 0.118 against 0.018 for cargo with cargo, which is the control on the import. **(b)** In the single-cell data, 5,617 of 22,553 cells carry one guide; of 85 targets with enough cells, 1 changes the parasite total transcriptional output. The assumption every pseudobulk host analysis rests on -- that the parasite is otherwise unchanged -- therefore largely holds.


**Figure 11. What the map can fill, and what it cannot**

*[panel not yet drawn; the result is measured and the numbers below are from it]*

1,047,574 held-out predictions over 12 targets and 15 strategies, under frozen partitions. 5 of 12 targets are called better than by guessing the commonest answer; the best is lopit_pf_location at 76.8% against 23.9%, answering 3.0% of the time. Accuracy, answer rate and baseline are reported together because any one of them alone can be made to look good.


**Figure 12. What nothing measures, and which experiment would close the most**

![](results/corpus_inventory_2026_11_18/panel_coverage.pdf){width=6.5in}

Across 8,140 genes and 13 categories of knowledge, built from 434 measurement columns. The largest gap is host_response: 7,707 genes (94.7%) have never been measured that way, so one assay of that kind would give them a first measurement. A column with no missing value at all is not a measurement of every gene but a count or a flag in which absence was written as zero, and such columns are excluded.

