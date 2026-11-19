% ParaScape: one gene space for the apicomplexan literature% 18 November 2026
## Results

### A corpus that is nearly complete in transcription and nearly empty elsewhere

ParaScape holds 177 datasets in one gene space of 8,140 genes
(Figure 1a). Read as coverage rather than as effort, the corpus is lopsided: transcript reaches
99.9% of genes and chromatin
98.5%, while function reaches
16.1%, phenotype 6.9% and
host response 5.3% (Figure 1b). What is missing from the
literature is not what a gene is transcribed at, but what it does and what it does to the host.

Coverage here counts only measurements. 51 columns are excluded because
they carry no missing value at all: a count that is zero where nobody looked, or a flag that is
False for an unassayed gene, would otherwise report a category as measured for every gene in the
genome. The distinction recurs throughout and is not pedantic -- it changed four categories from
fully covered to half covered.

### Studies cannot be compared until their untreated arms agree

Two perturbations measured in different laboratories can only be set against each other if the two
untreated arms describe the same starting parasite. Over 12 studies with a control
group that is never also a treated arm, agreement between those arms rises from
0.832 to 0.887 when the cell-cycle axis
is removed. A random axis of the same size and shape gains +0.018
(sd 0.005, 200 draws), so the cell-cycle excess is
+0.038 at p = 0.0050 (Figure 2). Cell cycle is a
substantial part of what separates one laboratory's control from another's.

That this matters at all is shown directly by the corpus. Across 20,100 pairs of stored
within-study contrasts, two from the same study agree at a median absolute Spearman of
0.284 and two from different studies at
0.085 -- a factor of 3.3 (Figure 4). Sharing the
KIND of perturbation across studies buys almost nothing. The laboratory, not the perturbation, is
the strongest thing structuring the published record, which makes comparability a prerequisite
rather than a precaution. Of 178 parsed contrasts over 19 studies,
1 perturbation is stated by two independent
laboratories, so whether perturbation effects replicate cannot be asked of this literature at all.

### Genes that move together are compartments, and transcription places the ones built in waves

Cut on co-expression components from 57 independent studies, 8,018 genes form
24 modules, of which 14 hold a compartment above a baseline matched for
size and for how much is already known about their members (Figure 3). One module is 78.6 per cent
ribosome against 3.8 per cent matched. The apical secretory machinery splits into two modules
rather than one.

The same coupling is visible from the other direction. Of the 13 coarse hyperLOPIT compartments,
five are recovered as well or better from RNA co-expression alone as from all allowed evidence
together -- ribosome, rhoptries, apical, micronemes and the IMC -- while the plasma membrane (0.67)
and the ER (0.78) are not (Figure 6a). Every one of the five is built in a wave during division, so
its genes share a schedule that transcription records; the two that resist are made continuously.

### What sets protein abundance, and what the answer depends on

In Plasmodium falciparum, transcript abundance explains
48.6% of the variance in protein abundance on the
439 genes carrying every layer, and adding turnover, ribosome occupancy and
codon usage takes it only to 50.5% (Figure 6b). The
qualification belongs beside the claim: on the wider 695-gene set transcript
explains 24.3%, so protein is roughly twice as predictable in
the genes everybody could measure. Compared with the human literature in its own currency, one
averaged mRNA measure against protein gives r2 = 0.396
over 2,611 genes, at the uncorrected end of published estimates.

Across every pair of knowledge categories, 11 of 14
recover more than five points of variance above their own shuffled baseline (Figure 6c). Protein
abundance explains fitness at R2 = 0.373: what the parasite makes most of is what it can least
afford to lose.

### Evidence of many kinds, weighted by what each is worth

Recovering a hidden compartment from a gene's neighbours, the best single edge layer is co-mention
at 71.3% on 157 genes, while every layer at once reaches
3,837 genes at 49.2% (Figure 7) -- twenty-four times the reach at
two thirds the accuracy. The naive union is dragged down by its largest layer: co-fitness carries
88,997 edges and performs below the commonest-answer baseline. Weighting each layer by its own
measured worth, fitted on half the labelled genes and scored on the other half, gives
54.8% against
49.1% unweighted.

### What varies between strains is what the parasite puts into the host

Nonsynonymous over synonymous SNPs across the sequenced strains, against a baseline matched gene by
gene on coding length: dense granules reach
2.06 and apical secretory
1.73, against
1.36, while the nucleus
(0.89), the cytosol (0.72) and the most required
genes (0.89) fall at or below it
(Figure 8). The matched baseline is necessary rather than decorative: longer genes carry more SNPs
of both kinds and secreted proteins are not a random draw from the length distribution.

### Which layers survive the species boundary

On 2,221 one-to-one orthologue pairs, where a protein is transfers
between Toxoplasma and Plasmodium at 0.630 against a 0.374
baseline, whether losing it costs fitness at 0.314, and how its message is written
not at all (Figure 9). Protein length transfers at 0.907, which is the control on
the pairing: a sequence property must transfer between true orthologues, so a layer that does not is
reporting biology rather than a mapping fault. 305
orthogroups are set aside as many-to-many rather than averaged.

### The host response belongs mostly to single effectors

Of 7,625 host genes moved by at least one of
21 parasite effectors, 3,435 are moved by exactly
one and 68 by half or more (Figure 10a). The experiment carries its own
control and it passes: the translocon members overlap each other at
0.118 against 0.018 for cargo
with cargo. In the single-cell data, of 85 targets with enough cells
1 changes the parasite's own transcriptional output, so
the assumption on which every pseudobulk host attribution rests largely holds (Figure 10b).

### A treatment reaches some bradyzoite subtypes and not others

Tested against all 201 stored within-study contrasts with 1,000
expression-matched random gene sets behind each, 135 contrasts reach at least one
of the five in vivo bradyzoite subtypes. Subtype C is reached by
4 of them, against 86 for the most reached
(Figure 5). Subtype C is testable and its matched sets are centred in every contrast checked,
so this is not a failure of the comparison: nothing in the published catalogue of drugs, knockouts,
pH switches and time points moves it. If a compound must reach every subtype to clear a tissue cyst,
subtype C is the gap.

### What the map can fill, and what it cannot

Over 1,047,574 held-out predictions under frozen partitions,
5 of 12 targets are called better
than by guessing the commonest answer (Figure 11). The best is lopit_pf_location at
76.8% against 23.9%, answering
3.0% of the time. Accuracy is never shown without its baseline: `is_exported`
is 95.9 per cent right against a 96.6 per cent baseline, which is worse than guessing.

### What nothing has measured

Across 8,140 genes and 13 categories built from
434 measurement columns, the largest gap is host_response:
7,707 genes (94.7%) have never been measured
that way (Figure 12). One assay of that kind would give them a first measurement, which is the most
any single experiment in this catalogue could add.

## Methods

### Baselines are matched, never prior category size

Every enrichment is divided by what the same procedure returns on random gene sets drawn to match
the real set in size and in expression or in another stated property, not by the prior size of the
category. A small category wins a comparison by shape alone, and dividing by its prior size
measures that shape. Where a matched set cannot be centred on the real one -- a set at the extreme
of the expression range has nothing more extreme to be compared with -- the result is reported as
not separable rather than as an effect.

### Absence is unknown, not zero

A gene with no value in an assay has not been shown to lack the quantity; the assay did not reach
it. Columns with no missing value are therefore excluded from any coverage statement, because a
count that is zero where nobody looked and a flag that is False for an unassayed gene both read as
measurements of every gene. Where a cluster or a class is too thinly measured to test, it is
reported as unknown and never as unaffected.

### A negative carries the smallest effect it could have seen

Every comparison that returns no effect also reports the smallest effect its own baseline could
have resolved, so "nothing there" and "could not have seen it" are distinguishable. Probabilities
are corrected across the comparisons made within each analysis.

### Held-out scoring

Predictions are scored only on facts hidden from the method, under frozen partitions, and reported
with three numbers together: how often the call is right when it answers, how often it answers, and
what the commonest answer would give. Any one of the three alone can be made to look good.

### Reproduction

Each analysis writes its result beside a record of the inputs and settings it was made from, and
each is runnable from within the application itself.

## Figures

**Figure 1. What ParaScape is, and how much of the genome it reaches**

![](results/figure1_schematic_2026_11_18/panel_schematic.pdf){width=6.5in}

**(a)** 177 datasets imported into one gene space of 8,140 genes and 362 slots, made comparable across studies, and read through 24 modules. Counts are read from the dataset catalogue at draw time. **(b)** The share of genes anything has ever measured, per category of knowledge. Transcript reaches 99.9% and chromatin 98.5%, while function reaches 16.1%, phenotype 6.9% and host response 5.3%. 51 columns are set aside as derivations, flags or counts in which absence was written as zero.


**Figure 2. Which studies can be compared at all**

![](results/paper_panels_2026_11_18/figure02_comparability.pdf){width=6.5in}

Agreement between the untreated arms of 12 independent studies, over 8,157 genes. Removing the cell-cycle axis raises the median agreement from 0.832 to 0.887. A random axis of the same size and shape gains +0.018 (sd 0.005, 200 draws), so the cell-cycle excess is +0.038 at p = 0.0050. A study whose control profile carries one value is refused rather than scored.


**Figure 3. The co-expression modules, and what each one is**

![](results/paper_panels_2026_11_18/figure03_modules.pdf){width=6.5in}

8,018 genes cut into 24 modules on co-expression components from 57 independent studies. 14 modules hold a compartment above their own matched baseline; the rest hold none, which is reported rather than forced. The baseline is a random gene set matched to each module in size and in how much is already known about its members (1000 draws), never the prior size of the compartment.


**Figure 4. The laboratory structures the corpus, not the perturbation**

![](results/paper_panels_2026_11_18/figure04_what_structures.pdf){width=6.5in}

20,100 pairs of stored contrasts. Two from the same study agree at a median absolute Spearman of 0.284, two from different studies at 0.085 -- 3.3 times. Of 178 parsed contrasts over 19 studies, 1 perturbation is stated by two independent laboratories, so cross-laboratory replication cannot be asked of this corpus.


**Figure 5. The bradyzoite subtype nothing published reaches**

![](results/bradyzoite_reach_2026_11_18/panel_reach_by_subtype.pdf){width=6.5in}

Five subtypes of in vivo bradyzoite tested against all 201 stored within-study contrasts, 1,000 expression-matched random gene sets behind each. 135 contrasts reach at least one subtype. Subtype C is reached by 4, against 86 for the most reached. Each subtype is marked by 150 genes.


**Figure 6. Which measurement explains which other measurement**

![](results/coexpression_components_localization_2026_11_06/panel_which_evidence_places_a_protein.pdf){width=6.5in}

**(a)** The 13 coarse hyperLOPIT compartments recovered from RNA co-expression alone against all allowed evidence, as cross-validated AUC. Five are placed as well or better by transcription alone, and each is built in a wave during division. **(b)** What sets protein abundance in *Plasmodium*: transcript abundance explains 48.6% of the variance on the 439 genes carrying every layer, and every layer together 50.5%. On the wider 695-gene set transcript explains 24.3%. Comparable with the human literature, one averaged mRNA measure against protein gives r2 = 0.396 over 2,611 genes. **(c)** Every pair of categories: 11 of 14 recover more than five points above their own shuffled baseline. Pairs in which the predictors determine the target by arithmetic are excluded by name (translation->transcript).


**Figure 7. The whole network against its best single layer**

![](results/paper_panels_2026_11_18/figure07_network_layers.pdf){width=6.5in}

Compartment recovered from a gene's neighbours, layer by layer, over 3,855 labelled genes against a commonest-answer baseline of 32.7%. The best single layer is comention at 71.3% on 157 genes; every layer at once is 49.2% on 3,837 genes. Weighting each layer by what it is worth, fitted on half the labelled genes and scored on the other half, gives 54.8% against 49.1% unweighted. The `compartment` layer is excluded because it joins genes that share the label being recovered.


**Figure 8. What differs between strains is what goes into the host**

![](results/paper_panels_2026_11_18/figure08_strain_variation.pdf){width=6.5in}

Nonsynonymous over synonymous SNPs across the sequenced strains, by class of gene, over 7,939 genes carrying a ratio. Dense granules reach 2.06 and apical secretory 1.73, against a length-matched baseline near 1.36; nucleus 0.89 and cytosol 0.72 fall below it. The null is matched gene by gene on coding length (2,000 draws), because a longer gene carries more SNPs of both kinds. The ratio is an index and not dN/dS: the counts are not normalised by the number of sites at which each change could occur.


**Figure 9. Which layers transfer between Toxoplasma and Plasmodium**

![](results/paper_panels_2026_11_18/figure09_species_transfer.pdf){width=6.5in}

2,221 one-to-one orthologue pairs; 305 orthogroups are set aside as many-to-many rather than averaged. Where the protein is transfers at 0.630 against a 0.374 baseline; whether losing it costs fitness at 0.314. Protein length transfers at 0.907, which is the control on the pairing: a sequence property must transfer between true orthologues, so a layer that does not is reporting biology rather than a mapping fault.


**Figure 10. The host response is mostly private to one effector**

![](results/paper_panels_2026_11_18/figure10_host_shared.pdf){width=6.5in}

**(a)** 7,625 host genes moved by at least one of 21 effectors: 3,435 are moved by exactly one and 68 by 10 or more. The translocon is the shared part, overlapping itself at 0.118 against 0.018 for cargo with cargo, which is the control on the import. **(b)** In the single-cell data, 5,617 of 22,553 cells carry one guide; of 85 targets with enough cells, 1 changes the parasite total transcriptional output. The assumption every pseudobulk host analysis rests on -- that the parasite is otherwise unchanged -- therefore largely holds.


**Figure 11. What the map can fill, and what it cannot**

![](results/paper_panels_2026_11_18/figure11_filling_the_map.pdf){width=6.5in}

1,047,574 held-out predictions over 12 targets and 15 strategies, under frozen partitions. 5 of 12 targets are called better than by guessing the commonest answer; the best is lopit_pf_location at 76.8% against 23.9%, answering 3.0% of the time. Accuracy, answer rate and baseline are reported together because any one of them alone can be made to look good.


**Figure 12. What nothing measures, and which experiment would close the most**

![](results/corpus_inventory_2026_11_18/panel_coverage.pdf){width=6.5in}

Across 8,140 genes and 13 categories of knowledge, built from 434 measurement columns. The largest gap is host_response: 7,707 genes (94.7%) have never been measured that way, so one assay of that kind would give them a first measurement. A column with no missing value at all is not a measurement of every gene but a count or a flag in which absence was written as zero, and such columns are excluded.

