# Which Toxoplasma genes are specifically important under mouse immune pressure?

Built by `scripts/build_mouse_immune_dependence_notebook.py`; the executed notebook is
`notebooks/mouse_immune_dependence_2026_11_04.ipynb`. Asked by Einar 2026-11-04 as a first, concrete
inference target (easier than host-facing) **with a built-in score against known biology**.

## The question
A mouse macrophage activated by IFN-gamma attacks the parasitophorous vacuole with IRG GTPases. Mice have a
large IRG family; humans have essentially none. So: which Toxoplasma genes does the parasite need *only* when
the mouse macrophage is activated?

## The answer in one line
`ifng_dependence_mouse_macrophage = fit_ifng - fit_naive_bmdm` (Wang 2020, PMID 33067458), ranked over
7402 genes measured in both conditions. Rank 1 = most IFN-gamma-dependent. Negative = the mutant loses
ground specifically when the macrophage is activated.

## Bottom line before the detail
The contrast is computable, is free of the scale artefact, and is **not strong enough to answer the question**.
GRA12 and ROP18 land near the top, ROP17 and GRA7 do not, the rank-sum p over the four independent known genes
is about 0.06, the top 30 contains almost no secreted effectors, and the one compartment that survives
multiple testing is micronemes rather than dense granules or rhoptries. The top of the difference is extreme in
both screen arms at once, which is how a difference of two noisy scores behaves, and `fit_ifng` alone scores
the known genes better than the difference does. Treat the ranking as a hypothesis list, not as an answer;
nothing was filled into the knowledge matrix.

## The scale trap, and how it was avoided
The fitness columns in `parascape/data/nodes.parquet` are **not on a common scale and must not be subtracted
across studies**. Re-verified in the notebook with the 158 ribosomal proteins (median, ribosomal vs all genes):

| column | study | ribosomal median | all genes | meaning |
|---|---|---|---|---|
| `fit_invitro_hff` | Sidik 2016, human fibroblasts | -4.25 | -0.81 | absolute growth; negative = important |
| `fit_naive_bmdm` | Wang 2020, naive mouse BMDM | 0.49 | 0.35 | growth already normalised out |
| `fit_ifng` | Wang 2020, IFN-gamma BMDM | -0.04 | -0.20 | growth already normalised out |
| `fit_invivo_*` | Giuliano 2024, six mouse tissues | — | — | differential fitness x significance: a **ranking statistic**, not an effect size (GRA12 reads -41.56) |

Therefore the only subtraction performed anywhere in this analysis is **within Wang 2020** (same cells, same
library, same control-relative baseline). `fit_invitro_hff` is reported beside the result, never subtracted,
and `fit_invivo_*` is used only through Spearman rank correlation.

## Validation (a check, never a training set)
Positive set: the IRG-pathway effectors, mouse-specific because humans lack IRGs, plus GRA45 from the Wang
abstract itself (fetched from PubMed; `wang2020_pubmed_abstract.txt`). No threshold or parameter was chosen
using them.

| gene | ME49 id | ifng_dependence | rank | percentile | why expected |
|---|---|---|---|---|---|
| ROP18 | TGME49_205250 | -1.900 | 447 of 7402 | top 6.04 % | IRG pathway: kinase phosphorylating IRGs |
| ROP17 | TGME49_258580 | -0.750 | 2635 of 7402 | top 35.60 % | IRG pathway: second kinase of the ROP5 complex |
| GRA7 | TGME49_203310 | -0.379 | 4823 of 7402 | top 65.16 % | IRG pathway: adaptor anchoring the complex at the PVM |
| GRA12 | TGME49_288650 | -2.533 | 161 of 7402 | top 2.18 % | IFN-gamma resistance, separate from the ROP5 complex |
| ROP5 | TGME49_308090 | not measured | — | — | IRG pathway scaffold; tandem array, untargetable by CRISPR |
| GRA45 | TGME49_316250 | -4.299 | 21 of 7402 | top 0.28 % | named and validated in the Wang 2020 abstract (not independent) |

* ROP5 is **not measured** in the Wang screen, as expected: it is a tandem array of near-identical copies that
  a CRISPR library cannot target cleanly. That is a property of the technology, not a negative result.
* GRA45 comes from the same screen, so it is a positive control for the screen, not independent literature.

Test, independent literature genes only (n = 4 measured, background 7398):
Mann-Whitney one-sided **p = 0.057**, ROC-AUC **0.728**, median percentile
20.82 %; 1 in the top 5 % (expected
0.2, hypergeometric p = 0.19),
0 in the top 1 % (expected 0.04,
p = 1).
With GRA45 added: p = 0.015, AUC 0.782.

**Read honestly: this is a weak, not a passing, result.** GRA12 (top 2 %) and ROP18 (top 6 %) land where the
biology says they should; ROP17 (36th percentile) and GRA7 (65th percentile) do not, and GRA7 sits in the
middle of the ranking, i.e. the contrast says nothing about it. The one-sided rank-sum p of about 0.06 is a
shift in the right direction that does **not** reach conventional significance. Adding GRA45 takes p below
0.02, but GRA45 was named by this very screen, so that is internal consistency, not recovery. With four
independent positives this cannot measure precision, and it is not a held-out benchmark in the instruction 70
sense -- it is a sanity score.

A likely reason GRA7 is flat is biological rather than technical: GRA7 is also needed for normal vacuole
architecture in unactivated cells, so it costs fitness in *both* Wang conditions and the subtraction cancels
it. A contrast can only see a gene whose requirement is conditional.

## Is the subtraction helping? (`score_comparison.tsv`)
The top of the difference is extreme in **both** arms at once, in opposite directions (top 30: median
`fit_naive_bmdm` +1.99 and median `fit_ifng` -2.64, against
+0.35 and -0.20 for all screened genes; all 30 have
`fit_naive_bmdm` > 0 and `fit_ifng` < 0). That is the signature of a difference of two noisy scores, not of a
gene simply needed under IFN-gamma. Scoring three alternatives on the four known genes:

| score | rank-sum p | ROC-AUC | Spearman vs brain in vivo |
|---|---|---|---|
| `fit_ifng - fit_naive_bmdm` (the contrast) | 0.057 | 0.728 | 0.106 |
| `fit_ifng` alone | 0.012 | 0.827 | 0.053 |
| `fit_naive_bmdm` alone (opposite direction expected) | 0.670 | 0.436 | -0.118 |

Because both Wang columns are already control-normalised, `fit_ifng` alone is itself a relative quantity and
may not need the naive arm subtracted; on the known genes it scores better, and the naive arm alone carries no
signal in its own expected direction. **But the evidence is mixed**: the contrast agrees better with the in
vivo tissue rankings and gives the stronger compartment enrichment, and with four positives the difference
between p = 0.012 and p = 0.057 is not real. **Neither score is demonstrated to be better, and neither should
be deployed.** Deciding needs a genuine positive set -- the seven genes Wang validated individually, or an
independent IFN-gamma-resistance list -- scored the instruction 70 way with positives neither score has seen
(comparing scores on this small set is already a use of it for selection).

## The top of the ranking
| rank | gene | product | ifng_dependence | hyperLOPIT | fit_invitro_hff |
|---|---|---|---|---|---|
| 1 | TGME49_269620 | hypothetical protein | -5.514 | — | -0.34 |
| 2 | TGME49_313600 | DDHD domain-containing protein | -5.450 | — | -1.04 |
| 3 | TGME49_237210 | tyrosine kinase-like protein TKL4 | -5.390 | — | -1.26 |
| 4 | TGME49_251440 | troponin c, isotype gamma, putative | -5.248 | — | -0.57 |
| 5 | TGME49_212770 | inner membrane complex protein IMC34 | -5.100 | — | -1.80 |
| 6 | TGME49_232800 | hypothetical protein | -5.032 | — | -1.98 |
| 7 | TGME49_321530 | cathepsin CPL | -5.002 | endomembrane vesicles | 0.68 |
| 8 | TGME49_295340 | UV excision repair protein Rad23 protein | -4.985 | — | -0.37 |
| 9 | TGME49_312420 | dense granule protein GRA38 | -4.861 | dense granules | -1.15 |
| 10 | TGME49_305010 | pre-mRNA branch site protein p14, putative | -4.849 | nucleus | 0.41 |

**This top has no effector signature**, and that is the main negative result of the analysis: hypothetical
proteins, a splicing factor, a DNA-repair protein, an IMC protein, a troponin. Only GRA38 and GRA45 (rank 21)
are dense granule proteins. The expected answer -- a list dominated by secreted vacuole-facing effectors --
did not appear. See the next section for why.

Full table: `ranking.tsv`; the top 30 with signal peptide, host_facing and publication count: `top30.tsv`.
Artefact check: 0 ribosomal proteins in the top 5 %
(expected 7.0), and the ribosomal median
ifng_dependence is -0.531
against -0.563 for all screened genes -- no sign of a growth/scaling
artefact driving the top.

## Compartment enrichment (Fisher, one-sided, top 5 % vs the rest of the screened hyperLOPIT markers)
705 screened genes carry a hyperLOPIT compartment; 34 of them are in the top 5 %.

| compartment | markers screened | in top 5 % | expected | odds ratio | Fisher p | BH q |
|---|---|---|---|---|---|---|
| micronemes | 28 | 6 | 1.35 | 6.32 | 0.0014 | 0.024 |
| tubulin cytoskeleton | 8 | 2 | 0.39 | 6.93 | 0.053 | 0.39 |
| nucleus | 94 | 8 | 4.53 | 2.09 | 0.07 | 0.39 |
| endomembrane vesicles | 13 | 2 | 0.63 | 3.75 | 0.13 | 0.54 |
| dense granules | 37 | 3 | 1.78 | 1.81 | 0.26 | 0.78 |
| IMC | 38 | 3 | 1.83 | 1.76 | 0.28 | 0.78 |
| Golgi | 20 | 1 | 0.96 | 1.04 | 0.63 | 1 |
| rhoptries | 45 | 2 | 2.17 | 0.91 | 0.65 | 1 |
| plasma membrane | 64 | 2 | 3.09 | 0.61 | 0.83 | 1 |
| mitochondrion | 68 | 2 | 3.28 | 0.57 | 0.86 | 1 |
| ER | 39 | 1 | 1.88 | 0.50 | 0.86 | 1 |
| apical | 41 | 1 | 1.98 | 0.48 | 0.88 | 1 |
| apicoplast | 46 | 1 | 2.22 | 0.42 | 0.9 | 1 |
| cytosol | 62 | 0 | 2.99 | 0.00 | 1 | 1 |
| proteasome | 28 | 0 | 1.35 | 0.00 | 1 | 1 |
| nucleolus | 15 | 0 | 0.72 | 0.00 | 1 | 1 |
| ribosome | 59 | 0 | 2.85 | 0.00 | 1 | 1 |

**The expected enrichment is not the one that appears.** The prediction was dense granules and rhoptries --
the secreted effectors. What comes out instead is **micronemes** (the only compartment surviving multiple
testing), with dense granules and rhoptries near their expected counts (odds ratios around 1.8 and 0.9).
Micronemes are the invasion organelle, so an IFN-gamma-specific requirement for them is not an obvious
artefact -- a parasite that invades or egresses poorly spends longer exposed to an activated macrophage -- but
it is a hypothesis generated after the fact from 6 genes out of 28 markers, not a finding. The reassuring half
is negative: cytosol, ribosome, nucleolus and proteasome contribute **zero** genes to the top 5 %, and
apicoplast and mitochondrion are below expectation, so the top is not a growth or scaling artefact.

Caveat: hyperLOPIT has **no PVM class** -- the vacuole membrane is not a fraction in that experiment, so
"dense granules" and "rhoptries" are the only available proxies for vacuole-facing, and several compartments
have very few screened markers, where one gene moves the p-value. Only 705 of the 7402 screened genes carry a
compartment at all.

## In vivo agreement (ranking only)
| tissue | genes in common | Spearman rho | p |
|---|---|---|---|
| heart | 7383 | 0.109 | 5e-21 |
| brain | 7383 | 0.106 | 5.9e-20 |
| lung | 7383 | 0.069 | 3e-09 |
| spleen | 7383 | 0.066 | 1.7e-08 |
| liver | 7383 | 0.063 | 5e-08 |
| PE | 7383 | 0.038 | 0.0013 |
| young2019 | 115 | -0.118 | 0.21 |

Agreement is **real but very weak**: the strongest tissue (heart, brain) reaches Spearman rho of only about
0.11, which is highly significant only because 7,383 genes are being correlated. It means the contrast is not
pure noise; it does not mean the two screens measure the same thing. `fit_invivo_young2019` covers only 115
genes and is not interpretable here. Both screens are in mice and use overlapping library designs, so
agreement means "not specific to one screen", not independent confirmation.

## The human comparison, stated honestly
Of the top 30 IFN-gamma-dependent genes, 20 have
`fit_invitro_hff > -1`, i.e. they are dispensable for growth in human fibroblasts (median
-0.34 against -0.81 for
all screened genes). The supported claim is **"needed under mouse immune pressure but not needed for growth in
human fibroblasts"**, which is *not* "more important in mice than in humans": the two screens differ in host
species, in immune activation and in scale all at once, and a resting fibroblast is not an activated immune
cell. The missing measurement is a genome-wide Toxoplasma CRISPR screen in **IFN-gamma-activated human
macrophages against matched unactivated human cells** -- the human twin of the Wang design -- so that the two
immune-dependence contrasts are each internally normalised and comparable.

## Every other caveat
* The contrast is a difference of two noisy screen scores, so its variance exceeds either column and its
  extremes are biased towards genes that were noisy in both arms; read single genes as a rank, not as an effect
  size, and do not quote a gene from the top 30 as "IFN-gamma-dependent" without an individual experiment.
* Wang reported 353 significant genes with their own statistics; `nodes.parquet` carries only the two fitness
  columns, so no per-gene significance or replicate variance is available here. Having the published FDR or
  the per-guide counts would let the contrast be tested properly instead of ranked raw, and is the first thing
  worth acquiring.
* Genes absent from the Wang library (poor guide coverage, tandem arrays like ROP5) are missing, not negative.
* IFN-gamma dependence measured in bone-marrow macrophages is one cell type of one mouse strain; it does not
  cover IRG-independent restriction, other cytokines, or other host cells.
* hyperLOPIT compartments cover only a few hundred of the screened genes, so compartment enrichment speaks
  about marker proteins, not about the whole top of the list.
* `host_facing` in the information map is a curated flag over 73 positives and 494 negatives; most of the
  ranking is simply empty there, and empty is not "no".
* Nothing here was filled into the information map. The derived column is registered as tier-2 evidence
  (`family='crispr_macrophage'`, `normalization='relative_to_control'`, `role='measurement'`) and the
  registration line is printed in the notebook for review; `parascape/information_map.py` was not edited.

## Files
| file | contents |
|---|---|
| `ranking.tsv` | all 7402 screened genes: product, ifng_dependence, rank, percentile, both Wang columns, fit_invitro_hff, compartment, host_facing, n_publications |
| `known_genes.tsv` | the validation set with ranks; the test result is appended as comment lines |
| `top30.tsv` | the 30 most IFN-gamma-dependent genes with signal peptide and host_facing |
| `compartment_enrichment.tsv` | Fisher test per coarse hyperLOPIT compartment |
| `score_comparison.tsv` | the contrast against `fit_ifng` alone and `fit_naive_bmdm` alone on the known genes |
| `invivo_agreement.tsv` | Spearman of the ranking against each mouse tissue column |
| `known_genes_rank_curve.png` | the positive set drawn on the full ranking |
| `wang2020_pubmed_abstract.txt` | the fetched abstract the GRA45 entry came from |
