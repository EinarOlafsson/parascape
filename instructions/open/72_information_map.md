# 72 · The information map: gene × trait values, a matching probability matrix, and calibrated gap filling

**Status: open (requested 2026-11-04).** Builds on 62 (track record), 63 (claims with measured certainty),
70 (scoreboard, deployment rule) and 71 (evidence layers, the three target traits).

## The request (Einar, 2026-11-04)
* "Create an information map where each gene is a row and each column is a trait. There should be an
  accompanying matrix with the same dimensions that holds the probability score for each data point. Then the
  goal would be to fill in the gaps in some systematic way with very probable entries."
* Priority traits: **host-facing**, **binding partners** (especially of host-facing proteins), **function**.
* **Fill bar: 95 % precision (α = 0.05)**, measured on held-out direct-evidence genes. A predicted entry is
  written as *filled* only at or above the calibrated probability where held-out precision reaches 0.95;
  everything below stays visible as *exploratory*, never as a fact.

## Shape
Three aligned tables over the same genes (ME49 IDs, the 8,140 rows of `nodes.parquet`) and trait columns:

| Table | Cell content |
|---|---|
| `values` | the fact (yes/no, a category, a set of terms, a number) |
| `probability` | calibrated probability that the fact is true (0–1) |
| `source` | `observed-direct` · `observed-curated` · `observed-screen` · `observed-proximity` · `filled:<strategy>@<version>` · `exploratory:<strategy>` · empty; plus the reference (PMID, dataset, run) in a long-form evidence table keyed by gene × trait |

* Pair-valued traits (binding partners) live in a companion **pairs** table (gene_a, gene_b, side
  parasite/host, value, probability, source, evidence); the main matrix holds summaries per gene: best observed
  partner, its probability, partners above 0.95, partners above 0.5.
* Multi-valued traits (functions, GO terms) are stored long-form (gene, term, probability, source); the matrix
  holds the top term and the count of terms at or above the fill bar.
* Missing is empty, never zero or "no". Unlabeled is not negative.

## Probabilities for observed facts
* v0 (provisional, recorded as such): tier defaults -- direct experiment 0.95, curated marker 0.90, single
  high-throughput screen call 0.80, proximity labelling 0.70 -- except where the evidence carries its own
  error rate (StarPath crosslinks: 1 − local FDR; Sage identifications: 1 − q-value).
* v1: replace tier defaults by measured agreement: for each source, how often its calls are confirmed by an
  independent direct source on the genes both cover (with a confidence interval). A source never gets a
  probability above what its replication supports.

## Filling the gaps
1. Per trait, every strategy (evidence layer) predicts held-out known cells under the 70 rules: orthogroup-grouped
   folds, whole-source hold-out, cold-node splits for pairs, baselines in the same table.
2. Out-of-fold scores are calibrated to probabilities (isotonic, `claims.CertaintyModel`), with the reliability
   curve and expected calibration error reported.
3. The fill threshold per strategy × trait is the lowest calibrated probability whose held-out precision on
   direct-evidence genes is ≥ 0.95 with its lower 90 % confidence bound reported; a strategy with too few direct
   genes to show that cannot fill.
4. Filled cells are features for other traits only with their probability attached, and **never become training
   labels**.
5. A joint multi-trait model (ToxoNet v2 design, masked-cell training) is compared on the same scoreboard and
   replaces per-trait strategies only where it beats them.
6. Conflicts between sources are kept side by side and flagged, never averaged away.
7. Every release of the filled map is frozen with a date for a prospective check against later publications.

## v0 build (first step)
`notebooks/information_map_v0_2026_11_04.ipynb` assembles observed facts only:
* host-facing: the 18 curated host interactions, `pvm_proximity_positive`, experimental GO host-cell terms;
  the literature list (`results/host_facing_truth_candidates_2026_11_04/`) enters after Einar reviews it.
* compartment: the 718 hyperLOPIT markers (curated) and experimental GO cellular component (direct); TAGM
  calls are predictions and are not observed facts.
* partners: StarPath crosslinks (1 − local FDR), parsed interaction studies (provisional until the CRAPome
  background exists), curated host interactions.
* function: experimental GO biological process / molecular function; screen phenotypes in `nodes.parquet`.
* Output: `values`, `probability`, `source` parquet files plus the long evidence table, and a coverage table
  (how many genes have an observed fact per trait).
