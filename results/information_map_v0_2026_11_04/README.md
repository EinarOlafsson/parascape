# Information map v0 (2026-11-04)

Instruction 72, v0 build. Observed facts only; no predictions, no gap filling. Built by
`scripts/build_information_map_v0_notebook.py`, executed notebook `notebooks/information_map_v0_2026_11_04.ipynb`.

## Files
* `values.parquet`, `probability.parquet`, `source.parquet` -- 8140 ME49 genes (index `gene_id`, the rows of
  `parascape/data/nodes.parquet`) x 20 trait columns, identical shape. Empty (NaN) = no observed fact. Never zero, never "no".
  `host_facing` is only ever `True`; screen columns hold both True and False because a screened gene with no phenotype is a measurement.
* `evidence.parquet` -- long form, one row per fact: gene_id, trait, value, probability, source_tier, reference, dataset (3726 rows).
* `pairs.parquet` -- binding partners: gene_a (parasite), gene_b (parasite gene, or host protein name when partner_side = host), partner_side, probability, source_tier, reference (36507 rows).
* `coverage.tsv` -- genes with an observed fact per trait and tier.
* `sources/` -- raw QuickGO pages with retrieval time, GO term names, unmapped UniProt accessions, hyperLOPIT markers outside the gene universe.

## Provisional probabilities (v0)
Tier defaults, NOT measured: observed-direct 0.95, observed-curated 0.90, observed-screen 0.80, observed-proximity 0.70.
StarPath crosslinks: 1 - local_fdr. Interaction-study members: 0.5, provisional until a background table exists.
v1 replaces these by measured agreement with independent sources.

## Rules followed
* TAGM / `lopit_*` node columns are classifier predictions and are excluded; only the 718 hyperLOPIT markers (curated) are used.
* PVM proximity "likely negatives" are NOT written as `host_facing = no` (not near the PVM is not proof of not exported).
* GO annotations with a NOT qualifier are skipped. Experimental evidence codes: ECO:0000269, 314, 353, 315, 316, 270 (with descendants).
* When several sources speak about one gene, the matrix cell takes the highest-probability one; all remain in `evidence.parquet`.
  Hyperlopit and GO CC are separate columns so disagreement stays visible.
* UniProt to ME49 mapping via UniProt VEuPathDB cross-references, restricted to the 8,140-gene universe; ambiguous (>1 gene): 1; unmapped accessions: 11 (listed in `sources/unmapped_uniprot.tsv`).
* The literature host-facing candidate list is not included; it enters after Einar reviews it.

## Coverage (genes with an observed fact)
```
source_tier                   any tier  observed-curated  observed-direct  observed-proximity  observed-screen
trait                                                                                                         
compartment_go_cc                   28                 0               28                   0                0
compartment_hyperlopit             716               716                0                   0                0
function_go_bp                      24                 0               24                   0                0
function_go_mf                      27                 0               27                   0                0
host_facing                        585               548               50                  73                0
partners                          7828                14             1630                   0             7822
screen_actin_phenotype             318                 0                0                   0              318
screen_any_phenotype               318                 0                0                   0              318
screen_apicoplast_phenotype        318                 0                0                   0              318
screen_egress_phenotype            318                 0                0                   0              318
screen_replication_phenotype       318                 0                0                   0              318
screen_scorers_agree                74                 0                0                   0               74
```
