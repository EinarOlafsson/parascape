# Held-out member recovery for seed expansion (instruction 64.12)

Builder `scripts/build_seed_expansion_benchmark_notebook.py` (one niche per process); notebook
`notebooks/seed_expansion_benchmark_2026_11_07.ipynb`. Each niche's curation-grade hyperLOPIT markers split by
orthogroup, half seeds and half hidden (seed 0); seed expansion run from the seeds with the hyperLOPIT
compartment excluded (checked: no candidate links to the seeds through it); hidden members counted in the
ranked genome at depths fixed before the run. Universe = every gene but the seeds.

| niche | seeds | hidden | in top 50 | in top 100 | recall at 100 | lift at 100 |
|---|---|---|---|---|---|---|
| rhoptries | 21 | 26 | 14 | 21 | 0.808 | 65.6 |
| micronemes | 15 | 14 | 11 | 12 | 0.857 | 69.6 |
| dense granules | 19 | 21 | 9 | 12 | 0.571 | 46.4 |
| IMC | 19 | 19 | 9 | 10 | 0.526 | 42.7 |
| apicoplast | 24 | 24 | 4 | 7 | 0.292 | 23.7 |
| mitochondrion | 33 | 35 | 7 | 12 | 0.343 | 27.8 |

The secretory and division organelles are found strongly; apicoplast and mitochondrion more weakly. Genes that
are not markers are not verified negatives, so precision is a lower bound. One split per niche; a repeat over
several seeds would give intervals.

## Second adapter: positive-unlabelled learning, same splits

Recall of hidden members at depth 100 (`recovery.tsv` has every depth):

| niche | positive-unlabelled | seed expansion |
|---|---|---|
| rhoptries | 0.769 | 0.808 |
| micronemes | 0.786 | 0.857 |
| dense granules | **0.762** | 0.571 |
| IMC | 0.474 | 0.526 |
| apicoplast | **0.458** | 0.292 |
| mitochondrion | 0.371 | 0.343 |

Neither wins everywhere: the network walk is better for the invasion organelles, the classifier for dense
granules and the apicoplast. Each finished niche is cached as a row file and reused by the notebook.

## Third adapter: set enrichment (profile resemblance), same splits

Recall at 100: rhoptries 0.692, micronemes 0.643, dense granules 0.524, IMC 0.526, apicoplast 0.333, mitochondrion
0.229 -- never the best of the three. Its profile draws on every feature category the seeds differ in, including
publication counts (study effort), which the other two do not lean on; read it with that in mind.

## Fourth: gene-set search, which is set retrieval, not a ranking

`geneset_hunt` returns ONE cluster (about 16 candidates), the best of 16 map and clustering settings, so depth-100
recall is not its measure. Hidden members inside its candidates: dense granules 3, mitochondrion 3, rhoptries 2,
IMC 2, micronemes 0, apicoplast 0 -- roughly 0.1-0.2 precision against a base rate of 0.002-0.004 where it finds
any, and nothing for two niches.
