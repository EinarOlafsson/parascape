# 72 · The information map: gene × trait values, a matching probability matrix, and calibrated gap filling

**Status: open (requested 2026-11-04).** Builds on 62 (track record), 63 (claims with measured certainty),
70 (scoreboard, deployment rule) and 71 (evidence layers, the three target traits).

## The request (Einar, 2026-11-04)
* "Create an information map where each gene is a row and each column is a trait. There should be an
  accompanying matrix with the same dimensions that holds the probability score for each data point. Then the
  goal would be to fill in the gaps in some systematic way with very probable entries."
* Priority traits: **host-facing**, **binding partners** (especially of host-facing proteins), **function**.
* **Fill bar: 90 % precision** (Einar 2026-11-04, changed from 95 %), measured on held-out direct-evidence genes. A
  predicted entry is written as *filled* only at or above the calibrated probability where held-out precision reaches 0.90;
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
   direct-evidence genes is ≥ 0.90 with its lower 90 % confidence bound reported; a strategy with too few direct
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

## Adding a column (Einar 2026-11-04: "make sure it is easy to add new columns")
* `parascape/information_map.py`: one `register(Trait(name, kind, description))` line, then evidence rows with
  that `trait` name from any source; `build()` places the column in all three tables. Kinds: flag, category,
  terms, measure. `tests/test_information_map.py` checks that adding a column leaves every other column unchanged.

## Trait tiers (Einar, 2026-11-04)
* **Tier 1 = conclusions** that depend on many factors: host_facing, binding partners, function, compartment.
  These are what the map predicts.
* **Tier 2 = measurements**: per-condition relative expression (RH, PRU, ...), protein abundance, fitness,
  screen results, sequence and structure features. These are the evidence predictions are learned from.
* Every trait carries `tier`, `family` (held out together), `normalization`, `role` and `depends_on`
  (`parascape/information_map.py`). Learning is run per trait order: tier 2 only; tier 2 + other observed
  tier-1 traits; and each tier-2 family alone, so it is visible which data carry the signal.
* **Normalization preference (Einar):** `relative_to_control` > `internal_reference` (housekeeping or
  within-sample ratio) > `library_size` > `raw`. Where both a relative and an absolute form of one quantity
  exist, learners take the relative one; absolute values are used only where no control exists.
* Leakage: predicting a trait excludes every family in its `depends_on` (e.g. PVM proximity labelling is
  excluded when predicting host_facing).
* Each expression sample carries a quality score (`results/rna_sample_quality_*`), low-quality samples are
  down-weighted when a column is built, and a column's probability inherits that quality.

## Sequence and structure traits (Einar, 2026-11-04)
Tier 2, one row per gene, from DeepTMHMM topology (`datasets/deeptmhmm_toxo.csv`, all 8,141 genes, per-domain
boundaries) and the existing Foldseek 3Di databases (`~/foldseek/naive_db`, `cleaned_db`; ~6,700 Toxoplasma
entries): `seq_full`, `seq_mature` (signal peptide and known propeptide removed), `seq_nterm_flank`
(N-terminus to the first TM), `seq_cterm_flank` (last TM to the C-terminus), and the same four in the 3Di
alphabet (`tdi_*`). **A protein with no TM domain has identical N- and C-terminal flank columns, equal to its
mature sequence** -- intended, stated in the notebook and tested. Each flank also records the side DeepTMHMM
assigns (inside or outside), since only an outward-facing region can reach the host.

## Definite host-facing negatives (Einar, 2026-11-04: "you also know which proteins are definitely not host facing")
`host_facing` is three-valued: True (shown beyond the parasite), **False (shown to stay inside)**, empty
(untested). Negatives come from two sources, both approved:
* the 16 curated gold negatives (`datasets/hostfacing_goldsets.json` `gold_negative`, e.g. GRA9 by immunogold
  EM in the dense granules then the IVN), tier `observed-curated`;
* hyperLOPIT markers of niches physically inside the parasite: apicoplast, mitochondrion, nucleus, nucleolus,
  40S/60S ribosome, 19S/20S proteasome, Golgi, ER, cytosol, tubulin cytoskeleton, IMC. Probability 0.80 and
  tier `observed-curated`, because this is an inference from a curated resident marker, not a direct export test.
  **Excluded from the negatives:** dense granules, rhoptries, micronemes, apical, plasma membrane and
  endomembrane vesicles -- secretory or boundary compartments whose proteins can be exported.
* The 1,201 PVM-proximity "likely negatives" stay empty: a labelling screen missing a protein is not proof.

## Two cautions found while building the first columns (2026-11-04)
* **DeepTMHMM misses some signal peptides.** GRA16 (TGME49_208830), a secreted effector, is called GLOBULAR
  with no signal peptide, so its `seq_mature` equals `seq_full`. The column is still informative -- a signal
  peptide is predicted for 55 of the 73 host-facing positives and 30 of the 41 dense granule markers, against
  25 of the 494 negatives -- but `seq_mature` is wrong wherever the prediction misses. SignalP or a curated
  cleavage table would fix it; recorded as a gap.
* **The negative set makes a trivial baseline look strong.** Because most negatives are markers of
  compartments inside the parasite (cytosol, nucleus, ribosome), "has a signal peptide" alone separates the
  classes well. Any host-facing strategy must therefore be scored against that baseline, and the question that
  matters is the harder one: **among secretory-pathway proteins, which are exported beyond the parasite?**
  Report that restricted comparison beside the overall one, and score no-knowledge genes separately (the CAFA
  rule) so easy genes do not flatter the result.

## Inference plan (Einar, 2026-11-04)
Stages, each scored on the 70 scoreboard before the next is built:
0. **Baselines**: majority, prevalence-weighted random, the single best feature (signal peptide), ortholog copy.
   Every strategy is reported as improvement over the best of these.
1. **One gradient-boosted-tree model per tier-1 trait** on the tabular tier-2 columns (missing values used
   natively, no imputation); out-of-fold predictions calibrated with isotonic regression.
2. **Sequence and structure as frozen feature blocks** (ESM embeddings, 3Di k-mer profiles), never a
   fine-tuned protein language model: with 73 positives it would memorise them.
3. **Guilt by association**: label propagation over co-expression, crosslinks, co-fitness and structural
   similarity, with the cofitness rule (high correlation AND top rank AND several independent studies).
4. **Joint model over the whole row**: mask known cells, predict them; replaces per-trait models only where it
   beats them on the same folds.
5. **Pairs** with whole proteins held out; the StarPath finding that the CRISPR-fitness *difference* of a pair
   is highly predictive is a feature to reuse.
6. **Fill** only cells above the measured threshold; record strategy, version and the top contributing evidence.

**Fill bound (Einar):** the **lower 90 % confidence bound** of held-out precision must reach 0.90, not the point
estimate. With few positives this may refuse to fill anything at first; that is the intended behaviour.

**Two risks to test for, not assume away:**
* **Too few positives to prove the bound.** 73 host-facing positives over 5 folds is ~15 per fold. Report the
  bound, and treat "cannot yet fill" as a result.
* **Publication-bias shortcut.** A well-studied gene has more filled columns, and attention columns
  (`n_publications`, `lit_tier`, `attention_depth`) could stand in for "famous effector". Keep them in, then
  repeat without them, and report genes with no prior annotation separately (the CAFA no-knowledge rule).

**Leakage specific to host_facing:** the negatives are hyperLOPIT markers, so every `lopit_*`, `compartment*`
and `ortholopit_*` column is excluded; so are `pvm_proximity_positive` and `n_host_targets`, which define the
positives. The eight conflicting genes are excluded from both training and scoring.

## First host-facing scoring run (2026-11-04): the bar is not honestly cleared
Results in `results/host_facing_inference_2026_11_04/` (73 positives, 486 negatives after dropping the 8
conflicts; 409 features; 5 orthogroup-grouped folds; calibration ECE 0.007).
* **Overall the model looks excellent and the appearance is misleading.** ROC-AUC 0.995, precision 0.93,
  recall 0.90, beating every baseline (best baseline: signal peptide alone, precision 0.75). The lowest
  threshold whose **lower 90 % bound** reaches 0.90 is 0.90 itself: 55 held-out calls, 53 correct, bound
  0.9061, which would fill 93 empty genes. But the bound ranges 0.877-0.906 across five reasonable feature
  sets and falls **below** the bar in two of them. With 73 positives this is at the edge of what the labels
  can demonstrate, not a cleared bar.
* **On the question that matters it fails.** Restricted to proteins with a predicted signal peptide (55
  positives, 19 negatives), the best lower bound is **0.868**: among secretory-pathway proteins we cannot yet
  say at 90 % which are exported.
* **The model separates the wrong thing.** Against the 477 hyperLOPIT internal markers its AUC is 0.997;
  against the 9 **curated gold negatives** (proteins shown to stay in the vacuole or the IVN) it collapses to
  **0.855 with a 33 % false-positive rate**. It distinguishes "inside the parasite" from "secreted", not
  "stays in the vacuole" from "exported into the host" -- which is the biological question.
* **Two features do much of the work, and neither is independent evidence.** The top feature is the
  host-transcription effector screen (`hosttx_padj`, `hosttx_T2`; PMID 37827122) -- a screen for genes whose
  loss changes host transcription, which is nearly the conclusion itself. The second is a publication count
  (`n_papers_incidental`). Attention columns alone reach AUC 0.902. Therefore `hosttx_effectors`,
  `organelle_surface_turboid` and `pubmed` are now declared in `host_facing.depends_on` and held out.
* **Decision: fill nothing for host_facing yet.** The next run scores it from genuinely independent evidence,
  and the limiting resource is not the model but **negatives**: 9 curated vacuole-resident proteins cannot
  support a 90 % claim. Expanding them, and Einar's review of the 140 literature candidates, come first.
