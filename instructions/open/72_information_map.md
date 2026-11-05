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

## Mouse immune-dependence attempt (2026-11-04): the contrast does not work
`results/mouse_immune_dependence_2026_11_04/`. Question: which genes matter specifically under mouse immune
pressure, from `fit_ifng` − `fit_naive_bmdm` (Wang et al., PMID 33067458; same cells, same study, so the only
safe subtraction among the fitness columns).
* **Known genes are not recovered significantly.** GRA12 rank 161 of 7,402 (top 2.2 %) and ROP18 rank 447
  (top 6 %), but ROP17 is at 35.6 % and GRA7 at 65.2 %: Mann-Whitney p = 0.057, AUC 0.73 on four positives.
  ROP5 is not measured at all (tandem array, not targetable), recorded as not measured rather than a miss.
  GRA45 ranks 21st but is named in the Wang abstract, so it tests the screen, not our ranking.
* **Subtracting two already-relative scores amplifies noise.** All 30 top-ranked genes are extreme in *both*
  arms in opposite directions (median naive +1.99, IFN-gamma −2.64, against +0.35 and −0.20 overall) -- the
  signature of a difference of two noisy numbers. `fit_ifng` **alone** does better on the known genes
  (p = 0.012, AUC 0.83) than the difference (p = 0.057, AUC 0.73), and `fit_naive_bmdm` alone carries nothing
  (p = 0.33). Einar's rule "relative, normalised to a control, is better" holds; the corollary learned here is
  that **subtracting two relative scores is not better still**.
* **Not a scaling artefact, though:** zero ribosomal proteins are in the top 5 % where 7 are expected, and
  cytosol, ribosome, nucleolus and proteasome contribute no top genes. The one compartment enriched after
  correction is **micronemes** (6 of 28, odds ratio 6.3, q = 0.024) -- not dense granules or rhoptries, and
  hyperLOPIT has no PVM class, so a PVM test is impossible here.
* **In vivo agreement is real but tiny** (brain rho 0.11, heart 0.11; significant only because n is large).
* **Decision: deploy neither score.** The blockers are four known positives and the absence of Wang's
  per-gene significance in `nodes.parquet`; acquiring the published FDR table is the first thing worth doing.
* **The human comparison cannot be made.** 20 of the top 30 are dispensable for growth in human fibroblasts,
  which supports "needed under mouse immune pressure but not for growth in resting human fibroblasts" and
  **not** "more important in mice than in humans": species, immune activation and scale all differ at once.
  The missing measurement is a genome-wide screen in IFN-gamma-activated **human** macrophages against matched
  unactivated human cells -- the human twin of the Wang design.
* **3Di completed to 92.2 % by alignment transfer (2026-11-04).** 520 of the 1,157 missing genes were filled by
  aligning the structure's own sequence to ours and copying the 3Di letter of each aligned residue (semi-global
  affine-gap alignment, identity scoring; accepted only at >= 95 % identity per structure and >= 50 % coverage of
  our protein; uncovered positions are `X`). Median identity 1.00, median coverage 0.93; 58 genes merged
  several structures. `tdi_source` ('aligned', 'aligned_multi') and the new `tdi_align_identity` /
  `tdi_align_coverage` columns must travel with any model that uses the 3Di columns: an aligned 3Di is weaker
  evidence than a measured one. 637 genes remain empty -- 373 whose best structure fails the rule (270 below
  95 % identity, i.e. a different isoform or gene) and 264 with no usable structure, which would need folding.
  Known gap: direct-source 3Di is still accepted on matching length alone, and two naive_db records have our
  length but not our sequence, so that kind of isoform can slip through.

## Einar's review lands (2026-11-04): the hard negatives now exist
`results/host_facing_truth_2026_11_04/`. He ruled on all 140 claims: 107 host-facing, 33 not, over 98 genes;
collapsed per gene by "strongest evidence wins" plus three explicit overrides (GRA60, GRA44, GRA47 -> not).
* The map now holds **65 host-facing** (48 at 0.95, 17 at 0.90) and **520 not**: 477 hyperLOPIT markers of
  compartments inside the parasite (0.80) and **43 vacuole-resident or curated negatives** (0.90-0.95), against
  9 before.
* **He rejected 22 dense granule proteins that the PVM proximity screen had called positive** (GRA1, GRA4,
  GRA22, GRA30, GRA36, GRA37, GRA40, GRA43, GRA46, GRA48, GRA49, GRA51-GRA59, GRA12B, ROP34), while accepting
  46 other proximity claims. That is exactly the distinction the first model failed: in the vacuole versus
  exported into the host. The raw proximity rows stay in `evidence.parquet` as the record and lose on
  probability.
* `conflicts.tsv` now lists 33 genes, up from 8, because a rejected proximity call counts as a disagreement.
  Where Einar's review is one of the two sources the conflict is **resolved**, not open.
* Next: rescore host-facing on these labels with `hosttx_effectors`, `organelle_surface_turboid`, `pubmed` and
  `pvm_proximity` held out, and report the secretory-restricted question as the headline.

## Second host-facing scoring run (2026-11-04): the bar fails, and run 1's pass was the leak
`results/host_facing_inference_v2_2026_11_04/`. 585 labels (65 positives, 520 negatives, of which **43 hard
negatives**); 398 features; the declared dependencies and the attention columns held out (47 columns); all 33
conflicted genes kept, because Einar's review resolves each one.

| Question | genes (pos/neg) | best precision | **lower 90 % bound** | reaches 0.90 | genes fillable |
|---|---|---|---|---|---|
| **Hard: secretory-pathway proteins** | 105 (58/47) | 0.822 | **0.727** | no | 0 |
| **Hardest: Einar's 98 reviewed genes** | 98 (66/32) | 0.673 | **0.605** | no | 0 |
| Easy: all labels (run 1's question) | 585 (65/520) | 0.776 | **0.681** | no | 0 |

* **Run 1's 0.906 does not reproduce.** The only changes are the reviewed labels and the held-out families, so
  that pass came from `hosttx_effectors` and the publication counts -- the near-circular features. Holding them
  out drops the same question from 0.906 to 0.681.
* **Measured against the hard negatives, the model is near chance.** Positives versus the 477 hyperLOPIT
  internal markers: AUC 0.984, false-positive rate 0.000. Versus the 43 hard negatives: **AUC 0.667, FPR 0.349**.
  Versus Einar's own 30 reviewed negatives: **AUC 0.586, FPR 0.467**. It tells inside-the-parasite from secreted
  and cannot tell stays-in-the-vacuole from exported.
* **No concentrated signal:** every permutation importance is <= 0.02 with an SD larger than its mean, spread
  across unrelated families. The calibrated [0.9, 1.0) bin holds 4 genes, 2 of them right. 0 of 14
  feature-set x training variants reach the bar.
* **More labels cannot fix this.** The Clopper-Pearson bound converges upward to the true precision, and the
  best precision reachable at any threshold is 0.822, below the bar itself. Extra labels would only measure
  0.822 more precisely. The binding constraint is the **evidence**, not the label count.
* **What would change the answer:** (1) a differential genome-wide measurement separating "in the vacuole" from
  "in the host cell" in the same cells -- selective permeabilisation, or proximity labelling anchored in the host
  cytosol against a PV-lumen-anchored control, with a per-gene FDR. The two existing genome-wide host-side
  assays cannot serve: the proximity screen defines the positives and the host-transcription screen is
  downstream of the conclusion. (2) SignalP or curated cleavage sites, so the secretory set is not defined by a
  predictor that calls GRA16 globular. (3) Only after a model clears 0.90 does the label count bind: at a true
  precision of 0.95 the bound needs ~65 held-out calls, about 150 labelled secretory genes against 105 today.
* **Decision: fill nothing for host_facing**, now for a measured reason rather than a suspicion.
* **Open for Einar:** GRA9 (TGME49_251540) and NTPase II (TGME49_277270). He accepted their PVM-proximity
  claims, while the curated gold-negative set (immunogold EM placing them in the IVN or PV lumen) says the
  opposite; the map resolves the 0.90-against-0.90 tie to "not host-facing". Flipping both still fails the bar
  (best bound 0.754), but the contested direction is the one that decides the question.

## Domain columns (2026-11-04), and the easy-contrast trap confirmed as a general pattern
`results/domain_traits_2026_11_04/`: **313 tier-2 columns** over the 8,140 genes -- 299 family flags present in
at least 10 genes (221 InterPro, 78 Pfam; the cutoff is visible at 5/10/20/50 in the notebook), 10 named classes
(protein kinase 173, WD40 129, SRS/SAG 110, RNA-binding 104, EF-hand 82, rhoptry kinase 50, AP2 45, ankyrin 30,
proteasome/ubiquitin 199, transmembrane-rich 668) and 4 counts. CDPK1, ROP18, AP2XI-5 and SAG1 are all caught by
the class they should be. 3,173 genes (39 %) carry no domain at all. Source: the ToxoDB InterPro gene-attribute
download, whose own registry note warns "domain annotation partly records study effort, not conserved
architecture"; a False flag means not found, not proven absent.

**`domain_novel` (no InterPro and no Pfam) is the clearest demonstration yet of the easy-contrast trap:**

| Comparison | host-facing with no domain | the other group | odds ratio | p |
|---|---|---|---|---|
| host-facing vs all genes | 46 of 65 (71 %) | 3,173 of 8,140 (39 %) | 3.79 | 2.9e-07 |
| host-facing vs all 520 negatives | 71 % | 72 of 520 (14 %) | 15.06 | 2e-21 |
| **host-facing vs Einar's 32 vacuole-resident negatives** | 71 % | **27 of 32 (84 %)** | **0.45** | **0.21** |

Against the whole negative set it looks spectacular; against the negatives that answer the question there is no
difference at all, and the vacuole-resident proteins are if anything *more* often domain-free. Structural
novelty separates secreted lineage-specific proteins from conserved housekeeping proteins -- not exported from
vacuole-resident. **Rule for every future column: report it against the hard negatives, never against "all
genes" or "all negatives".** This does not contradict the earlier foldseek result that host-facing domains are
structurally novel; it says novelty is a property of the secretory repertoire, so it cannot locate a protein
within it.

### Correction (2026-11-04, same evening): run 1's pass was the LABELS, not the leak
The leakage-attribution diagnostic (`leakage_attribution.tsv`) refits the easy question on the **new** labels
with the declared dependencies **put back**. Result: the bound moves from 0.681 to **0.694** -- nowhere near
run 1's 0.906.

| question | features | lower 90 % bound |
|---|---|---|
| easy | dependencies held out (run 2) | 0.681 |
| easy | dependencies put back (run 1 style) | **0.694** |
| hard | dependencies held out | 0.685 |
| hard | dependencies put back | 0.685 |

So the near-circular features were worth about **+0.013**, not +0.22. The earlier claim in this file and in the
commit message "the bar fails and run 1's pass was the leak" is **wrong** and is corrected here: run 1 passed
because **its labels were wrong**. It counted the 22 vacuole-resident dense granule proteins as host-facing
positives -- the easiest possible positives, since they look like every other GRA -- and scoring confident calls
on them as correct is what produced 0.906. Einar's review turned them into negatives, and the same confident
calls became errors.
`hosttx_padj` really was run 1's top feature, so the leak was real; it simply was not what passed the bar.
The practical lesson is sharper than the one first drawn: **a wrong label set flatters a model far more than a
leaky feature does**, and no feature audit would have caught it. Only the human review did.

## Function truth set (2026-11-04): 242 genes labelled, and the constraint is labels, not evidence
`results/function_truth_2026_11_04/`. The map held ~25 genes with an experimental function; it now holds **242**
with a coarse class, 145 of them at the two non-screen tiers.

| Tier | Probability | Genes | Source |
|---|---|---|---|
| observed-direct | 0.95 | 32 | experimental GO (QuickGO taxon 5811, 156 annotations, 85 terms) |
| observed-curated | 0.90 | 113 | reviewed Swiss-Prot curator FUNCTION comments (154 entries, all of Swiss-Prot Toxoplasma) |
| observed-screen | 0.80 | 97 | arrayed splitCas9 imaging screen (PMID 35538310), process claims only |

* 15 classes, all populated; 138 of 242 genes (57 %) carry more than one. Information content per label is
  recorded (`protein binding` 3.54 against `actin filament severing activity` 13.29), which is what stops a vague
  answer scoring like a specific one.
* Excluded on purpose, and this is why the set is small: InterPro/Pfam, ToxoDB product names, IEA/ISS/ISO GO,
  UniProt keywords and domain features, and CRISPR fitness (essentiality is not a function).
* Two guards on the text rules, both auditable: a negation cue in the same sentence (ASP1 "dispensable for
  bradyzoite viability" does not become a differentiation gene) and a host-context guard (ROP17 "promotes
  autophagy in human cells" is host modulation). 20 discards recorded.
* 183 PMIDs cited; a 30-id sample verified against NCBI, 30/30 resolve.
* **Unlike host-facing, the constraint here is the label count, not the evidence.** The lower 90 % bound needs
  ~65 held-out calls even at a true precision of 0.95, and the largest non-screen class has 37 genes, so **no
  class can be scored on its own yet**. What it can support is the pooled multi-label question plus the CAFA
  no-knowledge split. More curation therefore does help, which is why the review queue matters.
* `review_queue.tsv`: 254 genes, one row per gene, worst-first; 156 score >= 3, so the top 100 is a meaningful
  cut. Headed by genes with real evidence that no rule could classify (RON4, MIC21, RASP1/3, the enolases, GRA1
  — reviewed entries whose FUNCTION comment is empty).

### The strain-suffix mapping, checked
21 reviewed entries are the GT1/VEG/RH88 copy with no TGME49 cross-reference and were placed on
`TGME49_<same numeric suffix>`. `instructions/00_INDEX.txt` warns that **TGRH88 numbering does not correspond to
ME49**, so this needed checking. All 21 were cross-checked against the independent ToxoDB product name and all
21 agree: GC -> guanylyl cyclase, UGO -> guanylate cyclase organizer, POFUT2 -> O-fucosyltransferase, VIT ->
vacuolar iron transporter, ZNT -> zinc transporter, AMA1/AMA3, ATP4 -> P-type ATPase4, PIT/PT2 -> phosphate
transporters, DCX -> doublecortin, eIF6 -> binds the 60S subunit, CIAPIN1 -> cytosolic iron-sulfur assembly,
MIC2, RASP2 -> rhoptry apical surface protein, CDC50.1 -> traffics GC and UGO to the membrane (and GC and UGO
are themselves in the set). **Rule for future use: accept a strain-suffix mapping only when the product name
independently corroborates it**; the index warning stands for the general case.

## Binding partners scored for the first time (2026-11-04): 0.891, a near miss
`results/partner_inference_2026_11_04/`. Positives: the 2,842 StarPath XL-MS crosslink pairs over 1,630
proteins (probability 1 - local FDR). The 33,596 bait-member pull-down pairs were kept **out** of the headline,
because without a CRAPome background a sticky contaminant is indistinguishable from a real partner. No verified
non-interaction set exists, so every negative is a decoy and three constructions were scored.

| Question | Pairs | Lower 90 % bound | Reaches 0.90 |
|---|---|---|---|
| **Cold-node, degree-matched (headline)** | 1,248 | **0.891** | **no** |
| Cold-node, same-compartment decoys | 1,199 | 0.919 | yes |
| Cold-node, random decoys (easy) | 1,185 | 0.902 | yes |
| Known-node, degree-matched (easy) | 1,959 | 0.912 | yes |

* **The decoy construction decides the answer, exactly as for host-facing.** The "both proteins are popular"
  baseline scores AUC **0.995 against random decoys and 0.512 against degree-matched ones**: the popularity
  shortcut is real and degree-matching removes it. Model AUC falls from 0.933 to 0.788 the same way.
* **This is a near miss limited by held-out pairs, not by attainable precision.** Best precision with >= 20
  calls is 0.957, and a 0.90 bound at that precision needs ~52 held-out calls against the 47 available. Unlike
  host-facing, where attainable precision (0.822) was itself below the bar, here more pairs would settle it.
  Cold-node scoring keeps only 22 % of labelled pairs, which is inherent to holding out whole proteins.
* **What the model actually knows is co-localisation plus shared domains.** Family AUC alone: compartment 0.705,
  sequence 0.660, domain 0.581, fitness 0.548, co-expression 0.540, 3Di 0.502. There is **no per-pair
  measurement in the feature table at all**.
* **The StarPath fitness-difference finding does not reproduce**: the fitness family gives AUC 0.548 and
  removing it costs 0.0001 - consistent with the mouse-immune lesson that subtracting two relative scores
  amplifies noise.
* **AMA1-RON2 is a flat miss** (0.5th percentile): AMA1 is micronemal and RON2 rhoptry-neck, so a
  co-localisation predictor cannot find it. The MyoA glideosome does well (6 of 12 scorable pairs above the
  90th percentile). This is the clearest statement of the model's blind spot: the biologically interesting
  partners are the ones that cross compartments.
* **Prior shift:** the labelled set is 1:1 while true pair prevalence is tiny, so transferred probabilities are
  ranking scores, not posteriors. The 223,097 candidates above 0.5 must **not** be read as likely partners.
* **Decision: fill nothing for partners.** What would change it: the CRAPome background (it would turn 33,596
  provisional pairs into graded evidence and, more importantly, provide the independent assay this run lacks
  for held-out truth); reciprocal pull-downs or solved complexes for the first verified negatives; AF3 pair
  confidence on a shortlist as genuine per-pair evidence.

### The circularity worry, settled
The agent could not determine whether StarPath's pair FDR used localisation agreement, which would have made
the strongest feature partly circular. Checked here: same-compartment crosslink pairs do have a much lower FDR
(median 0.00000 against 0.00529, p = 6e-29) -- but **within bins of identification score the difference
disappears** (every bin p > 0.07), and identification score correlates -0.92 with local FDR. So the FDR is a
function of the mass-spectrometric evidence, and the association with compartment is explained by
same-compartment pairs yielding more crosslinks (median 2 against 1), which is what genuine interaction
predicts. **The crosslink labels are not circular with respect to compartment**, so the co-localisation feature
is legitimately predictive rather than leaked.

## Function scored (2026-11-04): fails the bar, but the learning curve is still climbing
`results/function_inference_2026_11_04/`. 242 genes x 15 classes = 3,630 (gene, class) pairs, 469 positive;
732 features in 86 families; 5 orthogroup-grouped folds; isotonic calibration (ECE 0.010).

* **No threshold from 0.20 to 0.99 reaches the 0.90 bound.** Verified against the tables: the best lower bound
  is **0.617** (at threshold 0.70) and the highest precision at any threshold is **0.760**, or 0.800 taking the
  top k of the raw ranking. (The agent's own summary quoted 0.525 and 0.680; the shipped tables give the higher
  figures, so its verdict was if anything conservative. The verdict is unchanged: the ceiling is below the bar.)
* **The model is real, just not precise enough**: +0.202 average precision over the best baseline (the class
  prior), beating it in 5 of 5 folds. Ortholog copy has coverage 0 under grouped folds, which is correct.
* **No class can be scored on its own.** Only 1 of 15 reaches even 22 held-out calls, the minimum at which the
  bound could reach 0.90 with perfect precision; 5 classes make no calls at all. `lipid_membrane_biology` is
  below chance (AUC 0.456).
* **The easy contrast appears here too.** Screen-only genes score AUC 0.893 against 0.701 for the non-screen
  genes, with all six `screen_*` columns dropped throughout -- so it is the screen's own gene selection, which
  no hold-out can remove. The hardest honest cut (non-screen genes, the 11 classes the screen does not support)
  gives 7 of 7 correct: a bound of 0.720 capped by arithmetic, not by errors.
* **The best news in the run:** on the 20 CAFA no-knowledge genes (no domain, no informative product name) the
  pooled AUC is **0.898** against 0.771 on the 222 annotated ones. The model is *not* worse on the blanks, which
  is the opposite of what the annotation-agreement worry predicted. Only 10 calls, so it cannot be filled on.
* **Not annotation agreement:** an annotation-free feature set (expression, fitness, proteome only) still
  reaches AUC 0.776 against 0.790 for everything. Domains alone 0.706, co-expression alone 0.663.
* **Labels bind here, unlike host-facing.** The learning curve, refitting on 40/60/80/100 % of labels dropped
  by orthogroup, is monotone and still rising: AUC 0.70 -> 0.79, average precision 0.26 -> 0.37, top-25
  precision 0.44 -> 0.68. **Actionable: a 65-call demonstration needs about 583 labelled genes against 242
  today, roughly 341 more.** Einar's 254-gene review queue is larger than that, and its top 100 alone takes the
  set to ~342 genes and the non-screen set from 145 to 245. Per class the need is far worse
  (`per_class_label_need.tsv`: 477 genes for `apicoplast_mitochondrion` to reach 65 calls).
* **Decision: fill nothing.** 22 unlabelled cells sit above a calibrated 0.90 and 2,240 above 0.50; all stay
  exploratory.
* Flagged and left in, each measured rather than assumed: `has_ec` is the second feature and an EC number is
  nearly the `metabolism` class; the protein-kinase domain is fourth while `signalling` is assigned by a
  `kinase` text rule; `n_fulltext` (attention) is eighth. Removing them costs little, and the no-knowledge and
  annotation-free results show they are not what carries the model.
* A negative pair means "not recorded", not "not true", so open-world annotation biases precision **downward**.
  That makes the failure somewhat pessimistic, but 0.76 is too far from 0.90 for it to matter.
