# 71 · What the aggregated omics are for: four products on the uniform RNA and MS data

**Status: open, planning (requested 2026-11-04).** Work can start on the runs already verified while
downloads continue (instruction 69); every product is rebuilt as the queue fills, never hand-patched.

## The request (Einar, 2026-11-04)
* "While all this data is downloading let's discuss what useful things we could do with this data."
* Chosen, in this order: (1) co-expression → localization, (2) perturbation atlas, (3) host response
  atlas, (4) Toxoplasma background frequency. All four.
* Context: ToxoDB hosts a curated subset of RNA-seq (genome-aligned, parasite-only, per-dataset gene
  values) and peptide evidence from submitted results. It has no cross-study matrix, no host side and no
  uniform proteomics re-search; its bulk downloads have needed a login since release 71. Overlapping
  datasets are a cross-check for our values, not a substitute.

## The end goal: one information map, three target traits (Einar, 2026-11-04)
* "In the end what I am most interested in is generating an information map where all current information is
  taken into account; this map I want to use to predict all the missing facts in the table."
* Target traits, in Einar's order. Localization (product 1) is a means for trait 1, not the end:

| Trait | Truth to score against (graded as in 70) | Evidence layers that should predict it |
|---|---|---|
| **1. Host-facing** (secreted beyond the parasite: PVM, PV lumen exported, injected into host cytosol/nucleus) | Direct: published host-cell/PVM localization by IFA/IEM, host-side BioID/TurboID hits, exported-protein screens; curation: dense granule / rhoptry hyperLOPIT markers (organelle is not the same as host-facing) | signal peptide and TEXEL-like motifs; co-expression with known secreted proteins (product 1); in vivo vs in vitro fitness; host response atlas (product 3: deleting it changes the host); fold novelty (host-facing domains are structurally novel, foldseek screen) |
| **2. Binding partners** (parasite–parasite and parasite–host, especially for host-facing proteins) | Direct: crosslinks (StarPath XL-MS), reciprocal pull-downs, solved complexes; curation: background frequency-filtered co-purification | co-purification across all uniformly searched pull-downs with background removed (product 4; host proteins are already in the search database as `HOST_`); co-expression and co-perturbation (products 1–2); co-localization; AF3 pair confidence for shortlisted pairs; ortholog interactions (Plasmodium) |
| **3. Function** | Direct: experimental GO (EXP/IDA…), phenotypes; curation: reviewed annotation | perturbation atlas (product 2: what changes when it is lost); co-expression modules; fold-based annotation (319 proteins annotated by structure); essentiality and condition-specific fitness |

* The map is one model over all layers (the ToxoNet v2 multi-task design: observed labels only, inferred evidence
  down-weighted), but every layer first earns its place on the 70 scoreboard on its own, so the combined model can be
  compared with its best single layer.
* Pairs (trait 2) need pair-level folds: hold out whole proteins (cold-node split, `splits.cold_node_pairs`) so a
  partner is never predicted from the same protein seen in training.

## Rules shared by all four (from 70, plus omics-specific)
* Scoring follows 70: graded truth, whole-source hold-out, orthogroup-grouped frozen folds, best-baseline
  comparison, the deployment rule. A product is "deployed" only through that rule.
* **Studies, not runs, are the unit of independence.** Hold-outs are by whole study; a gene module or a
  signature must replicate in studies it was not built from.
* Values come from the uniform reprocessing only (salmon/kallisto, Sage). Author-deposited matrices are a
  check, never mixed into the same feature.
* Host and parasite are never normalised together; each side is scaled within its own reads.
* Single-cell collections (instruction 69 exclusions) and the unreleased ToxoHub harvest may inform
  hypotheses and checks; the harvest never ships in a release.
* Normalization stays Einar's open decision (instruction 69 plan doc); each product records which
  normalization it used so it can be rebuilt.

## 1 · Co-expression → localization and function (first; feeds 70 directly)
* **Why:** genes of one compartment are made together in the cell cycle and respond together to
  perturbation (rhoptry, microneme, IMC waves). Robust cross-study modules are a localization feature
  that is independent of hyperLOPIT, which is what the 70 whole-source hold-out needs.
* **Build:** per study, within-study gene–gene correlation on log-TPM (parasite-only, studies with ≥ 6
  samples); combine across studies by rank aggregation (mean rank of correlation per pair), so no single
  large study dominates. Record per pair how many studies support it.
* **Evidence blocks registered in 70:** k-nearest co-expression neighbours' labels; module membership;
  peak phase from the ToxoHub cell-cycle clock (check-only while unreleased).
* **Scoring:** localization target of 70 (718 curation markers, direct tier as check), study-grouped
  plus orthogroup-grouped folds; baseline = ortholog copy and single-block kNN. Also function (GO BP
  experimental) as a second target.
* **Deliverable:** module table, per-gene neighbour list with support count, scoreboard rows.

## 2 · Perturbation atlas ("which known perturbation does my mutant look like")
* **Build:** for every study with a design contrast (sample registry, instruction 69: knockout,
  knockdown, drug, strain, stage), a within-study effect-size signature (log2 fold change with
  shrinkage, plus its standard error), never pooled raw values across studies.
* **Query:** a user signature (or a study held out) ranked against all signatures by weighted rank
  correlation; report the top matches with their study and contrast.
* **Validation:** replicate contrasts across independent studies (same gene knocked out in two labs, the
  same drug) must find each other; report the recovery rank. Regulator targets cross-checked against the
  Gubbels ToxoHub knockdowns (check-only).
* **Deliverable:** signature store, query function, recovery benchmark.

## 3 · Host response atlas
* **Build:** host-side counts from the competitive references (ME49_host and the host-specific ones),
  restricted to infected-vs-mock or mutant-vs-parent contrasts in the same host cell type; signatures as
  in 2, per host species, with human–mouse orthologs for comparison.
* **Questions it answers:** which parasite genes, strains or effectors change which host pathways;
  whether a host response replicates across cell types and labs.
* **Validation:** known effector phenotypes (e.g. GRA15 → NF-κB, ROP16 → STAT3/6, GRA24 → p38,
  TgIST → STAT1 repression) must appear in the matching contrasts; report which are recovered.
* **Deliverable:** host signature store, effector × host-pathway table with replication counts.

## 4 · Toxoplasma background frequency (background for pull-down MS)
* **Build:** from the uniform Sage results, every affinity-purification, BioID/TurboID or APEX
  project (classified from PRIDE metadata and sample names, checked by hand where unclear); per protein,
  the fraction of independent baits/projects in which it is detected, and its typical abundance.
* **Use:** a specificity score for any new IP-MS hit list (e.g. the EAF1 and ATG2 pull-downs): frequent
  flyers down-weighted, bait-specific hits kept. Comparable to the human background frequency (Mellacheruvu et al. 2013), which does not exist
  for Toxoplasma.
* **Validation:** known complexes (e.g. MyoA motor, AMA1–RON, MIC2–M2AP) must score as specific in
  their own pull-downs and not be removed as background.
* **Deliverable:** background frequency table with project counts, scoring function.

## Order and dependencies
1. Sample/context registry coverage for verified runs (instruction 69) → contrasts for 2 and 3.
2. Product 1 on the verified parasite-side runs; register its blocks in 70 and run the scoreboard for
   localization. This also restarts 70.
3. Products 2 and 3 share the signature code; 3 also needs per-host contrasts.
4. Product 4 needs the MS queue past its AP/BioID projects; classification can start now.

## background frequency first pass (2026-11-05): the frequency filter fails, the specificity score works
`results/pulldown_background_2026_11_05/`. 63 of 203 projects finished; 61 classified (25 affinity purification,
9 proximity labelling, 10 whole proteome, 7 PTM enrichment, 2 fractionation, 8 other, which are kept **out** of
the background because one misclassified whole-proteome run would raise every frequency at once). 23 of 61
classifications hinge on a judgement call and are flagged with the cue that decided them. Background: 4,449
parasite proteins over 34 experiments; detection = >= 2 PSMs at 1 % FDR mapping uniquely to one gene
(Spearman 0.95 against both a looser and a stricter rule, so nothing rests on the threshold).

**The frequency column must NOT be used as a filter.** Verified: the known complexes sit inside the
frequent-flyer range, not below it.

| | frequency range |
|---|---|
| Frequent flyers: actin 0.97, HSP70/90 0.76-0.94, cytosolic ribosome 0.79-0.91, GAPDH 0.76-0.88 | high, as expected |
| **Known complexes: MyoA glideosome 0.38-0.88, AMA1-RON 0.76-0.85, MIC2-M2AP 0.53-0.65, MORC/HDAC3 0.50-0.71** | **the same range** |

Separation AUC is only 0.743, and at a 50 % cutoff **12 of 13 complex members would be deleted**. No threshold
removes half the flyers while keeping 85 % of the complexes. The cause is biological, not a bug: 34 experiments
over 11 distinct baits, nearly all pull-downs from whole parasite lysate, so the abundant invasion and motility
machinery is in every list. The human CRAPome (Mellacheruvu et al. 2013), the published resource ours imitates, had over 400 experiments. **Honest statement: not enough
independent baits yet.** Enolase is additionally underestimated (0.24) by the unique-mapping rule on close
paralogues, which is recorded rather than patched.

**The specificity score (log2 abundance fold + log2 rarity) is usable for ranking today**, and the
abundance-fold term does the work. Bait recovery in its own pull-down, verified:

| project | bait | specificity rank | percentile | rank by raw abundance |
|---|---|---|---|---|
| PXD016383 | MYR1 | **1** of 363 | 100 | 1 |
| PXD038848 | TgBDP1 | **1** of 244 | 100 | 1 |
| PXD012103 | HCE1 | **1** of 334 | 100 | 11 |
| PXD039896 | C3 | 6 of 1,167 | 99.6 | 48 |
| PXD032360 | GRA28 | 14 of 428 | 97.0 | 4 |
| PXD032319 | ROP1 | 14 of 38 | 65.8 | 3 |

The tagged protein is in the top 5 % of its own hit list in 8 of 11 affinity purifications (median rank 14 of
~580, against 28 by abundance alone). The misses are correct answers: a TurboID ligase labels neighbours rather
than enriching itself, and PXD039922 is a SPARK auxin-degron IP in which the bait is deliberately degraded.
CDPK1 (PXD039432) is not recovered in its own pull-down at all and is worth a look.

**Practical guidance:** rank a new hit list by `specificity`, read `frequency` as context only, and treat
`median_abundance` as the useful column. For the 33,596 bait-member pairs this does not yet give a graded
probability by frequency, but it does give each pair an abundance-fold number it did not have. **Baits are the
limiting factor: 13 of 34 experiments have a confident bait, 10 more have a symbol absent from every ToxoDB
product name, and 7 have a complex or non-protein bait.** Identifying more baits is worth more than finishing
the remaining 140 projects.

## 5 · Condition fingerprints: annotate every dataset with metadata nobody published (Einar, 2026-11-05)
**The idea, in his words:** take a dataset that samples a condition axis in distinct samples (cell-cycle stages,
for example), build a fingerprint from it, then score every other dataset against that fingerprint. "That way we
not only compare, we annotate new metadata not known before to the datasets." On both the RNA and the protein
level. **Queued as a future product; nothing built yet.**

### Why it is sound, and the one refinement that matters
It is signature projection plus deconvolution, and it has precedent in this organism: ToxoHub transferred
cell-cycle phases to its own cells from the Smart-seq2 dataset GSE145080 by joint embedding and weighted vote.
**The refinement: a bulk sample is a mixture, so the output should be a proportion, not a label.** An asynchronous
tachyzoite culture contains every cell-cycle stage at once, and an alkaline-stressed culture is only partly
converted. "68 % tachyzoite, 32 % bradyzoite" is both more honest and more useful than "bradyzoite", and
deconvolution gives it directly. Estimated stage composition also becomes a covariate that explains differences
between studies which are currently attributed to biology.

### Rules it must follow (these are the traps already paid for in this project)
* **Score by within-sample ranks, not absolute values** (ssGSEA-style). A fingerprint learned in one study
  otherwise scores that study's batch, and rank scoring is what makes one fingerprint usable across platforms
  and between RNA and protein at all.
* **A fingerprint annotation is inferred, never observed.** It enters the map as `filled:<fingerprint>@<version>`
  or stays exploratory, and the fingerprint's source study must be declared in `depends_on`. Otherwise the
  obvious failure follows: annotate a dataset with an inferred stage, then use that stage as evidence for a
  trait, and the model learns from its own guess.
* **Validate by held-out studies that report the condition themselves.** Build the fingerprint on one study,
  score a study that states its own stages, and measure recovery — the instruction-70 pattern. Report accuracy
  per axis, and expect RNA -> protein transfer to be weaker than RNA -> RNA; score the two separately rather
  than assuming.
* Studies are the unit of independence, and sample quality weights the fingerprint as it weights everything else.

### Candidate axes and the data already on disk
| Axis | What we have | Worth |
|---|---|---|
| **Cell cycle** | 28 columns: a synchronised 12-hour series with async and blocked controls, plus ~5,377 periodic genes and the ToxoHub 24-window clock (unreleased, check-only) | Strongest. Phase-resolved co-expression is far more specific than overall co-expression, which would help the partner trait directly |
| **Bradyzoite differentiation** | 33 columns over several studies: tachyzoite/bradyzoite, bradyzoite subtypes A-E, alkaline stress, in vivo 28/90/120 days, brain time course, plus BFD1 and MORC knockdown series (18 more) | Strong, and multi-study, so it can be validated across studies rather than within one |
| **Oocyst and sporulation** | 6 columns: unsporulated, sporulating, sporulated | A clean three-point axis |
| **Merozoite and enteroepithelial** | 10 columns, including an enteroepithelial-versus-tachyzoite contrast | Covers the cat-stage biology nothing else reaches |
| **Extracellular versus intracellular** | 12 columns, RNA and Ribo-seq | An egress-state fingerprint |
| **Translational state** | 46 Ribo-seq and translational-efficiency columns | Orthogonal to transcript level: a protein can be made without its mRNA changing |
| **Nutrient and pH restriction** | 24 columns: pH 7.4/8.2, 2 and 5 weeks, serum and lipid restriction, glucose and glutamine withdrawal | Culture-condition covariates that currently masquerade as biology |
| **Host cell density** | 18 columns, confluent versus subconfluent | Same |
| **Strain type** | GT1, ME49 and VEG references already quantified per run | Type I/II/III without new data |
| **Protein level** | 12-column proteome time course, 8-column iTRAQ oocyst set, and the uniform Sage LFQ tables as the 203 projects finish | Lets the same axis be scored on protein, and tests RNA -> protein transfer |

### What it buys that nothing else does
The ~2,200 RNA runs and 203 proteomics projects would each gain estimated cell-cycle composition, differentiation
state and culture condition — **metadata the depositors never published**. That is new annotation rather than a
re-description, and it is the one product here that improves every other one: contrasts become comparable once
stage composition is a covariate, and a pull-down's lysate composition may explain the background frequency variability that
made the frequency filter fail.

## Product 1 through the harness (2026-11-06): co-expression alone places the division organelles

`results/coexpression_components_localization_2026_11_06/README.md`. Label-free co-expression components from 57
studies (735-run snapshot), scored per coarse niche through the harness: median cross-validated AUC **0.897**
against 0.934 for all allowed evidence; rhoptries 0.975, apical 0.958, micronemes 0.936 and IMC 0.913 match or
beat everything else together; plasma membrane 0.669 and ER 0.780 are where co-expression has little to say.
This is the block 70 asked product 1 to register; it ranks well, and like every localization question it cannot
fill until the frozen partitions hold more positives.

## Product 2 built (2026-11-07): the perturbation atlas finds the same perturbation in other labs

`parascape/perturbation_atlas.py` + `results/perturbation_atlas_2026_11_07/README.md`. Signature store = the
within-study contrast columns of the expression traits release; query = Spearman, own study held out; recovery
benchmark on four curated cross-study perturbations: median best-partner rank 1-3 against chance 59-99 (BFD1
stabilised: rank 1 both ways). Tests in `tests/test_perturbation_atlas.py`. Not done: weighting by the standard
errors, the ToxoHub knockdown cross-check (unreleased), and a module in the application.

## Product 3 built (2026-11-07): the host response atlas

`results/host_response_atlas_2026_11_07/README.md`: 106 host-side signatures over 18 studies from the same
contrasts as the parasite side. The MYR1 check holds: without MYR1 the infected-cell host response is half the
size (173 against 344 genes changed 2-fold) and only partly the same (Spearman 0.536). Open: orthologue mapping,
a query, mock-infected contrasts. Product 4 (background frequency) is instruction 80.

## Product 3, the query (2026-11-08)

`perturbation_atlas.host_store()` serves the host response atlas to the same `query` and `recovery` as the parasite
signatures (human and mouse never meet: they share no gene identifiers); `replicate_groups` finds perturbations
stated by at least two studies in the same host species (missing or baseline labels name nothing). Tests in
`tests/test_perturbation_atlas.py`.

**Result, stated plainly:** on the host side only one perturbation is stated by two independent studies -- deltaBFD1
in human cells (GSE132238, GSE139537). Each signature best other-lab partner ranks 8th-66th (median 21) against 16-40
expected by chance, never in the top 5: the host response to a bradyzoite-differentiation knockout is reproduced
weakly across labs, unlike the parasite side where perturbations find their partner at rank 1-3. The query works;
the host atlas needs more replicated perturbations (mock-infected and IFN-gamma contrasts) before its recovery can
be judged. Open: orthologue mapping between human and mouse signatures.

## Product 3, crossing the species (2026-11-08)

The last open point of product 3. `perturbation_atlas.host_symbol_map()` keys human and mouse host genes by the
gene symbol they share (human from `hs_genes.parquet`, mouse from the mouse protein-gene links joined to the
shipped mouse host table), and `cross_species(store)` re-indexes the whole atlas on that symbol.
`scripts/build_host_cross_species_notebook.py` -> `notebooks/host_cross_species_2026_11_08.ipynb`,
`results/host_cross_species_2026_11_08/`.

* **What the crossing costs, counted:** 31,151 symbols, **10,881 measured in both species**; 3.1 M gene-signature
  values dropped for want of a one-to-one symbol. It is a symbol mapping, not tested orthology, and a gene family
  whose member count differs between the species is dropped rather than guessed.
* **One perturbation is stated in both species: the BFD1 knockout.** gra18, LPS and BFD2 appear in the mouse
  signatures only, so they cannot be crossed.
* **Result: the host response to the BFD1 knockout crosses the species.** The mouse deltaBFD1 signature
  (GSE223621) finds a human deltaBFD1 signature as its **best partner of 99**, against a chance rank of 50.
  Spearman is low (0.203) -- the shared part is a minority of each response -- but the ordering is unambiguous.
* Contrast with the same-species result above: within human cells the deltaBFD1 partner ranks 8th-66th. Across
  species it ranks first. The host response has a reproducible core that a mouse cell and a human cell share,
  and a much larger species- or lab-specific part.

Tests: `test_the_symbol_map_is_one_to_one_within_each_species`,
`test_crossing_the_species_puts_human_and_mouse_on_one_index` (6 pass).

Product 3 has no open points left. Products 1-4 of this instruction are built; product 4 continues as 80.
