# 71 · What the aggregated omics are for: four products on the uniform RNA and MS data

**Status: open, planning (requested 2026-11-04).** Work can start on the runs already verified while
downloads continue (instruction 69); every product is rebuilt as the queue fills, never hand-patched.

## The request (Einar, 2026-11-04)
* "While all this data is downloading let's discuss what useful things we could do with this data."
* Chosen, in this order: (1) co-expression → localization, (2) perturbation atlas, (3) host response
  atlas, (4) Toxoplasma CRAPome. All four.
* Context: ToxoDB hosts a curated subset of RNA-seq (genome-aligned, parasite-only, per-dataset gene
  values) and peptide evidence from submitted results. It has no cross-study matrix, no host side and no
  uniform proteomics re-search; its bulk downloads have needed a login since release 71. Overlapping
  datasets are a cross-check for our values, not a substitute.

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

## 4 · Toxoplasma CRAPome (background for pull-down MS)
* **Build:** from the uniform Sage results, every affinity-purification, BioID/TurboID or APEX
  project (classified from PRIDE metadata and sample names, checked by hand where unclear); per protein,
  the fraction of independent baits/projects in which it is detected, and its typical abundance.
* **Use:** a specificity score for any new IP-MS hit list (e.g. the EAF1 and ATG2 pull-downs): frequent
  flyers down-weighted, bait-specific hits kept. Comparable to the human CRAPome, which does not exist
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
