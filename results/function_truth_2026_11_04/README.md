# The function truth set (2026-11-04)

Built by `scripts/build_function_truth_notebook.py`; raw responses fetched by
`scripts/fetch_function_sources.py` and pinned in `sources/` with their retrieval times and sha256.
The executed notebook is `notebooks/function_truth_2026_11_04.ipynb`.

Instruction 72, target trait 3. The map previously held experimental GO function annotations for about
25 genes. This set labels **242 genes** with at least one coarse function class.

| tier | probability | genes | what it is |
|---|---|---|---|
| `observed-direct` | 0.95 | 32 | experimental GO (EXP/IDA/IPI/IMP/IGI/IEP), BP + MF, QuickGO taxon 5811 |
| `observed-curated` | 0.90 | 113 | reviewed UniProt (Swiss-Prot) curator `FUNCTION` comment |
| `observed-screen` | 0.80 | 97 | arrayed splitCas9 imaging screen phenotype (PMID 35538310); a process claim only |

**These three probabilities are provisional tier defaults** (instruction 72 v0), recorded as such. They
are not measured error rates. v1 replaces each one by how often that source is confirmed by an
independent direct source on the genes both cover, with a confidence interval, and no source keeps a
probability above what its replication supports.

138 of 242 genes carry more than one class (57%),
because proteins do several things; the classes are a set, never one winner.

## Excluded, and this is why the set is small
InterPro/Pfam domain assignments, ToxoDB `product` names, and every GO annotation with an IEA, ISS, ISO
or other electronic/similarity evidence code are **not used**, for either labels or classes. They are
predictions from sequence -- the thing the information map is supposed to predict. Scoring against them
would measure agreement with a predictor rather than with biology (instruction 70), which is the
failure both host-facing scoring runs diagnosed: run 1 cleared the bar only because its strongest
features were downstream of the conclusion. UniProt keywords and domain features and the CRISPR fitness
scores are recorded in `context_not_used_for_classes.tsv` and assign no class either; essentiality is
not a function.

## The label space
15 coarse classes, derived from the curated text and GO terms by the rules printed verbatim in
`class_rules.md`. The coarse classes are a **scoring convenience**; the truth is the fine-grained
evidence in `function_evidence.tsv` (GO terms with their codes and PMIDs, and the curator function text).
Information content is reported per label in `label_information_content.tsv`: `ic_corpus` for how common
a label is in this set and `ic_structural` for how specific a GO term is in the ontology, so that a vague
answer ("protein binding", structural IC 3.54)
cannot score like a specific one.

Largest classes: apicoplast_mitochondrion (80), cell_cycle_division (62), cytoskeleton (58), egress (42), invasion_motility (37), protein_trafficking_secretion (32).

## Two things a reader should know before trusting a row
* **Gene ids.** Many reviewed entries are the GT1, VEG or RH88 copy of a gene and carry no TGME49
  cross-reference. Those are placed on ME49 through the ToxoDB convention that the numeric suffix of a
  gene id is shared across strains (`TGGT1_254370` -> `TGME49_254370`). That is an inference from an id
  convention, labelled `D (ToxoDB strain-suffix, inferred)` in every table and penalised in the review
  queue. 1 reviewed entry and 1 accession remain unplaced or ambiguous (`unmapped_accessions.tsv`).
* **Two guards on the text rules**, both auditable in `negated_matches.tsv`. A match is discarded when a
  negation cue precedes it in the same sentence (ASP1 is "dispensable for ... bradyzoite viability", so it
  is not a differentiation gene), and a match of `transcription_chromatin`, `translation` or
  `proteostasis` is discarded when its sentence names the host (ROP17 "Promotes autophagy in human cells"
  is host modulation, which it already has, not parasite proteostasis).

Genes whose only function evidence is uninformative -- an experimental `protein binding` annotation, or a
reviewed entry with no FUNCTION comment at all (RON4, MIC21, the enolases, GRA1) -- get no class and head
the review queue, because that is a gap in the rules or in the evidence, not an absence of knowledge.

## Is this enough to score function prediction at 90 % precision?
Not per class. A Clopper-Pearson lower 90 % bound needs 65 held-out calls to reach 0.90 even if the
true precision is 0.95, and 368 if it is 0.92. At the two non-screen tiers
145 genes are labelled in total and 0 of
15 classes reach that count on their own -- before splitting into folds. What the set *can*
support is the pooled multi-label question (all classes at once, a few hundred held-out class calls), plus
the CAFA no-knowledge rule: genes with no prior annotation scored separately.

## Files
* `function_truth.tsv` -- one row per gene: classes, tier, probability, deciding source, reference.
* `function_evidence.tsv` -- long form: gene x source x GO term or curator text, evidence code, PMID.
* `function_classes_long.tsv` -- gene x class x source, the per-class tiers behind the per-gene row.
* `class_rules.md` -- the rules verbatim (GO seeds with resolved ids, regular expressions, screen columns).
* `coverage.tsv` -- genes per class per tier, with each class information content.
* `label_information_content.tsv` -- per-label specificity, coarse classes and fine GO terms.
* `review_queue.tsv` -- **one row per gene**, worst-first, for Einar's ruling on the top 100.
* `unmapped_accessions.tsv` -- reviewed entries with no single ME49 gene; curated function the map cannot place.
* `pubmed_verification.tsv` -- a 30-id sample of the cited PMIDs checked against NCBI E-utilities.
* `negated_matches.tsv` -- every text match the two guards discarded, with the cue and the sentence.
* `context_not_used_for_classes.tsv` -- UniProt keywords/domains and CRISPR fitness, recorded only.
* `evidence.parquet` -- information-map rows for the `function` trait.
* `sources/` -- every raw API response with URL, status, retrieval time and sha256 (`manifest.json`).
