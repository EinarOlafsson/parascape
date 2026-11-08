# 65 · Audit dataset choices and retain a reproducible selection rule

Status: OPEN, 2026-10-07. User-requested priority alongside instruction 64.

User rule: when multiple datasets answer the same information slot, prefer the
publication with the highest citation rate per year, with a slight preference
for newer datasets and a preference for more comprehensive datasets.

Biological admission precedes ranking: organism, host/cell type, stage, assay,
quantity, units, measured versus inferred evidence, mapping, quality checks and
data access must fit the question. A highly cited review is not a dataset.
Incompatible cohorts remain separate. Popularity cannot establish ground-truth
accuracy. Unknown dates/counts/coverage remain unknown, not zero or guessed.

## Controlled action items

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| [65.01](../done/65_01_dataset_selection_policy.md) | ✅ | 0 h | Implement and test a versioned selection policy |
| [65.02](../done/65_02_publication_identity_audit.md) | ✅ | 0 h | Audit publication identity and citation rates for all 162 registered sources |
| [65.03](../done/65_03_slot_literature_discovery.md) | ✅ | 0 h | Search every current slot for literature challengers and expose comparison gaps |
| 65.04 | 12% | 14–28 h + compute | Review biological eligibility and deeper literature; validate/promote replacements and rerun affected checks |

65.02 depends on 65.01; 65.03 depends on 65.01/65.02; 65.04 depends on all three.
Each completion requires recorded evidence, relevant checks, a commit and nightly
publication. Repost the full progress table (instruction 64 plus these four rows)
when an item completes. Keep uncertain comparisons and unsuccessful replacements
visible rather than declaring every incumbent the best in the literature.

## Frozen audit scope

All 162 registry entries, all 297 slot declarations, current installed storage
coverage and candidate/refusal records. Host status is measured separately from
the pending host gene-space builders: seven human and three mouse dataset sources,
20,989/15,590 protein reference rows, 8/26 Tg and 4/24 Pf host slot views populated.
These are storage/slot counts, not host inference or benchmark completion.

Publication metadata comes from recorded primary-service queries and frozen
responses. Only identifiers verified by those responses may enter the audit's
resolved publication fields. Dataset accession links remain distinct from papers.
No wholesale data replacement follows from citation ranking alone; each proposed
replacement must pass its biological admission and affected validation checks.

## Completion evidence

65.01: `parascape/dataset_selection.py` implements the versioned, admission-gated
policy; `docs/dataset_selection.md` documents factors and unknowns. **389 focused
checks passed** for policy, docstrings and organism invariants. AGENTS.md and the
session handoff persist the user rule. Publication metadata and literature search
are next; no dataset replacements have been made by this policy-only step.

65.02/65.03: the frozen audit covers every source and slot; 95 recorded publication
identities resolved, 2,032 distinct literature publications retained, six diagnostic
quantity/context traps triaged, **576 relevant checks passed** (one optional pdoc
module skipped locally). See `results/dataset_selection_2026_10_07/README.md` and
the linked completion cards. **No dataset has been certified literature-best or
promoted**. 231 truncated slot searches, 28 zero-hit views and all source/deposit
lineage/assay-coverage decisions remain explicit in the review queue.

## Next controlled review batches (65.04)

1. Verify journal/preprint supplementary equivalence for the host rhoptry screen;
   retain first-publication/version citation lineage and original source provenance.
2. Review the 2017 quantitative RBC deposit and map its absolute copy numbers;
   retain the installed 2026 fraction PSMs as a distinct measured quantity.
3. Verify GTEx v10 release/publication lineage and actual cell-type fields, then
   FANTOM5 tissue/sample lineage and the host accession-only infection series.
4. Work through `source_review_queue.csv`, expanding truncated/zero-hit literature
   searches with accession-aware primary repositories and explicit matched assays.
5. For each admitted comparable alternative, record common measured coverage,
   mapping/QC, access, policy factors and a keep/add/replace decision. A replacement
   gets its own bounded action card and affected leakage/benchmark/release checks.

Popularity alone cannot close these review batches. The audit/search cards are
complete; the user's overarching best-dataset validation request remains open.

Host review update: 66.02 retains complementary RBC fraction/surface evidence and
keeps absolute copy counts as a separate candidate with measured group/mapping gaps.
66.03 verifies the rhoptry journal/preprint table exactly and upgrades the citation,
retaining first-publication lineage and exposing 39 ambiguous symbol projections.
66.04 now corrects the 39 ambiguous projections with original gene evidence intact;
expression-source replacement and host pack admission remain open.

Bounded truth-source partition **GT-SPATIAL-01** (65.04 with 64.10, depending on
64.03/64.04): review the two already registered primary spatial papers,
`lopit_tgon` (PMID 33053376) and `pf_spatial_proteome` (PMID 42218142).
Verify marker-field counts, exact installed-ID mapping and assay/label origins.
Retrieve only the primary named marker metadata, microscopy-method companions
and figures addressing experimental validation; retain original names and
primary PMC version/MD5/SHA receipts. The first preflight retrieved three
companions. The bounded extension uses figure addresses discovered in verified
article XML: Toxoplasma Figures 1/2 and Table S9, Plasmodium Figure 3.
No new paper, inferred identifier, replacement or biological admission follows
from this partition. Missing content becomes a recorded gap and the next source
proceeds. Final classifier marker sets can include the same paper's microscopy
and profile-based selections; they need row-level independence review.

Bounded host expression partition **GT-HOST-TX-01**: the existing
`host_gtex_transcriptome` and `host_mouse_tissue_transcriptome` records. Inspect
verified v10 tissue fields, exact installed-value reproduction and original
mapping receipts; check primary GTEx release metadata and the publicly listed
v11 median summary as a version candidate. Discover object addresses through
the public primary storage listing, verify provider checksum and retain original
filenames; no bulk per-sample/instrument files. Review E-MTAB-3579's original
assay/sample/age fields separately, preserving CAGE versus RNA-seq and adult
brain versus juvenile muscle distinctions. Acquisition and comparison do not
promote a runtime replacement; unavailable content is a gap while the other
source proceeds. Numerical replacement requires the existing affected checks.

GT-HOST-TX-01 evidence:
`results/host_expression_source_review_2026_10_08/` retrieves provider-checksummed
GTEx v11 summary/LCM README plus original Atlas `tpmss.tsv`; all Atlas parsed
values match the archived table. Full gene rows preserve PAR_Y. v11 has 74,628
genes versus 59,033, but only 11 additional native protein projections. Installed
replay exposes mapping/cutoff gaps; the documented 0.5 per-region filter reproduces
every overlapping legacy brain value. Mapping/coverage differences remain.
534 relevant checks passed; executed verification confirms rows/counts/checksums
and pending-admission policy refusal. No source declared literature-best or
promoted. Four additional source addresses reviewed across spatial/host partitions;
deeper literature/admission/promotion remains open (4%).

**GT-HOST-TX-02** depends on 66.04/64.21/64.22: original installed mapping lineage,
PAR_Y and explicit CAGE detection filter must be resolved, with raw gene evidence
intact, before corrected expression projections or v11 admission.

## Host documentary eligibility partition, 2026-10-08

The canonical `results/human_reference_terms_2026_10_08_v5/` independently pins
conditional GTEx public-data reuse terms, UniProt CC BY4/disclaimers and actual
GTEx V10→GENCODE39/GRCh38 reference association; all25 foundation/source receipts
remain unchanged. Future older-version packs require visible version notices and
source/date attribution. This is documentary eligibility evidence, not a new
citation-rate comparison, a source promotion, V11 association or an installed
human gene space. Release-to-publication lineage, measured version/quantity
mapping and best-source selection remain unresolved; 65.04 stays4%.

## 65.04, 2026-11-08: measured comprehensiveness, and batch 3 closed by 64.21

* **Comprehensiveness measured** (`scripts/build_dataset_comprehensiveness.py` ->
  `results/dataset_selection_2026_11_08/source_review_queue.csv`, a new queue version; the 2026-10-07 one is kept):
  for 140 of the 162 sources, the genes of their organism table with a value in any registered column, and that
  share -- median 0.80 (Tg 92 sources, median 0.80; Pf 48, median 0.85). 22 stay "unknown" with the reason: 12 have
  no registered column in the shipped tables, 10 are not gene tables of a shipped organism (host protein tables,
  source-only, metabolites). This is the "more comprehensive" term of Einar's selection rule, now a number.
* **Batch 3 (GTEx lineage)** is closed by instruction 64.21: the v10 release-to-GENCODE39/GRCh38 lineage and the GTEx
  public-data terms were verified from authoritative captures (`results/human_reference_terms_2026_10_08_v5/`), and
  the human gene space replays every GTEx value from source text.
* Batches 1-2 were closed under instruction 66 (rhoptry screen supplement, RBC deposit).

* **Batch 4, unresolved publications (2026-11-08, part):** of the 25, 10 sources now have their publication from
  the PRIDE project reference (`results/dataset_selection_2026_11_08/publication_resolution.tsv`; queue columns
  `publication_status`/`publication_id` updated). Two need a check before acceptance: morc_depletion (PXD058095 is
  linked to a bradyzoite-marker paper, PMID 40389644 -- confirm the RNA-seq workbook is from it) and
  pride_ubiquitination (PXD042937 linked by DOI only, 10.1128/MSPHERE.00511-23). Four corpus sources
  (bioid/ipms/pulldown corpora) are reclassified: many publications, cited per row. 11 remain unresolved: lopit_cpar,
  phosphosites, five ToxoDB dataset pages, plasmodb_pf3d7_expression, pf_febrile_stress, pf_sir2_perturbation,
  antisense_level -- the ToxoDB/PlasmoDB dataset pages sit behind the VEuPathDB login (see 38).

* **Batch 5, the slot decisions (2026-11-08):** `scripts/build_slot_decisions_notebook.py` ->
  `notebooks/slot_decisions_2026_11_08.ipynb`, `results/dataset_selection_2026_11_08/slot_decisions.tsv`,
  `slot_candidates_scored.tsv`, `slot_proposals.tsv`. 65.03 could rank **nothing** (`rankable_groups` 0 for all
  297 slots) because every candidate lacked a publication identity; the publications resolved in batch 4 make
  Einar's rule computable. Citation rates for the 10 newly resolved publications were fetched from Europe PMC
  (`citation_rates_resolved.tsv`); 105 of 162 sources now have one.

  Of the 60 slots with more than one candidate dataset: **27 separate cohorts** (the slot policy keeps them
  apart, so there is nothing to choose), **10 keep** (the installed default ranks first or is not beaten by the
  25 % margin), **7 additions proposed** (the slot has no installed default), **1 replacement proposed**
  (Pf lactylation: `pf_lactylome:has_lactyl` 6.5 against 5.0), **15 unknown** with the reason named per slot --
  mostly reference resources and the ToxoDB/PlasmoDB pages behind the login (38).

  Nothing is promoted automatically: the 8 proposals carry their margin and still owe the biological admission
  check (organism, stage, assay, units, mapping). That check is what remains of 65.04.

Still open: batch 4 (origin-paper verification for all 162; the 11 unresolved and 2 to confirm above) and the
biological admission check on the 8 batch-5 proposals.
