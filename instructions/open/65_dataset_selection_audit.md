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
  `slot_candidates_scored.tsv`, `slot_proposals.tsv`, `slot_proposal_admission.tsv`. 65.03 could rank **nothing**
  (`rankable_groups` 0 for all 297 slots) because every candidate lacked a publication identity; the publications
  resolved in batch 4 make Einar's rule computable. Citation rates for the 10 newly resolved publications were
  fetched from Europe PMC (`citation_rates_resolved.tsv`); 105 of 162 sources now have one.

  **Only a `one` slot has a default that can be chosen or replaced.** `audit_dataset_selection` fills
  `default_patterns` for policy `one` alone; a `fill` or `average` slot uses every candidate it has and a
  `separate` slot keeps them apart. A first version of this pass read those empty default lists as missing
  defaults and called nine slots additions -- wrong, and corrected here before anything was installed.

  Of the 60 slots with more than one candidate dataset:

  | Decision | Slots |
  |---|---:|
  | nothing to choose: separate cohorts | 27 |
  | nothing to choose: every candidate used (fill 5, average 4) | 9 |
  | keep the installed default | 10 |
  | **replacement proposed** | **1** |
  | unknown, reason named per slot | 13 |

  The one proposal is Pf lactylation: `pf_lactylome:has_lactyl` scores 6.5 against 5.0 for the installed default.
  It passes the organism check and its name matches the slot context and unit (`slot_proposal_admission.tsv`), and
  Einar chose on 2026-11-08 to leave it until he has looked at it. Nothing was installed.

  The 13 unknown slots are mostly reference resources and the ToxoDB/PlasmoDB pages behind the login (38).

Still open: batch 4 (origin-paper verification for all 162; the 11 unresolved and 2 to confirm above); Einar's
look at the Pf lactylation replacement.

## 65.04 batch 4, origin-paper verification (2026-11-08)

`scripts/verify_origin_papers.py` -> `results/dataset_selection_2026_11_08/origin_paper_verification.tsv`; the
queue column `origin_paper_verified` is written from it. The check that can be automated: does the paper recorded
for a source actually name the deposit the data came from? Europe PMC is searched for the accession and the
recorded publication must be among the papers that name it.

| Verdict | Sources |
|---|---:|
| **verified**: the recorded paper names the deposit | 25 |
| supplementary table of a named paper: the recorded location IS the paper, nothing to cross-check | 92 |
| no accession recorded, so this check cannot apply | 30 |
| no publication recorded (the 11 of batch 4) | 11 |
| accession named in other papers, not the recorded one -- read before accepting | 2 |
| real deposit that no indexed paper names (PXD058095, GSE302107) | 2 |

**PXD058095 is the morc_depletion link flagged earlier**, and this check does not support it: no indexed paper
names that deposit, so the PRIDE reference is the only evidence. It stays flagged for Einar.

25 of 162 sources are verified by deposit. For the 92 supplementary-table sources the paper is the location, so
verification means reading the table against the installed columns -- that is reading work, not a query, and it
is what remains of batch 4.

## 65.04 batch 4, 2026-11-10: the check was asking the wrong authority

The verification ran the wrong way round. It searched Europe PMC full text for an accession and required the
recorded paper to be among the hits. That fails whenever the origin paper does not spell its accession in
indexed text, which is common in older papers and in consortium data -- so a correct record reads as suspicious.
Both sources flagged "read before accepting" were correct:

* **`E-MTAB-3579`** is the FANTOM5 mouse CAGE deposit, and the recorded 24670764 IS the FANTOM5 flagship
  ("A promoter-level mammalian expression atlas", Nature 2014). The four papers that name the accession are
  re-users of the atlas, not its origin.
* **`GSE132248`** records exactly the citation GEO lists for it (31955846).

`scripts/verify_origin_papers.py` now asks the **deposit first** -- GEO carries the citation a submitter
attached to a series, BioStudies carries it for ArrayExpress and PRIDE -- and only falls back to full text when
the deposit names nobody. A source the deposit and the record disagree on is reported as a disagreement with
both named, rather than as a missing paper. Where neither route can settle it, `READ_AND_DECIDED` holds the
answer WITH its evidence, so the verdict can be checked rather than trusted.

| Verdict | Sources |
|---|---:|
| verified by the deposit | 14 |
| verified by the paper text | 12 |
| verified by reading, evidence recorded | 1 |
| supplementary table of a named paper | 63 |
| no publication recorded in the queue | 40 |
| no accession recorded | 30 |
| no paper names the accession (`PXD058095`, `GSE302107`) | 2 |

**The committed verification table was also stale.** It described 121 sources as having a publication where the
queue holds 105; the queue moved on 2026-11-08 and the verification was never re-run. They agree again now.

**Publications found for the queue gaps**, from the deposits themselves, ready to fill `publication_id`:

| accession | paper | source |
|---|---|---|
| GSE108740 | 30728393 | GEO citation |
| GSE206344 | 41761089 | Europe PMC |
| GSE22258 | 20709297 | Europe PMC |
| GSE200962 | 35652638 | Europe PMC |
| GSE223620 | 37081202 | GEO citation |
| GSE313048 | 42778561 | GEO citation |
| GSE132237 | 31955846 | GEO citation |
| PXD042937 | 37975677 | Europe PMC, confirms the PRIDE reference |
| GSE277553 | none | no paper names it; the deposit declares none |
| GSE335016 | none | no paper names it; the deposit declares none |

Still open in batch 4: filling those publications into the queue, and the 63 supplementary-table sources, where
verification means reading each table against the installed columns -- reading work, not a query.

## 65.04 batch 4, 2026-11-10: the 63 supplementary-table sources, read

`scripts/verify_supplementary_tables.py` ->
`results/dataset_selection_2026_11_08/supplementary_table_verification.tsv`.

A supplementary-table source records a LOCATION, not a deposit, so the origin check could only say "the
location is the paper itself". Verifying one means opening the table and asking whether the genes we carry a
value for are the genes that table is about. The script does NOT re-derive the column -- reimplementing 39
derivations would test the script rather than the data -- it checks the identity link, which is where this kind
of source goes wrong.

| Verdict | Sources |
|---|---:|
| the table names the genes we ship | 33 |
| read: keyed on something other than gene accessions | 3 |
| not keyed on genes, so the check does not apply | 3 |
| the file was never kept on this machine | 21 |
| no file recorded, only a paper | 3 |

**Two of my own mistakes were in the first run, and both looked like data problems.** It picked the first
shipped table that merely CONTAINED the column name, so `pf_meltome` was compared against Toxoplasma genes
because both tables carry `melting_temperature`; and it compared strain-prefixed accessions literally, so four
CRISPR screens whose tables are entirely GT1 reported 0 % against ME49 genes that are the same genes. Ten of
the eleven flags disappeared once the table was chosen by organism and the strain taken off the accession.

The three that remain are explained in the script with the evidence inside each file, not left as flags:

* **`pf_hsp90_chemoproteome`** is keyed on a "UniProt ID" column, so the link to our gene-keyed column runs
  through the identity layer and cannot be seen by matching accessions.
* **`host_k562_rhoptry_screen`** is MAGeCK output keyed on gene SYMBOLS while the host table is keyed on
  UniProt accessions; its first sheet is a README, which is why so few identifiers appear at all.
* **`pf_resistome`** is one row per SNV in compound-selected clones, and two of its six installed columns come
  from the Pf6 field data instead, so this table is expected to cover only part of the column set.

**What cannot be finished here: 21 of the 63 source files are not on this machine.** They are recorded with
their paths and were presumably kept where they were acquired. Those 21 need the files before the same check
can run; nothing else about them is unresolved.

## 65.04, 2026-11-10: "the file is not on this machine" was an assumption, and it hid something

The 21 unverifiable supplementary sources all record a URL, so not having the file was a statement about work
not done rather than about the data. `scripts/fetch_missing_supplementary.py` fetches them from what the
registry already says: most URLs are a Europe PMC `supplementaryFiles` endpoint, which serves a ZIP of
everything a paper deposited, so the archive is opened and the member matching the recorded basename is taken
out. Where that basename is absent the members are listed rather than a guess being made -- picking "the
spreadsheet that looks right" is how the wrong table gets installed.

**What the first run found is more useful than the files.** `os.makedirs` raised `FileNotFoundError` on a path
it should have been able to create, and the reason is that **part of this dataset tree is symlinked to the
other machine**:

```
datasets/post_translation/acetylome -> /mnt/firecuda2/Claude/repo/starplast/datasets/post_translation/acetylome
```

**10 of the 21 sit behind such a link**, including the whole `datasets/quarantine` subtree. They are not
missing: they live on the home machine, and here the links dangle. Two consequences worth writing down:

* a directory that `ls -ld` shows as present can resolve nowhere, so "the file is absent" and "the path is
  broken" look identical from a script;
* had `/mnt/firecuda2` been mounted, creating those directories would have written **into the home machine's
  own checkout of this repository**. The fetcher now refuses to write through a dangling link, keeps the file
  beside its download and records the real location in `supplementary_fetched_paths.json`, so verification
  never depends on repairing a tree that points somewhere else.

Reached before the run was stopped (the box was needed idle for an unrelated GPU measurement):

| source | outcome |
|---|---|
| `host_erythrocyte_proteome` | placed, 822,887 bytes |
| `host_erythrocyte_surface` | placed, 725,839 bytes |
| `host_gtex_transcriptome` | placed, 8,846,936 bytes |
| `cyst_wall_interactome` | the archive does not contain `cyst_wall` |
| `host_macrophage_surfaceome` | the archive does not contain `S1_File.xlsx` |

Those last two are a **provenance finding, not a download failure**: the recorded filename is not in the
supplementary bundle the paper actually deposited, so what is installed came from somewhere the record does
not name. Both need a person to look before anything is re-fetched. (`host_macrophage_surfaceome` is the CSPA
source whose download links were already found to serve the site HTML rather than a workbook.)

**Three fetched, two verified at 100 %.** `host_erythrocyte_proteome` (5,264 of 5,264 installed genes in the
table) and `host_erythrocyte_surface` (267 of 267). The GTEx matrix is a fourth "keyed differently" case: it
is keyed on ENSG gene identifiers while the human host table is keyed on UniProt accessions, so the strings
never match and the link is the identity layer. The reader also did not recognise `ENSG` identifiers at all
until now, which is why it first reported the file as unreadable rather than as keyed on something else.

Running total: **35 of 63 verified by reading, 4 keyed differently, 3 not gene-keyed, 3 recording no file, and
18 whose file is still not here** -- 10 of those behind the symlinks to the other machine. One Europe PMC
download (`oxidative_stress_screen`, PMC8216390) stalled rather than failed and was stopped by PID after ten
minutes; it is worth a retry when the box is free.

## 65.04, 2026-11-10 (later): fetching completed, and the tally now

Two tooling gaps were found by using the tools rather than reading them. The fetcher recorded where it had to
leave a file it could not place, but the verifier never read that mapping, so a successfully fetched file was
still reported as absent; `resolve` now consults `supplementary_fetched_paths.json`, which makes the sources
behind dangling symlinks checkable **where they sit**, with the broken tree left alone. And the fetcher had a
socket timeout but no deadline -- one Europe PMC endpoint trickled bytes for ten minutes without tripping it
and stalled the pass until it was killed by PID. Reads are chunked against a 240-second deadline now.

Fetch outcome over all 21: **2 placed, 2 fetched and left beside their download** (the symlinked ones),
**5 where the archive does not contain the recorded filename**, and **9 downloads that failed at the source**.

| Verdict | Sources |
|---|---:|
| the table names the genes we ship | 36 |
| the file is still not here | 14 |
| read: keyed on something other than gene accessions | 5 |
| not keyed on genes, so the check does not apply | 3 |
| no file recorded, only a paper | 3 |
| read but no identifier found -- needs a person | 2 |

**`pf_complexes` joins the keyed-differently group on direct evidence:** its `Clusters` sheet has a column
literally named `uniprot` and carries host and parasite accessions together (P02549, human spectrin, beside
A0A5K1K8E0), because the table is host-pathogen complexes. Our column is keyed on PF3D7 gene accessions.

**Five archives do not contain the filename the record names** -- `cyst_wall_interactome`,
`host_macrophage_surfaceome`, `pf_isoforms`, `pf_secretome`, `splitcas9_imaging_screen`,
`thermal_shift_cetsa`. That is a provenance finding rather than a download failure: what is installed came
from somewhere the record does not name. It is the single most useful thing this pass produced and it needs a
person, not a retry.

**Two read but yielded no identifier** -- `pf_cdpk1_dependent_sites` (an .xls whose real header sits below a
title row) and `spacr_escrt_screen` (the record names a directory, not a file). Both need the file named
before a check means anything.

The 14 still absent now wait on Europe PMC answering: `pf_crosslink_ms`, `pf_ip_ms`, `pf_lactylome`,
`pf_palmitome`, `cyst_wall_interactome` and `oxidative_stress_screen` all time out, `pf_acetylome` returns an
HTTP error. Not a tooling problem any more.

## 65.04, 2026-11-10: correcting my own finding -- three of the five "provenance defects" were my bug

I wrote above that five archives do not contain the file the record names, and called it "the single most
useful thing this pass produced", needing a person rather than a retry. **That was wrong for three of the
five, and the error was mine.**

Looking at what the archives actually held, rather than at the fact that a lookup failed:

* `host_macrophage_surfaceome` and `splitcas9_imaging_screen` held `[Content_Types].xml`, `_rels/.rels` and
  `xl/worksheets/sheet1.xml`. **That is the spreadsheet itself.** An `.xlsx` IS a ZIP, so `is_zipfile` says
  yes to one, and the code was looking for the file inside itself.
* `pf_isoforms` held the article's figures and `12936_2025_5376_MOESM1_ESM.zip` -- the supplementary data one
  level further in. `pf_secretome` the same shape with Wellcome Open Research `.tgz` members.

`is_office_file` now tells an Office document from a bundle by looking for the OOXML part, and one level of
nesting is followed. Three of the five recovered immediately: `host_macrophage_surfaceome` placed (2.7 MB),
`pf_isoforms` and `splitcas9_imaging_screen` fetched.

**38 of 63 verified now**, up from 36, and the absent count is 13.

What is left of that finding is one case, not five: `pf_secretome`, whose article bundle carries `.tgz`
members rather than `.zip` and is not followed. `thermal_shift_cetsa` times out and its record names a
directory rather than a file, which is a separate problem.

**The lesson is the one that keeps recurring in this item.** A lookup that fails is not evidence about the
data until you look at what it found instead. The strain prefixes, the ENSG identifiers, the organism-wrong
table and now the Office files were all the same mistake: I read my own tool's silence as a statement about
the world. Each time, looking at the actual content turned a "finding" into a bug.

## 65.04, 2026-11-10 (final): 39 of 63, and a map that stopped forgetting

Two more fixes, both found by watching the number go the wrong way.

**Nested tarballs.** `pf_secretome` and others bundle the supplementary archive inside the article archive as
`.tgz`, not `.zip`. One level of tar is now followed, with the same check against a member escaping its
directory that the zip path has.

**The fetched-file map was rewritten per run, so it forgot.** `supplementary_fetched_paths.json` is the only
record of where a file ended up when the tree is symlinked away, and rebuilding it from scratch each pass
discarded everything an earlier pass had fetched. `pf_isoforms` and `splitcas9_imaging_screen` had verified an
hour earlier and went back to "file not on this machine" because a later run happened to time out on them --
while one of the two files was still sitting on disk. The map now MERGES, and an entry is dropped only when
the file it points at is gone, which is the rule `add_deposits.py` already applies to a losing merge.

That regression is the useful part. **A count that goes down is information**; had I only looked at the
headline it would have read as "38, near enough 39" and two sources would have quietly reverted.

| Verdict | Sources |
|---|---:|
| the table names the genes we ship | 39 |
| the file is still not here | 9 |
| read: keyed on something other than gene accessions | 5 |
| not keyed on genes, so the check does not apply | 3 |
| no file recorded, only a paper | 3 |
| read but no identifier found | 3 |
| no installed column to compare | 1 |

From 33 verified and 21 absent at the start of this pass, to **39 verified and 9 absent**. The nine that
remain fail at the source, and `pf_ip_ms` proved those failures are transient rather than permanent: it timed
out twice and then fetched 8.3 MB on the third attempt. A later retry is worth more than a diagnosis.

## 65.04, 2026-11-10: none of the five "provenance defects" was real

`pf_secretome` was the last survivor of the five archives I reported as not containing the file the record
names -- the one I said was "a real record problem rather than five". The tar support found it one level in,
as a `.tgz` member of the Wellcome Open Research article bundle, 38,581 bytes. **So all five were tooling, and
the provenance finding I reported does not exist.**

That is worth stating as plainly as the original claim. I wrote it up as "the single most useful thing this
pass produced", and acting on it would have meant looking for a defect in five sources that are intact.

**40 of 63 verified.** Seven files remain absent and those failures are persistent, not transient: eight
downloads failed again on a third pass (`cyst_wall_interactome`, `oxidative_stress_screen`, `pf_acetylome`,
both `pf_berghei_*`, `pf_lactylome`, `pf_palmitome`, `thermal_shift_cetsa`), so the earlier note that a retry
beats a diagnosis was also too broad -- it was true of `pf_ip_ms` and of `pf_secretome`, and false of these.

Final tally for the pass: from 0 verified and no reading done, to **40 verified, 5 keyed on something other
than gene accessions, 3 not gene-keyed, 4 read without yielding an identifier, 3 recording no file, 1 with no
installed column, and 7 whose file this machine cannot fetch.**

## 65.04, 2026-11-12: "slow" and "failed" were the same word, and that hid one of seven

I called the seven unfetchable sources persistent failures at the source, on the evidence that they failed
three passes in an hour and again a day later. Probing the endpoints directly instead of trusting my own
error message: `PMC7002340` and `PMC4728587` both answer **HTTP 200 and stream throughout**, at roughly
40 KB/s. A 10 MB bundle therefore needs four minutes of perfectly healthy transfer.

`STALL_SECONDS` was a **total-time deadline wearing the name of a stall detector**. It killed a transfer that
was merely slow, raised `TimeoutError`, and I read that as the source being broken -- for three days, across
four passes. It is now a real stall detector: give up when nothing ARRIVES for 120 s, with a 30-minute ceiling
so one endpoint still cannot hang a pass.

With that fix `cyst_wall_interactome` downloads and gets far enough to open its archive. Six others still fail
and now fail *informatively*: five stall with no bytes for two minutes, and `thermal_shift_cetsa` raises
`IncompleteRead`, the connection dropping mid-transfer. Those six are genuine source failures. The
distinction between "slow" and "stalled" is the thing that was missing, and no amount of retrying would have
produced it.

**The pattern, for the record, since it is now five for five.** Every "finding" I reported in this item that
later proved false had the same shape: my tool failed, and I described the world instead of the tool. Strain
prefixes, ENSG identifiers, the wrong organism's table, Office files that are ZIPs, and now a timeout that
measured the wrong thing. In each case the correction came from looking at what the tool actually returned.

## 65.04, 2026-11-12: the cyst-wall source is identified, and the record names a directory

Fixing the timeout let `cyst_wall_interactome` download, and its archive answers the question the record could
not. `PMC7002340` (mBio 02699-19) holds 19 members, of which one is the data: **`mBio.02699-19-sd001.xls`**,
566 KB, sheets `Group 1`, `Group 1_raw`, `Group 2`, `Group 2_raw`.

Checked rather than assumed to be the right file: **all 56 genes carrying `cyst_wall_max_spectral` or
`cyst_wall_n_baits` appear in it, 56 of 56.** It is the source.

The registry records `datasets/quarantine/2026_08_16_unverified/Tg/cyst_wall/` -- a directory, with no file
named -- which is why this source could never be checked. The path should name
`mBio.02699-19-sd001.xls`; it cannot be written there yet because the whole `datasets/quarantine` subtree is
one of the symlinks pointing at the other machine. The verified file is recorded in
`supplementary_fetched_paths.json` so the check runs from where it actually sits.

**41 of 63 verified.** Six files remain unfetchable, and those six now fail for a reason rather than for a
timeout that measured the wrong thing: five stall with no bytes arriving for two minutes, one drops the
connection mid-transfer.
