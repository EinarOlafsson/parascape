# Next session: start here

**Selected original-measurement context (2026-10-30):** Current dataset goal
active/incomplete; previous goal turn was progress (row-wise admission and native
repeated-key QC). New `Inspect selected measurement…` button/double-click view in
`parascape/deposited_values_panel.py` shows the exact original row/column, literal
source token beside numeric value, all row annotations, column/submitted-factor
metadata, qualified identifier evidence and source/treatment context. One quantity
column need not be one biological sample. No sample/control/replicate/gene/scale
identity or biological quality is inferred. Missing metadata stays unavailable.

The selection is keyed by original result position stored with its table item,
so sorting cannot substitute a different repeated label or column. Sorting is
paused while rows populate to prevent cell mixing. Detail trees are read-only;
old dialogs close/delete on selection changes or query/source/namespace/rule/
refresh invalidation. Values, source hashes, treatment records and exports stay
unchanged. Existing whole-source QC remains unfiltered. This is a nightly source
follow-up, outside the public0.0.6.0 wheel/asset; new local139-table data also
awaits a future portable release. Do not modify frozen6.

Focused reader/QC/source-filter/flag/identifier/UI checks plus actual whole-window
tooltip audit passed46 checks. New selected-context tests cover repeated labels,
row/column sorting identity, original metadata/leading-zero replicate labels,
missing numeric versus empty token, duplicate factor/evidence occurrences,
no-selection/missing-annotation states and six stale-context invalidations.
Independent actualQt proof visits all139 locally available quantity tables with
NAS/network blocked and writes to input measurements forbidden. All selected
rows/columns/tokens/annotations and every original output hash rechecked; values
and treatment records unchanged. Proof and screenshot at
`/media/carruthers/mnt1/parascape-measurement-context-validation-20261030-v2`.
Source implementation SHA861f0378 matches the final validated proof. Actual
rendered screenshot inspected; original source_fields are expanded by default so
run/condition/replicate declarations are visible. Prototype proof/codeecb37bea
and source60b6458 retained,followup validates the visible-field improvement. The proof uses
hardlinks to immutable local data copies; never write through those input paths.

Download3386467,PRIDE3533750,quant3614034/kallisto512367 confirmed alive atturn
start. Latest325/1011RNAverified,324snapshot; native jobs continue,24CPUmax.
Continue acquiring/importing remaining formats, original row-design associations
and qualified accession/gene/scale compatibility. Full goal is not complete.

**Current row-wise admission (2026-10-30):** Previous goal turn was progress
(log2 quantities and expanded source-column links), not a wait. Dataset goal
remains active/incomplete; historical inference work remains paused. Source
f9df520 and prior followups are pushed to nightly; public0.0.6.0 unchanged.

New native `process_protein_rowwise_quantities.py`48635f1f admits only the explicit
ten-field protein-accession PG.Quantity/R.Condition/R.FileName/R.Replicate layout
into `prepared_deposited_values/protein_rowwise_quantities`. Peptide/site/fragment
extensions are excluded from this adapter, not mislabelled as protein reports.
PXD084269 source65535rows /1quantitycolumn /29564missing /5filename literals.
Original row ordinates, labels, quality declarations and all run/condition/replicate
strings retained. No pivot, extra structural cells, filtering, normalization,
imputation, duplicate collapse, averaging or assumed log/linear scale. One source
quantity column does not imply one biological sample. Group identifier modes and
source-plus group flag semantics do not interpret this accession assay.

Twelve synthetic corruptions refused. Every source cell/token/binary64/missing
cell/row/annotation and repeated-key count independently reconciled. All138 prior
quantity packs/manifests/outputs unchanged. Full139quantity catalogue now1384744
originalcells /71153missing /1097source columns (132protein/7Workbench tables).
Actual offlineQt source filtering, original row-metadata tooltips, column context,
CSV/treatment exports and unfiltered QC verified under
`/media/carruthers/mnt1/parascape-rowwise-validation-20261030`.
Native `audit_protein_rowwise_repeated_keys.py` fda9f15b also runs every5minutes
with1GiB/24CPUaffinity/BLAS2. It independently reparses and checksum-verifies the
original source and annotations, records missing/observed key categories and
field differences without changing rows, and archives code/operations/immutable
reports under `processing_audit/rowwise_repeated_key_qc`. Collection4cee5856… is
current after synthetic checks for identical/conflicting/mixed/missing cases and
case-sensitive keys. Expanded source-column collection now132packs/857columns;
all131 previous report files byte-identical. Its new scalar PG.Quantity column
correctly stays filename-role-unreviewed: actual filenames are row-specific.

New5minute native timer `parascape-protein-rowwise-values` enabled24CPU/2GiB,
BLAS2; existing download/quantification jobs remain alive. Salmon480572 completed and
quant supervisor3614034 advanced to active kallisto512367 on the same run;
325/1011entries stillverified. PRIDE3533750 has two open partial downloads;
confirmed transfer growth is distinct from completed/hash-verified file counts. Future portable release
must include new accession/log2/row-wise namespaces and updated source links;
don't modify public frozen6asset.

Repeated-key audit finds24820 exact accession+filename keys with2rows each. None
has every original field identical;24041mix observed+missing quantities,779all
missing,0with two observed quantities. Names/descriptions/qvalues can differ.
Do not deduplicate or pivot on this key, infer biological independence or impute
from another row. Exact65535data-row count alone is not proof of truncation.
Audit JSON retained in new namespace/runtime-validation. Row_layout metadata
records source rows, filename literal count and repeated-key multiplicities.
Continue explicit row-specific source-design associations and accession/gene/scale
compatibility work; full all-dataset/fair-comparison goal remains incomplete.

**Current dataset continuation (2026-10-30, supersedes older pending-release notes):**
Previous goal turn made verified progress: accession importer enabled and verified;
0.0.6.0 publication completed. Current dataset goal is active and incomplete.
Historical inference goals64/65/66 remain paused. Official6 publisher state at
`~/.local/state/parascape-release-0.0.6.0/status.json` is complete at38ad27e;
main publish run41144027558 succeeded. Live checkout/cache synchronization also
completed. Deferred source-flag followupb079984 finished normally; source and
documentation11fc841 are on nightly. PyPI6 remains frozen324RNA/136quantities;
source-flag view feature and new local admissions are outside that wheel/asset.

New native `process_protein_log2_quantities.py` preserves original accession rows
and explicitly labelled log2 filename columns in independent namespace
`prepared_deposited_values/protein_log2_quantities`. First table PXD084269:
7576rows ×5columns =37880cells /1909missing. Processor6c6a3f4c retains all
candidate/completeness/CV/average-log2 declarations, repeated accessions and exact
filenames. No exponentiation, imputation, normalization, averaging or filtering.
Ten synthetic corruptions refused; actual full-cell/token/binary64 and offline
Qt/CSV/context/treatment proof at
`/media/carruthers/mnt1/parascape-log2-validation-20261030`.
All137 prior packs/manifests/outputs unchanged. Native5minute timer enabled,
24CPUaffinity/2GiB, BLAS2. Unsupported layouts retain explicit review gaps.

The previous-turn accession abundance source PXD041907 has2496rows /6groups /
14976cells /1550missing. Source labels F1..F6 link only exact same-file presence
header candidates; no sample/replicate identity inferred. Its processor c53bef53,
8-corruption synthetic checks and actual offlineQt/CSV proof are archived under
`processing_audit/protein_accession_values`; local proof
`/media/carruthers/mnt1/parascape-accession-validation-20261030`.
Both accession assays remain distinct from protein-group quantities, do not enter
the group-qualified identifier index, and do not inherit group flag meanings.

Native `export_protein_column_design_links.py` d14bba33 now discovers all active
protein-group/accession namespaces with catalogue and manifest pins/race refusal.
Exact PG.Log2Quantity filename extraction added without changing old PG.Quantity
rules. Collection ed06f241… covers131protein packs /856columns:16exact primary
filematches,12unresolved,828unreviewed roles,10submitted rowoccurrences,
357experiment-label candidatecolumns. Every65 old report file is byte-identical.
New log2 source supplies5exact primary files+5submitted rowoccurrences; biological
sample/control/replicate independence is not inferred. Full131reader proof plus
actual offlineQt source tree/CSV are at
`/media/carruthers/mnt1/parascape-all-column-links-validation-20261030`.
Original sources are independently reparsed; ambiguity/duplicate-row/extension
conversion refusal checked. Existing5minute source-link timer stays enabled.

Live quantity catalogue138tables (131protein/7Workbench), 1319209originalcells /
41589missing /1096source columns. Counts repeat quantity families/columns and are
not independent samples or unique proteins. New source tables and expanded links
await future portable release; don't mutate frozen6. Downloads3386467,
PRIDE3533750, quant3614034 and active salmon480572 confirmed alive this turn.
Latest325/1011verified RNA entries,324integrated; nextbatch349. Resource guard
fast mode reserves12-15/28-31,24CPUmax; internal mapper may use8 effectively.
Continue acquisition, remaining original quantity formats and explicit sample/
identifier/scale compatibility work; do not mark full goal complete.

**Original source flag views (2026-10-30, after frozen6):** Independent follow-up
checkout `/media/carruthers/mnt1/parascape-source-flags-followup-checkout` starts
at38ad27e. New source_flags.apply annotates literal source fields and optionally
excludes exact+ in Reverse/Contaminant/Potentialcontaminant/Onlyidentifiedbysite
only for declared protein-group quantities. Default remains unfiltered. Missing
or unfamiliar values and other assays remain unreviewed/retained. Identity is
manifestSHA+source_row; repeated gene/group labels and source files cannot share
flags accidentally. Inconsistent row metadata is refused. No biological score,
normalization, imputation, pooling, sample/control/replicate inference or source
file mutation. This is an explicit view operation, not biological validation.

Other measurements adds visible flag declarations and an opt-in checkbox.
Changing its rule clears stale results and treatment state. CSV exports carry
uniquely named treatment JSON sidecars with CSV digest, source hashes, exact query,
identifier evidence, policy/code/version and excluded source rows; earlier records
are preserved. Separate treatment export works for all-excluded/no-match results.
All-excluded matches are distinguished from missing expression. Whole-matrix QC
remains unfiltered and says so explicitly.

Offline independent vectorized flag audit covers all136 tables /177781 original
source-row occurrences:12021 marked and22571 with no declared fields (counts
repeat quantity families, not independent samples or unique genes). Retained
numeric bits on each first original column and every original output-file hash
are unchanged. RealQt/CSV sidecar, all-excluded log and unfiltered whole-source QC
passed with NAS/network forbidden. Proof is retained under
`/media/carruthers/mnt1/parascape-source-flags-validation-20261030-v2`.
Keep frozen6checkout/source38ad27e unchanged while exactCI is pending. New follow-up
must wait for official6 publication and live-checkout synchronization before
nightly push, so it cannot cancel that CI. Full dataset goal remains active and
incomplete; historical inference goal paused unchanged.

**Expanded0.0.6.0 frozen preparation (2026-10-30):** Isolated source checkout
`/media/carruthers/mnt1/parascape-release-0.0.6.0-checkout` starts atccce694;
version bumped with scripts/release.py. Task-owned stage
`/media/carruthers/mnt1/parascape-release-data-0.0.6.0` freezes5690 files:
324 RNA entries /147 original matrices /54 platform indices /136 other quantity
tables (129protein/7Workbench) /601605 identifier evidence occurrences /52 FASTA
occurrences in36 content indices, source QC/designs and full processing history.
Archive985049339 bytes, SHA
`4f6c1a686a253200b025eefbeb4b220c9bcb84f8fca437d640f3809d43b6a17f`.
Strict archive installation verified every file. All72 previously published
quantity manifests/outputs and1747 pinned non-RNA files remain byte-identical;
all40891575 earlier RNA cells across120 matrices and their axes remain unchanged
while RNA matrices add verified entries. Frozen identifier and reference catalogue
contexts match. These facts are independently replayed in portable-validation.json.

Native prep/build/strict validation units006 completed successfully. Expanded
installed smoke now checks all136 quantities,129 identifier packs, source filtering
and stale-export refusal, with NAS/network forbidden. It awaits clean source
freeze/build/install-wheel validation and exact-commit CI. Initial629 focused checks
passed but one old packaging assertion wrongly demanded data inside the source
package; related path check had the same obsolete assumption. These are corrected
to require present installer icons/descriptor and writable, versioned external
cache resolution, consistent with the requested GitHub storage. Scientific data
and runtime path logic are unchanged. Retain failed/recovery test logs. The
corrected combined release/reader/UI/path/download/documentation check set now
passes662 tests (one existing runpy warning). Native installed-wheel verification
and exact-commit publication gates remain required before PyPI6 can be claimed.

**0.0.5.0 officially published (2026-10-30):** corrected929114d source CI and
manual publication run41137156827 succeeded. Official PyPI JSON independently
verifies wheel1840263 bytes, SHA
`c9f862a2d75692cb557e6b557ed9285820c3fe2a1e13e407e4358a661e4fc707`;
external872109077-byte data SHA44f3db3f…65522 unchanged. Published source tag
v0.0.5.0 remains90c2c40; workflow-only retry929114d leaves package/data unchanged.
Native publisher verified official wheel/descriptor/public gzip and synchronized
live checkout/cache. Deferred checked follow-upccce694 then pushed normally to
nightly and updated the live source checkout. Its completed service was disabled.
Native6 publication dependency must record previous version5 commit929114d as in
its completed state (published immutable tag remains90c2c40). Current live source
already has the136-table searchable/qualified browser; PyPI5 snapshot remains307
RNA entries and72 quantity tables. Do not claim PyPI6 before all publication gates.
Full all-dataset goal remains active/incomplete; historical inference goal paused.

**Searchable original-measurement browser (2026-10-30):** The independent
follow-up now filters Other measurements by case-insensitive text in study,
source analysis/file, assay and units. Feature queries stay case-sensitive.
Visible/total source counts remain on screen; refresh preserves the filter and a
still-visible selected source. Changing filter/source/identifier class/feature
query clears prior result/export state. Empty matches disable source actions;
feature edits cannot silently export the preceding query. Source choices and
numeric data remain separate; filtering never pools or changes quantities.
Four integration tests verify source/unit provenance, stale-export refusal,
selection retention, empty states and feature-edit invalidation. The all136-table
actualQt/CSV offline proof is retained under
`/media/carruthers/mnt1/parascape-source-filter-validation-20261030`.
The existing whole-window tooltip audit is included in final validation. Source
commit/push stays deferred behind exact929114d publication and checkout sync;
full all-dataset goal remains active and incomplete.

**Source-labelled scalar protein quantities (2026-10-30):** Native producer
2efa199f admits exact aggregate/channel/ratio field families in an independent
`prepared_deposited_values/protein_reported_quantities` namespace. The first
62 tables from29 acquired files /17 projects preserve85899 source cells and8662
missing cells:Intensity29, iBAQ17, IntensityL4, IntensityH4, RatioH/L4,
RatioH/Lnormalized4. Exact LFQ-intensity aggregate and IntensityM labels are also
supported when explicitly present. Scalar header slashes are escaped only in
internal pack names; original labels remain literal. No normalization, correction,
ratio conversion, imputation, filtering, group/channel pooling or inferred
sample/control/replicate identity. Full source annotations and flags remain.
Quantity-role metadata distinguishes source aggregate/unspecified scope, labelled
channels and reported ratios; sample assignment/orientation remain unresolved.

Every source cell/token/bit pattern/row/annotation and scalar role is independently
reconciled. Changed roles and source flags with rewritten checksums are refused.
All74 prior quantity manifests/outputs remain byte-identical and cached new outputs
stable. Native5minute CPU24/2GiB timer remains active. Offline all62-reader/QC checks
and actualQt/CSV checks for all6 currently present families are at
`/media/carruthers/mnt1/parascape-reported-quantities-validation-20261030`.

Live catalogue now136 quantity tables (129protein/7Workbench),1266353 cells across
19 protein projects. Automatic identifier discovery includes all129 protein packs:
175988 source-row occurrences and601605 evidence occurrences (468358 original
components,67055 submitted labels,19674 explicit reference gene declarations,
46518 GN-name candidates). Existing67 association packs remain byte-identical;
counts repeat quantity families and ambiguous components, not independent samples
or unique genes. Expanded offline129-pack API/Qt/CSV proof is retained under
`/media/carruthers/mnt1/parascape-expanded-protein-identifier-validation-20261030`.
These additions await a later data release; frozen5 is unchanged. Full dataset
goal remains active and incomplete; biological/gene/scale/sample comparability
requires further work, and historical inference work stays paused.

**All-namespace protein lookup and live RNA324 (2026-10-30):** Native producer
47f2e9d1 discovers every reconciled original protein quantity pack across root and
one-level namespaces, excluding other assays. It pins exact source catalogue
bytes/hashes, selected inventory and each source manifest; catalogue races and
unsafe/symlink namespace escapes are refused. The collection now covers67 tables
/90089 original source-row occurrences /239365 evidence occurrences:197471 original
components,21025 submitted labels,13116 explicit reference gene declarations and
7753 GN-name candidates. The reporter families add10104 component occurrences;
no new gene/name candidates. All65 prior association outputs and index entries
remain byte-identical. Reader/UI code is unchanged. Full67-pack offline API and
actual Qt/CSV proof (including component queries in both reporter families) is at
`/media/carruthers/mnt1/parascape-all-protein-identifier-validation-20261030`.

The immutable324-entry RNA snapshot finished all three build/extension/enrichment
steps and standard independent validation; its live NAS pointer now selects
`app_snapshots/20261030T214036Z-324/rna`.24 studies,120 matrices,1724 quantified
sample conventions and242 output checksums verified;1011 entries remain planned.
Original source samples, quantity conventions and processing/design metadata stay
separate. Full preexisting307-entry value/axis preservation is independently
replayed and retained as `prior-rna-preservation.json` in the new build folder.
All40891575 preexisting RNA cells across120 matrices are bit-identical after
307→324 expansion, with old feature/sample axes retained. Native snapshot worker
5de3c79e now enforces this preservation before every future live-pointer update
and archives its own exact source with each build. Synthetic tests refuse changed
values, imputation, signed-zero changes, missing samples/features/matrices and
duplicate identities while admitting new samples/features and safe row reordering.
The real new preservation gate replays the same307→324 snapshot independently;
proofs/code are retained under `processing_audit/rna_snapshot_preservation`.
This is a live NAS update; frozen0.0.5.0 still contains307 RNA entries and72 other
measurement packs. Release929114d Tests41132730010 remains in progress,
Documentation41132730719 succeeded; native publisher and deferred follow-up remain
active. Full acquisition/fair-comparison goal stays active and incomplete.

**Original reporter-channel admission (2026-10-30):** The separate native
`parascape-protein-reporter-values.timer` now processes exact numbered corrected
and uncorrected Reporter intensity fields from acquired protein-group tables.
First source PXD039148:4940 original group rows ×6 channels per family;2 packs /
59280 cells /12 family-channel columns /0 missing. Original count fields,
aggregate intensity, quality flags and group identities stay intact. Channel
ordinals are explicit source metadata, never inferred samples/labels/replicates.
No source correction reapplied, filtering, imputation, channel pooling or scale
harmonization. Existing Other measurements browser discovers the independent
`prepared_deposited_values/protein_reporter_channels` namespace automatically.
Live other-measurement catalogue now74 packs (67protein /7Workbench),1180454 cells;
all72 prior packs unchanged. These2 new tables are outside frozen0.0.5.0.

Producer7ad58236 independently reconciles every literal token, observed binary64
value, original row/annotation and channel/family identity via pandas versus CSV.
Synthetic changed channel ordinals and changed source flags with rewritten
checksums are refused. A first prototype shared the earlier worker's status path;
it produced valid isolated output but was corrected to its own lock/status state.
Earlier prototype packs/code remain retained outside the active catalogue;
previous protein worker status was restored from its own validated cached cycle.
Native5minute CPU24/2GiB timer and producer/selfcheck/proofs are retained. Offline
actualQt/CSV/whole-columnQC proof under
`/media/carruthers/mnt1/parascape-reporter-validation-20261030` confirms both
families and all12 original channel metadata with NAS/network blocked. Qualified identifier
lookup now covers all67 protein packs; reporter tables add original components,
with no new submitted gene/reference candidates. Source-design links still cover
the earlier65 packs; reporter sample/label assignments remain unreviewed. Future data
release preparation must include this independent namespace, not just protein_groups.
Full dataset acquisition, sample-design reconciliation and fair comparison remain
active and incomplete; historical inference goal stays paused.

**Publication environment retry (2026-10-30):** Release5 passed corrected source
CI and published its GitHub data/tag at90c2c40, but PyPI publication run41084621531
failed its installed-app CSV smoke because PARASCAPE_SMOKE_OUTPUT was absent.
Workflow-only commit929114d supplies `/tmp/parascape-release-smoke`; all package
sources/scripts/version/data descriptor remain byte-identical to the published
90c2c40 tag, which stays unchanged.31 release/packaging checks pass; actual installed
full72-pack/Qt/CSV smoke passes with the corrected environment and NAS/network off.
Native `parascape-publish-005-workflow-retry.service` pins this one-file change,
requires exact929114d nightly Tests+Documentation, promotes main without force,
explicitly requests the authorized same-version publication retry and verifies
its workflow event plus official PyPI wheel/digest/descriptor/public data URL.
Its plan records both the retry commit and immutable published source tag commit.
Old publisher is disabled; failed logs/proofs/plans remain archived. Follow-up
source push now waits for this corrected publication and live-checkout sync.

The native RNA queue has324 verified study/run entries; the current live snapshot
still307. A separate `parascape-snapshot-catchup-324.service` invokes the normal
immutable builder with `--force` to bypass only its25-entry batching threshold;
all receipt, sample identity, matrix, scale/convention and checksum checks remain.
Only a validated build may replace the live NAS pointer. Published snapshots stay
frozen. Normal timers and ongoing acquisition/processing continue independently.

**Qualified protein identifier lookup (2026-10-30, after frozen 0.0.5.0):**
The independent follow-up checkout at
`/media/carruthers/mnt1/parascape-identifier-followup-checkout` adds exact original
protein-component, submitted-label, reference gene/gene_id and qualified UniProt
GN-name candidate lookup to Compare → Other measurements. Whole original group
values and every ambiguous evidence occurrence remain intact. Duplicate evidence
never duplicates measurement records; source row identity prevents propagation to
another row with the same group label. CSV retains all supporting declarations.
No gene-version conversion, group collapse, inferred database use or unique-gene
abundance. Reference candidates exist in only five packs / two projects; submitted
labels in twelve packs / five projects. These are occurrence counts, not genes or
independent biological replicates.

Native five-minute `parascape-protein-identifier-evidence.timer` reconciles65 packs
/80209 source feature-row occurrences /229261 evidence records:187367 original
components,21025 submitted labels,13116 explicit reference gene declarations and
7753 GN-name candidates. Source policy and processor revisions, every immutable
Parquet association, collection receipts and operations are retained under NAS
`processing_audit/protein_identifier_evidence`. A revised processor uses4096-row
reference batches and retains only exact measured components within their declared
project. This fixes a real2GiB service OOM without raising its memory ceiling;
independent pandas eligibility and association reconstruction preserve all records.
Earlier verified collections and the OOM event remain retained. Original quantity
packs are unchanged. Runtime association data follows the frozen5snapshot and is
not included in that release.

Offline original-group/API/actualQt/CSV proof is retained under
`/media/carruthers/mnt1/parascape-protein-identifier-validation-20261030`;
572 focused reader/UI/QC/documentation/release-builder checks pass. Before nightly
push, wait for both official release5 completion at exact929114d AND its native
live-checkout synchronization receipt. Do not cancel frozen CI with a follow-up
push. Full acquisition and fair-comparison goal remains active and incomplete;
historical inference goal remains paused.

Release5 originally failed one GitHub UI audit (5987 tests passed): five reference
browser controls lacked tooltips. The reviewed correction90c2c40 passes608 focused
checks /1skip and is being independently rebuilt and checked before the same strict
publication gates resume. Old failed source and proof records remain archived;
the data archive is unchanged. Recovery unit:parascape-release-005-tooltip-recovery.


**0.0.5.0 frozen release preparation (2026-10-29):** clean detached checkout
`/media/carruthers/mnt1/parascape-release-0.0.5.0-checkout` starts from90e71e7;
version bumped with scripts/release.py. Task-owned data stage
`/media/carruthers/mnt1/parascape-release-data-0.0.5.0` contains4957 files,307 RNA
entries,147 original-value matrices,54 platform indexes,72 other measurements
(65 protein /7 Workbench), source QC, submitted designs, column associations and
46 FASTA source occurrences /30 content indexes. Archive872109077 bytes, SHA
`44f3db3f2ceb57620d3fd09b5c2f7486298649bf6f2f8c64c7bb818fec251522`.
All1380 prior non-audit files, including120 RNA matrix files, remain byte-identical
to published0.0.4.0. Strict archive install verifies all4957 files. Acquisition
counters are explicitly native processing observations, not bundled raw-file counts.
Only active measurement/index collections are bundled; historical source evidence
and exact processor revisions remain retained. No normalization/imputation or
identifier/probe/protein-group collapse. Native prep/build/validation sources and
operation records retained; downloads/quantification stay independent.

624 focused reader/UI/QC/download/release/path/documentation checks pass before
freezing the source. The installed-wheel smoke now exercises all72 measurement
packs, all reference indexes,209 PRIDE projects/52 designs and actual Compare tabs
and CSV exports with NAS/network blocked. Root nightly stays at90e71e7 during the
frozen release; after successful publication it can fast-forward to this commit
only if tracked files remain clean. Do not push followups that cancel this exact CI.

Release publication is still pending installed-wheel validation, exact-commit
Tests+Documentation and final native publication gates. Do not claim PyPI0.0.5.0
until the official version API/wheel/descriptor/public asset checks pass. Full
all-dataset goal remains active and incomplete after this bounded release.

**0.0.4.0 publication completed (2026-10-29):** frozen bb8ebca passed exact nightly
Tests40967543108 and Documentation, then main/publication gates. Official PyPI
0.0.4.0 JSON verified: wheel1756202 bytes, SHA
`a9234a6e40b1519aebc42544371be0ccc3d4cf7285410de346b477ab986eaa40`;
GitHub tag v0.0.4.0 data SHA remains9fa0662f…0dc9c09. Publisher status is complete.
The reviewed follow-up push automatically advanced nightly through35fc377 and
completed. Earlier frozen-release-wait notes are historical: coherent new checked
changes can now be pushed normally to nightly. Public0.0.4.0 includes307 RNA
entries,147 original-value matrices and54 platform indexes; later Workbench,
protein/metadata/source-link/reference-header additions await the next release.

**Original reference identifiers (2026-10-29):** Compare → Reference identifiers
now supports exact original primary-ID and explicit source gene-field lookup with
CSV provenance, preserving repeated rows, case, version suffixes and header text.
Native `parascape-reference-identifiers.timer` runs every5min atCPU24/1GiB.
Current54 planned FASTAs:46 acquired occurrences,8 pending;30 distinct content
indexes /297621 header records /75597 headers with explicit gene declarations.
Fourteen source occurrences contain such declarations. These are not unique
current genes or proof that a reference database was used for a measurement.
GFF/GTF files are outside this FASTA header workflow.

NAS `processing_audit/reference_identifiers` retains current/historical catalogues,
source signature/SHA verification, exact processor revisions, immutable Parquet
header packs, manifests, independent receipts and operation logs. Processor
fa5c6303 preserves all original headers and strict gene/gene_id/transcript/
transcript_id declarations; binary-line lexical validation independently replays
every text/regex-produced cell. No sequence transformation, symbol conversion,
canonical gene-version mapping or biological inference. Synthetic changed fields
with updated output checksums are refused; BOM/CRLF/trailing whitespace and
repeated IDs preserved. Cached inputs use recorded source signatures, with full
SHA verification on new/changed sources; output hashes remain checked.

Offline proof `/media/carruthers/mnt1/parascape-reference-identifiers-validation-20261029`
checks all46 source occurrences /30 content indexes, exact primary/gene queries
and actual Qt CSV preservation with NAS/network blocked.550 focused reader/UI/
docs/Compare/release-builder checks pass. These new data/browser additions are
unreleased. Next larger steps include linking source reference declarations to
measurement identifiers without collapsing groups or assuming cross-version
identity, then freezing the expanded data for a subsequent release. Full dataset
goal remains active/incomplete; historical inference goal remains paused.

**Submitted experiment-label candidates (2026-10-29):** Source columns and
submitted factors now also shows exact `Experiment` label candidates for original
Intensity/LFQ/iBAQ column suffixes within the same declared project. All matching
tables and fraction rows remain candidates, not confirmed quantification input
linkage, control pairing or independent samples. Existing exact filename/SDRF
links remain separate and unchanged. Original source columns/values stay intact.

Native column-design exporter revisioned160699 now retains20 explicit Name /
Fraction /Experiment tables (328rows) from the acquired schema inventory.
All245 rows from16 rectangular tables are independently reparsed with pandas,
including complete row coverage. Four irregular tables /83rows retain every
original cell/width as review-only JSON and are excluded from candidate links;
no trimming, padding, truncation or record deletion. Eligible sources produce357
column/row candidate occurrences across6 projects. These repeat across quantity
families and are not independent sample counts. Current filename counts remain11
exact /12 unresolved /748 unreviewed roles /5 submitted SDRF row occurrences.

Native five-minute timer stays CPU24/512MiB; immutable receipts bind all20 source
hashes and schema inventory. Source code and operation records are retained in
NAS `processing_audit/protein_column_design_links`. Synthetic checks retain
leading zeros, exact case, multiple matching templates and fractions, and reject
silently selected templates or dropped source rows. Earlier filename ambiguity/
repeated-row/factor refusal checks still pass. Offline proof
`/media/carruthers/mnt1/parascape-column-design-validation-20261029-v3` checks
all72 packs /771 protein contexts and actual Qt experiment-candidate rows /CSV
exports with all original fields and unreviewed linkage scope.528 relevant checks
pass, then11 reader/dialog checks including candidate-versus-filename separation.

These changes remain local/unreleased behind the frozen0.0.4.0 exact CI gate.
Full dataset goal active/incomplete; historical inference goal paused. Remaining
work includes source input-table linkage, irregular designs, other quantitative
formats, gene/scale harmonization, raw acquisition and fair cross-study comparison.

**24-thread handover verified complete (2026-10-29):** native completion-boundary
handover exited successfully. Quantifier parent3614034 has actual environment
`PARASCAPE_THREADS=24`; its new fasterq-dump child3614170 uses `-e24` for
SRR23248136. The earlier eight-thread run was preserved through completion. The
resource guard currently uses protected CPUs0-7 under observed I/O pressure;
quiet ceiling remains CPUs0-11,16-27. Audit proof is NAS
`processing_audit/resource-policy-24-thread-handover-complete.json`. Earlier notes
about the pending handover describe the prior state, not a further required restart.

**Exact quantity-column source links (2026-10-29):** Other measurements now has
**Source columns and submitted factors…**, available without a feature query.
Original sample/column records remain unchanged; separately pinned filename links
and submitted rows appear in the dialog and value CSV `sample_metadata` under
`source_design_metadata`. Numbered fields preserve duplicate SDRF headings and
all repeated row/channel/fraction occurrences. Selection/refresh clears old dialogs.

Native `parascape-protein-column-design-links.timer` refreshes every5min at CPU24 /
512MiB, independently reparsing original SDRF cells and reconciling every exact
primary filename candidate. Immutable collection, receipt, source bindings, exact
exporter code and operation records live in NAS
`processing_audit/protein_column_design_links`. All65 tables /771 columns covered:
11 exact unique primary filename matches,12 unresolved explicit filenames,748
columns with unreviewed filename roles,5 matched submitted-row occurrences.
Filename tokens are extracted only from explicit numbered PG.Quantity headers;
no case/basename/extension conversion or experiment-label guessing. No biological
independence, paired controls or downloaded raw payloads inferred. Missing/changed
measurement manifests do not acquire stale associations.

Synthetic checks preserve ambiguous primary candidates and repeated/duplicate SDRF
fields; dropped repeated rows, silently resolved ambiguity and changed source
factors are refused. All72 packs and all771 protein-column contexts, actual Qt
no-query dialog and CSV exports validated with NAS/network reads blocked; proof
`/media/carruthers/mnt1/parascape-column-design-validation-20261029-v2`.
555 focused reader/UI/QC/docstring/docs/Compare/release-builder checks pass; final
10 reader/dialog checks pass after readable association captions/provenance grouping.
All72 measurement manifests and288 original runtime outputs remain byte-unchanged;
cached association export preserves the immutable collection pointer and receipt.
New metadata is an unreleased local follow-up. Current dataset goal remains active
and incomplete; historical inference goal remains paused. Next broader gaps include
other quantitative formats, unsupported filename roles, gene/scale harmonization,
raw acquisition and validated cross-study comparison.

**PG.Quantity original pivots (2026-10-29):** native protein adapter now also
supports one explicit `PG.ProteinGroups` axis with strict `[number] file.PG.Quantity`
columns, retaining all non-quantity annotations and quality fields. Two additional
sources add135116 cells /5075 missing /23 columns. Current local totals are65
protein tables from28 source files /16 projects,1063706 cells /7529 missing /771
source-column occurrences. Six other layouts remain explicit review gaps; long
and derived-statistic tables are not forced into abundance matrices.

Processor revision4776d479 is archived by SHA on NAS; five-minute timer remains
active atCPU24/2GiB. Independent stdlib-versus-pandas cell/annotation/QC checks
pass; separate synthetic pivot fixture checks repeated groups, original precision,
signed zero, missingness, strict quantity-column selection and refusal of changed
or ambiguous feature axes. All63 earlier matrices retain identical values/tokens/
column metadata (189 files), with identical row annotations and flags. Offline
all72 deposited-pack API/QC and actual Qt PG.Quantity CSV exports pass; proof
`/media/carruthers/mnt1/parascape-protein-values-validation-20261029-v2`.
Cached replay preserves all65 receipts and the catalogue; coverage16 projects /
65 tables /0 inferred biological samples confirmed. Source processor, preservation
proofs, tests and operation logs remain in NAS `prepared_deposited_values/protein_groups`.

Next sample-design task: resolve explicitly filename-bearing quantity headers
against exact primary project file identities and validated SDRF occurrences,
retaining all channel/fraction/repeated-row associations rather than guessing
independent samples or paired controls. Gene/scale reconciliation and other
quantitative formats remain unfinished. Frozen0.0.4.0 exact Tests40967543108 has
started on GitHub (job123776343818); the full regression step is running after
successful packaged-layout and real-OpenGL checks. Publisher remains live; no
follow-up push may cancel that frozen CI. Dataset goal active/incomplete, separate
historical inference goal paused.

**PRIDE submitted-design browsing (2026-10-29):** Dataset coverage now exposes
primary metadata for all209 PRIDE projects and52 original SDRF designs /1392rows
across49 projects. Numbered fields preserve duplicate headings; repeated labels,
empty cells, literal missingness words, fractions and channels stay source records,
not inferred biological replicates or paired controls. All1387 exact filename
occurrences and5 unresolved references remain visible. Community49/archive3
metadata origins and non-independent biological truth status are disclosed.
Projects without represented SDRF designs retain metadata plus an explicit gap.

Native `parascape-pride-sample-design-export.timer` runs every5min, CPU24/512MiB,
reparsing every original SDRF cell with pandas independently of the original
stdlib adapter and checking all primary filename candidates. Immutable portable
collections, source/code bindings, versions and operation receipts live at NAS
`processing_audit/pride_sample_designs`. No raw payload or measurement-column
pairing is inferred. Cached export preserves collection pointer/receipt; synthetic
changed factors and silently resolved ambiguous filename links are refused.
Actual Qt selections for all209 projects and all52 designs checked with NAS/network
blocked; proof `/media/carruthers/mnt1/parascape-pride-design-validation-20261029-v2`.
526 relevant checks pass, followed by11 final metadata/UI checks after adding full
field tooltips. The frozen0.0.4.0 release and its exact CI gates remain unchanged;
these additions await a later release. Dataset goal remains active/incomplete;
historical inference goal remains paused. Remaining broad gaps include quantitative
formats, source-to-measurement sample joins, raw acquisition, gene identity/scale
reconciliation and independent cross-study comparability.

**Protein-group original values (2026-10-29):** local Other measurements now lists
63 independently reconciled quantity tables from 26 source files / 14 projects,
alongside the existing seven Workbench tables. Intensity26/LFQ20/iBAQ17 remain
separate; 928,590 original cells, 2,454 missing cells and 748 submitted quantity
column occurrences checked. All original group identities, duplicate rows, literal
tokens, source annotations and quality flags retained without filtering, row
collapse, imputation, normalization or gene conversion. Source `NaN`/empty tokens
become missing numeric cells only; zero semantics require review. Exact flags
appear in tooltips and CSV `source_row_metadata`. Quantification columns are not
inferred biological samples or independent replicates.

Native `parascape-protein-group-values.timer` runs the source-hashed processor and
independent pandas-versus-stdlib full-cell/annotation/QC validator every five
minutes, CPU0-11,16-27 / MemoryMax2G / BLAS2. Data, manifests, receipts, original
source bindings, code and operation logs live in NAS
`prepared_deposited_values/protein_groups`; namespace catalogues coexist with the
Workbench root catalogue. Three Protein IDs sources lack unambiguous supported
quantities and remain explicit review gaps (PXD002380, PXD042937, PXD039148).
The earlier 27-table candidate audit included one ambiguous-only source; two
additional Protein IDs tables have no eligible family headers. Other schemas,
source confidence/filter decisions, enrichment context, sample/fraction/channel
roles and gene/scale reconciliation remain unfinished. No whole-cell RNA/protein
concordance or biological quality prediction is admitted.

Portable proof `/media/carruthers/mnt1/parascape-protein-values-validation-20261029`
queries all70 packs, checks full QC and actual Qt CSV export with NAS/network reads
blocked. Independent validation refuses changed source quality flags even with
updated output checksums. 544 reader/UI/QC/docstring/docs/Compare checks and two
portable release-builder checks pass. These additions remain unreleased; frozen
0.0.4.0 still waits exact bb8ebca CI. Native downloads and processing continue;
full all-dataset fair-comparison goal remains active and incomplete, historical
inference goal remains paused.

**Latest resource authorization (2026-10-29):** user permits up to 24 of 32 logical
threads, superseding the earlier eight-thread ceiling. Owned worker units now use
CPUs `0-11,16-27`, reserving four full physical cores/eight logical threads.
Quantifier environment is `PARASCAPE_THREADS=24`; the live older worker keeps its
loaded eight-thread configuration until native
`parascape-thread-limit-handover.service` observes a verified quiescent completed-run
boundary. It never restarts an active quantifier child or unfinished source run.
The resource guard uses a 24-CPU quiet ceiling, eight protected CPUs and 30 seconds
of sustained quiet before restoring capacity. Cache/storage caps remain unchanged.
Audit policy, previous unit/runtime sources and handover state are retained in
`~/.local/state/parascape-thread-limit` and NAS processing audit. Future quantifier
completion records include configured threads and the source hash captured at
process startup; the existing process and tutorials remain untouched.
Queue/cache regressions: 26 checks and five resource-guard checks pass; the handover gate rejects unfinished
runs or leftover scratch and requires both complete tool outputs. The first full
protein-group source audit inspected 27 acquired tables (~59.7MB), preserving
headers, rows, original group identities, quality flags and numeric-token review
counts in `published_pride/protein_group_candidate_audit.json`. The preliminary audit marked literal `NaN` tokens as nonfinite; the adapter now
handles only exact empty/NaN tokens as missing and preserves their original text.
Header-only selection must not treat `iBAQ peptides` or aggregate/channel fields
as ordinary sample abundance. See the newer protein-group import above for the subsequent original-value
workflow; source quantities, zero semantics, decoys/contaminants and assay context
remain explicit review requirements before RNA/protein comparisons.

**Workbench original measurements (2026-10-29):** local **Other measurements** now
reads seven independently reconciled original tables in
NAS `prepared_deposited_values`. All 57,468 source cells and literal tokens checked
via stdlib JSON generation versus independent pandas/ujson parsing and exact
binary64 comparisons; 21,939 missing cells remain missing. Source rows, units,
annotations and submitted factor associations are retained without normalization,
imputation, feature collapse or gene mapping. Full QC covers 240 sample/analysis
records; reused samples across analyses are not independent replicates. One
analysis has metadata only. The native Workbench pipeline now sequences source
acquisition, metadata reconciliation, preparation, independent cell validation/QC
and coverage refresh without model calls. Runtime catalogue exposes only completed
reconciled packs. Live launcher sets `PARASCAPE_DEPOSITED_VALUES`; the external
asset builder supports `deposited_values`. 564 relevant checks pass. All seven
packs queried and QC checked in a local portable copy with NAS/network reads
blocked; actual Qt value/QC CSV exports verified. A changed numeric cell remains
rejected even after updating file/manifest checksums. New additions await a later
portable release; frozen 0.0.4.0 still requires exact `bb8ebca` CI/publication.
Raw payload acquisition, cross-study harmonization, protein data admission,
current gene identity reconciliation and the full dataset goal remain unfinished.

**Workbench source gap (2026-10-29):** all seven catalogue studies now have acquired
primary responses: 90 API responses, eight analyses and 276 submitted factor rows.
NAS `metabolomics_workbench` retains exact source bytes, HTTP/hash receipts, failed
attempts, superseded reports and versioned collector sources. Persistent
`parascape-workbench-acquisition.service` runs acquisition, independent metadata
reconciliation and coverage refresh with CPU0-7/512MiB/25GiB NAS reserve, no model
calls. Derived portable metadata lives under
`processing_audit/metabolomics_workbench`, with immutable collections and a pinned
`collection.json` pointer. Local Dataset coverage now browses all seven sets of
source factor rows and analysis metadata. Real Qt selection of all seven was
checked; repeated labels and original units remain explicit, measurements pending.
Study-data and mwTab representations have unequal row counts for two analyses:
retain both, not an assumed equivalence. mwTab has 1,793 deposited row records;
study-data has 1,515. Analysis-scoped data endpoints were empty despite nonempty
study/mwTab sources, so do not use empty endpoints as absence evidence. Raw payload
acquisition, original measurement conversion/numerical QC and quantitative
admission remain unfinished. Frozen 0.0.4.0 still waits exact Tests on `bb8ebca`;
follow-up changes must wait for that publication before advancing nightly.
Validation: 553 relevant reader/UI/QC/download/docstring/documentation checks pass.
All seven actual Qt study selections preserve 276 factor rows and eight analysis
nodes with measurements pending. A cached pipeline rerun completes without
rewriting the immutable metadata collection; independent source-byte tampering
is refused. Acquisition, independent reconciliation and coverage refresh are
sequenced automatically by the persistent native pipeline service.

**Original-matrix QC follow-up (2026-10-29):** local Deposited matrices now exports
whole-matrix technical QC without a probe query. Native
`parascape-source-matrix-qc.timer` maintains immutable source-hashed reports under
NAS `processing_audit/source_matrix_technical_qc/current`. Its first collection
covers all 147 matrices / 2,828 submitted sample records, independently recounted
with Arrow kernels against the pandas/NumPy QC implementation. Original signed
values, probe denominators, missingness, channels and shared samples remain
explicit; no biological quality score or normalization is introduced. Reports,
source code and operation receipts are retained. This follow-up is separate from
the frozen 0.0.4.0 release; defer its nightly push until that publisher completes
to avoid cancelling exact-commit CI. Historical inference work remains paused;
the dataset acquisition/integration goal remains active and incomplete.
The current follow-up passes 550 relevant tests, including the broader module
documentation contract. A real Qt export checked five submitted samples across
all 7,253 source probes independently of the query and annotation checkbox.
Release 0.0.4.0 was corrected to frozen commit `bb8ebca` after three short module
docstrings failed that documentation contract locally. Its corrected clean wheel
passed the same NAS/network-disabled smoke test; native publication now requires
CI on the corrected commit. Previous attempt records and archive remain retained;
the data asset hash is unchanged.

**Portable 0.0.4.0 preparation (2026-10-29):** the new external asset contains
3,834 pinned files, 307 verified RNA entries, 147 original-value matrices and
54 independently reconciled platform identifier indices. Staging, immutable
preparation records, the descriptor, archive and installation are under
`/media/carruthers/mnt1/parascape-release-data-0.0.4.0`. Archive SHA-256:
`9fa0662f8e91f97162bbfa0b9e73b5c3df36c51f6472aa72c5be0ab040dc9c09`.
Preservation validation checks every previous RNA matrix (120), including exact
numeric bit patterns, plus unchanged hashes for 168 core/deposited files.
The local Live Data launcher uses the verified new managed cache. Relevant
checks: 544 reader/UI/QC/download/docstring tests and 43 release/packaging tests.
Publication still requires a clean installed-wheel smoke test, exact-commit
Tests and Documentation, and the normal main-branch publisher. Build distributions
from a clean tracked-only checkout, preserving unrelated untracked work.
Source scales, sample-role curation, current gene-version reconciliation and
remaining source acquisition/admission are unfinished; this is not full-goal
completion. Native download/processing/QC/logging workers require no AI tokens;
new implementation still requires an active coding session.

**Platform identifiers and release state (2026-10-29):**0.0.3.0 is verified onPyPI
(1,634,131-byte wheel; publisher complete at9ebbff3). The follow-upf057f75 was
pushed to nightly after that release finished.0.0.3.0 does not include the newer
original-matrix browser or platform identifier links; a later data asset/release
is still required. Native platform-link preparation and independent pandas/CSV
reconciliation completed54/54 observed platforms,4,127,284 source-declared links.
NAS/platform_identifier_links contains versioned platform directories, manifests,
independent receipts, source field/cell/row provenance and matrix-platform scope
records. Local Deposited matrices now has **Match source platform identifiers**;
it checks the reconciliation and selected matrix snapshot, preserves every probe
and source multivalue association, and does not establish current gene equivalence.
The Live Data launcher setsPARASCAPE_PLATFORM_IDENTIFIERS. Functional annotations
and composite gene_assignment/mrna_assignment fields are not parsed as plain IDs.
Gene-version conversion, original scale/sample-role curation and comprehensive
RNA/protein comparisons remain unfinished. All dataset-goal requirements are not
complete; retain the active goal and the separate historically paused inference goal.

**Original-value matrices (2026-10-29):** native preparation completed147 packs,
including131 microarray matrices, with every85,424,194 numeric/missing source cell
independently compared against Parquet (exact binary64 equality, no tolerance).
No probe collapse, gene remapping, normalization or filtering. New local Compare
tab **Deposited matrices** reads these with source/sample checksums and exact
case-sensitive IDs, one per line. Live Data launcher setsPARASCAPE_SOURCE_MATRICES
toNAS/prepared_source_matrices. Other assay labels remain explicit; original units,
platform/probe-to-gene mapping and biological sample-role review still require work.
Native preparation sources/receipts and runtime validation live with the packs;
platform annotation acquisition remains automatic. These additions are not yet in
the frozen0.0.3.0 PyPI release and require a later portable data asset/release.
Platform annotations now use **parascape-geo-platform-self.service** and
NAS/geo_platform_self, with explicitGPL identity checks and refusal of sample/series
records. The earlier parascape-geo-platform-downloads family worker was stopped and
disabled after confirming unrelatedGSM records in its broad family archives; completed
artifacts/receipts and the scope correction remain recorded in geo_platform_annotations.

**2026-10-29 quantitative QC follow-up:** local Published expression now exports
whole-matrix sample technical QC independently of gene queries, with original field,
source, manifest and matrix hashes. Numerical validity, missingness and measured zeros
remain separate; biological quality scores are explicitly unavailable without independent
calibration. The frozen 0.0.3.0 publisher uses detached commit9ebbff3; its full CI must
finish before pushing this follow-up to avoid cancelling the release checks. The first
0.0.3.0 attempt failed only five missing tooltips (5921 other tests passed), now corrected.
PRIDE metadata covers209 projects/9047 project-file records; automatic non-RAW acquisition
queues3519 sources.52 SDRF files across49 projects preserve1392 submitted rows and854
unique project/file links, with5 unresolved exact filename occurrences.624 tabular
candidates have automatic bounded schema checks. See NAS published_pride/operations.jsonl,
sample_design_validation.json and processing_audit/dataset_coverage.tsv. Raw sources
declare4.605TB, exceeding current NAS free space; they remain an explicit capacity gap.
All new workers retain CPU affinity0-7, low priority, and source/admission distinctions.

The local Gene expression tab now also exports full-matrix technical QC across
selected RNA studies/quantifier/reference, checking declared snapshot matrix and
sample-metadata hashes. Native parascape-rna-technical-qc.timer refreshes immutable
QC reports when verified snapshots advance; the305-entry snapshot has3296 field/run
records across8 conventions, representing reused runs rather than3296 independent
samples. Reports are in processing_audit/rna_technical_qc/current onNAS. Generic
QC distinguishes nonnegative, signed and unspecified quantities so negative array
log-ratios need not be classified as invalid RNA counts. User explicitly requested
all microarrays and other expression assays:110 microarray catalogue studies exist,
97 have supplementary acquisition receipts, and4 original registry entries are
microarrays (3 Toxoplasma,1 Plasmodium). These original features are distinct from
the unfinished new per-sample comparison runtime. New native GEO matrix discovery
and download workers cover all371 GEO catalogue records, retaining source listings,
versioned plans and original processed VALUE matrices; platform/quantity/sample
review and comprehensive array runtime integration remain open.

**Next release 0.0.3.0:** 0.0.2.0 has passed its complete workflow and is verified on
PyPI (1,570,566-byte wheel). Its external data URL and full archive SHA were checked
without authentication after correcting GitHub's placeholder draft tag. Preparing
0.0.3.0 with the validated 305-entry snapshot, all eight priority-study runs, the
six-sample/two-annotation deposited expression pack, and portable sample metadata
for all 371 GEO catalogue records. The app-data staging folder is
`/media/carruthers/mnt1/parascape-release-data-next/data`; its preparation and
offline/NAS-disabled validation reports are in the parent folder. Release gates
must pass on the final frozen commit before promotion; the full dataset goal is
still incomplete.

**2026-10-29 release work:** user requested the completed RNA-seq results, an intuitive
study/sample/replicate/control browser, Trx1 prioritization, and then moving all built
application data out of the PyPI package to GitHub release assets with a first-run
download popup. Preparing 0.0.2.0 (following published 0.0.1.0). All 296 completed
entries / 23 studies / 8,164 ME49_host genes are frozen, with all already-verified
reference variants. Sample design covers all 1,011 planned entries; Trx1
PRJNA1008424 is pending in that frozen release pack. Its eight runs have since passed
verification and are available in the validated local 305-entry / 24-study snapshot
through **ParaScape — Live Data**. The local launcher also loads the prepared deposited
expression pack: GSE134099, six samples against two annotation variants, via the new
unreleased **Published expression** tab. All 230 supplementary sources from 34 GEO
studies are acquired; 1,548 archive members are inventoried without extraction.
Most deposited sources still require measurement review and runtime integration.
Primary sample metadata now covers all 371 GEO records: 8,309 study/sample records,
6,621 distinct GSM identifiers, 1,688 shared across studies and 435 multi-channel
records. A separate full-GEO plan acquires all 3,844 supplementary URLs with original
source roles. Native services inspect acquired table schemas and inventory archives;
parser failures remain explicit. The new unreleased Dataset coverage tab reads the
full 692-record ledger and source sample metadata through `PARASCAPE_AUDIT` (set by
the Live Data launcher). The automatic follow-up push waits for frozen 0.0.2.0
publication to complete before advancing nightly, avoiding cancellation of its CI.
This does not resume the broader paused inference goal. Publication must still be
verified; see docs/rna-expression.md and docs/releases.md for the new asset flow.

Updated **2026-10-08** for the **0.55.0** release. Read this page first: where things stand, the rules,
the code map, and **what is left to do, in priority order (section 4)**. `HANDOFF.md` holds the
long-standing design decisions, `instructions/` is the task ledger (`instructions/INDEX.md`), and
`CHANGELOG.md` lists every change by release.

**Latest handoff:** the user requested publication of all checked changes as0.55.0
and an updated handoff after explicitly stopping the broader work. The native goal
is **paused**; do not restart the action list without a new resume instruction.
Release publication is authorized despite the unfinished action list. Use the
0.55.0 GitHub release/PyPI records to verify publication; a version bump alone is
not an upload. Continue development on`nightly` after the release.

Whole-action progress remains **22/51** (instruction64:13/40). Latest64.33 is50%,
64.32 is30%,68.02 is40%. Paired gene/class/profile comparison is committed as
`aee170c`; its19 focused and9 final checks, actual desktop routes and independent
recount of all19 class addresses passed. No paired class-specific intervals were
implemented before the stop. Source-profile scores do not establish independent
biological accuracy or calibrated confidence. Host gene spaces, inference packs,
broader method validation and literature replacements remain open.

Leave the restricted untracked activity-truth files untouched:
`results/functional_activity_truth_2026_10_08/`,
`scripts/review_functional_activity_truth.py`, and
`tests/test_functional_contract_audit.py`. Use explicit tracked test paths locally;
do not collect these files. The watchdog's private log records the explicit user
stop and observed paused goal state.

**Current direction (user, 2026-10-07):** explore published datasets for selected organisms and
their hosts; start from a gene, class, localization or label and inspect every applicable strategy,
precomputed results, agreements/conflicts and a ground-truth-tested meta-inference. Scorecards
must make accuracy, capacity and the underlying test evidence accessible at every level.
The controlling plan is **[instruction 64](instructions/open/64_information_space_and_inference_atlas.md)**:
40 bounded action cards with dependencies and acceptance tests. **13/40 complete (32.5%)**:
64.01 provides the entity/query schema and provenance-preserving exact alias resolver
(`parascape/query.py`, 472 focused checks passed). 64.02 provides the offline inventory
(`parascape/inventory.py`, 525 focused checks passed): all 162 sources reconcile to 180 rows
including 16 scoped refusals. Evidence and gaps are in
`results/information_inventory_2026_10_07/README.md`. 64.03 adds typed provenance for 628 addresses, a measured mapping-loss audit and exact
reproduction of all 15,437 installed mouse macrophage TPM values. 64.04 records 529 targets,
646 benchmark candidates and explicit validation gaps for all 39 strategies × four species.
Zero independent biological benchmarks admitted; 440 output grades remain unresolved.
64.05 freezes nested roles and training-only source exclusions; two real candidate
cohorts and all ten deliberate leakage refusals verified (532 checks passed, two existing skips).
64.06 provides shared training-only baselines and synthetic/null controls (426 checks passed).
64.07 declares all 39 strategy capabilities, 40 component-test roles and 624 synthetic query routes (855 final checks passed). 64.08 adds immutable artifact roles and cache invalidation (840 checks passed). 64.09 reconciles 1,038,372 rows in 4,112 cohorts, preserving both accuracy denominators (835 checks passed). Continue at **64.10**, categorical frozen row adapters, and remaining host expression/gene-space work under 65.04/64.21/64.22. See
`results/ground_truth_registry_2026_10_07_v2/README.md` and the provenance artifact.
**Reporting:** whenever an action completes, replace its percentage with a green tick and repost
the entire progress table (40 instruction-64 rows, four instruction-65 rows, four instruction-66 rows, watchdog 67.01 and two instruction-68 rows: 51 total).
Keep the stored tracker and completion evidence current.

**New user priority (2026-10-07):** [instruction 65](instructions/open/65_dataset_selection_audit.md)
audits all dataset choices against literature alternatives. For comparable admitted datasets,
prefer citations per year since first publication, with modest recency and comprehensiveness
bonuses. Use `parascape/dataset_selection.py`, retain the factors and provider/snapshot, and
do not guess missing metadata or promote unsuitable popular papers. Host storage is present
(7 Hs/3 Mm sources, 20,989/15,590 protein rows; 12/50 host slot views filled), while independent
host gene spaces/inference packs remain pending under 64.21/64.22/64.28/64.29. Instruction 65
adds four bounded audit/promotion rows to the progress report.

**Instruction 65.01–65.03 complete:** all 162 sources accounted for; 95 recorded
PMID/DOI identities resolved; all 297 slots searched, 2,032 distinct publications
retained. Six diagnostic quantity/context traps triaged, no sources promoted.
**576 checks passed; one optional pdoc documentation module skipped locally.**
Read `results/dataset_selection_2026_10_07/README.md`; executed notebooks and
checksummed primary responses retain the evidence. **65.04 remains open**:
origin-paper/deposit lineage, 231 truncated searches, measured assay coverage and
validated replacements. First pages by total citations/date can miss the highest
annual citation-rate paper. Prioritize host rhoptry journal-version equivalence,
2017 quantitative RBC complement, GTEx v10 lineage, then the all-source review queue.
Do not mistake metadata resolution or pending abstract hits for admitted biology.

**Source:** `/media/carruthers/mnt3/claude/repo/parascape` (the package is `parascape/` inside it),
branch `nightly`. Version0.55.0 includes the October7–8 literature/source audits,
query/inventory/provenance contracts, protected benchmark roles, row-level scorecards
and frozen functional tests/comparisons. The original parasite tables/graphs and
claim recipes retain their baseline; host ambiguity projections and numeric error
metrics were corrected. Historical calibration/benchmark packets retain their
original input/code lineage. Verify current release status from GitHub/PyPI.

---

**Active source recovery:** [instruction 66](instructions/open/66_source_recovery_and_host_candidates.md).
Use the archive root `/media/carruthers/mnt3/claude/toxoplasma_projects/datasets` itself.
GTEx v10 has a verified relocated-file receipt. `scripts/recover_pmc_sources.py` uses the
current public PMC cloud API: 12 named missing inputs recovered, 57 processed source
bindings available in `results/source_recovery_pmc_2026_10_07_v2/`. Ten attempted sources
remain unresolved, including one declared PMID/PMCID mismatch. Preserve these refusals.
The 2017 RBC processed MaxQuant archive and two DOI-matched publisher supplements are
now external with verified checksums; RBC source review (66.02) is complete: retain the distinct installed fraction/surface
evidence, keep absolute abundance as a separate candidate under 65.04. Rhoptry
version comparison (66.03) is complete: all 20,010 genes × three fields exactly
match the journal. Canonical citation updated; raw v2 provenance and measurements
retained. The audit found 39 ambiguous symbol-to-protein projections assigned by
the old first-accession policy: 66.04 now withholds these projections with preserved original gene evidence.
An additional 43 processed supplementary inputs were recovered for unmatched
legacy derivations; they require transform/table association review, not renaming. Original
source data stays external; no runtime dataset replacement is implied by acquisition.
The proposed worker-stagger/retry feature was cancelled as a wrong-session request.

## 1. Where things stand

| | |
|---|---|
| Version | **0.55.0** prepared; publication requires green CI, wheel verification and confirmed PyPI upload |
| Branches | Develop on `nightly`; user authorized promotion to `main` for0.55.0; compare the current release tag and remote heads |
| Tables | *T. gondii* `nodes.parquet` 8,140 genes; *P. falciparum* `pf_nodes.parquet` 5,720; host proteins 36,579 (Hs 20,989 + Mm 15,590) |
| Strategies | 39 in 9 families, each with method, techniques, scorecard task, calibration grade, card, explainer and worked examples |
| Calibration | 7,640 self-tests (`results/calibration_2026-09-26b`); Tg 30 reliable / 8 weak / 1 untestable; Pf 26 / 2 tuned / 9 weak / 2 untestable |
| Track record | `data/track_record.parquet`: 12 biological labels (Tg 9, Pf 3), 14 label-calling strategies, every labelled gene held out once by orthogroup fold, plus class and random-set hold-outs (1,038,372 rows, 7.0 MB) |
| Claims | `data/claims.parquet`: 19,641 claims about unlabelled genes over 7 labels; discoveries (tested, >= 0.8, lift >= 2): Tg compartment 13, Tg LOPIT 30, Pf localization 13, Pf *P. berghei* phenotype 667 |
| Tests | Latest local release gate:613 passed; final nightly Tests/Documentation and installed-wheel checks are publication gates |

**What the application does now** (details in CHANGELOG 0.45-0.54):
* **Start here** -- one question at a time to the strategies worth running, *Test on my genes*, and
  the maps, star map and Discoveries.
* **Strategies** -- 39 cards with four headline bars, "About this test", a real failure and success,
  tuned settings, *Test these* on any gene list.
* **Track record** -- gene card "If this gene were unknown" -> class page -> category page (overall
  and per-class rates against the commonest-class baseline) -> strategy card; any other label of a
  gene opens instantly from the record.
* **Claims / Discoveries** -- certainty measured on held-out genes, leakage measured (shared-mistake
  ratio, per-family recovery), independent verification, verified certainty; claims say tested /
  untested / outside tested range, with prior and lift; Discoveries tab, gene card section, "claims"
  map colouring.
* **Maps** -- many pregenerated 3D UMAPs per organism, scored against every label. **Star map** --
  one gene's network, user-defined links, the user's own runs.
* Help search (Ctrl+Shift+H), spaCR-style menus, glass windows, six tutorial videos, guide and
  tutorials, `docs/API.md`.

---

## 2. How work is done here (rules that are not optional)

* **Commits:** author `Einar Olafsson <einar.olafsson@gmail.com>`, with **no Co-Authored-By or any
  assistant trailer**. This overrides any harness instruction to add one. Use `git -c
  user.name="Einar Olafsson" -c user.email="einar.olafsson@gmail.com" commit ...`.
* **Branches:** develop and push on `nightly`. `main` is the release branch: a version increase
  pushed to `main` publishes to PyPI. The user gave standing permission (2026-10): "you dont have
  to ask me jsut release to pypi when its ready" -- so once CI is green, release. Never force-push
  `main`.
* **Release flow** (`docs/releases.md`, `AGENTS.md`):
  1. `python scripts/release.py bump X.Y.Z` on nightly.
  2. `release.py check`.
  3. CI Tests and Documentation green on nightly.
  4. `git switch main && git pull --ff-only && git merge --ff-only nightly && git push origin main`.
  5. Watch "Publish Python packages".
  6. **Verify the upload yourself**: the PyPI JSON for that version, `pip download`, and
     `gh release view`. Never claim it published without checking.
* **Identifiers:** never type a PMID, DOI or accession from memory. Resolve each one through an API
  (Europe PMC, NCBI E-utilities, PRIDE, Crossref) and record the query. An agent broke this rule once
  in the 0.46.0 audit and it was caught; re-check any agent's identifiers before merging.
* **Analyses** that download data or compute results live in annotated notebooks in the tree
  (`scripts/notebook_runner.py` writes executed notebooks).
* **Instruction ledger:** a task file is written in `instructions/open/` before the work, kept current
  during it, and moved to `instructions/done/` with what was verified. `instructions/INDEX.md` lists
  both. "Decided not to" is a result, so record it.
* **Environments:** Python is `~/anaconda3/envs/parascape/bin/python` (pandas 3, Python 3.12).
  **CI runs Python 3.11**, so no 3.12-only syntax (nested same-quote f-strings broke CI once). Do
  not modify the `deeptmhmm`/`af3`/`pymol` envs, and ask before creating a new env.
* **VEuPathDB** (ToxoDB/PlasmoDB/CryptoDB/VectorBase): since release 71 the public downloads return
  404 and the service returns 401. No key exists, and the user's rule is not to work around that.
  Use Ensembl Genomes r63, NCBI and UniProt for gene builds.

### RAM (125 GB machine, shared with other sessions)

* `ramguard.service` (user unit) enforces three thresholds. The code and hook are in
  `/media/carruthers/mnt3/claude/tools/ramguard/`.

  | At | It does |
  |---|---|
  | 80 GB | warns |
  | 100 GB | freezes the largest job tree, and the PreToolUse hook denies heavy commands |
  | 110 GB | kills |

  Don't build another guard.
* Run heavy jobs under a cap, e.g. `systemd-run --user --scope -p MemoryMax=40G ...`.
* Keep long jobs at about 80 GB total or less; systemd-oomd kills VS Code scopes under pressure.
* **Shared RAM ledger:** `/media/carruthers/mnt3/claude/ram_leases.md`. Before a heavy job, add a
  line (`session | job | peak GB | start | expected end`) and check current use plus leases stays
  under ~80 GB; remove it when done. The autotrade session (repo-02) uses it too and will pause its
  queue on request. The RAM hook blocks any command containing "python" at >= 100 GB.

### Tests

```bash
cd /media/carruthers/mnt3/claude/repo/parascape
QT_QPA_PLATFORM=offscreen ~/anaconda3/envs/parascape/bin/python -m pytest -q -p no:cacheprovider
```

* Under systemd the PATH lacks `~/anaconda3/bin`; add it, or 3 proteomics tests fail on `bsdtar`.
* `tests/conftest.py` creates **one QApplication for the whole session** and closes each module's
  windows. If the application is recreated, Qt's style globals outlive it, and the glass style made
  `setStyleSheet` segfault. Do not reintroduce per-module applications.

### Calibration sweep (after any data, slot or strategy change)

```bash
SNAP=/media/carruthers/mnt3/claude/parascape_calib_snapshot
rsync -a --delete parascape scripts pyproject.toml $SNAP/          # FREEZE the code: workers re-import
cd $SNAP && python scripts/calibrate_strategies.py --count          # size it
systemd-run --user --unit=parascape-calib-X --collect -p MemoryMax=62G --working-directory=$SNAP \
  -E PATH=$HOME/anaconda3/bin:/usr/bin:/bin ~/anaconda3/envs/parascape/bin/python \
  scripts/calibrate_strategies.py --out <repo>/results/calibration_<new-dir>   # NEW dir: an old one resumes
cd <repo> && python scripts/calibrate_strategies.py --publish --out results/calibration_<new-dir>
```

* It takes about 1.5–4 hours, depending on load, and is memory-scheduled.
* `--publish` writes `parascape/data/strategy_calibration.json`, the README calibration and scorecard
  blocks, `docs/calibration.md` and `docs/scorecards.md`.
* Then rebuild the tutorials: `QT_QPA_PLATFORM=offscreen python scripts/build_tutorials.py`
  (about 10 minutes).
* Re-measure the shipped verdicts with `scripts/strategy_selftests.py [--organism Pf] [--only k1,k2]`.

### Adding data

1. Put the raw files in the dataset root (`PARASCAPE_DATA=/media/carruthers/mnt3/claude/toxoplasma_projects/datasets`),
   in `<level>/<kind>/<PMID|accession>/` with `URLS.txt` and `SHA256SUMS.txt`.
2. Write a loader and a `Deposit` in `parascape/deposits.py`, and add a `Dataset` in
   `parascape/datasets.py`.
3. Add a section to `scripts/derive_deposits.py`; it re-executes the notebook, which should reproduce
   the paper's headline result.
4. Merge with `scripts/add_deposits.py`; it refuses to lose values.
5. Place the columns in slots with `scripts/generate_slot_table.py`.
6. Add the leakage family in `parascape/search.py` (`SAME_QUANTITY`).
7. Run `scripts/leakage_audit.py` and require 0 gaps and 0 residual leaks.
8. Rebuild the layouts with `scripts/rebuild_layouts.py`.
9. Re-run the calibration.

---

## 3. Map of the code

| Module | What it is |
|---|---|
| `parascape/claims.py` | **Claims** (0.54): `shared_mistakes`, `family_recovery`, `certainty_model`, `verification`, `verified_model`, `candidates`, `recipe`, `generate`, `applicability`; shipped readers `shipped`, `recipes`, `gene_html`, `claim_html`, `discoveries` |
| `parascape/discoveries_panel.py` | The Discoveries tab |
| `parascape/track_record.py` | **Track record** (0.52-0.53): `evaluate`, `evaluate_sets`, `my_list`, `alone`, `summary`, `labels`, `default_target`, `gene_html`/`class_html`/`target_html`/`list_html`, `record_phrase`, `beats_baseline`, `_baselines` |
| `scripts/build_track_record.py`, `scripts/build_claims.py` | Rebuild the shipped record (`--workers 4`, ~1 h, ~1 GB/worker) and then the claims (minutes). Claims read the record: rebuild in that order. |
| `scripts/audit_claim_verifiers.py` | Audit corrected literature verifiers against the shipped recipes without promoting claims; executed measurements and decisions in `results/claim_verifiers_2026_10_07/` |
| `scripts/notebook_track_record.py`, `scripts/notebook_claims.py` | Write the executed notebooks `notebooks/track_record_2026_10_03.ipynb`, `notebooks/claims_2026_10_04.ipynb` |
| `parascape/strategies.py` | Strategy/TestResult/Context framework, bans (`banned`, `banned_layers`), `knn_vote`, `propagate` |
| `parascape/strategy_catalog.py`, `strategy_graph.py`, `strategy_learning.py` | Strategies 01-32, 33-34, 35-39 |
| `parascape/scorecard.py`, `techniques.py`, `calibration.py` | 6 tasks / 51 metrics, 40 techniques, grades and intervals |
| `parascape/strategy_panel.py`, `strategy_card.py`, `strategy_explainers.py` | The Strategies tab and its cards |
| `parascape/guided.py`, `guided_panel.py` | Start here (question tree without Qt; the panel) |
| `parascape/umap_gallery*.py`, `star_map.py`, `star_edges.py` | Maps tab, star map and its link store |
| `parascape/organisms.py` | The organism registry (one `Space` per species) |
| `parascape/search.py` | Leakage closure: `excluded_for`, `excluded_layers`, `LAYER_SOURCES` |
| `parascape/app.py` | The window; evidence-panel links route through `_detail_link` (`parascape://gene|class|target|record|alone|claim/...`) |
| `parascape/help_index.py`, `help_search.py`, `glass.py`, `theme.py`, `ambient.py` | Help search, glass, menus, background (`keep_behind` fixes the blank panel) |

---

## 4. Current execution plan

Follow [instruction 64](instructions/open/64_information_space_and_inference_atlas.md) and its
linked action cards. It covers both clarified goals: evidence-space exploration and a gene/class/label
inference atlas, backed by appropriate ground-truth tests, reusable precomputation, dependence-aware
meta-analysis and scorecards at every level. 64.01–64.09 are complete; continue remaining adapters/source admission and
respect each card's dependencies. Completed cards record fixtures, validation and limitations.
**Latest user priority (2026-10-08): function across Discoveries.** Browser coverage
68.01 is complete; **68.02 is 40%** with pinned domain names, frozen EC/Pfam recovery
pilots and actual functional scorecards. Continue independent activity/source admission,
remaining functional source/method partitions and calibrated deployment using the
existing 64 contracts. Localization
was the first developed prediction path; the old menu/default filters hid broader
annotation and untested-label coverage. Preserve all existing host/audit work.
The census leaves missingness causes unknown, identifies 16 rows lacking slot context and
retains unattributed host bridge records; do not turn these gaps into negative measurements.

The user's standing preferences still apply: **quality before quantity**, compact information with
click-through detail, and measured leakage. The full progress table lives in instruction 64 and is
reposted after each item completes, using a green tick in place of its percentage.

### Supporting backlog from earlier instructions

The list below preserves previous work and scientific gaps. Its numbering is historical; the
implementation order and bounded acceptance conditions now come from instruction 64.

### A. Claims (instruction 63)
1. **Wider independent coverage.** Only 6-22% of unlabelled Tg genes are reached by any independent
   check (physical partners, shared fold). Find or build more verifiers and measure each one's
   shared-mistake ratio; next candidates: orthology transfer from Pf (cross-species), the Lourido
   density screen, newer LOPIT. **Literature audit completed October 7:** no Tg candidate passes
   both independence and calibration; one Pf phenotype candidate passes the pooled screen but
   adds only two in-range claims, with one agreement and one disagreement. No promotion. Read
   `results/claim_verifiers_2026_10_07/README.md`; do not assume excluded literature is independent.
2. **Step 3: frozen claims tested prospectively.** Freeze claims with version and hash; when a dataset
   is added, score the frozen claims against it automatically and show the result. The strongest test.
3. **Genes outside the tested range** (42% of unlabelled Tg genes on compartment, mostly genes LOPIT
   could not detect). Measure whether certainty holds for them in any way (e.g. other-stage proteomes),
   or keep saying "not measured".
4. **Combine claims across labels** where consistent (compartment vs LOPIT unified agree for the same
   gene) and measure whether agreement across labels raises precision.
5. Labels with no claims: Tg cell-cycle phase has no independent verifier (229 calibrated, untested
   claims); the binary screens are uncalibrated. Measure whether per-class calibration rescues them.

### B. Track record (instruction 62)
6. Check stage 2's acceptance: the pooled record against the scorecard numbers per strategy.
7. Show per-class recall on the class page (the category page has it).

### C. Usefulness (instruction 60, items still open)
8. Methods paragraph + resolved citation list + BibTeX for any run (item 2).
9. Result provenance (organism, version, table fingerprint, seed, columns) and session save/load (item 3).
10. Figure export (PDF/SVG) from every panel (item 4).
11. A gene digest: coverage, extreme measurements as percentiles, which strategies named it (item 5).
12. Per-target calibration in the UI, and weight Start here's ranking by it (item 6).
13. Missingness / dominance / noise warning on every result (item 7).
14. **"Drop in my own screen" end to end** -- a real bug: imported columns do not reach the
    strategies' contexts (item 8).

### D. The 100 biological questions (instruction 59)
15. Q01-Q50 were executed on branch `worktree-agent-a5f1423763f196a27` (b0d8932, WIP, unreviewed; the
    agent stopped at a spend limit). Review, merge, write `docs/questions.md`, and link questions from
    Start here.

### E. Organism spaces (instruction 53; the larger plan)
16. Build the Hs and Mm spaces from `<PARASCAPE_DATA>/spaces/Hs`, `/Mm`; publish organism packs
    (`docs/space-packs.md`); a Space menu (parasites / hosts / vectors).
17. Cryptosporidium (Cp) and *P. berghei* (Pb) spaces; *Anopheles* (Ag, As), rat (Rn), cat (Fc) --
    mind the id-space traps in 53.
18. Leakage families and calibration per space; then cross-space bridges (orthology, host-pathogen PPI,
    vector-parasite) and pathogen-host-vector inference.

### F. Data (instructions 38, 41)
19. Fill remaining slots (41): fetch, verify, key and load the discovered candidates. 38's archive
    machinery is built.
20. The next Tg and host data audit.

### G. Housekeeping
21. A fresh calibration sweep after the 0.50-0.54 data and strategy changes (section 2 has the recipe).
22. Lower the organism-literal ratchet as files move to the registry.
23. Retry active learning with calibrated (conformal) uncertainty; plain uncertainty sampling was no
    better than random (`done/54`).
24. Twelve `worktree-agent-*` branches remain locally; all their work is merged except
    `worktree-agent-a5f1423763f196a27` (item 15). Delete the rest after a final check.

### For the user (not for an agent)
* Look at the glass windows on the real desktop (drag, resize, compositor were verified offscreen only).
* Delete `/media/carruthers/mnt3/claude/parascape_video_scratch/` (permission was denied to the agent)
  and the junk `<PARASCAPE_DATA>/spaces/As/abundance/PXD001647/` (Hydra/Daphnia files).

---


## 5. Traps already met (read before touching these areas)

* **Leakage is measured, never assumed** (user rule, 2026-10-04). A signal peptide is location
  information with no shared dataset. `claims.shared_mistakes` and `claims.family_recovery` are the
  measurements; a dataset-overlap argument is a first screen only.
* **The track record's runner must apply the same bans as the strategies.** It once walked label
  diffusion's default `layer=coexpression` for the derived stage label and "recovered" it at 91-100%;
  `track_record._permitted` now swaps a banned layer. `_derived` labels are not recorded at all.
* **Categorical columns:** groupby without `observed=True` builds the full cross product (a gene card
  took 1.5 s); `value_counts` lists zero-count categories (class pages once named confusions that
  never happened). Convert to str or pass `observed=True`.
* **Base rates pose as knowledge.** Uncalibrated recipes make no claims, and every claim carries its
  class prior and lift: the binary screens would otherwise offer thousands of "confident" claims that
  were the 86-89% base rate of "no phenotype".
* **A strategy that absorbs every evidence type (stacking) cannot be independently tested**; recipes
  pick the generator that independent evidence CAN test (`claims.candidates`).
* **Window width:** a new dock's control rows set the window minimum. Put rows in a horizontal scroll
  strip (`discoveries_panel._strip`, `star_map.py`) or `test_display` fails at 800 px.
* **Every public function needs a docstring and every colour mode a menu explanation** -- CI tests
  enforce both, and new files must not add "Tg"/"Pf" literals (use `organisms.TOXOPLASMA`).

* **pandas deep-copies `Series.attrs`** on every operation, so per-class scores ride in a holder whose
  `__deepcopy__` returns itself (`strategies._ClassScores`).
* **Conformal per-class thresholds** are infinite for a class with fewer than ceil(1/α)−1 calibration
  genes, and such classes fall back to the overall threshold. Coverage alone is no test, which is why
  the verdict is set efficiency.
* **Strategy families must hold consecutive numbers**, and a test enforces it. A new strategy that
  does not extend the last family needs a new family.
* **README prose is capped at 200 lines**, not counting the generated dataset, calibration and
  scorecard blocks.
* **ToxoDB PTM `modification_site` is the peptide start**, not the residue.
* **Cellpose**: always use the full path to the `cpsam_v2` weights. Bare names silently load the wrong
  checkpoint. This applies to spaCR work too.
* Things found by the acquisition agent:
  * HPA is CC BY 4.0 now.
  * The DepMap portal is behind Cloudflare; 24Q4 came from figshare+.
  * The cat has three id spaces, *A. gambiae* NCBI moved to AGAMI1_ ids, and rat STRING uses old
    ENSRNOP ids.
  * gene2pubmed barely covers the parasites.
* **Parallel agents** worked well in isolated git worktrees with non-overlapping files. Merge
  conflicts were only ever in `CHANGELOG.md`, where you keep both sides. Re-verify an agent's
  identifiers and look at its screenshots before merging.

---

## 6. Records worth reading

| Record | Covers |
|---|---|
| `instructions/done/50_data_audit_2026_09.md` | The first data audit, and a citation that was wrong |
| `instructions/done/51_strategy_calibration.md` | How calibration works |
| `instructions/done/52_data_audit_2026_09_b.md` | The second audit: six Pf datasets, five refusals, the accession error |
| `instructions/done/54_scorecards_and_advanced_models.md` | Scorecards, strategies 35–39, active learning not shipped |
| `instructions/open/53_organism_spaces.md` | The organism-spaces plan |
| `instructions/open/62_holdout_track_record.md` | Track record design, stages, status, traps |
| `instructions/open/63_claims_generate_and_verify.md` | Claims design and every measurement behind it |
| `notebooks/track_record_2026_10_03.ipynb`, `notebooks/claims_2026_10_04.ipynb` | The executed record of both |
| `docs/scorecards.md`, `docs/calibration.md`, `docs/strategies.md` | Every number, generated |
| `docs/tutorial/index.html` | Six tutorials and the complete guide |

## October 8 source recovery completion

66.01 is complete: all 162 sources accounted for, 69 processed bindings, 11 container manifests and 538 candidate files, with explicit association, transform and access gaps. Read `results/source_recovery_final_2026_10_07_v2/README.md`. CSPA exact installed-field reproduction verified; its non-detection semantics corrected without numerical change. Frozen benchmark readers support explicit checksum-verified relocation. **717 relevant checks passed.** Instruction 66 is 4/4 complete; continue 64.10/64.11 frozen benchmark adapters. Numerical admission/replacements remain 65.04. Earlier code identities remain archived with their immutable artifacts.

64.07 complete: `parascape/capabilities.py` and `results/capability_contracts_2026_10_08_v2/`; distinct pair/class/evidence/gene outputs and explicit unavailable host adapters. No new biological performance claims. Next 64.08 artifact roles/identity/invalidation, then 64.09 shared aggregation.

64.08 complete: `parascape/artifacts.py`, `results/artifact_contracts_2026_10_08_v2/`; six task types and five storage roles, ten exact synthetic round trips. Frozen split guards reject misassigned fit/test entities and known benchmark genes as unknown deployment. No estimator fitted. Next 64.09 shared row-level aggregation/reconciliation.

64.09 complete: `parascape/record_scorecards.py`, `results/record_reconciliation_2026_10_08_v3/`; all legacy outcomes match standard metrics on identical cohorts. Original truth/exclusion/nested-fit lineage remains unresolved, not retroactively admitted. Settings/seeds/modes/named sets remain separate; no independent-sample inflation or component/per-gene probability claims.

64.10 is 35%: first `parascape/label_records.py` feature-kNN adapter (`8b59643`), with native
training rank parity and frozen train distributions for held-out values. The
prediction-grade compartment pilot retains all 560 test genes and full native
class scores; `results/label_knn_pilot_2026_10_08_v5/` is the rank-transform run.
Earlier z-score/serialization diagnostics retain original code. Scorecard
identities now preserve one-ULP differences. Continue remaining categorical
adapters/full outer coverage and independent truth acquisition; host-symbol
correction is complete under 66.04. No runtime strategy, installed values or calibration
has changed. The missing categorical ortholog-transfer adapter is also verified:
`parascape/transfer_records.py`, `results/ortholog_transfer_pilot_2026_10_08/`.
Same 560 test genes: 123 calls/437 abstentions, 65 correct/58 wrong; all-hidden
surrogate agreement 0.116071 versus called-only 0.528455. Native method has no
class-score/support output; no scores invented. Both target/source are stored
spatial predictions with context/independence gaps. 464 checks passed and exact
projection/call/card/baseline replay verified (commit `445e5f1`).
`parascape/conformal_records.py` adds kNN conformal sets with verified base model
and separate calibration labels, exact native parity on frozen roles. Shared
cards retain set coverage/size/singleton/empty share/native efficiency and class
metrics. `results/conformal_label_pilot_2026_10_08_v2/`: empirical prediction-grade
coverage 0.973214, mean set size 24.455357/26, efficiency 0.061786, zero singleton
calls. All-training-classes control has coverage 1, size 26, efficiency 0. 908
checks passed; optional pdoc skipped. No biological/exchangeability guarantee.
Continue missing `holdout_search`/`multiplex_modules`, logistic conformal variant,
full outer coverage and independent truth. The former two need explicit label
benchmark-task declarations alongside their native cluster tests; graph/map
transductive construction and learned representation fitting must follow 64.05.

GT-SPATIAL-01 truth review: `results/spatial_truth_source_review_2026_10_08/`
retrieves/verifies seven primary spatial companions, complete marker fields and
Plasmodium Figure 3 outcomes. Toxoplasma 718 final markers include its 62 new IFA
results; those cannot independently test the final classifier. Nine selected
Plasmodium IFA outcomes map exactly and are outside both marker fields; six
other attempted targets remain unknown. Coarse taxonomy, assay/context, source
dependence and target-selection gaps remain. 426 relevant checks passed; exact
complete-row/media/mapping replay. No biological admission or runtime changes.
Continue host GTEx/FANTOM source lineage (65.04) and remaining categorical
adapters independently of unavailable content.

Host expression review **GT-HOST-TX-01**:
`results/host_expression_source_review_2026_10_08/` downloads primary GTEx v11
summary/LCM README and original Atlas `tpmss.tsv`, with verified associations and
provider/acquisition hashes. Atlas parsed values exactly match the archived table.
v11 remains a candidate with no new citation-age clock. Full gene rows preserve
PAR_Y, which native split-at-dot collapses. Current replay has three missing/two
extra GTEx projections and Atlas mapping/cutoff gaps. The documented 0.5 filter
per region reproduces every overlapping legacy value; small mapping/coverage gaps
remain. All rows/gaps/refusal diagnostics retained. 534 relevant checks passed;
counts/values/hashes/policy refusal replayed. 65.04 is 4%; 66.04 is complete.
Continue **GT-HOST-TX-02**, host mapping correction and original gene evidence.
No installed values or calibration changed.

Historical 66.04 first builder milestone (then 25%): `host.uniprot_index` now withholds ambiguous
symbols; Ensembl mappings/schema unchanged. `deposits.k562_rhoptry_evidence`
retains all original source rows before protein projection.
`results/host_symbol_mapping_2026_10_08/`: 20,010 exact gene records, 39 ambiguous
and 1,271 unmapped retained. The original builder compared 18,700 unambiguous
candidate protein rows using pandas' default tolerance; this did not prove
binary identity. 574 relevant checks passed, two existing skips. The later final
migration below completed withdrawal/gene preservation and strict precision
review, preserving pre-existing source/legacy encoding differences explicitly.

## Latest controlled completion: 66.04, 2026-10-08

[66.04](instructions/done/66_04_host_symbol_projections.md) completes instruction 66
(4/4): 39 ambiguous rhoptry protein projections withheld (117 cells), every other
installed human value exactly preserved, all 20,010 original gene records shipped.
Verified withdrawal ledger, exact external backups and executed migration/replay
are in `results/host_symbol_mapping_2026_10_08/migration/`. 575 follow-up checks
passed, two existing skips; two initial failures resolved. Strict precision review
corrects earlier tolerance-based exactness wording and retains pre-existing source/
legacy encoding gaps. Mouse/parasite feature/graph/track-record/claim hashes unchanged.
Total completed actions: 16/48; instruction 64 remains 9/40, categorical records
35%, literature-wide admission/replacement audit 4%. Continue frozen benchmark
adapters; host expression lineage and host gene spaces remain separate open work.

Latest 64.10 partition: native logistic conformal pilot completed and exact native
replayed in `results/conformal_logistic_pilot_2026_10_08/`. Fixed C=0.5/alpha=0.1,
train-only rank/model state; separate 569 calibration/560 test genes. Stored-label
set coverage 0.905357, mean size 6.975/26, zero singleton calls. 1,088 relevant
checks passed. Categorical action now 35%; no new completion tick, independent
biology/full outer coverage/missing adapters remain pending.

64.11 metric prerequisite now 10%: numeric cards retain MAE/RMSE for constant
baselines and small answered cohorts, with matched answered-row baseline skill.
`results/numeric_metric_review_2026_10_08/` contains executed analytic controls.
544 final checks pass; original historical calibration/data remain unchanged.
Next bounded numeric adapter: training-only native ridge on the existing
`fit_invitro_hff` candidate, original direct-experiment grade with unresolved
source/units/admission gaps; freeze scope before expensive fitting.

Latest numeric partition: 64.11 now 25%, native ridge on the existing HFF
fitness candidate, canonical `results/ridge_value_pilot_2026_10_08_v2/`.
Training-only exclusions/ranks/model; 4,015 train/1,106 tune/1,103 calibration/
1,101 test genes, 331 features. All training matrices/native predictions replay
exactly; MAE 1.269871, RMSE 1.590402, Spearman 0.701225, matched training-mean
MAE/MSE skill 0.335439/0.474868, coverage 1. Source/units/mapping/independent
biological admission remain gaps. Constant-float baseline decile diagnostic
preserved; revised cards retain identical values/models and correct code hashes.
1,123 final checks passed. Continue other numeric/categorical adapters and
biological truth; no additional full action complete (still 16/48).

## Active continuation/watchdog, 2026-10-08

The user explicitly requests a watchdog to keep the list moving and log stops.
Native goal continuation is active for unfinished instructions 64/65/66.
`parascape-work-watchdog.timer` is enabled/active, polling once/minute under
128 MB. `.parascape-watchdog/events.jsonl` retains private stop/progress events;
read it and state.json at resume, get actual goal status, record current item and
reason through `scripts/work_watchdog.py record`. Before a turn ends, log its
actual reason; progress reports are automatic-continuation boundaries, not
completion. Four earlier runtime invalid_prompt/content-restriction errors were
observed and logged, alongside eight other historical turn boundaries. Do not
bypass restrictions or invent a cause after an abrupt loss of runtime activity.
403 checks pass; completed [67.01](instructions/done/67_01_work_continuation_watchdog.md).
Total now 17/49 complete; 64.10 stays 35%, 64.11 25%, 65.04 4%. Repost all 49
rows whenever a complete item earns a tick. Details: [watchdog](docs/work-watchdog.md).

Latest numeric interval partition: **64.11 now 35%**, verified train-only ridge
reference plus separate 1,103 calibration genes and all 1,101 test genes in
`results/value_conformal_pilot_2026_10_08/`. Native point predictions, bounds,
absolute-error quantile and row/baseline cards replay exactly. Fixed nominal
coverage 0.9; stored-value coverage 0.921889, half-width 2.762314, finite
availability 1, mean width 5.524628 (2.517286 eligible-truth SD). Explicit
finite/unbounded/unavailable interval statuses and glossary metrics prevent
coverage being mistaken for capacity. Point-only cards gain null interval fields;
historical snapshots retain original values/code. 1,009 relevant checks pass,
10 initial fixture module-reference failures resolved. No biological admission
or exchangeability guarantee; units remain unresolved. Full model state is
verified upstream, not duplicated as a standalone deployment bundle.
No full action completed (17/49); 64.10 remains 35%, 65.04 4%.
Continue controlled remaining numeric adapters/variants and biological truth;
full outer coverage and shared user-facing scorecard component remain open.

Latest host partition: **64.21 human gene space now 30%**. Executed
`results/human_gene_space_foundation_2026_10_08/` builds 58,988 actual canonical
source-gene rows from verified GTEx v10, preserving all 59,033 original source
records and all 118,066 selected source cells exactly. The 45 qualified PAR_Y
records keep their 90 cells separately; they never overwrite ordinary X-gene
measurements. Reviewed UniProt crossrefs retain all 22,320 pairs, with 39,630
unmapped, 19,208 one-protein and 150 multi-protein genes; 1,713 source genes have
a shared protein. Counts/source positions/annotation versions are separate
metadata and never biological features. No protein feature is projected.
Per-column metadata plus typed measurement provenance retain exact source
hashes, median-TPM units, cultured-adult-fibroblast versus HFF context gap, reference
scope and redistribution/benchmark gaps. 618 checks passed; two wrapper source-role
errors before output writes recorded in `results/human_gene_space_preflight_2026_10_08/`.
All inspected installed host/parasite table/graph/track-record hashes unchanged.
No host space registered or distributable pack claimed. Continue controlled
reference/license/graph/independent-opening gates for this source set, then mouse
gene foundation and remaining sources. GTEx v11 still pending 65.04 comparison;
no claim that v10 is the best available replacement. Still 17/49 complete.

Latest Discoveries partition: **68.01 complete; 18/51 full actions complete**.
`parascape/discovery_labels.py` plus the actual Discoveries panel now browse 41
Toxoplasma and 30 Plasmodium labels, with function first. InterPro/Pfam/EC classes,
phenotypes, stages, structural labels and annotation flags no longer require a
generated claim to appear. Search terms/IDs, select overlapping class members,
open genes and existing legacy label/class scorecards. Claim filters/export/map
behaviour preserved; verifier-column rendering and empty-recipe lookup corrected.
Executed `results/discoveries_label_coverage_2026_10_08_v4/` independently replays
every variable and all 21,787/18,879 functional memberships exactly. Source tables,
claims, recipes and 1,038,372 held-out rows have unchanged hashes. **413 final
checks pass, one optional pdoc skip**; earlier broader app/navigation/help checks
pass with one corrected new fixture assertion. Three audit diagnostics are kept:
EC replacement codes mentioned in descriptions are not extra assigned classes,
empty InterPro delimiters and internal EC semicolons preserve source semantics.
No new functional claims or held-out rows: **68.02 is the next priority**, with
declared curation/domain-prediction/orthology lineage, current ontology and
unknown-negative semantics before fitting or reporting biological accuracy.
Instruction 64 categorical/numeric actions remain 35%; human gene space 30%;
dataset promotion/admission 4%. No native-goal completion or release implied.

Latest functional partition: **68.02 now 25%; still 18/51 whole actions complete**.
Official ENZYME/InterPro/Pfam nomenclature is pinned. All 34,983 domain memberships
and original descriptions preserved; 7,926 of 8,043 IDs have current names, 112
InterPro IDs unresolved and five Pfam IDs withdrawn with no reassignment. Current
names do not establish original assignment releases or independently measured activity.
Strict EC resolution gives 1,226/1,050 complete curated-source profiles in the two
parasites, plus 1,055 separate Plasmodium orthology-derived profiles.

Frozen FN-EC-01 seed23/k15/min-share0.3 fits training-only numeric kNN: 679 train,
175 tune,190 calibration,182 test; tuning/calibration unused. 349 inputs explicitly
withhold EC/domain/homology/attention summaries; graph empty. Exact recovery71/182
(39.0%), coverage171/182 (94.0%), correct71/171 among calls (41.5%); majority66/182
(36.3%). All native scores/model/card outputs replay exactly, including11abstentions
and an unseen full profile. This is weak source-profile recovery; **zero independently
admitted activity benchmarks, calibrated probabilities or deployment claims**.
Canonical pilot `results/functional_ec_knn_pilot_2026_10_08_v4/`, identity
`e7cf3b3933dbe4314c842098aa38c53a828765c351a5dd1a765d3ed1da60d461`.
Source reviewv2/domain metadatav4/UI bundlev4 and
`results/functional_discoveries_audit_2026_10_08/` retain executed evidence.
Actual Discoveries offers label/strategy/profile/seven major-class/control cards,
gene outcomes and selected-domain provenance. Both organisms' kinase searches
reach InterPro/Pfam; Plasmodium benchmark availability remains explicitly zero.
Original node/claim/recipe/record hashes unchanged. 653 distinct focused checks
verified (652 broad pass, then19 functional UI/reader checks after correcting
one HTML-tooltip fixture). Earlier Qt tuple lookup and method-docstring failures fixed.

Next bounded function work: protected-group Plasmodium EC recovery; complete
InterPro/Pfam profiles; source admission and independent activity truth; additional
adapters/calibration/deployment. Domain-content and domain-edge families differ,
and search.LAYER_SOURCES omits domain: explicitly exclude domain graph edges and
derived maps/representations with sentinel tests before graph benchmarks. Attention
aliases need a general exclusion follow-up; already explicitly withheld here.
Three worker agents were explicitly authorized by the user and delivered metadata,
UI scorecards and the ordinary precompute framework. Keep concurrent small jobs
within the combined2GB lease and root Qt/app checks serial; defer GPU/>=8GB jobs
during the plaque priority. Two actual worker invalid_prompt stops were logged;
root continued independent tasks without retrying restricted content. Partial
`scripts/review_functional_activity_truth.py`, `results/functional_activity_truth_2026_10_08/`
and `tests/test_functional_contract_audit.py` remain untouched/unreviewed/uncommitted,
not admitted evidence. Private watchdog logs retain actual causes; timer remains active.

Latest precompute partition: **64.25 now35%**, serial typed job graph, complete
parent content IDs, input snapshots, validated/corruption-refusing resume,
interruption recovery, failures/blocked descendants, atomic journal/writer lock
and cooperative attempt/time/currentRSS/disk/package budgets. Executed
`results/precompute_ec_views_2026_10_08/` retains all182 real frozen EC test rows,
recomputes exact profile/seven-member cards, stops after one partition, resumes
only two unfinished builders and replays with zero builders. Stages0.350/0.957/
0.328s; process lifetime peak187.3MiB, typed artifacts166,244bytes, recorded evidence
688,151bytes before addedREADME.83 pure framework/artifact/split checks pass;
34 runner checks also pass within the final broad app run using deterministic
fixtureRSS. Remaining production feature/map/model/deployment builders and adapter
coverage keep this item open. No new tick (18/51). Functional work pushed as682a399;
no original runtime strategy, dataset promotion, native-goal completion or release.

Latest integration partition: **68.02 30%, 64.17 50%; still 18/51 whole actions**.
Three authorized workers delivered complete recorded-domain profiles, the exact
functional coverage matrix and shared scorecard presentation; they continue bounded
profile-target, production-feature and class-card software partitions. Domain census
retains all34,983 memberships, including withdrawn/missing-current identifiers.
Coverage has880 addresses/11installedlabels/oneverifiedrecovery artifact and zero
independent biology/calibrated deployment. Discoveries exposes current-organism
240/320address coverage with exact recorded-result routes. Exact installed-node
hash/table binding prevents imported/altered contexts from borrowing archived accuracy.
Functional source exclusions close registered/declared derived inputs and reject
cycles, unresolved provenance and missing encoded questions. Existing artifacts
are not retroactively rewritten.

Shared immutable cards expose definitions/full record/export, source/split/control
gaps and readable evidence counts; actual human source status retains58,988canonical
genes without host inference/registration. Task-aware set metadata preserves native
algorithm schemas. Profile/major class populations remain distinct. The final
executed UI source replay is recorded in
`results/functional_coverage_ui_2026_10_08_v3/`; earlier prototypes remain immutable.
FN-EC-02 stopped with an actual runtime content-access restriction; do not retry or
reconstruct that task. Restricted partial activity-truth files remain unreviewed.
Watchdog logs the actual stop and ongoing independent work.

GitHubCI for b5b854b failed before regressions because default pytest collected
archived immutable test copies with duplicate module names. Commit043a2ad limits
default discovery to live tests and is verified pushednightly; normal collection
now succeeds without archive paths. Preserve archivedcopies. The replacement CI
run is pending; do not claim full-suite success or release readiness.

Final integration validation: **738 focused checks pass**, package-version check
passes, V3 actual UI/source audit passes, all51manifest-covered outputfiles verified.
Prepared-profile target software has25pure checks; numeric intermediate feature
builder/runner has59pure checks under400MB. The latter pins source/split/exclusion
and imported guardcode, preserves native train ranks/frozen queryECDF, and rejects
complex/attention inputs. No actual production feature pilot yet;64.25stays35%.
Root combined test RSS fixtures are deterministic; actual external limits remain.
Continue a bounded executed feature-operator pilot and domain prepared-target
capacity census, then remaining categorical/numeric/source admission partitions.

Latest bounded audit partition (2026-10-08): three worker jobs are terminal.
Prepared complete-Pfam targets are canonical in
`results/functional_profile_targets_2026_10_08_v3/`; earlier serialization failures
remain preserved. Source/target/split/capacity replay retains all 8,140 genes and
333 unsupported test profiles without fitting. The numeric operator pilot in
`results/precompute_features_2026_10_08/` exactly replays all 349 columns and 1,226
ordered EC entities, including child-only resume and zero-builder replay.
Its 401.6 MiB observed process peak differs from external 400 MiB cgroup accounting;
do not describe process peak as below 400 MiB. The six-task presentation audit in
`results/scorecard_task_views_2026_10_08/` passes 539 software checks, not biology.
Output receipts (20/33/20 files), input/code hashes and three script parsers pass.
**64.25 is now 45% (4–6 h); 64.17 stays 50%; 68.02 stays 30%; still 18/51 complete.**
User now requests minimum tokens and smallest tasks first: integrate ready evidence
and choose bounded existing acceptance gaps before expanding models or sweeps.
Restricted activity-truth/FN-EC-02 work remains unreviewed and must not be retried.


Latest acceptance: **64.17 complete; 19/51 full actions, 10/40 instruction-64 actions**.
Dedicated calibration navigation and explicit unavailable date/version indicators
preserve raw source snapshots and native metric schemas. The canonical desktop
packet is `results/desktop_scorecards_2026_10_08_v2/`: 1,573 checks across 14
actual Qt cards, 40 verified receipts, 0.68 seconds/150.0 MiB process peak.
The first attempt preserves an incidental JSON metric-order diagnostic, corrected
by comparing full metric records by identity with exact values/definitions.
Fresh Discoveries UI/source replay is canonical in
`results/functional_coverage_ui_2026_10_08_v4/`; all 51 nested receipts/input hashes
pass, as do 90 focused view/browser/host/class/Discoveries tests and version check.
This completes reusable presentation only; no biological truth, calibrated
functional deployment, host installation or scientific adapter admission added.
Choose the next small browser/navigation acceptance gap while larger biological
source and model work remains open. Restricted partials stay untouched.

Latest acceptance: **64.18 complete; 20/51 full actions, 11/40 instruction-64 actions**.
Tools → Datasets shares the Slots window. All 162 sources reconcile to 180 scoped
records, with explicit host/entity qualification, original stored values, source
cards, bounded paging and separate session imports. Processed availability is
separate from raw-assay provenance; stored slots do not imply experimental grade.
Canonical `results/dataset_browser_2026_10_08_v8/`: 3,900 notebook checks,
27 current input hashes, 24 verified receipts, 98.05 seconds/1,082.8 MiB peak.
505 focused regressions pass, one existing GL-context skip; six prior CI failures
and native interval checks pass (24). Full CI for old 9d7752b failed; its six
specific regressions are fixed, replacement nightly CI must still run.
Two native Qt test failures remain recorded with unproven causes; after explicit
new slot-test iterator cleanup, the same combined scope passes. A worker's 400 MiB
explanation audit was OOM-killed, then passed serially. Preserve all diagnostics.
No biological truth, strategy fits, host installation or source promotions added.
Next smallest bounded action: 64.19 gene evidence entry, reusing the query,
provenance, scorecard and dataset-browser contracts. Human foundation stays 30%,
mouse 0%; 65.04 stays 4%, 68.02 stays 30%. Restricted partials remain untouched.

Latest acceptance: **64.19 complete; 21/51 full actions, 12/40 instruction-64 actions**.
Exact organism-qualified canonical and recorded alias lookup exposes collisions,
including index targets absent from the current table. Product matches require
an explicit candidate click. Unavailable alias targets clear stale selection.
The gene card opens condensed/all-column evidence grouped by biological question,
with exact original values, source cards and dataset-address navigation. Label
and class routes preserve complete legacy ledger results; changed tables cannot
borrow that accuracy. Replacing the dialog releases its predecessor.
Canonical `results/gene_evidence_2026_10_08_v5/`: 1,246 executed checks over
442 Tg/167 Pf columns, GRA16/AMA1 aliases and two recorded label/class routes;
22 current input hashes, 20 output receipts, 102.68 seconds/804.1 MiB peak.
473 relevant regressions pass with one existing GL-context skip; final widget
width checks pass (4). Initial audit diagnostics preserve a wrong Pf column and
an unscored AMA1 class; no annotation or scorecard was fabricated to pass.
Next small action64.20: labels/classes as entry points, reusing Discoveries,
query/provenance/dataset/gene contracts. Scientific coverage, hosts and dataset
promotions remain open. Preserve restricted partials and nightly/full-CI gates.

Latest acceptance: **64.20 complete; 22/51 full actions, 13/40 instruction-64 actions**.
Tools → Labels and classes opens all 41 Tg/30 Pf native labels, function first,
with membership/unknown counts and separate records for all 39 strategies.
Original complete protein profiles, overlapping classes and known False values
remain intact. Class precision/recall/confusion and label overall/macro metrics
retain exact legacy settings/seed/cohort; missing biology/hierarchy stays explicit.
Actual class→gene→label→source navigation passes for both organisms. Ambiguous
source questions clear stale cards and records. Changed contexts refuse accuracy.
Final immutable V4 Tg/Pf packets: 418+459 checks, 18 input hashes/20 receipts each;
178.64/181.66 s, 1849.4/1373.2 MiB peaks, serial under the same 1900 MiB cap.
494 relevant checks pass; version/whitespace pass. Preserve first categorical
failure, successful V2 and combined V3 OOM diagnostics. No biological truth,
source promotion, strategy fitting, deployment or host installation added.
Latest small partition: 68.02 training-only complete-Pfam baseline controls on
verified prepared-target V3. The exact control partition is accepted: 16 focused checks, 35 input hashes and
22 receipts verify; the executed 400 MiB-capped notebook finishes in 6.676 s.
All 635 test genes/333 unsupported profiles remain; majority/prevalence recover
19/2 profiles. This is annotation-control arithmetic, not biological accuracy.
68.02 stays 30%; classifier fitting/source admission/deployment remain open.
Full goal remains active; nightly CI/release gates and restricted partials remain.

Latest documentary gate: human V5 terms/reference review is accepted, with 25
foundation/source receipts verified twice, five authoritative captures, five
GTEx V10→GENCODE39/GRCh38 checks and 19 output receipts. Conditional GTEx public
redistribution and UniProt CC BY4 documented; future pinned v10 must show older
version/source/date attribution. GENCODE39 is GRCh38.p13. V11 linkage, pack notices,
reference-gene mapping/graph/opening, publication and biological admission remain
open. Human stays30%, mouse0%,65.04 stays4%; no source/foundation promoted.
Preserve V1–V4 diagnostics/captures; actual V5 review111.2MiB/0.44s under400MiB.
Functional attention-alias closure64 checks pushed37f6770. Complete-Pfam native
kNN pilot is frozen at k15/share0.3 with train-only source-closed numeric inputs.
Initial categorical parity failure is NaN/None representation; only missingness
normalization changed,13 synthetic checks pass. Root verifies actual full matrix
and native outputs under800MiB before canonical replay. Two signal9 worker
partials remain diagnostic, without asserting an unproven OOM cause. No new tick.

Latest functional acceptance: canonical complete-Pfam kNN V4 at
`results/functional_profile_knn_2026_10_08_v4/`; fixed k15/share0.3, 348 inputs,
2,381 training/635 test genes, 333 unsupported test profiles. Exact native ranks,
votes/support/scores, refit, serialization and typed artifact replay pass.
12 correct/eight wrong/615 abstentions: all-eligible recovery1.89%, coverage3.15%,
correctness among calls60%; majority control recovers19/635 (2.99%). This weak
annotation recovery is not admitted biology, calibrated confidence or deployment.
47 input/47 output/12 artifact receipts independently verified in the acceptance
notebook;124 focused software checks pass. Runtime46.61s/975.92MiB process peak
serial under1900MiB, all workers terminal. Earlier unit/population metadata refusal
and authoritative800MiB systemd OOM preserved; only categorical missingness,
metadata and retained-copy lifetime changed. No relaxed numeric checks/settings.
Whole actions22/51,64.10 stays35%,68.02 stays30%. Next small partition: generalize
the existing EC-only functional result reader/packager to complete domain profiles,
then expose this exact Pfam benchmark and class/gene cards in Discoveries. Independent
truth admission, other mechanisms, human/mouse packs and dataset promotions remain.

Latest UI integration accepted: Discoveries now ships exact EC+complete-Pfam
benchmarks, with namespace-specific profiles,479 Pfam profile cards,583 domain
membership cards and all635 held-out genes/333 unsupported profiles. Recorded
presence/complement never claims biological negatives. Native control scopes and
full split identity/role counts retained. Candidate/merge/source/UI receipts pass;
118 relevant regressions pass. Actual UI V3 audit10.17s/1,067,978,752byte process
peak under1900MiB, all workers terminal. Preserve V1 authoritative800MiB OOM and
V2 incorrect export-envelope key diagnostic. Shipped bundle externally pinned;
no fit, biological admission or deployment. 68.02 now35%, whole actions22/51.
Next smallest bounded function work: protected-group Plasmodium EC preparation
using existing source/exclusion/controls contracts, preserving missingness and
source-grade gaps. Restricted activity-truth partials remain untouched; native
goal/watchdog active, host spaces and dataset promotions still open.

Latest Pf preparation/control acceptance:
`results/pf_functional_ec_controls_2026_10_08/`. Whole5,720 genes; direct
EC1,050 complete/170 unresolved/4,500 unannotated,19 atomic major-class profiles.
Orthology EC1,055 complete/451 unresolved/4,214 unannotated kept separate.
Seed20261008, native groups, roles600/151/147/152; zero missing/fallback collisions,
four unsupported test genes retained. Train-only majority/prevalence recover43/35
of152; these are control arithmetic, not feature_knn performance or biological
accuracy. All rows/class/major cards and exact source/group/serialization replay
saved. 128 relevant checks pass,40 outputs/38 inputs independently verified.
Runtime4.13s/261,208KiB under400MiB. Corrected pf_enzyme_classification's
unsupported direct='curated' prose, regenerated catalogue; no source bytes changed.
Next bounded step: fixed native Pf EC kNN on this exact protected split, with
full-table registered source-closed training-only feature selection, native replay
and matched controls. Do not reuse projected-control exclusions as sufficient for
new feature columns. 68.02 stays35%,22/51 complete; independent biology, calibrated
deployment, host spaces and dataset promotion remain open. Restricted files untouched.

Latest accepted function partition: native Pf EC kNN canonical V2 under
`results/pf_functional_ec_knn_2026_10_08_v2/`. All152 test genes/four unsupported
retained:40 correct,103 wrong,nine abstentions; majority43/prevalence35, so the
candidate does not beat majority. First packet remains diagnostic because two
registered homology summaries escaped source exclusions. Closure policyv3
withholds ortholog_number/has_paralog and descendants; original split/settings
unchanged. Canonical runtime5.42s/280,400KiB under400MiB. Independent acceptance
V2 verifies49 outputs,66 inputs,13 artifact receipts, original native controls
and installed context; candidate reader1.24s/195,960KiB. 98 focused checks pass
under800MiB; prior resource/path failures preserved with original code.
Actual biological accuracy remains unknown. Next small partition: install and
verify the accepted Pf candidate in Discoveries. Shipped Tg EC/Pfam bundle still
unchanged; function35%, whole actions22/51, goal/watchdog active.

Latest integration: **68.02 now40%;22/51 whole actions**. Shipped bundle
SHAb5fa39aa80200e19560e4600502b6decee1a2ec84ba38c03e33b11f8c1fe6eb3
contains original Tg EC/Pfam entries plus accepted Pf direct EC test, all exact.
Actual Pf desktop audit verifies original controls, membership/profile exports,
gene/coverage routes and altered-source/wrong-organism refusals; screenshot
inspected.81 focused checks pass; independent acceptance verifies18 outputs and
15 current input hashes. Audit5.03s/717,492,224bytes; merge2.20s/519,588KiB under
serial1900MiB with all workers terminal. Initial installer400MiB OOM preserved
with original bundle/code; no fitting or biological/calibration/deployment
admission. Guide updated. Continue remaining functional source/method partitions,
then other instruction64/65/66 gaps; retain cancelled/restricted-task boundaries.

Latest native functional strategy: accepted serial random forest on unchanged
Pf EC target/cohort.59/152correct,93wrong,noabstentions; kNN40/152,majority43/152,
prevalence35/152; all four unsupported test genes retained. Native defaults
300trees/leaf2; only n_jobs4→1 for deterministic probability addition. Exact
native/refit/treeJSON/parquet/artifact replay; root serial packet V2 succeeds
11.44s/570,120KiB under1900MiB. Worker signal9 preserved, cause unproven.
52 outputs/70inputs/12artifact receipts independently verified. Shipped bundle
SHA7f917b83f211c7d619b90aefb2d6cf01d75efdded703f1b71cbbcfd5934c2134
now2Tg/2Pf; previous entries/native controls exact. Actual two-strategy desktop
routes and both screenshots verified.101 combined checks pass; prior Qt styling
segfault/84isolated passes retained, cause unresolved. Final receipt audit112outputs
and15current UIinputs. Source recovery only, no biological/calibration/deployment
or importance admission. Function40%,categorical35%,whole22/51. Next small task:
descriptive paired method agreement/error dependence on same cohort, shared-source
limits explicit; other host/literature/strategy gaps remain. Watchdog/goal active.

Latest paired partition:64.32/64.33 now20%, whole22/51 unchanged. Canonical
`results/pf_functional_agreement_2026_10_08_v2/`, report SHA
91e32f04ab83bfdab2f6e115be8456554b74e0f3cc5653143fc5aa1789253e89.
152genes/four unsupported/145groups retained,143joint calls:94agreements,
33correct/61same wrong;49conflicts (7kNN correct/25forest correct/17bothwrong).
Nine forest-only calls (one correct),78bothwrong overall. Shared labels/training/
exact native input/retained columns and dependency hashes disclosed; no independence,
confidence, biological accuracy, ensemble-selection or group-aware interval claim.
Shared scorecard HTML/JSON export, no new app button.61checks pass serial1900MiB,
workers terminal. Report1.186s under800MiB;101source outputs/25artifact receipts;
independent acceptanceV2 verifies13packet outputs/14current inputs and every count/rate.
Preserve800MiB test OOM/renderer diagnostics. OriginalV1 retained; V2 adds strict
fit-role/feature-access/admission guards, same values. Next smallest step: expose
this card from existing Pf Discoveries tests with context refusal and original
method cards preserved. Broader uncertainty/calibration/method/host/dataset gaps
remain. Restricted activity-truth partials untouched; native goal/watchdog active.

Latest paired desktop partition:64.33 now30%,64.32 stays20%,whole22/51 unchanged.
Discoveries → Functional tests → Comparison: kNN and random forest opens the exact
pinned report from either Pf strategy. All152 genes/four unsupported remain, both
calls visible with named-method outcomes; rates show fractions/percentages, exports
retain exact values. Gene click/paired outcome popup and return to original method
cards/controls verified. Missing/corrupt reports or changed current source/organism
clear stale comparison cards/rows/exports. Shipped report SHA940da900…8e98b matches
canonical bytes. Every functional selection rechecks context, including return
to a method after source mutation. Canonical desktop V3 packet under
`results/pf_functional_agreement_ui_2026_10_08_v3/`; actual7.43s/810,123,264-byte
process peak under serial1900MiB, workers terminal. Final68 focused checks pass;
Tg desktop regression also passes. Independent acceptanceV2 verifies54 output
receipts/18 current inputs. Preserve initial deleted-dialog
test teardown and acceptance-script import failure; both fixed, original code/logs
retained. No original nodes/claims/benchmark bundle, fits, calibration or biological
admission changed. Next smallest partition: group-aware uncertainty for this
matched source-profile pair, explicit denominators/unsupported genes retained.
Other strategy/host/literature gaps remain; restricted partials untouched.

Latest group-uncertainty partition:64.32 now30%,64.33 stays30%,whole22/51 unchanged.
`parascape/group_ratios.py` implements shared whole-group ratio resampling,500draws,
seed20261008,minimum5groups,linear2.5/97.5 percentiles,explicit undefined draws and
exact saved-selection replay. Canonical V2 packet
`results/functional_group_uncertainty_2026_10_08_v2/`:all152genes/four unsupported/
145recorded groups retained. Agreement recovery35.1%, descriptive interval25.0–45.8%;
bothwrong54.5%,46.1–63.0%; right-minus-left recovery12.5percentage points,5.3–20.5.
12focused checks pass under800MiB; independent audit expands all500draws to original
rows and reproduces every count/ratio/six intervals exactly,10outputs/six current
inputs verified. Shared scorecard export removes the obsolete unavailable-interval
statement while keeping biological independence/calibration unresolved. Initial11
pure checks retained; no failures, fits, source changes or deployment. These are
separate exported intervals; existing app report is unchanged. Next smallest task:
attach the pinned interval result to the existing paired Discoveries card, preserving
exact original rates/source-context refusals. Other method/host/literature gaps remain.

Latest interval integration:64.33 now40%,64.32 stays30%,whole22/51 unchanged.
Existing Pf comparison shows five exact rates with95% descriptive whole-group
percentile bounds and forest-minus-kNN recovery12.5points (5.3–20.5). Compact text
shows500resamples/145groups; exchangeability/source-independence assumptions are
explicit. Shipped interval JSON SHA0bc3e966…061f9 is identical to the numerical
packet; parent report SHA940da900…8e98b and benchmark bundle remain unchanged.
Missing/corrupt intervals retain valid point comparisons; incompatible parent/
source contexts refuse data. All152genes/four unsupported, gene outcomes and exact
exports preserved; no per-gene probabilities.48 focused checks pass;18 final
presentation checks pass. Canonical desktopV2 under
`results/functional_group_uncertainty_ui_2026_10_08_v2/`,6.88s/864,231,424byte
process peak serial1900MiB, workers terminal. Independent48outputs/20current
inputs verified. No failures, fits or biological/calibration/source admission.
Next smallest bounded step: class/profile comparison navigation on this paired
cohort with explicit annotation-negative semantics, then other method/host/literature
gaps. Restricted partials remain untouched; native goal/watchdog active.

Latest class navigation:64.33 now50%,whole22/51 unchanged. Both Pf methods expose
seven major-class/twelve complete-profile comparisons. Full152genes/four unsupported
remain in metrics/exports; displayed union includes false positives. Nonmembership
is source-recorded, not biological absence; abstentions unknown. Native cards/raw
profiles retained; class-specific intervals/calibration/consensus controls unavailable.
19 focused+9 final checks pass. Desktop9.49s/948936704byte peak serial1900MiB;
`results/functional_pair_class_ui_2026_10_08/`. Independent acceptanceV2 exactly
recounts all19 addresses and verifies21 outputs/21 inputs/unchanged frozen bytes.
Preserve initial missing-fixture-metadata failures and compact/expanded overlap
acceptance assertion. Next smallest work: remaining bounded class-comparison
uncertainty or an independent source/host gap; avoid widening into unvalidated fits.
Other method/host/literature work remains. Restricted activity-truth partials untouched.
