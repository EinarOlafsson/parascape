# Browse datasets

Open **Tools → Datasets**. The browser shares the catalog window with **Slots**.
Filter by organism or host, entity unit, evidence family, measured context and
availability. Search also matches source names, columns, labels and recorded gaps.

Selecting a source opens its coverage and provenance card alongside its stored
values. **Previous** and **Next** reach every original row in pages of 200.
Hover over a shortened cell or heading to read its full text. Double-click a gene
identifier to open that gene in the corresponding organism's map. Host proteins,
metabolites and pair records retain their own identifiers and units.

Coverage counts nonmissing stored values in the declared table. Zero and `False`
are preserved; missing values retain unknown interpretation. Pair counts retain
an unknown assay denominator where none was recorded. Source cards show declared
derivations, publication links and supplied lineage; missing units, mapping,
source dates, versions and calibration are explicit. Biological accuracy requires
its own benchmark.

Session imports appear as separate sources with their preprocessing record. The
original source-file identity remains unavailable when the import record lacks it.
Importing values refreshes an open browser without changing installed cache files.

For Python callers, `parascape.dataset_space.DatasetSpace` accepts an explicitly
qualified inventory and local tables. `filter_rows`, `facets`, `source_card` and
`entities(row, start=0, limit=200)` expose the same scopes, immutable evidence cards
and original values used by the desktop browser.

Acceptance is recorded by `python -m scripts.audit_dataset_browser --out NEW_DIRECTORY`.
The executed notebook freezes local input/code hashes and checks source cards,
filters, exact stored values, pagination, session imports and gene navigation.
Use a fresh directory for each replay.

Gene search resolves exact recorded aliases within the selected organism. Multiple
alias or text matches show a candidate list; select an accession to open its gene.
The gene card’s **All evidence by biological question** link opens a condensed
view. **Show all evidence** includes every stored column and missing value.
Each row shows its question, context, quantity unit, declared origin and gaps;
source selection opens the original coverage/provenance card. Label and class
buttons reach existing recorded results when the current table matches the
installed context. Changed tables require their own validation.

Python callers can use `gene_evidence.build_resolver` for exact scoped aliases
and `gene_evidence.evidence_rows` for original per-column records. Gene entry
acceptance runs with `python -m scripts.audit_gene_evidence --out NEW_DIRECTORY`.


## Original protein identifier lookup

In **Compare → Other measurements**, use the source filter to find a study,
source file, assay or quantity label. Matching source metadata is case-insensitive;
feature identifiers retain their exact case. The source count shows how many
analyses match. Refresh keeps the filter and the selected source when it is still
available. Changing the filter, selected source, identifier class or feature query
clears earlier values and exports so their source context cannot be confused.

Select a protein quantity table, then choose
an identifier class above the query box. Enter one exact identifier per line and
select **Show original measurements**. The default searches the complete original feature label.
The other choices search original semicolon-separated protein-group components,
submitted gene labels, explicit reference gene declarations, or reference GN names.
These classes remain separate. Matching preserves case, version/isoform suffixes,
repeated source rows and ambiguous associations.

Results retain each whole original protein-group quantity, its missing values,
source tokens, units and row/column metadata. Multiple matching declarations
annotate the same measurement; they do not duplicate its value. Hover over cells
or export CSV to inspect all supporting declarations, reference headers, source
hashes and original row positions. A submitted label or same-project reference
candidate does not establish a unique gene abundance, the database actually used
for quantification, or identity with a current RNA identifier. No match means the
selected evidence class has no matching declaration, not that expression is absent.

The automatic identifier index includes protein quantity tables from each admitted
measurement namespace, including both reporter-channel families. Its current
601,605 evidence occurrences cover 129 tables; repeated components/declarations are
preserved and are not unique genes or independent observations. The reporter
addition supplies exact component lookup without introducing new reference gene
or submitted-name claims. Each collection retains the selected source inventory,
original catalogue hashes and exact measurement-manifest identities.

Identifier evidence is installed separately under
`processing_audit/protein_identifier_evidence` alongside the original measurement
packs. `PARASCAPE_AUDIT` selects an audit root; `PARASCAPE_PROTEIN_IDENTIFIERS`
can select this particular collection. If it is absent, complete original-label
lookup remains available. This addition follows the frozen 0.0.5.0 snapshot and
requires a later data release for portable installation.

Python callers can use `protein_identifier_evidence.records(pack, queries,
namespace=...)` for every supporting association and `profiles(...)` for original
measurements with evidence attached. Collection receipts bind both the source
measurement manifest and the association files. Gene declarations use only the
literal `gene`/`gene_id` header fields; UniProt GN values stay name candidates,
following the [official FASTA header format](https://rest.uniprot.org/help/fasta-headers).


### Reporter-channel quantities

Reporter tables appear in **Compare → Other measurements** as separate
**Reporter intensity** and **Reporter intensity corrected** analyses. Select the
source family, query a complete original protein-group label, and inspect the
source columns or export measurements/QC as usual. The source's corrected values
are preserved as deposited; ParaScape does not apply that correction again.

Each column retains its original numbered channel heading and ordinal in sample
metadata. The ordinal alone does not establish an isotopic label, biological
sample, control or replicate. Reporter count fields and aggregate intensity remain
original row annotations. No averaging, channel pooling or normalization is
applied. The first locally admitted source has 4,940 rows and six channels in each
family; their shared rows do not constitute independent biological observations.
The automatic importer revisits acquired tables every five minutes, preserving
immutable packs, full-cell reconciliation, descriptive QC and operation logs in
`deposited_values/protein_reporter_channels`. The corrected and uncorrected reporter collections are included in portable
data release 0.0.6.0.


### Aggregate, labelled-channel and reported-ratio fields

Additional protein tables expose exact source-labelled scalar fields such as
**Intensity**, **iBAQ**, **Intensity L/H**, **Ratio H/L** and
**Ratio H/L normalized**. Each appears as its own analysis in Other measurements,
with one original quantity column. These are source quantities, not inferred
independent samples. The column metadata states whether its label denotes an
aggregate/unspecified scope, a labelled channel or a reported ratio; experimental
sample assignments, aggregation scope and ratio orientation still require source
review. A normalized source field stays as deposited without normalization again.

The first admission reconciles 62 tables from 29 source files, preserving all
85,899 cells, including 8,662 missing cells. Original protein groups, source flags,
per-sample fields, peptide/count fields and variability annotations remain intact.
The importer runs every five minutes under the independent
`deposited_values/protein_reported_quantities` namespace. These values are not
approved for direct whole-cell RNA/protein concordance. These 62 tables are
included in portable data release 0.0.6.0. The prepared exact identifier index includes these tables while
keeping every ambiguous declaration and whole-group quantity separate.


### Original source flags and view treatment records

Other measurements shows exact original flags beside each queried value. The
optional **Exclude rows with source "+" flags** checkbox applies only to deposited
protein-group quantities and literal `+` values in `Reverse`, `Contaminant`,
`Potential contaminant` or `Only identified by site`. It starts unchecked. Unknown
values, whitespace variants and missing declarations are retained as unreviewed;
a blank or missing flag is not a biological quality score or proof of a clean row.
Other assays do not inherit protein flag meanings.

Exclusions use the original measurement manifest and row ordinal. Repeated labels
and different sources therefore keep their own flags. Original values, units,
missing cells and files stay unchanged; the selected view can omit marked rows.
A query whose matches are all excluded is reported as filtered, not absent
expression. **Export column technical QC** still describes every original source
row and is not restricted by this view filter.

**Export displayed measurements** writes the retained original values plus a
separate uniquely named `.treatment-….json` file beside the CSV. The record pins
the CSV digest, source manifest/matrix/token hashes, exact query and identifier
context, policy/version and excluded original rows. Existing treatment records
are preserved when a CSV is exported again. **Export view treatment record** also
works when every query match was excluded or no exact match was found. Changing
the source, query or rule clears prior results and treatment state.

This rule uses source declarations. It does not calibrate biological accuracy,
validate sample identity, infer independent replicates or replace a held-out
quality benchmark. The source flag policy is available to Python callers through
`parascape.source_flags.apply`, which returns the annotated view and its operation
record while retaining unfamiliar or missing flags.


### Original accession abundance groups

The local acquisition loop also admits supported accession tables under
`deposited_values/protein_accession_values`. Search **Other measurements** for
`PXD041907`, select **Proteins.txt :: Abundances (Grouped)**, and enter exact
original `Accession` labels. This first admission contains 2,496 source rows, six
source-labelled groups and 14,976 original cells, including 1,550 missing cells.
It was prepared after the frozen 0.0.6.0 snapshot and is not in that public asset.

Every original UID, `Master`, FDR confidence, q-value and per-row presence
declaration remains available in row metadata and CSV exports. No master-only
selection or confidence cutoff is applied. Repeated accession labels retain
separate original row identities; accessions are not assumed to identify unique
proteins, protein groups or current genes.

**Source columns** shows exact same-file group-label matches between abundance
headers such as `Abundances (Grouped): F1` and presence headers such as
`Found in Sample: [S1] F1: Sample`. These are source-header candidates, not verified
sample, control or replicate assignments. Presence declarations never replace
missing quantitative cells. The original protein-group component class and literal-plus
group flag rules do not interpret this different assay.

The enabled five-minute importer preserves original tokens, binary64 values,
missing cells and all row annotations. Each immutable pack records the source,
processor and dependency hashes, software versions, full-cell reconciliation,
unfiltered technical QC and operation history. These quantities retain their
source grouping and scale; cross-study scale compatibility remains unreviewed.


### Source-labelled log2 protein quantities

The local `deposited_values/protein_log2_quantities` collection preserves exact
`PG.ProteinAccessions` rows and numbered filename columns ending in
`.PG.Log2Quantity`. Search Other measurements for `PXD084269` and choose the
`BGS_Analysis_Grid_View_Report-Pivot.tsv :: PG.Log2Quantity` analysis. This first
table contains 7,576 rows, five source columns and 37,880 original cells, including
1,909 missing cells. It was admitted after frozen data release 0.0.6.0.

Values stay on the source-labelled log2 scale. No exponentiation, normalization,
averaging, row selection or imputation is applied. Original candidate,
completeness, variability and average-log2 fields remain annotations. A source
log2 zero is retained as a numeric value; it does not establish absence of
expression. Repeated accession labels and repeated filenames keep their original
rows and numbered column identities. The five-minute importer retains source and
processor hashes, literal tokens, independent full-cell reconciliation,
unfiltered technical QC and operation records.

The source-column association loop now discovers protein quantities across all
active namespaces. It pins catalogue bytes and manifests, refuses changes during
reconciliation, and preserves every earlier association report. Exact numbered
`PG.Log2Quantity` filename headers can link to primary file records and original
submitted design rows without trimming filenames or changing case. The first
log2 table adds five exact primary-file matches and five submitted row
occurrences. Source columns and CSV metadata retain all matching occurrences.
These links expose the deposited metadata; they do not establish control pairing
or biological independence. Other source columns retain explicit unresolved
roles. Both new local collections await a later portable data release.


### Row-wise protein quantities and repeated source records

The local `deposited_values/protein_rowwise_quantities` collection admits the
explicit row-wise protein report containing `PG.Quantity`, `PG.ProteinAccessions`
and original `R.Condition`, `R.FileName`, `R.Replicate` fields. Search Other
measurements for `PXD084269` and choose `ProteinQuant-Normal.tsv :: PG.Quantity`.
Each original row retains its own ordinal and annotations. In the nightly browser,
select a measurement and choose **Inspect selected measurement…** to inspect its
original run, condition, replicate and quality declarations. You can also hover
the measurement or export CSV to inspect the full `source_row_metadata`. These source declarations
do not establish independently verified biological sample identities.

The first table preserves 65,535 original quantity cells, including 29,564 missing
cells, with five distinct filename literals. Its single quantity column represents
a row-wise report, not one biological sample. No pivot creates additional cells,
and no repeated row is averaged or removed. Source scale and absolute units remain
unreviewed; the importer does not infer a logarithmic or linear scale.

There are 24,820 repeated exact accession/filename keys. None has identical values
in every original field: descriptions, names and quality declarations can differ.
Of these keys, 24,041 contain both a reported quantity and a missing quantity;
779 have only missing quantities. Their repeated rows remain distinct. These
patterns do not establish duplicate observations, independent replicates, or
permission to fill missing cells. A row count of 65,535 alone does not establish
source completeness or prove truncation.

The enabled five-minute importer accepts this explicit protein layout while
keeping peptide, modification-site and fragment reports separate. It preserves
source and implementation hashes, original tokens, all row annotations, repeated
key counts, independent full-cell reconciliation, unfiltered technical QC and
operation logs. All earlier 138 quantity tables remain unchanged. This local
admission is outside frozen portable data release 0.0.6.0.

The recurring repeated-key audit also reparses every admitted row-wise source and
checks its original annotations. Immutable, source-bound reports distinguish
all-missing, mixed observed/missing, and all-observed key groups, flag differing
original fields, and retain operation history under
`processing_audit/rowwise_repeated_key_qc`. These are descriptive ambiguity checks,
not biological quality scores or permission to remove or combine rows.


### Inspect one original measurement

The nightly browser offers **Inspect selected measurement…**, also opened by
double-clicking a measurement. It shows that exact source row and quantity column,
the numeric value beside its original token, original row annotations, column
metadata and submitted factor occurrences. Original row fields are expanded when
the view opens. Qualified identifier queries retain
every matching evidence occurrence. Source declarations and missing values remain
distinct from inferred sample roles or biological quality.

Sorting the table preserves measurement identity even when feature labels repeat.
Changing the selection, query, source, identifier mode or flag rule closes earlier
details. Inspection does not change values, the current treatment record, source
files or exports. Missing row annotations are explicitly unavailable rather than
filled from another row. This view is a source follow-up after public 0.0.6.0 and
awaits a later software release.


### Row-specific filename associations

Row-wise reports can put filename declarations on each original row rather than
in quantity-column headers. The nightly browser attaches separately reconciled
`row_design_metadata` to these measurements and includes the association receipt
and reader implementation hash in CSV treatment records. Inspect selected
measurement → **Row-specific filename links and review-only candidates** shows
the exact source label and all retained associations. Submitted design fields
appear as ordered header/value pairs, retaining duplicate headings and repeated
source-row occurrences.

The first row-wise source, PXD084269, has five `R.FileName` literals without the
`.raw` extension present in archive filenames. There are **zero exact filename
matches** and five separately labelled filename-stem candidates, each with a
submitted design-row candidate. The candidate rule appends one literal terminal
`.raw` and performs no trimming, case change, basename extraction or matching of
other extensions. Empty filename declarations never produce links. These remain
**review required**; no input-table linkage, sample equivalence, control pairing
or biological replicate independence is inferred. File metadata associations do
not prove raw payloads were downloaded.

The enabled five-minute reconciliation loop independently reparses all original
row annotations and submitted designs, pins their source/manifest/metadata hashes,
and retains every original row membership. It stores immutable reports and
operation history in `processing_audit/protein_row_design_links`. Missing or stale
association collections stay explicitly unavailable. All 139 earlier quantity
pack outputs remain unchanged. This collection and reader await a future portable
data/software release after 0.0.6.0.

Python callers can use `parascape.row_design_links.enrich(pack, frame,
audit_root=...)` with a frame from `deposited_values.profiles`. It returns the
original measurements with separately encoded row associations and a provenance
record. It does not convert scales, fill missing cells, merge repeated rows or
assign biological quality scores.


### Literal accession-cell components

The nightly identifier selector has a separate **Original accession-cell
component** class for accession-based assays. It matches exact semicolon-separated
components of original `Accession` or `PG.ProteinAccessions` cells. Source case,
whitespace, leading zeros, suffixes and repeated component occurrences remain
unchanged. Matching multiple components annotates the same original row rather
than creating extra measurements. Whole original row quantities and all ambiguous
cell membership remain intact; components do not establish unique proteins,
protein groups, gene equivalence or independent biological observations.

The automatic five-minute index currently covers three accession tables: 75,607
original source-row occurrences and 76,260 component occurrences. It stores the
original field, full cell and component position under
`processing_audit/accession_identifier_evidence`. The existing 129-table
protein-group index stays separate, with its output files unchanged. The browser
enables component classes for their declared assay and rejects applying the group
component class to accession cells or the accession class to protein groups.

Explicit submitted gene fields and same-project FASTA declarations remain
separately qualified namespaces when available. The two newly covered projects
have no acquired same-project FASTA references and currently contribute no such
reference gene/name mappings. The browser reports **no indexed declarations in
the selected identifier class** when the class has no evidence. This is not
evidence of absent expression. This condition, query class, evidence occurrences
and source receipts remain in treatment records even when no row matches.

Python callers use `protein_identifier_evidence.profiles(pack, queries,
namespace='original_accession_cell_component')`. Assay-specific roots are selected
automatically; `PARASCAPE_ACCESSION_IDENTIFIERS` can override the accession index,
while `PARASCAPE_PROTEIN_IDENTIFIERS` continues to select the separate group index.
The new index and browser behavior await a later portable release after 0.0.6.0.
