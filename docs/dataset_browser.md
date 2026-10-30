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
`deposited_values/protein_reporter_channels`. This collection is outside the frozen
0.0.5.0 snapshot and awaits a later portable data release.


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
approved for direct whole-cell RNA/protein concordance and are outside the frozen
0.0.5.0 snapshot. The prepared exact identifier index includes these tables while
keeping every ambiguous declaration and whole-group quantity separate.
