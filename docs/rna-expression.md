# Compare RNA-seq expression

The GitHub data download includes an offline snapshot of verified study/run entries,
selected from the ongoing processing queue. Use `snapshot()` below for the exact
coverage and provenance of your installed data release. This is a
partial release; a study may not yet contain all replicates or its control group.
All eight runs of PRJNA1008424 are included.

Open **compare → Gene expression**, enter one or more exact Toxoplasma gene accessions,
then select **Compare genes**. Filter by study, choose Salmon or kallisto, and export
the selected values with their sample metadata. Genes absent from the selected
measurements are reported explicitly. A measured zero is retained as zero.

**Condition comparison** provides the existing within-study group comparison.
Check the sample titles and assign the control and treatment groups yourself;
the guessed condition labels are only a convenience. A study lacking adequate
replicates or controls cannot support the intended contrast.

Selections must contain distinct runs in disjoint groups. A requested run that
is absent from the matrix is an error. Known repeated sequencing of a submitted
sample cannot count as extra replicates. Unknown sample identities remain marked
unknown; distinct submitted accessions do not prove biological independence.
Changing a study, normalization or group selection clears the old comparison.

Saving a comparison writes the CSV/TSV and a matching `.json` treatment record.
The record retains groups, submitted sample metadata, unknown identities, ordered
model-input hashes, filtering thresholds and library sizes, normalization factors,
software/source hashes, and the exported table checksum. The displayed call
thresholds are recorded separately from the model fit. Exports use fresh filenames.
The API accepts `sample_metadata=` in `expression_contrast.contrast()` and returns
the record in `result.processing_record`. Neither the statistical result nor its
record establishes a biological quality score or verified control pairing.

Open ParaScape once and accept the data download, or run `parascape-data --download`.
After installation the command line works without a display or NAS:

```bash
parascape-expression GENE_ID ANOTHER_GENE_ID --output expression.csv
parascape-expression GENE_ID --study STUDY_ACCESSION --quantifier kallisto
```

Replace the placeholders with actual reference accessions from your study.
The Python interface returns a pandas table:

```python
from parascape.rna_profiles import profiles, snapshot
print(snapshot())
values = profiles(["GENE_ID", "ANOTHER_GENE_ID"])
```

## Measurements and limitations

In **compare → Other measurements**, filter a source transcript table, choose
**Reported parent gene ID → original transcript rows**, and enter the exact
`gene_id` reported by that source. This returns every matching original transcript
row, including repeated labels. Reported parent IDs can differ between reference
annotations. The query does not remove transcript suffixes, translate gene versions
or add transcripts together. Exports retain each transcript's original values,
source annotations and the parent-ID query record.

The bundled reference is **ME49_host**: ME49 parasite transcripts compete with human
and mouse transcripts; the parasite genome supplies Salmon decoys. Both Salmon
and kallisto results are included. All already-verified ME49, GT1 and VEG alternative-reference quantifications for
these completed runs are included in the GitHub pack. Availability varies by study.
Alternative strain builds retain their own exact reference gene identifiers.
Set `PARASCAPE_OMICS` before starting ParaScape to select another matrix pack.

Transcript estimated counts and TPM are summed by the exact reference transcript-to-gene
mapping. Host genes are excluded from parasite tables; the host count fraction is
retained per sample, and parasite TPM is reclosed to one million. The GitHub data pack includes counts, TPM and effective lengths for both quantifiers
for every completed run. Lengths are transcript-TPM weighted, with a mean-length
fallback for unexpressed genes.
Matrices use float32 storage; full-precision quantification sources remain in the
processing archive. This storage rounding does not supply additional precision.

Cross-study TPM is descriptive expression, **not a batch-corrected fold change**.
Protocols, library composition, sample quality and biological context differ. These
tables do not establish gene function or independent biological validation. They
are separate from the prediction features, maps and existing calibration records.

`parascape/data/rnaseq/snapshot.json` records the frozen plan hash, per-study planned
and completed runs, quantifier versions, reference source/checksum lineage, input
hashes, matrix hashes, transformation and precision. Sample sheets preserve source
accessions, sample titles and original run receipts. Cite the original studies;
their titles and accessions are also listed in the bundled omics catalogue.

## Rebuild from the processing archive

```bash
python scripts/build_rna_snapshot.py --root /path/to/omics \
  --verified /path/to/automation/verified.json \
  --plan /path/to/automation/plan.jsonl --out /path/to/new-snapshot
```

Only verified entries from the frozen plan are admitted. Every selected compressed
quantification is checked against its audit hash and read fully, including gzip CRC
validation. Unmapped transcripts or changing run receipts abort the build. The
manifest is written only after all output files have been built and hashed.


## Study design and repeats

The **Studies & samples** tab includes every planned run, marking verified versus
pending measurements. Browse study → condition → submitted sample → sequencing run;
search by study, control label, replicate label, sample or run and export the design.
Source ENA sample attributes are cached for completed runs and the prioritized Trx1
study. Other pending samples retain run metadata and explicitly lack enriched attributes.

Grouping preserves strain, treatment, genotype, stage and other submitted factors.
Title-derived group labels remain suggestions requiring review. A source-labelled
control is distinguished from a candidate baseline; matching replicate numbers do not
prove paired controls. Repeated runs of one BioSample are linked, and the condition
comparison refuses to count them as independent replicates. Shared raw runs across
studies are explicitly labelled. Biological versus technical replicate type is left
unspecified unless the source declares it.

## Original deposited expression and catalogue coverage

The **Other measurements** tab browses independently reconciled source analyses
outside the gene-expression matrix format. Choose an analysis and enter exact
original feature labels, one per line. Every matching source row remains separate;
the table and CSV export include its row number, source-column label, binary64 value,
literal source token, units and provenance. Source feature annotations are retained
without inferring gene equivalence. **Export column technical QC…** uses all source
rows, independently of the query, and keeps missingness, measured zeros, signed
values and uncalibrated biological quality distinct. Local collections use
`PARASCAPE_DEPOSITED_VALUES`; portable data uses `deposited_values`.

The first Workbench collection contains seven tables and 57,468 source cells,
including 21,939 missing cells. Only empty strings or JSON null are converted to
missing numeric values, and their literal tokens are retained. Independent JSON
parsers check every token and binary64 value without tolerance. QC has 240
sample/analysis records; repeated analyses of the same submitted samples do not
create independent replicates. One additional analysis has metadata but no embedded
measurement table. No scaling, imputation or cross-study harmonization is added.

The local protein-group collection adds 63 quantity tables from 26 source files
across 14 projects: 26 Intensity, 20 LFQ intensity and 17 iBAQ tables. All 928,590
source cells were independently reconciled, including 2,454 missing cells. Literal
empty strings and `NaN` become numeric missingness; the original strings remain
available. Intensity, LFQ and iBAQ stay in separate tables. Exact `Protein IDs`
cells, repeated protein-group rows and every original annotation/quality flag are
retained. No reverse, contaminant, site-only or confidence filtering is applied.
Flags appear in row tooltips and the CSV `source_row_metadata` field.

The 748 submitted quantity-column occurrences are not inferred biological samples
or independent replicates. Zero may encode nonquantification; numeric coverage is
not biological detection coverage. Three additional protein-group sources have no
unambiguous supported quantity columns and remain recorded for review. Other
proteomics formats remain pending. Sample/channel/fraction roles, enrichment
context, source confidence thresholds, gene equivalence and compatible scales
must be reviewed before pooling or RNA/protein concordance. Native preparation
and independent validation run every five minutes, with input/output hashes,
software versions, exact adapter source and operation records retained. Independent
source namespaces preserve the existing Workbench catalogue. These collections are
included in the 0.0.5.0 external data snapshot.

The protein importer also accepts the explicit wide layout `PG.ProteinGroups`
with columns of the form `[number] original-file.PG.Quantity`. It retains the
complete original headers and does not reinterpret quality/count/statistic fields
as abundance. The two current tables add 135,116 cells, including 5,075 missing
cells, over 23 quantity columns. The combined local protein collection now has
65 tables / 28 source files / 16 projects, 1,063,706 cells and 7,529 missing cells.
All 63 earlier tables retain byte-identical numeric data, literal tokens and
column metadata, with unchanged original row annotations. The native validator
checks the original feature-column identity as well as every original cell.
Six source layouts now remain explicit review gaps, including long tables and
derived-statistic tables. Filename-bearing headers do not by themselves establish
sample-design pairing; exact submitted factor associations remain a separate step.

Use **Source columns and submitted factors…** in Other measurements to inspect
column context without entering a feature label. The dialog retains the original
column record and separately shows filename associations, complete source rows and
provenance. Numbered fields preserve duplicate SDRF headings. The original column
metadata is not overwritten by later associations; value CSV exports retain both
under `sample_metadata`, with associations in `source_design_metadata`.

The first association collection covers all 65 protein tables / 771 quantity
columns. Of the 23 explicit filename-bearing `PG.Quantity` columns, 11 match a
primary file identity exactly and 12 remain unresolved. Five columns also match
five submitted SDRF row occurrences. The other 748 quantity columns have no
reviewed filename role; similar experiment names are not substituted for filenames.
Repeated rows, channels, fractions and ambiguous primary candidates remain
separate. Exact filename equality does not establish independent biological
samples, paired controls or downloaded raw payloads. A native five-minute task
independently reparses the original SDRF cells and checks every filename candidate;
immutable reports bind each association to the exact measurement manifest and
original source column. Missing or changed contexts remain explicit gaps.

For original Intensity, LFQ intensity and iBAQ columns, the same dialog also shows
**Experiment row candidate** records when the literal quantity-column suffix equals
an explicit `Experiment` cell in an acquired `Name` / `Fraction` / `Experiment`
design table from the same declared project. Leading zeros, case, every matching
table and every fraction row remain unchanged. These are label candidates: they
do not establish which design table was used for a particular quantification,
paired controls or independent biological samples. Candidate metadata is separate
from filename-based SDRF associations and is retained separately in CSV context.

The current source audit retains all 20 explicit design tables / 328 rows. All
245 rows from 16 rectangular tables are independently reparsed and checked; four
irregular tables / 83 rows remain preserved but excluded from candidate links.
No short row is padded, long row truncated or blank record discarded. The eligible
tables provide 357 quantity-column candidates across six projects. Counts repeat
source records across quantity families and are not independent sample counts.
Missing candidates apply only to the acquired eligible source tables. Input-table
linkage, irregular layouts and other unreviewed column roles remain explicit gaps.

The **Reference identifiers** tab helps inspect identifier declarations before
combining gene and protein measurements. Select an acquired FASTA source, choose
**Original primary identifier** or **Explicit source gene declaration**, enter an
exact case-sensitive literal and select **Show original identifier records**.
Repeated primary identifiers and gene declarations retain distinct source rows.
CSV exports include exact original headers, declaration fields, source project
associations, input SHA-256 and header-index SHA-256.

The current header collection covers 46 acquired source occurrences from 54
planned FASTAs. Identical bytes shared across deposits produce 30 content indexes
with 297,621 header records; 75,597 records have explicit `gene` / `gene_id`
declarations. These are header records, not unique current genes or independent
reference versions. Eight files still await acquisition; GFF/GTF annotation files
are outside this FASTA header workflow. Native processing independently compares
text/regex parsing with binary/lexical parsing for every row and field, preserves
original header whitespace and records source/code/software hashes and operations.
No sequence transformation is performed.

Supported literal declarations are `gene`, `gene_id`, `transcript` and
`transcript_id` with `=` or `:` separators. Gene lookup uses only `gene` and
`gene_id`; unrecognised fields and gene symbols such as `GN` are not converted to
canonical identifiers. A source FASTA can be a contaminant list or small sequence
set; its presence does not prove it was the search database for a measurement
file. Source-to-measurement database linkage and current gene-version equivalence
remain unreviewed. These header data and the browser are included in
the 0.0.5.0 external data snapshot.

For PRIDE records, **Dataset coverage** shows primary project metadata for all
209 catalogue projects. The submitted SDRF collection contains 52 design files and
1,392 source rows across 49 projects. Open a project, expand **Submitted SDRF**, and
inspect each numbered source row. Column numbers preserve duplicate headings;
repeated sample labels, empty cells and literal missingness words remain unchanged.
Hover over a field to read its complete label or value. File associations report
exact primary file identifiers, including ambiguous or unresolved matches, without
claiming downloaded raw payloads or measurement-column pairing. The current
collection has 1,387 exact filename occurrences and five unresolved occurrences;
these counts are not independent samples. Source metadata origin and biological
truth status remain explicit. Projects without an SDRF in this collection retain
their project metadata and a visible design gap.

A native five-minute timer independently reparses all original SDRF cells with a
second parser, reconciles exact filename candidates against primary project file
records, and refreshes an immutable portable collection. Original source hashes,
validated project context, software versions, adapter source and operation logs
are retained. Controls, biological/technical replicate roles and assay compatibility
must follow submitted evidence; they are not inferred from shared filenames or
repeated channels. This metadata browser and its source records are included in 0.0.5.0.

For Metabolomics Workbench records, **Dataset coverage** shows original submitted
sample factors and analysis metadata, retaining source units and discrepancies
between deposited representations. Repeated sample labels remain separate rows;
raw filenames are declarations rather than evidence that payloads were downloaded.
These metadata records do not establish cross-assay compatibility or independent
biological replicates. Availability of original measurements is reported separately.

The **Dataset coverage** tab searches the complete local catalogue by accession,
title, source or assay. Its acquisition counts and local measurement status remain
separate: a downloaded source file or schema preview is not an admitted matrix.
Select a GEO study to inspect original sample labels, channels, submitted factors,
source relations and links to the same submitted sample in other studies. Controls,
paired samples and biological independence are not inferred from these links.
Use **Refresh coverage** to read updated acquisition progress. Set
`PARASCAPE_AUDIT` to a separately installed audit folder, containing the coverage
ledger and `all_geo_samples` metadata, when launching another installation.
The Live Data launcher selects the latest local audit automatically. After the
version 0.0.3.0 data download, other installations use the installed audit pack.
This view represents all 692 catalogue records; primary sample metadata covers
371 GEO records. Other sources and unavailable measurements retain explicit gaps.
Source-file counts describe the processing archive as of the displayed snapshot;
they do not mean the original source files are installed on this computer.

The **Published expression** tab reads separately installed, validated source packs.
Select the study and annotation variant, enter exact source gene identifiers, then
choose **Show deposited measurements**. CSV export includes the original values,
sample titles and source metadata, software version, source URL and source checksum.
Identifiers are case-sensitive and are not remapped to the current reference.

The local Live Data launcher selects the prepared published-expression folder.
After the version 0.0.3.0 data download, other installations use the included source
pack automatically. To select a different pack collection, set `PARASCAPE_PUBLISHED`
to the directory containing validated pack folders before starting the app.

The initial deposited RSEM pack represents GSE134099: six submitted samples against
ToxoDB release 38 and a StringTie annotation derived from that reference. These are
two representations of the same six samples. The two annotation variants remain
separate from each other and from the common-pipeline Salmon/kallisto matrices.
Original expected counts, TPM, FPKM, effective length and length retain their source
units; no normalization, filtering, rounding, imputation or gene remapping is added.
Source titles and attributes are preserved; biological independence and control
pairing are not inferred. Reported units alone do not establish comparability.

The source is complementary evidence retained with its original processing. No
citation-based replacement ranking is claimed; bibliometric comparison remains an
explicit unresolved gap. These measurements are not added to prediction features
or calibration. Pack manifests link archived inputs, software, transformations,
sample identities and output checksums; readers reject changed outputs or samples.
The original pack's independent validation checks every numeric source field.

For command-line CSV export:

```bash
python -m parascape.published_expression GENE_ID --pack /path/to/validated-pack \
  --reference ToxoDB_release38 --output deposited_expression.csv
```

The **Published expression → Export sample technical QC…** control exports metrics
for the complete selected annotation, independently of the gene query. Each sample
has separate TPM, FPKM and expected-count records with feature count, measured zeros,
positive values, missingness and numerical validity, plus source, matrix and manifest
hashes. Annotation variants remain separate representations of the same samples.
The generic numerical report preserves missing values and does not normalize or
impute measurements. Validity among observed values is a technical integrity metric;
it is not calibrated biological quality. Biological quality scores remain unavailable.

The **Gene expression → Export sample technical QC…** control provides the same
numerical metrics for all genes in the selected studies, quantifier and reference.
It retains submitted run entries separately and includes matrix, sample-metadata
and snapshot hashes. Declared snapshot output hashes must match before export;
external matrices without declared hashes are marked explicitly. Source metadata
and the snapshot's original normalization remain visible. The export adds no
normalization or imputation and does not treat repeated runs as biological replicates.

The **Deposited matrices** tab reads separately prepared GEO Series Matrix packs.
Choose one original study/platform-specific matrix and enter exact original probe
or submitted gene identifiers, one per line. Values retain their original sign,
missingness and source sample/channel metadata. The CSV export includes source and
matrix hashes, assay, quantity and comparison limits. Gene remapping, probe collapse
and cross-study normalization are not applied. A quantitative matrix can describe
an assay other than expression; its catalogue assay label remains visible.
Platform annotations and original scale/units need review before gene-level or
cross-assay comparisons. Local collections can be selected with
`PARASCAPE_SOURCE_MATRICES`; packaged collections live in external data under
`source_matrices`, keeping scientific data outside the application wheel.

Enable **Match source platform identifiers** to search explicit identifier fields
from the original platform annotation. Every matching probe remains separate;
measurements are not averaged by gene or alias. Hover over a value or export CSV
to inspect source fields, original annotation cells and multivalue associations.
Identifier matching is case-sensitive and requires a reconciled annotation index
prepared for the selected matrix. It reproduces a source association and does not
establish current gene equivalence. Use `PARASCAPE_PLATFORM_IDENTIFIERS` for a local
collection; portable collections live under external data `platform_identifiers`.

The **Deposited matrices → Export sample technical QC…** control checks every
original probe for each submitted sample, independently of the identifier query
or annotation checkbox. The report separates missing values, measured zeros,
negative finite values and nonfinite observed values. Source VALUE scales remain
unspecified, so negative values do not automatically fail numerical QC. Coverage
uses the original probe count; multiple probes are not expanded or collapsed
through gene identifiers. Original sample/channel metadata, assay and scale
limits, plus source, matrix, sample and manifest hashes accompany every row.
These reports describe technical integrity and coverage. They do not establish
biological quality, independent replicates or cross-study scale compatibility.
