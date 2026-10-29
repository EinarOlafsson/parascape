# Compare RNA-seq expression

The GitHub data download includes an offline snapshot of **305 verified study/run entries
from 24 studies**, selected from the ongoing 1,011-entry processing queue. This is a
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
