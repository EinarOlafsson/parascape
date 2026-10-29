# Compare RNA-seq expression

The GitHub data download includes an offline snapshot of **296 verified study/run entries
from 23 studies**, selected from the ongoing 1,011-entry processing queue. This is a
partial release; a study may not yet contain all replicates or its control group.

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
