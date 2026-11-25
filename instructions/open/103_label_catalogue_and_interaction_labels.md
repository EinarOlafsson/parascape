# 103 · Label catalogue with metadata; STRING and CORUM labels for Toxoplasma

**Einar, 2026-11-25** (his words): "CORUM and STING, etc, should both have toxo annotations i think, these
could be added to the Toxoplasma labels." Then, choosing to keep measured and copied evidence apart: "i agree
with the not mixing just make sure the label name indicates what type of evidense it is. also there should be
assoceated metadata with the label catalog detailing where each label column comes from and how it was treated
and generated and the labels it contains, etc."

## What the sources hold for Toxoplasma (checked 2026-11-25)

- **STRING v12.0** has T. gondii VEG (taxon 432359, a "periphery" genome). Measured in Toxoplasma itself:
  1,814 experimental pairs over 119 proteins and 909 curated-database pairs over 307 proteins. Everything else
  is transferred from other organisms (about 1.1 M experimental and 0.67 M database pairs, both directions
  counted). VEG locus numbers equal ME49 ones (`TGVEG_n` -> `TGME49_n`).
- **CORUM** is mammalian only; Toxoplasma labels come only by orthology (OrthoMCL 6.21 orthogroups).

## Built

- `scripts/labels/build_interaction_labels.py`: three labels, each name stating its evidence:
  - `interaction module (STRING, measured in Toxoplasma)`: 199 genes, 8 modules (observed);
  - `interaction module (STRING, copied by orthology)`: 2,469 genes, 14 modules (copied; scores >= 700 on
    transferred experiments or databases);
  - `complex (CORUM human, copied by orthology)`: 414 genes, 26 complexes (each gene in its most
    specific complex that keeps 10 genes).
  STRING co-expression and text mining are not used (co-expression is circular with the transcript view).
  Pair tables for the binding-partner benchmark: 2,555 measured, 28,651 copied STRING pairs; 18,477 CORUM
  pairs. Source files and their SHA-256 are listed in `sources.json`. Build output:
  `mnt2/parascape/labels/interaction_2026_11_25/`; raw downloads in `mnt2/parascape/labels/raw/`.
- `parascape/label_catalogue.py`: one entry per label (evidence type, sources resolved through the dataset
  registry, source column, how it was made, what is left out of the features); `write()` adds genes
  labelled, the classes with counts and the columns actually left out. E1 and E8 write
  `label_catalogue.tsv/.json` into every run folder. A label without an entry raises.
- Every benchmark script reads the new labels through `e1_information.added_labels` (the build output
  path is `PARASCAPE_INTERACTION_LABELS`) until a release joins them to the node table.
- Tests: `tests/test_label_catalogue.py`.

## Still to do

- [ ] Dataset registry entries for the three sources and the join at the next data release (0.0.14.0).
- [ ] Score the new labels in E1 (they are in E8 now) and use the pair tables in the binding-partner benchmark.
- [ ] Older labels keep their names (frozen results use them); their evidence is in the catalogue. Rename
      only if Einar asks.
