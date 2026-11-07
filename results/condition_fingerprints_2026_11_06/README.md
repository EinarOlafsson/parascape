# Condition fingerprints v1 (instruction 75)

Notebook: `notebooks/condition_fingerprints_2026_11_06.ipynb` (executed). Builder: `scripts/build_condition_fingerprints_notebook.py`. Snapshot 20261107T012052Z-735, 668 bulk runs scored.

## Validation on held-out studies

| axis | built_on | tested_on | value | chance | p | detail | ships |
|---|---|---|---|---|---|---|---|
| cell cycle | Xue 2020 single-cell phases | synchronised series cellcycle19092 (hours 1-12) | 1.0 | 0.502 | 0.001 | 1.44 turns of the cycle over 12 h; dominant phase by hour: S M M C G1a G1a G1b G1b S M G1a G1a | True |
| bradyzoite response | GSE223869 alkaline vs standard | 7 other studies that state the condition | 0.934 | 0.5 | nan | 6 of 7 studies above 0.7; converter removed vs present AUC 1.00, 1.00 | True |
| bradyzoite conversion | GSE223869 alkaline vs standard and without BFD2 | 7 other studies that state the condition | 0.922 | 0.5 | nan | 6 of 7 studies above 0.7; converter removed vs present AUC 1.00, 1.00 | True |

## Removing the converter under stress

| study | removed | fingerprint | present_runs | removed_runs | fraction_present | fraction_removed | auc |
|---|---|---|---|---|---|---|---|
| GSE132238 | BFD1 knocked out | response | 5 | 5 | 0.97 | 0.58 | 1.0 |
| GSE132238 | BFD1 knocked out | conversion | 5 | 5 | 1.0 | 0.51 | 1.0 |
| GSE223819 | BFD1 depleted by auxin | response | 2 | 2 | 0.98 | 0.65 | 1.0 |
| GSE223819 | BFD1 depleted by auxin | conversion | 2 | 2 | 1.0 | 0.54 | 1.0 |

Both bradyzoite fingerprints fall when BFD1 is removed; the conversion-specific one falls only slightly further, so most of the alkaline response already depends on BFD1.

## GSE223819: genotype field and sample title disagree

| genotype_field | sample_title | auxin | vehicle |
|---|---|---|---|
| 200385 | 200385 | 0.54 | 1.0 |
| 208020 | 311100 | 0.32 | 0.29 |
| 224630 | 224630 | 1.0 | 1.0 |
| 253790 | 306620 | 1.0 | 1.0 |
| 306620 | 253790 | 1.0 | 1.0 |
| 311100 | 208020 | 1.0 | 1.0 |
| parental | parental | 1.0 | 1.0 |

The line titled 311100 (BFD2) converts poorly even on vehicle; the fingerprint sides with the titles. Evidence only; the registry keeps the deposit text.

Files: `sample_annotations.tsv` (every run: phase mixture and angle, both bradyzoite fractions, the stated condition, `filled:`/`exploratory:` source and `depends_on`), `fingerprints.tsv`, `validation.tsv`, `bradyzoite_per_study.tsv`, `converter_removal.tsv`, `unlabelled_studies.tsv`. Annotations are inferred and are never evidence for a gene trait.
