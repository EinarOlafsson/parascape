# Per-sample RNA-seq quality (v0, 2026-11-06)

Built by scripts/build_rna_sample_quality_notebook.py (notebook: notebooks/rna_sample_quality_2026_11_06.ipynb) from snapshot 20261107T012052Z-735.
717 verified runs scored; quality quartiles 0.794 / 0.916 / 0.967; 176 flagged.

## Formula
Each metric is mapped linearly to [0,1] between fixed anchors (same for every study), clipped. Quality = geometric mean of the available
component scores (each floored at 0.01); `n_components` is how many entered (1 to 6). `quality_no_replicate` omits replicate agreement and is used for the leave-out check.

| component | score 0 | score 1 |
|---|---|---|
| parasite reads (sum of counts in salmon ME49_host matrix, log10) | 1e5 | 1e7 |
| parasite genes with TPM >= 1 (of 8157) | 3000 | 7000 |
| mapping rate (salmon percent_mapped / 100, combined parasite + host reference) | 0.30 | 0.80 |
| compatible fragment ratio (salmon lib_format_counts) | 0.80 | 0.98 |
| salmon vs kallisto Spearman of parasite gene TPM | 0.80 | 0.98 |
| replicate Spearman (log2 TPM+1 vs median of other runs of study x condition_group, >= 2 others, genes TPM >= 1 in either) | 0.80 | 0.98 |

Reported but not scored: total reads, host share, parasite share of mapped reads (1 - host_share).
Flags: parasite reads < 1e6; genes < 4000; mapping < 0.5; compatible ratio < 0.9; salmon-kallisto < 0.9; replicate < 0.9; parasite share < 2%; fewer than 3 components.

## Leave-out check
Quality without the replicate component vs replicate agreement (n=442 runs in 41 studies): Spearman 0.496 overall, 0.373 after removing study means.

## Caveats (v0)
* Anchors are a first guess, not fitted; revisit once more studies are verified and the score distribution is known.
* Replicate agreement depends on `condition_group` being a true replicate grouping; the design marks many groupings as needing review, and a wrong grouping looks like a bad replicate.
* Mapping rate is on the combined parasite + host reference, so it does not say how well the parasite alone was captured; parasite depth and share cover that.
* A run in a low-replicate study has fewer components and its score is less constrained; check n_components.
* Biological differences between conditions inside a group (time courses grouped together) lower replicate agreement without being technical problems.
