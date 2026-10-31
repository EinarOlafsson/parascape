# Coexpression

Open **Compare → Coexpression**, select **RNA** or **Protein**, select a source,
and click **Load selected source**. Missing optional protein or deposited RNA
files use the existing verified component download dialog. Review the submitted
column metadata, select columns and enter an exact original feature identifier.
Choose an original row identity when a label occurs more than once.

The tab computes seed-to-feature correlations across selected columns of one
source and quantity. RNA quantifier/reference combinations remain distinct;
deposited gene, transcript, target and protein-group rows retain their source
identity. Ratios and isoform percentages are excluded from the abundance selector.
Other original quantity roles, including spectral counts, remain explicitly named.
Some sources have too few columns for this analysis.

Spearman is the default. Pearson and explicit log2(x + 1) are available; the latter
requires nonnegative finite values. Already logged quantities are not automatically
logged again. Missing and nonfinite observations are excluded separately for each
pair; measured zeros remain measured zeros. Ranks are recomputed on each pair's
overlap. Constant or insufficient pairs retain a reason instead of a coefficient.
Exports include overlap and nonzero counts, every candidate row, source hashes,
selected metadata, numerical input hash, transformation and implementation hash.

Results are descriptive. Column identities do not establish independent biological
samples. Known repeated sequencing of the same submitted sample is refused, while
unknown relationships require study-design review. Small cohorts have a visible
warning. This tab provides no significance tests, biological quality score,
physical-interaction inference or cross-study pooling. A strong correlation can
reflect stage mixtures, batch, host fraction or shared non-detection patterns.

Click a result to inspect its source-scale scatter plot. The heatmap displays the
seed and up to twelve top profiles with a **display-only** within-row z score.
Missing cells are masked. This visualization never replaces original measurements
or changes exported coefficients. Navigation controls support zooming and saving
figures. Editing the source, selected columns, seed or settings invalidates results
and disables stale exports. File verification and correlation run outside the UI
thread; cancellation retains checked source files.

This analysis differs from RNA–protein agreement **across genes**, which needs
separately reviewed gene mappings and matched biological contexts. It also differs
from a consensus network across studies or a full module-discovery method. Those
require additional design, replication and evaluation.

`tests/test_coexpression.py` checks numerical results against independent SciPy
calculations, missingness, extreme scales, repeated labels/samples, source integrity,
and RNA/protein Qt load/plot/export/invalidation. `scripts/smoke_coexpression.py`
replays actual checked source tables and records a software-only verification.
