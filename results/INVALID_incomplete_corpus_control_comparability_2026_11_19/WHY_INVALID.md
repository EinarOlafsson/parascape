# Not a result: run on an incomplete corpus

Re-run 2026-11-19 only to save the 200 null draws for figure 2C. It used **8 studies**
against the **12** of the valid run of 2026-11-18, with the same 8,157 genes, and its gain
fell from +0.055 (p = 0.005) to +0.031 (p = 0.090).

The four studies missing were absent for pipeline reasons, not biological ones:

| study | why it dropped |
|---|---|
| GSE223819 | 0 control runs, fewer than 2 |
| GSE223869 | no built matrix |
| GSE223877 | no built matrix |
| PRJNA1008424 | no built matrix |

Their expression matrices were being rebuilt by the quantification queue (the
`--single-overhang` re-quantification) at the time of the run. So this run measured an
incomplete corpus and its numbers must not be used.

Renamed OUT of the `control_comparability_*` pattern (to `INVALID_incomplete_corpus_...`) so the figure builder, which takes the NEWEST
`control_comparability_*` directory, cannot pick it up. A first rename to `control_comparability_INCOMPLETE_...` still matched the glob and sorted AFTER the valid run, so it would have redrawn
figure 2 silently on 8 studies. The valid result is `control_comparability_2026_11_18`.

Figure 2C needs the null draws; re-run once the queue has rebuilt all twelve matrices, and
check the study count is 12 before trusting it.
