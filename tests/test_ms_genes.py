"""Gene-level protein abundance from Sage label-free tables (Einar 2026-11-26)."""
import numpy as np
import pandas as pd

from parascape import ms_genes as M


def _lfq():
    rows = []
    for i in range(5):            # gene A: 5 unique peptides
        rows.append(("PEPA%d" % i, 2, "TGME49_000001-t26_1-p1", 0.001, 2 ** (20 - i), 2 ** (21 - i)))
    rows.append(("SHARED", 2, "TGME49_000002-t26_1-p1;HOST_ENSP00000000001.1", 0.001, 2 ** 25, 2 ** 25))
    rows.append(("TWOGENES", 2, "TGME49_000001-t26_1-p1;TGME49_000003-t26_1-p1", 0.001, 2 ** 25, 2 ** 25))
    rows.append(("BADQ", 2, "TGME49_000004-t26_1-p1", 0.2, 2 ** 25, 2 ** 25))
    rows.append(("PEPB", 2, "TGME49_000005-t26_1-p1", 0.001, 2 ** 10, 0.0))
    t = pd.DataFrame(rows, columns=["peptide", "charge", "proteins", "q_value", "f1.mzML.gz", "f2.mzML.gz"])
    t["score"] = 1.0
    t["spectral_angle"] = 0.9
    return t


def test_only_confident_peptides_unique_to_one_parasite_gene_count():
    G, peptides = M.gene_table(_lfq())
    assert set(G.index) == {"TGME49_000001", "TGME49_000005"}
    assert peptides["TGME49_000001"] == 5


def test_top3_and_an_unmeasured_file_is_missing_not_zero():
    lfq = _lfq()
    G, _ = M.gene_table(lfq)
    raw = np.log2(np.mean([2 ** 20, 2 ** 19, 2 ** 18]))       # top 3 of gene A in file 1
    # median-centred: compare differences between genes within a file
    assert np.isnan(G.loc["TGME49_000005", "f2.mzML.gz"])
    assert np.isclose(G.loc["TGME49_000001", "f1.mzML.gz"] - G.loc["TGME49_000005", "f1.mzML.gz"], raw - 10)


def test_the_summary_counts_what_a_project_contributes():
    lfq = _lfq()
    G, _ = M.gene_table(lfq)
    s = M.summary("PXDTEST", lfq, G, labelled=False)
    assert s["genes_quantified"] == 2 and s["files"] == 2 and s["genes_in_every_file"] == 1
