# Sequence and 3Di traits (instruction 72, 2026-11-04)

One row per gene of parascape/data/nodes.parquet (8,140 genes), file sequence_traits.parquet, index gene_id. Built by scripts/build_sequence_traits_notebook.py (notebook notebooks/sequence_traits_2026_11_04.ipynb).

## Rules
* seq_mature: DeepTMHMM signal peptide removed (sequence[sp_end:]); otherwise equal to seq_full.
* seq_nterm_flank: start of the mature sequence to the residue before the FIRST TM helix; seq_cterm_flank: residue after the LAST TM helix to the C-terminus.
* NO-TM RULE: a protein with no TM domain has BOTH flanks equal to the whole mature sequence, so the two flank columns are identical (intended; asserted).
* tdi_*: the same four regions with the same coordinates, cut from the Foldseek 3Di string, only where its length equals the protein length; empty otherwise.
* nterm_side / cterm_side: inside or outside when >= 90% of the flank residues carry that I/O call in the DeepTMHMM per-residue topology string; unknown if mixed, empty or absent.
* tdi_source: naive_db preferred, then cleaned_db (suffixes such as __sol stripped; the record matching the protein length is chosen, since cleaned_db also holds domain fragments).

## Coordinate convention check
1-based inclusive confirmed on 6909 of 6909 domains (signal peptides and TM helices, all genes). sequence[start-1:end] equals the recorded tm_<i>_sequence / sp_sequence.

## Gaps
* No propeptide / TEXEL cleavage: no reliable per-gene table; only the signal peptide is removed.
* 3Di is missing for 2,330 genes: 1,651 in neither Foldseek database and 679 whose AlphaFold-derived 3Di length differs from the annotated protein (excluded, not trimmed).
* DeepTMHMM file has 8140 rows (the 8,141 quoted in the spec counts the header line), all identifiers unique, all present in nodes.parquet.
* Long strings are registered with kind 'category' (no better-fitting kind).

## Coverage

```
                                         item  genes  fraction
                    column seq_full non-empty   8140     1.000
                  column seq_mature non-empty   8140     1.000
             column seq_nterm_flank non-empty   8140     1.000
             column seq_cterm_flank non-empty   8140     1.000
                    column tdi_full non-empty   5810     0.714
                  column tdi_mature non-empty   5810     0.714
             column tdi_nterm_flank non-empty   5810     0.714
             column tdi_cterm_flank non-empty   5810     0.714
                        column n_tm non-empty   8140     1.000
          column has_signal_peptide non-empty   8140     1.000
                  column nterm_side non-empty   8137     1.000
                  column cterm_side non-empty   8137     1.000
                     column seq_len non-empty   8140     1.000
                  column mature_len non-empty   8140     1.000
                   column nterm_len non-empty   8140     1.000
                   column cterm_len non-empty   8140     1.000
                  column tdi_source non-empty   5810     0.714
            column flanks_identical non-empty   8140     1.000
                     no TM domain (n_tm == 0)   6729     0.827
                               with TM domain   1411     0.173
                               signal peptide    945     0.116
                        signal peptide and TM    269     0.033
                            3Di from naive_db   5478     0.673
                          3Di from cleaned_db    332     0.041
no 3Di: record of different length (excluded)    679     0.083
               no 3Di: not in either database   1651     0.203
                             flanks identical   6729     0.827
               TM genes with nterm_side known   1408     0.173
               TM genes with cterm_side known   1408     0.173
```
