# Sequence and 3Di traits (instruction 72, 2026-11-04)

One row per gene of parascape/data/nodes.parquet (8,140 genes), file sequence_traits.parquet, index gene_id. Built by scripts/build_sequence_traits_notebook.py (notebook notebooks/sequence_traits_2026_11_04.ipynb).

## Rules
* seq_mature: DeepTMHMM signal peptide removed (sequence[sp_end:]); otherwise equal to seq_full.
* seq_nterm_flank: start of the mature sequence to the residue before the FIRST TM helix; seq_cterm_flank: residue after the LAST TM helix to the C-terminus.
* NO-TM RULE: a protein with no TM domain has BOTH flanks equal to the whole mature sequence, so the two flank columns are identical (intended; asserted).
* tdi_*: the same four regions with the same coordinates, cut from the Foldseek 3Di string (full protein length); empty only for genes with no 3Di.
* nterm_side / cterm_side: inside or outside when >= 90% of the flank residues carry that I/O call in the DeepTMHMM per-residue topology string; unknown if mixed, empty or absent.
* tdi_source: naive_db preferred, then cleaned_db (suffixes such as __sol stripped; the record matching the protein length is chosen), then af3, then af3_fragments, then aligned / aligned_multi (section below), else none.
* aligned / aligned_multi: for genes whose structure(s) have a DIFFERENT sequence from the ToxoDB gene model (ToxoDB model merged or extended relative to the UniProt entry AlphaFold folded). The structure own sequence is aligned to ours (semi-global, free end gaps, identity scoring match +5 / mismatch -4 / gap 12+2(k-1); Biopython not installed, own numba implementation) and the 3Di letter of each aligned residue is copied to the residue it aligns to; positions without an aligned residue are X. Accepted only when a structure has >= 95% identity (matches over alignment columns, internal gaps counted) and the merged result covers >= 50% of the protein; several structures are merged, higher identity first, no overwriting. Candidates: naive_db and cleaned_db records of the gene, AF3 models of its accessions, AF3 models of accessions ambiguous between genes or unknown to DeepTMHMM (screened by shared 6-mers, then aligned). tdi_align_identity / tdi_align_coverage hold the identity (over the structures that contributed) and coverage; NaN for other sources; tdi_n_fragments = number of structures merged. Per-candidate numbers: tdi_alignment_candidates.tsv; accepted genes: tdi_aligned_genes.tsv. AN ALIGNED 3Di IS WEAKER EVIDENCE THAN A DIRECT ONE: any model using tdi_* must carry tdi_source.
* af3: Foldseek 3Di computed here (version 7f4b94644d6b210af4530e04a9fdc9e4bb64ea66) from the whole-protein AlphaFold3 model af3_<UniProt>_model.cif in /media/carruthers/mnt1/af3/af_output, UniProt mapped to the gene through the accession column of deeptmhmm_toxo.csv; used only where no database 3Di of the right length exists.
* af3_fragments: folders af3_<UniProt>_f<i>of<n>_<start>-<end> stitched at their 1-based inclusive ranges into a string of the full protein length; residues in no fragment are 'X' (not a 3Di letter); overlaps are resolved by the higher mean pLDDT (atom_plddts of the confidences json), else the earlier fragment. tdi_coverage = fraction of residues with a real letter; tdi_n_fragments = models used.
* ACCEPTANCE (direct sources): a 3Di string is kept only when its length equals the nodes.parquet protein length (stitched: every fragment length equals its range and the ranges fit inside the protein). Otherwise empty and listed in tdi_length_mismatch.tsv; nothing padded or trimmed.
* Stitched strings contain 'X'; downstream code must treat X as missing.

## Coordinate convention check
1-based inclusive confirmed on 6909 of 6909 domains (signal peptides and TM helices, all genes). sequence[start-1:end] equals the recorded tm_<i>_sequence / sp_sequence.

## Gaps
* No propeptide / TEXEL cleavage: no reliable per-gene table; only the signal peptide is removed.
* 3Di is missing for 637 genes: 264 with no usable model anywhere (no database record and no af3_<UniProt> folder that maps unambiguously to the gene) and 373 whose model(s) differ in length from the annotated protein and could not be aligned at >= 95% identity with >= 50% coverage (tdi_length_mismatch.tsv, with the best alignment identity; excluded, not trimmed). 520 genes of the former mismatch/no-model group were filled by alignment (tdi_source aligned / aligned_multi); at 90%/30% 186 more would pass.
* The AF3 tree has 1,536 distinct af3_<UniProt> ids, all UniProt (none TGME49-named); the TGME49-named folders (bait_, pep_, xl_, ...) are domain or complex jobs and are not used. 123 accessions are not in the DeepTMHMM accession column and 16 map to more than one gene (attached to none; af3_ambiguous_accessions.tsv).
* Whole-length direct sources (naive_db, cleaned_db, af3) are accepted on length only: the sequence was not compared with the annotated one, so an isoform of identical length would pass. Aligned sources are checked by identity.
* DeepTMHMM file has 8140 rows (the 8,141 quoted in the spec counts the header line), all identifiers unique, all present in nodes.parquet.
* Long strings are registered with kind 'category' (no better-fitting kind).

## Coverage

```
                                                       item  genes  fraction
                                  column seq_full non-empty   8140     1.000
                                column seq_mature non-empty   8140     1.000
                           column seq_nterm_flank non-empty   8140     1.000
                           column seq_cterm_flank non-empty   8140     1.000
                                  column tdi_full non-empty   7503     0.922
                                column tdi_mature non-empty   7503     0.922
                           column tdi_nterm_flank non-empty   7503     0.922
                           column tdi_cterm_flank non-empty   7503     0.922
                                      column n_tm non-empty   8140     1.000
                        column has_signal_peptide non-empty   8140     1.000
                                column nterm_side non-empty   8137     1.000
                                column cterm_side non-empty   8137     1.000
                                   column seq_len non-empty   8140     1.000
                                column mature_len non-empty   8140     1.000
                                 column nterm_len non-empty   8140     1.000
                                 column cterm_len non-empty   8140     1.000
                                column tdi_source non-empty   7503     0.922
                              column tdi_coverage non-empty   8140     1.000
                           column tdi_n_fragments non-empty   8140     1.000
                        column tdi_align_identity non-empty    520     0.064
                        column tdi_align_coverage non-empty    520     0.064
                          column flanks_identical non-empty   8140     1.000
                                   no TM domain (n_tm == 0)   6729     0.827
                                             with TM domain   1411     0.173
                                             signal peptide    945     0.116
                                      signal peptide and TM    269     0.033
                                          3Di from naive_db   5478     0.673
                                        3Di from cleaned_db    332     0.041
                                 3Di from af3 (whole model)    872     0.107
                          3Di from af3_fragments (stitched)    301     0.037
                           3Di from aligned (one structure)    462     0.057
                3Di from aligned_multi (several structures)     58     0.007
                           aligned with gaps (coverage < 1)    407     0.050
  would also pass at 90% identity / 30% coverage (not used)    186     0.023
                          stitched with gaps (coverage < 1)     60     0.007
no 3Di: model(s) of different length, alignment fails 95/50    373     0.046
                           no 3Di: no usable model anywhere    264     0.032
                                           flanks identical   6729     0.827
                             TM genes with nterm_side known   1408     0.173
                             TM genes with cterm_side known   1408     0.173
```
