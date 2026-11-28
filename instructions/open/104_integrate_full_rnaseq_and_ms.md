# 104 · Integrate all RNA-seq and the uniform MS re-search into the app

**Einar, 2026-11-26** (asked "Are all the rnaseq and ms datasets downloaded/ redownloaded and integrated?", then chose
"Both, RNA first").

## State found (2026-11-26 ~08:30)

- RNA-seq: 2,535 runs planned (priority plan), 2,498 quantified, 2,129 verified by the queue. The app's
  snapshot (`app_snapshots/current`) was 20261121T120206Z-2002: 2,002 runs, 135 studies. The snapshot service
  was killed by the 2026-11-21 09:49 reboot mid-build and its timer never fired again.
- MS: 207 PRIDE projects selected; 156 searched (Sage, uniform), 13 unsupported formats, 38 left (3 held:
  PXD000297, PXD078724, PXD052656). Nothing from the uniform search reaches the gene table.
- PRIDE deposited tables: 3,474 of 3,519 files acquired; measurements not admitted.

## Done

- [x] Snapshot service restarted (timer + one build of 2,129 verified runs); later batches follow automatically.
- [x] `parascape/ms_genes.py` + `scripts/ms/build_ms_gene_tables.py` + `tests/test_ms_genes.py`: per project a
      Toxoplasma gene x raw-file table (q <= 0.01, peptides unique to one parasite gene, Top3, log2,
      median-centred) and `prepared_ms_genes/summary.csv`.

## Result of the first MS build (2026-11-26, sched 5053, 3 min 46 s)

- 163 project/instrument folders (155 projects); Toxoplasma genes quantified: >= 100 in 101, >= 500 in 74,
  >= 1,000 in 35, >= 2,000 in 10 (top: PXD044081 3,042 genes, PXD014133 2,736 over 72 files). 1,458 raw files in
  the 101 usable folders; 11 folders are isobaric (iTRAQ/TMT; values per raw file). Median parasite share of
  confident peptides 0.64.
- Positive control: per-gene mean abundance against the proteome already in the app (`proteome_log2_normalized*`):
  median Spearman 0.38 (IQR 0.31-0.53, max 0.81) over the 74 projects with >= 500 genes; against tachyzoite RNA
  median 0.37 (the usual protein-mRNA level). Near 0 only for PXD004426, PXD031649, PXD044091, PXD049694, which
  are candidates for enrichment / pull-down designs rather than whole proteomes - check before admission.

## To do

- [x] Snapshot 20261126T133001Z-2129 validated (snapshot_ready, 2026-11-26 ~09:54) and is `app_snapshots/current`:
      2,129 runs, 146 studies, 632 matrices (was 2,002 runs, 135 studies). The timer builds the next batches as
      the queue verifies the remaining quantified runs (2,498 quantified, 2,535 planned).
- [ ] Release the new snapshot to users (data release 0.0.14.0, with the STRING/CORUM and MS columns).
- [x] MS admission, per project (Einar 2026-11-26: "per-project now, designs later"): 93 of 158 folders admitted
      (>= 100 genes and Spearman >= 0.15 with the existing proteome; 57 too few genes, 8 enrichment-like), 5,422
      genes. `parascape/data/ms_uniform/` (table + per-folder reasons), registry entry `ms_uniform_research`
      (JOINED_AT_NEXT_RELEASE 0.0.14.0), joined in `build_graph.load_nodes` as `msu_<project>` columns.
- [x] CORRECTION (2026-11-26, after curation): the 91 admitted projects were curated by 4 Sonnet agents from
      PRIDE descriptions, protocols and file names (`parascape/data/ms_uniform/ms_curated_projects.tsv`,
      `ms_curated_files.tsv`; rules in mnt2/parascape/ms_design_curation/RULES.md; 40 low / 38 medium / 13 high
      confidence). Only 6 measure abundance (2 whole proteome, 3 secretome, 1 subcellular fraction); 54 IP, 16
      proximity labelling, 10 PTM enrichment, 2 crosslinking, 3 other. The abundance filter (Spearman >= 0.15) let
      enrichment projects through, so admission now also requires the curated type: 6 folders, 2,626 genes.
- [x] (2026-11-28) Enrichment projects as interaction evidence (`scripts/ms/build_ms_enrichment_pairs.py`; output
      mnt2/parascape/ms_enrichment/): 67 projects with >= 2 curated conditions, 352 conditions, per-condition
      log2 enrichment and z against the project's other files; top gene = putative bait, z >= 2 = prey: 7,438
      pairs, 283 baits. Pairs are 9.3x enriched for STRING-measured pairs and 3.9x for StarPath crosslinks over the
      same baits with random prey (crosslinks are the independent check; STRING's Toxoplasma evidence may include
      some of these IPs). Next: join as an edge layer / partner feature.
- [ ] (old item, superseded above) The 80 enrichment projects as INTERACTION evidence: bait (from the curated design) x prey enrichment over the
      project's controls, into the binding-partner benchmark (E12) and the network layers.
- [ ] Per-condition values for the 31 usable projects with an SDRF design
      (processing_audit/pride_sample_designs), then curated designs for the other 68 (Sonnet agents).
- [ ] Retry PXD000297 at 48 GB; re-download the corrupt PXD052656 file.
- [x] (2026-11-28) MS queue stalled on the largest projects: Sage OOM-killed at the unit's 16 GB (twice) and one
      truncated zip download. Unit MemoryMax raised to 32 GB (drop-in 70-memory-32g.conf; lease line in
      mnt3/claude/ram_leases.md; revert when the queue completes).
