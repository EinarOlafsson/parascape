# 102. One DuckDB + Parquet store for the release data

**Einar, 2026-11-24:** "what type of database is parascape using? I have recently overhauled spacrs database structure to use duckDB and parquett. what database structure would be best for parascape? ... is the data stored in a way that makes all parascape functions work as fast as possible?" He chose **"Start-up fixes, then DuckDB"**.

## Done first: the start-up fixes

Commit `7a9b1bcd` cut opening from 17.8 s to 6.7 s by building docks and panels when they are first shown. Reading the data took only 1.8 s of the 17.8.

## What a release holds (data-0.0.13.0, measured 2026-11-24)

| Files | Count | Size | Read by |
|---|---|---|---|
| JSON | 8,378, of which 5,513 in `processing_audit/` and 1,736 in `deposited_values/` | | the audit and coverage panels, on demand |
| Parquet | 998, of which 634 in `processing_audit/`, 147 in `source_matrices/` and 124 in `rnaseq/` | | the gene table, the track record, the omics |
| CSV / TSV | 1,113, of which 915 in `processing_audit/` and 125 in `rnaseq/` | | the same |
| all | 10,934 files | 1.2 GB | |

## Read costs (cold disk cache not forced; WS, mnt1)

| What | Time |
|---|---|
| `nodes.parquet` (8,140 x 443) | 0.38 s |
| `track_record.shipped()` (1,047,574 rows), first call | 0.37 s |
| all 5,513 audit JSON files | 2.2 s |
| all 124 RNA-seq Parquet files | 2.3 s |
| 906 audit CSV files | 0.5 s |

**Reading is not what makes ParaScape slow.** The case for DuckDB is structure, not speed:
- one queryable store instead of about 10,900 files;
- a schema per table;
- queries across studies without loading every table;
- the same layout as spaCR (DuckDB + Parquet).

## Plan, if Einar confirms it

1. **Store:** keep Parquet as the storage format. Write one `parascape.duckdb` per release, holding views over the Parquet files plus small tables converted from the JSON records. The release stays readable without DuckDB.
2. **Convert first:**
   - the 5,513 audit JSON records, as one Parquet table per record kind;
   - the 1,736 `deposited_values` JSON files;
   - the RNA-seq per-study tables, as one long table (study, run, gene, value) partitioned by study.
3. **Query layer:** `parascape.store`, the one place that opens the store. Loaders move to it one at a time, each with a test that the old and new readers give the same frame (within dtype tolerance, as spaCR does).
4. **Dependency:** `duckdb` (2.0.0, as in the spacr env) added to `pyproject.toml` and to the starplast env.

## Decision (Einar, 2026-11-25)

Asked whether the structural gain was worth a new dependency, given that reading is not slow, he chose **"Parquet first, no DuckDB yet"**:
- convert the 7,249 JSON records and the per-study RNA-seq tables into a few Parquet tables, with no new dependency;
- add DuckDB later, only if a query needs it.

Steps 1 and 4 of the plan above wait. Steps 2 and 3 go ahead with pandas/pyarrow as the reader.

## Done when

- the Parquet tables are built by the release build;
- every loader of the converted records reads through `parascape.store`;
- the file count of a release is under 500;
- the full suite passes.

## Progress 2026-11-25

- **`parascape/store.py`:** a packed folder is one `<folder>.pack.parquet` (path, exact bytes, SHA-256). The exact bytes are kept because the audit records are content-addressed: their names, receipts and manifests carry digests of the bytes. `store.read_bytes(path)` serves a file from disk, or else from its folder's pack.
- **Release builder:** it packs every audit folder that no module names (worked out from the code at build time, so a folder that gains a reader ships loose again), plus the folders in `store.READ_THROUGH_STORE`. Those are now `all_geo_samples` and `pride_projects`, whose only reader, the coverage panel, moved to the store.
- **Measured on data-0.0.13.0:** 10,934 files become **5,219**. The archive stays 945 MB. One pack restored byte for byte (1,894 of 1,894 files).
- **Next:**
  - `deposited_values` (1,736 files, 8 reader modules);
  - `protein_column_design_links` (611, read by `source_design_links`);
  - `source_matrices` (741);
  - `rnaseq` (251).
