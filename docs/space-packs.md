# Space packs

Large organism spaces can live outside the application wheel. The pack framework is implemented, and
so are the human and mouse gene packs: each carries its gene table, its protein-gene link table and
two lens layouts (`tissue`, `parasite_facing`), and opens as a map window. Published pack URLs are
pending: all ParaScape data hosting is moving to Hugging Face (Einar, 2026-11-13). The human/mouse
protein reference tables shipped in the wheel are separate spaces from these gene packs -- a protein
is not a gene -- and `parascape.host_links` carries a selection between the two.

A pack is a ZIP containing flat data files and `manifest.json`. The manifest records the organism,
version, node/graph filenames, each file's SHA256 and byte size, and a license code for every node
column. Three codes may be distributed:

| code | what it means |
|---|---|
| `O` | a verified open licence |
| `SA` | share-alike: it may travel, and the licence travels with it |
| `G` | a public deposit, redistributed with its accession and the paper behind it named |

`V` (verify the terms first) and `X` (local build only) are refused. A `G` column is accepted **only
when both an accession and a citation are recorded for it** (`column_sources` in the manifest, and
`column_sources=` on `packs.build` / `spaces.build_pack`): a deposit states no licence of its own,
so the permission to redistribute it is a permission to redistribute it WITH its attribution. The
check is in the manifest reader, not only in the builder, because the manifest is what an installed
pack is read back from. Einar gave that permission on 2026-11-13.

The builder's license declarations must come from the source records; hashes do not establish
redistribution rights or publisher identity.

Validate and install a local pack:

```bash
python -m parascape.packs verify /path/to/pack.zip --sha256 TRUSTED_SHA256
python -m parascape.packs install /path/to/pack.zip --sha256 TRUSTED_SHA256 --organism CODE
python -m parascape.packs installed CODE
```

For a published HTTPS URL, use `download URL --sha256 TRUSTED_SHA256 --organism CODE`.
Download requires a checksum obtained independently from the archive. Local installation can omit
the archive checksum when the local input is trusted; member checksums are always verified.

Installed versions live at `user_cache_dir()/spaces/CODE/VERSION/`. `PARASCAPE_STATE` overrides the
user cache root. Installation checks the entire archive before activation, refuses unsafe members,
and replaces `active.json` atomically. Previous versions remain available. Reinstalling identical
data is harmless; different data under the same version is refused. Reinstall an older archive to
select it again. No imports or path lookups download anything.

Prepare a pack for a registered organism:

```bash
python scripts/build_space.py --organism CODE --data /path/to/prepared/tables \
  --version VERSION --licenses /path/to/column-licenses.json --output /path/to/pack.zip
```

The license file is a JSON mapping of every node column to its license code. The shared builder
requires unique canonical gene IDs, matching graph order, finite coordinates and valid edge
indices. It compares with the installed table by default, or with `--previous TABLE.parquet`, and
refuses lost identifiers, columns or previously measured cells. This packages prepared data; it
does not acquire or derive new measurements.

Python callers use `parascape.spaces.build_pack`, `parascape.packs.install` and
`parascape.organisms.nodes_path` / `graph_path`. A registry entry with `distribution="pack"` resolves
only its active installed pack. Wheel spaces retain the existing `PARASCAPE_CACHE` resolution.
