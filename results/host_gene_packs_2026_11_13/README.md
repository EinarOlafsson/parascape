# Host gene packs 2026.11.13

`parascape_Hs_2026.11.13.zip` (19.5 MB) and `parascape_Mm_2026.11.13.zip` (3.9 MB) are **not in
git**. Their SHA256 are in `SHA256` and what they contain is in `packs.json`.

They are left out because the move of ParaScape data hosting to Hugging Face (instruction 99,
Einar 2026-11-13) exists precisely to keep archives out of the repository, and that work deleted
1,205 component assets from GitHub on the same day. Adding 23 MB of new ones beside it would undo
what it was for.

Rebuild them from the deposits already on disk:

    python scripts/build_host_gene_lenses.py     # ~1 h on 32 shared cores
    python scripts/build_host_gene_packs.py      # seconds

The archives then appear here and `space_catalogue` installs them. **Publishing them belongs to the
Hugging Face hosting work**, not to instruction 53: both should be uploaded to the
`parascape-data` dataset and `space_catalogue` pointed at that URL, with the SHA256 above.

Until then both gene spaces are listed as installable only in a checkout that has built them, and
every test that needs one skips with that reason rather than failing.
