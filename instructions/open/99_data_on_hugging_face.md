================================================================================
Host all ParaScape data on Hugging Face, not on GitHub releases
================================================================================

THE DECISION

    Einar, 2026-11-13: "all of the data that parascape uses should be hosted on
    huggingface in my account. otherwise the space limits on github will fill
    quickly." The account is `einarolafsson`.

    Also decided the same day, and recorded here because it belongs to the same
    publishing pass: a column derived from a GEO series ships with that accession
    and the paper it came from, so a reader of a value can reach its source.

STATE

    Every artefact the application or its tooling fetches, measured on
    2026-11-13 from `parascape/data_manifest.json` and from the component index
    `data_components.json` inside the installed 0.0.12.0 data.

    | artefact | host | size | checksum | code path |
    |---|---|---|---|---|
    | core data archive `parascape-data-0.0.12.0.tar.gz` | **Hugging Face** `einarolafsson/parascape-data` since today; GitHub release `v0.0.12.0` kept | 944.3 MB, 10,929 files | archive SHA-256 plus a SHA-256 and size for every file, in `parascape/data_manifest.json` in the wheel | `parascape/data_manager.py` -- `download()`, `install_archive()`; `parascape-data --download`; the first-launch dialog `parascape/data_download_dialog.py` |
    | optional component packs, 1,195 of them in four kinds: `deposited_values` 574 (8.0 GB), `sparse_axes` 358 (876.6 MB), `sparse_coordinates` 180 (33.3 GB), `workbook_cells` 83 (1.1 GB) | **Hugging Face** `einarolafsson/parascape-data`, under `components/`, since today; the GitHub release assets on tags v0.0.8.0, v0.0.9.0 and v0.0.9.0-quantities-data are kept | 43.2 GB in 1,205 parts, each at most 1 GiB | per-part SHA-256 and size, whole-archive SHA-256, and a SHA-256 for every extracted file, all in the `data_components.json` that ships inside the verified core archive | `parascape/data_components.py` -- `read_index()`, `validate()`, `download()`; chosen through `parascape/component_download_dialog.py` and the resolvers `optional_sparse_data.py`, `optional_workbook_data.py`, `optional_deposited_data.py` |
    | organism and host gene/protein packs (instruction 53) | nothing published yet | a pack is at most 512 MB packed, 1 GB unpacked | manifest inside the zip, with a SHA-256 for every member; the publisher's SHA-256 has to be supplied separately to `download()` | `parascape/packs.py` -- `build()`, `verify()`, `download()`, `install()`; resolved by `paths.space_file()` |
    | raw source datasets (GEO, PRIDE, ToxoDB and the rest) | their own archives | tens of GB | per-source, in the dataset registry | `scripts/download_apicomplexan_archive.py`, `parascape/datasets.py`; only `build_graph.py` needs them, and they are never fetched at run time |
    | protein structures | AlphaFold DB and the PDB | on request | none; cached by identifier | fetched when a structure is opened, cached in the user cache |

    So the 43.2 GB of component packs, not the 944.3 MB core archive, is what
    fills the GitHub quota. The core archive moved first because it is the one
    every installation fetches; the packs followed the same day, and both are
    now served from Hugging Face with the GitHub copies kept behind them.

WHAT MOVED TODAY

    `https://huggingface.co/datasets/einarolafsson/parascape-data`

        README.md                                     what it is, how to verify it, provenance
        SHA256SUMS                                    the archive checksum, one line per release
        releases/0.0.12.0/parascape-data-0.0.12.0.tar.gz
        releases/0.0.12.0/data_manifest.json          the descriptor, 10,929 file entries
        releases/0.0.12.0/SHA256SUMS                  those 10,929 checksums, one per line

    One dataset repository with a version directory, not one repository per
    release: the client builds its address from the host, the version and the
    asset name, so a new release needs no new repository and no new setting, and
    a reader sees every release in one place with one set of terms.

    The uploaded archive is byte-identical to the GitHub release asset for
    v0.0.12.0 -- same 944,323,064 bytes, same SHA-256
    `6e676191021fec9a5afd8e7762ae7f5ba042de6791b4c4eb238bcd25be958fdb`, checked
    against the descriptor in the wheel before the upload and not re-downloaded
    from GitHub. Nothing unpublished was uploaded: only the asset that was
    already public on the GitHub release.

    The archive itself is uploaded rather than its 10,929 files individually.
    That keeps the checksum contract and the atomic install exactly as they
    were, and keeps a first install at one request instead of eleven thousand.
    The per-file checksums are published beside it, so a client can still verify
    every file it extracted.

WHAT CHANGED IN THE CLIENT

    `parascape/data_manager.py`

        HOST             the Hugging Face base address, `$PARASCAPE_DATA_HOST`
        ENV_URL          `$PARASCAPE_DATA_URL`, the whole address of this archive
        sources()        the addresses to try, in order: an explicit URL alone,
                         else the configured host, then the descriptor's `url`
                         and `fallback_url`
        _fetch()          one address, resuming the part file, split out of
                         `_download_locked` so a second host can be tried
        host_name()       what to call the host where a person reads it

      `_download_locked` tries each address in turn and reports every failure
      when all of them fail. The pinned checksums are unchanged and still
      checked after the download and again per file during the install, so
      trying a second host cannot widen what gets installed.

    `scripts/build_data_release.py` records both addresses in the descriptor,
    Hugging Face as `url` and the GitHub release as `fallback_url`.
    `parascape/data_manifest.json` was updated in place for 0.0.12.0 -- the two
    address fields only; size, archive checksum and all 10,929 file entries are
    untouched, so every installed copy still verifies.

    `scripts/check_wheel.py` now requires both addresses in a published wheel.

    `parascape/data_components.py` gained `hosts()` and
    `$PARASCAPE_COMPONENT_HOST`. The parts of every published component stay on
    their GitHub releases and both bases are accepted, so moving them is a
    setting plus a rebuilt index rather than another edit here.

WHAT MOVED ON 2026-11-13, SECOND PASS: THE COMPONENT PARTS

    All of it. 1,205 parts, 43.2 GB, every part of all four kinds, into the same
    dataset repository: `sparse_axes` 358 parts (0.9 GB), `workbook_cells` 83
    (1.1 GB), `deposited_values` 574 (8.0 GB), `sparse_coordinates` 190
    (33.3 GB). Smallest kind first, each kind finished before the next began.

        components/<kind>/<part name>     one part of one pack, at most 1 GiB
        components/<kind>/SHA256SUMS      the SHA-256 of every part of that kind
        components/README.md              what is here, per kind, and how to check it

    Part names are already unique across every kind and every release tag --
    1,205 distinct names for 1,205 parts -- so the name alone addresses a part
    and no tag has to appear in the layout.

    Each part was fetched from the public GitHub release asset the published
    index names, checked against the SHA-256 in that index before being
    uploaded, uploaded, and then checked again from what Hugging Face stored:
    the size and, for a part held as a large file, the SHA-256 the repository
    itself records, and for a small one the SHA-256 of a re-download. Only then
    was the staged copy deleted, so peak disk stayed around 3 GB. A part whose
    checksum did not match would not have been uploaded; none failed. Nothing
    unpublished was uploaded: every byte was already a public release asset.

NO NEW CORE RELEASE IS NEEDED

    The earlier note said moving the parts "rebuilds the index inside the core
    archive, i.e. a new release". It does not have to. The index does record a
    GitHub address per part, and the index lives inside an already-published
    core archive that cannot be rewritten -- but nothing forces the client to
    read the address from the index alone.

    `data_components.part_sources(component, part)` now builds the first address
    from three things the client already has: the configured host, the kind,
    which is the first segment of the component directory, and the part name,
    which is the last segment of the recorded address. The recorded address is
    the second one tried. So an installation of 0.0.12.0, with the index it
    already has, fetches a moved part from Hugging Face and an unmoved part from
    its GitHub release, with no new data release and no new application
    release. Moving a part is a setting, and `$PARASCAPE_COMPONENT_HOST` moves
    all of them at once.

WHAT CHANGED IN THE CLIENT, SECOND PASS

    `parascape/data_components.py`

        part_sources()        the addresses for one part, configured host first
        _download_part_once() tries each address in turn, resuming the same
                              part file, and checks the part checksum once,
                              after the loop
        _fetch_part()         one address, split out of the old function
        PartDownloadError     every address failed, naming each; `retryable` is
                              false when all of them refused the request, so a
                              part nobody has is not retried four times

      `PartDownloadError` is both an `OSError` and a `ValueError`, because
      callers already caught either one from the single-address version.

      The checksum contract is untouched. The part checksum is checked after
      the hosts have been tried, not between them, so wrong bytes from the first
      host end the install instead of being quietly replaced from the second.
      The joined archive checksum and every extracted file checksum are checked
      as before.

      The fetch is still `urllib` with a `Range` header and a resumed part file.
      `huggingface_hub` is not a dependency and must not become one.

    `scripts/build_data_components.py` still records a GitHub release address
    for a newly built part, which is now the fallback rather than the only
    address. A component published only on Hugging Face would carry a recorded
    address that answers 404 and still install from the first address; giving
    the builder the host explicitly is tidier and is left for the pass that next
    publishes a component.

STILL ON GITHUB

    - Nothing is only on GitHub any more. The identical copy of every component
      part stays published on tags v0.0.8.0, v0.0.9.0 and
      v0.0.9.0-quantities-data, deliberately, because an installation of an
      earlier ParaScape reads the address recorded in its own index. Deleting
      those assets is what actually frees the GitHub quota, and it is safe only
      once every installation in use builds the Hugging Face address first,
      which means after 0.0.13.0 ships. Not done, and not to be done quietly.
    - The identical core archive for every version, deliberately: an
      installation of 0.0.11.0 or earlier reads only its own descriptor, which
      names GitHub, and must keep working.
    - The release workflow still reads the core asset from the GitHub release
      with `gh release download`, so it needs no change today. It will need one
      when the GitHub copy stops being published; `.github/workflows/` was not
      touched.

WHAT TO DO NEXT

    1. Publish the instruction-53 organism and host gene and protein packs the
       same way, under `packs/<organism>/<version>/` in the same dataset
       repository, with the pack SHA-256 in the application rather than supplied
       by hand. `packs.download()` already takes any https address, so it needs
       a default host and the same two settings.
    2. Decide, after 0.0.13.0 has shipped and is in use, whether to delete the
       component assets from the GitHub releases. That is the step that frees
       the quota; until then both copies exist and the fallback is what keeps an
       older installation working.
    3. Add the Hugging Face upload to the release process so a release cannot
       ship a descriptor whose primary address is empty. Until then the manual
       step is in `docs/releases.md`.
    4. GEO-derived columns: ship the accession and the paper with the column, as
       decided today. Instruction 65 holds the source table; this is the
       publishing half of it.

HOW TO KNOW IT WORKED

    - `tests/test_data_download.py` and `tests/test_data_components.py` on the
      `starplast` environment. The new tests stub the fetch and touch no
      network: Hugging Face first and GitHub second; an unanswering first host
      falls through and still verifies; a failure on every host names every
      host; an address off the known hosts is refused before any request; the
      wrong bytes from the right host are still refused; `$PARASCAPE_DATA_HOST`
      and `$PARASCAPE_DATA_URL` move the address.
    - Measured against the live repository once, outside the suite: a ranged
      request for the last 1,024 bytes of the published archive returns 206 with
      `Content-Range: bytes 944322040-944323063/944323064`, through the redirect
      to the content host, so resuming an interrupted download works.
    - `python -c "from parascape import data_manager as D; print(D.sources())"`
      prints the Hugging Face address first and the GitHub one second.

    For the component parts, in `tests/test_data_components.py`, all stubbed and
    touching no network: the configured host is tried first and the recorded
    address second; a part on Hugging Face is taken from there and GitHub is
    never asked; a part still only on GitHub is taken from there after the first
    host answers 404; a failure on every host raises `PartDownloadError` naming
    every host and installs nothing; wrong bytes from the first host are refused
    with a checksum error and are not fetched again from the second;
    `$PARASCAPE_COMPONENT_HOST` moves where every part is fetched from.

    Measured once against the live repository, outside the suite, with the real
    client code and no stubs: one component of each of the four kinds installed
    and verified end to end, all four from Hugging Face. The same check, run
    before `sparse_coordinates` had moved, installed that one from its GitHub
    release asset after Hugging Face answered 404, so both paths are known to
    work. A ranged request for the last 1,024 bytes of a
    moved part returns 206 with the expected `Content-Range`, so a resumed part
    download works.

TRAPS

    - `huggingface_hub` is NOT in the `starplast` environment and is not a
      dependency of the application. It must not become one: the client fetches
      over `urllib`, which is what gives it resumable ranged downloads with no
      new dependency. The upload was run with the `spacr` environment python,
      which has the library and the token.
    - A Hugging Face `resolve` address redirects to a content host, and the
      final response comes from a hostname nobody can predict. Check the host of
      the address being requested, never the host that answered.
    - Do not add `?download=true` to a resolve address. `data_components.py`
      refuses an address with a query, and the plain address already serves the
      file.
    - The descriptor's `size` and `sha256` are the contract with every installed
      copy. Editing them means rebuilding and republishing the archive. The
      address fields are the only part of it that may be edited in place.
    - Uploading publishes under Einar's own account. Upload only what is already
      public, which for the core archive means the asset on the GitHub release
      of the same version. Nothing from the NAS and nothing from `results/`.
    - A part under about 10 MB is stored as an ordinary file rather than a large
      one, and then the repository records no SHA-256 for it. The uploader
      re-downloads those and hashes them, so every uploaded part is checked
      against the published index either way.
    - The copy of a component already extracted under
      `/media/carruthers/mnt1/parascape_data/components` is the unpacked pack,
      not the parts. The parts only exist as release assets, so moving them
      means fetching all of them from GitHub first, whatever is on this machine.
    - Do not point `$PARASCAPE_COMPONENT_HOST` at the GitHub releases base. The
      first address is then built as `<base>/<kind>/<part name>`, which is not
      where a release asset lives.
