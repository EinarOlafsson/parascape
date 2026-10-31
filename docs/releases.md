# Development and releases

## Branches and commits

Develop on `nightly`; `main` is the release branch. Both branches run the regression
checks and build the documentation. Only `main` deploys the public documentation
and can publish packages. Merge completed, checked changes into `main` promptly.
A merge without a version increase does not publish a package.

Commit as Einar Olafsson using `einar.olafsson@gmail.com`, without additional
co-author trailers. Use a clear subject and a body describing the reason for the
change, resulting behaviour, and relevant validation. Automated publishing creates
no commits and does not add a code contributor.

## Local checks

Use a virtual environment and install the development and documentation tools:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,docs]"
python scripts/release.py check
QT_QPA_PLATFORM=offscreen OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 pytest -q
python scripts/build_docs.py
```

On Windows, activate `.venv\Scripts\activate` instead. On Linux, Qt may need
`libegl1`, `libgl1`, `libdbus-1-3`, `libxkbcommon-x11-0`, and `libxcb-cursor0` from
the system package manager. `QT_QPA_PLATFORM=offscreen` is for tests, not normal use.

The full test suite also exercises source-dataset ingestion, rendering, GPU
backends, and long searches. Some checks skip without their source data, optional
dependencies or hardware. GitHub checks run all non-slow tests and a real OpenGL
species-switch test. On Linux, install `xvfb` and `xauth`, then run the rendering
check separately:

```bash
QT_QPA_PLATFORM=xcb xvfb-run -a -s '-screen 0 1920x1200x24' pytest tests/test_app_contexts.py -q
```

These are software checks. The [0.43 benchmark](benchmark-0.43.md) separately
reports the scientific comparisons and their limitations. Before a version bump,
require passing checks on the final implementation on `nightly`, build the docs,
and verify a wheel installation outside the source checkout.

## Package layout

There is one PyPI distribution: **`parascape`**. It contains the Python modules,
console commands, a small data release descriptor, and icons. `parascape[gpu]` enables optional
CUDA dependencies and `parascape[ingest]` enables coverage-file imports; both are
extras of the same project.

`parascape[sequence]` adds the encoder dependencies for regenerating ESM features;
`parascape[structures]` adds local structure parsing. Reading the bundled feature
tables requires neither extra.

Scientific data is hosted as a versioned GitHub release asset. First launch offers
a download when files are missing; installation verifies archive and per-file SHA-256
checksums, stages extraction, and preserves existing data on failure. Once installed,
browsing is offline. `parascape-data --download` supports headless installations.
Original raw source datasets are not mirrored; see the [dataset catalogue](datasets.md).

## One-time PyPI setup

Sign in to **Einar's own PyPI account** at
[PyPI publishing](https://pypi.org/manage/account/publishing/) and create
a pending Trusted Publisher for **`parascape`** with these values:

| Field | Value |
|---|---|
| PyPI project name | `parascape` |
| Owner | `EinarOlafsson` |
| Repository | `parascape` |
| Workflow filename | `release.yml` |
| Environment | `pypi` |

For an existing project, add the same publisher in that project's Publishing
settings. See [PyPI's setup instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).
The workflow exchanges its GitHub identity for short-lived upload credentials;
no permanent PyPI token is required in GitHub secrets.
The account that registers a pending publisher owns the resulting PyPI project.
The `Owner` field above identifies the GitHub repository owner; it is not a PyPI
username. Do not create the project under a different account.

Create the GitHub environment `pypi` and restrict deployments to the `main` branch.
To publish automatically, leave required reviewers unset. Under repository
**Settings → Pages**, choose **GitHub Actions** as the source. The documentation
workflow builds pull requests and deploys pushes to `main`.

## Publish a version

```bash
git switch nightly
git pull --ff-only origin nightly
python scripts/release.py bump 0.43.0
python scripts/release.py check
```

Use a version greater than the current one. The bump command updates `pyproject.toml` and `parascape.__version__`.
Update `CHANGELOG.md`, run the release checks, commit with a descriptive message,
and push to `nightly`. Merge `nightly` into `main` through a pull request or locally:

```bash
git push origin nightly
git switch main
git pull --ff-only origin main
git merge --ff-only nightly
git push origin main
git switch nightly
```

If the branches have diverged, resolve the merge on `nightly` or use a pull request;
do not force-push `main`. A separate version-bump commit on `main` also triggers
the release, but preparing it on `nightly` keeps the normal development flow intact.

The release workflow compares the version with the repository state before the
push. An unchanged version does not publish; a downgrade fails. A version increase
runs the checks, builds a wheel and source distribution, verifies their contents,
and publishes `parascape`. After a successful PyPI upload, it creates
a GitHub release and `v<version>` tag at the exact commit that was built, with
generated release notes and the distribution files attached. Prerelease versions
are marked as prereleases on GitHub. Existing releases and their assets are left
intact on retries. Tags do not need to be created manually.

For the first upload or to retry a partially completed upload, manually run
**Publish Python packages** on `main`. Existing files are skipped; PyPI does not
allow an uploaded version to be replaced. Use a new version for changed code.
If publication failed before any upload, rerunning the original job preserves the
same tested artifacts.

## README on PyPI

`README.md` is the source for both the GitHub README and the PyPI description.
Package builds include it directly through the `readme` setting in `pyproject.toml`;
there is no separate PyPI README to maintain. Use absolute URLs for images and
repository files so links work on both sites. The linked slide cover uses a
standard image that renders on both GitHub and PyPI.

PyPI receives the README with each release. A README edit on GitHub appears on
PyPI when the next version is published; it does not change existing releases.

## Inspect release artifacts

```bash
python -m build
python scripts/check_wheel.py dist
python -m twine check --strict dist/*
```

Start with an empty `dist/` directory so files from older versions are not included.
The wheel check requires code, icons and the versioned GitHub data descriptor.
It rejects bundled scientific datasets, saved embeddings, unconditional CUDA
dependencies and files above PyPI’s default upload limit. The installed-wheel
smoke test verifies and installs the pinned GitHub asset before loading gene data.

The desktop installers in `packaging/` are separate from PyPI. Their existing
workflow runs manually or on a `v*` tag.


## Versioned GitHub data assets

After the version bump, build the data asset with `python scripts/build_data_release.py`.
Commit the generated `parascape/data_manifest.json` with the source changes on nightly.
After checks pass, publish `data-dist/parascape-data-VERSION.tar.gz` on the GitHub
release `vVERSION` at that exact commit **before** promoting to main. The release
workflow downloads this existing asset and verifies it against the descriptor while
smoke-testing the installed wheel. It does not rebuild an archive whose compression
bytes could vary across platforms. PyPI publication happens only after this check;
the final job attaches wheels and source distributions to the existing release.

## Optional large-data components (development)

The development installer in `parascape.data_components` supports independently
installed sparse matrix and axis components. This is the storage foundation for
the next release; the current 0.0.7.0 release continues using its existing core
archive. Optional assets are not public merely because they have been prepared.
Development builds connect this index to the sparse browser's missing-dataset
prompt. Publication still requires verified hosted assets and the next core bundle.

A core-data `data_components.json` index has `schema_version: 1`, a version,
`components`, and a dataset-to-component dependency map. Component keys are
`sparse_coordinates/<pack-sha>` or `sparse_axes/<pack-sha>`. Each component declares
its original files with exact sizes and SHA-256 checksums, one complete archive
checksum/size, and an ordered list of GitHub release URLs with individual part
sizes/checksums. Every part is at most 1 GiB. Shared axes use the same component
identity across datasets.

Prepare an existing independently validated pack without changing its contents:

```bash
python scripts/build_data_components.py /path/to/original/pack \
  sparse_coordinates/ORIGINAL_PACK_SHA /path/to/staging --tag v0.0.8.0
```

The builder retains original manifests, validation receipts and Parquet bytes.
It creates deterministic tar headers, streams archive parts without allocating a
whole archive in memory, verifies every source file while copying, and records
builder/source signatures. Repeated preparation reuses checked immutable parts.
This command does not upload files or publish a release. Upload gates must verify
all assets against the index before including that index in a public core bundle.

`data_components.read_index()` reads an index from verified core data.
`data_components.download(component)` installs one component in a shared user
cache, overridable with `PARASCAPE_COMPONENTS`. The client verifies each part,
the reconstructed archive and every extracted file before atomically exposing
the component directory. Cancellation retains partial/completed download pieces
for a later verified retry. Traversal, links, unexpected members, bad ranges and
checksum mismatches are refused. Existing unmanaged targets are preserved;
replaced managed components are retained for recovery. Core data and other
component directories are not replaced.

Contributors can read an installed component explicitly with
`SparseDataset(..., coordinate_root=storage_root() / 'sparse_coordinates',
axis_root=storage_root() / 'sparse_axes')`. The sparse browser resolves these roots automatically when a published index
declares the dataset, while retaining complete local collections and explicit
root overrides. A provisional index cannot trigger remote downloads. Keep download/installation receipts
separate from scientific normalization, identifier mapping and biological QC;
copying a verified component does not establish cross-study comparability.


A public optional index must set `publication_status: published_verified` only
after all of its hosted assets have passed the release verification gates. Local
staging indices deliberately retain a provisional publication status. Tests may
use a publication-marker fixture with local transport; that is not evidence of a
GitHub upload. Preserve that distinction in release receipts.
