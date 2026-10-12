#!/usr/bin/env bash
# Build a .deb. Run from the repo root:  bash packaging/build_debian.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(python3 -c "import tomllib,pathlib;print(tomllib.loads(pathlib.Path('$ROOT/pyproject.toml').read_text())['project']['version'])")"
STAGE="${STAGE_DIR:-$ROOT/build/deb}"          # stage off the root disk if it is small
PKG="$STAGE/parascape_${VERSION}_amd64"
rm -rf "$STAGE"; mkdir -p "$PKG/DEBIAN" "$PKG/opt" "$PKG/usr/share/applications"

pyinstaller --noconfirm --distpath "$STAGE/dist" --workpath "$STAGE/work" "$ROOT/packaging/parascape.spec"
cp -r "$STAGE/dist/parascape" "$PKG/opt/parascape"

cat > "$PKG/DEBIAN/control" <<EOF
Package: parascape
Version: ${VERSION}
Section: science
Priority: optional
Architecture: amd64
Depends: libgl1, libegl1, libxcb-cursor0
Maintainer: Einar Olafsson <einar.olafsson@gmail.com>
Description: A 3D browser for the Toxoplasma gondii knowledge map
 8,140 genes positioned by a UMAP embedding, twelve relation types kept separate,
 and one panel per gene showing what is actually known about it. Ships its own
 data cache, so it needs no network and no dataset at runtime.
EOF

cat > "$PKG/usr/share/applications/parascape.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=parascape
Comment=3D Toxoplasma gondii knowledge map
Exec=/opt/parascape/parascape
Categories=Science;Biology;
Terminal=false
EOF

dpkg-deb --build --root-owner-group "$PKG"
echo "built: ${PKG}.deb"
