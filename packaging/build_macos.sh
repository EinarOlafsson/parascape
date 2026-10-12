#!/usr/bin/env bash
# Build parascape.app and a .dmg.  bash packaging/build_macos.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(python3 -c "import tomllib,pathlib;print(tomllib.loads(pathlib.Path('$ROOT/pyproject.toml').read_text())['project']['version'])")"
rm -rf "$ROOT/build/mac"; mkdir -p "$ROOT/build/mac"
pyinstaller --noconfirm --distpath "$ROOT/build/mac/dist" --workpath "$ROOT/build/mac/work" \
            "$ROOT/packaging/parascape.spec"
APP="$ROOT/build/mac/dist/parascape.app"
[ -d "$APP" ] || { echo "no .app produced"; exit 1; }
hdiutil create -volname "parascape ${VERSION}" -srcfolder "$APP" -ov -format UDZO \
        "$ROOT/build/mac/parascape-${VERSION}.dmg"
echo "built: build/mac/parascape-${VERSION}.dmg"
