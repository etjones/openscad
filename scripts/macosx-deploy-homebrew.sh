#!/usr/bin/env bash
#
# Turn a CMake build of OpenSCAD.app made against Homebrew dependencies
# into a self-contained, ad-hoc signed bundle and a .dmg.
#
# Usage: scripts/macosx-deploy-homebrew.sh <build-dir> [<output-dir>]
#
# release-common.sh covers the from-source dependency prefix; this is the
# equivalent for a Homebrew (single-architecture) build: macdeployqt pulls
# in Qt and every Homebrew dylib, then leftover Homebrew rpaths and install
# names are scrubbed so nothing can resolve outside the bundle.
#
set -euo pipefail

BUILDDIR=${1:?build directory}
OUT=${2:-$PWD/release-mac}
VERSION=$(date +%Y.%m.%d)
QTBIN=$(brew --prefix qt)/bin

rm -rf "$OUT"
mkdir -p "$OUT"
cp -R "$BUILDDIR/OpenSCAD.app" "$OUT/OpenSCAD.app"
cd "$OUT"
/usr/libexec/PlistBuddy -c "Set :CFBundleVersion $VERSION" OpenSCAD.app/Contents/Info.plist

echo "== macdeployqt"
"$QTBIN/macdeployqt" OpenSCAD.app -no-strip 2>&1 | grep -iE 'error' || true

echo "== scrubbing Homebrew paths"
for f in $(find OpenSCAD.app -type f \( -perm -u+x -o -name '*.dylib' \)); do
  file "$f" | grep -q Mach-O || continue
  for rp in $(otool -l "$f" | awk '/LC_RPATH/{r=1} r&&/path /{print $2; r=0}' | grep '^/opt/homebrew' || true); do
    install_name_tool -delete_rpath "$rp" "$f" 2>/dev/null || true
  done
  case "$f" in
    *.dylib)
      id=$(otool -D "$f" | tail -1)
      case "$id" in /opt/homebrew*) install_name_tool -id "@executable_path/../Frameworks/$(basename "$f")" "$f";; esac
      ;;
  esac
done
LEFT=$(find OpenSCAD.app -type f \( -perm -u+x -o -name '*.dylib' \) -exec otool -L {} + 2>/dev/null | grep -c '/opt/homebrew' || true)
echo "references to /opt/homebrew left: $LEFT"

echo "== ad-hoc signing"
codesign --force --deep --sign - OpenSCAD.app
codesign --verify --deep --strict OpenSCAD.app && echo "signature ok"

echo "== dmg"
rm -f "OpenSCAD-$VERSION.dmg"
hdiutil create -volname "OpenSCAD $VERSION" -srcfolder OpenSCAD.app -ov -format UDZO "OpenSCAD-$VERSION.dmg" | tail -1
du -sh OpenSCAD.app "OpenSCAD-$VERSION.dmg"
