#!/bin/sh
# Assemble the SWEVO reproducibility package (code + run records + generated outputs + README) as a zip file.
#
# Usage (from anywhere):   sh repro_package_build/make_package.sh [VERSION]        (default VERSION: 1.0.0-rc)
# Environment:  MANUSCRIPT_DIR=path   optionally add the LaTeX sources of the paper as manuscript/ (default: not added;
#                                      the CHECK-tag scan and the theory text check of regenerate_all.sh need them)
# Steps: build_inventory.py classifies every candidate file of the repository (inventory.csv) and lists the included
# ones (package_files.txt; analysis/rev3_* files present at this moment are included automatically); the files are
# copied with their repository paths into dist/stage/SWEVO_repro_<VERSION>/, together with the package-level files
# of this folder; MANIFEST.sha256 (sha256 of every file of the package except itself) is written into the package
# and copied to this folder; the zip goes to dist/SWEVO_repro_<VERSION>.zip with its sha256 in a side file.
# Caches (__pycache__, *.pyc), regen_work/ folders and the dist/ folder itself are never packaged.
# The repository is only read; this script writes only inside repro_package_build/.
set -eu
B=$(cd "$(dirname "$0")" && pwd)
ROOT=$(dirname "$B")
VERSION=${1:-1.0.0-rc}
NAME="SWEVO_repro_$VERSION"
python3 "$B/build_inventory.py" --root "$ROOT" --out "$B"
STAGE="$B/dist/stage/$NAME"
rm -rf "$B/dist/stage"
mkdir -p "$STAGE"
while IFS= read -r p; do
  [ -n "$p" ] || continue
  case "$p" in *__pycache__*|*.pyc) continue ;; esac
  mkdir -p "$STAGE/$(dirname "$p")"
  cp -p "$ROOT/$p" "$STAGE/$p"
done < "$B/package_files.txt"
for f in README.md requirements.txt requirements-pywake.txt pip-freeze-full.txt Dockerfile CITATION.cff \
         LICENSE_SUGGESTION.md regenerate_all.sh compare_outputs.py verify_determinism.py count_records.py \
         build_inventory.py make_package.sh inventory.csv; do
  cp -p "$B/$f" "$STAGE/$f"
done
[ -d "$B/verification" ] && cp -R "$B/verification" "$STAGE/verification"
if [ -n "${MANUSCRIPT_DIR:-}" ]; then
  cp -R "$MANUSCRIPT_DIR" "$STAGE/manuscript"
  find "$STAGE/manuscript" \( -name '*.aux' -o -name '*.log' -o -name '*.out' -o -name '*.synctex.gz' \) -delete
fi
find "$STAGE" \( -name __pycache__ -o -name regen_work \) -prune -exec rm -rf {} +
( cd "$STAGE" && find . -type f ! -name MANIFEST.sha256 | sed 's|^\./||' | LC_ALL=C sort | \
  while IFS= read -r f; do sha256sum "$f"; done > MANIFEST.sha256 )
cp "$STAGE/MANIFEST.sha256" "$B/MANIFEST.sha256"
rm -f "$B/dist/$NAME.zip" "$B/dist/$NAME.zip.sha256"
( cd "$B/dist/stage" && find "$NAME" -type f | LC_ALL=C sort | zip -X -q -9 "$B/dist/$NAME.zip" -@ )
( cd "$B/dist" && sha256sum "$NAME.zip" > "$NAME.zip.sha256" )
NF=$(wc -l < "$STAGE/MANIFEST.sha256")
echo "package: $B/dist/$NAME.zip"
echo "  files: $NF (+ MANIFEST.sha256), unpacked $(du -sh "$STAGE" | cut -f1), zip $(du -h "$B/dist/$NAME.zip" | cut -f1)"
echo "  sha256: $(cut -d' ' -f1 "$B/dist/$NAME.zip.sha256")"
echo "  sensitivity-study files (rev3_*) included: $(grep -c '^analysis/rev3_' "$B/package_files.txt" || true)"
rm -rf "$B/dist/stage"
