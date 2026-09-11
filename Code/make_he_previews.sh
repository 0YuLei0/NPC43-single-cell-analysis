#!/usr/bin/env bash
# Build 2000 px PNG thumbnails for AICL FFPE_*.HE.tif files.
# ImageMagick CLI applies TIFF tiles/predictor correctly; do not read
# these 20k+ TIFFs into R.
#
# Usage (on a compute node, from or pointing at 03.HE):
#   bash make_he_previews.sh
#   bash make_he_previews.sh /home/ly385/project_pi_rf273/shared/ly385/AICL/03.HE 2000

set -euo pipefail

ROOT="${1:-.}"
MAX_PX="${2:-2000}"

if [[ ! -d "$ROOT" ]]; then
  echo "Not a directory: $ROOT" >&2
  exit 1
fi

if command -v magick >/dev/null 2>&1; then
  IM=(magick)
elif command -v convert >/dev/null 2>&1; then
  IM=(convert)
else
  echo "Need ImageMagick (magick or convert) on PATH" >&2
  exit 1
fi

mkdir -p "$ROOT/previews"
shopt -s nullglob
tifs=("$ROOT"/FFPE_*.HE.tif "$ROOT"/FFPE_*.HE.tiff "$ROOT"/FFPE_ALCL_*.HE.tif)
if [[ ${#tifs[@]} -eq 0 ]]; then
  echo "No FFPE_*.HE.tif files in $ROOT" >&2
  exit 1
fi

ok=0
skip=0
fail=0
for tif in "${tifs[@]}"; do
  base="$(basename "$tif")"
  stem="${base%.*}"
  out="$ROOT/previews/${stem}.png"
  if [[ -s "$out" ]]; then
    echo "SKIP $out"
    skip=$((skip + 1))
    continue
  fi
  echo "THUMB $base -> $out"
  if "${IM[@]}" "$tif" -resize "${MAX_PX}x" "$out"; then
    ok=$((ok + 1))
    ls -lh "$out"
  else
    echo "FAIL $tif" >&2
    fail=$((fail + 1))
  fi
done

echo "ok=$ok skip=$skip fail=$fail out=$ROOT/previews"
if [[ "$fail" -gt 0 ]]; then
  exit 1
fi
