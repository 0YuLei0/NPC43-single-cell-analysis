#!/usr/bin/env bash
# Rename each FFPE_LN sample FASTQ pair to <folder>_1.fq.gz and <folder>_2.fq.gz.
#
# Usage (from the directory that contains the FFPE_LN* folders):
#   bash rename_ffpe_ln_fastq.sh            # dry-run
#   bash rename_ffpe_ln_fastq.sh --apply    # actually rename
#
# Or pass the parent directory:
#   bash rename_ffpe_ln_fastq.sh /path/to/fastq_root --apply

set -euo pipefail

usage() {
  cat <<'EOF'
Rename FFPE_LN FASTQ pairs to <folder_name>_1.fq.gz and <folder_name>_2.fq.gz.

Usage:
  rename_ffpe_ln_fastq.sh [PARENT_DIR] [--apply]

Default is a dry-run. Pass --apply to perform the renames.
PARENT_DIR defaults to the current directory.
EOF
}

APPLY=0
PARENT="."

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --apply)
      APPLY=1
      shift
      ;;
    --dry-run)
      APPLY=0
      shift
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      PARENT="$1"
      shift
      ;;
  esac
done

if [[ ! -d "$PARENT" ]]; then
  echo "Parent directory does not exist: $PARENT" >&2
  exit 1
fi

# folder|original_r1|original_r2
MAPPINGS=(
  "FFPE_LN14734|LN14734_L1_1.fq.gz|LN14734_L1_2.fq.gz"
  "FFPE_LN18427|LN18427_CKDL250013248-1A_22T5VTLT4_L5_1.fq.gz|LN18427_CKDL250013248-1A_22T5VTLT4_L5_2.fq.gz"
  "FFPE_LN21720|LN21720_CKDL250001098-1A_22M5YFLT4_L2_1.fq.gz|LN21720_CKDL250001098-1A_22M5YFLT4_L2_2.fq.gz"
  "FFPE_LN22630|LN22630_CKDL250013247-1A_22T5WYLT4_L4_1.fq.gz|LN22630_CKDL250013247-1A_22T5WYLT4_L4_2.fq.gz"
  "FFPE_LN32298|NA_CKDL250010655-1A_22THYKLT4_L5_1.fq.gz|NA_CKDL250010655-1A_22THYKLT4_L5_2.fq.gz"
  "FFPE_LN3470|LN3470_L1_1.fq.gz|LN3470_L1_2.fq.gz"
  "FFPE_LN3524|NA_2_CKDL250010656-1A_22THYKLT4_L2_1.fq.gz|NA_2_CKDL250010656-1A_22THYKLT4_L2_2.fq.gz"
  "FFPE_LN4737|LN4737_L1_1.fq.gz|LN4737_L1_2.fq.gz"
  "FFPE_LN5034|LN5034_CKDL250001099-1A_22M5YFLT4_L3_1.fq.gz|LN5034_CKDL250001099-1A_22M5YFLT4_L3_2.fq.gz"
)

errors=0
renamed=0
skipped=0

rename_one() {
  local src="$1"
  local dst="$2"

  if [[ "$src" == "$dst" ]]; then
    echo "SKIP already named: $dst"
    skipped=$((skipped + 1))
    return 0
  fi

  if [[ ! -e "$src" ]]; then
    if [[ -e "$dst" ]]; then
      echo "SKIP already renamed: $dst"
      skipped=$((skipped + 1))
      return 0
    fi
    echo "ERROR missing source: $src" >&2
    errors=$((errors + 1))
    return 1
  fi

  if [[ -e "$dst" ]]; then
    echo "ERROR destination exists: $dst" >&2
    errors=$((errors + 1))
    return 1
  fi

  if [[ "$APPLY" -eq 1 ]]; then
    mv -- "$src" "$dst"
    echo "RENAMED $src -> $dst"
  else
    echo "DRY-RUN mv -- $src $dst"
  fi
  renamed=$((renamed + 1))
}

for row in "${MAPPINGS[@]}"; do
  IFS='|' read -r folder r1 r2 <<<"$row"
  dir="$PARENT/$folder"
  if [[ ! -d "$dir" ]]; then
    echo "ERROR missing folder: $dir" >&2
    errors=$((errors + 1))
    continue
  fi
  rename_one "$dir/$r1" "$dir/${folder}_1.fq.gz" || true
  rename_one "$dir/$r2" "$dir/${folder}_2.fq.gz" || true
done

echo
echo "renamed=$renamed skipped=$skipped errors=$errors apply=$APPLY"

if [[ "$errors" -gt 0 ]]; then
  exit 1
fi
