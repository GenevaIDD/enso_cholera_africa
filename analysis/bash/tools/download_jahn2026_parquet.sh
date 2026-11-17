#!/bin/bash
# tools/download_jahn2026_parquet.sh — fetch the Jahn et al. (2026) admin-level
# CHIRPSv3 / ERA5-Land parquet observations for Africa from the public Box share.
#
# Usage:
#   bash analysis/bash/tools/download_jahn2026_parquet.sh                     # all countries
#   bash analysis/bash/tools/download_jahn2026_parquet.sh --iso3 LSO,STP,COM  # subset
#   bash analysis/bash/tools/download_jahn2026_parquet.sh --dry-run           # plan + bytes
#   bash analysis/bash/tools/download_jahn2026_parquet.sh --outdir /scratch/cp
#
# One parquet file per admin unit, 1981-01-01..2024-12-31; it covers gids the
# climate_stats.sqlite lacks. Lands in data/climate_parquet/ (not
# data/climate_data/), so 02_preprocess_climate.R can read either source
# (--climate_source csv|parquet).
#
# Box rejects HTTP Range requests: archives download whole to a .part file and
# are promoted only when the byte count matches the manifest, so there is no
# --resume (a country is extracted or refetched). The share is a browser
# endpoint; the User-Agent below matters.
#
# Steps:
#   1. Read analysis/data_manifests/jahn2026_africa_zips.tsv (iso3, file_id, bytes).
#   2. Skip extracted countries; download, size-check and unzip the rest.
#   3. Append iso3 / bytes / sha256 / time to <outdir>/DOWNLOAD_MANIFEST.tsv,
#      which 02_preprocess_climate.R stamps into climate_db_provenance.yml.

set -euo pipefail

MANIFEST="analysis/data_manifests/jahn2026_africa_zips.tsv"
OUTDIR="data/climate_parquet"
SHARED_NAME="y4f0ywor4gvcrdxkewlzerfxw4lik4ts"
BASE_URL="https://imperialcollegelondon.app.box.com/index.php"
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"

ISO3_FILTER=""
DRY_RUN=false
REFETCH=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --iso3)     ISO3_FILTER="$2"; shift 2 ;;
    --outdir)   OUTDIR="$2";      shift 2 ;;
    --manifest) MANIFEST="$2";    shift 2 ;;
    --dry-run)  DRY_RUN=true;     shift   ;;
    --refetch)  REFETCH=true;     shift   ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 [--iso3 A,B,C] [--outdir DIR] [--manifest TSV] [--dry-run] [--refetch]" >&2
      exit 1 ;;
  esac
done

if [[ ! -f "$MANIFEST" ]]; then
  echo "Error: manifest not found: $MANIFEST" >&2
  exit 1
fi

command -v curl  >/dev/null || { echo "Error: curl not found."  >&2; exit 1; }
command -v unzip >/dev/null || { echo "Error: unzip not found." >&2; exit 1; }

# --iso3 LSO,STP -> a newline-delimited set we can grep against.
WANTED=""
if [[ -n "$ISO3_FILTER" ]]; then
  WANTED=$(echo "$ISO3_FILTER" | tr ',' '\n' | tr '[:lower:]' '[:upper:]' | sed '/^$/d')
fi

# ---------------------------------------------------------------------------
# 1. Parse the manifest into parallel arrays, applying --iso3.
# ---------------------------------------------------------------------------
ISO3S=(); FILE_IDS=(); SIZES=()
TOTAL_BYTES=0

while IFS=$'\t' read -r iso3 file_id bytes; do
  # Skip comments, blanks and the column-header row.
  [[ -z "${iso3:-}" || "$iso3" == \#* || "$iso3" == "iso3" ]] && continue
  if [[ -n "$WANTED" ]] && ! grep -qxF "$iso3" <<< "$WANTED"; then
    continue
  fi
  ISO3S+=("$iso3"); FILE_IDS+=("$file_id"); SIZES+=("$bytes")
  TOTAL_BYTES=$(( TOTAL_BYTES + bytes ))
done < "$MANIFEST"

if [[ ${#ISO3S[@]} -eq 0 ]]; then
  echo "Error: no manifest rows selected (--iso3 '$ISO3_FILTER' matched nothing)." >&2
  exit 1
fi

# Warn about an --iso3 value that is not on Box at all (GAB, BWA, GNQ, LBY, MUS, SYC)
# rather than silently downloading fewer countries than asked for.
if [[ -n "$WANTED" ]]; then
  while read -r w; do
    [[ -z "$w" ]] && continue
    if ! printf '%s\n' "${ISO3S[@]}" | grep -qxF "$w"; then
      echo "WARNING: $w is not in $MANIFEST (not published on Box) — skipping." >&2
    fi
  done <<< "$WANTED"
fi

echo "--- Jahn et al. (2026) parquet download"
echo "Manifest  : $MANIFEST"
echo "Outdir    : $OUTDIR"
echo "Countries : ${#ISO3S[@]}"
printf 'Total size: %.2f GB (%d bytes)\n' "$(echo "$TOTAL_BYTES" | awk '{print $1/1e9}')" "$TOTAL_BYTES"

if [[ "$DRY_RUN" == true ]]; then
  echo "--- --dry-run: nothing downloaded. Would fetch:"
  for i in "${!ISO3S[@]}"; do
    state="download"
    [[ -d "$OUTDIR/${ISO3S[$i]}" ]] && state="SKIP (already extracted)"
    printf '  %-4s %12s bytes  %s\n' "${ISO3S[$i]}" "${SIZES[$i]}" "$state"
  done
  exit 0
fi

# ---------------------------------------------------------------------------
# 2. Download, verify, extract.
# ---------------------------------------------------------------------------
mkdir -p "$OUTDIR/zips"
DL_MANIFEST="$OUTDIR/DOWNLOAD_MANIFEST.tsv"
if [[ ! -f "$DL_MANIFEST" ]]; then
  printf 'iso3\tbytes\tsha256\tn_parquet\tfetched_at\n' > "$DL_MANIFEST"
fi

n_ok=0; n_skip=0; FAILED=()

for i in "${!ISO3S[@]}"; do
  iso3="${ISO3S[$i]}"; file_id="${FILE_IDS[$i]}"; want_bytes="${SIZES[$i]}"
  dest_dir="$OUTDIR/$iso3"

  # Presence-checking IS the resume mechanism (see CONSTRAINTS above).
  if [[ -d "$dest_dir" && "$REFETCH" == false ]]; then
    n_parquet=$(find "$dest_dir" -name '*.parquet.gzip' -type f | wc -l)
    if [[ "$n_parquet" -gt 0 ]]; then
      printf '[%2d/%2d] %-4s already extracted (%s parquet files) — skipping\n' \
        "$((i+1))" "${#ISO3S[@]}" "$iso3" "$n_parquet"
      n_skip=$((n_skip + 1))
      continue
    fi
    echo "WARNING: $dest_dir exists but holds no parquet files — refetching." >&2
    rm -rf "$dest_dir"
  fi
  [[ "$REFETCH" == true ]] && rm -rf "$dest_dir"

  zip_path="$OUTDIR/zips/$iso3.zip"
  part_path="$zip_path.part"
  url="$BASE_URL?rm=box_download_shared_file&shared_name=$SHARED_NAME&file_id=f_$file_id"

  printf '[%2d/%2d] %-4s downloading %s bytes ...\n' \
    "$((i+1))" "${#ISO3S[@]}" "$iso3" "$want_bytes"

  rm -f "$part_path"
  if ! curl -sSL --fail --retry 3 --retry-delay 5 -A "$UA" -o "$part_path" "$url"; then
    echo "  ERROR: download failed for $iso3" >&2
    rm -f "$part_path"
    FAILED+=("$iso3")
    continue
  fi

  got_bytes=$(stat -c %s "$part_path")
  if [[ "$got_bytes" != "$want_bytes" ]]; then
    # Box serves an HTML error page with HTTP 200 when a share link goes stale,
    # so a size mismatch is the only reliable detector. Do NOT unzip it.
    echo "  ERROR: size mismatch for $iso3 — expected $want_bytes, got $got_bytes." >&2
    echo "         Box may have re-uploaded the archive; regenerate the manifest." >&2
    rm -f "$part_path"
    FAILED+=("$iso3")
    continue
  fi

  sha=$(sha256sum "$part_path" | cut -d' ' -f1)
  mv "$part_path" "$zip_path"

  if ! unzip -q -o "$zip_path" -d "$OUTDIR"; then
    echo "  ERROR: unzip failed for $iso3" >&2
    FAILED+=("$iso3")
    continue
  fi
  rm -f "$zip_path"

  if [[ ! -d "$dest_dir" ]]; then
    # Every archive inspected so far wraps its files in an <ISO3>/ folder; if one
    # does not, the reader's directory-per-country assumption breaks silently.
    echo "  ERROR: $iso3.zip did not produce $dest_dir — unexpected archive layout." >&2
    FAILED+=("$iso3")
    continue
  fi

  n_parquet=$(find "$dest_dir" -name '*.parquet.gzip' -type f | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$iso3" "$want_bytes" "$sha" "$n_parquet" "$(date -Is)" >> "$DL_MANIFEST"
  printf '         %-4s OK — %s parquet files\n' "$iso3" "$n_parquet"
  n_ok=$((n_ok + 1))
done

rmdir "$OUTDIR/zips" 2>/dev/null || true

echo "--- Done: $n_ok downloaded, $n_skip already present, ${#FAILED[@]} failed"
echo "--- Manifest: $DL_MANIFEST"

if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "--- FAILED: ${FAILED[*]}" >&2
  echo "    Do not run 02_preprocess_climate.R --climate_source parquet until these" >&2
  echo "    resolve: a country with no parquet directory is reported as a missing" >&2
  echo "    gid in generated_data/qc_climate_gids_missing.csv, which is" >&2
  echo "    indistinguishable from a country that is genuinely not published." >&2
  exit 1
fi
