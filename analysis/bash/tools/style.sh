#!/bin/bash
# tools/style.sh — format the repo's R code with air (settings: air.toml at the root).
# Run from the project root. air must be on PATH (https://posit-dev.github.io/air/).
#
# Usage:
#   bash analysis/bash/tools/style.sh --check [PATH...]   # report files air would change, change nothing
#   bash analysis/bash/tools/style.sh --fix   [PATH...]   # rewrite them in place
# PATH defaults to the whole repo (air.toml's exclude list applies).

set -euo pipefail

mode="${1:-}"
case "$mode" in
  --check) shift; args=(--check) ;;
  --fix)   shift; args=() ;;
  *) echo "usage: $0 --check|--fix [PATH...]" >&2; exit 2 ;;
esac

if ! command -v air >/dev/null 2>&1; then
  echo "style.sh: air is not on PATH (https://posit-dev.github.io/air/)" >&2
  exit 127
fi
[[ -f air.toml ]] || { echo "style.sh: run from the project root (no air.toml here)" >&2; exit 2; }

(( $# )) || set -- .
air format "${args[@]}" "$@"
