#!/bin/bash
# slurm/submit_manuscript.sh — SBATCH script for the manuscript
# figures (phase 10) and the key-numbers report (phase 11) as one SLURM job.
#
# Runs, in order:
#   0. package check (R packages used by 10 / 11 and their sourced utils
#      closure; taxdat and the taxdat internals 10 calls; sf's GDAL/GEOS/PROJ
#      linkage; quarto on PATH). Any missing hard requirement stops the job.
#   1. Rscript analysis/10_manuscript_figures.R -C <config> [--figures] [--supplementary]
#   2. quarto render analysis/11_manuscript_numbers.qmd -P config_path:<config>
#
# Both read only the phase-09 artefacts (results_require()), so submit this
# after stages/09_postprocess_results.sh. Normally submitted by the stage
# driver stages/10_manuscript_figures.sh (--after <RESULTS_JOB_IDS>); direct
# sbatch still works.
#
# Usage (from the repo root):
#   sbatch analysis/bash/slurm/submit_manuscript.sh
#   sbatch analysis/bash/slurm/submit_manuscript.sh --config analysis/configs/results/manuscript.yml
#   sbatch analysis/bash/slurm/submit_manuscript.sh --figures 4,H --supplementary false
#   sbatch analysis/bash/slurm/submit_manuscript.sh --skip-10      # report only
#   sbatch analysis/bash/slurm/submit_manuscript.sh --skip-11      # figures only
#   bash   analysis/bash/slurm/submit_manuscript.sh --check-only   # login node: package check, no run
#
#SBATCH --job-name=manuscript_10_11
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH -c 1

set -euo pipefail

# The MANUSCRIPT stack of the cluster profile (lib/cluster.sh): R with sf / GDAL (as OutbreakExtractR and slurm/submit_results.sh; 10
# needs sf: GADM / giscoR basemaps, taxdat lakes and rivers) plus quarto for 11
# (non-fatal: the quarto check below reports a miss, and --skip-11 does not
# need it). A profile without modules makes this a no-op, so --check-only also
# runs anywhere.
source analysis/bash/lib/cluster.sh   # cwd = repo root (SLURM runs a copy of this script)
cholera_env --purge MANUSCRIPT

[[ "${RESULTS_ALLOW_STALE_INPUTS:-false}" == "true" ]] && \
  echo "RESULTS_ALLOW_STALE_INPUTS=true: 10 / 11 read artefacts even if their inputs changed since 09 (stamp + decision checks still enforced)" >&2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

CONFIG="analysis/configs/results/manuscript.yml"
FIGURES=""
SUPPLEMENTARY=""
SKIP_10=false
SKIP_11=false
CHECK_ONLY=false
QUARTO="${QUARTO:-quarto}"   # override with QUARTO=/path/to/quarto if it is not on PATH

usage() {
  echo "Usage: $0 [--config <results cfg.yml>] [--figures 1,2,3,4,H] [--supplementary true|false] [--skip-10] [--skip-11] [--check-only]" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config|-C)     CONFIG="$2";        shift 2 ;;
    --figures)       FIGURES="$2";       shift 2 ;;
    --supplementary) SUPPLEMENTARY="$2"; shift 2 ;;
    --skip-10)       SKIP_10=true;       shift   ;;
    --skip-11)       SKIP_11=true;       shift   ;;
    --check-only)    CHECK_ONLY=true;    shift   ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ ! -f "$CONFIG" ]]; then echo "Error: config not found: $CONFIG (run from the repo root)" >&2; exit 1; fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found on PATH (the profile's MANUSCRIPT stack failed?)." >&2; exit 1
fi

echo "***** Beginning of script *****"
date
echo "Config : ${CONFIG}"
echo "Rscript: $(command -v Rscript)"

# ---------------------------------------------------------------------------
# 0. Package check
# ---------------------------------------------------------------------------
echo "--- Checking required R packages ..."
NEED_QUARTO=true
[[ "$SKIP_11" == "true" ]] && NEED_QUARTO=false

Rscript --vanilla - <<'RS' || { echo "Error: package check failed (see above)." >&2; exit 1; }
cat("R ", R.version.string, "\n", sep = "")
cat(".libPaths():\n"); cat(paste0("  ", .libPaths()), sep = "\n")

# Hard: attached by 10 / 11, or called via `::` on their code paths
# (10 + utils/{results_core,figures_cb,figures_io,deploy_core,figures_hindcast},
#  11 + utils/results_report, and the closure results_core / deploy_core source).
hard <- c(
  # 10_manuscript_figures.R
  "optparse", "ggplot2", "tibble", "tidyr", "readr", "purrr", "dplyr",
  "stringr", "lubridate", "qs2", "sf", "scales", "patchwork", "ragg",
  "taxdat", "systemfonts", "cowplot", "forcats", "ggh4x", "ggnewscale",
  "ggspatial", "gtable", "rmapshaper", "janitor", "giscoR",
  # sourced utils closure
  "yaml", "glue",
  # 11_manuscript_numbers.qmd (+ quarto's knitr engine)
  "knitr", "kableExtra", "rmarkdown"
)
# Soft: referenced inside sourced utils functions that 10 / 11 do not call
# today, or optional in 11 (gtsummary: Section 1's descriptive tables and
# tests, skipped when missing); a miss is reported but does not stop the job.
soft <- c("brms", "dlnm", "ncdf4", "nngeo", "digest", "zoo", "png", "nanoparquet", "gtsummary")

have <- function(p) requireNamespace(p, quietly = TRUE)
ver  <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)

miss_hard <- hard[!vapply(hard, have, logical(1))]
miss_soft <- soft[!vapply(soft, have, logical(1))]
for (p in setdiff(hard, miss_hard)) cat(sprintf("  [ok]   %-12s %s\n", p, ver(p)))
for (p in miss_hard)                cat(sprintf("  [MISS] %-12s (required)\n", p))
for (p in miss_soft)                cat(sprintf("  [warn] %-12s (soft; not on the 10/11 code path)\n", p))

ok <- length(miss_hard) == 0L

# taxdat: installed from GitHub (HopkinsIDD/taxdat), not CRAN; 10 also calls
# three of its internals (:::), which an older / different build may lack.
if (have("taxdat")) {
  ns   <- asNamespace("taxdat")
  fns  <- c("map_theme", "output_plot_map", "get_AFRO_region", "get_AFRO_region_levels",
            "get_lakes", "get_rivers", "color_rivers", "color_lake_fill", "color_lake_border")
  nofn <- fns[!vapply(fns, exists, logical(1), envir = ns, inherits = FALSE)]
  e <- new.env()
  has_afr <- isTRUE(tryCatch({ utils::data("afr_sf", package = "taxdat", envir = e); exists("afr_sf", envir = e) },
                             error = function(err) FALSE, warning = function(w) FALSE))
  src <- tryCatch(utils::packageDescription("taxdat")$RemoteSha, error = function(e) NULL)
  cat(sprintf("  taxdat %s at %s%s\n", ver("taxdat"), find.package("taxdat"),
              if (length(src)) paste0(" (RemoteSha ", substr(src, 1, 10), ")") else ""))
  if (length(nofn)) { cat("  [MISS] taxdat functions:", paste(nofn, collapse = ", "), "\n"); ok <- FALSE }
  if (!has_afr)     { cat("  [MISS] taxdat dataset afr_sf\n"); ok <- FALSE }
} else {
  cat("  taxdat is not installed in any of the .libPaths() above. Install it from\n",
      "  the HopkinsIDD source into your R 4.3.2 library (with the module line of\n",
      "  this script loaded, so sf/GDAL match), e.g. R CMD INSTALL <taxdat source dir>\n", sep = "")
}

# sf must link the GDAL stack loaded by the module line.
if (have("sf")) {
  v <- tryCatch(sf::sf_extSoftVersion(), error = function(e) NULL)
  if (is.null(v)) { cat("  [MISS] sf loads but sf_extSoftVersion() failed\n"); ok <- FALSE }
  else cat("  sf links:", paste(names(v), v, sep = " ", collapse = ", "), "\n")
}

if (!ok) { cat("\nMissing required packages:", paste(miss_hard, collapse = ", "), "\n"); quit(status = 1) }
cat("All required R packages present.\n")
RS

if [[ "$NEED_QUARTO" == "true" ]]; then
  if command -v "$QUARTO" >/dev/null 2>&1; then
    echo "quarto : $(command -v "$QUARTO") ($("$QUARTO" --version))"
  else
    echo "Error: quarto not found (set QUARTO=/path/to/quarto, or pass --skip-11)." >&2
    exit 1
  fi
fi

if [[ "$CHECK_ONLY" == "true" ]]; then
  echo "--- Check only: nothing run."
  exit 0
fi

# ---------------------------------------------------------------------------
# 1. Phase 10 — manuscript figures
# ---------------------------------------------------------------------------
if [[ "$SKIP_10" != "true" ]]; then
  ARGS=(-C "$CONFIG")
  [[ -n "$FIGURES" ]]       && ARGS+=(--figures "$FIGURES")
  [[ -n "$SUPPLEMENTARY" ]] && ARGS+=(--supplementary "$SUPPLEMENTARY")
  echo "--- Rscript analysis/10_manuscript_figures.R ${ARGS[*]}"
  Rscript analysis/10_manuscript_figures.R "${ARGS[@]}"
fi

# ---------------------------------------------------------------------------
# 2. Phase 11 — key-numbers report
# ---------------------------------------------------------------------------
if [[ "$SKIP_11" != "true" ]]; then
  echo "--- ${QUARTO} render analysis/11_manuscript_numbers.qmd -P config_path:${CONFIG}"
  "$QUARTO" render analysis/11_manuscript_numbers.qmd -P "config_path:${CONFIG}"
fi

echo "***** End of script *****"
date
