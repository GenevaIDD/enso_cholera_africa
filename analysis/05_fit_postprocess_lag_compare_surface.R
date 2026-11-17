# 05_fit_postprocess_lag_compare_surface.R — ENSO exposure-response surface
#   at the best single lag, overall + by country, per ENSO group
#
# Usage:
#   Rscript analysis/05_fit_postprocess_lag_compare_surface.R \
#     -d analysis/configs/<occ_singlelag_batch>/
#   Rscript analysis/05_fit_postprocess_lag_compare_surface.R \
#     -d analysis/configs/<batch>/ --se_mult 1 --group enso_west
#
# Details:
#   Reads a batch dir of single-lag ENSO configs and each group's cached
#   loo_compare_lags_<group>_<date>.rds ranking (05_fit_postprocess_loo_lags.R),
#   picks the best lag, and extracts that lag's brms_cb ENSO surface.
#   Writes generated_data/surface_best_lag_<batch>.rds, a list keyed by ENSO
#   group; read by 06_fit_figures_lag_compare.R.

# 1. Preamble ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/loo_lags.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-d", "--batch_dir"),
    default = NULL,
    type = "character",
    help = "Directory of single-lag YAML configs (required)"
  ),
  make_option(
    c("--se_mult"),
    default = 2,
    type = "double",
    help = paste(
      "SE multiplier for the best-lag rule (max single_lag within",
      "--se_mult SEs of elpd_diff from the rank-1 fit); passed to",
      "select_best_lag() [default: 2]"
    )
  ),
  make_option(
    c("--lag_rule"),
    default = "longest_within_se",
    type = "character",
    help = "Best-lag rule: longest_within_se | best_elpd (utils/loo_lags.R; deploy.lag_selection.rule) [default: %default]"
  ),
  make_option(
    c("--group"),
    default = NULL,
    type = "character",
    help = "Filter to this grouping group (e.g. enso_east) [default: all groups found]"
  ),
  make_option(
    c("--var"),
    default = "enso",
    type = "character",
    help = "Crossbasis variable name to extract [default: enso]"
  ),
  make_option(
    c("--n_at"),
    default = 50L,
    type = "integer",
    help = "Number of exposure grid points for the full surface [default: 50]"
  )
)
script_opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(script_opt$batch_dir) || !dir.exists(script_opt$batch_dir)) {
  stop("Provide a valid batch directory via -d <dir/>.")
}

se_mult <- as.numeric(script_opt$se_mult %||% 2)
var_nm <- script_opt$var %||% "enso"
n_at <- as.integer(script_opt$n_at %||% 50L)

config_paths <- list.files(
  script_opt$batch_dir,
  pattern = "\\.yml$",
  full.names = TRUE
)
config_paths <- config_paths[basename(config_paths) != "manifest.yml"]
if (length(config_paths) == 0) {
  stop("No .yml configs found in: ", script_opt$batch_dir)
}

batch <- basename(sub("/+$", "", script_opt$batch_dir))

# 3. Helpers -------------------------------------------------------------------

# Fill an opt list read from YAML with config_defaults (mirrors sibling scripts).
.fill_defaults <- function(opt, model) {
  defaults <- get_default_options()
  for (nm in setdiff(
    names(defaults),
    c("occurrence", "size", "duration", "monthly", "climate")
  )) {
    if (is.null(opt[[nm]])) opt[[nm]] <- defaults[[nm]]
  }
  model_defaults <- defaults[[model]] %||% list()
  for (nm in names(model_defaults)) {
    if (is.null(opt[[model]][[nm]])) opt[[model]][[nm]] <- model_defaults[[nm]]
  }
  opt
}

# Build a tibble from an extract_cb_surface()-shaped list (single-lag: 1 lag col).
.surface_to_tibble <- function(surf) {
  logRR <- surf$matfit[, 1]
  tibble::tibble(
    exposure = surf$at,
    logRR = logRR,
    logRR_lo = surf$matlo[, 1],
    logRR_hi = surf$mathi[, 1],
    RR = exp(logRR),
    RR_lo = exp(surf$matlo[, 1]),
    RR_hi = exp(surf$mathi[, 1])
  )
}

# Build a tibble (one row per country x exposure point) from
# extract_cb_country_surface()'s named-list-of-lists output.
.country_surface_to_tibble <- function(ctry_surf) {
  dplyr::bind_rows(lapply(names(ctry_surf), function(ctry) {
    s <- ctry_surf[[ctry]]
    tibble::tibble(
      country = ctry,
      exposure = s$at,
      logRR = s$matfit,
      logRR_lo = s$matlo,
      logRR_hi = s$mathi,
      RR = exp(s$matfit),
      RR_lo = exp(s$matlo),
      RR_hi = exp(s$mathi)
    )
  }))
}

.add_nina_nino_condition <- function(df) {
  dplyr::mutate(
    df,
    condition = ifelse(exposure < 0, "nina", "nino"),
    .after = exposure
  )
}

# 4. Locate the ranking file(s) ------------------------------------------------
#
# One loo_compare_lags_<group>_<date>.rds per ENSO group; if several dated
# files exist for the same group, the latest (filename date, then mtime) is
# used. utils/loo_lags.R::loo_compare_ranking_files() resolves the location
# (<outdir>/loo_compare/<batch>/, legacy in-config-dir fallback) and matches
# the --group token verbatim (ranking filenames keep the underscores, unlike
# the YAML config filenames' _grp-ensoeast form).

rank_tbl <- loo_compare_ranking_files(
  script_opt$batch_dir,
  group = script_opt$group
) %>%
  dplyr::filter(!is.na(group)) # grouped rankings only
if (nrow(rank_tbl) == 0) {
  stop(
    "No loo_compare_lags ranking found for ",
    script_opt$batch_dir,
    " — run 05_fit_postprocess_loo_lags.R first."
  )
}

message(sprintf(
  "Found %d ENSO group ranking(s) in %s: %s",
  nrow(rank_tbl),
  script_opt$batch_dir,
  paste(rank_tbl$group, collapse = ", ")
))

# 5. Per-group extraction ------------------------------------------------------
# Per ENSO group: pick the lag by the shared rule with the occurrence floor
# off (select_best_lag(..., enforce_occ_floor = FALSE)), load that lag's
# brms_cb fit + cb_meta, and extract the ERF overall (extract_cb_surface())
# and by member country (extract_cb_country_surface(), shared with 06). Two
# grids: full_surface (n_at points over the range) and nina_nino (ENSO = -2,
# -1, +1, +2, computed directly, not interpolated).

out <- list()

for (i in seq_len(nrow(rank_tbl))) {
  group_token <- rank_tbl$group[i]
  rank_file <- rank_tbl$file[i]
  message("\n=== Group: ", group_token, " (", rank_file, ") ===")

  rk <- readRDS(rank_file)
  ranking <- rk$ranking
  track <- rk$model

  ls_rule <- lag_selection_cfg(list(
    rule = script_opt$lag_rule,
    se_mult = se_mult
  ))
  sel <- select_best_lag(
    ranking,
    track = track,
    se_mult = se_mult,
    enforce_occ_floor = FALSE,
    min_lag = ls_rule$min_lag,
    rule = ls_rule$rule
  )
  best_lag <- sel$lag
  message(sprintf(
    "Track: %-10s  best lag (no occurrence floor, se_mult=%.1f): %d (rank %d)",
    track,
    se_mult,
    best_lag,
    sel$rank
  ))

  # Locate the matching config (exact "_sl-<lag>_grp-<grp_fname>" token; prefer
  # the canonical "<batch>_sl-<lag>_grp-<grp_fname>.yml" filename if more than
  # one candidate matches, e.g. a "_prevcases" sibling variant).
  grp_fname <- gsub("[^A-Za-z0-9]", "", group_token)
  sl_grp_token <- sprintf("_sl-%s_grp-%s", sl_token(best_lag), grp_fname)
  candidates <- config_paths[grepl(
    sl_grp_token,
    basename(config_paths),
    fixed = TRUE
  )]
  if (length(candidates) == 0) {
    warning(
      "No config matches '",
      sl_grp_token,
      "' for group '",
      group_token,
      "' — skipping."
    )
    next
  }
  if (length(candidates) > 1) {
    canon <- file.path(
      script_opt$batch_dir,
      paste0(batch, sl_grp_token, ".yml")
    )
    if (canon %in% candidates) {
      candidates <- canon
    } else {
      stop(
        "Ambiguous config match for group '",
        group_token,
        "', lag ",
        best_lag,
        " — candidates: ",
        paste(basename(candidates), collapse = ", ")
      )
    }
  }
  config_path <- candidates[1]

  opt <- yaml::read_yaml(config_path)
  opt <- .fill_defaults(opt, track)

  kpfx <- if (track == "size") {
    "size_"
  } else if (track == "duration") {
    "dur_"
  } else if (track == "monthly") {
    "mon_"
  } else {
    ""
  }
  fit_file <- make_run_file(opt, kind = paste0(kpfx, "fit"))
  meta_file <- make_run_file(opt, kind = paste0(kpfx, "cb_meta"))
  if (!file.exists(fit_file) || !file.exists(meta_file)) {
    warning(
      "Fit/meta missing for group '",
      group_token,
      "', lag ",
      best_lag,
      " — skipping.\n  fit:  ",
      fit_file,
      "\n  meta: ",
      meta_file
    )
    next
  }

  message("Loading fit: ", fit_file)
  brms_fit <- readRDS(fit_file)
  cb_meta <- readRDS(meta_file)
  vm <- cb_meta[[var_nm]]
  if (is.null(vm)) {
    warning(
      "Variable '",
      var_nm,
      "' not found in cb_meta for group '",
      group_token,
      "' — skipping."
    )
    rm(brms_fit)
    next
  }
  if (is.null(vm$single_lag)) {
    warning(
      "Variable '",
      var_nm,
      "' is not a single-lag crossbasis for group '",
      group_token,
      "' — this script only supports single-lag sweeps. Skipping."
    )
    rm(brms_fit)
    next
  }

  # Group member countries (defensive restriction; grouped fits are already
  # trained on only these countries, but this pins the country set to the
  # registry rather than trusting whatever ranef() happens to carry).
  reg <- load_grouping(
    opt$grouping$source %||% "enso_grouping",
    "analysis/groupings"
  )
  group_isos <- reg$iso3[reg$group == (opt$grouping$group %||% group_token)]
  if (length(group_isos) == 0) {
    warning(
      "No countries found in registry '",
      opt$grouping$source,
      "' for group '",
      opt$grouping$group,
      "'."
    )
  }

  # Full surface: n_at points over the variable's trained range
  at_full <- seq(vm$range[1], vm$range[2], length.out = n_at)
  overall_full <- .surface_to_tibble(
    extract_cb_surface(
      fit = brms_fit,
      meta = cb_meta,
      vars = var_nm,
      at = at_full
    )[[var_nm]]
  )
  ctry_full <- .country_surface_to_tibble(
    extract_cb_country_surface(
      fit = brms_fit,
      meta = vm,
      at = at_full,
      countries = group_isos
    )
  )

  # Nina (-2, -1) / Nino (+1, +2): recomputed exactly at these 4 points
  at_nn <- c(-2, -1, 1, 2)
  overall_nn <- .add_nina_nino_condition(.surface_to_tibble(
    extract_cb_surface(
      fit = brms_fit,
      meta = cb_meta,
      vars = var_nm,
      at = at_nn
    )[[var_nm]]
  ))
  ctry_nn <- .add_nina_nino_condition(.country_surface_to_tibble(
    extract_cb_country_surface(
      fit = brms_fit,
      meta = vm,
      at = at_nn,
      countries = group_isos
    )
  ))

  out[[group_token]] <- list(
    track = track,
    group = group_token,
    best_lag = best_lag,
    best_lag_rank = sel$rank,
    se_mult = se_mult,
    n_countries = length(unique(ctry_full$country)),
    full_surface = list(overall = overall_full, by_country = ctry_full),
    nina_nino = list(overall = overall_nn, by_country = ctry_nn)
  )

  rm(brms_fit)
  message(sprintf(
    "Group '%s': done (lag %d, %d countries).",
    group_token,
    best_lag,
    length(unique(ctry_full$country))
  ))
}

if (length(out) == 0) {
  stop("No group produced a surface — see warnings above.")
}

# 6. Write outputs -------------------------------------------------------------
# Output: a named list keyed by ENSO group (e.g. "enso_east"), each entry:
#   list(track=, group=, best_lag=, best_lag_rank=, se_mult=, n_countries=,
#        full_surface = list(overall = <tibble>, by_country = <tibble>),
#        nina_nino    = list(overall = <tibble>, by_country = <tibble>))
# Tibble columns: exposure [+ condition for nina_nino], logRR/logRR_lo/logRR_hi
# (posterior median + 95% CI, log-RR scale), RR/RR_lo/RR_hi (exponentiated);
# by_country tibbles additionally carry a `country` (ISO3) column.
out_file <- here::here(
  "generated_data",
  stringr::str_glue("surface_best_lag_{batch}.rds")
)
saveRDS(out, out_file)
message("\nSaved -> ", out_file)
