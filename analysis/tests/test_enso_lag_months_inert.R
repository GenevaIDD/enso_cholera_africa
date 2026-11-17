# test_enso_lag_months_inert.R — HARD: the legacy linear-ENSO keys are inert on
#   the deployed fits.
#
# Usage:
#   Rscript analysis/tests/test_enso_lag_months_inert.R
#   Rscript analysis/tests/test_enso_lag_months_inert.R --all
#
# Details:
#   `enso.lag_months` and `<track>.random_effects.re_adm1_enso_lag5_slope` are
#   read only under `<track>.covariates.enso` (the linear `enso_lag<N>` term,
#   its RE slope and the run-file `_enso-<src>lag<N>` token:
#   formula.R::.resolve_fixed_covariates() / build_formula_terms(),
#   paths.R::make_occ_run_file() and its siblings), which the deployed configs
#   have off, so a deploy may strip both keys without changing a fit. Quick,
#   no fit, no data/. --all checks every config instead of a sample.

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(brms)
})

source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/formula.R")
source("analysis/utils/priors.R")
source("analysis/utils/paths.R") # make_run_file(); nbr_covariate_name() (call-time dependency of formula.R)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

ALL <- "--all" %in% commandArgs(trailingOnly = TRUE)
N_PER_BATCH <- 12L
RESULTS_PATH <- "analysis/configs/results/manuscript.yml"
TRACK_OF <- c(occ_dir = "occurrence", dur_dir = "duration", mon_dir = "monthly")
KIND_PFX <- c(occurrence = "", duration = "dur_", monthly = "mon_")
SLOPE_KEY <- "re_adm1_enso_lag5_slope"

cat(
  "# enso.lag_months / re_adm1_enso_lag5_slope inert on the deployed fit batches\n"
)

# --- Synthetic frame (test_formula_builder_golden.R's AMPLE_DAT) -------------
.make_dat <- function(n, region, country, adm1, adm2_gid, ttlo, mso_band) {
  d <- data.frame(
    AFRO_region = region,
    country = country,
    ADM1 = adm1,
    ADM2 = adm2_gid,
    gid = adm2_gid,
    time_to_last_outbreak = factor(ttlo),
    mso_band = mso_band,
    stringsAsFactors = FALSE
  )
  ms_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(ms_cols) <- paste0("msoms", 1:5)
  seas_cols <- as.data.frame(matrix(seq_len(n * 5) / n, nrow = n, ncol = 5))
  names(seas_cols) <- paste0("seasms", 1:5)
  dplyr::bind_cols(d, ms_cols, seas_cols)
}
AMPLE_DAT <- .make_dat(
  n = 16,
  region = rep(c("East", "West"), each = 8),
  country = rep(c("KEN", "TZA", "NGA", "GHA"), each = 4),
  adm1 = rep(paste0("A", 1:8), each = 2),
  adm2_gid = paste0("G", 1:16),
  ttlo = rep(c("[0,3)", "[3,6)"), length.out = 16),
  mso_band = rep(c(1, 2, 3), length.out = 16)
)

# --- Variants and their record ---------------------------------------------
VARIANT_DIR <- file.path(tempdir(), "enso_lag_months_inert")
unlink(VARIANT_DIR, recursive = TRUE)

drop_slope_key <- function(y) {
  for (tr in TRACK_OF) {
    re <- y[[tr]][["random_effects"]]
    if (!is.null(re) && SLOPE_KEY %in% names(re)) {
      y[[tr]][["random_effects"]][[SLOPE_KEY]] <- NULL
    }
  }
  y
}
set_lag <- function(y, v) {
  y[["enso"]] <- y[["enso"]] %||% list()
  y[["enso"]][["lag_months"]] <- v
  y
}

VARIANTS <- list(
  as_is = function(y) y,
  lag3 = function(y) set_lag(y, 3L),
  lag5 = function(y) set_lag(y, 5L),
  no_lag = function(y) {
    if (!is.null(y[["enso"]])) {
      y[["enso"]][["lag_months"]] <- NULL
    }
    y
  },
  no_slopekey = drop_slope_key
)

# Same basename in a per-variant dir, in case anything reads the file name.
resolve_variant <- function(path, vname, edit) {
  d <- file.path(VARIANT_DIR, vname)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(d, basename(path))
  yaml::write_yaml(edit(yaml::read_yaml(path)), f)
  opt <- NULL
  invisible(capture.output(
    opt <- suppressWarnings(suppressMessages(
      make_options_from_config(
        option_list = list(),
        cli_argv = c("--config", f)
      )
    ))
  ))
  opt
}

record_of <- function(opt, track) {
  blk <- opt[[track]]
  brms_cfg <- blk$brms_cb %||% blk$brms %||% list()
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  if (track == "occurrence") {
    parts <- q(build_formula_terms(
      opt,
      n_country = 4L,
      n_region = 2L,
      dat = AMPLE_DAT
    ))
    priors <- q(build_occ_priors(opt, parts, brms_cfg))
  } else if (track == "duration") {
    parts <- q(build_duration_formula_parts(opt, AMPLE_DAT))
    priors <- q(build_duration_priors(opt, parts, brms_cfg, dat = AMPLE_DAT))
  } else {
    parts <- q(build_monthly_formula_parts(opt, AMPLE_DAT))
    priors <- q(build_monthly_priors(opt, brms_cfg, parts))
  }
  kinds <- paste0(
    KIND_PFX[[track]],
    c("fit", "re_stats", "coef_summary", "loo")
  )
  files <- vapply(
    kinds,
    function(k) {
      tryCatch(as.character(make_run_file(opt, kind = k)), error = function(e) {
        paste("ERROR:", conditionMessage(e))
      })
    },
    character(1)
  )
  list(
    terms = parts,
    dlnm = blk$dlnm,
    priors = as.data.frame(priors),
    files = files
  )
}

# --- [1] follow the results config to the fit batches -----------------------
# Follows analysis/configs/results/manuscript.yml -> results.deploy_config ->
# deploy.fits.{occ,dur,mon}_dir (as test_deployed_model_terms.R does) and
# takes a deterministic sample of each batch's configs (every k-th sorted
# config, first and last included; --all for every config).
cat("\n[1] results config -> deploy config -> fit batches\n")
deploy_path <- yaml::read_yaml(RESULTS_PATH)$results$deploy_config %||% ""
expect(
  nzchar(deploy_path) && file.exists(deploy_path),
  sprintf("results.deploy_config resolves (%s)", deploy_path)
)
if (!nzchar(deploy_path) || !file.exists(deploy_path)) {
  quit(status = 1)
}
fits <- yaml::read_yaml(deploy_path)$deploy$fits

# --- [2] the variants agree, config by config -------------------------------
# For each config: five raw-YAML variants (as is, lag_months 3, lag_months 5,
# lag_months removed, re_adm1_enso_lag5_slope removed from every track block),
# each resolved with make_options_from_config(), must give identical formula
# terms (the track's builder on a synthetic frame), the resolved dlnm block,
# the brms priors (build_*_priors()) and the run-file names (make_run_file():
# fit, re_stats, coef_summary, loo). (A stripped lag_months is filled with 11
# by the runtime, config.R.)
ctrl_cfg <- NULL
for (key in names(TRACK_OF)) {
  track <- TRACK_OF[[key]]
  bdir <- file.path("analysis/configs", fits[[key]] %||% "")
  cfgs <- sort(setdiff(
    list.files(bdir, pattern = "\\.yml$", full.names = TRUE),
    file.path(bdir, "manifest.yml")
  ))
  cat(sprintf(
    "\n[2] %s: %s (%d configs)\n",
    track,
    fits[[key]] %||% "?",
    length(cfgs)
  ))
  expect(
    length(cfgs) > 0L,
    sprintf("batch %s has configs", fits[[key]] %||% "?")
  )
  if (!length(cfgs)) {
    next
  }
  pick <- if (ALL || length(cfgs) <= N_PER_BATCH) {
    seq_along(cfgs)
  } else {
    unique(round(seq(1, length(cfgs), length.out = N_PER_BATCH)))
  }
  cfgs <- cfgs[pick]
  if (track == "occurrence") {
    ctrl_cfg <- cfgs[[1L]]
  }

  raw <- lapply(cfgs, yaml::read_yaml)
  expect(
    all(vapply(
      raw,
      function(y) !isTRUE(y[[track]][["covariates"]][["enso"]]),
      logical(1)
    )),
    sprintf(
      "%s.covariates.enso is off in all %d sampled configs (the gate every read sits behind)",
      track,
      length(cfgs)
    )
  )
  has_keys <- vapply(
    raw,
    function(y) !is.null(y[["enso"]][["lag_months"]]),
    logical(1)
  )
  cat(sprintf(
    "  %d / %d sampled configs carry enso.lag_months; %d carry %s\n",
    sum(has_keys),
    length(cfgs),
    sum(vapply(
      raw,
      function(y) {
        any(vapply(
          TRACK_OF,
          function(tr) {
            SLOPE_KEY %in% names(y[[tr]][["random_effects"]])
          },
          logical(1)
        ))
      },
      logical(1)
    )),
    SLOPE_KEY
  ))

  bad <- character(0)
  errs <- character(0)
  for (cfg in cfgs) {
    recs <- tryCatch(
      lapply(names(VARIANTS), function(v) {
        record_of(resolve_variant(cfg, v, VARIANTS[[v]]), track)
      }),
      error = function(e) {
        errs <<- c(errs, sprintf("%s: %s", basename(cfg), conditionMessage(e)))
        NULL
      }
    )
    if (is.null(recs)) {
      next
    }
    names(recs) <- names(VARIANTS)
    for (v in setdiff(names(VARIANTS), "as_is")) {
      for (part in names(recs$as_is)) {
        if (!identical(recs$as_is[[part]], recs[[v]][[part]])) {
          bad <- c(bad, sprintf("%s [%s: %s]", basename(cfg), v, part))
        }
      }
    }
  }
  expect(
    !length(errs),
    sprintf(
      "every variant resolves and builds (errors: %s)",
      paste(head(errs, 3), collapse = " / ")
    )
  )
  expect(
    !length(bad),
    sprintf(
      "%d configs x %d variants: identical terms, dlnm block, priors and run-file names%s",
      length(cfgs),
      length(VARIANTS) - 1L,
      if (length(bad)) {
        paste0(" (differ: ", paste(head(bad, 5), collapse = ", "), ")")
      } else {
        ""
      }
    )
  )
}

# --- [3] positive control ---------------------------------------------------
# With the track's covariates.enso switched on, lag 3 vs 5 does change the
# occurrence run-file name, so the identity above is the gate and not a blind
# spot of the comparison.
cat(
  "\n[3] positive control: covariates.enso on -> lag_months reaches the run-file name\n"
)
if (!is.null(ctrl_cfg)) {
  enso_on <- function(y) {
    y[["occurrence"]][["covariates"]][["enso"]] <- TRUE
    y
  }
  f3 <- tryCatch(
    as.character(make_run_file(
      resolve_variant(ctrl_cfg, "ctrl3", function(y) set_lag(enso_on(y), 3L)),
      kind = "fit"
    )),
    error = function(e) NA_character_
  )
  f5 <- tryCatch(
    as.character(make_run_file(
      resolve_variant(ctrl_cfg, "ctrl5", function(y) set_lag(enso_on(y), 5L)),
      kind = "fit"
    )),
    error = function(e) NA_character_
  )
  expect(
    !is.na(f3) &&
      !is.na(f5) &&
      !identical(f3, f5) &&
      grepl("lag3", f3, fixed = TRUE) &&
      grepl("lag5", f5, fixed = TRUE),
    sprintf(
      "with covariates.enso on, lag 3 vs 5 give different fit names (%s vs %s)",
      basename(f3 %||% "NA"),
      basename(f5 %||% "NA")
    )
  )
}
unlink(VARIANT_DIR, recursive = TRUE)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
