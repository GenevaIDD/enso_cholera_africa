# dlnm.R — DLNM cross-basis construction, replay and surfaces for the brms_cb
#   fits and the simulator kernels.
#
# Provides:
#   make_lag_matrix()               season-aware lag matrix of one variable
#   build_cb_block()                fit-time cross-basis columns + cb_meta
#   build_cb_columns_from_Q(), build_cb_design()   replay the stored basis
#   build_mso_basis()               months-since-onset lookup table (mon)
#   extract_cb_surface()            exposure- / lag-response surfaces
#   extract_cb_country_surface()    country-level single-lag surfaces
#
# Details:
#   Replay reuses the fit-time argvar / arglag (resolved knots) stored in
#   cb_meta, so a simulated or postprocessed basis is bit-identical to the fit.

#' .LATENT_CB_VARS
#'
#' Cross-basis variables of the monthly-size (`mon`) track that are NOT panel
#' columns: they are reconstructed from filtered/simulated latent state inside
#' the kernel month loop, so they differ per particle and cannot be precomputed
#' per (unit, month) the way every occ/siz/dur cross-basis variable can.
#'
#'   mso                       -- months since onset; maintained by the duration
#'                                promotion block (spatpomp_core.R state slot 11)
#'   cases_prev, log_cases_prev -- previous month's case count; a function of the
#'                                mon draw made one month earlier
#'
#' Consumers must skip these wherever a panel lag column is otherwise required:
#'   * `.lags_required_from_draws()` (utils/simulations.R) -- must not demand a
#'     panel column that by construction cannot exist;
#'   * the smoothing cross-basis rebuild (utils/smoothing_core.R) -- must not
#'     warn-and-drop them as "unsupported non-ENSO variables";
#'   * `build_cb_mats()` (utils/simulator_core.R) -- skips them; the kernel calls
#'     build_cb_columns_from_Q() directly on a latent-value lookup table instead.
.LATENT_CB_VARS <- c("mso", "cases_prev", "log_cases_prev")

#' .onebasis_replay_args
#'
#' Sanitise a stored `onebasis_attrs` list into the argument list for a replay
#' `do.call(dlnm::onebasis, c(list(x), args))`. Two things are stripped:
#'
#'   * `range` -- a recorded attribute, not a onebasis() argument.
#'   * `degree` WHEN `fun == "ns"` -- dlnm::onebasis() records a spurious
#'     `degree` attribute even for a natural-spline (`ns`) basis, but splines::ns()
#'     has no `degree` parameter and errors ("unused argument (degree = ...)")
#'     when it is forwarded on replay. `bs()` accepts `degree`, which is why this
#'     only surfaced once the ENSO argvar switched bs -> ns. Stripping
#'     it is a no-op for the basis values (ns ignores degree) and is required for
#'     already-fitted ns cb_meta sidecars, whose stored attrs are frozen.
.onebasis_replay_args <- function(ob_attrs) {
  args <- (ob_attrs %||% list())[setdiff(names(ob_attrs %||% list()), "range")]
  if (identical(args$fun, "ns")) {
    args$degree <- NULL
  }
  args
}

#' make_lag_matrix
#'
#' Build an n × (max_lag + 1) matrix of lagged exposures for a single variable,
#' suitable for passing as pre-lagged input to dlnm::crossbasis().
#' Columns correspond to lags 0, 1, …, max_lag (lag0 = concurrent value).
#'
#' Two dispatching branches:
#'
#'   var == "enso":
#'     Reuses the pre-built enso_lag0 … enso_lag{max_lag} columns written by
#'     join_enso_lags() (calendar-aligned global index; no season-year reset).
#'     max_lag must not exceed 12.  Any NA values (coverage gaps in the ONI
#'     series) are imputed to 0 (neutral ENSO) with a warning.
#'
#'   all other vars (climate covariates):
#'     Generalises make_climate_lags() to lags 0..max_lag with season-aware
#'     grouping.  Lags are computed within each (uid × season-year) group so
#'     they cannot cross cholera-season boundaries.  Season-start boundary
#'     cells (dplyr::lag() → NA) are set to 0, matching the ttlo_mat /
#'     season-boundary convention already used in the pipeline.
#'     Data must be sorted by uid, t_start (guaranteed by 03_prepare_occ.R).
#'
#' @param df      survival dataset (output of load_surv_data()): must contain
#'                the column `var` (climate) or enso_lag0..max_lag (ENSO).
#' @param var     variable stem: "mean_precip", "dev_temp_daymean", "enso", …
#' @param max_lag maximum lag in months (inclusive).  ENSO max: 12.
#' @return numeric matrix, n × (max_lag + 1), colnames "lag0" … "lag{max_lag}".
make_lag_matrix <- function(df, var, max_lag) {
  max_lag <- as.integer(max_lag)

  if (identical(var, "enso")) {
    # --- ENSO branch: reuse pre-built enso_lag* columns ---
    if (max_lag > 12L) {
      stop(
        "make_lag_matrix: ENSO max_lag cannot exceed 12 ",
        "(join_enso_lags() creates columns up to enso_lag12)"
      )
    }
    lag_cols <- sprintf("enso_lag%d", seq(0L, max_lag))
    missing <- setdiff(lag_cols, colnames(df))
    if (length(missing) > 0) {
      stop(
        "make_lag_matrix: required ENSO columns not found in df: ",
        paste(missing, collapse = ", "),
        ". Ensure join_enso_lags() was called with max_lag >= ",
        max_lag
      )
    }
    Q <- as.matrix(df[, lag_cols, drop = FALSE])
    na_count <- sum(is.na(Q))
    if (na_count > 0) {
      warning(
        "make_lag_matrix: ",
        na_count,
        " NA(s) in ENSO lag columns; ",
        "imputing to 0 (neutral ENSO conditions)."
      )
      Q[is.na(Q)] <- 0
    }
    colnames(Q) <- sprintf("lag%d", seq(0L, max_lag))
    return(Q)
  }

  # --- Climate branch: read pre-computed lag columns ----------------------
  # Lag columns (<var>_lag1 ... <var>_lag{max_lag}) must have been written by
  # make_climate_lags() during data preparation (03_prepare_occ.R /
  # 03b_prepare_siz.R).  If any are missing, a clear error is raised so the
  # user knows to re-run prep with the required max_lag.
  if (!var %in% colnames(df)) {
    stop(
      "make_lag_matrix: column '",
      var,
      "' not found in df. ",
      "Available columns: ",
      paste(head(colnames(df), 20), collapse = ", ")
    )
  }

  needed <- sprintf("%s_lag%d", var, seq(1L, max_lag))
  missing_cols <- setdiff(needed, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "make_lag_matrix: pre-computed lag columns not found: ",
      paste(missing_cols, collapse = ", "),
      ". Re-run 03_prepare_occ.R / 03b_prepare_siz.R with max_lag >= ",
      max_lag
    )
  }

  lag_col_names <- c(var, needed) # lag0 = <var>; lag1..max_lag = <var>_lagk
  Q <- as.matrix(df[, lag_col_names, drop = FALSE])
  Q[is.na(Q)] <- 0 # boundary NAs → 0 (ttlo_mat / season-boundary convention)
  colnames(Q) <- sprintf("lag%d", seq(0L, max_lag))
  Q
}


#' extract_cb_surface
#'
#' Exposure-response / lag-response surfaces from a brms fit with
#' dlnm::crossbasis() predictors: dlnm::crosspred() per posterior draw,
#' summarised as posterior median + credible intervals per variable.
#'
#' @param fit    brms fit object from a `04_fit_*_brms_cb.R` script.
#' @param meta   Named list of crossbasis metadata (from _cb_meta.rds sidecar).
#'               Each entry must contain: argvar, arglag, lag, df, range,
#'               cb_col_names (column names in the brms fit).
#' @param vars   Character vector of variable names to process.
#'               NULL (default) processes all names(meta).
#' @param at     Numeric vector of exposure values for the prediction grid.
#'               NULL (default): 50 equally-spaced values over the stored range.
#' @param cumul  Logical: also compute cumulative exposure-response?  Default FALSE.
#' @param probs  Quantile levels for the outer (95%) credible interval.
#'               Default c(0.025, 0.5, 0.975).
#' @param probs50 Quantile levels for the inner (50%) credible interval.
#'               Default c(0.25, 0.75).
#' @return Named list (one entry per variable) containing:
#'   $at     — exposure prediction grid
#'   $lag    — lag axis (from crosspred)
#'   $matfit — median surface (length(at) × n_lags matrix, log-RR scale)
#'   $matlo  / $mathi    — lower/upper 95% CI surfaces
#'   $matlo50 / $mathi50 — lower/upper 50% CrI surfaces
#'   $allfit — full posterior surface array (at × lags × draws)
#'   $cumfit / $cumlo / $cumhi       — cumulative curves + 95% CI (if cumul = TRUE)
#'   $cumlo50 / $cumhi50             — cumulative 50% CrI (if cumul = TRUE)
#'   $meta   — the original metadata entry
extract_cb_surface <- function(
  fit,
  meta,
  vars = NULL,
  at = NULL,
  cumul = FALSE,
  probs = c(0.025, 0.5, 0.975),
  probs50 = c(0.25, 0.75)
) {
  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stop("extract_cb_surface: package 'dlnm' is required.")
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("extract_cb_surface: package 'brms' is required.")
  }

  if (is.null(vars)) {
    vars <- names(meta)
  }
  vars <- intersect(vars, names(meta))
  if (length(vars) == 0) {
    warning("extract_cb_surface: no matching variables found in meta.")
    return(list())
  }

  results <- list()

  for (v in vars) {
    vm <- meta[[v]]
    cb_col_names <- vm$cb_col_names
    if (is.null(cb_col_names) || length(cb_col_names) == 0) {
      warning(
        "extract_cb_surface: no cb_col_names in meta for '",
        v,
        "'. Skipping."
      )
      next
    }

    # brms coefficient names for the CB columns (prefix "b_")
    b_names <- sprintf("b_%s", cb_col_names)

    draws <- tryCatch(
      brms::as_draws_matrix(fit, variable = b_names),
      error = function(e) {
        warning(
          "extract_cb_surface: could not extract draws for '",
          v,
          "': ",
          conditionMessage(e)
        )
        NULL
      }
    )
    if (is.null(draws) || nrow(draws) == 0) {
      next
    }
    draws <- unclass(draws) # plain numeric matrix (S × K)

    n_cols <- ncol(draws)

    # Reconstruct a basis stub for dlnm::crosspred()
    # (the stub carries only the stored attributes; data values irrelevant).
    # Single-lag fits use dlnm::onebasis (no lag axis); the distributed-lag
    # fits use dlnm::crossbasis.  crosspred() and the downstream matfit
    # stacking handle the single-column (one-lag) case transparently.
    if (!is.null(vm$single_lag)) {
      cb_stub <- matrix(0, nrow = 1L, ncol = n_cols)
      ob_attrs <- vm$onebasis_attrs %||% list()
      for (nm in names(ob_attrs)) {
        attr(cb_stub, nm) <- ob_attrs[[nm]]
      }
      class(cb_stub) <- c("onebasis", "matrix")
    } else {
      cb_stub <- structure(
        matrix(0, nrow = 1L, ncol = n_cols),
        class = c("crossbasis", "matrix"),
        argvar = vm$argvar,
        arglag = vm$arglag,
        lag = vm$lag,
        df = vm$df,
        range = vm$range
      )
    }

    # Prediction grid
    rng <- vm$range %||% c(NA_real_, NA_real_)
    if (is.null(at)) {
      if (anyNA(rng)) {
        warning(
          "extract_cb_surface: unknown range for '",
          v,
          "' — provide `at` explicitly. Skipping."
        )
        next
      }
      at_vals <- seq(rng[1], rng[2], length.out = 50L)
    } else {
      at_vals <- at
    }

    # Per-draw crosspred loop
    n_draws <- nrow(draws)
    draw_preds <- vector("list", n_draws)
    ok_count <- 0L
    for (i in seq_len(n_draws)) {
      b_i <- as.numeric(draws[i, ])
      p_i <- tryCatch(
        dlnm::crosspred(
          cb_stub,
          coef = b_i,
          vcov = diag(0, length(b_i)),
          at = at_vals,
          cen = 0,
          cumul = cumul
        ),
        error = function(e) NULL
      )
      if (!is.null(p_i)) {
        draw_preds[[i]] <- p_i
        ok_count <- ok_count + 1L
      }
    }
    draw_preds <- Filter(Negate(is.null), draw_preds)

    if (length(draw_preds) == 0) {
      warning(
        "extract_cb_surface: all crosspred draws failed for '",
        v,
        "'. Skipping."
      )
      next
    }
    message(sprintf(
      "extract_cb_surface: %d/%d successful draws for '%s'",
      ok_count,
      n_draws,
      v
    ))

    # Stack matfit into a 3D array: length(at_vals) × n_lags × n_draws
    # (dlnm 2.4+ stores log-RR in $matfit; $matRRfit does not exist)
    mat_list <- lapply(draw_preds, function(x) x$matfit)
    mat_arr <- simplify2array(mat_list)

    # Lag sequence derived from matfit column names ("lag0", "lag1", ...)
    lag_seq <- as.numeric(sub("^lag", "", colnames(draw_preds[[1]]$matfit)))

    out <- list(
      at = at_vals,
      lag = lag_seq,
      single_lag = vm$single_lag, # NULL for distributed-lag; integer for onebasis
      matfit = apply(mat_arr, c(1, 2), median),
      matlo = apply(mat_arr, c(1, 2), quantile, probs = probs[1]),
      mathi = apply(mat_arr, c(1, 2), quantile, probs = probs[3]),
      matlo50 = apply(mat_arr, c(1, 2), quantile, probs = probs50[1]),
      mathi50 = apply(mat_arr, c(1, 2), quantile, probs = probs50[2]),
      allfit = mat_arr,
      meta = vm
    )

    if (cumul && !is.null(draw_preds[[1]]$allfit)) {
      cum_list <- lapply(draw_preds, function(x) x$allfit)
      cum_arr <- do.call(cbind, lapply(cum_list, as.numeric))
      out$cumfit <- apply(cum_arr, 1, median)
      out$cumlo <- apply(cum_arr, 1, quantile, probs = probs[1])
      out$cumhi <- apply(cum_arr, 1, quantile, probs = probs[3])
      out$cumlo50 <- apply(cum_arr, 1, quantile, probs = probs50[1])
      out$cumhi50 <- apply(cum_arr, 1, quantile, probs = probs50[2])
    }

    results[[v]] <- out
  }

  results
}


#' extract_cb_country_surface
#'
#' Country-level exposure-response curves for a single-lag (onebasis)
#' crossbasis variable: population CB draws plus the "country" random-effect
#' draws, evaluated at `at` through the stored onebasis attributes (centred at
#' exposure 0, as in extract_cb_surface()). Draw-based, because ranef draws
#' have no crosspred entry point.
#'
#' @param fit brms fit object (single-lag CB fit; ranef() must contain a
#'   "country" group).
#' @param meta one variable's cb_meta entry (must have single_lag,
#'   onebasis_attrs, cb_col_names, range).
#' @param at numeric vector of exposure grid points. NULL (default): 50
#'   equally-spaced values over meta$range.
#' @param countries optional character vector restricting output to these
#'   ISO3 codes (matched against the country token parsed from the ranef
#'   dimname, `"<Region>_<ISO3>"` -> `ISO3`). NULL (default) = every country
#'   with a ranef estimate in the fit.
#' @param n_draws number of posterior draws to sample for the RE combination
#'   [default 200].
#' @param seed RNG seed for draw subsampling [default 1234L].
#' @param probs two quantile levels (lo, hi) for the outer (95%) credible
#'   interval [default c(0.025, 0.975)].
#' @param probs50 two quantile levels (lo, hi) for the inner (50%) credible
#'   interval [default c(0.25, 0.75)].
#' @return named list, one entry per country ISO3 present (after `countries`
#'   filtering), each entry: list(at=, matfit=, matlo=, mathi=, matlo50=,
#'   mathi50=) — numeric vectors of length(at), log-RR scale, posterior median
#'   + 95% and 50% credible intervals.
extract_cb_country_surface <- function(
  fit,
  meta,
  at = NULL,
  countries = NULL,
  n_draws = 200L,
  seed = 1234L,
  probs = c(0.025, 0.975),
  probs50 = c(0.25, 0.75)
) {
  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stop("extract_cb_country_surface: package 'dlnm' is required.")
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("extract_cb_country_surface: package 'brms' is required.")
  }
  if (is.null(meta$single_lag)) {
    stop(
      "extract_cb_country_surface: `meta` has no single_lag — this helper ",
      "only supports single-lag (onebasis) crossbasis variables."
    )
  }

  cb_names <- meta$cb_col_names
  K <- length(cb_names)
  b_pop_raw <- unclass(brms::as_draws_matrix(
    fit,
    variable = paste0("b_", cb_names)
  ))
  colnames(b_pop_raw) <- cb_names

  set.seed(seed)
  draw_idx <- sample.int(nrow(b_pop_raw), min(n_draws, nrow(b_pop_raw)))
  b_pop_draws <- b_pop_raw[draw_idx, , drop = FALSE] # [S x K]

  re_all <- brms::ranef(fit, summary = FALSE)
  if (!"country" %in% names(re_all)) {
    stop("extract_cb_country_surface: fit has no 'country' ranef group.")
  }
  re_country_draws <- re_all[["country"]][draw_idx, , cb_names, drop = FALSE]
  country_lvls <- dimnames(re_country_draws)[[2]] # e.g. "Central Africa_AGO"
  iso3 <- sub(".*_", "", country_lvls)

  keep <- if (is.null(countries)) {
    rep(TRUE, length(iso3))
  } else {
    iso3 %in% countries
  }
  if (!any(keep)) {
    stop(
      "extract_cb_country_surface: none of the fit's country ranef levels ",
      "(",
      paste(unique(iso3), collapse = ", "),
      ") match `countries`."
    )
  }

  rng <- meta$range %||% c(NA_real_, NA_real_)
  at_vals <- if (is.null(at)) {
    if (anyNA(rng)) {
      stop(
        "extract_cb_country_surface: unknown range in `meta` — provide `at` explicitly."
      )
    }
    seq(rng[1], rng[2], length.out = 50L)
  } else {
    at
  }

  ob_args <- .onebasis_replay_args(meta$onebasis_attrs)
  Z <- do.call(dlnm::onebasis, c(list(at_vals), ob_args))
  # Center on the basis evaluated AT exposure = 0 (matches extract_cb_surface's
  # crosspred(cen = 0)) -- NOT on whichever `at_vals` entry happens to be
  # nearest to 0. The two coincide when `at` is a dense grid spanning 0 (the
  # original use case: 50 points over the full range), but diverge for a
  # sparse custom grid with no point near 0 (e.g. Nina/Nino at = c(-2,-1,1,2)),
  # where the nearest-point approach would silently center on -1 instead.
  Z0 <- do.call(dlnm::onebasis, c(list(0), ob_args))
  Z_cen <- Z - matrix(Z0, nrow(Z), K, byrow = TRUE) # [n_at x K]

  results <- list()
  for (i in which(keep)) {
    r_ctry <- matrix(
      re_country_draws[, i, ],
      nrow = nrow(b_pop_draws),
      dimnames = list(NULL, cb_names)
    )
    log_rr_mat <- Z_cen %*% t(b_pop_draws + r_ctry) # [n_at x S]
    results[[iso3[i]]] <- list(
      at = at_vals,
      matfit = apply(log_rr_mat, 1, median),
      matlo = apply(log_rr_mat, 1, quantile, probs = probs[1]),
      mathi = apply(log_rr_mat, 1, quantile, probs = probs[2]),
      matlo50 = apply(log_rr_mat, 1, quantile, probs = probs50[1]),
      mathi50 = apply(log_rr_mat, 1, quantile, probs = probs50[2])
    )
  }
  results
}


#' .cb_lag_exposure
#'
#' Exposure vector of variable `v` at lag `L` (negative = ENSO lead), read from
#' the prep-written lag column (`enso_lag_col_name(L)` for ENSO, `<v>` at lag 0,
#' `<v>_lag<L>` otherwise), with boundary/coverage NAs set to 0 (the
#' make_lag_matrix() convention). Shared by the single-lag and multi-lag (`lags`)
#' modes of build_cb_block() and by build_cb_design(), so fit time and replay
#' read the same column.
#' @param who prefix for the error message (caller context)
#' @noRd
.cb_lag_exposure_col <- function(v, L) {
  if (identical(v, "enso")) {
    enso_lag_col_name(L)
  } else if (L == 0L) {
    v
  } else {
    sprintf("%s_lag%d", v, L)
  }
}
.cb_lag_exposure <- function(dat, v, L, who = "single_lag") {
  L <- as.integer(L)
  xcol <- .cb_lag_exposure_col(v, L)
  if (!xcol %in% colnames(dat)) {
    stop(
      who,
      ": exposure column '",
      xcol,
      "' not found for variable '",
      v,
      "'. Ensure prep wrote the lag column (",
      if (L < 0L) {
        paste0("a lead: prep with --enso_max_lead >= ", -L)
      } else {
        paste0("max_lag >= ", L)
      },
      ")."
    )
  }
  x <- as.numeric(dat[[xcol]])
  x[is.na(x)] <- 0 # boundary/coverage NAs -> 0 (make_lag_matrix convention)
  x
}

#' .cb_lags_spec
#'
#' Resolve a `lags:` spec to an integer lag vector: a length-2 value is an
#' inclusive range (`[-11, 11]` -> -11:11), anything longer is taken as the
#' explicit lag set (sorted, unique). NULL -> NULL (mode not requested).
#' @noRd
.cb_lags_spec <- function(lags) {
  if (is.null(lags)) {
    return(NULL)
  }
  lags <- as.integer(unlist(lags))
  if (length(lags) == 2L) {
    lags <- seq.int(lags[1L], lags[2L])
  }
  lags <- sort(unique(lags))
  if (length(lags) < 2L) {
    stop("dlnm `lags:` needs at least two lags (use `single_lag:` for one).")
  }
  lags
}

#' .cb_lag_col_names
#'
#' Column names of a multi-lag (`lags`) basis, lag-major:
#' `<v>_cb<j>_l<sl_token(lag)>`, j = 1..J within each lag. Also the column order
#' of the J x n_lag `col_matrix` stored in cb_meta (by column).
#' @noRd
.cb_lag_col_names <- function(v, J, lags) {
  as.vector(vapply(
    lags,
    function(l) sprintf("%s_cb%d_l%s", v, seq_len(J), sl_token(l)),
    character(J)
  ))
}

#' .replay_onebasis_lags
#'
#' Evaluate a stored onebasis (`ob_attrs`, resolved fit-time knots) on every
#' column of the exposure matrix `X` (n x n_lag) and bind the results lag-major
#' (n x J*n_lag). Fit time and replay (build_cb_columns_from_Q()) share it.
#' @noRd
.replay_onebasis_lags <- function(X, ob_attrs) {
  replay <- .onebasis_replay_args(ob_attrs)
  do.call(
    cbind,
    lapply(seq_len(ncol(X)), function(i) {
      m <- unclass(do.call(dlnm::onebasis, c(list(X[, i]), replay)))
      attributes(m) <- list(dim = dim(m))
      m
    })
  )
}

#' build_cb_block
#'
#' Fit-time construction of DLNM cross-basis predictors: one
#' dlnm::crossbasis() (distributed-lag) or dlnm::onebasis() (single-lag) per
#' variable in `dlnm_cfg$variables`, bound onto `dat` as flat
#' `<var>_cb1 … cbK` columns. The one loop of 04_fit_occ_brms_cb.R,
#' 04_fit_siz_brms_cb.R and 04_fit_dur_brms_cb.R; byte-identical to the
#' inline copies it replaces (guarded by test_cb_block_extraction.R).
#'
#' @param dat      data frame carrying the columns make_lag_matrix()/single_lag
#'                 lookup need for every variable in dlnm_cfg$variables
#'                 (enso_lag* for ENSO; `<var>`, `<var>_lag1..max_lag` for
#'                 climate-style vars).
#' @param dlnm_cfg one track's `dlnm:` config sub-block (`$enable`, `$defaults`,
#'                 `$variables`).
#' @param label    appended after "[dlnm]" in status messages/warnings for
#'                 caller context (e.g. " (monthly model)"); default "" matches
#'                 the pre-extraction message text exactly.
#' @return list(dat = <dat + bound basis columns>, cb_cols = character(),
#'              cb_meta = list()) — cb_meta has one entry per variable with
#'              the attributes extract_cb_surface()/build_cb_columns_from_Q()
#'              need for crosspred reconstruction.
build_cb_block <- function(dat, dlnm_cfg, label = "") {
  cb_cols <- character(0)
  cb_meta <- list()

  if (isTRUE(dlnm_cfg$enable) && length(dlnm_cfg$variables) > 0) {
    dlnm_defaults <- dlnm_cfg$defaults %||% list()

    for (v in names(dlnm_cfg$variables)) {
      var_spec <- dlnm_cfg$variables[[v]] %||% list()

      # Merge per-variable spec with shared defaults (per-var wins)
      spec <- dlnm_defaults
      for (k in names(var_spec)) {
        spec[[k]] <- var_spec[[k]]
      }

      max_lag <- as.integer(spec$max_lag %||% 3L)
      argvar <- spec$argvar %||% list(fun = "bs", degree = 3L, df = 4L)
      arglag <- spec$arglag %||% list(fun = "ns", df = 3L)
      single_lag <- spec$single_lag # NULL = full crossbasis; integer = onebasis at one lag

      lags <- .cb_lags_spec(spec$lags) # NULL unless the multi-lag mode is requested

      if (!is.null(lags)) {
        # --- Multi-lag mode: one onebasis per lag on a SHARED basis ---
        # The exposure is pooled over every lag column to fix knots/boundary
        # knots once, then the same basis is replayed on each lag column, so
        # column j means the same exposure function at every lag. Used by the
        # brms_kl backend (dlnm_kernel.R), which weights the lags with a
        # Gaussian kernel. Columns are lag-major (.cb_lag_col_names()).
        if (!is.null(single_lag)) {
          stop(
            "dlnm variable '",
            v,
            "': set either `single_lag:` or `lags:`, not both."
          )
        }
        message(sprintf(
          "[dlnm]%s var = %-20s  lags = %s..%s (%d)  argvar = %s  (shared onebasis per lag)",
          label,
          v,
          min(lags),
          max(lags),
          length(lags),
          paste(names(argvar), argvar, sep = "=", collapse = ",")
        ))
        X <- matrix(
          vapply(
            lags,
            function(l) .cb_lag_exposure(dat, v, l, who = "lags"),
            numeric(nrow(dat))
          ),
          nrow = nrow(dat)
        )
        ob0 <- do.call(dlnm::onebasis, c(list(as.vector(X)), argvar))
        ob_attrs <- attributes(ob0)
        ob_attrs[c("dim", "dimnames", "class")] <- NULL
        if (identical(ob_attrs$fun, "ns")) {
          ob_attrs$degree <- NULL
        }
        J <- ncol(ob0)
        B <- .replay_onebasis_lags(X, ob_attrs)
        cb_names <- .cb_lag_col_names(v, J, lags)
        colnames(B) <- cb_names
        dat <- dplyr::bind_cols(dat, tibble::as_tibble(as.data.frame(B)))
        cb_cols <- c(cb_cols, cb_names)
        cb_meta[[v]] <- list(
          lags = lags,
          onebasis_attrs = ob_attrs,
          argvar = argvar,
          range = range(X, na.rm = TRUE),
          df = J,
          cb_col_names = cb_names,
          col_matrix = matrix(
            cb_names,
            nrow = J,
            dimnames = list(NULL, sl_token(lags))
          )
        )
      } else if (!is.null(single_lag)) {
        # --- Single-lag mode: collapse the lag dimension via dlnm::onebasis() ---
        # Keeps the non-linear exposure-response (argvar); arglag/max_lag unused.
        L <- as.integer(single_lag)
        message(sprintf(
          "[dlnm]%s var = %-20s  single_lag = %d  argvar = %s  (onebasis; lag dimension collapsed)",
          label,
          v,
          L,
          paste(names(argvar), argvar, sep = "=", collapse = ",")
        ))

        x <- .cb_lag_exposure(dat, v, L, who = "single_lag")

        ob <- do.call(dlnm::onebasis, c(list(x), argvar))

        cb_names <- sprintf("%s_cb%d", v, seq_len(ncol(ob)))
        colnames(ob) <- cb_names
        dat <- dplyr::bind_cols(
          dat,
          tibble::as_tibble(as.data.frame(unclass(ob)))
        )
        cb_cols <- c(cb_cols, cb_names)

        # Store onebasis attributes (resolved knots etc.) for crosspred reconstruction.
        # Drop the spurious `degree` dlnm records for an `ns` basis so the stored
        # sidecar replays cleanly (see .onebasis_replay_args); harmless for `bs`.
        ob_attrs <- attributes(ob)
        ob_attrs[c("dim", "dimnames", "class")] <- NULL
        if (identical(ob_attrs$fun, "ns")) {
          ob_attrs$degree <- NULL
        }
        cb_meta[[v]] <- list(
          single_lag = L,
          onebasis_attrs = ob_attrs,
          argvar = argvar,
          range = range(x, na.rm = TRUE),
          df = ncol(ob),
          cb_col_names = cb_names
        )
      } else {
        message(sprintf(
          "[dlnm]%s var = %-20s  max_lag = %d  argvar = %s  arglag = %s",
          label,
          v,
          max_lag,
          paste(names(argvar), argvar, sep = "=", collapse = ","),
          paste(names(arglag), arglag, sep = "=", collapse = ",")
        ))

        # Build lag matrix from pre-computed columns (lag0 = var, lag1..max_lag)
        Q <- make_lag_matrix(dat, v, max_lag)

        # Compute crossbasis matrix
        cb <- dlnm::crossbasis(
          Q,
          lag = c(0L, max_lag),
          argvar = argvar,
          arglag = arglag
        )

        # Unique per-variable column prefix: <var>_cb1, <var>_cb2, ...
        cb_names <- sprintf("%s_cb%d", v, seq_len(ncol(cb)))
        colnames(cb) <- cb_names

        # Bind CB columns into dat
        dat <- dplyr::bind_cols(dat, tibble::as_tibble(as.data.frame(cb)))
        cb_cols <- c(cb_cols, cb_names)

        # Store metadata for postprocessing (attributes needed by crosspred)
        meta_attrs <- attributes(cb)
        cb_meta[[v]] <- list(
          argvar = meta_attrs$argvar,
          arglag = meta_attrs$arglag,
          lag = meta_attrs$lag,
          df = meta_attrs$df,
          range = meta_attrs$range,
          cb_col_names = cb_names
        )
      }
    }

    message(sprintf(
      "[dlnm]%s %d crossbasis column(s) across %d variable(s): %s",
      label,
      length(cb_cols),
      length(names(dlnm_cfg$variables)),
      paste(names(dlnm_cfg$variables), collapse = ", ")
    ))
  } else {
    if (!isTRUE(dlnm_cfg$enable)) {
      message(sprintf(
        "[dlnm]%s DLNM disabled (dlnm.enable: false). Running as standard model.",
        label
      ))
    } else {
      message(sprintf(
        "[dlnm]%s no variables configured; running as standard model.",
        label
      ))
    }
  }

  list(dat = dat, cb_cols = cb_cols, cb_meta = cb_meta)
}

#' build_cb_columns_from_Q
#'
#' Transform a single variable's lag matrix Q (n_obs × (max_lag + 1)) into its
#' DLNM cross-basis columns, reusing the stored fit-time basis exactly. The
#' resolved `argvar`/`arglag` in `meta` carry the knots dlnm computed at fit
#' time, so re-calling dlnm::crossbasis() reproduces the basis bit-for-bit (the
#' same guarantee extract_cb_surface() relies on for crosspred()).
#'
#' @param Q    numeric matrix n_obs × (max_lag + 1); columns lag0 … lag{max_lag}
#'             (output of make_lag_matrix()).
#' @param meta one entry of a cb_meta list: must contain argvar, arglag, lag,
#'             cb_col_names.
#' @return numeric matrix n_obs × length(meta$cb_col_names); colnames set to
#'         meta$cb_col_names.
#' .clamp_cb_range
#'
#' Hold values inside a fitted crossbasis' boundary knots (`cb_meta$<var>$range`)
#' before the basis is evaluated, so the cubic B-spline is never extrapolated:
#' effects beyond the training extremes are taken as equal to the response at
#' the extreme. Used by the smoothing phase (smoothing_core.R, enso_month_vec)
#' and, when `simulator.enso.clamp_to_knots` is set, by the forward simulation
#' (simulator_core.R, build_cb_mats). `range = NULL` is a no-op.
#'
#' `margin` widens each bound by that fraction of its own magnitude before
#' clamping (`lo - margin*|lo|`, `hi + margin*|hi|`), so the basis may
#' extrapolate a little beyond the training extremes (a natural-spline basis
#' continues linearly) and is held only beyond that. margin = 0 clamps at the
#' knots (the smoothing phase); the forward simulation uses
#' `simulator.enso.clamp_margin` (default 0.25, i.e. up to 125% of the max).
.clamp_cb_range <- function(x, range, margin = 0) {
  if (is.null(range)) {
    return(x)
  }
  lo <- range[1L] - margin * abs(range[1L])
  hi <- range[2L] + margin * abs(range[2L])
  pmin(pmax(x, lo), hi)
}

build_cb_columns_from_Q <- function(Q, meta) {
  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stop("build_cb_columns_from_Q: package 'dlnm' is required.")
  }
  if (!is.null(meta$lags)) {
    # Multi-lag (`lags`) basis: Q holds one exposure column per lag, in
    # meta$lags order (see build_cb_design()).
    Q <- if (is.matrix(Q)) Q else matrix(as.numeric(Q), ncol = 1L)
    if (ncol(Q) != length(meta$lags)) {
      stop(
        "build_cb_columns_from_Q: Q has ",
        ncol(Q),
        " column(s) but meta$lags has ",
        length(meta$lags),
        "."
      )
    }
    cb <- .replay_onebasis_lags(Q, meta$onebasis_attrs)
  } else if (!is.null(meta$single_lag)) {
    x <- if (is.matrix(Q)) Q[, 1L] else as.numeric(Q)
    ob_args <- .onebasis_replay_args(meta$onebasis_attrs)
    cb <- as.matrix(do.call(dlnm::onebasis, c(list(x = x), ob_args)))
  } else {
    cb <- as.matrix(dlnm::crossbasis(
      Q,
      lag = meta$lag,
      argvar = meta$argvar,
      arglag = meta$arglag
    ))
  }
  if (ncol(cb) != length(meta$cb_col_names)) {
    stop(
      "build_cb_columns_from_Q: reconstructed ",
      ncol(cb),
      " cross-basis columns but meta expects ",
      length(meta$cb_col_names),
      " (",
      paste(meta$cb_col_names, collapse = ", "),
      "). ",
      "argvar/arglag in cb_meta may be inconsistent with the fitted model."
    )
  }
  colnames(cb) <- meta$cb_col_names
  cb
}

#' build_mso_basis
#'
#' Build the months-since-onset cross-basis lookup table for the monthly
#' (`mon`) kernels. `mso` is latent (created per particle in the month loop),
#' so its basis is a `(mso_max + 1) x K` table indexed at run time by the
#' clamped value. Built with `build_cb_columns_from_Q()` from
#' `meta$onebasis_attrs`, so it is bit-identical to fit time.
#'
#' Clamping is required: the `mso` basis is bs(degree = 3), an unbounded cubic
#' past its boundary knots. The table stops at `ceiling(max(meta$range))` and
#' is indexed by `pmin(mso, mso_max) + 1L`, so mu flattens at the last trained
#' value. Do not cap it at the horizon length: a warm-started replicate can
#' enter the horizon with a large `mso`.
#'
#' @param cb_meta_mon the monthly track's cb_meta sidecar (artefact kind
#'   `mon_cb_meta`); must carry an `mso` entry, else NULL is returned.
#' @return `NULL` when the fit carries no `mso` cross-basis, else
#'   `list(B = <(mso_max + 1) x K matrix, colnames = meta$cb_col_names>,
#'         mso_max = <integer>)`. Row `j` of `B` is the basis at `mso = j - 1`.
build_mso_basis <- function(cb_meta_mon) {
  meta <- (cb_meta_mon %||% list())[["mso"]]
  if (is.null(meta)) {
    return(NULL)
  }
  rng <- meta$range
  if (is.null(rng) || !all(is.finite(rng))) {
    stop(
      "build_mso_basis: cb_meta$mso carries no finite `range`; the ",
      "boundary-knot clamp cannot be established. Re-run ",
      "analysis/04_fit_mon_brms_cb.R to regenerate the mon_cb_meta sidecar."
    )
  }
  mso_max <- max(1L, as.integer(ceiling(max(rng))))
  B <- build_cb_columns_from_Q(
    matrix(as.numeric(seq.int(0L, mso_max)), ncol = 1L),
    meta
  )
  list(B = B, mso_max = mso_max)
}

#' build_cb_design
#'
#' Reconstruct all DLNM cross-basis design columns (`<var>_cb1 … cbK`) for a
#' data frame of observations, from a cb_meta sidecar (written by
#' 04b_run_occ_brms_cb.R / 04d_run_siz_brms_cb.R). For each variable it builds
#' the lag matrix via make_lag_matrix() — identical to the fit-time code path — and
#' applies the stored basis. The variable name drives make_lag_matrix() dispatch:
#' "enso" reads the pre-built enso_lag* columns; any other name reads the
#' `<var>` / `<var>_lag*` climate lag columns (which must be present in `df`).
#'
#' @param df      data frame carrying the lag columns required by make_lag_matrix()
#'                for every variable in cb_meta (enso_lag* for ENSO;
#'                `<var>`, `<var>_lag1`…`<var>_lag{max_lag}` for climate vars).
#' @param cb_meta named list (one entry per cross-basis variable) or NULL.
#' @param vars    optional character vector restricting which variables to build
#'                (e.g. only "enso"); NULL = all names(cb_meta).
#' @return data frame with one column per cross-basis column (named by
#'         cb_col_names), nrow == nrow(df); NULL if cb_meta is empty / no vars.
build_cb_design <- function(df, cb_meta, vars = NULL) {
  if (is.null(cb_meta) || length(cb_meta) == 0) {
    return(NULL)
  }
  if (is.null(vars)) {
    vars <- names(cb_meta)
  }
  vars <- intersect(vars, names(cb_meta))
  if (length(vars) == 0) {
    return(NULL)
  }

  cols <- list()
  for (v in vars) {
    meta <- cb_meta[[v]]
    if (!is.null(meta$lags)) {
      # Multi-lag (`lags`) fit: one exposure column per lag, in meta$lags order.
      Q <- matrix(
        vapply(
          meta$lags,
          function(l) .cb_lag_exposure(df, v, l, who = "build_cb_design"),
          numeric(nrow(df))
        ),
        nrow = nrow(df)
      )
    } else if (!is.null(meta$single_lag)) {
      # Single-lag (onebasis) fit: the lag dimension was collapsed at fit time
      # onto a single exposure column (04_fit_occ_brms_cb.R). Reconstruct
      # from that exact column — NOT lag0 — mirroring build_cb_mats()
      # (simulator_core.R). build_cb_columns_from_Q() then reads
      # Q[,1L], so a one-column Q carrying the correct lag is required.
      Q <- matrix(
        .cb_lag_exposure(df, v, meta$single_lag, who = "build_cb_design"),
        ncol = 1L
      )
    } else {
      max_lag <- as.integer(meta$lag[2])
      Q <- make_lag_matrix(df, v, max_lag)
    }
    cb <- build_cb_columns_from_Q(Q, meta)
    cols[[v]] <- as.data.frame(cb)
  }
  dplyr::bind_cols(cols)
}
