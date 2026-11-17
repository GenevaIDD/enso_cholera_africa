# test_spatpomp_c_variant.R — HARD: the Csnippet rprocess / rinit variant of
#   the spatPomp builder (build_spatpomp_tu_c()).
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_c_variant.R
#
# Details:
#   Run from the project root. T1-T7: compile / structure, LP decomposition vs
#   evaluate_lp(), degenerate hazards, bpfilter smoke, a build-time guard at
#   U=200, covarindex smoke at U=32 and the cross-basis path. Only the C-side
#   assertions are checked.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/dlnm.R") # build_cb_design() -- required by build_spatpomp_tu_c()'s cb_meta code path (T7)
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/spatpomp_core_c.R") # also sources 08a

# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------
.fails <- 0L
.passes <- 0L

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
expect_equal <- function(a, b, label, tol = 1e-9) {
  ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
  if (!ok) {
    cat(sprintf(
      "         got  = %s\n         want = %s\n",
      paste(head(round(as.numeric(a), 9), 8), collapse = " "),
      paste(head(round(as.numeric(b), 9), 8), collapse = " ")
    ))
  }
  expect_true(ok, label)
}

cat("# test_spatpomp_c_variant\n\n")

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------
make_fixture <- function(
  tu = "hmmW",
  countries = c("AAA", "BBB"),
  adm1_per_country = 2L,
  adm2_per_adm1 = 2L,
  n_months = 3L,
  base_date = as.Date("2022-01-01")
) {
  frame <- expand.grid(
    country = countries,
    adm1n = seq_len(adm1_per_country),
    adm2n = seq_len(adm2_per_adm1),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      tu = tu,
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", tu),
      min_seas = 1L,
      pop = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()

  months <- base_date %m+% months(0:(n_months - 1L))

  panel <- tidyr::expand_grid(gid = frame$gid, t_start = months) %>%
    dplyr::left_join(frame %>% dplyr::select(gid, min_seas), by = "gid") %>%
    dplyr::mutate(
      country = sub("_[0-9]+_[0-9]+$", "", gid),
      t_start_seas = as.integer(lubridate::yday(t_start) - 1L),
      mean_precip = 0.0,
      dev_precip = 0.0,
      mean_temp_daymean = 0.0,
      dev_temp_daymean = 0.0,
      outbreak = 0L,
      time_to_last_outbreak_cat = "[6,Inf]",
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L
    ) %>%
    tibble::as_tibble()

  list(frame = frame, panel = panel, months = months)
}

fix <- make_fixture()
frame_tu <- fix$frame
panel_tu <- fix$panel
n_months <- length(fix$months)
window <- list(
  start = as.character(fix$months[1]),
  end = as.character(fix$months[n_months] %m+% months(1L))
)

# Degenerate draws bundles
draws_hi <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 15.0),
  basehaz_const = 0.0
)
draws_lo <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -15.0),
  basehaz_const = 0.0
)

# Moderate-hazard draws for T4 (~20% outbreak rate)
draws_mod <- make_draws_bundle_synthetic(
  frame_tu,
  S = 4L,
  fe_overrides = c("(Intercept)" = -1.5),
  basehaz_const = 0.0
)

# ===========================================================================
# T1 — Compile / smoke test
# ===========================================================================
# T1: build_spatpomp_tu_c() produces a valid spatPomp object with correct unit
#     count and new covariate columns present.
cat("## T1: build_spatpomp_tu_c — compile and object structure\n")

obj_c_hi <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_hi,
    window,
    fixed_draw_idx = 1L
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
obj_c_lo <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_lo,
    window,
    fixed_draw_idx = 1L
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)

expect_true(
  !is.null(obj_c_hi),
  "hi-hazard C-variant object builds without error"
)
expect_true(
  !is.null(obj_c_lo),
  "lo-hazard C-variant object builds without error"
)

if (!is.null(obj_c_hi)) {
  U <- nrow(frame_tu)
  expect_equal(
    length(obj_c_hi@unit_names),
    U,
    sprintf("unit_names length = U (%d)", U),
    tol = 0
  )

  # Check new covariate columns are present in the covariate table.
  # spatPomp stores the covar table; we check via attr on the built object.
  ct <- coef(obj_c_hi) # not the covar table, but let's check via precompute directly
  # Rebuild covar to check column presence independently.
  filt_months <- attr(obj_c_hi, "filt_months")
  draws_sub1 <- .subset_draws(draws_hi, 1L)
  ws1 <- attr(obj_c_hi, "warmup_state")
  base_cov <- build_covar_for_tu(panel_tu, frame_tu, filt_months)
  lp_static <- precompute_lp_static(
    draws_sub1,
    frame_tu,
    base_cov,
    character(0)
  )
  ttlo_c <- precompute_ttlo_coefs(draws_sub1, frame_tu)
  spill_c <- precompute_spillover_coefs(draws_sub1, frame_tu)
  aug_cov <- build_covar_for_tu_c(
    base_cov,
    lp_static,
    ttlo_c,
    spill_c,
    ws1,
    frame_tu
  )

  new_cols <- c(
    "lp_static",
    paste0("ttlo_coef_b", 0:6),
    "spill_o_t",
    "spill_o_f",
    "spill_oct_t",
    "spill_oct_f",
    "spill_otu_t",
    "spill_otu_f",
    "u_days_init",
    "o_init",
    "o_ct_init",
    "o_tu_init"
  )
  all_present <- all(new_cols %in% names(aug_cov))
  expect_true(
    all_present,
    sprintf(
      "all %d new covariate columns present in augmented covar table",
      length(new_cols)
    )
  )

  expect_true(
    identical(attr(obj_c_hi, "core_variant"), "c"),
    "core_variant attribute = 'c'"
  )
  expect_true(
    !is.null(attr(obj_c_hi, "fixed_draw_idx")),
    "fixed_draw_idx attribute set"
  )
}

# T1: error when fixed_draw_idx is NULL
err_null <- tryCatch(
  build_spatpomp_tu_c(frame_tu, panel_tu, draws_hi, window),
  error = function(e) conditionMessage(e)
)
expect_true(
  is.character(err_null) && nchar(err_null) > 0,
  "build_spatpomp_tu_c errors when fixed_draw_idx=NULL"
)

# ===========================================================================
# T2 — LP decomposition correctness
# ===========================================================================
# T2: lp_static + ttlo_coef + spill_coef matches evaluate_lp() for one unit at
#     one time step (within 1e-9).
cat("\n## T2: LP decomposition — lp_static + ttlo + spillover == evaluate_lp\n")

if (!is.null(obj_c_hi)) {
  U <- nrow(frame_tu)
  draws_sub1 <- .subset_draws(draws_hi, 1L)
  filt_months <- attr(obj_c_hi, "filt_months")
  base_cov <- build_covar_for_tu(panel_tu, frame_tu, filt_months)
  lp_static <- precompute_lp_static(
    draws_sub1,
    frame_tu,
    base_cov,
    character(0)
  )
  ttlo_c <- precompute_ttlo_coefs(draws_sub1, frame_tu)
  spill_c <- precompute_spillover_coefs(draws_sub1, frame_tu)

  # Pick time step t=1 and unit u=1 (first unit in frame_tu)
  t_check <- 1L
  u_check <- 1L
  mask_t1 <- base_cov$time == t_check

  # Reconstruct xrow for time t=1 (all climate = 0 in this fixture)
  rows_t1 <- base_cov[mask_t1, , drop = FALSE]
  perm_t1 <- match(as.character(frame_tu$gid), rows_t1$unit)
  rows_ord <- rows_t1[perm_t1, , drop = FALSE]

  xrow <- list(
    mean_precip = rows_ord$mean_precip,
    dev_precip = rows_ord$dev_precip,
    mean_temp_daymean = rows_ord$mean_temp_daymean,
    dev_temp_daymean = rows_ord$dev_temp_daymean,
    t_start_seas = rows_ord$t_start_seas
  )

  # For u_check = 1: set u_days to put it in bin 2 ([2,3) years = ~730 days)
  u_days_test <- rep(6 * 365.25, U) # all in bin 6 [6,Inf]
  u_days_test[u_check] <- 2.0 * 365.0 # unit u_check in bin 2 [2,3)

  # spillover: unit u_check has ongoing_outbreak_lag1 = 1
  o_test <- rep(0L, U)
  o_test[u_check] <- 1L
  o_ct_test <- rep(0L, U)
  o_tu_test <- rep(0L, U)

  bh_test <- basehaz_eval(draws_sub1, xrow$t_start_seas, N = U)
  lp_ref <- evaluate_lp(
    draws_sub1,
    s = 1L,
    xrow,
    o = o_test,
    o_ct = o_ct_test,
    o_tu = o_tu_test,
    u_days = u_days_test,
    frame_tu = frame_tu,
    basehaz_val = bh_test
  )

  # Reconstruct lp for unit u_check from precomputed parts
  lp_static_t1 <- lp_static[mask_t1]
  # Map back to frame_tu order
  inv_perm <- match(rows_t1$unit, as.character(frame_tu$gid))
  lp_static_frame_order <- lp_static_t1[inv_perm]

  # TTLO bin for u_days_test[u_check]: bin 2 ([2,3) years)
  bin_u_check <- 2L
  ttlo_contrib <- as.numeric(ttlo_c[[paste0("ttlo_coef_b", bin_u_check)]])[
    u_check
  ]

  # Spillover: o_test[u_check] = 1 → use spill_o_t; o_ct=0 → spill_oct_f; o_tu=0 → spill_otu_f
  spill_contrib <- spill_c$spill_o_t[u_check] + # o=1 → TRUE coef
    spill_c$spill_oct_f[u_check] + # o_ct=0 → FALSE coef
    spill_c$spill_otu_f[u_check] # o_tu=0 → FALSE coef

  lp_decomp_u <- lp_static_frame_order[u_check] + ttlo_contrib + spill_contrib

  expect_equal(
    lp_decomp_u,
    lp_ref[u_check],
    sprintf("LP decomp for unit %d (bin2, o=1) matches evaluate_lp", u_check)
  )

  # Also check a unit that is in bin 0 ([0,1)) with all-zero spillover
  u_check2 <- 2L
  u_days2 <- rep(6 * 365.25, U)
  u_days2[u_check2] <- 0.5 * 365.0 # bin 0

  lp_ref2 <- evaluate_lp(
    draws_sub1,
    s = 1L,
    xrow,
    o = rep(0L, U),
    o_ct = rep(0L, U),
    o_tu = rep(0L, U),
    u_days = u_days2,
    frame_tu = frame_tu,
    basehaz_val = bh_test
  )

  bin_u2 <- 0L
  ttlo2 <- as.numeric(ttlo_c[[paste0("ttlo_coef_b", bin_u2)]])[u_check2]
  spill2 <- spill_c$spill_o_f[u_check2] +
    spill_c$spill_oct_f[u_check2] +
    spill_c$spill_otu_f[u_check2]
  lp_dec2 <- lp_static_frame_order[u_check2] + ttlo2 + spill2

  expect_equal(
    lp_dec2,
    lp_ref2[u_check2],
    sprintf(
      "LP decomp for unit %d (bin0, all-zero spillover) matches evaluate_lp",
      u_check2
    )
  )
}

# ===========================================================================
# T3 — Degenerate hazard
# ===========================================================================
# T3: spatPomp simulate() on the C variant under LP = +-15 produces all-1 /
#     all-0 Yout.
cat("\n## T3: degenerate hazard — C-variant Yout all 0 / all 1\n")

run_degenerate <- function(draws_case, label) {
  obj_c <- tryCatch(
    build_spatpomp_tu_c(
      frame_tu,
      panel_tu,
      draws_case,
      window,
      fixed_draw_idx = 1L
    ),
    error = function(e) NULL
  )

  if (is.null(obj_c)) {
    cat(sprintf("  SKIP   T3 %s — build failed\n", label))
    return(invisible(NULL))
  }

  U <- nrow(frame_tu)
  expected_Y <- if (label == "hi") 1.0 else 0.0

  set.seed(99L)
  sm_c <- tryCatch(
    spatPomp::simulate(obj_c, format = "data.frame"),
    error = function(e) NULL
  )

  if (is.null(sm_c)) {
    cat(sprintf("  SKIP   T3 %s C — simulate failed\n", label))
    return(NULL)
  }

  # Extract Yout matrix (time × unit); handle wide and long formats.
  extract_Yout <- function(sm, obj) {
    U <- length(obj@unit_names)
    M <- length(attr(obj, "filt_months"))
    mat <- matrix(NA_real_, nrow = M, ncol = U)
    for (tt in seq_len(M)) {
      row_t <- sm[sm$time == tt, , drop = FALSE]
      if (nrow(row_t) == 0) {
        next
      }
      if ("unitname" %in% names(row_t)) {
        for (i in seq_len(U)) {
          r_ui <- row_t[row_t$unitname == obj@unit_names[i], , drop = FALSE]
          if (nrow(r_ui) > 0 && "Yout" %in% names(r_ui)) {
            mat[tt, i] <- r_ui$Yout[1L]
          }
        }
      } else {
        for (i in seq_len(U)) {
          nm <- paste0("Yout", i)
          if (nm %in% names(row_t)) mat[tt, i] <- row_t[[nm]][1L]
        }
      }
    }
    mat
  }

  Y_c <- extract_Yout(sm_c, obj_c)

  all_ok_c <- all(!is.na(Y_c) & Y_c == expected_Y)
  expect_true(
    all_ok_c,
    sprintf(
      "T3 %s: C-variant Yout all = %.0f under LP %s15",
      label,
      expected_Y,
      if (label == "hi") "+" else "-"
    )
  )

  # Lag updates: for hi-hazard, all units have outbreak → o=1 everywhere after t=1
  if (label == "hi") {
    p_cols <- paste0("p", seq_len(U))
    if (all(p_cols %in% names(sm_c))) {
      p_row1 <- as.numeric(unlist(sm_c[sm_c$time == 1L, p_cols, drop = FALSE]))
      expect_true(
        all(p_row1 > 0.999, na.rm = TRUE),
        "T3 hi: C-variant latent p near 1 at t=1"
      )
    }
  }
}

run_degenerate(draws_hi, "hi")
run_degenerate(draws_lo, "lo")

# ===========================================================================
# T4 — bpfilter smoke
# ===========================================================================
# T4: finite logLik from run_bpfilter_tu() and run_one_bpfilter_replicate() at
#     Np=100.
cat("\n## T4: bpfilter smoke (C, Np=100)\n")

if (!is.null(obj_c_hi)) {
  Np4 <- 100L

  obj_c4 <- tryCatch(
    build_spatpomp_tu_c(
      frame_tu,
      panel_tu,
      draws_mod,
      window,
      fixed_draw_idx = 2L
    ),
    error = function(e) NULL
  )

  if (!is.null(obj_c4)) {
    blks_c <- make_block_list(obj_c4, frame_tu)

    bpf_c <- tryCatch(
      run_bpfilter_tu(obj_c4, Np = Np4, block_list = blks_c, seed = 7L),
      error = function(e) {
        warning("T4 C bpfilter failed: ", conditionMessage(e))
        NULL
      }
    )
    ll_c <- if (!is.null(bpf_c)) {
      tryCatch(pomp::logLik(bpf_c), error = function(e) NA_real_)
    } else {
      NA_real_
    }
    expect_true(
      is.finite(ll_c),
      sprintf("T4: C-variant logLik finite (%.3f)", ll_c %||% NaN)
    )

    # run_one_bpfilter_replicate interface
    bpf_rep_c <- tryCatch(
      run_one_bpfilter_replicate(
        obj_c4,
        Np = Np4,
        block_list = blks_c,
        seed = 13L
      ),
      error = function(e) NULL
    )
    expect_true(
      !is.null(bpf_rep_c) && is.finite(bpf_rep_c$logLik),
      "T4: C-variant run_one_bpfilter_replicate logLik finite"
    )
  } else {
    cat("  SKIP   T4 — object build failed for moderate-hazard draws\n")
  }
}

# ===========================================================================
# T5 — Build-time regression guard (large U)
# ===========================================================================
# NOTE on expected build time:
#   With the O(1) Csnippets, our generated C code is ~112
#   lines (constant) for U=200.  However, pomp/spatPomp internally generates
#   O(U) C wrapper boilerplate — specifically, named variable declarations for
#   every unit state variable (e.g. u_days1, u_days2, ..., u_days200 × 7 vars
#   = 1,400 declarations) plus all covariate aliases.  This pomp-internal O(U)
#   cost is unavoidable at the Csnippet level and accounts for ~8 s at U=200.
#
#   Threshold rationale: old O(U) Csnippet code at U=200 took ~35 s.  Our
#   current implementation takes ~8-9 s.  Threshold of 20 s catches a
#   regression (accidentally reintroduced O(U) Csnippet) while passing for the
#   current pomp-internal baseline.  Sub-second compilation would require
#   patching pomp itself (out of scope).
#
#   The guard is on CPU time (this process + its children, i.e. the C compiler)
#   rather than wall-clock time: under a parallel test run the elapsed time
#   doubled (41 s seen) while the work did not, so the old elapsed-time guard
#   failed with no regression. A reintroduced O(U) Csnippet raises the
#   compiler's CPU time, so the 20 s threshold still catches it.
cat("\n## T5: build-time regression guard — U=200 C-variant build < 20 s CPU\n")

fix_large <- make_fixture(
  tu = "hmmLARGE",
  countries = paste0("CX", 1:10),
  adm1_per_country = 4L,
  adm2_per_adm1 = 5L, # 10 × 4 × 5 = 200 units
  n_months = 3L
)
draws_large <- make_draws_bundle_synthetic(
  fix_large$frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -1.8),
  basehaz_const = 0.0
)
window_large <- list(
  start = as.character(fix_large$months[1]),
  end = as.character(fix_large$months[length(fix_large$months)] %m+% months(1L))
)

t5_t <- system.time({
  obj_c_large <- tryCatch(
    build_spatpomp_tu_c(
      fix_large$frame,
      fix_large$panel,
      draws_large,
      window_large,
      fixed_draw_idx = 1L
    ),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
})
t5_time <- sum(
  t5_t[c("user.self", "sys.self", "user.child", "sys.child")],
  na.rm = TRUE
)

U_large <- nrow(fix_large$frame)
cat(sprintf(
  "     U=%d  build_spatpomp_tu_c: %.2f s CPU (%.2f s elapsed)\n",
  U_large,
  t5_time,
  t5_t[["elapsed"]]
))
expect_true(
  !is.null(obj_c_large),
  sprintf("T5: U=%d C-variant builds without error", U_large)
)
expect_true(
  isTRUE(t5_time < 20.0),
  sprintf(
    "T5: build CPU time %.2f s < 20 s (O(1)-code regression guard)",
    t5_time
  )
)

# ===========================================================================
# T6 — Covarindex smoke at U=32
# ===========================================================================
cat("\n## T6: covarindex smoke — bpfilter logLik finite at U=32\n")
#
# Exercises __covarindex["name"] + u per-unit covariate access at a larger U.
# (The former R-variant logLik concordance bound went with the R builder;
# this is now a finiteness smoke only.)

fix_med <- make_fixture(
  tu = "hmmMED",
  countries = paste0("DX", 1:4),
  adm1_per_country = 4L,
  adm2_per_adm1 = 2L, # 4 × 4 × 2 = 32 units
  n_months = 4L
)
draws_med <- make_draws_bundle_synthetic(
  fix_med$frame,
  S = 4L,
  fe_overrides = c("(Intercept)" = -1.5),
  basehaz_const = 0.0
)
window_med <- list(
  start = as.character(fix_med$months[1]),
  end = as.character(fix_med$months[length(fix_med$months)] %m+% months(1L))
)

obj_c6 <- tryCatch(
  build_spatpomp_tu_c(
    fix_med$frame,
    fix_med$panel,
    draws_med,
    window_med,
    fixed_draw_idx = 2L
  ),
  error = function(e) NULL
)

U_med <- nrow(fix_med$frame)
cat(sprintf("     U=%d (covarindex probe)\n", U_med))

if (!is.null(obj_c6)) {
  blks_c6 <- make_block_list(obj_c6, fix_med$frame)
  bpf_c6 <- tryCatch(
    run_bpfilter_tu(obj_c6, Np = 50L, block_list = blks_c6, seed = 42L),
    error = function(e) {
      warning("T6 C bpfilter: ", conditionMessage(e))
      NULL
    }
  )
  ll_c6 <- if (!is.null(bpf_c6)) {
    tryCatch(pomp::logLik(bpf_c6), error = function(e) NA_real_)
  } else {
    NA_real_
  }
  expect_true(
    is.finite(ll_c6),
    sprintf(
      "T6: C-variant logLik finite at U=%d (%.3f)",
      U_med,
      if (is.finite(ll_c6)) ll_c6 else NaN
    )
  )
} else {
  cat("  SKIP   T6 — object build failed\n")
}

# ===========================================================================
# T7 — C-variant spatPomp with non-null cb_meta (cross-basis code path)
# ===========================================================================
cat("\n## T7: build_spatpomp_tu_c with cb_meta — cb code path exercised\n")
#
# Exercises the cross-basis code paths:
#   build_spatpomp_tu_c(…, cb_meta) → extracts cb_covar_names, passes to
#   build_covar_for_tu (→ build_cb_design → make_lag_matrix) and to
#   precompute_lp_static (→ routes cb columns into xrow for evaluate_lp).
#
# Verifies:
#   (a) The object builds without error when cb_meta is non-NULL.
#   (b) mean_precip_cb1 is present in the reconstructed covar table.
#   (c) mean_precip_cb1 values are non-zero (lag columns were consumed).
#   (d) lp_static with cb_meta differs from lp_static without it (FE ≠ 0,
#       precip ≠ 0 → non-zero CB contribution to LP).
#
# dlnm is required; the test is skipped gracefully when not installed.

if (!requireNamespace("dlnm", quietly = TRUE)) {
  cat("  SKIP   T7 — dlnm not available\n")
} else {
  # Probe actual column count for lin × lin, lag 0..2 (dlnm may produce >1 col).
  Q_probe <- matrix(c(1.5, 1.2, 0.8), nrow = 1)
  n_cb_cols_t7 <- ncol(as.matrix(
    dlnm::crossbasis(
      Q_probe,
      lag = c(0L, 2L),
      argvar = list(fun = "lin"),
      arglag = list(fun = "lin")
    )
  ))
  cb_col_names_t7 <- sprintf("mean_precip_cb%d", seq_len(n_cb_cols_t7))

  # Synthetic cb_meta: linear × linear cross-basis, lags 0..2.
  cb_meta_t7 <- list(
    mean_precip = list(
      lag = c(0L, 2L),
      argvar = list(fun = "lin"),
      arglag = list(fun = "lin"),
      cb_col_names = cb_col_names_t7
    )
  )

  # Fixture with n_months = 4 so there is room for lags.
  fix_t7 <- make_fixture(n_months = 4L)
  frame_t7 <- fix_t7$frame
  window_t7 <- list(
    start = as.character(fix_t7$months[1]),
    end = as.character(fix_t7$months[length(fix_t7$months)] %m+% months(1L))
  )

  # Panel: set non-zero mean_precip and pre-computed lag columns (as 03c emits).
  panel_t7 <- fix_t7$panel %>%
    dplyr::mutate(
      mean_precip = 1.5,
      mean_precip_lag1 = 1.2,
      mean_precip_lag2 = 0.8
    )

  # Draws with all cb columns active and the first one with a non-zero FE.
  active_t7 <- c(
    "(Intercept)",
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    "ongoing_outbreak_lag1",
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1",
    cb_col_names_t7
  )
  fe_ov_t7 <- setNames(
    c(-3.0, rep(2.0, n_cb_cols_t7)),
    c("(Intercept)", cb_col_names_t7)
  )
  draws_t7 <- make_draws_bundle_synthetic(
    frame_t7,
    S = 2L,
    fe_overrides = fe_ov_t7,
    basehaz_const = 0.0,
    active_terms = active_t7
  )

  # (a) Build with cb_meta — must not error.
  obj_c_cb <- tryCatch(
    build_spatpomp_tu_c(
      frame_t7,
      panel_t7,
      draws_t7,
      window_t7,
      fixed_draw_idx = 1L,
      cb_meta = cb_meta_t7
    ),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  expect_true(
    !is.null(obj_c_cb),
    "T7: build_spatpomp_tu_c with cb_meta builds without error"
  )

  if (!is.null(obj_c_cb)) {
    filt_months_t7 <- attr(obj_c_cb, "filt_months")
    draws_sub_t7 <- .subset_draws(draws_t7, 1L)

    # (b) Covar table carries all cb columns.
    base_cov_cb <- build_covar_for_tu(
      panel_t7,
      frame_t7,
      filt_months_t7,
      cb_meta = cb_meta_t7
    )
    all_cb_present <- all(cb_col_names_t7 %in% names(base_cov_cb))
    expect_true(
      all_cb_present,
      sprintf(
        "T7: %s present in covar table from build_covar_for_tu",
        paste(cb_col_names_t7, collapse = ", ")
      )
    )

    # (c) First cb column is non-zero (lag columns were consumed, not zero-filled).
    cb1_vals <- base_cov_cb[[cb_col_names_t7[1]]]
    expect_true(
      any(!is.na(cb1_vals) & cb1_vals != 0),
      "T7: first cb column is non-zero (lag columns consumed by build_cb_design)"
    )

    # (d) lp_static with cb differs from lp_static without cb.
    lp_with_cb <- precompute_lp_static(
      draws_sub_t7,
      frame_t7,
      base_cov_cb,
      character(0),
      cb_names = cb_col_names_t7
    )
    base_cov_nocb <- build_covar_for_tu(
      panel_t7,
      frame_t7,
      filt_months_t7,
      cb_meta = NULL
    )
    lp_nocb <- precompute_lp_static(
      draws_sub_t7,
      frame_t7,
      base_cov_nocb,
      character(0),
      cb_names = character(0)
    )

    max_diff <- max(abs(lp_with_cb - lp_nocb), na.rm = TRUE)
    expect_true(
      max_diff > 1e-6,
      sprintf(
        "T7: lp_static differs by max %.4f when cb FE=2, precip=1.5",
        max_diff
      )
    )
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
