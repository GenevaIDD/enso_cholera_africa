# test_spatpomp_c_variant_mon.R — HARD: smoke / structural checks of the
#   monthly-cases (`mon`) spatPomp Csnippet kernel (spatpomp_core_c.R).
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_c_variant_mon.R
#
# Details:
#   Run from the project root. mon hard-requires a duration fit (`active` /
#   `mso` are duration-model latent states), so every fixture carries
#   draws_dur alongside draws_mon. T1 compile / structure, T2 degenerate
#   trajectories, T3-T4 bpfilter smoke (ungrouped, grouped), T5 the
#   covariate-name collision fix.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/dlnm.R") # build_mso_basis()
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R

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
      paste(head(round(as.numeric(a), 9), 12), collapse = " "),
      paste(head(round(as.numeric(b), 9), 12), collapse = " ")
    ))
  }
  expect_true(ok, label)
}

cat(
  "# test_spatpomp_c_variant_mon — monthly-size (mon) C kernel smoke/structural checks\n\n"
)

# ---------------------------------------------------------------------------
# Shared fixture (same helper as test_spatpomp_c_variant_dur.R)
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

#' Extract a state-variable matrix (time x unit) from spatPomp::simulate()
#' output, handling both long (unitname) and wide (var<i>) formats.
extract_state <- function(sm, obj, varname) {
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
        if (nrow(r_ui) > 0 && varname %in% names(r_ui)) {
          mat[tt, i] <- r_ui[[varname]][1L]
        }
      }
    } else {
      for (i in seq_len(U)) {
        nm <- paste0(varname, i)
        if (nm %in% names(row_t)) mat[tt, i] <- row_t[[nm]][1L]
      }
    }
  }
  mat
}

fix <- make_fixture()
frame_tu <- fix$frame
panel_tu <- fix$panel
n_months <- length(fix$months)
window <- list(
  start = as.character(fix$months[1]),
  end = as.character(fix$months[n_months] %m+% months(1L))
)

# Degenerate onset draws (occurrence): near-certain onset every month.
draws_occ_hi <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 15.0),
  basehaz_const = 0.0
)

# mso cross-basis fixture: dlnm::onebasis(fun="lin") -> exactly one column
# ("mso_cb1"), matching draws_mon$fe below. Same fixture shape as
# test_spatpomp_recovery_mon.R -- sufficient to exercise the mso lookup table
# end-to-end without pulling in spline boundary-knot machinery this test does
# not need to re-verify (that's test_monthly_kernel.R's job).
cb_meta_mon <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, n_months + 5L),
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)

make_mon_bundle <- function(S = 2L) {
  b <- make_draws_bundle_synthetic(
    frame_tu,
    S = S,
    active_terms = c("(Intercept)", "mso_cb1")
  )
  b$fe[, "(Intercept)"] <- log(20) # ~20 cases/month at onset (mso=0)
  b$fe[, "mso_cb1"] <- -0.15 # declining mean as the episode ages
  b$phi <- rep(8.0, S)
  b$offset_pop <- FALSE
  b$offset_logoffset <- FALSE
  b
}
draws_mon <- make_mon_bundle(S = 2L)

# ===========================================================================
# T1 — Compile / smoke test
# ===========================================================================
# T1: build_spatpomp_tu_c(..., draws_mon=, cb_meta_mon=) produces a valid
#     spatPomp object with draws_mon_sub/mso_basis attributes set and all
#     mon_* covariate columns present; backward compat (draws_dur without
#     draws_mon) is unaffected.
cat(
  "## T1: build_spatpomp_tu_c(..., draws_mon=, cb_meta_mon=) — compile and object structure\n"
)

draws_dur_never <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -30.0),
  active_terms = c("(Intercept)")
)

obj_c_mon <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_occ_hi,
    window,
    fixed_draw_idx = 1L,
    draws_dur = draws_dur_never,
    draws_mon = draws_mon,
    cb_meta_mon = cb_meta_mon,
    measurement = "occurrence_duration_mon_NB"
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)

expect_true(
  !is.null(obj_c_mon),
  "C-variant object with draws_mon builds without error"
)

if (!is.null(obj_c_mon)) {
  expect_true(
    !is.null(attr(obj_c_mon, "draws_mon_sub")),
    "draws_mon_sub attribute set"
  )
  expect_true(!is.null(attr(obj_c_mon, "mso_basis")), "mso_basis attribute set")

  U <- nrow(frame_tu)
  filt_months <- attr(obj_c_mon, "filt_months")
  ws1 <- attr(obj_c_mon, "warmup_state")
  draws_sub1 <- .subset_draws(draws_occ_hi, 1L)
  draws_dur_sub1 <- .subset_draws(draws_dur_never, 1L)
  draws_mon_sub1 <- .subset_draws(draws_mon, 1L)
  mso_basis <- attr(obj_c_mon, "mso_basis")
  mon_covar_map <- attr(obj_c_mon, "mon_covar_map")

  base_cov <- build_covar_for_tu(panel_tu, frame_tu, filt_months)
  lp_static <- precompute_lp_static(
    draws_sub1,
    frame_tu,
    base_cov,
    character(0)
  )
  ttlo_c <- precompute_ttlo_coefs(draws_sub1, frame_tu)
  spill_c <- precompute_spillover_coefs(draws_sub1, frame_tu)
  dur_lp <- precompute_dur_lp_static(
    draws_dur_sub1,
    frame_tu,
    base_cov,
    character(0)
  )
  dur_ttlo <- precompute_dur_ttlo_coefs(draws_dur_sub1, frame_tu)
  dur_spill <- precompute_dur_spillover_coefs(draws_dur_sub1, frame_tu)
  mon_comp <- precompute_mon_lp_components(
    frame_tu,
    base_cov,
    mon_covar_map,
    draws_mon_sub1,
    mso_basis = mso_basis
  )

  aug_cov <- build_covar_for_tu_c(
    base_cov,
    lp_static,
    ttlo_c,
    spill_c,
    ws1,
    frame_tu,
    dur_lp_static_vec = dur_lp,
    dur_ttlo_coefs_df = dur_ttlo,
    dur_spill_coefs_df = dur_spill,
    mon_lp_static_vec = mon_comp$lp,
    mon_ttlo_coefs_df = mon_comp$ttlo,
    mon_spill_coefs_df = mon_comp$spill,
    mon_mso_coefs_df = mon_comp$mso,
    size_phi_vec = mon_comp$phi
  )

  mon_cols <- c(
    "mon_lp_static",
    paste0("mon_ttlo_coef_b", 0:6),
    "mon_spill_o_t",
    "mon_spill_o_f",
    "mon_spill_oct_t",
    "mon_spill_oct_f",
    "mon_spill_otu_t",
    "mon_spill_otu_f",
    "mon_mso_coef_j0",
    "mon_mso_coef_j1",
    "size_phi_init"
  )
  all_present <- all(mon_cols %in% names(aug_cov))
  expect_true(
    all_present,
    sprintf(
      "all %d mon_* / mso-coefficient covariate columns present",
      length(mon_cols)
    )
  )
}

# T1: backward compat — build_spatpomp_tu_c() WITH draws_dur but WITHOUT
# draws_mon is unaffected.
obj_c_nomon <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_occ_hi,
    window,
    fixed_draw_idx = 1L,
    draws_dur = draws_dur_never
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(
  !is.null(obj_c_nomon),
  "C-variant object without draws_mon still builds"
)
expect_true(
  is.null(attr(obj_c_nomon, "draws_mon_sub")),
  "draws_mon_sub attribute NULL when draws_mon omitted"
)

# ===========================================================================
# T2 — Degenerate state trajectories (ungrouped)
# ===========================================================================
# T2: degenerate state trajectories (ungrouped) under the active/mso scenarios
#     of test_spatpomp_c_variant_dur.R, extended to the mon `mu` state:
#       A) onset certain, cessation impossible -> active stays 1, mso
#          increments every month -> mu > 0 every step (mon_mu_block's mso
#          lookup fed by a live, incrementing latent state).
#       B) onset certain, cessation certain (same-step) -> mu is read from the
#          PRE-cessation mso (mon_mu_block runs before the cessation block),
#          so mu > 0 even though active ends every step at 0.
cat("\n## T2: degenerate active/mso/mu trajectories — C-variant\n")

run_scenario <- function(draws_dur_case, draws_mon_case, label) {
  obj_c <- tryCatch(
    build_spatpomp_tu_c(
      frame_tu,
      panel_tu,
      draws_occ_hi,
      window,
      fixed_draw_idx = 1L,
      draws_dur = draws_dur_case,
      draws_mon = draws_mon_case,
      cb_meta_mon = cb_meta_mon,
      measurement = "occurrence_duration_mon_NB"
    ),
    error = function(e) {
      cat("  ERROR (C):", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(obj_c)) {
    cat(sprintf("  SKIP   T2 %s — build failed\n", label))
    return(invisible(NULL))
  }

  set.seed(2024L)
  sm_c <- tryCatch(
    spatPomp::simulate(obj_c, format = "data.frame"),
    error = function(e) NULL
  )

  if (is.null(sm_c)) {
    cat(sprintf("  SKIP   T2 %s C — simulate failed\n", label))
    return(NULL)
  }

  invisible(list(sm_c = sm_c, obj_c = obj_c))
}

# Scenario A: onset certain, cessation impossible -> active=1 mso=1,2,3...,
# mu > 0 every month and tracks the (declining) mso_cb1 slope.
draws_dur_A <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -30.0),
  active_terms = c("(Intercept)")
)
res_A <- run_scenario(draws_dur_A, draws_mon, "A(never-cease)")
expect_true(!is.null(res_A), "T2 A: C-variant builds and simulates")

if (!is.null(res_A)) {
  mat_active_c <- extract_state(res_A$sm_c, res_A$obj_c, "active")
  mat_mso_c <- extract_state(res_A$sm_c, res_A$obj_c, "mso")
  mat_mu_c <- extract_state(res_A$sm_c, res_A$obj_c, "mu")
  expect_true(
    all(mat_active_c == 1),
    "T2 A: C-variant active = 1 for all (unit, month)"
  )
  expect_equal(
    mat_mso_c[, 1],
    seq_len(n_months),
    "T2 A: C-variant mso increments 1,2,3,... (no re-onset, always active)"
  )
  expect_true(
    all(mat_mu_c > 0),
    "T2 A: C-variant mu > 0 on every active (unit, month)"
  )
}

# Scenario B: onset certain, cessation certain (within onset month) -> active
# pulses to 1 then back to 0 every step; mu is read from the PRE-cessation
# mso (mon_mu_block runs before dur_cessation_block), so mu > 0 despite
# active ending every step at 0 (the "mu>0 exactly matches active-DURING-the-
# month, not active-at-end-of-step" contract already pinned by
# test_spatpomp_recovery_mon.R's Step-1 assertions).
draws_dur_B <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 30.0),
  active_terms = c("(Intercept)")
)
res_B <- run_scenario(draws_dur_B, draws_mon, "B(always-cease)")
expect_true(!is.null(res_B), "T2 B: C-variant builds and simulates")

if (!is.null(res_B)) {
  mat_yout_c <- extract_state(res_B$sm_c, res_B$obj_c, "Yout")
  mat_active_c <- extract_state(res_B$sm_c, res_B$obj_c, "active")
  mat_mso_c <- extract_state(res_B$sm_c, res_B$obj_c, "mso")
  mat_mu_c <- extract_state(res_B$sm_c, res_B$obj_c, "mu")
  expect_true(
    all(mat_yout_c == 1),
    "T2 B: C-variant Yout = 1 every month (re-onset every step)"
  )
  expect_true(
    all(mat_active_c == 0),
    "T2 B: C-variant active = 0 at end of every step (cessation same-step)"
  )
  expect_true(
    all(mat_mso_c == 0),
    "T2 B: C-variant mso = 0 at end of every step (reset on cessation)"
  )
  expect_true(
    all(mat_mu_c > 0),
    "T2 B: C-variant mu > 0 every step (mon_mu_block reads PRE-cessation mso)"
  )
}

# ===========================================================================
# T3 — bpfilter smoke (ungrouped)
# ===========================================================================
# T3: bpfilter smoke (ungrouped): finite logLik through the mon lp_static/
#     ttlo/spillover covariate wiring plus the mso coefficient table
#     (precompute_mon_mso_coefs / mon_mso_coef_j*) under moderate hazards,
#     unconditional (no obs_data).
cat("\n## T3: bpfilter smoke, ungrouped\n")

draws_occ_mod <- make_draws_bundle_synthetic(
  frame_tu,
  S = 4L,
  fe_overrides = c("(Intercept)" = -1.0),
  basehaz_const = 0.0
)

draws_dur_mod <- make_draws_bundle_synthetic(
  frame_tu,
  S = 4L,
  active_terms = c(
    "(Intercept)",
    "mso_band2",
    "mso_band3",
    "mso_band4M5",
    "mso_band6M8",
    "mso_band9P"
  )
)
draws_dur_mod$fe[, "(Intercept)"] <- -0.8
draws_dur_mod$fe[, "mso_band2"] <- 0.4
draws_dur_mod$fe[, "mso_band3"] <- -0.3
draws_dur_mod$fe[, "mso_band4M5"] <- 0.6
draws_dur_mod$fe[, "mso_band6M8"] <- -0.5
draws_dur_mod$fe[, "mso_band9P"] <- 0.2
draws_dur_mod$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)

draws_mon_mod <- make_mon_bundle(S = 4L)

fix4 <- make_fixture(n_months = 4L)
frame_tu4 <- fix4$frame
panel_tu4 <- fix4$panel
window4 <- list(
  start = as.character(fix4$months[1]),
  end = as.character(fix4$months[length(fix4$months)] %m+% months(1L))
)
cb_meta_mon4 <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, 4L + 5L),
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)

obj_c3 <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu4,
    panel_tu4,
    draws_occ_mod,
    window4,
    fixed_draw_idx = 2L,
    draws_dur = draws_dur_mod,
    draws_mon = draws_mon_mod,
    cb_meta_mon = cb_meta_mon4,
    measurement = "occurrence_duration_mon_NB"
  ),
  error = function(e) {
    cat("  ERROR (C):", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(!is.null(obj_c3), "T3: C-variant (mon) builds")

if (!is.null(obj_c3)) {
  blks_c3 <- make_block_list(obj_c3, frame_tu4)
  bpf_c3 <- tryCatch(
    run_bpfilter_tu(obj_c3, Np = 100L, block_list = blks_c3, seed = 7L),
    error = function(e) {
      warning("T3 C bpfilter: ", conditionMessage(e))
      NULL
    }
  )
  ll_c3 <- if (!is.null(bpf_c3)) {
    tryCatch(pomp::logLik(bpf_c3), error = function(e) NA_real_)
  } else {
    NA_real_
  }
  expect_true(
    is.finite(ll_c3),
    sprintf("T3: C-variant logLik finite (%.3f)", ll_c3 %||% NaN)
  )
} else {
  cat("  SKIP   T3 — object build failed\n")
}

# ===========================================================================
# T4 — bpfilter smoke (GROUPED)
# ===========================================================================
# T4: bpfilter smoke (GROUPED): grouping_mon_sub set and finite logLik through
#     grouping$draws_mon_by_group / grouping$cb_meta_mon_by_group
#     (precompute_mon_lp_components()'s per-group dispatch and the
#     padded-to-common-width mso table stitching).
cat("\n## T4: bpfilter smoke, GROUPED\n")
# Exercises grouping$draws_mon_by_group / grouping$cb_meta_mon_by_group ->
# precompute_mon_lp_components()'s per-group dispatch and the mso table
# stitching (each group padded to the shared max width by repeating its own
# clamped last column).

make_2group_frame <- function() {
  tidyr::expand_grid(country = c("AAAW", "BBBW"), adm1n = 1:2, adm2n = 1:2) %>%
    dplyr::mutate(
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = "REG_W",
      tu = "hmmW",
      min_seas = 1L,
      pop = 1e5,
      group = dplyr::if_else(country == "AAAW", "ga", "gb")
    ) %>%
    dplyr::select(
      gid,
      country,
      AFRO_region,
      tu,
      min_seas,
      pop,
      adm1_unit,
      group
    ) %>%
    tibble::as_tibble()
}

make_region_dropped_occ_bundle <- function(frame_g, intercept = 0, S = 2L) {
  fe_nm <- c(
    "(Intercept)",
    "mean_precip",
    "dev_precip",
    "mean_temp_daymean",
    "dev_temp_daymean",
    "ongoing_outbreak_lag1",
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1"
  )
  zm <- function(cols) {
    m <- matrix(0, nrow = S, ncol = length(cols))
    colnames(m) <- cols
    m
  }
  fe <- zm(fe_nm)
  fe[, "(Intercept)"] <- intercept
  list(
    S = S,
    fe = fe,
    re_intercept = list(
      region = NULL,
      country = zm(unique(as.character(frame_g$country))),
      adm1 = zm(unique(paste(frame_g$country, frame_g$adm1_unit, sep = "_"))),
      adm2 = NULL
    ),
    re_slope = list(),
    basehaz = list(kind = "const", const = 0),
    phi = rep(1.0, S),
    active_terms = fe_nm
  )
}

make_region_dropped_dur_bundle <- function(frame_g, fe_vals, S = 2L) {
  fe_nm <- names(fe_vals)
  zm <- function(cols) {
    m <- matrix(0, nrow = S, ncol = length(cols))
    colnames(m) <- cols
    m
  }
  fe <- zm(fe_nm)
  for (nm in fe_nm) {
    fe[, nm] <- fe_vals[[nm]]
  }
  list(
    S = S,
    fe = fe,
    re_intercept = list(
      region = NULL,
      country = zm(unique(as.character(frame_g$country))),
      adm1 = zm(unique(paste(frame_g$country, frame_g$adm1_unit, sep = "_"))),
      adm2 = NULL
    ),
    re_slope = list(),
    basehaz = list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS),
    phi = rep(1.0, S),
    active_terms = fe_nm
  )
}

make_region_dropped_mon_bundle <- function(frame_g, fe_vals, S = 2L) {
  fe_nm <- names(fe_vals)
  zm <- function(cols) {
    m <- matrix(0, nrow = S, ncol = length(cols))
    colnames(m) <- cols
    m
  }
  fe <- zm(fe_nm)
  for (nm in fe_nm) {
    fe[, nm] <- fe_vals[[nm]]
  }
  list(
    S = S,
    fe = fe,
    re_intercept = list(
      region = NULL,
      country = zm(unique(as.character(frame_g$country))),
      adm1 = zm(unique(paste(frame_g$country, frame_g$adm1_unit, sep = "_"))),
      adm2 = NULL
    ),
    re_slope = list(),
    basehaz = list(kind = "const", const = 0),
    phi = rep(8.0, S),
    active_terms = fe_nm,
    offset_pop = FALSE,
    offset_logoffset = FALSE
  )
}

ft_grp <- make_2group_frame()
frame_ga <- ft_grp[ft_grp$group == "ga", , drop = FALSE]
frame_gb <- ft_grp[ft_grp$group == "gb", , drop = FALSE]

draws_occ_ga <- make_region_dropped_occ_bundle(
  frame_ga,
  intercept = -1.0,
  S = 2L
)
draws_occ_gb <- make_region_dropped_occ_bundle(
  frame_gb,
  intercept = -1.0,
  S = 2L
)

dur_band_nm <- c(
  "(Intercept)",
  "mso_band2",
  "mso_band3",
  "mso_band4M5",
  "mso_band6M8",
  "mso_band9P"
)
draws_dur_ga <- make_region_dropped_dur_bundle(
  frame_ga,
  setNames(c(-0.5, 0.8, -0.6, 0.9, -0.7, 0.5), dur_band_nm),
  S = 2L
)
draws_dur_gb <- make_region_dropped_dur_bundle(
  frame_gb,
  setNames(c(0.3, -0.9, 0.7, -0.4, 0.6, -0.8), dur_band_nm),
  S = 2L
)
draws_dur_ga$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)
draws_dur_gb$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)

# Distinct per-group mon coefficients.
mon_nm <- c("(Intercept)", "mso_cb1")
draws_mon_ga <- make_region_dropped_mon_bundle(
  frame_ga,
  setNames(c(log(20), -0.15), mon_nm),
  S = 2L
)
draws_mon_gb <- make_region_dropped_mon_bundle(
  frame_gb,
  setNames(c(log(8), 0.10), mon_nm),
  S = 2L
)

months4 <- as.Date("2022-01-01") %m+% months(0:3)
cb_meta_mon_ga <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, 9L),
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)
cb_meta_mon_gb <- cb_meta_mon_ga

grouping4 <- list(
  draws_occ_by_group = list(ga = draws_occ_ga, gb = draws_occ_gb),
  draws_dur_by_group = list(ga = draws_dur_ga, gb = draws_dur_gb),
  draws_mon_by_group = list(ga = draws_mon_ga, gb = draws_mon_gb),
  cb_meta_mon_by_group = list(ga = cb_meta_mon_ga, gb = cb_meta_mon_gb)
)

panel4 <- tidyr::expand_grid(gid = ft_grp$gid, t_start = months4) %>%
  dplyr::left_join(ft_grp %>% dplyr::select(gid, min_seas), by = "gid") %>%
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
window_grp <- list(
  start = as.character(months4[1]),
  end = as.character(months4[length(months4)] %m+% months(1L))
)

obj_c4 <- tryCatch(
  build_spatpomp_tu_c(
    ft_grp,
    panel4,
    draws_occ_ga,
    window_grp,
    fixed_draw_idx = 1L,
    grouping = grouping4,
    measurement = "occurrence_duration_mon_NB"
  ),
  error = function(e) {
    cat("  ERROR (C):", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(!is.null(obj_c4), "T4: grouped C-variant (mon) builds")

if (!is.null(obj_c4)) {
  expect_true(
    !is.null(attr(obj_c4, "grouping_mon_sub")),
    "T4: grouped C-variant grouping_mon_sub set"
  )

  blks_c4 <- make_block_list(obj_c4, ft_grp)
  bpf_c4 <- tryCatch(
    run_bpfilter_tu(obj_c4, Np = 100L, block_list = blks_c4, seed = 11L),
    error = function(e) {
      warning("T4 C bpfilter: ", conditionMessage(e))
      NULL
    }
  )
  ll_c4 <- if (!is.null(bpf_c4)) {
    tryCatch(pomp::logLik(bpf_c4), error = function(e) NA_real_)
  } else {
    NA_real_
  }
  expect_true(
    is.finite(ll_c4),
    sprintf("T4: C-variant (grouped) logLik finite (%.3f)", ll_c4 %||% NaN)
  )
} else {
  cat("  SKIP   T4 — object build failed\n")
}

# ===========================================================================
# T5 — Covariate-name collision fix: .sp_dedigit_covar_names()
# ===========================================================================
# spatPomp v1.1.0's construct_spatPomp() flattens the wide covariate table via
# paste0(covariate_name, unit_index), where unit_index is the bare positional
# integer 1..U (see .sp_dedigit_covar_names()'s docstring, spatpomp_core_c.R,
# for the full derivation). That concatenation collides whenever one
# covariate name is another name plus digits -- exactly what
# precompute_mon_mso_coefs() produces (mon_mso_coef_j0..j<mso_max>, deployed
# mso_max 12-26): "mon_mso_coef_j1"+11 == "mon_mso_coef_j11"+1 ==
# "mon_mso_coef_j111". Collision threshold is U >= 11 (no unit index below 11
# has a leading "1" long enough to alias the j1/j11 pair).
# T5: covariate-name collision fix (.sp_dedigit_covar_names(),
#     spatpomp_core_c.R): spatPomp's internal covariate flatten
#     (paste0(name, unit_index)) collides whenever one covariate name is
#     another name plus digits -- exactly what mon_mso_coef_j0..j<mso_max>
#     produces (deployed mso_max 12-26) once a TU has U >= 11 units. T5a is a
#     deterministic replay of spatPomp's own flatten logic on a minimal
#     covar_tbl at and below the U=11 threshold; T5b builds a real spatPomp
#     object at U=12 with mso_max=12, a config whose names collide without the
#     de-digiting.
cat("\n## T5: .sp_dedigit_covar_names() — covariate-name collision fix\n")

# T5a — direct, deterministic replay of spatPomp's own flatten logic
# (construct_spatPomp()'s covar branch) on a minimal covar_tbl, reproducing
# its exact error at U=11 and pinning that U=10 does NOT collide.
.replicate_spatpomp_covar_flatten <- function(covar_tbl, unit_covarnames) {
  unit_names <- unique(covar_tbl$unit)
  tmp <- seq_along(unit_names)
  names(tmp) <- unit_names
  pc <- dplyr::mutate(covar_tbl, ui = match(covar_tbl$unit, names(tmp)))
  pc <- tidyr::gather(
    pc,
    dplyr::all_of(unit_covarnames),
    key = "covname",
    value = "val"
  )
  pc <- dplyr::mutate(pc, covname = paste0(covname, ui))
  pc <- dplyr::select(pc, -dplyr::all_of("unit"), -ui)
  tidyr::spread(pc, key = "covname", value = "val")
}

make_collision_tbl <- function(U) {
  tidyr::expand_grid(time = 1:2, unit = sprintf("u%02d", seq_len(U))) %>%
    dplyr::mutate(mon_mso_coef_j1 = 1.0, mon_mso_coef_j11 = 2.0)
}

collision_cols <- c("mon_mso_coef_j1", "mon_mso_coef_j11")

err_u11 <- tryCatch(
  {
    .replicate_spatpomp_covar_flatten(make_collision_tbl(11L), collision_cols)
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect_true(
  !is.null(err_u11) && grepl("unique combination of keys", err_u11),
  "T5a: original mon_mso_coef_j1/j11 names collide at U=11 (reproduces spatPomp's own error)"
)

err_u10 <- tryCatch(
  {
    .replicate_spatpomp_covar_flatten(make_collision_tbl(10L), collision_cols)
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect_true(
  is.null(err_u10),
  "T5a: original names do NOT collide at U=10 (pins the U=11 threshold)"
)

tbl_u11_fixed <- .sp_dedigit_covar_names(make_collision_tbl(11L))
fixed_cols <- paste0(collision_cols, "_")
expect_true(
  all(fixed_cols %in% names(tbl_u11_fixed)) &&
    !anyDuplicated(names(tbl_u11_fixed)),
  "T5a: .sp_dedigit_covar_names() renames both collision columns, no duplicates"
)

flat_u11_fixed <- tryCatch(
  .replicate_spatpomp_covar_flatten(tbl_u11_fixed, fixed_cols),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(
  !is.null(flat_u11_fixed) &&
    !anyDuplicated(names(flat_u11_fixed)) &&
    nrow(flat_u11_fixed) == 2L,
  "T5a: de-digited names flatten without collision at U=11"
)

# T5b — integration test through the real builder: U=12 (3 countries x 2 adm1
# x 2 adm2), mso_max=12 (range upper bound 12), so mon_mso_coef_j1/j11/j12
# genuinely collide without .sp_dedigit_covar_names().
fix_t5 <- make_fixture(
  tu = "hmmT5",
  countries = c("AAA", "BBB", "CCC"),
  adm1_per_country = 2L,
  adm2_per_adm1 = 2L,
  n_months = 3L
)
U_t5 <- nrow(fix_t5$frame)
window_t5 <- list(
  start = as.character(fix_t5$months[1]),
  end = as.character(fix_t5$months[length(fix_t5$months)] %m+% months(1L))
)

draws_occ_t5 <- make_draws_bundle_synthetic(
  fix_t5$frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -1.0),
  basehaz_const = 0.0
)
draws_dur_t5 <- make_draws_bundle_synthetic(
  fix_t5$frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -0.8),
  active_terms = c("(Intercept)")
)

cb_meta_mon_t5 <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, 12L), # mso_max = 12 -> mon_mso_coef_j0..j12
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)
draws_mon_t5 <- make_draws_bundle_synthetic(
  fix_t5$frame,
  S = 2L,
  active_terms = c("(Intercept)", "mso_cb1")
)
draws_mon_t5$fe[, "(Intercept)"] <- log(20)
draws_mon_t5$fe[, "mso_cb1"] <- -0.15
draws_mon_t5$phi <- rep(8.0, 2L)
draws_mon_t5$offset_pop <- FALSE
draws_mon_t5$offset_logoffset <- FALSE

obj_t5 <- tryCatch(
  build_spatpomp_tu_c(
    fix_t5$frame,
    fix_t5$panel,
    draws_occ_t5,
    window_t5,
    fixed_draw_idx = 1L,
    draws_dur = draws_dur_t5,
    draws_mon = draws_mon_t5,
    cb_meta_mon = cb_meta_mon_t5,
    measurement = "occurrence_duration_mon_NB"
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(
  !is.null(obj_t5),
  sprintf(
    "T5b: U=%d, mso_max=12 (mon_mso_coef_j1/j11/j12 collision live) builds without error",
    U_t5
  )
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
