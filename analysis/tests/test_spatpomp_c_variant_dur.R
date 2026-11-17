# test_spatpomp_c_variant_dur.R — HARD: smoke / structural checks of the
#   duration-aware (cessation-hazard) spatPomp Csnippet kernel.
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_c_variant_dur.R
#
# Details:
#   Run from the project root. build_spatpomp_tu_c(..., draws_dur=): T1
#   compile / structure, T2 degenerate active / mso trajectories, T3-T4
#   bpfilter smoke (ungrouped and grouped mso-band hazards), T5 the duration
#   calendar harmonics. Only C-side assertions.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
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
  "# test_spatpomp_c_variant_dur — duration C kernel smoke/structural checks\n\n"
)

# ---------------------------------------------------------------------------
# Shared fixture (same helper as test_spatpomp_c_variant.R)
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

# Degenerate onset draws (occurrence): near-certain onset every month
# (mirrors test_spatpomp_c_variant.R's T3 fixture and test_duration_kernel.R's
# convention of intercept=+-15/30 for degenerate Bernoulli/hazard draws).
draws_occ_hi <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 15.0),
  basehaz_const = 0.0
)

# ===========================================================================
# T1 — Compile / smoke test
# ===========================================================================
# T1: build_spatpomp_tu_c(..., draws_dur=) produces a valid spatPomp object
#     with use_dur=TRUE and the new dur_* covariate columns present.
cat(
  "## T1: build_spatpomp_tu_c(..., draws_dur=) — compile and object structure\n"
)

draws_dur_never <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -30.0),
  active_terms = c("(Intercept)")
)

obj_c_dur <- tryCatch(
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
  !is.null(obj_c_dur),
  "C-variant object with draws_dur builds without error"
)

if (!is.null(obj_c_dur)) {
  expect_true(isTRUE(attr(obj_c_dur, "use_dur")), "use_dur attribute is TRUE")
  expect_true(
    !is.null(attr(obj_c_dur, "draws_dur_sub")),
    "draws_dur_sub attribute set"
  )

  U <- nrow(frame_tu)
  filt_months <- attr(obj_c_dur, "filt_months")
  ws1 <- attr(obj_c_dur, "warmup_state")
  draws_sub1 <- .subset_draws(draws_occ_hi, 1L)
  draws_dur_sub1 <- .subset_draws(draws_dur_never, 1L)

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

  aug_cov <- build_covar_for_tu_c(
    base_cov,
    lp_static,
    ttlo_c,
    spill_c,
    ws1,
    frame_tu,
    dur_lp_static_vec = dur_lp,
    dur_ttlo_coefs_df = dur_ttlo,
    dur_spill_coefs_df = dur_spill
  )

  dur_cols <- c(
    "dur_lp_static",
    paste0("dur_ttlo_coef_b", 0:6),
    "dur_spill_o_t",
    "dur_spill_o_f",
    "dur_spill_oct_t",
    "dur_spill_oct_f",
    "dur_spill_otu_t",
    "dur_spill_otu_f",
    "active_init",
    "mso_init"
  )
  all_present <- all(dur_cols %in% names(aug_cov))
  expect_true(
    all_present,
    sprintf(
      "all %d dur_* / active-mso covariate columns present",
      length(dur_cols)
    )
  )
}

# T1: backward compat — build_spatpomp_tu_c() WITHOUT draws_dur is unaffected.
obj_c_nodur <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_occ_hi,
    window,
    fixed_draw_idx = 1L
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(
  !is.null(obj_c_nodur),
  "C-variant object without draws_dur still builds"
)
expect_true(
  !isTRUE(attr(obj_c_nodur, "use_dur")),
  "use_dur attribute FALSE/NULL when draws_dur omitted"
)

# ===========================================================================
# T2 — Degenerate state trajectories (ungrouped, basehaz kind="const")
# ===========================================================================
# T2: degenerate state trajectories (ungrouped, basehaz kind="const") under
#     two forced scenarios:
#       A) onset certain, cessation impossible -> active stays 1, mso
#          increments every month, no re-onset (tests the onset gate +
#          update_lags() fed by `active` + mso-increment-while-active path).
#       B) onset certain, cessation certain -> active pulses to 1 and back to
#          0 within the same step every month (tests that the cessation hazard
#          reads the post-promotion mso, i.e. cessation-in-onset-month, and
#          that update_lags() still sees `active` correctly).
cat("\n## T2: degenerate active/mso trajectories — C-variant\n")

run_scenario <- function(draws_dur_case, label) {
  obj_c <- tryCatch(
    build_spatpomp_tu_c(
      frame_tu,
      panel_tu,
      draws_occ_hi,
      window,
      fixed_draw_idx = 1L,
      draws_dur = draws_dur_case
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

# Scenario A: onset certain, cessation impossible -> active=1 mso=1,2,3...
draws_dur_A <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -30.0),
  active_terms = c("(Intercept)")
)
res_A <- run_scenario(draws_dur_A, "A(never-cease)")
expect_true(!is.null(res_A), "T2 A: C-variant builds and simulates")

if (!is.null(res_A)) {
  U <- nrow(frame_tu)
  mat_active_c <- extract_state(res_A$sm_c, res_A$obj_c, "active")
  mat_mso_c <- extract_state(res_A$sm_c, res_A$obj_c, "mso")
  expect_true(
    all(mat_active_c == 1),
    "T2 A: C-variant active = 1 for all (unit, month)"
  )
  expect_equal(
    mat_mso_c[, 1],
    seq_len(n_months),
    "T2 A: C-variant mso increments 1,2,3,... (no re-onset, always active)"
  )
}

# Scenario B: onset certain, cessation certain (within onset month) -> active
# pulses to 1 then back to 0 every step; mso ends each step at 0.
draws_dur_B <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 30.0),
  active_terms = c("(Intercept)")
)
res_B <- run_scenario(draws_dur_B, "B(always-cease)")
expect_true(!is.null(res_B), "T2 B: C-variant builds and simulates")

if (!is.null(res_B)) {
  mat_yout_c <- extract_state(res_B$sm_c, res_B$obj_c, "Yout")
  mat_active_c <- extract_state(res_B$sm_c, res_B$obj_c, "active")
  mat_mso_c <- extract_state(res_B$sm_c, res_B$obj_c, "mso")
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
}

# ===========================================================================
# T3 — bpfilter smoke (ungrouped, basehaz kind="factor")
# ===========================================================================
# T3: bpfilter smoke (ungrouped, basehaz kind="factor"): finite logLik under
#     moderate (non-degenerate) hazards through the mso-band lookup table
#     (precompute_dur_basehaz_table / DBH[]) and the dur_ttlo/dur_spill
#     covariate wiring.
cat("\n## T3: bpfilter smoke, ungrouped factor mso-band hazard\n")

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

fix4 <- make_fixture(n_months = 4L)
frame_tu4 <- fix4$frame
panel_tu4 <- fix4$panel
window4 <- list(
  start = as.character(fix4$months[1]),
  end = as.character(fix4$months[length(fix4$months)] %m+% months(1L))
)

obj_c3 <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu4,
    panel_tu4,
    draws_occ_mod,
    window4,
    fixed_draw_idx = 2L,
    draws_dur = draws_dur_mod
  ),
  error = function(e) {
    cat("  ERROR (C):", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(!is.null(obj_c3), "T3: C-variant (factor mso-band) builds")

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
# T4 — bpfilter smoke (GROUPED, basehaz kind="factor")
# ===========================================================================
# T4: bpfilter smoke (GROUPED, basehaz kind="factor", distinct per-group band
#     coefficients): use_dur and finite logLik through dur_grp_id_sorted and
#     the per-group DBH[N_GRP][K] table.
cat("\n## T4: bpfilter smoke, GROUPED factor mso-band hazard\n")
# Exercises dur_grp_id_sorted (frame_sorted-order group topology) and the
# per-group DBH[N_GRP][K] table -- the part of the C duration path most at risk of a
# silent unit<->group mismatch bug (e.g. frame_tu- vs frame_sorted-order).

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
# Clearly distinct per-group band coefficients.
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

grouping4 <- list(
  draws_occ_by_group = list(ga = draws_occ_ga, gb = draws_occ_gb),
  draws_dur_by_group = list(ga = draws_dur_ga, gb = draws_dur_gb)
)

months4 <- as.Date("2022-01-01") %m+% months(0:3)
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
    grouping = grouping4
  ),
  error = function(e) {
    cat("  ERROR (C):", conditionMessage(e), "\n")
    NULL
  }
)
expect_true(!is.null(obj_c4), "T4: grouped C-variant builds")

if (!is.null(obj_c4)) {
  expect_true(
    isTRUE(attr(obj_c4, "use_dur")),
    "T4: grouped C-variant use_dur=TRUE"
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
# T5 — duration.seasonal_term harmonic2 in the static duration LP
# ===========================================================================
# T5: duration.seasonal_term harmonic2: a duration fit carrying the four
#     calendar harmonics gets them in its static LP (dur_lp_static shifts by
#     b . harmonics(calendar month) at every pomp time, month k at time k-1); a
#     fit without them leaves the covariate table and dur_covar_map as before
#     (no harmonic columns without the mon track).
cat("\n## T5: duration calendar harmonics fold into dur_lp_static\n")

fix12 <- make_fixture(n_months = 12L)
window12 <- list(
  start = as.character(fix12$months[1]),
  end = as.character(fix12$months[12] %m+% months(1L))
)
occ12 <- make_draws_bundle_synthetic(
  fix12$frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -1.0),
  basehaz_const = 0.0
)
harm5 <- c(
  sin_month = 0.7,
  cos_month = -0.4,
  sin2_month = 0.3,
  cos2_month = -0.2
)
dur_noh <- make_draws_bundle_synthetic(
  fix12$frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -1.5),
  active_terms = "(Intercept)"
)
dur_harm5 <- make_draws_bundle_synthetic(
  fix12$frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -1.5, harm5),
  active_terms = c("(Intercept)", names(harm5))
)
build12 <- function(d) {
  tryCatch(
    build_spatpomp_tu_c(
      fix12$frame,
      fix12$panel,
      occ12,
      window12,
      fixed_draw_idx = 1L,
      draws_dur = d
    ),
    error = function(e) {
      cat("  ERROR (T5):", conditionMessage(e), "\n")
      NULL
    }
  )
}
o_noh <- build12(dur_noh)
o_h <- build12(dur_harm5)
expect_true(
  !is.null(o_noh) && !is.null(o_h),
  "T5: both C-variant objects build"
)

if (!is.null(o_noh) && !is.null(o_h)) {
  tb_noh <- o_noh@covar@table
  tb_h <- o_h@covar@table
  expect_true(
    length(attr(o_noh, "dur_covar_map")) == 0L &&
      !any(grepl("^sin_month", rownames(tb_noh))),
    "T5: no harmonics -> empty dur_covar_map and no harmonic covariate columns"
  )
  expect_true(
    identical(unname(attr(o_h, "dur_covar_map")), names(harm5)),
    "T5: harmonics -> dur_covar_map maps the four terms to themselves"
  )
  lp_rows <- function(tb) {
    tb[grep("^dur_lp_static", rownames(tb)), , drop = FALSE]
  }
  d_lp <- lp_rows(tb_h) - lp_rows(tb_noh)
  m12 <- 1:12
  Xh <- cbind(
    sin(2 * pi * m12 / 12),
    cos(2 * pi * m12 / 12),
    sin(4 * pi * m12 / 12),
    cos(4 * pi * m12 / 12)
  )
  want <- drop(Xh %*% harm5) # calendar month k sits at pomp time k-1
  got <- d_lp[, seq_len(12L), drop = FALSE]
  expect_true(
    max(abs(sweep(got, 2, want))) < 1e-12,
    sprintf(
      "T5: dur_lp_static shift = b . harmonics(month) in every unit and month (max err %.2g)",
      max(abs(sweep(got, 2, want)))
    )
  )
  other <- setdiff(
    rownames(tb_noh),
    grep("^dur_lp_static", rownames(tb_noh), value = TRUE)
  )
  expect_true(
    identical(tb_h[other, , drop = FALSE], tb_noh[other, , drop = FALSE]),
    "T5: every other covariate column is unchanged"
  )
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
