# test_spatpomp_recovery_mon.R — simulation-recovery smoke check of the
#   monthly-cases (`mon`) measurement on the spatPomp Csnippet kernel.
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_recovery_mon.R
#
# Details:
#   Run from the project root; not registered in analysis/tests/run_all.R.
#   Companion to test_spatpomp_recovery.R: does mon_mu_block
#   (spatpomp_core_c.R, after duration promotion and before cessation) build,
#   simulate and filter without error, and does its `mu` state track the true
#   per-month NegBin mean once the filter conditions on observed occurrence +
#   monthly case counts? mon hard-requires a duration fit (see
#   generate_rprocess_csnippet()), so draws_dur always rides along.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
.fails <- 0L
.passes <- 0L

cat(
  "# test_spatpomp_recovery_mon — bpfilter recovers a known occurrence+duration+mon trajectory\n\n"
)

# ---------------------------------------------------------------------------
# Fixture (identical location/panel scaffold to test_spatpomp_recovery.R)
# ---------------------------------------------------------------------------
make_fixture <- function(
  tu = "hmmW",
  countries = c("AAA", "BBB"),
  adm1_per_country = 2L,
  adm2_per_adm1 = 1L,
  n_months = 24L,
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

#' Mean-across-particles (time x unit) matrix for a state variable, read from
#' a bpfilterd_spatPomp's bf@saved.states (list of (statevars x Np) matrices,
#' one per filter time step). Aligns on the LAST M entries so an optional
#' leading t0 slice (if bpfilter ever emits one) does not shift the index.
extract_filtered_mean <- function(bpf, obj, varname, M) {
  U <- length(obj@unit_names)
  saved <- bpf@saved.states
  offset <- max(0L, length(saved) - M)
  mat <- matrix(NA_real_, nrow = M, ncol = U)
  for (tt in seq_len(M)) {
    mt <- saved[[offset + tt]]
    for (i in seq_len(U)) {
      rn <- paste0(varname, i)
      if (rn %in% rownames(mt)) mat[tt, i] <- mean(mt[rn, ])
    }
  }
  mat
}

# ---------------------------------------------------------------------------
# Fixture + known parameters
# ---------------------------------------------------------------------------
fix <- make_fixture()
frame_tu <- fix$frame
panel_tu <- fix$panel
n_months <- length(fix$months)
window <- list(
  start = as.character(fix$months[1]),
  end = as.character(fix$months[n_months] %m+% months(1L))
)

draws_occ <- make_draws_bundle_synthetic(
  frame_tu,
  S = 1L,
  fe_overrides = c("(Intercept)" = -4.65),
  basehaz_const = 0.0
)

draws_dur <- make_draws_bundle_synthetic(
  frame_tu,
  S = 1L,
  active_terms = c(
    "(Intercept)",
    "mso_band2",
    "mso_band3",
    "mso_band4M5",
    "mso_band6M8",
    "mso_band9P"
  )
)
draws_dur$fe[, "(Intercept)"] <- -4.24
draws_dur$fe[, "mso_band2"] <- 0.4
draws_dur$fe[, "mso_band3"] <- -0.3
draws_dur$fe[, "mso_band4M5"] <- 0.6
draws_dur$fe[, "mso_band6M8"] <- -0.5
draws_dur$fe[, "mso_band9P"] <- 0.2
draws_dur$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)

# Monthly-size (`mon`) draws bundle: intercept + a linear months-since-onset
# slope. Offsets explicitly OFF (pop/log_offset) so mu is driven
# purely by intercept + mso, independent of the fixture's pop scale -- the
# offset wiring itself is covered by test_monthly_offset_wellformed.R, not
# this test.
draws_mon <- make_draws_bundle_synthetic(
  frame_tu,
  S = 1L,
  active_terms = c("(Intercept)", "mso_cb1")
)
draws_mon$fe[, "(Intercept)"] <- log(20) # ~20 cases/month at onset (mso=0)
draws_mon$fe[, "mso_cb1"] <- -0.15 # declining mean as the episode ages
draws_mon$phi <- 8.0 # NegBin dispersion
draws_mon$offset_pop <- FALSE
draws_mon$offset_logoffset <- FALSE

# The single linear column (rather than the shipped bs(degree=3) default) is
# sufficient to exercise build_mso_basis()/precompute_mon_mso_coefs()
# end-to-end without the spline boundary-knot machinery this test does not
# need to re-verify.
# mso cross-basis fixture: dlnm::onebasis(fun="lin") -> exactly one column
# ("mso_cb1"), matching draws_mon$fe above. range must cover the horizon the
# episodes can reach so the clamp in build_mso_basis()/precompute_mon_mso_coefs()
# is not immediately saturated.
cb_meta_mon <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, n_months),
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)

# ===========================================================================
# Step 1 — simulate ONE ground-truth realization (occurrence + duration + mon)
# ===========================================================================
# Design (mirrors test_spatpomp_recovery.R's N-independent-runs pattern):
#   1. Build a spatPomp object (measurement="occurrence_duration_mon_NB") under
#      known occurrence + duration + mon parameters and simulate() ONE
#      ground-truth realization -> true Yout/active/mso/mu trajectories.
#   2. Draw a monthly case count at every active month from
#      rnbinom(mu = true mu, size = size_phi) -- the mon Csnippet's own hurdle
#      shape -- to stand in for an "observed" outbreak-size panel.
#   3. Feed Yout + the drawn case counts back as obs_data (obs_Y, obs_CM via
#      mon_case_col) and rebuild under the state-conditioning sibling
#      "occurrence_duration_mon_NB_state" (near-delta on Yout, NegBin term on
#      obs_CM gated on active -- see .SP_DUNIT_MON_NB_STATE).
#   4. Run bpfilter (Np particles, save_states=TRUE); check:
#        - the object builds and bpfilter completes with a finite logLik
#          (the actual "does the C filter work" question),
#        - Yout recovery is near-exact (directly observed, as in the parent
#          recovery test),
#        - the filtered `mu` state is positive on active unit-months and
#          correlates with the true mu trajectory (sanity that mon_mu_block's
#          mso lookup table is actually being read, not zeroed/NaN).
cat("## Step 1: simulate ground-truth occurrence+duration+mon trajectory\n")

obj_truth <- tryCatch(
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_occ,
    window,
    fixed_draw_idx = 1L,
    draws_dur = draws_dur,
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
  !is.null(obj_truth),
  "ground-truth spatPomp object builds without error (mon path)"
)

sm_truth <- NULL
if (!is.null(obj_truth)) {
  set.seed(321L)
  sm_truth <- tryCatch(
    spatPomp::simulate(obj_truth, format = "data.frame"),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  expect_true(!is.null(sm_truth), "ground-truth simulate() succeeds (mon path)")
}

if (!is.null(obj_truth) && !is.null(sm_truth)) {
  mat_yout_truth <- extract_state(sm_truth, obj_truth, "Yout")
  mat_active_truth <- extract_state(sm_truth, obj_truth, "active")
  mat_mso_truth <- extract_state(sm_truth, obj_truth, "mso")
  mat_mu_truth <- extract_state(sm_truth, obj_truth, "mu")

  # `mat_active_truth` (the RECORDED active state) reflects status heading
  # INTO the next month, not "active during this month": dur_cessation_block
  # can reset active/mso to 0 in the same rprocess step that mon_mu_block
  # used their pre-cessation values to compute mu. So a genuinely active
  # FINAL month of an episode reads active==0 here. Reconstruct "active
  # DURING month t" independently -- onset this month, or already active
  # entering this month (one-month-lagged recorded active, cold-started at
  # row 1) -- mirroring smoothing_core.R's act_prev_mat derivation exactly.
  active_prev_truth <- rbind(
    rep(0, ncol(mat_active_truth)),
    mat_active_truth[-n_months, , drop = FALSE]
  )
  active_thismonth_truth <- (mat_yout_truth > 0.5) | (active_prev_truth > 0.5)

  n_onsets <- sum(mat_yout_truth)
  n_active <- sum(active_thismonth_truth)
  n_unitmonths <- length(active_thismonth_truth)
  max_mso <- suppressWarnings(max(mat_mso_truth, na.rm = TRUE))
  mu_active <- mat_mu_truth[active_thismonth_truth]

  cat(sprintf(
    "     onsets=%d  active unit-months=%d/%d  max mso=%d\n",
    n_onsets,
    n_active,
    n_unitmonths,
    max_mso
  ))
  cat(sprintf(
    "     mu on active unit-months: mean=%.2f range=[%.2f, %.2f]\n",
    mean(mu_active),
    min(mu_active),
    max(mu_active)
  ))

  expect_true(n_onsets >= 1L, "ground truth has at least one onset event")
  expect_true(
    is.finite(max_mso) && max_mso >= 2L,
    "ground truth: at least one episode persists beyond its onset month"
  )
  expect_true(
    length(mu_active) > 0L && all(is.finite(mu_active)) && all(mu_active > 0),
    "mu is finite and positive on every active unit-month (mon_mu_block wired correctly)"
  )
  expect_true(
    all(mat_mu_truth[!active_thismonth_truth] == 0),
    "mu is exactly 0 on inactive unit-months (mon_mu_block's av_new gate)"
  )
  expect_true(
    all((mat_mu_truth > 0) == active_thismonth_truth),
    "mu>0 exactly matches the independently-reconstructed active-this-month mask, INCLUDING an episode's final month (cessation-order fix)"
  )

  # ===========================================================================
  # Step 2 — draw monthly case counts from the mon hurdle's own NegBin shape,
  # standing in for an observed outbreak-size panel.
  # ===========================================================================
  cat("\n## Step 2: draw monthly case counts at every active unit-month\n")

  set.seed(654L)
  mat_cases_truth <- matrix(NA_real_, nrow = n_months, ncol = nrow(frame_tu))
  active_idx <- which(active_thismonth_truth)
  mat_cases_truth[active_idx] <- stats::rnbinom(
    length(active_idx),
    size = draws_mon$phi[1],
    mu = mat_mu_truth[active_idx]
  )

  cat(sprintf(
    "     drew %d monthly case observations (mean=%.2f)\n",
    length(active_idx),
    mean(mat_cases_truth[active_idx])
  ))

  # ===========================================================================
  # Step 3 — feed Yout + drawn cases back as obs_data; rebuild under the
  # state-conditioning ("_state") mon measurement.
  # ===========================================================================
  cat(
    "\n## Step 3: rebuild spatPomp object conditioned on the simulated observations\n"
  )

  unit_names_truth <- obj_truth@unit_names
  obs_data <- tibble::tibble(
    gid = rep(unit_names_truth, each = n_months),
    t_start = rep(fix$months, times = length(unit_names_truth)),
    obs_Y = as.integer(as.vector(mat_yout_truth)),
    cases_month = as.vector(mat_cases_truth)
  )

  obj_filt <- tryCatch(
    build_spatpomp_tu_c(
      frame_tu,
      panel_tu,
      draws_occ,
      window,
      fixed_draw_idx = 1L,
      draws_dur = draws_dur,
      draws_mon = draws_mon,
      cb_meta_mon = cb_meta_mon,
      measurement = "occurrence_duration_mon_NB_state",
      obs_data = obs_data,
      mon_case_col = "cases_month"
    ),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  expect_true(
    !is.null(obj_filt),
    "obs-conditioned spatPomp object builds without error (mon path)"
  )

  # ===========================================================================
  # Step 4 — run bpfilter (the actual "does the C filter work" question) and
  # read the filtered particle swarm.
  # ===========================================================================
  if (!is.null(obj_filt)) {
    cat("\n## Step 4: run bpfilter (Np=500) and extract filtered state means\n")

    blks <- make_block_list(obj_filt, frame_tu)
    Np <- 500L
    bpf <- tryCatch(
      run_bpfilter_tu(obj_filt, Np = Np, block_list = blks, seed = 123L),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        NULL
      }
    )
    expect_true(!is.null(bpf), "bpfilter runs without error (mon path)")

    if (!is.null(bpf)) {
      ll <- tryCatch(pomp::logLik(bpf), error = function(e) NA_real_)
      expect_true(
        is.finite(ll),
        sprintf("bpfilter logLik finite (%.3f)", ll %||% NaN)
      )

      # =========================================================================
      # Step 5 — compare filtered means against ground truth
      # =========================================================================
      cat("\n## Step 5: compare filtered trajectory against ground truth\n")

      mat_yout_filt <- extract_filtered_mean(bpf, obj_filt, "Yout", n_months)
      mat_mu_filt <- extract_filtered_mean(bpf, obj_filt, "mu", n_months)

      # Yout is directly observed (near-delta measurement) -> recovery should
      # be essentially exact, same bar as the parent recovery test.
      yout_acc <- mean(round(mat_yout_filt) == mat_yout_truth)
      cat(sprintf(
        "     Yout recovery accuracy (directly observed)  = %.3f\n",
        yout_acc
      ))
      expect_true(
        yout_acc >= 0.95,
        sprintf(
          "filtered Yout matches ground truth onset months (%.3f >= 0.95)",
          yout_acc
        )
      )

      # mu is latent-driven (mso lookup) and only indirectly informed by the
      # NegBin case term. Like the parent test's active/mso bar (active is
      # "never directly signalled" there either -- test_spatpomp_recovery.R,
      # majority-baseline + correlation), the right standard here is "not
      # degenerate", not exact recovery: a bootstrap/block filter at Np=500
      # can legitimately lose track of a long, weakly-observed episode for a
      # stretch of months (a whole block's particles agree cessation already
      # happened when truth says otherwise) without that being a wiring bug --
      # confirmed by inspection: the zero-filt cells cluster into a single
      # long run for one unit (its longest episode, mso up to 11), not
      # scattered across every episode's tail, which is what the
      # cessation-order bug (now fixed, see Step 1's assertions above) would
      # have produced.
      active_mask <- active_thismonth_truth
      mu_filt_act <- mat_mu_filt[active_mask]
      mu_true_act <- mat_mu_truth[active_mask]
      mu_cor <- suppressWarnings(cor(mu_filt_act, mu_true_act))
      frac_pos <- mean(mu_filt_act > 0)

      cat(sprintf(
        "     mu recovery: mean(filt)=%.2f mean(true)=%.2f correlation=%.3f\n",
        mean(mu_filt_act),
        mean(mu_true_act),
        mu_cor
      ))
      cat(sprintf(
        "     filtered mu > 0 on %.3f of active unit-months\n",
        frac_pos
      ))

      expect_true(
        all(is.finite(mu_filt_act)),
        "filtered mu is finite on every active unit-month"
      )
      expect_true(
        frac_pos >= 0.75,
        sprintf(
          "filtered mu is positive on most active unit-months (%.3f >= 0.75)",
          frac_pos
        )
      )
      expect_true(
        is.finite(mu_cor) && mu_cor > 0.3,
        sprintf(
          "filtered mu is positively associated with ground-truth mu (cor=%.3f > 0.3)",
          mu_cor
        )
      )
    }
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
