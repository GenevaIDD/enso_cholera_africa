# test_spatpomp_recovery.R — simulation-recovery check for the duration-aware
# (occurrence + cessation-hazard) spatPomp Csnippet kernel: does bpfilter
# recover the TRUE latent trajectory of a small synthetic outbreak process
# from its own simulated observations?
#
# Design (N-independent-runs, C-variant kernel):
#   1. Build one spatPomp object under known, moderate occurrence + duration
#      hazard parameters (fixed_draw_idx = 1L) and simulate() ONE ground-truth
#      realization -> true Yout (onset) / active / mso trajectories.
#   2. Feed the realized Yout back in as `obs_data` (gid, t_start, obs_Y) and
#      rebuild the spatPomp object under measurement="occurrence_bernoulli_state"
#      (the near-delta measurement that conditions on the LATENT Yout state,
#      not just the deterministic hazard `p` -- see .SP_STATE_MEAS in
#      spatpomp_core.R) with the SAME occurrence/duration parameters.
#   3. Run bpfilter (Np particles, save_states=TRUE) and read the filtered
#      particle swarm's mean state at every (time, unit) from
#      bpf@saved.states.
#   4. Compare the filtered means against the KNOWN ground truth:
#        - Yout is directly observed (near-delta) -> recovery should be ~exact.
#        - active/mso are NOT observed directly (cessation is never signalled
#          in the data -- only re-onset is) -> recovery is a genuine inference
#          problem. The bar is "clearly better than the always-predict-
#          majority-class baseline" and "positive correlation with the true
#          trajectory", not near-exact recovery.
#
# NOTE on log_offset: both the occurrence and duration hazards add
# log_offset = log(days_in_month) (~3.4) to the linear predictor (see
# spatpomp_core.R:417,715,829). The T3/T4 parity-test fixtures in
# test_spatpomp_c_variant_dur.R use Intercept=-1.0/-0.8, which — once this
# offset is included — are near-certain (~100%) hazards; fine for a bpfilter
# smoke check, useless for a recovery test
# (no variability to recover). This file's Intercepts are chosen to
# counteract the offset and land near a realistic ~25-35%/month hazard.
#
# Run: Rscript analysis/tests/test_spatpomp_recovery.R
#      (must be called from the project root)

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core_c.R")   # also sources spatpomp_core.R

source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal
.fails <- 0L; .passes <- 0L

cat("# test_spatpomp_recovery — bpfilter recovers a known occurrence+duration trajectory\n\n")

# ---------------------------------------------------------------------------
# Fixture (same helper convention as test_spatpomp_c_variant_dur.R)
# ---------------------------------------------------------------------------
make_fixture <- function(tu = "hmmW",
                          countries = c("AAA", "BBB"),
                          adm1_per_country = 2L,
                          adm2_per_adm1 = 1L,
                          n_months = 24L,
                          base_date = as.Date("2022-01-01")) {
  frame <- expand.grid(
    country = countries,
    adm1n   = seq_len(adm1_per_country),
    adm2n   = seq_len(adm2_per_adm1),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      tu          = tu,
      gid         = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit   = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", tu),
      min_seas    = 1L,
      pop         = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop,
                  adm1_unit) %>%
    tibble::as_tibble()

  months <- base_date %m+% months(0:(n_months - 1L))

  panel <- tidyr::expand_grid(gid = frame$gid, t_start = months) %>%
    dplyr::left_join(frame %>% dplyr::select(gid, min_seas), by = "gid") %>%
    dplyr::mutate(
      country           = sub("_[0-9]+_[0-9]+$", "", gid),
      t_start_seas      = as.integer(lubridate::yday(t_start) - 1L),
      mean_precip       = 0.0,
      dev_precip        = 0.0,
      mean_temp_daymean = 0.0,
      dev_temp_daymean  = 0.0,
      outbreak          = 0L,
      time_to_last_outbreak_cat     = "[6,Inf]",
      ongoing_outbreak_lag1         = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1      = 0L
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
    if (nrow(row_t) == 0) next
    if ("unitname" %in% names(row_t)) {
      for (i in seq_len(U)) {
        r_ui <- row_t[row_t$unitname == obj@unit_names[i], , drop = FALSE]
        if (nrow(r_ui) > 0 && varname %in% names(r_ui))
          mat[tt, i] <- r_ui[[varname]][1L]
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
  U      <- length(obj@unit_names)
  saved  <- bpf@saved.states
  offset <- max(0L, length(saved) - M)
  mat    <- matrix(NA_real_, nrow = M, ncol = U)
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
# Fixture + known parameters (moderate, non-degenerate occurrence/duration
# hazards -- realistic ~25-35%/month once log_offset is accounted for)
# ---------------------------------------------------------------------------
fix      <- make_fixture()
frame_tu <- fix$frame
panel_tu <- fix$panel
n_months <- length(fix$months)
window   <- list(start = as.character(fix$months[1]),
                 end   = as.character(fix$months[n_months] %m+% months(1L)))

draws_occ <- make_draws_bundle_synthetic(frame_tu, S = 1L,
               fe_overrides = c("(Intercept)" = -4.65), basehaz_const = 0.0)

draws_dur <- make_draws_bundle_synthetic(
  frame_tu, S = 1L,
  active_terms = c("(Intercept)", "mso_band2", "mso_band3",
                    "mso_band4M5", "mso_band6M8", "mso_band9P"))
draws_dur$fe[, "(Intercept)"] <- -4.24
draws_dur$fe[, "mso_band2"]   <-  0.4
draws_dur$fe[, "mso_band3"]   <- -0.3
draws_dur$fe[, "mso_band4M5"] <-  0.6
draws_dur$fe[, "mso_band6M8"] <- -0.5
draws_dur$fe[, "mso_band9P"]  <-  0.2
draws_dur$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)

# ===========================================================================
# Step 1 — simulate ONE ground-truth realization
# ===========================================================================
cat("## Step 1: simulate ground-truth occurrence+duration trajectory\n")

obj_truth <- tryCatch(
  build_spatpomp_tu_c(frame_tu, panel_tu, draws_occ, window,
                      fixed_draw_idx = 1L, draws_dur = draws_dur),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
expect_true(!is.null(obj_truth), "ground-truth spatPomp object builds without error")

if (!is.null(obj_truth)) {
  set.seed(321L)
  sm_truth <- tryCatch(spatPomp::simulate(obj_truth, format = "data.frame"),
                        error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  expect_true(!is.null(sm_truth), "ground-truth simulate() succeeds")
}

if (!is.null(obj_truth) && !is.null(sm_truth)) {
  mat_yout_truth   <- extract_state(sm_truth, obj_truth, "Yout")
  mat_active_truth <- extract_state(sm_truth, obj_truth, "active")
  mat_mso_truth    <- extract_state(sm_truth, obj_truth, "mso")

  n_onsets     <- sum(mat_yout_truth)
  n_active     <- sum(mat_active_truth)
  n_unitmonths <- length(mat_active_truth)
  max_mso      <- suppressWarnings(max(mat_mso_truth, na.rm = TRUE))

  cat(sprintf("     onsets=%d  active unit-months=%d/%d  max mso=%d\n",
              n_onsets, n_active, n_unitmonths, max_mso))

  expect_true(n_onsets >= 1L, "ground truth has at least one onset event")
  expect_true(is.finite(max_mso) && max_mso >= 2L,
              "ground truth: at least one episode persists beyond its onset month")

  # ===========================================================================
  # Step 2 — feed the realized Yout back as obs_data; rebuild under the
  # state-conditioning ("_state") measurement so bpfilter is actually pulled
  # toward the ground-truth trajectory.
  # ===========================================================================
  cat("\n## Step 2: rebuild spatPomp object conditioned on the simulated observations\n")

  unit_names_truth <- obj_truth@unit_names
  obs_data <- tibble::tibble(
    gid     = rep(unit_names_truth, each = n_months),
    t_start = rep(fix$months, times = length(unit_names_truth)),
    obs_Y   = as.integer(as.vector(mat_yout_truth))
  )

  obj_filt <- tryCatch(
    build_spatpomp_tu_c(frame_tu, panel_tu, draws_occ, window,
                        fixed_draw_idx = 1L, draws_dur = draws_dur,
                        measurement = "occurrence_bernoulli_state",
                        obs_data = obs_data),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  expect_true(!is.null(obj_filt), "obs-conditioned spatPomp object builds without error")

  # ===========================================================================
  # Step 3 — run bpfilter and read the filtered particle swarm
  # ===========================================================================
  if (!is.null(obj_filt)) {
    cat("\n## Step 3: run bpfilter (Np=500) and extract filtered state means\n")

    blks <- make_block_list(obj_filt, frame_tu)
    Np   <- 500L
    bpf  <- tryCatch(run_bpfilter_tu(obj_filt, Np = Np, block_list = blks, seed = 123L),
                      error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
    expect_true(!is.null(bpf), "bpfilter runs without error")

    if (!is.null(bpf)) {
      ll <- tryCatch(pomp::logLik(bpf), error = function(e) NA_real_)
      expect_true(is.finite(ll), sprintf("bpfilter logLik finite (%.3f)", ll %||% NaN))

      # =========================================================================
      # Step 4 — compare filtered means against ground truth
      # =========================================================================
      cat("\n## Step 4: compare filtered trajectory against ground truth\n")

      mat_yout_filt   <- extract_filtered_mean(bpf, obj_filt, "Yout",   n_months)
      mat_active_filt <- extract_filtered_mean(bpf, obj_filt, "active", n_months)

      # Yout is directly observed (near-delta measurement) -> recovery should
      # be essentially exact.
      yout_acc <- mean(round(mat_yout_filt) == mat_yout_truth)
      cat(sprintf("     Yout recovery accuracy (directly observed)  = %.3f\n", yout_acc))
      expect_true(yout_acc >= 0.95,
                  sprintf("filtered Yout matches ground truth onset months (%.3f >= 0.95)", yout_acc))

      # active/mso are latent (cessation is never directly signalled in the
      # data) -- the bar is "clearly better than always predicting the
      # majority class" and "positive association with ground truth", not
      # near-exact recovery.
      active_acc      <- mean(round(mat_active_filt) == mat_active_truth)
      active_frac      <- mean(mat_active_truth)
      majority_baseline <- max(active_frac, 1 - active_frac)
      active_cor       <- suppressWarnings(cor(as.vector(mat_active_filt), as.vector(mat_active_truth)))

      cat(sprintf("     active recovery accuracy = %.3f  (majority-class baseline = %.3f)\n",
                  active_acc, majority_baseline))
      cat(sprintf("     active recovery correlation (filtered mean vs. true active) = %.3f\n",
                  active_cor))

      expect_true(active_acc > majority_baseline + 0.05,
                  sprintf("filtered active-state recovery (%.3f) beats majority-class baseline (%.3f) by >0.05",
                          active_acc, majority_baseline))
      expect_true(is.finite(active_cor) && active_cor > 0.3,
                  sprintf("filtered active-state mean is positively associated with ground truth (cor=%.3f > 0.3)",
                          active_cor))
    }
  }
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", .passes, .fails))
if (.fails > 0L) quit(status = 1L)
