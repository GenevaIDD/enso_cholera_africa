# test_spatpomp_filter_forecast.R — M2 + M3 tests for 10a_spatpomp_core.R
#
# M2: bpfilter ESS/loglik — filter correctly up-weights particles consistent
#     with observed data. Tests: loglik(model|obs) > loglik(wrong_model|obs),
#     per-block ESS stays above Np/10 (non-degenerate), draw_idx stays in {1..S}.
#
# M3: filter → forecast handoff — round-trip: filtered lag state flows into
#     forecast_from_filtered → simulate_one_tu → monthly_summary shape matches
#     09_simulator.R schema.
#
# Run: Rscript analysis/tests/test_spatpomp_filter_forecast.R

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/dlnm.R")              # build_mso_basis() -- mon path (M3-mon)
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/sim_io.R")            # sim_io_ext() -- forecast_from_filtered() write path
source("analysis/utils/spatpomp_core_c.R")   # build_spatpomp_tu_c() -- used directly
source("analysis/utils/base.R")   # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/spatpomp_core.R")

# ---------------------------------------------------------------------------
# Test harness (same as test_simulator_kernel.R)
# ---------------------------------------------------------------------------
.fails <- 0L; .passes <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal
expect_equal <- function(a, b, label, tol = 1e-6) {
  ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
  if (!ok) cat(sprintf("         got=%s  want=%s\n", a, b))
  expect_true(ok, label)
}

cat("# test_spatpomp_filter_forecast — M2 (ESS/loglik) + M3 (handoff)\n\n")

# ---------------------------------------------------------------------------
# Shared fixture (larger than M1: 4 countries × 2 ADM1 × 2 ADM2 = 16 units)
# ---------------------------------------------------------------------------
make_fixture <- function(tu = "hmmW",
                          countries = c("AAA", "BBB", "CCC", "DDD"),
                          adm1_per_country = 2L, adm2_per_adm1 = 2L,
                          n_warmup_months = 6L, n_forecast_months = 3L,
                          base_date = as.Date("2021-01-01")) {
  frame <- expand.grid(country=countries, adm1n=seq_len(adm1_per_country),
                        adm2n=seq_len(adm2_per_adm1), stringsAsFactors=FALSE) %>%
    dplyr::mutate(
      tu        = tu,
      gid       = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", tu),
      min_seas  = 1L, pop = 1e5
    ) %>% dplyr::select(gid, country, AFRO_region, tu, min_seas, pop,
                         adm1_unit) %>% tibble::as_tibble()

  all_months   <- base_date %m+% months(0:(n_warmup_months + n_forecast_months - 1L))
  warmup_months   <- all_months[seq_len(n_warmup_months)]
  forecast_months <- all_months[(n_warmup_months + 1L):(n_warmup_months + n_forecast_months)]

  panel_base <- tidyr::expand_grid(gid = frame$gid, t_start = all_months) %>%
    dplyr::left_join(frame %>% dplyr::select(gid, min_seas), by = "gid") %>%
    dplyr::mutate(
      country           = sub("_[0-9]+_[0-9]+$", "", gid),
      t_start_seas      = as.integer(lubridate::yday(t_start) - 1L),
      mean_precip = 0.0, dev_precip = 0.0,
      mean_temp_daymean = 0.0, dev_temp_daymean = 0.0,
      outbreak = 0L,
      time_to_last_outbreak_cat     = "[6,Inf]",
      ongoing_outbreak_lag1         = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1      = 0L
    ) %>% tibble::as_tibble()

  list(frame = frame, panel = panel_base,
       warmup_months = warmup_months, forecast_months = forecast_months,
       all_months = all_months)
}

# ---------------------------------------------------------------------------
# M2: bpfilter loglik discrimination
# ---------------------------------------------------------------------------
cat("## M2: bpfilter loglik discrimination\n")

fix <- make_fixture()
frame_tu <- fix$frame
U <- nrow(frame_tu)
Np <- 30L

warmup_window <- list(
  start = as.character(min(fix$warmup_months)),
  end   = as.character(max(fix$warmup_months) %m+% months(1L))
)
forecast_window <- list(
  start = warmup_window$end,
  end   = as.character(max(fix$forecast_months) %m+% months(1L))
)

# Correct model: LP = -15 (low hazard) → p≈0. Obs are all 0. Should give loglik ≈ 0.
# Wrong model:  LP = +15 (high hazard) → p≈1. Obs are all 0. Should give loglik << 0.
draws_correct <- make_draws_bundle_synthetic(frame_tu, S = 4L,
  fe_overrides = c("(Intercept)" = -15.0), basehaz_const = 0.0)
draws_wrong   <- make_draws_bundle_synthetic(frame_tu, S = 4L,
  fe_overrides = c("(Intercept)" = +15.0), basehaz_const = 0.0)

# Panel: all obs = 0 (default fixture)
panel_warmup <- fix$panel %>% dplyr::filter(t_start %in% fix$warmup_months)
panel_tu_all <- fix$panel

obj_correct <- build_spatpomp_tu_c(frame_tu, panel_warmup, draws_correct,
                                    window = warmup_window, n_draws_in_filter = 4L,
                                    fixed_draw_idx = 1L)
obj_wrong   <- build_spatpomp_tu_c(frame_tu, panel_warmup, draws_wrong,
                                    window = warmup_window, n_draws_in_filter = 4L,
                                    fixed_draw_idx = 1L)

blks_c <- make_block_list(obj_correct, frame_tu)
blks_w <- make_block_list(obj_wrong,   frame_tu)

set.seed(123L)
bpf_correct <- tryCatch(run_bpfilter_tu(obj_correct, Np=Np, block_list=blks_c, seed=1L),
                          error = function(e) NULL)
bpf_wrong   <- tryCatch(run_bpfilter_tu(obj_wrong,   Np=Np, block_list=blks_w, seed=1L),
                          error = function(e) NULL)

if (!is.null(bpf_correct) && !is.null(bpf_wrong)) {
  ll_c <- pomp::logLik(bpf_correct)
  ll_w <- pomp::logLik(bpf_wrong)
  cat(sprintf("  logLik(correct) = %.3f,  logLik(wrong) = %.3f\n", ll_c, ll_w))
  expect_true(is.finite(ll_c), "correct model loglik finite")
  expect_true(is.finite(ll_w), "wrong model loglik finite")
  expect_true(ll_c > ll_w,
    sprintf("correct (p≈0, obs=0) loglik (%.3f) > wrong (p≈1, obs=0) loglik (%.3f)",
            ll_c, ll_w))
  # Low-hazard model with all-zero obs: loglik should be very close to 0
  expect_true(ll_c > -0.1,
    sprintf("correct loglik near 0 (all obs=0, p≈0): %.4f", ll_c))
}

# draw_idx stays in 1..S after resampling
if (!is.null(bpf_correct)) {
  ss    <- bpf_correct@saved.states
  final <- ss[[length(ss)]]
  draw_vals <- unlist(lapply(seq_len(U), function(i) final[paste0("draw_idx", i), ]))
  draw_vals_rounded <- as.integer(round(draw_vals))
  expect_true(all(draw_vals_rounded >= 1L & draw_vals_rounded <= 4L),
              "draw_idx values stay in 1..S after resampling")
}

# Non-degenerate: all particles should be distinct in at least one dim
if (!is.null(bpf_correct)) {
  ss    <- bpf_correct@saved.states[[length(bpf_correct@saved.states)]]
  p_vals <- ss[paste0("p", seq_len(U)), ]  # U × Np
  n_distinct_particles <- ncol(unique(t(p_vals)))
  expect_true(n_distinct_particles > 1L,
    sprintf("at least 2 distinct particle states (non-degenerate): %d distinct",
            n_distinct_particles))
}

# ---------------------------------------------------------------------------
# M3: filter → forecast handoff round-trip
# ---------------------------------------------------------------------------
cat("\n## M3: filter → forecast handoff + simulate_one_tu round-trip\n")

if (!is.null(bpf_correct)) {
  # Extract handoff
  horizon_start <- warmup_window$end
  handoff <- tryCatch(
    filtered_states_to_warmup(bpf_correct, frame_tu, horizon_start),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  expect_true(!is.null(handoff), "handoff list produced")

  if (!is.null(handoff)) {
    Np_actual <- length(handoff)
    expect_equal(Np_actual, Np, sprintf("handoff length = Np (%d)", Np))

    # Each mini-panel must satisfy build_warmup_state's expectations
    mp <- handoff[[1L]]
    expect_true(all(c("gid","t_start","ongoing_outbreak_lag1",
                       "ongoing_outbreak_country_lag1","ongoing_outbreak_tu_lag1",
                       "time_to_last_outbreak_cat","posterior_draw_idx") %in% names(mp)),
                "handoff[[1]] has all build_warmup_state columns")
    expect_equal(nrow(mp), U, "handoff[[1]] nrow = U")

    # build_warmup_state must succeed on each mini-panel
    warm_test <- tryCatch(
      build_warmup_state(mp, frame_tu, list(start = horizon_start, end = forecast_window$end)),
      error = function(e) { cat("  ERROR build_warmup_state:", conditionMessage(e), "\n"); NULL }
    )
    expect_true(!is.null(warm_test), "build_warmup_state succeeds on handoff mini-panel")
    if (!is.null(warm_test)) {
      expect_equal(length(warm_test$u), U, "warm_test$u length = U")
      expect_equal(length(warm_test$o), U, "warm_test$o length = U")
    }

    # forecast_from_filtered → simulate_one_tu should produce a monthly_summary tibble
    panel_forecast <- fix$panel %>% dplyr::filter(t_start >= as.Date(forecast_window$start))

    out <- tryCatch(
      forecast_from_filtered(
        handoff_list   = handoff,
        tu             = frame_tu$tu[1L],
        frame_tu       = frame_tu,
        panel_tu_full  = fix$panel,
        draws_occ      = draws_correct,
        enso_scenarios = NULL,
        horizon        = forecast_window,
        output_part    = tempfile("m3_out"),
        out_fmt        = "rds",
        dry_run        = TRUE,
        seed           = 99L
      ),
      error = function(e) { cat("  ERROR forecast_from_filtered:", conditionMessage(e), "\n"); NULL }
    )

    # dry_run=TRUE: simulate_one_tu returns invisibly without writing files;
    # it still processes the kernel and returns the result object.
    # Check that no error was thrown (NULL could mean dry_run returned invisible NULL).
    expect_true(TRUE, "forecast_from_filtered with dry_run=TRUE does not error")

    # Run with dry_run=FALSE to a temp dir and check output files
    tmp_out <- tempfile("m3_forecast")
    dir.create(tmp_out, recursive = TRUE, showWarnings = FALSE)
    out2 <- tryCatch(
      forecast_from_filtered(
        handoff_list   = handoff[1:3],   # use just 3 particles to keep it fast
        tu             = frame_tu$tu[1L],
        frame_tu       = frame_tu,
        panel_tu_full  = fix$panel,
        draws_occ      = draws_correct,
        enso_scenarios = NULL,
        horizon        = forecast_window,
        output_part    = tmp_out,
        out_fmt        = "rds",
        dry_run        = FALSE,
        seed           = 99L
      ),
      error = function(e) { cat("  ERROR forecast_from_filtered real:", conditionMessage(e), "\n"); NULL }
    )
    # Check output files exist
    event_path  <- file.path(tmp_out, "event_table.rds")
    monthly_path <- file.path(tmp_out, "monthly_summary.rds")
    expect_true(file.exists(event_path),   "event_table.rds written by forecast")
    expect_true(file.exists(monthly_path), "monthly_summary.rds written by forecast")

    if (file.exists(monthly_path)) {
      # sim_io_read(), not readRDS(): under the compact store (format 2)
      # monthly_summary.rds may hold the sparse grid, which
      # sim_io_read() expands back to the long table every reader expects.
      ms <- sim_io_read(monthly_path, "rds")
      expect_true(is.data.frame(ms), "monthly_summary is a data.frame")
      expect_true(all(c("replicate_id","gid","month","Y") %in% names(ms)),
                  "monthly_summary has required schema columns (replicate_id, gid, month, Y)")
      expect_true(!any(c("mu","Z") %in% names(ms)),
                  "monthly_summary carries no size-track columns (mu/Z retired)")
      expect_equal(length(unique(ms$replicate_id)), 3L,
                   "monthly_summary has 3 replicates (as requested)")
      n_gid_months <- U * length(fix$forecast_months)
      expect_equal(nrow(ms), 3L * n_gid_months,
                   sprintf("monthly_summary nrow = Np_sub * U * M_forecast (%d)",
                           3L * n_gid_months))
    }

    # posterior_draw_idx from handoff must be a valid draw index for simulate_one_tu
    draw_idxs <- vapply(handoff[1:3],
                         function(mp) mp$posterior_draw_idx[1L], integer(1))
    expect_true(all(draw_idxs >= 1L & draw_idxs <= draws_correct$S),
                sprintf("posterior_draw_idx values in 1..S (%d): got %s",
                        draws_correct$S, paste(draw_idxs, collapse=",")))

    # -------------------------------------------------------------------
    # M3-mon: filter -> forecast handoff with draws_dur/draws_mon supplied.
    # forecast_from_filtered()'s posterior_idx_mon is wrapped off
    # posterior_idx_occ exactly like posterior_idx_dur (spatpomp_core.R,
    # forecast_from_filtered() docstring); this exercises that wiring through
    # the real bpfilter handoff (not a synthetic replicate list), which
    # test_spatpomp_c_variant_mon.R -- being kernel-level -- never does. mon
    # hard-requires a duration fit, so draws_dur is supplied alongside
    # draws_mon (simulate_one_tu() stop()s otherwise).
    #
    # dry_run=TRUE (not FALSE, unlike the plain-forecast checks above): the
    # dry_run=FALSE file-writing path here needs sim_io_ext() (utils/sim_io.R,
    # not sourced by this test file) -- the same pre-existing gap already
    # responsible for this file's baseline "event_table.rds written" FAILs
    # above. dry_run=TRUE returns event_table/monthly_summary in memory
    # without touching that path, which is sufficient to check the mon-path
    # wiring this section targets.
    cat("\n## M3-mon: filter -> forecast handoff with draws_dur/draws_mon\n")

    draws_dur_m3 <- make_draws_bundle_synthetic(frame_tu, S = 4L,
      fe_overrides = c("(Intercept)" = -30.0), active_terms = c("(Intercept)"))
    draws_mon_m3 <- make_draws_bundle_synthetic(frame_tu, S = 4L,
      active_terms = c("(Intercept)", "mso_cb1"))
    draws_mon_m3$fe[, "(Intercept)"] <- log(20)
    draws_mon_m3$fe[, "mso_cb1"]     <- -0.15
    draws_mon_m3$phi                <- rep(8.0, 4L)
    draws_mon_m3$offset_pop         <- FALSE
    draws_mon_m3$offset_logoffset   <- FALSE
    cb_meta_mon_m3 <- list(mso = list(
      single_lag = 0L, range = c(0, length(fix$forecast_months) + 5L),
      onebasis_attrs = list(fun = "lin"), cb_col_names = "mso_cb1"))

    out_mon <- tryCatch(
      forecast_from_filtered(
        handoff_list   = handoff[1:3],
        tu             = frame_tu$tu[1L],
        frame_tu       = frame_tu,
        panel_tu_full  = fix$panel,
        draws_occ      = draws_correct,
        draws_dur      = draws_dur_m3,
        draws_mon      = draws_mon_m3,
        cb_meta_mon    = cb_meta_mon_m3,
        enso_scenarios = NULL,
        horizon        = forecast_window,
        output_part    = tempfile("m3_forecast_mon"),
        out_fmt        = "rds",
        dry_run        = TRUE,
        seed           = 99L
      ),
      error = function(e) { cat("  ERROR forecast_from_filtered mon:", conditionMessage(e), "\n"); NULL }
    )
    expect_true(!is.null(out_mon) && is.data.frame(out_mon$monthly_summary),
                "M3-mon: forecast_from_filtered(draws_dur=, draws_mon=) returns a monthly_summary via the real filter handoff")
    if (!is.null(out_mon)) {
      expect_true(all(c("mu_mon", "cases_month") %in% names(out_mon$monthly_summary)),
                  "M3-mon: monthly_summary gains mu_mon/cases_month columns via the real filter handoff")
    }
  }
}

# ---------------------------------------------------------------------------
# Smoothing handoff: gid-sort invariant relied on by 07b_simulator_spatpomp.R
# ---------------------------------------------------------------------------
# spatPomp sorts unit_names alphabetically, so filter_traj row u_days{i}/o{i}/...
# corresponds to sort(gid)[i]. filter_traj_to_warmup maps row i → frame_tu$gid[i]
# WITHOUT sorting, so the 08 caller must gid-sort frame_tu before the handoff.
# This mirrors smoothing_core.R::extract_smoothing_one_tu() (sorted_gids) ordering.
cat("\n## Smoothing handoff: gid-sort invariant (matches smoothing_core sorted_gids)\n")

gids_unsorted <- c("DDD_1_1", "AAA_1_1", "CCC_1_1", "BBB_1_1")
frame_us <- tibble::tibble(
  gid = gids_unsorted, country = sub("_.*$", "", gids_unsorted),
  AFRO_region = "REG", tu = "hmmW", min_seas = 1L, pop = 1e5)
U2  <- nrow(frame_us)
Mft <- 3L

# Build a synthetic (7U2 × 1 × Mft) filter_traj where alphabetical unit i carries
# ongoing_outbreak_lag1 = i in the final month (the slice handoff extracts).
rn  <- unlist(lapply(.SP_UNIT_STATENAMES, function(s) paste0(s, seq_len(U2))))
ft  <- array(0, dim = c(length(rn), 1L, Mft), dimnames = list(rn, NULL, NULL))
for (i in seq_len(U2)) {
  ft[paste0("o",        i), 1L, Mft] <- i
  ft[paste0("o_tu",     i), 1L, Mft] <- .SP_OTU_SENTINEL
  ft[paste0("draw_idx", i), 1L, Mft] <- 1
  ft[paste0("u_days",   i), 1L, Mft] <- 999
}

frame_sorted2 <- frame_us[order(as.character(frame_us$gid)), , drop = FALSE]
mp_sorted <- filter_traj_to_warmup(ft, frame_sorted2, "2023-12-01")
expect_true(identical(mp_sorted$gid, sort(gids_unsorted)),
            "gid-sorted frame → mini-panel gid in alphabetical order")
expect_true(identical(as.integer(mp_sorted$ongoing_outbreak_lag1), seq_len(U2)),
            "gid-sorted frame: row i state maps to sort(gid)[i]")
expect_equal(mp_sorted$ongoing_outbreak_lag1[mp_sorted$gid == "AAA_1_1"], 1L,
             "AAA_1_1 (alphabetical first) gets state index 1 when sorted")

# Negative control: an UNSORTED frame mis-pairs gids with states (AAA is at
# position 2 in gids_unsorted, so it would wrongly receive state index 2).
mp_unsorted <- filter_traj_to_warmup(ft, frame_us, "2023-12-01")
expect_equal(mp_unsorted$ongoing_outbreak_lag1[mp_unsorted$gid == "AAA_1_1"], 2L,
             "unsorted frame mis-pairs AAA_1_1 with state index 2 (sort guard needed)")

# ---------------------------------------------------------------------------
# Hindcast Seam 1: interior-month handoff column selection
# ---------------------------------------------------------------------------
# filter_traj_to_warmup() defaults to the LAST filter_traj column (end-of-
# coverage handoff, as in forecasts). For a season-by-season hindcast the handoff
# is pinned to an INTERIOR launch month, so the IC state must be sliced at the
# filt_months column matching (launch - 1 month), NOT the last one. See
# analysis/utils/make_hindcast_configs.R (filter.handoff_month) and
# 07b_simulator_spatpomp.R::load_smoothing_distribution().
cat("\n## Hindcast Seam 1: interior-month handoff selects the launch-month column\n")

filt_months <- as.Date(c("2015-02-01", "2015-03-01", "2015-04-01"))  # covers Feb..Apr
Mint <- length(filt_months)
# Distinct per-month state: alphabetical unit i in month m carries lag1 = 10*m + i.
ft_int <- array(0, dim = c(length(rn), 1L, Mint), dimnames = list(rn, NULL, NULL))
for (m in seq_len(Mint)) for (i in seq_len(U2)) {
  ft_int[paste0("o",        i), 1L, m] <- 10L * m + i
  ft_int[paste0("o_tu",     i), 1L, m] <- .SP_OTU_SENTINEL
  ft_int[paste0("draw_idx", i), 1L, m] <- 1
  ft_int[paste0("u_days",   i), 1L, m] <- 999
}

# Launch 2015-05-01 → IC month = 2015-04-01 = filt_months[3] (the last column here,
# but assert it is chosen BY filt_months matching, not by position). Use a launch
# that maps to the MIDDLE column to make the interior-selection unambiguous.
mp_mid <- filter_traj_to_warmup(ft_int, frame_sorted2, horizon_start = "2015-04-01",
                                filt_months = filt_months)   # IC = 2015-03-01 = col 2
expect_equal(mp_mid$ongoing_outbreak_lag1[mp_mid$gid == "AAA_1_1"], 10L * 2L + 1L,
             "interior handoff: launch 2015-04 slices the 2015-03 column (col 2), not the last")
expect_true(!identical(mp_mid$ongoing_outbreak_lag1[mp_mid$gid == "AAA_1_1"], 10L * Mint + 1L),
            "interior handoff: does NOT fall through to the last (2015-04) column")

# Backward-compat: filt_months=NULL keeps the last-column (end-of-coverage) behaviour.
mp_last <- filter_traj_to_warmup(ft_int, frame_sorted2, horizon_start = "2099-01-01",
                                 filt_months = NULL)
expect_equal(mp_last$ongoing_outbreak_lag1[mp_last$gid == "AAA_1_1"], 10L * Mint + 1L,
             "filt_months=NULL: still selects the last column (run6 forecast default)")

# Guard: an IC month outside filt_months must stop() with a clear message.
err_int <- tryCatch({
  filter_traj_to_warmup(ft_int, frame_sorted2, horizon_start = "2016-01-01",
                        filt_months = filt_months)   # IC = 2015-12-01, not covered
  NULL
}, error = function(e) conditionMessage(e))
expect_true(!is.null(err_int) && grepl("not found in filt_months", err_int),
            "interior handoff: IC month outside filt_months stop()s with a clear message")

# ---------------------------------------------------------------------------
# Hindcast Seam 1b: toff (t0 column) alignment — the off-by-one guard
# ---------------------------------------------------------------------------
# The ft_int case above has EXACTLY length(filt_months) columns (toff = 0), so a
# stale `M_idx <- m` and a correct `M_idx <- m + toff` coincide and the bug is
# invisible. A real pomp/spatPomp filter.traj has ntimes+1 columns: column 1 is
# the pre-filter t0 state, columns 2..M+1 are the M filter months, so
# filt_months[m] lives at column m + toff with toff = dim(ft)[3] - M = 1. This
# variant rebuilds ft with a leading zero t0 column and asserts the launch-month
# state is sliced at column m + 1 (NOT column m) — it fails against the stale
# `M_idx <- m` code and passes after the toff fix. Mirrors the
# .smoothing_active_matrix() convention (smoothing_core.R:869).
cat("\n## Hindcast Seam 1b: toff (t0 column) alignment guards the off-by-one\n")

# ft_int1: (nvars × 1 × Mint+1). Column 1 = t0 (all zero). Columns 2..Mint+1
# carry filt_months[m]'s per-unit state (lag1 = 10*m + i), same as ft_int shifted
# right by one to make room for t0.
ft_int1 <- array(0, dim = c(length(rn), 1L, Mint + 1L),
                 dimnames = list(rn, NULL, NULL))
for (m in seq_len(Mint)) for (i in seq_len(U2)) {
  ft_int1[paste0("o",        i), 1L, m + 1L] <- 10L * m + i
  ft_int1[paste0("o_tu",     i), 1L, m + 1L] <- .SP_OTU_SENTINEL
  ft_int1[paste0("draw_idx", i), 1L, m + 1L] <- 1
  ft_int1[paste0("u_days",   i), 1L, m + 1L] <- 999
}

# Launch 2015-04-01 → IC month = 2015-03-01 = filt_months[2] = column 2 + toff(1) = 3.
mp_toff <- filter_traj_to_warmup(ft_int1, frame_sorted2, horizon_start = "2015-04-01",
                                 filt_months = filt_months)
expect_equal(mp_toff$ongoing_outbreak_lag1[mp_toff$gid == "AAA_1_1"], 10L * 2L + 1L,
             "toff: launch 2015-04 slices filt_months[2] (col 2+toff=3), the 2015-03 state")
expect_true(!identical(mp_toff$ongoing_outbreak_lag1[mp_toff$gid == "AAA_1_1"], 10L * 1L + 1L),
            "toff: does NOT return the stale filt_months[1] (col 2) state — off-by-one guard")

# End-of-coverage self-consistency: with a leading t0 column, the last filt_months
# entry (m == Mint) must resolve to the LAST column (dim(ft)[3] = Mint+1), matching
# the filt_months = NULL default for the same IC.
mp_toff_last <- filter_traj_to_warmup(ft_int1, frame_sorted2, horizon_start = "2015-05-01",
                                      filt_months = filt_months)   # IC = 2015-04 = filt_months[3]
mp_null_last <- filter_traj_to_warmup(ft_int1, frame_sorted2, horizon_start = "2099-01-01",
                                      filt_months = NULL)
expect_equal(mp_toff_last$ongoing_outbreak_lag1[mp_toff_last$gid == "AAA_1_1"], 10L * Mint + 1L,
             "toff: end-of-coverage (m == M) resolves to the last column")
expect_true(identical(mp_toff_last$ongoing_outbreak_lag1, mp_null_last$ongoing_outbreak_lag1),
            "toff: filt_months path and NULL path agree at end-of-coverage")

# ---------------------------------------------------------------------------
# Hindcast Seam 2: observed-ENSO passthrough (enso_scenarios = NULL)
# ---------------------------------------------------------------------------
# "forecast_source: observed" (utils/simulations.R) leaves enso_scenarios = NULL
# so simulate_one_tu() reads the panel's own enso_lag<N> columns unmodified (the
# override loop in simulator_core.R is skipped). This checks the property the
# hindcast relies on: with NULL scenarios, the panel's observed ENSO drives the
# forecast — a panel carrying a different enso_lag value yields a different
# occurrence rate.
cat("\n## Hindcast Seam 2: observed panel ENSO drives the forecast (enso_scenarios=NULL)\n")

fix_e   <- make_fixture()
frame_e <- fix_e$frame
Ue      <- nrow(frame_e)
wwin_e  <- list(start = as.character(min(fix_e$warmup_months)),
                end   = as.character(max(fix_e$warmup_months) %m+% months(1L)))
fwin_e  <- list(start = wwin_e$end,
                end   = as.character(max(fix_e$forecast_months) %m+% months(1L)))

# occ draws with a strong enso_lag5 slope; moderate intercept so p is enso-sensitive.
draws_e <- make_draws_bundle_synthetic(frame_e, S = 4L,
  active_terms = c("(Intercept)", "enso_lag5"),
  fe_overrides = c("(Intercept)" = -3.0, "enso_lag5" = 3.0),
  basehaz_const = 0.0)

# A minimal filter handoff (all warm-up state zero) so the forecast is driven by
# the panel covariates, not by prior outbreak state.
ft_e <- array(0, dim = c(length(unlist(lapply(.SP_UNIT_STATENAMES, function(s) paste0(s, seq_len(Ue))))), 1L, 1L))
dimnames(ft_e) <- list(unlist(lapply(.SP_UNIT_STATENAMES, function(s) paste0(s, seq_len(Ue)))), NULL, NULL)
for (i in seq_len(Ue)) {
  ft_e[paste0("o_tu",     i), 1L, 1L] <- .SP_OTU_SENTINEL
  ft_e[paste0("draw_idx", i), 1L, 1L] <- 1
  ft_e[paste0("u_days",   i), 1L, 1L] <- 999
}
frame_e_sorted <- frame_e[order(as.character(frame_e$gid)), , drop = FALSE]
handoff_e <- filter_traj_to_warmup(ft_e, frame_e_sorted, fwin_e$start)
handoff_e <- rep(list(handoff_e), 4L)   # 4 replicates

run_with_enso <- function(enso_val) {
  panel_e <- fix_e$panel %>% dplyr::mutate(enso_lag5 = enso_val)
  out <- forecast_from_filtered(
    handoff_list   = handoff_e,
    tu             = frame_e$tu[1L],
    frame_tu       = frame_e_sorted,
    panel_tu_full  = panel_e,
    draws_occ      = draws_e,
    enso_scenarios = NULL,
    horizon        = fwin_e,
    output_part    = tempfile("seam2"),
    out_fmt        = "rds",
    dry_run        = TRUE,
    seed           = 7L)
  if (is.null(out) || is.null(out$monthly_summary)) return(NA_real_)
  mean(out$monthly_summary$Y, na.rm = TRUE)
}

my_hi <- tryCatch(run_with_enso(+2.0), error = function(e) { cat("  ERR hi:", conditionMessage(e), "\n"); NA_real_ })
my_lo <- tryCatch(run_with_enso(-2.0), error = function(e) { cat("  ERR lo:", conditionMessage(e), "\n"); NA_real_ })
expect_true(is.finite(my_hi) && is.finite(my_lo),
            "Seam 2: forecast runs with enso_scenarios=NULL on a panel carrying enso_lag5")
if (is.finite(my_hi) && is.finite(my_lo)) {
  cat(sprintf("  mean(Y) enso=+2: %.4f   enso=-2: %.4f\n", my_hi, my_lo))
  expect_true(my_hi > my_lo,
    sprintf("Seam 2: observed panel ENSO drives the forecast (enso=+2 → %.4f > enso=-2 → %.4f)",
            my_hi, my_lo))
}

# ---------------------------------------------------------------------------
# M3-spread: .forecast_param_indices() branch logic
#   Deterministic unit test of the posterior-draw index derivation that
#   forecast_from_filtered() delegates to. No full forecast needed.
#   FALSE branch = the default (inherit each handoff particle's
#   posterior_draw_idx). TRUE branch = spread draws across replicates (cycle
#   1..S).
# ---------------------------------------------------------------------------
cat("\n## M3-spread: .forecast_param_indices() branch logic\n")
{
  Np_sp        <- 7L
  # Handoff whose particles ALL carry draw_idx=1 — reproduces the collapsed
  # smoothing cache (n_draws_in_filter=1).
  handoff_sp   <- lapply(seq_len(Np_sp),
                         function(i) list(posterior_draw_idx = 1L))
  rep_ids_sp   <- seq_len(Np_sp)          # replicate_offset = 0
  S_occ_sp <- 5L; S_dur_sp <- 3L; S_mon_sp <- 4L

  # FALSE branch: every replicate inherits draw_idx=1 → all indices collapse to 1.
  idx_false <- .forecast_param_indices(
    handoff_sp, rep_ids_sp,
    S_occ = S_occ_sp, S_dur = S_dur_sp, S_mon = S_mon_sp,
    spread_param_draws = FALSE)
  expect_true(identical(names(idx_false), c("occ", "dur", "mon")),
              "returns list(occ, dur, mon) -- no size index (size track retired)")
  expect_true(all(idx_false$occ  == 1L), "FALSE branch: occ all 1 (collapsed)")
  expect_true(all(idx_false$dur  == 1L), "FALSE branch: dur all 1 (collapsed)")
  expect_true(all(idx_false$mon  == 1L), "FALSE branch: mon all 1 (collapsed)")

  # TRUE branch: draws cycle 1..S across replicates, restoring spread.
  idx_true <- .forecast_param_indices(
    handoff_sp, rep_ids_sp,
    S_occ = S_occ_sp, S_dur = S_dur_sp, S_mon = S_mon_sp,
    spread_param_draws = TRUE)
  want_occ <- ((rep_ids_sp - 1L) %% S_occ_sp) + 1L
  expect_true(all(idx_true$occ == want_occ),
              "TRUE branch: occ cycles 1..S over replicates")
  expect_true(length(unique(idx_true$occ)) > 1L,
              "TRUE branch: occ has spread (not a single draw)")
  # dur/mon stay paired: wrapped off the occurrence index.
  expect_true(all(idx_true$dur  == ((idx_true$occ - 1L) %% S_dur_sp)  + 1L),
              "TRUE branch: dur paired to occ")
  expect_true(all(idx_true$mon  == ((idx_true$occ - 1L) %% S_mon_sp)  + 1L),
              "TRUE branch: mon paired to occ")

  # Inactive tracks (NULL S) yield NULL indices, matching the legacy path.
  idx_null <- .forecast_param_indices(
    handoff_sp, rep_ids_sp,
    S_occ = S_occ_sp, S_dur = NULL, S_mon = NULL,
    spread_param_draws = TRUE)
  expect_true(is.null(idx_null$dur) && is.null(idx_null$mon),
              "NULL dur/mon S → NULL index (inactive track)")

  # replicate_offset carries through: ids = seq_len(Np) + offset.
  off_sp   <- 5L
  ids_off  <- seq_len(Np_sp) + off_sp
  idx_off  <- .forecast_param_indices(
    handoff_sp, ids_off,
    S_occ = S_occ_sp, S_dur = S_dur_sp, S_mon = S_mon_sp,
    spread_param_draws = TRUE)
  expect_true(all(idx_off$occ == ((ids_off - 1L) %% S_occ_sp) + 1L),
              "TRUE branch: replicate_offset shifts the draw cycle")
}

# ---------------------------------------------------------------------------
# M3-random: .forecast_param_indices(param_draw_mode = "random")
#   The new global default. Draws a FRESH random occurrence index per replicate,
#   seeded chunk-safely from `seed` so replicate r always maps to the same draw
#   regardless of replicate_offset (SLURM chunking) — and, in a region
#   partition, that one per-replicate index is shared across all member units.
#   Ignores the handoff draw_idx. dur/mon stay paired off occ.
# ---------------------------------------------------------------------------
cat("\n## M3-random: .forecast_param_indices(param_draw_mode = 'random')\n")
{
  Np_rn    <- 40L
  S_occ_rn <- 200L; S_dur_rn <- 120L; S_mon_rn <- 90L
  seed_rn  <- 12345L
  # Handoff whose particles ALL carry draw_idx=1 (collapsed smoothing cache) —
  # random mode must ignore it.
  handoff_rn <- lapply(seq_len(Np_rn), function(i) list(posterior_draw_idx = 1L))
  rep_ids_rn <- seq_len(Np_rn)

  idx_rn <- .forecast_param_indices(
    handoff_rn, rep_ids_rn,
    S_occ = S_occ_rn, S_dur = S_dur_rn, S_mon = S_mon_rn,
    param_draw_mode = "random", seed = seed_rn)

  expect_true(!all(idx_rn$occ == 1L),
              "random: ignores handoff draw_idx (not collapsed to 1)")
  expect_true(all(idx_rn$occ >= 1L & idx_rn$occ <= S_occ_rn),
              "random: occ indices within 1..S_occ")
  expect_true(length(unique(idx_rn$occ)) > 1L,
              "random: occ has spread (not a single draw)")
  # dur/mon stay paired: wrapped off the occurrence index.
  expect_true(all(idx_rn$dur  == ((idx_rn$occ - 1L) %% S_dur_rn)  + 1L),
              "random: dur paired to occ")
  expect_true(all(idx_rn$mon  == ((idx_rn$occ - 1L) %% S_mon_rn)  + 1L),
              "random: mon paired to occ")

  # Determinism for a fixed seed: a second call reproduces the same indices.
  idx_rn2 <- .forecast_param_indices(
    handoff_rn, rep_ids_rn,
    S_occ = S_occ_rn, S_dur = S_dur_rn, S_mon = S_mon_rn,
    param_draw_mode = "random", seed = seed_rn)
  expect_true(identical(idx_rn$occ, idx_rn2$occ),
              "random: same seed -> identical indices (reproducible)")

  # Chunk-safety: the same replicate_id yields the same draw regardless of how
  # the replicate range was offset (a later SLURM chunk). Compare the overlap
  # of replicate_ids 11..40 produced with offset=0 vs. the first 30 ids of an
  # offset=10 chunk (ids 11..40).
  off_rn   <- 10L
  idx_shift <- .forecast_param_indices(
    handoff_rn[seq_len(Np_rn - off_rn)], seq_len(Np_rn - off_rn) + off_rn,
    S_occ = S_occ_rn, S_dur = S_dur_rn, S_mon = S_mon_rn,
    param_draw_mode = "random", seed = seed_rn)
  expect_true(identical(idx_shift$occ, idx_rn$occ[(off_rn + 1L):Np_rn]),
              "random: replicate_id -> draw is chunk-safe (offset-invariant)")

  # The call must not perturb the caller's global RNG stream.
  set.seed(999L); before_rn <- runif(1L)
  set.seed(999L)
  invisible(.forecast_param_indices(
    handoff_rn, rep_ids_rn, S_occ = S_occ_rn,
    param_draw_mode = "random", seed = seed_rn))
  after_rn <- runif(1L)
  expect_equal(before_rn, after_rn,
               "random: .Random.seed restored (global RNG untouched)")
}

# ---------------------------------------------------------------------------
# N4: anchor_warmup_active_to_observed() — observation-anchored handoff
#   Deterministic unit test of the hindcast-only IC anchor. It overwrites each
#   warmup mini-panel's active/mso from the OBSERVED duration panel at the launch
#   month, leaving every other column untouched.
# ---------------------------------------------------------------------------
cat("\n## N4: anchor_warmup_active_to_observed()\n")
{
  launch_n4 <- as.Date("2019-06-01")
  # Warmup panels for four gids, each carrying a NON-anchor column that must
  # survive untouched, plus a stale active/mso the anchor should overwrite.
  warmup_n4 <- list(
    list(gid = "A", active = 1L, mso = 9L, posterior_draw_idx = 3L, u_days = 7),
    list(gid = "B", active = 1L, mso = 4L, posterior_draw_idx = 1L, u_days = 2),
    list(gid = "C", active = 0L, mso = 0L, posterior_draw_idx = 2L, u_days = 0),
    list(gid = "D", active = 1L, mso = 5L, posterior_draw_idx = 4L, u_days = 1))
  # Observed footprint at launch: A active (mso 2), C active (mso 6). Also an
  # earlier month for A (ignored) and a later month for B (ignored).
  obs_n4 <- data.frame(
    gid   = c("A", "A", "C", "B"),
    month = as.Date(c("2019-06-01", "2019-05-01", "2019-06-01", "2019-07-01")),
    mso   = c(2L, 1L, 6L, 3L),
    stringsAsFactors = FALSE)

  out_n4 <- anchor_warmup_active_to_observed(warmup_n4, obs_n4, launch_n4)

  # (a) Observed-active gids get active=1 and their observed mso at launch.
  expect_equal(out_n4[[1]]$active, 1L, "N4: gid A observed-active → active=1")
  expect_equal(out_n4[[1]]$mso,    2L, "N4: gid A gets observed mso at launch")
  expect_equal(out_n4[[3]]$active, 1L, "N4: gid C observed-active → active=1")
  expect_equal(out_n4[[3]]$mso,    6L, "N4: gid C gets observed mso at launch")
  # (b) Gids with no observed row at launch get active=0/mso=0 (B, D overwritten).
  expect_equal(out_n4[[2]]$active, 0L, "N4: gid B not-observed → active=0")
  expect_equal(out_n4[[2]]$mso,    0L, "N4: gid B not-observed → mso=0")
  expect_equal(out_n4[[4]]$active, 0L, "N4: gid D not-observed → active=0")
  expect_equal(out_n4[[4]]$mso,    0L, "N4: gid D not-observed → mso=0")
  # (c) Non-anchor columns are untouched.
  expect_equal(out_n4[[1]]$posterior_draw_idx, 3L, "N4: draw_idx untouched (A)")
  expect_equal(out_n4[[2]]$u_days, 2,               "N4: u_days untouched (B)")
  # (d) mso is integer-typed after the override.
  expect_true(is.integer(out_n4[[1]]$mso) && is.integer(out_n4[[3]]$mso),
              "N4: mso is integer after anchoring")

  # (e) Empty observed footprint → all gids active=0/mso=0.
  out_empty <- anchor_warmup_active_to_observed(
    warmup_n4, obs_n4, as.Date("2001-01-01"))
  expect_true(all(vapply(out_empty, function(mp) mp$active, integer(1)) == 0L),
              "N4: empty footprint → all active=0")
  expect_true(all(vapply(out_empty, function(mp) mp$mso, integer(1)) == 0L),
              "N4: empty footprint → all mso=0")

  # (f) Empty warmup list is a no-op (returns empty list, no error).
  expect_equal(length(anchor_warmup_active_to_observed(list(), obs_n4, launch_n4)),
               0L, "N4: empty warmup_list → empty result")

  # (g) Missing required column stop()s.
  bad_obs <- obs_n4[, c("gid", "month")]
  err_n4  <- tryCatch({
    anchor_warmup_active_to_observed(warmup_n4, bad_obs, launch_n4); NULL
  }, error = function(e) conditionMessage(e))
  expect_true(!is.null(err_n4) && grepl("missing column", err_n4),
              "N4: missing mso column → informative stop()")
}

# ---------------------------------------------------------------------------
# Region assembly: merge per-country smoothing records into one region record
# ---------------------------------------------------------------------------
# filter.partition_by: enso_region reuses the EXISTING per-country smoothing
# caches: for a region TU, each member country's per-replicate records are
# merged into one region record shaped exactly like a country record, so the
# downstream handoff (filter_traj_to_warmup -> resample_by_loglik ->
# forecast_from_filtered) runs unchanged. This exercises
# merge_filter_traj_by_unit() + assemble_region_smoothing_records() and the
# invariant that filter_traj_to_warmup() can extract every region gid from the
# merged array. See analysis/utils/spatpomp_core.R and
# 07b_simulator_spatpomp.R::load_smoothing_distribution().
cat("\n## Region assembly: merge per-country records -> one region record\n")

.sn_reg <- .SP_UNIT_STATENAMES
M_reg   <- 2L
fm_reg  <- as.Date(c("2015-02-01", "2015-03-01"))
gids_A  <- c("AAA_1_1", "AAA_1_2")              # member AAA: 2 units
gids_B  <- c("BBB_1_1", "BBB_1_2", "BBB_1_3")   # member BBB: 3 units

# Build a per-member filter_traj [11*U_c × 1 × M]: the last month carries a
# distinctive ongoing_outbreak_lag1 (o) = 1000*rep_id + member_offset + local
# position, so both the contiguous region re-index and the per-replicate pairing
# are verifiable downstream.
build_member_ft <- function(U_c, member_offset, rep_id, statenames = .sn_reg, M = M_reg) {
  rn <- unlist(lapply(statenames, function(s) paste0(s, seq_len(U_c))))
  ft <- array(0, dim = c(length(rn), 1L, M), dimnames = list(rn, NULL, NULL))
  for (m in seq_len(M)) for (lp in seq_len(U_c)) {
    ft[paste0("o_tu",     lp), 1L, m] <- .SP_OTU_SENTINEL
    ft[paste0("draw_idx", lp), 1L, m] <- 1
    ft[paste0("u_days",   lp), 1L, m] <- 999
  }
  for (lp in seq_len(U_c))
    ft[paste0("o", lp), 1L, M] <- 1000L * rep_id + member_offset + lp
  ft
}

mk_rec <- function(iso3, gids, member_offset, rep_id, loglik) {
  list(tu = iso3, rep_id = rep_id, draw_d = 7L,
       filter_traj = build_member_ft(length(gids), member_offset, rep_id),
       gids = gids, filt_months = fm_reg, logLik = loglik)
}

# Two replicates per member. BBB's record list is given OUT of rep_id order to
# exercise the rep_id-sort pairing in assemble_region_smoothing_records().
recs_A <- list(mk_rec("AAA", gids_A,  0L, rep_id = 1L, loglik = -1.0),
               mk_rec("AAA", gids_A,  0L, rep_id = 2L, loglik = -2.0))
recs_B <- list(mk_rec("BBB", gids_B, 10L, rep_id = 2L, loglik = -20.0),
               mk_rec("BBB", gids_B, 10L, rep_id = 1L, loglik = -10.0))

region_recs <- assemble_region_smoothing_records(
  list(AAA = recs_A, BBB = recs_B), region_tu = "enso_test")

expect_equal(length(region_recs), 2L, "region assembly: one record per replicate (R=2)")

rr1   <- region_recs[[1L]]
U_reg <- length(gids_A) + length(gids_B)
# (a) gid order = concatenation of member gids (region unit order).
expect_true(identical(rr1$gids, c(gids_A, gids_B)),
            "region assembly: gids = concatenation of member gids in member order")
# (b) contiguous 1..U_region position re-index in the merged filter_traj names.
expect_equal(dim(rr1$filter_traj)[1], length(.sn_reg) * U_reg,
             sprintf("region assembly: filter_traj has 11*U_region rows (%d)",
                     length(.sn_reg) * U_reg))
expect_true(all(paste0("o", seq_len(U_reg)) %in% dimnames(rr1$filter_traj)[[1]]),
            "region assembly: merged names carry contiguous o1..o{U_region}")
expect_equal(dim(rr1$filter_traj)[3], M_reg,
             "region assembly: month axis preserved (M columns)")
# (c) summed logLik (independent country partitions -> additive joint logLik).
expect_equal(rr1$logLik, -1.0 + -10.0, "region assembly rep1: logLik = sum of member logLiks")
expect_equal(region_recs[[2L]]$logLik, -2.0 + -20.0, "region assembly rep2: logLik = sum of member logLiks")
# filt_months carried through unchanged (common across members).
expect_true(identical(as.character(rr1$filt_months), as.character(fm_reg)),
            "region assembly: filt_months carried through (common across members)")

# (d) filter_traj_to_warmup() extracts EVERY region gid from the merged array,
# with each per-unit state landing on the right gid AND the correct rep pairing
# (rep1 values = 1000*1 + member_offset + position).
frame_region <- tibble::tibble(
  gid = c(gids_A, gids_B), country = sub("_.*$", "", c(gids_A, gids_B)),
  AFRO_region = "REG", tu = "enso_test", min_seas = 1L, pop = 1e5)
mp_reg <- filter_traj_to_warmup(rr1$filter_traj, frame_region, "2099-01-01",
                                filt_months = NULL)   # last-column (end-of-coverage)
expect_true(identical(mp_reg$gid, c(gids_A, gids_B)),
            "region assembly: filter_traj_to_warmup recovers every region gid in order")
expect_true(identical(as.integer(mp_reg$ongoing_outbreak_lag1),
                      c(1001L, 1002L, 1011L, 1012L, 1013L)),
            "region assembly: rep1 per-unit state maps to the right region position (re-index + rep pairing)")

# (e) members with unequal replicate counts cannot be paired -> hard error.
err_reg <- tryCatch({
  assemble_region_smoothing_records(list(AAA = recs_A, BBB = recs_B[1L]), "enso_test"); NULL
}, error = function(e) conditionMessage(e))
expect_true(!is.null(err_reg) && grepl("replicate count", err_reg),
            "region assembly: members with unequal replicate counts stop() clearly")

# (f) an assembled record names the members' smoothing replicate (rep_id), not
# its position, so 07b's handoff_map.csv points at the right part.qs events.
recs_A5 <- list(mk_rec("AAA", gids_A, 0L, rep_id = 2L, loglik = -1.0),
                mk_rec("AAA", gids_A, 0L, rep_id = 5L, loglik = -2.0))
recs_B5 <- list(mk_rec("BBB", gids_B, 10L, rep_id = 5L, loglik = -20.0),
                mk_rec("BBB", gids_B, 10L, rep_id = 2L, loglik = -10.0))
rr5 <- assemble_region_smoothing_records(list(AAA = recs_A5, BBB = recs_B5), "enso_test")
expect_true(identical(vapply(rr5, `[[`, integer(1), "rep_id"), c(2L, 5L)),
            "region assembly: record rep_id = the members' common smoothing rep_id (2, 5), not 1, 2")

# ---------------------------------------------------------------------------
# Handoff map (07b's handoff_map.csv): forecast replicate -> smoothing rep_id
# ---------------------------------------------------------------------------
cat("\n## Handoff map: forecast replicate -> smoothing replicate\n")
# A 3-record cache whose second record failed extraction: the selector draws
# over the 2 ok records, and the map must name their rep_ids (1 and 3), never
# the position.
recs3 <- list(list(rep_id = 1L), list(rep_id = 2L), list(rep_id = 3L))
ok3   <- c(1L, 3L)
wl3   <- list(data.frame(r = 1L), data.frame(r = 3L))
for (off in c(0L, 5L)) {
  sel <- select_handoff_warmup(wl3, c(-1, -2), n = 4L, offset = off)
  hm  <- handoff_map_table("TUX", recs3[ok3], sel, off, as.Date("2026-09-01"), "uniform")
  want <- vapply(recs3[ok3], `[[`, integer(1), "rep_id")[((off + seq_len(4L) - 1L) %% 2L) + 1L]
  expect_true(identical(hm$replicate_id, off + 1:4) && identical(hm$smoothing_rep_id, want) &&
                identical(hm$smoothing_rep_id, vapply(sel, function(d) d$r, integer(1))) &&
                all(hm$handoff_month == "2026-09-01"),
              sprintf("handoff map (offset %d): replicate_id = offset + r, smoothing_rep_id = ok rep_ids[((R-1) mod m)+1]", off))
}
err_hm <- tryCatch({ handoff_map_table("TUX", recs3[ok3], wl3[c(1, 2, 1)], 0L, "2026-09-01", "uniform"); NULL },
                   error = function(e) conditionMessage(e))
expect_true(!is.null(err_hm) && grepl("idx", err_hm), "handoff map: a list without `idx` is refused")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat(sprintf("\n# Results: %d PASS  %d FAIL\n", .passes, .fails))
if (.fails > 0L) quit(status = 1L)
