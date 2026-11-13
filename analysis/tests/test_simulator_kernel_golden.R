# test_simulator_kernel_golden.R -- golden characterization of the forward
# kernel simulate_one_tu() and its step helpers (analysis/utils/simulator_core.R).
#
# HARD gate (see analysis/tests/run_all.R). Guards the lookups inside
# evaluate_lp(), update_lags(), advance_ttlo(), the per-step covariate
# assembly and the event bookkeeping. Any change must leave the output
# bit-identical at a fixed seed (same RNG draws, same order), so this test runs
# deterministic synthetic worlds that between them reach every kernel branch the
# deployed runs use:
#   S1  grouped (4 fit groups, region dropped), occ + dur + mon, ENSO ns(df=4)
#       single-lag cross-basis on every track, latent mso basis, cyclic seasonal
#       spline, country/ADM1/ADM2 intercept REs, country + ADM1 RE slopes on
#       every cb column, own-o + FALSE slope, 7 TTLO levels incl. "none",
#       log_pop_c, nested phi REs, per-replicate warm-ups and ENSO scenarios,
#       covariate_feed "detected", force_o_tu, incidence cap, neighbour
#       contagion (synthetic ADM1 kernel; GADM is not read);
#   S2  ungrouped occurrence-only legacy path: AFRO-region-prefixed RE keys with
#       some matrices keyed the alternate way (dual-key fallback), climate
#       linear terms, o / o_ct / o_tu spillovers with FALSE slopes, NA-TU rows,
#       enso_lag<n> linear term, constant basehaz;
#   S3  ungrouped occ + dur, latent feed, shared warm-up, drop_country_tu_spillover;
#   S4  ungrouped occ + dur + mon, pop_year mode (pop_y<YYYY> columns), constant
#       phi, no cap;
# plus direct calls of evaluate_lp(), update_lags(), advance_ttlo(),
# compute_phi_vec(), basehaz_dur_eval() and monthlyLogMu() (shared with the
# spatPomp R builder, the Csnippet static precompute and smoothing).
# The captured list is compared with identical() against
# analysis/tests/golden/simulator_kernel_golden.rds. Each world is also run
# through the frozen pre-speed-up kernel (utils/simulator_core_legacy.R,
# simulate_one_tu_legacy()) on the same inputs, and the rewritten helpers are
# compared with their legacy copies on edge cases the worlds do not reach
# (NA / non-0/1 lags, warm-up ties, gids absent before the horizon, NA TTLO).
#
#   Rscript analysis/tests/test_simulator_kernel_golden.R --capture   # once, pre-refactor
#   Rscript analysis/tests/test_simulator_kernel_golden.R             # compare
#
# Synthetic data only -- no fit, no generated_data, no network.
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tibble); library(tidyr); library(lubridate)
})
cli <- parse_args(OptionParser(option_list = list(
  make_option("--capture", action = "store_true", default = FALSE,
              help = "write the golden file instead of comparing"))))
.fails <- 0L; .passes <- 0L
source("analysis/tests/helpers.R")
for (f in c("base.R", "gid.R", "paths.R", "clim_compile.R", "clim_db.R", "occ_data.R",
            "siz_data.R", "provenance.R", "sim_io.R", "dlnm.R", "population.R",
            "composite_lp.R", "simulator_core.R", "anchor.R", "spatial_contagion.R",
            "simulator_core_legacy.R"))
  source(file.path("analysis/utils", f))

GOLDEN_PATH <- "analysis/tests/golden/simulator_kernel_golden.rds"
cat("# simulator kernel golden (mode:", if (cli$capture) "capture" else "compare", ")\n")

# Synthetic neighbour kernel: the real one reads GADM polygons. Deterministic,
# dense, zero diagonal, ADM1 units in unique(frame_tu$adm1_unit) order.
nbr_kernel_for_frame <- function(frame_tu, spec) {
  a1 <- unique(as.character(frame_tu$adm1_unit))
  n  <- length(a1)
  K  <- outer(seq_len(n), seq_len(n), function(i, j) exp(-abs(i - j) / 3))
  diag(K) <- 0
  K <- K / rowSums(K)
  dimnames(K) <- list(a1, a1)
  K
}

TTLO_TERMS <- c("time_to_last_outbreak12", "time_to_last_outbreak23", "time_to_last_outbreak34",
                "time_to_last_outbreak45", "time_to_last_outbreak56", "time_to_last_outbreak6Inf")

mk_frame <- function(n_ctry, n_adm1, n_adm2, groups = NULL, regions = "REG", na_tu_ctry = NULL) {
  fr <- expand.grid(ci = seq_len(n_ctry), a1 = seq_len(n_adm1), a2 = seq_len(n_adm2)) %>%
    mutate(country = sprintf("C%02d", ci), adm1_unit = sprintf("%s_%d", country, a1),
           gid = sprintf("%s_%d", adm1_unit, a2),
           AFRO_region = regions[(ci - 1L) %% length(regions) + 1L],
           tu = ifelse(country %in% na_tu_ctry, NA_character_, "hmmX"),
           min_seas = ((ci - 1L) %% 4L) + 1L, pop = 1e5 * (1 + (a2 %% 5)),
           admin_level = "ADM2",
           first_onset = as.Date("2015-01-01") + 30L * ((ci + a1) %% 3L))
  fr$first_onset[fr$a2 == 1L & fr$ci == 1L] <- as.Date(NA)   # never-onset units
  if (!is.null(groups)) fr$group <- groups[(fr$ci - 1L) %% length(groups) + 1L]
  fr <- fr[order(fr$gid), , drop = FALSE]
  tibble::as_tibble(fr[, setdiff(names(fr), c("ci", "a1", "a2"))])
}

mk_panel <- function(frame, months, enso = sin(seq_along(months) / 2), ttlo = "[6,Inf]",
                     extra = NULL) {
  p <- tidyr::expand_grid(gid = frame$gid, t_start = months)
  n <- nrow(p)
  p <- p %>% mutate(
    t_start_seas = lubridate::yday(t_start),
    mean_precip = cos(seq_len(n) / 7), dev_precip = sin(seq_len(n) / 5),
    mean_temp_daymean = cos(seq_len(n) / 11), dev_temp_daymean = sin(seq_len(n) / 13),
    ongoing_outbreak_lag1 = as.integer(seq_len(n) %% 7L == 0L),
    ongoing_outbreak_country_lag1 = as.integer(seq_len(n) %% 5L == 0L),
    ongoing_outbreak_tu_lag1 = as.integer(seq_len(n) %% 3L == 0L),
    time_to_last_outbreak_cat = rep_len(c(ttlo, "[0,1)", "[2,3)", NA), n),
    enso_lag0 = rep(enso, times = nrow(frame)),
    enso_lag3 = rep(enso * 0.7, times = nrow(frame)))
  if (!is.null(extra)) p <- extra(p)
  p
}

# Random-effect matrix keyed the way the kernel looks it up; `alt` keys a share of
# the columns the other way to exercise match_keys()' fallback.
re_mat <- function(S, keys, sd, alt = NULL) {
  m <- matrix(rnorm(S * length(keys), 0, sd), S)
  colnames(m) <- if (is.null(alt)) keys else ifelse(seq_along(keys) %% 3L == 0L, alt, keys)
  m
}
re_keys <- function(frame, bucket, drop_region) unique(re_level_key(frame, bucket, drop_region))
mk_re <- function(frame, S, drop_region, sd = 0.2, levels = c("country", "adm1"), alt = FALSE) {
  out <- list()
  for (b in levels) {
    k <- re_keys(frame, b, drop_region)
    out[[b]] <- re_mat(S, k, sd, if (alt) re_keys(frame, b, !drop_region) else NULL)
  }
  out
}

# ---------------------------------------------------------------------------
# Bundles
# ---------------------------------------------------------------------------
mk_cb_meta <- function(panel, with_mso = FALSE) {
  vars <- list(enso = list(single_lag = 0L))
  if (with_mso) vars$mso <- list(single_lag = 0L, argvar = list(fun = "bs", degree = 3L, df = 4L))
  dat <- panel
  if (with_mso) dat$mso <- rep_len(1:9, nrow(dat))
  suppressMessages(build_cb_block(dat, list(enable = TRUE,
    defaults = list(argvar = list(fun = "ns", df = 4L)), variables = vars), ""))$cb_meta
}

mk_occ_grouped <- function(frame, S, cbn, pred_fn, K) {
  terms <- c("(Intercept)", "ongoing_outbreak_lag1", TTLO_TERMS, "time_to_last_outbreaknone",
             "nbr_lag1_r500_rn", cbn)
  slopes <- c(setNames(lapply(cbn, function(x) mk_re(frame, S, TRUE)), cbn),
              list(ongoing_outbreak_lag1 = mk_re(frame, S, TRUE),
                   ongoing_outbreak_lag1FALSE = mk_re(frame, S, TRUE),
                   time_to_last_outbreak01 = mk_re(frame, S, TRUE),
                   nbr_lag1_r500_rn = mk_re(frame, S, TRUE)),
              setNames(lapply(TTLO_TERMS, function(x) mk_re(frame, S, TRUE)), TTLO_TERMS))
  b <- make_draws_bundle_synthetic(frame, S = S, active_terms = terms,
         fe_overrides = c(`(Intercept)` = -7.0, ongoing_outbreak_lag1 = 1.0, nbr_lag1_r500_rn = 0.8),
         re_slope_overrides = slopes)
  b$fe[, cbn] <- matrix(rnorm(S * length(cbn), 0, 0.2), S)
  b$fe[, c(TTLO_TERMS, "time_to_last_outbreaknone")] <- matrix(rnorm(S * 7, -0.3, 0.1), S)
  b$re_intercept <- list(region = NULL,
    country = re_mat(S, re_keys(frame, "country", TRUE), 0.3),
    adm1    = re_mat(S, re_keys(frame, "adm1", TRUE), 0.3),
    adm2    = re_mat(S, re_keys(frame, "adm2", TRUE), 0.3))
  b$basehaz <- list(kind = "brms_seasonal_in_lp")
  b$smooth_pred_fn <- list(s1 = pred_fn)
  b$smooth_coef    <- list(s1 = matrix(rnorm(S * K, 0, 0.3), S))
  b
}
mk_dur <- function(frame, S, cbn, drop_region) {
  terms <- c("(Intercept)", "mso_band2", "mso_band3", "mso_band4M5", "mso_band6M8", "mso_band9P", cbn)
  b <- make_draws_bundle_synthetic(frame, S = S, active_terms = terms,
         fe_overrides = c(`(Intercept)` = -3.2),
         re_slope_overrides = setNames(lapply(cbn, function(x) mk_re(frame, S, drop_region)), cbn))
  if (length(cbn)) b$fe[, cbn] <- matrix(rnorm(S * length(cbn), 0, 0.2), S)
  b$fe[, c("mso_band2", "mso_band3", "mso_band4M5", "mso_band6M8", "mso_band9P")] <-
    matrix(rnorm(S * 5, 0.2, 0.2), S)
  b$re_intercept <- list(region = NULL,
    country = re_mat(S, re_keys(frame, "country", drop_region), 0.3),
    adm1    = re_mat(S, re_keys(frame, "adm1", drop_region), 0.3), adm2 = NULL)
  b$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)
  b
}
mk_mon <- function(frame, S, cbn, drop_region, phi_re = "region_country_adm1") {
  terms <- c("(Intercept)", "log_pop_c", "sin_month", "cos_month", "sin2_month", "cos2_month",
             "admin_levelADM2", cbn)
  b <- make_draws_bundle_synthetic(frame, S = S, active_terms = terms,
         fe_overrides = c(`(Intercept)` = -8.5, log_pop_c = 0.1, sin_month = 0.3, cos_month = -0.2,
                          admin_levelADM2 = 0.1),
         re_slope_overrides = setNames(lapply(cbn, function(x) mk_re(frame, S, drop_region)), cbn))
  if (length(cbn)) b$fe[, cbn] <- matrix(rnorm(S * length(cbn), 0, 0.2), S)
  b$re_intercept <- list(region = NULL,
    country = re_mat(S, re_keys(frame, "country", drop_region), 0.3),
    adm1    = re_mat(S, re_keys(frame, "adm1", drop_region), 0.3), adm2 = NULL)
  b$offset_pop <- TRUE; b$offset_logoffset <- TRUE
  b$log_pop_centres <- c(ADM2 = log(2e5))
  b$phi_re <- phi_re
  b$phi <- if (identical(phi_re, "constant")) rep(0.6, S) else
    list(fe = rnorm(S, -0.5, 0.1),
         region = if (!drop_region) re_mat(S, re_keys(frame, "region", FALSE), 0.2) else NULL,
         country = re_mat(S, re_keys(frame, "country", drop_region), 0.2),
         adm1 = re_mat(S, re_keys(frame, "adm1", drop_region), 0.2))
  b
}

cap <- list()

# Run the kernel and its frozen legacy copy on the same arguments; they must agree.
sim_both <- function(label, ...) {
  a <- simulate_one_tu(...)
  b <- simulate_one_tu_legacy(...)
  same <- identical(a, b)
  if (!same) print(all.equal(a, b))
  expect_true(same, sprintf("%s: kernel identical to simulate_one_tu_legacy()", label))
  a
}

# ===========================================================================
# S1 -- grouped, every deployed feature
# ===========================================================================
set.seed(101)
{
  groups <- c("east", "horn", "south", "west")
  fr  <- mk_frame(5, 4, 4, groups = groups)
  M   <- 13L; R <- 5L; S <- 12L
  mo  <- seq(as.Date("2026-02-01"), by = "month", length.out = M + 3L)
  pan <- mk_panel(fr, mo)
  hz  <- list(start = "2026-05-01", end = format(max(mo)))
  cbm_occ <- mk_cb_meta(pan)
  cbm_mon <- mk_cb_meta(pan, with_mso = TRUE)
  cbn     <- cbm_occ$enso$cb_col_names
  sm_raw  <- mgcv::smoothCon(mgcv::s(t_start_seas, bs = "cc", k = 8),
                             data = data.frame(t_start_seas = 0:364), absorb.cons = FALSE)[[1]]
  K       <- ncol(sm_raw$X)
  pred_fn <- function(new_t) mgcv::Predict.matrix(sm_raw, data.frame(t_start_seas = as.numeric(new_t)))
  grouping <- list(
    draws_occ_by_group   = setNames(lapply(groups, function(g) mk_occ_grouped(fr, S, cbn, pred_fn, K)), groups),
    cb_meta_occ_by_group = setNames(rep(list(cbm_occ), 4), groups),
    draws_dur_by_group   = setNames(lapply(groups, function(g) mk_dur(fr, S, cbn, TRUE)), groups),
    cb_meta_dur_by_group = setNames(rep(list(cbm_occ), 4), groups),
    draws_mon_by_group   = setNames(lapply(groups, function(g)
      mk_mon(fr, S, c(cbn, cbm_mon$mso$cb_col_names), TRUE)), groups),
    cb_meta_mon_by_group = setNames(rep(list(cbm_mon), 4), groups))
  # different draw counts per group exercise the per-context draw clamp
  grouping$draws_occ_by_group$horn$S <- 7L
  ens <- list(enso_lag0 = matrix(rnorm((M + 2L) * (R + 2L)), M + 2L, R + 2L))
  wm  <- mo[mo < as.Date(hz$start)]
  wl  <- lapply(seq_len(R + 2L), function(r) {
    p <- mk_panel(fr, wm)
    n <- nrow(p)
    p$active <- as.integer((seq_len(n) + r) %% 9L == 0L)
    p$mso    <- p$active * ((seq_len(n) %% 4L) + 1L)
    p$det    <- p$active * as.integer(seq_len(n) %% 2L == 0L)
    p
  })
  rid <- c(2L, 3L, 4L, 6L, 7L)
  res <- sim_both("S1", "hmmX", fr, pan,
    draws_occ = grouping$draws_occ_by_group[[1]], draws_dur = grouping$draws_dur_by_group[[1]],
    draws_mon = grouping$draws_mon_by_group[[1]],
    cb_meta_occ = cbm_occ, cb_meta_dur = cbm_occ, cb_meta_mon = cbm_mon, grouping = grouping,
    replicate_ids = rid, posterior_idx_occ = c(1L, 5L, 9L, 12L, 3L),
    posterior_idx_dur = c(2L, 2L, 8L, 11L, 4L), posterior_idx_mon = c(3L, 6L, 1L, 10L, 12L),
    seed = 777L, output_part = tempfile(), horizon = hz, dry_run = TRUE,
    enso_scenarios = ens, warmup_panels_list = wl, covariate_feed = "detected",
    delta_fn = 0.2, delta_fp = 0.03, force_o_tu = TRUE, use_nbr = TRUE,
    nbr_spec = list(col = "nbr_lag1_r500_rn", rho_km = 500L),
    clamp_enso_to_knots = TRUE, clamp_enso_margin = 0.25, incidence_cap = 0.02)
  cap$S1 <- res
  expect_true(nrow(res$event_table) > 0L, "S1: events drawn")
  expect_true(any(!is.na(res$monthly_summary$cases_month)), "S1: monthly cases drawn")
  expect_true(any(res$event_table$censored %in% 0L) && any(res$event_table$censored %in% 1L),
              "S1: both closed and censored episodes")
  cap$S1_inputs <- list(fr = fr, grouping = grouping)
}

# ===========================================================================
# S2 -- ungrouped occurrence-only legacy path
# ===========================================================================
set.seed(202)
{
  fr  <- mk_frame(4, 3, 3, regions = c("West Africa", "Central Africa"), na_tu_ctry = "C04")
  M   <- 10L; R <- 4L; S <- 6L
  mo  <- seq(as.Date("2025-10-01"), by = "month", length.out = M + 2L)
  pan <- mk_panel(fr, mo)
  hz  <- list(start = "2025-12-01", end = format(max(mo)))
  terms <- c("(Intercept)", "mean_precip", "dev_precip", "mean_temp_daymean", "dev_temp_daymean",
             "enso_lag3", "ongoing_outbreak_lag1", "ongoing_outbreak_country_lag1",
             "ongoing_outbreak_tu_lag1", TTLO_TERMS)
  sl <- list(enso_lag3 = mk_re(fr, S, FALSE, levels = c("region", "country", "adm1"), alt = TRUE),
             mean_precip = mk_re(fr, S, FALSE, levels = "adm1"),
             ongoing_outbreak_lag1 = mk_re(fr, S, FALSE, levels = c("region", "country")),
             ongoing_outbreak_lag1FALSE = mk_re(fr, S, FALSE, levels = c("region", "country")),
             ongoing_outbreak_country_lag1 = mk_re(fr, S, FALSE, levels = "country"),
             ongoing_outbreak_country_lag1FALSE = mk_re(fr, S, FALSE, levels = "country"),
             ongoing_outbreak_tu_lag1 = mk_re(fr, S, FALSE, levels = "country", alt = TRUE),
             time_to_last_outbreak01 = mk_re(fr, S, FALSE, levels = "country"),
             time_to_last_outbreak23 = mk_re(fr, S, FALSE, levels = "country"))
  d_occ <- make_draws_bundle_synthetic(fr, S = S, active_terms = terms, basehaz_const = -0.4,
             fe_overrides = c(`(Intercept)` = -5.5, mean_precip = 0.2, dev_precip = -0.1,
                              mean_temp_daymean = 0.05, dev_temp_daymean = 0.1, enso_lag3 = 0.3,
                              ongoing_outbreak_lag1 = 0.9, ongoing_outbreak_country_lag1 = 0.4,
                              ongoing_outbreak_tu_lag1 = 0.2),
             re_slope_overrides = sl)
  d_occ$fe[, TTLO_TERMS] <- matrix(rnorm(S * 6, -0.2, 0.1), S)
  d_occ$re_intercept$region  <- re_mat(S, re_keys(fr, "region", FALSE), 0.2)
  d_occ$re_intercept$country <- re_mat(S, re_keys(fr, "country", FALSE), 0.2, re_keys(fr, "country", TRUE))
  d_occ$re_intercept$adm1    <- re_mat(S, re_keys(fr, "adm1", FALSE), 0.2)
  d_occ$re_intercept$adm2    <- re_mat(S, re_keys(fr, "adm2", FALSE), 0.2)
  res <- sim_both("S2", "hmmX", fr, pan, draws_occ = d_occ, replicate_ids = 1:R,
                         seed = 31L, output_part = tempfile(), horizon = hz, dry_run = TRUE)
  cap$S2 <- res
  cap$S2_inputs <- list(fr = fr, d_occ = d_occ, pan = pan)
  expect_true(nrow(res$event_table) > 0L, "S2: events drawn")
}

# ===========================================================================
# S3 -- ungrouped occ + dur, latent feed, shared warm-up, drop ct/tu spillover
# ===========================================================================
set.seed(303)
{
  fr  <- mk_frame(3, 3, 4, regions = c("West Africa", "Central Africa"))
  M   <- 12L; R <- 4L; S <- 5L
  mo  <- seq(as.Date("2025-01-01"), by = "month", length.out = M + 2L)
  pan <- mk_panel(fr, mo, extra = function(p) { p$active <- as.integer(seq_len(nrow(p)) %% 6L == 0L)
                                                p$mso <- p$active * 2L; p })
  hz  <- list(start = "2025-03-01", end = format(max(mo)))
  d_occ <- make_draws_bundle_synthetic(fr, S = S, fe_overrides = c(`(Intercept)` = -4.5,
             ongoing_outbreak_lag1 = 0.7, ongoing_outbreak_country_lag1 = 0.5,
             ongoing_outbreak_tu_lag1 = 0.3, mean_precip = 0.1))
  d_occ$re_intercept$adm1 <- re_mat(S, re_keys(fr, "adm1", FALSE), 0.3)
  d_dur <- mk_dur(fr, S, character(0), FALSE)
  res <- sim_both("S3", "hmmX", fr, pan, draws_occ = d_occ, draws_dur = d_dur,
                         replicate_ids = c(3L, 5L, 8L, 9L), seed = 4242L, output_part = tempfile(),
                         horizon = hz, dry_run = TRUE, drop_country_tu_spillover = TRUE,
                         force_o_tu = TRUE)
  cap$S3 <- res
  expect_true(any(res$monthly_summary$active == 1L), "S3: active episodes")
}

# ===========================================================================
# S4 -- ungrouped occ + dur + mon, pop_year mode, constant phi
# ===========================================================================
set.seed(404)
{
  fr  <- mk_frame(3, 2, 3)
  fr$pop_y2025 <- fr$pop; fr$pop_y2026 <- fr$pop * 1.03; fr$pop_y2027 <- fr$pop * 1.06
  M   <- 14L; R <- 3L; S <- 4L
  mo  <- seq(as.Date("2025-11-01"), by = "month", length.out = M + 1L)
  pan <- mk_panel(fr, mo)
  hz  <- list(start = "2025-12-01", end = format(max(mo)))
  d_occ <- make_draws_bundle_synthetic(fr, S = S, fe_overrides = c(`(Intercept)` = -3.8,
             ongoing_outbreak_lag1 = 0.6))
  d_dur <- mk_dur(fr, S, character(0), FALSE)
  d_mon <- mk_mon(fr, S, character(0), FALSE, phi_re = "constant")
  d_mon$log_pop_centres <- c(ADM2 = log(3e5))
  res <- sim_both("S4", "hmmX", fr, pan, draws_occ = d_occ, draws_dur = d_dur, draws_mon = d_mon,
                         replicate_ids = 1:R, seed = 99L, output_part = tempfile(),
                         horizon = hz, dry_run = TRUE)
  cap$S4 <- res
  expect_true(any(!is.na(res$monthly_summary$cases_month)), "S4: monthly cases drawn")
}

# ===========================================================================
# Direct helper calls
# ===========================================================================
set.seed(505)
{
  g1  <- cap$S1_inputs; fr <- g1$fr[g1$fr$group == "east", , drop = FALSE]; N <- nrow(fr)
  d   <- g1$grouping$draws_occ_by_group$east
  cbn <- grep("_cb[0-9]+$", colnames(d$fe), value = TRUE)
  xr  <- c(setNames(lapply(cbn, function(x) rnorm(N)), cbn), list(t_start_seas = (seq_len(N) * 17) %% 365))
  o   <- as.integer(seq_len(N) %% 3L == 0L); o_ct <- as.integer(seq_len(N) %% 4L == 0L)
  o_tu <- o_ct; o_tu[1:2] <- NA_integer_
  u   <- c(NA, 0, 100, 400, 800, 1200, 1600, 2000, 2400, 1e9, rep(50, N - 10))
  nbr <- seq(0, 1, length.out = N)
  cap$lp_grouped <- lapply(c(1L, 4L), function(s)
    evaluate_lp(d, s, xr, o, o_ct, o_tu, u, fr, basehaz_val = 0.1, drop_region = TRUE, nbr = nbr))
  s2 <- cap$S2_inputs; fr2 <- s2$fr; N2 <- nrow(fr2)
  xr2 <- list(mean_precip = rnorm(N2), dev_precip = c(NA, rnorm(N2 - 1)),
              mean_temp_daymean = rnorm(N2), dev_temp_daymean = rnorm(N2), enso_lag3 = rep(0.4, N2))
  o2 <- as.integer(seq_len(N2) %% 2L == 0L); otu2 <- o2; otu2[fr2$country == "C04"] <- NA_integer_
  u2 <- rep(c(10, 500, 900, 3000), length.out = N2)
  cap$lp_ungrouped <- lapply(1:3, function(s) evaluate_lp(s2$d_occ, s, xr2, o2, rev(o2), otu2, u2, fr2,
                                                          basehaz_val = basehaz_eval(s2$d_occ, 0, N2)))
  cap$lp_ungrouped_dropct <- evaluate_lp(s2$d_occ, 2L, xr2, o2, rev(o2), otu2, u2, fr2, drop_ct_tu = TRUE)
  Y <- as.integer(seq_len(N2) %% 5L == 0L)
  cap$lags <- update_lags(Y, fr2)
  cap$lags0 <- update_lags(integer(N2), fr2)
  cap$ttlo <- advance_ttlo(c(NA, 0, 10, 1e9, 30), c(0L, 1L, 0L, 0L, 1L))
  dm <- g1$grouping$draws_mon_by_group$west
  frw <- g1$fr[g1$fr$group == "west", , drop = FALSE]
  cap$phi <- list(compute_phi_vec(dm, 3L, frw, drop_region = TRUE),
                  compute_phi_vec(mk_mon(fr2, 3, character(0), FALSE, "region_country"), 2L, fr2))
  dd <- g1$grouping$draws_dur_by_group$south
  cap$bh_dur <- basehaz_dur_eval(dd, 2L, c(0L, 1L, 2L, 3L, 5L, 7L, 12L, NA), 8L)
  Nw <- nrow(frw)
  xm <- c(setNames(lapply(cbn, function(x) rnorm(Nw)), cbn),
          list(sin_month = rep(0.5, Nw), cos_month = rep(0.86, Nw), sin2_month = rep(0.86, Nw),
               cos2_month = rep(-0.5, Nw), admin_levelADM2 = rep(1, Nw)))
  mso_tab <- build_mso_basis(g1$grouping$cb_meta_mon_by_group$west)
  cap$mon_mu <- monthlyLogMu(dm, 5L, xm, integer(Nw), integer(Nw), integer(Nw), rep(100, Nw), frw,
                             log_offset = log(31), mso = rep_len(0:20, Nw), mso_tab = mso_tab,
                             drop_region = TRUE)
}

# ===========================================================================
# Rewritten helpers vs their legacy copies, on edge cases (not in the golden)
# ===========================================================================
set.seed(606)
{
  s2 <- cap$S2_inputs; fr2 <- s2$fr; N2 <- nrow(fr2); mp <- lagMaps(fr2)
  for (k in 1:20) {
    Y <- rbinom(N2, 1, runif(1))
    expect_true(identical(update_lags(Y, fr2, maps = mp), update_lags_legacy(Y, fr2)) &&
                identical(update_lags(Y, fr2), update_lags_legacy(Y, fr2)),
                sprintf("update_lags coded == legacy (draw %d)", k))
  }
  for (Y in list(replace(integer(N2), 3, NA), replace(integer(N2), 5, 2L), as.numeric(rbinom(N2, 1, 0.3))))
    expect_true(identical(update_lags(Y, fr2, maps = mp), update_lags_legacy(Y, fr2)),
                "update_lags with NA / non-0/1 / double Y == legacy")
  u <- c(NA, 0, 1e9, 31, 400.5); Yt <- c(1L, NA, 0L, 1L, 0L)
  expect_true(identical(advance_ttlo(u, Yt), advance_ttlo_legacy(u, Yt)) &&
              identical(advance_ttlo(stats::setNames(u, letters[1:5]), Yt),
                        advance_ttlo_legacy(stats::setNames(u, letters[1:5]), Yt)),
              "advance_ttlo == legacy (NA Y, NA days, named days)")
  d <- c(NA, -5, 0, 364, 365, 366, 2190, 2191, 5000, 1e9, Inf, NaN)
  expect_true(identical(bin_ttlo(d), bin_ttlo_legacy(d)) &&
              identical(bin_ttlo(d, TRUE), bin_ttlo_legacy(d, TRUE)), "bin_ttlo == legacy")
  ms <- c(NA, -1, 0, 1, 2, 3, 4, 5, 6, 8, 9, 50, Inf, 2.5)
  expect_true(identical(.dur_bin_mso(ms), .dur_bin_mso_legacy(ms)) &&
              identical(.onset_log_offset_adj(ms), .onset_log_offset_adj_legacy(ms)),
              ".dur_bin_mso / .onset_log_offset_adj == legacy")
  # warm-up row selection: ties at the latest pre-horizon month, a gid with
  # only post-horizon rows, a frame gid absent from the panel, an extra panel
  # gid, factor gids, rows out of order.
  pan <- s2$pan
  hz  <- list(start = "2025-12-01")
  p1  <- rbind(pan, pan[pan$t_start == as.Date("2025-11-01"), ])          # ties
  p1$ongoing_outbreak_lag1[nrow(pan) + 1:3] <- 1L - p1$ongoing_outbreak_lag1[nrow(pan) + 1:3]
  p1  <- p1[!(p1$gid == fr2$gid[2] & p1$t_start < as.Date(hz$start)), ]   # post-only gid
  p1  <- p1[p1$gid != fr2$gid[3], ]                                       # absent gid
  p1  <- p1[sample(nrow(p1)), ]
  p2  <- p1; p2$gid <- factor(p2$gid)
  p3  <- rbind(p1, transform(p1[1:4, ], gid = "ZZZ_1_1"))                 # extra gid
  for (pp in list(p1, p2, p3, pan[0, ], pan[pan$t_start >= as.Date(hz$start), ]))
    expect_true(identical(build_warmup_state(pp, fr2, hz), build_warmup_state_legacy(pp, fr2, hz)),
                "build_warmup_state == legacy (ties / post-only / absent / extra / empty)")
}
cap$S1_inputs <- NULL; cap$S2_inputs <- NULL

# ---------------------------------------------------------------------------
if (cli$capture) {
  dir.create(dirname(GOLDEN_PATH), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cap, GOLDEN_PATH)
  cat(sprintf("captured %d entries -> %s\n", length(cap), GOLDEN_PATH))
} else {
  if (!file.exists(GOLDEN_PATH)) fail(paste("golden missing:", GOLDEN_PATH))
  gold <- readRDS(GOLDEN_PATH)
  expect_true(identical(sort(names(cap)), sort(names(gold))), "same captured entries")
  for (nm in names(gold)) {
    same <- identical(cap[[nm]], gold[[nm]])
    if (!same) print(all.equal(cap[[nm]], gold[[nm]]))
    expect_true(same, sprintf("%s identical to golden", nm))
  }
}
cat(sprintf("\n%d passed, %d failed\n", .passes, .fails))
if (.fails > 0L) quit(status = 1)
