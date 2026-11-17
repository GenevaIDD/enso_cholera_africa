# test_spatpomp_detected_feed.R — HARD: covariate_feed = "detected" and the
#   activity-scored measurement (..._state_detect_act).
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_detected_feed.R
#
# Details:
#   Pure R/C fixtures, run from the project root; no fits, no MCMC. T1-T2
#   statename layouts, T3 the C kernel under the detected feed, T4 semantics,
#   T5 the `_act` measurement density, T6 the config guard.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")
source("analysis/utils/dlnm.R") # build_mso_basis() -- the mon path's mso cross-basis
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R
source("analysis/tests/helpers.R")
.passes <- 0L
.fails <- 0L

cat("# test_spatpomp_detected_feed\n\n")

make_fixture <- function(n_months = 6L, base_date = as.Date("2022-01-01")) {
  frame <- expand.grid(
    country = c("AAA", "BBB"),
    adm1n = 1:2,
    adm2n = 1:2,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      tu = "hmmW",
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = "REG_hmmW",
      min_seas = 1L,
      pop = 1e5
    ) %>%
    select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    as_tibble()
  # Frame pre-sorted in the pomp unit order (C-locale gid order). The builder
  # sorts internally (row-order independence is pinned separately by
  # test_spatpomp_unit_order.R); kept sorted here so this test isolates the
  # detected-feed behaviour.
  frame <- frame[order(frame$gid, method = "radix"), ]
  months <- base_date %m+% months(0:(n_months - 1L))
  panel <- tidyr::expand_grid(gid = frame$gid, t_start = months) %>%
    mutate(
      min_seas = 1L,
      country = sub("_[0-9]+_[0-9]+$", "", gid),
      t_start_seas = as.integer(yday(t_start) - 1L),
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0,
      outbreak = 0L,
      time_to_last_outbreak_cat = "[6,Inf]",
      ongoing_outbreak_lag1 = 0L,
      ongoing_outbreak_country_lag1 = 0L,
      ongoing_outbreak_tu_lag1 = 0L
    ) %>%
    as_tibble()
  list(
    frame = frame,
    panel = panel,
    months = months,
    window = list(
      start = as.character(months[1]),
      end = as.character(months[n_months] %m+% months(1L))
    )
  )
}
extract_state <- function(sm, obj, varname) {
  U <- length(obj@unit_names)
  M <- length(attr(obj, "filt_months"))
  mat <- matrix(NA_real_, nrow = M, ncol = U)
  for (tt in seq_len(M)) {
    row_t <- sm[sm$time == tt, , drop = FALSE]
    for (i in seq_len(U)) {
      r_ui <- row_t[row_t$unitname == obj@unit_names[i], , drop = FALSE]
      if (nrow(r_ui) > 0 && varname %in% names(r_ui)) {
        mat[tt, i] <- r_ui[[varname]][1L]
      }
    }
  }
  mat
}

fix <- make_fixture()
frame_tu <- fix$frame
panel_tu <- fix$panel
window <- fix$window
draws_occ_hi <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 15.0),
  basehaz_const = 0.0
)
draws_occ_mod <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -3.0),
  basehaz_const = 0.0
)
draws_dur_mod <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -3.6),
  active_terms = "(Intercept)"
)
draws_dur_never <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -30),
  active_terms = "(Intercept)"
)

build_c <- function(draws_occ, draws_dur, ...) {
  build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    draws_occ,
    window,
    fixed_draw_idx = 1L,
    draws_dur = draws_dur,
    ...
  )
}
sim <- function(obj, seed) {
  set.seed(seed)
  spatPomp::simulate(obj, format = "data.frame")
}

# --------------------------------------------------------------------------
# T1: latent default is untouched: covariate_feed omitted == "latent" ->
#     identical unit statenames (11) and identical rprocess C.
# T2: the detected feed adds exactly one appended slot (`det`, 12th).
cat("## T1-T2: layouts\n")
o_def <- build_c(draws_occ_mod, draws_dur_mod)
o_lat <- build_c(draws_occ_mod, draws_dur_mod, covariate_feed = "latent")
o_det <- build_c(
  draws_occ_mod,
  draws_dur_mod,
  covariate_feed = "detected",
  delta_fn = 0.3,
  delta_fp = 0.2
)
expect_true(
  identical(o_def@unit_statenames, .SP_UNIT_STATENAMES) &&
    identical(o_lat@unit_statenames, .SP_UNIT_STATENAMES),
  "T1: default / latent keep the 11-slot layout"
)
# The generated C lives in the object's compiled-library sources; drop the
# per-build library name line so two builds of identical C compare equal.
rp <- function(o) {
  src <- unlist(strsplit(
    paste(vapply(o@solibs, function(s) s$src, ""), collapse = "\n"),
    "\n"
  ))
  src <- src[!grepl("^/\\* (Time|Salt): ", src)] # per-build header lines
  paste(gsub("pomp_[0-9a-f]+", "pomp_X", src), collapse = "\n")
}
expect_true(
  identical(rp(o_def), rp(o_lat)),
  "T1: latent rprocess C identical to default"
)
expect_true(
  !grepl("dn\\[", rp(o_lat)) && grepl("dn\\[", rp(o_det)),
  "T1: detection block emitted only under the detected feed"
)
expect_true(
  identical(o_det@unit_statenames, c(.SP_UNIT_STATENAMES, "det")),
  "T2: detected feed appends exactly `det`"
)

# --------------------------------------------------------------------------
# T3: C kernel under the detected feed with NON-degenerate detection draws
#     (delta_fn = 0.3, delta_fp = 0.2) and a non-degenerate cessation hazard:
#     simulates, and `det` takes both values.
cat("\n## T3: C kernel under the detected feed (non-degenerate)\n")
for (seed in c(11L, 2024L)) {
  sc <- sim(o_det, seed)
  d <- extract_state(sc, o_det, "det")
  expect_true(
    length(unique(as.vector(d))) == 2L,
    sprintf("T3 seed %d: det non-degenerate", seed)
  )
}

# --------------------------------------------------------------------------
# T4: det is constant within an episode; with detection near-impossible the
#     spillover lag `o` stays 0 although every unit is active, while the
#     latent feed reports o = 1.
cat("\n## T4: semantics\n")
o_never <- build_c(
  draws_occ_hi,
  draws_dur_never,
  covariate_feed = "detected",
  delta_fn = 0.5,
  delta_fp = 1e-9
)
sc <- sim(o_never, 7L)
act <- extract_state(sc, o_never, "active")
det <- extract_state(sc, o_never, "det")
expect_true(all(act == 1), "T4: never-cease fixture keeps every unit active")
expect_true(
  all(apply(det, 2, function(x) length(unique(x)) == 1L)),
  "T4: det is constant within an episode (drawn once at onset)"
)
o_blind <- build_c(
  draws_occ_hi,
  draws_dur_never,
  covariate_feed = "detected",
  delta_fn = 1 - 1e-12,
  delta_fp = 1e-12
)
o_lat2 <- build_c(draws_occ_hi, draws_dur_never)
ob <- extract_state(sim(o_blind, 3L), o_blind, "o")
ol <- extract_state(sim(o_lat2, 3L), o_lat2, "o")
expect_true(
  all(ob == 0) && all(ol == 1),
  "T4: undetected episodes leave o = 0 (detected feed); latent feed sees o = 1"
)
ud <- extract_state(sim(o_blind, 3L), o_blind, "u_days")
expect_true(
  all(diff(ud[, 1]) > 0),
  "T4: TTLO keeps advancing when onsets go undetected"
)

# --------------------------------------------------------------------------
# T5: `_act` measurement: dunit_measure() equals the documented zero-inflated
#     NegBin activity term at every branch, including the obs_Y = NA
#     (post-onset, censored) month the old guard skipped.
cat("\n## T5: activity measurement density\n")
cb_meta_mon <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, 11),
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)
draws_mon <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  active_terms = c("(Intercept)", "mso_cb1")
)
draws_mon$fe[, "(Intercept)"] <- log(20)
draws_mon$fe[, "mso_cb1"] <- -0.15
draws_mon$phi <- rep(8, 2L)
draws_mon$offset_pop <- draws_mon$offset_logoffset <- FALSE
obs_data <- tidyr::expand_grid(gid = frame_tu$gid, t_start = fix$months) %>%
  mutate(obs_Y = 0L, cases_month = 0)
dfn <- 0.25
obj_act <- build_spatpomp_tu_c(
  frame_tu,
  panel_tu,
  draws_occ_mod,
  window,
  fixed_draw_idx = 1L,
  draws_dur = draws_dur_mod,
  draws_mon = draws_mon,
  cb_meta_mon = cb_meta_mon,
  measurement = "occurrence_duration_mon_NB_state_detect_act",
  obs_data = obs_data,
  mon_case_col = "cases_month",
  delta_fn = dfn,
  delta_fp = 0.01
)
expect_true(
  all(c("obs_Y1", "obs_CM1") %in% rownames(obj_act@data)),
  "T5: obs_CM carried"
)
eval_d <- function(Yout, obs_Y, obs_CM, mu, phi = 8) {
  U <- length(obj_act@unit_names)
  sn <- as.vector(t(outer(obj_act@unit_statenames, seq_len(U), paste0)))
  x <- matrix(0, length(sn), 1L, dimnames = list(sn, NULL))
  x["Yout1", 1L] <- Yout
  x["mu1", 1L] <- mu
  x["size_phi1", 1L] <- phi
  on <- rownames(obj_act@data)
  y <- matrix(0, length(on), 1L, dimnames = list(on, NULL))
  y["obs_Y1", 1L] <- obs_Y
  y["obs_CM1", 1L] <- obs_CM
  as.numeric(spatPomp::dunit_measure(
    obj_act,
    y = y,
    x = x,
    unit = 1L,
    time = 1,
    params = pomp::coef(obj_act),
    log = TRUE
  ))
}
onset_ll <- function(Yout, obs_Y) {
  if (is.na(obs_Y)) {
    0
  } else if (obs_Y == 1) {
    (if (Yout == 1) log(1 - dfn) else log(0.01))
  } else {
    (if (Yout == 1) log(dfn) else log(1 - 0.01))
  }
}
lnb <- function(k, mu, r = 8) dnbinom(k, size = r, mu = mu, log = TRUE)
cases <- list(
  list("NA obs_CM", 0, 0, NA, 40, function() 0),
  list("active, zero observed", 0, 0, 0, 40, function() {
    log(dfn + (1 - dfn) * exp(lnb(0, 40)))
  }),
  list("active, positive observed", 0, 0, 12, 40, function() {
    log(1 - dfn) + lnb(12, 40)
  }),
  list("inactive, positive (floor)", 0, 0, 12, 0, function() {
    lnb(12, .SP_MON_MU_FLOOR)
  }),
  list("inactive, zero observed", 0, 0, 0, 0, function() 0),
  list("post-onset obs_Y NA, active 0", 0, NA, 0, 40, function() {
    log(dfn + (1 - dfn) * exp(lnb(0, 40)))
  }),
  list("post-onset obs_Y NA, active k", 0, NA, 30, 40, function() {
    log(1 - dfn) + lnb(30, 40)
  }),
  list("onset observed + cases", 1, 1, 25, 40, function() {
    log(1 - dfn) + lnb(25, 40)
  })
)
for (cs in cases) {
  got <- eval_d(cs[[2]], cs[[3]], cs[[4]], cs[[5]])
  want <- onset_ll(cs[[2]], cs[[3]]) + cs[[6]]()
  expect_true(
    abs(got - want) < 1e-9,
    sprintf("T5 %-30s got %.6f want %.6f", cs[[1]], got, want)
  )
}
expect_true(
  eval_d(1, 0, 0, 40) < eval_d(0, 0, 0, 0),
  "T5: a still-active particle in an empty observed month is penalised"
)

# T5b: the `_actn` sibling through the C builder: case_fp_rate /
# case_fp_mean reach the compiled dunit; only the inactive branches move.
obj_act_keep <- obj_act
eps_n <- 0.05
mu0_n <- 0.3
obj_act <- build_spatpomp_tu_c(
  frame_tu,
  panel_tu,
  draws_occ_mod,
  window,
  fixed_draw_idx = 1L,
  draws_dur = draws_dur_mod,
  draws_mon = draws_mon,
  cb_meta_mon = cb_meta_mon,
  measurement = "occurrence_duration_mon_NB_state_detect_actn",
  obs_data = obs_data,
  mon_case_col = "cases_month",
  delta_fn = dfn,
  delta_fp = 0.01,
  case_fp_rate = eps_n,
  case_fp_mean = mu0_n
)
cases_n <- list(
  list("actn active, zero observed", 0, 0, 0, 40, function() {
    log(dfn + (1 - dfn) * exp(lnb(0, 40)))
  }),
  list("actn active, positive observed", 0, 0, 12, 40, function() {
    log(1 - dfn) + lnb(12, 40)
  }),
  list("actn inactive, positive", 0, 0, 12, 0, function() {
    log(eps_n) + lnb(12, mu0_n)
  }),
  list("actn inactive, zero observed", 0, 0, 0, 0, function() {
    log((1 - eps_n) + eps_n * exp(lnb(0, mu0_n)))
  }),
  list("actn post-onset NA, inactive k", 0, NA, 30, 0, function() {
    log(eps_n) + lnb(30, mu0_n)
  })
)
for (cs in cases_n) {
  got <- eval_d(cs[[2]], cs[[3]], cs[[4]], cs[[5]])
  want <- onset_ll(cs[[2]], cs[[3]]) + cs[[6]]()
  expect_true(
    abs(got - want) < 1e-9,
    sprintf("T5b %-30s got %.6f want %.6f", cs[[1]], got, want)
  )
}
obj_act <- obj_act_keep

# --------------------------------------------------------------------------
# T6: `_act` without monthly_case_obs is an error.
cat("\n## T6: config guard\n")
expect_true(
  inherits(
    try(
      resolve_filter_cfg(list(
        measurement = "occurrence_duration_mon_NB_state_detect_act"
      )),
      silent = TRUE
    ),
    "try-error"
  ),
  "T6: _act without monthly_case_obs -> error"
)
fc <- resolve_filter_cfg(list(
  measurement = "occurrence_duration_mon_NB_state_detect_act",
  monthly_case_obs = "zero_filled",
  covariate_feed = "detected"
))
expect_true(
  identical(fc$mon_case_col, "cases_month") && fc$covariate_feed == "detected",
  "T6: resolved config carries mon_case_col + feed"
)
fc0 <- resolve_filter_cfg(list())
expect_true(
  is.null(fc0$mon_case_col) &&
    fc0$covariate_feed == "latent" &&
    fc0$measurement == "occurrence_duration_mon_NB",
  "T6: defaults unchanged"
)

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
