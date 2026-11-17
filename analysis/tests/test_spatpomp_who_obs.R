# test_spatpomp_who_obs.R — HARD: WHO dashboard national monthly counts as a
#   country-level observation in the C spatPomp kernel (..._act_who).
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_who_obs.R
#
# Details:
#   Pure R/C fixtures, run from the project root; no fits, no MCMC. W1-W1b the
#   weekly / cumulative -> monthly conversion and merge, W2 config guards, W3
#   layout, W4 simulated cases, W5 the dunit, W6 bpfilter smoke.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/who_dashboard.R")
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R
source("analysis/tests/helpers.R")
.passes <- 0L
.fails <- 0L
cat("# test_spatpomp_who_obs\n\n")

WHO_MEAS <- "occurrence_duration_mon_NB_state_detect_act_who"
ACT_MEAS <- "occurrence_duration_mon_NB_state_detect_act"
is_err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# --------------------------------------------------------------------------
# W1: who_weekly_to_monthly(): day-overlap split across a month boundary,
#     partial months, NA weeks not covered, duplicate handling;
#     who_obs_table() keeps complete months in the window; who_block_share()
#     sums to 1 per country; who_coverage_fraction() clamps and errors on a
#     missing ISO3.
cat("## W1: weekly -> monthly and per-unit constants\n")
wk <- tibble(
  iso_3_code = c("AAA", "AAA", "AAA", "BBB", "BBB"),
  # AAA: week of Mon 2023-01-30 straddles Jan/Feb (2 Jan days, 5 Feb days)
  date_wk = as.Date(c(
    "2023-01-30",
    "2023-02-06",
    "2023-02-13",
    "2023-02-06",
    "2023-02-13"
  )),
  cases = c(70, 14, NA, 7, 7)
)
m <- who_weekly_to_monthly(wk)
a_jan <- m[m$iso3 == "AAA" & m$t_start == as.Date("2023-01-01"), ]
a_feb <- m[m$iso3 == "AAA" & m$t_start == as.Date("2023-02-01"), ]
expect_true(
  abs(a_jan$cases - 20) < 1e-9 && a_jan$n_days == 2L && a_jan$partial,
  "W1: straddling week puts 2/7 in January (partial month)"
)
expect_true(
  abs(a_feb$cases - 64) < 1e-9 && a_feb$n_days == 12L && a_feb$partial,
  "W1: 5/7 + a full week in February; the NA week is not covered"
)
expect_true(
  is_err(who_weekly_to_monthly(rbind(
    wk,
    tibble(iso_3_code = "BBB", date_wk = as.Date("2023-02-06"), cases = 9)
  ))),
  "W1: conflicting duplicate week is an error"
)
expect_true(
  nrow(who_weekly_to_monthly(rbind(wk, wk[1, ]))) == nrow(m),
  "W1: exact duplicate collapsed"
)
full <- tibble(
  iso_3_code = "CCC",
  date_wk = as.Date("2023-02-27") + 7 * (0:5),
  cases = 7
)
mf <- who_weekly_to_monthly(full)
expect_true(
  identical(mf$partial[mf$t_start == as.Date("2023-03-01")], FALSE) &&
    abs(mf$cases[mf$t_start == as.Date("2023-03-01")] - 31) < 1e-9,
  "W1: fully covered month is complete with cases = days"
)
ot <- who_obs_table(mf, "2023-03-01", "2023-04-01")
expect_true(
  nrow(ot) == 1L && ot$country == "CCC",
  "W1: who_obs_table keeps complete months in window"
)
fr <- tibble(country = c("AAA", "AAA", "BBB"), pop = c(1, 3, 0))
expect_true(
  isTRUE(all.equal(who_block_share(fr, "pop"), c(.25, .75, 1))) &&
    isTRUE(all.equal(who_block_share(fr, "equal"), c(.5, .5, 1))),
  "W1: who_block_share sums to 1 per country (zero-pop country -> equal)"
)
wpp <- data.frame(iso3 = c("AAA", "BBB"), year = 2020L, pop_total = c(8, 1))
fr2 <- tibble(country = c("AAA", "AAA", "BBB"), pop = c(1, 3, 2))
cv <- who_coverage_fraction(fr2, wpp, 2020L)
expect_true(
  isTRUE(all.equal(cv$who_f, c(0.5, 1))),
  "W1: who_f = covered / national, clamped at 1"
)
expect_true(
  is_err(who_coverage_fraction(tibble(country = "ZZZ", pop = 1), wpp, 2020L)),
  "W1: country missing from WPP is an error"
)

# --------------------------------------------------------------------------
# W1b: who_cumulative_to_monthly(): year-to-date differencing with the
#      1 January reset, day-overlap split of a cut interval across a month
#      boundary, unchanged cumulative -> zero only for the allowlist or an
#      explicit last-28-days 0, downward revision dropped with a warning;
#      who_merge_monthly(): weekly rows (complete and partial) win, cumulative
#      rows fill only absent months, `source` labels; the weekly path through
#      the shared .who_intervals_to_monthly() reproduces the on-disk table.
cat("\n## W1b: cumulative snapshots -> monthly, and the gap-fill merge\n")
# DDD: year-to-date cuts 2025-12-28 (100), 2026-01-25 (50), 2026-02-22 (50),
# 2026-03-29 (80): Jan 1-25 = 50; Jan 26-Feb 22 = 0 (unchanged); Feb 23-Mar 29
# = 30. EEE: an explicit last-28-days 0 at the second cut. FFF: revised down.
snap <- tibble(
  iso3 = c("DDD", "DDD", "DDD", "DDD", "EEE", "EEE", "FFF", "FFF"),
  table_as_of = as.Date(c(
    "2025-12-28",
    "2026-01-25",
    "2026-02-22",
    "2026-03-29",
    "2026-01-25",
    "2026-02-22",
    "2026-01-25",
    "2026-02-22"
  )),
  cum_cases = c(100, 50, 50, 80, 10, 10, 40, 30),
  last28_cases = c(NA, NA, NA, 30, NA, 0, NA, NA)
)
cm0 <- suppressWarnings(who_cumulative_to_monthly(snap))
d_jan <- cm0[cm0$iso3 == "DDD" & cm0$t_start == as.Date("2026-01-01"), ]
expect_true(
  nrow(d_jan) == 1L &&
    abs(d_jan$cases - 50) < 1e-9 &&
    d_jan$n_days == 25L &&
    d_jan$partial,
  "W1b: first cut of the year runs from 1 January (year-to-date reset), 25 covered days"
)
d_feb0 <- cm0[cm0$iso3 == "DDD" & cm0$t_start == as.Date("2026-02-01"), ]
expect_true(
  nrow(d_feb0) == 1L &&
    d_feb0$n_days == 6L &&
    d_feb0$partial &&
    abs(d_feb0$cases - 30 * 6 / 35) < 1e-9,
  "W1b: unchanged cumulative without allowlist -> Jan 26-Feb 22 dropped; February keeps only the next interval's 6 days"
)
d_mar <- cm0[cm0$iso3 == "DDD" & cm0$t_start == as.Date("2026-03-01"), ]
expect_true(
  abs(d_mar$cases - 30 * 29 / 35) < 1e-9 &&
    d_mar$n_days == 29L &&
    d_mar$partial,
  "W1b: Feb 23-Mar 29 interval (35 days) puts 29/35 of its 30 cases in March"
)
expect_true(
  !any(cm0$iso3 == "DDD" & cm0$t_start < as.Date("2026-01-01")),
  "W1b: a lone December year-to-date cut (362 days) exceeds max_interval_days and is dropped"
)
cm_long <- suppressWarnings(who_cumulative_to_monthly(
  snap,
  max_interval_days = 400L
))
d_dec <- cm_long[
  cm_long$iso3 == "DDD" & cm_long$t_start == as.Date("2025-12-01"),
]
expect_true(
  nrow(d_dec) == 1L &&
    abs(d_dec$cases - 100 * 28 / 362) < 1e-9 &&
    d_dec$partial,
  "W1b: with a long max_interval_days the lone cut is spread over 1 Jan-cut (December partial)"
)
e_feb <- cm0[cm0$iso3 == "EEE" & cm0$t_start == as.Date("2026-02-01"), ]
expect_true(
  nrow(e_feb) == 1L && e_feb$cases == 0 && e_feb$n_days == 22L,
  "W1b: explicit last-28-days 0 -> reported zero interval"
)
expect_true(
  !any(cm0$iso3 == "FFF" & cm0$t_start == as.Date("2026-02-01")) &&
    any(cm0$iso3 == "FFF" & cm0$t_start == as.Date("2026-01-01")),
  "W1b: downward revision drops the interval, keeps the earlier one"
)
expect_true(
  inherits(
    tryCatch(who_cumulative_to_monthly(snap), warning = function(w) w),
    "warning"
  ),
  "W1b: downward revision warns"
)
cm1 <- suppressWarnings(who_cumulative_to_monthly(snap, zero_iso3 = "DDD"))
d_feb1 <- cm1[cm1$iso3 == "DDD" & cm1$t_start == as.Date("2026-02-01"), ]
expect_true(
  nrow(d_feb1) == 1L &&
    d_feb1$n_days == 28L &&
    !d_feb1$partial &&
    abs(d_feb1$cases - 30 * 6 / 35) < 1e-9,
  "W1b: allowlisted country: unchanged cumulative -> zero over Jan 26-Feb 22, February complete"
)
d_jan1 <- cm1[cm1$iso3 == "DDD" & cm1$t_start == as.Date("2026-01-01"), ]
expect_true(
  abs(d_jan1$cases - 50) < 1e-9 && d_jan1$n_days == 31L && !d_jan1$partial,
  "W1b: allowlisted zero interval completes January too"
)
expect_true(
  is_err(who_cumulative_to_monthly(rbind(
    snap,
    tibble(
      iso3 = "DDD",
      table_as_of = as.Date("2026-01-25"),
      cum_cases = 51,
      last28_cases = NA
    )
  ))),
  "W1b: conflicting duplicate cut is an error"
)
# merge: weekly has DDD Jan (partial) and Mar (complete) -> both win; Feb filled.
wk_m <- tibble(
  iso3 = "DDD",
  t_start = as.Date(c("2026-01-01", "2026-03-01")),
  cases = c(5, 7),
  n_days = c(10L, 31L),
  partial = c(TRUE, FALSE)
)
mg <- who_merge_monthly(wk_m, cm1)
g <- mg[mg$iso3 == "DDD" & mg$t_start >= as.Date("2026-01-01"), ]
expect_true(
  identical(
    g$source,
    c("dashboard_weekly", "sitrep_cumulative", "dashboard_weekly")
  ) &&
    abs(g$cases[1] - 5) < 1e-9 &&
    g$partial[1] &&
    abs(g$cases[3] - 7) < 1e-9,
  "W1b: weekly rows (partial January, March) win; February filled from the cumulative"
)
expect_true(
  all(who_merge_monthly(wk_m, NULL)$source == "dashboard_weekly") &&
    nrow(who_merge_monthly(wk_m, NULL)) == 2L,
  "W1b: no cumulative table -> weekly rows only, source = dashboard_weekly"
)
ot2 <- who_obs_table(mg, "2026-01-01", "2026-04-01")
ot2 <- ot2[ot2$country == "DDD", ]
expect_true(
  nrow(ot2) == 2L && all(ot2$t_start == as.Date(c("2026-02-01", "2026-03-01"))),
  "W1b: who_obs_table keeps the filled complete month and ignores `source`"
)
who_rds <- "generated_data/who_dashboard_monthly_africa.rds"
who_wk <- "data/who_annual_report_africa/who_dashboard_weekly_africa.rds"
if (file.exists(who_rds) && file.exists(who_wk)) {
  disk <- readRDS(who_rds)
  rebuilt <- who_weekly_to_monthly(readRDS(who_wk))
  if ("source" %in% names(disk)) {
    disk <- disk[
      disk$source == "dashboard_weekly",
      setdiff(names(disk), "source")
    ]
  }
  expect_true(
    isTRUE(all.equal(as.data.frame(disk), as.data.frame(rebuilt))),
    "W1b: shared interval helper reproduces the on-disk weekly rows exactly"
  )
} else {
  cat("  (skip) on-disk WHO tables not present\n")
}

# --------------------------------------------------------------------------
# W2: measurement <-> who_obs.enable both ways, contribution shares refused,
#     ADM1 blocks refused unless allowed.
cat("\n## W2: config guards\n")
base <- list(monthly_case_obs = "zero_filled", covariate_feed = "detected")
expect_true(
  is_err(resolve_filter_cfg(c(base, list(measurement = WHO_MEAS)))),
  "W2: _who without who_obs.enable -> error"
)
expect_true(
  is_err(resolve_filter_cfg(c(
    base,
    list(measurement = ACT_MEAS, who_obs = list(enable = TRUE))
  ))),
  "W2: who_obs.enable under another measurement -> error"
)
expect_true(
  is_err(resolve_filter_cfg(c(
    base,
    list(
      measurement = WHO_MEAS,
      who_obs = list(enable = TRUE, block_share = "contribution")
    )
  ))),
  "W2: contribution shares refused"
)
fc <- resolve_filter_cfg(c(
  base,
  list(measurement = WHO_MEAS, who_obs = list(enable = TRUE, phi = 3))
))
expect_true(
  fc$who_obs$enable && fc$who_obs$phi == 3 && fc$who_obs$block_share == "pop",
  "W2: who_obs resolved (phi from config, defaults filled)"
)
expect_true(
  !resolve_filter_cfg(c(base, list(measurement = ACT_MEAS)))$who_obs$enable,
  "W2: default who_obs disabled"
)
expect_true(
  is_err(.check_who_obs_blocks(fc$who_obs, "adm1")),
  "W2: ADM1 blocks refused by default"
)
expect_true(
  isTRUE(.check_who_obs_blocks(fc$who_obs, "country")),
  "W2: country blocks accepted"
)
w_adm1 <- tryCatch(
  .check_who_obs_blocks(
    modifyList(fc$who_obs, list(allow_adm1_blocks = TRUE)),
    "adm1"
  ),
  warning = function(w) "warned"
)
expect_true(identical(w_adm1, "warned"), "W2: allow_adm1_blocks warns")

# --------------------------------------------------------------------------
# Fixture: 2 countries x 4 units (as test_spatpomp_detected_feed.R)
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
    pop = 1e5 * adm1n
  ) %>%
  select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
  as_tibble()
frame <- frame[order(frame$gid, method = "radix"), ]
months <- as.Date("2022-01-01") %m+% months(0:5)
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
window <- list(
  start = as.character(months[1]),
  end = as.character(months[6] %m+% months(1L))
)
draws_occ <- make_draws_bundle_synthetic(
  frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -1.0),
  basehaz_const = 0.0
)
draws_dur <- make_draws_bundle_synthetic(
  frame,
  S = 2L,
  fe_overrides = c("(Intercept)" = -2.0),
  active_terms = "(Intercept)"
)
cb_meta_mon <- list(
  mso = list(
    single_lag = 0L,
    range = c(0, 11),
    onebasis_attrs = list(fun = "lin"),
    cb_col_names = "mso_cb1"
  )
)
draws_mon <- make_draws_bundle_synthetic(
  frame,
  S = 2L,
  active_terms = c("(Intercept)", "mso_cb1")
)
draws_mon$fe[, "(Intercept)"] <- log(20)
draws_mon$fe[, "mso_cb1"] <- -0.15
draws_mon$phi <- rep(8, 2L)
draws_mon$offset_pop <- draws_mon$offset_logoffset <- FALSE
obs_data <- tidyr::expand_grid(gid = frame$gid, t_start = months) %>%
  mutate(obs_Y = 0L, cases_month = 0)
who_tbl <- tibble(
  country = c("AAA", "AAA", "BBB"),
  t_start = months[c(2, 3, 2)],
  cases_who = c(500, 0, 80)
)
f_ctry <- c(AAA = 0.5, BBB = 0.8)
who_in <- function(share = "pop", phi = 4, mu_floor = 1, scale = 1) {
  list(
    tbl = who_tbl,
    w = who_block_share(frame, share),
    f = unname(f_ctry[frame$country]),
    phi = phi,
    mu_floor = mu_floor,
    scale = scale
  )
}
build <- function(meas, who = NULL, dfn = 0.25) {
  build_spatpomp_tu_c(
    frame,
    panel,
    draws_occ,
    window,
    fixed_draw_idx = 1L,
    draws_dur = draws_dur,
    draws_mon = draws_mon,
    cb_meta_mon = cb_meta_mon,
    measurement = meas,
    obs_data = obs_data,
    mon_case_col = "cases_month",
    covariate_feed = "detected",
    delta_fn = dfn,
    delta_fp = 0.01,
    who_obs = who
  )
}
rp <- function(o) {
  src <- unlist(strsplit(
    paste(vapply(o@solibs, function(s) s$src, ""), collapse = "\n"),
    "\n"
  ))
  src <- src[!grepl("^/\\* (Time|Salt): ", src)]
  paste(gsub("pomp_[0-9a-f]+", "pomp_X", src), collapse = "\n")
}

# --------------------------------------------------------------------------
# W3: `_who` appends exactly cm, ccm; the `_act` build is unchanged
#     (statenames and rprocess C); obs_W replicated per country.
cat("\n## W3: layout\n")
o_act <- build(ACT_MEAS)
o_who <- build(WHO_MEAS, who_in())
expect_true(
  identical(o_who@unit_statenames, c(.SP_UNIT_STATENAMES, "det", "cm", "ccm")),
  "W3: _who appends cm, ccm after det"
)
expect_true(
  identical(o_act@unit_statenames, c(.SP_UNIT_STATENAMES, "det")),
  "W3: _act layout unchanged"
)
expect_true(
  !grepl("rnbinom_mu", rp(o_act)) && grepl("rnbinom_mu", rp(o_who)),
  "W3: simulated-cases block emitted only for _who"
)
expect_true(
  is_err(build(WHO_MEAS)) && is_err(build(ACT_MEAS, who_in())),
  "W3: measurement / who_obs mismatch refused by the builder"
)
dat <- o_who@data
U <- length(o_who@unit_names)
w_row <- function(u) dat[paste0("obs_W", u), ]
aaa <- which(startsWith(o_who@unit_names, "AAA"))
bbb <- which(startsWith(o_who@unit_names, "BBB"))
expect_true(
  all(vapply(
    aaa,
    function(u) identical(unname(w_row(u)[2:3]), c(500, 0)),
    logical(1)
  )) &&
    all(vapply(
      bbb,
      function(u) unname(w_row(u)[2]) == 80 && is.na(w_row(u)[3]),
      logical(1)
    )) &&
    all(is.na(dat[paste0("obs_W", seq_len(U)), c(1, 4:6)])),
  "W3: obs_W replicated on every unit of the country, NA elsewhere"
)

# --------------------------------------------------------------------------
# W4: cm == 0 wherever mu == 0, cm > 0 somewhere, and ccm is the country sum
#     of cm on every unit.
cat("\n## W4: simulated cases\n")
sim <- function(obj, seed) {
  set.seed(seed)
  spatPomp::simulate(obj, format = "data.frame")
}
st <- function(sm, v) {
  M <- max(sm$time)
  mat <- matrix(NA_real_, M, U)
  for (i in seq_len(U)) {
    r <- sm[sm$unitname == o_who@unit_names[i], ]
    mat[r$time, i] <- r[[v]]
  }
  mat
}
sm <- sim(o_who, 5L)
mu <- st(sm, "mu")
cm <- st(sm, "cm")
ccm <- st(sm, "ccm")
expect_true(
  all(cm[mu == 0] == 0) && any(cm > 0),
  "W4: cm = 0 where mu = 0, positive somewhere"
)
ct <- sub("_.*", "", o_who@unit_names)
ccm_ref <- t(apply(cm, 1, function(x) ave(x, ct, FUN = sum)))
expect_true(
  isTRUE(all.equal(ccm, ccm_ref)),
  "W4: ccm is the country sum of cm on every unit"
)
expect_true(all(cm == round(cm)), "W4: cm are counts")

# --------------------------------------------------------------------------
# W5: summed over a country's units = `_act` sum + log NB(W; ccm * scale / f,
#     phi) for pop and equal shares; obs_W NA -> `_act` density; ccm == 0 ->
#     mu_floor.
cat("\n## W5: dunit\n")
dfn <- 0.25
eval_u <- function(obj, u, Yout, obs_Y, obs_CM, mu, obs_W, ccm, phi_m = 8) {
  sn <- as.vector(t(outer(obj@unit_statenames, seq_len(U), paste0)))
  x <- matrix(0, length(sn), 1L, dimnames = list(sn, NULL))
  x[paste0("Yout", u), 1L] <- Yout
  x[paste0("mu", u), 1L] <- mu
  x[paste0("size_phi", u), 1L] <- phi_m
  if (paste0("ccm", u) %in% sn) {
    x[paste0("ccm", u), 1L] <- ccm
  }
  on <- rownames(obj@data)
  y <- matrix(NA_real_, length(on), 1L, dimnames = list(on, NULL))
  y[paste0("obs_Y", u), 1L] <- obs_Y
  y[paste0("obs_CM", u), 1L] <- obs_CM
  if (paste0("obs_W", u) %in% on) {
    y[paste0("obs_W", u), 1L] <- obs_W
  }
  as.numeric(spatPomp::dunit_measure(
    obj,
    y = y,
    x = x,
    unit = u,
    time = 2,
    params = pomp::coef(obj),
    log = TRUE
  ))
}
for (share in c("pop", "equal")) {
  o <- build(
    WHO_MEAS,
    who_in(share = share, phi = 4, mu_floor = 2, scale = 1.5)
  )
  for (cfg in list(
    list(W = 500, ccm = 120),
    list(W = 500, ccm = 0),
    list(W = NA, ccm = 120)
  )) {
    tot_who <- sum(vapply(
      aaa,
      function(u) eval_u(o, u, 1, 1, 25, 40, cfg$W, cfg$ccm),
      numeric(1)
    ))
    tot_act <- sum(vapply(
      aaa,
      function(u) eval_u(o_act, u, 1, 1, 25, 40, NA, 0),
      numeric(1)
    ))
    m_nat <- if (cfg$ccm > 0) cfg$ccm * 1.5 / f_ctry[["AAA"]] else 2
    want <- tot_act +
      if (is.na(cfg$W)) 0 else dnbinom(cfg$W, size = 4, mu = m_nat, log = TRUE)
    expect_true(
      abs(tot_who - want) < 1e-8,
      sprintf(
        "W5 [%s] W=%s ccm=%s: country sum = _act + NB term (got %.6f want %.6f)",
        share,
        cfg$W,
        cfg$ccm,
        tot_who,
        want
      )
    )
  }
}

# --------------------------------------------------------------------------
# W6: bpfilter smoke with country blocks; .sp_strip_scratch_slots() leaves the
#     sv11/sv12 layout.
cat("\n## W6: bpfilter smoke + strip\n")
bl <- make_block_list(o_who, frame, granularity = "country")
res <- run_one_bpfilter_replicate(o_who, Np = 50L, block_list = bl, seed = 1L)
expect_true(
  is.finite(res$logLik),
  sprintf("W6: bpfilter logLik finite (%.3f)", res$logLik)
)
rn <- dimnames(res$filter_traj)[[1]]
expect_true(
  any(startsWith(rn, "ccm")),
  "W6: filter_traj carries the scratch slots in memory"
)
ft <- .sp_strip_scratch_slots(res$filter_traj)
expect_true(
  !any(grepl("^c?cm[0-9]+$", dimnames(ft)[[1]])) &&
    dim(ft)[1] == length(.sp_unit_statenames("detected")) * U,
  "W6: stripped trajectory has the sv12 (detected) layout"
)
expect_true(
  identical(.sp_strip_scratch_slots(ft), ft),
  "W6: strip is the identity without scratch slots"
)
res_a <- run_one_bpfilter_replicate(o_act, Np = 50L, block_list = bl, seed = 1L)
expect_true(
  is.finite(res_a$logLik) && res$logLik < res_a$logLik,
  sprintf(
    "W6: WHO term lowers the logLik (%.2f < %.2f)",
    res$logLik,
    res_a$logLik
  )
)

# --------------------------------------------------------------------------
# Per-country estimated scale (who_obs$scale_u) and the
# back-compat golden of the fixed-scale build.
cat("\n## W7: per-unit estimated scale\n")
sc_ctry <- c(AAA = 2.5, BBB = 0.4)
who_u <- c(
  who_in(phi = 4, mu_floor = 2),
  list(scale_u = unname(sc_ctry[frame$country]))
)
o_u <- build(WHO_MEAS, who_u)
for (ct_i in c("AAA", "BBB")) {
  us <- if (ct_i == "AAA") aaa else bbb
  tot_who <- sum(vapply(
    us,
    function(u) eval_u(o_u, u, 1, 1, 25, 40, 300, 90),
    numeric(1)
  ))
  tot_act <- sum(vapply(
    us,
    function(u) eval_u(o_act, u, 1, 1, 25, 40, NA, 0),
    numeric(1)
  ))
  want <- tot_act +
    dnbinom(
      300,
      size = 4,
      mu = 90 * sc_ctry[[ct_i]] / f_ctry[[ct_i]],
      log = TRUE
    )
  expect_true(
    abs(tot_who - want) < 1e-8,
    sprintf(
      "W7 [%s]: country sum = _act + log NB(W; ccm x scale_c / f_c)",
      ct_i
    )
  )
}
expect_true(
  identical(
    attr(o_u, "who_obs")$scale_u,
    unname(sc_ctry[sub("_.*", "", o_u@unit_names)])
  ),
  "W7: the object records scale_u in unit order"
)
expect_true(
  is_err(build(WHO_MEAS, c(who_in(), list(scale_u = c(1, 2))))) &&
    is_err(build(WHO_MEAS, c(who_in(), list(scale_u = rep(-1, nrow(frame)))))),
  "W7: scale_u of the wrong length or <= 0 refused"
)
expect_true(
  isTRUE(all.equal(who_multiplier(f_ctry, sc_ctry), sc_ctry / f_ctry)),
  "W7: who_multiplier() = scale / who_f"
)

cat("\n## W8: fixed-scale build = the pre-item-10 golden (c032a12)\n")
gold <- readRDS("analysis/tests/golden/who_obs_csnippet_golden.rds")
for (sc in c(1, 1.5)) {
  g <- gold[[sprintf("scale_%s", sc)]]
  o <- build(WHO_MEAS, who_in(scale = sc))
  res_g <- run_one_bpfilter_replicate(
    o,
    Np = 50L,
    block_list = make_block_list(o, frame, granularity = "country"),
    seed = 1L
  )
  expect_true(
    identical(rp(o), g$csnippet),
    sprintf("W8 scale %s: Csnippet text identical to the golden", sc)
  )
  expect_true(
    identical(res_g$logLik, g$logLik),
    sprintf("W8 scale %s: seeded bpfilter logLik identical", sc)
  )
  expect_true(
    identical(attr(o, "who_obs"), g$who_attr),
    sprintf("W8 scale %s: who_obs attribute identical", sc)
  )
}

cat("\n## W9: config keys (scale_mode / scale_est / db_overlap)\n")
r <- .resolve_who_obs_cfg(
  list(enable = TRUE, scale_mode = "estimated"),
  WHO_MEAS
)
expect_true(
  identical(r$scale, 1) && is.null(list(scale_mode = "estimated")[["scale"]]),
  "W9: `scale` absent -> 1 (no partial match of scale_mode / scale_est)"
)
expect_true(
  identical(r$scale_est$fallback, "pooled") &&
    identical(r$scale_est$clamp, c(0.1, 10)) &&
    identical(r$db_overlap, "score_both"),
  "W9: scale_est / db_overlap defaults"
)
expect_true(
  is_err(.resolve_who_obs_cfg(
    list(enable = TRUE, scale_mode = "x"),
    WHO_MEAS
  )) &&
    is_err(.resolve_who_obs_cfg(
      list(enable = TRUE, db_overlap = "x"),
      WHO_MEAS
    )) &&
    is_err(.resolve_who_obs_cfg(
      list(
        enable = TRUE,
        scale_mode = "estimated",
        scale_est = list(clamp = c(1, 0.5))
      ),
      WHO_MEAS
    )),
  "W9: bad scale_mode / db_overlap / clamp refused"
)
expect_true(
  is.null(who_in()[["scale_u"]]),
  "W9: who_obs[['scale_u']] is NULL when only `scale` is set"
)

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
