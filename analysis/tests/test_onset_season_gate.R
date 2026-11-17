# test_onset_season_gate.R — HARD: simulator.onset_eligibility, the season
#   onset gate in the v1 and spatPomp C kernels.
#
# Usage:
#   Rscript analysis/tests/test_onset_season_gate.R
#
# Details:
#   Synthetic, no MCMC, no data; run from the project root. Covers the gate
#   in simulate_one_tu() and build_spatpomp_tu_c(): helpers (H), the frame
#   refusal (F), the covariate table and Csnippet (K), degenerate (D) and
#   moderate-hazard (M) kernels, v1 and C alike, and cache identity (P).

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
  library(tidyr)
})
source("analysis/utils/smoothing_core.R") # base, simulator_core, spatpomp_core{,_c}
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
err_of <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
cat("# test_onset_season_gate\n\n")

# ===========================================================================
# H  helpers: onset_eligibility_of() / check_onset_eligibility(),
#    season_month_index(), season_onset_open() on hand-worked u_days /
#    seas_k cases (onset in the season's first month blocks the rest of it;
#    an onset in the previous season's last month does not block month 0;
#    sentinel and NA eligible).
cat("## H: helpers\n")
expect(
  identical(onset_eligibility_of(NULL), "inactive") &&
    identical(onset_eligibility_of(list()), "inactive") &&
    identical(
      onset_eligibility_of(list(onset_eligibility = "season")),
      "season"
    ),
  "H1 onset_eligibility_of: default inactive, season accepted"
)
expect(
  all(vapply(
    list("Season", "", NA_character_, c("season", "inactive"), TRUE, 1),
    function(v) {
      !is.null(err_of(onset_eligibility_of(list(onset_eligibility = v))))
    },
    logical(1)
  )),
  "H2 onset_eligibility_of refuses anything else"
)
mo <- as.Date("2024-01-01") %m+% months(0:13)
sk <- season_month_index(c(1L, 4L, 12L), mo)
expect(
  identical(sk[1, 1:13], c(0:11, 0L)) &&
    identical(sk[2, 1:5], c(9L, 10L, 11L, 0L, 1L)) &&
    identical(sk[3, 1:2], c(1L, 2L)) &&
    identical(dim(sk), c(3L, 14L)),
  "H3 season_month_index = (month - min_seas) mod 12"
)
step <- 365.25 / 12
# onset at season month k0 leaves u_days = 0 going into k0 + 1, then + step per month
expect(
  identical(
    season_onset_open(c(0, step, 10 * step), c(1L, 2L, 11L)),
    c(FALSE, FALSE, FALSE)
  ),
  "H4 onset earlier in this season blocks the rest of the season"
)
expect(
  identical(season_onset_open(0, 0L), TRUE),
  "H5 onset in the previous season's last month does not block month 0"
)
expect(
  identical(
    season_onset_open(c(11 * step, 3 * step), c(11L, 3L)),
    c(TRUE, TRUE)
  ),
  "H6 onset in the previous season (n = k) is eligible"
)
expect(
  identical(
    season_onset_open(c(.TTLO_NONE_DAYS, NA, 6 * 365.25), c(11L, 5L, 11L)),
    c(TRUE, TRUE, TRUE)
  ),
  "H7 no-prior sentinel, NA and the [6,Inf] warm start are eligible"
)
# exact sums: 12 advance_ttlo steps from 0 land on 12 * step exactly
u <- 0
for (i in 1:12) {
  u <- advance_ttlo(u, 0L)
}
expect(
  identical(u, 12 * step) && identical(floor(u / 30.4375 + 0.5), 12),
  "H8 u_days is an exact multiple of 30.4375 (no rounding drift)"
)

# ===========================================================================
# F  check_onset_eligibility_frame(): no-op under inactive; refuses a frame
#    whose min_seas was defaulted (build_location_frame()'s
#    "min_seas_defaulted" attribute) under season, counting only gids still
#    in the frame; build_location_frame() sets the attribute only when it
#    defaulted a gid.
cat("\n## F: frame refusal\n")
mk_frame <- function(ms) {
  tibble(
    gid = paste0("AAA_1_", seq_along(ms)),
    country = "AAA",
    AFRO_region = "R",
    adm1_unit = "AAA_1",
    tu = "T",
    min_seas = ms,
    pop = 1e5
  )
}
f_ok <- mk_frame(c(1L, 4L))
f_df <- f_ok
attr(f_df, "min_seas_defaulted") <- "AAA_1_1"
expect(
  is.null(err_of(check_onset_eligibility_frame(f_df, "inactive"))),
  "F1 inactive: no check"
)
expect(
  is.null(err_of(check_onset_eligibility_frame(f_ok, "season"))),
  "F2 season: own min_seas passes"
)
expect(
  !is.null(err_of(check_onset_eligibility_frame(f_df, "season"))),
  "F3 season: defaulted min_seas refused"
)
expect(
  !is.null(err_of(check_onset_eligibility_frame(
    mk_frame(c(1L, NA)),
    "season"
  ))),
  "F4 season: NA min_seas refused"
)
expect(
  is.null(err_of(check_onset_eligibility_frame(
    f_ok,
    "season",
    defaulted = "AAA_1_9"
  ))) &&
    !is.null(err_of(check_onset_eligibility_frame(
      f_ok,
      "season",
      defaulted = "AAA_1_2"
    ))),
  "F4b season: only defaulted gids still in the (final) frame count"
)
pan <- tibble(
  gid = c("AAA_1_1", "AAA_1_2"),
  country = "AAA",
  AFRO_region = "R",
  tu = "T",
  pop = 1e5
)
fr_def <- suppressMessages(build_location_frame(pan))
fr_own <- suppressMessages(build_location_frame(mutate(
  pan,
  min_seas = c(3L, 5L)
)))
expect(
  identical(attr(fr_def, "min_seas_defaulted"), c("AAA_1_1", "AAA_1_2")) &&
    is.null(attr(fr_own, "min_seas_defaulted")),
  "F5 build_location_frame: min_seas_defaulted set only when a gid was defaulted"
)

# ===========================================================================
# Fixture: 4 units (2 countries x 2 ADM1, one ADM2 each), min_seas 1/4/7/11,
# 15 months from 2024-01, warm start [6,Inf] (no recent onset).
min_seas <- c(1L, 4L, 7L, 11L)
frame_tu <- tidyr::expand_grid(country = c("AAA", "BBB"), adm1n = 1:2) %>%
  mutate(
    gid = sprintf("%s_%d_1", country, adm1n),
    adm1_unit = sprintf("%s_%d", country, adm1n),
    AFRO_region = "REG_W",
    tu = "hmmW",
    min_seas = min_seas,
    pop = 1e5
  ) %>%
  select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
  as_tibble()
n_m <- 15L
months <- as.Date("2024-01-01") %m+% months(0:(n_m - 1L))
pre <- as.Date("2023-12-01")
panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = c(pre, months)) %>%
  left_join(select(frame_tu, gid, min_seas), by = "gid") %>%
  mutate(
    country = sub("_[0-9]+_[0-9]+$", "", gid),
    t_start_seas = as.integer(
      yday(t_start %m+% period(-(min_seas - 1L), "months")) - 1L
    ),
    mean_precip = 0,
    dev_precip = 0,
    mean_temp_daymean = 0,
    dev_temp_daymean = 0,
    outbreak = 0L,
    ongoing_outbreak_lag1 = 0L,
    ongoing_outbreak_country_lag1 = 0L,
    ongoing_outbreak_tu_lag1 = 0L,
    time_to_last_outbreak_cat = "[6,Inf]"
  ) %>%
  as_tibble()
horizon <- list(start = "2024-01-01", end = "2025-03-31")
window <- list(start = "2024-01-01", end = "2025-04-01")
sk_fix <- season_month_index(frame_tu$min_seas, months) # units x months
season_id <- function(ms, mon) {
  # season-year id per (unit, month)
  a <- as.integer(format(mon, "%Y")) * 12L + as.integer(format(mon, "%m")) - 1L
  (a - (ms - 1L)) %/% 12L
}

d_occ_hi <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 15),
  active_terms = "(Intercept)"
)
d_occ_mod <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -4),
  active_terms = "(Intercept)"
)
d_dur_end <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 30),
  active_terms = "(Intercept)"
)

# ===========================================================================
# K  covariate table + Csnippet: seas_k only under season, values
#    (month - min_seas) mod 12 on the one-month-shifted grid; the emitted
#    onset line carries the gate; season without a duration fit refused.
cat("\n## K: covariate table and Csnippet\n")
cv0 <- build_covar_for_tu(panel_tu, frame_tu, months)
cv1 <- build_covar_for_tu(panel_tu, frame_tu, months, season_k = TRUE)
expect(
  !"seas_k" %in% names(cv0) && identical(cv0, select(cv1, -seas_k)),
  "K1 season_k = FALSE: table unchanged; TRUE only appends seas_k"
)
# time t carries month t+1's covariates (one-month shift); time M duplicates month M
want_k <- sk_fix[cbind(match(cv1$unit, frame_tu$gid), pmin(cv1$time + 1L, n_m))]
expect(
  identical(cv1$seas_k, as.numeric(want_k)),
  "K2 seas_k = (month - min_seas) mod 12 on the shifted grid"
)
o_inact <- build_spatpomp_tu_c(
  frame_tu,
  panel_tu,
  d_occ_hi,
  window,
  fixed_draw_idx = 1L,
  draws_dur = d_dur_end
)
o_seas <- build_spatpomp_tu_c(
  frame_tu,
  panel_tu,
  d_occ_hi,
  window,
  fixed_draw_idx = 1L,
  draws_dur = d_dur_end,
  onset_eligibility = "season"
)
src_of <- function(o) paste(vapply(o@solibs, `[[`, "", "src"), collapse = "\n")
gate_c <- "if (av[u] < 0.5 && (ISNAN(ud[u]) || floor(ud[u] / 30.4375 + 0.5) >= sk[u]))"
expect(
  !grepl("sk[u]", src_of(o_inact), fixed = TRUE) &&
    grepl(gate_c, src_of(o_seas), fixed = TRUE),
  "K3 the emitted onset line carries the gate only under season"
)
expect(
  is.null(attr(o_inact, "onset_eligibility")) &&
    identical(attr(o_seas, "onset_eligibility"), "season"),
  "K4 object attribute only under season"
)
expect(
  !is.null(err_of(build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    d_occ_hi,
    window,
    fixed_draw_idx = 1L,
    onset_eligibility = "season"
  ))) &&
    !is.null(err_of(simulate_one_tu(
      tu = "hmmW",
      frame_tu = frame_tu,
      panel_tu = panel_tu,
      draws_occ = d_occ_hi,
      replicate_ids = 1L,
      seed = 1L,
      output_part = tempfile(),
      horizon = horizon,
      dry_run = TRUE,
      onset_eligibility = "season"
    ))),
  "K5 season without a duration fit is refused (C and v1)"
)

# ===========================================================================
# Kernel runners: Y as an array [replicate, unit (frame_tu order), month].
run_v1 <- function(d_occ, gate, R, seed = 11L) {
  ms <- suppressWarnings(
    simulate_one_tu(
      tu = "hmmW",
      frame_tu = frame_tu,
      panel_tu = panel_tu,
      draws_occ = d_occ,
      draws_dur = d_dur_end,
      replicate_ids = seq_len(R),
      seed = seed,
      output_part = tempfile("osg_"),
      horizon = horizon,
      dry_run = TRUE,
      onset_eligibility = gate
    )$monthly_summary
  )
  a <- array(NA_integer_, c(R, nrow(frame_tu), n_m))
  a[cbind(
    ms$replicate_id,
    match(ms$gid, frame_tu$gid),
    match(as.Date(ms$month), months)
  )] <- ms$Y
  a
}
run_c <- function(d_occ, gate, R, seed = 11L) {
  o <- build_spatpomp_tu_c(
    frame_tu,
    panel_tu,
    d_occ,
    window,
    fixed_draw_idx = 1L,
    draws_dur = d_dur_end,
    onset_eligibility = gate
  )
  a <- array(NA_integer_, c(R, nrow(frame_tu), n_m))
  for (r in seq_len(R)) {
    set.seed(seed + r)
    sm <- spatPomp::simulate(o, format = "data.frame")
    sm <- sm[sm$time >= 1 & sm$time <= n_m, ]
    a[cbind(r, match(sm$unitname, frame_tu$gid), sm$time)] <- as.integer(
      sm$Yout
    )
  }
  a
}
max_per_season <- function(a) {
  # max onsets in any (replicate, unit, season)
  sid <- t(vapply(
    seq_len(nrow(frame_tu)),
    function(i) season_id(frame_tu$min_seas[i], months),
    integer(n_m)
  ))
  mx <- 0L
  for (r in seq_len(dim(a)[1])) {
    for (i in seq_len(nrow(frame_tu))) {
      mx <- max(mx, tapply(a[r, i, ], sid[i, ], sum))
    }
  }
  mx
}

# ===========================================================================
# D  degenerate kernels (onset certain, cessation certain, 15 months, four
#    units with min_seas 1/4/7/11): inactive re-onsets every month; season
#    onsets exactly in month 1 and in every season-start month, v1 and C
#    alike.
cat("\n## D: degenerate kernels (onset and cessation certain)\n")
want_season <- matrix(0L, nrow(frame_tu), n_m)
want_season[, 1] <- 1L
want_season[sk_fix == 0L] <- 1L
for (k in c("v1", "C")) {
  run <- if (k == "v1") run_v1 else run_c
  a_in <- run(d_occ_hi, "inactive", 2L)
  a_se <- run(d_occ_hi, "season", 2L)
  expect(all(a_in == 1L), sprintf("D1 %s inactive: re-onset every month", k))
  expect(
    identical(a_se[1, , ], want_season) && identical(a_se[2, , ], want_season),
    sprintf(
      "D2 %s season: onsets in month 1 and each season-start month only",
      k
    )
  )
}

# ===========================================================================
# M  moderate onset hazard, certain cessation, 30 replicates: inactive gives
#    repeated onsets within a unit-season; season gives at most one per
#    unit-season, v1 and C alike.
cat("\n## M: moderate onset hazard, 30 replicates\n")
for (k in c("v1", "C")) {
  run <- if (k == "v1") run_v1 else run_c
  m_in <- max_per_season(run(d_occ_mod, "inactive", 30L))
  m_se <- max_per_season(run(d_occ_mod, "season", 30L))
  cat(sprintf(
    "  %s: max onsets per unit-season inactive %d, season %d\n",
    k,
    m_in,
    m_se
  ))
  expect(
    m_in >= 2L,
    sprintf("M1 %s inactive: repeated onsets within a unit-season", k)
  )
  expect(
    m_se <= 1L,
    sprintf("M2 %s season: at most one onset per unit-season", k)
  )
}

# ===========================================================================
# P  cache identity: .make_filter_cache_path() `_onset-season` token only
#    under season; smoothing meta.yml field written only when not inactive;
#    a cache built under the other value is refused.
cat("\n## P: cache identity\n")
ww <- list(start = "2023-01-01", end = "2024-01-01")
p0 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15)
p1 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15, "inactive")
p2 <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15, "season")
expect(
  identical(p0, p1) && identical(p2, sub("\\.qs$", "_onset-season.qs", p0)),
  "P1 filter cache path: `_onset-season` only under season"
)
expect(
  is.null(
    smoothing_cache_identity(
      list(NP = 10L, N_REPS = 2L),
      optional = smoothing_cache_optional_fields(list())
    )$onset_eligibility
  ) &&
    identical(
      smoothing_cache_identity(
        list(NP = 10L, N_REPS = 2L),
        optional = smoothing_cache_optional_fields(
          list(ONSET_ELIGIBILITY = "season")
        )
      )$onset_eligibility,
      "season"
    ),
  "P2 smoothing meta.yml: onset_eligibility written only when not inactive"
)
cd <- tempfile("osg_cache_")
dir.create(cd)
yaml::write_yaml(
  smoothing_cache_meta("transmission_unit"),
  smoothing_cache_meta_path(cd)
)
exp_season <- c(
  smoothing_cache_meta("transmission_unit"),
  list(onset_eligibility = "season")
)
exp_inact <- c(
  smoothing_cache_meta("transmission_unit"),
  list(onset_eligibility = "inactive")
)
expect(
  !is.null(err_of(check_smoothing_cache_meta(cd, exp_season))) &&
    is.null(err_of(check_smoothing_cache_meta(cd, exp_inact))),
  "P3 a legacy (inactive) cache is refused by a season run, accepted by an inactive one"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
