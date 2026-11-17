# test_who_scale_estimate.R — HARD: the per-country WHO term scale.
#
# Usage:
#   Rscript analysis/tests/test_who_scale_estimate.R
#
# Details:
#   Synthetic, no MCMC, no data; run from the project root. Sections E, S, O,
#   D, V, M: observed DB months, the scale estimate, the DB-overlap drop, the
#   scale digest, the WHO-vs-smoothing table and the smoothing cache identity.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/smoothing_core.R") # who_dashboard.R, spatpomp_core.R, base.R
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
cat("# test_who_scale_estimate\n\n")

mo <- as.Date("2023-01-01") %m+% months(0:11)
frame <- tibble(
  gid = c("AAA_1_1", "AAA_2_1", "BBB_1_1", "CCC_1_1"),
  country = c("AAA", "AAA", "BBB", "CCC")
)

# ===========================================================================
# E: who_db_observed_months(): the filter's own per-gid cases summed over frame
#    gids; "observed" = any non-NA gid; window; non-frame gids ignored.
cat("## E: who_db_observed_months()\n")
obs <- tidyr::expand_grid(gid = c(frame$gid, "ZZZ_1_1"), t_start = mo) %>%
  mutate(cases_month = 10)
obs$cases_month[obs$gid == "AAA_2_1" & obs$t_start == mo[1]] <- NA # one AAA gid NA
obs$cases_month[obs$gid == "CCC_1_1"] <- NA # CCC never observed
db <- who_db_observed_months(obs, frame, mo[1], mo[7])
expect(
  identical(sort(unique(db$country)), c("AAA", "BBB", "CCC")) &&
    all(db$t_start < mo[7]),
  "E1 frame countries only (ZZZ ignored), window [start, end)"
)
expect(
  db$db_cases[db$country == "AAA" & db$t_start == mo[1]] == 10 &&
    db$db_cases[db$country == "AAA" & db$t_start == mo[2]] == 20,
  "E2 cases summed over the country's frame gids (NA gid skipped)"
)
expect(
  all(db$db_observed[db$country %in% c("AAA", "BBB")]) &&
    !any(db$db_observed[db$country == "CCC"]),
  "E3 observed = at least one non-NA gid"
)
expect(
  !is.null(err_of(who_db_observed_months(
    dplyr::select(obs, -cases_month),
    frame,
    mo[1],
    mo[7]
  ))),
  "E4 no cases_month column (monthly_case_obs none) refused"
)

# ===========================================================================
# S: estimate_who_scale(): kappa = 0 gives W / E and the multiplier W / D
#    (who_f cancels); shrinkage formula; fallbacks (too few cases / months, no
#    overlap) pooled and fixed; clamp; partial-WHO / DB-unobserved months
#    excluded; no identified country with fallback pooled -> error; counts.
cat("\n## S: estimate_who_scale()\n")
cov <- tibble(
  country = c("AAA", "BBB", "CCC"),
  pop_covered = c(5, 8, 1),
  pop_national = c(10, 10, 10),
  who_f = c(0.5, 0.8, 0.1)
)
dbm <- tidyr::expand_grid(country = c("AAA", "BBB", "CCC"), t_start = mo) %>%
  mutate(
    db_cases = c(AAA = 50, BBB = 10, CCC = 0)[country],
    db_observed = country != "CCC"
  )
who <- tidyr::expand_grid(country = c("AAA", "BBB", "CCC"), t_start = mo) %>%
  mutate(cases_who = c(AAA = 200, BBB = 15, CCC = 30)[country])
st <- function(...) {
  estimate_who_scale(dbm, who, cov, mo[1], mo[12] %m+% months(1), ...)
}
s0 <- st(min_db_cases = 100, min_months = 6, shrink_cases = 0)
a <- s0[s0$country == "AAA", ]
expect(
  isTRUE(all.equal(a$scale, (12 * 200) / (12 * 50 / 0.5))) &&
    isTRUE(all.equal(a$multiplier, 200 / 50)) &&
    isTRUE(all.equal(a$ratio_raw, 4)),
  "S1 kappa = 0: scale = W / E, multiplier = W / D (who_f cancels)"
)
b <- s0[s0$country == "BBB", ]
expect(
  b$source == "estimated" && isTRUE(all.equal(b$multiplier, 15 / 10)),
  "S2 BBB identified (120 DB cases >= 100, 12 months)"
)
cc <- s0[s0$country == "CCC", ]
pooled <- (12 * 200 + 12 * 15) / (12 * 50 / 0.5 + 12 * 10 / 0.8)
expect(
  cc$source == "fallback_pooled" &&
    isTRUE(all.equal(cc$scale, pooled)) &&
    cc$est_months == 0L,
  "S3 CCC (DB never observes it) falls back to the pooled ratio"
)
s_k <- st(min_db_cases = 100, min_months = 6, shrink_cases = 500)
expect(
  isTRUE(all.equal(
    s_k$scale[s_k$country == "BBB"],
    (12 * 15 + 500 * pooled) / (12 * 10 / 0.8 + 500)
  )),
  "S4 shrinkage (W + kappa s0) / (E + kappa)"
)
s_m <- st(
  min_db_cases = 100,
  min_months = 13,
  shrink_cases = 0,
  fallback = "fixed",
  fixed_scale = 1.3
)
expect(
  all(s_m$source == "fallback_fixed") && all(s_m$scale == 1.3),
  "S5 too few overlap months -> fixed fallback = `scale`"
)
expect(
  !is.null(err_of(st(min_db_cases = 1e6, min_months = 1))),
  "S6 no identified country + pooled -> error"
)
s_c <- st(min_db_cases = 100, min_months = 6, clamp = c(0.5, 1.5))
expect(
  isTRUE(all.equal(s_c$scale[s_c$country == "AAA"], 1.5)) &&
    s_c$clamped[s_c$country == "AAA"],
  "S7 clamp, flagged"
)
dbm2 <- dbm
dbm2$db_observed[dbm2$country == "AAA" & dbm2$t_start %in% mo[1:2]] <- FALSE
who2 <- who[!(who$country == "AAA" & who$t_start == mo[3]), ] # an incomplete WHO month
s_x <- estimate_who_scale(
  dbm2,
  who2,
  cov,
  mo[1],
  mo[12] %m+% months(1),
  min_db_cases = 100,
  min_months = 6
)
expect(
  s_x$est_months[s_x$country == "AAA"] == 9L,
  "S8 DB-unobserved and WHO-absent months excluded from O_c"
)
dbm3 <- dbm
dbm3$db_cases[dbm3$country == "BBB" & dbm3$t_start == mo[1]] <- 0
who3 <- who
who3$cases_who[who3$country == "BBB" & who3$t_start == mo[2]] <- 0
s_d <- estimate_who_scale(
  dbm3,
  who3,
  cov,
  mo[1],
  mo[12] %m+% months(1),
  min_db_cases = 100,
  min_months = 6
)
expect(
  s_d$months_who_pos_db_zero[s_d$country == "BBB"] == 1L &&
    s_d$months_db_pos_who_zero[s_d$country == "BBB"] == 1L,
  "S9 D9 diagnostic counts"
)

# ===========================================================================
# O: who_drop_db_overlap() drops exactly the DB-observed country-months.
cat("\n## O: who_drop_db_overlap()\n")
kept <- who_drop_db_overlap(who, dbm2)
expect(
  nrow(kept) == 12L + 2L &&
    all(kept$country %in% c("CCC", "AAA")) &&
    all(kept$t_start[kept$country == "AAA"] %in% mo[1:2]),
  "O1 drops exactly the DB-observed country-months"
)

# ===========================================================================
# D: digest_who_scale(): stable, row-order independent, content-sensitive.
cat("\n## D: digest_who_scale()\n")
d1 <- digest_who_scale(s0)
d2 <- digest_who_scale(s0[c(3, 1, 2), ])
expect(
  is.character(d1) && nchar(d1) == 16L && identical(d1, d2),
  "D1 stable, row-order independent"
)
s0b <- s0
s0b$scale[1] <- s0b$scale[1] * 1.001
expect(!identical(digest_who_scale(s0b), d1), "D2 content-sensitive")

# ===========================================================================
# V: who_vs_smoothing(): who_scale NULL = the fixed-scale table (frozen copy);
#    a scale multiplies the bands by who_multiplier().
cat("\n## V: who_vs_smoothing()\n")
# Frozen copy of the fixed-scale band scaling (smoothing_core.R at c032a12).
old_scale <- function(rib, countries, who_f) {
  f <- unname(who_f[countries])
  f[is.na(f)] <- 1
  fc <- stats::setNames(f, countries)
  qcols <- grep("^n_cases_", names(rib), value = TRUE)
  rib[qcols] <- lapply(rib[qcols], function(v) v / fc[rib$country])
  rib$who_f <- unname(fc[rib$country])
  rib
}
long_mon <- tidyr::expand_grid(
  rep_id = 1:3,
  t_start = mo[1:4],
  country = c("AAA", "BBB")
) %>%
  mutate(cases_month = rep_id * 7 + as.integer(country == "AAA") * 5)
who_m <- tibble(
  iso3 = rep(c("AAA", "BBB"), each = 4),
  t_start = rep(mo[1:4], 2),
  cases = c(30, 40, 50, 60, 5, 6, 7, 8),
  partial = FALSE
)
rw <- tibble(rep_id = 1:3, w = c(0.2, 0.5, 0.3))
wf <- c(AAA = 0.5, BBB = 0.8)
v0 <- who_vs_smoothing(long_mon, who_m, rw, mo[1], mo[5], who_f = wf)
lm_ <- long_mon[, c("rep_id", "t_start", "country", "cases_month")]
rib_ref <- old_scale(
  .agg_and_ribbon_impl(
    lm_,
    "country",
    rw,
    mo[1:4],
    src_col = "cases_month",
    out_col = "n_cases"
  ),
  c("AAA", "BBB"),
  wf
)
expect(
  identical(v0$ribbon, rib_ref),
  "V1 who_scale NULL: ribbon identical to the pre-item-10 scaling"
)
v1 <- who_vs_smoothing(
  long_mon,
  who_m,
  rw,
  mo[1],
  mo[5],
  who_f = wf,
  who_scale = c(AAA = 2, BBB = 1)
)
expect(
  isTRUE(all.equal(
    v1$ribbon$n_cases_q50[v1$ribbon$country == "AAA"],
    v0$ribbon$n_cases_q50[v0$ribbon$country == "AAA"] * 2
  )) &&
    isTRUE(all.equal(
      v1$ribbon$n_cases_q50[v1$ribbon$country == "BBB"],
      v0$ribbon$n_cases_q50[v0$ribbon$country == "BBB"]
    )) &&
    "who_scale" %in% names(v1$summary),
  "V2 a per-country scale multiplies the bands (scale / who_f)"
)

# ===========================================================================
# M: smoothing cache identity: who_scale_digest / db_overlap written only when
#    not the legacy value; record-only on the forecast read.
cat("\n## M: smoothing cache identity\n")
S_fixed <- list(
  N_REPS = 2L,
  NP = 10L,
  PARTITION_BY = "country",
  run_tag = "x",
  WHO_OBS = list(enable = TRUE, db_overlap = "score_both")
)
expect(
  is.null(
    smoothing_cache_identity(
      S_fixed,
      "a",
      smoothing_cache_optional_fields(S_fixed)
    )$who_scale_digest
  ) &&
    is.null(
      smoothing_cache_identity(
        S_fixed,
        "a",
        smoothing_cache_optional_fields(S_fixed)
      )$db_overlap
    ),
  "M1 fixed scale + score_both: nothing written"
)
S_est <- modifyList(
  S_fixed,
  list(
    WHO_SCALE_DIGEST = d1,
    WHO_OBS = list(enable = TRUE, db_overlap = "skip_who")
  )
)
m_est <- smoothing_cache_identity(
  S_est,
  "a",
  smoothing_cache_optional_fields(S_est)
)
expect(
  identical(m_est$who_scale_digest, d1) &&
    identical(m_est$db_overlap, "skip_who"),
  "M2 estimated + skip_who: both recorded"
)
dm <- tempfile("who_cache_")
dir.create(dm)
yaml::write_yaml(m_est, smoothing_cache_meta_path(dm))
expect(
  is.null(err_of(check_smoothing_cache_meta(
    dm,
    smoothing_cache_meta("country")
  ))),
  "M3 a forecast (no expectation) reads an estimated-scale cache: record-only fields not checked"
)
expect(
  !is.null(err_of(check_smoothing_cache_meta(
    dm,
    c(smoothing_cache_meta("country"), list(who_scale_digest = "other"))
  ))),
  "M4 an explicit expectation is checked"
)
unlink(dm, recursive = TRUE)

cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
