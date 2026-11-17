# test_results_core_sim_aggregates.R — HARD: the simulator-partition helpers
#   of utils/results_core.R shared by 09b and 07d.
#
# Usage:
#   Rscript analysis/tests/test_results_core_sim_aggregates.R
#
# Details:
#   Quick-safe: no fit, no data, no MCMC. Builds a tiny fake
#   <sim_root>/<run_tag>/parts/tu=KEN/repoff=0/{monthly_summary.qs,_SUCCESS}
#   tree. load_sim_run_monthly() (cases = the mon track's cases_month)
#   equals a reference load_run() fixture, and a run without cases_month
#   stops; adm1_prob_stats(), summarise_reps() and the paired stats are
#   checked against hand-computed values.

suppressPackageStartupMessages({
  library(dplyr)
})
source("analysis/utils/results_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
# route expect_equal() through the PASS/FAIL counter this file uses
expect_equal <- function(a, b, label, tol = 1e-9) {
  expect(isTRUE(all.equal(a, b, tolerance = tol)), label)
}

cat(
  "# results_core simulator aggregates (load_sim_run_monthly / ADM1 stats / ribbons)\n"
)

# ---- fake partition tree -------------------------------------------------
root <- tempfile("sim_spatpomp_")
tag <- "fake_run"
tu <- file.path(root, tag, "parts", "tu=KEN", "repoff=0")
dir.create(tu, recursive = TRUE)
invisible(file.create(file.path(tu, "_SUCCESS")))
ms <- expand.grid(
  replicate_id = 1:3,
  gid = c("KEN.1_1", "KEN.2_1"),
  month = c("2025-01-01", "2025-02-01", "2025-03-01"),
  stringsAsFactors = FALSE
)
ms$Y <- as.integer((ms$replicate_id + seq_len(nrow(ms))) %% 2)
ms$cases_month <- ifelse(ms$Y == 1L, 3 * ms$replicate_id, NA_real_)
qs2::qs_save(ms, file.path(tu, "monthly_summary.qs"))

w0 <- as.Date("2025-01-01")
w1 <- as.Date("2025-02-01")

# ---- literal fixture of the original body ----------------------------------
fixture_old <- function(col = "cases_month") {
  m <- ms %>%
    transmute(
      replicate_id = as.integer(replicate_id),
      gid = as.character(gid),
      month = as.Date(month),
      Y = as.integer(Y),
      cases = ifelse(is.na(.data[[col]]), 0, as.numeric(.data[[col]]))
    )
  list(
    country_month = m %>%
      group_by(replicate_id, month) %>%
      summarise(
        n_cases = sum(cases),
        n_outbreaks = sum(Y),
        .groups = "drop"
      ) %>%
      mutate(country = "KEN", run = tag),
    adm1_month = m %>%
      filter(month >= w0, month <= w1) %>%
      mutate(
        adm1_unit = stringr::str_extract(
          stringr::str_replace_all(gid, "\\.", "_"),
          "^[A-Z]{3}_[0-9]+"
        )
      ) %>%
      group_by(replicate_id, adm1_unit, month) %>%
      summarise(cases = sum(cases), outbreaks = sum(Y), .groups = "drop") %>%
      mutate(run = tag)
  )
}
got <- load_sim_run_monthly(root, tag, w0, w1)
want <- fixture_old()
# max_gid_cases (per-gid explosive cap) is the one column the
# reference fixture does not have.
expect(
  isTRUE(all.equal(
    as.data.frame(dplyr::select(got$country_month, -max_gid_cases)),
    as.data.frame(want$country_month)
  )),
  "country_month (cases_month) equals the original body"
)
expect(
  isTRUE(all.equal(
    as.data.frame(got$adm1_month),
    as.data.frame(want$adm1_month)
  )),
  "adm1_month (cases_month) equals the original body"
)
expect(
  !"cases_col" %in% names(formals(load_sim_run_monthly)),
  "load_sim_run_monthly() has no cases_col switch"
)
expect(
  identical(names(got)[1:2], c("country_month", "adm1_month")),
  "load_sim_run_monthly() keeps country_month / adm1_month first (07d indexes them)"
)
expect(
  identical(names(got$gids), c("gid", "country", "run")) &&
    setequal(got$gids$gid, unique(ms$gid)) &&
    !anyDuplicated(got$gids$gid) &&
    all(got$gids$country == "KEN"),
  "gids = the distinct simulated gids with their ISO3 (population denominators)"
)

# ---- region-partitioned layout: tu=<enso group>, gids from several countries
# (partition_by = enso_region). The country must come from the gid
# ISO3, not the directory name (else every gid of the region would be filed
# under the "country" "tu=enso_east").
tag_reg <- "fake_run_region"
tu_reg <- file.path(root, tag_reg, "parts", "tu=enso_east", "repoff=0")
dir.create(tu_reg, recursive = TRUE)
invisible(file.create(file.path(tu_reg, "_SUCCESS")))
ms_reg <- expand.grid(
  replicate_id = 1:2,
  gid = c("KEN.1_1", "UGA.3_1", "UGA.4_1"),
  month = c("2025-01-01", "2025-02-01"),
  stringsAsFactors = FALSE
)
ms_reg$Y <- 1L
ms_reg$cases_month <- ifelse(grepl("^KEN", ms_reg$gid), 10, 1)
qs2::qs_save(ms_reg, file.path(tu_reg, "monthly_summary.qs"))
got_reg <- load_sim_run_monthly(root, tag_reg, w0, w1)$country_month
expect(
  setequal(unique(got_reg$country), c("KEN", "UGA")),
  "region partition: country rows keyed by gid ISO3 (KEN, UGA), not 'tu=enso_east'"
)
expect(
  all(got_reg$n_cases[got_reg$country == "KEN"] == 10) &&
    all(got_reg$n_cases[got_reg$country == "UGA"] == 2) &&
    all(got_reg$n_outbreaks[got_reg$country == "UGA"] == 2),
  "region partition: per-country sums are per country"
)
expect(
  identical(
    names(got_reg),
    c(
      "replicate_id",
      "month",
      "n_cases",
      "n_outbreaks",
      "country",
      "max_gid_cases",
      "run"
    )
  ),
  "country_month keeps its column order"
)
# A gid with no ISO3 prefix is refused, naming it.
ms_bad <- ms_reg
ms_bad$gid[1] <- "composite-1" # no leading ISO3
tag_bad <- "fake_run_badgid"
tu_bad <- file.path(root, tag_bad, "parts", "tu=enso_east", "repoff=0")
dir.create(tu_bad, recursive = TRUE)
invisible(file.create(file.path(tu_bad, "_SUCCESS")))
qs2::qs_save(ms_bad, file.path(tu_bad, "monthly_summary.qs"))
err_bad <- try(load_sim_run_monthly(root, tag_bad, w0, w1), silent = TRUE)
expect(
  inherits(err_bad, "try-error") && grepl("composite-1", err_bad),
  "a gid without an ISO3 prefix stops, naming the gid"
)
# A run without cases_month (not configured with grouping.monthly_enable) stops.
tag_nomon <- "fake_run_nomon"
tu_nomon <- file.path(root, tag_nomon, "parts", "tu=KEN", "repoff=0")
dir.create(tu_nomon, recursive = TRUE)
invisible(file.create(file.path(tu_nomon, "_SUCCESS")))
qs2::qs_save(
  ms[, setdiff(names(ms), "cases_month")],
  file.path(tu_nomon, "monthly_summary.qs")
)
err_nomon <- try(load_sim_run_monthly(root, tag_nomon, w0, w1), silent = TRUE)
expect(
  inherits(err_nomon, "try-error") && grepl("monthly_enable", err_nomon),
  "a run without cases_month stops, naming grouping.monthly_enable"
)
expect(
  inherits(
    try(
      load_sim_run_monthly(root, tag, w0, w1, series = "detected"),
      silent = TRUE
    ),
    "try-error"
  ),
  "series = 'detected' on a run without a detected column stops"
)
expect(
  inherits(
    try(load_sim_run_monthly(root, "nope", w0, w1), silent = TRUE),
    "try-error"
  ),
  "missing run dir stops"
)

# adm1_prob_stats() against hand-computed outer() fractions, keep_means =
# FALSE slim form.
# ---- adm1_prob_stats ------------------------------------------------------
f <- c(1, 5, 9)
cc <- c(2, 4, 6)
d <- bind_rows(
  tibble(
    run = "forecast",
    adm1_unit = "KEN_1",
    replicate_id = 1:3,
    month = 1,
    cases = f
  ),
  tibble(
    run = "counterfactual",
    adm1_unit = "KEN_1",
    replicate_id = 1:3,
    month = 1,
    cases = cc
  )
)
thr <- 0.5
ps <- adm1_prob_stats(d, "cases", thr)
expect_equal(
  ps$prob_larger,
  mean(outer(f, (1 + thr) * cc, `>`)),
  "prob_larger == outer() fraction"
)
expect_equal(
  ps$prob_smaller,
  mean(outer(f, (1 - thr) * cc, `<`)),
  "prob_smaller == outer() fraction"
)
expect_equal(
  c(ps$mean_forecast, ps$mean_counterfactual),
  c(mean(f), mean(cc)),
  "arm means"
)
slim <- adm1_prob_stats(d, "cases", thr, keep_means = FALSE)
expect(
  identical(names(slim), c("adm1_unit", "prob_larger", "prob_smaller")),
  "keep_means = FALSE returns only the two probabilities"
)
fin <- finish_adm1_prob(
  ps,
  low_signal_thresh = 2,
  sig_threshold = 0.5,
  country_region = tibble(country = "KEN", region = "East"),
  adm1_group_horizon = tibble(adm1_unit = "KEN_1", horizon = 1)
)
expect(
  fin$country == "KEN" && !fin$low_signal && fin$region == "East",
  "finish_adm1_prob attaches country / low_signal / crosswalks"
)
expect(
  fin$signed_prob ==
    ifelse(ps$prob_larger >= ps$prob_smaller, ps$prob_larger, -ps$prob_smaller),
  "signed_prob picks the dominant direction"
)

# paired_adm1_stats() / finish_adm1_ratio() / finish_adm1_prob() on a
# hand-built case.
# ---- paired stats / ratio ---------------------------------------------------
pr <- paired_adm1_stats(d, "cases", ci_level = 0.5)
expect_equal(pr$mean_diff, mean(f - cc), "paired mean_diff")
expect_equal(
  pr$prob_forecast_higher,
  mean(f - cc > 0),
  "paired prob_forecast_higher"
)
rt <- finish_adm1_ratio(
  pr,
  low_signal_thresh = 2,
  country_region = tibble(country = "KEN", region = "East"),
  adm1_group_horizon = tibble(adm1_unit = "KEN_1", horizon = 1),
  ratio_eps = 0.25
)
expect_equal(
  rt$ratio,
  (mean(f) + 0.25) / (mean(cc) + 0.25),
  "ratio uses the pseudocount"
)
rt_pop <- finish_adm1_ratio(
  pr,
  2,
  tibble(country = "KEN", region = "East"),
  tibble(adm1_unit = "KEN_1", horizon = 1),
  adm1_pop = tibble(adm1_unit = "KEN_1", pop_adm1 = 1e5)
)
expect_equal(
  rt_pop$incidence_forecast,
  mean(f),
  "adm1_pop adds per-1e5 incidence"
)

# summarise_reps() against stats::quantile().
# ---- summarise_reps ---------------------------------------------------------
set.seed(1)
big <- tibble(
  run = rep(c("a", "b"), each = 40),
  n_cases = rpois(80, 5),
  n_outbreaks = rpois(80, 1)
)
sr <- summarise_reps(big, "run")
expect_equal(
  sr$n_cases_q50[sr$run == "a"],
  unname(quantile(big$n_cases[big$run == "a"], 0.5)),
  "q50"
)
expect_equal(
  sr$n_outbreaks_q90[sr$run == "b"],
  unname(quantile(big$n_outbreaks[big$run == "b"], 0.9)),
  "q90"
)
expect(
  identical(names(sr)[2:3], c("n_cases_q10", "n_cases_q25")),
  "column naming n_cases_q10 ..."
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
