# test_results_attrib_metrics.R — HARD: report_attrib_metrics
#   (utils/results_derive.R) on synthetic paired runs.
#
# Usage:
#   Rscript analysis/tests/test_results_attrib_metrics.R
#
# Details:
#   Quick-safe: no fit, no data, no simulation output. Attributable number /
#   fraction, probability ratio with the probabilities of causation, and
#   signal-to-noise of paired forecast / ENSO-counterfactual runs.

suppressPackageStartupMessages({
  library(dplyr)
})
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# report_attrib_metrics\n")
near <- function(a, b, tol = 1e-9) isTRUE(all.equal(a, b, tolerance = tol))
I1 <- function(n) matrix(1L, n, 1L) # identity resample: every replicate once

# [1] kernels by hand: af_stats(), pr_stats() (Inf, exact bound when p0 = 0,
#     both zero, p0 = 1), sn_stats(), .boot_q() with Inf;
cat("# [1] kernels\n")
a <- af_stats(c(30, 30, 30, 0), rep(10, 4), I1(4), pop = 1e5)
expect(
  a$n_pairs == 4L &&
    a$n_f_zero == 1L &&
    near(a$mean_f, 22.5) &&
    near(a$an_mean, 12.5),
  "af_stats: means, AN of the means, zero-forecast count"
)
expect(
  near(a$af_of_means, 12.5 / 22.5) &&
    near(a$af_med, 2 / 3) &&
    near(a$an_per_100k, 12.5),
  "af_stats: AF of the means, per-replicate AF median (NA where f = 0 dropped), AN per 100k"
)
expect(
  near(a$an_lo, a$an_mean) &&
    near(a$an_hi, a$an_mean) &&
    near(a$af_lo, a$af_of_means),
  "af_stats: identity resample -> interval collapses on the point"
)
w2 <- matrix(c(2L, 1L, 1L, 0L), 4L, 1L) # replicate 4 dropped, 1 doubled
a2 <- af_stats(c(30, 30, 30, 0), rep(10, 4), w2)
expect(
  near(a2$an_lo, 20),
  "af_stats: resampled AN = mean over the resample (f = 30 x 4 draws)"
)
expect(
  is.na(af_stats(numeric(0), numeric(0), matrix(0L, 0L, 3L))$an_mean),
  "af_stats: no pairs -> NA"
)

p <- pr_stats(c(1, 1, 1), c(0, 0, 0), I1(3), alpha = 0.05)
expect(
  is.infinite(p$pr) &&
    p$p0_zero &&
    near(p$far, 1) &&
    near(p$pns, 1) &&
    near(p$ps, 1),
  "pr_stats: p0 = 0 -> PR Inf, FAR 1, PS 1, PNS 1"
)
expect(
  near(p$p0_upper, 1 - 0.05^(1 / 3)) && near(p$pr_lo, 1 / (1 - 0.05^(1 / 3))),
  "pr_stats: exact one-sided bound on p0 -> pr_lo = p1 / p0_upper"
)
p2 <- pr_stats(c(1, 1, 0, 0), c(1, 0, 0, 0), I1(4), alpha = 0.05)
expect(
  near(p2$pr, 2) &&
    near(p2$far, 0.5) &&
    near(p2$pns, 0.25) &&
    near(p2$ps, 0.25 / 0.75) &&
    !p2$p0_zero,
  "pr_stats: p1 0.5, p0 0.25 -> PR 2, FAR (PN) 0.5, PNS 0.25, PS 1/3"
)
p3 <- pr_stats(c(0, 0), c(0, 0), I1(2), alpha = 0.05)
expect(
  is.na(p3$pr) && is.na(p3$far) && !p3$p0_zero,
  "pr_stats: neither arm has the event -> no ratio"
)
p4 <- pr_stats(c(1, 1), c(1, 1), I1(2), alpha = 0.05)
expect(
  near(p4$pr, 1) && is.na(p4$ps) && near(p4$pns, 0) && near(p4$far, 0),
  "pr_stats: p0 = 1 -> PS undefined"
)
expect(
  identical(.boot_q(c(1, Inf, 2, Inf), c(0.025, 0.975)), c(1, Inf)) &&
    !anyNA(.boot_q(c(Inf, Inf), 0.5)),
  ".boot_q: Inf-safe (no NaN)"
)

s <- sn_stats(c(30, 30, 30, 0), rep(10, 4))
expect(
  near(s$mean_diff, 12.5) &&
    near(s$sd_diff, 15) &&
    near(s$sn, 12.5 / 15) &&
    near(s$p_sign, 0.75),
  "sn_stats: d = 20,20,20,-10 -> mean 12.5, sd 15, S/N 0.833, p_sign 0.75"
)
expect(
  is.na(sn_stats(rep(30, 3), rep(10, 3))$sn),
  "sn_stats: constant difference -> sd 0 -> S/N NA"
)
expect(
  sn_stats(rep(0, 3) + c(1, 2, 3), rep(10, 3))$sn < 0,
  "sn_stats: sign kept"
)

# [2] attrib_metrics_tables() end to end on the report_attrib_exposure
#     fixture (2 groups x 3 ADM1 x 4 replicates, one explosive replicate):
#     guard, scope sums, Africa common ids, population tables, determinism;
cat("# [2] attrib_metrics_tables() end to end\n")
a0 <- as.Date("2026-05-01")
mk_am <- function(ken2_f = 12) {
  am <- expand.grid(
    replicate_id = 1:4,
    adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
    run = c("forecast", "counterfactual"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      month = a0,
      analog_set = ifelse(run == "forecast", NA_character_, "all"),
      enso_group = ifelse(
        substr(adm1_unit, 1, 3) == "KEN",
        "enso_east",
        "enso_horn"
      ),
      group_horizon_end = as.Date("2026-06-01"),
      cases = dplyr::case_when(
        run == "counterfactual" ~ 10,
        adm1_unit == "KEN_1" ~ 30,
        adm1_unit == "KEN_2" ~ ken2_f[pmin(replicate_id, length(ken2_f))],
        TRUE ~ c(30, 30, 30, 0)[replicate_id]
      ),
      outbreaks = as.integer(cases > 10)
    ) %>%
    tibble::as_tibble()
  am$cases[
    am$run == "forecast" & am$replicate_id == 4 & am$adm1_unit == "KEN_1"
  ] <- 2e5
  am
}
mk_cm <- function(am) {
  am %>%
    dplyr::mutate(country = substr(adm1_unit, 1, 3)) %>%
    dplyr::group_by(
      replicate_id,
      month,
      country,
      run,
      analog_set,
      enso_group,
      group_horizon_end
    ) %>%
    dplyr::summarise(
      n_cases = sum(cases),
      n_outbreaks = sum(outbreaks),
      .groups = "drop"
    )
}
am3 <- mk_am()
cm3 <- mk_cm(am3)
pop3 <- tibble::tibble(
  adm1_unit = c("KEN_1", "KEN_2", "UGA_1"),
  pop = c(100, 40, NA)
)
gs <- c(enso_east = "East", enso_horn = "Horn")
out <- attrib_metrics_tables(
  am3,
  cm3,
  pop3,
  gs,
  "all",
  a0,
  cap = 2e5,
  pr_case_threshold = 20,
  pr_alpha = 0.05,
  sn_levels = c(1, 2),
  n_boot = 50L
)
expect(
  identical(out$explosive$replicate_id, 4L),
  "explosive replicate keyed to its group"
)
af <- dplyr::filter(out$af_adm1, metric == "cases")
expect(
  identical(af$n_pairs, c(3L, 3L, 4L)),
  "guard drops replicate 4 from the KEN ADM1s only"
)
expect(
  near(af$an_mean[af$adm1_unit == "KEN_1"], 20) &&
    near(af$af_of_means[af$adm1_unit == "KEN_1"], 2 / 3),
  "KEN_1: AN 20, AF 2/3"
)
afs <- dplyr::filter(out$af_scope, metric == "cases")
east <- dplyr::filter(afs, level == "region", enso_group == "enso_east")
expect(
  near(east$an_mean, sum(af$an_mean[af$enso_group == "enso_east"])),
  "region AN = sum of its ADM1 ANs (same paired ids)"
)
afr <- dplyr::filter(afs, level == "continent")
expect(
  afr$n_pairs == 3L &&
    near(afr$mean_f, 30 + 12 + 30) &&
    near(afr$an_mean, 72 - 30),
  "Africa: ids common to every group (1-3), totals summed over ADM1s"
)
expect(
  afr$an_lo <= afr$an_mean && afr$an_mean <= afr$an_hi,
  "Africa AN interval brackets the point"
)

pr <- dplyr::filter(out$pr_adm1, event == "onset_any")
expect(
  all(is.infinite(pr$pr)) && all(pr$p0_zero),
  "onset_any: no counterfactual onset (10 > 10 false) -> PR Inf"
)
expect(
  near(pr$pr_lo[pr$adm1_unit == "UGA_1"], 0.75 / (1 - 0.05^(1 / 4))),
  "UGA_1 PR lower bound, n = 4"
)
cg <- dplyr::filter(out$pr_adm1, event == "cases_ge")
expect(
  is.na(cg$pr[cg$adm1_unit == "KEN_2"]),
  "cases_ge (>= 20): KEN_2 12 vs 10 -> neither arm -> no ratio"
)
cgk <- dplyr::filter(
  out$pr_scope,
  event == "cases_ge",
  level == "country",
  scope == "KEN"
)
expect(
  near(cgk$pr, 1) && is.na(cgk$ps),
  "cases_ge at country scope saturates: 42 vs 20 -> PR 1, PS undefined"
)
pp <- dplyr::filter(out$pr_pop, event == "onset_any", rule == "lower")
expect(
  near(pp$pop[pp$level == "continent"], sum(pp$pop[pp$level == "region"])) &&
    pp$n_adm1[pp$level == "continent"] == 3L,
  "pr_pop: Africa = sum of regions; every ADM1 has pr_lo > 1"
)
expect(
  near(pp$pct[pp$level == "continent"], 1),
  "pr_pop: pct = pop / pop_universe"
)

am5 <- mk_am(ken2_f = c(12, 13, 18, 18))
cm5 <- mk_cm(am5)
o5 <- attrib_metrics_tables(
  am5,
  cm5,
  pop3,
  gs,
  "all",
  a0,
  cap = 2e5,
  pr_case_threshold = 20,
  pr_alpha = 0.05,
  sn_levels = c(1, 2),
  n_boot = 20L
)
sn5 <- dplyr::filter(o5$sn_adm1, metric == "cases", adm1_unit == "KEN_2")
expect(
  near(sn5$sn, mean(c(2, 3, 8)) / sd(c(2, 3, 8))),
  "KEN_2 S/N over the 3 guarded pairs (1.35)"
)
sp <- dplyr::filter(o5$sn_pop, metric == "cases", level == "continent")
at <- function(L, d) sp$n_adm1[sp$sn_level == L & sp$direction == d]
expect(
  at(1, "increase") == 1L &&
    at(2, "increase") == 0L &&
    at(1, "decrease") == 0L &&
    at(1, "either") == at(1, "increase") + at(1, "decrease"),
  "sn_pop: KEN_2 emerged at |S/N| >= 1 (1.35), not 2; either = increase + decrease"
)
expect(
  near(sp$pop[sp$sn_level == 1 & sp$direction == "increase"], 40),
  "sn_pop: population of KEN_2"
)

out2 <- attrib_metrics_tables(
  am3,
  cm3,
  pop3,
  gs,
  "all",
  a0,
  cap = 2e5,
  pr_case_threshold = 20,
  pr_alpha = 0.05,
  sn_levels = c(1, 2),
  n_boot = 50L
)
expect(identical(out, out2), "deterministic (fixed bootstrap seed)")

# [3] attrib_metrics_windows(): window x guard cells.
cat("# [3] attrib_metrics_windows()\n")
am4 <- dplyr::bind_rows(
  am3,
  dplyr::mutate(am3, month = as.Date("2026-06-01"), cases = 10, outbreaks = 0L)
)
cm4 <- mk_cm(am4)
ow <- attrib_metrics_windows(
  am4,
  cm4,
  pop3,
  gs,
  "all",
  attrib_cell_table(a0, as.Date("2026-06-01"), cap = 2e5),
  pr_case_threshold = 20,
  pr_alpha = 0.05,
  sn_levels = c(1, 2),
  n_boot = 20L
)
expect(
  all(vapply(
    ow[c("af_adm1", "pr_scope", "sn_pop")],
    function(d) {
      identical(names(d)[1:2], c("window", "explosive_guard"))
    },
    logical(1)
  )),
  "every table starts with window, explosive_guard"
)
expect(
  setequal(
    unique(paste(ow$af_scope$window, ow$af_scope$explosive_guard)),
    c(
      "attribution TRUE",
      "attribution FALSE",
      "post_fork TRUE",
      "post_fork FALSE"
    )
  ),
  "four window x guard cells"
)
pf <- dplyr::filter(
  ow$af_scope,
  window == "post_fork",
  explosive_guard,
  level == "continent",
  metric == "cases"
)
pfr <- dplyr::filter(
  ow$pr_scope,
  window == "post_fork",
  explosive_guard,
  level == "continent",
  event == "onset_any"
)
expect(
  near(pf$an_mean, 0) && is.na(pfr$pr),
  "post_fork with equal arms: AN 0, no onset in either arm -> no PR"
)
ng <- dplyr::filter(
  ow$af_adm1,
  window == "attribution",
  !explosive_guard,
  metric == "cases"
)
expect(all(ng$n_pairs == 4L), "explosive_guard FALSE keeps every replicate")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
