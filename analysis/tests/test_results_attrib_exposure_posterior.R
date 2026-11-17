# test_results_attrib_exposure_posterior.R — HARD: the posterior interval of
#   the Childs population.
#
# Usage:
#   Rscript analysis/tests/test_results_attrib_exposure_posterior.R
#
# Details:
#   Quick-safe: no fit, no data, no simulation output. Pins
#   results_derive.R::attrib_exposure_posterior_tables() and the
#   attrib_exposure_posterior step: recovery of known per-draw AFs,
#   unidentified ADM1s, windows and explosive guards, refusals.

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
cat("# posterior Childs population\n")
near <- function(a, b, tol) isTRUE(abs(a - b) <= tol)

#' A fixed-draw design: D draws x R replicates; group g on units u1 (AAA_1,
#' AAA_2) and u2 (BBB_1), group h on u3 (CCC_1). Per ADM1 and draw d, an
#' expected forecast total mu and a true AF af[d]; counterfactual mean
#' mu * (1 - af[d]); negative-binomial noise.
synth <- function(D = 99L, R = 50L, nb_size = 2, seed = 1L) {
  set.seed(seed)
  a0 <- as.Date("2026-05-01")
  spec <- list(
    AAA_1 = list(g = "g", u = "u1", mu = 3000, af = rnorm(D, 0.40, 0.15)),
    AAA_2 = list(g = "g", u = "u1", mu = 2000, af = rep(0.60, D)),
    BBB_1 = list(g = "g", u = "u2", mu = 50, af = rep(0, D)),
    CCC_1 = list(g = "h", u = "u3", mu = 1000, af = rep(0.70, D))
  )
  n <- D * R
  rid <- seq_len(n)
  draw <- ((rid - 1L) %% D) + 1L
  rows <- lapply(names(spec), function(a) {
    s <- spec[[a]]
    f <- rnbinom(n, mu = s$mu, size = nb_size)
    cf <- rnbinom(n, mu = s$mu * (1 - s$af[draw]), size = nb_size)
    dplyr::bind_rows(
      tibble::tibble(run = "forecast", analog_set = NA_character_, cases = f),
      tibble::tibble(run = "counterfactual", analog_set = "neutral", cases = cf)
    ) %>%
      dplyr::mutate(
        adm1_unit = a,
        enso_group = s$g,
        unit = s$u,
        replicate_id = rep(rid, 2L),
        month = a0,
        group_horizon_end = as.Date("2026-12-01"),
        outbreaks = as.numeric(cases > s$mu)
      )
  })
  am <- dplyr::bind_rows(rows)
  cm <- am %>%
    dplyr::mutate(country = substr(adm1_unit, 1, 3)) %>%
    dplyr::group_by(
      replicate_id,
      month,
      country,
      run,
      analog_set,
      enso_group,
      unit,
      group_horizon_end
    ) %>%
    dplyr::summarise(
      n_cases = sum(cases),
      n_outbreaks = sum(outbreaks),
      .groups = "drop"
    )
  draws <- tidyr::expand_grid(
    unit = c("u1", "u2", "u3"),
    replicate_id = rid
  ) %>%
    dplyr::mutate(draw = ((replicate_id - 1L) %% D) + 1L)
  list(am = am, cm = cm, draws = draws, spec = spec, a0 = a0)
}
sy <- synth()
pop <- tibble::tibble(
  adm1_unit = c("AAA_1", "AAA_2", "BBB_1", "CCC_1"),
  pop = c(100, 40, 60, 30)
)
gs <- c(g = "G", h = "H")
M <- 0.5 # Childs threshold: ratio >= 1.5, AF >= 1/3

# --- [1] recovery -------------------------------------------------------------
# [1] recovery: a synthetic fixed-draw design (99 draws x 50 replicates,
#     negative-binomial epidemic noise) with KNOWN per-draw AFs per ADM1: an
#     ADM1 whose AF straddles the Childs threshold passes on the share of
#     draws its true AF does; constant strong / null effects pass on every /
#     no draw (a noisy null keeps its pooled class); the population interval
#     follows (Africa, region); pop_point = report_attrib_exposure's Childs
#     population;
cat("# [1] recovery of known per-draw AFs\n")
out <- suppressMessages(attrib_exposure_posterior_tables(
  sy$am,
  sy$cm,
  pop,
  gs,
  "neutral",
  sy$a0,
  M,
  sy$draws,
  cap = 2e5,
  n_mc = 2000L
))
ad <- dplyr::filter(out$adm1, metric == "cases")
p_of <- function(a) ad$p_adj[ad$adm1_unit == a]
p_true <- mean(sy$spec$AAA_1$af >= M / (1 + M))
expect(
  near(p_of("AAA_1"), p_true, 0.10),
  sprintf(
    "straddling ADM1: share of draws passing %.2f ~ true %.2f (+/- 0.10)",
    p_of("AAA_1"),
    p_true
  )
)
expect(
  p_of("AAA_2") == 1 && p_of("CCC_1") == 1,
  "constant strong effects pass on every draw"
)
expect(
  p_of("BBB_1") == 0,
  sprintf(
    "noisy null (noise share %.2f) keeps its pooled class: no draw passes",
    ad$noise_share[ad$adm1_unit == "BBB_1"]
  )
)
expect(
  all(ad$n_draws == 99L) && out$n_min == 50L,
  "99 draws per ADM1, 50 replicates per (unit, draw)"
)
ex <- dplyr::filter(out$exposure, metric == "cases", variant == "adj")
af <- dplyr::filter(ex, level == "continent")
expect(
  af$pop_lo == 70 && af$pop_hi == 170 && af$pop_med == 170,
  sprintf(
    "Africa: 95%% [%g, %g] median %g = AAA_2 + CCC_1 (+ AAA_1 on its passing draws)",
    af$pop_lo,
    af$pop_hi,
    af$pop_med
  )
)
expect(
  near(af$pop_mean, 70 + 100 * p_of("AAA_1"), 3),
  "Africa mean = fixed part + AAA_1 x its pass share"
)
rg <- dplyr::filter(ex, level == "region", scope == "G")
expect(rg$pop_lo == 40 && rg$pop_hi == 140, "region G: 95% [40, 140]")
rh <- dplyr::filter(ex, level == "region", scope == "H")
expect(rh$pop_lo == 30 && rh$pop_hi == 30, "region H: no posterior spread")
x0 <- suppressMessages(attrib_exposure_tables(
  sy$am,
  sy$cm,
  pop,
  gs,
  "neutral",
  sy$a0,
  M,
  0.66,
  0.8,
  n_boot = 0L,
  cap = 2e5
))
ch <- dplyr::filter(x0$exposure, class == "childs")
cmp <- dplyr::inner_join(
  dplyr::select(
    dplyr::filter(out$exposure, variant == "adj"),
    metric,
    level,
    scope,
    analog_set,
    pop_post = pop_point,
    n_post = n_adm1_point
  ),
  dplyr::select(ch, metric, level, scope, analog_set, pop, n_adm1),
  by = c("metric", "level", "scope", "analog_set")
)
expect(
  nrow(cmp) == nrow(ch) &&
    nrow(cmp) > 0L &&
    all(cmp$pop_post == cmp$pop) &&
    all(cmp$n_post == cmp$n_adm1),
  "pop_point / n_adm1_point = report_attrib_exposure's Childs population"
)
raw_af <- dplyr::filter(
  out$exposure,
  metric == "cases",
  variant == "raw",
  level == "continent"
)
expect(
  raw_af$pop_lo <= af$pop_lo && raw_af$pop_hi >= af$pop_hi,
  "raw interval at least as wide as adj"
)

# --- [1b] unidentified ADM1s ----------------------------------
# [1b] at noise_share_max 0.8 the noisy null is
#      unidentified (NA hits, own population variant), the rest unchanged;
cat("# [1b] unidentified ADM1s\n")
o8 <- suppressMessages(attrib_exposure_posterior_tables(
  sy$am,
  sy$cm,
  pop,
  gs,
  "neutral",
  sy$a0,
  M,
  sy$draws,
  cap = 2e5,
  noise_share_max = 0.8
))
a8 <- dplyr::filter(o8$adm1, metric == "cases")
id8 <- stats::setNames(a8$identified, a8$adm1_unit)
expect(
  isTRUE(id8[["AAA_1"]]) && !any(id8[c("AAA_2", "BBB_1", "CCC_1")]),
  "threshold 0.8: only AAA_1 (AF varying across draws) is identified; constant effects and the null are noise"
)
expect(
  near(a8$p_adj[a8$adm1_unit == "AAA_1"], p_of("AAA_1"), 1e-12) &&
    all(is.na(a8$p_adj[a8$adm1_unit != "AAA_1"])),
  "identified ADM1: same adj hit share as without the threshold; unidentified: NA"
)
e8 <- dplyr::filter(o8$exposure, metric == "cases")
un <- function(v, lv, sc = NULL) {
  x <- dplyr::filter(e8, variant == v, level == lv)
  if (!is.null(sc)) {
    x <- dplyr::filter(x, scope == sc)
  }
  x
}
expect(
  setequal(
    unique(e8$variant),
    c("adj", "raw", "unidentified", "unidentified_pass")
  ),
  "variants: adj (identified ADM1s), raw, unidentified, unidentified_pass"
)
expect(
  un("unidentified", "continent")$pop_point == 130 &&
    un("unidentified_pass", "continent")$pop_point == 70 &&
    un("unidentified", "region", "G")$pop_lo == 100,
  "unidentified = AAA_2 + BBB_1 + CCC_1 (130); on the pooled estimate AAA_2 + CCC_1 pass (70)"
)
expect(
  un("adj", "continent")$pop_point == 100 &&
    un("adj", "continent")$pop_lo == 0 &&
    un("adj", "continent")$pop_hi == 100,
  "adj: point and interval over the identified ADM1s (AAA_1)"
)
expect(
  un("adj", "continent")$pop_point +
    un("unidentified_pass", "continent")$pop_point ==
    af$pop_point,
  "adj point + unidentified_pass point = the Childs point of every ADM1"
)
expect(
  !"identified" %in% names(out$adm1) &&
    !"unidentified" %in% out$exposure$variant,
  "noise_share_max = Inf (default of the core): no flag column, no unidentified rows"
)
oi <- suppressMessages(attrib_exposure_posterior_tables(
  sy$am,
  sy$cm,
  pop,
  gs,
  "neutral",
  sy$a0,
  M,
  sy$draws,
  cap = 2e5,
  africa = "independent",
  n_mc = 2000L
))
ai <- dplyr::filter(
  oi$exposure,
  metric == "cases",
  variant == "adj",
  level == "continent"
)
expect(
  ai$pop_lo == af$pop_lo && ai$pop_hi == af$pop_hi,
  "independent combination (kept side by side): same Africa interval here (one unit per group)"
)

# --- [2] windows and guards ---------------------------------------------------
# [2] windows and guards: rows per cell x explosive guard, an explosive
#     replicate dropped only from the guarded rows, same seed -> same output;
cat("# [2] windows and explosive guards\n")
am2 <- sy$am
cm2 <- sy$cm
hot <- am2$adm1_unit == "AAA_1" & am2$run == "forecast" & am2$replicate_id == 7L
am2$cases[hot] <- 3e5
cm2$n_cases[
  cm2$country == "AAA" & cm2$run == "forecast" & cm2$replicate_id == 7L
] <- 3e5
w <- suppressMessages(attrib_exposure_posterior_windows(
  am2,
  cm2,
  pop,
  gs,
  "neutral",
  attrib_cell_table(sy$a0, sy$a0, cap = 2e5),
  M,
  sy$draws,
  n_mc = 500L
))
expect(
  setequal(
    unique(paste(w$exposure$window, w$exposure$explosive_guard)),
    c(
      "attribution TRUE",
      "attribution FALSE",
      "post_fork TRUE",
      "post_fork FALSE"
    )
  ),
  "rows per window x explosive guard"
)
mf <- function(gd) {
  w$adm1$mean_f[
    w$adm1$adm1_unit == "AAA_1" &
      w$adm1$metric == "cases" &
      w$adm1$window == "attribution" &
      w$adm1$explosive_guard == gd
  ]
}
expect(
  mf(FALSE) > mf(TRUE) + 50,
  "the explosive replicate counts only in the unguarded rows"
)
w2 <- suppressMessages(attrib_exposure_posterior_windows(
  am2,
  cm2,
  pop,
  gs,
  "neutral",
  attrib_cell_table(sy$a0, sy$a0, cap = 2e5),
  M,
  sy$draws,
  n_mc = 500L
))
expect(identical(w, w2), "same seed: identical output")

# --- [3] refusals -------------------------------------------------------------
# [3] refusals: replicates without a posterior draw, too few replicates per
#     draw.
cat("# [3] refusals\n")
err <- function(e) {
  inherits(try(suppressMessages(e), silent = TRUE), "try-error")
}
expect(
  err(attrib_exposure_posterior_tables(
    sy$am,
    sy$cm,
    pop,
    gs,
    "neutral",
    sy$a0,
    M,
    cap = 2e5,
    dplyr::filter(sy$draws, replicate_id != 5L),
    n_mc = 10L
  )),
  "refused: a replicate without a posterior draw"
)
expect(
  err(attrib_exposure_posterior_tables(
    sy$am,
    sy$cm,
    pop,
    gs,
    "neutral",
    sy$a0,
    M,
    sy$draws,
    cap = 2e5,
    min_reps = 60L,
    n_mc = 10L
  )),
  "refused: fewer than min_reps replicates per draw"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
