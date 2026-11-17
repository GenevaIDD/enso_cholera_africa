# test_results_attrib_information.R — HARD: the information component of the
#   forecast attribution (results_derive.R: report_attrib_information).
#
# Usage:
#   Rscript analysis/tests/test_results_attrib_information.R
#
# Details:
#   Quick-safe: synthetic per-replicate totals and temp partitions only.
#   Pins the decomposition total = information + enso, the arm pairs, the
#   explosive guard's effect on the gap and the handoff-map check.

suppressPackageStartupMessages({
  library(dplyr)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")
cat("# information component\n")
near <- function(a, b, tol = 1e-6) {
  isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b)), na.rm = FALSE))
}

#' D draws x R replicates on units u1 (group g, AAA) and u2 (g, BBB), u3 (h, CCC);
#' per replicate three totals S, M, C with draw-dependent means.
synth <- function(D = 25L, R = 24L, seed = 3L) {
  set.seed(seed)
  n <- D * R
  rid <- seq_len(n)
  draw <- ((rid - 1L) %% D) + 1L
  mu <- exp(rnorm(D, log(2000), 0.3))
  units <- list(u1 = c("g", "AAA"), u2 = c("g", "BBB"), u3 = c("h", "CCC"))
  p <- bind_rows(lapply(names(units), function(u) {
    s <- rnbinom(n, mu = mu[draw] * 1.6, size = 3)
    m <- rnbinom(n, mu = mu[draw] * 1.3, size = 3)
    cc <- rnbinom(n, mu = mu[draw], size = 3)
    tibble::tibble(
      cell = "attribution",
      enso_group = units[[u]][1],
      unit = u,
      country = units[[u]][2],
      replicate_id = rid,
      draw = draw,
      S = s,
      M = m,
      C = cc
    )
  }))
  arm <- function(run, aset, v) {
    transmute(
      p,
      cell,
      enso_group,
      unit,
      country,
      run = run,
      analog_set = aset,
      replicate_id,
      cases = as.numeric(.data[[v]]),
      onsets = as.numeric(.data[[v]] > 1500)
    )
  }
  none <- tibble::tibble(
    cell = character(0),
    enso_group = character(0),
    run = character(0),
    analog_set = character(0),
    replicate_id = integer(0)
  )
  list(
    reps_s = bind_rows(
      arm("forecast", NA_character_, "S"),
      arm("counterfactual", "neutral", "C")
    ),
    reps_a = bind_rows(
      arm("forecast", NA_character_, "M"),
      arm("counterfactual", "neutral", "C")
    ),
    expl_s = none,
    expl_a = none,
    draws = distinct(p, cell, unit, replicate_id, draw),
    p = p
  )
}
gs <- c(g = "G", h = "H")
x <- synth()

# ---------------------------------------------------------------- [1]
# [1] decomposition: on a fixed-draw design with a displayed factual S, an
#     attribution factual M and a counterfactual C, the pooled attributable
#     numbers add (total = information + enso, gap 0) at region, country and
#     Africa level, and per draw (AN_total_d = AN_info_d + AN_enso_d);
cat("[1] decomposition\n")
t1 <- suppressWarnings(attrib_information_tables(
  x$reps_s,
  x$reps_a,
  x$expl_s,
  x$expl_a,
  x$draws,
  "neutral",
  gs,
  min_reps = 5L,
  n_mc = 200L
))
dec <- filter(t1$decomposition, explosive_guard, metric == "cases")
expect(
  nrow(dec) > 0L && all(c("continent", "region", "country") %in% dec$level),
  "decomposition rows at every level"
)
expect(
  near(dec$gap, rep(0, nrow(dec)), 1e-6),
  "pooled AN: total = information + enso (gap 0) at every scope"
)
reg <- filter(dec, level == "region", key == "g")
expect(
  near(
    reg$total,
    mean(
      tapply(
        x$p$S[x$p$enso_group == "g"] - x$p$C[x$p$enso_group == "g"],
        x$p$draw[x$p$enso_group == "g"],
        sum
      ) /
        24
    ),
    1e-9
  ),
  "region g total AN = mean over draws of the summed per-draw mean differences"
)
dr <- filter(
  t1$draws,
  explosive_guard,
  metric == "cases",
  level == "region",
  key == "g"
)
w <- tidyr::pivot_wider(
  select(dr, component, draw, an),
  names_from = component,
  values_from = an
)
expect(
  near(w$total, w$information + w$enso, 1e-9),
  "per draw: AN_total = AN_information + AN_enso"
)

# ---------------------------------------------------------------- [2]
# [2] the pairs: total = S vs C, information = S vs M (analog set
#     "information"), enso = M vs C; components and AF rows present;
cat("[2] pairs\n")
expect(
  setequal(unique(t1$scope$component), c("total", "information", "enso")),
  "three components"
)
expect(
  all(
    t1$scope$analog_set[t1$scope$component == "information"] == "information"
  ) &&
    all(t1$scope$analog_set[t1$scope$component != "information"] == "neutral"),
  "information rows carry analog set 'information', the others the counterfactual set"
)
pr <- attrib_information_pairs(
  x$reps_s,
  x$reps_a,
  x$expl_s,
  x$expl_a,
  "neutral"
)
s1 <- function(d, rn) sum(d$cases[d$run == rn])
expect(
  s1(pr$total$reps, "forecast") == sum(x$p$S) &&
    s1(pr$total$reps, "counterfactual") == sum(x$p$C),
  "total = S vs C"
)
expect(
  s1(pr$information$reps, "forecast") == sum(x$p$S) &&
    s1(pr$information$reps, "counterfactual") == sum(x$p$M),
  "information = S vs M"
)
expect(
  s1(pr$enso$reps, "forecast") == sum(x$p$M) &&
    s1(pr$enso$reps, "counterfactual") == sum(x$p$C),
  "enso = M vs C"
)
afg <- filter(
  t1$scope,
  explosive_guard,
  metric == "cases",
  level == "region",
  key == "g"
)
expect(
  all(is.finite(afg$af_adj_med)) && all(afg$af_adj_lo <= afg$af_adj_hi),
  "AF rows with finite intervals"
)

# ---------------------------------------------------------------- [3]
# [3] the explosive guard drops a replicate from the pairs it is explosive in
#     (S here: total and information, not enso), so the gap is no longer 0;
cat("[3] explosive guard\n")
ex <- tibble::tibble(
  cell = "attribution",
  enso_group = "g",
  run = "forecast",
  analog_set = NA_character_,
  replicate_id = 7L
)
t3 <- suppressWarnings(attrib_information_tables(
  x$reps_s,
  x$reps_a,
  ex,
  x$expl_a,
  x$draws,
  "neutral",
  gs,
  min_reps = 5L,
  n_mc = 200L
))
d3 <- filter(
  t3$decomposition,
  explosive_guard,
  metric == "cases",
  level == "region",
  key == "g"
)
d3u <- filter(
  t3$decomposition,
  !explosive_guard,
  metric == "cases",
  level == "region",
  key == "g"
)
expect(
  abs(d3$gap) > 0 && near(d3u$gap, 0, 1e-6),
  "a replicate explosive in S leaves the guarded gap non-zero, the unguarded one 0"
)

# ---------------------------------------------------------------- [4]
# [4] attrib_information_check_handoff(): equal handoff maps pass; a replicate
#     starting from another smoothing replicate, or a run without the map, is
#     refused.
cat("[4] handoff maps\n")
td <- tempfile("attinfo_")
dir.create(td)
cfg <- function(tag) {
  p <- file.path(td, paste0(tag, ".yml"))
  yaml::write_yaml(list(simulator = list(run_tag = tag)), p)
  p
}
part <- function(tag, map) {
  pd <- file.path(td, "sim", tag, "parts", "tu=u1", "repoff=000000")
  dir.create(pd, recursive = TRUE)
  file.create(file.path(pd, "_SUCCESS"))
  if (!is.null(map)) {
    utils::write.csv(map, file.path(pd, "handoff_map.csv"), row.names = FALSE)
  }
}
hm <- function(sm) {
  data.frame(
    tu = "u1",
    replicate_id = 1:4,
    smoothing_rep_id = sm,
    handoff_month = "2026-05-01",
    handoff_weighting = "uniform"
  )
}
mk_man <- function(fac) {
  list(
    units = list(
      u1 = list(
        groups = list("enso_g"),
        factual_config = fac,
        primary_analog_set = "neutral",
        counterfactual_sets = list(neutral = list(config = cfg("cf")))
      )
    )
  )
}
part("fs", hm(c(1L, 2L, 3L, 4L)))
part("fa", hm(c(1L, 2L, 3L, 4L)))
part("fb", hm(c(1L, 2L, 4L, 3L)))
part("fn", NULL)
ok <- tryCatch(
  {
    attrib_information_check_handoff(
      file.path(td, "sim"),
      mk_man(cfg("fs")),
      mk_man(cfg("fa")),
      "enso_g",
      "neutral"
    )
    TRUE
  },
  error = function(e) FALSE
)
expect(ok, "equal handoff maps pass")
e <- tryCatch(
  {
    attrib_information_check_handoff(
      file.path(td, "sim"),
      mk_man(cfg("fs")),
      mk_man(cfg("fb")),
      "enso_g",
      "neutral"
    )
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(e) && grepl("2 of 4 replicates differ", e, fixed = TRUE),
  "a replicate starting from another smoothing replicate is refused"
)
e <- tryCatch(
  {
    attrib_information_check_handoff(
      file.path(td, "sim"),
      mk_man(cfg("fs")),
      mk_man(cfg("fn")),
      "enso_g",
      "neutral"
    )
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(e) && grepl("no handoff_map.csv", e, fixed = TRUE),
  "a run without the handoff map is refused"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
