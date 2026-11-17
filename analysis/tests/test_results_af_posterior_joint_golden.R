# test_results_af_posterior_joint_golden.R — golden of the af_posterior /
#   attrib_exposure_posterior cores.
#
# Usage:
#   Rscript analysis/tests/test_results_af_posterior_joint_golden.R
#   Rscript analysis/tests/test_results_af_posterior_joint_golden.R --capture
#
# Details:
#   HARD, quick-safe: no fit, no data, no simulation output. The Africa row
#   is the replicate-level sum of every unit, and a scope whose Monte-Carlo
#   noise share reaches results.thresholds.af_noise_share_max is flagged
#   unidentified. The golden holds the reference cores' outputs on three
#   synthetic fixed-draw designs (written with --capture).

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})
cli <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--capture",
      action = "store_true",
      default = FALSE,
      help = "write the golden file instead of comparing"
    )
  )
))
source("analysis/utils/results_core.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/results_derive.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
GOLDEN_PATH <- "analysis/tests/golden/af_posterior_joint_golden.rds"
cat(
  "# af_posterior joint golden (mode:",
  if (cli$capture) "capture" else "compare",
  ")\n"
)

# Under the legacy settings (guard_scope = "group", noise_share_max = 1),
# every region and country row and every per-draw region / country row must
# be identical to the golden; only the Africa rows (and the Africa
# population rows of the exposure posterior) may differ. Group rows of the
# same unit under guard_scope = "unit" are checked by
# test_results_af_posterior.R.
#' Fixed-draw design of test_results_af_posterior.R, plus `shared`: groups g
#' and h on one unit (u1), so the Africa total has within-unit dependence.
synth_af <- function(D = 49L, R = 30L, nb_size = 2, seed = 1L, shared = FALSE) {
  set.seed(seed)
  af_true <- list(g = rnorm(D, 0.45, 0.12), h = rnorm(D, 0.10, 0.05))
  mu_f <- list(
    g = exp(rnorm(D, log(3000), 0.3)),
    h = exp(rnorm(D, log(800), 0.3))
  )
  n <- D * R
  rid <- seq_len(n)
  draw <- ((rid - 1L) %% D) + 1L
  units <- if (shared) {
    list(
      a = c("g", "u1", "AAA", 0.6),
      b = c("g", "u2", "BBB", 0.4),
      c = c("h", "u1", "CCC", 1)
    )
  } else {
    list(
      a = c("g", "u1", "AAA", 0.6),
      b = c("g", "u2", "BBB", 0.4),
      c = c("h", "u3", "CCC", 1)
    )
  }
  rows <- lapply(units, function(x) {
    grp <- x[1]
    s <- as.numeric(x[4])
    f <- rnbinom(n, mu = s * mu_f[[grp]][draw], size = nb_size)
    cf <- rnbinom(
      n,
      mu = s * mu_f[[grp]][draw] * (1 - af_true[[grp]][draw]),
      size = nb_size
    )
    tibble::tibble(
      cell = "w",
      enso_group = grp,
      unit = x[2],
      country = x[3],
      replicate_id = rid,
      draw = draw,
      f = f,
      cf = cf
    )
  })
  p <- dplyr::bind_rows(rows)
  arm <- function(run, aset, v) {
    dplyr::transmute(
      p,
      cell,
      enso_group,
      unit,
      country,
      run = run,
      analog_set = aset,
      replicate_id,
      cases = .data[[v]],
      onsets = as.numeric(.data[[v]] > 0)
    )
  }
  reps <- dplyr::bind_rows(
    arm("forecast", NA_character_, "f"),
    arm("counterfactual", "neutral", "cf")
  )
  # Replicates 1-3 of group g explosive in the forecast arm (exercises the guard).
  ex <- dplyr::distinct(
    dplyr::filter(
      reps,
      enso_group == "g",
      replicate_id %in% 1:3,
      run == "forecast"
    ),
    cell,
    enso_group,
    run,
    analog_set,
    replicate_id
  )
  list(
    reps = reps,
    explosive = ex,
    draws = dplyr::distinct(p, cell, unit, replicate_id, draw)
  )
}

#' Fixed-draw design of test_results_attrib_exposure_posterior.R.
synth_ax <- function(D = 49L, R = 30L, nb_size = 2, seed = 1L) {
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
  rows <- lapply(names(spec), function(a) {
    s <- spec[[a]]
    draw <- ((rid - 1L) %% D) + 1L
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
  list(am = am, cm = cm, draws = draws, a0 = a0)
}

gs <- c(g = "G", h = "H")
pop <- tibble::tibble(
  adm1_unit = c("AAA_1", "AAA_2", "BBB_1", "CCC_1"),
  pop = c(100, 40, 60, 30)
)

#' The cores under the legacy settings. Before the change the extra arguments
#' do not exist, so they are passed only when the function accepts them.
legacy_args <- function(fn) {
  a <- list()
  if ("africa" %in% names(formals(fn))) {
    a$africa <- "independent"
  }
  if ("guard_scope" %in% names(formals(fn))) {
    a$guard_scope <- "group"
  }
  if ("noise_share_max" %in% names(formals(fn))) {
    a$noise_share_max <- Inf
  }
  a
}
run_af <- function(sy) {
  do.call(
    af_posterior_tables,
    c(
      list(
        sy$reps,
        sy$explosive,
        sy$draws,
        "neutral",
        gs,
        guards = c(TRUE, FALSE)
      ),
      legacy_args(af_posterior_tables)
    )
  )
}
run_ax <- function(sy) {
  suppressMessages(do.call(
    attrib_exposure_posterior_tables,
    c(
      list(
        sy$am,
        sy$cm,
        pop,
        gs,
        "neutral",
        sy$a0,
        c(0.25, 0.5),
        sy$draws,
        cap = 2e5
      ),
      legacy_args(attrib_exposure_posterior_tables)
    )
  ))
}

cur <- list(
  af_split = run_af(synth_af()),
  af_shared = run_af(synth_af(shared = TRUE)),
  ax = run_ax(synth_ax())
)

if (cli$capture) {
  saveRDS(cur, GOLDEN_PATH)
  cat("captured", GOLDEN_PATH, "\n")
  quit(status = 0L)
}
gold <- readRDS(GOLDEN_PATH)

#' Rows of a scope / draws table outside Africa, the columns the golden had,
#' sorted, so new columns (af_identified, ...) do not break the comparison.
non_africa <- function(x, cols) {
  x <- as.data.frame(x)
  if ("level" %in% names(x)) {
    x <- x[x$level != "continent", , drop = FALSE]
  }
  x <- x[, intersect(cols, names(x)), drop = FALSE]
  x[do.call(order, unname(as.list(x))), , drop = FALSE] |> `rownames<-`(NULL)
}
for (k in c("af_split", "af_shared")) {
  for (tb in c("scope", "draws")) {
    g <- gold[[k]][[tb]]
    n <- cur[[k]][[tb]]
    keep <- setdiff(names(g), "n_mc")
    expect(
      isTRUE(all.equal(
        non_africa(n, keep),
        non_africa(g, keep),
        tolerance = 0
      )),
      sprintf(
        "%s$%s: region and country rows identical to the golden (legacy settings)",
        k,
        tb
      )
    )
  }
}
ga <- non_africa(gold$ax$adm1, names(gold$ax$adm1))
na <- non_africa(cur$ax$adm1, names(gold$ax$adm1))
expect(
  isTRUE(all.equal(na, ga, tolerance = 0)),
  "ax$adm1: per-ADM1 table identical to the golden (legacy settings)"
)
gx <- non_africa(gold$ax$exposure, names(gold$ax$exposure))
nx <- non_africa(cur$ax$exposure, names(gold$ax$exposure))
expect(
  nrow(gx) == nrow(nx) &&
    isTRUE(all.equal(nx$pop_point, gx$pop_point, tolerance = 0)) &&
    isTRUE(all.equal(nx$n_adm1_point, gx$n_adm1_point, tolerance = 0)),
  "ax$exposure: region / country point populations identical to the golden"
)
expect(identical(cur$ax$n_min, gold$ax$n_min), "ax: n_min unchanged")

# The legacy path ("independent", "group", Inf) is the old code kept side by
# side: every output, Africa included, must equal the golden exactly.
for (k in c("af_split", "af_shared")) {
  expect(
    identical(cur[[k]], gold[[k]]),
    sprintf("%s: legacy path identical to the golden, Africa included", k)
  )
}
expect(
  identical(cur$ax, gold$ax),
  "ax: legacy path identical to the golden, Africa included"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
