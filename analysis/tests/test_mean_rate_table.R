# test_mean_rate_table.R — unit tests for the country x admin-level
# mean_rate table (analysis/utils.R: buildMeanRateTable() / resolveMeanRate()
# / fitTrainingMeanRate()) and its fallback ladder.
#
# Synthetic-fixture tests only -- no real outbreak-size data is read.
#
# Run: Rscript analysis/tests/test_mean_rate_table.R

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R")   # transitive dependency of provenance.R (utils/README.md)
source("analysis/utils/paths.R")   # transitive dependency of provenance.R (utils/README.md)
source("analysis/utils/provenance.R")   # buildMeanRateTable()/fitTrainingMeanRate()/resolveMeanRate() -- used directly

.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

cat("# mean_rate table (buildMeanRateTable / resolveMeanRate) tests\n\n")

# --------------------------------------------------------------------------
# Shared fixture: raw outbreak-size-analysis-data-shaped tibble, using real
# ENSO-registry iso3 (analysis/groupings/enso_grouping.csv) so the T3 group
# path is genuinely exercised, not just a synthetic label. KEN/UGA/TZA/RWA
# are all `enso_east`.
#   KEN, UGA: training rows at both ADM1 and ADM2 -> own cells at T1; feed
#             the country x level_factor estimation (needs >=2 countries
#             spanning >1 level).
#   TZA:      training rows at ADM1 only -> its ADM2 query resolves at T2
#             (country present, level absent -> parent + level_factor).
#   BIGN:     a large (n=500) ADM1 cell -> used for the w -> 1 shrinkage
#             check (not in the ENSO registry; irrelevant to the group path).
#   RWA:      absent from training, present in the registry -> T3.
#   ZZZ:      absent from training AND the registry -> T4.
# --------------------------------------------------------------------------
make_raw_dat <- function() {
  mk_cell <- function(country, level, n, base_rate) {
    tibble::tibble(
      country         = country,
      admin_level     = level,
      pop             = round(stats::runif(n, 2e4, 8e4))
    ) %>%
      dplyr::mutate(suspected_cases = round(pop * base_rate * stats::runif(n, 0.5, 1.5)))
  }
  set.seed(42)
  dplyr::bind_rows(
    mk_cell("KEN",  "ADM1", 5,   5e-4),
    mk_cell("KEN",  "ADM2", 5,   1.2e-3),
    mk_cell("UGA",  "ADM1", 5,   4e-4),
    mk_cell("UGA",  "ADM2", 5,   1.0e-3),
    mk_cell("TZA",  "ADM1", 3,   6e-4),
    mk_cell("BIGN", "ADM1", 500, 5e-4)
  ) %>%
    dplyr::mutate(TL = as.Date("2015-01-01"), TR = as.Date("2015-06-01"))
}

raw_dat <- make_raw_dat()
raw_rds <- tempfile(fileext = ".rds")
saveRDS(raw_dat, raw_rds)

opt <- list(size = list(min_pop = 10, max_attack_rate = 0.10,
                          mean_rate = list(estimator = "geometric")))

mrt <- suppressWarnings(buildMeanRateTable(
  opt, path = raw_rds,
  enso_grouping_path = "analysis/groupings/enso_grouping.csv"
))

cat("## 1) buildMeanRateTable: schema + own-cell tiers (T1/T2)\n")
{
  expect_true(all(c("cell", "group", "global", "meta") %in% names(mrt)),
              "returns list(cell, group, global, meta)")
  expect_true(all(c("country", "admin_level", "mean_rate", "mr_tier", "n_cell",
                     "w_shrink", "mr_raw", "mr_parent", "estimator", "k",
                     "level_factor") %in% names(mrt$cell)),
              "cell table has the documented schema")

  ken1 <- dplyr::filter(mrt$cell, country == "KEN", admin_level == "ADM1")
  ken2 <- dplyr::filter(mrt$cell, country == "KEN", admin_level == "ADM2")
  tza2 <- dplyr::filter(mrt$cell, country == "TZA", admin_level == "ADM2")
  expect_equal(ken1$mr_tier, "T1", "KEN/ADM1 (own data) -> T1")
  expect_equal(ken2$mr_tier, "T1", "KEN/ADM2 (own data) -> T1")
  expect_equal(tza2$mr_tier, "T2", "TZA/ADM2 (country present, level absent) -> T2")
  expect_true(is.na(tza2$n_cell), "TZA/ADM2 has no own-cell rows (n_cell NA)")
}

cat("\n## 2) resolveMeanRate: full ladder T1/T2/T3/T4\n")
{
  q <- resolveMeanRate(
    country     = c("KEN",  "TZA",  "RWA",    "RWA",    "ZZZ"),
    admin_level = c("ADM1", "ADM2", "ADM1",   "ADM3",   "ADM1"),
    mrt
  )
  expect_equal(q$mr_tier[1], "T1", "KEN/ADM1 -> T1 (own cell)")
  expect_equal(q$mr_tier[2], "T2", "TZA/ADM2 -> T2 (country present, level absent)")
  expect_equal(q$mr_tier[3], "T3-own",
               "RWA/ADM1 -> T3-own (country absent, but enso_east has ADM1 rows)")
  expect_equal(q$mr_tier[4], "T3-parent",
               "RWA/ADM3 -> T3-parent (enso_east has no ADM3 rows -> group parent + level_factor)")
  expect_equal(q$mr_tier[5], "T4",
               "ZZZ/ADM1 -> T4 (absent from training AND the ENSO registry)")
  expect_true(all(is.finite(q$mean_rate) & q$mean_rate > 0),
              "every resolved mean_rate is a finite positive number")
}

cat("\n## 3) resolveMeanRate: T5 hard error on an unresolvable admin_level\n")
{
  err <- tryCatch(resolveMeanRate("RWA", "ADM4", mrt),
                   error = function(e) conditionMessage(e))
  expect_true(!is.null(err) && grepl("T5", err),
              "unrecognised admin_level (no cell/group/global row) -> T5 hard error")

  # A genuinely unresolvable pair should error even when the code path is
  # reached via the initial cell join rather than the fallback join.
  err2 <- tryCatch(resolveMeanRate("ZZZ", "ADM4", mrt),
                    error = function(e) conditionMessage(e))
  expect_true(!is.null(err2) && grepl("T5", err2),
              "country AND admin_level both unresolvable -> T5 hard error")
}

cat("\n## 4) EB shrinkage: w = n/(n+k) increases toward 1 as n grows\n")
{
  k    <- mrt$meta$k
  tza1 <- dplyr::filter(mrt$cell, country == "TZA",  admin_level == "ADM1")
  bign <- dplyr::filter(mrt$cell, country == "BIGN", admin_level == "ADM1")

  expect_true(is.finite(k) && k > 0, "k (EB shrinkage constant) is finite and positive")
  expect_equal(tza1$w_shrink, tza1$n_cell / (tza1$n_cell + k),
               "TZA/ADM1 w_shrink matches the n/(n+k) formula exactly")
  expect_equal(bign$w_shrink, bign$n_cell / (bign$n_cell + k),
               "BIGN/ADM1 w_shrink matches the n/(n+k) formula exactly")
  expect_true(bign$w_shrink > tza1$w_shrink,
              "larger n_cell (500 vs 3) shrinks less (w closer to 1)")
  expect_true(bign$w_shrink > 0.98,
              "w -> 1 as n -> infinity (n=500 essentially untouched by shrinkage)")
}

cat("\n## 5) Global rescaling of the table leaves fitted mu invariant (intercept absorption)\n")
{
  # mu = pop * mean_rate * exp(lp). A constant multiplicative rescale of every
  # mean_rate in the table is exactly offset by a matching additive shift of
  # -log(c) to lp -- the offset's level is unidentifiable on its own; only its
  # per-unit CONTRASTS matter (this is the design note's "global level is
  # unidentifiable" finding). We rescale a copy of `mrt` by a constant c and
  # confirm mu is unchanged when lp is shifted by -log(c) to compensate.
  c_scale <- 3.7
  mrt_scaled <- mrt
  mrt_scaled$cell$mean_rate   <- mrt$cell$mean_rate   * c_scale
  mrt_scaled$cell$mr_raw      <- mrt$cell$mr_raw      * c_scale
  mrt_scaled$cell$mr_parent   <- mrt$cell$mr_parent   * c_scale
  mrt_scaled$group$mean_rate  <- mrt$group$mean_rate  * c_scale
  mrt_scaled$global$mean_rate <- mrt$global$mean_rate * c_scale

  qry_country <- c("KEN",  "TZA",  "RWA",  "ZZZ")
  qry_level   <- c("ADM1", "ADM2", "ADM1", "ADM1")
  pop         <- c(5e4,    3e4,    2e4,    6e4)
  lp          <- c(0.10,  -0.20,   0.05,   0.30)

  before <- resolveMeanRate(qry_country, qry_level, mrt)
  after  <- resolveMeanRate(qry_country, qry_level, mrt_scaled)

  expect_true(all(abs(after$mean_rate / before$mean_rate - c_scale) < 1e-9),
              "every resolved mean_rate scaled by exactly c_scale")
  expect_equal(after$mr_tier, before$mr_tier,
               "tier assignment is unaffected by a global rescale")

  mu_before <- pop * before$mean_rate * exp(lp)
  mu_after  <- pop * after$mean_rate  * exp(lp - log(c_scale))
  expect_true(all(abs(mu_before - mu_after) < 1e-8 * abs(mu_before)),
              "mu = pop * mean_rate * exp(lp) is invariant under a global rescale + compensating intercept shift")
}

cat("\n## 6) fitTrainingMeanRate: hard error on a per-row (non-constant) mean_rate\n")
{
  fit_const <- list(data = data.frame(mean_rate = rep(2.5e-4, 10)))
  mr <- fitTrainingMeanRate(fit_const)
  expect_equal(mr, 2.5e-4, "constant mean_rate column -> returns that scalar unchanged")

  fit_varying <- list(data = data.frame(mean_rate = c(1e-4, 2e-4, 3e-4, 1e-4)))
  err <- tryCatch({ fitTrainingMeanRate(fit_varying); NULL },
                  error = function(e) conditionMessage(e))
  expect_true(!is.null(err),
              "per-row (non-constant) mean_rate -> fitTrainingMeanRate() now hard-errors (D5)")
  expect_true(!is.null(err) && grepl("not constant", err),
              "error message explains why (fit trained with a per-unit country x admin-level rate)")

  fit_no_data <- list(data = NULL)
  expect_true(is.na(fitTrainingMeanRate(fit_no_data)),
              "missing fit$data -> NA (unchanged legacy behaviour, not an error)")
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0) quit(status = 1)
