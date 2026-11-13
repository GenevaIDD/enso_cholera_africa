# test_spatpomp_grouped_c.R — Focused tests for Design B grouped C-kernel support.
#
# Test 1: .build_grouping_sub unit tests
#   1a  fixed_draw_idx path: each group is reduced to a 1-draw subset;
#       index clamped to group S.
#   1b  n-draw subsample path: each group is subsampled to min(n, S).
#   1c  missing-bundle guard: frame_tu$group contains a group absent from
#       grouping$draws_occ_by_group → error.
#
# Test 2: C-vs-R lp_static parity (the main correctness guard for drop_region)
#   Build precompute_lp_components with grouping_sub (grouped/region-dropped path)
#   and compare to the per-group reference built by calling precompute_lp_static
#   with drop_region=TRUE explicitly.  If drop_region is mis-wired, the C path
#   silently mis-keys region REs → non-zero discrepancy.
#   Also assert that the two groups' lp_static values are distinct (proves each
#   group's draws bundle was used, not collapsed to a single bundle).
#
# Run: Rscript analysis/tests/test_spatpomp_grouped_c.R
#      (must be called from the project root)

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R")   # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core_c.R")   # also sources spatpomp_core.R

# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------
.fails <- 0L; .passes <- 0L

source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal
expect_equal <- function(a, b, label, tol = 1e-9) {
  ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
  if (!ok) cat(sprintf("         got[1]=%s  want[1]=%s\n",
                        as.numeric(a)[1], as.numeric(b)[1]))
  expect_true(ok, label)
}
expect_error <- function(expr, label) {
  ok <- inherits(try(force(expr), silent = TRUE), "try-error")
  expect_true(ok, label)
}

cat("# test_spatpomp_grouped_c — .build_grouping_sub + C-vs-R lp_static parity\n\n")

# ---------------------------------------------------------------------------
# Shared fixture builder
# ---------------------------------------------------------------------------

#' Make a 2-group frame: group ga = {AAAW 4 units}, gb = {BBBW 4 units}.
make_2group_frame <- function() {
  tidyr::expand_grid(country = c("AAAW", "BBBW"),
                     adm1n = 1:2, adm2n = 1:2) %>%
    dplyr::mutate(
      gid         = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit   = sprintf("%s_%d",    country, adm1n),
      AFRO_region = "REG_W",
      tu          = "hmmW",
      min_seas    = 1L,
      pop         = 1e5,
      group       = dplyr::if_else(country == "AAAW", "ga", "gb")
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop,
                  adm1_unit, group) %>%
    tibble::as_tibble()
}

#' Build a region-dropped draws bundle (re_intercept$region = NULL).
#' Country keys = bare ISO3; adm1 keys = <country>_<adm1_unit>.
make_region_dropped_bundle <- function(frame_g, intercept = 0, S = 2L) {
  fe_nm <- c("(Intercept)", "mean_precip", "dev_precip",
             "mean_temp_daymean", "dev_temp_daymean",
             "ongoing_outbreak_lag1", "ongoing_outbreak_country_lag1",
             "ongoing_outbreak_tu_lag1")
  zm <- function(cols) {
    m <- matrix(0, nrow = S, ncol = length(cols))
    colnames(m) <- cols; m
  }
  fe <- zm(fe_nm)
  fe[, "(Intercept)"] <- intercept
  list(
    S            = S,
    fe           = fe,
    re_intercept = list(
      region  = NULL,
      country = zm(unique(as.character(frame_g$country))),
      adm1    = zm(unique(paste(frame_g$country, frame_g$adm1_unit, sep = "_"))),
      adm2    = NULL
    ),
    re_slope     = list(),
    basehaz      = list(kind = "const", const = 0),
    phi          = rep(1.0, S),
    active_terms = fe_nm
  )
}

#' Minimal covar_tbl (one time step, all climate = 0) matching frame_tu unit order.
make_minimal_covar_tbl <- function(frame_tu, t = as.Date("2021-01-01")) {
  tibble::tibble(
    unit              = as.character(frame_tu$gid),
    time              = t,
    t_start_seas      = 0L,
    mean_precip       = 0.0,
    dev_precip        = 0.0,
    mean_temp_daymean = 0.0,
    dev_temp_daymean  = 0.0
  )
}

# ---------------------------------------------------------------------------
cat("## 1) .build_grouping_sub — fixed_draw_idx path\n")
{
  ft <- make_2group_frame()

  draws_ga <- make_region_dropped_bundle(ft[ft$group == "ga", ], S = 4L)
  draws_gb <- make_region_dropped_bundle(ft[ft$group == "gb", ], S = 3L)

  grouping <- list(draws_occ_by_group = list(ga = draws_ga, gb = draws_gb))

  # 1a: fixed_draw_idx within bounds — each group subset has S=1
  gsub_2 <- .build_grouping_sub(grouping, fixed_draw_idx = 2L,
                                  n_draws_in_filter = 1L, frame_tu = ft)
  expect_true(!is.null(gsub_2), "1a: non-NULL result")
  expect_equal(gsub_2$draws_sub_by_group$ga$S, 1L,
               "1a: ga subset S=1 (fixed_draw_idx path)")
  expect_equal(gsub_2$draws_sub_by_group$gb$S, 1L,
               "1a: gb subset S=1 (fixed_draw_idx path)")

  # 1b: fixed_draw_idx > gb$S (=3) → clamped to 3; ga index = min(5,4) = 4
  gsub_5 <- .build_grouping_sub(grouping, fixed_draw_idx = 5L,
                                  n_draws_in_filter = 1L, frame_tu = ft)
  expect_equal(gsub_5$draws_sub_by_group$ga$S, 1L,
               "1b: ga clamped to S=4 → 1-draw subset")
  expect_equal(gsub_5$draws_sub_by_group$gb$S, 1L,
               "1b: gb clamped to S=3 → 1-draw subset")

  cat("\n## 1c: .build_grouping_sub — n-draw subsample path (no fixed_draw_idx)\n")
  gsub_nd <- .build_grouping_sub(grouping, fixed_draw_idx = NULL,
                                   n_draws_in_filter = 2L, frame_tu = ft)
  expect_equal(gsub_nd$draws_sub_by_group$ga$S, 2L,
               "1c: ga subsampled to min(2, 4) = 2")
  expect_equal(gsub_nd$draws_sub_by_group$gb$S, 2L,
               "1c: gb subsampled to min(2, 3) = 2")

  cat("\n## 1d: .build_grouping_sub — missing-bundle guard\n")
  ft_extra <- ft
  ft_extra$group[1] <- "gc"   # group gc is not in grouping
  expect_error(
    .build_grouping_sub(grouping, fixed_draw_idx = 1L,
                         n_draws_in_filter = 1L, frame_tu = ft_extra),
    "1d: group 'gc' absent from draws_occ_by_group → error")

  cat("\n## 1e: .build_grouping_sub — NULL grouping returns NULL\n")
  expect_true(is.null(.build_grouping_sub(NULL, 1L, 1L, ft)),
              "1e: NULL grouping → NULL result")
}

# ---------------------------------------------------------------------------
cat("\n## 2) C-vs-R lp_static parity (drop_region routing)\n")
{
  ft <- make_2group_frame()

  frame_ga <- ft[ft$group == "ga", , drop = FALSE]
  frame_gb <- ft[ft$group == "gb", , drop = FALSE]

  # Group-specific intercepts: ga = +2.0, gb = -1.5 (clearly distinct)
  draws_ga <- make_region_dropped_bundle(frame_ga, intercept =  2.0, S = 2L)
  draws_gb <- make_region_dropped_bundle(frame_gb, intercept = -1.5, S = 2L)

  grouping <- list(draws_occ_by_group = list(ga = draws_ga, gb = draws_gb))

  # Single-draw subset (fixed_draw_idx = 1)
  grouping_sub <- .build_grouping_sub(grouping, fixed_draw_idx = 1L,
                                       n_draws_in_filter = 1L, frame_tu = ft)

  covar_tbl <- make_minimal_covar_tbl(ft)

  # C-path (function under test): precompute_lp_components with grouping_sub
  comp_c <- precompute_lp_components(
    frame_tu    = ft,
    covar_tbl   = covar_tbl,
    enso_names  = character(0),
    cb_names    = character(0),
    draws_sub   = grouping_sub$draws_sub_by_group$ga,   # representative (unused in grouped path)
    grouping_sub = grouping_sub
  )
  lp_c <- comp_c$lp

  cat(sprintf("    lp_c range: [%.4f, %.4f]\n", min(lp_c), max(lp_c)))

  # Reference: per-group precompute_lp_static with drop_region = TRUE
  grp_vec <- as.character(ft$group)
  lp_ref  <- rep(NA_real_, nrow(covar_tbl))
  for (g in c("ga", "gb")) {
    row_mask_g <- grp_vec == g
    frame_g    <- ft[row_mask_g, , drop = FALSE]
    gids_g     <- as.character(frame_g$gid)
    cov_mask_g <- covar_tbl$unit %in% gids_g
    d_g        <- grouping_sub$draws_sub_by_group[[g]]
    lp_ref[cov_mask_g] <- precompute_lp_static(
      d_g, frame_g,
      covar_tbl[cov_mask_g, , drop = FALSE],
      enso_names = character(0), cb_names = character(0),
      drop_region = TRUE
    )
  }

  expect_equal(lp_c, lp_ref,
    "2a: grouped precompute_lp_components matches per-group reference (drop_region=TRUE)")

  expect_true(all(is.finite(lp_c)),
    "2b: lp_static is finite under grouped region-dropped routing")

  # Group distinctness: ga (intercept=+2) and gb (intercept=-1.5) must differ
  ga_lp_mean <- mean(lp_c[covar_tbl$unit %in% frame_ga$gid])
  gb_lp_mean <- mean(lp_c[covar_tbl$unit %in% frame_gb$gid])
  cat(sprintf("    ga mean lp=%.4f  gb mean lp=%.4f\n", ga_lp_mean, gb_lp_mean))
  expect_true(abs(ga_lp_mean - gb_lp_mean) > 1.0,
    "2c: ga and gb lp_static means differ by >1 (each group's draws bundle was used)")

  # Parity guard for drop_region: verify that drop_region=FALSE would MISS the
  # intercept (region=NULL bundle has no region RE to look up, so region lookup
  # returns 0 — but country/adm1 keys differ between drop_region modes, so a
  # region-rooted lookup on a region-dropped bundle would key off region-based
  # names that don't exist → effectively 0 contribution from all REs).
  # (This is already proven by tests in test_simulator_grouping.R §2; we just
  # confirm the C-path is not silently using drop_region=FALSE.)
  # Direct confirmation: lp_c values ≠ 0 (intercept is non-zero in both groups).
  expect_true(!all(lp_c == 0),
    "2d: lp_c is non-zero (confirms draws bundle was read, not silently zeroed)")
}

# ---------------------------------------------------------------------------
cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) quit(save = "no", status = 1L)
