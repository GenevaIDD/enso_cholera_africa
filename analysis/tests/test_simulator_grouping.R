# test_simulator_grouping.R — the country-grouping simulator building blocks
#   (drop_region routing).
#
# Usage:
#   Rscript analysis/tests/test_simulator_grouping.R
#
# Details:
#   Synthetic fixtures, no real fit object: region-dropped RE routing and the
#   per-group artefact-discovery helpers. Covers re_level_key() and
#   evaluate_lp() with drop_region = TRUE (simulator_core.R),
#   .lags_required_from_draws(), .discover_group_fit() and
#   check_frame_vs_draws() (utils/simulations.R). Not registered in
#   analysis/tests/run_all.R.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/dlnm.R") # .LATENT_CB_VARS -- required by .lags_required_from_draws()
source("analysis/utils/simulator_core.R")
source("analysis/utils/simulations.R")

# --------------------------------------------------------------------------
# Test harness
# --------------------------------------------------------------------------

.fails <- 0L
.passes <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal
expect_error <- function(expr, label) {
  ok <- inherits(try(force(expr), silent = TRUE), "try-error")
  expect_true(ok, label)
}

cat("# simulator country-grouping (C2 / drop_region) tests\n\n")

# --------------------------------------------------------------------------
# Fixture: 2 ADM1 × 2 ADM2 in country AAAW, region REG_W (one TU = hmmW).
# Mirrors fixture_frame() in test_simulator_kernel.R.
# --------------------------------------------------------------------------

fixture_frame <- function() {
  expand.grid(
    tu = c("hmmW", "hmmCE"),
    country = c("AAA", "BBB"),
    adm1n = c(1, 2),
    adm2n = c(1, 2),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      country = paste0(country, substr(tu, 4, 4)),
      gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
      adm1_unit = sprintf("%s_%d", country, adm1n),
      AFRO_region = paste0("REG_", substr(tu, 4, 5)),
      min_seas = 1L,
      pop = 1e5
    ) %>%
    dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
    tibble::as_tibble()
}

frame_tu <- fixture_frame() %>% dplyr::filter(tu == "hmmW")

# --------------------------------------------------------------------------
cat("## 1) re_level_key(drop_region=TRUE) drops the region prefix\n")
{
  # Ungrouped (region-rooted) keys paste region first.
  k_region_full <- re_level_key(frame_tu, "adm2", drop_region = FALSE)
  expect_true(
    all(grepl("^REG_W_", k_region_full)),
    "ungrouped adm2 keys are region-rooted (REG_W_...)"
  )

  # Grouped (region-dropped) keys root at country.
  expect_true(
    all(is.na(re_level_key(frame_tu, "region", drop_region = TRUE))),
    "region bucket -> NA under drop_region"
  )
  expect_equal(
    re_level_key(frame_tu, "country", drop_region = TRUE),
    frame_tu$country,
    "country bucket -> bare country"
  )
  expect_equal(
    re_level_key(frame_tu, "adm1", drop_region = TRUE),
    paste(frame_tu$country, frame_tu$adm1_unit, sep = "_"),
    "adm1 bucket -> country_adm1unit"
  )
  expect_equal(
    re_level_key(frame_tu, "adm2", drop_region = TRUE),
    paste(frame_tu$country, frame_tu$adm1_unit, frame_tu$gid, sep = "_"),
    "adm2 bucket -> country_adm1unit_gid"
  )
}

# --------------------------------------------------------------------------
cat(
  "\n## 2) evaluate_lp(drop_region=TRUE) finds REs in a region-dropped bundle\n"
)
{
  N <- nrow(frame_tu)
  adm2_keys <- re_level_key(frame_tu, "adm2", drop_region = TRUE)
  target_gid <- "AAAW_1_1"
  target_key <- adm2_keys[frame_tu$gid == target_gid]
  re_adm2 <- matrix(
    0,
    nrow = 2,
    ncol = length(adm2_keys),
    dimnames = list(NULL, adm2_keys)
  )
  re_adm2[, target_key] <- 2.0
  draws <- list(
    S = 2L,
    fe = matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "Intercept")),
    re_intercept = list(
      region = NULL,
      country = NULL,
      adm1 = NULL,
      adm2 = re_adm2
    ),
    re_slope = list(),
    basehaz = list(kind = "none"),
    phi = c(5, 5),
    active_terms = "Intercept"
  )
  xrow <- list(
    mean_precip = rep(0, N),
    dev_precip = rep(0, N),
    mean_temp_daymean = rep(0, N),
    dev_temp_daymean = rep(0, N)
  )
  lp <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu,
    drop_region = TRUE
  )
  exp_lp <- ifelse(frame_tu$gid == target_gid, 2.0, 0.0)
  expect_equal(
    lp,
    exp_lp,
    "region-dropped adm2 intercept RE adds to LP for its gid only"
  )

  # evaluate_lp() uses a dual-key lookup: primary keys per
  # drop_region, falling back to the alternate format for misses (grouped
  # brms_cb fits mix bare-country and AFRO-prefixed keys). So the same bundle
  # with drop_region=FALSE still finds the RE via the fallback.
  lp0 <- evaluate_lp(
    draws,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu,
    drop_region = FALSE
  )
  expect_equal(
    lp0,
    exp_lp,
    "region-rooted lookup falls back to the region-dropped colnames"
  )

  # Routing still depends on drop_region: when both key formats are present,
  # the primary format for the given drop_region wins.
  rooted_key <- re_level_key(frame_tu, "adm2", drop_region = FALSE)[
    frame_tu$gid == target_gid
  ]
  re_both <- cbind(
    re_adm2,
    matrix(5.0, nrow = 2, ncol = 1, dimnames = list(NULL, rooted_key))
  )
  draws_both <- draws
  draws_both$re_intercept$adm2 <- re_both
  lp_dr <- evaluate_lp(
    draws_both,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu,
    drop_region = TRUE
  )
  lp_rr <- evaluate_lp(
    draws_both,
    1L,
    xrow,
    o = rep(0L, N),
    o_ct = rep(0L, N),
    o_tu = rep(0L, N),
    u_days = rep(NA_real_, N),
    frame_tu = frame_tu,
    drop_region = FALSE
  )
  expect_equal(
    lp_dr,
    exp_lp,
    "both formats present: drop_region=TRUE prefers region-dropped key"
  )
  expect_equal(
    lp_rr,
    ifelse(frame_tu$gid == target_gid, 5.0, 0.0),
    "both formats present: drop_region=FALSE prefers region-rooted key"
  )
}

# --------------------------------------------------------------------------
cat("\n## 3) .lags_required_from_draws: enso FE + cb_meta lag expansion\n")
{
  draws <- list(
    fe = matrix(
      0,
      nrow = 1,
      ncol = 3,
      dimnames = list(NULL, c("Intercept", "enso_lag3", "mean_precip"))
    )
  )
  # No cb_meta: only the enso FE lag column is required.
  expect_equal(
    .lags_required_from_draws(draws, NULL),
    "enso_lag3",
    "enso FE lag column detected; no cb expansion without cb_meta"
  )

  # cb_meta with enso (max_lag 7) + a climate var (max_lag 2).
  cb_meta <- list(
    enso = list(lag = c(0, 7)),
    mean_precip = list(lag = c(0, 2))
  )
  got <- .lags_required_from_draws(draws, cb_meta)
  exp <- unique(c(
    "enso_lag3",
    sprintf("enso_lag%d", 0:7),
    "mean_precip",
    "mean_precip_lag1",
    "mean_precip_lag2"
  ))
  expect_equal(
    sort(got),
    sort(exp),
    "cb_meta expands enso_lag0..7 and climate var + _lag1.._lag2"
  )
}

# --------------------------------------------------------------------------
cat("\n## 4) .discover_group_fit: glob by _grp- token, exactly-one semantics\n")
{
  tdir <- tempfile("grpfits_")
  dir.create(tdir)
  on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
  src <- "afro_subregion"
  # Occurrence fit for group 'east' + sidecars that must NOT match.
  occ_east <- file.path(
    tdir,
    sprintf(
      "outbreaks_survival_fit_base-1990_2019_grp-%s-east_adm-adm2.rds",
      src
    )
  )
  file.create(c(
    occ_east,
    file.path(
      tdir,
      sprintf(
        "outbreaks_survival_re_stats_base-1990_2019_grp-%s-east_adm-adm2.rds",
        src
      )
    ),
    file.path(
      tdir,
      sprintf(
        "outbreaks_survival_coef_summary_base-1990_2019_grp-%s-east_adm-adm2.rds",
        src
      )
    ),
    file.path(
      tdir,
      sprintf(
        "outbreak_size_fit_base-1990_2019_grp-%s-west_fam-negbinomial_adm-ADM1.rds",
        src
      )
    )
  ))

  expect_equal(
    .discover_group_fit(tdir, src, "east", "occ"),
    occ_east,
    "occ glob matches exactly the _fit_ file (excludes re_stats/coef_summary)"
  )
  # The size track is no longer simulated: "size" is an unknown track and
  # errors even though an outbreak_size _fit_ file is present.
  expect_error(
    .discover_group_fit(tdir, src, "west", "size"),
    "which = 'size' -> unknown-track error (no outbreak_size fallback)"
  )
  expect_error(
    .discover_group_fit(tdir, src, "east", "bogus"),
    "unknown which -> error"
  )
  expect_error(
    .discover_group_fit(tdir, src, "central", "occ"),
    "no match -> error (group 'central' absent)"
  )

  # Two occ fits for the same group -> ambiguity error.
  file.create(file.path(
    tdir,
    sprintf(
      "outbreaks_survival_fit_base-1980_2009_grp-%s-east_adm-adm2.rds",
      src
    )
  ))
  expect_error(
    .discover_group_fit(tdir, src, "east", "occ"),
    ">1 match -> ambiguity error"
  )
}

# --------------------------------------------------------------------------
cat("\n## 5) check_frame_vs_draws(drop_region=TRUE) uses region-dropped keys\n")
{
  # Build an adm2 RE matrix that only contains the region-DROPPED key for one
  # gid. Under drop_region=TRUE the restrict policy must keep that gid and drop
  # the rest. Under drop_region=FALSE the (region-rooted) lookup would miss the
  # seen gid too and drop everything -> proves the flag is wired through.
  keep_gid <- "AAAW_1_1"
  seen_key <- paste(
    frame_tu$country[frame_tu$gid == keep_gid][1],
    frame_tu$adm1_unit[frame_tu$gid == keep_gid][1],
    keep_gid,
    sep = "_"
  )
  re_adm2 <- matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, seen_key))
  draws <- list(re_intercept = list(adm2 = re_adm2, adm1 = NULL))
  scope <- list(admin_level = "ADM2", unseen_adm1 = "restrict")

  kept <- suppressWarnings(suppressMessages(
    check_frame_vs_draws(frame_tu, draws, scope, drop_region = TRUE)
  ))
  expect_equal(
    sort(unique(kept$gid)),
    keep_gid,
    "restrict keeps only the gid present in the region-dropped RE matrix"
  )

  # Same matrix, wrong (region-rooted) key space -> everything unseen -> error
  # because restrict removes all rows.
  expect_error(
    suppressWarnings(suppressMessages(
      check_frame_vs_draws(frame_tu, draws, scope, drop_region = FALSE)
    )),
    "region-rooted lookup misses every gid -> restrict empties frame -> error"
  )
}

# --------------------------------------------------------------------------
cat(
  "\n## 6) simulate_one_tu grouped routing: high-LP group fires, low-LP group silent\n"
)
{
  # Frame: hmmW TU — AAAW → "ga" (high intercept) / BBBW → "gb" (low intercept).
  ft <- fixture_frame() %>%
    dplyr::filter(tu == "hmmW") %>%
    dplyr::mutate(group = dplyr::if_else(country == "AAAW", "ga", "gb"))

  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = ft$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = 0L,
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  horizon <- list(start = "2024-01-01", end = "2024-06-30")

  # Region-dropped draws bundle (region=NULL; country-rooted adm1 keys matching
  # re_level_key(frame_g, "adm1", drop_region=TRUE) = paste(country, adm1_unit)).
  make_grouped_bundle <- function(frame_g, intercept, phi = 1.0, S = 2) {
    fe_nm <- c(
      "(Intercept)",
      "mean_precip",
      "dev_precip",
      "mean_temp_daymean",
      "dev_temp_daymean",
      "ongoing_outbreak_lag1",
      "ongoing_outbreak_country_lag1",
      "ongoing_outbreak_tu_lag1"
    )
    zm <- function(cols) {
      m <- matrix(0, S, length(cols))
      colnames(m) <- cols
      m
    }
    fe <- zm(fe_nm)
    fe[, "(Intercept)"] <- intercept
    list(
      S = S,
      fe = fe,
      re_intercept = list(
        region = NULL,
        country = zm(unique(frame_g$country)),
        adm1 = zm(unique(paste(frame_g$country, frame_g$adm1_unit, sep = "_"))),
        adm2 = NULL
      ),
      re_slope = list(),
      basehaz = list(kind = "const", const = 0),
      phi = rep(phi, S),
      active_terms = fe_nm
    )
  }

  frame_ga <- ft %>% dplyr::filter(group == "ga")
  frame_gb <- ft %>% dplyr::filter(group == "gb")

  ga_occ <- make_grouped_bundle(frame_ga, intercept = 20, phi = 5)
  gb_occ <- make_grouped_bundle(frame_gb, intercept = -1e6, phi = 5)

  res <- simulate_one_tu(
    tu = "hmmW",
    frame_tu = ft,
    panel_tu = panel,
    draws_occ = ga_occ, # positional arg; per-group routing overrides via grouping=
    replicate_ids = 1L:3L,
    seed = 11L,
    output_part = tempfile("sim_grp_routing_"),
    horizon = horizon,
    dry_run = TRUE,
    grouping = list(
      draws_occ_by_group = list(ga = ga_occ, gb = gb_occ),
      cb_meta_occ_by_group = list(ga = NULL, gb = NULL)
    )
  )

  ga_gids <- unique(frame_ga$gid)
  gb_gids <- unique(frame_gb$gid)

  gb_monthly <- res$monthly_summary %>% dplyr::filter(gid %in% gb_gids)
  expect_true(
    all(gb_monthly$Y == 0L),
    "routing: gb-group (LP=-1e6) shows Y=0 in all months/reps"
  )

  ga_monthly <- res$monthly_summary %>% dplyr::filter(gid %in% ga_gids)
  expect_true(
    all(ga_monthly$Y == 1L),
    "routing: ga-group (LP=+20) shows Y=1 in all months/reps"
  )

  expect_true(
    nrow(res$event_table) > 0L && all(res$event_table$gid %in% ga_gids),
    "routing: event_table non-empty and every row from ga-group (AAAW)"
  )
}

# --------------------------------------------------------------------------
cat(
  "\n## 7) simulate_one_tu single-group equivalence: grouped == ungrouped (zero RE)\n"
)
{
  # One group "all" wraps the full-frame bundle. All RE values are zero, so the
  # region-dropped (grouped) and region-rooted (ungrouped) RE lookups both add 0
  # → LP = FE intercept in both cases → identical p → same Y with same seed.
  ft <- fixture_frame() %>%
    dplyr::filter(tu == "hmmW") %>%
    dplyr::mutate(group = "all")

  months <- as.Date(sprintf("2024-%02d-01", 1:6))
  panel <- tidyr::expand_grid(gid = ft$gid, t_start = months) %>%
    dplyr::mutate(
      t_start_seas = 0L,
      mean_precip = 0,
      dev_precip = 0,
      mean_temp_daymean = 0,
      dev_temp_daymean = 0
    )
  horizon <- list(start = "2024-01-01", end = "2024-06-30")

  bun_occ <- make_draws_bundle_synthetic(
    ft,
    S = 4,
    fe_overrides = c(`(Intercept)` = -2)
  )

  common_args <- list(
    tu = "hmmW",
    frame_tu = ft,
    panel_tu = panel,
    draws_occ = bun_occ,
    replicate_ids = 1L:4L,
    seed = 42L,
    output_part = tempfile("sim_equiv_"),
    horizon = horizon,
    dry_run = TRUE
  )

  res_ug <- do.call(simulate_one_tu, common_args)

  res_g <- do.call(
    simulate_one_tu,
    c(
      common_args,
      list(
        grouping = list(
          draws_occ_by_group = list(all = bun_occ),
          cb_meta_occ_by_group = list(all = NULL)
        )
      )
    )
  )

  expect_equal(
    res_ug$monthly_summary$Y,
    res_g$monthly_summary$Y,
    "equivalence: Y identical between ungrouped and single-group grouped (zero RE)"
  )
  expect_equal(
    nrow(res_ug$event_table),
    nrow(res_g$event_table),
    "equivalence: event_table row count identical (zero RE)"
  )
}

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(save = "no", status = 1L)
}
