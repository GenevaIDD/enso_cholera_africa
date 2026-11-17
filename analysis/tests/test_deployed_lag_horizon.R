# test_deployed_lag_horizon.R — HARD: the deployed best-lag lookup and the
#   forecast / hindcast horizon rule (occ / dur / mon tracks; no size track).
#
# Usage:
#   Rscript analysis/tests/test_deployed_lag_horizon.R
#
# Details:
#   Covers utils/paths.R::read_deployed_lag(), make_forecast_configs.R
#   generate_forecast_configs(), make_hindcast_configs.R
#   generate_hindcast_configs() and the event generator's season ends.
#   Tempdir fixtures only (empty RDS files named like deployed fits): no fit,
#   no data.

suppressPackageStartupMessages({
  library(yaml)
  library(lubridate)
})

source("analysis/utils/make_forecast_configs.R") # sources config/base/paths/enso/provenance/anchor
source("analysis/utils/make_hindcast_configs.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# deployed best-lag lookup + occ/dur/mon horizon rule\n")

errmsg <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

root <- tempfile("deployed_lag_")
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

GROUPING <- "enso_grouping"
fit_name <- function(prefix, kind, group, lag, extra = "", ep = "") {
  sprintf(
    "%s_%s_base-1990_2019_grp-%s-%s%s_sl%d%s.rds",
    prefix,
    kind,
    GROUPING,
    group,
    extra,
    lag,
    ep
  )
}
touch <- function(dir, name, opt = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  obj <- list()
  if (!is.null(opt)) {
    attr(obj, "opt") <- opt
  }
  saveRDS(obj, file.path(dir, name))
  invisible(file.path(dir, name))
}

# --- fixture dir: enso_east has all three tracks (+ noise), enso_west lacks mon
fits <- file.path(root, "fits")
touch(
  fits,
  fit_name(
    "outbreaks_survival_brms_cb",
    "fit",
    "enso_east",
    5L,
    "_bh-cyclic",
    ""
  ),
  opt = list(occurrence = list(offset_baserate = FALSE))
)
touch(fits, fit_name("outbreaks_survival_brms_cb", "cb_meta", "enso_east", 1L)) # not a fit
touch(
  fits,
  fit_name("outbreaks_survival_brms_cb", "fit", "enso_east", 2L, "", "_ep-gs")
) # gid_season
touch(
  fits,
  fit_name(
    "outbreak_duration_brms_cb",
    "fit",
    "enso_east",
    3L,
    "_bhd-factor_adm-ADM1"
  )
)
touch(fits, fit_name("outbreak_duration_brms_cb", "re_stats", "enso_east", 1L)) # not a fit
touch(
  fits,
  fit_name("outbreak_monthly_brms", "fit", "enso_east", 7L, "_fam-negbinomial")
) # plain brms backend
touch(fits, fit_name("outbreak_size_brms_cb", "fit", "enso_east", 1L)) # retired: ignored
touch(fits, fit_name("outbreaks_survival_brms_cb", "fit", "enso_west", 4L))
touch(fits, fit_name("outbreak_duration_brms_cb", "fit", "enso_west", 6L))

# --- [1] per kind ---------------------------------------------------------------
# Per kind (occ / dur / mon): parses the sl<N> tag, matches only `_fit_` files
# (cb_meta etc. ignored), is backend-agnostic, returns NULL when absent,
# refuses `kind = "size"`.
cat("\n[1] read_deployed_lag() per kind\n")
r_occ <- read_deployed_lag(fits, "enso_east", GROUPING, "occ")
r_dur <- read_deployed_lag(fits, "enso_east", GROUPING, "dur")
r_mon <- read_deployed_lag(fits, "enso_east", GROUPING, "mon")
expect(
  identical(r_occ$lag, 5L) && grepl("_fit_", basename(r_occ$file)),
  "occ: the lp_season _fit_ file's sl5 tag (cb_meta and _ep-gs ignored)"
)
expect(identical(r_dur$lag, 3L), "dur: sl3 (re_stats ignored)")
expect(
  identical(r_mon$lag, 7L) &&
    grepl("^outbreak_monthly_brms_fit_", basename(r_mon$file)),
  "mon: sl7, a plain-brms (non-cb) fit matches too (backend-agnostic prefix)"
)
expect(
  is.null(read_deployed_lag(fits, "enso_west", GROUPING, "mon")),
  "NULL when the track has no deployed fit for the group"
)
expect(
  is.null(read_deployed_lag(fits, "enso_horn", GROUPING, "occ")),
  "NULL for a group with no fits at all"
)
expect(
  !is.null(errmsg(read_deployed_lag(fits, "enso_east", GROUPING, "size"))),
  "kind = \"size\" is refused (size is retired downstream)"
)
expect(
  identical(names(.DEPLOYED_FIT_PREFIX), c("occ", "dur", "mon")),
  ".DEPLOYED_FIT_PREFIX covers exactly occ / dur / mon"
)

# --- [2] episode key filter ----------------------------------------------------
# The "_ep-gs" filter both ways (mirrors simulations.R::.discover_group_fit()):
# an lp_season lookup ignores a gid_season fit and vice versa.
cat("\n[2] the _ep-gs filter, both ways\n")
r_gs <- read_deployed_lag(
  fits,
  "enso_east",
  GROUPING,
  "occ",
  episode_key = "gid_season"
)
expect(
  identical(r_gs$lag, 2L) && grepl("_ep-gs", basename(r_gs$file), fixed = TRUE),
  "gid_season lookup returns the _ep-gs fit (sl2)"
)
expect(
  !grepl("_ep-gs", basename(r_occ$file), fixed = TRUE),
  "lp_season lookup ignores the _ep-gs fit"
)
expect(
  is.null(read_deployed_lag(
    fits,
    "enso_east",
    GROUPING,
    "dur",
    episode_key = "gid_season"
  )),
  "gid_season lookup ignores lp_season fits (NULL)"
)

# --- [3] duplicates ---------------------------------------------------------------
# More than one matching fit is an error.
cat("\n[3] more than one matching fit is an error\n")
dups <- file.path(root, "dups")
touch(dups, fit_name("outbreak_duration_brms_cb", "fit", "enso_horn", 2L))
touch(dups, fit_name("outbreak_duration_brms_cb", "fit", "enso_horn", 4L))
e <- errmsg(read_deployed_lag(dups, "enso_horn", GROUPING, "dur"))
expect(
  !is.null(e) && grepl("multiple dur fit files", e, fixed = TRUE),
  sprintf("two dur fits for one group -> error (%s)", e %||% "no error")
)
touch(
  dups,
  fit_name("outbreak_monthly_brms_cb", "fit", "enso_horn", 2L, "", "_ep-gs")
)
touch(dups, fit_name("outbreak_monthly_brms_cb", "fit", "enso_horn", 3L))
expect(
  identical(read_deployed_lag(dups, "enso_horn", GROUPING, "mon")$lag, 3L),
  "one lp_season + one gid_season fit is not a duplicate"
)

# --- [4] forecast generator: horizon = last CFSv2 + min(occ, dur, mon) --------
# horizon_end = last CFSv2 month + min(best_lag_occ, best_lag_dur,
# best_lag_mon); a deployed size fit with a smaller lag is ignored; a group
# missing a track is skipped (naming it); the manifest carries the occ/dur/mon
# fields and no size field; no generated config pins a mean_rate table.
cat(
  "\n[4] generate_forecast_configs(): min(best_lag_occ, best_lag_dur, best_lag_mon)\n"
)
mk_tpl <- function(path, tag, start = "2026-04-01", handoff = NULL) {
  flt <- list(
    smoothing_input = file.path(root, "cache"),
    forecast_window = list(start = start, end = "2027-12-01")
  )
  if (!is.null(handoff)) {
    flt$handoff_month <- handoff
  }
  yaml::write_yaml(
    list(
      baseline = "1990_2019",
      simulator = list(
        backend = "spatpomp",
        run_tag = tag,
        climate_panel = file.path(root, "no_panel.rds"),
        location_frame = list(size_data_path = file.path(root, "no_size.rds")),
        horizon = list(start = start, end = "2027-12-01"),
        filter = flt,
        output = list(format = "qs")
      )
    ),
    path
  )
  path
}
fac_tpl <- mk_tpl(file.path(root, "fac.yml"), "fac_tag")
cf_tpl <- mk_tpl(file.path(root, "cf.yml"), "cf_tag")
meta <- file.path(root, "cfsv2_meta.yml")
yaml::write_yaml(list(target_range = list("2026-05-01", "2026-10-01")), meta)

fopt <- list(
  fits_dir = fits,
  grouping = GROUPING,
  cfsv2_meta = meta,
  factual_template = fac_tpl,
  counterfactual_template = cf_tpl,
  config_outdir = file.path(root, "fc_configs"),
  manifest = file.path(root, "forecast_horizons.yml"),
  outdir = file.path(root, "out/"),
  episode_key = "",
  partition_by = "",
  analog_sets = "all",
  primary_analog_set = "all",
  prune_stale = FALSE
)
msgs <- character(0)
fm <- withCallingHandlers(
  tryCatch(generate_forecast_configs(fopt), error = function(e) {
    cat("  generate_forecast_configs() error: ", conditionMessage(e), "\n")
    NULL
  }),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
)
expect(!is.null(fm), "generate_forecast_configs() runs on the fixture dir")
east <- fm$groups$enso_east
expect(
  identical(east$min_lag, 3L) && identical(east$horizon_end, "2027-01-01"),
  sprintf(
    "enso_east: min(5, 3, 7) = 3 -> horizon_end 2026-10-01 + 3 mo = 2027-01-01 (got %s / %s)",
    format(east$min_lag),
    format(east$horizon_end)
  )
)
expect(
  identical(
    c(east$best_lag_occ, east$best_lag_dur, east$best_lag_mon),
    c(5L, 3L, 7L)
  ),
  "manifest carries best_lag_occ / best_lag_dur / best_lag_mon"
)
expect(
  all(c("occ_fit_file", "dur_fit_file", "mon_fit_file") %in% names(east)) &&
    !any(grepl("size", names(east))),
  "manifest carries {occ,dur,mon}_fit_file and no size field"
)
expect(
  identical(sort(names(fm$groups)), "enso_east") &&
    setequal(
      unlist(fm$skipped_groups),
      c("enso_central", "enso_horn", "enso_south", "enso_west")
    ),
  "only enso_east is generated; the others are skipped"
)
expect(
  any(grepl(
    "SKIP group 'enso_west' (no deployed mon fit)",
    msgs,
    fixed = TRUE
  )),
  "the skip message names the missing track (enso_west: mon)"
)
fac_cfg <- yaml::read_yaml(east$factual_config)
cf_cfg <- yaml::read_yaml(east$counterfactual_config)
expect(
  identical(fac_cfg$simulator$horizon$end, "2027-01-01") &&
    identical(fac_cfg$simulator$filter$forecast_window$end, "2027-01-01"),
  "factual config horizon.end / forecast_window.end = the min-lag horizon"
)
expect(
  is.null(fac_cfg$simulator$location_frame$mean_rate_table_path) &&
    is.null(cf_cfg$simulator$location_frame$mean_rate_table_path) &&
    is.null(fac_cfg$simulator$location_frame$base_rate_table_paths),
  "no generated config pins a mean_rate table (and no anchor: no fit is anchored)"
)
expect(
  any(grepl(
    "best_lag_dur=3, best_lag_mon=7 -> min_lag=3",
    readLines(east$factual_config, n = 6L),
    fixed = TRUE
  )),
  "config header records the three lags and the min"
)

# An empty fits dir generates nothing and says why.
e <- errmsg(suppressMessages(generate_forecast_configs(
  modifyList(
    fopt,
    list(
      fits_dir = file.path(root, "empty"),
      config_outdir = file.path(root, "fc_empty"),
      manifest = file.path(root, "fc_empty.yml")
    )
  )
)))
expect(
  !is.null(e) && grepl("occurrence, duration AND monthly", e, fixed = TRUE),
  "no group with all three tracks -> error naming the three tracks"
)

# gid_season: the deploy's key selects the _ep-gs occ fit -- which has no
# dur/mon sibling here, so every group is skipped.
e <- errmsg(suppressMessages(generate_forecast_configs(
  modifyList(
    fopt,
    list(
      episode_key = "gid_season",
      config_outdir = file.path(root, "fc_gs"),
      manifest = file.path(root, "fc_gs.yml")
    )
  )
)))
expect(
  !is.null(e) && grepl("episodes.key = gid_season", e, fixed = TRUE),
  "episode_key = gid_season looks fits up under _ep-gs (none complete here)"
)

# --- [4b] per-arm launch (handoff) months ---------------------------------------
cat("\n[4b] per-arm launches: manifest keys, attribution start, refusals\n")
expect(
  identical(east$factual_launch, "2026-04-01") &&
    identical(east$counterfactual_launch, "2026-04-01") &&
    identical(east$counterfactual_sets$all$launch, "2026-04-01") &&
    identical(fm$attribution_start, "2026-04-01"),
  "equal launches: factual_launch / counterfactual_launch / set launch / attribution_start recorded"
)
expect(is.null(fm$schema), "v1 manifest still carries no schema token")
al <- manifest_arm_launches(fm)
expect(
  nrow(al) == 2L &&
    all(al$source == "manifest") &&
    all(al$launch == as.Date("2026-04-01")),
  "manifest_arm_launches() reads both arms from the manifest keys"
)
expect(
  any(grepl(
    "Launch (smoothing handoff) 2026-04-01",
    readLines(east$factual_config, n = 8L),
    fixed = TRUE
  )),
  "config header states the launch (stale 2025-01-01 handoff text gone)"
)
# Counterfactual launching earlier (at the analog anchor): stitch window.
cf_early <- mk_tpl(
  file.path(root, "cf_early.yml"),
  "cf_tag",
  start = "2026-02-01",
  handoff = "2026-02-01"
)
fm2 <- suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    counterfactual_template = cf_early,
    config_outdir = file.path(root, "fc_early"),
    manifest = file.path(root, "fc_early.yml")
  )
)))
e2 <- fm2$groups$enso_east
expect(
  identical(e2$factual_launch, "2026-04-01") &&
    identical(e2$counterfactual_launch, "2026-02-01") &&
    identical(fm2$attribution_start, "2026-02-01"),
  "earlier counterfactual: per-arm launches kept, attribution_start = the earlier one"
)
expect(
  manifest_attrib_start(fm2) == as.Date("2026-02-01") &&
    manifest_attrib_start(fm2, alignment = "common_start") ==
      as.Date("2026-04-01") &&
    manifest_factual_start(fm2) == as.Date("2026-04-01"),
  "manifest_attrib_start(): stitch = earlier launch, common_start = later launch"
)
expect(
  identical(
    attribution_ratio_label("2026-04-01", "2026-02-01"),
    "ENSO + post-handoff information"
  ) &&
    identical(attribution_ratio_label("2026-04-01", "2026-04-01"), "ENSO"),
  "attribution_ratio_label() names the information gap only when the launches differ"
)
# Older manifests (no launch keys) fall back to each config's horizon.start.
old <- fm2
for (k in c("factual_launch", "counterfactual_launch")) {
  old$groups$enso_east[[k]] <- NULL
}
old$groups$enso_east$counterfactual_sets$all$launch <- NULL
al_old <- manifest_arm_launches(old)
expect(
  all(al_old$source == "config") &&
    identical(sort(format(al_old$launch)), c("2026-02-01", "2026-04-01")),
  "manifest without launch keys: launches read from the configs"
)
# Refusals.
cf_late <- mk_tpl(
  file.path(root, "cf_late.yml"),
  "cf_tag",
  start = "2026-06-01"
)
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    counterfactual_template = cf_late,
    config_outdir = file.path(root, "fc_late"),
    manifest = file.path(root, "fc_late.yml")
  )
))))
expect(
  !is.null(e) && grepl("after the factual template", e, fixed = TRUE),
  "a counterfactual launching after the factual arm is refused"
)
# Closed horizons: firm == launch gives one firm month
# (the launch month is simulated) and is accepted; only firm < launch is refused.
fac_eq <- mk_tpl(file.path(root, "fac_eq.yml"), "fac_tag", start = "2027-01-01")
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    factual_template = fac_eq,
    config_outdir = file.path(root, "fc_firm_eq"),
    manifest = file.path(root, "fc_firm_eq.yml")
  )
))))
expect(
  is.null(e),
  "firm horizon == factual launch (2027-01-01) is accepted: one firm month"
)
fac_late <- mk_tpl(
  file.path(root, "fac_late.yml"),
  "fac_tag",
  start = "2027-02-01"
)
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    factual_template = fac_late,
    config_outdir = file.path(root, "fc_firm"),
    manifest = file.path(root, "fc_firm.yml")
  )
))))
expect(
  !is.null(e) &&
    grepl("factual arm's launch 2027-02-01", e, fixed = TRUE) &&
    grepl("is before", e, fixed = TRUE),
  "firm horizon before the later (factual) launch is refused and names the arm"
)

# --- [4c] firm horizon tied to the ENSO seam -----------------------------------
cat("\n[4c] firm horizon = enso_seam_month(last CFSv2, seam) + min lag\n")
expect(
  is.null(fm$firm_oni_end) && is.null(fm$firm_horizon_rule),
  "legacy seam: no firm_oni_end / firm_horizon_rule in the manifest (an older manifest regenerates identically)"
)
expect(
  identical(manifest_last_firm_oni(fm), "2026-10-01"),
  "legacy: manifest_last_firm_oni() = last CFSv2 month"
)
fm_s <- suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    enso_continuation = list(method = "analog", seam = "last_full_oni"),
    config_outdir = file.path(root, "fc_seam"),
    manifest = file.path(root, "fc_seam.yml")
  )
)))
expect(
  identical(fm_s$groups$enso_east$horizon_end, "2026-12-01"),
  sprintf(
    "seam last_full_oni: firm horizon = 2026-09-01 + 3 mo = 2026-12-01 (got %s)",
    format(fm_s$groups$enso_east$horizon_end)
  )
)
expect(
  identical(fm_s$firm_oni_end, "2026-09-01") &&
    identical(fm_s$firm_horizon_rule, "last_full_oni + min lag"),
  "seam last_full_oni: firm_oni_end / firm_horizon_rule recorded"
)
expect(
  identical(manifest_last_firm_oni(fm_s), "2026-09-01"),
  "manifest_last_firm_oni() reads firm_oni_end"
)
fs_cfg <- yaml::read_yaml(fm_s$groups$enso_east$factual_config)
expect(
  identical(fs_cfg$simulator$horizon$end, "2026-12-01") &&
    identical(fs_cfg$simulator$filter$forecast_window$end, "2026-12-01"),
  "seam last_full_oni: horizon.end == forecast_window.end, one month earlier"
)
expect(
  any(grepl(
    "firm ONI end, last CFSv2 - 1 2026-09-01 + 3 mo",
    readLines(fm_s$groups$enso_east$factual_config, n = 8L),
    fixed = TRUE
  )),
  "seam last_full_oni: the config header names the firm ONI end"
)
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    enso_continuation = list(method = "analog", seam = "bogus"),
    config_outdir = file.path(root, "fc_seam_bad"),
    manifest = file.path(root, "fc_seam_bad.yml")
  )
))))
expect(!is.null(e), "an unknown seam is refused")
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    enso_continuation = list(method = "neutral", seam = "last_full_oni"),
    config_outdir = file.path(root, "fc_seam_neu"),
    manifest = file.path(root, "fc_seam_neu.yml")
  )
))))
expect(
  !is.null(e) &&
    grepl("needs enso_continuation.method 'analog'", e, fixed = TRUE),
  "seam last_full_oni without the analog continuation is refused"
)
expect(
  identical(fs_cfg$simulator$enso$continuation$seam, "last_full_oni"),
  "seam last_full_oni is written into the forward config's continuation block"
)
fm_d <- suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    enso_continuation = list(method = "analog", seam = "last_cfsv2"),
    config_outdir = file.path(root, "fc_seam_def"),
    manifest = file.path(root, "fc_seam_def.yml")
  )
)))
fd_cfg <- yaml::read_yaml(fm_d$groups$enso_east$factual_config)
expect(
  is.null(fd_cfg$simulator$enso$continuation$seam) &&
    is.null(fm_d$firm_oni_end),
  "a stated default seam (last_cfsv2) is not written: configs and manifest as without it"
)

# --- [4d] counterfactual analog cache, library mode -----------------------------
cat(
  "\n[4d] enso_counterfactual: legacy builds nothing; library builds the per-run cache in-process\n"
)
expect(
  is.null(fm$enso_counterfactual),
  "legacy (no block): no enso_counterfactual in the manifest"
)
fm_l <- suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    enso_counterfactual = list(match = "legacy"),
    config_outdir = file.path(root, "fc_cf_leg"),
    manifest = file.path(root, "fc_cf_leg.yml")
  )
)))
expect(
  is.null(fm_l$enso_counterfactual),
  "match: legacy stated: manifest unchanged"
)
# Library fixture: synthetic piControl library + ONI cache in the generator's outdir.
set.seed(81)
lib_series <- lapply(c(A = 1, B = 1.4), function(k) {
  as.numeric(stats::filter(
    rnorm(12L * 200L, 0, 0.25),
    0.93,
    method = "recursive"
  )) *
    k
})
lib_feats <- do.call(
  rbind,
  lapply(names(lib_series), function(m) {
    enso_continuation_features(
      lib_series[[m]],
      rep(1:12, 200L),
      model = m,
      horizon = 48L
    )
  })
)
lib_q <- list(
  series = lib_series,
  features = lib_feats,
  feature_sd = vapply(
    c("v", "tend", "phase"),
    function(f) sd(lib_feats[[f]][lib_feats$ok]),
    numeric(1)
  ),
  meta = list(
    horizon_months = 48L,
    rescale = list(window = "1991_2020", obs_sd = 0.8)
  )
)
lib_q$meta$digest <- enso_continuation_library_digest(lib_q)
qs2::qs_save(lib_q, file.path(root, "lib.qs"))
od <- seq(as.Date("2024-01-01"), as.Date("2026-08-01"), by = "1 month")
dir.create(file.path(root, "out"), showWarnings = FALSE)
saveRDS(
  data.frame(
    year = as.integer(format(od, "%Y")),
    month = as.integer(format(od, "%m")),
    enso_index = 0.8 * sin(seq_along(od) / 5),
    date = od
  ),
  file.path(root, "out", "enso_oni_cache.rds")
)
cf_tpl_l <- yaml::read_yaml(cf_tpl)
cf_tpl_l$simulator$enso <- list(
  counterfactual = list(
    scenario = "lenssen",
    cache_path = file.path(
      root,
      "lenssen_prod_run",
      "lenssen_nino34_analogs.rds"
    )
  )
)
yaml::write_yaml(cf_tpl_l, file.path(root, "cf_lib.yml"))
fm_lib <- suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    counterfactual_template = file.path(root, "cf_lib.yml"),
    enso_counterfactual = list(
      match = "library",
      n_analogs = 20L,
      winters = "all_read",
      min_set_size = 1L
    ),
    enso_continuation = list(
      method = "analog",
      library = file.path(root, "lib.qs")
    ),
    config_outdir = file.path(root, "fc_cf_lib"),
    manifest = file.path(root, "fc_cf_lib.yml")
  )
)))
ecf <- fm_lib$enso_counterfactual
expect(
  !is.null(ecf) &&
    file.exists(ecf$cache_path) &&
    identical(ecf$match, "library") &&
    nzchar(ecf$digest),
  "library: per-run cache built in-process, its digest in the manifest"
)
expect(
  identical(ecf$read_end, "2026-10-01"),
  sprintf(
    "read end = enso_east horizon 2027-01 - min lag 3 = 2026-10, the last CFSv2 month (got %s)",
    ecf$read_end %||% "NULL"
  )
)
expect(
  identical(as.integer(unlist(ecf$winters$primary)), 2026L) &&
    !length(unlist(ecf$winters$other)),
  "no DJF read in the window: sets on the first DJF after the launch, nothing else constrained"
)
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    enso_counterfactual = list(match = "library"),
    enso_continuation = list(
      method = "analog",
      library = file.path(root, "lib.qs")
    ),
    config_outdir = file.path(root, "fc_cf_shared"),
    manifest = file.path(root, "fc_cf_shared.yml")
  )
))))
expect(
  !is.null(e) && grepl("shared legacy cache", e),
  "library mode with the shared default cache path is refused"
)
bad_ho <- mk_tpl(
  file.path(root, "fac_badho.yml"),
  "fac_tag",
  handoff = "2026-03-01"
)
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  fopt,
  list(
    factual_template = bad_ho,
    config_outdir = file.path(root, "fc_badho"),
    manifest = file.path(root, "fc_badho.yml")
  )
))))
expect(
  !is.null(e) && grepl("handoff_month", e, fixed = TRUE),
  "a template whose handoff_month differs from horizon.start is refused at generation"
)

# --- [5] hindcast generator ------------------------------------------------------
# Same min-lag rule on top of the season end (--add_best_lag), same manifest
# fields.
cat("\n[5] generate_hindcast_configs(): season end + min(occ, dur, mon)\n")
hopt <- list(
  template = fac_tpl,
  fits_dir = fits,
  grouping = GROUPING,
  climate_panel = file.path(root, "no_hist_panel.rds"),
  launch_month = "05-01",
  launch_years = "2015",
  season_end_by_group = "",
  horizon_months = 12L,
  add_best_lag = TRUE,
  no_add_best_lag = FALSE,
  smoothing_input = "",
  episode_key = "",
  partition_by = "",
  run_tag_base = "hc_tag",
  config_outdir = file.path(root, "hc_configs"),
  manifest = file.path(root, "hindcast_horizons.yml"),
  outdir = file.path(root, "out/")
)
hm <- tryCatch(
  suppressMessages(generate_hindcast_configs(hopt)),
  error = function(e) {
    cat("  generate_hindcast_configs() error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(!is.null(hm), "generate_hindcast_configs() runs on the fixture dir")
he <- hm$groups$enso_east
expect(
  identical(he$min_lag, 3L) &&
    identical(he$years[["2015"]]$horizon_end, "2016-08-01"),
  sprintf(
    "enso_east 2015: launch 2015-05-01 + 12 mo + min_lag 3 = 2016-08-01 (got %s)",
    format(he$years[["2015"]]$horizon_end)
  )
)
expect(
  identical(
    c(he$best_lag_occ, he$best_lag_dur, he$best_lag_mon),
    c(5L, 3L, 7L)
  ) &&
    !any(grepl("size", names(he))),
  "hindcast manifest carries the three lags and no size field"
)
hc_cfg <- yaml::read_yaml(he$years[["2015"]]$config)
expect(
  is.null(hc_cfg$simulator$location_frame$mean_rate_table_path),
  "hindcast config pins no mean_rate table"
)
hm2 <- suppressMessages(generate_hindcast_configs(
  modifyList(
    hopt,
    list(
      add_best_lag = FALSE,
      no_add_best_lag = TRUE,
      config_outdir = file.path(root, "hc_nolag"),
      manifest = file.path(root, "hc_nolag.yml")
    )
  )
))
expect(
  identical(hm2$groups$enso_east$years[["2015"]]$horizon_end, "2016-05-01"),
  "--no_add_best_lag: horizon_end = launch + 12 mo"
)

# --- [6] partition_by: registry -- units mixing fit groups ------------------------
cat(
  "\n[6] registry units: max member firm horizon, per-country firm horizons\n"
)
# enso_horn gets all three tracks (min lag 8); a temporary partition registry
# merges enso_east + enso_horn into one unit 'eh' (other groups have no fits).
touch(fits, fit_name("outbreaks_survival_brms_cb", "fit", "enso_horn", 8L))
touch(fits, fit_name("outbreak_duration_brms_cb", "fit", "enso_horn", 9L))
touch(fits, fit_name("outbreak_monthly_brms_cb", "fit", "enso_horn", 10L))
preg <- load_grouping("enso_grouping")
preg$group <- ifelse(
  preg$group %in% c("enso_east", "enso_horn"),
  "eh",
  sub("^enso_", "", preg$group)
)
PSRC <- sprintf("zz_test_partition_%d", Sys.getpid())
ppath <- file.path("analysis/groupings", paste0(PSRC, ".csv"))
utils::write.csv(preg, ppath, row.names = FALSE, quote = FALSE)
ropt <- modifyList(
  fopt,
  list(
    partition_by = "registry",
    partition_source = PSRC,
    config_outdir = file.path(root, "fc_reg"),
    manifest = file.path(root, "fc_reg.yml")
  )
)
rm_ <- tryCatch(
  suppressMessages(generate_forecast_configs(ropt)),
  error = function(e) {
    cat(
      "  registry generate_forecast_configs() error: ",
      conditionMessage(e),
      "\n"
    )
    NULL
  }
)
expect(!is.null(rm_), "registry-mode forecast generation runs")
if (!is.null(rm_)) {
  u <- rm_$units$eh
  expect(
    identical(names(rm_$units), "eh"),
    "one unit (the only one with deployed fits)"
  )
  expect(
    identical(unlist(u$groups), c("enso_east", "enso_horn")),
    "unit 'eh' mixes enso_east and enso_horn"
  )
  expect(
    identical(u$horizon_end, "2027-06-01"),
    sprintf(
      "unit horizon = max(firm east 2027-01-01, horn 2027-06-01) (got %s)",
      u$horizon_end
    )
  )
  expect(
    identical(rm_$countries$KEN$firm_horizon_end, "2027-01-01") &&
      identical(rm_$countries$KEN$unit_horizon_end, "2027-06-01") &&
      identical(rm_$countries$ETH$firm_horizon_end, "2027-06-01"),
    "countries carry their own firm horizon and the unit horizon"
  )
  expect(
    identical(rm_$schema, 2L) &&
      identical(rm_$partition_source, PSRC) &&
      identical(rm_$partition_digest, registry_digest(PSRC)) &&
      identical(rm_$grouping_digest, registry_digest(GROUPING)),
    "v2 manifest identity keys"
  )
  expect(
    identical(rm_$africa_min_horizon, "2027-01-01") &&
      identical(rm_$africa_max_horizon, "2027-06-01"),
    "africa_min_horizon = min firm, africa_max_horizon = max unit"
  )
  expect(
    identical(rm_$groups$enso_east$horizon_end, "2027-01-01") &&
      identical(unlist(rm_$groups$enso_east$units), "eh") &&
      is.null(rm_$groups$enso_east$factual_config),
    "groups keep their firm horizon and list their units (no run fields)"
  )
  rc <- yaml::read_yaml(u$factual_config)
  expect(
    basename(u$factual_config) == "forecast_eh.yml" &&
      identical(rc$simulator$run_tag, "fac_tag_eh"),
    "config + run tag keyed by unit"
  )
  expect(
    identical(rc$simulator$filter$partition_by, "registry") &&
      identical(rc$simulator$filter$partition_source, PSRC) &&
      identical(rc$simulator$nbr_contagion$partition_source, PSRC),
    "config carries partition_by/source and the contagion registry (default = partition)"
  )
  expect(
    setequal(
      strsplit(rc$simulator$scope$countries, ",")[[1]],
      preg$iso3[preg$group == "eh"]
    ),
    "scope.countries = unit members"
  )
  expect(
    identical(rc$simulator$horizon$end, "2027-06-01"),
    "config horizon.end = unit horizon"
  )
}
# a partition registry missing a country with deployed fits is refused
preg2 <- preg[preg$iso3 != "KEN", ]
utils::write.csv(preg2, ppath, row.names = FALSE, quote = FALSE)
e <- errmsg(suppressMessages(generate_forecast_configs(modifyList(
  ropt,
  list(
    config_outdir = file.path(root, "fc_reg2"),
    manifest = file.path(root, "fc_reg2.yml")
  )
))))
expect(
  !is.null(e) && grepl("KEN", e),
  "partition registry not covering KEN refused"
)
utils::write.csv(preg, ppath, row.names = FALSE, quote = FALSE)
rh <- tryCatch(
  suppressMessages(generate_hindcast_configs(modifyList(
    hopt,
    list(
      partition_by = "registry",
      partition_source = PSRC,
      season_end_by_group = "enso_east=7",
      config_outdir = file.path(root, "hc_reg"),
      manifest = file.path(root, "hc_reg.yml")
    )
  ))),
  error = function(e) {
    cat("  registry hindcast error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(!is.null(rh), "registry-mode hindcast generation runs")
if (!is.null(rh)) {
  y <- rh$units$eh$years[["2015"]]
  # east: season end 2015-07 + 3 = 2015-10-01; horn: 2015-05 + 12 + 8 = 2017-01-01
  expect(
    identical(y$horizon_end, "2017-01-01") &&
      identical(rh$groups$enso_east$years[["2015"]]$horizon_end, "2015-10-01"),
    sprintf(
      "hindcast unit-year horizon = max of member firm windows (got %s)",
      y$horizon_end
    )
  )
  expect(
    basename(y$config) == "hindcast_eh_2015.yml" && identical(rh$schema, 2L),
    "hindcast configs keyed (unit, year); v2 manifest"
  )
}
unlink(ppath)


cat("\n[7] read_deployed_lag(): lead fits (slm<k>, 2026-10-17)\n")
dl <- file.path(tempdir(), "deployed_leads")
dir.create(dl, showWarnings = FALSE)
for (n in c(
  "outbreaks_survival_brms_cb_fit_x_grp-gg-east_cb-enso_slm6_ns3_ep-gs_nz.rds",
  "outbreaks_survival_brms_cb_fit_x_grp-gg-west_cb-enso_sl11_ns3_ep-gs_nz.rds"
)) {
  file.create(file.path(dl, n))
}
expect(
  identical(
    read_deployed_lag(dl, "east", "gg", "occ", "gid_season", TRUE)$lag,
    -6L
  ),
  "slm6 parses to lag -6"
)
expect(
  identical(
    read_deployed_lag(dl, "west", "gg", "occ", "gid_season", TRUE)$lag,
    11L
  ),
  "sl11 still parses to 11"
)

# With deployed ENSO leads (negative min lag) the window ends at the season
# end; the lead does not shorten it (observed ENSO).
cat(
  "\n[8] hindcast: a deployed lead does not shorten the window (2026-10-18)\n"
)
fl <- file.path(root, "fits_lead")
for (n in c(
  "outbreaks_survival_brms_cb_fit_base-1990_2019_grp-enso_grouping-enso_west_slm6.rds",
  "outbreak_duration_brms_cb_fit_base-1990_2019_grp-enso_grouping-enso_west_slm4.rds",
  "outbreak_monthly_brms_cb_fit_base-1990_2019_grp-enso_grouping-enso_west_slm3.rds",
  "outbreaks_survival_brms_cb_fit_base-1990_2019_grp-enso_grouping-enso_east_sl0.rds",
  "outbreak_duration_brms_cb_fit_base-1990_2019_grp-enso_grouping-enso_east_slm1.rds",
  "outbreak_monthly_brms_cb_fit_base-1990_2019_grp-enso_grouping-enso_east_sl3.rds"
)) {
  touch(fl, n)
}
hl <- tryCatch(
  suppressMessages(generate_hindcast_configs(modifyList(
    hopt,
    list(
      fits_dir = fl,
      season_end_by_group = "enso_west=11,enso_east=7",
      config_outdir = file.path(root, "hc_lead"),
      manifest = file.path(root, "hc_lead.yml")
    )
  ))),
  error = function(e) {
    cat("  lead hindcast error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(
  !is.null(hl),
  "lead fits (west min lag -6, season May-Nov) no longer stop the generator"
)
if (!is.null(hl)) {
  expect(
    identical(hl$groups$enso_west$min_lag, -6L) &&
      identical(hl$groups$enso_west$years[["2015"]]$horizon_end, "2015-11-01"),
    sprintf(
      "enso_west 2015: window ends at the season end 2015-11-01 (got %s)",
      hl$groups$enso_west$years[["2015"]]$horizon_end
    )
  )
  expect(
    identical(hl$groups$enso_east$years[["2015"]]$horizon_end, "2015-07-01"),
    sprintf(
      "enso_east 2015 (min lag -1): season end 2015-07-01, not 2015-06-01 (got %s)",
      hl$groups$enso_east$years[["2015"]]$horizon_end
    )
  )
}

# Hindcast and event generators on a synthetic panel calendar: per-country
# firm horizons (hindcast: first own end >= launch + 6 months + max(min lag,
# 0); event: first own end on/after the peak + max lag), group = max over
# members, unit = its latest member, countries.<iso> fields,
# season_end_source, the CSV; by_group unchanged plus the > 4-month-span
# warning.
cat(
  "\n[9] --season_end by_country_from_panel (per-country season ends, 2026-10-27)\n"
)
# Synthetic panel calendar: KEN min_seas 8 (end Jul), TZA 12 (end Nov), ETH 5
# (end Apr), every other east / horn country 1 (end Dec). Fits from [1]-[6]:
# enso_east lags 5/3/7 (min 3, max 7), enso_horn 8/9/10 (min 8, max 10);
# enso_west has no mon fit (skipped).
source("analysis/utils/make_event_attribution_configs.R")
ms_of <- c(KEN = 8L, TZA = 12L, ETH = 5L)
reg9 <- load_grouping(GROUPING)
iso9 <- reg9$iso3[reg9$group %in% c("enso_east", "enso_horn")]
pan9 <- tidyr::expand_grid(
  iso = iso9,
  k = 1:2,
  t_start = seq(as.Date("2014-01-01"), as.Date("2018-12-01"), by = "month")
)
pan9$gid <- paste0(pan9$iso, "_", pan9$k)
pan9$country <- pan9$iso
pan9$min_seas <- ifelse(pan9$iso %in% names(ms_of), ms_of[pan9$iso], 1L)
pan9 <- as.data.frame(pan9[, c("gid", "country", "min_seas", "t_start")])
ppan <- file.path(root, "hist_panel_calendar.rds")
saveRDS(pan9, ppan)

hc9 <- tryCatch(
  suppressWarnings(suppressMessages(generate_hindcast_configs(modifyList(
    hopt,
    list(
      climate_panel = ppan,
      season_end = "by_country_from_panel",
      season_end_by_group = "enso_east=7",
      config_outdir = file.path(root, "hc_bc"),
      manifest = file.path(root, "hc_bc.yml")
    )
  )))),
  error = function(e) {
    cat("  by_country hindcast error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(!is.null(hc9), "by_country_from_panel hindcast generation runs")
if (!is.null(hc9)) {
  c9 <- hc9$countries
  expect(
    identical(c9$KEN$years[["2015"]]$firm_horizon_end, "2016-10-01") &&
      identical(c9$TZA$years[["2015"]]$firm_horizon_end, "2016-02-01") &&
      identical(c9$BDI$years[["2015"]]$firm_horizon_end, "2016-03-01"),
    "east 2015: KEN Jul>=Nov15 -> 2016-07+3, TZA 2015-11+3, Dec-end 2015-12+3"
  )
  expect(
    identical(c9$ETH$years[["2015"]]$firm_horizon_end, "2016-12-01") &&
      identical(c9$SDN$years[["2015"]]$firm_horizon_end, "2016-08-01"),
    "horn 2015: ETH 2016-04+8, Dec-end 2015-12+8"
  )
  expect(
    identical(hc9$groups$enso_east$years[["2015"]]$horizon_end, "2016-10-01") &&
      identical(
        hc9$groups$enso_horn$years[["2015"]]$horizon_end,
        "2016-12-01"
      ) &&
      is.null(hc9$groups$enso_east$season_end_month),
    "group horizon = max over its members; no group season end"
  )
  expect(
    identical(c9$KEN$season_end_month, 7L) &&
      identical(c9$KEN$min_seas_mode, 8L) &&
      identical(c9$KEN$modal_share, 1) &&
      identical(c9$KEN$unit, "enso_east") &&
      identical(c9$TZA$years[["2015"]]$unit_horizon_end, "2016-10-01"),
    "countries.<iso>: season end, min_seas mode, modal share, unit and unit horizon"
  )
  expect(
    identical(hc9$season_end, "by_country_from_panel") &&
      identical(hc9$season_min_window_months, 6L) &&
      identical(hc9$season_end_source$panel, ppan) &&
      nzchar(hc9$season_end_source$min_seas_digest),
    "manifest records the season-end mode, window rule and source (panel + digest)"
  )
  csv9 <- file.path(root, "season_ends_hc_bc.csv")
  expect(
    file.exists(csv9) &&
      all(
        c(
          "iso3",
          "season_end_month",
          "firm_horizon_end",
          "unit_horizon_end",
          "group_season_end_month"
        ) %in%
          names(read.csv(csv9))
      ),
    "season_ends_<manifest stem>.csv next to the manifest"
  )
  hz9 <- manifest_country_horizons(hc9, "firm", year = 2015)
  expect(
    identical(hz9[["TZA"]], "2016-02-01") &&
      identical(hz9[["KEN"]], "2016-10-01"),
    "sim_io reads the per-country firm horizons back"
  )
  cfg9 <- yaml::read_yaml(hc9$groups$enso_east$years[["2015"]]$config)
  expect(
    identical(cfg9$simulator$horizon$end, "2016-10-01"),
    "the run goes to the latest member"
  )
}
# by_group on the same panel: output as without the panel, plus a diagnostic
hg9 <- withCallingHandlers(
  suppressMessages(generate_hindcast_configs(modifyList(
    hopt,
    list(
      climate_panel = ppan,
      config_outdir = file.path(root, "hc_bg"),
      manifest = file.path(root, "hc_bg.yml")
    )
  ))),
  warning = function(w) {
    assign("w9", conditionMessage(w), envir = globalenv())
    invokeRestart("muffleWarning")
  }
)
expect(
  identical(hg9$groups$enso_east$years[["2015"]]$horizon_end, "2016-08-01") &&
    is.null(hg9$countries) &&
    is.null(hg9$season_end) &&
    !file.exists(file.path(root, "season_ends_hc_bg.csv")),
  "by_group (default): windows unchanged, no countries / season_end keys, no CSV"
)
expect(
  exists("w9") && grepl("enso_east", w9) && grepl("span", w9),
  "by_group warns when a group's member season ends span > 4 months (Jul / Nov / Dec)"
)
expect(
  inherits(
    try(
      suppressMessages(generate_hindcast_configs(modifyList(
        hopt,
        list(
          season_end = "by_country",
          config_outdir = file.path(root, "hc_x"),
          manifest = file.path(root, "hc_x.yml")
        )
      ))),
      silent = TRUE
    ),
    "try-error"
  ),
  "an unknown --season_end is refused"
)

# event generator, registry units: east + horn in one unit 'eh'
utils::write.csv(preg, ppath, row.names = FALSE, quote = FALSE)
cf_ev <- file.path(root, "cf_ev.yml")
tt <- yaml::read_yaml(fac_tpl)
tt$simulator$enso$counterfactual <- list(scenario = "lenssen")
yaml::write_yaml(tt, cf_ev)
eopt9 <- list(
  template = fac_tpl,
  counterfactual_template = cf_ev,
  events = "elnino_2015",
  events_registry = "analysis/configs/enso_events_v2.yml",
  analog_set = "neutral",
  lenssen_dir = file.path(root, "lenssen"),
  fits_dir = fits,
  grouping = GROUPING,
  climate_panel = ppan,
  season_end = "by_country_from_panel",
  season_end_by_group = "",
  smoothing_input = "",
  episode_key = "",
  partition_by = "registry",
  partition_source = PSRC,
  nbr_partition_source = NULL,
  run_tag_base = "ev_tag",
  config_outdir = file.path(root, "ev_bc"),
  manifest = file.path(root, "ev_bc.yml"),
  outdir = file.path(root, "out/"),
  param_draw_mode = ""
)
ev9 <- tryCatch(
  suppressMessages(generate_event_attribution_configs(eopt9))$manifest,
  error = function(e) {
    cat("  by_country event error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(!is.null(ev9), "by_country_from_panel event generation runs")
if (!is.null(ev9)) {
  e9 <- ev9$events$elnino_2015
  expect(
    identical(e9$countries$KEN$firm_horizon_end, "2017-02-01") &&
      identical(e9$countries$TZA$firm_horizon_end, "2017-06-01") &&
      identical(e9$countries$BDI$firm_horizon_end, "2016-07-01") &&
      identical(e9$countries$ETH$firm_horizon_end, "2017-02-01") &&
      identical(e9$countries$SDN$firm_horizon_end, "2016-10-01"),
    "peak Dec 2015: first own season end on/after the peak + the group's max lag (7 east, 10 horn)"
  )
  expect(
    identical(e9$groups$enso_east$horizon_end, "2017-06-01") &&
      identical(e9$groups$enso_horn$horizon_end, "2017-02-01") &&
      is.null(e9$groups$enso_east$season_end_month),
    "group horizon = its latest member"
  )
  expect(
    identical(e9$units$eh$horizon_end, "2017-06-01") &&
      all(vapply(e9$countries, `[[`, "", "unit_horizon_end") == "2017-06-01"),
    "unit runs to its latest member country"
  )
  expect(
    identical(e9$countries$KEN$season_end_month, 7L) &&
      identical(e9$countries$KEN$season_end, "2016-07-01"),
    "countries.<iso> carries the season end"
  )
  expect(
    identical(ev9$season_end, "by_country_from_panel") &&
      identical(ev9$season_end_source$panel, ppan) &&
      file.exists(file.path(root, "season_ends_ev_bc.csv")),
    "event manifest records season_end + source; CSV next to it"
  )
  v9 <- event_manifest_view(ev9, "elnino_2015")
  expect(
    identical(manifest_country_horizons(v9, "firm")[["BDI"]], "2016-07-01"),
    "the 09 view censors each country at its own horizon"
  )
}
e_bg <- tryCatch(
  suppressMessages(generate_event_attribution_configs(modifyList(
    eopt9,
    list(
      season_end = "by_group",
      config_outdir = file.path(root, "ev_bg"),
      manifest = file.path(root, "ev_bg.yml")
    )
  ))),
  error = function(e) conditionMessage(e)
)
expect(
  is.character(e_bg) && grepl("lacks fit group", e_bg),
  "by_group still requires a season end for every fit group"
)
# --- [10] multi-winter event with per-winter sub-windows ------------------------
# A multi-winter event (enso_events_v2.yml shape): consecutive per-winter
# sub-windows per group and per country, schema-3 manifest, library cache
# labelled on every primary winter; legacy matcher refused; single-window
# events still give schema 2, unchanged.
cat(
  "\n[10] one multi-winter event: per-winter sub-windows, schema-3 manifest, library cache on every winter\n"
)
# Registry: a two-winter event on the [9] panel (2014-2018), launch 2014-09,
# peaks Dec 2015 / Dec 2016. Library fixture with a long horizon and an ONI
# cache over 2010-2019 in its own outdir.
reg10 <- file.path(root, "events_v2.yml")
yaml::write_yaml(
  list(
    events = list(
      elnino_2015 = list(
        phase = "elnino",
        label = "El Nino 2015-16",
        peak = "2015-12-01",
        peak_season_year = 2015L,
        launch = "2014-09-01",
        launch_rule = "crossing"
      ),
      multi_2015_17 = list(
        phase = "elnino",
        label = "Two winters",
        peaks = list("2015-12-01", "2016-12-01"),
        peak_season_years = list(2015L, 2016L),
        launch = "2014-09-01",
        launch_rule = "crossing",
        sub_windows = "per_winter"
      )
    )
  ),
  reg10
)
set.seed(810)
lib10_series <- lapply(c(A = 1, B = 1.3), function(k) {
  as.numeric(stats::filter(
    rnorm(12L * 300L, 0, 0.25),
    0.93,
    method = "recursive"
  )) *
    k
})
lib10_feats <- do.call(
  rbind,
  lapply(names(lib10_series), function(m) {
    enso_continuation_features(
      lib10_series[[m]],
      rep(1:12, 300L),
      model = m,
      horizon = 72L
    )
  })
)
lib10 <- list(
  series = lib10_series,
  features = lib10_feats,
  feature_sd = vapply(
    c("v", "tend", "phase"),
    function(f) sd(lib10_feats[[f]][lib10_feats$ok]),
    numeric(1)
  ),
  meta = list(
    horizon_months = 72L,
    rescale = list(window = "1991_2020", obs_sd = 0.8)
  )
)
lib10$meta$digest <- enso_continuation_library_digest(lib10)
qs2::qs_save(lib10, file.path(root, "lib10.qs"))
od10 <- seq(as.Date("2010-01-01"), as.Date("2019-12-01"), by = "1 month")
dir.create(file.path(root, "out10"), showWarnings = FALSE)
saveRDS(
  data.frame(
    year = as.integer(format(od10, "%Y")),
    month = as.integer(format(od10, "%m")),
    enso_index = 0.8 * sin(seq_along(od10) / 5),
    date = od10
  ),
  file.path(root, "out10", "enso_oni_cache.rds")
)
utils::write.csv(preg, ppath, row.names = FALSE, quote = FALSE)
eopt10 <- modifyList(
  eopt9,
  list(
    events = "elnino_2015,multi_2015_17",
    events_registry = reg10,
    lenssen_dir = file.path(root, "lenssen10"),
    outdir = file.path(root, "out10/"),
    config_outdir = file.path(root, "ev10"),
    manifest = file.path(root, "ev10.yml"),
    enso_counterfactual = list(
      match = "library",
      n_analogs = 40L,
      winters = "all_read",
      min_set_size = 1L
    ),
    enso_continuation = list(library = file.path(root, "lib10.qs"))
  )
)
ev10 <- tryCatch(
  suppressMessages(generate_event_attribution_configs(eopt10))$manifest,
  error = function(e) {
    cat("  multi-winter event error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(
  !is.null(ev10),
  "a registry with a multi-winter event generates (library matcher)"
)
if (!is.null(ev10)) {
  m10 <- ev10$events$multi_2015_17
  expect(
    identical(ev10$schema, 3L),
    "an event with sub-windows makes the manifest schema 3"
  )
  expect(
    identical(m10$peak, "2016-12-01") &&
      identical(m10$peak_season_year, 2015L) &&
      identical(unlist(m10$peak_season_years), c(2015L, 2016L)),
    "peak = the last peak (whole-event window), peak_season_year = the first winter, all winters listed"
  )
  sw <- m10$sub_windows
  expect(
    length(sw) == 2L &&
      identical(vapply(sw, `[[`, "", "id"), c("winter_2015", "winter_2016")),
    "two sub-windows, ids winter_<peak season year>"
  )
  # east lags 5/3/7 (max 7), horn 8/9/10 (max 10); KEN end Jul, TZA Nov, ETH Apr, others Dec.
  expect(
    identical(
      sw[[1]]$countries$KEN,
      list(start = "2014-09-01", end = "2017-02-01")
    ) &&
      identical(
        sw[[2]]$countries$KEN,
        list(start = "2017-03-01", end = "2018-02-01")
      ) &&
      identical(sw[[1]]$countries$TZA$end, "2017-06-01") &&
      identical(sw[[2]]$countries$TZA$end, "2018-06-01") &&
      identical(sw[[1]]$countries$BDI$end, "2016-07-01") &&
      identical(sw[[2]]$countries$BDI$start, "2016-08-01"),
    "per-country sub-windows: own season end on/after each peak + the group's max lag, consecutive"
  )
  expect(
    identical(
      sw[[1]]$groups$enso_east,
      list(start = "2014-09-01", end = "2017-06-01")
    ) &&
      identical(
        sw[[2]]$groups$enso_east,
        list(start = "2017-07-01", end = "2018-06-01")
      ) &&
      identical(sw[[2]]$groups$enso_horn$end, "2018-02-01"),
    "per-group sub-windows: the latest member's end, consecutive from the launch"
  )
  expect(
    identical(sw[[2]]$groups$enso_east$end, m10$groups$enso_east$horizon_end) &&
      identical(sw[[2]]$countries$KEN$end, m10$countries$KEN$firm_horizon_end),
    "the last sub-window ends at the whole-event window end (group and country)"
  )
  expect(
    is.null(ev10$events$elnino_2015$sub_windows) &&
      is.null(ev10$events$elnino_2015$peaks),
    "a single-window event in the same manifest carries no sub-window fields"
  )
  ca10 <- readRDS(m10$lenssen_cache)
  expect(
    identical(
      as.integer(unlist(ca10$meta$spec$primary_years)),
      c(2015L, 2016L)
    ) &&
      identical(colnames(ca10$stat_matrix)[1:2], c("2015", "2016")),
    "the counterfactual cache is labelled on both primary winters (every DJF of the window classified)"
  )
  expect(
    identical(unlist(m10$lenssen_winters$primary), c(2015L, 2016L)),
    "the manifest records both primary winters"
  )
  hdr <- readLines(m10$units$eh$factual_config, n = 30L)
  expect(
    any(grepl("sub-window winter_2016", hdr, fixed = TRUE)),
    "the config header lists the sub-windows"
  )
}
# Legacy matcher: a multi-winter event is refused (the legacy builder labels one DJF).
e10 <- errmsg(suppressMessages(generate_event_attribution_configs(modifyList(
  eopt10,
  list(
    events = "multi_2015_17",
    enso_counterfactual = NULL,
    config_outdir = file.path(root, "ev10l"),
    manifest = file.path(root, "ev10l.yml")
  )
))))
expect(
  !is.null(e10) &&
    grepl("needs enso_counterfactual.match 'library'", e10, fixed = TRUE),
  "legacy matcher + multi-winter event: refused"
)
# A v1-only selection from the v2 registry: schema 2, the [9] manifest's events.
ev10b <- tryCatch(
  suppressMessages(generate_event_attribution_configs(modifyList(
    eopt9,
    list(
      events_registry = reg10,
      config_outdir = file.path(root, "ev10b"),
      manifest = file.path(root, "ev10b.yml")
    )
  )))$manifest,
  error = function(e) NULL
)
expect(
  !is.null(ev10b) &&
    !is.null(ev9) &&
    identical(ev10b$schema, 2L) &&
    identical(
      gsub("/ev10b/", "/ev_bc/", yaml::as.yaml(ev10b$events), fixed = TRUE),
      yaml::as.yaml(ev9$events)
    ),
  "single-window events: schema 2, events identical to the v1 registry's ([9])"
)
# Event-level constraint scope: the same multi-winter event with
# other_winters: any in its registry entry. Its cache is built under "any" (the
# deploy's other settings unchanged), the later read winters still classified,
# and the manifest records the override; the event without it is unchanged.
reg10c <- file.path(root, "events_v2c.yml")
r10c <- yaml::read_yaml(reg10)
r10c$events$multi_2015_17$other_winters <- "any"
yaml::write_yaml(r10c, reg10c)
ev10c <- tryCatch(
  suppressMessages(generate_event_attribution_configs(modifyList(
    eopt10,
    list(
      events_registry = reg10c,
      lenssen_dir = file.path(root, "lenssen10c"),
      config_outdir = file.path(root, "ev10c"),
      manifest = file.path(root, "ev10c.yml")
    )
  )))$manifest,
  error = function(e) {
    cat("  other_winters any error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(
  !is.null(ev10c),
  "a multi-winter event with other_winters: any generates"
)
if (!is.null(ev10c) && !is.null(ev10)) {
  m10c <- ev10c$events$multi_2015_17
  ca10c <- readRDS(m10c$lenssen_cache)
  ca10 <- readRDS(ev10$events$multi_2015_17$lenssen_cache)
  expect(
    identical(ca10c$meta$spec$other_winters, "any") &&
      identical(ca10$meta$spec$other_winters, "neutral"),
    "the event's cache spec carries other_winters any; without the key the deploy's neutral"
  )
  expect(
    identical(colnames(ca10c$stat_matrix), colnames(ca10$stat_matrix)) &&
      identical(unname(ca10c$stat_matrix), unname(ca10$stat_matrix)) &&
      identical(ca10c$oni_trajectories, ca10$oni_trajectories),
    "same pool and every read winter still classified (stat_matrix identical)"
  )
  oth <- setdiff(colnames(ca10$stat_matrix), c("2015", "2016"))
  expect(
    length(oth) >= 1L,
    sprintf("the window reads later winters (%s)", paste(oth, collapse = ", "))
  )
  expect(
    sum(ca10c$set_labels == "mixed", na.rm = TRUE) <=
      sum(ca10$set_labels == "mixed", na.rm = TRUE) &&
      all(
        ca10c$set_labels[
          ca10$set_labels != "mixed" & !is.na(ca10$set_labels)
        ] ==
          ca10$set_labels[ca10$set_labels != "mixed" & !is.na(ca10$set_labels)]
      ),
    "any: no label set by the primary winters changes; only 'mixed' rows from later winters are released"
  )
  expect(
    identical(m10c$lenssen_winters$other_winters, "any") &&
      is.null(ev10$events$multi_2015_17$lenssen_winters$other_winters) &&
      is.null(ev10c$events$elnino_2015$lenssen_winters$other_winters),
    "the manifest records other_winters only for the event that sets it"
  )
  expect(
    identical(
      yaml::as.yaml(ev10c$events$elnino_2015$groups),
      yaml::as.yaml(ev10$events$elnino_2015$groups)
    ),
    "the event without the key: same windows as before"
  )
}

# --- [11] --season_window -------------------------------------------------------
cat(
  "\n[11] --season_window: per-country troughs, one run per (unit, launch month, season year)\n"
)
# Registry unit 'eh' (east + horn), [9] panel calendar (KEN min_seas 8, TZA 12,
# ETH 5, others 1). KEN's occurrence panel has onsets every Jun-Sep 2014-2018
# (20 onsets, onset-rate trough Jan); nobody else has onsets (min_seas fallback).
utils::write.csv(preg, ppath, row.names = FALSE, quote = FALSE)
sv11 <- do.call(
  rbind,
  lapply(iso9, function(iso) {
    ms <- seq(as.Date("2014-01-01"), as.Date("2018-12-01"), by = "month")
    data.frame(
      gid = paste0(iso, "_1"),
      country = iso,
      t_start = ms,
      outbreak = as.integer(
        iso == "KEN" & as.integer(format(ms, "%m")) %in% 6:9
      )
    )
  })
)
psv <- file.path(root, "surv11.rds")
saveRDS(sv11, psv)
sw_opt <- modifyList(
  hopt,
  list(
    climate_panel = ppan,
    partition_by = "registry",
    partition_source = PSRC,
    season_window = TRUE,
    launch_floor_month = 5L,
    trough_rule = "onset_rate",
    trough_min_onsets = 10L,
    trough_min_stability = 0.8,
    trough_years = "",
    trough_smooth_months = 3L,
    trough_bootstrap = 50L,
    trough_seed = 1L,
    survival_panel = psv,
    config_outdir = file.path(root, "hc_sw"),
    manifest = file.path(root, "hc_sw.yml")
  )
)
hs <- tryCatch(
  suppressMessages(generate_hindcast_configs(sw_opt)),
  error = function(e) {
    cat("  season-window hindcast error: ", conditionMessage(e), "\n")
    NULL
  }
)
expect(!is.null(hs), "season-window hindcast generation runs")
if (!is.null(hs)) {
  c11 <- hs$countries
  expect(
    identical(c11$KEN$trough_month, 1L) &&
      identical(c11$KEN$trough_source, "onset_rate") &&
      identical(c11$KEN$n_onsets, 20L) &&
      identical(c11$TZA$trough_source, "min_seas") &&
      identical(c11$TZA$trough_month, 12L) &&
      identical(c11$ETH$trough_month, 5L),
    "troughs: KEN onset rate (Jan), the others their modal min_seas"
  )
  expect(
    identical(names(hs$units$eh$years), c("2015-05", "2015-06")),
    "unit eh, season 2015: one run per launch month (May for the floored starts, Jun for ETH)"
  )
  y5 <- hs$units$eh$years[["2015-05"]]
  y6 <- hs$units$eh$years[["2015-06"]]
  expect(
    identical(unlist(y6$scored), "ETH") &&
      !"ETH" %in% unlist(y5$scored) &&
      setequal(c(unlist(y5$scored), unlist(y6$scored)), iso9) &&
      identical(y5$launch_year, 2015L) &&
      identical(y6$launch, "2015-06-01"),
    "each member scored in exactly one run; runs carry launch, launch_year, scored"
  )
  expect(
    identical(c11$KEN$years[["2015-05"]]$season_end, "2016-01-01") &&
      identical(c11$KEN$years[["2015-05"]]$firm_horizon_end, "2016-04-01") &&
      identical(c11$TZA$years[["2015-05"]]$season_end, "2015-12-01") &&
      identical(c11$ETH$years[["2015-06"]]$firm_horizon_end, "2017-01-01") &&
      identical(c11$SDN$years[["2015-05"]]$firm_horizon_end, "2016-09-01"),
    "windows: KEN May..Jan + 3, TZA May..Dec, ETH Jun..May + 8, SDN May..Jan + 8"
  )
  expect(
    identical(y5$horizon_end, "2016-09-01") &&
      identical(y6$horizon_end, "2017-01-01"),
    "each run goes to the latest firm horizon of the countries it scores"
  )
  expect(
    isTRUE(c11$KEN$floored) &&
      identical(c11$KEN$months_cut, 3L) &&
      identical(c11$KEN$onset_share_cut, 0) &&
      !isTRUE(c11$ETH$floored) &&
      identical(c11$TZA$months_cut, 4L) &&
      is.null(c11$TZA$onset_share_cut),
    "floor reporting: start month cut, onset share of the cut months (none for TZA: no onsets)"
  )
  expect(
    is.null(c11$KEN$years[["2015-06"]]) && is.null(c11$ETH$years[["2015-05"]]),
    "a country's season entry sits only under the key of the run that scores it"
  )
  expect(
    identical(hs$season_end, "season_window") &&
      identical(hs$season_window$launch_floor_month, 5L) &&
      identical(hs$season_window$survival_panel, psv) &&
      nzchar(hs$season_window$onsets_digest),
    "manifest records the rule and its inputs"
  )
  expect(
    identical(hs$groups$enso_horn$years[["2015"]]$launch, "2015-05-01") &&
      identical(hs$groups$enso_horn$years[["2015"]]$horizon_end, "2017-01-01"),
    "group season-year entry: earliest launch .. latest firm horizon"
  )
  hz11 <- manifest_country_horizons(hs, "firm", year = "2015-06")
  expect(
    identical(names(hz11), "ETH") && identical(hz11[["ETH"]], "2017-01-01"),
    "sim_io: a run key returns only the countries scored in that run"
  )
  cfg11 <- yaml::read_yaml(y6$config)
  expect(
    basename(y6$config) == "hindcast_eh_2015-06.yml" &&
      identical(cfg11$simulator$horizon$start, "2015-06-01") &&
      identical(cfg11$simulator$filter$handoff_month, "2015-06-01") &&
      identical(cfg11$simulator$horizon$end, "2017-01-01") &&
      identical(cfg11$simulator$run_tag, "hc_tag_eh_201506") &&
      setequal(strsplit(cfg11$simulator$scope$countries, ",")[[1]], iso9),
    "config keyed (unit, launch month): whole unit simulated, handoff at its launch"
  )
  csv11 <- read.csv(file.path(root, "season_ends_hc_sw.csv"))
  expect(
    nrow(csv11) == length(iso9) &&
      all(
        c(
          "trough_source",
          "floored",
          "months_cut",
          "onset_share_cut",
          "run_key"
        ) %in%
          names(csv11)
      ),
    "season_ends CSV: one row per country x season year, with the floor reporting"
  )
}
e <- errmsg(suppressMessages(generate_hindcast_configs(modifyList(
  sw_opt,
  list(
    partition_by = "",
    config_outdir = file.path(root, "hc_sw2"),
    manifest = file.path(root, "hc_sw2.yml")
  )
))))
expect(
  !is.null(e) && grepl("registry", e),
  "--season_window without registry units is refused"
)
hs_ms <- suppressMessages(generate_hindcast_configs(modifyList(
  sw_opt,
  list(
    trough_rule = "min_seas",
    survival_panel = file.path(root, "absent.rds"),
    config_outdir = file.path(root, "hc_sw3"),
    manifest = file.path(root, "hc_sw3.yml")
  )
)))
expect(
  identical(hs_ms$countries$KEN$trough_month, 8L) &&
    identical(hs_ms$countries$KEN$trough_source, "min_seas") &&
    identical(names(hs_ms$units$eh$years), c("2015-05", "2015-06", "2015-09")),
  "--trough_rule min_seas reads no occurrence panel; KEN (trough Aug) launches in Sep"
)

unlink(ppath)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
