# test_forecast_attrib_configs.R -- the attribution factual:
# utils/make_forecast_attrib_configs.R (deploy_prod_run.R --stage gen_attrib)
# and the results-side manifest checks (utils/results_core.R).
#
# HARD gate, quick-safe: synthetic temp configs and manifests only.
#   [1] one forecast_<unit>_attrib.yml per unit: the factual config with
#       horizon.start, filter.forecast_window.start, filter.handoff_month moved
#       to the counterfactual launch and `_attrib` appended to run_tag / run_id;
#       every other leaf unchanged;
#   [2] the attribution manifest: factual fields swapped, launches equal to the
#       counterfactual's, combined.factual suffixed, provenance keys, the same
#       counterfactual runs;
#   [3] rerun: equal configs left as they are; a differing existing config is
#       refused without --force and overwritten with it; dry-run writes nothing;
#   [4] counterfactual sets launching at different months are refused;
#   [5] aligned launches: no config, the manifest's own factual runs, stale
#       _attrib configs pruned;
#   [6] results_attrib_factual_differs() / results_attrib_manifest_check(): a
#       manifest generated from another forecast manifest, or naming other
#       counterfactual runs, is refused.
#
# Run: Rscript analysis/tests/test_forecast_attrib_configs.R

suppressPackageStartupMessages({ library(yaml) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/make_forecast_attrib_configs.R")
source("analysis/utils/results_core.R")
cat("# attribution factual configs\n")

td <- tempfile("attrib_"); dir.create(td)
cdir <- file.path(td, "configs"); dir.create(cdir)
sim_cfg <- function(tag, launch, extra = list()) {
  list(baseline = "1990_2019", seed = 20261002L,
       simulator = c(list(run_tag = tag, run_id = paste0("20261026_", tag),
                          horizon = list(start = launch, end = "2027-10-01"),
                          max_draws = 199L, n_replicates = 9950L,
                          forecast = list(param_draw_mode = "cycle"),
                          filter = list(Np = 10000L, smoothing_input = "generated_data/sim_spatpomp/hist/cache",
                                        forecast_window = list(start = launch, end = "2027-10-01"),
                                        handoff_month = launch)), extra))
}
wcfg <- function(cfg, name) { p <- file.path(cdir, name); yaml::write_yaml(cfg, p); p }
make_manifest <- function(fac_launch = "2026-09-01", cf_launch = c(u1 = "2026-05-01", u2 = "2026-05-01")) {
  units <- list()
  for (u in c("u1", "u2")) {
    fac <- wcfg(sim_cfg(paste0("fc_", u), fac_launch, list(enso = list(cfsv2_calibration = "cal.rds"))),
                paste0("forecast_", u, ".yml"))
    cf  <- wcfg(sim_cfg(paste0("cf_", u, "_set-neutral"), cf_launch[[u]]), paste0("counterfactual_", u, "_set-neutral.yml"))
    units[[u]] <- list(members = list(paste0("C", u)), groups = list(paste0("enso_", u)),
                       horizon_end = "2027-10-01", factual_config = fac,
                       factual_merged_dir = file.path(td, "sim", paste0("fc_", u), "merged"),
                       factual_launch = fac_launch, analog_sets = list("neutral"), primary_analog_set = "neutral",
                       counterfactual_sets = list(neutral = list(analog_set = "neutral", config = cf,
                                                                 merged_dir = "x", launch = cf_launch[[u]])),
                       counterfactual_config = cf, counterfactual_merged_dir = "x",
                       counterfactual_launch = cf_launch[[u]])
  }
  man <- list(schema = 2L, analog_sets = list("neutral"), primary_analog_set = "neutral",
              attribution_start = min(c(fac_launch, cf_launch)),
              combined = list(factual = file.path(td, "sim", "fc_combined", "merged")),
              groups = list(enso_u1 = list(horizon_end = "2027-10-01"), enso_u2 = list(horizon_end = "2027-10-01")),
              units = units)
  p <- file.path(td, "forecast_horizons.yml"); yaml::write_yaml(man, p); p
}
quiet <- function(expr) suppressMessages(expr)
err_of <- function(expr) tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
am_path <- file.path(td, "forecast_attribution.yml")

# ---------------------------------------------------------------- [1]
cat("[1] attribution configs\n")
mp <- make_manifest()
res <- quiet(generate_forecast_attrib_configs(list(manifest = mp, attrib_manifest = am_path)))
expect(identical(sort(basename(res$written)), c("forecast_u1_attrib.yml", "forecast_u2_attrib.yml")),
       "one forecast_<unit>_attrib.yml per unit, next to the factual configs")
fac <- yaml::read_yaml(file.path(cdir, "forecast_u1.yml"))
att <- yaml::read_yaml(file.path(cdir, "forecast_u1_attrib.yml"))
fa <- config_flatten_leaves(fac); fb <- config_flatten_leaves(att)
diff_leaves <- names(fa)[!vapply(names(fa), function(l) identical(fa[[l]], fb[[l]]), logical(1))]
expect(setequal(diff_leaves, c("simulator.run_tag", "simulator.run_id", "simulator.horizon.start",
                               "simulator.filter.forecast_window.start", "simulator.filter.handoff_month")) &&
         setequal(names(fa), names(fb)),
       "only run_tag, run_id and the three launch dates differ from the factual config")
expect(identical(att$simulator$run_tag, "fc_u1_attrib") && identical(att$simulator$run_id, "20261026_fc_u1_attrib"),
       "run_tag and run_id carry the _attrib suffix")
expect(identical(att$simulator$horizon$start, "2026-05-01") && identical(att$simulator$filter$handoff_month, "2026-05-01") &&
         identical(att$simulator$filter$forecast_window$start, "2026-05-01"),
       "the attribution config launches at the counterfactual handoff")
expect(identical(att$simulator$horizon$end, fac$simulator$horizon$end), "horizon end unchanged")

# ---------------------------------------------------------------- [2]
cat("[2] attribution manifest\n")
am <- yaml::read_yaml(am_path); man <- yaml::read_yaml(mp)
expect(identical(am$units$u1$factual_config, file.path(cdir, "forecast_u1_attrib.yml")) &&
         identical(am$units$u2$factual_launch, "2026-05-01"),
       "the units' factual fields point at the attribution configs")
expect(grepl("sim_spatpomp/fc_u1_attrib/merged", am$units$u1$factual_merged_dir, fixed = TRUE),
       "the attribution run's merged dir follows its run tag")
expect(identical(am$combined$factual, file.path(td, "sim", "fc_combined_attrib", "merged")),
       "combined.factual: the _attrib combined dir")
expect(identical(am$kind, "forecast_attribution") && identical(am$source_manifest, mp) &&
         identical(am$display_factual_launch, "2026-09-01") && identical(am$attribution_factual, "attrib") &&
         identical(am$attribution_start, "2026-05-01"),
       "provenance keys: kind, source manifest, displayed factual launch")
al <- attrib_arm_launches(am)
expect(al$factual == al$counterfactual && al$factual == as.Date("2026-05-01"),
       "manifest readers: factual and counterfactual launch together")
expect(identical(am$units$u1$counterfactual_sets, man$units$u1$counterfactual_sets),
       "the counterfactual runs are the forecast manifest's")

# ---------------------------------------------------------------- [3]
cat("[3] rerun policy\n")
mt <- file.mtime(file.path(cdir, "forecast_u1_attrib.yml"))
Sys.sleep(1.1)
res2 <- quiet(generate_forecast_attrib_configs(list(manifest = mp, attrib_manifest = am_path)))
expect(length(res2$written) == 0L && length(res2$unchanged) == 2L &&
         file.mtime(file.path(cdir, "forecast_u1_attrib.yml")) == mt,
       "equal configs are left as they are")
x <- yaml::read_yaml(file.path(cdir, "forecast_u1_attrib.yml")); x$simulator$filter$Np <- 500L
yaml::write_yaml(x, file.path(cdir, "forecast_u1_attrib.yml"))
e <- err_of(quiet(generate_forecast_attrib_configs(list(manifest = mp, attrib_manifest = am_path))))
expect(!is.null(e) && grepl("--force", e, fixed = TRUE), "a differing existing config is refused without --force")
expect(yaml::read_yaml(file.path(cdir, "forecast_u1_attrib.yml"))$simulator$filter$Np == 500L,
       "the refusal writes nothing")
quiet(generate_forecast_attrib_configs(list(manifest = mp, attrib_manifest = am_path, force = TRUE)))
expect(yaml::read_yaml(file.path(cdir, "forecast_u1_attrib.yml"))$simulator$filter$Np == 10000L,
       "--force overwrites it")
invisible(file.remove(file.path(cdir, "forecast_u2_attrib.yml")))
quiet(generate_forecast_attrib_configs(list(manifest = mp, attrib_manifest = am_path, dry_run = TRUE)))
expect(!file.exists(file.path(cdir, "forecast_u2_attrib.yml")), "dry-run writes nothing")

# ---------------------------------------------------------------- [4]
cat("[4] mixed counterfactual launches\n")
mp4 <- make_manifest(cf_launch = c(u1 = "2026-05-01", u2 = "2026-06-01"))
e <- err_of(quiet(generate_forecast_attrib_configs(list(manifest = mp4, attrib_manifest = am_path))))
expect(!is.null(e) && grepl("counterfactual arm", e, fixed = TRUE), "counterfactual sets at different months are refused")

# ---------------------------------------------------------------- [5]
cat("[5] aligned launches\n")
mp5 <- make_manifest(fac_launch = "2026-05-01")
stale <- file.path(cdir, "forecast_u1_attrib.yml")
expect(file.exists(stale), "(setup) an _attrib config from the earlier run exists")
res5 <- quiet(generate_forecast_attrib_configs(list(manifest = mp5, attrib_manifest = am_path)))
am5 <- yaml::read_yaml(am_path); man5 <- yaml::read_yaml(mp5)
expect(isTRUE(res5$aligned) && length(res5$written) == 0L, "aligned: no attribution config")
expect(identical(am5$units$u1$factual_config, man5$units$u1$factual_config) &&
         identical(am5$attribution_factual, "forecast"),
       "aligned: the manifest names the forecast manifest's own factual runs")
expect(!file.exists(stale), "aligned: stale _attrib configs pruned")

# ---------------------------------------------------------------- [6]
cat("[6] results-side manifest checks\n")
mp6 <- make_manifest()
quiet(generate_forecast_attrib_configs(list(manifest = mp6, attrib_manifest = am_path)))
man6 <- yaml::read_yaml(mp6); am6 <- yaml::read_yaml(am_path)
expect(isTRUE(results_attrib_factual_differs(man6, am6)) && !isTRUE(results_attrib_factual_differs(man6, man6)),
       "results_attrib_factual_differs(): TRUE for the attribution manifest, FALSE for the forecast one")
expect(is.null(err_of(results_attrib_manifest_check(man6, am6, mp6, am_path))), "a matching attribution manifest passes")
bad <- am6; bad$source_manifest <- "analysis/configs/other.yml"
e <- err_of(results_attrib_manifest_check(man6, bad, mp6, am_path))
expect(!is.null(e) && grepl("generated from", e, fixed = TRUE), "a manifest from another forecast manifest is refused")
bad <- am6; bad$units$u2$counterfactual_sets$neutral$config <- "other_cf.yml"; bad$units$u2$counterfactual_config <- "other_cf.yml"
e <- err_of(results_attrib_manifest_check(man6, bad, mp6, am_path))
expect(!is.null(e) && grepl("other counterfactual runs", e, fixed = TRUE), "other counterfactual runs are refused")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
