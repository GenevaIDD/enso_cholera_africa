# test_deploy_episode_key.R -- `deploy.episodes.key`, the rollout's single
# declaration of the cross-track episode definition.
#
# HARD gate (see analysis/tests/run_all.R).
#
# WHAT THIS GUARDS
# ----------------
# 11_deploy_prod_run.R threads deploy.episodes.key into (i) the fit collector,
# which refuses per-lag configs / fits / prep artefacts of another key and
# never prunes another episode definition's fits, (ii) the smoothing stage,
# which refuses a historical template declaring another key, and (iii) the
# forecast / hindcast generators, which write the key into every generated
# config (applyEpisodeKeyToSimConfig()) and pre-flight the prep artefacts it
# will read (preflightSimEpisodeInputs()). The default `lp_season` must leave
# every generated config byte-identical.
#
# Pure synthetic. No fits, no MCMC, no data/.
#
# Run: Rscript analysis/tests/test_deploy_episode_key.R

suppressPackageStartupMessages({ library(dplyr) })

source("analysis/utils/deploy_core.R")   # sources config/base/paths/provenance/loo_lags

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
throws <- function(expr, pattern = NULL) {
  r <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(r) && (is.null(pattern) || grepl(pattern, r))
}

cat("# deploy.episodes.key\n")

cat("\n[1] schema default\n")
expect(identical(get_default_options()$deploy$episodes$key, "lp_season"),
       "get_default_options()$deploy$episodes$key == 'lp_season'")
td <- tempfile("dep"); dir.create(td)
dy <- file.path(td, "d.yml")
yaml::write_yaml(list(deploy = list(fits = list(dest = "x"))), dy)
expect(identical(load_deploy_config(dy)$episodes$key, "lp_season"),
       "a deploy YAML without `episodes` inherits lp_season")
yaml::write_yaml(list(deploy = list(episodes = list(key = "gid_season"))), dy)
expect(identical(load_deploy_config(dy)$episodes$key, "gid_season"),
       "a deploy YAML can declare gid_season")

cat("\n[2] applyEpisodeKeyToSimConfig()\n")
tpl <- list(baseline = "1990_2019", simulator = list(run_tag = "t"))
expect(identical(applyEpisodeKeyToSimConfig(tpl, ""), tpl) &&
         identical(applyEpisodeKeyToSimConfig(tpl, NULL), tpl),
       "no key (standalone CLI) -> template untouched")
expect(identical(applyEpisodeKeyToSimConfig(tpl, "lp_season"), tpl),
       "lp_season on a template without `episodes` -> byte-identical config")
g <- applyEpisodeKeyToSimConfig(tpl, "gid_season")
expect(identical(g$episodes, list(key = "gid_season")), "gid_season is written")
tpl_gs <- c(tpl, list(episodes = list(key = "gid_season")))
expect(throws(applyEpisodeKeyToSimConfig(tpl_gs, "lp_season", "t"), "template declares"),
       "template declaring another key -> error")
expect(identical(applyEpisodeKeyToSimConfig(tpl_gs, "gid_season"), tpl_gs),
       "template already agreeing -> unchanged")
expect(throws(applyEpisodeKeyToSimConfig(tpl, "gid")), "unknown key -> error")

cat("\n[3] preflightSimEpisodeInputs()\n")
mk <- function(nm, key = NULL, digest = NA_character_) {
  p <- file.path(td, nm); saveRDS(1, p)
  if (!is.null(key))
    persistPrepPolicyStamp(p, prepPolicyStamp("crosswalk", "drop_collapse", .9, "derived",
                                              episode_key = key, episode_digest = digest))
  p
}
panel_lp <- mk("panel_lp.rds"); size_lp <- mk("size_lp.rds"); surv_lp <- mk("surv_lp.rds")
panel_gs <- mk("panel_gs.rds", "gid_season", "abc")
size_gs  <- mk("size_gs.rds",  "gid_season", "abc")
surv_gs  <- mk("surv_gs.rds",  "gid_season", "abc")
size_gs2 <- mk("size_gs2.rds", "gid_season", "zzz")
cfg <- function(panel, size, key = NULL)
  c(list(baseline = "1990_2019",
         simulator = list(climate_panel = panel,
                          location_frame = list(size_data_path = size))),
    if (!is.null(key)) list(episodes = list(key = key)))
expect(!throws(preflightSimEpisodeInputs(cfg(panel_lp, size_lp), surv_path = surv_lp)),
       "historical config + unstamped artefacts -> silent")
expect(!throws(preflightSimEpisodeInputs(cfg(panel_gs, size_gs, "gid_season"), surv_path = surv_gs)),
       "gid_season config + gid_season artefacts sharing a digest -> passes")
expect(throws(preflightSimEpisodeInputs(cfg(panel_lp, size_gs, "gid_season"), surv_path = surv_gs),
              "climate_panel"),
       "gid_season config + an unrebuilt climate panel -> error naming it")
expect(throws(preflightSimEpisodeInputs(cfg(panel_gs, size_gs2, "gid_season"), surv_path = surv_gs)),
       "digest mismatch between the size frame and the panels -> error")
expect(throws(preflightSimEpisodeInputs(cfg(panel_gs, size_gs, NULL), surv_path = surv_gs)),
       "historical config pointed at gid_season artefacts -> error")

cat("\n[4] the fit collector never prunes the other episode definition\n")
files <- c("outbreaks_survival_poisson_fit_base-1990_2019_grp-enso_grouping-enso_east_sl5.rds",
           "outbreaks_survival_poisson_fit_base-1990_2019_grp-enso_grouping-enso_east_sl5_ep-gs.rds")
expect(identical(.deploy_foreign_episode_artefacts(files, "lp_season"), files[2]),
       "lp_season deploy: the _ep-gs fit is foreign")
expect(identical(.deploy_foreign_episode_artefacts(files, "gid_season"), files[1]),
       "gid_season deploy: the untokenised fit is foreign")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
