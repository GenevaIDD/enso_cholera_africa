# test_deploy_config_coverage.R — guard the prod_run deploy config.
#
# Hard-fail gate. The `deploy:` block is the authoritative record of a production
# deployment's decisions (which occurrence/duration/monthly batch feeds the
# shared best-lag archive, the neighbour-contagion naming + kernel gate, the
# spatPomp templates/run_tags, the forecast/hindcast config outdirs +
# shared manifests, and the hindcast launch grid). It mirrors the
# manuscript-results contract (test_results_config_coverage.R), so the same two
# failure modes apply:
#
#   1. The shipped deploy config drifts out of sync with the defaults — a key
#      added to get_default_options()$deploy but never restated in
#      analysis/configs/deploy/prod_run.yml would be silently supplied from
#      config_defaults.yml, so reading prod_run.yml would no longer tell you the
#      whole deployment decision set. The reverse (an orphan key) means a
#      decision nobody consumes.
#   2. The deploy provenance manifest path stops tracking the config it was
#      derived from.
#
# Tests:
#   1. Leaf-path parity between defaults and every deploy config, BOTH directions
#   2. each deploy config is genuinely self-contained (zero default fills occurred)
#   3. makeDeployRunFile() derives the run tag from the config stem
#   4. The registry entry actually merges recursively at depth >= 2
#   5. The season-end keys (hindcast.season_end, season_min_window_months,
#      event_attribution.season_end) and event_attribution.n_replicates reach
#      the generator opts
#
# Pure config resolution — no fit, no data, no network.
#
# Run: Rscript analysis/tests/test_deploy_config_coverage.R

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R")    # %||%, used by paths.R
source("analysis/utils/config.R")
source("analysis/utils/paths.R")   # makeDeployRunFile()

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

# Every shipped deploy config is guarded (prod_run, synth_run1, ...): sections
# [1]/[2] loop over the whole directory so a new deploy rollout cannot drift out
# of sync with the defaults. Section [3]'s stem assertions stay pinned to
# prod_run.yml.
CONFIG_PATHS <- sort(list.files("analysis/configs/deploy",
                                pattern = "\\.yml$", full.names = TRUE))
CONFIG_PATH  <- "analysis/configs/deploy/prod_run.yml"

# Every leaf path of a nested list, as dot-separated strings (a character
# vector, an unnamed list and an empty mapping each stay one leaf).
leaf_paths <- config_leaf_paths   # config.R

# deploy.sim.<arm> is a free-form overlay of dotted template leaves
# (--stage gen_templates): its four arm keys are required, its contents are
# opaque. Collapse every sim.<arm>.* path to sim.<arm> on both sides, the same
# way test_config_defaults_parity.R::drop_deploy() treats the deploy subtree.
collapse_opaque <- function(p) {
  unique(sub("^(sim\\.(shared|historical|forecast|counterfactual))(\\..*)?$", "\\1", p))
}

defaults  <- get_default_options()

cat("# prod_run deploy config coverage\n")
cat(sprintf("# guarding %d deploy config(s): %s\n",
            length(CONFIG_PATHS),
            paste(basename(CONFIG_PATHS), collapse = ", ")))

expect(!is.null(defaults$deploy), "get_default_options() has a deploy block")
expect(length(CONFIG_PATHS) >= 1L, "at least one analysis/configs/deploy/*.yml exists")

for (cfg_path in CONFIG_PATHS) {
  cfg_name <- basename(cfg_path)
  shipped  <- yaml::read_yaml(cfg_path)

  # --- 1. Leaf-path parity, both directions ----------------------------------
  cat(sprintf("\n[1] defaults <-> %s leaf-path parity\n", cfg_name))

  expect(!is.null(shipped$deploy), sprintf("%s has a deploy block", cfg_name))

  def_leaves  <- sort(collapse_opaque(leaf_paths(defaults$deploy)))
  ship_leaves <- sort(collapse_opaque(leaf_paths(shipped$deploy)))

  missing <- setdiff(def_leaves, ship_leaves)
  orphan  <- setdiff(ship_leaves, def_leaves)

  expect(length(missing) == 0L,
         sprintf("%s restates every default key (missing: %s)", cfg_name,
                 if (length(missing)) paste(missing, collapse = ", ") else "none"))
  expect(length(orphan) == 0L,
         sprintf("%s has no orphan key (orphan: %s)", cfg_name,
                 if (length(orphan)) paste(orphan, collapse = ", ") else "none"))
  expect(setequal(names(defaults$deploy), names(shipped$deploy)),
         sprintf("%s: same sub-blocks (defaults: %s | shipped: %s)", cfg_name,
                 paste(sort(names(defaults$deploy)), collapse = ","),
                 paste(sort(names(shipped$deploy)), collapse = ",")))
  # YAML 1.1 trap: a bare off / on / yes / no reads as a
  # logical, so a string-valued default leaf restated unquoted (e.g.
  # convergence.action: off) arrives as FALSE. Every leaf whose default is a
  # single string must be shipped as a string.
  str_leaves <- Filter(function(p) {
    v <- config_get_leaf(defaults$deploy, p)
    is.character(v) && length(v) == 1L
  }, setdiff(def_leaves, c("sim.shared", "sim.historical", "sim.forecast", "sim.counterfactual")))
  bad_type <- Filter(function(p) is.logical(config_get_leaf(shipped$deploy, p)), str_leaves)
  expect(length(bad_type) == 0L,
         sprintf("%s: no string-valued key shipped as a YAML boolean (bare off/on/yes/no: %s)", cfg_name,
                 if (length(bad_type)) paste(bad_type, collapse = ", ") else "none"))

  # --- 2. config is self-contained -------------------------------------------
  cat(sprintf("\n[2] resolving %s performs zero default fills\n", cfg_name))

  opt <- make_options_from_config(option_list = list(),
                                  cli_argv = c("--config", cfg_path))

  expect(identical(opt$deploy, shipped$deploy),
         sprintf("%s: resolved opt$deploy is identical to the raw YAML (nothing filled in)",
                 cfg_name))
  expect(is.null(opt$occurrence) && is.null(opt$size),
         sprintf("%s: no model block is materialised (detected_model is NULL for a deploy config)",
                 cfg_name))
}

# --- 3. makeDeployRunFile() tracks the config stem ---------------------------
cat("\n[3] makeDeployRunFile() derives the run tag from the config stem\n")

# Re-resolve against the pinned prod_run config (the loop above leaves `opt` at the
# last config). The run tag is derived from config_path's stem, not opt's content.
opt <- make_options_from_config(option_list = list(),
                                cli_argv = c("--config", CONFIG_PATH))

run_dir  <- makeDeployRunFile(opt, "deploy",          config_path = CONFIG_PATH)
manifest <- makeDeployRunFile(opt, "deploy_manifest", config_path = CONFIG_PATH)

expect(grepl("^generated_data/+deploy/prod_run/?$", run_dir),
       sprintf("deploy run dir lands at <outdir>/deploy/<config run tag> (got %s)", run_dir))
expect(grepl("^generated_data/+deploy/prod_run/deploy_manifest\\.qs$", manifest),
       sprintf("deploy manifest is deploy_manifest.qs under the run dir (got %s)", manifest))

other <- makeDeployRunFile(opt, "deploy_manifest",
                           config_path = "analysis/configs/deploy/synth_run1.yml")
expect(!identical(other, manifest),
       "a different config stem changes the deploy manifest path")

# --- 4. Registry entry merges recursively ------------------------------------
cat("\n[4] deploy is registered as a recursive, always-filled block\n")

expect(identical(.CONFIG_TOPLEVEL_BLOCKS$deploy$merge, "recursive"),
       "deploy merges recursively (the block is >=2 levels deep)")
expect(identical(.CONFIG_TOPLEVEL_BLOCKS$deploy$fill_when_absent, TRUE),
       "deploy is materialised even when a config omits it")

tmp <- tempfile(fileext = ".yml")
on.exit(unlink(tmp), add = TRUE)
writeLines(c("deploy:", "  contagion:", "    rho_km: 999"), tmp)
opt_min <- make_options_from_config(option_list = list(),
                                    cli_argv = c("--config", tmp))

expect(identical(opt_min$deploy$contagion$rho_km, 999L) ||
         identical(opt_min$deploy$contagion$rho_km, 999),
       "explicit nested deploy.contagion.rho_km overrides the default")
expect(identical(opt_min$deploy$contagion$scale,
                 defaults$deploy$contagion$scale),
       "untouched nested sibling deploy.contagion.scale still fills (.merge_recursive ran at depth 2)")
expect(identical(opt_min$deploy$fits, defaults$deploy$fits),
       "untouched sibling sub-block deploy.fits still fills")

# --- 5. Season-end keys reach the generators ----------------------------------
cat("\n[5] deploy_prod_run.R .hindcast_opt() / .event_opt() pass the season-end keys\n")

source("analysis/utils/deploy_core.R")   # load_deploy_config()
# The two builders are lifted out of the CLI script (it parses argv at top
# level) and evaluated against a deploy block, as the script does.
dpr <- parse("analysis/utils/deploy_prod_run.R", keep.source = FALSE)
fn_env <- new.env(parent = globalenv())
for (e in dpr) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c(".hindcast_opt", ".event_opt", ".deploy_nbr_partition_source"))
    eval(e, fn_env)
}
opts_for <- function(cfg_path) {
  fn_env$DEP <- load_deploy_config(cfg_path)
  fn_env$CONFIG_PATH <- cfg_path
  fn_env$opt <- list(outdir = "generated_data/")
  ev <- if (nzchar(as.character(fn_env$DEP$event_attribution$events %||% "")))
    fn_env$.event_opt() else NULL
  list(h = fn_env$.hindcast_opt(), e = ev)
}
# Fixtures: prod_run.yml with the keys a case needs edited, written to a temp
# file (the generators read nothing but the deploy block).
td5 <- tempfile("dep5"); dir.create(td5)
deploy_fixture <- function(name, edit) {
  y <- yaml::read_yaml(CONFIG_PATH)
  y$deploy <- edit(y$deploy)
  p <- file.path(td5, paste0(name, ".yml")); yaml::write_yaml(y, p)
  p
}
op <- opts_for(CONFIG_PATH)
expect(identical(op$h$season_end, "by_country_from_panel") &&
         identical(op$h$season_min_window_months, 6L),
       "season_end by_country_from_panel: hindcast generator gets it, min window 6")
expect(identical(op$e$season_end, "by_country_from_panel") &&
         identical(op$e$season_min_window_months, 6L),
       "season_end by_country_from_panel: event generator inherits the hindcast's (event_attribution.season_end \"\")")
p_grp <- deploy_fixture("by_group", function(d) {
  d$hindcast$season_end <- "by_group"; d$event_attribution$season_end <- ""; d
})
og <- opts_for(p_grp)
expect(identical(og$h$season_end, "by_group") && identical(og$e$season_end, "by_group"),
       "season_end by_group: by_group for both generators")
p5 <- deploy_fixture("by_group_event_override", function(d) {
  d$hindcast$season_end <- "by_group"
  d$event_attribution$season_end <- "by_country_from_panel"; d
})
o5 <- opts_for(p5)
expect(identical(o5$h$season_end, "by_group") && identical(o5$e$season_end, "by_country_from_panel"),
       "an explicit event_attribution.season_end wins over the hindcast's")
expect(identical(defaults$deploy$hindcast[["season_end"]], "by_group") &&
         identical(defaults$deploy$hindcast$season_min_window_months, 6L) &&
         identical(defaults$deploy$event_attribution[["season_end"]], ""),
       "defaults: hindcast.season_end by_group, season_min_window_months 6, event_attribution.season_end \"\"")
p0 <- deploy_fixture("event_reps_template", function(d) {
  d$event_attribution$n_replicates <- 0L; d
})
o0 <- opts_for(p0)
expect(identical(op$e$n_replicates, 10000L) && identical(o0$e$n_replicates, 0L) &&
         identical(defaults$deploy$event_attribution$n_replicates, 0L),
       "event_attribution.n_replicates reaches the event generator (10000 = 100 x 100; 0 = template's; default 0)")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
