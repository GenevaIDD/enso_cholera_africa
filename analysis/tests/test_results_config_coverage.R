# test_results_config_coverage.R — HARD: guard the manuscript results
#   config.
#
# Usage:
#   Rscript analysis/tests/test_results_config_coverage.R
#
# Details:
#   The `results:` block is the authoritative record of the output decisions
#   behind the manuscript figures. Guards two failures: manuscript.yml
#   drifting from get_default_options()$results (a key silently supplied
#   from config_defaults.yml, or an orphan key nobody consumes), and two
#   derived tables colliding on one artefact path. Pure config resolution --
#   no fit, no data, no network.

suppressPackageStartupMessages({
  library(yaml)
})

source("analysis/utils/base.R") # %||%, used by paths.R
source("analysis/utils/config.R")
source("analysis/utils/paths.R") # make_results_run_file()
# make_options_from_config() hard-stops if results.deploy_config is set but
# resolve_results_inputs_from_deploy() is absent (see config.R). manuscript.yml
# now points at a deploy config, so source deploy_core.R to make the resolver
# available — exactly as 09/10 do. Inert when deploy_config is "".
source("analysis/utils/deploy_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

CONFIG_PATH <- "analysis/configs/results/manuscript.yml"

# Every leaf path of a nested list, as dot-separated strings (a character
# vector, an unnamed list and an empty mapping each stay one leaf).
leaf_paths <- config_leaf_paths # config.R

defaults <- get_default_options()
shipped <- yaml::read_yaml(CONFIG_PATH)

cat("# manuscript results config coverage\n")

# --- 1. Leaf-path parity, both directions ------------------------------------
# Failure 1: the shipped run config drifts out of sync with the defaults — a
# key added to get_default_options()$results but never restated in
# analysis/configs/results/manuscript.yml would be silently supplied from
# config_defaults.yml, so reading manuscript.yml would no longer tell you
# the whole decision set (main analog set, excluded groups, attribution
# thresholds, IRR subsample, Figure-1 lag rule). The reverse (an orphan key
# in manuscript.yml) means a decision nobody consumes.
cat("\n[1] defaults <-> manuscript.yml leaf-path parity\n")

expect(!is.null(defaults$results), "get_default_options() has a results block")
expect(!is.null(shipped$results), "manuscript.yml has a results block")

def_leaves <- sort(leaf_paths(defaults$results))
ship_leaves <- sort(leaf_paths(shipped$results))

missing <- setdiff(def_leaves, ship_leaves)
orphan <- setdiff(ship_leaves, def_leaves)

expect(
  length(missing) == 0L,
  sprintf(
    "manuscript.yml restates every default key (missing: %s)",
    if (length(missing)) paste(missing, collapse = ", ") else "none"
  )
)
expect(
  length(orphan) == 0L,
  sprintf(
    "manuscript.yml has no orphan key (orphan: %s)",
    if (length(orphan)) paste(orphan, collapse = ", ") else "none"
  )
)
expect(
  setequal(names(defaults$results), names(shipped$results)),
  sprintf(
    "same sub-blocks (defaults: %s | shipped: %s)",
    paste(sort(names(defaults$results)), collapse = ","),
    paste(sort(names(shipped$results)), collapse = ",")
  )
)

# `redo_cache` is an operational override, not an output decision. It must stay
# a CLI-only flag on 09 — a YAML key would imply it is part of the record.
expect(
  !("redo_cache" %in% unlist(lapply(defaults$results, names))),
  "redo_cache is NOT a results config key (CLI-only operational override)"
)

# --- 2. manuscript.yml is self-contained -------------------------------------
cat("\n[2] resolving manuscript.yml performs zero default fills\n")

opt <- make_options_from_config(
  option_list = list(),
  cli_argv = c("--config", CONFIG_PATH)
)

expect(
  identical(opt$results, shipped$results),
  "resolved opt$results is identical to the raw YAML (nothing filled in)"
)
expect(
  is.null(opt$occurrence) && is.null(opt$size),
  "no model block is materialised (detected_model is NULL for a results config)"
)

# --- 2b. manuscript_synth.yml: the same contract -----------------------------
# The synthetic-run record is aligned by
# deploy_prod_run.R --stage align_results --results-config <it>, so it must
# restate every key (an absent inputs line is never inserted) and resolve
# with zero default fills, exactly as manuscript.yml.
SYNTH_PATH <- "analysis/configs/results/manuscript_synth.yml"
cat("\n[2b] manuscript_synth.yml: leaf parity + zero default fills\n")
ship_s <- yaml::read_yaml(SYNTH_PATH)
l_s <- sort(leaf_paths(ship_s$results))
expect(
  identical(sort(leaf_paths(defaults$results)), l_s),
  sprintf(
    "manuscript_synth.yml restates every default key, no orphan (missing: %s | orphan: %s)",
    paste(setdiff(def_leaves, l_s), collapse = ", "),
    paste(setdiff(l_s, def_leaves), collapse = ", ")
  )
)
opt_s <- make_options_from_config(
  option_list = list(),
  cli_argv = c("--config", SYNTH_PATH)
)
expect(
  identical(opt_s$results, ship_s$results),
  "manuscript_synth.yml resolves to the raw YAML (nothing filled in)"
)
expect(
  identical(
    ship_s$results$deploy_config,
    "analysis/configs/deploy/synth_run1.yml"
  ),
  "manuscript_synth.yml points at synth_run1"
)
expect(
  identical(ship_s$results$derivation$enso_phase$smoothing, "episode") &&
    identical(ship_s$results$derivation$rate_pop_year, "window") &&
    identical(ship_s$results$inputs$episode_national_zero_rule, "r3") &&
    identical(ship_s$results$inputs$episode_at_risk_start, "after_episode") &&
    identical(ship_s$results$inputs$episode_adm1_overlap, "censor_orphan"),
  "manuscript_synth.yml carries the manuscript choices (enso_phase episode, rate_pop_year window, r3 / ae / a1c)"
)

# --- 3. make_results_run_file() path distinctness -------------------------------
# Failure 2: two different derived tables collide on one artefact path, so a
# re-derivation for one track/set list silently overwrites another's.
cat("\n[3] make_results_run_file() paths are pairwise distinct\n")

SCHEMA <- c(
  enso_relation = "v1",
  attribution = "v5",
  smoothing = "v1",
  oni = "v1"
)
tracks <- shipped$results$scope$fig1_tracks
set_lists <- list(
  shipped$results$scope$main_set,
  c(shipped$results$scope$main_set, shipped$results$scope$supp_sets)
)

paths <- c(
  vapply(
    tracks,
    function(tr) {
      make_results_run_file(
        opt,
        "results_enso_relation",
        track = tr,
        schema = SCHEMA[["enso_relation"]],
        config_path = CONFIG_PATH
      )
    },
    character(1)
  ),
  vapply(
    set_lists,
    function(s) {
      make_results_run_file(
        opt,
        "results_attribution",
        sets = s,
        schema = SCHEMA[["attribution"]],
        config_path = CONFIG_PATH
      )
    },
    character(1)
  ),
  vapply(
    set_lists,
    function(s) {
      make_results_run_file(
        opt,
        "results_oni",
        sets = s,
        schema = SCHEMA[["oni"]],
        config_path = CONFIG_PATH
      )
    },
    character(1)
  ),
  # the other sets-keyed kinds: report_forecast and report_attrib_exposure share
  # their runs and set keys, so their names must stay distinct
  vapply(
    set_lists,
    function(s) {
      make_results_run_file(
        opt,
        "results_report_forecast",
        sets = s,
        schema = "v4",
        config_path = CONFIG_PATH
      )
    },
    character(1)
  ),
  vapply(
    set_lists,
    function(s) {
      make_results_run_file(
        opt,
        "results_report_attrib_exposure",
        sets = s,
        schema = "v1",
        config_path = CONFIG_PATH
      )
    },
    character(1)
  ),
  vapply(
    set_lists,
    function(s) {
      make_results_run_file(
        opt,
        "results_report_attrib_metrics",
        sets = s,
        schema = "v1",
        config_path = CONFIG_PATH
      )
    },
    character(1)
  ),
  make_results_run_file(
    opt,
    "results_smoothing",
    schema = SCHEMA[["smoothing"]],
    config_path = CONFIG_PATH
  ),
  # hindcast skill kinds (unkeyed bundles)
  vapply(
    c("hindcast_runs", "hindcast_occ", "hindcast_cases", "report_hindcast"),
    function(k) {
      make_results_run_file(
        opt,
        paste0("results_", k),
        schema = "v1",
        config_path = CONFIG_PATH
      )
    },
    character(1)
  )
)

expect(
  length(unique(paths)) == length(paths),
  sprintf("all %d track x sets x kind paths are distinct", length(paths))
)
# `outdir` defaults to "generated_data/" WITH a trailing slash, so file.path()
# yields "generated_data//results/...". That doubled separator is the house
# behaviour of every sibling builder (make_sim_run_file, paths.R:1011-1013) and is
# inert on disk; normalise rather than diverge from the siblings here.
expect(
  all(grepl("^generated_data/+results/manuscript/", paths)),
  "artefacts land in <outdir>/results/<config run tag>/"
)

bumped <- make_results_run_file(
  opt,
  "results_attribution",
  sets = set_lists[[1]],
  schema = "v99",
  config_path = CONFIG_PATH
)
expect(
  !identical(bumped, paths[length(tracks) + 1L]),
  "bumping the schema token changes the artefact path"
)

# A sets key must not depend on the order the sets were requested in, or the
# same derived table would be written twice under two names.
expect(
  identical(
    make_results_run_file(
      opt,
      "results_oni",
      sets = c("nina", "all"),
      schema = "v1",
      config_path = CONFIG_PATH
    ),
    make_results_run_file(
      opt,
      "results_oni",
      sets = c("all", "nina"),
      schema = "v1",
      config_path = CONFIG_PATH
    )
  ),
  "sets key is order-invariant"
)

expect(
  inherits(
    try(
      make_results_run_file(
        opt,
        "results_enso_relation",
        track = "occurrence",
        config_path = CONFIG_PATH
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  "omitting schema is an error (an unversioned filename cannot be invalidated)"
)

# --- 4. Registry entry merges recursively ------------------------------------
cat("\n[4] results is registered as a recursive, always-filled block\n")

expect(
  identical(.CONFIG_TOPLEVEL_BLOCKS$results$merge, "recursive"),
  "results merges recursively (the block is 2 levels deep)"
)
expect(
  identical(.CONFIG_TOPLEVEL_BLOCKS$results$fill_when_absent, TRUE),
  "results is materialised even when a config omits it"
)

tmp <- tempfile(fileext = ".yml")
on.exit(unlink(tmp), add = TRUE)
writeLines(c("results:", "  thresholds:", "    low_signal: 9"), tmp)
opt_min <- make_options_from_config(
  option_list = list(),
  cli_argv = c("--config", tmp)
)

expect(
  identical(opt_min$results$thresholds$low_signal, 9L) ||
    identical(opt_min$results$thresholds$low_signal, 9),
  "explicit nested results.thresholds.low_signal overrides the default"
)
expect(
  identical(
    opt_min$results$thresholds$prob_threshold,
    defaults$results$thresholds$prob_threshold
  ),
  "untouched nested sibling results.thresholds.prob_threshold still fills (.merge_recursive ran at depth 2)"
)
expect(
  identical(opt_min$results$scope$main_set, defaults$results$scope$main_set),
  "untouched sibling sub-block results.scope still fills"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
