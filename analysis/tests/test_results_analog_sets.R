# test_results_analog_sets.R — HARD: 09 side of the counterfactual analogs.
#
# Usage:
#   Rscript analysis/tests/test_results_analog_sets.R
#
# Details:
#   Quick, synthetic caches only. Pins analog_set_tables(), the
#   counterfactual cache identity assert, that the oni step reads the
#   counterfactual config's cache, and the run design table.

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/base.R")
source("analysis/utils/enso.R")
source("analysis/utils/sim_io.R")
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
errmsg <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
env <- new.env()
sys.source("analysis/utils/results_derive.R", envir = env, keep.source = FALSE)

oni_lkp <- c("2026-05-01" = 0.3)

# [1] analog_set_tables(): D29 rows of a legacy cache (one winter) and of a
#     library cache (several winters, `mixed`): n, model shares, seam error,
#     per-winter means, the min_set_size flag
cat("# [1] analog_set_tables\n")
set.seed(35)
traj <- matrix(rnorm(40 * 12, 0, 0.6), nrow = 40)
traj[, 1] <- 0.3 + rnorm(40, 0, 0.05)
leg <- list(
  oni_trajectories = traj,
  set_labels = rep(c("nina", "neutral", "elnino", "neutral"), 10),
  set_stat = rnorm(40),
  model = rep(c("A", "B"), each = 20),
  meta = list(
    horizon_start = "2026-05-01",
    set_definition = list(peak_season_year = 2026L)
  )
)
t1 <- env$analog_set_tables(leg, "forecast", "x.rds", oni_lkp)
s1 <- t1$sets
expect(
  identical(s1$set, c("all", "nina", "neutral", "elnino", "mixed")) &&
    identical(s1$n, c(40L, 10L, 20L, 10L, 0L)),
  "legacy: counts per set (mixed empty)"
)
expect(
  isTRUE(all.equal(
    s1$seam_error[s1$set == "neutral"],
    mean(traj[leg$set_labels == "neutral", 1]) - 0.3
  )),
  "seam error = set mean at launch - observed ONI"
)
expect(
  identical(s1$match[1], "legacy") &&
    identical(s1$primary_winter[1], "2026/27") &&
    all(s1$below_min_set_size[s1$set %in% c("nina", "elnino", "neutral")]) &&
    !s1$below_min_set_size[1],
  "legacy: primary DJF 2026/27, sets of 10 / 20 flagged under the default 30, `all` never flagged"
)
expect(
  nrow(t1$winters) == 4L && all(t1$winters$role == "primary"),
  "legacy: one winter per non-empty set"
)
expect(
  isTRUE(all.equal(
    sum(t1$model_shares$share[t1$model_shares$set == "all"]),
    1
  )),
  "model shares sum to 1"
)

lib <- leg
lib$set_labels[1:5] <- "mixed"
lib$stat_matrix <- cbind(`2026` = rnorm(40), `2027` = rnorm(40))
lib$meta$match <- "library"
lib$meta$digest <- "abc"
lib$meta$spec <- list(min_set_size = 5L)
lib$meta$set_definition$winters <- list(
  primary = list(2026L),
  other = list(2027L)
)
t2 <- env$analog_set_tables(lib, "elnino_2015", "y.rds", oni_lkp)
expect(
  identical(t2$sets$match[1], "library") &&
    identical(t2$sets$constrained_winters[1], "2026/27;2027/28") &&
    t2$sets$n[t2$sets$set == "mixed"] == 5L &&
    identical(t2$sets$cache_digest[1], "abc"),
  "library: two classified winters, mixed rows counted, digest carried"
)
expect(
  !any(t2$sets$below_min_set_size[t2$sets$set != "mixed"]),
  "library: min_set_size from the cache spec (5)"
)
w <- t2$winters[t2$winters$set == "all", ]
expect(
  identical(w$role, c("primary", "constrained")) &&
    isTRUE(all.equal(w$mean_djf_oni, unname(colMeans(lib$stat_matrix)))),
  "per-winter means and roles"
)
# Event-level constraint scope: other_winters "any"
# records the later read winter as read_unconstrained, with its mean; an explicit
# "neutral" is the default (identical tables).
lib_n <- lib
lib_n$meta$set_definition$winters$other_winters <- "neutral"
expect(
  identical(env$analog_set_tables(lib_n, "elnino_2015", "y.rds", oni_lkp), t2),
  "other_winters neutral (explicit) == absent: identical tables"
)
lib_a <- lib
lib_a$meta$set_definition$winters$other_winters <- "any"
t3 <- env$analog_set_tables(lib_a, "lanina_2020_23", "z.rds", oni_lkp)
wa <- t3$winters[t3$winters$set == "all", ]
expect(
  identical(wa$role, c("primary", "read_unconstrained")) &&
    isTRUE(all.equal(wa$mean_djf_oni, unname(colMeans(lib$stat_matrix)))),
  "other_winters any: the later winter is read_unconstrained, its mean still recorded"
)
expect(
  identical(t3$sets$constrained_winters[1], "2026/27") &&
    identical(t3$sets$n, t2$sets$n),
  "other_winters any: constrained_winters lists the primary only; set counts from the labels"
)

# [2] sim_io.R::assert_counterfactual_cache_identity(): no-op for a legacy
#     manifest; stops on a rebuilt cache and on a run that drew from another
#     (or a legacy) cache; runs without partition_meta.yml only reported
cat("# [2] assert_counterfactual_cache_identity\n")
tmp <- tempfile("cfid_")
dir.create(tmp)
cpath <- file.path(tmp, "cache.rds")
saveRDS(lib, cpath)
run_dir <- file.path(tmp, "sim_spatpomp", "cf_tag")
dir.create(file.path(run_dir, "merged"), recursive = TRUE)
runs <- list(list(config = "c.yml", merged_dir = file.path(run_dir, "merged")))
expect(
  identical(assert_counterfactual_cache_identity(NULL, runs), 0L),
  "legacy manifest (no record): no-op"
)
yaml::write_yaml(
  list(enso_counterfactual = list(digest = "abc")),
  file.path(run_dir, "partition_meta.yml")
)
expect(
  identical(
    assert_counterfactual_cache_identity(
      list(cache_path = cpath, digest = "abc"),
      runs
    ),
    1L
  ),
  "matching cache and run: checked"
)
e <- errmsg(assert_counterfactual_cache_identity(
  list(cache_path = cpath, digest = "zzz"),
  runs
))
expect(
  !is.null(e) && grepl("rebuilt after the runs", e),
  "manifest digest != cache digest: refused"
)
yaml::write_yaml(
  list(enso_continuation = list(method = "analog")),
  file.path(run_dir, "partition_meta.yml")
)
e <- errmsg(assert_counterfactual_cache_identity(
  list(cache_path = cpath, digest = "abc"),
  runs
))
expect(
  !is.null(e) && grepl("legacy cache", e),
  "a run without a library-cache record: refused"
)
unlink(file.path(run_dir, "partition_meta.yml"))
expect(
  identical(
    suppressMessages(assert_counterfactual_cache_identity(
      list(cache_path = cpath, digest = "abc"),
      runs
    )),
    0L
  ),
  "run without partition_meta.yml here: reported, not failed"
)
e <- errmsg(assert_counterfactual_cache_identity(
  list(cache_path = file.path(tmp, "gone.rds"), digest = "abc"),
  runs
))
expect(!is.null(e) && grepl("missing", e), "missing cache file: refused")
unlink(tmp, recursive = TRUE)

# [3] the oni step reads the counterfactual config's cache_path (no
#     hard-coded LENSSEN_CACHE read left) and lists the caches as inputs
cat("# [3] the oni step reads the config's cache\n")
oni_src <- paste(deparse(env$derive_oni_envelope), collapse = "\n")
expect(
  !grepl("ctx$LENSSEN_CACHE", oni_src, fixed = TRUE) &&
    grepl("lcf$cache_path", oni_src, fixed = TRUE),
  "derive_oni_envelope samples from lenssen_counterfactual_args()$cache_path"
)
ens_src <- paste(deparse(env$ensure_oni), collapse = "\n")
expect(
  grepl("caches", ens_src, fixed = TRUE),
  "ensure_oni lists the Lenssen caches as inputs"
)

cat("# [4] run design (report_analog_sets v2 run_design, table D34)\n")
cfg <- list(
  simulator = list(
    n_replicates = 5000L,
    max_draws = 200L,
    forecast = list(param_draw_mode = "cycle", handoff_weighting = "shuffled"),
    enso = list(clamp_to_knots = TRUE, clamp_margin = 0.25)
  )
)
lib <- list(
  oni_trajectories = matrix(0, 1000, 3),
  model = rep(LETTERS[1:10], each = 100),
  meta = list(
    match = "library",
    n_analogs_per_model = 100L,
    n_models = 10L,
    n_total_analogs = 1000L
  )
)
r <- env$run_design_row("counterfactual", cfg, lib)
expect(
  identical(r$n_replicates, 5000L) &&
    identical(r$max_draws, 200L) &&
    identical(r$reps_per_draw, 25) &&
    identical(r$draw_mode, "cycle") &&
    identical(r$handoff, "shuffled") &&
    identical(r$n_analogs_per_model, 100L) &&
    identical(r$n_models, 10L) &&
    identical(r$clamp_margin, 0.25) &&
    identical(r$analog_match, "library"),
  "cycle 200 x 25, shuffled handoff, library pool 100 per model x 10 models, clamp 0.25"
)
h <- env$run_design_row(
  "hindcast",
  list(
    simulator = list(
      n_replicates = 5000L,
      max_draws = 200L,
      forecast = list(param_draw_mode = "random")
    )
  )
)
expect(
  is.na(h$reps_per_draw) &&
    identical(h$draw_mode, "random") &&
    identical(h$handoff, "uniform") &&
    is.na(h$n_models) &&
    is.na(h$clamp_margin) &&
    identical(h$clamp_to_knots, FALSE),
  "random draws: no fixed reps per draw; config defaults (uniform handoff, clamp off); no pool"
)
sm <- env$run_design_row(
  "smoothing",
  list(simulator = list(n_replicates = 100L)),
  handoff = FALSE
)
expect(
  identical(sm$max_draws, 200L) &&
    is.na(sm$handoff) &&
    identical(sm$draw_mode, "random"),
  "smoothing: max_draws default 200, no handoff"
)
bl <- function(spread) {
  env$run_design_row(
    "forecast",
    list(
      simulator = list(
        n_replicates = 400L,
        max_draws = 200L,
        forecast = list(param_draw_mode = "", spread_param_draws = spread)
      )
    )
  )
}
expect(
  identical(bl(TRUE)$draw_mode, "cycle") &&
    identical(bl(TRUE)$reps_per_draw, 2) &&
    identical(bl(FALSE)$draw_mode, "inherited"),
  "blank param_draw_mode as the simulator: spread_param_draws ? cycle : inherited"
)
na <- env$run_design_row("events", NULL)
expect(
  is.na(na$n_replicates) && is.na(na$draw_mode) && is.na(na$clamp_to_knots),
  "absent arm: NA row"
)
legc <- leg
legc$meta$n_models <- NULL
expect(
  identical(env$run_design_row("events", cfg, legc)$n_models, 2L),
  "legacy cache: models counted from $model"
)
# The fixture manifests (analysis/tests/fixtures/) and the run configs of their
# first unit (fixtures/run_configs/): one config per arm, caches absent here.
suppressPackageStartupMessages(source("analysis/utils/results_core.R"))
FX <- "analysis/tests/fixtures"
dep_fx <- tempfile("deploy_fx_", fileext = ".yml")
yaml::write_yaml(
  list(
    deploy = list(
      templates = list(
        historical = file.path(
          FX,
          "run_configs",
          "sim_historical_smoothing.yml"
        )
      )
    )
  ),
  dep_fx
)
ctx <- list2env(
  list(
    MAN = load_manifest(file.path(FX, "forecast_horizons_fixture.yml")),
    MAN_ATTRIB = load_manifest(file.path(
      FX,
      "forecast_attribution_fixture.yml"
    )),
    EXCLUDE_GROUPS = character(0),
    MAIN_SET = "neutral",
    EVENT_MANIFEST_PATH = file.path(
      FX,
      "event_attribution_horizons_fixture.yml"
    ),
    HINDCAST_MANIFEST_PATH = file.path(FX, "hindcast_horizons_fixture.yml"),
    RES = list(deploy_config = dep_fx)
  ),
  parent = emptyenv()
)
rd <- env$run_design_configs(ctx)
expect(
  identical(rd$arm, env$RUN_DESIGN_ARMS) &&
    all(!is.na(rd$config)) &&
    all(file.exists(rd$config)),
  "fixture: a run config for every arm"
)
expect(
  grepl("set-neutral", rd$config[rd$arm == "counterfactual"]) &&
    grepl("_attrib", rd$config[rd$arm == "attribution_factual"]) &&
    grepl("hindcast", rd$config[rd$arm == "hindcast"]) &&
    grepl("historical", rd$config[rd$arm == "smoothing"]),
  "main-set counterfactual, attribution factual, hindcast and smoothing configs"
)
td <- env$run_design_table(rd)
expect(
  nrow(td) == 6L &&
    td$n_replicates[td$arm == "forecast"] == 9950L &&
    td$max_draws[td$arm == "forecast"] == 199L &&
    td$reps_per_draw[td$arm == "counterfactual"] == 50 &&
    td$n_replicates[td$arm == "smoothing"] == 100L &&
    is.na(td$handoff[td$arm == "smoothing"]),
  "fixture: 199 x 50 = 9,950 forward replicates, smoothing 100"
)

cat(sprintf("\n[test_results_analog_sets] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
