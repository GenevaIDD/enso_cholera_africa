# test_config_size_retired.R — HARD: the config surface enforces that the size
#   track is retired downstream of fitting.
#
# Usage:
#   Rscript analysis/tests/test_config_size_retired.R
#
# Details:
#   The simulator, smoothing, deploy and results phases model case counts with
#   the monthly (`mon`) track only. Checks the defaults [1], the refusals [2]
#   and what still resolves [3]. Pure config resolution: no fit, no data.

suppressPackageStartupMessages(library(yaml))

source("analysis/utils/base.R")
source("analysis/utils/config.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# size track retired downstream: config surface\n")

errmsg <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
tmpdir <- tempfile("size_retired_")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

# Resolve a config written from `cfg` (or given as YAML text) through the same
# entry point every simulator / smoothing driver uses.
resolve <- function(cfg = NULL, text = NULL, model = NULL) {
  p <- tempfile("cfg_", tmpdir, fileext = ".yml")
  if (!is.null(text)) {
    writeLines(text, p)
  } else {
    yaml::write_yaml(cfg, p)
  }
  errmsg(invisible(capture.output(
    make_options_from_config(list(), cli_argv = c("-C", p), model = model)
  )))
}
sim_cfg <- function(...) {
  list(
    simulator = utils::modifyList(
      list(backend = "spatpomp", occurrence_fit = "occ.rds"),
      list(...)
    )
  )
}

# --- [1] defaults ----------------------------------------------------------------
# get_default_options() / config_defaults.yml carry no simulator.size_fit,
# simulator.size_cb_meta, simulator.size, location_frame.mean_rate_table_path
# or deploy.fits.siz_dir; results.scope.fig1_tracks has no "size";
# results.inputs.lag_loo_batch is {default, occurrence, duration, monthly}.
cat("\n[1] defaults carry no size key\n")
d <- get_default_options()
yml <- yaml::read_yaml("analysis/config_defaults.yml")
for (src in list(
  list(name = "get_default_options()", x = d),
  list(name = "config_defaults.yml", x = yml)
)) {
  s <- src$x$simulator
  expect(
    !any(c("size_fit", "size_cb_meta", "size") %in% names(s)) &&
      !"mean_rate_table_path" %in% names(s$location_frame),
    sprintf(
      "%s: no simulator.size_fit / size_cb_meta / size / location_frame.mean_rate_table_path",
      src$name
    )
  )
  r <- src$x$results
  expect(
    !"size" %in% unlist(r$scope$fig1_tracks) &&
      setequal(
        unlist(r$scope$fig1_tracks),
        c("occurrence", "duration", "monthly")
      ),
    sprintf(
      "%s: results.scope.fig1_tracks = occurrence/duration/monthly",
      src$name
    )
  )
  expect(
    identical(
      sort(names(r$inputs$lag_loo_batch)),
      sort(c("default", "occurrence", "duration", "monthly"))
    ),
    sprintf(
      "%s: results.inputs.lag_loo_batch = {default, occurrence, duration, monthly}",
      src$name
    )
  )
}
expect(
  !"siz_dir" %in% names(d$deploy$fits) &&
    all(c("occ_dir", "dur_dir", "mon_dir") %in% names(d$deploy$fits)),
  "get_default_options()$deploy$fits: occ_dir / dur_dir / mon_dir, no siz_dir"
)
expect(
  identical(d$simulator$posterior_strategy, "paired"),
  "simulator.posterior_strategy defaults to \"paired\""
)

# --- [2] refusals ----------------------------------------------------------------
# make_options_from_config() REFUSES a simulator config that sets
# simulator.size_fit / size_cb_meta / size / location_frame.
# mean_rate_table_path, or a posterior_strategy other than "paired", naming
# the key.
cat("\n[2] retired simulator keys are refused, naming the key\n")
cases <- list(
  "simulator.size_fit" = sim_cfg(size_fit = "size.rds"),
  "simulator.size_cb_meta" = sim_cfg(size_cb_meta = "size_cb_meta.rds"),
  "simulator.size" = sim_cfg(size = list(admin_level = "ADM1-ADM2")),
  "simulator.location_frame.mean_rate_table_path" = sim_cfg(
    location_frame = list(mean_rate_table_path = "mean_rate_table.rds")
  )
)
for (k in names(cases)) {
  e <- resolve(cases[[k]])
  expect(
    !is.null(e) &&
      grepl(k, e, fixed = TRUE) &&
      grepl("retired", e, fixed = TRUE),
    sprintf("%s set -> refused, naming it", k)
  )
}
e <- resolve(sim_cfg(size_fit = "size.rds", size = list(admin_level = "ADM1")))
expect(
  !is.null(e) && grepl("simulator.size_fit, simulator.size", e, fixed = TRUE),
  "several retired keys are all named in one error"
)
e <- resolve(sim_cfg(posterior_strategy = "independent"))
expect(
  !is.null(e) &&
    grepl("posterior_strategy", e, fixed = TRUE) &&
    grepl("paired", e, fixed = TRUE),
  "posterior_strategy: independent -> refused (only \"paired\" remains)"
)

# --- [3] what still resolves -------------------------------------------------------
# It does NOT refuse a null (`~`) size_fit, `posterior_strategy: paired`, a
# minimal simulator config, or a FIT-side size config (top-level `size:` model
# block, 04_fit_siz_*), which stays supported.
cat(
  "\n[3] null keys, paired strategy and fit-side size configs still resolve\n"
)
e <- resolve(
  text = c(
    "simulator:",
    "  backend: spatpomp",
    "  occurrence_fit: occ.rds",
    "  size_fit: ~",
    "  size_cb_meta: ~",
    "  location_frame:",
    "    mean_rate_table_path: ~"
  )
)
expect(
  is.null(e),
  sprintf(
    "YAML `~` for the retired keys reads as NULL and passes (%s)",
    e %||% "ok"
  )
)
e <- resolve(sim_cfg(posterior_strategy = "paired"))
expect(
  is.null(e),
  sprintf("posterior_strategy: paired passes (%s)", e %||% "ok")
)
e <- resolve(sim_cfg(
  monthly_fit = "mon.rds",
  duration_fit = "dur.rds",
  location_frame = list(size_data_path = "size_data.rds")
))
expect(
  is.null(e),
  sprintf(
    "a mon/dur simulator config with a size_data_path passes (%s)",
    e %||% "ok"
  )
)

siz_cfgs <- Sys.glob(
  "analysis/configs/*size_grouped_enso_singlelag/*_sl-*_grp-*.yml"
)
# The size batches are deploy-generated configs: a code-release checkout
# (analysis/RELEASE) ships only the hand-written ones, so there is nothing to check.
if (!length(siz_cfgs) && file.exists("analysis/RELEASE")) {
  cat("  SKIP  release checkout: no fit-side size batch config\n")
}
if (!length(siz_cfgs) && !file.exists("analysis/RELEASE")) {
  expect(FALSE, "fit-side size batch configs exist to check")
}
if (length(siz_cfgs)) {
  e <- errmsg(invisible(capture.output(
    make_options_from_config(
      list(),
      cli_argv = c("-C", siz_cfgs[[1L]]),
      model = "size"
    )
  )))
  expect(
    is.null(e),
    sprintf(
      "a fit-side size config (%s) still resolves (%s)",
      basename(siz_cfgs[[1L]]),
      e %||% "ok"
    )
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
