# test_utils_module_purity.R — HARD: utils modules never source() each other.
#
# Usage:
#   Rscript analysis/tests/test_utils_module_purity.R
#
# Details:
#   Pure parse + sys.source of definition-only files; no fit, no data, no
#   network. analysis/utils/README.md's rule: sourced-only libraries are pure
#   definitions and never source() another module; a caller sources a module
#   plus its transitive closure (README dependency table). Second check: the
#   config modules each load alone after base.R.

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/tests/helpers_parse.R") # parse_data_of()

cat(
  "# analysis/utils modules never source() each other (documented exceptions only)\n"
)

# The documented exceptions are the `_core` kernels that source their own
# closure (deploy_core, results_core, smoothing_core), spatpomp_core_c.R's
# conditional self-source of spatpomp_core.R, and the non-step CLI tools
# (README section 3 / 5), which are run with Rscript and load their own
# dependencies.
# What counts: any source() / sys.source() CALL token in R's parse data,
# anywhere in the file (top level or inside a function). Comments and strings
# do not count. Ratchet: ALLOWED_SOURCERS is exactly today's set; a new
# sourcer fails, and an entry that stops sourcing must be deleted (a stale
# entry also fails).
ALLOWED_SOURCERS <- c(
  # _core kernels that source their own closure (README §2 / §5)
  "deploy_core.R",
  "results_core.R",
  "smoothing_core.R",
  # Csnippet kernel: conditional self-source of spatpomp_core.R
  "spatpomp_core_c.R",
  # non-step CLI tools (README §3 / §5)
  "deploy_prod_run.R",
  "make_batch_configs.R",
  "make_forecast_configs.R",
  "make_hindcast_configs.R",
  "make_event_attribution_configs.R",
  "make_forecast_attrib_configs.R",
  "make_pipeline_docs.R",
  "make_kl_grouping.R",
  "make_kl_sweep_groupings.R",
  "preprocess_paths.R",
  "compact_sim_outputs.R"
)

files <- list.files("analysis/utils", pattern = "\\.R$", full.names = TRUE)
expect(
  length(files) > 30L,
  sprintf("scanned analysis/utils (%d files)", length(files))
)

sourcers <- character(0)
for (f in files) {
  pd <- parse_data_of(f)
  if (is.null(pd)) {
    sourcers <- c(sourcers, sprintf("%s:unparseable", basename(f)))
    next
  }
  calls <- pd$text[pd$token == "SYMBOL_FUNCTION_CALL"]
  if (any(calls %in% c("source", "sys.source"))) {
    sourcers <- c(sourcers, basename(f))
  }
}

new_sourcers <- setdiff(sourcers, ALLOWED_SOURCERS)
expect(
  length(new_sourcers) == 0L,
  sprintf(
    "no utils module sources another outside the documented exceptions (offenders: %s)",
    if (length(new_sourcers)) paste(new_sourcers, collapse = ", ") else "none"
  )
)
stale <- setdiff(ALLOWED_SOURCERS, sourcers)
expect(
  length(stale) == 0L,
  sprintf(
    "no stale ALLOWED_SOURCERS entry (stale: %s)",
    if (length(stale)) paste(stale, collapse = ", ") else "none"
  )
)

# The config.R split modules: definitions only, and each loads alone after base.R ----
# The four config modules (config.R, formula.R, formula_legacy.R, priors.R)
# each load alone after base.R in a fresh environment: no hidden load-time
# dependency on a sibling.
for (m in c("config", "formula", "formula_legacy", "priors")) {
  f <- file.path("analysis/utils", paste0(m, ".R"))
  ex <- parse(f, keep.source = FALSE)
  defs_only <- all(vapply(
    ex,
    function(e) is.call(e) && as.character(e[[1]]) %in% c("<-", "="),
    TRUE
  ))
  expect(defs_only, sprintf("%s.R holds top-level assignments only", m))
  env <- new.env(parent = globalenv())
  sys.source("analysis/utils/base.R", env)
  loaded <- tryCatch(
    {
      sys.source(f, env)
      TRUE
    },
    error = function(e) conditionMessage(e)
  )
  expect(
    isTRUE(loaded),
    sprintf(
      "%s.R loads alone after base.R%s",
      m,
      if (isTRUE(loaded)) "" else paste0(" (", loaded, ")")
    )
  )
}

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
