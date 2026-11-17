# test_size_batch_config_determinism.R — HARD: make_batch_configs.R produces
#   deterministic, idempotent size-batch outputs.
#
# Usage:
#   Rscript analysis/tests/test_size_batch_config_determinism.R
#
# Details:
#   The BATCH_REGISTRY entries drive every multi-config sensitivity batch. A
#   regression that reorders expand.grid columns, randomises hash / file
#   order, changes the filename token convention (DIM_ABBREV /
#   make_config_filename) or silently mutates defaults would break downstream
#   batch postprocessing (keyed on filenames) and force every cached fit to be
#   re-run.

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
})

source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/make_batch_configs.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

# generate_batch() writes into get_config_dir() / <batch>. We need to
# redirect it to a temp dir for the determinism check. The cleanest way is
# to monkey-patch get_config_dir() in the running session.
.with_tmp_config_dir <- function(expr) {
  tmpdir <- tempfile("cfgdir_")
  dir.create(tmpdir, recursive = TRUE)
  saved <- environment(generate_batch)$get_config_dir
  unlockBinding("get_config_dir", as.environment(asNamespace("base")))
  # Override the helper in the script's namespace (sourced, so global).
  assign("get_config_dir", function() tmpdir, envir = .GlobalEnv)
  on.exit(
    {
      assign(
        "get_config_dir",
        saved %||% function() here::here("analysis/configs/"),
        envir = .GlobalEnv
      )
      unlink(tmpdir, recursive = TRUE)
    },
    add = TRUE
  )
  force(expr)
  tmpdir
}

# Robust pure-tempdir generator: regenerate into a fresh tempdir each call.
generate_batch_into_tmp <- function(batch_name) {
  tmpdir <- tempfile("cfgdir_")
  dir.create(tmpdir, recursive = TRUE)
  # Stash & override
  prev_fun <- get("get_config_dir", envir = .GlobalEnv, inherits = FALSE)
  assign("get_config_dir", function() tmpdir, envir = .GlobalEnv)
  on.exit(assign("get_config_dir", prev_fun, envir = .GlobalEnv), add = FALSE)
  generate_batch(batch_name)
  assign("get_config_dir", prev_fun, envir = .GlobalEnv)
  file.path(tmpdir, batch_name)
}

# Read every YAML in a dir and return list(filename -> parsed list).
read_dir_yamls <- function(d) {
  fls <- list.files(d, pattern = "\\.yml$", full.names = TRUE)
  out <- lapply(fls, function(f) yaml::read_yaml(f))
  names(out) <- basename(fls)
  out
}

# Compare two YAML maps key-by-key. Returns TRUE if all keys / values match.
# Optional `ignore_keys` for fields that intentionally differ (e.g.,
# manifest$created).
yaml_equal <- function(a, b, ignore_keys = character()) {
  a_keys <- setdiff(names(a), ignore_keys)
  b_keys <- setdiff(names(b), ignore_keys)
  if (!identical(sort(a_keys), sort(b_keys))) {
    return(FALSE)
  }
  for (k in a_keys) {
    if (!identical(a[[k]], b[[k]])) return(FALSE)
  }
  TRUE
}

SIZE_BATCHES <- c(
  "size_run1_benchmark",
  "size_seasonal_sensitivity",
  "size_enso_lag_sensitivity",
  "size_re_structure"
)

# --- 1. Determinism: two runs of generate_batch produce identical files --
# Strategy: generate each size batch into a TEMP directory, then call
# generate_batch() a second time into another TEMP directory and verify
# byte-identical YAML content per filename.
cat("\n[1] determinism — two runs produce identical files\n")

for (b in SIZE_BATCHES) {
  d1 <- generate_batch_into_tmp(b)
  d2 <- generate_batch_into_tmp(b)

  # Same filename set.
  f1 <- sort(list.files(d1, pattern = "\\.yml$"))
  f2 <- sort(list.files(d2, pattern = "\\.yml$"))
  expect(identical(f1, f2), sprintf("%s: same filename set across runs", b))

  # Per-config byte-identical content (ignore manifest$created — depends
  # on first-seen mtime, never on grid).
  y1 <- read_dir_yamls(d1)
  y2 <- read_dir_yamls(d2)
  ok <- TRUE
  for (nm in names(y1)) {
    ignore <- if (nm == "manifest.yml") "created" else character()
    if (!yaml_equal(y1[[nm]], y2[[nm]], ignore_keys = ignore)) {
      ok <- FALSE
      cat(sprintf("    DIFF in %s/%s\n", b, nm))
    }
  }
  expect(ok, sprintf("%s: per-config YAML content stable across runs", b))

  unlink(c(dirname(d1), dirname(d2)), recursive = TRUE)
}

# --- 2. Idempotent regeneration over existing dir keeps content stable ---
cat("\n[2] idempotence — regenerating in-place yields the same content\n")

for (b in SIZE_BATCHES) {
  d <- generate_batch_into_tmp(b)
  y1 <- read_dir_yamls(d)

  # Now regenerate IN PLACE using the same temp dir.
  prev_fun <- get("get_config_dir", envir = .GlobalEnv, inherits = FALSE)
  assign("get_config_dir", function() dirname(d), envir = .GlobalEnv)
  generate_batch(b)
  assign("get_config_dir", prev_fun, envir = .GlobalEnv)

  y2 <- read_dir_yamls(d)
  ok <- identical(sort(names(y1)), sort(names(y2)))
  for (nm in intersect(names(y1), names(y2))) {
    ignore <- if (nm == "manifest.yml") "created" else character()
    if (!yaml_equal(y1[[nm]], y2[[nm]], ignore_keys = ignore)) ok <- FALSE
  }
  expect(
    ok,
    sprintf("%s: in-place regeneration is a no-op (modulo timestamp)", b)
  )

  unlink(dirname(d), recursive = TRUE)
}

# --- 3. Filename token convention sanity --------------------------------
# Catch a regression in DIM_ABBREV / make_config_filename. Each size batch
# should produce filenames that:
#   - start with the batch name
#   - end with .yml
#   - contain exactly the abbrev tokens for varying dims (one per varying
#     grid column)
cat("\n[3] filename token convention\n")

for (b in SIZE_BATCHES) {
  d <- generate_batch_into_tmp(b)
  fls <- sort(setdiff(list.files(d, pattern = "\\.yml$"), "manifest.yml"))
  # All start with batch name
  expect(
    all(startsWith(fls, b)),
    sprintf("%s: all filenames begin with batch name", b)
  )
  # All end with .yml
  expect(
    all(grepl("\\.yml$", fls)),
    sprintf("%s: all filenames end in .yml", b)
  )

  # Manifest claims n_configs == len(fls)
  mf <- yaml::read_yaml(file.path(d, "manifest.yml"))
  expect(
    identical(mf$n_configs, length(fls)),
    sprintf("%s: manifest$n_configs matches actual file count", b)
  )
  expect(
    identical(sort(unlist(mf$configs)), fls),
    sprintf("%s: manifest$configs lists every generated YAML", b)
  )

  unlink(dirname(d), recursive = TRUE)
}

# --- 4. Grid expansion — count of YAMLs == nrow(expand.grid) -------------
cat("\n[4] grid expansion produces exactly n_configs YAMLs\n")

for (b in SIZE_BATCHES) {
  spec <- BATCH_REGISTRY[[b]]
  expected_n <- nrow(do.call(
    expand.grid,
    c(spec$grid, list(stringsAsFactors = FALSE))
  ))
  d <- generate_batch_into_tmp(b)
  actual_n <- length(setdiff(
    list.files(d, pattern = "\\.yml$"),
    "manifest.yml"
  ))
  expect(
    identical(actual_n, expected_n),
    sprintf("%s: %d YAMLs == nrow(expand.grid(grid))", b, expected_n)
  )
  unlink(dirname(d), recursive = TRUE)
}

# --- 5. Committed batches still match the registry output ----------------
# Cross-check against the committed batches (analysis/configs/<batch>/) for
# the same filename set and key-by-key content equality (allowing for
# timestamp-only drift in manifest.yml).
# Drift guard: if BATCH_REGISTRY is edited without re-running generate_batch,
# the committed YAMLs go stale. This is a "soft" check — we expect equality
# for stable batches but accept that the test author may have just edited
# the registry and not yet regenerated.
cat("\n[5] committed-vs-registry parity (size_run1_benchmark)\n")

committed <- file.path(here::here("analysis/configs"), "size_run1_benchmark")
if (dir.exists(committed)) {
  d <- generate_batch_into_tmp("size_run1_benchmark")
  c1 <- read_dir_yamls(committed)
  c2 <- read_dir_yamls(d)
  # Same per-config keys (ignore manifest$created drift only)
  same_names <- identical(sort(names(c1)), sort(names(c2)))
  expect(
    same_names,
    "size_run1_benchmark: committed filename set == registry-generated set"
  )
  if (same_names) {
    n_ok <- 0L
    for (nm in names(c1)) {
      ignore <- if (nm == "manifest.yml") "created" else character()
      if (yaml_equal(c1[[nm]], c2[[nm]], ignore_keys = ignore)) {
        n_ok <- n_ok + 1L
      }
    }
    expect(
      n_ok == length(c1),
      sprintf(
        "size_run1_benchmark: %d/%d configs byte-equal to committed",
        n_ok,
        length(c1)
      )
    )
  }
  unlink(dirname(d), recursive = TRUE)
} else {
  cat("  SKIP  committed size_run1_benchmark dir not found\n")
}

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) {
  quit(save = "no", status = 1)
}
message("PASS  test_size_batch_config_determinism")
