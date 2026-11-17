# test_smoothing_cache_schema.R — HARD: the smoothing handoff-cache schema
#   versioning.
#
# Usage:
#   Rscript analysis/tests/test_smoothing_cache_schema.R
#
# Details:
#   Pure function tests: no data, no fit, only a tempdir to prove the glob
#   patterns against real files. The cache (hist_filter_tu=<tu>_*.qs, written
#   by 08c_smoothing_compile.R, read by load_smoothing_distribution() in
#   07b_simulator_spatpomp.R) carries no statename manifest: a cache written
#   before the .SP_UNIT_STATENAMES 9->11 change silently misaligns the
#   named-vector active/mso lookups in filter_traj_to_warmup() (NA, not an
#   error) unless callers check for it explicitly.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})

source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R") # bin_ttlo() -- used by filter_traj_to_warmup()
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/spatpomp_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

# 1. sm_cache_schema_token() / sm_cache_schema_token_parse(): the `_sv<NN>`
#    filename token and its round-trip parse, including the NA-on-legacy-
#    filename case.
cat("\n## sm_cache_schema_token() / sm_cache_schema_token_parse()\n")

expect(
  identical(
    sm_cache_schema_token(),
    sprintf("sv%02d", length(.SP_UNIT_STATENAMES))
  ),
  "token format is sv<NN>, NN = length(.SP_UNIT_STATENAMES)"
)
expect(
  .SM_CACHE_SCHEMA_VERSION == length(.SP_UNIT_STATENAMES),
  ".SM_CACHE_SCHEMA_VERSION tracks .SP_UNIT_STATENAMES length"
)

fn_new <- sprintf(
  "hist_filter_tu=USA_Np=100_%s_nrep=50_20180101_20191201.qs",
  sm_cache_schema_token()
)
fn_old <- "hist_filter_tu=USA_Np=100_nrep=50_20180101_20191201.qs" # pre-versioning (legacy)
fn_bad <- "hist_filter_tu=USA_Np=100_sv09_nrep=50_20180101_20191201.qs" # stale schema

expect(
  sm_cache_schema_token_parse(fn_new) == .SM_CACHE_SCHEMA_VERSION,
  "parse round-trips the current schema token"
)
expect(
  is.na(sm_cache_schema_token_parse(fn_old)),
  "legacy (token-less) filename parses to NA, not a false-compatible version"
)
expect(
  sm_cache_schema_token_parse(fn_bad) == 9L,
  "parse extracts an old/mismatched version verbatim (sv09)"
)

# The token's filename position (after `_Np=<Np>_`, not between `tu=` and
# `_Np=`) does not break the existing TU-name extraction regex in
# 07b_simulator_spatpomp.R's dry run or either discovery glob.
cat("\n## Filename-position regression: TU-name regex + discovery globs\n")

# 07b_simulator_spatpomp.R (dry-run TU discovery) -- non-greedy capture up to the first "_Np="
tu_extract <- function(fn) sub("^hist_filter_tu=(.*?)_Np=.*$", "\\1", fn)
expect(
  identical(tu_extract(fn_new), "USA"),
  "TU-name regex still extracts bare TU from a token-bearing filename"
)
expect(
  identical(tu_extract(fn_old), "USA"),
  "TU-name regex extracts bare TU from a legacy filename (baseline)"
)

expect(
  grepl("^hist_filter_tu=USA_.*\\.qs$", fn_new),
  "resolve-cache glob (hist_filter_tu=<tu>_*.qs) matches token-bearing filename"
)
expect(
  grepl("^hist_filter_tu=.*_Np=.*\\.qs$", fn_new),
  "discovery glob (hist_filter_tu=*_Np=*.qs) matches token-bearing filename"
)

# Prove the actual Sys.glob() patterns used in 07b against real files, not
# just string matching, in an isolated tempdir.
tmp <- tempfile("sm_cache_schema_glob_")
dir.create(tmp)
withr_like_cleanup <- function() unlink(tmp, recursive = TRUE, force = TRUE)
on.exit(withr_like_cleanup(), add = TRUE)
invisible(file.create(file.path(tmp, fn_new)))
invisible(file.create(file.path(tmp, fn_old)))

hits_by_tu <- Sys.glob(file.path(tmp, "hist_filter_tu=USA_*.qs"))
hits_by_np <- Sys.glob(file.path(tmp, "hist_filter_tu=*_Np=*.qs"))
expect(
  length(hits_by_tu) == 2L,
  "Sys.glob(hist_filter_tu=<tu>_*.qs) finds both old and new files"
)
expect(
  length(hits_by_np) == 2L,
  "Sys.glob(hist_filter_tu=*_Np=*.qs) finds both old and new files"
)

# 2. check_smoothing_cache_schema(): the loud "rebuild the smoothing cache"
#    gate that load_smoothing_distribution() calls before reading a cache.
cat("\n## check_smoothing_cache_schema() -- the load-time gate\n")

ok_res <- tryCatch(
  {
    check_smoothing_cache_schema(fn_new, "USA")
    "no_error"
  },
  error = function(e) conditionMessage(e)
)
expect(
  identical(ok_res, "no_error"),
  "current-schema cache passes the gate silently"
)

err_old <- tryCatch(
  {
    check_smoothing_cache_schema(fn_old, "USA")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(err_old) &&
    grepl("predates schema versioning", err_old) &&
    grepl("rebuild the smoothing cache", err_old),
  "legacy (token-less) cache fails loudly with a 'rebuild the smoothing cache' instruction"
)

err_bad <- tryCatch(
  {
    check_smoothing_cache_schema(fn_bad, "USA")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(err_bad) &&
    grepl("sv09", err_bad) &&
    grepl(sprintf("sv%02d", .SM_CACHE_SCHEMA_VERSION), err_bad) &&
    grepl("rebuild the smoothing cache", err_bad),
  "mismatched-version cache fails loudly naming both the found and expected versions"
)

# 3. filter_traj_to_warmup()'s explicit active/mso name-presence guard,
#    exercised against a synthetic stale (9-state) trajectory array so a
#    schema regression can't silently return NA again.
cat("\n## filter_traj_to_warmup(): active/mso name-presence guard\n")

frame_tu <- tibble::tibble(gid = c("A_1", "B_1"))
U <- nrow(frame_tu)

.mk_ft_array <- function(statenames) {
  nm <- as.vector(outer(statenames, seq_len(U), paste0))
  array(0, dim = c(length(nm), 1, 1), dimnames = list(nm, NULL, NULL))
}

# Fresh (11-state, current schema) array: guard must NOT fire.
ft_fresh <- .mk_ft_array(.SP_UNIT_STATENAMES)
err_fresh <- tryCatch(
  {
    filter_traj_to_warmup(ft_fresh, frame_tu, "2020-01-01")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  is.null(err_fresh) || !grepl("no 'active'/'mso' entries", err_fresh %||% ""),
  "fresh 11-state trajectory does not trip the active/mso guard"
)

# Stale (9-state) array: guard MUST fire, not silently return NA.
statenames_9 <- setdiff(.SP_UNIT_STATENAMES, c("active", "mso"))
ft_stale <- .mk_ft_array(statenames_9)
err_stale <- tryCatch(
  {
    filter_traj_to_warmup(ft_stale, frame_tu, "2020-01-01")
    NULL
  },
  error = function(e) conditionMessage(e)
)
expect(
  !is.null(err_stale) &&
    grepl("no 'active'/'mso' entries", err_stale) &&
    grepl("08a/08b/08c_smoothing", err_stale),
  "stale 9-state trajectory errors loudly instead of silently returning NA"
)

cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1)
}
