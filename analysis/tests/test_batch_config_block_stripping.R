# test_batch_config_block_stripping.R — verify generate_batch() emits only
# the fields relevant to the computation each generated config describes.
#
# Hard-fail gate (companion to test_size_batch_config_determinism.R, which
# covers determinism/idempotence but not structural content). Covers:
#
#   1. `simulator:` is never emitted (make_batch_configs.R never generates a
#      simulator batch; every simulator config in the repo is hand-authored).
#   2. Within the active model block (occurrence/size/duration/monthly), only
#      the active backend's own sub-block (brms/brms_cb/brms_kl)
#      survives — including for `_cb` backends, where the plain-backend
#      sibling (e.g. `brms` for `brms_cb`) must be ABSENT, not retained.
#   3. `dlnm:` is present iff the backend is a `_cb` variant.
#   4. Regression guard for the root-cause fix: when a `_cb`-backend batch
#      sets a field (re_uncorrelated) ONLY under the plain-backend override,
#      strip_backend_siblings() must merge it into the `_cb` block BEFORE
#      dropping the plain sibling from the written file — and the value must
#      still resolve correctly when the written file is round-tripped through
#      make_options_from_config() (config.R: .fill_model_block() sibling
#      inheritance for `_cb` backends).

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
})

source("analysis/utils/base.R")    # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/make_batch_configs.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

# generate_batch() writes into get_config_dir() / <batch>. Redirect to a
# tempdir for every call so nothing lands in the tracked analysis/configs/
# tree (same pattern as test_size_batch_config_determinism.R).
generate_batch_into_tmp <- function(batch_name) {
  tmpdir <- tempfile("cfgdir_")
  dir.create(tmpdir, recursive = TRUE)
  prev_fun <- get("get_config_dir", envir = .GlobalEnv, inherits = FALSE)
  assign("get_config_dir", function() tmpdir, envir = .GlobalEnv)
  on.exit(assign("get_config_dir", prev_fun, envir = .GlobalEnv), add = FALSE)
  generate_batch(batch_name)
  assign("get_config_dir", prev_fun, envir = .GlobalEnv)
  file.path(tmpdir, batch_name)
}

BACKEND_SUBBLOCKS <- c("brms", "brms_cb", "brms_kl")
TRACKS <- c("occurrence", "size", "duration", "monthly")

# Structural assertions shared by every fixture: simulator stripped, only the
# active backend's sub-block survives, dlnm present iff backend is `_cb`.
check_config_structure <- function(cfg, track, label) {
  expect(is.null(cfg$simulator),
         sprintf("%s: simulator: block stripped", label))

  blk <- cfg[[track]]
  expect(!is.null(blk), sprintf("%s: %s: block present", label, track))
  if (is.null(blk)) return(invisible())

  be <- blk$backend
  expect(!is.null(be), sprintf("%s: %s$backend resolved", label, track))
  if (is.null(be)) return(invisible())

  present_subblocks <- intersect(BACKEND_SUBBLOCKS, names(blk))
  expect(setequal(present_subblocks, be),
         sprintf("%s: only active backend sub-block '%s' present (found: %s)",
                 label, be, paste(present_subblocks, collapse = ", ")))
  expect(!is.null(blk[[be]]),
         sprintf("%s: %s$%s retained", label, track, be))

  expect(identical(!is.null(blk$dlnm), grepl("_(cb|kl)$", be)),
         sprintf("%s: dlnm present iff backend is a _cb/_kl variant (backend=%s)",
                 label, be))
  expect(identical(!is.null(blk$kernel), endsWith(be, "_kl")),
         sprintf("%s: kernel present iff backend is a _kl variant (backend=%s)",
                 label, be))

  # No other model block leaks through.
  for (other in setdiff(TRACKS, track)) {
    expect(is.null(cfg[[other]]),
           sprintf("%s: sibling model block '%s' absent", label, other))
  }
}

# --- 1. Real BATCH_REGISTRY fixtures, one per track/backend combination ----
cat("\n[1] structural stripping across real BATCH_REGISTRY fixtures\n")

FIXTURES <- list(
  list(batch = "run1_benchmark",                       track = "occurrence"),
  list(batch = "enso_cb_singlelag_occ",                 track = "occurrence"),
  list(batch = "size_run1_benchmark",                   track = "size"),
  list(batch = "size_enso_cb_singlelag",                track = "size"),
  list(batch = "prod_run_dur_grouped_enso_singlelag",   track = "duration"),
  list(batch = "mon_structure_compare",                 track = "monthly")
)

for (fx in FIXTURES) {
  d   <- generate_batch_into_tmp(fx$batch)
  fls <- setdiff(list.files(d, pattern = "\\.yml$", full.names = TRUE),
                 file.path(d, "manifest.yml"))
  expect(length(fls) > 0L, sprintf("%s: at least one config generated", fx$batch))
  for (f in fls) {
    cfg <- yaml::read_yaml(f)
    check_config_structure(cfg, fx$track, sprintf("%s/%s", fx$batch, basename(f)))
  }
  unlink(dirname(d), recursive = TRUE)
}

# --- 2. Synthetic duration/brms fixture (no real registry entry exists) ----
# grep-confirmed: the only duration BATCH_REGISTRY entry
# (prod_run_dur_grouped_enso_singlelag) is brms_cb; there is no plain-brms
# duration batch on record. Add a minimal synthetic entry (in-memory only,
# not written back to make_batch_configs.R) to cover this combination.
cat("\n[2] synthetic duration/brms fixture\n")

BATCH_REGISTRY[["test_duration_brms_synth"]] <- list(
  description = "Synthetic (test-only): duration/brms, no real registry entry exists for plain brms",
  grid        = list(`duration.backend` = "brms"),
  overrides   = list(country = "all", seed = 1234L)
)

d   <- generate_batch_into_tmp("test_duration_brms_synth")
fls <- setdiff(list.files(d, pattern = "\\.yml$", full.names = TRUE),
               file.path(d, "manifest.yml"))
expect(length(fls) > 0L, "test_duration_brms_synth: at least one config generated")
for (f in fls) {
  cfg <- yaml::read_yaml(f)
  check_config_structure(cfg, "duration", sprintf("test_duration_brms_synth/%s", basename(f)))
  expect(identical(cfg$duration$backend, "brms"),
         "test_duration_brms_synth: backend is brms")
}
unlink(dirname(d), recursive = TRUE)

# --- 3. Regression guard: plain-sibling-only field survives strip+merge ----
# Root-cause scenario: a `_cb`-backend batch
# sets re_uncorrelated ONLY via the plain-backend override key
# (occurrence.brms.re_uncorrelated), never under occurrence.brms_cb directly.
# strip_backend_siblings() must merge it into the brms_cb block BEFORE
# dropping brms from the written file -- otherwise the value is lost forever
# (a config file with the plain sibling already gone gives
# .fill_model_block()'s sibling-inheritance nothing to inherit from on the
# NEXT load).
cat("\n[3] regression guard: re_uncorrelated set only on plain sibling survives strip+merge\n")

BATCH_REGISTRY[["test_cb_sibling_reuncorrelated_synth"]] <- list(
  description = "Synthetic (test-only): occurrence/brms_cb, re_uncorrelated set only on the plain brms sibling",
  grid        = list(`occurrence.backend` = "brms_cb"),
  overrides   = list(
    country = "all",
    seed    = 1234L,
    `occurrence.brms.re_uncorrelated` = TRUE
  )
)

d   <- generate_batch_into_tmp("test_cb_sibling_reuncorrelated_synth")
fls <- setdiff(list.files(d, pattern = "\\.yml$", full.names = TRUE),
               file.path(d, "manifest.yml"))
expect(length(fls) == 1L, "test_cb_sibling_reuncorrelated_synth: exactly one config generated")

if (length(fls) == 1L) {
  f   <- fls[[1]]
  cfg <- yaml::read_yaml(f)
  check_config_structure(cfg, "occurrence", "test_cb_sibling_reuncorrelated_synth (written file)")

  # 3a. Written-file check: value merged onto brms_cb, plain sibling absent.
  # NB: use names() rather than `$brms` -- `$` partial-matches list names, so
  # `cfg$occurrence$brms` would silently resolve to `brms_cb` and mask a
  # regression here.
  expect(!("brms" %in% names(cfg$occurrence)),
         "written file: plain brms sibling absent")
  expect(isTRUE(cfg$occurrence$brms_cb$re_uncorrelated),
         "written file: re_uncorrelated merged onto brms_cb despite living only on the plain sibling in overrides")

  # 3b. Round-trip through the real merge path: the value must still resolve
  # after make_options_from_config() re-loads the (already-stripped) file.
  opt <- make_options_from_config(option_list = list(), cli_argv = c("--config", f))
  expect(isTRUE(opt$occurrence$brms_cb$re_uncorrelated),
         "round-trip: opt$occurrence$brms_cb$re_uncorrelated survives make_options_from_config() despite stripped plain sibling")
}
unlink(dirname(d), recursive = TRUE)

# 5. spec$drop_superseded_prior_sd drops the legacy prior_sd from the
# ACTIVE backend block only when both tier keys are set there; without the
# flag the written configs are unchanged.
cat("\n[E9] drop_superseded_prior_sd\n")
e9_spec <- function(tiers, flag) {
  bc <- list(sampling = 10L, prior_sd = "normal(0, 2.5)")
  if ("rc" %in% tiers)   bc$prior_sd_region_country <- "normal(0, 0.35)"
  if ("adm1" %in% tiers) bc$prior_sd_adm1 <- "normal(0, 0.35)"
  s <- list(description = "Synthetic (test-only): E9 prior_sd",
            grid = list(`duration.backend` = "brms_cb"),
            overrides = list(country = "all", seed = 1234L, `duration.brms_cb` = bc))
  if (!is.null(flag)) s$drop_superseded_prior_sd <- flag
  s
}
e9_cfg <- function(tiers, flag) {
  BATCH_REGISTRY[["test_e9_prior_sd_synth"]] <<- e9_spec(tiers, flag)
  d <- generate_batch_into_tmp("test_e9_prior_sd_synth")
  f <- setdiff(list.files(d, pattern = "\\.yml$", full.names = TRUE), file.path(d, "manifest.yml"))
  out <- list(cfg = yaml::read_yaml(f[[1]]), text = readLines(f[[1]]))
  unlink(dirname(d), recursive = TRUE)
  out
}
both_on   <- e9_cfg(c("rc", "adm1"), TRUE)
one_on    <- e9_cfg("rc", TRUE)
none_on   <- e9_cfg(character(0), TRUE)
both_off  <- e9_cfg(c("rc", "adm1"), NULL)
both_fls  <- e9_cfg(c("rc", "adm1"), FALSE)
expect(!("prior_sd" %in% names(both_on$cfg$duration$brms_cb)),
       "flag + both tiers set: prior_sd dropped from the active backend block")
expect(identical(both_on$cfg$duration$brms_cb$prior_sd_region_country, "normal(0, 0.35)") &&
         identical(both_on$cfg$duration$brms_cb$prior_sd_adm1, "normal(0, 0.35)"),
       "flag + both tiers set: the tiers stay")
expect(identical(one_on$cfg$duration$brms_cb$prior_sd, "normal(0, 2.5)"),
       "flag + one tier only: prior_sd kept")
expect(identical(none_on$cfg$duration$brms_cb$prior_sd, "normal(0, 2.5)"),
       "flag + no tier: prior_sd kept")
expect(identical(both_off$cfg$duration$brms_cb$prior_sd, "normal(0, 2.5)"),
       "no flag: prior_sd kept (both tiers set)")
expect(identical(both_off$text, both_fls$text),
       "flag FALSE writes the same file as no flag")
expect(identical(setdiff(both_off$text, both_on$text), "    prior_sd: normal(0, 2.5)") &&
         length(both_off$text) == length(both_on$text) + 1L,
       "the flag removes exactly the prior_sd line, nothing else")
BATCH_REGISTRY[["test_e9_prior_sd_synth"]] <- NULL

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
message("PASS  test_batch_config_block_stripping")
