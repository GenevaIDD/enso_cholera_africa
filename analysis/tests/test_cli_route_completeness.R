# test_cli_route_completeness.R — HARD: every enforce_options flag of a fit
#   worker routes through .CLI_ROUTE or is a shared top-level key.
#
# Usage:
#   Rscript analysis/tests/test_cli_route_completeness.R
#
# Details:
#   .route_cli_to_nested() routes only intersect(enforce_options, names(route))
#   (analysis/utils/config.R), so a flag with no route slot is silently
#   dropped: `--sampling 50` once ran 1000 draws with no error anywhere. The
#   shared .CLI_ROUTE registry closes that gap; this test guards it for every
#   (model, backend) pair. Pure config resolution: no fit, no data, no network.

source("analysis/utils/base.R") # `%||%`
source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# CLI route completeness (enforce_options -> .CLI_ROUTE)\n")

# Shared top-level keys applied directly by make_options_from_config()
# (analysis/utils/config.R, the `for (nm in c("baseline", "country", ...))`
# loop) — these are never routed through .CLI_ROUTE and need no entry there.
SHARED_TOPLEVEL_KEYS <- c(
  "baseline",
  "country",
  "TL",
  "TR",
  "seed",
  "outdir",
  "overwrite"
)

# Every 04_fit_<track>_<backend>.R worker script, with its (model, backend)
# pair as understood by .CLI_ROUTE and its literal enforce_options vector.
# Parsed by hand from each script's make_options_from_config() call, so this
# table itself drifts if a script's enforce_options changes without a matching
# update here — the intersect-based checks below still catch a missing route
# even if this table falls a little stale, but keep it in sync on edits.
WORKER_SCRIPTS <- list(
  list(
    script = "analysis/04_fit_occ_brms.R",
    model = "occurrence",
    backend = "brms",
    enforce = c(
      "country",
      "admin_level",
      "basehaz_poisson",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "brms_init",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_occ_brms_cb.R",
    model = "occurrence",
    backend = "brms_cb",
    enforce = c(
      "country",
      "admin_level",
      "basehaz_poisson",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "brms_init",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_occ_brms_kl.R",
    model = "occurrence",
    backend = "brms_kl",
    enforce = c(
      "country",
      "basehaz_poisson",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "brms_init",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_siz_brms.R",
    model = "size",
    backend = "brms",
    enforce = c("country", "baseline", "admin_level")
  ),
  list(
    script = "analysis/04_fit_siz_brms_cb.R",
    model = "size",
    backend = "brms_cb",
    enforce = c(
      "country",
      "baseline",
      "admin_level",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_dur_brms.R",
    model = "duration",
    backend = "brms",
    enforce = c(
      "country",
      "admin_level",
      "basehaz_duration",
      "basehaz_poisson",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "brms_init",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_dur_brms_cb.R",
    model = "duration",
    backend = "brms_cb",
    enforce = c(
      "country",
      "admin_level",
      "basehaz_duration",
      "basehaz_poisson",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "brms_init",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_mon_brms.R",
    model = "monthly",
    backend = "brms",
    enforce = c(
      "country",
      "baseline",
      "admin_level",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "overwrite"
    )
  ),
  list(
    script = "analysis/04_fit_mon_brms_cb.R",
    model = "monthly",
    backend = "brms_cb",
    enforce = c(
      "country",
      "baseline",
      "admin_level",
      "sampling",
      "warmup",
      "chains",
      "cores",
      "seed",
      "overwrite"
    )
  )
)

# Deliberate exclusions: a flag that is routable in principle but is known and
# accepted to have no .CLI_ROUTE entry. Empty today — every flag below routes
# cleanly once .CLI_ROUTE is in place. Keep the mechanism so a future
# deliberate exclusion doesn't have to bypass the gate outright.
KNOWN_UNROUTED <- character(0) # format: "model.backend.flag"

# --- 1. every WORKER_SCRIPTS entry matches the file on disk ------------------
# For each 04_fit_<track>_<backend>.R worker script, the enforce_options =
# c(...) call is parsed literally (no sourcing, no fitting). Every flag must be
# either (a) one of the shared top-level keys make_options_from_config()
# applies directly (these never need a .CLI_ROUTE entry), or (b) present in
# names(.CLI_ROUTE[[model]][[backend]]), including an explicit NULL mapping
# (e.g. admin_level), which list() retains as a name even though the value is
# NULL. A flag failing both is a silent-drop gap: routable, not shared, but
# with no route entry.
cat("\n[1] enforce_options literal in this table matches the script on disk\n")

.extract_enforce_options <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  m <- regmatches(
    txt,
    regexpr("enforce_options\\s*=\\s*c\\([^)]*\\)", txt, perl = TRUE)
  )
  if (length(m) == 0L || !nzchar(m)) {
    return(character(0))
  }
  inner <- sub("^enforce_options\\s*=\\s*c\\((.*)\\)$", "\\1", m, perl = TRUE)
  vals <- regmatches(inner, gregexpr('"[^"]*"', inner))[[1]]
  gsub('"', "", vals)
}

for (w in WORKER_SCRIPTS) {
  on_disk <- .extract_enforce_options(w$script)
  expect(
    setequal(on_disk, w$enforce),
    sprintf(
      "%s: table's enforce list matches script (on disk: %s | table: %s)",
      w$script,
      paste(sort(on_disk), collapse = ","),
      paste(sort(w$enforce), collapse = ",")
    )
  )
}

# --- 2. every non-shared flag has a .CLI_ROUTE entry -------------------------
cat("\n[2] every enforce_options flag routes or is a shared top-level key\n")

for (w in WORKER_SCRIPTS) {
  route <- .CLI_ROUTE[[w$model]][[w$backend]]
  expect(
    !is.null(route),
    sprintf(
      "%s.%s: .CLI_ROUTE has an entry for this (model, backend) pair",
      w$model,
      w$backend
    )
  )
  if (is.null(route)) {
    next
  }

  routable <- setdiff(w$enforce, SHARED_TOPLEVEL_KEYS)
  for (flag in routable) {
    key <- sprintf("%s.%s.%s", w$model, w$backend, flag)
    if (key %in% KNOWN_UNROUTED) {
      next
    }
    expect(
      flag %in% names(route),
      sprintf(
        "%s: '%s' has a .CLI_ROUTE entry (would otherwise be silently dropped)",
        key,
        flag
      )
    )
  }
}

# --- 3. KNOWN_UNROUTED is not stale -------------------------------------------
cat("\n[3] KNOWN_UNROUTED allowlist is not stale\n")

still_missing <- character(0)
for (w in WORKER_SCRIPTS) {
  route <- .CLI_ROUTE[[w$model]][[w$backend]]
  if (is.null(route)) {
    next
  }
  routable <- setdiff(w$enforce, SHARED_TOPLEVEL_KEYS)
  for (flag in routable) {
    key <- sprintf("%s.%s.%s", w$model, w$backend, flag)
    if (key %in% KNOWN_UNROUTED && flag %in% names(route)) {
      still_missing <- c(still_missing, key)
    }
  }
}
expect(
  length(still_missing) == 0L,
  sprintf(
    "no KNOWN_UNROUTED entry now routes cleanly (stale entries: %s)",
    paste(still_missing, collapse = ", ")
  )
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
