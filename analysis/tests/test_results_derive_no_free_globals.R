# test_results_derive_no_free_globals.R — HARD: results_derive.R reads the 09
#   script state only through `ctx`.
#
# Usage:
#   Rscript analysis/tests/test_results_derive_no_free_globals.R
#
# Details:
#   The derivers in utils/results_derive.R read the ~40 run variables (MAN,
#   RES, SETS, the input paths, ...) from the read-only context built by
#   results_context(), as ctx$NAME. A bare reference would not fail at load
#   time; it would fail (or, worse, silently pick up a same-named global)
#   only when its branch runs -- e.g. a track not in
#   results.scope.fig1_tracks, which the artefact-equivalence gate never
#   exercises. codetools sees every branch. Pure parse -- no data, no fit.

suppressPackageStartupMessages(library(codetools))
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# results_derive.R: former 09 globals are read only via ctx\n")
env <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env)
sys.source("analysis/utils/results_derive.R", env)
ctx_vars <- get(".RESULTS_CTX_VARS", env)
expect(
  length(ctx_vars) >= 40L,
  sprintf(
    "results_context() exports the former globals (%d names)",
    length(ctx_vars)
  )
)

fns <- Filter(function(n) is.function(get(n, env)), ls(env, all.names = TRUE))
fns <- setdiff(fns, c("%||%", "results_context")) # results_context() DEFINES them (locals)
offenders <- character()
for (f in fns) {
  free <- intersect(findGlobals(get(f, env), merge = TRUE), ctx_vars)
  if (length(free)) {
    offenders <- c(offenders, sprintf("%s: %s", f, paste(free, collapse = ",")))
  }
}
expect(
  length(fns) >= 40L,
  sprintf("scanned every function in results_derive.R (%d)", length(fns))
)
expect(
  length(offenders) == 0L,
  sprintf(
    "no function reads a former 09 global directly (offenders: %s)",
    if (length(offenders)) paste(offenders, collapse = "; ") else "none"
  )
)

# results_context() itself returns exactly the advertised names.
ctx_fn_locals <- all.names(body(get("results_context", env)))
expect(
  all(setdiff(ctx_vars, "opt") %in% ctx_fn_locals),
  "every exported name is assigned inside results_context()"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
