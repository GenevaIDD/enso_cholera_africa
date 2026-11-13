# test_results_exclude_groups.R — the stale enso_central exclusion.
#
# HARD gate, quick, no data (tempdir manifests and registries).
#   [1] defaults and the shipped results records exclude nothing
#   [2] results_derive.R::.resolve_exclude_groups(): an unknown name is an
#       error naming the registry; a registry group absent from the manifest is
#       a no-op with a message
#   [3] results_core.R::.RESULTS_STAMP_NORMALISE on scope.exclude_groups: a
#       stamp ["enso_central"] is fresh under [] when the manifest lacks the
#       group and stale when it has it; ["enso_west"] (a manifest group) stale
#       under []; an unreadable manifest compares strictly
#   [4] the derivation.enso_phase.smoothing back-compat entry keeps an
#       artefact stamped before the key fresh
#
# Run: Rscript analysis/tests/test_results_exclude_groups.R

suppressPackageStartupMessages({ library(yaml); library(dplyr) })
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/paths.R")
source("analysis/utils/results_core.R")
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
errmsg <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))

cat("# [1] configs\n")
expect(!length(get_default_options()$results$scope$exclude_groups), "default exclude_groups is empty")
for (f in sort(list.files("analysis/configs/results", pattern = "\\.yml$"))) {
  y <- yaml::read_yaml(file.path("analysis/configs/results", f))
  expect(!length(y$results$scope$exclude_groups), sprintf("%s excludes nothing", f))
}

# Fixtures: a kl-like registry (no enso_central) and two manifests.
root <- tempfile("excl_"); dir.create(root)
reg <- file.path(root, "enso_grouping_kl.csv")
writeLines(c("iso3,group", "KEN,enso_east", "NGA,enso_west", "ZAF,enso_south"), reg)
man_kl <- file.path(root, "man_kl.yml")
yaml::write_yaml(list(groups = list(enso_east = list(), enso_west = list()),
                      skipped_groups = list("enso_south")), man_kl)
man_leg <- file.path(root, "man_leg.yml")
yaml::write_yaml(list(groups = list(enso_east = list(), enso_central = list())), man_leg)

cat("# [2] registry guard\n")
env <- new.env()
sys.source("analysis/utils/results_derive.R", envir = env, keep.source = FALSE)
env$load_enso_grouping <- load_enso_grouping
man <- yaml::read_yaml(man_kl)
e <- errmsg(env$.resolve_exclude_groups("enso_central", man, reg))
expect(!is.null(e) && grepl("enso_central", e) && grepl(reg, e, fixed = TRUE),
       "a name that is no registry group is an error naming the registry")
expect(is.null(errmsg(env$.resolve_exclude_groups(character(0), man, file.path(root, "absent.csv")))),
       "nothing excluded: the registry is not read")
msg <- NULL
withCallingHandlers(env$.resolve_exclude_groups("enso_south", man, reg),
                    message = function(m) { msg <<- conditionMessage(m); invokeRestart("muffleMessage") })
expect(is.null(msg), "a skipped manifest group is accepted silently")
y2 <- yaml::read_yaml(man_kl); y2$skipped_groups <- NULL
withCallingHandlers(env$.resolve_exclude_groups("enso_south", y2, reg),
                    message = function(m) { msg <<- conditionMessage(m); invokeRestart("muffleMessage") })
expect(!is.null(msg) && grepl("no-op", msg), "a registry group absent from the manifest is a no-op with a message")

cat("# [3] manifest-aware stamp normaliser\n")
cfg_res <- function(excl, manifest)
  list(inputs = list(manifest = manifest), scope = list(exclude_groups = excl),
       derivation = list(n_irr_draws = 1L), thresholds = list(prob_threshold = 0.9))
stamp_of <- function(excl, manifest) list(results = cfg_res(excl, manifest))
d <- .results_decision_diff(stamp_of("enso_central", man_kl), cfg_res(list(), man_kl))
expect(!length(d), "stamp [enso_central] fresh under [] when the manifest lacks the group (kl registry)")
d <- .results_decision_diff(stamp_of("enso_central", man_leg), cfg_res(list(), man_leg))
expect(length(d) == 1L && grepl("scope.exclude_groups", d), "stale when the manifest has the group (legacy registry)")
d <- .results_decision_diff(stamp_of("enso_west", man_kl), cfg_res(list(), man_kl))
expect(length(d) == 1L && grepl("scope.exclude_groups", d), "stamp [enso_west] (a manifest group) stale under []")
d <- .results_decision_diff(stamp_of(c("enso_west", "enso_central"), man_kl), cfg_res("enso_west", man_kl))
expect(!length(d), "only the manifest groups of the exclusion are compared")
bad <- file.path(root, "absent.yml")
d <- .results_decision_diff(stamp_of("enso_central", bad), cfg_res(list(), bad))
expect(length(d) == 1L && grepl("scope.exclude_groups", d), "unreadable manifest: strict comparison")
expect(is.list(.RESULTS_STAMP_NORMALISE) && all(vapply(.RESULTS_STAMP_NORMALISE, is.function, logical(1))) &&
         identical(names(formals(.RESULTS_STAMP_NORMALISE[[1]])), c("x", "res")),
       ".RESULTS_STAMP_NORMALISE: named list of function(x, res)")

cat("# [4] derivation.enso_phase.smoothing back-compat\n")
want <- cfg_res(list(), man_kl); want$derivation$enso_phase <- list(smoothing = "threshold")
d <- .results_decision_diff(stamp_of(list(), man_kl), want)
expect(!length(d), "artefact stamped before the key is fresh under threshold")
want$derivation$enso_phase$smoothing <- "episode"
d <- .results_decision_diff(stamp_of(list(), man_kl), want)
expect(length(d) == 1L && grepl("enso_phase.smoothing", d), "and re-derives under episode")

unlink(root, recursive = TRUE)
cat(sprintf("\n[test_results_exclude_groups] %d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
