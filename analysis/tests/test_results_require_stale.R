# test_results_require_stale.R -- the opt-in RESULTS_ALLOW_STALE_INPUTS override
# of results_require() (utils/results_core.R): an artefact whose recorded input
# changed or disappeared since derivation stops by default and only warns with
# the variable set; a decision mismatch still stops either way.
#
# HARD gate, quick-safe: a throwaway artefact in tempdir(), no pipeline data.
#
# Run: Rscript analysis/tests/test_results_require_stale.R

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/utils/results_core.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
cat("# results_require() stale-input override\n")

td  <- tempfile("rr_stale_"); dir.create(td)
inp <- file.path(td, "panel.rds"); saveRDS(1:3, inp)
art <- file.path(td, "artefact.qs")
cfg <- file.path(td, "results.yml"); writeLines("results: {}", cfg)
res <- list(inputs = list(), scope = list(a = 1), derivation = list(), thresholds = list(m = 0.5))
results_write(list(x = 1), art, res, cfg, "v1", inputs = inp)

errs <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
Sys.unsetenv("RESULTS_ALLOW_STALE_INPUTS")
expect(!errs(results_require(art, cfg, res)), "fresh artefact loads")

Sys.setFileTime(inp, Sys.time() + 3600); saveRDS(1:4, inp)   # the input changes
expect(errs(results_require(art, cfg, res)), "changed input: stops by default")

Sys.setenv(RESULTS_ALLOW_STALE_INPUTS = "true")
w <- NULL
obj <- withCallingHandlers(results_require(art, cfg, res),
                           warning = function(c) { w <<- conditionMessage(c); invokeRestart("muffleWarning") })
expect(identical(obj$x, 1) && grepl("RESULTS_ALLOW_STALE_INPUTS", w %||% "") && grepl("panel.rds", w %||% ""),
       "changed input + override: loads, warns naming the input")

unlink(inp)
expect(!errs(suppressWarnings(results_require(art, cfg, res))), "missing input + override: loads")

res2 <- res; res2$thresholds$m <- 0.9
expect(errs(suppressWarnings(results_require(art, cfg, res2))), "decision mismatch still stops under the override")

Sys.setenv(RESULTS_ALLOW_STALE_INPUTS = "false")
expect(errs(results_require(art, cfg, res)), "RESULTS_ALLOW_STALE_INPUTS=false: stops")
Sys.unsetenv("RESULTS_ALLOW_STALE_INPUTS")
unlink(td, recursive = TRUE)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
