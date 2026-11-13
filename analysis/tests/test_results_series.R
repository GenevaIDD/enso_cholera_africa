# test_results_series.R -- results_select_series() (results_core.R): the
# results.series switch between the latent and the detected simulated series.
# Run: Rscript analysis/tests/test_results_series.R   (project root)
suppressPackageStartupMessages(library(dplyr))
source("analysis/utils/results_core.R")
source("analysis/tests/helpers.R")
PASS <- 0L; FAIL <- 0L

ev  <- tibble::tibble(rep_id = 1:4, cases = c(5, 6, 7, 8), detected = c(1L, 0L, 1L, 0L))
grd <- tibble::tibble(Y = c(1L, 1L, 0L), cases_month = c(10, 20, NA), detected = c(1L, 0L, 0L))
expect(identical(results_select_series(ev, "latent", "x"), ev), "latent: events unchanged")
expect(identical(results_select_series(ev, NULL, "x"), ev), "NULL series = latent")
expect(identical(results_select_series(ev, "detected", "x")$rep_id, c(1L, 3L)),
       "detected: event rows filtered to detected episodes")
g <- results_select_series(grd, "detected", "x", kind = "grid")
expect(nrow(g) == 3L && identical(g$Y, c(1L, 0L, 0L)) && identical(g$cases_month, c(10, NA, NA)),
       "detected: grid kept complete, undetected Y/cases zeroed")
expect(inherits(try(results_select_series(dplyr::select(ev, -detected), "detected", "run X"),
                    silent = TRUE), "try-error"),
       "detected on a run without the column -> error (no silent fall-back)")
expect(inherits(try(results_select_series(ev, "observed", "x"), silent = TRUE), "try-error"),
       "unknown series -> error")
cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
