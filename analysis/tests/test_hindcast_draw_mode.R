# test_hindcast_draw_mode.R -- the hindcasts' posterior-draw mode and replicate
# count (utils/make_hindcast_configs.R: hindcast_param_draw_mode(),
# hindcast_n_replicates()).
#
# HARD gate, quick-safe: no fit, no data, no simulation output.
#   [1] defaults: "" and 0 leave a patched config unchanged (patch_hindcast_config()'s
#       random at the template's count, the behaviour before the keys);
#   [2] joint: written with a uniform handoff and max_draws = the smoothing's
#       replicates; refused with a shuffled handoff, a max_draws that differs, or
#       an unknown replicate count; random is written as is; unknown modes refused;
#   [3] replicates: a positive count replaces the template's, negatives refused.
#
# Run: Rscript analysis/tests/test_hindcast_draw_mode.R

suppressPackageStartupMessages({ library(optparse) })
source("analysis/utils/base.R")
source("analysis/utils/make_hindcast_configs.R")   # CLI guarded (sys.nframe() == 0L)

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
cat("# hindcast posterior-draw mode and replicates\n")
err <- function(e) inherits(try(e, silent = TRUE), "try-error")

cfg0 <- list(simulator = list(n_replicates = 10000L, max_draws = 100L,
                              forecast = list(param_draw_mode = "random", handoff_weighting = "uniform")))

cat("# [1] defaults\n")
expect(identical(hindcast_param_draw_mode(cfg0, ""), cfg0) && identical(hindcast_param_draw_mode(cfg0, NULL), cfg0),
       "mode '' / NULL: config unchanged")
expect(identical(hindcast_n_replicates(cfg0, 0L), cfg0) && identical(hindcast_n_replicates(cfg0, NULL), cfg0),
       "n_replicates 0 / NULL: config unchanged")

cat("# [2] joint\n")
j <- hindcast_param_draw_mode(cfg0, "joint", n_smoothing = 100L)
expect(identical(j$simulator$forecast$param_draw_mode, "joint"), "joint written with uniform handoff, max_draws = nrep")
sh <- cfg0; sh$simulator$forecast$handoff_weighting <- "shuffled"
expect(err(hindcast_param_draw_mode(sh, "joint", 100L)), "joint refused with a shuffled handoff")
expect(err(hindcast_param_draw_mode(cfg0, "joint", 200L)), "joint refused when max_draws != smoothing replicates")
expect(err(hindcast_param_draw_mode(cfg0, "joint", NA_integer_)) && err(hindcast_param_draw_mode(cfg0, "joint", NULL)),
       "joint refused without the smoothing's replicate count")
expect(identical(hindcast_param_draw_mode(sh, "random")$simulator$forecast$param_draw_mode, "random"),
       "random written as is (any handoff)")
expect(err(hindcast_param_draw_mode(cfg0, "cycle")), "unknown mode refused")

cat("# [3] replicates\n")
expect(hindcast_n_replicates(cfg0, 5000L)$simulator$n_replicates == 5000L, "5000 replaces the template's 10000")
expect(err(hindcast_n_replicates(cfg0, -1L)), "a negative count refused")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
