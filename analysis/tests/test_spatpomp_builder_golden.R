# test_spatpomp_builder_golden.R -- golden characterization of the spatPomp
# builder build_spatpomp_tu_c().
#
# HARD gate (see analysis/tests/run_all.R). Every builder variant below must
# produce the same object as the golden after a refactor: its builder attributes
# (frame, filter months, warm-up
# state, draw subset, ENSO / cross-basis names, covariate maps, ...), unit
# names, times, data, params, covariate table and the generated C source
# (pomp's library hash normalised). Builds only -- no particle filter runs.
#
# Fixture: the 2-unit x 6-month one of test_spatpomp_dunit_detect.R.
#
# Season golden (simulator.onset_eligibility): a second file,
# spatpomp_builder_season_golden.rds, holds the duration-path build of the same
# fixture (unit min_seas 1 and 4) under onset_eligibility inactive and season,
# so the emitted season gate and the seas_k covariate are pinned too. It has its
# own capture flag, so capturing it never rewrites the main golden.
#
#   Rscript analysis/tests/test_spatpomp_builder_golden.R --capture          # once, pre-refactor
#   Rscript analysis/tests/test_spatpomp_builder_golden.R --capture-season   # season golden only
#   Rscript analysis/tests/test_spatpomp_builder_golden.R                    # compare both

suppressPackageStartupMessages({
  library(optparse); library(tibble); library(dplyr); library(lubridate); library(tidyr)
})
cli <- parse_args(OptionParser(option_list = list(
  make_option("--capture", action = "store_true", default = FALSE,
              help = "write the golden file instead of comparing"),
  make_option("--capture-season", action = "store_true", default = FALSE,
              dest = "capture_season",
              help = "write the season golden file instead of comparing it"))))

source("analysis/utils.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")
source("analysis/utils/spatpomp_core_c.R")   # also sources spatpomp_core.R
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

GOLDEN_PATH <- "analysis/tests/golden/spatpomp_builder_golden.rds"
cat("# spatPomp builder golden (mode:", if (cli$capture) "capture" else "compare", ")\n")

# Fixture: reuse test_spatpomp_dunit_detect.R's make_fixture() + synthetic draws.
src <- readLines("analysis/tests/test_spatpomp_dunit_detect.R")
i1 <- grep("^make_fixture <- function", src); i2 <- grep("^build_fix <- function", src) - 2L
eval(parse(text = src[i1:i2]))

snapshot <- function(o) {
  extra <- setdiff(names(attributes(o)), c(methods::slotNames(o), "class"))
  srcs <- vapply(o@solibs, function(s)   # mask pomp's library hash + generation time
    gsub("/\\* Time: [^*]*\\*/", "/* Time: X */", gsub("pomp_[0-9a-f]+", "pomp_X", s$src)), "")
  list(attrs = lapply(setNames(extra, extra), function(a) attr(o, a)),
       unit_names = o@unit_names, times = o@times, data = o@data,
       params = o@params, covar = o@covar, solib_src = srcs)
}
variants <- list(
  c_bern        = list(builder = build_spatpomp_tu_c, measurement = "occurrence_bernoulli", fixed = 1L),
  c_state       = list(builder = build_spatpomp_tu_c, measurement = "occurrence_bernoulli_state", fixed = 1L))
if (exists(".SP_DETECT_MEAS")) for (m in intersect(.SP_DETECT_MEAS, c("occurrence_bernoulli_state_detect"))) {
  variants[[paste0("c_", m)]] <- list(builder = build_spatpomp_tu_c, measurement = m, fixed = 1L)
}
out <- list()
for (nm in names(variants)) {
  v <- variants[[nm]]
  set.seed(20261006)
  o <- suppressMessages(v$builder(fix$frame, fix$panel, draws_occ, fix$window,
                                  fixed_draw_idx = v$fixed, measurement = v$measurement,
                                  obs_data = fix$obs))
  out[[nm]] <- snapshot(o)
}
expect(length(out) >= 2L, sprintf("built %d builder variants", length(out)))

if (cli$capture) {
  saveRDS(out, GOLDEN_PATH)
  cat("  captured", length(out), "variants ->", GOLDEN_PATH, "\n")
} else {
  expect(file.exists(GOLDEN_PATH), paste("golden file exists:", GOLDEN_PATH))
  if (file.exists(GOLDEN_PATH)) {
    gold <- readRDS(GOLDEN_PATH)
    expect(identical(sort(names(gold)), sort(names(out))), "same builder variants as the golden")
    for (nm in intersect(names(gold), names(out))) for (k in names(gold[[nm]])) {
      ok <- identical(gold[[nm]][[k]], out[[nm]][[k]])
      expect(ok, sprintf("identical to golden: %s / %s", nm, k))
      if (!ok) print(all.equal(gold[[nm]][[k]], out[[nm]][[k]]))
    }
  }
}
# --- season golden (onset_eligibility, duration path) ------------------------
SEASON_GOLDEN_PATH <- "analysis/tests/golden/spatpomp_builder_season_golden.rds"
fr_s <- fix$frame; fr_s$min_seas <- c(1L, 4L)
pa_s <- fix$panel; pa_s$min_seas <- fr_s$min_seas[match(pa_s$gid, fr_s$gid)]
draws_dur_s <- make_draws_bundle_synthetic(fr_s, S = 1L, fe_overrides = c("(Intercept)" = -1.0),
                                           active_terms = "(Intercept)")
out_s <- list()
for (oe in c("inactive", "season")) {
  set.seed(20261107)
  o <- suppressMessages(build_spatpomp_tu_c(fr_s, pa_s, draws_occ, fix$window, fixed_draw_idx = 1L,
                                            measurement = "occurrence_bernoulli_state",
                                            obs_data = fix$obs, draws_dur = draws_dur_s,
                                            onset_eligibility = oe))
  out_s[[paste0("c_state_dur_", oe)]] <- snapshot(o)
}
if (isTRUE(cli$capture_season)) {
  saveRDS(out_s, SEASON_GOLDEN_PATH)
  cat("  captured", length(out_s), "season variants ->", SEASON_GOLDEN_PATH, "\n")
} else {
  expect(file.exists(SEASON_GOLDEN_PATH), paste("season golden file exists:", SEASON_GOLDEN_PATH))
  if (file.exists(SEASON_GOLDEN_PATH)) {
    gold_s <- readRDS(SEASON_GOLDEN_PATH)
    expect(identical(sort(names(gold_s)), sort(names(out_s))), "same season variants as the golden")
    for (nm in intersect(names(gold_s), names(out_s))) for (k in names(gold_s[[nm]])) {
      ok <- identical(gold_s[[nm]][[k]], out_s[[nm]][[k]])
      expect(ok, sprintf("identical to season golden: %s / %s", nm, k))
      if (!ok) print(all.equal(gold_s[[nm]][[k]], out_s[[nm]][[k]]))
    }
  }
}
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
