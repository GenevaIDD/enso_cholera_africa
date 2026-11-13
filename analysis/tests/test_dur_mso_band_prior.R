# test_dur_mso_band_prior.R — duration prior_mso_band.
#
# Hard-fail gate. The brms-block key prior_mso_band gives the episode-age band
# effects their own prior (one class = "b", coef = "mso_band<level>" row per
# occupied non-reference band). Checked here:
#   1. unset: build_duration_priors() output is identical with and without dat;
#   2. set: one row per occupied non-reference band, brms-named (4-5 -> 4M5,
#      9+ -> 9P), with the set prior; the class-b default row is kept;
#   3. an unoccupied band gets no row (brms rejects a prior on a missing coef);
#   4. set without dat errors;
#   5. makeDurationRunFile(): no token when unset, "_mbp1p5" for normal(0, 1.5).

suppressPackageStartupMessages(library(brms))
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/formula.R")
source("analysis/utils/priors.R")
source("analysis/utils/paths.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

opt <- get_default_options()
opt$duration$backend <- "brms_cb"
opt$duration$brms_cb <- opt$duration$brms
parts <- list(fixed_terms = c("1", "enso_cb1"), bh_terms = "mso_band", seas_terms = character(0),
              re_terms = "(1 | country)", off_terms = character(0))
lv  <- c("1", "2", "3", "4-5", "6-8", "9+")
dat <- data.frame(mso_band = factor(rep(lv, 2), levels = lv))

cat("\n[1] unset: unchanged\n")
p0 <- build_duration_priors(opt, parts, opt$duration$brms_cb)
p1 <- build_duration_priors(opt, parts, opt$duration$brms_cb, dat = dat)
expect(identical(p0, p1), "prior_mso_band unset: dat is not read, output identical")
expect(!any(nzchar(as.character(p0$coef))), "prior_mso_band unset: no coef-level rows")

cat("\n[2] set: one row per non-reference band\n")
cfg <- opt$duration$brms_cb; cfg$prior_mso_band <- "normal(0, 1.5)"
p2 <- build_duration_priors(opt, parts, cfg, dat = dat)
rows <- p2[nzchar(as.character(p2$coef)), ]
expect(setequal(rows$coef, c("mso_band2", "mso_band3", "mso_band4M5", "mso_band6M8", "mso_band9P")),
       "band rows named as brms names the coefficients")
expect(all(rows$prior == "normal(0, 1.5)") && all(rows$class == "b"), "band rows carry the set prior, class b")
expect(any(p2$class == "b" & !nzchar(as.character(p2$coef)) & p2$prior == "normal(0, 0.5)"),
       "class-b default row kept for the other coefficients")

cat("\n[3] unoccupied band: no row\n")
dat3 <- data.frame(mso_band = factor(c("1", "2", "3", "4-5"), levels = lv))
r3 <- build_duration_priors(opt, parts, cfg, dat = dat3)
expect(setequal(r3$coef[nzchar(as.character(r3$coef))], c("mso_band2", "mso_band3", "mso_band4M5")),
       "only occupied non-reference bands get a row")

cat("\n[4] set without dat: error\n")
expect(inherits(tryCatch(build_duration_priors(opt, parts, cfg), error = function(e) e), "error"),
       "prior_mso_band set and dat missing -> error")

cat("\n[5] fit-name token\n")
f0 <- basename(makeDurationRunFile(opt, kind = "dur_fit"))
opt5 <- opt; opt5$duration$brms_cb$prior_mso_band <- "normal(0, 1.5)"
f5 <- basename(makeDurationRunFile(opt5, kind = "dur_fit"))
expect(!grepl("_mbp", f0), "no token when unset")
expect(grepl("_mbp1p5_", f5) && identical(sub("_mbp1p5", "", f5), f0), "token _mbp1p5 is the only difference")

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
message("PASS  test_dur_mso_band_prior")
