# test_loo_compare_paths.R -- LOO ranking artefacts live under
# generated_data/loo_compare/<batch>/, and every reader goes through ONE locator.
#
# HARD gate (see analysis/tests/run_all.R). No fit, no data, no MCMC.
#   [1] makeLooCompareRunFile(): dir/batch-name equivalence, five distinct kinds;
#       the lag-sidecar name never matches the ranking regex.
#   [2] loo_compare_ranking_files(): new dir beats legacy; legacy fallback;
#       mtime tiebreak at equal date; group = "" / NULL semantics.
#   [3] Static scan: no non-comment line outside utils/loo_lags.R and
#       utils/paths.R spells the `loo_compare_lags_` filename pattern, i.e. no
#       reader re-implements the glob.
#
# Run: Rscript analysis/tests/test_loo_compare_paths.R

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/loo_lags.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# loo_compare paths + locator\n")

# [1] path helper ----------------------------------------------------------
cat("\n[1] makeLooCompareRunFile\n")
d_path <- makeLooCompareRunFile(list(outdir = "gd"), "loo_compare_dir", batch = "analysis/configs/x/")
d_name <- makeLooCompareRunFile(list(outdir = "gd"), "loo_compare_dir", batch = "x")
expect(identical(d_path, d_name), "batch dir path and bare batch name resolve to the same dir")
expect(identical(unname(d_name), file.path("gd", "loo_compare", "x")), "dir is <outdir>/loo_compare/<batch>")
kinds <- c("loo_compare_dir", "loo_compare_ranking", "loo_compare_png", "loo_compare_coef",
           "loo_compare_lag_sidecar")
paths <- vapply(kinds, function(k) makeLooCompareRunFile(list(outdir = "gd"), k, batch = "x",
                                                         group = "enso_east", date = as.Date("2026-01-02")),
                character(1))
expect(!anyDuplicated(paths), "the five kinds give pairwise distinct paths")
expect(grepl("lag_loo_sidecar_enso_east_20260102\\.rds$", paths[["loo_compare_lag_sidecar"]]),
       "sidecar filename = lag_loo_sidecar_<group>_<YYYYMMDD>.rds")
side_names <- basename(c(paths[["loo_compare_lag_sidecar"]],
                         makeLooCompareRunFile(list(), "loo_compare_lag_sidecar", batch = "x",
                                               date = as.Date("2026-01-02"))))
expect(!any(grepl(LOO_COMPARE_RANKING_RE, side_names)) && all(grepl(LAG_LOO_SIDECAR_RE, side_names)),
       "sidecar filenames (grouped / ungrouped) never match the ranking regex")
expect(grepl(LOO_COMPARE_RANKING_RE, basename(paths[["loo_compare_ranking"]])) &&
         !grepl(LAG_LOO_SIDECAR_RE, basename(paths[["loo_compare_ranking"]])),
       "and the ranking filename never matches the sidecar regex")
expect(grepl("loo_compare_lags_enso_east_20260102\\.rds$", paths[["loo_compare_ranking"]]),
       "ranking filename = loo_compare_lags_<group>_<YYYYMMDD>.rds")
expect(grepl("/loo_compare_lags\\.rds$|loo_compare_lags_20260102\\.rds$",
             makeLooCompareRunFile(list(), "loo_compare_ranking", batch = "x", date = as.Date("2026-01-02"))),
       "ungrouped ranking filename omits the group token")

# [2] locator ---------------------------------------------------------------
cat("\n[2] loo_compare_ranking_files\n")
tmp    <- tempfile("loo_paths_"); dir.create(tmp)
opt    <- list(outdir = file.path(tmp, "gd"))
batch  <- file.path(tmp, "configs", "b_occ_grouped_enso_singlelag")   # legacy dir = batch dir
newdir <- makeLooCompareRunFile(opt, "loo_compare_dir", batch = batch)
dir.create(batch, recursive = TRUE); dir.create(newdir, recursive = TRUE)
touch <- function(dir, name, age_days = 0) {
  f <- file.path(dir, name); saveRDS(list(ranking = data.frame(single_lag = 1L)), f)
  Sys.setFileTime(f, Sys.time() - 86400 * age_days); f
}

f_leg <- touch(batch, "loo_compare_lags_enso_east_20260901.rds")
lf <- loo_compare_ranking_files(batch, "enso_east", opt = opt)
expect(nrow(lf) == 1L && lf$legacy && identical(lf$file, f_leg),
       "new dir empty -> legacy fallback is read (flagged legacy = TRUE)")

f_new <- touch(newdir, "loo_compare_lags_enso_east_20260101.rds")
lf <- loo_compare_ranking_files(batch, "enso_east", opt = opt)
expect(nrow(lf) == 1L && !lf$legacy && identical(lf$file, f_new),
       "new dir beats legacy even when the legacy file has a newer date")

f_new2 <- touch(newdir, "loo_compare_lags_enso_east_20260301.rds", age_days = 30)
invisible(touch(newdir, "loo_compare_lags_enso_west_20260301.rds"))
invisible(touch(newdir, "loo_compare_lags_20260301.rds"))            # ungrouped
expect(identical(latest_ranking_file(batch, "enso_east", opt = opt), f_new2),
       "latest = filename date descending (older mtime, newer date wins)")
all_rows <- loo_compare_ranking_files(batch, opt = opt)
expect(nrow(all_rows) == 3L, "group = NULL, latest_only: one row per group (east, west, ungrouped)")
invisible(touch(newdir, basename(makeLooCompareRunFile(opt, "loo_compare_lag_sidecar", batch = batch,
                                                       group = "enso_east", date = as.Date("2026-12-31")))))
expect(identical(latest_ranking_file(batch, "enso_east", opt = opt), f_new2) &&
         nrow(loo_compare_ranking_files(batch, opt = opt, latest_only = FALSE)) == 4L,
       "a (newer-dated) lag sidecar in the dir is invisible to the ranking locator")
expect(nrow(loo_compare_ranking_files(batch, "", opt = opt)) == 1L &&
         is.na(loo_compare_ranking_files(batch, "", opt = opt)$group),
       'group = "" matches the ungrouped file only')
expect(nrow(loo_compare_ranking_files(batch, "enso_east", opt = opt, latest_only = FALSE)) == 2L,
       "latest_only = FALSE returns every dated file for the group")
expect(is.null(latest_ranking_file(batch, "enso_nowhere", opt = opt)), "no match -> NULL")
expect(is.null(read_latest_ranking(batch, "enso_nowhere", opt = opt, required = FALSE)),
       "read_latest_ranking(required = FALSE) -> NULL")
expect(inherits(try(read_latest_ranking(batch, "enso_nowhere", opt = opt), silent = TRUE), "try-error"),
       "read_latest_ranking(required = TRUE) stops")

# [3] static scan -------------------------------------------------------------
cat("\n[3] no reader re-implements the filename glob\n")
allowed <- c("analysis/utils/loo_lags.R", "analysis/utils/paths.R",
             "analysis/tests/test_loo_compare_paths.R", "analysis/tests/test_loo_lags_helpers.R")
files <- c(list.files("analysis", pattern = "\\.(R|sh)$", full.names = TRUE),
           list.files("analysis/utils", pattern = "\\.R$", full.names = TRUE),
           list.files("analysis/scratchpads", pattern = "\\.R$", full.names = TRUE),
           list.files("analysis/bash", pattern = "\\.sh$", full.names = TRUE, recursive = TRUE),
           list.files("notebooks", pattern = "\\.Rmd$", full.names = TRUE))
files <- setdiff(unique(files[!nzchar(Sys.readlink(files))]), allowed)   # skip bash/ back-compat symlinks
# Code lines only: drop `#` comments (R/sh); for .Rmd keep just the fenced
# chunk bodies, since markdown prose may legitimately name the file.
code_lines <- function(f) {
  txt <- readLines(f, warn = FALSE)
  if (grepl("\\.Rmd$", f)) {
    in_chunk <- cumsum(grepl("^```\\{", txt)) > cumsum(grepl("^```\\s*$", txt)) &
                !grepl("^```", txt)
    txt <- txt[in_chunk]
  }
  txt[!grepl("^\\s*#", txt)]
}
offenders <- Filter(function(f) any(grepl("loo_compare_lags_", code_lines(f), fixed = TRUE)), files)
expect(!length(offenders),
       if (length(offenders)) paste("loo_compare_lags_ pattern re-implemented in:", paste(offenders, collapse = ", "))
       else "no non-comment `loo_compare_lags_` outside the locator + path helper")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
