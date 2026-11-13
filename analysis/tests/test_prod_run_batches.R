# test_prod_run_batches.R — the prod_run fit batches and kl class sweep
# (integration).
#
# make_batch_configs.R's prod_run entries, generated into a tempdir:
#   [1] occ / dur / mon: 115 configs each = ONI lags -11..11 x the 5 groups of
#       enso_grouping_prod_run, one per (lag, group) cell; brms_cb, ADM1-ADM2,
#       episodes.key gid_season; the occurrence panel spec = deploy.episodes
#   [2] duration carries seasonal_term harmonic2 and prior_mso_band
#       normal(0, 1.5) on every config
#   [3] fit names are distinct across lags and groups on every track
#   [4] the kl class sweep: 8 cells (ns4 / ns5 x K 2-5), kernel.classes = K,
#       ENSO basis ns df 4 / 5, outdir generated_data/prod_run_kl_sweep_classes/
#   [5] deferral: a fit batch whose registry is absent is deferred
#       (generate_batch() stops with the reason); the sweep never is
#   [6] prod_run.yml agrees with the registry: fits dirs and kl_sweep.batch
#       are the registered batch names, grouping = the batches' registry,
#       fits.groups = the registry's groups
#
# Pure config generation: no fit, no data. HARD, quick.
# Run: Rscript analysis/tests/test_prod_run_batches.R

suppressPackageStartupMessages({ library(yaml) })
suppressMessages(source("analysis/utils/make_batch_configs.R"))
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

REG    <- "enso_grouping_prod_run"
TRACKS <- c(occ = "occurrence", dur = "duration", mon = "monthly")
BATCH  <- setNames(sprintf("prod_run_%s_grouped_enso_singlelag", names(TRACKS)), names(TRACKS))
KL     <- "prod_run_kl_sweep_classes"
LAGS   <- -11:11
DEP    <- yaml::read_yaml("analysis/configs/deploy/prod_run.yml")$deploy

out <- tempfile("prod_run_batches_"); dir.create(out)
get_config_dir <- function() out
gen <- function(b) invisible(capture.output(generate_batch(b)))
cfgs_of <- function(b) sort(setdiff(list.files(file.path(out, b), full.names = TRUE),
                                    file.path(out, b, "manifest.yml")))
opt_of <- function(f) suppressMessages(make_options_from_config(option_list = list(), cli_argv = c("--config", f)))

reg <- load_grouping(REG)
groups <- sort(unique(as.character(reg$group)))
expect(length(groups) == 5L, sprintf("%s has 5 groups (%s)", REG, paste(groups, collapse = ", ")))

# --- [1] the three fit batches ------------------------------------------------
cat("\n[1] prod_run fit batches: 23 lags x 5 groups\n")
cfg <- list()
for (tr in names(TRACKS)) {
  b <- BATCH[[tr]]
  expect(identical(BATCH_REGISTRY[[b]]$grouping, REG), sprintf("%s grouped by %s", b, REG))
  gen(b)
  fs <- cfgs_of(b)
  ys <- lapply(fs, yaml::read_yaml)
  cfg[[tr]] <- ys
  blk <- lapply(ys, `[[`, TRACKS[[tr]])
  cells <- data.frame(lag = vapply(blk, function(x) as.integer(x$dlnm$defaults$single_lag), 1L),
                      group = vapply(ys, function(x) x$grouping$group, ""))
  expect(length(fs) == 115L, sprintf("%s: 115 configs (got %d)", tr, length(fs)))
  expect(setequal(cells$lag, LAGS) && setequal(cells$group, groups) &&
           !anyDuplicated(cells) && all(table(cells$group) == length(LAGS)),
         sprintf("%s: one config per (lag -11..11, group) cell", tr))
  expect(all(vapply(ys, function(x) identical(x$grouping$source, REG) && isTRUE(x$grouping$enable), TRUE)),
         sprintf("%s: grouping.source = %s on every config", tr, REG))
  expect(all(vapply(blk, function(x) identical(x$backend, "brms_cb") && identical(x$admin_level, "ADM1-ADM2"), TRUE)),
         sprintf("%s: brms_cb on ADM1-ADM2", tr))
  expect(all(vapply(ys, function(x) identical(x$episodes$key, "gid_season"), TRUE)),
         sprintf("%s: episodes.key gid_season", tr))
  expect(all(vapply(blk, function(x) identical(x$dlnm$defaults$argvar$fun, "ns") &&
                      identical(as.integer(x$dlnm$defaults$argvar$df), 4L), TRUE)),
         sprintf("%s: ENSO basis ns df 4", tr))
}
ep_dep <- DEP$episodes; ep_dep$national_zero_threshold <- as.integer(ep_dep$national_zero_threshold)
expect(all(vapply(cfg$occ, function(x) identical(x$episodes[names(ep_dep)], ep_dep) &&
                    setequal(names(x$episodes), names(ep_dep)), TRUE)),
       "occurrence episodes = deploy.episodes (the occurrence panel spec)")

# --- [2] duration --------------------------------------------------------------
cat("\n[2] duration: harmonic2 cessation seasonality, months-since-onset band prior\n")
expect(all(vapply(cfg$dur, function(x) identical(x$duration$seasonal_term, "harmonic2"), TRUE)),
       "duration.seasonal_term harmonic2 on every config")
expect(all(vapply(cfg$dur, function(x) identical(x$duration$brms_cb$prior_mso_band, "normal(0, 1.5)"), TRUE)),
       "duration.brms_cb.prior_mso_band normal(0, 1.5) on every config")

# --- [3] fit names ---------------------------------------------------------------
cat("\n[3] fit names distinct across lags and groups\n")
fit_name <- function(tr, f) {
  invisible(capture.output(o <- opt_of(f)))
  basename(switch(tr, occ = makeOccRunFile(o, "fit"), dur = makeDurationRunFile(o, "dur_fit"),
                  mon = makeMonthlyRunFile(o, "mon_fit")))
}
for (tr in names(TRACKS)) {
  fs <- cfgs_of(BATCH[[tr]])
  # every lag of one group, and every group at one lag
  pick <- fs[grepl("_grp-ensoeast\\.yml$", fs) | grepl("_sl-3_grp-", fs)]
  nm <- vapply(pick, function(f) fit_name(tr, f), "")
  expect(length(pick) == 27L && !anyDuplicated(nm),
         sprintf("%s: %d distinct fit names (23 lags of East + 4 other groups at lag 3)", tr, length(unique(nm))))
}

# --- [4] kl class sweep ----------------------------------------------------------
cat("\n[4] the kl class sweep\n")
gen(KL)
kf <- cfgs_of(KL)
ky <- lapply(kf, yaml::read_yaml)
cell <- sub("^.*_basis-(ns[45])_K-([0-9])\\.yml$", "\\1_K\\2", basename(kf))
expect(length(kf) == 8L && setequal(cell, as.vector(outer(c("ns4", "ns5"), 2:5, paste, sep = "_K"))),
       sprintf("8 cells ns4 / ns5 x K 2-5 (%s)", paste(cell, collapse = ", ")))
ok_cell <- mapply(function(y, cl) {
  k <- as.integer(sub("^ns[45]_K", "", cl)); df <- as.integer(substr(cl, 3, 3))
  identical(as.integer(y$occurrence$kernel$classes), k) &&
    identical(as.integer(y$occurrence$dlnm$defaults$argvar$df), df) &&
    identical(y$occurrence$backend, "brms_kl")
}, ky, cell)
expect(all(ok_cell), "each cell: brms_kl, kernel.classes = K, ENSO basis ns df = the basis")
expect(all(vapply(ky, function(y) identical(y$outdir, "generated_data/prod_run_kl_sweep_classes/"), TRUE)),
       "outdir generated_data/prod_run_kl_sweep_classes/ on every cell")
expect(is.null(ky[[1]]$kl_sweep), "the synthetic kl_sweep.* grid keys are removed")

# --- [5] deferral ------------------------------------------------------------------
cat("\n[5] deferral while a fit batch's registry is absent\n")
spec_absent <- BATCH_REGISTRY[[BATCH[["dur"]]]]; spec_absent$grouping <- "zz_registry_not_written"
expect(grepl("not written yet", .batch_deferred(spec_absent) %||% ""),
       "a fit batch whose registry is absent is deferred")
BATCH_REGISTRY$zz_deferred_batch <- spec_absent
e <- tryCatch({ generate_batch("zz_deferred_batch"); NULL }, error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("is deferred", e), "generate_batch() stops with the deferral reason")
BATCH_REGISTRY$zz_deferred_batch <- NULL
expect(all(vapply(BATCH, function(b) is.null(.batch_deferred(BATCH_REGISTRY[[b]])), TRUE)),
       "the prod_run fit batches are not deferred (their registry exists)")
expect(is.null(.batch_deferred(BATCH_REGISTRY[[KL]])), "the kl sweep is never deferred (it derives the registry)")

# --- [6] prod_run.yml --------------------------------------------------------------
cat("\n[6] prod_run.yml agrees with the registry\n")
expect(identical(c(DEP$fits$occ_dir, DEP$fits$dur_dir, DEP$fits$mon_dir), unname(BATCH)),
       "deploy.fits.{occ,dur,mon}_dir = the registered batch names")
expect(identical(DEP$kl_sweep$batch, KL) && KL %in% names(BATCH_REGISTRY),
       "deploy.kl_sweep.batch = prod_run_kl_sweep_classes (registered)")
expect(identical(DEP$grouping, REG), "deploy.grouping = the fit batches' registry")
expect(setequal(unlist(DEP$fits$groups), sub("^enso_", "", groups)),
       "deploy.fits.groups = the registry's groups")
expect(identical(DEP$kl_sweep$cell, "ns4_K5") && "ns4_K5" %in% cell,
       "the pinned kl cell is one of the sweep's cells")

unlink(out, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
