# test_basename_length.R — artefact and config names fit NAME_MAX (255 bytes).
#
# A file name over 255 bytes cannot be created, and the failure can come late
# (writing cb_meta after the fit is saved, or inside 05_fit_postprocess.R's
# tryCatch). paths.R::assertBasenameLength() stops in every name builder,
# and checkRunFileNames() builds every kind of a config's track up front (the
# fit workers and generate_batch() call it), so a long name fails before a fit.
#
#   [1] assertBasenameLength(): passes at the limit, stops one byte over with the
#       builder, the length and the longest tokens; checks every path component
#   [2] every name every builder emits for every config in analysis/configs/
#       (each track builder x kind, the data-file builders, the sim / smooth /
#       spatpomp, results, deploy, loo_compare and kl_sweep builders, and the
#       sidecars derived from them by sub / paste) fits the limit
#   [3] checkRunFileNames() passes on every fit config and stops on a config
#       whose longest kind is over the limit
#   [4] every file and directory name under analysis/configs/ fits the limit
#
# Pure path algebra: no fit, no data. HARD, quick.
# Run: Rscript analysis/tests/test_basename_length.R
#      Rscript analysis/tests/test_basename_length.R --dump <file.tsv>
#        (writes every enumerated name; the golden for a naming refactor)

suppressPackageStartupMessages({ library(yaml) })
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")   # currentPopDigestToken(), fitPopStampPath(), prepPolicyStampPath()
source("analysis/utils/enso.R")         # ENSO_ANALOG_SETS
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
dump <- if ("--dump" %in% args) args[[which(args == "--dump") + 1L]] else NULL
LIMIT <- if (exists("BASENAME_MAX_BYTES")) BASENAME_MAX_BYTES else 255L

CONFIG_ROOT <- "analysis/configs"
TRACK_BUILDERS <- list(occurrence = "makeOccRunFile", duration = "makeDurationRunFile",
                       monthly = "makeMonthlyRunFile", size = "makeSizeRunFile")
DATA_BUILDERS <- c("makeSurvDatFile", "makeDurationDatFile",
                   "makeMonthlyDatFile", "makeSizeDatFile", "makeSimPanelFile",
                   "makeSimPanelProjectionFile", "makeFootprintMembersFile")
kinds_of <- function(fn) eval(formals(get(fn))$kind)

# --- Enumerate every name -----------------------------------------------------
rows <- list()
add <- function(src, builder, kind, path)
  rows[[length(rows) + 1L]] <<- data.frame(source = src, builder = builder, kind = kind,
                                           path = as.character(path), stringsAsFactors = FALSE)
# A builder call; an error is recorded as a row so that it fails [2] below.
try_add <- function(src, builder, kind, expr) {
  p <- tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(p, "err")) add(src, builder, kind, paste0("<error> ", unclass(p)))
  else for (q in p) add(src, builder, kind, q)
}

yml <- sort(list.files(CONFIG_ROOT, "\\.ya?ml$", recursive = TRUE, full.names = TRUE))
yml <- yml[!basename(yml) %in% c("manifest.yml", "LAG_SELECTION.yml")]
is_results <- startsWith(yml, file.path(CONFIG_ROOT, "results/"))
is_deploy  <- startsWith(yml, file.path(CONFIG_ROOT, "deploy/"))
opt_errors <- character()
batch_groups <- list()
fit_opts <- list()

for (f in yml[!is_results & !is_deploy]) {
  raw <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
  if (!is.list(raw)) next
  o <- tryCatch({
         out <- NULL
         invisible(capture.output(out <- suppressMessages(suppressWarnings(
           make_options_from_config(option_list = list(), cli_argv = c("--config", f))))))
         out
       }, error = function(e) e)
  # A file without a track or simulator block (the generators' *_horizons
  # parameter files) is not a pipeline config; one with either must resolve.
  if (inherits(o, "error")) {
    if (length(intersect(c(names(TRACK_BUILDERS), "simulator"), names(raw))))
      opt_errors[f] <- conditionMessage(o)
    next
  }
  tracks <- intersect(names(TRACK_BUILDERS), names(raw))
  if (!length(tracks) && is.null(raw[["simulator"]])) tracks <- "occurrence"
  for (tr in tracks) {
    fn <- TRACK_BUILDERS[[tr]]
    for (k in kinds_of(fn)) try_add(f, fn, k, get(fn)(opt = o, kind = k))
    fit_opts[[length(fit_opts) + 1L]] <- list(file = f, track = tr, opt = o)
  }
  for (fn in DATA_BUILDERS) try_add(f, fn, "", get(fn)(o))
  if (!is.null(raw[["simulator"]]))
    for (fn in c("makeSimRunFile", "makeSmoothRunFile", "makeSpatpompRunFile"))
      for (k in kinds_of(fn)) try_add(f, fn, k, get(fn)(o, kind = k))
  # Sidecars built outside paths.R by sub / paste on a builder's name.
  if ("occurrence" %in% tracks) {
    fit <- tryCatch(makeOccRunFile(o, "fit"), error = function(e) NULL)
    # 07b_simulator_v1.R / simulations.R / smoothing_core.R: cb_meta from the fit
    if (!is.null(fit)) add(f, "derived", "fit_to_cb_meta", sub("_fit_base-", "_cb_meta_base-", fit, fixed = TRUE))
    if (identical(o$occurrence$backend, "brms_kl")) {   # sidecars of the kl fits only
      pin  <- tryCatch(makeOccRunFile(o, "kl_pin"), error = function(e) NULL)
      part <- tryCatch(makeOccRunFile(o, "kl_partition"), error = function(e) NULL)
      ch   <- paste(seq_len(max(4L, as.integer(o$chains %||% 4L))), collapse = "-")
      # make_kl_grouping.R --write_fit: <fit>_ch<chains>.rds
      if (!is.null(fit))  add(f, "derived", "fit_ch", sub("\\.rds$", sprintf("_ch%s.rds", ch), fit))
      if (!is.null(pin))  add(f, "derived", "kl_pin_arr", sub("\\.rds$", "_arr99.rds", pin))   # dlnm_kernel.R
      if (!is.null(part)) add(f, "derived", "kl_partition_csv", sub("\\.rds$", ".csv", part))  # 04_fit_occ_brms_kl.R
    }
    surv <- tryCatch(as.character(makeSurvDatFile(o)), error = function(e) NULL)
    if (!is.null(surv)) {
      add(f, "derived", "surv_sidecars", c(unlist(occPanelSidecarPaths(surv)), occPanelLpMembersPath(surv),
                                           contagionCoveragePath(surv), prepPolicyStampPath(surv),
                                           sub("\\.rds$", ".min_seas_source.csv", surv)))
    }
  }
  for (tr in intersect(c("duration", "monthly"), tracks)) {
    fit <- tryCatch(makeRunFile(o, if (tr == "duration") "dur_fit" else "mon_fit"), error = function(e) NULL)
    if (!is.null(fit)) add(f, "derived", "fit_to_cb_meta", sub("_fit_base-", "_cb_meta_base-", fit, fixed = TRUE))
  }
  if ("size" %in% tracks) {
    fit <- tryCatch(makeSizeRunFile(o, "size_fit"), error = function(e) NULL)
    mrt <- tryCatch(makeSizeRunFile(o, "size_mean_rate"), error = function(e) NULL)
    if (!is.null(fit)) add(f, "derived", "size_fit_sidecars",
                           c(fitPopStampPath(fit), paste0(tools::file_path_sans_ext(fit), ".stan")))
    if (!is.null(mrt)) add(f, "derived", "mean_rate_csv",   # 05_fit_postprocess.R
                           c(sub("\\.rds$", ".csv", mrt), sub("\\.rds$", "_tier_summary.csv", mrt)))
  }
  bd <- dirname(f)
  if (!identical(bd, CONFIG_ROOT))
    batch_groups[[bd]] <- unique(c(batch_groups[[bd]], o$grouping$group %||% ""))
}

# Batch-keyed builders: one per batch dir, every group seen in it.
for (bd in names(batch_groups)) {
  for (k in kinds_of("makeLooCompareRunFile"))
    for (g in batch_groups[[bd]])
      try_add(bd, "makeLooCompareRunFile", k,
              makeLooCompareRunFile(list(), k, batch = bd, group = g, date = as.Date("2026-11-08")))
  for (k in kinds_of("makeKlSweepRunFile"))
    try_add(bd, "makeKlSweepRunFile", k, makeKlSweepRunFile(list(), k, batch = bd))
}

# Results configs: every kind, the longest schema token and every analog set.
schema_max <- "v99"
for (f in yml[is_results]) {
  for (k in kinds_of("makeResultsRunFile"))
    try_add(f, "makeResultsRunFile", k,
            makeResultsRunFile(list(), k, track = "occurrence", sets = ENSO_ANALOG_SETS,
                               schema = schema_max, config_path = f))
  try_add(f, "makeFigureDir", "", makeFigureDir(list(), track = "cb", config_path = f, subdir = "adm1"))
}

# Deploy YAMLs: the run root, and every run tag / batch name they state.
for (f in yml[is_deploy]) {
  for (k in kinds_of("makeDeployRunFile")) try_add(f, "makeDeployRunFile", k, makeDeployRunFile(list(), k, config_path = f))
  d <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
  dep <- d$deploy %||% d
  tags <- unlist(c(dep$run_tags, dep$fits), use.names = FALSE)
  tags <- tags[is.character(tags)]
  if (length(tags)) add(f, "deploy_tags", "", file.path("generated_data", basename(tags)))
}

names_df <- do.call(rbind, rows)
names_df$bytes <- nchar(basename(names_df$path), type = "bytes")
names_df$max_component <- vapply(strsplit(names_df$path, "/", fixed = TRUE),
                                 function(p) max(nchar(p, type = "bytes")), 0L)

if (!is.null(dump)) {
  write.table(names_df[, c("source", "builder", "kind", "path")], dump, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat(sprintf("wrote %d names (%d configs did not resolve) to %s\n",
              nrow(names_df), length(opt_errors), dump))
  quit(save = "no", status = 0L)
}

# --- [1] the guard ---------------------------------------------------------------
cat("\n[1] assertBasenameLength()\n")
at  <- paste0(strrep("a", LIMIT - 4L), ".rds")
over <- paste0("outbreaks_x_", strrep("b", LIMIT), "_ep-gs.rds")
expect(identical(assertBasenameLength(file.path("generated_data", at), "makeTest"),
                 file.path("generated_data", at)),
       sprintf("a %d-byte name passes and is returned unchanged", LIMIT))
e <- tryCatch({ assertBasenameLength(file.path("generated_data", over), "makeTest"); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("makeTest", e) &&
       grepl(sprintf("%d bytes", nchar(over, type = "bytes")), e) &&
       grepl(sprintf("%d-byte limit", LIMIT), e) && grepl(strrep("b", 20L), e),
       "one over the limit stops, naming the builder, the length, the limit and the longest token")
e2 <- tryCatch({ assertBasenameLength(file.path("generated_data", strrep("d", LIMIT + 1L), "x.rds"), "makeDir"); NULL },
               error = function(e) conditionMessage(e))
expect(!is.null(e2) && grepl("makeDir", e2), "a directory component over the limit stops too")
expect(identical(assertBasenameLength(c("a/b.rds", "c.qs"), "makeTest"), c("a/b.rds", "c.qs")),
       "vectorised: every element checked, input returned")
multi <- paste0(strrep("é", 128L), ".rds")   # 260 bytes, 132 characters
expect(inherits(tryCatch(assertBasenameLength(multi, "makeTest"), error = identity), "error"),
       "the limit is in bytes, not characters")
expect(identical(.capFileName(paste0(strrep("a", 300L), ".rds")),
                 assertBasenameLength(.capFileName(paste0(strrep("a", 300L), ".rds")), "makeTest")),
       ".capFileName() output (occurrence) always passes the guard")

# --- [2] every enumerated name ------------------------------------------------
cat("\n[2] every name the builders emit for analysis/configs/\n")
errs <- names_df[startsWith(names_df$path, "<error>"), ]
expect(nrow(errs) == 0L,
       sprintf("no builder errors (%d)%s", nrow(errs),
               if (nrow(errs)) paste0(": ", paste(head(unique(paste(errs$builder, errs$kind, errs$path)), 5L),
                                              collapse = " | ")) else ""))
ok_df <- names_df[!startsWith(names_df$path, "<error>"), ]
bad <- ok_df[ok_df$max_component > LIMIT, ]
expect(nrow(bad) == 0L,
       sprintf("%d names, every path component <= %d bytes%s", nrow(ok_df), LIMIT,
               if (nrow(bad)) paste0(": ", paste(head(basename(bad$path), 3L), collapse = " | ")) else ""))
expect(length(opt_errors) == 0L,
       sprintf("every fit / simulator config resolves%s",
               if (length(opt_errors)) paste0(": ", names(opt_errors)[[1L]], " -- ", opt_errors[[1L]]) else ""))
summ <- aggregate(bytes ~ builder + kind, ok_df, max)
summ <- summ[order(-summ$bytes), ]
cat("  longest per builder x kind (bytes, margin):\n")
for (i in seq_len(min(12L, nrow(summ))))
  cat(sprintf("    %-24s %-26s %3d  %3d\n", summ$builder[i], summ$kind[i], summ$bytes[i], LIMIT - summ$bytes[i]))

# --- [3] checkRunFileNames() -----------------------------------------------------
cat("\n[3] checkRunFileNames()\n")
chk <- vapply(fit_opts, function(x)
  tryCatch({ checkRunFileNames(x$opt, x$track); "" }, error = function(e) conditionMessage(e)), "")
expect(all(!nzchar(chk)),
       sprintf("passes on all %d fit configs%s", length(chk),
               if (any(nzchar(chk))) paste0(": ", chk[nzchar(chk)][[1L]]) else ""))
mon <- fit_opts[vapply(fit_opts, function(x) identical(x$track, "monthly"), TRUE)][[1L]]$opt
mon$grouping$enable <- TRUE; mon$grouping$source <- strrep("g", 240L); mon$grouping$group <- "east"
e3 <- tryCatch({ checkRunFileNames(mon, "monthly"); NULL }, error = function(e) conditionMessage(e))
expect(!is.null(e3) && grepl("makeMonthlyRunFile", e3),
       "stops on a monthly config whose names are over the limit, naming the builder")
expect(inherits(tryCatch(checkRunFileNames(mon, "nope"), error = identity), "error"),
       "an unknown track stops")

# --- [4] the config tree itself -------------------------------------------------
cat("\n[4] file and directory names under analysis/configs/\n")
tree <- list.files(CONFIG_ROOT, recursive = TRUE, include.dirs = TRUE)
comp <- unique(unlist(strsplit(tree, "/", fixed = TRUE)))
lc <- nchar(comp, type = "bytes")
expect(all(lc <= LIMIT), sprintf("%d names, longest %d bytes (%s)", length(comp), max(lc),
                                 comp[which.max(lc)]))

cat(sprintf("\nPASS: %d   FAIL: %d\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
