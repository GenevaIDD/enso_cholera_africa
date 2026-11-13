# test_data_vintage_default.R -- the data leg reads one outbreak-database vintage by default.
#
# HARD gate (see analysis/tests/run_all.R). PIPELINE.md runs the data-leg
# scripts without `--indir`, so per-script defaults could mix two database
# vintages (LP ids of a newer refresh unmatched in an older crosswalk are
# dropped silently by 03). utils/paths.R::defaultOutbreakDataDir() is the one
# default, and 00 / 01 / 02
# crosswalk stamp the directory they read into generated_data/prep_counts/<step>.yml.
# What this pins:
#   1. defaultOutbreakDataDir() is ~/secure/cholera_outbreak_data/outbreak_data_Nov_2026;
#   1b. CHOLERA_OUTBREAK_DATA_DIR replaces that root in paths.R and config.R;
#   2. no numbered script (analysis/NN*.R) and no maintained bash driver holds an
#      `outbreak_data_<Mon>_<YYYY>` string literal (comments are free);
#   3. every numbered script with an `--indir` option defaults it to
#      defaultOutbreakDataDir(), and sources utils/paths.R before its option list;
#   4. every writeStepCounts() call in a numbered script passes `indir =`;
#   5. writeStepCounts(indir =) records indir / data_vintage, readStepCounts()
#      returns the counts-only shape with or without them, readStepProvenance()
#      reads them (NA in an old file), the reserved names are refused, and 09's
#      report_cascade() gives the same table with or without the fields.
#
# Pure parse + temp files -- no data/, no secure files, no fit.
# Run: Rscript analysis/tests/test_data_vintage_default.R

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/tests/helpers_parse.R")   # parse_data_of()

VINTAGE_RX <- "outbreak_data_[A-Z][a-z]{2}_[0-9]{4}"

.env_root <- Sys.getenv("CHOLERA_OUTBREAK_DATA_DIR", NA)
Sys.unsetenv("CHOLERA_OUTBREAK_DATA_DIR")   # [1] pins the default root
env <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env)
sys.source("analysis/utils/paths.R", env)
sys.source("analysis/utils/provenance.R", env)

# 1. the default ---------------------------------------------------------------
cat("\n[1] defaultOutbreakDataDir()\n")
d <- env$defaultOutbreakDataDir()
expect(identical(d, "~/secure/cholera_outbreak_data/outbreak_data_Nov_2026"),
       sprintf("defaultOutbreakDataDir() is the Nov_2026 snapshot (got %s)", d))
expect(identical(basename(d), env$OUTBREAK_DATA_VINTAGE), "its basename is OUTBREAK_DATA_VINTAGE")
expect(identical(normalizePath(dirname(path.expand(d)), mustWork = FALSE),
                 normalizePath(sub("/+$", "", env$SECURE_OUTBREAK_DIR), mustWork = FALSE)),
       "it lies directly under SECURE_OUTBREAK_DIR")

# 1b. CHOLERA_OUTBREAK_DATA_DIR replaces the root (a synthetic world) ---------
cat("\n[1b] CHOLERA_OUTBREAK_DATA_DIR\n")
.alt <- file.path(tempdir(), "synthetic_secure")
Sys.setenv(CHOLERA_OUTBREAK_DATA_DIR = paste0(.alt, "/"))
env_alt <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env_alt)
sys.source("analysis/utils/paths.R", env_alt)
expect(identical(env_alt$defaultOutbreakDataDir(), file.path(.alt, env$OUTBREAK_DATA_VINTAGE)),
       "defaultOutbreakDataDir() sits under the override (trailing slash dropped)")
expect(identical(env_alt$SECURE_OUTBREAK_DIR, paste0(.alt, "/")),
       "SECURE_OUTBREAK_DIR is the override with one trailing slash")
expect(identical(as.character(env_alt$makeOutbreakMetadataFile()), file.path(.alt, "outbreak_metadata.rds")),
       "the root-level builders follow it (makeOutbreakMetadataFile())")
env_cfg_alt <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env_cfg_alt)
sys.source("analysis/utils/config.R", env_cfg_alt)
.d_alt <- env_cfg_alt$get_default_options()
expect(identical(.d_alt$simulator$location_frame$outbreak_metadata_path, file.path(.alt, "outbreak_metadata.rds")) &&
         identical(dirname(.d_alt$population$lp_shapefiles), env_alt$defaultOutbreakDataDir()),
       "config.R's two defaults under the root follow it too")
Sys.unsetenv("CHOLERA_OUTBREAK_DATA_DIR")
if (!is.na(.env_root)) Sys.setenv(CHOLERA_OUTBREAK_DATA_DIR = .env_root)

# 2. no vintage literal in numbered scripts / maintained bash --------------------
cat("\n[2] no per-script vintage literal\n")
numbered <- list.files("analysis", pattern = "^[0-9][0-9].*\\.R$", full.names = TRUE)
expect(length(numbered) > 40L, sprintf("scanned %d numbered scripts", length(numbered)))
str_hits <- character()
for (f in numbered) {
  pd <- parse_data_of(f)
  if (is.null(pd)) { str_hits <- c(str_hits, sprintf("%s: unparseable", f)); next }
  s  <- pd[pd$token == "STR_CONST" & grepl(VINTAGE_RX, pd$text), , drop = FALSE]
  if (nrow(s)) str_hits <- c(str_hits, sprintf("%s:%d %s", f, s$line1, s$text))
}
expect(length(str_hits) == 0L,
       sprintf("no numbered script holds an outbreak_data_<Mon>_<YYYY> string literal (%s)",
               if (length(str_hits)) paste(str_hits, collapse = "; ") else "none"))

bash_files <- list.files(file.path("analysis/bash", c("stages", "campaigns", "tools", "slurm")),
                         pattern = "\\.sh$", full.names = TRUE)
bash_hits <- character()
for (f in bash_files) {
  l <- readLines(f, warn = FALSE)
  k <- which(grepl(VINTAGE_RX, l) & !grepl("^\\s*#", l))
  if (length(k)) bash_hits <- c(bash_hits, sprintf("%s:%d", f, k))
}
expect(length(bash_hits) == 0L,
       sprintf("no maintained bash driver names a vintage outside comments (%s)",
               if (length(bash_hits)) paste(bash_hits, collapse = "; ") else "none"))

# 3. --indir defaults ---------------------------------------------------------
cat("\n[3] every --indir defaults to defaultOutbreakDataDir()\n")
EXPECTED_INDIR <- c("00_ingest_outbreakextractr.R", "01_preprocess_outbreaks.R",
                    "01b_preprocess_outbreaks_monthly.R", "01_preprocess_national_zero_months.R",
                    "01_preprocess_figures_data_overview.R", "02_preprocess_climate.R",
                    "02_preprocess_crosswalk.R")
walk_calls <- function(e, fn, visit) {
  kids <- NULL
  if (is.call(e)) {
    if (is.name(e[[1]]) && identical(as.character(e[[1]]), fn)) visit(e)
    kids <- as.list(e)[-1]
  } else if (is.expression(e) || is.list(e) || is.pairlist(e)) {
    kids <- as.list(e)
  }
  for (i in seq_along(kids))   # skip empty arguments (e.g. function formals)
    if (!identical(kids[[i]], quote(expr = ))) walk_calls(kids[[i]], fn, visit)
}
with_indir <- character()
for (f in numbered) {
  ex <- parse(f, keep.source = FALSE)
  found <- list()
  walk_calls(ex, "make_option", function(cl) {
    args  <- as.list(cl)[-1]
    nm    <- names(args) %||% rep("", length(args))
    os    <- args[nm == "opt_str"]
    if (!length(os)) os <- args[nm == ""][1]
    # opt_str is "--x" or c("-x", "--x")
    flags <- tryCatch(as.character(eval(os[[1]], baseenv())), error = function(e) character())
    if ("--indir" %in% flags) found[[length(found) + 1L]] <<- cl
  })
  if (!length(found)) next
  with_indir <- c(with_indir, basename(f))
  for (cl in found) {
    def <- as.list(cl)$default
    expect(is.call(def) && identical(deparse(def), "defaultOutbreakDataDir()"),
           sprintf("%s: --indir default is defaultOutbreakDataDir() (got %s)",
                   basename(f), paste(deparse(def), collapse = " ")))
  }
  # paths.R is sourced (top level) before the option list is built.
  top <- as.list(ex)
  src_paths <- which(vapply(top, function(e) is.call(e) && identical(e[[1]], as.name("source")) &&
                              identical(e[[2]], "analysis/utils/paths.R"), logical(1)))
  opt_list  <- which(vapply(top, function(e) {
    hit <- FALSE; walk_calls(e, "make_option", function(cl) hit <<- TRUE); hit }, logical(1)))
  expect(length(src_paths) && length(opt_list) && min(src_paths) < min(opt_list),
         sprintf("%s sources utils/paths.R before its options", basename(f)))
}
expect(setequal(with_indir, EXPECTED_INDIR),
       sprintf("the --indir scripts are the seven data-leg scripts (found: %s)",
               paste(sort(with_indir), collapse = ", ")))

# 4. step counts carry the indir --------------------------------------------------
cat("\n[4] writeStepCounts() callers pass indir =\n")
n_calls <- 0L
for (f in numbered) {
  walk_calls(parse(f, keep.source = FALSE), "writeStepCounts", function(cl) {
    n_calls <<- n_calls + 1L
    expect("indir" %in% names(as.list(cl)),
           sprintf("%s: writeStepCounts(%s, ...) passes indir =", basename(f), deparse(cl[[2]])))
  })
}
expect(n_calls >= 3L, sprintf("found %d writeStepCounts() calls (00, 01, 02 crosswalk)", n_calls))

# 5. round trip + backward compatibility --------------------------------------------
cat("\n[5] prep_counts provenance round trip\n")
tmp <- tempfile("prep_counts_"); dir.create(tmp)
cnt <- list(n_lp_in = 10L, n_lp_matched = 8L)
suppressMessages(env$writeStepCounts("02_preprocess_crosswalk", cnt, dir = tmp,
                                     indir = "~/secure/cholera_outbreak_data/outbreak_data_Aug_2026/"))
raw <- yaml::read_yaml(env$stepCountsPath("02_preprocess_crosswalk", tmp))
expect(identical(raw$data_vintage, "outbreak_data_Aug_2026"), "data_vintage is the indir's basename")
expect(identical(raw$indir, normalizePath(path.expand("~/secure/cholera_outbreak_data/outbreak_data_Aug_2026"),
                                          mustWork = FALSE)),
       "indir is tilde-expanded, without a trailing slash")
rd <- env$readStepCounts("02_preprocess_crosswalk", tmp)
expect(identical(names(rd), c("step", "written", names(cnt))),
       "readStepCounts() returns step, written, counts (provenance left out)")
pv <- env$readStepProvenance("02_preprocess_crosswalk", tmp)
expect(identical(pv$data_vintage, "outbreak_data_Aug_2026") && identical(pv$indir, raw$indir),
       "readStepProvenance() returns indir / data_vintage")

# an old file (no provenance) and a file written without indir =
suppressMessages(env$writeStepCounts("00_ingest_outbreakextractr", cnt, dir = tmp))
old <- yaml::read_yaml(env$stepCountsPath("00_ingest_outbreakextractr", tmp))
expect(!any(c("indir", "data_vintage") %in% names(old)), "indir = NULL writes the pre-2026-10-27 shape")
expect(identical(env$readStepCounts("00_ingest_outbreakextractr", tmp), old),
       "readStepCounts() of an old file is the file itself")
pv0 <- env$readStepProvenance("00_ingest_outbreakextractr", tmp)
expect(is.na(pv0$indir) && is.na(pv0$data_vintage), "readStepProvenance() of an old file: NA fields")
expect(is.null(env$readStepProvenance("01_preprocess_outbreaks", tmp)), "absent step: NULL")
expect(inherits(try(suppressMessages(env$writeStepCounts("x", list(indir = 1L), dir = tmp)),
                    silent = TRUE), "try-error"), "a count named indir is refused")
expect(inherits(try(suppressMessages(env$writeStepCounts("x", cnt, dir = tmp, indir = c("a", "b"))),
                    silent = TRUE), "try-error"), "a non-scalar indir is refused")

# report_cascade(): same table whether the files carry provenance or not.
renv <- new.env(parent = env)
sys.source("analysis/utils/results_derive.R", renv)
ctx_of <- function(dir) list(PREP_COUNTS_DIR = dir,
                             SIZE_DATA_PATH = file.path(dir, "none_siz.rds"),
                             SURV_PATH      = file.path(dir, "none_occ.rds"),
                             DUR_DATA_PATH  = file.path(dir, "none_dur.rds"),
                             MON_DATA_PATH  = file.path(dir, "none_mon.rds"))
with_pv <- tempfile("pc_with_"); dir.create(with_pv)
no_pv   <- tempfile("pc_without_"); dir.create(no_pv)
for (st in renv$.CASCADE_STEPS) {
  suppressMessages(env$writeStepCounts(st, cnt, dir = with_pv, indir = env$defaultOutbreakDataDir()))
  suppressMessages(env$writeStepCounts(st, cnt, dir = no_pv))
}
quiet_cascade <- function(dir) {
  w <- character()
  out <- withCallingHandlers(renv$report_cascade(ctx_of(dir)),
                             warning = function(cond) { w <<- c(w, conditionMessage(cond))
                                                        invokeRestart("muffleWarning") })
  list(tbl = out, warnings = w)
}
a <- quiet_cascade(with_pv); b <- quiet_cascade(no_pv)
expect(identical(a$tbl, b$tbl), "report_cascade(): identical table with and without the provenance fields")
expect(nrow(a$tbl) == length(renv$.CASCADE_STEPS) * length(cnt) && !anyNA(a$tbl$value),
       "report_cascade(): only the counts, no NA metric from indir / data_vintage")
expect(!any(grepl("coercion", a$warnings)), "report_cascade(): no coercion warning")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
