# test_hindcast_iterator.R -- utils/sim_io.R's hindcast manifest helpers.
#
# HARD gate (see analysis/tests/run_all.R). read_hindcast_manifest() and
# for_each_hindcast_run() carry the group x launch-year loop shared by 07c's
# run_hindcast_skill() and 07d's run_hindcast_figures(). A miscounted skip or a run
# silently not scored would only show up
# as a thinner hindcast table on the HPC, so the skip paths are pinned here on a
# synthetic manifest: missing config, ready merged run, run built by
# `ensure_summary`, run with no summary at all.
#
# Synthetic temp files only -- no fit, no generated_data, no network.
# Run: Rscript analysis/tests/test_hindcast_iterator.R

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(qs2) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "config", "paths", "provenance", "gid", "occ_data", "sim_io"))
  source(sprintf("analysis/utils/%s.R", m))

cat("# hindcast manifest iteration (read_hindcast_manifest / for_each_hindcast_run)\n")

tmp <- tempfile("hcit_"); dir.create(tmp)
cfg_for <- function(tag) {
  p <- file.path(tmp, paste0(tag, ".yml"))
  yaml::write_yaml(list(outdir = file.path(tmp, "gd/"),
                        simulator = list(run_tag = tag, output = list(format = "qs"))), p)
  p
}
ms <- tibble(replicate_id = 1L, gid = "A", month = "2018-01-01", Y = 1L)
md_ready <- file.path(tmp, "ready"); dir.create(md_ready)
qs_save(ms, file.path(md_ready, "monthly_summary.qs"))
man <- list(groups = list(
  g1 = list(years = list(
    `2017` = list(config = file.path(tmp, "missing.yml")),
    `2018` = list(config = cfg_for("ready"), merged_dir = md_ready))),
  g2 = list(years = list(
    `2019` = list(config = cfg_for("built")),
    `2020` = list(config = cfg_for("empty"))))))
man_path <- file.path(tmp, "m.yml"); yaml::write_yaml(man, man_path)
panel <- file.path(tmp, "panel.rds"); saveRDS(tibble(gid = "A"), panel)
ol <- list(optparse::make_option(c("-C", "--config"), type = "character"))

# read_hindcast_manifest: argument checks
expect(inherits(try(read_hindcast_manifest(file.path(tmp, "x.yml"), panel), silent = TRUE), "try-error"),
       "read_hindcast_manifest stops on a missing manifest")
expect(inherits(try(read_hindcast_manifest(man_path, file.path(tmp, "x.rds")), silent = TRUE), "try-error"),
       "read_hindcast_manifest stops on a missing obs panel")
m2 <- file.path(tmp, "m2.yml"); yaml::write_yaml(list(groups = list()), m2)
expect(inherits(try(read_hindcast_manifest(m2, panel), silent = TRUE), "try-error"),
       "read_hindcast_manifest stops on a manifest without groups")
man_r <- read_hindcast_manifest(man_path, panel)

built <- character()
ensure <- function(opt, md, fmt) {
  if (identical(opt$simulator$run_tag, "built")) {
    dir.create(md, recursive = TRUE, showWarnings = FALSE)
    qs_save(ms, file.path(md, paste0("monthly_summary", sim_io_ext(fmt))))
  }
  built <<- c(built, opt$simulator$run_tag)
}
warns <- character()
it <- withCallingHandlers(
  suppressMessages(for_each_hindcast_run(man_r, ol, "[t]",
    function(run) paste(run$group, run$year, basename(run$merged_dir)),
    announce = "before", ensure_summary = ensure)),
  warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })

expect(identical(it$n_ok, 2L) && identical(it$n_skip, 2L), "2 runs scored, 2 skipped")
expect(identical(unlist(it$results), c("g1 2018 ready", "g2 2019 merged")),
       "fn sees the manifest merged_dir, else the run's spatpomp merged dir")
expect(identical(built, c("built", "empty")), "ensure_summary is called only where the summary is missing")
expect(length(warns) == 2L && grepl("config missing", warns[1]) && grepl("no monthly_summary", warns[2]),
       "one warning per skip, naming the reason")

it2 <- suppressWarnings(suppressMessages(for_each_hindcast_run(man_r, ol, "[t]", function(run) NULL,
                                                                announce = "after")))
expect(identical(it2$n_ok, 2L) && identical(length(it2$results), 0L),
       "NULL returns are counted but not collected")

# Parallel runs (forked workers): same results, counts and order as serial, and
# the workers' skip warnings re-raised in the parent. The "built" summary exists
# by now, so ensure_summary is a no-op here.
warns_p <- character()
itp <- withCallingHandlers(
  suppressMessages(for_each_hindcast_run(man_r, ol, "[t]",
    function(run) paste(run$group, run$year, basename(run$merged_dir)),
    announce = "before", ensure_summary = function(opt, md, fmt) NULL, workers = 2L)),
  warning = function(w) { warns_p <<- c(warns_p, conditionMessage(w)); invokeRestart("muffleWarning") })
expect(identical(itp$n_ok, it$n_ok) && identical(itp$n_skip, it$n_skip) &&
         identical(itp$results, it$results),
       "workers = 2: same results, counts and order as the serial loop")
expect(length(warns_p) == 2L && grepl("config missing", warns_p[1]) && grepl("no monthly_summary", warns_p[2]),
       "workers = 2: skip warnings re-raised in the parent")

# Parallel dispatch order: largest merged monthly_summary first, results still
# in manifest order. Three runs whose summaries grow down the manifest.
md_sz <- vapply(c("s1", "s2", "s3"), function(tag) {
  d <- file.path(tmp, paste0("md_", tag)); dir.create(d)
  qs_save(tibble(replicate_id = seq_len(10^match(tag, c("s1", "s2", "s3"))), gid = "A",
                 month = "2018-01-01", Y = runif(10^match(tag, c("s1", "s2", "s3")))),
          file.path(d, "monthly_summary.qs"))
  d
}, character(1))
ents <- lapply(names(md_sz), function(tag) list(entry = list(merged_dir = md_sz[[tag]])))
ents[[4]] <- list(entry = list())   # no merged_dir in the manifest
sz <- hindcast_entry_size(ents)
expect(sz[3] > sz[2] && sz[2] > sz[1] && sz[4] == 0,
       "hindcast_entry_size(): summary bytes per entry, 0 without a merged_dir")
man_s <- list(groups = list(g = list(years = stats::setNames(lapply(names(md_sz), function(tag)
  list(config = cfg_for(tag), merged_dir = md_sz[[tag]])), c("2018", "2019", "2020")))))
# Submission order = the order mapWorkers() receives (it dispatches dynamically
# in input order); recorded through a wrapper so no worker timing is involved.
mw_orig <- mapWorkers; submitted <- NULL
mapWorkers <- function(X, FUN, workers = 1L) {
  submitted <<- vapply(X, function(e) e$entry$merged_dir, character(1))
  mw_orig(X, FUN, workers = workers)
}
its <- suppressMessages(for_each_hindcast_run(man_s, ol, "[t]", function(run) run$year, workers = 2L))
mapWorkers <- mw_orig
expect(identical(unlist(its$results), c("2018", "2019", "2020")),
       "workers = 2: results in manifest order after largest-first dispatch")
expect(identical(unname(submitted), unname(md_sz[c("s3", "s2", "s1")])),
       "workers = 2: runs submitted largest first")

# v2 (registry) manifest: runs keyed by unit; fn sees the unit's fit groups and
# every member country's firm window end for the launch year.
man_u <- list(groups = list(
  enso_east = list(iso3 = list("KEN"), years = list(`2018` = list(horizon_end = "2018-12-01"))),
  enso_horn = list(iso3 = list("ETH"), years = list(`2018` = list(horizon_end = "2019-06-01")))),
  units = list(eh = list(members = list("KEN", "ETH"), groups = list("enso_east", "enso_horn"),
                         years = list(`2018` = list(config = cfg_for("ready"), merged_dir = md_ready,
                                                    horizon_end = "2019-06-01")))),
  countries = list(KEN = list(group = "enso_east", unit = "eh"),
                   ETH = list(group = "enso_horn", unit = "eh")))
itu <- suppressMessages(for_each_hindcast_run(man_u, ol, "[t]", function(run)
  list(g = run$group, ug = run$unit_groups, fh = run$firm_horizon_by_country)))
r1 <- itu$results[[1L]]
expect(identical(itu$n_ok, 1L) && identical(r1$g, "eh"), "v2: one run keyed by its unit")
expect(identical(r1$ug, c("enso_east", "enso_horn")), "v2: run carries its fit groups")
expect(identical(r1$fh[c("KEN", "ETH")], c(KEN = "2018-12-01", ETH = "2019-06-01")),
       "v2: per-country firm windows for the launch year")
itv1 <- suppressWarnings(suppressMessages(for_each_hindcast_run(man_r, ol, "[t]",
  function(run) list(fh = run$firm_horizon_by_country, ug = run$unit_groups))))
expect(is.null(itv1$results[[1L]]$fh) && identical(itv1$results[[1L]]$ug, "g1"),
       "v1 runs carry no per-country firm map; unit_groups = the group")
expect(identical(r1$fh, itu$results[[1L]]$fh) &&
         identical(suppressMessages(for_each_hindcast_run(man_u, ol, "[t]", function(run)
           list(ly = run$launch_year, sc = run$scored_countries)))$results[[1L]],
           list(ly = 2018L, sc = NULL)),
       "v2: launch_year = the integer key; no scored_countries outside season windows")

# Season-window manifest (make_hindcast_configs.R --season_window): run keys
# "<year>-<MM>", only the scored countries read; the unit's other members get
# the month before the launch as horizon, so every censor drops them.
man_s <- man_u
man_s$season_end <- "season_window"
man_s$units$eh$years <- list(`2018-06` = list(config = cfg_for("ready"), merged_dir = md_ready,
                                              launch = "2018-06-01", horizon_end = "2019-06-01",
                                              launch_year = 2018L, scored = list("ETH")))
man_s$countries$ETH$years <- list(`2018-06` = list(firm_horizon_end = "2019-06-01"))
man_s$countries$KEN$years <- list(`2018-05` = list(firm_horizon_end = "2019-01-01"))
its <- suppressMessages(for_each_hindcast_run(man_s, ol, "[t]", function(run)
  list(y = run$year, ly = run$launch_year, sc = run$scored_countries, fh = run$firm_horizon_by_country)))
rs <- its$results[[1L]]
expect(identical(rs$y, "2018-06") && identical(rs$ly, 2018L) && identical(rs$sc, "ETH"),
       "season windows: run key, integer season year, scored countries")
expect(identical(rs$fh[c("ETH", "KEN")], c(ETH = "2019-06-01", KEN = "2018-05-01")),
       "season windows: the scored country's firm end; the unscored member padded to launch - 1 month")
expect(identical(manifest_country_horizons(man_s, "firm", year = "2018-05"), c(KEN = "2019-01-01")),
       "manifest_country_horizons(): a run key returns only its own countries (no group fallback)")
df <- data.frame(gid = c("KEN.1_1", "ETH.1_1"), month = as.Date(c("2018-06-01", "2018-06-01")))
expect(identical(censor_sim_by_country(df, rs$fh)$gid, "ETH.1_1"),
       "season windows: censor_sim_by_country() drops the unscored member's rows")

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
