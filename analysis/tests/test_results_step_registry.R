# test_results_step_registry.R -- the phase-09 step registry is complete,
# acyclic and in the historical dispatch order.
#
# HARD gate (see analysis/tests/run_all.R). results_dispatch() runs
# .RESULTS_STEP_TABLE (utils/results_derive.R), expanded by results_steps(),
# which 09a turns into one SLURM task per step and the stage driver runs as wave-chained arrays.
# What this pins:
#   1. every RESULTS_SCHEMA kind is written by exactly one step (an artefact
#      no step owns would silently never be derived on the HPC path);
#   2. step ids are unique, every `after` names a step, every step's wave is
#      above its dependencies' (so a wave's jobs never derive each other's
#      artefacts concurrently);
#   3. the expanded order is the former hard-coded dispatch order (a full run
#      is call-for-call identical to before the registry);
#   4. each step's runner calls its own ensure_*() (and only that one);
#   5. the hindcast step disappears when results.scope.hindcast.enable is off;
#   6. results_select_steps() rejects an unknown id and keeps registry order.
#
# Pure parse -- no data, no fit.
# Run: Rscript analysis/tests/test_results_step_registry.R

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

env <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env)
sys.source("analysis/utils/paths.R", env)
sys.source("analysis/utils/results_derive.R", env)
schema_kinds <- local({
  # RESULTS_SCHEMA without sourcing results_core.R's closure: parse its literal.
  ex <- parse("analysis/utils/results_core.R")
  hit <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name("RESULTS_SCHEMA")), as.list(ex))
  env$RESULTS_SCHEMA <- eval(hit[[1]][[3]], baseenv())
  names(env$RESULTS_SCHEMA)
})

fake_ctx <- function(tracks = c("occurrence", "duration", "monthly"), hindcast = TRUE,
                     kl_sweep = TRUE, events = TRUE, af_posterior = TRUE, attrib = TRUE,
                     lag_avg = TRUE, panel_sidecars = TRUE)
  list2env(list(FIG1_TRACKS = tracks, HC_SCOPE = list(enable = hindcast), AF_POSTERIOR = af_posterior,
                NZ_AUDIT_PATH = if (panel_sidecars) "x.audit.csv" else "",
                OCC_EXCL_PATH = if (panel_sidecars) "x.adm1_overlap.csv" else "",
                CONTAGION_COV_PATH = if (panel_sidecars) "x.contagion_coverage.csv" else "",
                LAG_AVG = list(enable = lag_avg), KL_SWEEP_CELL = if (kl_sweep) "ns5_K5" else "",
                ATTRIB_FACTUAL = if (attrib) "attrib" else "forecast",
                KL_SWEEP_DIR = if (kl_sweep) "generated_data/kl_sweep_compare/x" else "",
                EVENT_MANIFEST_PATH = if (events) "analysis/configs/event_attribution_horizons_x.yml" else "",
                SETS = c("neutral", "all"), opt = list(outdir = "generated_data/"),
                # report_who_scale is on with a smoothing run
                SMOOTHING_PARTS_DIR = "generated_data/sim_spatpomp/x/parts",
                # report_who_enso_ratio (fig. WR) is on with a smoothing run and groups
                WER_DERIV = list(groups = "enso_west"),
                CONFIG_PATH = "analysis/configs/results/manuscript.yml"), parent = emptyenv())

ctx   <- fake_ctx()
steps <- env$results_steps(ctx)
ids   <- names(steps)

cat("# [1] kind ownership\n")
owned <- unlist(lapply(steps, `[[`, "kinds"), use.names = FALSE)
table_kinds <- unlist(lapply(env$.RESULTS_STEP_TABLE, `[[`, "kinds"))
expect(!anyDuplicated(table_kinds), "each kind appears in one registry row")
expect(setequal(table_kinds, schema_kinds),
       sprintf("registry kinds == RESULTS_SCHEMA kinds (missing: %s; extra: %s)",
               paste(setdiff(schema_kinds, table_kinds), collapse = ","),
               paste(setdiff(table_kinds, schema_kinds), collapse = ",")))
expect(setequal(unique(owned), schema_kinds), "the expanded steps write every kind")

cat("# [2] ids, dependencies, waves\n")
expect(!anyDuplicated(ids), sprintf("%d unique step ids", length(ids)))
deps <- unlist(lapply(steps, `[[`, "after"))
expect(all(deps %in% ids), "every `after` names an existing step")
ok_wave <- vapply(steps, function(s)
  !length(s$after) || all(s$wave > vapply(steps[s$after], `[[`, integer(1), "wave")),
  logical(1))
expect(all(ok_wave), "every step's wave is above its dependencies'")
expect(identical(sort(ids[vapply(steps, `[[`, integer(1), "wave") == 2L]),
                 sort(c(paste0("report_effects:", ctx$FIG1_TRACKS),
                        "report_enso_exposure", "report_smoothing",
                        "report_convergence", "report_software", "report_lags", "modelled_units",
                        "seasonality", "report_who_enso_ratio", "oni_composite"))),
       "wave 2 = report_effects:<track>, report_enso_exposure, report_smoothing (after oni), report_convergence, report_software, report_lags, modelled_units (after onset_rate), seasonality (after rain_enso, occ_seasonality), report_who_enso_ratio (after group_labels), oni_composite (after oni)")
expect(inherits(try(env$results_step_waves(list(
  a = list(after = "b"), b = list(after = "a"))), silent = TRUE), "try-error"),
  "a dependency cycle stops")

cat("# [3] historical dispatch order\n")
tr <- ctx$FIG1_TRACKS
legacy <- c(paste0("enso_relation:", tr), paste0("track_params:", tr), "onset_rate",
            "attribution", "oni", "smoothing", "report_dataset",
            paste0("report_effects:", tr), "report_smoothing", "report_forecast",
            "report_enso_exposure", "hindcast",   # report_sensitive_pop unregistered
            # supplementary tables, appended after the historical order
            paste0("convergence:", tr), "report_convergence", "report_lags",
            "report_anchor", "report_software", "report_kl_sweep",
            # forecast vs counterfactual exposure, appended
            "report_attrib_exposure",
            # attribution metrics of the same runs, appended
            "report_attrib_metrics", "report_event_attrib", "report_event_af_model",
            # posterior AF / AN of fixed-draw runs, appended
            "af_posterior", "event_af_posterior",
            # posterior Childs population of the same runs, appended
            "attrib_exposure_posterior",
            # information component of the attribution factual, appended
            "report_attrib_information",
            # Figure 1d's rainfall tables, appended
            "rain_enso",
            # group display lags + per-country season ends, appended
            "group_labels",
            # lag averaging over the single-lag fits, appended
            paste0("lag_averaging:", tr),
            # membership uncertainty of the deployed classes, appended
            "report_membership",
            # per-country WHO term scale of the smoothing, appended
            "report_who_scale",
            # counterfactual analog sets, D29, appended
            "report_analog_sets",
            # occurrence-panel tables, appended
            "report_national_zero_audit", "report_occ_panel_exclusions", "report_contagion_coverage",
            # supplement maps MU / GM / KC: modelled units + fit membership, appended
            "modelled_units",
            # Figure 1c's onset-season column: occurrence seasonal spline, appended
            "occ_seasonality",
            # supplement fig. SE: rainfall vs the three tracks by month, appended
            "seasonality",
            # supplement fig. WR: observed vs modelled ENSO ratio of 2026, appended
            "report_who_enso_ratio",

            # Figure 1c's composite ONI row, appended
            "oni_composite",
            # forecast vs counterfactual incidence by country, appended
            "report_attrib_incidence",
            # ENSO contrasts, all vs confirmed outbreaks, appended
            "report_confirmed_agreement")
expect(identical(ids, legacy), "expanded order == former results_dispatch() order + supplementary steps")

cat("# [4] each runner calls its own ensure_*()\n")
for (s in steps) {
  called <- grep("^ensure_", all.names(body(s$run)), value = TRUE)
  expect(identical(called, paste0("ensure_", s$stem)) &&
           exists(called, envir = env, inherits = FALSE),
         sprintf("%s runs ensure_%s()", s$id, s$stem))
}

cat("# [5] hindcast toggle and track scope\n")
expect(!"hindcast" %in% names(env$results_steps(fake_ctx(hindcast = FALSE))),
       "hindcast step absent when results.scope.hindcast.enable is false")
expect(!any(grepl("^lag_averaging:", names(env$results_steps(fake_ctx(lag_avg = FALSE))))),
       "lag_averaging steps absent when results.derivation.lag_averaging.enable is false")
expect(!"report_kl_sweep" %in% names(env$results_steps(fake_ctx(kl_sweep = FALSE))),
       "report_kl_sweep step absent when results.kl_sweep_dir is empty")
expect(!"report_membership" %in% names(env$results_steps(fake_ctx(kl_sweep = FALSE))),
       "report_membership step absent when results.kl_sweep_dir / kl_sweep_cell are empty")
expect(!any(c("report_national_zero_audit", "report_occ_panel_exclusions", "report_contagion_coverage") %in% names(env$results_steps(fake_ctx(panel_sidecars = FALSE)))),
       "the occurrence-panel table steps are absent when their sidecars do not exist (an older run)")
expect(!any(c("report_event_attrib", "report_event_af_model") %in% names(env$results_steps(fake_ctx(events = FALSE)))),
       "the event steps are absent when results.event_manifest is empty")
expect(!any(c("af_posterior", "event_af_posterior", "attrib_exposure_posterior") %in%
              names(env$results_steps(fake_ctx(af_posterior = FALSE)))),
       "the posterior AF steps are absent when results.af_posterior is false")
expect(!"report_attrib_information" %in% names(env$results_steps(fake_ctx(attrib = FALSE))) &&
         !"report_attrib_information" %in% names(env$results_steps(fake_ctx(af_posterior = FALSE))),
       "report_attrib_information needs an attribution factual and results.af_posterior")
no_ev <- names(env$results_steps(fake_ctx(events = FALSE)))
expect("af_posterior" %in% no_ev && !"event_af_posterior" %in% no_ev,
       "event_af_posterior also needs results.event_manifest; af_posterior does not")
one <- env$results_steps(fake_ctx(tracks = "occurrence"))
expect(identical(one[["report_enso_exposure"]]$after, "enso_relation:occurrence"),
       "`enso_relation:*` follows fig1_tracks")

cat("# [6] selection\n")
expect(inherits(try(env$results_select_steps(ctx, "bogus"), silent = TRUE), "try-error"),
       "an unknown step id stops")
sel <- env$results_select_steps(ctx, c("oni", "enso_relation:duration"))
expect(identical(names(sel), c("enso_relation:duration", "oni")),
       "a selection keeps registry order")

cat("# [7] artefact paths\n")
p <- env$results_step_paths(ctx, steps[["hindcast"]])
expect(length(p) == 4L && all(startsWith(p, "generated_data//results/manuscript/")),
       "hindcast writes four artefacts under the config's results dir")
expect(identical(env$makeResultsRunFile(ctx$opt, "results_tasks_dir", config_path = ctx$CONFIG_PATH),
                 "generated_data//results/manuscript/tasks"),
       "results_tasks_dir sits in the config's results dir")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
