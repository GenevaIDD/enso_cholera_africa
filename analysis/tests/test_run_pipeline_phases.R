# test_run_pipeline_phases.R -- tools/run_pipeline.sh plans every artefact 09 reads.
#
# HARD gate (see analysis/tests/run_all.R). The SLURM relay chain is the
# whole rollout; a step it does not chain is silently missing from a launch:
#   [1] past ENSO event attribution: the sim phase runs `deploy_prod_run.R
#       --stage gen_events` and submits stages/07_event_attribution.sh (with a
#       --list preflight) when deploy.event_attribution.events is set, and the
#       results phase waits on those arrays. --event-replicates N reaches the
#       stage driver as --n_replicates N; without it --plan says the configs'
#       simulator.n_replicates is used.
#   [2] a deploy without events (the block absent, or events "") plans no
#       event step and says so.
#   [3b] --after queues the first phase as a relay job; bad ids and --phase refused.
#   [4] --plan submits nothing: a stub
#       `sbatch` first on PATH records any call and fails loudly.
#   [5] stages/07_event_attribution.sh skips a config whose run has a _SUCCESS
#       sentinel in every partition (as 07_forecast_groups.sh), before the
#       Lenssen cache check; --list marks it [skip]; an incomplete run (a
#       missing country dir or repoff chunk) is submitted; --overwrite lists
#       every config. Temp configs + a recording sbatch stub that prints a fake
#       job id (no real submission).
#   [6] the attrib phase with no attribution array (fake drivers in a temp
#       copy of the repo layout): the empty id capture reaches the no-array
#       branch (stops under --until attrib, runs results inline under --until
#       results) and exits 0; with an array id it relays results afterok on it.
#       Every `ids=$(...)` / `id=$(...)` capture in run_pipeline.sh ends in
#       `|| true` (an empty grep under pipefail would abort the script).
#
# Committed deploy YAMLs are only read; the no-event variants are temp copies.
# No fit, no data, no cluster. Shells out to bash (each plan ~2 s).
# Run: Rscript analysis/tests/test_run_pipeline_phases.R

suppressPackageStartupMessages(library(yaml))

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# run_pipeline.sh: event attribution chained\n")

DRIVER  <- "analysis/bash/tools/run_pipeline.sh"
DEP     <- "analysis/configs/deploy/prod_run.yml"
RES     <- "analysis/configs/results/manuscript.yml"

# --- stub sbatch: any call is recorded and fails ------------------------------
STUB_DIR <- tempfile("stub_bin_")
dir.create(STUB_DIR)
MARKER <- file.path(STUB_DIR, "sbatch_called")
writeLines(c("#!/bin/bash",
             sprintf("echo \"$*\" >> '%s'", MARKER),
             "echo 'STUB sbatch: a --plan / --dry-run must never submit' >&2",
             "exit 97"),
           file.path(STUB_DIR, "sbatch"))
Sys.chmod(file.path(STUB_DIR, "sbatch"), "0755")
TMP <- tempfile("run_pipeline_test_")
dir.create(TMP)
# (no top-level on.exit: it fires at once outside a function; R removes its
# tempdir() -- STUB_DIR, TMP -- when the session ends)

run_bash <- function(script, args) {
  out <- suppressWarnings(system2("bash", c(script, args), stdout = TRUE, stderr = TRUE,
                                  env = sprintf("PATH=%s:%s", STUB_DIR, Sys.getenv("PATH"))))
  st <- attr(out, "status")
  list(out = out, status = if (is.null(st)) 0L else st)
}
plan <- function(...) run_bash(DRIVER, c(..., "--plan"))
has  <- function(r, pattern) any(grepl(pattern, r$out, fixed = TRUE))
line <- function(r, prefix) grep(prefix, r$out, fixed = TRUE, value = TRUE)

# --- [1] prod_run: events in the sim phase -------------------------------------
cat("\n[1] prod_run --plan: gen_events + 07_event_attribution.sh in the sim phase\n")
r19 <- plan("--deploy", DEP, "--results-config", RES)
expect(r19$status == 0L, sprintf("--plan on %s exits 0", basename(DEP)))
sim19 <- line(r19, "  sim     :")
expect(length(sim19) == 1L && grepl("gen_hindcast, gen_events (-> analysis/configs/event_attribution_horizons_prod_run.yml)",
                                    sim19, fixed = TRUE),
       "sim runs gen_events after gen_hindcast, naming the deploy's event manifest")
expect(has(r19, sprintf("stages/07_event_attribution.sh --deploy %s --list (preflight)", DEP)),
       "sim preflights the event configs with 07_event_attribution.sh --list")
expect(has(r19, sprintf("then stages/07_event_attribution.sh --deploy %s", DEP)) &&
         has(r19, "events elnino_2015,lanina_2020_23,elnino_2023"),
       "sim submits 07_event_attribution.sh on the deploy's events")
expect(has(r19, "event replicates: each event config's simulator.n_replicates"),
       "without --event-replicates the plan says the configs' simulator.n_replicates is used")
expect(!any(grepl("--n_replicates", r19$out, fixed = TRUE)),
       "without --event-replicates no --n_replicates is passed")
expect(has(r19, "results waits afterok on the forecast arrays, the hindcast postprocess and the event arrays"),
       "results waits on the event arrays too")

cat("\n[1b] --event-replicates / --concurrency reach the event stage driver\n")
r18 <- plan("--deploy", DEP, "--from", "sim", "--event-replicates", "9950", "--concurrency", "7")
expect(r18$status == 0L, sprintf("--plan on %s --from sim exits 0", basename(DEP)))
ev18 <- grep("then stages/07_event_attribution.sh", r18$out, fixed = TRUE, value = TRUE)
expect(length(ev18) == 1L && grepl("--concurrency 7 --n_replicates 9950", ev18, fixed = TRUE),
       "--event-replicates 9950 --concurrency 7 -> 07_event_attribution.sh --concurrency 7 --n_replicates 9950")
expect(has(r18, "event replicates: 9950 (--event-replicates)"), "the plan names the explicit replicate count")
rbad <- plan("--deploy", DEP, "--event-replicates", "0")
expect(rbad$status != 0L && has(rbad, "--event-replicates must be a positive integer"),
       "--event-replicates 0 is refused")

# Static: the executed sim phase submits what the plan prints, and the relay
# carries --event-replicates to the next phases.
src <- readLines(DRIVER, warn = FALSE)
expect(any(grepl('07_event_attribution.sh --deploy "$DEPLOY" "${EVENT_PASS[@]}"', src, fixed = TRUE)) &&
         any(grepl('07_event_attribution.sh --deploy "$DEPLOY" --list', src, fixed = TRUE)),
       "run_phase sim calls 07_event_attribution.sh (--list preflight, then the submit with EVENT_PASS)")
expect(any(grepl('echo "$ev" | grep -E', src, fixed = TRUE)),
       "run_phase sim collects the event array ids for the results relay")
expect(any(grepl('args+=(--event-replicates "$EVENT_REPLICATES")', src, fixed = TRUE)),
       "relay() forwards --event-replicates")

# --- [2] no events: skipped ----------------------------------------------------
cat("\n[2] a deploy without events plans no event step\n")
d19_lines <- readLines(DEP, warn = FALSE)
ev_idx <- grep("^\\s+events:\\s", d19_lines)
expect(length(ev_idx) == 1L, "prod_run.yml has one event_attribution.events line")
off_cfg <- file.path(TMP, "prod_runX_events_off.yml")
writeLines(sub("^(\\s+events:\\s+).*$", "\\1\"\"", d19_lines), off_cfg)
absent_cfg <- file.path(TMP, "prod_runX_events_absent.yml")
raw <- yaml::read_yaml(DEP)
raw$deploy$event_attribution <- NULL
yaml::write_yaml(raw, absent_cfg)
for (cfg in c(off_cfg, absent_cfg)) {
  r <- plan("--deploy", cfg, "--from", "sim", "--results-config", RES)
  lbl <- basename(cfg)
  expect(r$status == 0L, sprintf("%s: --plan exits 0", lbl))
  expect(has(r, "events: (skipped: deploy.event_attribution.events is unset)"),
         sprintf("%s: the plan reports the event step skipped", lbl))
  expect(!any(grepl("gen_events|07_event_attribution", r$out)),
         sprintf("%s: no gen_events, no 07_event_attribution.sh", lbl))
}

# --- recording sbatch stub: logs the call, prints a fake job id ---------------
OK_DIR <- tempfile("stub_ok_bin_")
dir.create(OK_DIR)
OK_MARKER <- file.path(OK_DIR, "sbatch_calls")
writeLines(c("#!/bin/bash",
             sprintf("echo \"$*\" >> '%s'", OK_MARKER),
             "echo 4242"),
           file.path(OK_DIR, "sbatch"))
Sys.chmod(file.path(OK_DIR, "sbatch"), "0755")
ok_calls <- function() if (file.exists(OK_MARKER)) readLines(OK_MARKER) else character(0)
# Run `bash <script> <args>` from `root`, the recording stub first on PATH.
run_in <- function(root, script, args, extra_env = character()) {
  cmd <- sprintf("cd %s && bash %s %s", shQuote(root), shQuote(script),
                 paste(shQuote(args), collapse = " "))
  out <- suppressWarnings(system2("bash", c("-c", shQuote(cmd)), stdout = TRUE, stderr = TRUE,
                                  env = c(sprintf("PATH=%s:%s", OK_DIR, Sys.getenv("PATH")), extra_env)))
  st <- attr(out, "status")
  list(out = out, status = if (is.null(st)) 0L else st)
}

# --- [5] 07_event_attribution.sh skip-complete guard ---------------------------
cat("\n[5] stages/07_event_attribution.sh skips complete runs\n")
EV_DRIVER <- normalizePath("analysis/bash/stages/07_event_attribution.sh")
EVR <- file.path(TMP, "ev")
dir.create(file.path(EVR, "configs"), recursive = TRUE)
writeLines("cache", file.path(EVR, "cache.rds"))
ev_cfg <- function(name, tag, cache = NULL) {
  writeLines(c(sprintf("outdir: \"%s/out/\"", EVR),
               "simulator:",
               sprintf("  run_tag: %s", tag),
               "  scope:",
               "    countries: KEN,TZA",
               if (!is.null(cache)) c("  enso:", "    counterfactual:", sprintf("      cache_path: %s", cache))),
             file.path(EVR, "configs", paste0(name, ".yml")))
}
sentinel <- function(tag, rel) {
  d <- file.path(EVR, "out", "sim_spatpomp", tag, "parts", rel)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(d, "_SUCCESS"))
}
# ev1 factual: complete (both countries).
ev_cfg("event_ev1_uA", "tagA"); sentinel("tagA", "tu=KEN"); sentinel("tagA", "tu=TZA")
# ev1 counterfactual: incomplete (TZA chunk repoff=50 has no _SUCCESS).
ev_cfg("event_ev1_uA_cf-neutral", "tagB", file.path(EVR, "cache.rds"))
sentinel("tagB", "tu=KEN"); sentinel("tagB", "tu=TZA/repoff=0")
dir.create(file.path(EVR, "out/sim_spatpomp/tagB/parts/tu=TZA/repoff=50"), recursive = TRUE)
# ev2 counterfactual: complete, its Lenssen cache missing (the guard runs first).
ev_cfg("event_ev2_uA_cf-neutral", "tagC", file.path(EVR, "no_such_cache.rds"))
sentinel("tagC", "tu=KEN"); sentinel("tagC", "tu=TZA")
# ev3 factual: one of two country partitions written -> incomplete.
ev_cfg("event_ev3_uA", "tagD"); sentinel("tagD", "tu=KEN")

rl <- run_in(EVR, EV_DRIVER, c("--config-dir", "configs", "--list"))
expect(rl$status == 0L, "--list exits 0")
expect(has(rl, "[skip] spev_ev1_uA: run complete") && has(rl, "[skip] spev_ev2_uA_cf-neutral: run complete"),
       "--list marks the two complete runs [skip] (the complete counterfactual despite its missing cache)")
expect(has(rl, "[list] spev_ev1_uA_cf-neutral") && has(rl, "[list] spev_ev3_uA"),
       "--list lists the incomplete runs (a missing repoff chunk; a missing country dir)")
expect(has(rl, "[list] 2 complete run(s) would be skipped.") && length(ok_calls()) == 0L,
       "--list reports the skip count and calls no sbatch")

rs <- run_in(EVR, EV_DRIVER, c("--config-dir", "configs"))
calls <- ok_calls()
expect(rs$status == 0L && has(rs, "Submitted 2 array job(s); skipped 2 complete run(s)."),
       "submit: 2 incomplete runs submitted, 2 complete runs skipped")
expect(length(calls) == 2L && any(grepl("event_ev1_uA_cf-neutral.yml", calls, fixed = TRUE)) &&
         any(grepl("event_ev3_uA.yml", calls, fixed = TRUE)) &&
         !any(grepl("event_ev1_uA.yml|event_ev2_uA_cf-neutral.yml", calls)),
       "sbatch (stub) called only for the two incomplete configs")
expect(sum(grepl("array job id: 4242", rs$out, fixed = TRUE)) == 2L,
       "each submission prints its array job id (what run_pipeline.sh collects)")

ro <- run_in(EVR, EV_DRIVER, c("--config-dir", "configs", "--event", "ev1", "--overwrite", "--list"))
expect(ro$status == 0L && has(ro, "[list] spev_ev1_uA ") && has(ro, "[list] spev_ev1_uA_cf-neutral") &&
         has(ro, "[list] 0 complete run(s) would be skipped."),
       "--overwrite lists the complete run too (nothing skipped)")
expect(length(ok_calls()) == 2L, "--list / --overwrite --list add no sbatch call")

# --- [6] attrib phase: empty id capture reaches the no-array branch ------------
cat("\n[6] attrib phase without an attribution array\n")
FR <- file.path(TMP, "fake_repo")
for (d in c("analysis/bash/tools", "analysis/bash/stages", "analysis/bash/lib", "analysis/utils"))
  dir.create(file.path(FR, d), recursive = TRUE)
file.copy(DRIVER, file.path(FR, DRIVER))
file.copy("analysis/bash/lib/cluster.sh", file.path(FR, "analysis/bash/lib/cluster.sh"))   # no profile: plain sbatch
writeLines(c("a <- commandArgs(TRUE); st <- a[which(a == '--stage') + 1L]",
             "if (identical(st, 'resolve')) cat('DEPLOY_TAG=prod_runT', 'HISTORICAL_TEMPLATE=hist.yml',",
             "  'FORECAST_ATTRIB_MANIFEST=attrib_manifest.yml', 'KL_SWEEP_BATCH=', 'KL_SWEEP_CELL=',",
             "  'EVENTS=', 'EVENT_CONFIG_DIR=', 'EVENT_MANIFEST=', sep = '\\n') else",
             "  cat('[fake deploy] stage', st, '\\n')"),
           file.path(FR, "analysis/utils/deploy_prod_run.R"))
writeLines(c("#!/bin/bash",
             "echo \"[skip] spfc_fa_x: run complete; --overwrite resubmits it\"",
             "[[ -n \"${FAKE_ARRAY:-}\" ]] && echo \"    array job id: ${FAKE_ARRAY}\"",
             "echo \"Submitted 0 array job(s); skipped 1 complete run(s).\"",
             "exit 0"),
           file.path(FR, "analysis/bash/stages/07_forecast_groups.sh"))
writeLines(c("#!/bin/bash", "echo 'RESULTS_JOB_IDS=11:12'"),
           file.path(FR, "analysis/bash/stages/09_postprocess_results.sh"))
writeLines("deploy: {}", file.path(FR, "deploy.yml"))
writeLines("results: {}", file.path(FR, "res.yml"))
fake_args <- function(until) c("--deploy", "deploy.yml", "--from", "attrib", "--until", until,
                               "--results-config", "res.yml")
FR_DRIVER <- file.path(FR, DRIVER)
n0 <- length(ok_calls())

ra <- run_in(FR, FR_DRIVER, fake_args("attrib"))
expect(ra$status == 0L && has(ra, "chain stops here (no attribution array was needed)"),
       "no array, --until attrib: reaches the no-array branch and exits 0")
rr <- run_in(FR, FR_DRIVER, fake_args("results"))
expect(rr$status == 0L && has(rr, "no attribution array to wait on; running results now") &&
         has(rr, "[fake deploy] stage align_results") && has(rr, "--until results: chain stops here."),
       "no array, --until results: runs results inline (align_results + 09) and exits 0")
expect(length(ok_calls()) == n0, "no array: no sbatch call")
rq <- run_in(FR, FR_DRIVER, fake_args("results"), extra_env = "FAKE_ARRAY=555")
new_calls <- setdiff(ok_calls(), ok_calls()[seq_len(n0)])
expect(rq$status == 0L && has(rq, "phase 'results' queued as job 4242 (afterok:555)") &&
         length(new_calls) == 1L && grepl("--dependency=afterok:555", new_calls, fixed = TRUE),
       "with an array id: results is relayed afterok on it")
src2 <- readLines(DRIVER, warn = FALSE)
caps <- grep("^\\s+(ids|id)=\\$\\(", src2, value = TRUE)
multi <- grep("collect_ids \\|\\| true\\)$", src2, value = TRUE)
expect(length(caps) >= 8L && all(grepl("\\|\\| true\\)$", caps) |
                                   grepl("^\\s+ids=\\$\\( \\{ echo \"\\$fc\"", caps)) &&
         length(multi) >= 6L,
       sprintf("every id capture (%d) ends in || true", length(caps)))

# --- [3b] --after ----------------------------------------------------------------
cat("\n[3b] --after: the first phase is queued on the given jobs\n")
ra <- plan("--deploy", DEP, "--results-config", RES, "--after", "123:456")
expect(ra$status == 0L && any(grepl("^  start   : phase '(classes|fits)' queued as a relay job, afterok:123:456$", ra$out)),
       "--after plans the first phase as a relay job on those ids")
rb <- plan("--deploy", DEP, "--results-config", RES, "--after", "12,34")
expect(rb$status != 0L, "--after refuses ids that are not colon-joined integers")
rc <- plan("--deploy", DEP, "--results-config", RES, "--after", "12", "--phase", "fits")
expect(rc$status != 0L, "--after and --phase are refused together")
src_a <- readLines(DRIVER, warn = FALSE)
expect(any(grepl('^elif \\[\\[ -n "\\$AFTER" \\]\\]; then$', src_a)) &&
         any(grepl('^  relay "\\$FROM" "\\$AFTER"$', src_a)),
       "with --after the driver relays the first phase instead of running it")

# --- [4] nothing submitted ------------------------------------------------------
cat("\n[4] nothing submitted\n")
expect(file.exists(file.path(STUB_DIR, "sbatch")), "the failing sbatch stub was in place throughout")
expect(!file.exists(MARKER),
       if (file.exists(MARKER)) paste("sbatch was called:", paste(readLines(MARKER), collapse = " | "))
       else "no sbatch call across every --plan / --dry-run above")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
