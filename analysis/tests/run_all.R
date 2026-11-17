# run_all.R — runs the test harness; any HARD failure makes the exit status 1.
#
# Usage:
#   Rscript analysis/tests/run_all.R
#   Rscript analysis/tests/run_all.R --quick --jobs 4
#   Rscript analysis/tests/run_all.R --baseline 1990_2019
#
# Details:
#   Every test runs; a hard failure makes the exit status 1 at the end, a soft
#   failure only shows in the report. Writes SUMMARY.md and timings.csv.
#   --jobs 1 (default) runs the tests in registration order with their output
#   streamed; --jobs N > 1 runs them on N worker processes, longest first by
#   the last run's timings.csv, each test's output in
#   test_reports/logs/<test>.log (a failing test's log tail is printed at the
#   end). The tests write only to their own temp dirs and report files, so the
#   order does not matter.

suppressPackageStartupMessages({
  library(optparse)
})

opt <- parse_args(OptionParser(
  option_list = list(
    make_option(
      c("-b", "--baseline"),
      default = "1990_2019",
      type = "character"
    ),
    make_option(
      c("--quick"),
      action = "store_true",
      default = FALSE,
      help = "skip slow tests (diagnostics)"
    ),
    make_option(
      c("-j", "--jobs"),
      default = 1L,
      type = "integer",
      help = "tests run in parallel; 1 = in order, output streamed [default: %default]"
    )
  )
))

report_dir <- "generated_data/test_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(report_dir, "SUMMARY.md")
timings_path <- file.path(report_dir, "timings.csv")

# Registers a test; the tests run after the last registration (section "Run").
run_one <- function(script, args = character(), hard_fail = FALSE) {
  list(script = script, args = args, hard_fail = hard_fail)
}

# Runs one registered test: output streamed, or into `log` when given.
exec_one <- function(job, log = NULL) {
  cmd <- paste(c("Rscript", job$script, job$args), collapse = " ")
  t0 <- proc.time()[["elapsed"]]
  status <- if (is.null(log)) {
    message("\n>>> ", cmd)
    system(cmd)
  } else {
    system(paste(cmd, ">", shQuote(log), "2>&1"))
  }
  list(
    script = job$script,
    args = paste(job$args, collapse = " "),
    status = status,
    hard_fail = job$hard_fail,
    elapsed = proc.time()[["elapsed"]] - t0
  )
}

# Precondition: the deployed run's fit batches (prod_run) are deploy-generated
# configs; a checkout that does not version them (the public release)
# regenerates them first, so the config-scanning gates (formula / prior goldens,
# deployed model terms) see the configs their goldens were captured on.
if (!dir.exists("analysis/configs/prod_run_occ_grouped_enso_singlelag")) {
  message(
    "run_all: generating the prod_run fit batches (deploy_prod_run.R --stage gen_batches)"
  )
  if (
    system(
      "Rscript analysis/utils/deploy_prod_run.R -C analysis/configs/deploy/prod_run.yml --stage gen_batches > /dev/null 2>&1"
    ) !=
      0
  ) {
    stop("run_all: gen_batches of prod_run failed")
  }
}

results <- list()

# Hard-fail tests --------------------------------------------------------
results$occurrence_config_merge <- run_one(
  "analysis/tests/test_occurrence_config_merge.R",
  hard_fail = TRUE
)

results$occurrence_dispatch <- run_one(
  "analysis/tests/test_occurrence_dispatch.R",
  hard_fail = TRUE
)

# Shared top-level config blocks are merged from .CONFIG_TOPLEVEL_BLOCKS
# (analysis/utils/config.R). HARD: a block registered in get_default_options() but
# missing from the registry resolves to NULL for every config that omits it,
# with no error anywhere.
results$config_block_registry <- run_one(
  "analysis/tests/test_config_block_registry.R",
  hard_fail = TRUE
)

# analysis/configs/results/manuscript.yml is the authoritative record of the
# manuscript's output decisions, and phases 09/10 both read it. HARD: if the
# shipped config stops restating every get_default_options()$results key, the
# missing decision is silently supplied from config_defaults.yml and the file
# no longer documents what produced a figure. Also gates artefact-path
# collisions between tracks/set-lists in make_results_run_file().
results$results_config_coverage <- run_one(
  "analysis/tests/test_results_config_coverage.R",
  hard_fail = TRUE
)

# Every enforce_options flag in a 04_fit_*.R worker script must have a
# .CLI_ROUTE entry (analysis/utils/config.R) or it is silently dropped by
# .route_cli_to_nested()'s intersect(). HARD for the same
# reason as config_block_registry: the failure mode is silent.
results$cli_route_completeness <- run_one(
  "analysis/tests/test_cli_route_completeness.R",
  hard_fail = TRUE
)

# Consumer-side hardcoded literals (paths.R DEFAULT_* constants and every
# %||% "literal" fallback, 04_fit_run.R's dispatcher default, and every %||%
# literal of the simulator / smoothing modules) vs get_default_options() drift.
# HARD for the same reason as cli_route_completeness: the "config computed but
# not threaded to its consumer" failure mode is silent. Its own allowlist is
# empty, same as cli_route_completeness's KNOWN_UNROUTED.
results$hardcoded_default_parity <- run_one(
  "analysis/tests/test_hardcoded_default_parity.R",
  hard_fail = TRUE
)

# One outbreak-database vintage for the data leg: every numbered script's --indir
# defaults to utils/paths.R::default_outbreak_data_dir() (outbreak_data_Nov_2026), no
# numbered script or maintained bash driver holds its own outbreak_data_<Mon>_<YYYY>
# literal, and 00/01/02 stamp the indir into prep_counts (report_cascade() unchanged
# with or without it). HARD, quick-safe: two vintages mixed by default drop
# unmatched location periods silently.
results$data_vintage_default <- run_one(
  "analysis/tests/test_data_vintage_default.R",
  hard_fail = TRUE
)

# population.lp_shapefiles (read by 02a/02b) defaults to the same vintage
# directory as default_outbreak_data_dir(), in get_default_options() and
# config_defaults.yml. HARD, quick-safe: config.R cannot source paths.R, so
# the default is a literal that would otherwise drift silently.
results$lp_shapefiles_vintage <- run_one(
  "analysis/tests/test_lp_shapefiles_vintage.R",
  hard_fail = TRUE
)

# Every artefact and config name the builders emit for analysis/configs/ fits
# NAME_MAX (255 bytes; paths.R::assert_basename_length() / check_run_file_names()).
# HARD, quick-safe (pure path algebra, ~100 s): a name over the limit cannot be
# written, and the post-fit kinds are the longest, so it fails after the MCMC.
results$basename_length <- run_one(
  "analysis/tests/test_basename_length.R",
  hard_fail = TRUE
)

# Nothing numbered, in utils/, or in bash/ may invoke analysis/scratchpads/ (the
# "a number means it runs" rule). HARD: a regression is silent --
# the HPC driver just keeps running an unnumbered script. Empty allowlist, same
# ratchet shape as cli_route_completeness's KNOWN_UNROUTED.
results$no_scratchpad_callers <- run_one(
  "analysis/tests/test_no_scratchpad_callers.R",
  hard_fail = TRUE
)

# utils/ modules never source() each other outside the documented exceptions
# (_core kernels, spatpomp_core_c.R, the non-step CLI tools); the config.R split
# modules each load alone after base.R. HARD: a hidden cross-source breaks the
# "caller sources the transitive closure" contract silently. Ratchet allowlist.
results$utils_module_purity <- run_one(
  "analysis/tests/test_utils_module_purity.R",
  hard_fail = TRUE
)

# Ownership table of the former config.R symbols across config.R / formula.R /
# formula_legacy.R / priors.R. HARD: a symbol moving or
# appearing twice changes which source() lines a caller needs.
results$config_split_symbols <- run_one(
  "analysis/tests/test_config_split_symbols.R",
  hard_fail = TRUE
)

# utils/results_derive.R (09's derivers) reads the 09 globals
# only through the results_context() environment. HARD: a bare reference would
# only fail -- or silently bind a same-named global -- when its branch runs.
results$results_derive_no_free_globals <- run_one(
  "analysis/tests/test_results_derive_no_free_globals.R",
  hard_fail = TRUE
)

# The phase-09 step registry (utils/results_derive.R .RESULTS_STEP_TABLE) that
# 09a turns into SLURM tasks and 09b runs. HARD: a kind no step owns is never
# derived on the HPC path, and a wrong wave lets two array jobs derive the same
# artefact at once.
results$results_step_registry <- run_one(
  "analysis/tests/test_results_step_registry.R",
  hard_fail = TRUE
)

# Shared best-lag / grouped-config helpers the figure scripts and drivers use
# (utils/loo_lags.R::best_lags_by_group(), paths.R::grouped_lag_config_path(),
# config.R::detect_model_from_yaml()). HARD: a wrong lag pick does not error, it
# silently pictures a different fit than the deployed/manuscript one.
results$loo_lags_helpers <- run_one(
  "analysis/tests/test_loo_lags_helpers.R",
  hard_fail = TRUE
)

# Cross-track canonical episodes (episodes.key: gid_season) -- the shared
# builder (utils/footprint.R), the gid-season occurrence history (occ_data.R),
# the episode digest and its loader/consumer guards (provenance.R), and the
# default lp_season artefact names. HARD: a drift here makes the four tracks
# disagree on what one outbreak is, which no downstream number reveals.
results$episode_key <- run_one(
  "analysis/tests/test_episode_key.R",
  hard_fail = TRUE
)

# National-zero months in the coverage mask (utils/national_zero.R,
# 01_preprocess_national_zero_months.R): the national / WHO / sub-national /
# active-outbreak / exclusion rule and build_coverage_mask()'s union (on by
# default). HARD: a wrong month counts simulated cholera against an observed
# zero that never was.
results$national_zero_months <- run_one(
  "analysis/tests/test_national_zero_months.R",
  hard_fail = TRUE
)

# National zero-fill of the occurrence panel (episodes.national_zero_fill,
# 03_prepare_occ.R --national_zero_fill): which rows national_zero_surv_rows()
# adds (never at/after an onset, never a duplicate), the "_nz" artefact token,
# the loader stamp guard and the sim-config flag. HARD: a wrong row adds
# exposure no report supports; the default panel must stay unchanged.
results$national_zero_fill <- run_one(
  "analysis/tests/test_national_zero_fill.R",
  hard_fail = TRUE
)

# One occurrence-panel spec / token / exact whole-suffix matcher:
# names, blocks and generated configs byte-identical at the defaults; a "_nz60"
# lookup never matches "_nz60-r2" / "_nz60_ae_a1c" files. HARD: a loose match
# silently deploys or simulates the wrong panel's fits.
results$occ_panel_spec <- run_one(
  "analysis/tests/test_occ_panel_spec.R",
  hard_fail = TRUE
)

# Overwrite guards of the other fixed-name 03 artefacts (siz, dur,
# mon, the three sim panels, the WHO monthly table): stamp and guard before
# every write, --force threaded, identical rebuilds pass. HARD: a data-vintage
# refresh must not replace a deployed run's prep inputs silently.
results$prep_overwrite_guards <- run_one(
  "analysis/tests/test_prep_overwrite_guards.R",
  hard_fail = TRUE
)

# Occurrence-panel masks (at_risk_start / adm1_overlap), the
# orphan ADM1 classification, discarded onsets, the TTLO level "none" and the
# D31 table. HARD: the default panel must stay identical; a wrong mask removes
# exposure or onsets the data support.
results$occ_panel_masks <- run_one(
  "analysis/tests/test_occ_panel_masks.R",
  hard_fail = TRUE
)

# Contagion coverage (stamp only): the contagion-coverage share columns never match
# the covariate pattern, the coverage definition and its D33 table. HARD: a
# share column read as a model term would enter a fit.
results$contagion_coverage_report <- run_one(
  "analysis/tests/test_contagion_coverage_report.R",
  hard_fail = TRUE
)

# deploy.episodes.key: the rollout's single episode declaration -- written into
# every generated forecast/hindcast config, pre-flighted against the prep
# artefacts, and never pruning another episode definition's fits. HARD: a
# rollout silently on the wrong definition only shows up in simulated numbers.
results$deploy_episode_key <- run_one(
  "analysis/tests/test_deploy_episode_key.R",
  hard_fail = TRUE
)

# Every analysis/configs/deploy/*.yml restates get_default_options()$deploy leaf
# for leaf and resolves with zero default fills. HARD: a key added to the
# defaults but missing from a shipped deploy YAML is filled silently.
results$deploy_config_coverage <- run_one(
  "analysis/tests/test_deploy_config_coverage.R",
  hard_fail = TRUE
)

# The deploy YAML is the single source of run identity: run-specific deploy
# defaults are "" and refused, `deploy_prod_run.R --stage resolve` prints the
# KEY=value contract bash drivers read, and --stage gen_batches plans/generates
# fit-batch dirs. HARD: a stale run inherited from a default goes unnoticed.
results$deploy_resolve <- run_one(
  "analysis/tests/test_deploy_resolve.R",
  hard_fail = TRUE
)

# Past-ENSO-event attribution (utils/make_event_attribution_configs.R, deploy
# --stage gen_events): the event registry, its launch rule against the cached
# ONI, the per-group windows and the per-event Lenssen cache checks. HARD,
# quick-safe: a wrong window or a stale cache silently attributes the wrong months.
results$event_attribution_configs <- run_one(
  "analysis/tests/test_event_attribution_configs.R",
  hard_fail = TRUE
)

# Per-country cholera-season ends (utils/season_calendar.R; the hindcast and
# event generators' --season_end by_country_from_panel): the
# min_seas -> season-end map, the gid -> country mode and tie rule, the
# hindcast 6-month and event peak rules. HARD, quick-safe: a wrong end month
# silently moves every scored or attributed window.
results$season_calendar <- run_one(
  "analysis/tests/test_season_calendar.R",
  hard_fail = TRUE
)

# Group display labels: results_core.R::group_display_meta()
# / relabel_groups() (letters by the deployed occurrence lag, ids and colours
# unchanged, identity under "registry"), the 09 group_labels step on a
# synthetic best-lag archive (lags read from the fit files, a disagreeing
# archive or manifest stops), the per-country season-end tables of S30 and the
# per-country event windows. HARD, quick-safe: a wrong letter mislabels every
# figure and table.
results$group_display <- run_one(
  "analysis/tests/test_group_display.R",
  hard_fail = TRUE
)

# The supplementary hindcast skill-vs-climatology figure (10 --figures H
# --supplementary true): hindcast_clim_skill_frame() and its
# builder (utils/figures_supp.R) on synthetic skill tables. HARD, quick-safe: a
# wrong filter plots another reference, rule or support under the same caption.
results$supp_hindcast_clim_skill <- run_one(
  "analysis/tests/test_supp_hindcast_clim_skill.R",
  hard_fail = TRUE
)

# Phase-11 supplement helpers: the WHO-AFRO subregion
# column (with_afro), the data release split into supp_data/ and github/ with
# its manifest, and the supplement table builders (D34, D35, T1-T5), on
# synthetic frames. HARD, quick-safe: a wrong release map ships a table to the
# wrong place.
results$report11_supplement_tables <- run_one(
  "analysis/tests/test_report11_supplement_tables.R",
  hard_fail = TRUE
)

# The supplement figures (10 --figures DC,SV,FA,AA): their
# phase-10 builders and frame helpers (utils/figures_supp.R) and 09's smoothing
# validation tables, on synthetic artefacts. HARD, quick-safe: a builder that
# breaks on a real artefact stops the figure run.
results$supp_missing_figs <- run_one(
  "analysis/tests/test_supp_missing_figs.R",
  hard_fail = TRUE
)


# Membership uncertainty of the deployed classes: the
# 09 report_membership step on a synthetic kl-sweep comparison -- MAP ==
# registry assertion, elpd differences = loo::loo_compare(), carried-over
# countries NA, label-free sensitivity, the M2 membership-weighted ERF from
# lag sidecars. HARD, quick-safe: a wrong p_max misstates every group's
# membership certainty.
results$results_occ_panel_tables <- run_one(
  "analysis/tests/test_results_occ_panel_tables.R",
  hard_fail = TRUE
)

results$results_membership <- run_one(
  "analysis/tests/test_results_membership.R",
  hard_fail = TRUE
)

# The fit registry derived from a pinned cell of a brms_kl class sweep
# (utils/kl_grouping.R, deploy_prod_run.R --stage gen_grouping; relay phase
# `classes` of bash/tools/run_pipeline.sh), on fake kl_fits: same bytes as the
# reference make_kl_grouping.R body, overwrite policy, fits.groups check, resolve
# before the registry exists. HARD: the registry decides every fit's groups.
results$deploy_kl_grouping <- run_one(
  "analysis/tests/test_deploy_kl_grouping.R",
  hard_fail = TRUE
)

# Grouping-registry sidecars (analysis/groupings/<source>.yml), content digests
# and label helpers (utils/config.R). HARD: the
# enso_grouping sidecar must reproduce the GROUP_* tables
# that phase-09 artefacts store as data, and partition labels become paths.
results$grouping_meta <- run_one(
  "analysis/tests/test_grouping_meta.R",
  hard_fail = TRUE
)

# filter.partition_by: registry: the partition
# (TU) registry decoupled from the fit grouping -- stamping, smoothing branch,
# exact cache resolver, cache-dir meta.yml identity, native-or-stitched
# handoff, and the consolidated pre-flight (contagion must refine the TUs).
results$smoothing_partition_registry <- run_one(
  "analysis/tests/test_smoothing_partition_registry.R",
  hard_fail = TRUE
)

# registry with partition_source == the fit grouping reproduces enso_region bit
# for bit: structural always, plus a real 2-replicate 07b run in both modes when
# the real-data inputs are present (~25 s). HARD: the new mode must not perturb the
# deployed one.
# Neighbour contagion over an arbitrary registry: nbr_lag1_r<rho>[_rn]
# _part<token> names, parse round trip, token -> registry uniqueness, prep
# columns equal to the simulator kernel. HARD: the fit -> simulator round trip
# rests on the covariate name alone.
results$nbr_partition_token <- run_one(
  "analysis/tests/test_nbr_partition_token.R",
  hard_fail = TRUE
)

# Grouped fits of ONE country (e.g. enso_nga) keep their country-level REs
# under a normal(0, 0.05) sd prior instead of having them trimmed: the trimmed
# group names (ADM1, ADM1:gid) are not mapped by the simulator / smoothing
# extractor, which zeroed those REs silently. HARD: replay correctness.
results$single_country_re <- run_one(
  "analysis/tests/test_single_country_re.R",
  hard_fail = TRUE
)

# Per-country firm-horizon censoring + forecast/hindcast manifest v2 readers:
# the one place a unit's raw output is cut per country; v1 manifests
# read as units = groups.
results$sim_censor_helper <- run_one(
  "analysis/tests/test_sim_censor_helper.R",
  hard_fail = TRUE
)

# ENSO past the CFSv2 span: neutral or AR(1) continuation; legacy_seam
# reproduces the legacy seam sampler exactly. HARD: the
# horizon month of every group's min-lag track reads the last-span ONI.
results$enso_ar1_continuation <- run_one(
  "analysis/tests/test_enso_ar1_continuation.R",
  hard_fail = TRUE
)

# ENSO analog continuation past the CFSv2 span (k-NN from the piControl
# library): neighbours, faded join, neutral / ar1 identical to the pre-analog
# sampler, leads covered, refusals, provenance. HARD.
results$enso_analog_continuation <- run_one(
  "analysis/tests/test_enso_analog_continuation.R",
  hard_fail = TRUE
)

# Skill by lead of the ENSO continuation (03_prepare_enso_continuation_library.R
# --validate; utils/enso_continuation_skill.R): fair CRPS,
# lead indexing, references vs themselves, analog = the deployed sampler path,
# output schema, settings overlay. HARD: 09 reports the continuation CRPSS of
# each group from this table.
results$enso_continuation_skill <- run_one(
  "analysis/tests/test_enso_continuation_skill.R",
  hard_fail = TRUE
)

# One ENSO phase function: classify_enso_phase() equals the five frozen
# reference code paths bit for bit; only the episode rule is "CPC"; no +/-0.5
# literal outside the constants; report_smoothing's phase rule. HARD, quick.
results$enso_phase <- run_one(
  "analysis/tests/test_enso_phase.R",
  hard_fail = TRUE
)

# One horizon-end convention: closed forward horizons,
# half-open filter windows; build_climate_panel_for_tu(end_mode=) identical to
# the frozen slice; assert_arm_launch end check; config-time panel coverage;
# enso_seam_month(). HARD, quick.
results$horizon_convention <- run_one(
  "analysis/tests/test_horizon_convention.R",
  hard_fail = TRUE
)

# The stale enso_central exclusion: defaults / shipped results records
# exclude nothing, registry guard, manifest-aware stamp normaliser
# (.RESULTS_STAMP_NORMALISE), enso_phase back-compat. HARD, quick.
results$results_exclude_groups <- run_one(
  "analysis/tests/test_results_exclude_groups.R",
  hard_fail = TRUE
)

# Counterfactual ENSO analogs: the legacy Lenssen builder's
# helpers moved verbatim (frozen copies), the library matcher, multi-winter
# sets, read winters, cache spec / digest / status, analog seam. HARD, quick.
results$enso_counterfactual_match <- run_one(
  "analysis/tests/test_enso_counterfactual_match.R",
  hard_fail = TRUE
)

# 09 side of the counterfactual analogs: report_analog_sets' D29 rows (legacy + library
# caches), the counterfactual cache identity assert, the oni step reading the
# config's cache_path. HARD, quick.
results$results_analog_sets <- run_one(
  "analysis/tests/test_results_analog_sets.R",
  hard_fail = TRUE
)

# Replicate chunking of the ENSO scenarios (--chunk_size / --n_replicates_total):
# chunks bound together equal the unchunked sample; both loaders slice. HARD.
results$enso_chunk_invariance <- run_one(
  "analysis/tests/test_enso_chunk_invariance.R",
  hard_fail = TRUE
)

# Fast paths of the results loaders (data.table grouping in
# load_sim_run_monthly(), sort-and-count in adm1_prob_stats()): identical to
# the dplyr / outer() versions they replace. HARD.
results$results_loader_fast <- run_one(
  "analysis/tests/test_results_loader_fast.R",
  hard_fail = TRUE
)

# 09 derivation speed-ups (data.table grouping / anti-join, one-call ribbon
# quantiles, one run load per job): identical to the frozen reference helpers. HARD.
results$results_derive_fast <- run_one(
  "analysis/tests/test_results_derive_fast.R",
  hard_fail = TRUE
)

# report_attrib_metrics v2's scope-level IPCC class (ipcc_scope): equal to a
# reference pd_scope_class() on synthetic pairs; v1 tables unchanged. HARD.
results$results_ipcc_scope <- run_one(
  "analysis/tests/test_results_ipcc_scope.R",
  hard_fail = TRUE
)

# Forecast handoff states among the smoothing replicates: uniform (default,
# balanced by global replicate id, chunk-safe) or loglik (SIR). HARD.
results$handoff_weighting <- run_one(
  "analysis/tests/test_handoff_weighting.R",
  hard_fail = TRUE
)

# Golden of the forward draw indices / handoff map under the
# pre-existing modes (unchanged by the joint design). HARD.
results$forward_draw_golden <- run_one(
  "analysis/tests/test_forward_draw_golden.R",
  hard_fail = TRUE
)

# Joint state-parameter draws (shared_cycle smoothing,
# param_draw_mode joint, meta.yml identity, 09 design check, event generator). HARD.
results$joint_draws <- run_one(
  "analysis/tests/test_joint_draws.R",
  hard_fail = TRUE
)

# The hindcasts' posterior-draw mode (joint)
# and replicate count, applied after patch_hindcast_config(). HARD.
results$hindcast_draw_mode <- run_one(
  "analysis/tests/test_hindcast_draw_mode.R",
  hard_fail = TRUE
)

# The kl program's cmdstanr compile is serialised across array tasks sharing
# stan_dir (guards against concurrent compiles); its .stan is written atomically. HARD.
results$kl_compile_lock <- run_one(
  "analysis/tests/test_kl_compile_lock.R",
  hard_fail = TRUE
)

# Results phase on a v2 (registry) manifest: unit runs relabelled by fit
# group and censored per country; v1 lookups and config lists unchanged; a
# stamp predating beyond_firm_horizon still reads as current.
# Golden characterization of the 09 loader with a smoothing prefix (stitched,
# smoothing-as-factual, per-country horizons): the per-smoothing-replicate
# aggregation reproduces the forecast-replicate grid exactly.
results$results_stitched_load_golden <- run_one(
  "analysis/tests/test_results_stitched_load_golden.R",
  hard_fail = TRUE
)
# Compact (format-2) monthly summary: every writer / reader pair
# sees the same long table (round trips, merged / combined files, the
# 09 loader against the stitched-load golden on compact parts).
results$sim_monthly_compact <- run_one(
  "analysis/tests/test_sim_monthly_compact.R",
  hard_fail = TRUE
)
results$results_unit_runs <- run_one(
  "analysis/tests/test_results_unit_runs.R",
  hard_fail = TRUE
)


# manuscript.yml follows the deploy YAML and says so: `--stage align_results`
# edits deploy_config + the eight deploy-derived inputs.* as TEXT (comments kept,
# re-read verified), and an explicit inputs leaf that disagrees with the deploy
# is refused unless named under results.allow_overrides. HARD: otherwise
# manuscript.yml can point at one run while reading another run's batches, or
# lose the decision record's comments silently.
results$deploy_align_results <- run_one(
  "analysis/tests/test_deploy_align_results.R",
  hard_fail = TRUE
)

# The spatPomp templates are generated (--stage gen_templates) from the
# analysis/configs/sim_templates/ skeletons + the deploy `sim:` overlay: every
# deploy YAML with templates on disk reproduces them leaf for leaf (bar run_id),
# derived leaves cannot be overlaid, and a differing template is refused
# before anything is written. HARD: a template drifting from its deploy YAML
# would smooth / forecast a run under settings nobody recorded.
results$deploy_sim_templates <- run_one(
  "analysis/tests/test_deploy_sim_templates.R",
  hard_fail = TRUE
)

# Forecast / hindcast horizon = min(best_lag_occ, best_lag_dur, best_lag_mon)
# over the DEPLOYED fits (utils/paths.R::read_deployed_lag(); _ep-gs filter,
# duplicate refusal; a deployed size fit is ignored). HARD: a wrong horizon
# forecasts months whose ENSO covariate lies beyond the observed/forecast ENSO.
results$deployed_lag_horizon <- run_one(
  "analysis/tests/test_deployed_lag_horizon.R",
  hard_fail = TRUE
)

# Size retired downstream: no size key in the defaults, and a
# simulator config setting simulator.size_fit / size_cb_meta / size /
# location_frame.mean_rate_table_path or posterior_strategy != paired is
# refused (fit-side size configs still resolve). HARD: a retired key would
# otherwise be a silent no-op.
results$config_size_retired <- run_one(
  "analysis/tests/test_config_size_retired.R",
  hard_fail = TRUE
)

# No maintained bash driver (stages/ campaigns/ tools/ slurm/) hardcodes a
# prod_run<N>; deploy-aware drivers resolve via --deploy. HARD: a hardcoded run
# default drifts silently. Empty allowlist, same ratchet as
# no_scratchpad_callers.
results$no_hardcoded_run_in_bash <- run_one(
  "analysis/tests/test_no_hardcoded_run_in_bash.R",
  hard_fail = TRUE
)

# tools/run_pipeline.sh --plan on prod_run (and temp no-event
# copies): the sim phase chains gen_events + stages/07_event_attribution.sh
# (results waits on the event arrays; --event-replicates -> --n_replicates),
# a deploy without events skips them, and the reports phase reads the chain's
# --results-config. Also stages/07_event_attribution.sh's skip-complete guard
# (temp _SUCCESS sentinels) and the attrib phase's no-array branch (fake
# drivers). Stub sbatch: nothing real is submitted. HARD,
# quick-safe: a step the relay does not chain is silently missing from a launch.
results$run_pipeline_phases <- run_one(
  "analysis/tests/test_run_pipeline_phases.R",
  hard_fail = TRUE
)

# analysis/bash/lib/cluster.sh and the local backend (lib/local/sbatch):
# cholera_sbatch's flags under a temp profile (class partition, account / qos,
# CHOLERA_RES_<KEY> > --defaults > class, the caller's options win, array log
# names, a payload's key), cholera_env's stack
# fallback, and the shim (arrays with %cap and :step, --parsable, --wrap,
# --export, the payload's #SBATCH job name / cpus, a failed task -> non-zero
# exit). Temp payloads only. HARD, quick-safe: a wrong flag sends every job of
# a class to the wrong partition; a shim that swallows a failure lets the
# local chain run on missing artefacts.
results$cluster_backend <- run_one(
  "analysis/tests/test_cluster_backend.R",
  hard_fail = TRUE
)


# campaigns/preprocess.sh and the 00-02 stage drivers --plan (stub sbatch):
# 00 -> 01 (+01b, 01_nz on the deploy's rule) and 00 -> 02 population -> 02c,
# crosswalk afterok 01 + 02c, climate afterok the crosswalk; --from / --until
# slices; utils/preprocess_paths.R names the vintage and the outputs. HARD,
# quick-safe: a wrong edge lets a step read the previous vintage's artefact.
results$preprocess_stages_plan <- run_one(
  "analysis/tests/test_preprocess_stages_plan.R",
  hard_fail = TRUE
)

# analysis/pipeline_steps.yml is the only run-order record (no executed entry
# point): every numbered script, report, stage driver and deploy stage has an
# entry, the list order is a valid run order, and PIPELINE.md's generated table
# and notebooks/pipeline_run_graph.dot are fresh. HARD: a missing step is a
# silent hole in the only place the order is written down.
results$pipeline_doc_coverage <- run_one(
  "analysis/tests/test_pipeline_doc_coverage.R",
  hard_fail = TRUE
)

# Numbered scripts open with their own filename, and every `Rscript analysis/NN*.R`
# command in numbered scripts, utils/ and maintained bash names an existing file.
# HARD: a rename can leave headers (and their usage lines) on stale names.
results$script_header_names <- run_one(
  "analysis/tests/test_script_header_names.R",
  hard_fail = TRUE
)

# Every R file's header follows the README.md "Code style" template, and
# numbered scripts use numbered sections starting with 1. Preamble.
# HARD: a new or edited file would drift from the documented layout.
results$script_header_layout <- run_one(
  "analysis/tests/test_script_header_layout.R",
  hard_fail = TRUE
)

# The R code is formatted as air formats it (air.toml; skips without air).
# HARD: an unformatted file would reformat on the next style.sh --fix.
results$style_air <- run_one(
  "analysis/tests/test_style_air.R",
  hard_fail = TRUE
)

# random_effects.re_adm2_leaf_scope (size/monthly): the ADM2-only leaf term and
# the simulator's read-back into its adm2 bucket. HARD: a mis-parsed leaf is a
# silent RE of 0 at replay.
results$re_leaf_scope <- run_one(
  "analysis/tests/test_re_leaf_scope.R",
  hard_fail = TRUE
)

# Simulator-partition aggregates shared by 09b_postprocess_results.R and
# 07d_simulator_figures_counterfactual.R (utils/results_core.R). HARD: the Z vs
# cases_month switch and the ADM1 shift statistics change numbers, not behaviour.
results$results_core_sim_aggregates <- run_one(
  "analysis/tests/test_results_core_sim_aggregates.R",
  hard_fail = TRUE
)

# Phase-11 per-draw population helpers + incidence denominators
# (utils/results_derive.R, figures_cb.R::enso_conds_named()). HARD, quick-safe.
results$results_sensitive_pop <- run_one(
  "analysis/tests/test_results_sensitive_pop.R",
  hard_fail = TRUE
)

# report_attrib_exposure (forecast vs counterfactual exposed population) and the
# helpers extracted for it, each against the code it replaced. HARD, quick-safe.
results$results_attrib_exposure <- run_one(
  "analysis/tests/test_results_attrib_exposure.R",
  hard_fail = TRUE
)

# report_attrib_metrics (attributable fraction / number, probability ratio +
# probabilities of causation, signal-to-noise): kernels by hand and the tables
# end to end on synthetic paired runs. HARD, quick-safe.
results$results_attrib_metrics <- run_one(
  "analysis/tests/test_results_attrib_metrics.R",
  hard_fail = TRUE
)

# report_event_attrib (past strong ENSO events): smoothing-as-factual
# reading through the handoff map, the per-event manifest view, the per-event
# cells and the observed-anchored AN, on synthetic runs. HARD, quick-safe.
results$results_event_attrib <- run_one(
  "analysis/tests/test_results_event_attrib.R",
  hard_fail = TRUE
)

# One-event La Nina 2020-23: enso_events_v2.yml parse, the
# per-winter sub-windows, the schema-3 manifest table and the 09 sub-window
# cells, with frozen-copy equivalence of the schema-2 path. HARD, quick-safe.
results$event_sub_windows <- run_one(
  "analysis/tests/test_event_sub_windows.R",
  hard_fail = TRUE
)

# af_posterior / event_af_posterior: the posterior interval of the
# attributable fraction from fixed-draw runs (recovery of a synthetic known
# posterior, Africa / AN definitions, design refusals) and the bootstrap switch
# (attrib_boot_b = 0). HARD, quick-safe.
results$results_af_posterior <- run_one(
  "analysis/tests/test_results_af_posterior.R",
  hard_fail = TRUE
)

# attrib_exposure_posterior: the posterior interval of the Childs
# population (recovery of known per-ADM1 AFs, the point = report_attrib_exposure's
# Childs population, guards, refusals). HARD, quick-safe.
results$results_attrib_exposure_posterior <- run_one(
  "analysis/tests/test_results_attrib_exposure_posterior.R",
  hard_fail = TRUE
)

# Attribution factual: the generator of the factual launched at
# the counterfactual handoff (make_forecast_attrib_configs.R, --stage
# gen_attrib), the attribution manifest and its results-side checks; and the
# information component (report_attrib_information: total = information +
# ENSO, pairs, guard, handoff-map refusal). HARD, quick-safe.
results$forecast_attrib_configs <- run_one(
  "analysis/tests/test_forecast_attrib_configs.R",
  hard_fail = TRUE
)
results$results_attrib_information <- run_one(
  "analysis/tests/test_results_attrib_information.R",
  hard_fail = TRUE
)
# One Africa window for the attribution steps:
# results_africa_window(), the cell table, table S28's continuation leads and
# the cells on a synthetic run with continuation months. HARD, quick-safe.
results$results_africa_window <- run_one(
  "analysis/tests/test_results_africa_window.R",
  hard_fail = TRUE
)
# The explosive-replicate guard as a results decision: the S29
# dropped-pair summary by hand, cap and modes from the context (no constant),
# group-level keys for the country ratios, stamp back-compat. HARD, quick-safe.
results$results_explosive_guard <- run_one(
  "analysis/tests/test_results_explosive_guard.R",
  hard_fail = TRUE
)
# Fit-level model specification of the deployed fits (report_model_spec, S34):
# extraction from fit-like objects built from
# model_terms_deployed.yml, the check against it, registry and paths. HARD, quick-safe.
results$results_model_spec <- run_one(
  "analysis/tests/test_results_model_spec.R",
  hard_fail = TRUE
)

# rain_enso: Figure 1d's rainfall tables -- the ADM1 series, the
# rainfall-deviation / lagged-ONI correlation with within-group agreement and
# the rainy seasons on a synthetic panel, plus (when the inputs are present)
# the reference CSVs reproduced. HARD, quick-safe.
results$results_rain_enso <- run_one(
  "analysis/tests/test_results_rain_enso.R",
  hard_fail = TRUE
)

# Climate column plumbing: dev_sd_* never rescaled, .clim_dev_wide(), the opt-in
# get_stats_gid(). HARD, quick-safe.
results$climate_column_plumbing <- run_one(
  "analysis/tests/test_climate_column_plumbing.R",
  hard_fail = TRUE
)

# seasonality: supplementary fig. SE -- the three tracks' observed
# monthly rates (admin-level rule, gap months, sparse flag, Poisson interval),
# the monthly model's harmonic curve, the peak month and its lag from the
# rainy season, and the 4-row figure on synthetic frames. HARD, quick-safe.
results$results_seasonality <- run_one(
  "analysis/tests/test_results_seasonality.R",
  hard_fail = TRUE
)

# report_who_enso_ratio: supplementary fig. WR -- the WHO annual
# cases (absence NA, epiweek years, censoring), ENSO states of calendar years,
# the observed and modelled ratios of means (pseudo-count, complete neutral
# years, handoff-paired spread), the tables and the 2-panel figure on
# synthetic inputs. HARD, quick-safe.
results$results_who_enso_ratio <- run_one(
  "analysis/tests/test_results_who_enso_ratio.R",
  hard_fail = TRUE
)

# results_require()'s opt-in RESULTS_ALLOW_STALE_INPUTS override: warns instead
# of stopping on changed inputs, never on a decision mismatch. HARD, quick-safe.
results$results_require_stale <- run_one(
  "analysis/tests/test_results_require_stale.R",
  hard_fail = TRUE
)

# Supplementary-table producers: fit_diagnostics() on a fake kl_fit
# and, when present locally, a saved deployed brmsfit (no sampling); the kl sweep
# helpers on Stan-free fake fits; and the registered sweep batch resolving leaf
# for leaf as the hand-written configs (so it reuses their fits). HARD, quick-safe.
results$fit_diagnostics <- run_one(
  "analysis/tests/test_fit_diagnostics.R",
  hard_fail = TRUE
)
# Convergence gate of deploy_prod_run.R --stage deploy_fits / gen_grouping:
# synthetic diagnostics, injected read / diagnostics
# functions, refuse-before-prune on a mocked deploy_collect_fits(). HARD, no
# sampling.
results$deploy_convergence_gate <- run_one(
  "analysis/tests/test_deploy_convergence_gate.R",
  hard_fail = TRUE
)
# report_convergence v2 (D27 / S8): each deployed fit's deploy-gate verdict from
# the stored convergence tables equals fit_diagnostics_gate() on the fit; the
# deploy's convergence block (off / warn / refuse, allow). HARD, quick.
results$results_convergence_gate <- run_one(
  "analysis/tests/test_results_convergence_gate.R",
  hard_fail = TRUE
)
results$kl_sweep_helpers <- run_one(
  "analysis/tests/test_kl_sweep_helpers.R",
  hard_fail = TRUE
)

# prod_run fit batches: 115 configs per track (23 lags x 5 groups of
# enso_grouping_prod_run), duration harmonic2 + prior_mso_band, distinct fit
# names, the 8-cell kl class sweep, the deferral while a registry is absent,
# prod_run.yml vs the registry.
# HARD, quick-safe: a wrong leaf refits the wrong model on 345 fits.
results$prod_run_batches <- run_one(
  "analysis/tests/test_prod_run_batches.R",
  hard_fail = TRUE
)

# LOO ranking caches live under generated_data/loo_compare/<batch>/
# and every reader goes through utils/loo_lags.R's one locator. HARD because the
# failure mode is silent: a reader that re-implements the glob against an
# in-config-dir location just finds nothing (06 skips its elpd panel without error).
results$loo_compare_paths <- run_one(
  "analysis/tests/test_loo_compare_paths.R",
  hard_fail = TRUE
)

# Pareto-k of each lag's PSIS-LOO
# (utils/loo_lags.R::pareto_k_summary() / ranking_pareto()) and the per-(batch,
# group) lag sidecar of 05_fit_postprocess_loo_lags.R --lag_sidecar (pointwise
# elpd, Pareto-k, per-lag crossbasis draws, row-key digests). HARD, quick, no
# MCMC: a misaligned sidecar would silently stack the wrong rows across lags.
# Parts (b)/(c): the deploy floor, report_lags v2 rows, paired / k-robust elpd
# differences, stacking and pseudo-BMA+ weights, the common ERF grid and the
# systematic lag mixture, and the 09 lag_averaging step on synthetic sidecars
# (skipped with a message when none exists).
results$lag_averaging <- run_one(
  "analysis/tests/test_lag_averaging.R",
  hard_fail = TRUE
)


results$response_invariant <- run_one(
  "analysis/tests/test_response_invariant.R",
  hard_fail = TRUE
)

results$simulator_kernel <- run_one(
  "analysis/tests/test_simulator_kernel.R",
  hard_fail = TRUE
)

# Golden characterization of simulate_one_tu() and its step helpers
# (evaluate_lp, update_lags, advance_ttlo, compute_phi_vec, basehaz_dur_eval,
# monthly_log_mu) over four synthetic worlds that reach every deployed kernel
# branch. HARD: the kernel speed-up rewrites these internals and must stay
# bit-identical at a fixed seed (same draws, same order).
results$simulator_kernel_golden <- run_one(
  "analysis/tests/test_simulator_kernel_golden.R",
  hard_fail = TRUE
)

results$mean_rate_table <- run_one(
  "analysis/tests/test_mean_rate_table.R",
  hard_fail = TRUE
)

# The four-track generalisation of the above (analysis/utils/anchor.R). HARD for
# the same reason: every failure mode is a silent rescale of mu, never an error.
# A wrong anchor is an offset, so it multiplies every prediction of that track;
# a zero anchor is offset(log(0)) = -Inf, i.e. mu = 0 for a whole country.
results$base_rate_table <- run_one(
  "analysis/tests/test_base_rate_table.R",
  hard_fail = TRUE
)

# Duration-track data prep (analysis/utils/dur_data.R). HARD: every failure mode
# here is silent and rescales the cessation hazard rather than erroring — an
# off-by-one in the partial-month exposure, an `ended` flag on a censored
# episode, or a censoring rule wired to the wrong series all produce a panel
# that fits perfectly well and is wrong. The synthetic half needs no data; the
# TZA-2019 corpus gate self-skips when the prepared panel is absent.
results$duration_censoring <- run_one(
  "analysis/tests/test_duration_censoring.R",
  hard_fail = TRUE
)

# `duration:` block merge + formula construction. HARD: duration is a MODEL
# block (.fill_model_block), so a missing model-detection arm resolves every
# unspecified key to NULL and the fit then runs with basehaz_duration = NULL --
# no baseline hazard at all, and no error. The formula half fails the same way:
# a dropped mso_band gives a constant hazard, a dropped offset turns a per-day
# hazard into a per-row count. Both fit perfectly well and are wrong.
results$duration_config_merge <- run_one(
  "analysis/tests/test_duration_config_merge.R",
  hard_fail = TRUE
)

# make_duration_run_file() collision gate. HARD: every fit script early-exits when
# its output path already exists, so two configs that differ in a load-bearing
# field but resolve to one filename neither overwrite nor error -- the second
# run reports success and leaves the FIRST one's fit on disk under its name.
# This has already happened on the size track (see .size_estimand_path_token()).
results$duration_filename_roundtrip <- run_one(
  "analysis/tests/test_duration_filename_roundtrip.R",
  hard_fail = TRUE
)

# Phase-5 size<->duration dependence check (analysis/utils/dependence.R /
# 05_fit_postprocess_dependence.R). HARD: a wrong PIT transform doesn't error,
# it silently feeds a plausible-but-wrong Kendall's tau into the
# pre-registered decision rule for whether the simulator needs to
# couple size and duration draws. Fully synthetic (known discrete-hazard /
# NB generative process, known injected copula tau) -- no data/ access, no
# fitted model.
results$duration_dependence <- run_one(
  "analysis/tests/test_duration_dependence.R",
  hard_fail = TRUE
)

# Duration-aware simulator kernel (analysis/utils/simulator_core.R):
# extract_duration_draws()/basehaz_dur_eval()/mso_band helpers, and the
# has_dur = TRUE branch of simulate_one_tu() -- onset gating, active/mso
# bookkeeping, cessation, end-of-horizon censoring, the duty-cycle fix in
# update_lags()'s input, and the log(size) covariate hook. HARD: this is
# new modelling logic with no other test coverage; a silent regression here
# produces plausible-looking simulated outbreak durations that are simply
# wrong (e.g. re-onset every month, or contagion lags silently switching off
# the month after onset).
results$duration_kernel <- run_one(
  "analysis/tests/test_duration_kernel.R",
  hard_fail = TRUE
)

# Duration-aware simulator optionality / backward compatibility. HARD: every
# run without duration draws (has_dur = FALSE / draws_dur = NULL / omitted)
# must stay byte-identical to the onset-only kernel, including the
# exact RNG draw sequence. A regression here would silently change simulator
# output for every existing config that has never heard of the duration track.
results$duration_optional_backcompat <- run_one(
  "analysis/tests/test_duration_optional_backcompat.R",
  hard_fail = TRUE
)

# Smoothing handoff-cache schema versioning (analysis/utils/spatpomp_core.R:
# sm_cache_schema_token()/check_smoothing_cache_schema(), filter_traj_to_warmup()'s
# active/mso guard). HARD: without this, a stale (9-state) smoothing
# cache does not error at load time -- filter_traj_to_warmup()'s named-vector
# active/mso lookup silently returns NA, seeding every forecast location as
# inactive with no warning. Pure function tests, no data/ access.
results$smoothing_cache_schema <- run_one(
  "analysis/tests/test_smoothing_cache_schema.R",
  hard_fail = TRUE
)

# build_adm1_dist_km() cache coverage: a country requested at build time but with
# no ADM1 rows is covered via the sidecar's iso3 (only a sidecar written with
# gid_to_adm1() keying), so the simulator never rebuilds from the geopackage (needs
# `sf`, absent on some clusters). HARD: an uncoverable country made every enso_west
# smoothing task fail on the cluster.
results$adm1_dist_cache <- run_one(
  "analysis/tests/test_adm1_dist_cache.R",
  hard_fail = TRUE
)

# ADM1 key from a gid (gid.R::gid_to_adm1() + its inlined copies): GADM 4.1's
# dot-less Ghana gids ("GHA5_2", "GHA14_13_2") map to GHA_<n>, and an
# unparseable gid is an error. HARD: the NA ADM1 dropped Ghana from every fit.
results$gid_adm1 <- run_one(
  "analysis/tests/test_gid_adm1.R",
  hard_fail = TRUE
)

# Region-level (enso_region) bpfilter smoothing partition: the shared
# partition_frame_column()/assert_partition_column() helpers (spatpomp_core.R),
# smoothing_apply_partition()'s enso_region branch (part=group, tu left intact),
# per-partition frames spanning multiple member countries, and the forward
# direct-vs-assemble selector probe in load_smoothing_distribution() (07b). HARD:
# a silently-missing enso_region branch falls through to the TU default (one
# partition per HMM TU, no cross-country coupling) with no error; a wrong
# direct/assemble probe either re-assembles a native region cache or loads a
# per-country cache under a region key. Pure — no fit / no data / no network.
results$smoothing_partition_region <- run_one(
  "analysis/tests/test_smoothing_partition_region.R",
  hard_fail = TRUE
)

# Population-denominator QC gates (analysis/utils/population_qc.R).
# Synthetic fixtures only — passes on
# both a pre-patch corpus (provenance columns absent) and a post-patch one.
results$population_gates <- run_one(
  "analysis/tests/test_population_gates.R",
  hard_fail = TRUE
)

# LP x gid crosswalk algebra (analysis/utils/population_core.R) with a stubbed
# extractor. HARD: the crosswalk is what stops a composite LP dumping all its
# member states' population onto one arbitrary gid, and every failure mode here
# (silently dropped GEOMETRYCOLLECTIONs, NA pieces, unbound cross-border LPs)
# produces plausible-looking numbers rather than an error. No raster, no
# network, no data/ access, no OutbreakExtractR install.
results$population_crosswalk <- run_one(
  "analysis/tests/test_population_crosswalk.R",
  hard_fail = TRUE
)

# compute_gid_cover() in isolation from geometry -- composite detection and
# cover_threshold boundary behaviour. Complementary to population_crosswalk
# above, which drives the same derivation through the full geometry pipeline.
results$composite_detection <- run_one(
  "analysis/tests/test_composite_detection.R",
  hard_fail = TRUE
)

# build_lp_observation_operator() / expected_lp_from_gid_rates() (analysis/utils.R) --
# the coverage-weighted operator that maps simulated gid cases into LP space
# for validation (Stage 6). HARD: a wrong operator silently mis-attributes
# simulated cases rather than erroring.
results$lp_observation_operator <- run_one(
  "analysis/tests/test_lp_observation_operator.R",
  hard_fail = TRUE
)

# ---------------------------------------------------------------------------
# Population crosswalk WIRING (as distinct from the crosswalk algebra above).
# The four gates above test the artefact; these three test what the pipeline
# does with it. All three are HARD, and for the same reason: every failure mode
# they cover produces numbers that look entirely reasonable.
# ---------------------------------------------------------------------------

# apportion_composite_cases() -- splitting a composite LP's cases across its
# member gids by population share w = H/P_l. HARD: the failure is CONSERVATION.
# If the split does not sum back to the original, cases are created or
# destroyed silently, and the totals stay plausible either way. Includes the
# apportion -> collapse -> reconcile composition the pipeline actually runs,
# which neither helper's own tests cover.
results$case_apportionment <- run_one(
  "analysis/tests/test_case_apportionment.R",
  hard_fail = TRUE
)

# lp_composite_status() / collapse_composite_lps() -- the REDUCING counterpart to
# apportion_composite_cases(), selected by --composite_policy. HARD, and for a
# different reason than conservation: the two coexisting definitions of
# "composite" (DOMINANCE gid_cover < 0.90 vs MULTIPLICITY n_members > 1) select
# 22% of LPs and 97% of rows respectively, so keying `drop` on the wrong one
# discards nearly the whole corpus without erroring. Also pins that collapse
# scales NOTHING and denominates by P_l, so the attack rate is unchanged.
results$composite_policy <- run_one(
  "analysis/tests/test_composite_policy.R",
  hard_fail = TRUE
)

# reconcile_footprint_pop() -- the pop = min(sum_l H, Q_g) denominator replacing
# `pop = sum(pop, na.rm = TRUE)` (03_prepare_siz.R). HARD: the sum rule
# aggregated a STOCK with a flow rule; it doubled NGA_19_1 to 209M against
# Nigeria's ~220M and returned pop = 0 for the empty group COD_4_6_1. A wrong
# denominator does not error, it just rescales every rate under it.
results$footprint_pop_reconciliation <- run_one(
  "analysis/tests/test_footprint_pop_reconciliation.R",
  hard_fail = TRUE
)

# pop_provenance_stamp() / assert_pop_provenance() -- the mean_rate ordering guard.
# HARD: the mean_rate table is suspected_cases / pop, so changing pop
# invalidates every cell of it, and the simulator reads that table as a frozen
# sidecar (mu = pop * mean_rate * exp(LP), the two factors multiply). New
# populations paired with old rates is invisible to every other check here --
# both artefacts remain individually plausible and only their provenance
# distinguishes a correct pairing from an incorrect one.
results$pop_provenance_guard <- run_one(
  "analysis/tests/test_pop_provenance_guard.R",
  hard_fail = TRUE
)

# RMD duration-ratio estimand (analysis/utils/dur_post.R:
# duration_ratio_draws()/duration_ratio_summary()). HARD: the manuscript "XX times
# longer" claim is a ratio of restricted mean durations, NOT 1/HR, and the two
# agree only in the K -> inf limit. A regression that quietly returned 1/HR at a
# finite horizon would produce a plausible-but-wrong duration ratio with no
# error. Pins the constant-hazard closed form: finite-K RMD exact, convergence to
# 1/HR only as tail mass vanishes, disagreement with 1/HR at a short horizon, and
# the honest restricted-mean warning. Fully synthetic (stub fit), no data/ access.
results$duration_ratio_estimand <- run_one(
  "analysis/tests/test_duration_ratio_estimand.R",
  hard_fail = TRUE
)

results$simulator_validation <- run_one(
  "analysis/tests/test_simulator_validation.R",
  hard_fail = TRUE
)

# Golden characterization of every simulator validation metric (run_validation()
# end to end under both total-case bases, partition aggregation/combination,
# compute_hindcast_series()) on one synthetic world. HARD: guards
# the metrics in utils/sim_metrics.R and its CRPS/geography dedup -- a metric
# that drifts by a rounding error fails here. Regenerate only on a deliberate
# metric change: --capture.
results$sim_metrics_golden <- run_one(
  "analysis/tests/test_sim_metrics_golden.R",
  hard_fail = TRUE
)

# Golden characterization of resolve_gid_population() (all four tiers) and
# censor_sim_by_country() on one synthetic fixture. HARD: guards the vectorised
# T2 rescale and the per-unique-gid ISO3 of the hindcast-postprocess speed-up --
# a population or a censored row that drifts fails here. Regenerate only on a
# deliberate change: --capture.
results$pop_censor_golden <- run_one(
  "analysis/tests/test_pop_censor_golden.R",
  hard_fail = TRUE
)

# Hindcast manifest iteration shared by 07c --hindcast and 07d --hindcast
# (utils/sim_io.R::for_each_hindcast_run()). HARD: a miscounted skip or a
# silently unscored run only shows up as a thinner hindcast table on the HPC.
results$hindcast_iterator <- run_one(
  "analysis/tests/test_hindcast_iterator.R",
  hard_fail = TRUE
)

# Hindcast forecast-skill scoring (utils/hindcast_skill.R, derived by 09). HARD:
# the primitives pin the CRPS estimators, the leave-one-window-out reference,
# the case supports and the episode-key / PIT-support guards on synthetic data;
# the equivalence test pins the 09 artefacts to the notebook's reference
# numbers (SKIPs when the artefacts are absent).
results$hindcast_skill_primitives <- run_one(
  "analysis/tests/test_hindcast_skill_primitives.R",
  hard_fail = TRUE
)

# Golden characterization of build_spatpomp_tu_c() (the only spatPomp builder) on the
# dunit fixture (attributes, data, params, covariate table, generated C).
# HARD: guards the builders' shared helper blocks.
results$spatpomp_builder_golden <- run_one(
  "analysis/tests/test_spatpomp_builder_golden.R",
  hard_fail = TRUE
)

# make_batch_configs.R: deterministic, idempotent batch output and the
# backend-sibling / simulator / dlnm block stripping. Hard-fail gates.
results$smoothing_aggregates_golden <- run_one(
  "analysis/tests/test_smoothing_aggregates_golden.R",
  hard_fail = TRUE
)

# 08c's observed monthly-cases series keeps one admin level per (gid, month),
# the filter's footprint rule (smoothing_obs_cases_mon()). HARD: summing every
# level double-counts nested reports on the observed side of the cases ribbons.
results$smoothing_obs_cases_footprint <- run_one(
  "analysis/tests/test_smoothing_obs_cases_footprint.R",
  hard_fail = TRUE
)

results$size_batch_config_determinism <- run_one(
  "analysis/tests/test_size_batch_config_determinism.R",
  hard_fail = TRUE
)
results$batch_config_block_stripping <- run_one(
  "analysis/tests/test_batch_config_block_stripping.R",
  hard_fail = TRUE
)

results$simulator_dlnm <- run_one(
  "analysis/tests/test_simulator_dlnm.R",
  hard_fail = TRUE
)

results$simulator_enso_cb_grouping <- run_one(
  "analysis/tests/test_simulator_enso_cb_grouping.R",
  hard_fail = TRUE
)

results$offset_wellformed <- run_one(
  "analysis/tests/test_offset_wellformed.R",
  args = sprintf("--baseline %s", opt$baseline),
  hard_fail = TRUE
)

# Formula-builder golden file (analysis/utils/formula.R). HARD: this is the
# byte-identity guard on the shared-helper consolidation across all four
# tracks (occ/size/dur/monthly) -- a regression here means a config resolves
# to a formula that has silently changed without anyone touching that track's
# code.
results$formula_builder_golden <- run_one(
  "analysis/tests/test_formula_builder_golden.R",
  hard_fail = TRUE
)

# Deployed model terms: every config of the deployed fit
# batches (results.deploy_config -> fits.*_dir) builds one term set per track,
# equal to analysis/model_terms_deployed.yml. HARD, quick-safe: a deploy that
# turns a covariate on or off would otherwise leave CLAUDE.md and the methods
# describing a model that is not the fitted one.
results$deployed_model_terms <- run_one(
  "analysis/tests/test_deployed_model_terms.R",
  hard_fail = TRUE
)

# build_size_priors() prior-class assembly (has_slopes 'b' gating, phi_re
# 'shape'-dpar priors, custom prior string round-tripping). HARD for the same
# reason as the formula gates above: a wrong has_slopes/phi_re branch does not
# error, brms just silently fits with (or without) a class of priors it
# should not have.
results$size_priors <- run_one(
  "analysis/tests/test_size_priors.R",
  hard_fail = TRUE
)

# duration prior_mso_band: own prior on the episode-age band effects, one row
# per occupied non-reference band, and its fit-name token. HARD: a misnamed coef
# is rejected by brms only at fit time, on the cluster.
results$dur_mso_band_prior <- run_one(
  "analysis/tests/test_dur_mso_band_prior.R",
  hard_fail = TRUE
)

# Prior-builder golden file (analysis/utils/priors.R). HARD: the byte-identity
# guard on the four-track prior-construction consolidation (build_occ_priors/
# build_duration_priors/build_size_priors/build_monthly_priors sharing
# .PRIOR_DEFAULTS + private helpers) -- a regression here means a config
# resolves to brms priors that have silently changed. Also exercises the
# brms_cb vs brms prior_* split (build_size_priors()/build_monthly_priors()
# read the brms_cb sub-block when called from a _cb fit script).
results$prior_builder_golden <- run_one(
  "analysis/tests/test_prior_builder_golden.R",
  hard_fail = TRUE
)

# build_cb_block() (analysis/utils/dlnm.R) byte-identity vs the frozen
# per-track inline DLNM crossbasis loops. HARD: a silent divergence here changes
# every _cb backend's fixed/RE terms across all four tracks at once.
results$cb_block_extraction <- run_one(
  "analysis/tests/test_cb_block_extraction.R",
  hard_fail = TRUE
)

# brms_kl building blocks (dlnm `lags:` basis, kernel formula, Stan patch
# landmarks, stanc syntax, data names, pinned program golden, adjacency /
# prior-partition helpers, artefact names, config refusals). No sampling.
results$kl_stan_patch <- run_one(
  "analysis/tests/test_kl_stan_patch.R",
  hard_fail = TRUE
)

# Monthly-track (mon) make_monthly_run_file() collision gate -- sibling of
# duration_filename_roundtrip above, same failure mode (silent second-run
# success that leaves the FIRST config's fit on disk under the second one's
# name).
results$mon_filename_roundtrip <- run_one(
  "analysis/tests/test_mon_filename_roundtrip.R",
  hard_fail = TRUE
)

# Monthly-track formula wiring: cases_month LHS, pop/log_offset offsets,
# intercept RE nesting, phi/shape sub-formula, cyclic+seasonal-RE-slope
# guard, and CB slope-RE assembly (mirrors 04_fit_mon_brms_cb.R inline logic).
results$mon_formula <- run_one(
  "analysis/tests/test_mon_formula.R",
  hard_fail = TRUE
)

# Monthly-track load_monthly_data() well-formedness: pop > 0 (log(pop)
# finite), cases_gap rows never dropped, row count matches the panel's own
# country/admin_level filter. HARD for the same reason as offset_wellformed
# above -- every failure mode here is silent and produces a fit that runs to
# completion on a subtly wrong panel.
results$mon_offset_wellformed <- run_one(
  "analysis/tests/test_mon_offset_wellformed.R",
  args = sprintf("--baseline %s", opt$baseline),
  hard_fail = TRUE
)

# Monthly-track (mon) SIMULATOR-side offset contract: extract_monthly_draws()
# unconditionally forces offset_logoffset = TRUE (there is no offset_meanrate
# flag) regardless of what the inner extractor returns, and simulate_one_tu()'s mon
# mu block actually honours that (mean_rate never read; pop/log_offset do
# move mu). Distinct from mon_offset_wellformed above, which gates the
# fit-side load_monthly_data() panel, not the simulator kernel. HARD for the
# same reason: a wrong offset here is a silent multiplicative bias with no
# visible failure mode.
results$monthly_offset_wellformed <- run_one(
  "analysis/tests/test_monthly_offset_wellformed.R",
  hard_fail = TRUE
)

# Monthly-track (mon) config-validation guards: unknown `measurement` string,
# mon-without-duration hard-require (spatPomp C generator), the
# retired occurrence_negbin* names being rejected, monthly.prev_cases: TRUE.
# HARD -- an unguarded predicate chain here would silently degrade a mis-wired
# mon config to the Bernoulli path instead of erroring.
results$mon_measurement_enum <- run_one(
  "analysis/tests/test_mon_measurement_enum.R",
  hard_fail = TRUE
)

# Monthly-track (mon) P0-P4 denominator/elasticity apparatus: estimand
# resolution + validation, resolve_size_pop_column()'s asymmetric fallback (gid
# must HARD-ERROR rather than silently hand back a footprint population),
# .add_log_pop_centred()'s within-level nesting invariant, the `mon_estimand`
# stamp round-trip into .recover_log_pop_centres(), and the _estnd- filename
# separation. HARD: this is the apparatus that fixes the CMR ~76x monthly
# overestimate, and every one of its failure modes is a silent multiplicative
# rescale -- a fit that runs to completion and a forecast with no diagnostic
# that could detect the wrong population denominator or the wrong P4 pivot.
results$mon_pop_elasticity <- run_one(
  "analysis/tests/test_mon_pop_elasticity.R",
  hard_fail = TRUE
)

# compute_phi_vec(): the NegBin dispersion broadcast (constant and nested
# region/country/adm1 shape REs, grouping-mode re-rooting, unseen levels) used
# by the monthly track, under a track-neutral name. HARD: a wrong phi silently skews
# every mon case draw.
results$nb_phi_re_simulator <- run_one(
  "analysis/tests/test_nb_phi_re_simulator.R",
  hard_fail = TRUE
)

# Detection-error (`*_state_detect`) dunit family: enum wiring, the frozen
# golden emitted-C strings for the 5 state-conditioning dunit variants left
# (no NegBin size family), builder dispatch
# (C builder), delta_fn/delta_fp config plumbing + .check_deltas()
# validation, and a bpfilter smoke run. HARD -- a moved golden string means the
# emitted C changed.
results$spatpomp_dunit_detect <- run_one(
  "analysis/tests/test_spatpomp_dunit_detect.R",
  hard_fail = TRUE
)

# The normalized inactive-state case term of the `_actn` measurements
# (mass 1 at mu = 0 for every phi, compiled from the emitted C) and the frozen
# `_act` C golden (side-by-side rule). HARD -- a moved `_act` golden means the
# deployed measurement changed.
results$spatpomp_case_norm <- run_one(
  "analysis/tests/test_spatpomp_case_norm.R",
  hard_fail = TRUE
)

# Monthly-track (mon) optionality: has_mon = FALSE (draws_mon NULL/omitted)
# reproduces pre-mon simulate_one_tu() output byte-for-byte, mu_mon/
# cases_month columns are strictly opt-in, and a plain-intercept draws_mon
# bundle (no cb_meta_mon) needs no `dlnm` package. Sibling of
# duration_optional_backcompat one phase earlier.
results$monthly_optional_backcompat <- run_one(
  "analysis/tests/test_monthly_optional_backcompat.R",
  hard_fail = TRUE
)

# Monthly-track (mon) kernel machinery: build_mso_basis() bit-identity vs. a
# direct dlnm::onebasis() re-evaluation, .inject_mso()'s out-of-range clamp
# (positive control against the documented CFSv2-style unbounded-cubic
# extrapolation failure), and simulate_one_tu()'s active-gated mu_mon/
# cases_month emission (populated iff active == 1). HARD -- an unclamped or
# mis-gated mon kernel would silently corrupt simulated case counts.
results$monthly_kernel <- run_one(
  "analysis/tests/test_monthly_kernel.R",
  hard_fail = TRUE
)

# Monthly incidence cap (simulator.incidence_cap, default 0.15): monthly cases
# per unit <= cap x pop on the mon mean and the drawn count, in the v1 /
# forecast kernel, the spatPomp C filter kernel and the smoothing
# post-hoc draw; NULL/Inf reproduce the uncapped kernels exactly. HARD: the
# failure mode is silent (a smoothing can draw most of a unit's population in
# one month).
results$incidence_cap <- run_one(
  "analysis/tests/test_incidence_cap.R",
  hard_fail = TRUE
)

# Onset eligibility (simulator.onset_eligibility): the season
# gate (no second onset in a unit-season) in the v1 and C kernels, the seas_k
# covariate, the defaulted-min_seas refusal and the cache identity. HARD: a
# wrong gate silently changes every simulated onset under season.
results$onset_season_gate <- run_one(
  "analysis/tests/test_onset_season_gate.R",
  hard_fail = TRUE
)

# Cessation-month exposure (simulator.cessation_exposure): the
# default `full` is byte-identical to a pre-change golden in the v1 and C
# kernels; `half` halves exactly the ceased unit-months' mean, the C kernel's
# RNG stream unchanged; cache identity. HARD: a wrong order or a misplaced
# halving silently biases every simulated case count.
results$cessation_exposure <- run_one(
  "analysis/tests/test_cessation_exposure.R",
  hard_fail = TRUE
)

# Golden for the shared monthly-mean / monthly-draw helpers (monthly_log_mu(),
# draw_monthly_cases(), simulator_core.R): v1 grouped/ungrouped and the
# smoothing post-hoc draw on C-kernel bpfilter trajectories, all uncapped, captured
# before the four inline copies were routed through them. HARD: a drift is a
# silent change of every simulated case count.
results$monthly_mean_golden <- run_one(
  "analysis/tests/test_monthly_mean_golden.R",
  hard_fail = TRUE
)

# Zero-filled climate guard (simulator_core.R::
# .assert_no_zero_filled_climate()): a fit with a local-climate term must not
# be replayed on the panels' zero-filled post-DB months. HARD, quick-safe: the
# failure is a forecast at rainfall / temperature 0 with no error.
results$climate_zero_fill_guard <- run_one(
  "analysis/tests/test_climate_zero_fill_guard.R",
  hard_fail = TRUE
)

# Shared infrastructure: base.R::content_digest() (and
# registry_digest() delegating to it, values pinned), the smoothing-cache
# identity helpers, manifest_add_provenance() and the .TTLO_* constants. HARD,
# quick-safe: any drift here changes a stamped
# identity, a cache meta, a manifest or a panel without an error.
results$t2_shared_infra <- run_one(
  "analysis/tests/test_t2_shared_infra.R",
  hard_fail = TRUE
)

# Population year of the simulator / smoothing: HARD,
# quick-safe (synthetic; short bpfilter runs, no MCMC). The default 2020L is
# the reference frame; a mode with a constant population reproduces every
# kernel's output exactly; smoothing threading; cache / deploy identity.
results$sim_pop_year <- run_one(
  "analysis/tests/test_sim_pop_year.R",
  hard_fail = TRUE
)

# Compact bpfilter trajectory store (simulator.filter.traj_storage,
# spatpomp_core.R::.bpfilter_traj_compact()): HARD, quick-safe (synthetic,
# Np <= 300, no MCMC). "full" == the pre-change spatPomp::bpfilter golden;
# "pruned" identical() to "full"; "pruned_float32" exact on integer states and
# within 1e-6 elsewhere, logLik identical. A drift changes every smoothing
# trajectory without an error.
results$bpfilter_traj_compact <- run_one(
  "analysis/tests/test_bpfilter_traj_compact.R",
  hard_fail = TRUE
)

# Per-country WHO term scale: HARD, quick-safe (synthetic).
# Estimator, DB-overlap drop, digest, who_vs_smoothing scale fix, cache identity.
# The Csnippet golden of the fixed-scale build is in test_spatpomp_who_obs.R W8.
results$who_scale_estimate <- run_one(
  "analysis/tests/test_who_scale_estimate.R",
  hard_fail = TRUE
)

# TTLO level "none" in the simulator: HARD,
# quick-safe (synthetic). Default bins / chains identical to the code without
# the level; with the level the v1 and C kernels seed and score it alike.
results$ttlo_none_level <- run_one(
  "analysis/tests/test_ttlo_none_level.R",
  hard_fail = TRUE
)

# The duration policy stamp's episode-span fields
# and guard counts, 03_prepare_dur.R --restamp on a synthetic panel
# in a temp dir, loo_lags.R::check_occ_floor_offset(), and no caller of
# guard_outbreak_episodes() outside scratchpads. HARD, quick-safe (synthetic
# only, ~4 s): a wrong span misstates the lag floor's reference silently.
results$duration_span_stamp <- run_one(
  "analysis/tests/test_duration_span_stamp.R",
  hard_fail = TRUE
)

# enso.lag_months and re_adm1_enso_lag5_slope
# are inert on the deployed fit batches (results.deploy_config -> fits.*_dir;
# 12 configs per batch, --all for all 345): lag 3 / 5 / removed and the slope
# key dropped give identical formula terms, priors and run-file names. HARD,
# quick-safe (~4 s): the prod_run batches strip both keys on this premise.
results$enso_lag_months_inert <- run_one(
  "analysis/tests/test_enso_lag_months_inert.R",
  hard_fail = TRUE
)

if (!opt$quick) {
  # Synthetic world (00_ingest_synthetic.R): the data phase of
  # tools/run_synthetic_pipeline.sh in a temp sandbox, checked with the
  # pipeline's own readers (no real country, r3 digest, one episode definition,
  # panels complete, offline). ~2 min, no fit. SOFT: needs taxdat / terra / arrow.
  results$synthetic_world <- run_one(
    "analysis/tests/test_synthetic_world.R",
    hard_fail = FALSE
  )
  # brms_kl: the patched program has brms's own log density at K = 1
  # (deterministic; compiles two programs, ~4 min first time, cached after).
  results$kl_lp_equivalence <- run_one(
    "analysis/tests/test_kl_lp_equivalence.R",
    hard_fail = TRUE
  )
  # brms_kl MCMC tests (minutes; SOFT -- MCMC noise): synthetic-truth recovery
  # (single kernel, latent classes, prior partition, ICAR/DAGAR membership) and
  # the end-to-end smoke on the real panel (skips without the prep artefact).
  results$kl_synthetic_recovery <- run_one(
    "analysis/tests/test_kl_synthetic_recovery.R",
    hard_fail = FALSE
  )
  # The ENSO crossbasis extractors lifted out of 09 into utils/figures_cb.R reproduce
  # 09's stored manuscript artefact (loads one deployed brms fit; no sampler). SKIPs
  # (passes) when the artefact / deployed-fit archive is absent.
  results$figures_cb_extractors_equivalence <- run_one(
    "analysis/tests/test_figures_cb_extractors_equivalence.R",
    hard_fail = TRUE
  )
  # 09's manuscript artefacts vs a baseline copy. Compares the artefacts as they stand;
  # SKIPs (passes) without a
  # baseline. Re-deriving first is opt-in (`--run`, tens of minutes).
  results$results_artefacts_equivalence <- run_one(
    "analysis/tests/test_results_artefacts_equivalence.R",
    hard_fail = TRUE
  )
  # C-kernel smoke/structural checks, unit-order invariance and the detected feed. Each
  # compiles Csnippets (slow), so full suite only.
  results$spatpomp_c_variant <- run_one(
    "analysis/tests/test_spatpomp_c_variant.R",
    hard_fail = TRUE
  )
  results$spatpomp_c_variant_dur <- run_one(
    "analysis/tests/test_spatpomp_c_variant_dur.R",
    hard_fail = TRUE
  )
  results$spatpomp_unit_order <- run_one(
    "analysis/tests/test_spatpomp_unit_order.R",
    hard_fail = TRUE
  )
  results$spatpomp_detected_feed <- run_one(
    "analysis/tests/test_spatpomp_detected_feed.R",
    hard_fail = TRUE
  )
  results$spatpomp_who_obs <- run_one(
    "analysis/tests/test_spatpomp_who_obs.R",
    hard_fail = TRUE
  )
  results$spatpomp_filter_forecast <- run_one(
    "analysis/tests/test_spatpomp_filter_forecast.R",
    hard_fail = TRUE
  )
  results$spatpomp_grouped_c <- run_one(
    "analysis/tests/test_spatpomp_grouped_c.R",
    hard_fail = TRUE
  )
  # Smoke/structural checks for the monthly-size (mon) spatPomp Csnippet kernel
  # (mon_mu_block): degenerate active/mso/mu trajectories plus ungrouped and
  # grouped bpfilter logLik finiteness. Sibling of spatpomp_grouped_c,
  # one phase later (mon builds on duration's active/mso machinery).
  results$spatpomp_c_variant_mon <- run_one(
    "analysis/tests/test_spatpomp_c_variant_mon.R",
    hard_fail = TRUE
  )
}

# config_defaults.yml vs get_default_options(): same key paths and same values.
# HARD: the YAML-only simulator.filter.* /
# cfsv2_target_peak keys are defaulted in R at their consumers' fallbacks and
# the YAML values follow R, so the allowlists are empty; a YAML value that is
# not the runtime default is silent documentation drift (and the batch
# generator copies it into every fit-batch config).
results$config_defaults_parity <- run_one(
  "analysis/tests/test_config_defaults_parity.R",
  hard_fail = TRUE
)

# Single-lag DLNM (onebasis) machinery — soft gate.
results$cb_singlelag <- run_one(
  "analysis/tests/test_cb_singlelag.R",
  hard_fail = FALSE
)

# ENSO lead columns for the -11..11 single-lag sweep (join_enso_lags max_lead) — soft gate.
results$enso_leads <- run_one(
  "analysis/tests/test_enso_leads.R",
  hard_fail = FALSE
)

# Deployed ENSO leads in the simulator: kernel lead
# crossbasis + scenario override, sampler enso_lead<k> matrices — hard gate.
results$simulator_enso_leads <- run_one(
  "analysis/tests/test_simulator_enso_leads.R",
  hard_fail = TRUE
)

# make_figure_dir()/06_fit_figures.R collision ratchet. SOFT: 0
# collisions across all configs; a regression here is a
# silent figure overwrite, not a fit-breaking error, so it's tracked but does
# not halt the suite.
results$figure_path_collisions <- run_one(
  "analysis/tests/test_figure_path_collisions.R",
  hard_fail = FALSE
)

# Lenssen analog-set subsetting (counterfactual ENSO-state scenarios).
# HARD gate: case 7 is the only check that the row filter actually reaches
# sample.int() (a misplaced filter still returns plausible numbers), and case 10
# is what licenses renaming the pre-existing `_set-all` output trees instead of
# recomputing them. Needs no fitted models or real cache — fully synthetic.
results$enso_analog_sets <- run_one(
  "analysis/tests/test_enso_analog_sets.R",
  hard_fail = TRUE
)

# CFSv2 seam when the ONI record ends after the raw Niño 3.4 record (CPC's
# normal state). HARD gate: a broken seam puts a one-month ONI = 0 notch
# in the forecast ENSO input (manuscript ONI envelope AND simulator forecast
# arm). Synthetic records only.
results$enso_cfsv2_seam <- run_one(
  "analysis/tests/test_enso_cfsv2_seam.R",
  hard_fail = TRUE
)

# CFSv2 hindcast calibration (CPC PDF + spread correction replica) and the
# forecast-phase crossbasis knot clamp helper. HARD gate: synthetic records only.
results$enso_cfsv2_calibration <- run_one(
  "analysis/tests/test_enso_cfsv2_calibration.R",
  hard_fail = TRUE
)

# Run --------------------------------------------------------------------

t_start <- proc.time()[["elapsed"]]
if (opt$jobs <= 1L) {
  results <- lapply(results, exec_one)
} else {
  log_dir <- file.path(report_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  # Longest first, so a slow test does not start last; a test without a
  # previous timing counts as the longest.
  prev <- if (file.exists(timings_path)) {
    tryCatch(
      utils::read.csv(timings_path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }
  est <- if (is.null(prev)) {
    rep(NA_real_, length(results))
  } else {
    prev$elapsed[match(names(results), prev$test)]
  }
  est[is.na(est)] <- Inf
  ord <- names(results)[order(-est)]
  message(sprintf(
    "run_all: %d test(s) on %d worker(s); per-test logs in %s",
    length(results),
    opt$jobs,
    log_dir
  ))
  done <- parallel::mclapply(
    ord,
    function(nm) {
      exec_one(results[[nm]], file.path(log_dir, paste0(nm, ".log")))
    },
    mc.cores = opt$jobs,
    mc.preschedule = FALSE
  )
  names(done) <- ord
  for (nm in ord) {
    if (!is.list(done[[nm]])) {
      # a worker that died returns an error object
      done[[nm]] <- list(
        script = results[[nm]]$script,
        args = paste(results[[nm]]$args, collapse = " "),
        status = -1L,
        hard_fail = results[[nm]]$hard_fail,
        elapsed = NA_real_
      )
    }
  }
  results <- done[names(results)]
  for (nm in names(results)) {
    if (results[[nm]]$status != 0) {
      log <- file.path(log_dir, paste0(nm, ".log"))
      message("\n--- FAIL ", nm, " (", log, "), last lines:")
      if (file.exists(log)) {
        message(paste(
          utils::tail(readLines(log, warn = FALSE), 30),
          collapse = "\n"
        ))
      }
    }
  }
}
t_wall <- proc.time()[["elapsed"]] - t_start
utils::write.csv(
  data.frame(
    test = names(results),
    elapsed = vapply(results, function(r) r$elapsed, numeric(1))
  ),
  timings_path,
  row.names = FALSE
)

# Summary ----------------------------------------------------------------

con <- file(summary_path, "w")
on.exit(close(con), add = TRUE)
cat("# Test harness summary\n\n", file = con)
cat(
  sprintf(
    "Generated: %s (jobs %d, wall %.0f s)\n\n",
    Sys.time(),
    opt$jobs,
    t_wall
  ),
  file = con
)
cat("| test | status | hard_fail | result | elapsed_s |\n", file = con)
cat("|------|--------|-----------|--------|-----------|\n", file = con)
hard_failed <- FALSE
for (nm in names(results)) {
  r <- results[[nm]]
  res <- if (r$status == 0) "PASS" else "FAIL"
  if (r$status != 0 && r$hard_fail) {
    hard_failed <- TRUE
  }
  cat(
    sprintf(
      "| %s | %d | %s | %s | %.1f |\n",
      nm,
      r$status,
      r$hard_fail,
      res,
      r$elapsed
    ),
    file = con
  )
}

if (hard_failed) {
  message("\nrun_all: HARD FAILURE — see ", summary_path)
  quit(status = 1)
} else {
  message("\nrun_all: all hard gates PASSED — see ", summary_path)
  quit(status = 0)
}
