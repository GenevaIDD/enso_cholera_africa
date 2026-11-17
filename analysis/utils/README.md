# `analysis/utils/` — no numeric prefix

Everything in this directory shares one property: **no numeric prefix, not a
pipeline step**. Most of it is also a sourced-only library (never invoked
directly, only `source()`-d by the numbered pipeline scripts in `analysis/`),
but §3 is the explicit exception — three standalone CLI tools that live here
too. See `CLAUDE.md` → "Pipeline file naming convention" for the
runnable-step / library / non-step-tool three-way rule.

1. **The 17 `utils.R` modules** (§1 below) — domain modules; every top-level
   definition lives in exactly one module. `analysis/utils.R` is a facade that
   `source()`s all 17 in dependency order.
2. **9 single-purpose sourced-only libraries** (§2 below) — no numeric prefix,
   no CLI, sourced by name.
3. **Standalone CLI tools** (§3 below) — `make_batch_configs.R`,
   `make_forecast_configs.R`, `install_occ_brms.R` and others. Unlike groups
   1–2, these are **not** sourced-only — they are invoked directly via
   `Rscript`. They live here because they are not pipeline steps.

## 1. The 17 `utils.R` modules

| Module | Covers | Top-level fns |
|---|---|---|
| `base.R` | `%||%`, `wquantile()`, the ENSO lag naming (`enso_lag_col_name()`, `sl_token()`, `sl_tag()`), `map_workers()`, and the shared infrastructure of §7 ( `content_digest()`, the `.TTLO_*` constants) | 1 |
| `paths.R` | `make_run_file` and the whole `make_*_run_file` family — the universal artefact-naming hub (incl. `make_loo_compare_run_file()` for the per-batch LOO lag-ranking outputs, `grouped_lag_config_path()` for the `<batch>_sl-<lag>_grp-<group>.yml` config of a grouped single-lag batch, and `read_deployed_lag()` — the single reader of a deployed occ/dur/mon best-lag fit's `sl<N>` tag, `_ep-gs`-aware, shared by both forecast/hindcast generators; `make_enso_continuation_library_file()` for the ENSO analog continuation library and `make_enso_continuation_skill_file()` for its skill-by-lead table; `default_outbreak_data_dir()` / `OUTBREAK_DATA_VINTAGE`, the data leg's one `--indir` default, outbreak_data_Aug_2026, used by 00, 01, 01b, 01_nz, the data-overview figure and both 02 scripts; `BASENAME_MAX_BYTES` / `assert_basename_length()`, the 255-byte name guard every builder returns through, and `check_run_file_names()`, every kind of a fit config's track built up front by the 04 workers and `generate_batch()`) | 14 (+`make_occ_run_file`, see below) |
| `gid.R` | GADM id manipulation (`flatten_gid`, ADM1 filtering, `mark_tza_mainland_adm0()` -- TZA mainland is ADM0 throughout) | 9 |
| `clim_compile.R` | CSV/parquet climate file discovery, per-gid stats/extremes compilation | 34 |
| `rain_enso.R` | ENSO and rainfall, descriptive (Figure 1c rain row, T1a / T1b): rainy-season climatology (`is_rainy_month()`, `circular_month_runs()`, `rain_seasons_by_group()`: weighted-vote main rainy season(s) of groups of units) and the rainfall tables of 09's `rain_enso` step and `03_prepare_rain_adm1.R` (`adm1_monthly_precip()`, `precip_dev_series()`, `precip_oni_corr()` / `precip_oni_group_corr()`, `adm1_population_from_gid()`, `adm1_rain_seasons()`). No model term reads it. | 34 |
| `clim_db.R` | SQLite accessors for the compiled climate stats/extremes DBs | 5 |
| `occ_data.R` | `make_surv_data`, `load_surv_data`, scaling, climate lags (occurrence track); `lp_coverage_season_years()` (extracted from `add_no_outbreak_years()`), and the gid-season builder `gid_at_risk_years()` / `gid_season_surv_input()` (`episodes.key: gid_season`); `national_zero_surv_rows()` (no-onset rows for national-zero months, `03_prepare_occ.R --national_zero_fill`) | 13 |
| `siz_data.R` | `make_size_data`, `read_size_data` (size track) | 4 |
| `enso.R` | ENSO phase (`classify_enso_phase(x, rule = "threshold" or "episode")`, `enso_phase_definition()`, `ENSO_PHASE_THRESHOLDS` / `ENSO_EPISODE_MIN_RUN` / `ENSO_PHASE_LABELS`: every ±0.5 phase label), ENSO index fetch (`latest_oni_month()` = the cached ONI's last month; ERSSTv6 / OISST Niño3.4), plume/analog/CFSv2 sampling, CFSv2 hindcast fetch + CPC-style PDF/spread calibration (`fit_cfsv2_calibration()`, `apply_cfsv2_calibration()`), analog-set classification, ENSO past the CFSv2 span (`continue_enso_ar1()`; k-NN analog continuation `load_enso_continuation_library()`, `enso_analog_neighbours()`, `continue_enso_analog()`, `enso_continuation_provenance()`; `qs2` package) | 11 |
| `population.R` | Population crosswalk, WPP rescaling, LP resolution, composite-case apportionment; `canonical_level_lookup()` (geometry-derived level for the canonical episode table) | 12 |
| `composite_lp.R` | `lp_composite_status`, `collapse_composite_lps`, `.match_lp_members` | 3 |
| `provenance.R` | Population/fit provenance stamps, mean-rate table, prep-policy stamps; the prep_counts files (`write_step_counts(indir =)` records the outbreak-data dir + `data_vintage`, `read_step_counts()` leaves them out, `read_step_provenance()` reads them); the shared P3 gate `canonical_admin_level_enabled()` / `.apply_p3()` (all four loaders + anchor views; calls `population.R::apply_canonical_admin_level()` only when the re-key runs); the episode fingerprint (`digest_episodes()`, `episode_stamp_of()`, `fit_episode_stamp_of()`, `assert_episode_key()`, `assert_episode_coherence()`, `sim_episode_check()`, `prep_stamp_inherit()`, and for deploy-generated configs `apply_episode_key_to_sim_config()` / `preflight_sim_episode_inputs()`) and the shared 03_prepare_* CLI surface (`prep_cli_options()`, `validate_prep_cli()`, `resolve_inherited_episode_key()`) | 33 |
| `occ_post.R` | Occurrence RE-stat adapter (brms) + `canonicalize_term()` | 5 |
| `siz_post.R` | Size RE-stat adapter (brms) + pp_check table builder | 3 |
| `dlnm.R` | DLNM crossbasis design matrix / surface extraction; `.clamp_cb_range()` boundary-knot hold (smoothing + forecast) | 4 |
| `plots.R` | `plot_surv_dat`, `plot_clim_data`, `decorate_covar`, covariate label/type helpers | 5 |
| `sim_io.R` | Simulator IO: warmup panels, task-plan build/write (the surveillance-coverage mask is in `coverage.R`); the compact (format-2) monthly summary `sim_monthly_save()` / `sim_monthly_compact()` / `sim_monthly_expand()` / `sim_monthly_bind()` / `sim_monthly_cells()`, expanded by `sim_io_read()` so long-table readers see the long table (gate `tests/test_sim_monthly_compact.R`) | 6 |
| `coverage.R` | The surveillance-coverage definition `build_coverage_mask()` (+ national-zero months from `make_national_zero_months_file()`, see `national_zero.R`). Sourced by 07c, `smoothing_core.R` (08c), 09a/09b and the `utils.R` facade, next to `sim_io.R`. | 1 |

`make_run_file` is a router that dispatches on the `kind` prefix (`size_` →
`make_size_run_file`, `sim`/`sim_` → `make_sim_run_file`, `smooth`/`smooth_` →
`make_smooth_run_file`, `spatpomp`/`spatpomp_` → `make_spatpomp_run_file`, `pop_` →
`make_pop_run_file`, else → `make_occ_run_file`, the occurrence path builder in
`paths.R`).

### Cross-module dependency table

Modules are **pure definitions — they never `source()` each other** (that
would re-trigger `.POP_XWALK_CACHE`'s reset and create load cycles). A script
that sources a subset of modules must source **that module plus its full
transitive closure**, in any order (R resolves free variables at call time, so
intra-closure ordering never matters). `base.R` is the universal floor —
always source it.

| Module | Depends on (direct) |
|---|---|
| `base.R` | — |
| `gid.R` | `base` |
| `paths.R` | `base`, `provenance`, `config` *(call-time: `.flag_val()` in two CLI-flag readers)*; `digest` via `::` (`.cap_file_name()`, occurrence names over `BASENAME_MAX_BYTES`) |
| `clim_compile.R` | `gid`, `paths` |
| `rain_enso.R` | none (base R, `dplyr::`) |
| `clim_db.R` | `gid` |
| `occ_data.R` | `clim_compile`, `clim_db`, `gid`, `paths`, `siz_data`, `provenance` *(call-time: `assert_contagion_partition()` / `assert_episode_key()` read the prep policy stamp)*, `population` *(call-time, only under `occurrence.canonical_admin_level`: `.apply_p3()` -> `apply_canonical_admin_level()`)* |
| `siz_data.R` | `base`, `clim_compile`, `clim_db`, `gid`, `provenance`, `anchor` *(§5 module: `read_size_data()` builds `mean_rate` through the shared `build_base_rate_table()`)*, `formula` *(§5 module; call-time, only when the ADM2 leaf scope is not "all": `.leaf_columns()`)* |
| `enso.R` | — |
| `population.R` | `composite_lp` |
| `composite_lp.R` | `population` *(mutual dependency with `population.R` — source both together)* |
| `provenance.R` | `base`, `paths`, `optparse` *(call-time, `prep_cli_options()` only)*, `population` *(call-time, `build_mean_rate_table()`'s P3 via `.apply_p3()`)* |
| `occ_post.R` | `paths`, `plots` |
| `siz_post.R` | `base` |
| `dlnm.R` | `base` |
| `plots.R` | `clim_compile`, `clim_db` |
| `sim_io.R` | `base`, `occ_data`, `paths`, `spatpomp_core` *(call-time, only when `build_obs_panel()` gets `mon_dat`: `.aggregate_footprint_counts()`)*, `config` *(call-time: `for_each_hindcast_run()` -> `make_options_from_config()`)*, `gid` *(call-time: `gid_iso3()` -> `flatten_gid()`)* |
| `coverage.R` | `base` (`%>%`), `paths` *(call-time: `make_national_zero_months_file()`)* |
| `spatial_contagion.R` | `base`, `gid`, `paths`, `config`, `coverage` / `occ_data` *(the contagion coverage report only: `build_coverage_mask()`, `active_last_month()`)* *(all call-time only — see §2 "Other domain modules"; `config` only for a partition-restricted term: `load_grouping()`)* |

This table was derived from `codetools::findGlobals()` run on each module's
function bodies (AST-based, so it is blind to symbols referenced only inside
non-standard-evaluation contexts, e.g. a purrr formula lambda; none of the
module-internal call sites hit that limitation).

`paths.R` ↔ `provenance.R` is not circular despite both appearing in each
other's neighbourhood at first glance: `paths.R` depends on `provenance.R`
(some `make_*_run_file` variants read provenance-stamped metadata), and
`provenance.R` depends on `paths.R` (stamp file paths are themselves built via
`make_*_run_file`). Both are always sourced together in practice by every
pipeline script that needs either.

## 2. Single-purpose libraries

| File | Original path | Covers |
|---|---|---|
| `config.R` | `analysis/config.R` | Config-driven option parsing (the **options core**): `make_options_from_config()`, `get_default_options()`, CLI routing, grouping, `detect_model_from_yaml()` (which model block a raw config carries). Three-layer merge: `config_defaults.yml` → YAML config → CLI flags. Its formula and prior builders are in `formula.R` / `formula_legacy.R` / `priors.R` (§5). Depends on `base`; `deploy_core` at call time only (`exists()`-guarded). |
| `simulations.R` | `analysis/utils_simulations.R` | Simulator shared helpers: `load_artefacts(opt)`, `check_frame_vs_draws()`. |
| `figures_cb.R` | `analysis/utils_figures_cb.R` | Pure helpers for the DLNM crossbasis figure scripts (`build_onebasis_grid()` etc.), plus the ENSO crossbasis posterior extractors shared with the manuscript phase (`extract_adm1_irr()`, `extract_country_irr()`, `extract_country_erf()`, `cb_target_matrix()`, `country_irr_from_adm1()`, and the per-draw `extract_adm1_logrr_draws()` that 09 caches), `plot_irr_map()`, `COMPOSITE_CONDS` and its label rule `enso_conds_named()` (09 builds the configured `results.derivation.enso_conds` with it), the taxdat basemap loaders (`load_taxdat_afr_sf()`, `load_africa_basemap()`), `drop_small_islands()` (drops Mauritius / Réunion from the Africa outlines), the map hatch pattern `mask_hatch()` (diagonal / crossed) with its guard `pattern_fill_supported()` (ggplot2 >= 3.5), and `ps_diverging_fill_scale()` (the manuscript maps' horizontal diverging colourbar, title on the left). Call-time dependency: `build_onebasis_grid()`/`brms_group_re_scaffold()` only; all packages `::`-qualified. |
| `population_core.R` | `analysis/02a_population_core.R` | Population-extraction-array kernel shared by the `02a`/`02b`/`02c` population sub-pipeline. |
| `population_qc.R` | `analysis/population_qc.R` | Population-conservation QC gates. |
| `simulator_core.R` | `analysis/07a_simulator_core.R` | v1 simulator kernel: `simulate_one_tu()`, `build_climate_panel_for_tu()`, `derive_min_seas_from_panel()`. Refactors reuse its primitives rather than change it. Also owns the monthly (`mon`) case-count recipe shared by every kernel: `monthly_log_mu()` (evaluate_lp + mso basis via `.inject_mso()` + the pop / log_offset / anchor offsets), `draw_monthly_cases()` (the NB draw, truncated at the cap), and the monthly incidence cap `simulator.incidence_cap` (default 0.15): `incidence_cap_of()` reads and validates it, `monthly_mu_cap()` turns it into the per-unit ceiling cap x pop used by `simulate_one_tu()`, the C kernel's `mon_mu_cap` covariate (`spatpomp_core_c.R`) and `extract_smoothing_one_tu()`. Gates: `tests/test_monthly_mean_golden.R` (uncapped output of every site), `tests/test_incidence_cap.R`. Also owns the zero-filled climate guard: `.assert_no_zero_filled_climate()` refuses a run whose fits carry a local-climate fixed effect when a panel row in the run window has all four climate columns exactly 0 (the post-DB zero fill of `03_prepare_sim_panel_projection.R`); with `.draws_fe_names()` and `.sim_panel_window()`, called by both `simulations.R` loaders, `07b_simulator_v1.R`'s own loader and `smoothing_core.R::smoothing_load_artefacts()` (no extra module dependency: all of them source this file). Gate: `tests/test_climate_zero_fill_guard.R`. Kernel speed-up (bit-identical): `lp_plan()` resolves `evaluate_lp()`'s random-effect key matches once per (bundle, rows, `drop_region`) and `.lp_coefs()` caches the per-location coefficients per posterior draw (`evaluate_lp(plan =, smooth_X =)`, `monthly_log_mu(plan =, smooth_X =)`; no plan = a throwaway one); `build_lag_maps()` + `update_lags(maps =)` aggregate on integer codes (keyed `tapply()` path for NA / non-0/1); `.warm_row_index()` selects `build_warmup_state()`'s rows without dplyr; `simulate_one_tu()` keeps a plan per context and track, the spline design per month, the NegBin shape per replicate and the event rows as vectors. The pre-speed-up bodies are frozen in `simulator_core_legacy.R`. Gate: `tests/test_simulator_kernel_golden.R` (HARD). |
| `spatpomp_core.R` | `analysis/07a_spatpomp_core.R` | Shared library of the spatPomp filter: observation/covariate tables, measurement Csnippets, block particle filter runners. Also the builder's shared blocks (`.sp_filter_window_months()`, `.sp_enso_covar_names()`, `.sp_cb_covar_names()`), the warm-start row builder `.state_vec_to_warmup_row()` behind both handoff paths, and the smoothing-cache helpers (`.resolve_smoothing_cache()`, `.handoff_from_cache()`, `.cache_handoff_label()`, `get_tu_list()`, `.make_filter_cache_path()`, `load_smoothing_distribution()`). Forecast handoff selection among the smoothing replicates: `select_handoff_warmup()` (`simulator.forecast.handoff_weighting`: `uniform` (default), balanced by global replicate id; `loglik` = `resample_by_loglik()`) and `handoff_loglik_ess()`. Compact trajectory store: `run_one_bpfilter_replicate()` dispatches on `simulator.filter.traj_storage` (`full` = `spatPomp::bpfilter`; `pruned` (the default) / `pruned_float32` = `.bpfilter_traj_compact()`, a copy of spatPomp 1.1.0's `bpfilter.internal`, run in the spatPomp namespace, with ancestry pruning every `traj_prune_every` steps and optional float32 storage via the `float` package, needed only by `pruned_float32`). Gate: `tests/test_bpfilter_traj_compact.R`. |
| `spatpomp_core_c.R` | `analysis/07a_spatpomp_core_c.R` | The spatPomp builder, `build_spatpomp_tu_c()` (Csnippet rprocess / rinit; the only builder). Conditionally `source()`s `spatpomp_core.R` (a documented exception to the no-cross-sourcing rule). Its three grouped precompute wrappers share `.stitch_grouped_lp()`. |
| `smoothing_core.R` | `analysis/08a_smoothing_core.R` | Historical-smoothing sub-pipeline kernel: `smoothing_make_config`, `smoothing_filter_chunk`, ribbon-plot helpers. A `_core` kernel: sources its own 15-module closure. The event-table and ribbon aggregators are one implementation each (`.smoothing_events_long_impl()`, `.agg_and_ribbon_impl()`). 08c's observed monthly-cases series: `smoothing_obs_cases_mon()`, one admin level per (gid, month) by the filter's footprint rule via `footprint.R::finest_level_rows()` (gate `tests/test_smoothing_obs_cases_footprint.R`). |

These are single-purpose libraries; they live in `utils/` per the
library-vs-runnable-step rule in `CLAUDE.md`.

### Other domain modules

| File | Covers |
|---|---|
| `spatial_contagion.R` | Distance-weighted neighbour-contagion covariate for the occurrence track: `build_adm1_dist_km()` (cached GADM ADM1 centroid distance matrix, km), `kernel_from_dist()`, `neighbour_pressure()`, `build_contagion_table()`, and `build_nbr_kernel()` -- the single kernel definition (scale raw/rownorm x partition none/group) shared by prep and every simulator kernel; simulator side: `resolve_nbr_spec()` (scale/partition parsed from the fit's FE name), `build_partition_nbr_kernel()`, `nbr_kernel_for_frame()`. Covariate names come from `paths.R::nbr_covariate_name()` / `nbr_covariate_pattern()`. Sourced by `03_prepare_occ.R` and the simulator drivers. |

A domain module, named by domain with no phase semantics per the `CLAUDE.md`
naming convention. Its two intra-`utils/` dependencies (`gid.R::flatten_gid()`,
`paths.R::make_adm1_dist_file()`) are **call-time, not source-time**, following
the same pattern `gid.R` already documents for its own
`population_core.R::load_gadm_units()` dependency: the module does not
`source()` them, so callers must source all three.

## 3. Standalone CLI tools

| File | Original path | Covers |
|---|---|---|
| `make_batch_configs.R` | `analysis/make_batch_configs.R` | Batch config generator: `expand.grid()` over parameter sets → YAML files + `manifest.yml`. Handles grouped batches (expand over `source × group`). Usage: `Rscript analysis/utils/make_batch_configs.R --batch <name> \| --all \| --list`. |
| `make_forecast_configs.R` | `analysis/make_forecast_configs.R` | Per-ENSO-group forecast/counterfactual config generator. Derives each group's horizon from the deployed best-lag fits' `sl<N>` tag and cached CFSv2 metadata (`last_cfsv2_month + min(best_lag_occ, best_lag_dur, best_lag_mon)`); writes `configs/forecast_country_by_group/*.yml` + `configs/forecast_horizons.yml`. Usage: `Rscript analysis/utils/make_forecast_configs.R [--fits_dir <dir>] [--analog_sets ...]`. |
| `install_occ_brms.R` | `analysis/install_occ_brms.R` | One-off helper to install/build the occurrence brms backend dependencies. Usage: `Rscript analysis/utils/install_occ_brms.R`. |
| `compact_sim_outputs.R` | *(new)* | Rewrites existing forward-sim `monthly_summary.qs` files (parts, then merged) in the compact format 2 of `sim_io.R`; a file is replaced only when it reads back identical, mtime and mode kept, plus a sidecar `<file>.format2` (original size / mtime, new size) that `results_core.R::.results_stamp_stale_reason()` reads through `sim_io.R::sim_monthly_converted_from()`, so stamps that recorded the file before conversion stay fresh. Usage: `Rscript analysis/utils/compact_sim_outputs.R --root generated_data/sim_spatpomp --pattern '<glob>' [--dry-run] [--log <csv>]`; run by hand. Gate: `tests/test_sim_monthly_compact.R` [6]. |

Unlike §1–2, these are not library-shaped; they live here because they are
not pipeline steps. Each keeps its CLI entrypoint (`optparse` option list),
guarded by `sys.nframe() == 0L`, so tests and other tools can `source()` it to
use its internals (`test_size_batch_config_determinism.R` and
`test_size_seasonal_term.R` both `source("analysis/utils/make_batch_configs.R")`).

## 4. Shared helpers

- **`build_coverage_mask()`** — one definition (`coverage.R`, §1 table above),
  used by `utils/smoothing_core.R` and `07c_simulator_postprocess.R` directly
  or via the `utils.R` facade.
- **Test-assertion micro-harness helpers** (`expect`, `fail`, `expect_true`,
  `expect_equal`) — in `analysis/tests/helpers.R`, sourced by the test files.
  A few files keep a local override (defined after the `source()` call, so it
  shadows the shared one) where the behaviour genuinely differs: looser
  floating-point tolerance for cross-kernel (R vs Csnippet)
  numerical-equivalence checks, or extra failure-label tracking. `helpers.R`'s
  header comment lists them.

## 5. Further modules

| Module | Covers | Dependencies |
|---|---|---|
| `dur_data.R` | Duration track: `make_duration_panel`, `derive_series_end`, `flag_right_censored`, `bin_mso` / `.DUR_MSO_BREAKS`, `load_duration_data`. Expands the size track's episode table into one row per **active month** and marks the month of observed cessation, so outbreak duration is fitted as a discrete-time cessation hazard (Poisson with a log-days offset, mirroring the occurrence backend). Also the duration stamp's summaries: `duration_span_summary()` (whole-month episode span median / IQR and median elapsed months, stamped by `03_prepare_dur.R` as `median_span_months`, `span_months_q25` / `_q75`, `median_elapsed_months`) and `outbreak_episode_guard_flags()` (per-episode flags of what `occ_data.R::guard_outbreak_episodes()` would drop, stamped as `n_guard_*`; nothing is dropped). Gate: `tests/test_duration_span_stamp.R`. | `clim_db`, `occ_data`, `paths`, `siz_data`, `gid` *(call-time, only under `duration.admin_level: ADM1-ADM2`: `apply_admin_level_filter()` in `duration_episode_rows()`)*, `formula` *(call-time, only when the ADM2 leaf scope is not "all": `.leaf_columns()`)* |
| `footprint.R` | The episode-footprint collapse, shared between `03_prepare_siz.R` (producer) and `03_prepare_mon.R` (consumer): `FOOTPRINT_KEY`, `unify_footprint_ids`, `footprint_member_rows`, `collapse_to_footprint`, `aggregate_monthly_to_footprint` (with a `key` argument), so the monthly track re-derives the *same* grouping instead of guessing at it from `(uid, location_period_id)`. Also the cross-track canonical episode table (`episodes.key: gid_season`): `build_canonical_episodes()` -- called with identical arguments by `03_prepare_occ.R` and `03_prepare_siz.R` -- over `canonical_season_calendar()`, `canonical_episodes()`, `finest_level_totals()` (parity-tested against `spatpomp_core.R::.aggregate_footprint_counts()`), `finest_level_rows()` (the rows of the level `finest_level_totals()` counts, order kept; the observed monthly cases of 08c's `smoothing_obs_cases_mon()` and 09's `derive_smoothing()`), `representative_ids_by_gid()`, `footprint_key()`. The geometry-derived level comes in as a caller-supplied `level_of` (`population.R::canonical_level_lookup()`), so the module stays base-only. | `base` |
| `dur_post.R` | Duration-track postprocessing adapter. Sibling of `occ_post.R` / `siz_post.R`: the duration track reuses the occurrence RE-stats adapter unchanged, so this module supplies only what duration needs of its own — `make_duration_pp_check_tbl()`, a two-panel posterior predictive check (monthly cessation-hazard calibration by `mso_band`, and the implied duration distribution derived by pushing per-month hazards through to a duration) rather than a predicted-vs-observed scatter, since the duration response (`y = as.integer(ended)`, one row per active month) is not itself the endpoint of interest. | `base` |
| `anchor.R` | Empirical **baseline-rate anchors**, shared by all four brms tracks: `.ANCHOR_SPEC` / `anchor_spec()` / `anchor_enabled()`, the four per-track views, `build_base_rate_table()`, `resolve_base_rate()`, `attach_base_rate()`, `persist_base_rate_table()`. Generalises the size track's `offset(log(mean_rate))` to occurrence / duration / monthly, so every track satisfies `denominator(anchor) == product of that track's OTHER offsets` and a single `Intercept` prior at residual scale is correctly calibrated on all four. | `base`, `paths`, `provenance`, `occ_data` |
| `figures_io.R` | Figure-saving helpers: `write_figure_files()` / `figure_saver()` (`dir.create(recursive = TRUE)` + `ggsave()` + `message("Saved: ...")`), used by the figure scripts (`06_fit_figures.R`, `_cb.R`, `_cb_adm1.R`, `_lag_compare.R`, `07d_simulator_figures.R`); plus `true_minus()` / `ps_true_minus()` (true minus sign in numeric labels / in a plot's scales, for `10_manuscript_figures.R`). Deliberately not folded into `make_run_file()` (`paths.R`'s `kind` is a closed `match.arg()` set) — path algebra (`make_figure_dir()`, `paths.R`) stays separate from writing (`write_figure_files()`, here) so `paths.R` never needs a ggplot2 dependency. | `base` (callers separately source `paths.R` for `make_figure_dir()`) |
| `mon_data.R` | Monthly track: fit-time loader for the monthly case-count panel — country/admin filtering, `cases_prev`/`log_cases_prev` derivation from the full panel before filtering, harmonic seasonal columns, and integer-rounding of `cases_month` for the NegBin response. Sibling of `occ_data.R` / `siz_data.R` / `dur_data.R`. | `base`, `paths`, `provenance`, `siz_data`, `formula` *(call-time, only when the ADM2 leaf scope is not "all": `.leaf_columns()`)* |
| `loo_lags.R` | Shared ENSO single-lag selection rule used by `05_fit_postprocess_loo_lags.R`, `05_fit_postprocess_lag_compare_surface.R`, and the deploy path: `select_best_lag()` picks the max lag within `se_mult` SEs of the rank-1 elpd (occurrence), with an occurrence-floor-constrained variant for the size/duration/monthly tracks. Also owns the one ranking-file locator for the `loo_compare_lags_*.rds` caches (`loo_compare_ranking_files()` / `latest_ranking_file()` / `read_latest_ranking()`: `generated_data/loo_compare/<batch>/` first, legacy `analysis/configs/<batch>/` fallback, latest by filename date then mtime) and `best_lags_by_group()` — the per-group best lag (latest ranking -> `select_best_lag()`, with overrides) the grouped figure scripts use, so they pick the same lag as the deployed fits and the manuscript. Also the Pareto-k summary of each lag's PSIS-LOO (`pareto_k_summary()`, `ranking_pareto()`: the additive `PARETO_K_COLS` 05 appends to every ranking) and the per-(batch, group) lag sidecar of `05_fit_postprocess_loo_lags.R --lag_sidecar` (`lag_sidecar_entry()` / `build_lag_loo_sidecar()` / `check_lag_sidecar_rows()`, row-key digest `fit_row_digest()`, population crossbasis draws `fit_cb_pop_draws()`, reader `latest_lag_sidecar_file()` / `read_lag_loo_sidecar()`; file shape `LAG_LOO_SIDECAR_RE`, which never matches the ranking regex). Gate: `tests/test_lag_averaging.R`. Also `check_occ_floor_offset()`: the deploy floor `occ_floor_offset` against the duration stamp's `median_span_months` (warns on a mismatch, messages when unstamped), called by `deploy_core.R::deploy_select_lags()`; gate `tests/test_duration_span_stamp.R`. | `base`, `paths` (call-time, for `make_loo_compare_run_file()` / `derive_config_run_tag()`); `loo` / `digest` / `posterior` / `brms` via `::` |
| `dependence.R` | Phase-5 size↔duration dependence recipe behind `05_fit_postprocess_dependence.R`: randomized PIT for the NegBin size fit and censoring-aware chained-hazard PIT for the duration fit, per-draw Kendall's τ, and ADM1-intercept correlation decomposition. | `siz_data` |
| `obs_sim_compare.R` | Honest observed-vs-simulated case-count building blocks for the hindcast skill: population-downweighted ADM1 totals, strict gid intersection (`obs_sim_bases()`, which also returns the per-key `obs_rk` totals), the scalar scoring rules `crps_ens()` / `crps_fair()` (unbiased, fair across ensemble sizes) / `pinball_loss()`, the observed-frame builder `LVL_RANK` / `build_mon_obs_frame()`, and `map_distinct()` (a function evaluated once per distinct value of a replicate-frame column). Its `crps_ens()` / `region_key_of()` are also the single copies `sim_metrics.R` uses. | `gid`, `population` |
| `sim_metrics.R` | Simulator validation and hindcast-skill metrics, used by `07c_simulator_postprocess.R` and `07d_simulator_figures.R`: `load_observed_from_panel()`, `compute_pit_binom()`, `compute_monthly_freq()`, `compute_crps_annual()`, `compute_annual_count_by_geo()`, `compute_ar_metrics()`, `reconstruct_lag_series()`, `compute_contagion_conditionals()`, `compute_episode_case_totals()` / `attach_episode_case_totals()`, `compute_crps_total_cases()`, `.km_estimate()`, `compute_duration_dist()`, `compute_active_fraction()`, `compute_sim_kendall_tau()`, `compute_hindcast_series()`, plus the shared geography scorers `.sim_frame_geo()` / `.geo_crps_level()`. Golden gate: `tests/test_sim_metrics_golden.R`. | `base`, `gid`, `obs_sim_compare`, `sim_io`, `simulator_core` (`update_lags()`), `dependence` (`paired_kendall_tau()`) |
| `results_core.R` | Shared derive/plot module for the manuscript-results phase (sourced by both `09b_postprocess_results.R` and `10_manuscript_figures.R`): manifest/grouping loaders, the `ENSO_ANALOG_SETS`-derived registries, the `make_results_run_file()` freshness/provenance wrappers `results_read_if_fresh()` / `results_write()` / `results_require()`, the deployed-fit archive readers (`load_best_lag_manifest()` / `track_fit_files()`), the `TRACK_ADVERSE_SIGN` direction registry, and the simulator-partition aggregates shared with `07d_simulator_figures_counterfactual.R` (`load_sim_run_monthly()`, `summarise_reps()`, `paired_adm1_stats()`, `finish_adm1_ratio()`, `adm1_prob_stats()`, `finish_adm1_prob()`), and the group display helpers `group_display_meta()` / `relabel_groups()` (results.style.group_labels, `GROUP_LABEL_MODES`: `registry` = the sidecar's names, `lag` = class letters by the deployed occurrence lag, `centroid` = the sidecar's display-only `centroid_label` / `centroid_short` in `centroid_order`, e.g. "West [GW]"), and the fixed ISO3 -> WHO-AFRO subregion map `AFRO_SUBREGION` / `afro_subregion_of()` (11's `with_afro()`). | `base`, `paths`, `provenance`, `enso`, `loo_lags`, `gid`, `figures_cb`, `sim_io` (per-country censoring + manifest readers), `footprint` (`finest_level_totals()`, observed cases of the event attribution) (a `_core` kernel: sources its own closure) |
| `results_derive.R` | The manuscript-results derivers of phase 09 (`09b_postprocess_results.R`): `results_context(opt)` (results blocks, manifest, sets, thresholds, input paths, resolved once into a read-only environment, in 09's order), `results_dispatch(ctx)` over the step registry (`.RESULTS_STEP_TABLE` expanded by `results_steps()`, with `results_step_waves()` / `results_select_steps()` / `results_run_steps()` / `results_step_paths()` / `results_steps_table()`: one step per artefact group, track-keyed steps once per `fig1_tracks` entry, each with its dependency wave; 09a writes one task per step and 09b runs all, `--steps` or one `--task`), every `derive_*()` / `ensure_*()` pair (thin over the shared `ensure_artefact()` freshness/derive/write tail), `load_group_arm_runs()` (forecast + counterfactual partitions per group), `load_smoothing_scope()` / `smoothing_scope_gids()` and the report helpers (the per-draw population sums `pop_mass_by_draw()` / `dir_pop_ci()` over the draw array cached in `enso_relation`, the direction / magnitude / certainty exposure tables `derive_report_enso_exposure()` (directions via `exposure_mask()`: adverse / protective / either) and the superseded, unregistered ENSO-sensitive population `derive_report_sensitive_pop()` (kept side by side; the `either` direction reproduces it), both thin over the shared per-draw machinery `enso_draw_universe()` / `enso_track_masks()` / `enso_combine_masks()` / `enso_scopes()` / `mask_pop_rows()` / `prob_grid_rows()`, the incidence denominators `build_sim_scope_pop()` / `per_rep_totals()` / `window_n_months()`), and the forecast vs counterfactual exposed population `derive_report_attrib_exposure()` (Childs ratio-of-means point + paired bootstrap, IPCC robustness classes; pure core `attrib_exposure_tables()`) over helpers shared with attribution / report_forecast (`window_cum_per_rep()` / `drop_explosive_replicates()` / `explosive_replicate_keys()` with `EXPLOSIVE_CASE_CAP`, `pair_adm1_cum()`, `build_forecast_adm1_pop()`) and with the ENSO tables (`.point_scope_row()`, `mask_pop_rows(cond_cols =)`). Gate: `tests/test_results_attrib_exposure.R`. Beside it, the attribution metrics of the same paired runs `derive_report_attrib_metrics()` (attributable number / fraction, probability ratio with FAR / PS / PNS for `onset_any` and `cases_ge`, signal-to-noise with the emerged population; pure core `attrib_metrics_tables()` over the kernels `af_stats()` / `pr_stats()` / `sn_stats()` and `scope_pair_vectors()`; inputs from `attrib_arm_inputs()`, shared with `derive_report_attrib_exposure()`). Gate: `tests/test_results_attrib_metrics.R`. Also `derive_group_labels()` (each group's deployed lags read from the fit files via `paths.R::read_deployed_lag()`, asserted against the archive and the forecast manifest, plus table S30's `season_ends_country_table()` / `season_windows_table()`; gate `tests/test_group_display.R`), `report_lags` v2 rows `report_lags_selected_row()` over the deploy floor `lag_floor_of()` / `results_lag_floor()`, the lag-averaging step `derive_lag_averaging()` over the pure `lag_averaging_group()`, `lag_elpd_diff()`, `lag_krobust()`, `lag_weights()`, `lag_erf_logrr()`, `lag_common_grid()`, `lag_mixture_draws()` (gate `tests/test_lag_averaging.R`), and the membership step `derive_report_membership()` over `membership_map_partition()`, `membership_check_map()`, `membership_elpd_diff()`, `membership_sensitivity()`, `membership_cell_partitions()`, `membership_country_table()`, `membership_group_table()`, `membership_m2()` (gate `tests/test_results_membership.R`). Every function that reads the context takes `ctx` first. Not a `_core`: sources nothing. Confirmed-subset agreement: `derive_report_confirmed_agreement()` / `ensure_report_confirmed_agreement()` over the pure `confirmed_agreement_table()` (ENSO contrasts of occurrence, duration and monthly cases, all episodes vs those with a confirmed case; `.agree_rate_ratio()`, `.agree_boot_ratio()`, `.agree_boot_cluster()`) and `confirmed_characteristics_table()` (within-country ratios, confirmed / other, of cases, duration, attack rate and population), on the analysis set (the size dataset, ADM0 excluded, as report_dataset); method in `results.derivation.confirmed_agreement`. Smoothing v4: `smoothing_validation_tables()` (fig. SV's PIT histograms, seasonality, durations, active fraction and contagion conditionals; gate `tests/test_supp_missing_figs.R`). Gates: `tests/test_results_artefacts_equivalence.R` (every artefact identical to a reference baseline), `tests/test_results_derive_no_free_globals.R`, `tests/test_results_step_registry.R`. | the closure 09 sources: `results_core` (+ its own closure), `gid`, `sim_io`, `simulator_core`, `figures_cb`, `loo_lags`, `population`, `dur_post`, `dlnm`, `plots`, `clim_compile`, `clim_db`, `deploy_core`, `obs_sim_compare`, `sim_metrics`, `hindcast_skill` (the hindcast-skill section `ensure_hindcast()`), `anchor`, `diagnostics`, `season_calendar` (group_labels' season ends), `kl_sweep` (report_membership's `adjusted_rand_index()` / `kl_match_classes()`), `who_dashboard` (report_who_scale), `who_annual` (report_who_enso_ratio's `build_who_annual_cases()` / `who_window_totals()`) |
| `hindcast_skill.R` | Forecast-skill scoring of the spatPomp hindcast experiment, with fair variants: run index `hindcast_run_index()` (over `for_each_hindcast_run()`, refuses skipped runs and mixed episode keys), `hindcast_panel_paths()` (the runs' own episode-key panels), the per-run loader `hindcast_load_run()` (asserts the PIT cells equal the panel's at-risk cells), the occurrence scorers (`onset_window_counts()` at window / calendar-year grain, `score_occurrence_model()`, `score_occurrence_references()` over in-sample / leave-one-window-out / full-history / regional / zero / persistence references, `score_cells()`), the case scorers (`case_support()` over five supports, `score_cases_run()`, `score_case_references()`, `case_monthly_views()`), the notebook's metric primitives (`brier()`, `logscore()`, `auc_mw()`, `reliability_curve()`, `roc_points()`, `coverage_at_width()`, `randomised_rank()`, `interval_score()`, `crps_binom()`), and the orchestration `hindcast_skill_bundle()` / `hindcast_report_tables()` behind 09's `ensure_hindcast()`. Every ensemble score goes through one vectorised kernel (`sort_rows()`, `crps_sorted()`, `quantile7_sorted()`, `score_sorted()`), every case score through one unit builder (`case_sim_index()`, `case_unit_arrays()`); the replicate rows of a run are read once into (gid, month) cell x replicate matrices (`sim_cells()`, shared by `onset_window_counts()` and `case_sim_index()`), so no scorer passes over the replicate rows and every skill table through `.skill_rollup()` (optional season-year bootstrap `hindcast_boot_spec()` and per-area summaries); section 5b adds the monthly and season case skill (`case_obs_history()`, `case_skill_units()`, reference registry `.CASE_SKILL_REFERENCES()`). Gates: `tests/test_hindcast_skill_primitives.R`, `tests/test_hindcast_skill_golden.R` (the bundle against a golden capture), `tests/test_hindcast_skill_equivalence.R` (reproduces `notebooks/hindcast_skill.qmd`'s figure 9 and basis-B/C tables exactly). | `base`, `gid`, `obs_sim_compare`, `sim_io`, `sim_metrics`; *call-time: `config`, `paths`, `provenance`, `population` (`resolve_gid_population()`, case skill with `case_pop_scale: wpp`)* |
| `diagnostics.R` | MCMC convergence diagnostics of one saved fit: `fit_diagnostics(fit, chains, exclude, thresholds)` for a `brmsfit` (any backend; sampler table from `brms::nuts_params()`, E-BFMI from `energy__`, settings from the fit's `stan_args`, versions from `fit$version`) or a brms_kl `kl_fit` (sampler table from its `diagnostics`; a chain subset), returning per-parameter R-hat / bulk / tail ESS, per-chain divergences / max-treedepth hits / E-BFMI and scalar summaries against `DIAG_THRESHOLDS` (R-hat 1.01, ESS 400); `fit_diagnostics_row()`. The deploy convergence gate: `DIAG_SCOPES` (`all`; `population` = every parameter but `^(r_|z_|zs_)`), `DIAG_GATE` (R-hat 1.05 on every parameter, E-BFMI 0.3), `diag_scope_summary()` and the pure verdict `fit_diagnostics_gate()` (one row: scoped R-hat / ESS against `DIAG_THRESHOLDS`, the gross bound, divergences, E-BFMI; treedepth hits reported only; reuses `.diag_scalars()`; `fit_diagnostics()` output unchanged). Used by 09's `convergence` step, `kl_sweep.R`, `deploy_core.R` and `tests/test_diagnostics.R`. Gates: `tests/test_fit_diagnostics.R`, `tests/test_deploy_convergence_gate.R`. | `base`; `posterior` / `brms` / `tibble` via `::` |
| `kl_sweep.R` | Comparison of a brms_kl latent-class sweep behind `05_fit_postprocess_kl_sweep.R`: `kl_sweep_inventory()` (configs -> fits; K from kernel.classes under identify repel, else from the seeds, basis from argvar, prior from contiguity), `kl_chain_modes()` / `kl_agreeing_chains()` (the chain mode rule), `kl_fit_row()`, `kl_fit_diagnostics_tables()`, `kl_class_rows()`, `kl_partition_rows()`, membership uncertainty (`kl_responsibility_rows()` -> `responsibilities.csv`, `kl_country_mixture_rows()` -> `country_mixture.csv` (M1: per draw sum_k r_ck kcon_k / kmu_k vs the hard-MAP class), `kl_partition_stability()` -> `stability.csv` (partitions sampled from Cat(r_c), ARI vs the MAP), `kl_sweep_prior_grouping()`, pure `adjusted_rand_index()` / `kl_match_classes()`), `kl_sweep_check_comparable()`, `kl_sweep_loo()` (country-block LOO / WAIC, loo_compare, stacking weights) and the driver `kl_sweep_compare()`. Gate: `tests/test_kl_sweep_helpers.R` (Stan-free fake fits). | `base`, `config`, `paths`, `dlnm_kernel`, `diagnostics`; `dplyr` / `posterior` / `loo` via `::` |
| `kl_grouping.R` | The fit registry from a brms_kl latent-class fit: `kl_grouping_from_fit()` (argmax partition + membership, chain subset, `fill_from`; the `make_kl_grouping.R` body), `kl_grouping_sidecar()` (labels / colours / order from `meta_from` for shared groups, an unused fallback colour for new ones), `kl_grouping_status()` (absent / identical / differs + per-country diff), `kl_grouping_header()` (provenance lines), `write_grouping_registry()` (CSV + `.yml` + membership CSV, read back) and `kl_grouping_plan_from_sweep()` (the pinned cell of a 05 sweep comparison: cell + agreeing chains from `loo.rds`, `fits.groups` check), and `kl_grouping_convergence()` (the pinned cell's all-chain / agreeing-chain diagnostics from `fits.csv` + `diagnostics_chains.csv`; under `deploy.convergence.action: refuse` it refuses when fewer than `KL_GROUPING_MIN_AGREE` (2) chains agree or the agreeing R-hat exceeds 1.05, else warns beyond 1.01 / 400). The sidecar carries the centroid display names of `style.group_labels: centroid` (`kl_centroid_fields()`: the registry's previous sidecar, then `meta_from`, then `KL_CENTROID_NAMES` by id). Used by `deploy_prod_run.R --stage gen_grouping` and `make_kl_grouping.R`. Gates: `tests/test_deploy_kl_grouping.R`, `tests/test_deploy_convergence_gate.R`. | `base`, `config`, `dlnm_kernel` (+ `dlnm`), `diagnostics` (`kl_grouping_convergence()` only); `dplyr` attached; `tidyr` / `purrr` / `yaml` / `digest` via `::` |
| `figures_kl.R` | ggplot builders of `06_fit_figures_kl_sweep.R`: batch figures from the 05 tables (`plot_kl_sweep_elpd()`, `_weights()`, `_centres()`, `_convergence()`, `_partition()`) and the per-fit panels (`plot_kl_traces()`, `plot_kl_centres_by_chain()`, `plot_kl_membership_maps()`, `plot_kl_class_erf()`, `plot_kl_re_sds()`); never save. | `base`, `config`, `dlnm_kernel`, `kl_sweep`, `figures_cb` (basemap); `ggplot2` / `dplyr` / `tidyr` / `sf` via `::` |
| `figures_supp.R` | ggplot builders of supplementary figures S5 (`build_figS5_convergence()`, sampler diagnostics of the deployed models), S6 (`build_figS6_lag_loo()`, LOO lag rankings), the by-country forecast / smoothing figures (`build_figS_forecast_country()`, `build_figS_smoothing_country()`) AE (`build_figS_attrib_exposure()`, ratio-of-means map with IPCC hatching + population by class; map layers and hatch fills passed in; its map built by `attrib_ratio_map_sf()` / `attrib_ratio_map_layers()`, shared with Figure 4b), and MB (`build_figS_membership()`, responsibility matrix + class paths across the sweep cells) for `10_manuscript_figures.R`, from phase-09 artefacts; never save. The supplement maps `build_figS_modelled_units()` (unit status from 09's `modelled_units`; display labels / colours `MODELLED_UNIT_STATUS` / `MODELLED_UNIT_COL`), `build_figS_group_map()`, `kl_classes_by_K_frame()` / `build_figS_kl_classes_by_K()` (10 `--figures MU,GM,KC`; gate `tests/test_report11_supplement_tables.R`). Further supplement figures (10 `--figures DC,SV,FA,AA`; gate `tests/test_supp_missing_figs.R`): `data_coverage_frames()` / `build_figS_data_coverage()` (report_dataset's `outbreaks`), `pit_hist_rebin()` / `build_figS_smoothing_validation()` (smoothing v4 `validation`), `forecast_analogue_frame()` / `build_figS_forecast_analogues()` (attribution ribbons), `attribution_analogue_map_sf()` / `build_figS_attribution_analogues()` (report_attrib_exposure's `adm1`). The processing cascade as a flow diagram (10 `--figures PC`; same gate): `processing_cascade_frames()` / `build_figS_processing_cascade()` (report_dataset's `cascade`). | `ggplot2` / `dplyr` / `tidyr` / `patchwork` / `sf` / `scales` / `grid` / `grDevices` via `::` |
| `figures_hindcast.R` | ggplot builders for the hindcast-skill artefacts (`hc_plot_*()`, `hc_scope_labels()`, `hc_scope_palette()`): skill bars by reference, skill by lead, cell scores, case support, per-country lollipop, case intervals, monthly trajectories, PIT / rank histograms, reliability, ROC, coverage curve, interval score, seasonal cycle. Take tidy artefact tables plus labels / palette / theme; never save. Shared by `10_manuscript_figures.R` (figure H) and `notebooks/hindcast_skill.qmd`. | *call-time: `ggplot2`, `dplyr`, `tidyr`, `scales` via `::`* |
| `figures_manuscript_fig1.R` | Figure 1 builders of `10_manuscript_figures.R` (`build_fig1()`, the coverage panel + inset map, the lag-weight / rainfall-correlation / onset-season panel, the ERF grid, their legend block), 10's Section 9a. **10-only**: definitions load alone, but the functions read 10's globals when called (style `FONT` / `BASE_SIZE` / `theme_ms()` / `PS_*`, `GROUP_*`, `RES` / `opt`, the `FIG1_*` constants, `read_results()` and 10's other helpers). | *call-time: 10's globals; `figures_cb`, `figures_io`, `figures_manuscript_fig34` (`fig3_oni_strip()`), `paths`, `provenance`* |
| `figures_manuscript_fig2.R` | Figure 2 builders of `10_manuscript_figures.R` (`build_fig2()`, the bar panel, the gtable edits: header rules, phase fills, row keys, bar rows, tag, height fit), 10's Section 9b. **10-only**, as above (`FIG2_*` constants). | *call-time: 10's globals* |
| `figures_manuscript_fig34.R` | Figure 3 and 4 builders of `10_manuscript_figures.R` (`build_fig3()`, `build_fig4()`, the forecast panels `fig3_panel_a/b/c()` and ONI strip `fig3_oni_strip()` they share, the event and ratio / AF panels, the legends), 10's Sections 9c and 9d. **10-only**, as above (`F3_*`, `F3_EV_*`, `FIG4_*` constants). | *call-time: 10's globals; `figures_cb`, `figures_io`, `figures_supp` (`attrib_ratio_map_sf()` / `_layers()`), `figures_manuscript_fig2` (`fig2_group_map_inset()`, the group-map insets of Figures 3b / 4c)* |
| `results_report.R` | Display helpers of the phase-11 key-numbers report (`11_manuscript_numbers.qmd`): the `fmt_*()` formatters (incl. `fmt_rate()` and `fmt_by_quantity()`), `short_enso()`, `scenario_lab()`, `report_levels()` (region / ONI-condition / scope / scenario orderings from the results config), `QUANTITY_LABELS`, `report_wide()` (the long-to-wide table shaping), `report_with_africa()`, `report_table_counter()` / `report_emit_html()` (kable emission) and `report_loader()` (`results_require()` of any `results_*` kind); the supplement helpers: `with_afro()` (WHO-AFRO subregion column), `report_data_id()` / `report_release_split()` (the `style.data_release` split into `supp_data/` / `github/` with its manifest), `report_run_design_table()` (D34), `report_analysis_set_table()` (D35), `report_t1a_mechanisms()` (from `rain_enso$mechanisms`), `report_t1b_rain_lags()` (09's `rain_enso$group_lags` + the deployed lags looked up; T2 is `report_lags$loo_summary`, printed as is), `report_outside_range()` / `RANGE_MARK` (T1-T3; gate `tests/test_report11_supplement_tables.R`). Formats only; never computes a quantile. Gate: `tests/test_results_sensitive_pop.R` covers the 09 helpers it displays. | call time: `results_core` registries (`GROUP_SHORT`, `GROUP_ORDER`, `RESULTS_SCHEMA`, `make_results_run_file()`, `results_require()`, `afro_subregion_of()`, `group_lag_signed()`), `figures_cb::enso_conds_named()`; `knitr` / `kableExtra` / `tidyr` via `::` |
| `deploy_core.R` | Domain logic for the deploy orchestrator (`deploy_prod_run.R`, §5 tools below): best-lag fit selection (`select_best_lag` via `loo_lags`), canonical artefact naming, fit collection/`MANIFEST` into `generated_data/enso_best_lag_fits/`, deploy-vs-results config comparison/alignment (`compare_deploy_results()`, `deploy_plan_results_alignment()` / `align_results_to_deploy()` — a line-level text edit of `manuscript.yml` that keeps its comments and verifies by re-read; `resolve_results_inputs_from_deploy()` refuses an explicit `inputs.*` leaf that disagrees with the deploy unless `results.allow_overrides` names it), and provenance stamping. `deploy.episodes.key` is the rollout's single episode declaration: `deploy_collect_fits()` refuses per-lag configs of another key, checks (gid_season) every source fit's episode record against the prep artefacts BEFORE pruning/copying, and never prunes the other definition's fits (`.deploy_foreign_episode_artefacts()`); the generators write the key into every config via `provenance.R::apply_episode_key_to_sim_config()` and pre-flight its inputs with `preflight_sim_episode_inputs()`. **Sim-template generation** (`--stage gen_templates`): `deploy_build_sim_template()` (skeleton `analysis/configs/sim_templates/<arm>.yml` + deploy `sim.shared` / `sim.<arm>` dotted-leaf overlay + `.deploy_sim_derived_leaves()`), `deploy_plan_sim_templates()` (absent / unchanged / differs, ignoring `simulator.run_id`), `deploy_write_sim_template()`, `deploy_check_sim_template()` (the on-disk template agrees with the deploy: run_tag, nbr gate, episode key, fits dir, forward `smoothing_input`), `deploy_sim_run_id()`, `deploy_smoothing_dir()` (the tracks are occ/dur/mon); and the generic `deploy_leaf_diff()` / `deploy_print_leaf_diff()` (also behind `compare_deploy_results()`). **Convergence gate**: `deploy_convergence_cfg()` (deploy.convergence: action off / warn / refuse, scope, the `"default"` threshold sentinels resolved from `diagnostics.R::DIAG_THRESHOLDS`, max_divergent, `allow` ids `<track>:<short group>`; a missing block = off), `deploy_check_convergence()` (every source fit through `fit_diagnostics_gate()` BEFORE the prune; refuse stops, `--dry-run` prints "WOULD refuse"; verdicts to `dest/CONVERGENCE.csv`), `.deploy_read_source_fits()` (one `readRDS` per source fit, shared with the episode check). `deploy_assert_sim_pop_year()`: the three arms must resolve one `simulator.location_frame.pop_year` (in `deploy_plan_sim_templates()` and, against the historical template on disk, in `deploy_check_sim_template()`). `.deploy_leaf_same()` is typed (a deploy-derived number must be a YAML number, a logical a boolean; also behind the post-align re-read and `compare_deploy_results()`; `deploy_leaf_diff()` stays string-based), `.deploy_leaf_type_hint()` gives the refusal's unquoting hint, `.yaml_set_scalar_line()` writes numbers bare; `deploy_select_lags()` checks `occ_floor_offset` against the duration panel's stamped episode span (`loo_lags.R::check_occ_floor_offset()`, `provenance.R::read_prep_policy_stamp()`). | `base`, `config`, `paths`, `provenance`, `loo_lags`, `diagnostics`, `simulator_core` (`sim_pop_year_spec()`) |
| `formula.R` | Model-term assembly for all four tracks: `build_formula_terms()` (occurrence), `build_{size,duration,monthly}_formula_parts()` + `build_*_formula()`, the RE/fixed-term helpers the `_cb` fit scripts also call (`re_group_by_level()`, `nested_slope_re_terms()`, `filter_re_terms_for_data()`, `filter_fixed_factor_terms_for_data()`), slope resolution (`validate_size_options()`, `resolve_slopes()`), the offset spec + base-rate anchor gate (`.OFFSET_SPEC`, `.offset_terms()`, `.offset_spec()`, `.anchor_gate_on()`), `.occ_basehaz()`, the occurrence baseline term `.occ_baseline_block()` and `integrate_cb_terms()` (dlnm basis columns into the formula terms; shared by `04_fit_occ_brms_cb.R` and `04_fit_occ_brms_kl.R`). Golden gate: `tests/test_formula_builder_golden.R`. | `base`, `paths` *(call-time: `nbr_covariate_name()`; `.leaf_scope_tag()` only when the leaf scope is not "all")*, `config` *(call-time: `get_default_options()`, only when a size / monthly block has no backend)* |
| `formula_legacy.R` | `build_formula_terms_legacy()` (+ its `.occ_group_by_level()`), the frozen pre-`resolve_slopes()` occurrence builder, kept side-by-side; only `tests/test_occurrence_config_merge.R` calls it. | `base` |
| `simulator_core_legacy.R` | `simulate_one_tu_legacy()` and the nine helpers it calls (`evaluate_lp`, `update_lags`, `advance_ttlo`, `bin_ttlo`, `.dur_bin_mso`, `basehaz_dur_eval`, `.onset_log_offset_adj`, `build_warmup_state`, `monthly_log_mu`, each `_legacy`), the frozen pre-speed-up forward kernel, kept side by side and never modified; only `tests/test_simulator_kernel_golden.R` calls it. | `simulator_core` (+ its closure) |
| `priors.R` | brms prior construction for all four tracks: `build_occ_priors()`, `build_duration_priors()`, `build_size_priors()`, `build_monthly_priors()`, the shared `.PRIOR_DEFAULTS` table (also read by `make_batch_configs.R`) and the tiered-`sd` helpers (`.re_group_tiers()`, `.tiered_sd_priors()`). `build_occ_priors()` / `.build_hazard_track_priors()` / `.tiered_sd_priors()` take an optional `nlpar` (brms_kl: `"lin"`); `NULL`, every other caller, is byte-identical. Golden gate: `tests/test_prior_builder_golden.R`. | `base`, `formula` (`build_size_formula_parts()` / `build_monthly_formula_parts()` when `parts` is NULL, `resolve_slopes()`, `.anchor_gate_on()`, `.occ_basehaz()`) |
| `who_annual.R` | WHO annual cholera reports: `build_who_annual_zero_mask()` combines the WHO annual report (`sCh`, 1970–2024) and the WHO dashboard annual totals (2023 onwards) from `data/who_annual_report_africa/` into one row per (country ISO3, year) with `suspected_cases_max` (max of the sources, tested by the national-zero threshold) and `suspected_cases` (mean of the two sources where both report) and `is_zero`. Consumed by `01_preprocess_national_zero_months.R`. `build_who_annual_cases()`: annual cases per country-year for comparing years, absence NA (never zero), the dashboard weekly series summed by epiweek year, report + dashboard averaged, a target-year series stopping before the cut-off censored; `who_window_totals()`: cases over a set of calendar months from the WHO monthly table, complete windows only. Both consumed by 09's `report_who_enso_ratio` (supplement fig. WR). | `dplyr`, `tidyr`, `tibble` |
| `who_dashboard.R` | WHO dashboard national weekly counts (`data/who_annual_report_africa/who_dashboard_weekly_africa.rds`, 2023 onwards) and the year-to-date cumulative snapshots of the WHO multi-country updates (`who_sitrep_cumulative_africa.csv`): `who_weekly_to_monthly()` (day-overlap split into calendar months; unreported weeks are gaps, `partial` flags incomplete months), `who_cumulative_to_monthly()` (cuts differenced within the year, same split; an unchanged cumulative is a zero only for `zero_iso3` or an explicit last-28-days 0), `who_merge_monthly()` (cumulative rows fill only months with no weekly row; `source` column), `who_obs_table()` (complete months in a window), `who_coverage_fraction()` (covered population share of each country vs WPP), `who_block_share()` (per-unit share of the country term). Run by `03_prepare_who_monthly.R`; consumed by `smoothing_core.R` for the `_who` measurement (`simulator.filter.who_obs`). | `dplyr` |
| `dlnm_kernel.R` | Gaussian ENSO lag kernel and latent country classes for the occurrence `brms_kl` backend (`04_fit_occ_brms_kl.R`, sandbox): `resolve_kernel_cfg()` / `kl_engine()` / `kl_spatial_prior()` (validated `occurrence.kernel` block), the brms non-linear kernel formula and priors (`build_kl_expr()`, `build_kl_formula()`, `build_kl_priors()` over `priors.R::build_occ_priors(nlpar = "lin")`, `assemble_kl_model()` -- basis from `dlnm.R::build_cb_block()` `lags:` mode, `lin` terms via `formula.R::build_formula_terms()` + `integrate_cb_terms()`), the patcher `patch_kl_latent_class_stancode()` (brms's program for that formula turned into a vectorised K-class country mixture; landmarks `kl_stan_landmarks()`, each asserted once) with the Stan fragments `.KL_STAN_*` (incl. a DAGAR field), data/inits/fit (`kl_prior_partition()`, `kl_adjacency_data()`, `build_kl_latent_standata()`, `kl_sort_by_country()`, `make_kl_class_inits()`, `fit_kl_stan()` -> class `kl_fit`), and summaries (`kl_weights()`, `summarise_kernel_draws()`, `summarise_kl_fit()`, `summarise_kl_latent()`, `kl_responsibilities()`, `kl_comembership()`, `kl_chain_agreement()`, `kl_map_partition()`, `write_kl_partition()`, `kl_loo_country()`). Two class forms (`kernel.form`): `gaussian` (kernel centre/width x exposure curve; identification `kernel.identify` kmu / anchor / erf via `kl_identify()`, `kl_erf_direction()`, `kl_name_classes()`, `kl_init_strategy()` / `kl_multistart_optimize()`) and `dlnm` (B-spline lag basis `kl_lag_basis()` x exposure basis per class under a regularised group horseshoe over lag windows; exchangeable classes relabelled after sampling by `kl_relabel_draws()`; summaries `summarise_kl_dlnm_class()`); a third form, `discrete`, puts the Gaussian kernel centres on the lag grid (widths from `kernel.discrete.widths`, shared, summed out) with class centres ordered, ties allowed, and summed out exactly over the choose(NC + K - 1, K) combinations (`kl_disc_combos()`, `kl_disc_spec()`; identification `lagorder`, ties ordered by kcon in `kl_relabel_draws(lag = "klag")`; drawn centres exposed as `kmu` by `kl_disc_rename_draws()`); `kernel.pin` (form discrete) fixes which class sits at which lag rank by seed countries (`kl_pin_seeds()`, `kl_pin_fix()`, Stan `kl_fix`), one run per candidate order of the prior groups (`kl_pin_arrangements()`, from per-group pilot fits `kl_pin_pilots()` or an existing sweep's LOO rankings `kl_pin_sweep_orders()` via `loo_lags.R::read_latest_ranking()`), runs in parallel (`kl_pin_fit_arrangements()`) and combined by stacking on country-level LOO (`kl_pin_combine()`; driver `kl_pin_run()`, or the SLURM split `kl_pin_plan()` / `kl_pin_fit_one()` / `kl_pin_finish()` behind `04_fit_occ_brms_kl.R --pin_stage`); `kernel.re_slope` (forms discrete / gaussian) adds country and ADM1 deviations of the class exposure coefficients (Stan `kl_disc_re_term()` / `kl_ctry_logp()`, the analogue of brms_cb's dlnm `re_adm1_slope`). Gaussian form extras: `identify: seed` (seed countries fix class k = prior group k, centres unordered, one run; `kl_pin_seeds()` / `kl_pin_fix()`), `ksd_fixed`, `init: sweep` (`kl_sweep_modes()`), and `reduce_sum` over countries (`kl_partial_ctry()`) when `brms_kl.threads_per_chain > 1` (country path `kl_rs`; the discrete form cannot be threaded: its summed-out lag combinations couple the countries). MCMC tests: synthetic cases A-J in `tests/test_kl_synthetic_recovery.R`. Gates: `tests/test_kl_stan_patch.R` (+ pinned program `tests/golden/kl_latent_program.stan`), `tests/test_kl_lp_equivalence.R` (log density equals brms's at K = 1; form discrete equals a labelled brute-force sum); MCMC: `tests/test_kl_synthetic_recovery.R`, `tests/test_kl_smoke_e2e.R`. | `base`, `dlnm`, `formula`, `priors`, `config` (`load_grouping()`), `paths` (`read_deployed_lag()`, only with `kernel.prior_grouping.deploy_config`), `loo_lags` (only with `kernel.pin.order_source: sweep`); `brms` / `cmdstanr` / `posterior` / `dlnm` / `loo` / `digest` via `::` |
| `national_zero.R` | National-zero months for the surveillance-coverage mask: `national_month_status()` (national weekly series → per-month zero / positive), `subnational_positive_months()`, `active_episode_months()`, `parse_national_zero_exclude()`, `db_who_consistency()` (database deepest-level, year-apportioned episode cases vs the WHO total; flags under-reported WHO years), `build_national_zero_months()` (with `zero_threshold` > 0 also WHO country-years under the threshold, `who_below_threshold`, last `window` months only if the next year is also below) — the rule `((national 0 AND no WHO > 0) OR (WHO 0 AND no national > 0)) AND no sub-national cases AND no active outbreak AND not excluded`. Run by `01_preprocess_national_zero_months.R`; its output is read by `sim_io.R::build_coverage_mask()` (on by default). | `base` (`%||%`); `who_annual` output as input |
| `season_calendar.R` | Per-country cholera-season ends from the pipeline's own seasonality: `season_end_month()` (min_seas − 1, 1 → 12), `season_end_by_country()` (gid → country: modal min_seas, ties → the later end month; modal share, n_gid, n_distinct), `season_end_span()` (circular spread), `first_season_end()`, the window rules `event_window_end()` (first end on/after the peak + max positive lag; `make_event_attribution_configs.R::event_group_window()` delegates to it) and `hindcast_window_end()` (first end ≥ launch + `min_window_months` + max(min lag, 0)), `panel_season_calendar()` / `read_panel_season_calendar()` / `season_end_source()` (a sim panel's gid calendar + digest), `season_end_spread_warnings()` (the by_group > 4-month diagnostic), `season_ends_csv_path()` / `write_season_ends_csv()`, and `gid_min_seas_source()` (each gid's canonical min_seas and its source, estimate vs January default, for `03_prepare_occ.R`'s `.min_seas_source.csv` sidecar). Sourced by `make_hindcast_configs.R`, `make_event_attribution_configs.R`, `03_prepare_occ.R`, `09b_postprocess_results.R` (group_labels' table S30). Gate: `tests/test_season_calendar.R`. | `base` (`%||%`); `digest` optional |
| `enso_continuation_skill.R` | Skill by lead of the ENSO continuation past the CFSv2 span on observed seams (perfect seam), for `03_prepare_enso_continuation_library.R --validate` -> `generated_data/enso_continuation_skill.csv` (`paths.R::make_enso_continuation_skill_file()`), which 09 reads for each group's continuation CRPSS: seams `enso_skill_seams()` / `enso_skill_seam_state()`, method closures `enso_skill_method_neutral()` / `enso_skill_method_ar1()` / `enso_skill_method_clim()` (leave-out same-calendar-month climatology) / `enso_skill_method_analog()` (the deployed `continue_enso_analog()` with every `simulator.enso.continuation` knob), `enso_continuation_library_rehorizon()` (candidate rule for another horizon), the scorer `enso_continuation_hindcast()` (fair CRPS per seam x method x lead, one RNG stream per seed), `enso_skill_by_lead()` / `compute_enso_continuation_skill()` (CRPS, RMSE, spread, CRPSS vs climatology / neutral; seed means + seed range) and `enso_continuation_skill_settings()` (defaults' continuation block overlaid by a deploy YAML's `deploy.enso_continuation`). Gate: `tests/test_enso_continuation_skill.R`. | `base`, `enso` (`continue_enso_analog()`, `enso_continuation_features()`), `obs_sim_compare` (`crps_fair()`); `yaml` via `::` |

**`config.R`, `formula.R`, `formula_legacy.R`, `priors.R`.** `config.R` holds the
options core (defaults, YAML merge, CLI routes, grouping, dotted-leaf and
simulator-CLI helpers); the formula builders are in `formula.R`, the legacy
occurrence builder in `formula_legacy.R` and the prior builders in `priors.R`.
`%||%` lives in `base.R` only, so callers source `base.R` first and the formula /
prior modules next to `config.R`. Two HARD gates pin the layout:
`tests/test_config_split_symbols.R` (module → symbol ownership) and
`tests/test_utils_module_purity.R` (no utils module `source()`s another outside the
documented exceptions — the `_core` kernels, `spatpomp_core_c.R`, the CLI tools —
and the four config modules each load alone after `base.R`).

**07c / 07d.** `07c_simulator_postprocess.R` is a thin orchestrator over
`sim_metrics.R`; `aggregate_partitions()` / `combine_group_runs()` are in
`sim_io.R`. Shared helpers: `paths.R::make_sim_engine_dirs()` (the `--engine`
run-dir dispatch of 07c and 07d) and `partition_sentinels()` (the `_SUCCESS` scan of
`aggregate_partitions()` and `results_core.R::load_sim_run_monthly()`);
`sim_io.R::sim_metrics_fmt()` (the metrics-format rule) and
`read_hindcast_manifest()` / `for_each_hindcast_run()` (the hindcast loop of 07c
`--hindcast` and 07d `--hindcast`); `obs_sim_compare.R::crps_ens()` /
`region_key_of()`; and `.sim_frame_geo()` / `.geo_crps_level()` for the geography
scoring of both CRPS scorers. Gates: `tests/test_sim_metrics_golden.R`
(`run_validation()` end to end under both total-case bases, partition aggregation,
`compute_hindcast_series()`) and `tests/test_hindcast_iterator.R` (the hindcast
skip paths).

**spatPomp kernels and batch tools.** The builders share `spatpomp_core.R`
helpers (`tests/test_spatpomp_builder_golden.R`, seven R/C variants incl.
generated C); `.state_vec_to_warmup_row()` is behind `filtered_states_to_warmup()`
and `filter_traj_to_warmup()` (only the cache path checks the active/mso schema);
`spatpomp_core_c.R::.stitch_grouped_lp()` behind the three grouped precompute
wrappers; `.cache_handoff_label()` behind 07b's two dry-run branches;
`smoothing_core.R`'s event-table and ribbon families are one implementation each
(`tests/test_smoothing_aggregates_golden.R`); `paths.R::partition_sentinels()` is
used by 08c/02c; `paths.R::.FIT_FILE_STEM` is behind both `.DEPLOYED_FIT_PREFIX`
and `simulations.R::.discover_group_fit()`. Kept separate on purpose: the
builders' draw-subset/duration/monthly blocks (divergent), `run_bpfilter_tu()` vs
`run_one_bpfilter_replicate()` (different contracts), `smoothing_load_artefacts()`
vs `simulations.R::load_artefacts()` (different backend defaulting),
`.load_cb_sidecar()` vs the closure inside `load_artefacts()`, and the two
sentinel writers (different contents).

**09.** `09b_postprocess_results.R` is a guarded CLI: `results_context(opt)` then
`results_dispatch(ctx)`. Its functions live in `results_derive.R` and take `ctx`
first; `ctx` is an environment, so `ctx$MAN` can never partial-match. Shared
pieces: `ensure_artefact()` (the read-if-fresh / derive / write tail of the
`ensure_*()` wrappers); `load_group_arm_runs()` (the forecast/counterfactual load
loop of the attribution and report-forecast derivers); `load_smoothing_scope()` /
`smoothing_scope_gids()` (the set-up of both smoothing derivers);
`ribbon_from_events()` reuses `rep_month_totals()`; `base.R::wquantile()`,
`figures_cb.R::adm1_logrr_draws()` (the ADM1 draw combination of
`extract_adm1_irr()` and `dir_pop_ci()`) and `paths.R::partition_sentinels()` (the
smoothing-parts scan) are the single copies.

**Shared helpers in existing modules.**
`config.R` has `config_get_leaf()` / `config_set_leaf()` / `config_leaf_paths()` /
`config_flatten_leaves()`, the dotted-path nested-list helpers (used by
`make_batch_configs.R`, `deploy_core.R` and the two coverage tests).
`provenance.R` has `generated_yaml_header()` / `write_yaml_with_header()`, the
"GENERATED ... DO NOT EDIT BY HAND" banner and writer that every generated
config uses (forecast, hindcast and deploy sim templates).


A **fourth standalone CLI tool**, `make_hindcast_configs.R`, is the hindcast counterpart of §3's `make_forecast_configs.R`
(same `#!/usr/bin/env Rscript` guard + `optparse` entrypoint, not sourced-only). It
generates per-(ENSO group, launch year) hindcast configs that pin the smoothing handoff
to an interior historical month and drive observed-ENSO + historical-climate in-sample
forecast-skill scoring; writes `configs/hindcast_by_group*/*.yml` + `configs/hindcast_horizons.yml`.

A further standalone CLI tool, `make_event_attribution_configs.R`, builds
the past-ENSO-event attribution runs: per event of `analysis/configs/enso_events_v2.yml`
and simulation unit, a factual config (a hindcast, built by
`make_hindcast_configs.R::patch_hindcast_config()`)
and a counterfactual one (the same plus `simulator.enso.counterfactual`, Lenssen
analogs from a per-event cache), and the event manifest
(`configs/event_attribution_horizons_run<N>.yml`) read by 09's `report_event_attrib`.
It sources `make_hindcast_configs.R` (CLI guarded) and is run in-process by
`deploy_prod_run.R --stage gen_events`. The smoothing-cache filename parser both
generators use is `paths.R::smoothing_cache_coverage()`.

Both generators take `--season_end by_group | by_country_from_panel` and
`--season_min_window_months` (default 6; hindcast only). The default `by_group`
reads `--season_end_by_group` and warns when a group's member countries' cholera-year ends span more
than 4 months. `by_country_from_panel` takes each country's end from the gid
calendar of `--climate_panel` (`season_calendar.R`), writes per-country firm horizons
under `countries.<iso>` (hindcasts: `years.<y>`), a group `horizon_end` = max over its
members, a unit horizon = its latest member country, `season_end_source` (panel +
min_seas digest), and `season_ends_<manifest stem>.csv` next to the manifest.

Another standalone CLI tool, `make_forecast_attrib_configs.R`, builds the
attribution factual of a forecast rollout: per simulation unit, the factual config named
by the forecast manifest launched at the counterfactual handoff (only `horizon.start`,
`filter.forecast_window.start`, `filter.handoff_month`, `run_tag` and `run_id` change;
written as `forecast_<unit>_attrib.yml` next to the factual configs), and an attribution
manifest (the forecast manifest with its factual fields swapped) that 09's attribution
steps read through `results.inputs.attrib_manifest`. It sources `base.R`, `config.R`,
`paths.R`, `sim_io.R` and `provenance.R`, and is run in-process by
`deploy_prod_run.R --stage gen_attrib` (also by `bash/tools/run_pipeline.sh`'s sim and
attrib phases).

A **fifth standalone CLI tool**, `deploy_prod_run.R`, is the `prod_run<N>` deploy
orchestrator.
It wires configs and collects fits and transforms no data, and it runs at several
points between phases 05 and 07, so it is a non-step tool rather than a phase.
From one `configs/deploy/<stem>.yml` it collects and canonically names the best-lag
fits (`deploy_core.R`), generates forecast and hindcast configs **in-process** by
`source()`-ing the two generators above (their `sys.nframe()==0L` CLI guards make
that safe), and **prints** the SLURM commands, but never submits them and never
shells out. Stages are selected with `--stage` (`gen_batches` | `deploy_fits` |
`gen_templates` | `smoothing` | `gen_forecast` | `gen_hindcast` | `gen_events` | `align_results` |
`all` | `resolve`); `--dry-run`, `--align-results` and `--force` modify them.
`align_results` (also run by `--align-results` after `gen_forecast` /
`gen_hindcast`) points `analysis/configs/results/manuscript.yml` at the deploy:
`results.deploy_config` plus the six deploy-derived `inputs.*` leaves, as a
line-level **text** edit (`deploy_core.R::align_results_to_deploy()`) — every
comment kept, the result re-read and checked all-match before it reports success,
the pre-edit copy under `generated_data/deploy/<stem>/`, and comment lines still
naming another run listed for a hand edit; `--dry-run` prints the line diff. Its
counterpart on the read side is the guard in `resolve_results_inputs_from_deploy()`:
with `deploy_config` set, an explicit `inputs.*` leaf that disagrees with the
deploy is refused unless named under `results.allow_overrides`. `all` runs
`deploy_fits -> gen_templates -> smoothing`. `gen_templates`
generates the run's three spatPomp templates (`deploy.templates.*`) from the
skeletons in `analysis/configs/sim_templates/` and the deploy YAML's `sim:`
overlay. An existing template that matches is left untouched; one that differs
is refused with a leaf diff unless `--force`. `gen_batches` (run before fitting, so not part of
`all`) generates the fit-batch config dirs named by the YAML's `fits.*_dir` via
`make_batch_configs.R::generate_batch()`, loaded into an isolated environment
(`deploy_core.R::deploy_batch_env()`; that file's CLI guard is `sys.nframe() == 0L`
like the other generators); a dir that already holds a `manifest.yml` is skipped.
`resolve` prints the run's paths as `KEY=value` lines
(`deploy_core.R::deploy_resolved_keys()`), which every bash driver taking `--deploy`
parses; the deploy YAML is the only place a `prod_run<N>` is named, and
`deploy_require_keys()` refuses a YAML that leaves a run-specific key unset (their
defaults are `""`). Usage: `Rscript analysis/utils/deploy_prod_run.R -C analysis/configs/deploy/prod_run.yml --stage all --dry-run`.
Executed by `bash/tools/run_pipeline.sh` (stage `deploy_fits` and the others) and by
the `--deploy` drivers (stage `resolve`; see `bash/README.md`).

**brms_kl sandbox additions to the §1 modules** (behaviour of every
other caller unchanged): `occ_data.R::prepare_occ_brms_data()` (the occurrence brms
workers' shared data preamble), `dlnm.R` `lags:` mode in `build_cb_block()` /
`build_cb_columns_from_Q()` / `build_cb_design()` (one shared onebasis replayed
per lag; helpers `.cb_lag_exposure()`, `.cb_lags_spec()`, `.cb_lag_col_names()`,
`.replay_onebasis_lags()`), and in `config.R` `resolve_brms_init()` /
`resolve_brms_threads()` (occurrence workers' init and threading blocks),
`load_adjacency()` over the registry reader `.read_registry_csv()` shared with
`load_grouping()`, and the `brms_kl` backend in `.validate_basehaz()`,
`.fill_model_block()` (`_kl` inherits `brms` like `_cb`; `kernel:` filled
recursively) and `.CLI_ROUTE`.

A **seventh standalone CLI tool**, `make_country_adjacency.R` (brms_kl sandbox),
writes the country adjacency registry `analysis/groupings/country_adjacency.csv`
(`iso3_a,iso3_b`) from the `taxdat::afr_sf` polygons (`--tol_km`, default 1 km),
read back by `config.R::load_adjacency()` for the `brms_kl` spatial membership
term. Depends on nothing in `utils/`.

A further standalone CLI tool, `kl_compile_batch.R`, compiles
every Stan program of a brms_kl batch once before its fit array: it runs
`04_fit_run.R -C <config> --compile_only` per config in parallel (configs that
share a program meet at `dlnm_kernel.R::kl_with_compile_lock()`; the first
compiles, the others reuse the executable). Submitted by
`bash/stages/04_fit_kl.sh` via `slurm/submit_kl_compile.sh`, the fit
array afterok it. Sources nothing in `utils/` (it shells out to the fit script).

An **eighth standalone CLI tool**, `make_kl_grouping.R` (brms_kl sandbox),
writes a groupings registry `analysis/groupings/<out>.csv` from a saved
`brms_kl` latent-class fit: each fitted country goes to its class of maximum
mean posterior probability (`dlnm_kernel.R::kl_map_partition()`), optionally on
a subset of chains (`kl_subset_chains()`; `--write_fit` saves that subset beside
the fit) and with `--fill_from` carrying over a registry's unfitted countries.
Per-country probabilities go to `generated_data/kl_groupings/<out>_membership.csv`.
Sources `base`, `config`, `paths`, `dlnm`, `dlnm_kernel`, `kl_grouping` (its body is
`kl_grouping.R::kl_grouping_from_fit()`, shared with `deploy_prod_run.R --stage gen_grouping`,
which derives a deploy's registry, with its `.yml` sidecar, from a pinned sweep cell).

A **ninth standalone CLI tool**, `make_kl_sweep_groupings.R`, writes the prior
registries of the brms_kl class sweep, `analysis/groupings/enso_grouping_k{2,3,4}.csv`
(+ `.yml` sidecars), as whole-group merges of `enso_grouping` (K = 5 is the
existing `enso_grouping_nga`). Sources `base`, `config`. The sweep itself is the
registered batch `kl_sweep_classes_prodtu_all` (`make_batch_configs.R`), compared
by `05_fit_postprocess_kl_sweep.R`.

A **sixth standalone CLI tool**, `make_pipeline_docs.R`,
renders the machine-readable run graph `analysis/pipeline_steps.yml` into the
generated step table of `analysis/PIPELINE.md` and into
`notebooks/pipeline_run_graph.dot`, so the two cannot drift apart. `--check`
writes nothing and exits 1 if either is stale. Its functions
(`pipeline_steps_load()`, `pipeline_validate()`, `pipeline_render_*()`) are
sourced by `analysis/tests/test_pipeline_doc_coverage.R` behind the same
`sys.nframe() == 0L` CLI guard. Depends on `base`.

`dur_data.R` sits downstream of `siz_data.R` rather than beside it: it consumes
`outbreak_size_analysis_data.rds` instead of rebuilding the seasonality-match
cascade, which both avoids duplicating ~250 lines of `03_prepare_siz.R` and
makes the size and duration corpora identical by construction — a prerequisite
for the size–duration dependence analysis, which pairs the two fits' residuals
per episode.

It does **not** use `build_coverage_mask()` (`sim_io.R:504`) for its censoring
indicator, despite the surface similarity. That helper bridges
`location_period_id → gid` and filters to ADM2/ADM3 because the simulator works
in gid space; the censoring question is asked of the observation record, which
is in LP space. Measured on the current corpus, the bridge left 31 of 1582
episodes with no coverage window and produced 30 spurious censoring flags for
NGA (country×ADM2 fallback window ending 2020-01 against episodes running to
2022-08). Keying directly on `location_period_id` matches all 1582. See the
`derive_series_end()` roxygen for the full comparison.

`anchor.R` is a **sibling of, not a replacement for**, `provenance.R`'s
`build_mean_rate_table()` / `resolve_mean_rate()` / `persist_mean_rate_table()`, which
stay side by side. Its table's schema is deliberately
identical to theirs — same `cell` / `group` / `global` / `meta` shape, same
T1→T5 tier labels, same `mean_rate` column — which buys two things: `resolve_base_rate()`
delegates to `resolve_mean_rate()` rather than carrying a second copy of a
five-tier fallback ladder, and the size track can be built through *either*
function and compared cell-for-cell. It is: under `estimator = "geometric"`,
`build_base_rate_table(opt, "siz")` reproduces `build_mean_rate_table(opt)` at exactly
zero relative difference on `mean_rate` (and identically on group/global/`k`/
`level_factor`), with the sidecar filename unchanged. The only residual
difference is a `names` attribute the original carries on its *diagnostic*
`mr_parent` / `level_factor` columns. So `read_size_data()` (`siz_data.R`) builds
through the shared builder (the original hardcodes the geometric estimator, and
size's default is `pooled`).

The **estimator**: `geometric` (size only) estimates each cell by a geometric
mean of per-row `log(cases/pop)`; the other three tracks cannot, since occurrence and duration responses are 0/1
and monthly has structural zero months, so per-row `log(rate)` is `-Inf` on most
rows. They need a pooled `log(sum(num)/sum(den))` instead, with the EB shrinkage
weight rederived for it: the Poisson delta-method variance of `log(Y/E)` is
`1/Y`, so `tau2/(tau2 + 1/Y)` reduces to the *existing* `n/(n + k)` expression
under `n := Y_c` and `k := 1/tau2_between` — the shrinkage code is shared
verbatim, only those two inputs change. Zero-event cells fall back to the parent
exactly (tier `"T1-zero"`, a label distinct from `T2` so diagnostics can tell
"observed, saw nothing" from "never observed").

**All four tracks default to `pooled`**, so there is one estimator. For size,
pooled vs geometric moves the frame aggregate `sum(mu)` by +47.8% (median per-cell ratio 1.245, 24.4% of cells beyond
2× either way, mostly in the Jensen direction), but **zero cells and zero panel
rows change fallback tier** — values move, the information source does not — and
the aggregate is absorbed by a −0.391 intercept shift, comfortably inside the
anchored `normal(0, 1)`. `geometric` stays fully supported per-config
(`size.mean_rate.estimator: geometric`) and reproduces `build_mean_rate_table()`
cell-for-cell at exactly zero relative difference. The
sidecar filename carries an `est-<estimator>` token, so the two tables can never
collide on disk.

`footprint.R` exists because an episode in the size and duration frames is a
**`(gid, year_seas, admin_level)` footprint**, not a location period. The
collapse at the end of `03_prepare_siz.R` sums `suspected_cases` across every
LP resolving to the same gid in the same season-year and takes all other
columns — including `uid`, `location_period_id`, `TL`, `TR` — with
`dplyr::first()`. So those id columns are a *representative label*, not a key.
Measured on the current corpus, 314 of 1581 episodes (19.9%) absorb cases from
more than one LP, covering 29.7% of total case mass (238,350 cases; worst ratio
127×, concentrated in COD/ETH/CMR/NGA). Joining a monthly case series onto the
panel by `(uid, location_period_id, month)` therefore passes every row-count
and duplicate check while silently returning only the label LP's series. The
module's answer is `footprint_member_rows()`, which captures the full membership
before the collapse and is persisted by `03_prepare_siz.R` to
`make_footprint_members_file()`.

**Ordering matters**: `footprint_member_rows()` must be called *before*
`unify_footprint_ids()`. Under `--pop_source legacy`, unification overwrites
`uid`/`location_period_id` with a single representative across all season-years,
so a map built afterwards would attribute every season's cases to one LP and
reintroduce the double-count. The `footprint_member_rows()` roxygen states this and
`aggregate_monthly_to_footprint()` returns its conservation check rather than
enforcing it, so the caller decides fatality (`03_prepare_mon.R` gate M3 does).

## 6. Partition registry

The simulation partition (TU) is decoupled from the ENSO fit grouping
(`simulator.filter.partition_by: registry` + `partition_source`). Helpers, by
module:

| Module | Added | Notes |
|---|---|---|
| `config.R` | `registry_path()`, `registry_digest()` (xxhash64 of a registry CSV), `load_grouping_meta()` (registry + optional sidecar `analysis/groupings/<source>.yml`: label_prefix, order, label/short/short_fig/colour), `resolve_group_label()`, `grouping_label_token()`, `grouping_source_from_token()`, `assert_partition_labels()` | `registry_digest()` delegates to `base.R::content_digest()` on the file's raw bytes (call-time `base`; config.R sources nothing). Ownership pinned by `test_config_split_symbols.R`. |
| `spatpomp_core.R` | `load_partition_map()`, `stamp_partition_column()`, `check_partition_vs_grouping()` (pre-flight: fit-group coverage, mixed-TU duration basehaz, contagion blocks must refine the TUs), cache-dir identity `smoothing_cache_meta()` / `check_smoothing_cache_meta()` / `smoothing_cache_meta_path()`, exact resolver `.smoothing_cache_hits()`, `assemble_partition_smoothing_records` (alias; refuses mismatched replicate sets) | Call-time: `config` (`load_grouping()`, `registry_digest()`, `assert_partition_labels()`). |
| `smoothing_core.R` | registry branch of `smoothing_apply_partition()`, `smoothing_preflight_partition()`, `.smoothing_nbr_partition_source()` | |
| `paths.R` | `nbr_covariate_name()` / `parse_nbr_covariate()` / `nbr_artefact_token()` take `partition = "registry"` (`_part<token>`), `.nbr_partition_token()` | `_pg` and plain names byte-identical. |
| `spatial_contagion.R` | `build_contagion_table(partitions =)`; `resolve_nbr_spec()` registry terms; `nbr_cfg_from_sim()` defaults the contagion registry to `filter.partition_source` | Call-time: `config` (`grouping_source_from_token()`). |
| `sim_io.R` | `gid_iso3()`, `censor_sim_by_country()`, `manifest_units()`, `manifest_runs()`, `manifest_country_horizons()` (a country's own `countries.<iso>` horizon first), `manifest_has_country_horizons()`, `manifest_country_groups()`; `combine_group_runs(horizon_by_country=)`; `for_each_hindcast_run()` iterates v2 units | |
| `enso.R` | `group_colour_map(meta)` (`enso_group_colours()` reads the enso_grouping sidecar); `fit_oni_ar1()`, `continue_enso_ar1()`; `sample_enso_cfsv2_trajectories(beyond_cfsv2 =, ar1_phi =, ar1_fit_window =, legacy_seam =)` | Call-time: `config` (`load_grouping_meta()`); `base` (`enso_lag_col_name()`: `join_enso_lags()`, `add_enso_leads_from_lag0()` and the four samplers' `enso_lag<k>` / `enso_lead<k>` names). |
| `results_core.R` | `results_group_meta()`, `results_run_configs()`, `.manifest_run_entry()`, `.RESULTS_STAMP_BACKCOMPAT`; `load_sim_run_monthly(horizon_by_country =, keep_beyond =)`; `GROUP_*` tables read the enso_grouping sidecar | Sources `sim_io.R`. |
| `results_derive.R` | `.load_unit_arm_runs()` (v2 manifests) | |
| `deploy_core.R` | `deploy_partition_source()`, `deploy_nbr_partition_source()` (refinement check), `.deploy_group_meta()` / `.deploy_group_full()` | |

## 7. Shared infrastructure

Shared helpers used by sections 8-10, with no output change of their own. Gate: `tests/test_t2_shared_infra.R` (HARD, quick).

| Module | Added | Notes |
|---|---|---|
| `base.R` | `content_digest(x, digits = 15)` (xxhash64 of a canonical text serialisation: data.frame -> fixed-format CSV text in the given column / row order, numbers `%.15g`, strings quoted, NA bare; lists as length-prefixed blocks; a raw vector hashed as its bytes), internals `.content_canonical_text()` / `.content_canonical_values()` / `.content_lines()` | `digest` via `requireNamespace()` (NA when absent). `config.R::registry_digest()` delegates to it on the file's raw bytes (same value, pinned). Users: the Lenssen cache digest, the WHO-scale table digest, the contagion coverage digests. |
| `spatpomp_core.R` | `smoothing_cache_identity(S, tus, optional)` (the meta.yml list 08c writes, built in one place), `.SM_META_OPTIONAL_FIELDS` (field -> legacy value), `.cache_meta_add_optional()`, `.cache_meta_optional_check(meta, field, expected, legacy)`, `.cache_meta_same()` | Rule for every optional field: written only when not the legacy value, an absent field reads as the legacy value, a mismatch is an error naming the field and both values. `check_smoothing_cache_meta()` checks every registered field. 08c's meta.yml is byte-identical. Fields: `pop_year`, `who_scale_digest` / `db_overlap`. |
| `provenance.R` | `manifest_add_provenance(manifest, name, value, default)` | Adds `manifest[[name]] <- value` only when `!identical(value, default)` (a NULL value adds nothing), so a default run's manifest regenerates byte-identically. In provenance.R because all four manifest generators (`make_forecast_configs.R`, `make_hindcast_configs.R`, `make_event_attribution_configs.R`, `make_forecast_attrib_configs.R`) already source it. Users: `enso_counterfactual`, `pop_year`, `firm_oni_end`. |
| `base.R` | `.TTLO_ORIGIN` (2000-01-01), `.TTLO_LEVELS` (the seven `time_to_last_outbreak` levels, `"[0,1)"` .. `"[6,Inf]"`), `.TTLO_NONE_DAYS` (1e9, the sentinel of the level `none`, §8) | Used by (identical by test against frozen copies): `occ_data.R::make_time_to_last_outbreak_cat()` (`cut(labels =)` + the NA level), `03_prepare_occ.R`'s TTLO origin, `simulator_core.R::bin_ttlo()` / `ttlo_init_from_panel()` / `evaluate_lp()`'s reference bin. `simulator_core.R` thereby reads `base` at call time (as for `enso_lag_col_name()`), so its callers source `base.R` first. Left as literals: the Csnippet bin chains in `spatpomp_core_c.R` (`.ttlo_bin_c()`), `spatpomp_core.R`'s NA -> `"[6,Inf]"`, `siz_data.R`'s size-track copy, `occ_post.R`'s regexes and the size-track origin in `03_prepare_siz.R` (overwritten there). |

## 8. Simulator: population year, WHO scale, TTLO level

Population year of the simulator and the smoothing. Key
`simulator.location_frame.pop_year` (default `2020L` = a fixed 2020 population; `season_year` / `calendar_year` = a per-month population). Gate:
`tests/test_sim_pop_year.R` (HARD, quick).

| Module | Added | Notes |
|---|---|---|
| `simulator_core.R` | `.SIM_POP_YEAR_DEFAULT`, `.SIM_POP_YEAR_MODES`, `sim_pop_year_spec(sim)` (the one reader of the key; refuses `simulator.population.base_year` != 2020), `sim_pop_year_is_mode()`, `sim_frame_pop_args(sim, opt_paths, gid_pop_path, wpp_path)` (the frame-population block shared by 07a, 07b_v1, 07c, `simulations.R` and `smoothing_core.R`), `season_year_of(month, min_seas)`, `frame_pop_has_years()`, `frame_pop_at(frame, month)`, `frame_pop_matrix(frame, months)`, `sim_pop_window(sim)` (`.sim_panel_window()` widened to the filter windows); `build_location_frame(pop_window =)` (mode: `pop_y<YYYY>` + `pop_shift_m` columns); `monthly_log_mu(pop =)`; `simulate_one_tu()`'s per-month population / cap matrix | An integer pop_year adds no column and every kernel keeps its static path, so the default is bit-identical. |
| `spatpomp_core_c.R` | per-month population in `precompute_mon_lp_static()` (a temporary `.pop_t` covariate column set by `build_spatpomp_tu_c()`), `build_covar_for_tu_c(mon_mu_cap_by_row =)` | Absent by default: covariate table unchanged. |
| `spatpomp_core.R` | `.SM_META_OPTIONAL_FIELDS$pop_year` (legacy 2020L), `smoothing_cache_optional_fields(S)` | 08c writes `pop_year` to the cache `meta.yml` only when not 2020; 07b's forecast expects its own. |
| `smoothing_core.R` | `smoothing_frame_pop_args(S)`; `smoothing_make_config()` threads `location_frame.pop_year` / `gid_pop_path` / `wpp_path` as 07 does | `who_f`'s population year under a mode: year of the WHO window start. |
| `results_derive.R` | `rate_pop_weights(months)`, `.sim_pop_year_fields()`; `build_sim_scope_pop(pop_year_weights =, wpp =)`, `build_adm1_obs_pop(wpp =)`, `build_forecast_adm1_pop()`'s `pop_rate` column, `attrib_adm1_universe()`'s `pop_rate0` | Stamped `results.derivation.rate_pop_year` (`ref` default and back-compat; `window` = month-weighted mean population over the months read) moves only the per-capita rates (`inc_ann_100k`, `an_per_100k`). |
| generators | `make_forecast_configs.R`, `make_hindcast_configs.R`, `make_event_attribution_configs.R` record `pop_year` in the manifest via `manifest_add_provenance()` (only when not 2020) and source `simulator_core.R` for `sim_pop_year_spec()` | The attribution manifest inherits it from the forecast manifest. |

Per-country WHO term scale. Keys `simulator.filter.who_obs.{scale_mode,
scale_est.*, db_overlap}` (defaults `fixed` / `score_both` = a fixed scale).
Gates: `tests/test_who_scale_estimate.R` (HARD, quick) and `tests/test_spatpomp_who_obs.R`
W7-W9 (the fixed-scale Csnippet text / logLik equal `tests/golden/who_obs_csnippet_golden.rds`).

| Module | Added | Notes |
|---|---|---|
| `who_dashboard.R` | `who_db_observed_months()`, `estimate_who_scale()`, `who_multiplier()`, `who_drop_db_overlap()`, `digest_who_scale()` | Pure; `digest_who_scale()` calls `base.R::content_digest()` at call time. |
| `spatpomp_core.R` | `.resolve_who_obs_cfg()` resolves these keys (read with `[[ ]]`); `.SM_META_OPTIONAL_FIELDS` includes `who_scale_digest` (legacy `""`) / `db_overlap` (legacy `"score_both"`), record-only (`.SM_META_RECORD_ONLY`: checked only when the reader expects them) | |
| `spatpomp_core_c.R` | `build_spatpomp_tu_c()` takes `who_obs[["scale_u"]]` (per-unit estimated scale; validated) | NULL: the fixed-scale expression. |
| `smoothing_core.R` | `smoothing_who_scale_paths()`, `smoothing_write_who_scale()`; `smoothing_load_artefacts()` returns `who_scale` / `who_skipped`; `who_vs_smoothing(who_scale =)` | 08a writes `<run>/who_scale.csv` + `who_scale_meta.yml` and the digest into each task, 08b refuses a mismatch, 08c records the digest in the cache meta.yml. |
| `results_derive.R` | `derive_report_who_scale()` / `ensure_report_who_scale()`, `.who_scale_files()`, `.who_scale_frame_meta()`, `.who_db_consistency_path()` | Kind `report_who_scale` v1, table D30. 09b sources `who_dashboard.R`. |

TTLO level "none" in the simulator (contract with the occurrence panel: label "none" after "[6,Inf]", "[0,1)" the reference, brms column
`time_to_last_outbreaknone`, sentinel `base.R::.TTLO_NONE_DAYS`). Gate:
`tests/test_ttlo_none_level.R` (HARD, quick).

| Module | Added | Notes |
|---|---|---|
| `simulator_core.R` | `.TTLO_NONE_LEVEL`, `.TTLO_NONE_COL`, `bin_ttlo(has_none =)`, `.draws_has_ttlo_none()`, `frame_first_onset()`, `ttlo_no_prior_onset_at()`, `.ttlo_no_prior_for_frame()`, `build_warmup_state(no_prior_onset =)`; `ttlo_init_from_panel()` / `ttlo_col_to_level()` know "none"; `evaluate_lp()` bins with the level only when the bundle has it | Fits without the level: unchanged bins (the sentinel never exists). |
| `spatpomp_core_c.R` | `.ttlo_bin_c()` (replaces the three C bin chains), `precompute_ttlo_coefs()`'s `ttlo_coef_b7`, `generate_rprocess_csnippet(ttlo_none =)` | Without the level the C text is string-identical (W8 golden of test_spatpomp_who_obs.R covers the three chains). |
| `spatpomp_core.R` | the warm-start row hands a sentinel state off as "none" (`bin_ttlo(has_none = TRUE)`) | |
| loaders | `simulations.R` (both), `07b_simulator_v1.R`, `smoothing_core.R` put `first_onset` on the frame when the occurrence draws carry the level | |

## ENSO phase, horizon convention, counterfactual analogs

Gates: `tests/test_enso_phase.R`,
`tests/test_horizon_convention.R` (both HARD, quick).

| Module | Added | Notes |
|---|---|---|
| `enso.R` | `ENSO_PHASE_THRESHOLDS`, `ENSO_EPISODE_MIN_RUN`, `ENSO_PHASE_CODES`, `ENSO_PHASE_LABELS`, `ENSO_PHASE_STATISTICS`; `classify_enso_phase(x, rule, thresholds, min_run)`, `enso_phase_definition(rule, statistic)` | The six +/-0.5 sites (`.classify_enso_episodes()`, a wrapper; `classify_enso_analog_sets()`'s labels; `enso_event_launch()`; 09 `oni`, `report_dataset`, `report_smoothing`) go through `classify_enso_phase()`. Only the episode rule is "NOAA CPC". |
| `enso.R` | `ENSO_SEAMS`, `enso_seam_month(last_cfsv2_month, seam)` | Firm ONI end (`make_forecast_configs.R`: firm horizon = this + min lag) and the CFSv2 sampler's analog seam. `last_cfsv2` = today. |
| `base.R` | `horizon_month_seq(start, end)` (closed), `filter_window_months(start, end)` (half-open) | The one window convention: forward horizons closed, filter windows half-open with end = handoff. |
| `simulator_core.R` | `build_climate_panel_for_tu(end_mode = "legacy" / "closed" / "open")` | Same slice in every mode; `closed` warns when the end month is missing (07b spatPomp and v1), `open` for the smoothing's filter window. |
| `sim_io.R` | `assert_arm_launch()` also checks `horizon.end == filter.forecast_window.end`; `assert_panel_covers_horizon(panel, end, label)` | Config-time panel coverage (hard error in `make_forecast_configs.R`; the event generator's own check is closed). Reuses `manifest_last_firm_oni()`. |
| `enso.R` | `LENSSEN_MODELS`, `lenssen_extract_rda()`; `.lenssen_select_analogs_legacy()` / `.lenssen_extract_trajectories_legacy()` / `.lenssen_rescale_legacy()` | Both 03 builders (legacy cache, continuation library) share the model list and the Rda extraction. |
| `enso.R` | `enso_analog_neighbours(horizon =)`, `.enso_analog_distance()`, `enso_analog_paths()`, `enso_seam_state()` | `continue_enso_analog()` calls `enso_analog_paths()` (RNG-identical). `enso_continuation_skill.R::enso_continuation_library_rehorizon()` is untouched. |
| `enso.R` | `ENSO_COUNTERFACTUAL_DEFAULTS`, `enso_counterfactual_cfg()`, `build_lenssen_counterfactual_pool()`, `classify_enso_analog_sets_multi()`, `enso_read_winters()`, `enso_read_end()`, `.enso_shift_months()`, `lenssen_n_lead_months()`, `lenssen_build_command()`, `lenssen_cache_spec()`, `lenssen_cache_digest()` (`base.R::content_digest()`), `build_lenssen_cache()`, `lenssen_cache_status()`, `lenssen_cache_provenance()`, `lenssen_counterfactual_args()`, `LENSSEN_SHARED_CACHE`, `ensure_lenssen_counterfactual_cache()` | The library matcher (one matcher, two uses) and the per-run caches built in-process by `make_forecast_configs.R` / `make_event_attribution_configs.R`. Call-time: `base` (`content_digest()`). The event generator's `event_lenssen_lead_months()` / `event_lenssen_command()` are wrappers. |
| `enso.R` | `ENSO_COUNTERFACTUAL_SEAMS`, `.enso_pinned_paths()`, `enso_launch_oni()`, `.lenssen_cache_spec_base()`; `build_lenssen_counterfactual_pool(pin =)`, `lenssen_cache_spec(launch_oni =)` | Launch-month seam pin: `enso_counterfactual.seam_at` `launch` (library default) pins cache column 1 to the observed launch ONI and fades the offset from there; `launch_minus_1` = no pin (spec / cache / digest of the unpinned builder). `lenssen_cache_status()` compares the union of the spec fields. |
| `enso.R` | `sample_enso_cfsv2_trajectories(analog_seam =)`, `enso_continuation_provenance()` records `seam` off the default; `.rebase_offset()`, `cfsv2_correction_provenance()` | Runtime key `simulator.enso.continuation.seam`; the seam rebase (numeric-string rebase; projection-panel stamp fields only off the default). |
| `sim_io.R` | `assert_counterfactual_cache_identity(expected, runs, label)` | 09 (`results_context()`, `report_event_attrib`) checks manifest digest = cache digest = each run's `partition_meta.yml` digest; no-op for legacy manifests. |
| `results_derive.R` | `analog_set_tables()`, `analog_set_caches()`, `derive_report_analog_sets()` / `ensure_report_analog_sets()` (kind `report_analog_sets` v2, table D29, + `run_design`, table D34, see section 12); `smoothing_phase_rule()`, `smoothing_month_phases()`, `smoothing_phases_from_runs()`, `.oni_artefact_file()`, `.resolve_exclude_groups()` | The oni step reads each counterfactual config's `cache_path` (`ctx$LENSSEN_CACHE` stays a context variable) and lists the caches as inputs. |
| `make_event_attribution_configs.R` (tool) | `load_enso_events()` reads `peaks` / `peak_season_years` / `sub_windows`; `event_sub_windows()`, `.event_sub_window_frame()`, `event_sub_window_ids()` | `analysis/configs/enso_events_v2.yml` (one-event La Niña 2020–23); consecutive per-winter sub-windows per group and per country; event manifest schema 3 (`events.<id>.sub_windows[k]`) only when an event has them. |
| `sim_io.R` | `event_sub_windows_table(em, event_id)` | One row per sub-window x key (country bounds when the manifest has them, else group); NULL for a schema-2 event. Call-time: `dplyr`, `tibble`. |
| `results_derive.R` | `.cell_bounds_frame()` (`.cell_frame()` reads a cell's `bounds`), `event_sub_window_cells()`, `event_cells_all()`, `event_attrib_cells(cells =)`, `event_observed_adm1(country_first =)`, `event_sub_window_observed()`, `event_sub_window_keep()` | `report_event_attrib` v3, `event_af_posterior` v3, `report_event_af_model` v2 carry `sub_windows` (one cell per sub-window beside the whole event) only for schema-3 manifests; schema-2 payloads unchanged (`tests/test_event_sub_windows.R`, frozen copies). Call-time: `sim_io` (`event_sub_windows_table()`). |
| `results_core.R` | `.RESULTS_STAMP_NORMALISE` (named list keyed by dotted leaf path, `function(x, res)`), `.results_normalise_leaves()`, `.results_manifest_groups()` | `scope.exclude_groups` compared on the manifest's fit groups only; a failing normaliser compares strictly. |

## 8. Occurrence panel: masks, national zero-fill rule, contagion coverage

Occurrence-panel masks, the national zero-fill rule and the contagion coverage
stamp. Each behaviour sits behind a key or flag whose default keeps the base
panel, names and configs byte for byte.

| Module | Added | Notes |
|---|---|---|
| `paths.R` | the occurrence-panel spec `occ_panel_spec_of(opt)` / `.occ_panel_spec_from_episodes(ep)` / `.occ_panel_spec()` (fields `key`, `national_zero_fill`, `national_zero_threshold`, `national_zero_rule`, `at_risk_start`, `adm1_overlap`; defaults `.OCC_PANEL_SPEC_DEFAULTS`, enums `OCC_PANEL_AT_RISK_STARTS` / `OCC_PANEL_ADM1_OVERLAPS`, `.national_zero_rule_check()`), the one filename suffix `.occ_panel_variant_token()` (`_ep-gs` -> `_nz<thr>[-<rule>]` -> `_ae` -> `_a1c`; `.occ_panel_nz_token()` is the nz slot), the exact whole-suffix matcher `.occ_panel_variant_file_match()` (`.occ_panel_file_variant()`, `.OCC_PANEL_VARIANT_RE`) and `.episodes_block(spec, always)` (the `episodes` block every generated config / opt wrote by hand). `make_surv_dat_file()`, the occurrence base-rate path and `make_occ_run_file()` take the token; `read_deployed_lag(occ_panel =)`. | `.national_zero_fill_path_token()` and `.national_zero_fill_file_match()` stay as wrappers. Both match exactly (a `_nz60` lookup never matches `_nz60_ae_a1c` or `_nz60-r2` files). Gate: `tests/test_occ_panel_spec.R` (HARD). |
| `provenance.R` | `assert_occ_panel_variant()` (every spec field against the prep stamp; a stamp without a field reads as its default); `apply_episode_key_to_sim_config(occ_panel =)` (the panel fields under the same authority rule, written only when non-default) | `assert_national_zero_fill()` is a wrapper. |
| `national_zero.R` | the rule registry reader `national_zero_rule_params()` (`analysis/national_zero_rules.yml`, `NATIONAL_ZERO_RULES_PATH`; `legacy` = the 01 CLI defaults, `r2` = head 4, consistency scope all, ETH 2010-2016), `national_zero_rule_digest()`, `.national_zero_end_guard()` (no end inside the dashboard's partial year), `.year_edge_claim()` (one edge guard for tail and head), `build_national_zero_months(head_window =, consistency_scope =, edge_detail =)`, the D32 table `national_zero_audit()` | At `head_window = 0` and scope `threshold` the output, attributes included, is identical to a frozen copy (`tests/test_national_zero_months.R` [8]). 01 persists its table on the rule-free stem (`national_zero_months[_t<N>]_db_who_consistency.csv`, `paths.R::national_zero_consistency_path(threshold)`), the path `report_who_scale` reads. |
| `paths.R` | `make_national_zero_months_file(threshold, rule)` (`_t60_r2` / `_r2`), `national_zero_months_file_of(opt | spec)`, `national_zero_audit_path(months)` / `national_zero_consistency_path(threshold)` | `.occ_panel_spec()` carries `national_zero_rule`; token `_nz<thr>-<rule>`. |
| `occ_data.R` | `.own_episode_span_hits(surv, ep, rule = span | running)` (shared by the masks and the zero-fill rule), `national_zero_surv_rows(episodes =)` (a named rule keeps the fill out of own episode spans; attribute `n_in_episode`), `assert_nz_rows_outside_episodes()` (post-mask gate; asserted under a named rule, report-only under legacy) | `episodes = NULL` (legacy) gives the identical rows and attributes. |
| `deploy_core.R` | `deploy_check_national_zero_rule()` (rule in the registry; local months file built with its parameters); sources `national_zero.R` | Silent for the legacy rule. |
| `results_derive.R` / `results_core.R` | 09 kind `report_national_zero_audit` v1 (table D32; on only when the panel's months file has its audit) | 11 writes `D32_national_zero_audit.csv`. Gate `tests/test_results_occ_panel_tables.R` (HARD). |
| `gid.R` | `.overlap_adm1_set()` (the ADM1s with ADM2 rows; one definition behind `filter_overlapping_adm1()`, `drop_overlapping_adm1_from_frame()` and `adm1_overlap_parents()`, which classifies each ADM1 parent episode orphan / covered by month-span overlap with its ADM2 children's episodes) | Both filters identical to frozen copies (`tests/test_occ_panel_masks.R` [4]). |
| `occ_data.R` | `mask_occ_panel_rows(surv, ep, spec, overlap)` (at_risk_start after_episode via `.own_episode_span_hits("running")`; adm1_overlap censor_orphan; identical panel at the defaults; attributes n_* and `masked_rows`), `make_time_to_last_outbreak_cat(no_prior_level =)` + `ttlo_no_prior_onset()` (level `"none"` after `"[6,Inf]"`, load-time recode; `load_surv_data()` reads `occurrence.ttlo_no_prior_level`); `load_surv_data()` calls `assert_occ_panel_variant()` | Call-time `gid` (`gid_to_adm1()`). |
| `footprint.R` | `discarded_onset_table(ob_ts2, ob_ts3_first, ob_ts3_final, ep_res)` (reasons lp_second_onset / redate_duplicate / absorbed_member / absorbed_after_gap), `adm0_episode_counts()` | Read by 03 only. |
| `paths.R` / `occ_post.R` | `occ_panel_sidecar_paths()` (`<panel>.{adm1_overlap,discarded_onsets,adm0_dropped,masked_rows}.csv`), fit token `_ttn` (before the panel token) in `make_occ_run_file()`; `canonicalize_term()` maps `time_to_last_outbreaknone` -> `ttlo_none` | |
| `results_derive.R` / `results_core.R` | 09 kind `report_occ_panel_exclusions` v1 (table D31, per country; `occ_panel_exclusions_table()`, `OCC_PANEL_EXCLUSION_COLS`), on only when the sidecars exist | 11 writes `D31_occ_panel_exclusions.csv`. |
| `spatial_contagion.R` | `contagion_coverage_gid()` (coverage "reporting": `build_coverage_mask()` over the member bridge with ADM1 windows and the panel's national zeros, UNION the recorded episode months), `adm1_coverage_grid()` (at t-1), `adm1_activity_grid()` (global-month lag-1 activity), `contagion_coverage_shares()` (raw-kernel neighbour coverage share over the covariate's masked universe) | The covariate and `build_contagion_table()` are unchanged. |
| `occ_data.R` / `paths.R` | `active_last_month()`; `nbr_coverage_share_name()` (`nbr_cov_lag1_r<rho>[_pg|_part<tok>]`, never a model term), `contagion_coverage_csv_path()` | |
| `results_derive.R` / `results_core.R` | 09 kind `report_contagion_coverage` v1 (table D33; `contagion_coverage_table()`, `CONTAGION_COVERAGE_COLS`), on only when the panel's coverage sidecar exists | 11 writes `D33_contagion_coverage.csv`. Gate `tests/test_contagion_coverage_report.R` (HARD). |
| `provenance.R` | `assert_prep_stamp_overwrite(path, new_stamp, force)` -- 03_prepare_occ.R builds its whole stamp before writing anything and refuses to replace an existing panel whose stamp differs in any field but `stamped_at` (or that has none) unless `--force` | A columns-only change keeps the filename; this keeps an existing panel from being rewritten in place. Gate `tests/test_occ_panel_spec.R` [6]. |
| `coverage.R` | `build_coverage_mask()`, whose `national_zero_path` is REQUIRED while `national_zero` is on (the default threshold-0 file would silently mismatch a thresholded run); `coverage_lp_bridge(obs_panel, surv_path)`: the panel's full LP -> gid member bridge sidecar `<panel>.lp_members.csv` (`paths.R::occ_panel_lp_members_path()`, written by 03 under gid_season) when it exists, else the panel's representative LP per gid | Sourced next to `sim_io.R` by 07c, `smoothing_core.R`, 09a/09b and the `utils.R` facade. 07c, 08c and 09's smoothing step bridge through it; 09 lists the sidecar as an input only when it exists, so artefacts built without it stay fresh. Gate: `tests/test_national_zero_months.R` [3]. |
| `simulations.R`, `smoothing_core.R`, `deploy_core.R`, `hindcast_skill.R`, `results_derive.R`, `dlnm_kernel.R`, the three generators | the spec is threaded through `.discover_group_fit(occ_panel =)`, `S$OCC_PANEL`, `.deploy_foreign_episode_artefacts(occ_panel =)` + the per-lag / template checks, `hindcast_run_index(occ_panel =)` / `hindcast_panel_paths(occ_panel =)`, every `read_deployed_lag()` call | The `national_zero_*` pair arguments are kept (side by side). `config.R`'s results -> episodes copy keeps its own `modifyList()` (config.R sources nothing, so it cannot call `.episodes_block()`). |

## 9. Deploy defaults, deployed-run batches, convergence report

| Module | Added | Notes |
|---|---|---|
| `config.R` | `get_default_options()$deploy` has `convergence` (off, population, "default", "default", 0, []), `enso_counterfactual` (= `enso.R::ENSO_COUNTERFACTUAL_DEFAULTS`), `enso_continuation.seam` (last_cfsv2) and `episodes.{national_zero_rule, at_risk_start, adm1_overlap}` ("", season, drop) | The in-function defaults read with `[[ ]]`; restated in every deploy YAML. Pinned by `test_deploy_convergence_gate.R` [T5], `test_enso_counterfactual_match.R` [T5]; `test_deploy_config_coverage.R` refuses a string leaf shipped as a YAML boolean (a bare `off`). |
| `make_batch_configs.R` | `.drop_superseded_prior_sd()` (opt-in `spec$drop_superseded_prior_sd`), `.deployed_lag_rows()`, `.batch_deferred()`; `generate_batch()` honours the optional `spec$rows_fn` (row selection after the varying dims are fixed) and stops on a deferred batch; `generate_all_batches()` skips deferred batches; sources `paths.R` (`nbr_covariate_name()` for the kl sweep's coefficient-prior names) | Entries `prod_run_{occ,dur,mon}_grouped_enso_singlelag`, `prod_run_kl_sweep_classes`. Gates `tests/test_prod_run_batches.R`, `test_batch_config_block_stripping.R` [E9] (HARD). |
| `results_derive.R` / `results_core.R` | `cf_constrained_winter()` reads every classified winter of a library-mode cache (`set_definition$winters`: primary + the other read winters when held neutral), `.cf_winter_cell()` joins them in S28; a legacy cache keeps its one `peak_season_year` winter. `report_convergence` v2, `.convergence_gate_cfg()` (the results deploy's `deploy.convergence`, `deploy_convergence_cfg(validate_allow = FALSE)`; defaults without a deploy) and `.convergence_gate_rows()` (per fit `fit_diagnostics_gate()` on the summary / sampler tables the convergence artefact stores: `gate_verdict` pass / allowed / fail, `gate_reasons`, the scope's worst R-hat and parameter) | No fit is reloaded (the gate reads only what `convergence` v1 stores; `convergence` stays v1); the deploy YAML joins `report_convergence`'s inputs. Gates `tests/test_results_convergence_gate.R`, `test_results_africa_window.R` [3b] (HARD). |

## 10. Event constraint scope, prep overwrite guards

| Module | Added | Notes |
|---|---|---|
| `make_event_attribution_configs.R` (tool) | `load_enso_events()` reads an optional per-event `other_winters` (neutral / any); `event_counterfactual_cfg(cf, ev)` | Constraint scope: the event's registry `other_winters` overrides `deploy.enso_counterfactual.other_winters` for that event only (`enso_events_v2.yml`: `lanina_2020_23: other_winters: any`, so only DJF 2020/21-2022/23 constrain its sets; DJF 2023/24 is still classified and recorded). An event without the key is unaffected; the manifest's `lenssen_winters.other_winters` appears only for an event that sets it. Gates `tests/test_event_sub_windows.R` [1], `test_deployed_lag_horizon.R` [10]. |
| `results_derive.R` | `analog_set_tables()`: winter role `read_unconstrained` (cache `set_definition$winters$other_winters == "any"`); `constrained_winters` then lists the primary winters only | D29 tables of caches built under `neutral` are unaffected (`tests/test_results_analog_sets.R` [1]). 11 adds one sentence only when a `read_unconstrained` row exists. |
| `provenance.R` | `assert_prep_stamp_overwrite(ignore =, old_stamp =)` (defaults = the base guard), `prep_force_option(what)`, `who_monthly_guard_stamp(meta, data)` | The overwrite guard of `03_prepare_occ.R` also runs in `03_prepare_siz.R` (size frame + footprint map), `03_prepare_dur.R` (`ignore` adds `restamped_at`), `03_prepare_mon.R`, `03_prepare_sim_panel.R`, `03_prepare_sim_panel_projection.R`, `03_prepare_sim_panel_extended.R` and `03_prepare_who_monthly.R` (its `_meta.yml` less `written` / input mtimes + a content digest, `old_stamp`): each builds its stamp first, refuses an existing artefact whose stamp differs or that has none unless `--force`, then writes. An identical rebuild passes. `03_prepare_who_monthly.R` sources `provenance.R`; `who_monthly_guard_stamp()` calls `base.R::content_digest()` and `yaml` at call time. Gate `tests/test_prep_overwrite_guards.R` (HARD, quick). |

## 11. Data leg on a SLURM cluster, faster 00-02

Speed-ups of phases 00-02 that keep their outputs, the per-country population runner, and the
Slurm stages of the data leg (`bash/campaigns/preprocess.sh`). Equality: the per-country runner
against the per-year one in `tests/test_population_crosswalk.R` [12].

| Module / tool | Added / changed | Notes |
|---|---|---|
| `population_core.R` | `run_population_country_task()`; `run_population_task()` built from `.pop_country_geometry()`, `.pop_extract_roles()`, `.pop_assemble_task()`; `pop_part_fresh()`, `pop_gid_cache_key()`, `read_gid_pop_cache()` / `write_gid_pop_cache()`; `load_gadm_units()` reads with an SQL filter on `GID_0` | One country, all years: geometry once, one extraction per year (the per-year call is run_population_task's, so a year's tables are identical). Q_g (gadm1 / gadm2 roles) is cached per (country, year) keyed on the GADM file, raster / boundary dirs and the OutbreakExtractR version, not the outbreak vintage; with Q_g cached only the LP roles are extracted. A partition counts as done only when its manifest's `lp_md5` (or, for older manifests, path + mtime) matches the current LP corpus. `workers` (> 1) extracts the years at once on forked processes (`future::multicore` via `furrr::future_map`, no RNG), result identical to the serial loop. `qs2`, `yaml` via `::`. |
| `national_zero.R` | `.month_span()` | Vectorised month expansion (`rep` / `sequence`) in `national_month_status()`, `subnational_positive_months()`, `active_episode_months()`; (an empty sub-national input returns a Date `month` column). |
| `clim_db.R` | `get_stats_gid(path =)`, `get_extremes_gid(path =)` | `NULL` (default) = the default databases. |
| `clim_compile.R` | `read_clim_csv()` calls `readr::read_csv` | furrr workers do not attach readr. |
| `paths.R` | `make_pop_run_file(kind = "pop_gid_cache_dir")` | `generated_data/population/gid_pop_cache/`. |
| `preprocess_paths.R` (tool) | CLI | `KEY=value` for bash: `OUTBREAK_VINTAGE`, `VINTAGE_DIR`, `SECURE_ROOT`; `--deploy` adds `NZ_RULE` / `NZ_THRESHOLD` (deploy.episodes); `--outputs` (with `-C`) lists the 01-02 outputs as `SEC=` / `SECDIR=` / `GEN=` / `GENDIR=` to copy back from the cluster. Sources `config`, `base`, `paths`, `provenance`. |

## 12. Supplement tables and maps

| Module | Added / changed | Notes |
|---|---|---|
| `results_derive.R` | `RUN_DESIGN_ARMS`, `run_design_row()`, `run_design_configs()`, `run_design_table()`; `report_analog_sets` v2 (+ `run_design`, table D34) | One run config (and its analog cache) per simulation arm: replicates, max_draws (unset = 200, as the simulator), draw mode, handoff, analog pool, clamp margin. `.sim_param_draw_mode()` resolves the draw mode as 07b does (unset = `random`; blank = `spread_param_draws ? cycle : inherited`), shared with `af_posterior_check_design()`. Gate `tests/test_results_analog_sets.R`. |
| `results_derive.R` / `results_core.R` / `paths.R` | 09 kind `modelled_units` v1 (step `modelled_units`, after `onset_rate`): `modelled_units_tables()`, `.gadm_adm1_gids()`, `derive_modelled_units()` / `ensure_modelled_units()`; context `GADM_ADM1_PATH`; `make_results_run_file("results_modelled_units")` | Units of the survival panel by level (`gid_level()`) and ADM1 key (`gid.R::gid_to_adm1()`) with status adm1 / adm2 / no_onset, the ADM1s of the modelled countries without a recorded episode (GADM ADM_1 attributes), and the fit registry's `members`. 10's figures MU / GM / KC read it (Rule 10). Inputs: survival panel, registry CSV, GADM gpkg. |
| `results_derive.R` / `results_core.R` | `rain_enso` v2: `group_lags` (`rain_enso_lag_summary()`, table T1b's rainfall columns per weighting x group) and `mechanisms` (`rain_enso_mechanisms()`, a checked copy of `results.inputs.enso_rain_mechanisms`, table T1a); `report_lags` v3: `loo_summary` (`report_lags_loo_summary()`, table T2, at `rule$se_mult`) | 11 only displays them. |
| `config.R` / `config_defaults.yml` | `results.inputs.enso_rain_mechanisms` (default `analysis/tables/enso_rain_mechanisms.yml`; stamped, back-filled with the same path), `results.style.data_release` (unstamped, D / T table -> supp / github) | 11 falls back to `get_default_options()$results$style$data_release` when the style key is unset. |
| `results_core.R` | `GROUP_LABEL_MODES` includes `centroid`; `AFRO_SUBREGION`, `afro_subregion_of()` | Display only (style). |
| `kl_grouping.R` | `KL_CENTROID_NAMES`, `kl_centroid_fields()`; `kl_grouping_sidecar(centroid_from =)` | A regenerated fit-registry sidecar keeps `centroid_order` / `centroid_short` / `centroid_label`. Gate `tests/test_deploy_kl_grouping.R` [3b]. |
| `figures_supp.R` | `build_figS_modelled_units()`, `build_figS_group_map()`, `kl_classes_by_K_frame()`, `build_figS_kl_classes_by_K()` | Call-time packages: `sf`, `scales`, `grid`, `grDevices` beside ggplot2 / dplyr / tidyr / patchwork. |
| `results_report.R` | `with_afro()`, `report_data_id()`, `report_release_split()`, `report_run_design_table()`, `report_analysis_set_table()`, `report_t1a_mechanisms(mech, ...)`, `report_t1b_rain_lags()` (lookup only), `report_outside_range()` / `RANGE_MARK` | Formats only. Gate `tests/test_report11_supplement_tables.R`. |

## 13. Supplement figures

| Module | Added / changed | Notes |
|---|---|---|
| `figures_supp.R` | `data_coverage_frames()`, `build_figS_data_coverage()` | Fig. DC (`supplementary/figS_data_coverage`): episodes per country x cholera season year (ENSO-group blocks, display labels from 10's `group_display_meta()`), episodes per season year by admin level; reads 09's `report_dataset$outbreaks` (no schema change). |
| `figures_supp.R` | `processing_cascade_frames()`, `build_figS_processing_cascade()` | Fig. PC (`supplementary/figS_processing_cascade`): the processing cascade (Data S1 D1) as a flow diagram, spine + side boxes of what leaves it + fork into the three datasets; reads 09's `report_dataset$cascade` (no schema change). |
| `results_derive.R` / `results_core.R` | `report_event_attrib` v4: `.event_attrib_stamp()` renames the probability-ratio tables' outbreak event (`event`: onset_any / cases_ge) to `pr_event` before the past event's name is stamped, so the two stay distinguishable. 11's D14 table filters the headline `pr_event`. | Data S1 D14 has `pr_event`; the forecast tables (`report_attrib_metrics`) keep `event`. |
| `results_report.R` | `report_release_split()` copies through `.release_copy()`: columns empty in every row (e.g. bootstrap intervals with `attrib_boot_b: 0`) are left out of the released copies in `supp_data/` and `github/`; cells are copied as text; the numbered tables keep every column. | Gate `tests/test_report11_supplement_tables.R` [2]. |
| `results_derive.R` / `results_core.R` | `smoothing` v4: `derive_smoothing()` adds `validation` (`smoothing_validation_tables()`, pure, on the event tables + observed panel + location frame); `ensure_smoothing()` inputs include the aggregated series of the coverage mask and the HMM cluster file (when present) | Fig. SV (`supplementary/figS_smoothing_validation`). Randomised PIT of the monthly onsets (at-risk unit-months) and cases (covered unit-months, empirical replicate distribution), stored as 100-bin counts per scope x metric; onsets per at-risk unit-month by calendar month; episode-duration shares (runs of consecutive active months); fraction of covered unit-months active; onset frequency given an ongoing outbreak at m - 1 in the same ADM1 / another ADM1 of the country / another country of the TU (the kernel's `update_lags()` levels), among units not active themselves at m - 1. Every replicate, uniform weights. |
| `figures_supp.R` | `pit_hist_rebin()`, `build_figS_smoothing_validation()` | Fig. SV's builder; 10 relabels the scopes (`relabel_groups()`), display bins and duration cap are 10 constants. |
| `figures_supp.R` | `forecast_analogue_frame()`, `build_figS_forecast_analogues()` | Fig. FA (`supplementary/figS_forecast_analogues`): monthly onsets and cases, attribution factual (`attrib_ribbon` / `attrib_africa_ribbon`, else the forecast rows) vs the counterfactual of each analogue set (`region_ribbon` / `africa_ribbon`), Africa + each group, from the attribution start (to the common Africa window's end under a `_common` headline cell); reads the existing `attribution` artefact (no schema change). |
| `figures_supp.R` | `attribution_analogue_map_sf()`, `build_figS_attribution_analogues()` | Fig. AA (`supplementary/figS_attribution_analogues`): per metric, ADM1 maps of the ratio of means (Figure 4b's log scale and cap) and of `p_sign` (probability the factual exceeds the counterfactual over the paired replicates), one row per analogue set, at the headlined cell and explosive mode; reads the existing `report_attrib_exposure$adm1` (no schema change). Fill scales and base layers (a fresh grey layer per map, ggnewscale) are passed in by 10. |
| `figures_supp.R` | `hindcast_clim_skill_frame()`, `build_figS_hindcast_clim_skill()` | Hindcast skill against climatology (`hindcast/figH_supp_skill_vs_climatology`, drawn by 10's figure H with `--supplementary true`; gate `tests/test_supp_hindcast_clim_skill.R`): CRPSS with season-year bootstrap intervals by lead class, Africa + each group, for onsets per unit and season window (`hindcast_occ$skill`, window grain, `crps_fair`, `scope.hindcast.headline_occ_reference`) and monthly country cases, raw and log(1 + y) (`hindcast_cases$month_skill`, headline support and `headline_case_skill_reference`); no schema change. Display floor `HS_YLIM_LO` is a 10 constant. |
| `figures_supp.R` | `membership_heatmap()`, `membership_paths()`, `build_figS_enso_groups()`, `build_figS_membership_paths()` | Fig. EG (`supplementary/figS_enso_groups`, `figS_enso_groups_paths`; 09 `modelled_units`, `report_membership`, `report_kl_sweep`): the group map, membership probability per class, most probable class at K = 2..4 and, with `paths = TRUE`, the class paths under the map; `build_figS_membership()` is built from the two panel helpers. |
| `figures_supp.R` | `who_enso_ratio_frame()`, `build_figS_who_enso_ratio()` | Fig. WR (`supplementary/figS_who_enso_ratio`; 09 `report_who_enso_ratio`): per country of a group and the group, the observed WHO ratio of 2026 to the ENSO-neutral years against the smoothing / ENSO-counterfactual ratio of means, one panel per observed estimator; labels with `ggrepel::geom_text_repel()` (ggrepel is in renv.lock). Gate `tests/test_results_who_enso_ratio.R`. |

## 14. Figure 1 onset-season artefact

| Module | Added / changed | Notes |
|---|---|---|
| `results_derive.R` / `results_core.R` / `paths.R` | 09 kind `occ_seasonality` v1 (step `occ_seasonality`, appended, no `after`): `.OCC_SEAS_MONTH_STARTS`, `occ_gid_min_seas()`, `occ_seas_smooth_label()`, `occ_seasonality_fit_files()`, `derive_occ_seasonality()` / `ensure_occ_seasonality()`; `make_results_run_file("results_occ_seasonality")` | Figure 1c's onset-season column: per group, the deployed occurrence fit's cyclic `t_start_seas` spline at the 12 season-month starts, `exp(s - mean(s))` per draw, posterior mean, rotated to calendar months per gid by its `min_seas` (survival panel), mean over the group's gids. Loads the occurrence fits; inputs = the fit files, the survival panel, the grouping. |

## 15. Manuscript figure builders

| Module | Added / changed | Notes |
|---|---|---|
| `figures_manuscript_fig1.R`, `figures_manuscript_fig2.R`, `figures_manuscript_fig34.R` | the figure builders of `10_manuscript_figures.R`'s Sections 9a-9d, each with its roxygen block | 10 keeps the CLI, style, every load-time constant (`FIG1_*`, `FIG2_*`, `F3_*`, `FIG4_*`, the `HL_EXPLOSIVE` check), the shared helpers of Section 2c and the `fig_step_*()` dispatch, and sources the three modules in its module loop. Definitions only (they load alone, `tests/test_utils_module_purity.R`), but the functions read 10's globals at call time, so no other script may source them. `tests/test_enso_phase.R` [3] scans the modules for ±0.5 literals; `tests/test_hardcoded_default_parity.R` classifies `figures_manuscript_fig34.R`'s `f23` sites. |

## 16. Figure 1c seasonality design

| Module | Added / changed | Notes |
|---|---|---|
| `results_derive.R` / `results_core.R` / `paths.R` | `occ_seasonality` v2 (+ `draws`, `occ_seas_calendar_draws()`); `seasonality` v2 (+ the occurrence band); 09 kind `oni_composite` v1 (step `oni_composite`, appended, after `oni`): `oni_composite_cycle()`, `derive_oni_composite()` / `ensure_oni_composite()`; `make_results_run_file("results_oni_composite")` | occ_seasonality keeps each group's per-draw calendar curves (the gid mean per draw), so seasonality gives occurrence a band at its `interval`, as cessation and cases have. oni_composite: the observed ONI over Y, Y+1 averaged over the oni artefact's El Niño / La Niña episodes, aligned on their peak; its input is the oni artefact. |
| `figures_manuscript_fig1.R` | `fig1c_*()` builders and `build_fig1()`; the earlier layout's builder `build_fig1_v91()` (kept side by side, with its panel helpers) | Panel c: composite ONI, rainfall tiles, and per track (occurrence, cessation, monthly incidence) the ELPD lag weights over the seasonality relative to the annual mean; d as the last row of c's grid. Reads 09's oni_composite, rain_enso, seasonality v2 and lag_rank. Constants `FIG1C_*` in 10's Section 9a. |
