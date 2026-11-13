# Pipeline run graph

This file is the single record of the order in which the `analysis/` pipeline
runs. The repository has no Makefile and no `targets` pipeline, so the order is
written down here, and tests check it mechanically.
`analysis/bash/tools/run_pipeline.sh` chains the phases (see `analysis/bash/README.md`).

- The **step table** below is generated from `analysis/pipeline_steps.yml`.
  That file is the machine-readable run graph: one entry per numbered script,
  report and non-step tool stage, with dependencies, drivers and outputs.
- `notebooks/pipeline_run_graph.dot` is generated from the same file, so the
  table and the diagram cannot drift apart.
- Phase meanings and file naming are in the top-level `README.md`. The bash tiers are described in `analysis/bash/README.md`.

## Phase numbers are identifiers, not run order

Numbers name phases and are **not** renumbered to follow run order. Where the two differ:

- **08 runs before 07b on the spatPomp path.** `08c_smoothing_compile.R`
  writes the per-TU handoff cache (`cache/hist_filter_tu=<TU>_*.qs`) that the
  forecast and hindcast configs point to (`simulator.filter.smoothing_input`).
  `07b_simulator_spatpomp.R` cannot start a forecast without it. Only the v1
  simulator (`07a` → `07b_simulator_v1.R`) is independent of 08.
- **The deploy tool runs between phases** (`utils/deploy_prod_run.R`, which
  has no phase number). `gen_batches` runs before 04. `deploy_fits` and
  `smoothing` (together, `--stage all`) run between 05 and 08.
  `gen_forecast` and `gen_hindcast` run between 08 and 07b.
- **Phase 02 has an internal order across its unlettered files:**
  `02a` → `02b` → `02c` → `02_preprocess_crosswalk.R` → `02_preprocess_climate.R`.
  The climate step reads the crosswalk gids only under `--gid_source crosswalk`
  (its default). Crosswalk and climate also read 01's `outbreak_shapefiles.gpkg`,
  so they run after 01; the population steps read 00's LP shapefiles and run
  beside 01.
- **The simulator inputs from 03 belong to every rollout**, although no deploy
  stage names them: `03_prepare_sim_panel{,_projection}.R` and
  `03_prepare_lenssen_analogs.R`.
- **11 is a report, not a script.** `11_manuscript_numbers.qmd` is rendered
  after 09 and reads only phase-09 artefacts.
- **06 is not a gate.** Only `06_fit_figures_lag_compare.R` is part of the
  rollout's postfit campaign. The other 06 scripts are figures you run when
  you need them.

## A prod rollout, in order

Each block runs once the previous HPC jobs have finished. The run itself is
named only in the deploy YAML (`--deploy`; see `analysis/bash/README.md`
"Run identity").

```bash
# Data leg (when outbreak / climate / population inputs change)
# On a SLURM cluster, one campaign chains 00 -> 01 (+01b, 01_nz) -> 02 population -> 02 crosswalk -> 02 climate:
#   bash analysis/bash/campaigns/preprocess.sh --config analysis/configs/pop_extract/pop_extract.yml --deploy $D --plan
# Locally (the vintage folder is paths.R::defaultOutbreakDataDir(); 00 writes into it):
Rscript analysis/00_ingest_outbreakextractr.R
Rscript analysis/01_preprocess_outbreaks.R -r TRUE
Rscript analysis/01b_preprocess_outbreaks_monthly.R
Rscript analysis/01_preprocess_national_zero_months.R --rule <deploy.episodes.national_zero_rule> --zero_threshold <...threshold>
bash analysis/bash/stages/02_population.sh --config analysis/configs/pop_extract/pop_extract.yml   # 02a -> 02b per-country array -> 02c (afterok)
Rscript analysis/02_preprocess_crosswalk.R -C analysis/configs/pop_extract/pop_extract.yml
Rscript analysis/02_preprocess_climate.R -s TRUE   # defaults = the deployed build (parquet, crosswalk gids, dynamic)
Rscript analysis/03_prepare_occ.R -b 1990_2019
Rscript analysis/03_prepare_siz.R -b 1990_2019
Rscript analysis/03_prepare_dur.R -b 1990_2019
Rscript analysis/03_prepare_mon.R -b 1990_2019
Rscript analysis/03_prepare_sim_panel.R -b 1990_2019
Rscript analysis/03_prepare_sim_panel_projection.R -b 1990_2019 --cfsv2_correction calibration --enso_max_lead 11   # the deployed panel: CPC-style calibration (cfsv2_rebase is ignored under it) + ENSO leads
Rscript analysis/03_prepare_lenssen_analogs.R

# Model + simulator leg
D=analysis/configs/deploy/prod_run<N>.yml
Rscript analysis/utils/deploy_prod_run.R -C $D --stage gen_batches              # fit-batch config dirs
bash analysis/bash/campaigns/enso_grouped_singlelag.sh --deploy $D              # 04 (+05) occ; dur_/mon_ siblings
bash analysis/bash/campaigns/postfit.sh --deploy $D --tracks occ,dur,mon        # per track: 05 cb + loo_lags, 06 lag_compare
Rscript analysis/utils/deploy_prod_run.R -C $D --stage all                      # deploy_fits + gen_templates + smoothing check
bash analysis/bash/stages/08_smoothing.sh --config <printed historical template>  # 08a -> 08b -> 08c
Rscript analysis/utils/deploy_prod_run.R -C $D --stage gen_forecast
Rscript analysis/utils/deploy_prod_run.R -C $D --stage gen_attrib                       # attribution factual (deploy.outputs.forecast_attrib_manifest; no-op when unset)
Rscript analysis/utils/deploy_prod_run.R -C $D --stage gen_hindcast --align-results
Rscript analysis/utils/deploy_prod_run.R -C $D --stage align_results --dry-run   # or without --dry-run: point results/manuscript.yml at this run
bash analysis/bash/stages/07_forecast_groups.sh --deploy $D --arm all           # 07b spatPomp (+07c/07d): factual, attribution factual, counterfactual; complete runs skipped
bash analysis/bash/stages/07_hindcast.sh --deploy $D
# by hand: 07c_simulator_postprocess.R --combine-groups (07c across groups / arms), 07d_simulator_figures_counterfactual.R

# Manuscript leg
Rscript analysis/09b_postprocess_results.R -C analysis/configs/results/manuscript.yml   # local, serial; or on the HPC:
bash analysis/bash/stages/09_postprocess_results.sh --config analysis/configs/results/manuscript.yml   # 09a tasks -> 09b arrays, one per wave
Rscript analysis/10_manuscript_figures.R  -C analysis/configs/results/manuscript.yml
quarto render analysis/11_manuscript_numbers.qmd
bash analysis/bash/stages/10_manuscript_figures.sh --config analysis/configs/results/manuscript.yml --after <RESULTS_JOB_IDS>   # 10 + 11 on the HPC, one job
quarto render notebooks/hindcast_skill.qmd   # hindcast skill (09 artefacts; not a numbered step)
```

## Keeping this current

1. Add, rename or re-wire a numbered script, a tier-1 stage driver or a deploy
   stage → edit `analysis/pipeline_steps.yml`. List order is run order:
   anything a step needs must appear above it.
2. Regenerate the table below and the `.dot`:
   `Rscript analysis/utils/make_pipeline_docs.R` (`--check` only reports).
3. `analysis/tests/test_pipeline_doc_coverage.R` (HARD gate) fails when:
   - a numbered script, the phase-11 report or a `bash/stages/` driver is
     missing from `pipeline_steps.yml` or from this file;
   - an entry names a file or driver that does not exist;
   - the list order is not a valid run order;
   - this file's generated block or the `.dot` is stale.
4. `analysis/tests/test_script_header_names.R` (HARD gate) is the header
   convention:
   - the first comment line of every numbered script names its own filename
     (`# <file>.R — ...`), after any shebang or `# ====` rule;
   - every `Rscript analysis/NN*.R` command written in a numbered script, a
     `utils/` module or a maintained bash driver resolves to an existing file.
   This catches the `05b_loo_compare_lags.R` class of drift, where a renamed
   script kept describing itself under its old name.

## Step table

`Rollout`: `core` = run in every `prod_run<N>` rollout; `optional` =
benchmarks, sensitivity checks, extra figures, the v1 simulator. `After` =
the steps whose artefacts it reads (`dispatch:` = the dispatcher that routes
to it). Entries with phase `--` are stages of the non-step tools in `utils/`.

<!-- BEGIN GENERATED: pipeline_steps.yml via analysis/utils/make_pipeline_docs.R -- do not edit by hand -->

### Step table (75 entries, in run order)

| # | Step | Phase | Track | Rollout | Runs | After | Driver(s) | Writes |
|---|---|---|---|---|---|---|---|---|
| 1 | `00_ingest_outbreakextractr.R` | 00 | data | core | local+slurm | — | `stages/00_ingest.sh`, `campaigns/preprocess.sh` | ~/secure/cholera_outbreak_data/outbreak_data_&lt;Mon_YYYY&gt;/time_series_{preoutbreak,outbreak}_extraction.rds<br>~/secure/cholera_outbreak_data/outbreak_data_&lt;Mon_YYYY&gt;/outbreak_shapefiles.rds<br>generated_data/prep_counts/00_ingest_outbreakextractr.yml (cascade counts)<br>*input:* OutbreakExtractR outputs: combined_outbreaks_cholera.csv + stage1_geo_\*.geojson<br>*note:* Run when the outbreak data are refreshed; AFR/EMR filter (the &gt;=50-case filter is upstream, in OutbreakExtractR). Writes into the vintage folder (--outdir defaults to --indir = paths.R::defaultOutbreakDataDir()). |
| 2 | `00_ingest_synthetic.R` | 00 | data | optional | local | — | `tools/run_synthetic_pipeline.sh` | $CHOLERA_OUTBREAK_DATA_DIR/outbreak_data_&lt;Mon_YYYY&gt;/combined_outbreaks_cholera.csv + stage1_geo/\*.geojson (the input of 00_ingest_outbreakextractr.R)<br>data/geodata/{gadm_410-levels,gadm_valid}.gpkg, basemap shapefiles; data/climate_parquet/; data/worldpop/ + data/country_boundaries/<br>data/who_annual_report_africa/, data/best_seasonality_estimates.csv, data/hmm_country_clusters.csv, data/lenssen/extracted/, C3S plume stand-in<br>generated_data/enso_\*_cache.rds, generated_data/wpp2024_national_totals.csv, generated_data/synthetic_world/{truth.yml,epidemics.rds}<br>*input:* analysis/groupings/synthetic_countries.csv (fake countries, ISO 3166 user-assigned codes XA?)<br>*note:* The synthetic world: every raw input in the format of the real source, simulated with known ENSO effects (utils/synthetic_world.R), so the real 00 -&gt; 11 steps run without the restricted data. Runs only with CHOLERA_OUTBREAK_DATA_DIR set inside a sandbox (never ~/secure/); tools/run_synthetic_pipeline.sh chains every step on it with the deploy configs/deploy/synth_run1.yml. |
| 3 | `01_preprocess_outbreaks.R` | 01 | data | core | local+slurm | 00_ingest_outbreakextractr.R | `stages/01_preprocess.sh`, `campaigns/preprocess.sh` | ~/secure/cholera_outbreak_data/outbreak_metadata.rds<br>~/secure/cholera_outbreak_data/aggregated{,_admin_level}_ts_preoutbreak_extraction.rds<br>~/secure/cholera_outbreak_data/outbreak_data_&lt;Mon_YYYY&gt;/outbreak_shapefiles.gpkg<br>generated_data/prep_counts/01_preprocess_outbreaks.yml (cascade counts)<br>figures/outbreak_map.png (-f TRUE; the Slurm stage passes -f FALSE) |
| 4 | `01_preprocess_national_zero_months.R` | 01 | data | core | local+slurm | 01_preprocess_outbreaks.R | `stages/01_preprocess.sh`, `campaigns/preprocess.sh` | ~/secure/cholera_outbreak_data/national_zero_months.rds<br>~/secure/cholera_outbreak_data/national_zero_months_t&lt;N&gt;.rds<br>~/secure/cholera_outbreak_data/national_zero_months[_t&lt;N&gt;][_&lt;rule&gt;].audit.csv<br>~/secure/cholera_outbreak_data/national_zero_months[_t&lt;N&gt;]_db_who_consistency.csv<br>*note:* National (ADM0) + WHO annual zeros -&gt; country-months the coverage mask (coverage.R::build_coverage_mask) counts as observed 0 for every gid. --zero_threshold N (e.g. 60) adds WHO country-years under N cases (last --window months only if the next year is also under N; database-vs-WHO inconsistent years left out) as national_zero_months_t&lt;N&gt;.rds. --rule &lt;name&gt; (analysis/national_zero_rules.yml; r2 = head guard, consistency scope all, ETH 2010-2016) writes national_zero_months[_t&lt;N&gt;]_&lt;name&gt;.rds; every run writes the audit (&lt;months&gt;.audit.csv) and the database-vs-WHO table beside the months. The Slurm stage builds the deploy's deploy.episodes rule. |
| 5 | `01b_preprocess_outbreaks_monthly.R` | 01 | mon | core | local+slurm | 00_ingest_outbreakextractr.R, 01_preprocess_outbreaks.R | `stages/01_preprocess.sh`, `campaigns/preprocess.sh` | ~/secure/cholera_outbreak_data/outbreak_monthly_cases.rds<br>*note:* Month-grain re-aggregation of the same raw series; response of the monthly track. Gate B reads 01's outbreak metadata. |
| 6 | `01_preprocess_figures_data_overview.R` | 01 | data | optional | local | 00_ingest_outbreakextractr.R, 01_preprocess_outbreaks.R | — | figures/data_overview/<br>*note:* Manuscript Figure 1 (data overview). Reads 01's outbreak_shapefiles.gpkg. |
| 7 | `utils/preprocess_paths.R` | -- | data | core | local+slurm | — | `stages/01_preprocess.sh`, `campaigns/preprocess.sh` | KEY=value on stdout: OUTBREAK_VINTAGE / VINTAGE_DIR / SECURE_ROOT, NZ_RULE / NZ_THRESHOLD (--deploy), SEC= / GEN= output list (--outputs)<br>*note:* Bash reads the data-leg paths and the deploy's national-zero rule here; it never names the vintage. |
| 8 | `02a_preprocess_population_tasks.R` | 02 | pop | core | local | 00_ingest_outbreakextractr.R | `stages/02_population.sh` | task YAMLs, per country (--task_unit country) or per (country, year) (TASKS_DIR= / N_TASKS= on stdout)<br>*note:* Login node; the national-boundary cache and the WorldPop rasters are downloaded before the 02b array (stages/02_population.sh checks the boundary cache). |
| 9 | `02b_preprocess_population_run.R` | 02 | pop | optional | local+slurm | 02a_preprocess_population_tasks.R | `stages/02_population.sh` | generated_data/population/parts/iso3=&lt;ISO3&gt;/year=&lt;YYYY&gt;/<br>*input:* WorldPop constrained rasters; national boundaries<br>*note:* Per-(country, year) runner (stages/02_population.sh --task-unit country_year), kept beside the per-country one. |
| 10 | `02b_preprocess_population_run_country.R` | 02 | pop | core | local+slurm | 02a_preprocess_population_tasks.R | `stages/02_population.sh`, `campaigns/preprocess.sh` | generated_data/population/parts/iso3=&lt;ISO3&gt;/year=&lt;YYYY&gt;/ (same partitions as 02b_preprocess_population_run.R)<br>generated_data/population/gid_pop_cache/iso3=&lt;ISO3&gt;/year=&lt;YYYY&gt;.qs (Q_g, independent of the outbreak vintage)<br>*input:* WorldPop constrained rasters (prestaged); national boundaries<br>*note:* One task per country, all years: geometry once, Q_g reused across vintages. A partition is skipped only when built from the current LP corpus (popPartFresh). |
| 11 | `02c_preprocess_population_compile.R` | 02 | pop | core | local+slurm | 02b_preprocess_population_run_country.R | `stages/02_population.sh`, `campaigns/preprocess.sh` | generated_data/population_crosswalk_&lt;Y0&gt;-&lt;Y1&gt;.qs<br>generated_data/gid_population_&lt;Y0&gt;-&lt;Y1&gt;.qs<br>generated_data/population_qc/<br>*note:* The stage driver submits it afterok on the array. |
| 12 | `02_preprocess_crosswalk.R` | 02 | pop | core | local+slurm | 00_ingest_outbreakextractr.R, 01_preprocess_outbreaks.R, 02c_preprocess_population_compile.R | `stages/02_crosswalk.sh`, `campaigns/preprocess.sh` | generated_data/matched_ob_gids_crosswalk.rds<br>generated_data/u_gids_crosswalk.rds<br>generated_data/prep_counts/02_preprocess_crosswalk.yml (cascade counts)<br>*note:* Reads 01's outbreak_shapefiles.gpkg (attributes). |
| 13 | `02_preprocess_climate.R` | 02 | data | core | local+slurm | 00_ingest_outbreakextractr.R, 01_preprocess_outbreaks.R, 02_preprocess_crosswalk.R | `stages/02_climate.sh`, `campaigns/preprocess.sh` | generated_data/matched_ob_gids.rds (--stage match)<br>generated_data/climate_stats/ per-gid caches (--stage extract)<br>generated_data/climate_stats.sqlite, generated_data/climate_extremes.sqlite (--stage db)<br>*input:* data/climate_parquet/ (tools/download_jahn2026_parquet.sh) or data/climate_data/ (CHIRPS / ERA5-Land CSV, --climate_source csv)<br>*note:* Defaults = the deployed build (parquet, crosswalk gids, dynamic weighting); a database of other provenance is not overwritten without --force. Reads 01's outbreak_shapefiles.gpkg. |
| 14 | `03_prepare_occ.R` | 03 | occ | core | local+slurm | 01_preprocess_outbreaks.R, 01_preprocess_national_zero_months.R, 02_preprocess_climate.R, 02c_preprocess_population_compile.R | `stages/03_prepare.sh` | generated_data/outbreaks_survival_data_base-&lt;baseline&gt;.rds<br>generated_data/matched_ob_gids_seasonality_occ.rds<br>generated_data/outbreaks_survival_data_base-&lt;baseline&gt;&lt;variant&gt;.lp_members.csv<br>generated_data/outbreaks_survival_data_base-&lt;baseline&gt;&lt;variant&gt;.{adm1_overlap,discarded_onsets,adm0_dropped,masked_rows}.csv<br>generated_data/outbreaks_survival_data_base-&lt;baseline&gt;&lt;variant&gt;.contagion_coverage.csv<br>*note:* Refuses to overwrite an existing panel whose prep stamp differs (any field but stamped_at) unless --force. --national_zero_fill (gid_season only) reads 01_preprocess_national_zero_months.R's output and writes the _nz panel; with --national_zero_threshold N it reads the _t&lt;N&gt; months and writes the _nz&lt;N&gt; panel. |
| 15 | `03_prepare_siz.R` | 03 | siz | core | local+slurm | 01_preprocess_outbreaks.R, 02_preprocess_climate.R, 02c_preprocess_population_compile.R | `stages/03_prepare.sh` | generated_data/outbreak_size_analysis_data\*.rds<br>generated_data/outbreak_footprint_members\*.rds<br>generated_data/matched_ob_gids_seasonality_siz.rds<br>*note:* Refuses to overwrite an existing size frame or footprint map whose prep stamp differs (any field but stamped_at) or that has none, unless --force (the matched_ob_gids_seasonality_siz.rds cache is written before the guard). |
| 16 | `03_prepare_dur.R` | 03 | dur | core | local+slurm | 03_prepare_siz.R | `stages/03_prepare.sh` | generated_data/outbreaks_duration_data_base-&lt;baseline&gt;\*.rds<br>*note:* Refuses to overwrite an existing panel whose prep stamp differs (any field but stamped_at / restamped_at) or that has none, unless --force. |
| 17 | `03_prepare_mon.R` | 03 | mon | core | local+slurm | 03_prepare_dur.R, 01b_preprocess_outbreaks_monthly.R | `stages/03_prepare.sh` | generated_data/outbreaks_monthly_data_base-&lt;baseline&gt;\*.rds<br>*note:* Refuses to overwrite an existing panel whose prep stamp differs (any field but stamped_at) or that has none, unless --force. |
| 18 | `03_prepare_sim_panel.R` | 03 | sim | core | local+slurm | 03_prepare_occ.R | `stages/03_prepare.sh` | generated_data/sim_climate_panel_base-&lt;baseline&gt;.rds<br>*note:* Continuous per-(gid, month) historical panel (smoothing, hindcast, v1). Refuses to overwrite an existing panel whose prep stamp differs (any field but stamped_at) or that has none, unless --force; ENSO values are not in the stamp (back up and diff). |
| 19 | `03_prepare_who_monthly.R` | 03 | sim | core | local+slurm | — | `stages/03_prepare.sh` | generated_data/who_dashboard_monthly_africa.rds<br>generated_data/who_dashboard_monthly_africa_meta.yml<br>*input:* data/who_annual_report_africa/who_dashboard_weekly_africa.rds<br>*input:* data/who_annual_report_africa/who_sitrep_cumulative_africa.csv<br>*note:* WHO dashboard national weekly counts -&gt; monthly (day overlap), gap-filled from the WHO multi-country update cumulative snapshots (unchanged cumulative = zero only for --zero_iso3); scored by the smoothing's `_who` measurement (simulator.filter.who_obs). Refuses to overwrite an existing table whose _meta.yml (less `written` and the input mtimes) or content differs, or that has no meta, unless --force. |
| 20 | `03_prepare_sim_panel_projection.R` | 03 | sim | core | local+slurm | 03_prepare_occ.R | `stages/03_prepare.sh` | generated_data/sim_climate_panel_projection_base-&lt;baseline&gt;.rds<br>*input:* ONI cache + NOAA CPC CFSv2 ENSO forecast<br>*note:* Forecast-horizon panel past the climate DB's last month. Refuses to overwrite an existing panel whose prep stamp differs (any field but stamped_at) or that has none, unless --force; ENSO values are not in the stamp (back up and diff). |
| 21 | `03_prepare_sim_panel_extended.R` | 03 | sim | optional | local+slurm | 03_prepare_sim_panel.R, 03_prepare_sim_panel_projection.R | `stages/03_prepare.sh` | generated_data/sim_climate_panel_base-&lt;baseline&gt;[_ep-gs]_to-YYYY-MM.rds<br>*note:* Historical panel + projection rows up to the latest observed ONI month: the smoothing panel when the smoothing runs to the ENSO seam. Refuses to overwrite an existing panel whose prep stamp differs (any field but stamped_at) or that has none, unless --force. |
| 22 | `03_prepare_lenssen_analogs.R` | 03 | sim | core | local | — | — | generated_data/lenssen_nino34_analogs.rds<br>generated_data/lenssen_nino34_analogs_meta.yml<br>*input:* Lenssen et al. (2024) CMIP6 piControl Nino3.4, Zenodo 10045687<br>*note:* ENSO analog cache for the counterfactual forecast (enso.counterfactual.scenario: lenssen). --match legacy (default) = the original matcher on the Zenodo Rdas; --match library draws the pool from generated_data/enso_continuation_library.qs, which gen_forecast / gen_events also build in-process per run under deploy.enso_counterfactual.match: library. |
| 23 | `03_prepare_enso_cfsv2_calibration.R` | 03 | sim | optional | local | — | — | generated_data/enso_cfsv2_calibration.rds<br>generated_data/enso_cfsv2_calibration_meta.yml<br>*input:* NMME CFSv2 hindcast Nino3.4 (IRIDL PENTAD_SAMPLES_FULL, archived NetCDF)<br>*input:* CPC detrended ERSST and OISST Nino3.4 indices<br>*note:* CPC-style PDF + spread correction of the CFSv2 forecast (simulator.enso.cfsv2_calibration). |
| 24 | `03_prepare_enso_continuation_library.R` | 03 | sim | optional | local | — | — | generated_data/enso_continuation_library.qs<br>generated_data/enso_continuation_library_meta.yml<br>generated_data/enso_continuation_skill.csv (--validate)<br>generated_data/enso_continuation_skill_meta.yml (--validate)<br>*input:* Lenssen et al. (2024) CMIP6 piControl Nino3.4 (data/lenssen/extracted/ or the Zenodo 10045687 zips)<br>*input:* ONI cache (rescale sd; --validate: the observed seams)<br>*note:* Seam-agnostic k-NN analog library for ENSO past the CFSv2 span (simulator.enso.continuation.method: analog; also read by 03_prepare_sim_panel_projection.R --enso_continuation analog). --validate (no build) scores the library as deployed (--deploy's enso_continuation over the defaults) on observed seams: fair CRPS / CRPSS vs climatology and neutral by lead (utils/enso_continuation_skill.R), read by 09 for each group's continuation CRPSS. |
| 25 | `03_prepare_rain_adm1.R` | 03 | results | core | local+slurm | 02_preprocess_climate.R | `stages/03_prepare.sh` | generated_data/adm1_monthly_precip.qs<br>generated_data/adm1_monthly_precip_meta.yml<br>*note:* ADM1 x month rainfall reduced from climate_stats.sqlite (8.6 GB) for 09's rain_enso step (Figure 1c); built locally (copy the .qs to the cluster before 09 runs there) or on the cluster by stages/03_prepare.sh. |
| 26 | `utils/deploy_prod_run.R --stage gen_batches` | -- | all | core | local | — | — | analysis/configs/&lt;fits.\*_dir&gt;/\*.yml + manifest.yml<br>analysis/configs/&lt;kl_sweep.batch&gt;/\*.yml + manifest.yml (when the deploy derives its grouping from a kl sweep)<br>*note:* --stage gen_batches: the deploy YAML's fit-batch dirs (and kl class sweep dir), from make_batch_configs.R's BATCH_REGISTRY. A sweep rollout runs it twice: the sweep dir first (the fit dirs are 'deferred' while their grouping registry is absent), then the fit dirs after deploy:gen_grouping (run_pipeline.sh phases `classes` and `fits`). A fit batch grouped by another registry than deploy.grouping is refused. |
| 27 | `utils/make_batch_configs.R` | -- | all | optional | local | — | — | analysis/configs/&lt;batch&gt;/\*.yml + manifest.yml<br>*note:* Any registered batch outside a deploy (benchmarks, sensitivity). |
| 28 | `04_fit_run.R` | 04 | all | core | local+slurm | deploy:gen_batches, utils/make_batch_configs.R | `stages/04_fit.sh`, `campaigns/enso_grouped_singlelag.sh`, `campaigns/dur_enso_grouped_singlelag.sh`, `campaigns/mon_enso_grouped_singlelag.sh` | (routes to 04_fit_&lt;track&gt;_&lt;backend&gt;.R from the config's &lt;track&gt;.backend)<br>*note:* Dispatcher; slurm/submit_fit.sh runs it, then 05_fit_postprocess.R. |
| 29 | `04_fit_occ_brms_cb.R` | 04 | occ | core | local+slurm | 03_prepare_occ.R, (dispatch: 04_fit_run.R) | — | generated_data/ occurrence fit + cb_meta (makeRunFile) |
| 30 | `04_fit_occ_brms_kl.R` | 04 | occ | optional | local+slurm | 03_prepare_occ.R, (dispatch: 04_fit_run.R) | `stages/04_fit_kl.sh` | generated_data/ occurrence fit + kl_meta + kl_summary (+ kl_resp, kl_comembership, kl_partition for latent classes; kl_pin_plan, kl_pin + per-order _arrNN fits under kernel.pin) (makeRunFile)<br>*note:* ENSO lag kernel / latent country classes (utils/dlnm_kernel.R). A config batch (e.g. the class sweep kl_sweep_classes_prodtu_all, utils/make_batch_configs.R) on SLURM: stages/04_fit_kl.sh (array, -c 32 / 96G; --postprocess chains 05_fit_postprocess_kl_sweep.R; --deploy reads the batch from deploy.kl_sweep.batch and the comparison's chain selection from kl_sweep.n_agree / mode_tol). |
| 31 | `04_fit_siz_brms_cb.R` | 04 | siz | core | local+slurm | 03_prepare_siz.R, (dispatch: 04_fit_run.R) | — | generated_data/ size fit + cb_meta (makeSizeRunFile) |
| 32 | `04_fit_dur_brms_cb.R` | 04 | dur | core | local+slurm | 03_prepare_dur.R, (dispatch: 04_fit_run.R) | — | generated_data/ duration fit + cb_meta (makeDurationRunFile) |
| 33 | `04_fit_mon_brms_cb.R` | 04 | mon | core | local+slurm | 03_prepare_mon.R, (dispatch: 04_fit_run.R) | — | generated_data/ monthly fit + cb_meta |
| 34 | `04_fit_occ_brms.R` | 04 | occ | optional | local+slurm | 03_prepare_occ.R, (dispatch: 04_fit_run.R) | — | generated_data/outbreaks_survival_poisson_fit_\*.rds |
| 35 | `04_fit_siz_brms.R` | 04 | siz | optional | local+slurm | 03_prepare_siz.R, (dispatch: 04_fit_run.R) | — | generated_data/ size fit (makeSizeRunFile) |
| 36 | `04_fit_dur_brms.R` | 04 | dur | optional | local+slurm | 03_prepare_dur.R, (dispatch: 04_fit_run.R) | — | generated_data/ duration fit (makeDurationRunFile) |
| 37 | `04_fit_mon_brms.R` | 04 | mon | optional | local+slurm | 03_prepare_mon.R, (dispatch: 04_fit_run.R) | — | generated_data/ monthly fit |
| 38 | `05_fit_postprocess.R` | 05 | all | core | local+slurm | 04_fit_run.R | `stages/04_fit.sh` | generated_data/ RE stats, coef summaries, LOO, pp_check (make\*RunFile)<br>*note:* Runs in the same SLURM job as the fit (slurm/submit_fit.sh). |
| 39 | `05_fit_postprocess_cb.R` | 05 | all | core | local+slurm | 04_fit_occ_brms_cb.R, 04_fit_siz_brms_cb.R | `stages/05_postprocess_cb.sh`, `campaigns/postfit.sh`, `tools/run_batch_cb_local.sh` | generated_data/ crossbasis exposure-/lag-response surfaces |
| 40 | `05_fit_postprocess_loo_lags.R` | 05 | all | core | local+slurm | 04_fit_occ_brms_cb.R, 04_fit_siz_brms_cb.R, 04_fit_dur_brms_cb.R, 04_fit_mon_brms_cb.R | `stages/05_loo_compare.sh`, `campaigns/postfit.sh` | generated_data/loo_compare/&lt;batch&gt;/loo_compare_lags_&lt;group&gt;_&lt;YYYYMMDD&gt;.rds (+ png; ranking with per-lag Pareto-k columns)<br>generated_data/loo_compare/&lt;batch&gt;/lag_loo_sidecar_&lt;group&gt;_&lt;YYYYMMDD&gt;.rds (opt-in --lag_sidecar: pointwise elpd, Pareto-k, per-lag crossbasis draws)<br>*note:* Lag ranking per ENSO group; --sidecar_only writes only the sidecar (ranking untouched). |
| 41 | `05_fit_postprocess_lag_compare_surface.R` | 05 | all | optional | local | 05_fit_postprocess_loo_lags.R | — | generated_data/ best-lag exposure-response surface + Nina/Nino contrasts |
| 42 | `05_fit_postprocess_dependence.R` | 05 | all | optional | local | 04_fit_siz_brms_cb.R, 04_fit_dur_brms_cb.R | — | generated_data/ size-duration Kendall tau + joint/independent decision<br>*note:* One-off decision for the simulator's size/duration sampling. |
| 43 | `05_fit_postprocess_kl_sweep.R` | 05 | occ | optional | local+slurm | utils/make_batch_configs.R, 04_fit_occ_brms_kl.R | `stages/04_fit_kl.sh` | generated_data/kl_sweep_compare/&lt;batch&gt;/{fits,classes,partition,loo,loo_country,diagnostics,diagnostics_chains,responsibilities,country_mixture,stability}.csv + loo.rds<br>*input:* analysis/configs/&lt;sweep batch&gt;/\*.yml and the brms_kl fits they name<br>*note:* brms_kl class-sweep comparison: chain modes, agreeing chains, R-hat / ESS, country-block LOO + stacking weights, MAP partitions, membership uncertainty (responsibility matrix, M1 country mixture, partition stability). Read by 06_fit_figures_kl_sweep.R, deploy:gen_grouping and 09 report_kl_sweep (results.kl_sweep_dir). |
| 44 | `utils/deploy_prod_run.R --stage gen_grouping` | -- | occ | optional | local | 05_fit_postprocess_kl_sweep.R | — | analysis/groupings/&lt;deploy.grouping&gt;.csv + .yml (provenance header: deploy, batch, cell, chains, fit sha256)<br>generated_data/kl_groupings/&lt;deploy.grouping&gt;_membership.csv<br>*input:* generated_data/kl_sweep_compare/&lt;kl_sweep.batch&gt;/{fits.csv,partition.csv,loo.rds} + the fit of the pinned cell (+ diagnostics_chains.csv for the convergence report)<br>*note:* --stage gen_grouping (deploy.kl_sweep.batch set): the fit registry from the pinned cell deploy.kl_sweep.cell, each country in its class of maximum mean posterior probability (utils/kl_grouping.R, the make_kl_grouping.R body). The cell is never chosen automatically; its chains are the comparison's agreeing chains (kl_sweep.chains: agree; the comparison must have been run with the deploy's kl_sweep.n_agree / mode_tol: dominant mode, or the N best-agreeing chains) or a fixed list; deploy.fits.groups must name exactly the derived groups. Under deploy.convergence.action warn \| refuse it prints the cell's all-chain and agreeing-chain diagnostics (kl_grouping.R::kl_grouping_convergence()) and, under refuse, refuses when fewer than 2 chains agree or the agreeing-chain R-hat exceeds 1.05 (warns beyond 1.01 / 400). Absent -&gt; written, identical -&gt; left, different -&gt; refused unless --force. Consumed by deploy:gen_batches (the fit dirs) and every later deploy stage; first step of the `fits` phase of tools/run_pipeline.sh, after the `classes` phase (04_fit_kl.sh --deploy --postprocess). |
| 45 | `06_fit_figures_lag_compare.R` | 06 | all | core | local+slurm | 05_fit_postprocess_loo_lags.R, 05_fit_postprocess_cb.R | `stages/06_fit_figures_lag_compare.sh`, `campaigns/postfit.sh` | figures/ cross-lag comparison figures |
| 46 | `06_fit_figures.R` | 06 | all | optional | local | 05_fit_postprocess.R | — | figures/ per-config fit figures (makeFigureDir) |
| 47 | `06_fit_figures_cb.R` | 06 | all | optional | local | 05_fit_postprocess_cb.R | `tools/run_batch_cb_local.sh` | figures/cb/&lt;cfg-stem&gt;/ |
| 48 | `06_fit_figures_cb_adm1.R` | 06 | all | optional | local | 04_fit_occ_brms_cb.R, 04_fit_siz_brms_cb.R | — | figures/cb/ per-ADM1 crossbasis figures |
| 49 | `06_fit_figures_cb_country.R` | 06 | all | optional | local | 04_fit_occ_brms_cb.R, 04_fit_siz_brms_cb.R | — | figures/cb/ per-country crossbasis figures |
| 50 | `06_fit_figures_cb_grouped.R` | 06 | all | optional | local | 05_fit_postprocess_loo_lags.R | — | figures/cb/&lt;cfg-stem&gt;/grouped/ |
| 51 | `06_fit_figures_cb_grouped_composite.R` | 06 | all | optional | local | 05_fit_postprocess_loo_lags.R | — | figures/cb/&lt;batch-dir-name&gt;/composite/ |
| 52 | `06_fit_figures_kl_sweep.R` | 06 | occ | optional | local | 05_fit_postprocess_kl_sweep.R | — | figures/kl/&lt;batch&gt;/kl_sweep_\*.png (+ per_fit/&lt;config&gt;/ with --per_fit)<br>*note:* Figures of the class sweep from the 05 tables; --per_fit adds per-fit panels (traces, class centres, membership maps). |
| 53 | `utils/deploy_prod_run.R --stage deploy_fits` | -- | all | core | local | 05_fit_postprocess_loo_lags.R | — | generated_data/enso_best_lag_fits/ (canonically named best-lag fits + MANIFEST + CONVERGENCE.csv when the gate is on)<br>*note:* --stage deploy_fits (first part of --stage all). Convergence gate (deploy.convergence; absent = off): every source fit through diagnostics.R::fit_diagnostics_gate() (scope population: R-hat &lt;= 1.01 and bulk / tail ESS &gt;= 400 on all but r_/z_/zs_, R-hat &lt;= 1.05 on every parameter, max_divergent, E-BFMI &gt;= 0.3) BEFORE the prune; warn \| refuse, allow-listed &lt;track&gt;:&lt;group&gt; ids pass as allowed; --dry-run prints the verdicts and never stops. |
| 54 | `utils/deploy_prod_run.R --stage gen_templates` | -- | sim | core | local | deploy:deploy_fits | — | analysis/configs/&lt;run_tags.{historical,forecast,counterfactual}&gt;/&lt;run_tag&gt;.yml (spatPomp templates: analysis/configs/sim_templates/&lt;arm&gt;.yml + deploy sim: overlay)<br>*note:* --stage gen_templates (second part of --stage all); refuses to overwrite a template that differs unless --force. |
| 55 | `utils/deploy_prod_run.R --stage smoothing` | -- | sim | core | local | deploy:gen_templates | — | (prints the stages/08_smoothing.sh submit command for the historical template)<br>*note:* --stage smoothing (third part of --stage all); writes nothing. |
| 56 | `08a_smoothing_tasks.R` | 08 | sim | core | local | deploy:smoothing, 03_prepare_sim_panel.R, 03_prepare_who_monthly.R | `stages/08_smoothing.sh` | generated_data/sim_spatpomp/&lt;run_tag&gt;/tasks/task_NNNN.yml<br>generated_data/sim_spatpomp/&lt;run_tag&gt;/who_scale.csv + who_scale_meta.yml (who_obs.scale_mode: estimated only) |
| 57 | `08b_smoothing_spatpomp.R` | 08 | sim | core | slurm | 08a_smoothing_tasks.R | `stages/08_smoothing.sh` | generated_data/sim_spatpomp/&lt;run_tag&gt;/parts/ (one .qs + _SUCCESS per (TU, chunk)) |
| 58 | `08c_smoothing_compile.R` | 08 | sim | core | slurm | 08b_smoothing_spatpomp.R | `stages/08_smoothing.sh` | generated_data/sim_spatpomp/&lt;run_tag&gt;/cache/hist_filter_tu=&lt;TU&gt;_\*.qs (forecast handoff cache)<br>generated_data/sim_spatpomp/&lt;run_tag&gt;/merged/ (ribbons) |
| 59 | `utils/deploy_prod_run.R --stage gen_forecast` | -- | sim | core | local | 08c_smoothing_compile.R, deploy:deploy_fits | — | analysis/configs/forecast_country_by_group\*/ + counterfactual configs<br>analysis/configs/forecast_horizons.yml<br>per-run Lenssen counterfactual cache (simulator.enso.counterfactual.cache_path) when deploy.enso_counterfactual.match: library<br>*note:* --stage gen_forecast, via utils/make_forecast_configs.R; prints the submit command. Library-mode counterfactual analogs: builds / verifies the per-run cache in-process and records its digest in the manifest (enso_counterfactual). |
| 60 | `utils/deploy_prod_run.R --stage gen_attrib` | -- | sim | core | local | deploy:gen_forecast | `tools/run_pipeline.sh` | analysis/configs/forecast_country_by_group\*/forecast_&lt;unit&gt;_attrib.yml (factual launched at the counterfactual handoff)<br>analysis/configs/forecast_attribution_run\*.yml (attribution manifest, results.inputs.attrib_manifest)<br>*note:* --stage gen_attrib (deploy.outputs.forecast_attrib_manifest; no-op when unset), via utils/make_forecast_attrib_configs.R; prints the stages/07_forecast_groups.sh --arm attrib submit. run_pipeline.sh runs it in its sim phase, or alone with --from attrib. |
| 61 | `utils/deploy_prod_run.R --stage gen_hindcast` | -- | sim | core | local | 08c_smoothing_compile.R, deploy:deploy_fits | — | analysis/configs/hindcast_by_group\*/ + analysis/configs/hindcast_horizons.yml<br>*note:* --stage gen_hindcast, via utils/make_hindcast_configs.R; --align-results also runs deploy:align_results. |
| 62 | `utils/deploy_prod_run.R --stage gen_events` | -- | sim | optional | local | 08c_smoothing_compile.R, deploy:deploy_fits | `tools/run_pipeline.sh` | analysis/configs/event_attribution_run\*/ (factual + counterfactual per event x unit)<br>analysis/configs/event_attribution_horizons_run\*.yml (schema 3 with events.&lt;id&gt;.sub_windows when an event has them)<br>*input:* analysis/configs/enso_events_v2.yml (past-ENSO-event registry, with the one-event La Nina 2020-23 in per-winter sub-windows; deploy.event_attribution.registry)<br>*note:* --stage gen_events (deploy.event_attribution), via utils/make_event_attribution_configs.R; prints the per-event 03_prepare_lenssen_analogs.R builds and the stages/07_event_attribution.sh submit. run_pipeline.sh runs it in its sim phase when deploy.event_attribution.events is set (after gen_hindcast), preflights 07_event_attribution.sh --list and submits the event arrays with the forecast / hindcast ones; results waits on them (--event-replicates N -&gt; --n_replicates N, else each config's simulator.n_replicates). |
| 63 | `utils/deploy_prod_run.R --stage align_results` | -- | results | core | local | deploy:gen_forecast, deploy:gen_attrib, deploy:gen_hindcast | — | analysis/configs/results/manuscript.yml (or the --results-config file): deploy_config + the deploy-derived inputs.\* leaves (text edit, comments kept)<br>*note:* --stage align_results (or --align-results on gen_forecast / gen_hindcast); re-read + all-match verified; --dry-run prints the line diff and comments still naming another run. --results-config (also on tools/run_pipeline.sh) aligns another results config than manuscript.yml (a run reviewed before the headline moves to it). |
| 64 | `07b_simulator_spatpomp.R` | 07 | sim | core | slurm | deploy:gen_forecast, deploy:gen_attrib, deploy:gen_hindcast, deploy:gen_events, 08c_smoothing_compile.R, 03_prepare_sim_panel_projection.R, 03_prepare_lenssen_analogs.R, 03_prepare_enso_cfsv2_calibration.R, 03_prepare_enso_continuation_library.R | `stages/07_forecast_groups.sh`, `stages/07_hindcast.sh`, `stages/07_event_attribution.sh` | generated_data/sim_spatpomp/&lt;run_tag&gt;/ partitions<br>*note:* Forecast, attribution-factual (stages/07_forecast_groups.sh --arm attrib\|all), counterfactual and hindcast arrays; SLURM calls it directly (07_simulator_run.R is the local entry). |
| 65 | `07c_simulator_postprocess.R` | 07 | sim | core | local+slurm | 07b_simulator_spatpomp.R | `stages/07_forecast_groups.sh`, `stages/07_hindcast.sh` | &lt;sim outdir&gt;/&lt;run_tag&gt;/ aggregates + validation/metrics.rds |
| 66 | `07d_simulator_figures.R` | 07 | sim | core | local+slurm | 07c_simulator_postprocess.R | `stages/07_forecast_groups.sh`, `stages/07_hindcast.sh` | figures/sim/&lt;run_tag&gt;/ diagnostic figures |
| 67 | `07d_simulator_figures_counterfactual.R` | 07 | sim | core | local | 07c_simulator_postprocess.R | — | figures/sim/&lt;manifest-stem&gt;/counterfactual/ |
| 68 | `07_simulator_run.R` | 07 | sim | optional | local | 03_prepare_sim_panel.R | — | (routes to 07b_simulator_&lt;backend&gt;.R from simulator.backend)<br>*note:* Dispatcher; local entry point for both simulator backends. |
| 69 | `07a_simulator_tasks.R` | 07 | sim | optional | local | 03_prepare_sim_panel.R | — | generated_data/sim/&lt;run_tag&gt;/tasks/task_NNNN.yml<br>*note:* v1 path only; independent of 08. |
| 70 | `07b_simulator_v1.R` | 07 | sim | optional | local+slurm | 07a_simulator_tasks.R, 03_prepare_sim_panel.R, (dispatch: 07_simulator_run.R) | — | generated_data/sim/&lt;run_tag&gt;/ partitions<br>*note:* Pure-R historical replay; not part of the spatPomp rollout. |
| 71 | `09a_postprocess_results_tasks.R` | 09 | results | core | local+slurm | 07b_simulator_spatpomp.R, 07c_simulator_postprocess.R, 08b_smoothing_spatpomp.R, deploy:deploy_fits, deploy:gen_hindcast, deploy:align_results, 05_fit_postprocess_loo_lags.R, 05_fit_postprocess_kl_sweep.R, 03_prepare_rain_adm1.R | `stages/09_postprocess_results.sh` | generated_data/results/&lt;config stem&gt;/tasks/&lt;submission&gt;/task_NNNN.yml (one per results step, with its dependency wave; one folder per submission)<br>*note:* Builds the results context and writes one task per step of the registry in utils/results_derive.R; derives nothing. |
| 72 | `09b_postprocess_results.R` | 09 | results | core | local+slurm | 09a_postprocess_results_tasks.R | `stages/09_postprocess_results.sh` | generated_data/results/manuscript/\*.qs (stamped with the results: config)<br>generated_data/results/manuscript/{hindcast_runs,hindcast_occ,hindcast_cases,report_hindcast}_v1.qs (hindcast skill)<br>generated_data/results/manuscript/{convergence_v1_track-&lt;track&gt;,report_convergence,report_lags,report_anchor,report_software,report_kl_sweep}_v1.qs (supplementary tables)<br>generated_data/results/manuscript/{group_labels_v1,lag_averaging_v1_track-&lt;track&gt;,report_membership_v1}.qs (group display lags + S30 season ends; lag averaging from the 05 lag sidecars; membership uncertainty from the 05 kl-sweep comparison)<br>*note:* Run alone for a full serial derivation (no 09a needed), or per step (--steps / --task) as wave-chained arrays by the stage driver. Hindcast skill reads each manifest run's 07c --hindcast validation/metrics + merged summary. |
| 73 | `10_manuscript_figures.R` | 10 | results | core | local | 09b_postprocess_results.R | `stages/10_manuscript_figures.sh` | figures/manuscript/<br>figures/manuscript/hindcast/ (figure H)<br>figures/manuscript/supplementary/figS5_convergence, figS6_lag_loo (--figures S5,S6)<br>figures/manuscript/supplementary/figS_forecast_country_{cases,onsets}, figS_smoothing_country_{cases,onsets} (--figures FC,SC)<br>*note:* Plots only; refuses stale or decision-mismatched artefacts. |
| 74 | `10_manuscript_figures_legacy.R` | 10 | results | optional | local | 09b_postprocess_results.R | — | figures/manuscript_legacy/ (or $LEGACY_FIGURES_OUT_DIR; never figures/manuscript/)<br>*note:* Frozen copy of 10 (side-by-side rule): the earlier figure versions. Same inputs and CLI as 10; display only. |
| 75 | `11_manuscript_numbers.qmd` | 11 | results | core | local | 09b_postprocess_results.R | `stages/10_manuscript_figures.sh` | analysis/11_manuscript_numbers.html<br>figures/manuscript/tables/&lt;NN&gt;_&lt;slug&gt;.csv (one per table; results.style.tables_dir)<br>*note:* quarto render analysis/11_manuscript_numbers.qmd; display only. |

### Bash driver → steps

| Driver (`analysis/bash/`) | Runs |
|---|---|
| `campaigns/dur_enso_grouped_singlelag.sh` | 04_fit_run.R |
| `campaigns/enso_grouped_singlelag.sh` | 04_fit_run.R |
| `campaigns/mon_enso_grouped_singlelag.sh` | 04_fit_run.R |
| `campaigns/postfit.sh` | 05_fit_postprocess_cb.R, 05_fit_postprocess_loo_lags.R, 06_fit_figures_lag_compare.R |
| `campaigns/preprocess.sh` | 00_ingest_outbreakextractr.R, 01_preprocess_outbreaks.R, 01_preprocess_national_zero_months.R, 01b_preprocess_outbreaks_monthly.R, utils/preprocess_paths.R, 02b_preprocess_population_run_country.R, 02c_preprocess_population_compile.R, 02_preprocess_crosswalk.R, 02_preprocess_climate.R |
| `stages/00_ingest.sh` | 00_ingest_outbreakextractr.R |
| `stages/01_preprocess.sh` | 01_preprocess_outbreaks.R, 01_preprocess_national_zero_months.R, 01b_preprocess_outbreaks_monthly.R, utils/preprocess_paths.R |
| `stages/02_climate.sh` | 02_preprocess_climate.R |
| `stages/02_crosswalk.sh` | 02_preprocess_crosswalk.R |
| `stages/02_population.sh` | 02a_preprocess_population_tasks.R, 02b_preprocess_population_run.R, 02b_preprocess_population_run_country.R, 02c_preprocess_population_compile.R |
| `stages/03_prepare.sh` | 03_prepare_occ.R, 03_prepare_siz.R, 03_prepare_dur.R, 03_prepare_mon.R, 03_prepare_sim_panel.R, 03_prepare_who_monthly.R, 03_prepare_sim_panel_projection.R, 03_prepare_sim_panel_extended.R, 03_prepare_rain_adm1.R |
| `stages/04_fit_kl.sh` | 04_fit_occ_brms_kl.R, 05_fit_postprocess_kl_sweep.R |
| `stages/04_fit.sh` | 04_fit_run.R, 05_fit_postprocess.R |
| `stages/05_loo_compare.sh` | 05_fit_postprocess_loo_lags.R |
| `stages/05_postprocess_cb.sh` | 05_fit_postprocess_cb.R |
| `stages/06_fit_figures_lag_compare.sh` | 06_fit_figures_lag_compare.R |
| `stages/07_event_attribution.sh` | 07b_simulator_spatpomp.R |
| `stages/07_forecast_groups.sh` | 07b_simulator_spatpomp.R, 07c_simulator_postprocess.R, 07d_simulator_figures.R |
| `stages/07_hindcast.sh` | 07b_simulator_spatpomp.R, 07c_simulator_postprocess.R, 07d_simulator_figures.R |
| `stages/08_smoothing.sh` | 08a_smoothing_tasks.R, 08b_smoothing_spatpomp.R, 08c_smoothing_compile.R |
| `stages/09_postprocess_results.sh` | 09a_postprocess_results_tasks.R, 09b_postprocess_results.R |
| `stages/10_manuscript_figures.sh` | 10_manuscript_figures.R, 11_manuscript_numbers.qmd |
| `tools/run_batch_cb_local.sh` | 05_fit_postprocess_cb.R, 06_fit_figures_cb.R |
| `tools/run_pipeline.sh` | deploy:gen_attrib, deploy:gen_events |
| `tools/run_synthetic_pipeline.sh` | 00_ingest_synthetic.R |

<!-- END GENERATED -->
