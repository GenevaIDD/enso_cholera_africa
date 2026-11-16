# ENSO and cholera in Africa

Code for the analysis of:

Perez-Saez, J., Lee, E.C., Zheng, Q., Kamadjeu, R., Zaitchik, B.F. and Azman, A.S., 2026. El Nino Southern Oscillation shapes cholera in Africa. [Preprint] medRxiv. Available at: https://www.medrxiv.org/content/10.64898/2026.11.11.26367838v1.

The analysis estimates the effect of ENSO (the Oceanic Niño Index, ONI) on three outcomes: cholera outbreak occurrence, duration (cessation), and monthly incidence rates.

| Track | Outcome | Model |
|---|---|---|
| `occ` | outbreak occurrence (onset hazard) | Poisson GLMM per unit-month at risk |
| `dur` | outbreak duration (cessation hazard) | discrete-time Poisson GLMM per active episode-month |
| `mon` | monthly case counts within an outbreak | negative-binomial GLMM per active episode-month |

The deployed model terms and priors are listed in [`analysis/model_terms_deployed.yml`](analysis/model_terms_deployed.yml).
The fits feed a spatially coupled stochastic outbreak simulator used for the forecasts, the attribution and the hindcasts. Historical states are reconstructed by block particle filtering as implemented in
[spatPomp](https://cran.r-project.org/package=spatPomp).

## Supplementary data

[`generated_data/supplementary_data/`](generated_data/supplementary_data/) holds the
tables cited as "GitHub data folder" in the Supplementary Materials, one CSV per
table, with a [`manifest.csv`](generated_data/supplementary_data/manifest.csv) of
their contents. The outbreak episode records and the model inputs derived from
them are not distributed.

## Run it on synthetic data

The outbreak surveillance database is not public, so the repository ships a
**synthetic dataset**: four fake countries with simulated outbreaks, climate,
population, WHO series and ENSO, written in the formats of the real inputs.
Every step of the pipeline runs on it, from the raw-data ingest to the manuscript
figures and key numbers:

```bash
bash analysis/bash/tools/run_synthetic_pipeline.sh --cores 4
```

It builds a sandbox in `synthetic_world/` and takes about 25 minutes on a
four-core laptop. [`TUTORIAL.md`](TUTORIAL.md) walks through it phase by phase,
and then through the real-data run on a SLURM cluster or a workstation. `--until data` stops after the
model-ready panels; `--from` / `--until` run any range of phases (`data`, `fits`,
`postfit`, `deploy`, `sim`, `results`, `manuscript`). The synthetic data and their known ENSO
effects are described in [`analysis/utils/synthetic_world.R`](analysis/utils/synthetic_world.R);
the run is the deploy config [`analysis/configs/deploy/synth_run1.yml`](analysis/configs/deploy/synth_run1.yml)
and its results config [`analysis/configs/results/manuscript_synth.yml`](analysis/configs/results/manuscript_synth.yml).

## Repository layout

```
analysis/                  all code
  NN[letter]_<stem>_*.R    runnable pipeline steps, phases 00-11
  utils/                   sourced-only modules and non-step tools (see utils/README.md)
  configs/                 deploy/ run YAMLs, results/, sim_templates/, pop_extract/
  groupings/               ISO3 -> group registries (.csv) with .yml label / colour sidecars
  tables/                  hand-written tables read by the results phase
  data_manifests/          download manifests of public inputs
  bash/                    pipeline drivers: tools/run_pipeline.sh (relay), stages/, campaigns/,
                           slurm/ payloads, lib/ (cluster profile, local backend)
  tests/                   test harness (run_all.R)
  pipeline_steps.yml       machine-readable run graph
  PIPELINE.md              run order, generated from pipeline_steps.yml
notebooks/pipeline_run_graph.dot   file-level run graph (generated)
renv.lock, renv/, .Rprofile        package lockfile and opt-in renv activation
generated_data/supplementary_data/ the supplementary data tables
```

Not versioned, created by the pipeline: `data/` (raw inputs), `generated_data/`
(derived datasets, fits, simulations, results; all but `supplementary_data/`), `figures/`.

## Pipeline

Scripts follow one naming rule: **a number means it runs; no number means it is
sourced.** `NN` is the phase, an optional letter is an ordered substep, then the
phase verb and qualifiers (track, backend, variant).

| Phase | Stem | Concern |
|---|---|---|
| 00 | `ingest` | external tool outputs → repo-native RDS (or the synthetic data) |
| 01 | `preprocess` | outbreak time series |
| 02 | `preprocess` | climate, location-period crosswalk, population (`02a`→`02b`→`02c`) |
| 03 | `prepare` | model-ready datasets (occ / dur / mon panels, simulator panels, ENSO inputs) |
| 04 | `fit` | model fitting, through the `04_fit_run.R` dispatcher |
| 05 | `fit_postprocess` | postprocessing of the fits, LOO lag selection |
| 06 | `fit_figures` | figures from the fits |
| 07 | `simulator` | forward simulator (`07a` tasks → `07b` kernels → `07c` postprocess → `07d` figures) |
| 08 | `smoothing` | historical smoothing (`08a` tasks → `08b` filter → `08c` compile) |
| 09 | `postprocess_results` | derives the manuscript artefacts |
| 10 | `manuscript_figures` | plots the phase-09 artefacts |
| 11 | `manuscript_numbers` | key-numbers report (Quarto) |

The full run order, with dependencies, inputs and outputs, is in
[`analysis/PIPELINE.md`](analysis/PIPELINE.md). A run is named in one place, a
deploy YAML: the manuscript run is [`analysis/configs/deploy/prod_run.yml`](analysis/configs/deploy/prod_run.yml),
the synthetic run `synth_run1.yml`. Its fit batches, simulator configs and results
inputs are generated from it by `analysis/utils/deploy_prod_run.R` (`--stage gen_batches`,
`deploy_fits`, `gen_templates`, `gen_forecast`, `gen_attrib`, `gen_hindcast`,
`align_results`; `--stage resolve` prints every path of a run).
[`analysis/bash/tools/run_pipeline.sh`](analysis/bash/tools/run_pipeline.sh)
runs a deploy phase by phase on a SLURM cluster (described once in
`analysis/bash/cluster.profile`) or locally (`--local`).

## Running on the real data

1. **Outbreaks.** The episodes come from the Global Cholera Taxonomy Database,
   extracted with [OutbreakExtractR](https://github.com/HopkinsIDD/OutbreakExtractR)
   and aggregated with `taxdat` ([cholera-mapping-pipeline](https://github.com/HopkinsIDD/cholera-mapping-pipeline),
   `packages/taxdat`). Place the
   OutbreakExtractR stage-2 output (`combined_outbreaks_cholera.csv`,
   `stage1_geo/`) in `<root>/outbreak_data_<Mon_YYYY>/` and set
   `CHOLERA_OUTBREAK_DATA_DIR=<root>` (default `~/secure/cholera_outbreak_data`);
   the vintage is `OUTBREAK_DATA_VINTAGE` in `analysis/utils/paths.R`.
2. **Public inputs.** CHIRPS / ERA5-Land climate per GADM unit
   (`bash analysis/bash/tools/download_jahn2026_parquet.sh`), GADM 4.1
   (`data/geodata/gadm_410-levels.gpkg`), WorldPop constrained rasters
   (`data/worldpop/`), WHO annual and dashboard series (`data/who_annual_report_africa/`),
   the ONI and the CFSv2 ENSO forecast (downloaded by `analysis/utils/enso.R`),
   and the Lenssen et al. (2024) CMIP6 piControl Niño 3.4 archive
   (Zenodo 10045687, `data/lenssen/`).
3. **Seasonality estimates.** `data/best_seasonality_estimates.csv` (the ADM1
   seasonal profiles of a previous analysis), from
   [HopkinsIDD/cholera_seasonality_ssa](https://github.com/HopkinsIDD/cholera_seasonality_ssa/blob/main/generated_data/best_seasonality_estimates.csv).
4. Run the deploy `analysis/configs/deploy/prod_run.yml` (see [`TUTORIAL.md`](TUTORIAL.md)):

```bash
bash analysis/bash/tools/run_pipeline.sh --deploy analysis/configs/deploy/prod_run.yml --plan
bash analysis/bash/tools/run_pipeline.sh --deploy analysis/configs/deploy/prod_run.yml
```

## Software

R 4.6 with the packages pinned in [`renv.lock`](renv.lock); CmdStan through
`cmdstanr`; Quarto for phase 11. Two packages are not on CRAN:

```r
remotes::install_github("HopkinsIDD/cholera-mapping-pipeline", subdir = "packages/taxdat")
remotes::install_github("HopkinsIDD/OutbreakExtractR")
```

renv is opt-in: `.Rprofile` activates it only when `CHOLERA_USE_RENV=1`.

```bash
CHOLERA_USE_RENV=1 Rscript -e 'renv::restore()'                 # once: install the locked versions
CHOLERA_USE_RENV=1 Rscript analysis/tests/run_all.R --quick      # the test harness on the renv library
```

`Rscript analysis/tests/run_all.R --quick` runs the test harness (no model fit;
the hard gates must pass); without `--quick` it also runs the synthetic
data phase and the slower tests.

The deployed fits used the module set `GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1` on the
cluster; the results phase (09-11) `R/4.3.2` with GDAL 3.7.1. Each fit records the
versions of brms, rstan, cmdstanr and CmdStan it ran with, and phase 09 collects
them (`report_software`).

## License

GPL-3.0, see [LICENSE](LICENSE).

## Contacts

- javier.perez@unige.ch
- andrew.azman@unige.ch
