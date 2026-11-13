# Tutorial: running the pipeline end to end

This tutorial runs the whole analysis twice: first on the bundled **synthetic
data**, on a laptop, then on the **real data**, locally or on a SLURM cluster.
Every command is run from the repository root.

## 1. Install

You need R 4.6, CmdStan (through `cmdstanr`), Quarto (for the phase-11 report)
and the packages pinned in `renv.lock`. Two packages are not on CRAN:

```r
remotes::install_github("HopkinsIDD/cholera-mapping-pipeline", subdir = "packages/taxdat")
remotes::install_github("HopkinsIDD/OutbreakExtractR")
cmdstanr::install_cmdstan()
```

To install the exact locked versions instead, use renv (opt-in):

```bash
CHOLERA_USE_RENV=1 Rscript -e 'renv::restore()'
```

Check the installation with the test harness (it skips the slow model fits):

```bash
Rscript analysis/tests/run_all.R --quick
```

## 2. The synthetic data, end to end

The outbreak database is not public. The synthetic data replace it: four fake
countries (codes `XAA`–`XAD`, invented names, polygons over open ocean) with
simulated outbreaks whose ENSO effects are known (`analysis/utils/synthetic_world.R`).

```bash
bash analysis/bash/tools/run_synthetic_pipeline.sh --cores 4
```

The driver creates a sandbox in `synthetic_world/` (the repository's `analysis/`
is linked into it) and runs these phases in order:

| Phase | What runs | Main outputs (in the sandbox) |
|---|---|---|
| `data` | 00 synthetic ingest, 01 outbreaks, 02 climate / crosswalk / population, 03 panels | `generated_data/outbreaks_survival_data_*.rds`, simulator panels |
| `fits` | 04 fits: 4 ENSO lags × 2 groups × 3 tracks | `generated_data/brms_cb_*/` |
| `postfit` | 05 postprocessing, LOO lag selection | `generated_data/loo_compare/` |
| `deploy` | best-lag selection, simulator templates, 08 historical smoothing | `generated_data/enso_best_lag_fits_*`, `generated_data/sim_spatpomp/` |
| `sim` | 07 forecasts, ENSO counterfactuals, hindcasts, past-event attribution | `generated_data/sim_spatpomp/` |
| `results` | 09 manuscript artefacts | `generated_data/results/manuscript_synth/` |
| `manuscript` | 10 figures, 11 key-numbers report | `figures/manuscript/`, `analysis/11_manuscript_numbers.html` |

It takes about 25 minutes on four cores. `--from <phase>` and `--until <phase>`
run a range, e.g. `--until data` to stop after the panels. After the data
phase it runs the production chain locally,
`run_pipeline.sh --deploy analysis/configs/deploy/synth_run1.yml --local --from fits`
(section 3). Logs are in `synthetic_world/logs/synthetic/`, one per job in its
`jobs/` folder.

The run is described by two YAML files:
[`analysis/configs/deploy/synth_run1.yml`](analysis/configs/deploy/synth_run1.yml)
(the deploy: fit batches, lag selection, simulator settings, experiments) and
[`analysis/configs/results/manuscript_synth.yml`](analysis/configs/results/manuscript_synth.yml)
(the results phase). Everything else is generated from the deploy; to list
every path of the run:

```bash
Rscript analysis/utils/deploy_prod_run.R -C analysis/configs/deploy/synth_run1.yml --stage resolve
```

## 3. The real data

### Inputs

1. **Outbreaks.** The Global Cholera Taxonomy Database extract (OutbreakExtractR
   stage-2 output: `combined_outbreaks_cholera.csv`, `stage1_geo/`). Put it in
   `<root>/outbreak_data_<Mon_YYYY>/` and set
   `export CHOLERA_OUTBREAK_DATA_DIR=<root>`; the vintage name is
   `OUTBREAK_DATA_VINTAGE` in `analysis/utils/paths.R`.
2. **Public inputs.** CHIRPS / ERA5-Land climate per GADM unit
   (`bash analysis/bash/tools/download_jahn2026_parquet.sh`), GADM 4.1
   (`data/geodata/gadm_410-levels.gpkg`), WorldPop constrained rasters
   (`data/worldpop/`), the WHO annual and dashboard series
   (`data/who_annual_report_africa/`), and the Lenssen et al. CMIP6 piControl
   Niño 3.4 archive (Zenodo 10045687, `data/lenssen/`). The ONI and the CFSv2
   ENSO forecast are downloaded by the pipeline.
3. **Seasonality estimates.** `data/best_seasonality_estimates.csv`, from
   [HopkinsIDD/cholera_seasonality_ssa](https://github.com/HopkinsIDD/cholera_seasonality_ssa/blob/main/generated_data/best_seasonality_estimates.csv).

### The production run

The manuscript run is one deploy YAML,
[`analysis/configs/deploy/prod_run.yml`](analysis/configs/deploy/prod_run.yml),
and its results config, [`analysis/configs/results/manuscript.yml`](analysis/configs/results/manuscript.yml).
`deploy_prod_run.R` generates from it the fit batches (115 fits per track:
ENSO lags −11..11 × 5 groups), the class sweep that derives the country
grouping, the simulator templates and the forecast, attribution, hindcast and
past-event experiments.

### The data phases (00–03)

The data phases build the model panels from the inputs above. One campaign
chains 00 ingest, 01 outbreaks, 02 population, crosswalk and climate; the 03
panels follow:

```bash
D=analysis/configs/deploy/prod_run.yml
bash analysis/bash/campaigns/preprocess.sh --check-packages
bash analysis/bash/campaigns/preprocess.sh --config analysis/configs/pop_extract/pop_extract.yml --deploy $D --plan
bash analysis/bash/campaigns/preprocess.sh --config analysis/configs/pop_extract/pop_extract.yml --deploy $D
bash analysis/bash/stages/03_prepare.sh --deploy $D
```

`--plan` prints the jobs and submits nothing. On a SLURM cluster, set up the
cluster profile first (next section); on a workstation, prefix each command with
`CHOLERA_BACKEND=local` (see [Locally](#locally)).

### On a SLURM cluster

Every driver submits through one library, `analysis/bash/lib/cluster.sh`, which
reads a cluster profile. Describe your cluster once:

```bash
cp analysis/bash/cluster.profile.example analysis/bash/cluster.profile
# edit: software environment (modules or conda), partitions per job class
# (SHORT, LONG, BIGMEM), account / qos, per-job overrides
```

Every variable is optional: without a profile, jobs go to the default partition
with each job script's own `#SBATCH` memory and time. To see the `sbatch` line a
job would get:

```bash
bash -c 'source analysis/bash/lib/cluster.sh; cholera_sbatch_cmd SHORT --job-name x analysis/bash/slurm/submit_results.sh; echo "${CHOLERA_SBATCH_CMD[*]}"'
```

Then plan and run the model chain:

```bash
bash analysis/bash/tools/run_pipeline.sh --deploy analysis/configs/deploy/prod_run.yml --plan
bash analysis/bash/tools/run_pipeline.sh --deploy analysis/configs/deploy/prod_run.yml --from fits
```

`run_pipeline.sh` is a relay: each phase submits its jobs, then submits the
next phase as a small job that waits (`afterok`) on them. The phases are
`classes` (the brms_kl class sweep that derives the country grouping), `fits`,
`postfit` (LOO lag rankings), `deploy` (best-lag selection, simulator templates,
historical smoothing), `sim` (forecasts, ENSO counterfactuals, hindcasts,
past-event attribution), `results` (09) and `manuscript` (10 + 11). `--from` and
`--until` restart or stop at a phase. `squeue -u $USER` shows the jobs; logs are
in `logs/` (`CHOLERA_LOG_DIR`).

The shipped `prod_run.yml` pins the manuscript's class-sweep cell and ships its
grouping, so `--from fits` starts directly. To re-derive the grouping, run
`--until classes`, read the sweep comparison
(`generated_data/kl_sweep_compare/<batch>/`), set `kl_sweep.cell` in the deploy,
then continue with `--from fits` (the `fits` phase writes the grouping).

### Locally

The same drivers run without SLURM: `--local` replaces job submission with
local processes. Array tasks run in parallel (`--cores N`), phases one after the
other; no profile is needed.

```bash
bash analysis/bash/tools/run_pipeline.sh --deploy analysis/configs/deploy/prod_run.yml --local --cores 8 --from fits
```

The production run fits 345 models and runs thousands of forward simulations
per experiment; budget days of CPU time on a workstation. The synthetic run of
section 2 uses exactly this path.

## 4. Results

Phase 09 derives every manuscript artefact into `generated_data/results/manuscript/`,
phase 10 draws the figures into `figures/manuscript/`, and phase 11 renders the
key-numbers report:

```bash
Rscript analysis/09b_postprocess_results.R -C analysis/configs/results/manuscript.yml --list   # status of every step
Rscript analysis/09b_postprocess_results.R -C analysis/configs/results/manuscript.yml
Rscript analysis/10_manuscript_figures.R  -C analysis/configs/results/manuscript.yml
quarto render analysis/11_manuscript_numbers.qmd
```

Phase 09 re-derives only what is stale: an artefact records its inputs and the
results settings it was derived with, and `--list` shows which are out of date.

## 5. Troubleshooting

- **`run from a synthetic sandbox`**: `00_ingest_synthetic.R` writes fake inputs
  over `data/` and `generated_data/`, so it refuses to run outside a directory
  holding `.synthetic_sandbox`. Use `run_synthetic_pipeline.sh`, which creates one.
- **A fit stops on a missing panel**: run the 03 step named in the error; the
  panel name encodes the episode definition of the deploy (`episodes:` block).
- **Phase 09 says an input is stale**: an upstream file changed after the
  artefact was derived; re-run 09 (it re-derives the affected steps only).
- **Cluster jobs pending on a partition that does not exist**: set the
  partitions in your `analysis/bash/cluster.profile`.
