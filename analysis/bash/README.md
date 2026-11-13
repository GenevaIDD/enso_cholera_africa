# analysis/bash

Shell drivers of the pipeline (run order: `analysis/PIPELINE.md`). They run every phase as
SLURM jobs on a cluster described by a profile, or on this machine with the local backend.
Run every script from the repository root: `bash analysis/bash/<tier>/<script>.sh ...`.

## Cluster profile and backends

Every driver and payload sources `lib/cluster.sh`, which reads the cluster profile
`${CHOLERA_CLUSTER_PROFILE:-analysis/bash/cluster.profile}`:

```bash
cp analysis/bash/cluster.profile.example analysis/bash/cluster.profile   # then fill it in
```

- `cholera_env [--purge] [STACK]`: the software environment loaded before every `Rscript`
  (`CHOLERA_ENV_SETUP`, or `CHOLERA_ENV_SETUP_<STACK>` for the `GEO`, `MANUSCRIPT` and
  `REPORTS` stacks).
- `cholera_sbatch CLASS ...`: `sbatch` plus the profile's partition, account, qos and sizes
  for the job class: `SHORT` (<= 12 h), `LONG` (days), `BIGMEM` (smoothing). Per-payload
  overrides: `CHOLERA_RES_<KEY>`. Precedence: the driver's explicit options >
  `CHOLERA_RES_<KEY>` > the driver's defaults > the class > the payload's `#SBATCH` lines.

Without a profile, jobs go to the cluster's default partition with each payload's `#SBATCH`
memory, time and cores, and nothing is loaded before `Rscript`.

**Local backend.** `bash analysis/bash/tools/run_pipeline.sh --deploy D.yml --local [--cores N]`
(or `CHOLERA_BACKEND=local`) runs everything on this machine: `lib/local/sbatch` stands in for
`sbatch` (array tasks `N` at a time, logs in `logs/`), and the phases run one after the other.

## Tiers

| Tier | Directory | Role |
|---|---|---|
| 0 | `slurm/` | sbatch payload: `#SBATCH` header (mem / time / cpus), `cholera_env`, one `Rscript` call. |
| 1 | `stages/` | Driver of one R phase: builds its tasks, submits tier-0 payloads, chains them with `--dependency=afterok`. Prefix = the R phase. |
| 2 | `campaigns/` | Chains tier-1 drivers: the 00-02 data leg (`preprocess.sh`), the fits of the three tracks, the postfit steps. |
| — | `tools/` | `run_pipeline.sh` (the whole run, phase by phase), `run_synthetic_pipeline.sh`, local helpers. |
| — | `lib/` | Sourced only: `cluster.sh` and the local backend (`local/sbatch`, `local/squeue`). |

## Entry points

| Script | Purpose |
|---|---|
| `tools/run_pipeline.sh` | A production run from its deploy YAML: fits → postfit → deploy (smoothing) → sim (forecast, attribution, counterfactual, hindcast, past ENSO events) → results (09) → manuscript (10 + 11). On SLURM each phase submits its jobs and queues the next phase as a relay job; `--local` runs them inline. `--plan` prints the phases and submits nothing; `--from` / `--until` restart or stop at a phase. |
| `campaigns/preprocess.sh` | The 00-02 data leg (ingest, outbreaks, population, crosswalk, climate); then `stages/03_prepare.sh` builds the model panels. `--plan` prints the job chain. |
| `tools/run_synthetic_pipeline.sh` | The whole pipeline on the synthetic world, on this machine: its data phase (00 synthetic → 03), then `run_pipeline.sh --deploy analysis/configs/deploy/synth_run1.yml --local --from fits --until manuscript` (`--from` / `--until`, `--cores`). |
| `tools/run_batch_cb_local.sh` | Postprocess (crossbasis surfaces) and draw one batch of fits locally. |
| `tools/download_jahn2026_parquet.sh` | Download the CHIRPS / ERA5-Land admin-unit parquet corpus (Jahn et al. 2026) listed in `analysis/data_manifests/`. |

Rules: bash never builds artefact names (R prints them as `KEY=value` lines, e.g.
`deploy_prod_run.R --stage resolve`), and a driver that needs a run takes `--deploy`.
