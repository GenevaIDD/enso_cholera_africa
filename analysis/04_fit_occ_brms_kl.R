# 04_fit_occ_brms_kl.R — occurrence Poisson GLMM with a Gaussian ENSO lag kernel
#
# Usage:
#   Rscript analysis/04_fit_run.R \
#     -C analysis/configs/smoke_kl_occ/smoke_kl_occ.yml
#   Rscript analysis/04_fit_occ_brms_kl.R \
#     -C analysis/configs/smoke_kl_occ/smoke_kl_occ_k2.yml \
#     --chains 2 --sampling 100 --warmup 100
#
# Details:
#   Selected by occurrence.backend: brms_kl (dispatched by 04_fit_run.R). The
#   single-lag ENSO onebasis is replaced by a normalised Gaussian kernel over
#   the dlnm lags: range; the best lag is the posterior of the centre kmu.
#   kernel.classes >= 2 or kernel.prior_grouping: K latent country classes.
#   Writes (make_run_file kinds) fit, kl_meta, kl_summary; latent classes
#   also kl_resp, kl_comembership, kl_partition.
#
# See: analysis/utils/dlnm_kernel.R

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(brms)
  library(splines2)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")
source("analysis/utils/priors.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/composite_lp.R") # transitive dependency of population.R
source("analysis/utils/population.R") # .apply_p3() -- canonical_admin_level (load_surv_data)
source("analysis/utils/siz_data.R")
source("analysis/utils/enso.R")
source("analysis/utils/provenance.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/dlnm_kernel.R") # lag kernel + latent classes
source("analysis/utils/loo_lags.R") # kernel.pin.order_source sweep (read_latest_ranking())
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-a", "--admin_level"),
    default = NULL,
    type = "character",
    help = "Admin level (ADM1 | ADM2) [config default: ADM2]"
  ),
  make_option(
    c("-b", "--baseline"),
    default = NULL,
    type = "character",
    help = "Baseline period [config default: 1990_2019]"
  ),
  make_option(
    c("-c", "--country"),
    default = NULL,
    type = "character",
    help = "ISO3 code(s), comma-separated (e.g. 'ETH' or 'ETH,KEN'), or 'all'"
  ),
  make_option(
    c("--basehaz_poisson"),
    default = NULL,
    type = "character",
    help = "Baseline hazard: factor|cyclic|seasonal|ms [config default: seasonal]"
  ),
  make_option(
    c("--sampling"),
    default = NULL,
    type = "integer",
    help = "Post-warmup sampling iterations [config default: 1000]"
  ),
  make_option(
    c("--warmup"),
    default = NULL,
    type = "integer",
    help = "Warmup iterations [config default: 1000]"
  ),
  make_option(
    c("--chains"),
    default = NULL,
    type = "integer",
    help = "MCMC chains [config default: 4]"
  ),
  make_option(
    c("--cores"),
    default = NULL,
    type = "integer",
    help = "Parallel cores [config default: 4]"
  ),
  make_option(
    c("--seed"),
    default = NULL,
    type = "integer",
    help = "RNG seed [config default: 1234]"
  ),
  make_option(
    c("--brms_init"),
    default = NULL,
    type = "character",
    help = "Init strategy (stage 1): random | 0 | <number> [config default: random]"
  ),
  make_option(
    c("--overwrite"),
    default = NULL,
    action = "store_true",
    help = "Re-fit even if output file already exists [default: FALSE]"
  ),
  make_option(
    c("--pin_stage"),
    default = "all",
    type = "character",
    help = paste(
      "kernel.pin only: all (one process) | plan (seeds + orders -> kl_pin_plan) |",
      "run (one order, --pin_task) | combine (stacking + artefacts) [default: all]"
    )
  ),
  make_option(
    c("--pin_task"),
    default = NULL,
    type = "integer",
    help = "kernel.pin --pin_stage run: 1-based order index (the SLURM array task id)"
  ),
  make_option(
    c("--compile_only"),
    default = FALSE,
    action = "store_true",
    help = paste(
      "Engine stan: build and compile this config's program, then exit without",
      "sampling (the batch precompile, utils/kl_compile_batch.R) [default: FALSE]"
    )
  )
)

opt <- make_options_from_config(
  option_list = option_list,
  enforce_options = c(
    "country",
    "basehaz_poisson",
    "sampling",
    "warmup",
    "chains",
    "cores",
    "seed",
    "brms_init",
    "overwrite"
  ),
  model = "occurrence"
)
stopifnot(opt$occurrence$backend == "brms_kl")
print_options(opt)

# brms_kl falls back to the brms sub-block for MCMC settings
brms_cfg <- opt$occurrence$brms_kl %||% opt$occurrence$brms %||% list()

basehaz <- .occ_basehaz(brms_cfg)
stopifnot(basehaz %in% c("factor", "cyclic", "ms"))

kcfg <- resolve_kernel_cfg(opt)
latent <- kl_is_latent(kcfg)
engine <- kl_engine(kcfg)
if (latent && isTRUE(opt$grouping$enable)) {
  stop(
    "brms_kl: latent classes (kernel.classes >= 2 or kernel.prior_grouping) ",
    "infer the country partition -- disable grouping.enable and use `country:` ",
    "to restrict the fitted countries."
  )
}

# 3. Early exit if fit already exists ------------------------------------------

fit_file <- make_run_file(opt, kind = "fit")
check_run_file_names(opt, "occurrence") # every artefact name fits NAME_MAX before the fit (paths.R)
if (file.exists(fit_file) && !isTRUE(opt$overwrite)) {
  message("Fit already exists (use --overwrite to re-fit): ", fit_file)
  quit(save = "no", status = 0)
}

# 4. Load data + Poisson view --------------------------------------------------
# prepare_occ_brms_data() (utils/occ_data.R): load_surv_data() -> episode
# stamp -> base-rate anchor sidecar -> .occ_poisson_view() -> M-spline columns.

prep <- prepare_occ_brms_data(opt, basehaz)
opt <- prep$opt
pois_dat <- prep$pois_dat
ep_rec <- prep$ep_rec

if (nrow(pois_dat) == 0L) {
  message("No rows in pois_dat after country/grouping filter — skipping fit.")
  quit(save = "no", status = 0)
}

# 5. Model: basis, formula terms, nl formula, priors ---------------------------
# Kernel over the dlnm `lags:` range (negative = ENSO lead):
#   ENSO term = sum_l w_l(kmu, ksd) * f(ENSO at lag l),
# with f one natural-spline onebasis shared by every lag.
# kernel.classes <= 1 and no prior_grouping: one kernel (grouping mode allowed).
# Otherwise K latent country classes, each with its kernel, marginalised per
# country; optional prior partition (kernel.prior_grouping) and ICAR membership
# (kernel.spatial).
# kernel.engine: stan (default; brms's program patched and vectorised, any K,
# cmdstanr; fit = kl_fit) or brms (brm() on the nl formula, single kernel, a
# slow cross-check; fit = brmsfit).

m <- assemble_kl_model(
  opt,
  pois_dat,
  brms_cfg,
  .occ_baseline_block(basehaz, prep$ms_basis),
  kcfg
)
message(
  "Linear part:\n  lin ~ ",
  paste(
    c(
      "1",
      .occ_baseline_block(basehaz, prep$ms_basis),
      m$terms$fixed_terms,
      m$terms$re_terms,
      m$terms$off_terms
    ),
    collapse = " + "
  )
)
message(sprintf(
  "Kernel: lags %d..%d (%d), J = %d basis functions, %s, engine %s",
  kcfg$lag_min,
  kcfg$lag_max,
  length(kcfg$lags),
  m$meta$df,
  if (latent) "latent classes" else "single kernel",
  engine
))

kl_meta <- list(
  meta = m$meta,
  kcfg = kcfg,
  cb_meta = m$cb_meta,
  lin_terms = m$terms,
  latent = latent,
  engine = engine
)

# 6. Fit -----------------------------------------------------------------------
t0 <- Sys.time()
if (engine == "brms" && isTRUE(opt$compile_only)) {
  message(
    "--compile_only: kernel.engine brms compiles inside brms::brm(); nothing precompiled."
  )
  quit(save = "no", status = 0)
}
if (engine == "brms") {
  # brms non-linear formula (single kernel) ---------------------------------
  mcmc <- resolve_brms_mcmc(brms_cfg, label = "brms_kl")
  fit <- brms::brm(
    formula = m$formula,
    data = m$data,
    family = m$family,
    prior = m$priors,
    chains = brms_cfg$chains %||% 4L,
    iter = mcmc$iter,
    warmup = mcmc$warmup,
    cores = brms_cfg$cores %||% 4L,
    seed = opt$seed,
    init = resolve_brms_init(brms_cfg),
    refresh = mcmc$refresh,
    backend = "cmdstanr",
    threads = resolve_brms_threads(brms_cfg),
    control = list(
      adapt_delta = brms_cfg$adapt_delta %||% 0.95,
      max_treedepth = brms_cfg$max_treedepth %||% 12L
    )
  )
  summ <- summarise_kl_fit(fit, m$meta, cen = kcfg$cen)
} else {
  # Patched vectorised program (single kernel or latent classes) -----------
  if (
    !is.null(brms_cfg$threads_per_chain) &&
      as.integer(brms_cfg$threads_per_chain) > 1L
  ) {
    message(
      if (identical(kcfg$form, "discrete")) {
        "brms_kl engine stan: threads_per_chain ignored for kernel.form discrete (lag combinations couple the countries)."
      } else {
        sprintf(
          "brms_kl engine stan: reduce_sum over countries, %d threads per chain.",
          as.integer(brms_cfg$threads_per_chain)
        )
      }
    )
  }
  m <- kl_sort_by_country(m)
  ctry_levels <- levels(droplevels(as.factor(m$data$country)))
  pg <- kl_prior_partition(kcfg, ctry_levels)
  adj <- if (isTRUE(kcfg$spatial$enable)) {
    kl_adjacency_data(
      load_adjacency(kcfg$spatial$adjacency %||% "country_adjacency"),
      ctry_levels
    )
  } else {
    NULL
  }
  if (!is.null(adj) && length(adj$isolated)) {
    message(
      "ICAR: countries with no fitted neighbour (field pinned at 0): ",
      paste(adj$isolated, collapse = ", ")
    )
  }
  message(sprintf(
    "Classes: K = %d (%s), identify = %s, %d countries%s%s",
    pg$K,
    paste(pg$class_labels, collapse = ", "),
    kl_identify(kcfg),
    length(ctry_levels),
    if (pg$use) {
      sprintf(", prior partition '%s'", kcfg$prior_grouping$source)
    } else {
      ""
    },
    if (!is.null(adj)) {
      sprintf(", ICAR on %d edges / %d components", adj$n_edges, adj$n_comp)
    } else {
      ""
    }
  ))
  stan_dir <- file.path(opt$outdir %||% "generated_data/", "stan_kl")
  # Batch precompile (bash/stages/04_fit_kl.sh): compile the program the fit
  # would compile, then exit, so the array tasks find the executable built.
  if (isTRUE(opt$compile_only)) {
    if (kl_pin_cfg(kcfg)$enable) {
      message(
        "--compile_only: kernel.pin compiles per order; nothing precompiled."
      )
    } else {
      prog <- kl_compile_model(m, brms_cfg, stan_dir)
      message(
        "--compile_only: compiled ",
        prog$file,
        " (",
        prog$threads,
        " thread(s) per chain)"
      )
    }
    quit(save = "no", status = 0)
  }
  pin <- NULL
  pin_stage <- opt$pin_stage %||% "all"
  if (!pin_stage %in% c("all", "plan", "run", "combine")) {
    stop(
      "--pin_stage must be all | plan | run | combine (got '",
      pin_stage,
      "')."
    )
  }
  if (!identical(pin_stage, "all") && !kl_pin_cfg(kcfg)$enable) {
    stop("--pin_stage ", pin_stage, " needs occurrence.kernel.pin.enable.")
  }
  if (kl_pin_cfg(kcfg)$enable) {
    # kernel.pin: group lag distributions (pilots or a sweep) -> candidate group
    # orders -> one seed-pinned run per order -> stacking; the highest-weight
    # run is `fit`, every run and the weights go to the kl_pin artefact.
    # --pin_stage splits this for SLURM: plan (one job) -> run (array, one order
    # per task) -> combine (one job).
    bb <- .occ_baseline_block(basehaz, prep$ms_basis)
    pin_file <- make_run_file(opt, kind = "kl_pin")
    plan_file <- make_run_file(opt, kind = "kl_pin_plan")
    read_plan <- function() {
      if (!file.exists(plan_file)) {
        stop(
          "kernel.pin plan not found (run --pin_stage plan first): ",
          plan_file
        )
      }
      readRDS(plan_file)
    }
    if (identical(pin_stage, "plan")) {
      plan <- kl_pin_plan(
        opt,
        pois_dat,
        brms_cfg,
        bb,
        m,
        pg,
        seed = opt$seed,
        stan_dir = stan_dir,
        cores = as.integer(brms_cfg$cores %||% 4L)
      )
      saveRDS(plan, plan_file)
      message("Saved kernel.pin plan -> ", plan_file)
      cat(sprintf(
        "N_TASKS=%d\nPLAN_FILE=%s\n",
        nrow(plan$arrangements),
        plan_file
      ))
      quit(save = "no", status = 0)
    }
    if (identical(pin_stage, "run")) {
      plan <- read_plan()
      i <- opt$pin_task
      if (is.null(i) || is.na(i) || i < 1L) {
        stop("--pin_stage run needs --pin_task >= 1.")
      }
      if (i > nrow(plan$arrangements)) {
        message(sprintf(
          "kernel.pin task %d: the plan has %d order(s) -- nothing to do.",
          i,
          nrow(plan$arrangements)
        ))
        quit(save = "no", status = 0)
      }
      run_file <- kl_pin_run_file(pin_file, i)
      if (file.exists(run_file) && !isTRUE(opt$overwrite)) {
        message(
          "kernel.pin run already exists (use --overwrite to re-fit): ",
          run_file
        )
        quit(save = "no", status = 0)
      }
      message(sprintf(
        "kernel.pin task %d: order %s",
        i,
        plan$arrangements$arrangement[i]
      ))
      f <- kl_pin_fit_one(
        m,
        brms_cfg,
        pg,
        adj,
        plan$seeds,
        plan$arrangements,
        i,
        seed = opt$seed,
        stan_dir = stan_dir
      )
      attr(f, "opt") <- opt
      saveRDS(f, run_file)
      message("Saved kernel.pin run -> ", run_file)
      quit(save = "no", status = 0)
    }
    if (identical(pin_stage, "combine")) {
      plan <- read_plan()
      files <- kl_pin_run_file(pin_file, seq_len(nrow(plan$arrangements)))
      if (any(!file.exists(files))) {
        stop(
          "kernel.pin combine: missing run file(s):\n  ",
          paste(files[!file.exists(files)], collapse = "\n  ")
        )
      }
      pin <- kl_pin_finish(plan, lapply(files, readRDS))
    } else {
      pin <- kl_pin_run(
        opt,
        pois_dat,
        brms_cfg,
        bb,
        m,
        pg,
        adj,
        seed = opt$seed,
        stan_dir = stan_dir,
        cores = as.integer(brms_cfg$cores %||% 4L) *
          kl_pin_cfg(kcfg)$parallel_runs
      )
    }
    print(pin$table)
    fit <- pin$fits[[pin$best]]
    message("[kl_pin] reported fit: arrangement ", fit$arrangement)
  } else {
    # identify = seed: the seed countries of prior group k are fixed to class k
    seed_fix <- NULL
    seed_labels <- NULL
    if (identical(kl_identify(kcfg), "seed")) {
      reg <- load_grouping(kcfg$prior_grouping$source, "analysis/groupings")
      groups <- pg$class_labels[pg$class_labels %in% reg$group]
      kl_meta$seeds <- kl_pin_seeds(
        m$data,
        reg,
        groups,
        as.integer(kcfg$seeds$per_group %||% 1L),
        kcfg$seeds[["iso3"]]
      )
      seed_fix <- kl_pin_fix(ctry_levels, kl_meta$seeds, groups)
      seed_labels <- pg$class_labels
      message(
        "Seeds (fixed class): ",
        paste(
          sprintf(
            "%s=%s",
            groups,
            vapply(kl_meta$seeds, paste, "", collapse = "+")
          ),
          collapse = "; "
        )
      )
    }
    fit <- fit_kl_stan(
      m,
      brms_cfg,
      pg,
      adj,
      seed = opt$seed,
      stan_dir = stan_dir,
      fix = seed_fix,
      labels = seed_labels
    )
  }
  summ <- summarise_kl_latent(fit, m$meta, cen = kcfg$cen)
  kl_meta$classes <- fit$K
  kl_meta$class_labels <- fit$class_labels
  kl_meta$identify <- fit$identify
  kl_meta$form <- fit$form
  kl_meta$relabel <- fit$relabel
  kl_meta$naming <- fit$naming
  kl_meta$country_levels <- fit$country_levels
  kl_meta$prior_partition <- fit$prior_partition
  kl_meta$adjacency <- fit$adjacency
  kl_meta$program_sha1 <- fit$program_sha1
  kl_meta$program_file <- fit$program_file
  kl_meta$inits <- fit$inits
  kl_meta$init_strategy <- fit$init_strategy
  kl_meta$init_table <- fit$init_table
  kl_meta$chain_agreement <- kl_chain_agreement(fit)
  if (!is.null(pin)) {
    kl_meta$pin <- list(
      table = pin$table,
      seeds = pin$seeds,
      arrangements = pin$arrangements[, c("arrangement", "prob")],
      pilots = lapply(pin$pilots, function(p) {
        p[c("group", "countries", "lags", "lag_prob", "kcon")]
      }),
      reported = fit$arrangement
    )
    pin_file <- make_run_file(opt, kind = "kl_pin")
    saveRDS(pin, pin_file)
    message(
      "Saved pinned-arrangement runs (all fits, pilots, stacking table) -> ",
      pin_file
    )
  }
  kl_meta$loco_note <- paste(
    "kl_log_lik_ctry is the leave-one-country-out unit;",
    "country-block PSIS-LOO (kl_loo_country()) is descriptive only."
  )
}
message(sprintf(
  "Fit time: %.1f min",
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
))

# 7. Save ----------------------------------------------------------------------
# Output files (via make_run_file):
#   kind = "fit"        — kl_fit (engine stan) or brmsfit (engine brms)
#   kind = "kl_meta"    — kernel basis metadata (+ class/partition/adjacency info)
#   kind = "kl_summary" — kernel parameter / lag-weight / exposure-response tables
#   latent classes only:
#   kind = "kl_resp"         — class responsibilities per country
#   kind = "kl_comembership" — posterior co-membership matrix
#   kind = "kl_partition"    — MAP partition (rds; + .csv in the groupings format)

attr(fit, "opt") <- opt
attr(fit, "episode_stamp") <- ep_rec
saveRDS(fit, fit_file)
message("Saved fit -> ", fit_file)

meta_file <- make_run_file(opt, kind = "kl_meta")
saveRDS(kl_meta, meta_file)
message("Saved kernel meta -> ", meta_file)

summ_file <- make_run_file(opt, kind = "kl_summary")
saveRDS(summ, summ_file)
message("Saved kernel summary -> ", summ_file)

if (engine == "brms") {
  print(summ$params)
} else {
  print(summ$global)
  for (s in summ$classes) {
    print(s$params)
  }
}
message(
  "RE stats skipped for brms_kl: the brms adapter (compute_nested_surv_re_stats_brms()) ",
  "does not parse nlpar-prefixed / patched-program draws."
)
if (latent) {
  resp <- kl_responsibilities(fit)
  saveRDS(resp, make_run_file(opt, kind = "kl_resp"))
  saveRDS(kl_comembership(fit), make_run_file(opt, kind = "kl_comembership"))
  part <- kl_map_partition(fit, "argmax", resp = resp)
  part_file <- make_run_file(opt, kind = "kl_partition")
  saveRDS(part, part_file)
  write_kl_partition(part, sub("\\.rds$", ".csv", part_file))
  message(
    "Saved responsibilities / co-membership / MAP partition -> ",
    dirname(part_file)
  )
  print(as.data.frame(part))
  print(kl_meta$chain_agreement)
  if (!is.null(kcfg$partition_name) && nzchar(kcfg$partition_name)) {
    gfile <- file.path(
      "analysis/groupings",
      paste0(kcfg$partition_name, ".csv")
    )
    write_kl_partition(part, gfile)
    message("Wrote MAP partition to the groupings registry -> ", gfile)
  }
}
