# 05_fit_postprocess.R — unified postprocessing of phase-04 fits, all backends
#
# Usage:
#   Rscript analysis/05_fit_postprocess.R -C <config.yml>
#   Rscript analysis/05_fit_postprocess.R -d analysis/configs/<batch>/
#   Rscript analysis/05_fit_postprocess.R -d analysis/configs/<batch>/ --redo
#
# Details:
#   Takes one config (-C) or a batch directory (-d); for each config loads the
#   fit RDS and caches outputs under generated_data/ (skipped if present
#   unless --redo): *_re_stats.rds, *_coef_summary.rds (mean + 95% CI),
#   *_loo.rds and, per batch, *_compare.rds.

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  # Explicit, eager load: brms::ranef() on a reloaded (not freshly-fit)
  # brms_cb fit dispatches S4 methods (e.g. dimnames) on the "stanfit" class,
  # which triggers R's lazy auto-load of rstan via .requirePackage()/require().
  # That lazy path swallows an underlying rstan load failure (e.g. an
  # Rcpp/StanHeaders ABI mismatch) into a warning and leaves rstan half-loaded,
  # surfacing only as "unable to find required package 'rstan'" deep in the
  # S4 dispatch stack. Loading rstan cleanly upfront avoids that fragile path.
  library(rstan)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/clim_compile.R")
source("analysis/utils/clim_db.R")
source("analysis/utils/occ_data.R")
source("analysis/utils/siz_data.R")
source("analysis/utils/dur_data.R")
source("analysis/utils/population.R")
source("analysis/utils/composite_lp.R")
source("analysis/utils/provenance.R")
source("analysis/utils/occ_post.R")
source("analysis/utils/siz_post.R")
source("analysis/utils/dur_post.R")
source("analysis/utils/plots.R")
source("analysis/utils/anchor.R")

# 2. Options -------------------------------------------------------------------
# Script-level options (separate from per-run options inside each config)
script_options <- list(
  make_option(
    c("-d", "--batch_dir"),
    default = NULL,
    type = "character",
    help = "Directory of YAML config files to postprocess (batch mode)"
  ),
  make_option(
    c("--redo"),
    default = FALSE,
    action = "store_true",
    help = "Recompute and overwrite existing output files [default: FALSE]"
  ),
  make_option(
    c("--no_loo"),
    default = FALSE,
    action = "store_true",
    help = "Skip LOO computation (slow for large fits) [default: FALSE]"
  ),
  make_option(
    c("--no_compare"),
    default = FALSE,
    action = "store_true",
    help = "Skip cross-backend comparison table [default: FALSE]"
  ),
  make_option(
    c("--kfold_group"),
    default = FALSE,
    action = "store_true",
    help = paste(
      "Monthly model only: additionally compute",
      "brms::kfold(folds = 'grouped', group = 'episode_id')",
      "-- corrects row-wise LOO's within-episode",
      "correlation optimism [default: FALSE]"
    )
  ),
  make_option(
    c("--kfold_K"),
    default = 10L,
    type = "integer",
    help = "Number of folds for --kfold_group [default: 10]"
  )
)

# Inject --config / --update_config via make_options_from_config so we can
# accept a single-config path on -C as well.
script_opt <- make_options_from_config(option_list = script_options)
print_options(script_opt)

# Recover the config path: make_options_from_config intentionally drops the
# "config" key from opt (to avoid writing it into --update_config YAML dumps),
# so we re-extract it from the raw argv here.
local({
  argv <- commandArgs(trailingOnly = TRUE)
  for (i in seq_along(argv)) {
    if (argv[i] %in% c("-C", "--config") && i < length(argv)) {
      script_opt$config <<- argv[i + 1L]
      break
    }
  }
})

# 3. Collect config paths ------------------------------------------------------

config_paths <- character(0)

if (!is.null(script_opt$batch_dir)) {
  if (!dir.exists(script_opt$batch_dir)) {
    stop("Batch directory not found: ", script_opt$batch_dir)
  }
  config_paths <- list.files(
    script_opt$batch_dir,
    pattern = "\\.yml$",
    full.names = TRUE
  )
  config_paths <- config_paths[basename(config_paths) != "manifest.yml"]
  if (length(config_paths) == 0) {
    stop("No .yml files found in: ", script_opt$batch_dir)
  }
  message(sprintf(
    "Batch mode: %d config(s) in %s",
    length(config_paths),
    script_opt$batch_dir
  ))
} else if (!is.null(script_opt$config)) {
  config_paths <- script_opt$config
} else {
  stop("Provide either -C <config.yml> (single) or -d <dir/> (batch).")
}

# 4. Per-config processing -----------------------------------------------------

.extract_fixed_effects <- function(fit, backend) {
  switch(
    backend,
    brms = ,
    brms_cb = {
      # brms_cb (DLNM crossbasis) fits are still brms::brmsfit objects underneath
      # -- brms::fixef() dispatches on class, not on our config-level backend
      # string, so the crossbasis variant needs no separate branch. Same fix as
      # .compute_loo() below: before this, brms_cb fell through switch()'s
      # implicit no-match into stop("Unknown backend: brms_cb") here, which
      # tryCatch converted to a silent "Coef summary failed" warning -- no
      # coef_summary file was ever written for ANY brms_cb config (occurrence,
      # size, duration, or monthly), across every track, until this fix.
      #
      # brms: fixef() returns Estimate / Est.Error / Q2.5 / Q97.5
      fe <- brms::fixef(fit, probs = c(0.025, 0.975))
      tibble::tibble(
        term = rownames(fe),
        mean = fe[, "Estimate"],
        q025 = fe[, "Q2.5"],
        q975 = fe[, "Q97.5"]
      )
    },
    stop("Unknown backend: ", backend)
  )
}

.compute_loo <- function(fit, backend) {
  switch(
    backend,
    brms = ,
    brms_cb = {
      # brms_cb (DLNM crossbasis) fits are still brms::brmsfit objects underneath
      # -- brms::loo() dispatches on class, not on our config-level backend
      # string, so the crossbasis variant needs no separate branch. Before this
      # fix, brms_cb fell through switch()'s implicit no-match, which returns
      # NULL SILENTLY (no error, no warning): 05_fit_postprocess.R produced no
      # LOO file for any brms_cb config, and nothing in its output said so. The
      # gap was masked because 05_fit_postprocess_loo_lags.R computes and caches
      # LOO itself rather than calling this helper.
      tryCatch(brms::loo(fit), error = function(e) {
        warning("LOO failed for brms fit: ", conditionMessage(e))
        NULL
      })
    },
    stop("Unknown backend: ", backend)
  )
}

.re_adapter <- function(backend) {
  switch(
    backend,
    brms = compute_nested_surv_re_stats_brms,
    stop("Unknown backend: ", backend)
  )
}

process_one <- function(
  config_path,
  redo = FALSE,
  no_loo = FALSE,
  kfold_group = FALSE,
  kfold_K = 10L
) {
  message("\n=== ", config_path, " ===")

  opt <- yaml::read_yaml(config_path)
  model <- if (!is.null(opt$occurrence)) {
    "occurrence"
  } else if (!is.null(opt$size)) {
    "size"
  } else if (!is.null(opt$duration)) {
    "duration"
  } else if (!is.null(opt$monthly)) {
    "monthly"
  } else {
    stop(
      "config must have occurrence:, size:, duration: or monthly: block: ",
      config_path
    )
  }
  # The track's own default (monthly is brms_cb), as the fit scripts resolve it.
  backend <- opt[[model]]$backend %||% get_default_options()[[model]]$backend
  opt[[model]]$backend <- backend
  # brms_kl fits are kl_fit objects with their own summaries (written by
  # 04_fit_occ_brms_kl.R) and are compared across a sweep by
  # 05_fit_postprocess_kl_sweep.R; none of the adapters below applies.
  if (identical(backend, "brms_kl")) {
    message("  brms_kl: skipped (see 05_fit_postprocess_kl_sweep.R)")
    return(invisible(NULL))
  }

  # Fill top-level shared keys from defaults
  defaults <- get_default_options()
  for (nm in setdiff(names(defaults), c("occurrence", "size"))) {
    if (is.null(opt[[nm]])) opt[[nm]] <- defaults[[nm]]
  }
  # Fill model-block keys (one level deep) from defaults
  model_defaults <- defaults[[model]] %||% list()
  for (nm in names(model_defaults)) {
    if (is.null(opt[[model]][[nm]])) opt[[model]][[nm]] <- model_defaults[[nm]]
  }

  # Kind prefix: size model uses "size_fit", "size_re_stats", etc.; duration
  # uses "dur_fit", "dur_re_stats", etc.; monthly uses "mon_fit", "mon_re_stats",
  # etc. (routed by make_run_file() to make_duration_run_file() / make_monthly_run_file()
  # -- see utils/paths.R:655-666 and make_monthly_run_file()).
  kpfx <- if (model == "size") {
    "size_"
  } else if (model == "duration") {
    "dur_"
  } else if (model == "monthly") {
    "mon_"
  } else {
    ""
  }

  fit_file <- make_run_file(opt, kind = paste0(kpfx, "fit"))
  if (!file.exists(fit_file)) {
    warning("Fit file not found, skipping: ", fit_file)
    return(invisible(NULL))
  }

  message("Loading fit: ", fit_file)
  fit <- readRDS(fit_file)

  # Load matching data
  surv_dat_all <- tryCatch(
    {
      load_surv_data(opt)$dat # opt$country mutation discarded; make_run_file already called above
    },
    error = function(e) {
      warning("Could not load survival data: ", conditionMessage(e))
      NULL
    }
  )

  # Duration panel (one row per active month); only needed for the duration
  # track's RE-stats call below. Mirrors load_duration_data() as called at
  # 04_fit_dur_brms.R:135 -- plain load, no ms-spline basis reconstruction,
  # since RE-stats enumeration only needs the RE grouping columns and the
  # smoke config (basehaz_duration: factor) never touches the basis columns.
  dur_dat_all <- if (model == "duration") {
    tryCatch(
      {
        load_duration_data(opt)
      },
      error = function(e) {
        warning("Could not load duration data: ", conditionMessage(e))
        NULL
      }
    )
  } else {
    NULL
  }

  results <- list(config = config_path, backend = backend, opt = opt)

  # RE stats
  re_file <- make_run_file(opt, kind = paste0(kpfx, "re_stats"))
  if (!file.exists(re_file) || redo) {
    drop_region <- isTRUE(opt$grouping$enable)
    tryCatch(
      {
        re_stats <- if (model == "size") {
          if (identical(backend, "brms")) {
            compute_nested_size_re_stats_brms(
              fit,
              opt,
              drop_region = drop_region
            )
          } else {
            warning("RE stats: unknown size backend '", backend, "'; skipping")
            NULL
          }
        } else if (model == "monthly") {
          # Reused unchanged: the size adapter unpacks brms::ranef() purely from
          # the fit object (opt is accepted but not read -- see siz_post.R), so
          # it works for monthly's deeper (AFRO_region/country/ADM1/gid) nesting
          # exactly as it does for size's (AFRO_region/country/ADM1) -- the
          # colon-count-based label parser is nesting-depth-agnostic.
          if (identical(backend, "brms") || identical(backend, "brms_cb")) {
            compute_nested_size_re_stats_brms(
              fit,
              opt,
              drop_region = drop_region
            )
          } else {
            warning(
              "RE stats: unknown monthly backend '",
              backend,
              "'; skipping"
            )
            NULL
          }
        } else if (model == "duration") {
          # Reused, not forked: the duration cessation hazard is Poisson with
          # the same AFRO_region/country/ADM1 nesting as occurrence, so the
          # occurrence RE-stats adapter works unchanged.
          if (identical(backend, "brms") && !is.null(dur_dat_all)) {
            compute_nested_surv_re_stats_brms(
              fit,
              dur_dat_all,
              drop_region = drop_region
            )
          } else {
            warning(
              "RE stats: unknown duration backend '",
              backend,
              "' or missing duration data; skipping"
            )
            NULL
          }
        } else if (!is.null(surv_dat_all)) {
          if (identical(backend, "brms")) {
            compute_nested_surv_re_stats_brms(
              fit,
              surv_dat_all,
              drop_region = drop_region
            )
          } else {
            .re_adapter(backend)(fit, surv_dat_all)
          }
        } else {
          NULL
        }
        if (!is.null(re_stats)) {
          saveRDS(re_stats, re_file)
          message("Saved RE stats -> ", re_file)
        }
      },
      error = function(e) {
        warning("RE stats failed: ", conditionMessage(e))
      }
    )
  } else {
    message("RE stats already exist (--redo to recompute): ", re_file)
  }

  # Fixed-effect coefficient summary
  coef_file <- make_run_file(opt, kind = paste0(kpfx, "coef_summary"))
  if (!file.exists(coef_file) || redo) {
    tryCatch(
      {
        # Explicit per-track branch: each model's baseline-hazard knob lives at
        # a different schema location (occurrence nests per-backend; size and
        # duration are flat), so there is no shared fallback expression that
        # resolves correctly for all three. The fallbacks are the runtime
        # defaults (`defaults`, bound above), never literals: the flat size /
        # duration / monthly keys are already filled by the one-level loop above,
        # the per-backend occurrence keys are not (a backend block present
        # without the key).
        .basehaz_val <- if (model == "occurrence") {
          if (backend == "brms") {
            opt$occurrence[[backend]][["basehaz_poisson"]] %||%
              defaults$occurrence[[backend]][["basehaz_poisson"]]
          } else {
            NA_character_
          }
        } else if (model == "size") {
          opt$size[["seasonal_term"]] %||% defaults$size[["seasonal_term"]]
        } else if (model == "duration") {
          opt$duration[["basehaz_duration"]] %||%
            defaults$duration[["basehaz_duration"]]
        } else if (model == "monthly") {
          opt$monthly[["seasonal_term"]] %||%
            defaults$monthly[["seasonal_term"]]
        } else {
          NA_character_
        }
        coef_tbl <- .extract_fixed_effects(fit, backend) %>%
          dplyr::mutate(
            backend = backend,
            baseline = opt$baseline,
            country = opt$country,
            basehaz = .basehaz_val
          )
        saveRDS(coef_tbl, coef_file)
        message("Saved coef summary -> ", coef_file)
        results$coef_summary <- coef_tbl
      },
      error = function(e) {
        warning("Coef summary failed: ", conditionMessage(e))
      }
    )
  } else {
    message("Coef summary already exists (--redo to recompute): ", coef_file)
    results$coef_summary <- readRDS(coef_file)
  }

  # LOO
  if (!no_loo) {
    if (model == "monthly") {
      message(
        "Note: row-wise LOO treats each episode-month as independent; ",
        "within-episode correlation (repeated months of the same ",
        "outbreak) can make it optimistic. Use --kfold_group for a ",
        "leave-one-episode-out alternative."
      )
    }
    loo_file <- make_run_file(opt, kind = paste0(kpfx, "loo"))
    if (!file.exists(loo_file) || redo) {
      tryCatch(
        {
          loo_obj <- .compute_loo(fit, backend)
          if (!is.null(loo_obj)) {
            saveRDS(loo_obj, loo_file)
            message("Saved LOO -> ", loo_file)
          }
        },
        error = function(e) {
          warning("LOO failed: ", conditionMessage(e))
        }
      )
    } else {
      message("LOO already exists (--redo to recompute): ", loo_file)
    }
  }

  # K-fold (monthly model only, opt-in via --kfold_group) --------------------
  # Leave-one-episode-out alternative to row-wise LOO: groups folds by
  # episode_id so repeated within-episode months never split across the
  # train/test boundary, correcting the correlation-optimism noted above.
  if (isTRUE(kfold_group) && model == "monthly") {
    kfold_file <- make_run_file(opt, kind = "mon_kfold")
    if (!file.exists(kfold_file) || redo) {
      tryCatch(
        {
          if (identical(backend, "brms") || identical(backend, "brms_cb")) {
            kfold_obj <- brms::kfold(
              fit,
              folds = "grouped",
              group = "episode_id",
              K = kfold_K
            )
            saveRDS(kfold_obj, kfold_file)
            message(
              "Saved k-fold (grouped by episode_id, K = ",
              kfold_K,
              ") -> ",
              kfold_file
            )
          } else {
            warning("k-fold: unknown monthly backend '", backend, "'; skipping")
          }
        },
        error = function(e) {
          warning("k-fold failed: ", conditionMessage(e))
        }
      )
    } else {
      message("k-fold already exists (--redo to recompute): ", kfold_file)
    }
  }

  # PP check (duration model) ------------------------------------------------
  # Own builder rather than the size scatter: the duration response is a 0/1
  # cessation indicator per active month, so the checkable quantities are the
  # hazard-by-mso calibration and the IMPLIED duration distribution (see
  # utils/dur_post.R header). Needs the panel filtered to the fit's own rows,
  # which under grouping mode means resolving the group's countries into
  # opt$country first -- dur_dat_all above is deliberately the full panel
  # (the RE-stats adapter only enumerates RE labels and does not care).
  if (model == "duration") {
    pp_check_file <- make_run_file(opt, kind = "dur_pp_check")
    if (!file.exists(pp_check_file) || redo) {
      tryCatch(
        {
          if (!identical(backend, "brms") && !identical(backend, "brms_cb")) {
            warning(
              "pp_check: unknown duration backend '",
              backend,
              "'; skipping"
            )
          } else {
            dur_dat_fit <- load_duration_data(resolve_group_countries(opt))
            pp_check_obj <- make_duration_pp_check_tbl(fit, dur_dat_fit)
            saveRDS(pp_check_obj, pp_check_file)
            message("Saved pp_check -> ", pp_check_file)
          }
        },
        error = function(e) {
          warning("pp_check failed: ", conditionMessage(e))
        }
      )
    } else {
      message("pp_check already exists (--redo to recompute): ", pp_check_file)
    }
  }

  # PP check (size model only) -----------------------------------------------
  if (model == "size") {
    pp_check_file <- make_size_run_file(opt, kind = "size_pp_check")
    if (!file.exists(pp_check_file) || redo) {
      tryCatch(
        {
          if (identical(backend, "brms")) {
            pp_check_obj <- make_size_pp_check_tbl(fit)
          } else {
            warning("pp_check: unknown size backend '", backend, "'; skipping")
            pp_check_obj <- NULL
          }
          if (!is.null(pp_check_obj)) {
            saveRDS(pp_check_obj, pp_check_file)
            message("Saved pp_check -> ", pp_check_file)
          }
        },
        error = function(e) {
          warning("pp_check failed: ", conditionMessage(e))
        }
      )
    } else {
      message("pp_check already exists (--redo to recompute): ", pp_check_file)
    }

    # mean_rate table QC (S8): CSV export of the cell table + a tier summary,
    # emitted once per (min_pop, max_attack_rate, estimator) triple -- the same
    # triple every size config with matching values resolves to via
    # make_size_run_file(opt, "size_mean_rate"), so this is written next to
    # that sidecar rather than per-config.
    mrt_rds <- make_size_run_file(opt, kind = "size_mean_rate")
    mrt_csv <- sub("\\.rds$", ".csv", mrt_rds)
    tier_csv <- sub("\\.rds$", "_tier_summary.csv", mrt_rds)
    if (!file.exists(mrt_rds)) {
      message(
        "mean_rate table sidecar not found (written by persist_mean_rate_table() ",
        "at fit time): ",
        mrt_rds
      )
    } else if (!file.exists(mrt_csv) || redo) {
      tryCatch(
        {
          mrt <- readRDS(mrt_rds)
          utils::write.csv(mrt$cell, mrt_csv, row.names = FALSE)
          tier_summary <- mrt$cell %>%
            dplyr::count(mr_tier, name = "n_cells") %>%
            dplyr::mutate(pct = round(100 * n_cells / sum(n_cells), 1))
          utils::write.csv(tier_summary, tier_csv, row.names = FALSE)
          message("Saved mean_rate QC -> ", mrt_csv, " + ", tier_csv)
        },
        error = function(e) {
          warning("mean_rate QC export failed: ", conditionMessage(e))
        }
      )
    } else {
      message("mean_rate QC already exists (--redo to recompute): ", mrt_csv)
    }
  }

  # PP check (monthly model only) --------------------------------------------
  # Reuses make_size_pp_check_tbl() with response = "cases_month" (instead of
  # the size default "suspected_cases") -- no
  # mean_rate sidecar here, monthly has no analogous QC export.
  if (model == "monthly") {
    pp_check_file <- make_run_file(opt, kind = "mon_pp_check")
    if (!file.exists(pp_check_file) || redo) {
      tryCatch(
        {
          if (identical(backend, "brms") || identical(backend, "brms_cb")) {
            pp_check_obj <- make_size_pp_check_tbl(
              fit,
              response = "cases_month"
            )
          } else {
            warning(
              "pp_check: unknown monthly backend '",
              backend,
              "'; skipping"
            )
            pp_check_obj <- NULL
          }
          if (!is.null(pp_check_obj)) {
            saveRDS(pp_check_obj, pp_check_file)
            message("Saved pp_check -> ", pp_check_file)
          }
        },
        error = function(e) {
          warning("pp_check failed: ", conditionMessage(e))
        }
      )
    } else {
      message("pp_check already exists (--redo to recompute): ", pp_check_file)
    }
  }

  rm(fit)
  gc()

  invisible(results)
}

# 5. Process all configs -------------------------------------------------------

all_results <- lapply(
  config_paths,
  process_one,
  redo = isTRUE(script_opt$redo),
  no_loo = isTRUE(script_opt$no_loo),
  kfold_group = isTRUE(script_opt$kfold_group),
  kfold_K = script_opt$kfold_K %||% 10L
)

# 6. Cross-backend comparison table --------------------------------------------
# *_compare.rds: side-by-side coef table across all configs in the batch,
# produced only when the batch contains >= 2 backends or basehaz variants
# with matching baseline × country.
# Produced when batch contains >= 2 configs, or when --no_compare is not set.

if (!isTRUE(script_opt$no_compare) && length(all_results) >= 2) {
  coef_list <- purrr::keep(
    purrr::map(all_results, "coef_summary"),
    ~ !is.null(.)
  )

  if (length(coef_list) >= 2) {
    compare_tbl <- dplyr::bind_rows(coef_list) %>%
      dplyr::arrange(term, backend)

    # Determine output path from the first config's opt
    first_opt <- all_results[[1]]$opt
    compare_file <- make_run_file(
      modifyList(first_opt, list(backend = "compare")),
      kind = "compare"
    )
    # Batch mode: a per-batch name under <outdir>/loo_compare/<batch>/
    # (make_loo_compare_run_file(), utils/paths.R) -- the batch config dir itself
    # holds inputs only.
    if (!is.null(script_opt$batch_dir)) {
      compare_file <- make_loo_compare_run_file(
        first_opt,
        "loo_compare_coef",
        batch = script_opt$batch_dir
      )
      dir.create(dirname(compare_file), recursive = TRUE, showWarnings = FALSE)
    }
    saveRDS(compare_tbl, compare_file)
    message("\nCross-backend comparison table -> ", compare_file)
    message(sprintf(
      "  %d terms × %d backends",
      length(unique(compare_tbl$term)),
      length(unique(compare_tbl$backend))
    ))

    # Print a readable summary to console
    cat("\n--- Fixed-effect comparison (posterior mean [95% CI]) ---\n")
    compare_tbl %>%
      dplyr::select(term, backend, basehaz, country, mean, q025, q975) %>%
      dplyr::mutate(ci = sprintf("%.3f [%.3f, %.3f]", mean, q025, q975)) %>%
      dplyr::select(term, backend, ci) %>%
      tidyr::pivot_wider(names_from = backend, values_from = ci) %>%
      print(n = Inf)
  }
}

message("\nDone.")
