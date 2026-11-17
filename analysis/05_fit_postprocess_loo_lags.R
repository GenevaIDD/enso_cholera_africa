# 05_fit_postprocess_loo_lags.R — LOO/elpd comparison of single-lag ENSO fits
#
# Usage:
#   Rscript analysis/05_fit_postprocess_loo_lags.R \
#     -d analysis/configs/enso_cb_singlelag_occ/
#   Rscript analysis/05_fit_postprocess_loo_lags.R -d <occ_batch>/ \
#     --group enso_west --copy_best generated_data/enso_best_lag_fits/
#   Rscript analysis/05_fit_postprocess_loo_lags.R -d <size_batch>/ \
#     --group enso_west --copy_best <dir> --min_lag 8
#   Rscript analysis/05_fit_postprocess_loo_lags.R -d <batch>/ \
#     --group enso_west --lag_sidecar --sidecar_only
#
# Details:
#   Each config in the batch dir (-d, inputs only) is one ENSO lag; the track
#   is auto-detected per config. Loads each fit, caches brms::loo, and ranks
#   the lags by elpd (loo::loo_compare). Writes under
#   <outdir>/loo_compare/<batch>/ (make_loo_compare_run_file()) the ranking
#   rds + png and, with --lag_sidecar, the per-lag LOO sidecar; --copy_best
#   copies the winning fit's artefacts.
#
# See: analysis/utils/loo_lags.R (best-lag rule, ranking readers)

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(ggplot2)
  library(loo)
})

source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")
source("analysis/utils/loo_lags.R")

# 2. Options -------------------------------------------------------------------
script_options <- list(
  make_option(
    c("-d", "--batch_dir"),
    default = NULL,
    type = "character",
    help = "Directory of single-lag YAML configs to compare (required)"
  ),
  make_option(
    c("--group"),
    default = NULL,
    type = "character",
    help = "Filter configs to this grouping group (e.g. enso_east) [default: all]"
  ),
  make_option(
    c("--redo"),
    default = FALSE,
    action = "store_true",
    help = "Recompute and overwrite cached per-fit LOO [default: FALSE]"
  ),
  make_option(
    c("--copy_best"),
    default = NULL,
    type = "character",
    help = paste(
      "If set, copy the winning fit's artefacts",
      "(fit/re_stats/loo/cb_meta) into this directory"
    )
  ),
  make_option(
    c("--min_lag"),
    default = NULL,
    type = "integer",
    help = paste(
      "Occurrence lag for this group; required by --copy_best on",
      "size/duration/monthly tracks to enforce track-lag >= occ-lag"
    )
  ),
  make_option(
    c("--lag_rule"),
    default = "longest_within_se",
    type = "character",
    help = "Best-lag rule: longest_within_se | best_elpd (utils/loo_lags.R; deploy.lag_selection.rule) [default: %default]"
  ),
  make_option(
    c("--occ_floor_offset"),
    default = NA,
    type = "double",
    help = paste(
      "best_elpd only: a non-occurrence pick must be >= --min_lag minus this many",
      "months (NA = no occurrence floor) [default: %default]"
    )
  ),
  make_option(
    c("--se_mult"),
    default = 2,
    type = "double",
    help = paste(
      "SE multiplier for the occurrence best-lag rule (max single_lag",
      "within --se_mult SEs of elpd_diff from the rank-1 fit) [default: 2]"
    )
  ),
  make_option(
    c("--lag_sidecar"),
    default = FALSE,
    action = "store_true",
    help = paste(
      "Also write lag_loo_sidecar[_<group>]_<YYYYMMDD>.rds: N x L pointwise elpd",
      "and Pareto-k, per-lag b_enso_cb* draws + cb_meta, row-key digests.",
      "Loads every fit even when its LOO is cached [default: FALSE]"
    )
  ),
  make_option(
    c("--sidecar_only"),
    default = FALSE,
    action = "store_true",
    help = paste(
      "Write only the sidecar (implies --lag_sidecar): the ranking rds and",
      "png are not written [default: FALSE]"
    )
  )
)
opt_parser <- OptionParser(option_list = script_options)
script_opt <- parse_args(opt_parser)

if (is.null(script_opt$batch_dir) || !dir.exists(script_opt$batch_dir)) {
  stop("Provide a valid batch directory via -d <dir/>.")
}
if (isTRUE(script_opt$sidecar_only)) {
  script_opt$lag_sidecar <- TRUE
  if (!is.null(script_opt$copy_best)) {
    stop(
      "--sidecar_only writes the sidecar alone; drop --copy_best (or --sidecar_only)."
    )
  }
}
WRITE_SIDECAR <- isTRUE(script_opt$lag_sidecar)
WRITE_RANKING <- !isTRUE(script_opt$sidecar_only)

config_paths <- list.files(
  script_opt$batch_dir,
  pattern = "\\.yml$",
  full.names = TRUE
)
config_paths <- config_paths[basename(config_paths) != "manifest.yml"]
if (length(config_paths) == 0) {
  stop("No .yml configs found in: ", script_opt$batch_dir)
}

if (!is.null(script_opt$group)) {
  # Config filenames strip non-alphanumeric chars (make_batch_configs.R line ~1486)
  grp_fname <- gsub("[^A-Za-z0-9]", "", script_opt$group)
  keep <- grepl(paste0("_grp-", grp_fname, "\\.yml$"), basename(config_paths))
  config_paths <- config_paths[keep]
  if (length(config_paths) == 0) {
    stop(
      "No configs match group '",
      script_opt$group,
      "' (filename token: '",
      grp_fname,
      "') in: ",
      script_opt$batch_dir
    )
  }
  message(sprintf(
    "Group filter '%s': %d config(s) selected",
    script_opt$group,
    length(config_paths)
  ))
}

message(sprintf(
  "Comparing %d single-lag fit(s) in %s",
  length(config_paths),
  script_opt$batch_dir
))

# 3. Helpers -------------------------------------------------------------------

# Fill an opt list read from YAML with config_defaults (mirrors 05_fit_postprocess.R).
.fill_defaults <- function(opt, model) {
  defaults <- get_default_options()
  for (nm in setdiff(
    names(defaults),
    c("occurrence", "size", "duration", "monthly")
  )) {
    if (is.null(opt[[nm]])) opt[[nm]] <- defaults[[nm]]
  }
  model_defaults <- defaults[[model]] %||% list()
  for (nm in names(model_defaults)) {
    if (is.null(opt[[model]][[nm]])) opt[[model]][[nm]] <- model_defaults[[nm]]
  }
  opt
}

# Resolve the single_lag value from a model block (defaults, then per-variable).
.resolve_single_lag <- function(model_blk) {
  dlnm_cfg <- model_blk$dlnm %||% list()
  sl <- dlnm_cfg$defaults$single_lag
  if (is.null(sl)) {
    for (vv in names(dlnm_cfg$variables %||% list())) {
      if (!is.null(dlnm_cfg$variables[[vv]]$single_lag)) {
        sl <- dlnm_cfg$variables[[vv]]$single_lag
        break
      }
    }
  }
  if (is.null(sl)) NA_integer_ else as.integer(sl)
}

# Sidecar entry of one lag (utils/loo_lags.R::lag_sidecar_entry()): the
# population ENSO crossbasis draws and the cb_meta "enso" entry of the fit,
# its row-key digest, and the loo's pointwise elpd / Pareto-k. NULL (with a
# warning) when the fit is missing or has no ENSO crossbasis.
.sidecar_entry_for <- function(
  fit,
  loo_obj,
  opt,
  kpfx,
  sl,
  label,
  fit_file,
  loo_file,
  cb_var = "enso"
) {
  meta_file <- make_run_file(opt, kind = paste0(kpfx, "cb_meta"))
  vm <- if (file.exists(meta_file)) readRDS(meta_file)[[cb_var]] else NULL
  cb_names <- vm$cb_col_names
  if (is.null(cb_names)) {
    cb_names <- sub(
      "^b_",
      "",
      grep(cb_coef_pattern(cb_var), posterior::variables(fit), value = TRUE)
    )
    warning(sprintf(
      "No cb_meta '%s' entry for lag %d (%s); population draws from the fit's b_%s_cb* only.",
      cb_var,
      sl,
      basename(meta_file),
      cb_var
    ))
  }
  if (!length(cb_names)) {
    warning(
      "No ",
      cb_var,
      " crossbasis coefficients in ",
      basename(fit_file),
      "; no sidecar entry."
    )
    return(NULL)
  }
  lag_sidecar_entry(
    single_lag = sl,
    label = label,
    loo = loo_obj,
    fit_data = fit$data,
    b_draws = fit_cb_pop_draws(fit, cb_names),
    cb_meta_enso = vm,
    exclude_cols = cb_names,
    files = c(
      fit = fit_file,
      loo = loo_file,
      cb_meta = if (file.exists(meta_file)) meta_file else NA_character_
    )
  )
}

# Per-config: load fit, compute+cache LOO, return list(single_lag, loo, label).
# With sidecar = TRUE the fit is loaded even when its LOO is cached, and the
# result carries the lag's sidecar entry (NULL when it could not be built).
.loo_for_config <- function(config_path, redo = FALSE, sidecar = FALSE) {
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
    warning("No occurrence:/size:/duration:/monthly: block: ", config_path)
    return(NULL)
  }
  opt <- .fill_defaults(opt, model)

  sl <- .resolve_single_lag(opt[[model]])
  if (is.na(sl)) {
    warning("No single_lag found; skipping ", config_path)
    return(NULL)
  }

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
  loo_file <- make_run_file(opt, kind = paste0(kpfx, "loo"))

  fit <- NULL
  if (file.exists(loo_file) && !redo) {
    message("Using cached LOO: ", loo_file)
    loo_obj <- readRDS(loo_file)
  } else {
    if (!file.exists(fit_file)) {
      warning("Fit not found, skipping: ", fit_file)
      return(NULL)
    }
    message("Loading fit: ", fit_file)
    fit <- readRDS(fit_file)
    loo_obj <- tryCatch(brms::loo(fit), error = function(e) {
      warning("LOO failed for ", basename(fit_file), ": ", conditionMessage(e))
      NULL
    })
    if (is.null(loo_obj)) {
      return(NULL)
    }
    saveRDS(loo_obj, loo_file)
    message("Saved LOO -> ", loo_file)
  }

  label <- if (sl < 0L) paste0("lag", sl_token(sl)) else sprintf("lag%02d", sl)
  side <- NULL
  if (isTRUE(sidecar)) {
    if (is.null(fit) && file.exists(fit_file)) {
      message("Loading fit (sidecar): ", fit_file)
      fit <- readRDS(fit_file)
    }
    side <- if (is.null(fit)) {
      warning("Fit not found for the sidecar: ", fit_file)
      NULL
    } else {
      .sidecar_entry_for(fit, loo_obj, opt, kpfx, sl, label, fit_file, loo_file)
    }
  }
  rm(fit)
  invisible(gc(verbose = FALSE))

  list(
    single_lag = sl,
    model = model,
    loo = loo_obj,
    label = label,
    opt = opt,
    kpfx = kpfx,
    sidecar = side
  )
}

# 4. Collect LOO per lag -------------------------------------------------------

results <- lapply(
  config_paths,
  .loo_for_config,
  redo = isTRUE(script_opt$redo),
  sidecar = WRITE_SIDECAR
)
results <- Filter(Negate(is.null), results)
if (length(results) < 2L) {
  stop("Need >= 2 successful LOO fits to compare; got ", length(results))
}

# Order by lag for stable labelling
results <- results[order(vapply(results, function(r) r$single_lag, integer(1)))]
loo_list <- lapply(results, function(r) r$loo)
names(loo_list) <- vapply(results, function(r) r$label, character(1))
model <- results[[1]]$model
cb_label <- "ENSO"

# 5. loo_compare (rows sorted best→worst by elpd_diff) -------------------------

cmp <- loo::loo_compare(loo_list) # cols incl. elpd_diff, se_diff
cmp_df <- as.data.frame(cmp)
# Label column: loo >= 2.8 returns a data.frame carrying an explicit `model`
# column with integer rownames ("1".."N"); older loo returns a matrix whose
# rownames ARE the model labels. Reading rownames unconditionally silently
# yields "1".."N" on the newer shape, which then fails the join below and
# leaves elpd_diff / se_diff as NA.
cmp_df$label <- if ("model" %in% names(cmp_df)) {
  as.character(cmp_df$model)
} else {
  rownames(cmp_df)
}

# Per-fit elpd_loo + se + p_loo from each loo object's estimates table
est <- lapply(results, function(r) {
  e <- r$loo$estimates
  data.frame(
    label = r$label,
    single_lag = r$single_lag,
    elpd_loo = e["elpd_loo", "Estimate"],
    se_elpd = e["elpd_loo", "SE"],
    p_loo = e["p_loo", "Estimate"],
    stringsAsFactors = FALSE
  )
})
est_df <- dplyr::bind_rows(est)

ranking <- est_df %>%
  dplyr::left_join(
    dplyr::select(cmp_df, label, elpd_diff, se_diff),
    by = "label"
  ) %>%
  dplyr::arrange(dplyr::desc(elpd_loo)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::select(
    rank,
    single_lag,
    label,
    elpd_loo,
    se_elpd,
    elpd_diff,
    se_diff,
    p_loo
  )

# Pareto-k summary of each lag's PSIS-LOO, appended (utils/loo_lags.R::
# ranking_pareto(); additive -- no selection rule reads these columns).
ranking <- ranking_pareto(
  ranking,
  loo_files = stats::setNames(
    lapply(results, `[[`, "loo"),
    vapply(results, function(r) as.character(r$single_lag), character(1))
  )
)

message("\n", strrep("-", 60))
message(sprintf(
  "Single-lag %s comparison (%s track) — best lag: %d",
  cb_label,
  model,
  ranking$single_lag[1]
))
message(strrep("-", 60))
print(ranking, row.names = FALSE)

# 6. Write outputs -------------------------------------------------------------
# Outputs (under <outdir>/loo_compare/<batch>/ — make_loo_compare_run_file(),
# utils/paths.R; the batch config dir -d holds INPUTS only. Readers locate the
# latest ranking via utils/loo_lags.R::loo_compare_ranking_files(), which also
# falls back to the legacy in-config-dir location):
#   loo_compare_lags[_<group>]_<YYYYMMDD>.rds
#                                    — ranking tibble (best→worst), with elpd_loo,
#                                      se_elpd, elpd_diff, se_diff, p_loo, single_lag,
#                                      and the Pareto-k summary
#                                      columns of each lag's PSIS-LOO
#                                      (utils/loo_lags.R::PARETO_K_COLS; additive,
#                                      no selection rule reads them)
#   loo_compare_lags[_<group>].png   — elpd-by-lag plot (point + SE), best lag flagged
#   lag_loo_sidecar[_<group>]_<YYYYMMDD>.rds
#                                    — (opt-in, --lag_sidecar) N x L pointwise elpd and
#                                      Pareto-k matrices, per-lag population ENSO
#                                      crossbasis draws (b_enso_cb*) + the cb_meta
#                                      "enso" entry, per-lag row-key digests, track
#                                      and group (utils/loo_lags.R::build_lag_loo_sidecar()).
#                                      Loads every fit even when its LOO is cached.
#                                      --sidecar_only writes only this file: the
#                                      ranking rds and png are left untouched, so a
#                                      retrofit does not move 09's enso_relation
#                                      freshness inputs.
#   per-fit LOO cached via make_run_file(opt, "loo" | "size_loo")

# Outputs go under <outdir>/loo_compare/<batch>/ (make_loo_compare_run_file(),
# utils/paths.R) -- NOT into the batch config dir, which holds inputs only.
# outdir comes from the batch's own configs (all share one); readers resolve
# the same location via utils/loo_lags.R::loo_compare_ranking_files().
w_opt <- results[[1]]$opt %||% list()
out_dir <- make_loo_compare_run_file(
  w_opt,
  "loo_compare_dir",
  batch = script_opt$batch_dir
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Lag sidecar (--lag_sidecar): one file per (batch, group); every lag of the
# ranking must have its entry, and the rows must line up across lags
# (build_lag_loo_sidecar() stops otherwise). Written after the ranking / png /
# copy_best (so a sidecar failure never blocks them), or alone under
# --sidecar_only.
.write_lag_sidecar <- function() {
  missing_side <- vapply(results, function(r) is.null(r$sidecar), logical(1))
  if (any(missing_side)) {
    stop(
      "--lag_sidecar: no sidecar entry for ",
      paste(vapply(results[missing_side], `[[`, "", "label"), collapse = ", "),
      " (fit or ENSO crossbasis missing; see the warnings above). Sidecar not written."
    )
  }
  sidecar <- build_lag_loo_sidecar(
    lapply(results, `[[`, "sidecar"),
    track = model,
    group = script_opt$group %||% NA_character_,
    batch = script_opt$batch_dir
  )
  out_side <- make_loo_compare_run_file(
    w_opt,
    "loo_compare_lag_sidecar",
    batch = script_opt$batch_dir,
    group = script_opt$group
  )
  saveRDS(sidecar, out_side)
  message(sprintf(
    "\nSaved lag sidecar (%d lags x %d rows, %d draws) -> %s",
    length(sidecar$lags),
    sidecar$n_obs,
    nrow(sidecar$cb[[1]]$b_draws),
    out_side
  ))
  invisible(out_side)
}

if (!WRITE_RANKING) {
  .write_lag_sidecar()
  message("--sidecar_only: ranking rds and png not written.")
  quit(save = "no", status = 0)
}

out_rds <- make_loo_compare_run_file(
  w_opt,
  "loo_compare_ranking",
  batch = script_opt$batch_dir,
  group = script_opt$group
)
saveRDS(list(ranking = ranking, loo_compare = cmp, model = model), out_rds)
message("\nSaved ranking -> ", out_rds)

best_lag <- ranking$single_lag[1]
plot_df <- ranking %>% dplyr::mutate(is_best = single_lag == best_lag)

grp_label <- if (!is.null(script_opt$group)) {
  paste0(" — ", script_opt$group)
} else {
  ""
}
p <- ggplot(plot_df, aes(x = single_lag, y = elpd_loo)) +
  geom_errorbar(
    aes(ymin = elpd_loo - se_elpd, ymax = elpd_loo + se_elpd),
    width = 0.2,
    colour = "grey50"
  ) +
  geom_line(colour = "grey60") +
  geom_point(aes(colour = is_best), size = 3) +
  scale_colour_manual(
    values = c(`FALSE` = "black", `TRUE` = "firebrick"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = sort(unique(plot_df$single_lag))) +
  labs(
    title = sprintf(
      "%s single-lag elpd comparison (%s%s)",
      cb_label,
      model,
      grp_label
    ),
    subtitle = sprintf("Best lag = %d (highest elpd_loo)", best_lag),
    x = paste(
      cb_label,
      if (any(plot_df$single_lag < 0L)) {
        "lag (months; < 0 = lead)"
      } else {
        "lag (months)"
      }
    ),
    y = "elpd_loo (\u00b1 SE)"
  ) +
  theme_bw()

out_png <- make_loo_compare_run_file(
  w_opt,
  "loo_compare_png",
  batch = script_opt$batch_dir,
  group = script_opt$group
)
ggsave(out_png, p, width = 7, height = 4.5, dpi = 150)
message("Saved elpd-by-lag plot -> ", out_png)

# 7. Copy best-fit artefacts ---------------------------------------------------
# --copy_best <dir> copies the winning fit's artefacts (fit / re_stats / loo /
# cb_meta) into <dir>, printing each copied file.
# Best-lag selection (shared with bash/tools/collect_enso_best_lag_fits.sh and
# the loo_lag_results_by_group.Rmd diagnostic report — see utils/loo_lags.R):
#   occurrence        — the MAXIMUM single_lag among fits within --se_mult SEs of
#                        elpd_diff from the rank-1 (best-elpd) fit.
#   size/duration/monthly — the best (lowest-rank) single_lag >= --min_lag, which
#                        callers pass as the occurrence lag for the same group.
#                        Required for --copy_best on these tracks.
#
# Persist the winning fit's artefacts (fit / re_stats / loo / cb_meta) into a
# best-fit folder, per the shared (hybrid) selection rule in
# analysis/utils/loo_lags.R: occurrence — max single_lag within --se_mult SEs
# of elpd_diff from rank-1; other tracks — the same unconstrained SE-window
# pick if it already clears --min_lag (the occurrence lag), else the
# SE-window pick within the single_lag >= --min_lag eligible subset. The
# ranking/plot above are unchanged — they always reflect the raw LOO order;
# only the *copied* fit is subject to the selection rule.

if (!is.null(script_opt$copy_best)) {
  dest <- script_opt$copy_best
  if (!dir.exists(dest)) {
    dir.create(dest, recursive = TRUE)
  }

  ls_rule <- lag_selection_cfg(list(
    rule = script_opt$lag_rule,
    se_mult = script_opt$se_mult,
    occ_floor_offset = if (is.na(script_opt$occ_floor_offset)) {
      "none"
    } else {
      script_opt$occ_floor_offset
    }
  ))
  if (identical(ls_rule$rule, "best_elpd")) {
    occ_l <- if (model == "occurrence") NULL else script_opt$min_lag
    sel <- select_best_lag_best_elpd(
      ranking,
      min_lag = ls_rule$min_lag,
      occ_lag = occ_l,
      occ_floor_offset = ls_rule$occ_floor_offset
    )
    rule_desc <- sprintf(
      "best_elpd: best-ranked lag%s",
      if (is.null(occ_l) || is.null(ls_rule$occ_floor_offset)) {
        ""
      } else {
        sprintf(" >= occ_lag (%d) - %g", occ_l, ls_rule$occ_floor_offset)
      }
    )
  } else if (model == "occurrence") {
    sel <- select_best_lag_occurrence(ranking, se_mult = script_opt$se_mult)
    rule_desc <- sprintf(
      "max single_lag within %.1f SE elpd_diff of rank-1 fit",
      script_opt$se_mult
    )
  } else {
    if (is.null(script_opt$min_lag)) {
      stop(
        "--copy_best on a non-occurrence track requires --min_lag <occ_lag> ",
        "(the occurrence lag selected for this group)."
      )
    }
    sel <- select_best_lag_downstream(
      ranking,
      occ_lag = script_opt$min_lag,
      se_mult = script_opt$se_mult
    )
    rule_desc <- sprintf(
      "max single_lag within %.1f SE of rank-1 if >= occ_lag (%d), else within %.1f SE of best eligible fit",
      script_opt$se_mult,
      script_opt$min_lag,
      script_opt$se_mult
    )
  }
  copy_lag <- sel$lag
  win <- Find(function(r) r$single_lag == copy_lag, results)

  message("\n", strrep("-", 60))
  message(sprintf(
    "Copying best fit (lag %d, rank %d) -> %s",
    copy_lag,
    sel$rank,
    dest
  ))
  message(sprintf("Rule: %s", rule_desc))
  message(strrep("-", 60))

  kinds <- paste0(win$kpfx, c("fit", "re_stats", "loo", "cb_meta"))
  copied <- character(0)
  for (k in kinds) {
    src <- make_run_file(win$opt, kind = k)
    if (!file.exists(src)) {
      message(sprintf("  [skip] %-14s (not found): %s", k, basename(src)))
      next
    }
    tgt <- file.path(dest, basename(src))
    if (isTRUE(file.copy(src, tgt, overwrite = TRUE, copy.date = TRUE))) {
      copied <- c(copied, tgt)
      message(sprintf("  [ok]   %-14s -> %s", k, tgt))
    } else {
      message(sprintf("  [FAIL] %-14s -> %s", k, tgt))
    }
  }
  message(sprintf("Copied %d file(s) to %s", length(copied), dest))
}

if (WRITE_SIDECAR) {
  .write_lag_sidecar()
}
