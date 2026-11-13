# 06_fit_figures_kl_sweep.R — figures of a brms_kl latent-class sweep
#
# Batch-level figures from the tables 05_fit_postprocess_kl_sweep.R wrote
# (never re-derived here): country-block elpd difference by K / prior / basis,
# stacking weights, class lag centres, convergence (all vs agreeing chains)
# and the MAP partition per fit. With --per_fit, also per-fit panels for every
# fit on disk: traces, class centres by chain, membership maps, class
# exposure-response and deviation sds. Builders: utils/figures_kl.R.
#
# Inputs:  <outdir>/kl_sweep_compare/<tag>/{fits,classes,partition,loo}.csv
#          (utils/paths.R::makeKlSweepRunFile()); --per_fit: the fits and their
#          kl_meta sidecars named by the configs.
# Outputs: figures/kl/<tag>/kl_sweep_*.png and, --per_fit,
#          figures/kl/<tag>/per_fit/<config stem>/*.png
#
# Usage:
#   Rscript analysis/06_fit_figures_kl_sweep.R -d analysis/configs/kl_sweep_classes_prodtu_all [--per_fit]
#   Rscript analysis/06_fit_figures_kl_sweep.R -d <dir1>,<dir2> --tag <name>   # as given to 05

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(ggplot2) })

option_list <- list(
  make_option(c("-d", "--batch_dir"), type = "character", default = NULL,
              help = "sweep config dir(s), comma-separated, as given to 05 (required)"),
  make_option("--tag", type = "character", default = "",
              help = "tag used by 05 [default: basename of the first dir]"),
  make_option("--per_fit", action = "store_true", default = FALSE,
              help = "also draw the s65 panels for every fit on disk"),
  make_option("--mode_tol", type = "double", default = 1.5,
              help = "chain-mode tolerance of the per-fit panels (lags) [default: %default]"),
  make_option("--format", type = "character", default = "png", help = "png | pdf [default: %default]"),
  make_option("--dpi", type = "integer", default = 150L, help = "[default: %default]"),
  make_option("--outdir", type = "character", default = "generated_data/",
              help = "root of the 05 outputs [default: %default]")
)

main <- function() {
  args <- parse_args(OptionParser(option_list = option_list))
  if (is.null(args$batch_dir) || !nzchar(args$batch_dir)) stop("-d/--batch_dir is required.")
  dirs <- trimws(strsplit(args$batch_dir, ",")[[1]])
  for (m in c("base", "config", "paths", "dlnm", "dlnm_kernel", "diagnostics", "kl_sweep",
              "figures_io", "figures_cb", "figures_kl"))
    source(sprintf("analysis/utils/%s.R", m))
  tag <- if (nzchar(args$tag)) args$tag else dirs[1]
  P <- function(kind) makeKlSweepRunFile(list(outdir = args$outdir), kind, batch = tag)
  need <- c("kl_sweep_fits", "kl_sweep_classes", "kl_sweep_partition", "kl_sweep_loo")
  miss <- need[!file.exists(vapply(need, P, ""))]
  if (length(miss))
    stop("05 outputs missing under ", P("kl_sweep_dir"), " (", paste(miss, collapse = ", "),
         "): run Rscript analysis/05_fit_postprocess_kl_sweep.R -d ", args$batch_dir,
         if (nzchar(args$tag)) paste(" --tag", args$tag) else "")
  rd <- function(kind) utils::read.csv(P(kind), stringsAsFactors = FALSE)
  fits <- rd("kl_sweep_fits"); classes <- rd("kl_sweep_classes")
  part <- rd("kl_sweep_partition"); loo <- rd("kl_sweep_loo")

  fig_dir <- makeFigureDir(list(), track = "kl", config_path = tag)
  save <- figureSaver(fig_dir, dpi = args$dpi, formats = args$format)
  nb <- length(unique(fits$basis))
  save(plot_kl_sweep_elpd(loo, "agree"), "kl_sweep_elpd_agree", width = 3 + 3 * nb, height = 4)
  save(plot_kl_sweep_elpd(loo, "all"), "kl_sweep_elpd_all", width = 3 + 3 * nb, height = 4)
  save(plot_kl_sweep_weights(loo), "kl_sweep_stacking_weights", width = 3 + 3 * nb, height = 3.8)
  save(plot_kl_sweep_centres(classes), "kl_sweep_centres", width = 3 + 3.5 * nb, height = 4.5)
  save(plot_kl_sweep_convergence(fits), "kl_sweep_convergence", width = 3 + 3.5 * nb, height = 5.5)
  save(plot_kl_sweep_partition(part, fits), "kl_sweep_partition",
       width = 3 + 0.35 * nrow(fits), height = 1.5 + 0.2 * length(unique(part$iso3)))

  if (args$per_fit) {
    inv <- kl_sweep_inventory(dirs)
    inv <- inv[inv$present, ]
    afr <- load_taxdat_afr_sf()
    for (i in seq_len(nrow(inv))) {
      r <- inv[i, ]
      opt <- NULL
      invisible(utils::capture.output(opt <- suppressMessages(make_options_from_config(
        option_list = list(), cli_argv = c("--config", r$config), model = "occurrence"))))
      fit  <- readRDS(r$fit_file)
      meta <- readRDS(makeRunFile(opt, "kl_meta"))
      sv <- figureSaver(file.path(fig_dir, "per_fit", r$model), dpi = args$dpi, formats = args$format)
      K <- fit$K
      sv(plot_kl_traces(fit), "traces", width = 2 + 2.6 * K, height = 5)
      sv(plot_kl_centres_by_chain(fit), "centres_by_chain", width = 2 + 2.6 * K, height = 3.4)
      sv(plot_kl_membership_maps(fit, meta, afr, args$mode_tol), "membership_maps",
         width = 11, height = 3.4 * ceiling((3 + fit$chains) / 3))
      sv(plot_kl_class_erf(fit, meta, args$mode_tol), "class_erf", width = 2 + 2.6 * K, height = 3.4)
      re <- plot_kl_re_sds(fit, meta)
      if (!is.null(re)) sv(re, "re_sds", width = 2 + 2 * fit$data_dims$J, height = 3.2)
      cl <- plot_kl_country_lag(fit)                     # kernel.lag_re fits only
      if (!is.null(cl)) sv(cl, "country_lag", width = 7, height = 1.5 + 0.22 * length(fit$country_levels))
      message("[06 kl_sweep] per-fit panels: ", r$id)
    }
  }
  message("[06 kl_sweep] figures in ", fig_dir)
}

if (sys.nframe() == 0L) main()
