# 05_fit_postprocess_kl_sweep.R — comparison of a brms_kl latent-class sweep
#
# The logic lives in utils/kl_sweep.R.
#
# Per fit (brms_kl occurrence fits of one or more sweep config dirs; fits not
# on disk are listed and skipped): sampler diagnostics; chain modes (chains
# whose per-chain mean class centres kmu all lie within --mode_tol lags,
# complete linkage); the "agreeing" chains = the dominant mode, or with
# --n_agree n the n best-agreeing chains; R-hat and bulk / tail ESS on all and
# on the agreeing chains (utils/diagnostics.R::fit_diagnostics()); class
# centres, El Nino contrasts, weights and kappa; MAP partition; membership
# uncertainty: the country x class responsibility matrix (all and agreeing
# chains, with each country's prior registry group), the M1 country mixture
# (per draw sum_k r_ck * kcon_k and sum_k r_ck * kmu_k vs the hard-MAP class,
# agreeing chains) and partition stability (per draw a partition sampled from
# Cat(r_c), adjusted Rand index vs the MAP; all and agreeing chains).
# Across fits (same rows and countries, checked): country-block PSIS-LOO and
# WAIC on all and on the agreeing chains (descriptive: high Pareto k is
# expected), loo_compare differences and stacking weights, per-country elpd.
#
# Inputs:  analysis/configs/<sweep batch>/*.yml (manifest.yml ignored) and the
#          fits they name (makeRunFile(opt, "fit")).
# Outputs: <outdir>/kl_sweep_compare/<tag>/{fits,classes,partition,loo,
#          loo_country,diagnostics,diagnostics_chains,responsibilities,
#          country_mixture,stability}.csv + loo.rds
#          (utils/paths.R::makeKlSweepRunFile()); read by
#          06_fit_figures_kl_sweep.R and the 09 report_kl_sweep step.
#
# Usage:
#   Rscript analysis/05_fit_postprocess_kl_sweep.R -d analysis/configs/kl_sweep_prodtu_classes_basis
#   Rscript analysis/05_fit_postprocess_kl_sweep.R -d <dir1>,<dir2> --tag <name> [--n_agree 4]
# HPC: bash analysis/bash/stages/04_fit_kl.sh --batch <dir> --postprocess

suppressPackageStartupMessages({ library(optparse); library(dplyr) })

option_list <- list(
  make_option(c("-d", "--batch_dir"), type = "character", default = NULL,
              help = "sweep config dir(s), comma-separated (required)"),
  make_option("--tag", type = "character", default = "",
              help = "output tag [default: basename of the first dir]"),
  make_option("--mode_tol", type = "double", default = 1.5,
              help = "max |difference| of per-chain mean class centres (lags) within one mode [default: %default]"),
  make_option("--n_agree", type = "integer", default = 0L,
              help = "keep the n best-agreeing chains per fit instead of the dominant mode (0 = dominant mode) [default: %default]"),
  make_option("--outdir", type = "character", default = "generated_data/",
              help = "output root [default: %default]")
)

main <- function() {
  args <- parse_args(OptionParser(option_list = option_list))
  if (is.null(args$batch_dir) || !nzchar(args$batch_dir)) stop("-d/--batch_dir is required.")
  dirs <- trimws(strsplit(args$batch_dir, ",")[[1]])
  miss <- dirs[!dir.exists(dirs)]
  if (length(miss)) stop("batch dir(s) not found: ", paste(miss, collapse = ", "))
  for (m in c("base", "config", "paths", "dlnm", "dlnm_kernel", "diagnostics", "kl_sweep"))
    source(sprintf("analysis/utils/%s.R", m))

  tag <- if (nzchar(args$tag)) args$tag else dirs[1]
  out <- kl_sweep_compare(dirs, mode_tol = args$mode_tol, n_agree = args$n_agree)
  P   <- function(kind) makeKlSweepRunFile(list(outdir = args$outdir), kind, batch = tag)
  dir.create(P("kl_sweep_dir"), recursive = TRUE, showWarnings = FALSE)
  for (nm in c("fits", "classes", "partition", "loo", "loo_country", "diagnostics", "diagnostics_chains",
                "responsibilities", "country_mixture", "stability"))
    utils::write.csv(out[[nm]], P(paste0("kl_sweep_", nm)), row.names = FALSE)
  saveRDS(out[c("inv", "agree", "res", "settings")], P("kl_sweep_rds"))

  r2 <- function(d) mutate(d, across(where(is.numeric), ~ round(.x, 2)))
  cat("\n== Sampler and convergence (agreeing =",
      if (args$n_agree > 0) sprintf("the %d best-agreeing chains", args$n_agree) else "dominant mode",
      "; mode tolerance", args$mode_tol, "lags)\n")
  print(as.data.frame(r2(select(out$fits, id, chains, mode_sizes, agreeing, kmu_spread_agree,
                                divergences, treedepth_hits, min_ebfmi, max_rhat_all,
                                max_rhat_agree, n_rhat105_agree, min_ess_bulk_agree))), row.names = FALSE)
  for (s in c("all", "agree")) {
    cat(sprintf("\n== Country-block LOO / WAIC, %s chains (descriptive)\n", s))
    print(as.data.frame(r2(out$loo %>% filter(chains == s) %>% arrange(desc(elpd_diff_loo)) %>%
      select(id, elpd_diff_loo, se_diff_loo, stacking_weight, elpd_diff_waic, se_diff_waic,
             p_loo, n_k07))), row.names = FALSE)
  }
  cat("\n== MAP partition (agreeing chains; p = posterior probability of the class)\n")
  pw <- out$partition %>% mutate(v = sprintf("%s %.2f", group, p)) %>% select(id, iso3, v) %>%
    tidyr::pivot_wider(names_from = id, values_from = v)
  print(as.data.frame(pw), row.names = FALSE)
  cat("\nWrote", P("kl_sweep_dir"), "\n")
  invisible(out)
}

if (sys.nframe() == 0L) main()
