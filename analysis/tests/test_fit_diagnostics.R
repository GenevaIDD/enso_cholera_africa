# test_fit_diagnostics.R — HARD: the utils/diagnostics.R::fit_diagnostics()
#   contract.
#
# Usage:
#   Rscript analysis/tests/test_fit_diagnostics.R
#   Rscript analysis/tests/test_fit_diagnostics.R --fit <brmsfit.rds>
#
# Details:
#   On a Stan-free fake kl_fit (tests/helpers_synth.R::make_fake_kl_fit()) and,
#   when the deployed best-lag archive is present locally (or --fit is given),
#   on one real deployed brmsfit (no sampling).

suppressPackageStartupMessages(library(optparse))
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/diagnostics.R")
source("analysis/tests/helpers_synth.R")

args <- parse_args(OptionParser(
  option_list = list(
    make_option(
      "--fit",
      default = "",
      type = "character",
      help = "a saved brmsfit to check [default: first duration fit of the newest local best-lag archive]"
    )
  )
))

FIELDS <- c(
  "summary",
  "sampler",
  "settings",
  "versions",
  "thresholds",
  "n_params",
  "max_rhat",
  "min_ess_bulk",
  "min_ess_tail",
  "n_rhat_gt",
  "n_ess_bulk_lt",
  "n_ess_tail_lt",
  "divergences",
  "treedepth_hits",
  "min_ebfmi"
)

cat("# [1] kl_fit branch (fake fit, 4 chains, 2 modes)\n")
fk <- make_fake_kl_fit(
  rbind(c(0, 9), c(0, 9), c(-6, 4), c(0.1, 9)),
  divergent = 2L
)
d <- fit_diagnostics(fk)
expect(identical(names(d), FIELDS), "fields in contract order")
expect(
  nrow(d$sampler) == 4L && identical(d$sampler$chain, 1:4),
  "one sampler row per chain"
)
expect(identical(d$divergences, 8L), "divergences summed over chains")
expect(
  d$n_params == nrow(d$summary) &&
    !any(d$summary$variable %in% c("lp__", "lprior")),
  "n_params = summary rows; lp__ excluded"
)
expect(isTRUE(d$max_rhat > 1.1), "mixed modes give a large R-hat")
d2 <- fit_diagnostics(fk, chains = c(1, 2, 4))
expect(
  identical(d2$sampler$chain, c(1L, 2L, 4L)) && d2$settings$n_chains == 3L,
  "chain subset keeps the chosen chains"
)
expect(isTRUE(d2$max_rhat < d$max_rhat), "agreeing chains reduce R-hat")
expect(all(is.na(unlist(d$versions))), "kl_fit versions are NA")
r <- fit_diagnostics_row(d2)
expect(nrow(r) == 1L && r$n_chains == 3L, "fit_diagnostics_row: one row")
expect(
  inherits(tryCatch(fit_diagnostics(list()), error = identity), "error"),
  "unsupported class is an error"
)

cat("# [2] brmsfit branch (a saved deployed fit; skipped when none is local)\n")
f <- args$fit
if (!nzchar(f)) {
  arch <- Sys.glob("generated_data/enso_best_lag_fits*")
  arch <- arch[order(file.mtime(arch), decreasing = TRUE)]
  cand <- unlist(lapply(
    arch,
    list.files,
    pattern = "^outbreak_duration_brms_cb_fit_.*\\.rds$",
    full.names = TRUE
  ))
  f <- if (length(cand)) cand[1] else ""
}
if (!nzchar(f) || !file.exists(f)) {
  cat("  SKIP: no saved brmsfit found\n")
} else {
  fit <- readRDS(f)
  b <- fit_diagnostics(fit)
  expect(identical(names(b), FIELDS), "fields in contract order")
  expect(nrow(b$sampler) == brms::nchains(fit), "one sampler row per chain")
  expect(
    all(is.finite(b$sampler$ebfmi) & b$sampler$ebfmi > 0),
    "E-BFMI finite and > 0"
  )
  expect(
    identical(b$versions$brms, as.character(fit$version$brms)),
    "versions from fit$version"
  )
  expect(b$settings$n_draws == brms::ndraws(fit), "n_draws = ndraws(fit)")
  expect(
    is.finite(b$settings$max_treedepth),
    "max_treedepth read from the fit's control args"
  )
}

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
quit(status = if (FAIL > 0L) 1L else 0L)
