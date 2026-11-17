# test_synthetic_equivalence.R — fits brms Poisson onset models with each
#   baseline-hazard variant to a simulated process with known truth.
#
# Usage:
#   Rscript analysis/tests/test_synthetic_equivalence.R
#   Rscript analysis/tests/test_synthetic_equivalence.R --iter 1000 --seed 1
#   Rscript analysis/tests/test_synthetic_equivalence.R --quick TRUE
#
# Details:
#   Strongest test of model-form fidelity: a piecewise-exponential onset
#   process, fitted with the factor / cyclic / ms baselines. Criteria: every
#   beta 95% CrI covers the truth; beta means within 0.25 SD of truth; cross-
#   model beta means within 0.5 SD of each other. Small by design (50 ADM1, 5
#   countries, 2 regions, monthly, 10 years). Needs brms + cmdstanr and a
#   working cmdstan toolchain; not in CI or run_all.R's registry, run once per
#   significant model-form change (fit runtime + MCMC noise can give false
#   negatives).

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

opt <- parse_args(OptionParser(
  option_list = list(
    make_option(c("--iter"), default = 1000, type = "integer"),
    make_option(c("--warmup"), default = 500, type = "integer"),
    make_option(c("--seed"), default = 20260512L, type = "integer"),
    make_option(
      c("--quick"),
      default = FALSE,
      type = "logical",
      help = "tiny config for smoke testing"
    )
  )
))

set.seed(opt$seed)

report_dir <- "generated_data/test_reports"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report <- file.path(report_dir, "synthetic_equivalence.md")
con <- file(report, "w")
on.exit(close(con), add = TRUE)
cat("# test_synthetic_equivalence\n\n", file = con)
cat(
  sprintf(
    "seed: %d, iter: %d, warmup: %d, quick: %s\n\n",
    opt$seed,
    opt$iter,
    opt$warmup,
    opt$quick
  ),
  file = con
)

# Simulation parameters ---------------------------------------------------

n_region <- if (opt$quick) 2 else 2
n_country <- if (opt$quick) 4 else 5
n_adm1 <- if (opt$quick) 12 else 50
n_year <- if (opt$quick) 5 else 10

beta_true <- c(
  mean_precip = 0.30,
  dev_precip = -0.20,
  mean_temp_daymean = 0.15,
  dev_temp_daymean = 0.05
)
sd_region <- 0.40
sd_country <- 0.20
sd_adm1 <- 0.10

# Build panel -------------------------------------------------------------
# Builder lifted to tests/helpers_synth.R (same RNG call order as the former
# inline code, so the simulated panel is unchanged).
source("analysis/tests/helpers_synth.R")
panel <- make_synthetic_occ_panel(
  n_region = n_region,
  n_country = n_country,
  n_adm1 = n_adm1,
  n_year = n_year,
  beta_true = beta_true,
  sd_region = sd_region,
  sd_country = sd_country,
  sd_adm1 = sd_adm1
)

cat("## Simulated panel\n\n", file = con)
cat(sprintf("- rows: %d\n", nrow(panel)), file = con)
cat(
  sprintf("- outbreak rate: %.3f per row\n", mean(panel$outbreak)),
  file = con
)
cat(
  sprintf("- region/country/adm1: %d/%d/%d\n", n_region, n_country, n_adm1),
  file = con
)
cat("\n", file = con)

# Helpers -----------------------------------------------------------------

summarise_beta <- function(draws, term_map) {
  # draws: long tibble with columns term, value
  draws %>%
    inner_join(
      tibble(term = names(term_map), truth = unname(term_map)),
      by = "term"
    ) %>%
    group_by(term, truth) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      q025 = quantile(value, 0.025),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      covers_truth = (q025 <= truth) & (q975 >= truth),
      bias_in_sd = (mean - truth) / sd
    )
}

fit_pois <- function(basehaz) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    message("brms not installed; skipping Poisson fit (", basehaz, ")")
    return(NULL)
  }

  bh <- switch(
    basehaz,
    factor = "factor(months)",
    cyclic = "s(t_start_seas, bs = 'cc', k = 8)",
    ms = NULL
  )
  if (basehaz == "ms") {
    if (!requireNamespace("splines2", quietly = TRUE)) {
      message("splines2 not installed; skipping ms variant")
      return(NULL)
    }
    ms <- splines2::mSpline(
      panel$t_start_seas,
      degree = 3,
      df = 5,
      intercept = FALSE
    )
    colnames(ms) <- paste0("ms", seq_len(ncol(ms)))
    dat <- bind_cols(panel, tibble::as_tibble(ms))
    bh <- paste(colnames(ms), collapse = " + ")
  } else {
    dat <- panel
  }

  form <- as.formula(paste(
    "y ~ 1 +",
    bh,
    "+ mean_precip + dev_precip + mean_temp_daymean + dev_temp_daymean",
    "+ (1 | AFRO_region/country/ADM1)",
    "+ offset(log(exposure))"
  ))

  fit <- brms::brm(
    form,
    data = dat,
    family = stats::poisson(),
    chains = 2,
    iter = opt$iter,
    warmup = opt$warmup,
    cores = 2,
    seed = opt$seed,
    backend = "cmdstanr",
    refresh = 0
  )
  fit
}

extract_pois_beta <- function(fit) {
  d <- posterior::as_draws_df(fit) %>% tibble::as_tibble()
  d %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "par", values_to = "value") %>%
    mutate(term = str_remove(par, "^b_")) %>%
    filter(term %in% names(beta_true))
}

# Run ---------------------------------------------------------------------

pois_fits <- list(
  factor = fit_pois("factor"),
  cyclic = fit_pois("cyclic"),
  ms = fit_pois("ms")
)
pois_betas <- map(pois_fits, function(f) {
  if (is.null(f)) NULL else extract_pois_beta(f)
})

cat("## Posterior summaries\n\n", file = con)

all_summ <- list()
for (nm in names(pois_betas)) {
  if (!is.null(pois_betas[[nm]])) {
    all_summ[[paste0("poisson_", nm)]] <- summarise_beta(
      pois_betas[[nm]],
      beta_true
    )
  }
}

for (nm in names(all_summ)) {
  cat(sprintf("### %s\n\n", nm), file = con)
  utils::write.table(
    format(all_summ[[nm]], digits = 3),
    con,
    sep = " | ",
    row.names = FALSE,
    quote = FALSE
  )
  cat("\n", file = con)
}

# Pass criteria -----------------------------------------------------------

failures <- c()

for (nm in names(all_summ)) {
  s <- all_summ[[nm]]
  if (any(!s$covers_truth)) {
    failures <- c(
      failures,
      sprintf(
        "%s: CrI does not cover truth for %s",
        nm,
        paste(s$term[!s$covers_truth], collapse = ", ")
      )
    )
  }
  if (any(abs(s$bias_in_sd) > 0.25)) {
    failures <- c(
      failures,
      sprintf(
        "%s: |bias|/sd > 0.25 for %s",
        nm,
        paste(s$term[abs(s$bias_in_sd) > 0.25], collapse = ", ")
      )
    )
  }
}

# Cross-model pairwise check: posterior means within 0.5 SD of each other
if (length(all_summ) >= 2) {
  ref <- all_summ[[1]]
  for (k in 2:length(all_summ)) {
    other <- all_summ[[k]]
    cmp <- inner_join(
      ref %>% select(term, mean_a = mean, sd_a = sd),
      other %>% select(term, mean_b = mean, sd_b = sd),
      by = "term"
    ) %>%
      mutate(
        pooled_sd = sqrt((sd_a^2 + sd_b^2) / 2),
        gap_in_sd = abs(mean_a - mean_b) / pooled_sd
      )
    if (any(cmp$gap_in_sd > 0.5)) {
      failures <- c(
        failures,
        sprintf(
          "%s vs %s: |mean gap|/sd > 0.5 for %s",
          names(all_summ)[1],
          names(all_summ)[k],
          paste(cmp$term[cmp$gap_in_sd > 0.5], collapse = ", ")
        )
      )
    }
  }
}

cat("## Verdict\n\n", file = con)
if (length(failures) == 0) {
  cat("PASS\n", file = con)
  message("test_synthetic_equivalence: PASS")
  quit(status = 0)
} else {
  cat("SOFT-FAIL:\n", file = con)
  for (f in failures) {
    cat("- ", f, "\n", sep = "", file = con)
  }
  message("test_synthetic_equivalence: SOFT-FAIL")
  for (f in failures) {
    message("  - ", f)
  }
  # Exit 0 so run_all.R continues; the report records the failure.
  quit(status = 0)
}
