# test_results_loader_fast.R — HARD: the fast paths of the results loaders are
#   identical to the dplyr / outer() code they replace.
#
# Usage:
#   Rscript analysis/tests/test_results_loader_fast.R
#
# Details:
#   Quick, synthetic inputs only. Fast paths under test:
#     .sum_by_keys()                   data.table keyby, dplyr's row / col order
#     .prob_pairs_gt() / _lt()         sort-and-count for mean(outer(f, g, op))
#     gid -> country / ADM1 parsing    once per distinct gid
#   Why: at 10^4 replicates per arm, load_sim_run_monthly()'s ADM1 x month
#   group_by took 489 s per arm in dplyr (4 s in data.table), and
#   adm1_prob_stats()'s outer() built a 10^8-cell matrix per ADM1.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(qs2)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/results_core.R")

# [1] .sum_by_keys() == dplyr group_by + summarise (character / Date / factor /
#     logical keys, integer and double sums, an empty frame).
cat("\n[1] .sum_by_keys() == dplyr\n")
set.seed(1)
n <- 5000L
df <- tibble(
  replicate_id = sample(1:40, n, TRUE),
  month = as.Date("2026-05-01") + 31L * sample(0:6, n, TRUE),
  country = sample(c("KEN", "AGO", "GHA", "BDI"), n, TRUE),
  run = factor(
    sample(c("forecast", "counterfactual"), n, TRUE),
    levels = c("forecast", "counterfactual")
  ),
  flag = sample(c(TRUE, FALSE), n, TRUE),
  cases = as.numeric(rnbinom(n, mu = 20, size = 0.4)),
  Y = rbinom(n, 1, 0.2)
)
ref <- function(d, keys) {
  d %>%
    group_by(across(all_of(keys))) %>%
    summarise(n_cases = sum(cases), n_outbreaks = sum(Y), .groups = "drop")
}
for (keys in list(
  c("replicate_id", "month", "country"),
  c("replicate_id", "month", "country", "flag"),
  c("run", "country", "replicate_id")
)) {
  got <- .sum_by_keys(df, keys, c(n_cases = "cases", n_outbreaks = "Y"))
  expect(
    identical(got, ref(df, keys)),
    sprintf("keys %s", paste(keys, collapse = " x "))
  )
}
expect(
  identical(
    .sum_by_keys(
      df[0, ],
      c("replicate_id", "country"),
      c(n_cases = "cases", n_outbreaks = "Y")
    ),
    ref(df[0, ], c("replicate_id", "country"))
  ),
  "empty frame"
)

# [2] .prob_pairs_gt() / _lt() bit-identical to mean(outer()) on random counts,
#     ties and zeros included.
cat("\n[2] .prob_pairs_gt() / _lt() == mean(outer())\n")
bad <- 0L
nt <- 0L
for (it in 1:500) {
  f <- rnbinom(
    sample(c(1:9, 50, 301), 1),
    mu = runif(1, 0, 40),
    size = runif(1, .2, 3)
  )
  g <- rnbinom(
    sample(c(1:9, 50, 299), 1),
    mu = runif(1, 0, 40),
    size = runif(1, .2, 3)
  )
  for (t in c(0, .25, .5, 1)) {
    nt <- nt + 2L
    bad <- bad +
      !identical(
        .prob_pairs_gt(f, (1 + t) * g),
        mean(outer(f, (1 + t) * g, `>`))
      ) +
        !identical(
          .prob_pairs_lt(f, (1 - t) * g),
          mean(outer(f, (1 - t) * g, `<`))
        )
  }
}
expect(bad == 0L, sprintf("%d comparisons bit-identical (%d differ)", nt, bad))

# [3] load_sim_run_monthly() and adm1_prob_stats() identical to their dplyr /
#     outer() versions on synthetic partitions: two TUs, a replicate chunk, NA
#     cases, Ghana's dot-less gids, horizon_by_country with and without
#     keep_beyond, series detected.
cat("\n[3] loaders == the dplyr / outer() versions (3884439)\n")
head_src <- tryCatch(
  system2(
    "git",
    c("show", "3884439:analysis/utils/results_core.R"),
    stdout = TRUE,
    stderr = FALSE
  ),
  warning = function(w) NULL,
  error = function(e) NULL
)
if (is.null(head_src) || any(grepl("^\\.sum_by_keys <- function", head_src))) {
  cat("  SKIP  3884439's results_core.R not readable (git unavailable?)\n")
} else {
  old <- new.env(parent = globalenv())
  exprs <- parse(text = head_src)
  for (e in exprs) {
    if (
      is.call(e) &&
        identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("load_sim_run_monthly", "adm1_prob_stats")
    ) {
      eval(e, envir = old)
    }
  }
  legacy_names(old)
  root <- tempfile("simroot_")
  run_tag <- "run_x"
  gids <- c(
    "KEN.1.2_1",
    "KEN.1.3_1",
    "KEN.4_1",
    "GHA5_2",
    "GHA14.13_2",
    "AGO.1.2_1",
    "AGO.2_1"
  )
  months <- format(seq(
    as.Date("2026-03-01"),
    as.Date("2027-02-01"),
    by = "month"
  ))
  write_part <- function(tu, off, reps, g) {
    d <- file.path(
      root,
      run_tag,
      "parts",
      paste0("tu=", tu),
      sprintf("repoff=%06d", off)
    )
    dir.create(d, recursive = TRUE)
    ms <- tidyr::expand_grid(
      replicate_id = off + seq_len(reps),
      gid = g,
      month = months
    ) %>%
      mutate(
        Y = rbinom(n(), 1, .15),
        active = Y,
        mu_mon = runif(n()),
        cases_month = ifelse(
          runif(n()) < .1,
          NA,
          rnbinom(n(), mu = 30, size = .5)
        ),
        detected = rbinom(n(), 1, .7)
      )
    qs2::qs_save(ms, file.path(d, "monthly_summary.qs"))
    file.create(file.path(d, "_SUCCESS"))
  }
  write_part("hmmA", 0L, 12L, gids[1:5])
  write_part("hmmA", 12L, 8L, gids[1:5])
  write_part("hmmB", 0L, 20L, gids[6:7])
  hz <- c(
    KEN = as.Date("2026-11-01"),
    GHA = as.Date("2027-01-01"),
    AGO = as.Date("2026-09-01")
  )
  args <- list(
    plain = list(),
    horizon = list(horizon_by_country = hz),
    beyond = list(horizon_by_country = hz, keep_beyond = TRUE),
    detected = list(series = "detected")
  )
  for (k in names(args)) {
    a <- c(
      list(
        sim_root = root,
        run_tag = run_tag,
        window_start = as.Date("2026-05-01"),
        window_end = as.Date("2027-01-01")
      ),
      args[[k]]
    )
    new_out <- suppressMessages(do.call(load_sim_run_monthly, a))
    old_out <- suppressMessages(do.call(old$load_sim_run_monthly, a))
    # max_gid_cases (per-gid explosive cap) is the one added
    # column; everything else is unchanged.
    mg <- new_out$country_month$max_gid_cases
    new_out$country_month$max_gid_cases <- NULL
    expect(
      identical(new_out, old_out),
      sprintf("load_sim_run_monthly (%s): every table identical", k)
    )
    if (identical(k, "plain")) {
      raw <- bind_rows(lapply(
        list.files(
          file.path(root, run_tag, "parts"),
          "monthly_summary.qs",
          recursive = TRUE,
          full.names = TRUE
        ),
        qs2::qs_read
      )) %>%
        mutate(
          country = substr(gid, 1, 3),
          month = as.Date(month),
          cases = ifelse(is.na(cases_month), 0, cases_month)
        ) %>%
        group_by(replicate_id, month, country) %>%
        summarise(mx = max(cases), .groups = "drop")
      got <- new_out$country_month %>%
        mutate(max_gid_cases = mg) %>%
        inner_join(raw, by = c("replicate_id", "month", "country"))
      expect(
        nrow(got) == nrow(new_out$country_month) &&
          isTRUE(all.equal(got$max_gid_cases, got$mx)),
        "load_sim_run_monthly: max_gid_cases = the largest gid's cases of the country-month"
      )
    }
  }
  am <- suppressMessages(load_sim_run_monthly(
    root,
    run_tag,
    as.Date("2026-05-01"),
    as.Date("2027-01-01")
  ))$adm1_month
  both <- bind_rows(
    mutate(am, run = "forecast"),
    mutate(
      am,
      run = "counterfactual",
      cases = round(cases * runif(n(), .5, 1.5))
    )
  )
  for (kk in c(TRUE, FALSE)) {
    for (t in c(0.25, 0.5)) {
      expect(
        identical(
          adm1_prob_stats(both, "cases", t, keep_means = kk),
          old$adm1_prob_stats(both, "cases", t, keep_means = kk)
        ),
        sprintf(
          "adm1_prob_stats (threshold %g, keep_means %s) identical",
          t,
          kk
        )
      )
    }
  }
  unlink(root, recursive = TRUE)
}

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
