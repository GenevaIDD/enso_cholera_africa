# test_spatpomp_unit_order.R — HARD: spatPomp builders must not depend on the
#   ROW ORDER of frame_tu.
#
# Usage:
#   Rscript analysis/tests/test_spatpomp_unit_order.R
#
# Details:
#   No fits, no MCMC. spatPomp holds per-unit states in C-locale-sorted gid
#   order (build_obs_for_tu()'s order), so every positional per-unit vector a
#   builder hands its kernel (warm start, LP inputs, lag aggregation, size_phi,
#   the nbr ADM1 kernel) must be in that order. build_spatpomp_tu_c() sorts
#   (frame_sorted); filtered_states_to_warmup() must map states by unit name.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
})
source("analysis/utils/base.R") # .TTLO_LEVELS (simulator_core.R TTLO bins; base.R is the universal floor)
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")
source("analysis/utils/spatpomp_core_c.R") # also sources spatpomp_core.R
source("analysis/tests/helpers.R")
.passes <- 0L
.fails <- 0L
cat("# test_spatpomp_unit_order\n\n")

n_months <- 4L
months <- as.Date("2022-01-01") %m+% months(0:(n_months - 1L))
frame_unsorted <- expand.grid(
  country = c("AAA", "BBB"),
  adm1n = 1:2,
  adm2n = 1:2,
  stringsAsFactors = FALSE
) %>%
  mutate(
    tu = "hmmW",
    gid = sprintf("%s_%d_%d", country, adm1n, adm2n),
    adm1_unit = sprintf("%s_%d", country, adm1n),
    AFRO_region = "REG_hmmW",
    min_seas = 1L,
    pop = 1e5
  ) %>%
  select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
  as_tibble()
frame_sorted <- frame_unsorted[order(frame_unsorted$gid, method = "radix"), ]
stopifnot(!identical(frame_unsorted$gid, frame_sorted$gid)) # fixture must be unsorted
panel <- tidyr::expand_grid(gid = frame_unsorted$gid, t_start = months) %>%
  mutate(
    min_seas = 1L,
    country = sub("_[0-9]+_[0-9]+$", "", gid),
    t_start_seas = as.integer(yday(t_start) - 1L),
    mean_precip = 0,
    dev_precip = 0,
    mean_temp_daymean = 0,
    dev_temp_daymean = 0,
    outbreak = 0L,
    time_to_last_outbreak_cat = "[6,Inf]",
    ongoing_outbreak_lag1 = 0L,
    ongoing_outbreak_country_lag1 = 0L,
    ongoing_outbreak_tu_lag1 = 0L
  ) %>%
  as_tibble()
window <- list(
  start = as.character(months[1]),
  end = as.character(months[n_months] %m+% months(1L))
)

# Base LP -30 (never fires); +60 on ADM1 AAA_1 only (always fires). The RE
# matrix is keyed by re_level_key(): <region>_<country>_<adm1_unit>.
draws <- make_draws_bundle_synthetic(
  frame_unsorted,
  S = 2L,
  basehaz_const = 0,
  fe_overrides = c("(Intercept)" = -30)
)
adm1_keys <- unique(re_level_key(frame_unsorted, "adm1"))
draws$re_intercept$adm1 <- matrix(
  0,
  nrow = 2L,
  ncol = length(adm1_keys),
  dimnames = list(NULL, adm1_keys)
)
draws$re_intercept$adm1[, "REG_hmmW_AAA_AAA_1"] <- 60
fires <- c("AAA_1_1", "AAA_1_2")

sim_states <- function(obj, seed = 1L) {
  set.seed(seed)
  sm <- spatPomp::simulate(obj, format = "data.frame")
  get <- function(v) {
    m <- sapply(obj@unit_names, function(u) sm[[v]][sm$unitname == u])
    if (is.null(dim(m))) {
      m <- matrix(m, nrow = 1L, dimnames = list(NULL, obj@unit_names))
    }
    m
  }
  list(Yout = get("Yout"), o = get("o"), o_ct = get("o_ct"))
}
check_answer <- function(st, label) {
  want <- as.numeric(colnames(st$Yout) %in% fires)
  expect_true(
    all(sweep(st$Yout, 2, want, "==")),
    sprintf("%s: Yout fires exactly AAA_1_*", label)
  )
  expect_true(
    all(sweep(st$o, 2, want, "==")),
    sprintf("%s: o = 1 exactly for AAA_1_*", label)
  )
  # o_ct: leave-one-ADM1-out within country -> AAA_2_* see AAA_1 active; BBB none.
  want_ct <- as.numeric(colnames(st$o_ct) %in% c("AAA_2_1", "AAA_2_2"))
  expect_true(
    all(sweep(st$o_ct, 2, want_ct, "==")),
    sprintf("%s: o_ct = 1 exactly for AAA_2_*", label)
  )
}

# Fixture: 8 ADM2 units in 4 ADM1s, frame rows in expand.grid order (NOT
# sorted). ADM1 random intercepts make AAA_1 certain to onset every month and
# every other ADM1 impossible, so the correct answer is known exactly:
#   Yout = 1 and o = 1 exactly for AAA_1_1 and AAA_1_2, 0 elsewhere.
# T2: C builder, unsorted frame: Yout / o / o_ct match the known answer.
cat("## T2: known answer with an UNSORTED frame (C builder)\n")
obj_c_un <- build_spatpomp_tu_c(
  frame_unsorted,
  panel,
  draws,
  window,
  fixed_draw_idx = 1L
)
st_c_un <- sim_states(obj_c_un)
check_answer(st_c_un, "T2 C builder")

# T4: filtered_states_to_warmup(): per-gid warm start identical for sorted vs
#     unsorted frame input, and AAA_1_* carry ongoing_outbreak_lag1 = 1.
cat("\n## T4: filtered_states_to_warmup() maps states by unit name\n")
bpf <- run_bpfilter_tu(
  obj_c_un,
  Np = 20L,
  block_list = make_block_list(obj_c_un, frame_unsorted),
  seed = 3L
)
w_un <- filtered_states_to_warmup(
  bpf,
  frame_unsorted,
  months[n_months] %m+% months(1L)
)[[1]]
w_so <- filtered_states_to_warmup(
  bpf,
  frame_sorted,
  months[n_months] %m+% months(1L)
)[[1]]
expect_true(
  identical(arrange(w_un, gid), arrange(w_so, gid)),
  "T4 per-gid warm start identical for sorted vs unsorted frame"
)
expect_true(
  all(w_un$ongoing_outbreak_lag1[match(fires, w_un$gid)] == 1L) &&
    all(w_un$ongoing_outbreak_lag1[!(w_un$gid %in% fires)] == 0L),
  "T4 ongoing_outbreak_lag1 = 1 exactly for AAA_1_*"
)

cat(sprintf("\n%d pass, %d fail\n", .passes, .fails))
if (.fails > 0L) {
  quit(status = 1L)
}
