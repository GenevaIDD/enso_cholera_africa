# test_incidence_cap.R — HARD: simulator.incidence_cap bounds monthly cases
#   per unit by incidence_cap x pop, in every kernel.
#
# Usage:
#   Rscript analysis/tests/test_incidence_cap.R
#
# Details:
#   Synthetic, no MCMC; run from the project root. The cap applies to the
#   monthly NB mean and to the drawn monthly case count: helpers (H), the
#   v1 / spatPomp forecast simulate_one_tu() (V), the spatPomp C filter
#   kernel (K) and the smoothing post-hoc draw (S). NULL and Inf reproduce
#   the uncapped output exactly.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(lubridate)
  library(tidyr)
})
source("analysis/utils/smoothing_core.R") # sources simulator_core.R, spatpomp_core{,_c}.R, dlnm.R
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# test_incidence_cap\n\n")

# ===========================================================================
# H -- helpers
# ===========================================================================
# H  helpers: monthly_mu_cap(), draw_monthly_cases(), incidence_cap_of(),
#    .make_filter_cache_path()'s `_cap-` token.
cat("## H: helpers\n")
pop <- c(1e3, 5e4, 2e5)
expect(
  identical(monthly_mu_cap(pop, 0.15), pop * 0.15),
  "H1 monthly_mu_cap = cap x pop"
)
expect(
  identical(monthly_mu_cap(pop), rep(Inf, 3)),
  "H2 monthly_mu_cap(NULL) = Inf"
)
expect(
  identical(monthly_mu_cap(pop, Inf), rep(Inf, 3)),
  "H3 monthly_mu_cap(Inf) = Inf"
)
expect(
  identical(monthly_mu_cap(pop, 0.15, 1.0), pop * 0.15) &&
    identical(monthly_mu_cap(pop, 0.5, 0.2), pop * 0.2) &&
    identical(monthly_mu_cap(pop, NULL, 1.0), pop * 1.0),
  "H4 monthly_mu_cap takes min(incidence_cap, mu_pop_frac)"
)

mu <- c(5, 4e4, 30, NA, 12)
phi <- c(0.4, 0.4, 0.4, 0.4, -1)
cap <- c(1e9, 150, 1e9, 1e9, 1e9)
set.seed(7L)
ref <- stats::rnbinom(3L, mu = mu[1:3], size = phi[1:3])
set.seed(7L)
d_inf <- draw_monthly_cases(mu, phi, Inf)
expect(
  identical(d_inf[1:3], ref) && all(is.na(d_inf[4:5])),
  "H5 draw_monthly_cases(cap = Inf) = plain rnbinom stream; NA where not drawable"
)
set.seed(7L)
d_cap <- draw_monthly_cases(mu, phi, cap)
expect(
  identical(d_cap[c(1, 3)], ref[c(1, 3)]) &&
    d_cap[2] == min(ref[2], 150) &&
    all(is.na(d_cap[4:5])),
  "H6 binding cap truncates only the capped unit, same stream"
)
set.seed(7L)
d_frac <- draw_monthly_cases(c(1e6, 1e6), c(1, 1), c(10.7, 10.7))
expect(all(d_frac == 10), "H7 draw is truncated at floor(cap)")

expect(
  identical(incidence_cap_of(NULL), 0.15) &&
    identical(incidence_cap_of(list()), 0.15),
  "H8 incidence_cap_of default 0.15"
)
expect(
  identical(incidence_cap_of(list(incidence_cap = 0.3)), 0.3) &&
    identical(incidence_cap_of(list(incidence_cap = Inf)), Inf) &&
    identical(incidence_cap_of(list(incidence_cap = 1L)), 1),
  "H9 incidence_cap_of accepts numbers > 0 and Inf"
)
bad <- list(0, -0.1, "0.15", c(0.1, 0.2), NA_real_, TRUE)
expect(
  all(vapply(
    bad,
    function(v) {
      inherits(
        try(incidence_cap_of(list(incidence_cap = v)), silent = TRUE),
        "try-error"
      )
    },
    logical(1)
  )),
  "H10 incidence_cap_of refuses 0, negative, string, vector, NA, logical"
)

ww <- list(start = "2009-01-01", end = "2011-01-01")
p_null <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww)
p_inf <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, Inf)
p_cap <- .make_filter_cache_path("TUA", "c", "occ", 100L, 4L, ww, 0.15)
expect(
  identical(p_null, p_inf) && !grepl("_cap-", p_null),
  "H11 filter cache key: no token when uncapped (NULL / Inf)"
)
expect(
  grepl("_warmup=20090101_20110101_cap-0.15\\.qs$", p_cap),
  "H12 filter cache key: `_cap-0.15` token under a finite cap"
)

# ===========================================================================
# Fixture: 2 countries x 2 ADM1 x 1 ADM2, unequal pop, exploding mon mean
# ===========================================================================
frame_tu <- tidyr::expand_grid(country = c("AAA", "BBB"), adm1n = 1:2) %>%
  dplyr::mutate(
    gid = sprintf("%s_%d_1", country, adm1n),
    adm1_unit = sprintf("%s_%d", country, adm1n),
    AFRO_region = "REG_W",
    tu = "hmmW",
    min_seas = 1L,
    pop = c(2e4, 1.8e5, 6e4, 3e5)
  ) %>%
  dplyr::select(gid, country, AFRO_region, tu, min_seas, pop, adm1_unit) %>%
  tibble::as_tibble()
cap_v <- 0.15 * frame_tu$pop
months <- as.Date("2024-01-01") %m+% months(0:5)
panel_tu <- tidyr::expand_grid(gid = frame_tu$gid, t_start = months) %>%
  dplyr::left_join(frame_tu %>% dplyr::select(gid, min_seas), by = "gid") %>%
  dplyr::mutate(
    country = sub("_[0-9]+_[0-9]+$", "", gid),
    t_start_seas = as.integer(lubridate::yday(t_start) - 1L),
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
  tibble::as_tibble()
horizon <- list(start = "2024-01-01", end = "2024-06-30")
window <- list(start = "2024-01-01", end = "2024-07-01")

# Onset certain in month 1, cessation ~impossible: every unit active throughout.
d_occ <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = 5)
)
d_dur <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -30),
  active_terms = "(Intercept)"
)
d_dur$basehaz <- list(kind = "factor", breaks = .SIM_DUR_MSO_BREAKS)
# Per-person-day rate e^-3 with pop + log_offset offsets: mu ~ 1.5 x pop, far
# above the 0.15 cap but finite (so the uncapped path is a real comparison).
d_mon <- make_draws_bundle_synthetic(
  frame_tu,
  S = 2L,
  fe_overrides = c("(Intercept)" = -3),
  active_terms = "(Intercept)",
  phi = 0.45
)
d_mon$offset_pop <- TRUE
d_mon$offset_logoffset <- TRUE

# ===========================================================================
# V -- simulate_one_tu()
# ===========================================================================
# V  simulate_one_tu() (v1 and the spatPomp forecast phase) under an
#    exploding monthly intercept: mu_mon <= cap x pop and cases_month <=
#    floor(cap x pop); NULL and Inf reproduce the uncapped output exactly.
cat("\n## V: simulate_one_tu()\n")
run_v1 <- function(...) {
  simulate_one_tu(
    tu = "hmmW",
    frame_tu = frame_tu,
    panel_tu = panel_tu,
    draws_occ = d_occ,
    draws_dur = d_dur,
    draws_mon = d_mon,
    replicate_ids = 1:5,
    seed = 99L,
    output_part = tempfile("ic_"),
    horizon = horizon,
    dry_run = TRUE,
    ...
  )$monthly_summary
}
ms_def <- suppressWarnings(run_v1())
ms_nul <- suppressWarnings(run_v1(incidence_cap = NULL))
ms_inf <- suppressWarnings(run_v1(incidence_cap = Inf))
ms_cap <- suppressWarnings(run_v1(incidence_cap = 0.15))
pop_r <- frame_tu$pop[match(ms_cap$gid, frame_tu$gid)]
expect(all(ms_cap$active == 1L), "V0 fixture: every unit active every month")
expect(
  identical(ms_def, ms_nul) && identical(ms_def, ms_inf),
  "V1 incidence_cap NULL / Inf reproduce the default (uncapped) kernel exactly"
)
expect(
  any(ms_def$mu_mon > 0.15 * pop_r),
  "V2 fixture: uncapped mu exceeds 0.15 x pop"
)
expect(all(ms_cap$mu_mon <= 0.15 * pop_r), "V3 capped mu_mon <= 0.15 x pop")
expect(
  isTRUE(all.equal(ms_cap$mu_mon, 0.15 * pop_r)),
  "V4 capped mu_mon pinned at 0.15 x pop"
)
expect(
  all(ms_cap$cases_month <= floor(0.15 * pop_r)),
  "V5 capped cases_month <= floor(0.15 x pop)"
)
expect(
  any(ms_def$cases_month > floor(0.15 * pop_r)),
  "V6 fixture: uncapped draws exceed the cap"
)
expect(
  identical(ms_cap$Y, ms_def$Y) && identical(ms_cap$active, ms_def$active),
  "V7 the cap leaves onsets and activity untouched"
)

# ===========================================================================
# K -- spatPomp C filter kernel
# ===========================================================================
# K  spatPomp C filter kernel (mon_mu_block): the `mu` state is pinned at
#    cap x pop; NULL leaves it uncapped; the uncapped C build emits no
#    mon_mu_cap covariate.
cat("\n## K: spatPomp C filter kernel\n")
mu_state <- function(obj) {
  set.seed(2024L)
  sm <- spatPomp::simulate(obj, format = "data.frame")
  U <- length(obj@unit_names)
  M <- length(attr(obj, "filt_months"))
  mat <- matrix(NA_real_, M, U, dimnames = list(NULL, obj@unit_names))
  for (tt in seq_len(M)) {
    for (i in seq_len(U)) {
      mat[tt, i] <- sm$mu[sm$time == tt & sm$unitname == obj@unit_names[i]][1L]
    }
  }
  mat
}
build_k <- function(builder, ...) {
  builder(
    frame_tu,
    panel_tu,
    d_occ,
    window,
    n_draws_in_filter = 1L,
    fixed_draw_idx = 1L,
    draws_dur = d_dur,
    draws_mon = d_mon,
    measurement = "occurrence_duration_mon_NB",
    ...
  )
}
o_c_cap <- build_k(build_spatpomp_tu_c, incidence_cap = 0.15)
o_c_nul <- build_k(build_spatpomp_tu_c)
o_c_inf <- build_k(build_spatpomp_tu_c, incidence_cap = Inf)
pop_u <- frame_tu$pop[match(o_c_cap@unit_names, frame_tu$gid)]
cap_m <- matrix(
  0.15 * pop_u,
  nrow = length(months),
  ncol = length(pop_u),
  byrow = TRUE
)
mc_cap <- mu_state(o_c_cap)
mc_nul <- mu_state(o_c_nul)
expect(
  isTRUE(all.equal(unname(mc_cap), unname(cap_m))),
  "K2 C kernel: mu pinned at 0.15 x pop"
)
expect(all(mc_nul > cap_m), "K4 NULL: C kernel uncapped")
has_cap_cov <- function(o) any(grepl("^mon_mu_cap", rownames(o@covar@table)))
expect(
  !has_cap_cov(o_c_nul) && identical(o_c_nul@covar@table, o_c_inf@covar@table),
  "K6 C build: no mon_mu_cap covariate when uncapped (NULL == Inf)"
)
expect(
  has_cap_cov(o_c_cap),
  "K7 C build: mon_mu_cap covariate present under a cap"
)

# ===========================================================================
# S -- smoothing post-hoc draw
# ===========================================================================
# S  smoothing post-hoc draw (extract_smoothing_one_tu()) on a C-kernel
#    bpfilter trajectory: cases <= floor(cap x pop); NULL == Inf.
cat("\n## S: extract_smoothing_one_tu()\n")
results <- lapply(1:2, function(r) {
  obj <- build_k(build_spatpomp_tu_c, incidence_cap = 0.15)
  res <- run_one_bpfilter_replicate(
    obj,
    Np = 20L,
    block_list = make_block_list(obj, frame_tu),
    seed = 500L + r
  )
  list(
    tu = "hmmW",
    rep_id = r,
    draw_d = 1L,
    filt_months = attr(obj, "filt_months"),
    filter_traj = res$filter_traj,
    logLik = res$logLik,
    gids = sort(as.character(frame_tu$gid), method = "radix")
  )
})
run_sm <- function(...) {
  set.seed(31L)
  sm <- extract_smoothing_one_tu(
    results,
    frame_tu,
    panel_tu,
    FILTER_START = as.Date("2024-01-01"),
    FILTER_END = as.Date("2024-07-01"),
    draws_mon = d_mon,
    ...
  )
  lapply(sm$reps, `[[`, "cases_mon_mat")
}
sm_nul <- run_sm()
sm_inf <- run_sm(incidence_cap = Inf)
sm_cap <- run_sm(incidence_cap = 0.15)
pop_s <- frame_tu$pop[match(sort(frame_tu$gid, method = "radix"), frame_tu$gid)]
capped_ok <- all(vapply(
  sm_cap,
  function(m) all(m <= floor(0.15 * pop_s), na.rm = TRUE),
  logical(1)
))
expect(sum(!is.na(unlist(sm_cap))) > 10, "S0 fixture: smoothing draws cases")
expect(identical(sm_nul, sm_inf), "S1 NULL / Inf identical")
expect(capped_ok, "S2 capped smoothing cases <= floor(0.15 x pop)")
expect(
  any(unlist(sm_nul) > max(floor(0.15 * pop_s)), na.rm = TRUE),
  "S3 fixture: uncapped smoothing draws exceed the cap"
)

cat(sprintf("\n# Results: %d PASS  %d FAIL\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
