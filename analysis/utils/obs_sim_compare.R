# obs_sim_compare.R — obs-vs-sim total-case comparison on a fair denominator,
#   plus the scalar scoring rules shared by the skill and validation metrics
#
# Provides:
#   obs_sim_bases()          observed totals on the naive / B / C bases
#   cover_shares(), sim_replicate_totals()   footprint shares, sim totals
#   build_mon_obs_frame()    observed frame (finest level per region-month)
#   crps_ens(), crps_fair(), pinball_loss()  scalar scoring rules
#   region_key_of()          ADM1 region key of a gid
#
# Details:
#   The math notebooks/hindcast_skill.qmd needs (quantiles, CRPS, per-country,
#   monthly). crps_ens() and region_key_of() are also the one copy used by
#   sim_metrics.R and hindcast_skill.R, so this module is in 07c's and 07d's
#   sourced closure. Callers source gid.R and pass the frames loaded with
#   population.R::load_gid_population() and sim_io.R::sim_io_read().
#
# See: obs_sim_bases() for the B / C denominators

#' region_key_of
#'
#' ADM1 region key (`ISO3_<n>`) from a gid, after dot->underscore flattening.
region_key_of <- function(g) {
  gid_to_adm1(g, "region_key_of")
}

#' map_distinct
#'
#' f(x) evaluated once per distinct value of x (a vector, or a data frame
#' compared row-wise) and spread back to every element, in x's order. The
#' replicate frames repeat every gid / month / (country, month) once per
#' replicate, so string and date work done per row is done R times over.
map_distinct <- function(x, f) {
  u <- vctrs::vec_unique(x)
  vctrs::vec_slice(f(u), vctrs::vec_match(x, u))
}

#' crps_ens
#'
#' Empirical ensemble CRPS of a replicate sample `x` against a scalar `y`
#' (energy form). Moved verbatim from the hindcast notebook's lp-cov inline def.
crps_ens <- function(x, y) {
  xs <- sort(x)
  n <- length(xs)
  mean(abs(xs - y)) - sum((2 * seq_len(n) - n - 1) * xs) / (n^2)
}

#' crps_fair
#'
#' "Fair" ensemble CRPS (Ferro 2014): the same sorted energy form as crps_ens()
#' with the unbiased pairwise term (divisor n(n-1) instead of n^2), so ensembles
#' of different sizes -- e.g. a 300-member forecast vs a 10-member climatology --
#' are scored on an equal footing. NA members are dropped; a 1-member ensemble
#' reduces to |x - y|.
crps_fair <- function(x, y) {
  xs <- sort(x[!is.na(x)])
  n <- length(xs)
  if (!n) {
    return(NA_real_)
  }
  if (n == 1L) {
    return(abs(xs - y))
  }
  mean(abs(xs - y)) - sum((2 * seq_len(n) - n - 1) * xs) / (n * (n - 1))
}

#' pinball_loss
#'
#' Quantile (pinball) loss of the predictive `tau`-quantile `q` against `y`.
#' Proper for the tau-quantile and robust to heavy tails (unlike raw CRPS on
#' case totals). Vectorised over q/y/tau.
pinball_loss <- function(q, y, tau) {
  d <- y - q
  pmax(tau * d, (tau - 1) * d)
}

#' LVL_RANK
#'
#' Admin-level rank used by the finest-level dedup in obs_sim_bases()
#' (higher = finer). Rows whose admin_level is not a name here get lvl = NA.
LVL_RANK <- c(ADM1 = 1L, ADM2 = 2L, ADM3 = 3L)

#' build_mon_obs_frame
#'
#' The observed monthly cases frame obs_sim_bases() expects (`mon_all`), built
#' exactly as notebooks/hindcast_skill.qmd's lp-cov chunk did: gid, admin_level,
#' month (floored), cases_month, region_key, lvl. `mon` is the monthly-cases
#' model frame (03_prepare_mon.R output) or its path.
#'
#' Rows whose admin_level is outside `lvl_rank` get lvl = NA. obs_sim_bases()
#' drops them from basis B / naive by default but keeps them in C; their count
#' is returned as attr(, "n_na_lvl") and reported with message() so the drop is
#' never silent.
build_mon_obs_frame <- function(mon, lvl_rank = LVL_RANK) {
  if (is.character(mon)) {
    mon <- readRDS(mon)
  }
  out <- mon |>
    dplyr::transmute(
      gid = as.character(gid),
      admin_level = as.character(admin_level),
      month = lubridate::floor_date(as.Date(month), "month"),
      cases_month = as.numeric(cases_month)
    ) |>
    dplyr::mutate(region_key = region_key_of(gid), lvl = lvl_rank[admin_level])
  n_na <- sum(is.na(out$lvl))
  if (n_na > 0) {
    message(
      sprintf(
        "[build_mon_obs_frame] %d row(s) with admin_level outside {%s} ",
        n_na,
        paste(names(lvl_rank), collapse = ", ")
      ),
      "get lvl = NA (dropped from basis B/naive, kept in C)."
    )
  }
  attr(out, "n_na_lvl") <- n_na
  out
}

#' cover_shares
#'
#' Per-region simulated/ADM1 population share for a sim footprint. `pop_y` is a
#' load_gid_population() frame ALREADY FILTERED to the window-start year (columns
#' gid, Q, gadm_level). Returns tibble(region_key, sim_pop, adm1_pop, s) with
#' s = min(sim_pop / adm1_pop, 1). sim_pop sums own-level Q over the sim gids;
#' adm1_pop is the authoritative ADM1-level population for the region_key.
cover_shares <- function(sim_gids, pop_y) {
  ratio <- function(a, b) ifelse(b > 0, a / b, NA_real_)
  pop_y <- dplyr::mutate(
    pop_y,
    gid = as.character(gid),
    region_key = region_key_of(gid)
  )
  gidpop <- dplyr::distinct(pop_y, gid, Q)
  adm1pop <- pop_y |>
    dplyr::filter(gadm_level == "ADM1") |>
    dplyr::group_by(region_key) |>
    dplyr::summarise(adm1_pop = sum(Q), .groups = "drop")
  dplyr::tibble(gid = as.character(sim_gids)) |>
    dplyr::left_join(gidpop, by = "gid") |>
    dplyr::mutate(region_key = region_key_of(gid)) |>
    dplyr::group_by(region_key) |>
    dplyr::summarise(sim_pop = sum(Q, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(adm1pop, by = "region_key") |>
    dplyr::mutate(s = pmin(ratio(sim_pop, adm1_pop), 1))
}

#' sim_replicate_totals
#'
#' Per-replicate windowed case total over a gid set. `sim_in` is the sim frame
#' already restricted to the hindcast window (columns gid, replicate_id,
#' cases_month). Returns a numeric vector of per-replicate sums (callers take
#' whatever mean/quantiles they need).
sim_replicate_totals <- function(sim_in, keep_gids) {
  sim_in |>
    dplyr::filter(gid %in% keep_gids) |>
    dplyr::group_by(replicate_id) |>
    dplyr::summarise(v = sum(cases_month, na.rm = TRUE), .groups = "drop") |>
    dplyr::pull(v)
}

# Why: the sim footprint only partially covers (by population) each ADM1 it
# touches, so a naive full-ADM1 obs total is the wrong denominator. Two honest
# fixes side by side:
#
#   B  region_key (ADM1) scale, POPULATION-DOWN-WEIGHTED. Full-ADM1 obs scaled
#      by the simulated population share s(r)=min(sim_pop/adm1_pop, 1) at
#      window-start year.
#   C  strict exact-gid INTERSECTION. Keep only obs recorded at exactly a sim gid.
#   naive  full-ADM1 reference (the unfair denominator B/C correct).
#
# Common support on each side = obs-positive units within the sim footprint.
#' obs_sim_bases
#'
#' Orchestrate the naive / B / C building blocks for one run. Returns the raw frames +
#' support sets so callers can
#' format totals, quantiles, CRPS, per-country rollups, and monthly trajectories.
#'
#' @param ms      sim frame from sim_io_read() (gid, month, replicate_id, cases_month)
#' @param mon_all obs monthly frame (gid, admin_level, month, cases_month,
#'                region_key, lvl) -- caller builds it with region_key_of + LVL_RANK
#' @param gid_pop full load_gid_population() frame (filtered to window-start year inside)
#' @param na_lvl  rows of mon_all with lvl = NA: "drop" (default, historical
#'                behaviour) excludes them from B/naive; "finest" treats them as
#'                the finest level present (they are always kept in C)
#'
#' @return list with:
#'   HS, HE, yr                 window start/end month, window-start year
#'   sim_gids, R_sim            sim footprint gids, touched region_keys
#'   cover, cover_popwt         per-region cover_shares() frame + pop-wt coverage scalar
#'   obs_month                  tibble(region_key, country, month, cases_month, s)
#'                              finest-admin-deduped per (region_key, month), s joined
#'   obs_gid                    tibble(gid, obs) per-gid obs at exact sim gids
#'   S_B, simB_gids             B support regions; sim gids within them
#'   G_obs                      C support (sim gids with exact obs > 0)
#'   obs_naive, obs_B, obs_C    the three observed window totals
#'   obs_rk                     tibble(region_key, obs_adm1, obs_fp) per touched
#'                              region: full-ADM1 and pop-down-weighted (s(r))
#'                              observed totals, zeros not included
#'   sim_in                     ms windowed, with region_key + country columns
obs_sim_bases <- function(ms, mon_all, gid_pop, na_lvl = c("drop", "finest")) {
  na_lvl <- match.arg(na_lvl)
  if (na_lvl == "finest") {
    mon_all <- dplyr::mutate(
      mon_all,
      lvl = dplyr::coalesce(lvl, .Machine$integer.max)
    )
  }
  ms <- dplyr::mutate(
    ms,
    gid = as.character(gid),
    month = map_distinct(month, function(m) {
      lubridate::floor_date(as.Date(m), "month")
    })
  )
  HS <- min(ms$month)
  HE <- max(ms$month)
  yr <- lubridate::year(HS)
  sim_gids <- sort(unique(ms$gid))
  R_sim <- unique(stats::na.omit(region_key_of(sim_gids)))

  # Nearest available population year: gid_pop covers only 2015-2030, so a pre-2015
  # (or post-2030) launch would otherwise find no rows -> empty cover -> NA coverage
  # and s coalesced to 1 (B silently collapses to naive). s(r)=sim_pop/adm1_pop is a
  # within-year ratio over the same gids and is nearly year-invariant, so the nearest
  # year is a safe proxy. In-range years match exactly (no-op).
  avail_yr <- sort(unique(gid_pop$year))
  yr_pop <- avail_yr[which.min(abs(avail_yr - yr))]
  pop_y <- dplyr::filter(gid_pop, year == yr_pop)
  cover <- cover_shares(sim_gids, pop_y)
  cover_popwt <- {
    sp <- sum(cover$sim_pop, na.rm = TRUE)
    ap <- sum(cover$adm1_pop, na.rm = TRUE)
    ifelse(ap > 0, sp / ap, NA_real_)
  }

  # sim per-replicate monthly cases in window (shared across B/C; carries keys)
  # (keys per distinct gid; the window rows are not copied when the window
  # keeps every row, as it does unless a month is NA)
  keep <- ms$month >= HS & ms$month <= HE
  keep[is.na(keep)] <- FALSE
  sim_in <- (if (all(keep)) ms else dplyr::filter(ms, !!keep)) |>
    dplyr::mutate(
      region_key = map_distinct(gid, region_key_of),
      country = map_distinct(region_key, function(r) stringr::str_sub(r, 1, 3))
    )

  # ── B: region_key, finest-admin dedup per (region_key, month), pop share s ──
  obs_month <- mon_all |>
    dplyr::filter(
      month >= HS,
      month <= HE,
      !is.na(region_key),
      region_key %in% R_sim,
      !is.na(lvl)
    ) |>
    dplyr::group_by(region_key, month) |>
    dplyr::filter(lvl == max(lvl)) |>
    dplyr::summarise(
      cases_month = sum(cases_month, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(dplyr::select(cover, region_key, s), by = "region_key") |>
    dplyr::mutate(
      s = dplyr::coalesce(s, 1),
      country = stringr::str_sub(region_key, 1, 3)
    )

  obs_rk <- obs_month |>
    dplyr::group_by(region_key) |>
    dplyr::summarise(
      obs_adm1 = sum(cases_month, na.rm = TRUE),
      obs_fp = sum(cases_month * s, na.rm = TRUE),
      .groups = "drop"
    )
  S_B <- obs_rk$region_key[obs_rk$obs_adm1 > 0]
  simB_gids <- sim_gids[region_key_of(sim_gids) %in% S_B]
  obs_naive <- sum(obs_rk$obs_adm1[obs_rk$region_key %in% S_B])
  obs_B <- sum(obs_rk$obs_fp[obs_rk$region_key %in% S_B])

  # ── C: strict exact-gid intersection ────────────────────────────────────────
  obs_gid <- mon_all |>
    dplyr::filter(month >= HS, month <= HE, gid %in% sim_gids) |>
    dplyr::group_by(gid) |>
    dplyr::summarise(obs = sum(cases_month, na.rm = TRUE), .groups = "drop")
  G_obs <- obs_gid$gid[obs_gid$obs > 0]
  obs_C <- sum(obs_gid$obs[obs_gid$gid %in% G_obs])

  list(
    HS = HS,
    HE = HE,
    yr = yr,
    sim_gids = sim_gids,
    R_sim = R_sim,
    cover = cover,
    cover_popwt = cover_popwt,
    obs_month = obs_month,
    obs_gid = obs_gid,
    S_B = S_B,
    simB_gids = simB_gids,
    G_obs = G_obs,
    obs_naive = obs_naive,
    obs_B = obs_B,
    obs_C = obs_C,
    obs_rk = obs_rk,
    sim_in = sim_in
  )
}
