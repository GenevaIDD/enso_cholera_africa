# rain_enso.R — ENSO and rainfall, descriptive (Figure 1c rain row, tables T1a / T1b).
#
# Sourced-only module (sources nothing). The rainy-season climatology of a group
# of units and the rainfall-anomaly vs ONI correlations that 03_prepare_rain_adm1.R
# and 09's rain_enso / seasonality steps compute from climate_stats.sqlite.
# Descriptive only: no model term reads it (the deployed fits use ENSO alone).

# Rainy-season climatology -------------------------------------------------

#' rainyMonths
#' Rainy calendar months of one unit's monthly rainfall climatology.
#'
#' Simple rule: a calendar month is rainy when its climatological mean rainfall
#' exceeds the unit's own average month (annual total / 12). Two-season
#' (bimodal) regimes need no special handling: each wet spell above the average
#' month comes out as its own run of rainy months. A unit drier than
#' `min_annual_mm` per year has no rainy season (all FALSE), so hyper-arid
#' units do not report a "season" made of a few millimetres.
#'
#' @param mm numeric(12), climatological mean rainfall of January..December
#'   (e.g. mm/month over a baseline period).
#' @param min_annual_mm annual total below which the unit has no season.
#' @return logical(12), TRUE = rainy month; all NA when `mm` is not 12 finite values.
rainyMonths <- function(mm, min_annual_mm = 100) {
  if (length(mm) != 12 || any(!is.finite(mm))) return(rep(NA, 12))
  if (sum(mm) < min_annual_mm) return(rep(FALSE, 12))
  mm > mean(mm)
}

#' circularMonthRuns
#' Runs of consecutive TRUE calendar months, wrapping December -> January.
#'
#' @param on logical(12), January..December.
#' @return list of integer vectors of months (1..12), each in season order
#'   (a Nov-Mar run is c(11, 12, 1, 2, 3)); list() when no month is on,
#'   list(1:12) when all are.
circularMonthRuns <- function(on) {
  on <- as.logical(on) & !is.na(on)
  if (all(on)) return(list(1:12))
  if (!any(on)) return(list())
  starts <- which(on & !on[c(12, 1:11)])
  lapply(starts, function(s) { run <- s; m <- s %% 12 + 1
    while (on[m]) { run <- c(run, m); m <- m %% 12 + 1 }
    as.integer(run) })
}

#' rainSeasonsByGroup
#' Main rainy season(s) of groups of spatial units (e.g. the ENSO groups of
#' Figure 1d) from the units' monthly rainfall climatology.
#'
#' Method (simple rule, weighted vote):
#'   1. per unit, the rainy months by `rainyMonths()` (month mean > the unit's
#'      average month; units with < `min_annual_mm` per year have no season);
#'   2. per group and calendar month, the weighted share of the group's units in
#'      season, share = sum(w * rainy) / sum(w), weights `weight` (population
#'      by default, so sparsely populated arid units do not widen or
#'      narrow a group's season; weight = 1 gives the share of units). Units
#'      with no climatology or a missing / non-positive weight are left out;
#'      dry units count as not in season;
#'   3. a month is in the group's main season when share >= `share_cut`; each
#'      circular run of such months (`circularMonthRuns()`) is one season, so a
#'      bimodal group gets two and a season crossing the year end stays one.
#'
#' Pure: no I/O, base R only.
#'
#' @param clim data.frame with columns `unit`, `mo` (1..12) and `mm` (the
#'   climatological mean rainfall of that calendar month), one row per unit x month.
#' @param units data.frame with columns `unit`, `group` and `weight`.
#' @param share_cut weighted share of units in season at or above which a month
#'   belongs to the group's main season.
#' @param min_annual_mm passed to `rainyMonths()`.
#' @return list of
#'   `unit`: unit, group, mo, rainy (per unit x month);
#'   `month`: group, mo, share, n_units, weight_total, in_season, mm (the
#'     weighted mean climatology of the group);
#'   `seasons`: group, season (1, 2, ... in calendar order of the first month),
#'     months ("Nov-Dec-Jan-Feb-Mar"), first, last (calendar months), n_months.
rainSeasonsByGroup <- function(clim, units, share_cut = 0.5, min_annual_mm = 100) {
  stopifnot(all(c("unit", "mo", "mm") %in% names(clim)),
            all(c("unit", "group", "weight") %in% names(units)),
            !anyDuplicated(units$unit), share_cut > 0, share_cut <= 1)
  units <- units[is.finite(units$weight) & units$weight > 0 & !is.na(units$group), , drop = FALSE]
  clim  <- clim[clim$unit %in% units$unit, , drop = FALSE]
  per_unit <- lapply(split(clim, clim$unit), function(d) {
    mm <- rep(NA_real_, 12); mm[d$mo] <- d$mm
    data.frame(unit = d$unit[1], mo = 1:12, mm = mm, rainy = rainyMonths(mm, min_annual_mm))
  })
  u <- do.call(rbind, per_unit)
  u <- u[!is.na(u$rainy), , drop = FALSE]
  u <- merge(u, units[, c("unit", "group", "weight")], by = "unit")
  month <- do.call(rbind, lapply(split(u, list(u$group, u$mo), drop = TRUE), function(d)
    data.frame(group = d$group[1], mo = d$mo[1],
               share = sum(d$weight * d$rainy) / sum(d$weight),
               n_units = nrow(d), weight_total = sum(d$weight),
               mm = sum(d$weight * d$mm) / sum(d$weight))))
  month$in_season <- month$share >= share_cut
  month <- month[order(month$group, month$mo), , drop = FALSE]
  seasons <- do.call(rbind, lapply(split(month, month$group), function(d) {
    runs <- circularMonthRuns(d$in_season[order(d$mo)])
    if (!length(runs)) return(NULL)
    runs <- runs[order(vapply(runs, `[`, 1L, 1L))]
    data.frame(group = d$group[1], season = seq_along(runs),
               months = vapply(runs, function(m) paste(month.abb[m], collapse = "-"), ""),
               first = vapply(runs, `[`, 1L, 1L),
               last = vapply(runs, function(m) m[length(m)], 1L),
               n_months = lengths(runs))
  }))
  rownames(month) <- NULL; rownames(seasons) <- NULL
  list(unit = u[order(u$unit, u$mo), c("unit", "group", "mo", "rainy")],
       month = month[, c("group", "mo", "share", "n_units", "weight_total", "in_season", "mm")],
       seasons = seasons)
}


# Rainfall anomaly vs ONI (Figure 1c) ---------------------------------------
# The ADM1 series of 03_prepare_rain_adm1.R and the correlations, seasons and
# population weights of 09's rain_enso step.

#' adm1MonthlyPrecip
#' ADM1 x month rainfall from the rows of climate_stats.sqlite's
#' `climate_stats` table: an ADM1's own row where it has one, else the mean of
#' its ADM2 rows.
#'
#' Ghana's flattened gids lack the separator after the ISO code ("GHA14_2",
#' "GHA14_11_2"); it is inserted so the token count gives the admin level
#' everywhere (4 tokens = ADM2).
#'
#' @param pr data.frame with columns `gid`, `month` (days since 1970-01-01, as
#'   stored, or Date) and `precip`.
#' @return tibble adm1, month (Date), precip, src ("adm1_row" | "adm2_mean"),
#'   iso3, year, mo.
adm1MonthlyPrecip <- function(pr) {
  pr <- dplyr::as_tibble(pr) |>
    dplyr::mutate(gid = sub("^([A-Z]{3})([0-9])", "\\1_\\2", gid))
  tok <- strsplit(pr$gid, "_")
  pr <- pr |> dplyr::mutate(month = as.Date(month, origin = "1970-01-01"),
                            ntok  = lengths(tok),
                            adm1  = ifelse(ntok == 4, paste(sapply(tok, `[`, 1), sapply(tok, `[`, 2), sep = "_"),
                                           sub("_[0-9]+$", "", gid)),
                            level = ifelse(ntok == 4, "ADM2", "ADM1"))
  has_adm1 <- unique(pr$adm1[pr$level == "ADM1"])
  dplyr::bind_rows(
    pr |> dplyr::filter(level == "ADM1") |> dplyr::transmute(adm1, month, precip, src = "adm1_row"),
    pr |> dplyr::filter(level == "ADM2", !adm1 %in% has_adm1) |> dplyr::group_by(adm1, month) |>
      dplyr::summarise(precip = mean(precip, na.rm = TRUE), src = "adm2_mean", .groups = "drop")) |>
    dplyr::mutate(iso3 = substr(adm1, 1, 3), year = as.integer(format(month, "%Y")),
                  mo = as.integer(format(month, "%m")))
}

#' precipDevSeries
#' Standardised monthly rainfall deviations per ADM1, smoothed over 3 months.
#'
#' dev_sd = (precip - mean) / sd of the same calendar month over the baseline
#' years; calendar months with a baseline mean below `dry_mm_month` (or no
#' spread) are dropped (NA). D(c, y) is the centred 3-month running mean of
#' dev_sd, kept where at least 2 of the 3 months are present.
#'
#' @param adm1_m adm1MonthlyPrecip() output.
#' @param baseline integer(2), first and last baseline year.
#' @param dry_mm_month baseline mean (mm/month) below which a month is dropped.
#' @return list(D = months x ADM1 matrix, months = Date vector, units =
#'   ADM1 ids (the columns), mo, yr = calendar month / year of each row,
#'   dry_share = share of unit-months dropped).
precipDevSeries <- function(adm1_m, baseline = c(1990L, 2019L), dry_mm_month = 10) {
  clim <- adm1_m |> dplyr::filter(year >= baseline[1], year <= baseline[2]) |>
    dplyr::group_by(adm1, mo) |>
    dplyr::summarise(mu = mean(precip), sdv = sd(precip), .groups = "drop")
  dv <- adm1_m |> dplyr::inner_join(clim, by = c("adm1", "mo")) |>
    dplyr::mutate(dry = mu < dry_mm_month | !is.finite(sdv) | sdv <= 0,
                  dev_sd = dplyr::if_else(dry, NA_real_, (precip - mu) / sdv))
  months_all <- seq(min(dv$month), max(dv$month), by = "month")
  units <- sort(unique(dv$adm1))
  M <- matrix(NA_real_, length(months_all), length(units),
              dimnames = list(as.character(months_all), units))
  M[cbind(match(as.character(dv$month), rownames(M)), match(dv$adm1, units))] <- dv$dev_sd
  D <- M
  for (i in seq_len(nrow(M))) {
    X <- M[max(1, i - 1):min(nrow(M), i + 1), , drop = FALSE]
    k <- colSums(!is.na(X))
    D[i, ] <- ifelse(k >= 2, colMeans(X, na.rm = TRUE), NA)
  }
  list(D = D, months = months_all, units = units,
       mo = as.integer(format(months_all, "%m")), yr = as.integer(format(months_all, "%Y")),
       dry_share = mean(dv$dry))
}

#' precipOniCorr
#' Pearson correlation of series-mean rainfall deviations with lagged ONI, per
#' centre month and lag.
#'
#' The series of key k is the mean of D over the ADM1 columns whose key is k
#' (NA keys are left out). For centre month c and lag l (negative = ONI after
#' the season), r = cor(D(c, y), ONI at month c - l of year y) over all years
#' with both, when there are at least `min_n` of them; p_r from the t-test on
#' n - 2 df; r_lo / r_hi the Fisher-z 95% interval.
#'
#' @param dev precipDevSeries() output.
#' @param keys character, one key per column of dev$D (e.g. the ADM1's ENSO
#'   group, or its country).
#' @param oni data.frame year, month, oni.
#' @param lags integer lags.
#' @param min_n minimum number of years for r.
#' @return tibble series, centre, lag, n, r, p_r, r_lo, r_hi.
precipOniCorr <- function(dev, keys, oni, lags = -11:11, min_n = 30) {
  stopifnot(length(keys) == ncol(dev$D))
  ks <- sort(unique(stats::na.omit(keys)))
  S <- do.call(cbind, lapply(ks, function(k) {
    X <- dev$D[, which(keys == k), drop = FALSE]; n <- rowSums(!is.na(X))
    ifelse(n > 0, rowMeans(X, na.rm = TRUE), NA) }))
  colnames(S) <- ks
  oni_key <- paste(oni$year, oni$month)
  oni_at <- function(y, m) { yy <- y + (m - 1) %/% 12; mm <- (m - 1) %% 12 + 1
    oni$oni[match(paste(yy, mm), oni_key)] }
  dplyr::bind_rows(lapply(1:12, function(c) {
    rows <- which(dev$mo == c); yrs <- dev$yr[rows]
    dplyr::bind_rows(lapply(lags, function(l) {
      x <- oni_at(yrs, c - l)
      dplyr::bind_rows(lapply(colnames(S), function(s) {
        y <- S[rows, s]; ok <- is.finite(y) & is.finite(x)
        r <- if (sum(ok) >= min_n) cor(y[ok], x[ok]) else NA_real_
        dplyr::tibble(series = s, centre = c, lag = l, n = sum(ok), r = r,
                      p_r = if (is.na(r)) NA_real_ else
                        2 * pt(-abs(r * sqrt(sum(ok) - 2) / sqrt(1 - r^2)), sum(ok) - 2))
      }))
    }))
  })) |>
    dplyr::mutate(se_z = 1 / sqrt(replace(n - 3, n <= 3, NA)),   # no interval below 4 years
                  r_lo = tanh(atanh(r) - 1.96 * se_z), r_hi = tanh(atanh(r) + 1.96 * se_z),
                  se_z = NULL)
}

#' precipOniGroupCorr
#' precipOniCorr() per group, with the within-group agreement: the share of
#' the group's countries whose own series is significant (p < 0.05) with the
#' group's sign in that cell (agree_r, over n_ctry countries).
#'
#' @param dev precipDevSeries() output.
#' @param group named character, ISO3 -> group (the countries of the grouping).
#' @inheritParams precipOniCorr
#' @return tibble series (the group), centre, lag, n, r, p_r, r_lo, r_hi,
#'   n_ctry, agree_r.
precipOniGroupCorr <- function(dev, group, oni, lags = -11:11, min_n = 30) {
  iso_u <- substr(dev$units, 1, 3)
  st_g <- precipOniCorr(dev, unname(group[iso_u]), oni, lags, min_n)
  st_c <- precipOniCorr(dev, ifelse(iso_u %in% names(group), iso_u, NA), oni, lags, min_n) |>
    dplyr::mutate(group = unname(group[series]))
  agree <- st_c |> dplyr::select(series, group, centre, lag, r, p_r) |>
    dplyr::inner_join(st_g |> dplyr::select(group = series, centre, lag, gr = r),
                      by = c("group", "centre", "lag")) |>
    dplyr::group_by(group, centre, lag) |>
    dplyr::summarise(n_ctry = dplyr::n(), agree_r = mean(p_r < .05 & sign(r) == sign(gr), na.rm = TRUE),
                     .groups = "drop")
  st_g |> dplyr::left_join(agree, by = c("series" = "group", "centre", "lag"))
}

#' adm1PopulationFromGid
#' ADM1 population of one year from the gid_population table: the ADM1's own
#' row, else the sum of its ADM2 rows (the ADM1 row wins).
#'
#' @param pop gid_population table (gid, gadm_level, year, Q).
#' @param year population year.
#' @return tibble adm1, Q.
adm1PopulationFromGid <- function(pop, year) {
  pop <- pop |> dplyr::filter(.data$year == !!year) |>
    dplyr::mutate(gid = sub("^([A-Z]{3})([0-9])", "\\1_\\2", gid))
  dplyr::bind_rows(
    pop |> dplyr::filter(gadm_level == "ADM1") |> dplyr::transmute(adm1 = sub("_[0-9]+$", "", gid), Q),
    pop |> dplyr::filter(gadm_level == "ADM2") |> dplyr::mutate(adm1 = sub("_[0-9]+_[0-9]+$", "", gid)) |>
      dplyr::group_by(adm1) |> dplyr::summarise(Q = sum(Q), .groups = "drop")) |>
    dplyr::group_by(adm1) |> dplyr::slice(1) |> dplyr::ungroup()
}

#' adm1RainSeasons
#' Main rainy season(s) per group (rainSeasonsByGroup()) from the ADM1
#' climatology over the baseline years, weighted by population (`pop`) and
#' by ADM1 count (`adm1`).
#'
#' @param adm1_m adm1MonthlyPrecip() output.
#' @param group named character, ISO3 -> group.
#' @param pop_adm1 adm1PopulationFromGid() output.
#' @param baseline integer(2), first and last baseline year.
#' @param share_cut,min_annual_mm passed to rainSeasonsByGroup().
#' @return list(seasons = weighting, group, season, months, first, last,
#'   n_months; month = weighting, group, mo, share, n_units, weight_total,
#'   in_season, mm; n_units, n_pop_matched).
adm1RainSeasons <- function(adm1_m, group, pop_adm1, baseline = c(1990L, 2019L),
                            share_cut = 0.5, min_annual_mm = 100) {
  clim <- adm1_m |> dplyr::filter(year >= baseline[1], year <= baseline[2], iso3 %in% names(group)) |>
    dplyr::group_by(adm1, iso3, mo) |> dplyr::summarise(mu = mean(precip, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(group = unname(group[iso3])) |> dplyr::arrange(adm1, mo)
  units <- clim |> dplyr::distinct(adm1, group) |> dplyr::left_join(pop_adm1, by = "adm1")
  clim_in <- clim |> dplyr::transmute(unit = adm1, mo, mm = mu)
  res <- list(pop  = rainSeasonsByGroup(clim_in, units |> dplyr::transmute(unit = adm1, group, weight = Q),
                                        share_cut = share_cut, min_annual_mm = min_annual_mm),
              adm1 = rainSeasonsByGroup(clim_in, units |> dplyr::transmute(unit = adm1, group, weight = 1),
                                        share_cut = share_cut, min_annual_mm = min_annual_mm))
  list(seasons = dplyr::bind_rows(lapply(names(res), function(w)
         dplyr::mutate(res[[w]]$seasons, weighting = w, .before = 1))),
       month = dplyr::bind_rows(lapply(names(res), function(w)
         dplyr::mutate(res[[w]]$month, weighting = w, .before = 1))),
       n_units = nrow(units), n_pop_matched = sum(is.finite(units$Q)))
}
