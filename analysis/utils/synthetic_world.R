# synthetic_world.R — a small fake world whose raw inputs every pipeline step can read.
#
# Sourced-only module (sources nothing; callers source base, gid and enso first).
# Used by analysis/00_ingest_synthetic.R and tests/test_synthetic_world.R.
#
# The world is built from analysis/groupings/synthetic_countries.csv: fake
# countries with ISO 3166 user-assigned codes (XAA..XAD, never a real country)
# and invented ADM1 / ADM2 names, laid out as 1-degree squares over open ocean
# inside the Africa bounding box, so no real name, code or polygon appears in
# any synthetic file. It writes only RAW inputs, in the formats of the real
# sources, so the real 00 -> 02 steps run on them unmodified:
#
#   <root>/<vintage>/combined_outbreaks_cholera.csv, stage1_geo/*.geojson (OutbreakExtractR stage 2)
#   data/geodata/gadm_410-levels.gpkg, gadm_valid.gpkg                (GADM 4.1, ADM_1 / ADM_2)
#   data/climate_parquet/<ISO3>/Africa_<ISO3>_v410_..._observation.parquet.gzip (daily CHIRPS / ERA5-Land)
#   data/worldpop/<iso3>_pop_<year>_CN_100m_R2025A_v1.tif, data/country_boundaries/<ISO3>_adm0.geojson
#   data/who_annual_report_africa/*  (WHO annual report, dashboard annual / weekly, sitrep cumulative)
#   data/best_seasonality_estimates.csv, data/hmm_country_clusters.csv
#   data/lenssen/extracted/<model>/step01_InitVerifInfo.Rda            (piControl Nino3.4)
#   data/geodata/Africa_waterbody.shp, rivers_africa_37333.shp         (basemap layers)
#   generated_data/enso_{oni,nino34,nino34_oisst}_cache.rds, enso_cfsv2_cache.rds (+ download metas)
#   generated_data/wpp2024_national_totals.csv
#
# Generative model (monthly, per reporting unit u in country c of group g):
#   onset   P(onset | inactive) = 1 - exp(-h),
#           log h = a0 + a_c + a_u + s_g(month) + beta_g * ONI(t + L_g) + gamma * nbr_active(t - 1)
#   cessation P(end | active, m months in) = plogis(c0 + c1 * m - 0.3 * beta_g * ONI(t + L_g))
#   cases   NB(mean = pop_u * rate * exp(0.5 * s_g(month)), size 3), weekly split multinomial
# L_g and beta_g are the registry's enso_lag / enso_beta (negative L = ENSO leads),
# so the fits have a known answer to recover.

# ---------------------------------------------------------------------------
# Geography
# ---------------------------------------------------------------------------

#' synthRegistry
#' The synthetic country registry.
#' @param path analysis/groupings/synthetic_countries.csv
#' @return data.frame, one row per fake country.
synthRegistry <- function(path = "analysis/groupings/synthetic_countries.csv") {
  reg <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("iso3", "name", "group", "afro_region", "hmm_cluster", "lon0", "lat0",
            "enso_lag", "enso_beta", "season_peak", "rain_peak")
  miss <- setdiff(need, names(reg))
  if (length(miss)) stop("synthRegistry(): ", path, " lacks ", paste(miss, collapse = ", "))
  if (!all(grepl("^XA[A-Z]$", reg$iso3)))
    stop("synthRegistry(): synthetic codes must be ISO 3166 user-assigned XA? codes (got ",
         paste(reg$iso3, collapse = ", "), ")")
  reg
}

# One rectangle as an sf MULTIPOLYGON (GADM's geometry type).
.msquare <- function(x0, y0, w, h)
  sf::st_multipolygon(list(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0)))))

.ADM1_SUFFIX <- c("North", "South", "East", "West", "Central", "Upper", "Lower", "Coast")
.ADM2_SUFFIX <- c("A", "B", "C", "D")

#' synthGadm
#' GADM-like ADM1 / ADM2 units of the synthetic countries: each country a
#' 2 x 2 block of 1-degree ADM1 squares (lon0, lat0 = south-west corner), each
#' ADM1 split into `n_adm2` vertical strips.
#' @return list(adm1, adm2) of sf in GADM 4.1's layer schemas (EPSG:4326)
synthGadm <- function(reg, n_adm1 = 4L, n_adm2 = 2L) {
  stopifnot(n_adm1 <= length(.ADM1_SUFFIX), n_adm2 <= length(.ADM2_SUFFIX))
  ncol <- ceiling(sqrt(n_adm1))
  a1 <- list(); a2 <- list()
  for (k in seq_len(nrow(reg))) {
    r <- reg[k, ]
    for (i in seq_len(n_adm1)) {
      x0 <- r$lon0 + ((i - 1) %% ncol); y0 <- r$lat0 + ((i - 1) %/% ncol)
      gid1 <- sprintf("%s.%d_1", r$iso3, i); nm1 <- paste(r$name, .ADM1_SUFFIX[i])
      a1[[length(a1) + 1L]] <- data.frame(
        GID_0 = r$iso3, COUNTRY = r$name, GID_1 = gid1, NAME_1 = nm1, VARNAME_1 = NA_character_,
        NL_NAME_1 = NA_character_, TYPE_1 = "Province", ENGTYPE_1 = "Province", CC_1 = NA_character_,
        HASC_1 = NA_character_, ISO_1 = NA_character_, x0 = x0, y0 = y0, stringsAsFactors = FALSE)
      for (j in seq_len(n_adm2)) {
        a2[[length(a2) + 1L]] <- data.frame(
          GID_0 = r$iso3, COUNTRY = r$name, GID_1 = gid1, NAME_1 = nm1, NL_NAME_1 = NA_character_,
          GID_2 = sprintf("%s.%d.%d_1", r$iso3, i, j), NAME_2 = paste(nm1, .ADM2_SUFFIX[j]),
          VARNAME_2 = NA_character_, NL_NAME_2 = NA_character_, TYPE_2 = "District",
          ENGTYPE_2 = "District", CC_2 = NA_character_, HASC_2 = NA_character_,
          x0 = x0 + (j - 1) / n_adm2, y0 = y0, w = 1 / n_adm2, stringsAsFactors = FALSE)
      }
    }
  }
  a1 <- do.call(rbind, a1); a2 <- do.call(rbind, a2)
  g1 <- sf::st_sfc(lapply(seq_len(nrow(a1)), function(i) .msquare(a1$x0[i], a1$y0[i], 1, 1)), crs = 4326)
  g2 <- sf::st_sfc(lapply(seq_len(nrow(a2)), function(i) .msquare(a2$x0[i], a2$y0[i], a2$w[i], 1)), crs = 4326)
  list(adm1 = sf::st_sf(a1[, setdiff(names(a1), c("x0", "y0"))], geom = g1),
       adm2 = sf::st_sf(a2[, setdiff(names(a2), c("x0", "y0", "w"))], geom = g2))
}

#' synthUnits
#' The reporting units (location periods) of the synthetic database: the first
#' ADM1 of each country reports at ADM1, the others at ADM2, plus one national
#' unit per country. Location strings follow OutbreakExtractR
#' ("AFR::XAA::Astoria North::Astoria North A"; the "::" count gives the level).
#' @return data.frame(lpid, iso3, group, level, location, admin1, admin2, gid, pop)
synthUnits <- function(reg, gadm, seed = 1L) {
  set.seed(seed)
  rows <- list(); id <- 900000L
  for (k in seq_len(nrow(reg))) {
    r <- reg[k, ]
    a1 <- gadm$adm1[gadm$adm1$GID_0 == r$iso3, ]
    for (i in seq_len(nrow(a1))) {
      a2 <- gadm$adm2[gadm$adm2$GID_1 == a1$GID_1[i], ]
      pops <- round(stats::runif(nrow(a2), 5e4, 5e5))
      if (i == 1L) {
        id <- id + 1L
        rows[[length(rows) + 1L]] <- data.frame(
          lpid = as.character(id), iso3 = r$iso3, group = r$group, level = "ADM1",
          location = sprintf("AFR::%s::%s", r$iso3, a1$NAME_1[i]), admin1 = a1$NAME_1[i],
          admin2 = NA_character_, gid = a1$GID_1[i], pop = sum(pops), stringsAsFactors = FALSE)
      } else {
        for (j in seq_len(nrow(a2))) {
          id <- id + 1L
          rows[[length(rows) + 1L]] <- data.frame(
            lpid = as.character(id), iso3 = r$iso3, group = r$group, level = "ADM2",
            location = sprintf("AFR::%s::%s::%s", r$iso3, a1$NAME_1[i], a2$NAME_2[j]),
            admin1 = a1$NAME_1[i], admin2 = a2$NAME_2[j], gid = a2$GID_2[j], pop = pops[j],
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  u <- do.call(rbind, rows)
  nat <- do.call(rbind, lapply(split(u, u$iso3), function(d) {
    id <<- id + 1L
    data.frame(lpid = as.character(id), iso3 = d$iso3[1], group = d$group[1], level = "ADM0",
               location = sprintf("AFR::%s", d$iso3[1]), admin1 = NA_character_,
               admin2 = NA_character_, gid = d$iso3[1], pop = sum(d$pop), stringsAsFactors = FALSE)
  }))
  rbind(u, nat)
}

#' synthUnitGeometry
#' Geometry of each reporting unit (its GADM unit; the country union for ADM0).
synthUnitGeometry <- function(units, gadm) {
  g <- lapply(seq_len(nrow(units)), function(i) {
    u <- units[i, ]
    geo <- switch(u$level,
      ADM1 = sf::st_geometry(gadm$adm1[gadm$adm1$GID_1 == u$gid, ]),
      ADM2 = sf::st_geometry(gadm$adm2[gadm$adm2$GID_2 == u$gid, ]),
      ADM0 = sf::st_union(sf::st_geometry(gadm$adm1[gadm$adm1$GID_0 == u$iso3, ])))
    sf::st_cast(geo, "MULTIPOLYGON")[[1]]
  })
  sf::st_sf(location_period_id = units$lpid, geom = sf::st_sfc(g, crs = 4326))
}

# ---------------------------------------------------------------------------
# ENSO
# ---------------------------------------------------------------------------

#' synthNino34
#' Monthly Nino3.4 anomaly: an AR(2) with a ~4-year quasi-period, sd ~0.8.
#' @return data.frame(year, month, date, nino34)
synthNino34 <- function(start = "1950-01-01", end = "2026-09-01", seed = 7L) {
  set.seed(seed)
  dates <- seq(as.Date(start), as.Date(end), by = "month")
  n <- length(dates) + 24L
  x <- numeric(n); e <- stats::rnorm(n, sd = 0.28)
  for (t in 3:n) x[t] <- 1.6 * x[t - 1] - 0.68 * x[t - 2] + e[t]
  x <- x[-(1:24)]
  x <- x / stats::sd(x) * 0.8
  data.frame(year = as.integer(format(dates, "%Y")), month = as.integer(format(dates, "%m")),
             date = dates, nino34 = x)
}

#' synthOni
#' ONI = centred 3-month running mean of the Nino3.4 anomaly, in the shape of
#' utils/enso.R::fetchEnsoIndex(): year, month, enso_index, date, enso_phase
#' (CPC episode rule, ensoPhase(rule = "episode")).
synthOni <- function(n34) {
  x <- stats::filter(n34$nino34, rep(1 / 3, 3), sides = 2)
  ok <- is.finite(x)
  out <- data.frame(year = n34$year, month = n34$month, enso_index = round(as.numeric(x), 2),
                    date = n34$date)[ok, ]
  out$enso_phase <- ensoPhase(out$enso_index, rule = "episode")
  tibble::as_tibble(out)
}

#' synthEnsoEvents
#' Past ENSO events of the synthetic ONI in the shape of
#' configs/enso_events_v2.yml: the strongest El Nino and La Nina episodes
#' (CPC episode phase) peaking in boreal winter (Oct-Feb) inside [from, to],
#' each launched the month before its run of episode months starts, before that
#' winter (launch_rule crossing; the event attribution classifies the DJF season).
#' @return named list event id -> list(phase, label, peak, peak_season_year, launch, launch_rule)
synthEnsoEvents <- function(oni, from = as.Date("2013-01-01"), to = as.Date("2024-06-01"),
                            n = c(elnino = 2L, lanina = 1L)) {
  ph <- oni$enso_phase; ph[is.na(ph)] <- "neutral"
  r <- rle(ph); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  ev <- list()
  for (k in seq_along(r$values)) {
    v <- r$values[k]; if (!v %in% c("elnino", "nina")) next
    i <- starts[k]:ends[k]; pk <- i[which.max(abs(oni$enso_index[i]))]
    if (oni$date[pk] < from || oni$date[pk] > to || starts[k] < 2L) next
    yr <- as.integer(format(oni$date[pk], "%Y")); mo <- as.integer(format(oni$date[pk], "%m"))
    if (!mo %in% c(10:12, 1:2)) next                       # a boreal-winter (DJF season) peak
    psy <- if (mo >= 10L) yr else yr - 1L                   # the DJF season starting in December psy
    if (oni$date[starts[k] - 1L] >= as.Date(sprintf("%d-12-01", psy))) next   # launched before that winter
    phase <- if (v == "elnino") "elnino" else "lanina"
    ev[[length(ev) + 1L]] <- list(id = sprintf("%s_%d", phase, psy), phase = phase,
                                  strength = abs(oni$enso_index[pk]),
                                  label = sprintf("%s %d-%02d (synthetic)", if (phase == "elnino") "El Nino" else "La Nina",
                                                  psy, (psy + 1L) %% 100L),
                                  peak = format(oni$date[pk]), peak_season_year = psy,
                                  launch = format(oni$date[starts[k] - 1L]), launch_rule = "crossing")
  }
  out <- list()
  for (p in names(n)) {
    e <- Filter(function(x) x$phase == p, ev)
    e <- e[order(-vapply(e, `[[`, 0, "strength"))][seq_len(min(n[[p]], length(e)))]
    for (x in e) out[[x$id]] <- x[c("phase", "label", "peak", "peak_season_year", "launch", "launch_rule")]
  }
  out[order(vapply(out, `[[`, "", "peak"))]
}

# ---------------------------------------------------------------------------
# Epidemics
# ---------------------------------------------------------------------------

#' synthEpidemics
#' Monthly onsets, activity and cases per reporting unit (ADM1 / ADM2 units only).
#' @param units synthUnits()
#' @param reg synthRegistry()
#' @param oni synthOni()
#' @param months first-of-month Dates simulated (the database is cut later)
#' @return data.frame(lpid, month, active, onset, episode, cases)
synthEpidemics <- function(units, reg, oni, gadm,
                           months = seq(as.Date("2009-01-01"), as.Date("2026-08-01"), by = "month"),
                           a0 = -3.3, gamma = 0.8, rate = 6e-4, seed = 11L) {
  set.seed(seed)
  u <- units[units$level != "ADM0", ]
  n_u <- nrow(u); n_t <- length(months)
  r <- reg[match(u$iso3, reg$iso3), ]
  oni_at <- function(d) oni$enso_index[match(d, oni$date)]
  mo <- as.integer(format(months, "%m"))
  # neighbours: units of the same country, and of the other country of the group
  cent <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(synthUnitGeometry(u, gadm))))
  dkm <- as.matrix(stats::dist(cent)) * 111
  W <- (dkm > 0 & dkm <= 500) * 1
  W <- W / pmax(rowSums(W), 1)
  a_c <- stats::setNames(stats::rnorm(nrow(reg), 0, 0.3), reg$iso3)[u$iso3]
  a_u <- stats::rnorm(n_u, 0, 0.3)
  seas <- function(m, peak) cos(2 * pi * (m - peak) / 12)
  active <- onset <- episode <- matrix(0L, n_u, n_t); cases <- matrix(0, n_u, n_t)
  months_in <- integer(n_u); ep_id <- integer(n_u)
  for (t in seq_len(n_t)) {
    x <- vapply(seq_len(n_u), function(i) {
      d <- seq(months[t], by = paste(r$enso_lag[i], "months"), length.out = 2)[2]
      v <- oni_at(d); if (is.na(v)) 0 else v }, 0)
    nbr <- if (t > 1) as.numeric(W %*% active[, t - 1]) else rep(0, n_u)
    for (i in seq_len(n_u)) {
      s <- seas(mo[t], r$season_peak[i])
      if (t > 1 && active[i, t - 1] == 1L) {
        p_end <- stats::plogis(-1.2 + 0.35 * months_in[i] - 0.3 * r$enso_beta[i] * x[i])
        if (stats::runif(1) < p_end) { months_in[i] <- 0L } else {
          active[i, t] <- 1L; months_in[i] <- months_in[i] + 1L; episode[i, t] <- ep_id[i]
        }
      } else {
        h <- exp(a0 + a_c[i] + a_u[i] + 0.9 * s + r$enso_beta[i] * x[i] + gamma * nbr[i])
        if (stats::runif(1) < 1 - exp(-h)) {
          active[i, t] <- 1L; onset[i, t] <- 1L; months_in[i] <- 1L
          ep_id[i] <- ep_id[i] + 1L; episode[i, t] <- ep_id[i]
        }
      }
      if (active[i, t] == 1L) {
        mu <- u$pop[i] * rate * exp(0.5 * s + 0.2 * r$enso_beta[i] * x[i])
        cases[i, t] <- stats::rnbinom(1, size = 3, mu = mu) + 10   # >= 10: a recorded outbreak month
      }
    }
  }
  data.frame(lpid = rep(u$lpid, n_t), month = rep(months, each = n_u),
             active = as.vector(active), onset = as.vector(onset),
             episode = as.vector(episode), cases = as.vector(cases), stringsAsFactors = FALSE)
}

#' synthWeekly
#' Weekly OutbreakExtractR-style rows (Monday-start weeks, TR = TL + 6) of the
#' units, cases of a month split multinomially over the weeks starting in it;
#' outbreak_number numbers each unit's episodes (0 outside an outbreak). The
#' national units carry the weekly sum of their country's units; a third of the
#' units report confirmed cases (cCh) in their outbreak weeks.
synthWeekly <- function(units, epi, first = as.Date("2010-01-04"), last = as.Date("2024-12-30"), seed = 13L) {
  set.seed(seed)
  wk <- seq(first, last, by = "week")
  wmon <- as.Date(format(wk, "%Y-%m-01"))
  sub <- units[units$level != "ADM0", ]
  rows <- lapply(seq_len(nrow(sub)), function(i) {
    e <- epi[epi$lpid == sub$lpid[i], ]
    sCh <- numeric(length(wk)); ob <- integer(length(wk))
    for (k in which(e$cases > 0 & e$month >= min(wmon) & e$month <= max(wmon))) {
      idx <- which(wmon == e$month[k]); if (!length(idx)) next
      sCh[idx] <- as.vector(stats::rmultinom(1, e$cases[k], rep(1, length(idx))))
      ob[idx] <- e$episode[k]
    }
    data.frame(lpid = sub$lpid[i], TL = wk, TR = wk + 6, sCh = sCh, outbreak_number = ob)
  })
  w <- do.call(rbind, rows)
  # Confirmed cases: a third of the units test, confirming ~5 % of suspected
  # cases in their outbreak weeks (NA elsewhere, as an untested week).
  testing <- sub$lpid[stats::runif(nrow(sub)) < 1 / 3]
  w$cCh <- NA_real_
  k <- w$lpid %in% testing & w$outbreak_number > 0 & w$sCh > 0
  w$cCh[k] <- stats::rbinom(sum(k), w$sCh[k], 0.05)
  w$iso3 <- units$iso3[match(w$lpid, units$lpid)]
  nat <- stats::aggregate(sCh ~ iso3 + TL + TR, data = w, FUN = sum)
  nat$lpid <- units$lpid[units$level == "ADM0"][match(nat$iso3, units$iso3[units$level == "ADM0"])]
  nat$outbreak_number <- 0L; nat$cCh <- NA_real_
  rbind(w, nat[, names(w)])
}

#' synthOutbreakCsv
#' combined_outbreaks_cholera.csv in OutbreakExtractR's stage-2 shape (the
#' columns 00_ingest_outbreakextractr.R reads; the rest it fills with NA).
synthOutbreakCsv <- function(weekly, units) {
  u <- units[match(weekly$lpid, units$lpid), ]
  data.frame(
    location = ifelse(u$level == "ADM0", paste0(u$location, ".", u$iso3), u$location),
    TL = format(weekly$TL), TR = format(weekly$TR), sCh = weekly$sCh, cCh = weekly$cCh, deaths = NA,
    epiweek = format(weekly$TL, "%Y-W%V"), start_weekday = "Monday",
    spatial_scale = c(ADM0 = "country", ADM1 = "admin1", ADM2 = "admin2")[u$level],
    composite_loc = FALSE, date_range = 7, who_region = "AFR", country_iso3 = u$iso3,
    admin1 = u$admin1, admin2 = u$admin2, admin3 = NA, admin4 = NA, admin5 = NA, phantom = FALSE,
    original_location_name = u$location, observation_collection_id = "synthetic",
    location_period_id = u$lpid, pop = u$pop, adj_factors = 1,
    outbreak_number = weekly$outbreak_number, threshold = round(u$pop * 1e-4, 1),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Climate, population, WHO, static tables
# ---------------------------------------------------------------------------

#' synthClimateDaily
#' Daily climate of one unit in the Jahn et al. parquet schema (Date,
#' {pr,tas,tasmax,tasmin,hurs,huss}_{dynamic,static,simple}): a seasonal
#' rainfall cycle peaking in `rain_peak`, wetter in El Nino for a positive
#' `beta` (drier for a negative one), and a seasonal temperature cycle.
synthClimateDaily <- function(oni, rain_peak, beta, seed,
                              start = as.Date("1981-01-01"), end = as.Date("2024-12-31")) {
  set.seed(seed)
  d <- seq(start, end, by = "day")
  m <- as.integer(format(d, "%m"))
  o <- oni$enso_index[match(as.Date(format(d, "%Y-%m-01")), oni$date)]; o[is.na(o)] <- 0
  mu <- (0.5 + 6 * pmax(0, cos(2 * pi * (m - rain_peak) / 12))^2) * exp(0.25 * sign(beta) * o)
  pr <- ifelse(stats::runif(length(d)) < 0.55, stats::rgamma(length(d), shape = 0.8, scale = mu / 0.8 / 0.55), 0)
  tas <- 25 + 3 * cos(2 * pi * (m - rain_peak + 2) / 12) + stats::rnorm(length(d), 0, 1) + 0.2 * o
  vals <- list(pr = round(pr, 3), tas = round(tas, 2), tasmax = round(tas + 5 + stats::rnorm(length(d), 0, 0.5), 2),
               tasmin = round(tas - 5 + stats::rnorm(length(d), 0, 0.5), 2),
               hurs = round(pmin(99, 60 + 3 * pr), 1), huss = round(0.012 + 0.0002 * pr, 5))
  out <- data.frame(Date = d)
  for (w in c("dynamic", "static", "simple")) for (v in names(vals)) out[[paste0(v, "_", w)]] <- vals[[v]]
  out
}

#' synthWhoTables
#' The WHO national series of the synthetic countries, from the simulated
#' national weekly totals (so the database and WHO agree): annual report
#' 2010-2022, dashboard annual 2023-2025 and weekly 2023-01..2026-08, sitrep
#' cumulative rows. National counts are the database's times `who_factor`.
synthWhoTables <- function(units, epi, reg, who_factor = 1.15) {
  nat <- epi
  nat$iso3 <- units$iso3[match(nat$lpid, units$lpid)]
  nat$year <- as.integer(format(nat$month, "%Y"))
  nat$cases <- nat$cases * who_factor   # WHO's national count exceeds the subnational database
  yr <- stats::aggregate(cases ~ iso3 + year, data = nat, FUN = sum)
  yr$cases <- round(yr$cases)
  rep_yrs <- yr[yr$year >= 2010 & yr$year <= 2022, ]
  annual <- tibble::tibble(
    iso3 = rep_yrs$iso3, country_name = reg$name[match(rep_yrs$iso3, reg$iso3)], oc_who_region = "AFR",
    location = paste0("AFR::", rep_yrs$iso3),
    TL = as.Date(sprintf("%d-01-01", rep_yrs$year)), TR = as.Date(sprintf("%d-12-31", rep_yrs$year)),
    sCh = rep_yrs$cases, cCh = NA_real_, deaths = NA_real_,
    spatial_scale = factor("country", levels = c("country", "admin1", "admin2", "admin3", "admin4")),
    composite_loc = FALSE, date_range = as.difftime(365, units = "days"),
    temporal_scale = factor("yearly", levels = c("multiyear", "yearly", "monthly", "biweekly", "weekly", "daily", "other", "unknown")),
    start_weekday = factor("Monday", levels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"), ordered = TRUE),
    epiweek = NA_character_, who_region = "AFR", country = rep_yrs$iso3, admin1 = NA_character_,
    admin2 = NA_character_, admin3 = NA_character_, admin4 = NA_character_, admin5 = NA_character_,
    admin6 = NA_character_, location_period_id = NA_character_, observation_collection_id = "synthetic",
    original_location_name = reg$name[match(rep_yrs$iso3, reg$iso3)])
  # weekly dashboard: monthly national totals spread over the Mondays of the month
  mon <- stats::aggregate(cases ~ iso3 + month, data = nat, FUN = sum)
  mon <- mon[mon$month >= as.Date("2023-01-01"), ]
  wk <- do.call(rbind, lapply(seq_len(nrow(mon)), function(i) {
    days <- seq(mon$month[i], seq(mon$month[i], by = "month", length.out = 2)[2] - 1, by = "day")
    mondays <- days[format(days, "%u") == "1"]
    data.frame(iso_3_code = mon$iso3[i], date_wk = mondays,
               cases = as.integer(round(mon$cases[i] / length(mondays))))
  }))
  wk <- wk[order(wk$iso_3_code, wk$date_wk), ]
  wk$epiyr <- as.integer(format(wk$date_wk, "%G")); wk$epiwk <- as.integer(format(wk$date_wk, "%V"))
  wk$adm0_name <- toupper(reg$name[match(wk$iso_3_code, reg$iso3)])
  wk$deaths <- 0L; wk$who_region <- "AFRO"; wk$iso_2_code <- substr(wk$iso_3_code, 1, 2)
  wk$center_lon <- reg$lon0[match(wk$iso_3_code, reg$iso3)] + 1
  wk$center_lat <- reg$lat0[match(wk$iso_3_code, reg$iso3)] + 1
  wk$adm0_viz_n <- reg$name[match(wk$iso_3_code, reg$iso3)]
  wk$cum_case <- stats::ave(wk$cases, wk$iso_3_code, wk$epiyr, FUN = cumsum)
  wk$cum_death <- 0L; wk$ObjectId <- seq_len(nrow(wk))
  weekly <- wk[, c("adm0_name", "date_wk", "cases", "deaths", "who_region", "iso_2_code", "iso_3_code",
                   "center_lon", "center_lat", "adm0_viz_n", "cum_case", "cum_death", "epiwk", "epiyr", "ObjectId")]
  da <- stats::aggregate(cases ~ iso_3_code + epiyr, data = wk, FUN = sum)
  dash_annual <- data.frame(
    report_year = da$epiyr, adm0_name = toupper(reg$name[match(da$iso_3_code, reg$iso3)]),
    who_region = "AFRO", iso_3_code = da$iso_3_code,
    first_epiwk = as.Date(sprintf("%d-01-01", da$epiyr)), last_epiwk = as.Date(sprintf("%d-12-31", da$epiyr)),
    case_total = as.integer(da$cases), death_total = 0L)
  dash_annual$last_epiwk[dash_annual$report_year == max(dash_annual$report_year)] <- max(wk$date_wk)
  last28 <- do.call(rbind, lapply(split(wk, wk$iso_3_code), function(d)
    data.frame(iso3 = d$iso_3_code[1], table_as_of = format(max(d$date_wk) + 6),
               cum_cases = as.integer(sum(d$cases[d$epiyr == max(d$epiyr)])),
               last28_cases = as.integer(sum(utils::tail(d$cases, 4))), report_id = "synthetic",
               source_url = NA_character_, flag = NA_character_)))
  list(annual = annual, dashboard_annual = dash_annual, dashboard_weekly = weekly, sitrep = last28)
}

#' synthSeasonality
#' best_seasonality_estimates.csv rows of every ADM1 (12 months each, month
#' order): the normalised seasonal profile of the unit's group.
synthSeasonality <- function(gadm, reg) {
  a1 <- sf::st_drop_geometry(gadm$adm1)
  do.call(rbind, lapply(seq_len(nrow(a1)), function(i) {
    pk <- reg$season_peak[match(a1$GID_0[i], reg$iso3)]
    m <- 1:12; v <- (cos(2 * pi * (m - pk) / 12) + 1) / 2
    data.frame(country = a1$GID_0[i], gid = a1$GID_1[i], name = a1$NAME_1[i], run_level = 1L,
               model = "offset", month = m, mean = round(v, 4),
               q025 = round(pmax(0, v - 0.15), 4), q975 = round(pmin(1, v + 0.15), 4))
  }))
}

#' synthWorldPopRaster
#' A WorldPop-like population raster of one country and year at 0.02 degrees:
#' each ADM2 strip carries its unit's population (grown 2.5 % a year from 2020),
#' spread evenly over its cells.
synthWorldPopRaster <- function(units, gadm, iso3, year, res = 0.02) {
  a2 <- gadm$adm2[gadm$adm2$GID_0 == iso3, ]
  bb <- sf::st_bbox(a2)
  r <- terra::rast(xmin = bb[["xmin"]], xmax = bb[["xmax"]], ymin = bb[["ymin"]], ymax = bb[["ymax"]],
                   resolution = res, crs = "EPSG:4326")
  pop2 <- vapply(seq_len(nrow(a2)), function(j) {
    g2 <- a2$GID_2[j]; g1 <- a2$GID_1[j]
    p <- units$pop[units$gid == g2]
    if (!length(p)) {                                  # an ADM1-level unit: split evenly
      p1 <- units$pop[units$gid == g1]; p <- p1 / sum(gadm$adm2$GID_1 == g1)
    }
    p[1] }, 0)
  # ADM2 index of every cell, then each strip's population over its cell count
  a2$k <- seq_len(nrow(a2))
  cells <- terra::values(terra::rasterize(terra::vect(a2), r, field = "k"))[, 1]
  n_k <- tabulate(cells, nbins = nrow(a2))
  terra::values(r) <- pop2[cells] * 1.025^(year - 2020) / n_k[cells]
  r
}
