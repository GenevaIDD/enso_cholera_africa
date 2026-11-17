# enso.R — ENSO index data, phases, lags and forecast / counterfactual ENSO
#   trajectories (CFSv2, C3S, AR(1) / analog continuation, Lenssen piControl).
#
# Provides:
#   fetch_enso_index(), join_enso_lags(), add_enso_leads_from_lag0()   ONI
#   classify_enso_phase(), enso_phase_definition(), enso_year_states()
#   enso_group_colours(), group_colour_map()     ENSO-group palette
#   fetch_enso_cfsv2_ensemble(), fit_cfsv2_calibration(),
#   apply_cfsv2_calibration()                    CFSv2 forecasts
#   sample_enso_forecast_trajectories(), sample_enso_cfsv2_trajectories()
#   continue_enso_ar1(), continue_enso_analog()  continuation past CFSv2
#   sample_enso_analog_scenarios(), sample_enso_lenssen_scenarios()
#   ensure_lenssen_counterfactual_cache(), lenssen_cache_status()
#   enso_seam_month(), enso_event_launch(), latest_oni_month()
#
# Details:
#   One ENSO phase rule for the whole pipeline (classify_enso_phase()).
#   The Lenssen cache helpers back 03_prepare_lenssen_analogs.R.

# ENSO region palette -------------------------------------------------------

#' group_colour_map / enso_group_colours
#'
#' `group_colour_map(meta)`: the colour vector of a registry's sidecar
#' (utils/config.R::load_grouping_meta()$colour, keyed by group label).
#' `enso_group_colours()`: the canonical palette for the five ENSO teleconnection groups
#' (`analysis/groupings/enso_grouping.csv`: enso_east, enso_horn, enso_west,
#' enso_south, enso_central), shared by every figure that colours by ENSO
#' group so the same region always gets the same colour across the
#' manuscript. Single source of truth.
#'
#' Deliberately avoids red and blue: those two colours are reserved
#' elsewhere in the manuscript figures for the ENSO phases (El Niño red, La
#' Niña blue) and the diverging blue-white-red scales of the effect maps
#' (`figures_cb.R::ps_diverging_fill_scale()`), so an ENSO-region palette
#' that also used red/blue would visually collide with those unrelated
#' encodings whenever both appear in the same figure. Uses a red/blue-free slice of the
#' ColorBrewer "Dark2" qualitative palette instead.
#'
#' The values live in analysis/groupings/enso_grouping.yml;
#' tests/test_grouping_meta.R pins them.
#'
#' @return Named character vector of hex colours, keyed by group label.
group_colour_map <- function(meta) meta$colour

# Call-time dependency on utils/config.R::load_grouping_meta().
enso_group_colours <- function() {
  group_colour_map(load_grouping_meta("enso_grouping"))
}

# ENSO data ---------------------------------------------------------------

# ENSO phase: one function, two rules -------------------------------------
#
# Every ENSO phase label in the pipeline comes from classify_enso_phase(). The RULE
# (instantaneous threshold vs the NOAA CPC episode rule) and the STATISTIC it
# is applied to (monthly ONI, the DJF ONI of a season-year, a trajectory's DJF
# mean, ...) are separate: the statistic is the caller's input, recorded with
# enso_phase_definition(). Only the episode rule is "NOAA CPC" / "official"; the
# threshold rule is a "threshold phase" (no persistence requirement).

#' ENSO phase cut-offs: La Niña at or below `nina`, El Niño at or above
#' `elnino` (both inclusive; the CPC operational +/-0.5 degC).
ENSO_PHASE_THRESHOLDS <- c(nina = -0.5, elnino = 0.5)
#' CPC episode rule: minimum run of consecutive overlapping 3-month seasons.
ENSO_EPISODE_MIN_RUN <- 5L
#' Phase codes classify_enso_phase() returns, in display order.
ENSO_PHASE_CODES <- c("nina", "neutral", "elnino")
#' Display labels of the phase codes.
ENSO_PHASE_LABELS <- c(
  nina = "La Niña",
  neutral = "Neutral",
  elnino = "El Niño"
)

#' classify_enso_phase
#'
#' ENSO phase of a numeric ENSO-index statistic.
#'
#' \code{rule = "threshold"}: each value on its own, \code{"nina"} at or below
#' \code{thresholds[["nina"]]}, \code{"elnino"} at or above
#' \code{thresholds[["elnino"]]}, else \code{"neutral"}.
#'
#' \code{rule = "episode"}: the NOAA CPC operational episode rule. A value at or
#' beyond a threshold counts only when it is part of a run of at least
#' \code{min_run} consecutive same-sign values (consecutive overlapping 3-month
#' seasons when \code{x} is a monthly ONI series in date order).
#'
#' @param x          numeric vector (chronological order for the episode rule)
#' @param rule       "threshold" | "episode"
#' @param thresholds named numeric c(nina =, elnino =), nina < elnino
#' @param min_run    episode rule: minimum run length
#' @return character vector of ENSO_PHASE_CODES, NA where \code{x} is NA
classify_enso_phase <- function(
  x,
  rule = c("threshold", "episode"),
  thresholds = ENSO_PHASE_THRESHOLDS,
  min_run = ENSO_EPISODE_MIN_RUN
) {
  rule <- match.arg(rule)
  if (!all(c("nina", "elnino") %in% names(thresholds))) {
    stop(
      "classify_enso_phase: `thresholds` must be a named numeric with elements 'nina' and 'elnino'.",
      call. = FALSE
    )
  }
  lo <- thresholds[["nina"]]
  hi <- thresholds[["elnino"]]
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    stop(
      "classify_enso_phase: `thresholds` must be finite with nina < elnino.",
      call. = FALSE
    )
  }
  n <- length(x)
  if (identical(rule, "threshold")) {
    out <- rep(NA_character_, n)
    ok <- !is.na(x)
    out[ok & x <= lo] <- "nina"
    out[ok & x >= hi] <- "elnino"
    out[ok & is.na(out)] <- "neutral"
    return(out)
  }
  phase <- rep("neutral", n)
  phase[is.na(x)] <- NA_character_
  mark_runs <- function(flag, label) {
    flag[is.na(flag)] <- FALSE
    rl <- rle(flag)
    ends <- cumsum(rl$lengths)
    starts <- ends - rl$lengths + 1L
    hits <- which(rl$values & rl$lengths >= min_run)
    for (i in hits) {
      phase[starts[i]:ends[i]] <<- label
    }
  }
  mark_runs(x >= hi, "elnino")
  mark_runs(x <= lo, "nina")
  phase
}

#' The statistics classify_enso_phase() is applied to in this pipeline.
ENSO_PHASE_STATISTICS <- c(
  monthly_oni = "monthly ONI (3-month running mean of Niño 3.4)",
  monthly_oni_cfsv2_ext = "monthly ONI, extended by the CFSv2 ensemble mean past the observed record",
  djf_oni = "ONI of the DJF season (the January value)",
  djf_mean_oni = "mean of the monthly ONI over December-February (a 5-month triangular mean of Niño 3.4)",
  event_run = "monthly ONI, contiguous run back from the event peak",
  year_mean_oni = "mean monthly ONI of a calendar year (January-December)",
  year_lag_mean_oni = "mean monthly ONI of a calendar year shifted by a deployed ENSO lag"
)

#' enso_phase_definition
#'
#' A self-describing record of how a phase label was made: the rule, the
#' statistic, the thresholds and a display label. Only the episode rule is
#' labelled "NOAA CPC episode"; the threshold rule is "threshold phase".
#'
#' @param rule      "threshold" | "episode"
#' @param statistic one of names(ENSO_PHASE_STATISTICS)
#' @return list(rule, statistic, statistic_label, thresholds, min_run, label)
enso_phase_definition <- function(
  rule = c("threshold", "episode"),
  statistic = names(ENSO_PHASE_STATISTICS),
  thresholds = ENSO_PHASE_THRESHOLDS,
  min_run = ENSO_EPISODE_MIN_RUN
) {
  rule <- match.arg(rule)
  statistic <- match.arg(statistic)
  cut <- sprintf(
    "%s / %s",
    format(thresholds[["nina"]]),
    format(thresholds[["elnino"]])
  )
  label <- if (identical(rule, "episode")) {
    sprintf(
      "NOAA CPC episode (%s at or beyond %s for >= %d consecutive overlapping seasons)",
      ENSO_PHASE_STATISTICS[[statistic]],
      cut,
      as.integer(min_run)
    )
  } else {
    sprintf(
      "threshold phase (%s at or beyond %s, no persistence)",
      ENSO_PHASE_STATISTICS[[statistic]],
      cut
    )
  }
  list(
    rule = rule,
    statistic = statistic,
    statistic_label = unname(ENSO_PHASE_STATISTICS[[statistic]]),
    thresholds = as.list(thresholds),
    min_run = if (identical(rule, "episode")) {
      as.integer(min_run)
    } else {
      NA_integer_
    },
    label = label
  )
}

#' ENSO state of a calendar year (09 report_who_enso_ratio): the
#' four statistics a calendar-year case total can be labelled by, each cut by
#' classify_enso_phase(rule = "threshold"). Names are the config values of
#' results.derivation.who_enso_ratio.year_definitions; values are
#' names(ENSO_PHASE_STATISTICS).
#'   annual_mean  mean monthly ONI, January-December of year y
#'   lag_occ      the same twelve months shifted by the group's deployed
#'                occurrence lag L (L = -4: September y-1 to August y)
#'   djf_start    the DJF season centred on January y
#'   djf_end      the DJF season centred on January y+1 (the rule of 09's
#'                report_dataset v4, which labels a cholera season-year)
ENSO_YEAR_DEFINITIONS <- c(
  annual_mean = "year_mean_oni",
  lag_occ = "year_lag_mean_oni",
  djf_start = "djf_oni",
  djf_end = "djf_oni"
)

#' enso_year_states
#'
#' The ENSO statistics and threshold phases of calendar years under the
#' ENSO_YEAR_DEFINITIONS. A mean over a window with any month past the record
#' uses the months available (`n_months_*` says how many); a DJF past the record
#' is NA, and so is its phase.
#'
#' @param oni        monthly ONI (year, month, enso_index), as fetch_enso_index()
#' @param years      integer calendar years
#' @param lag_months the deployed occurrence lag L of `lag_occ` (months; the
#'                   window is months L .. L + 11 relative to January y)
#' @param thresholds classify_enso_phase() thresholds
#' @return tibble(year, oni_annual_mean, n_months_annual_mean, oni_lag_occ,
#'   n_months_lag_occ, oni_djf_start, oni_djf_end, and state_<definition> for
#'   every definition: ENSO_PHASE_CODES, NA without data)
enso_year_states <- function(
  oni,
  years,
  lag_months = 0L,
  thresholds = ENSO_PHASE_THRESHOLDS
) {
  o <- tibble::as_tibble(oni) |>
    dplyr::transmute(
      year = as.integer(year),
      month = as.integer(month),
      v = as.numeric(enso_index)
    ) |>
    dplyr::filter(!is.na(v)) |>
    dplyr::mutate(mi = year * 12L + month - 1L)
  if (anyDuplicated(o$mi)) {
    stop("enso_year_states: duplicate months in the ONI series.", call. = FALSE)
  }
  years <- as.integer(years)
  lag_months <- as.integer(lag_months)
  win <- function(y, shift) {
    lo <- y * 12L + shift
    v <- o$v[o$mi >= lo & o$mi <= lo + 11L]
    c(if (length(v)) mean(v) else NA_real_, length(v))
  }
  jan <- function(y) {
    v <- o$v[o$year == y & o$month == 1L]
    if (length(v)) v else NA_real_
  }
  am <- vapply(years, win, numeric(2), shift = 0L)
  lg <- vapply(years, win, numeric(2), shift = lag_months)
  out <- tibble::tibble(
    year = years,
    oni_annual_mean = am[1, ],
    n_months_annual_mean = as.integer(am[2, ]),
    oni_lag_occ = lg[1, ],
    n_months_lag_occ = as.integer(lg[2, ]),
    oni_djf_start = vapply(years, jan, numeric(1)),
    oni_djf_end = vapply(years + 1L, jan, numeric(1))
  )
  for (d in names(ENSO_YEAR_DEFINITIONS)) {
    out[[paste0("state_", d)]] <- classify_enso_phase(
      out[[paste0("oni_", d)]],
      "threshold",
      thresholds
    )
  }
  out
}

#' Classify an ONI/RONI series into official NOAA CPC ENSO episodes.
#'
#' Thin wrapper of \code{classify_enso_phase(index, "episode")}, kept for its callers
#' (fetch_enso_index()). An El Niño (La Niña) episode is
#' declared for every season where the index is at or above +threshold (at or
#' below -threshold) \emph{and} that crossing is part of a run of at least
#' \code{min_run} consecutive overlapping 3-month seasons. Each row of
#' \code{fetch_enso_index(source="oni")} already represents one such overlapping
#' season, so the runs are computed directly on the monthly-indexed series.
#'
#' @param index     Numeric vector of ONI/RONI values, in chronological order.
#' @param threshold Threshold magnitude (default ENSO_PHASE_THRESHOLDS: 0.5).
#' @param min_run   Minimum run length (default ENSO_EPISODE_MIN_RUN: 5).
#' @return Character vector, same length as \code{index}, one of
#'   \code{"elnino"}, \code{"nina"}, \code{"neutral"} (or \code{NA} where
#'   \code{index} is \code{NA}).
.classify_enso_episodes <- function(
  index,
  threshold = ENSO_PHASE_THRESHOLDS[["elnino"]],
  min_run = ENSO_EPISODE_MIN_RUN
) {
  classify_enso_phase(
    index,
    "episode",
    thresholds = c(nina = -threshold, elnino = threshold),
    min_run = min_run
  )
}

#' fetch_enso_index
#'
#' Download the Oceanic Niño Index (ONI) from NOAA CPC and cache locally.
#' ONI is the 3-month running mean of ERSSTv6 SST anomalies in the
#' Niño 3.4 region (5°N-5°S, 120°-170°W) and is the official NOAA ENSO
#' event definition.
#'
#' @param source    Index to download. One of:
#'   "oni"    — Oceanic Niño Index (primary; 3-month running mean of ERSSTv6
#'              SST anomalies in Niño 3.4, fixed 30-yr baseline).
#'   "roni"   — Relative ONI (official NOAA index as of Feb 2026; same as ONI
#'              but with tropical-mean SST anomaly subtracted, removing the
#'              global-warming trend; preferred for projection work).
#'   "nino34" — Raw monthly Niño 3.4 SST anomaly, ERSSTv6, centered base
#'              periods (CPC `data/indices/detrend.nino34.ascii.txt`; its
#'              3-month mean is ONI). Seam/history series of the CFSv2 sampler
#'              and the "ersst" target of the CFSv2 calibration.
#'   "nino34_oisst" — Monthly OISST v2.1 Niño 3.4 anomaly (fixed 1991-2020
#'              base; CPC `sstoi.indices`). The observational reference of
#'              CPC's PDF-corrected CFSv2 Niño 3.4 forecast, used as one
#'              target of the hindcast calibration
#'              (`03_prepare_enso_cfsv2_calibration.R`).
#' @param cache_dir Directory to store the cached RDS and metadata sidecar.
#' @param refresh   If TRUE, re-download even when a cache file exists.
#' @return tibble with columns: year (int), month (int), enso_index (dbl),
#'   date (Date, first day of the centre month), enso_phase (chr; one of
#'   "elnino"/"nina"/"neutral" for \code{source \%in\% c("oni", "roni")} per
#'   the official NOAA CPC definition -- ONI at or beyond \eqn{\pm}0.5 for at
#'   least 5 consecutive overlapping 3-month seasons; \code{NA} for
#'   \code{source = "nino34"}, an unsmoothed sensitivity series the official
#'   definition does not apply to).
fetch_enso_index <- function(
  source = "oni",
  cache_dir = "generated_data/",
  refresh = FALSE
) {
  source <- match.arg(source, c("oni", "roni", "nino34", "nino34_oisst"))

  urls <- c(
    oni = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
    roni = "https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt",
    # CPC's current ERSSTv6 file. The former URL
    # (products/analysis_monitoring/ensostuff/detrend.nino34.ascii.txt) is a copy
    # frozen on 5 Aug 2026: older ERSST version, ends Jun 2026, differs by up to
    # 0.71 from v6.
    nino34 = "https://www.cpc.ncep.noaa.gov/data/indices/detrend.nino34.ascii.txt",
    nino34_oisst = "https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices"
  )
  url <- urls[[source]]
  cache_file <- file.path(cache_dir, sprintf("enso_%s_cache.rds", source))
  meta_file <- file.path(
    cache_dir,
    sprintf("enso_%s_download_meta.yml", source)
  )

  if (!refresh && file.exists(cache_file)) {
    message("fetch_enso_index: loading cached ", source, " from ", cache_file)
    return(readRDS(cache_file))
  }

  message("fetch_enso_index: downloading ", source, " from\n  ", url)
  raw <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) {
      stop("fetch_enso_index: download failed — ", conditionMessage(e))
    }
  )

  if (source %in% c("oni", "roni")) {
    # ONI:  SEAS YR TOTAL ANOM  (4 columns; anom is last)
    # RONI: SEAS YR ANOM        (3 columns; anom is last)
    # SEAS is a 3-letter season code; YR is the year of the centre month.
    # Centre-month mapping (DJF centre = January of YR, etc.):
    seas_to_month <- c(
      DJF = 1L,
      JFM = 2L,
      FMA = 3L,
      MAM = 4L,
      AMJ = 5L,
      MJJ = 6L,
      JJA = 7L,
      JAS = 8L,
      ASO = 9L,
      SON = 10L,
      OND = 11L,
      NDJ = 12L
    )
    data_lines <- raw[!grepl("^\\s*(SEAS|#|$)", raw)]
    # Detect column count from first data line
    ncols <- length(strsplit(trimws(data_lines[1]), "\\s+")[[1]])
    col_names <- if (ncols == 4L) {
      c("seas", "year", "total", "anom")
    } else {
      c("seas", "year", "anom")
    }
    dat <- utils::read.table(
      text = paste(data_lines, collapse = "\n"),
      header = FALSE,
      col.names = col_names,
      stringsAsFactors = FALSE
    )
    out <- dplyr::tibble(
      year = as.integer(dat$year),
      month = as.integer(seas_to_month[dat$seas]),
      enso_index = suppressWarnings(as.numeric(dat$anom))
    )
  } else if (source == "nino34_oisst") {
    # Format: YR MON NINO1+2 ANOM NINO3 ANOM NINO4 ANOM NINO3.4 ANOM
    data_lines <- raw[!grepl("^\\s*(YR|#|$)", raw)]
    dat <- utils::read.table(
      text = paste(data_lines, collapse = "\n"),
      header = FALSE
    )
    if (ncol(dat) != 10L) {
      stop(
        "fetch_enso_index: unexpected sstoi.indices layout (",
        ncol(dat),
        " columns)"
      )
    }
    out <- dplyr::tibble(
      year = as.integer(dat[[1]]),
      month = as.integer(dat[[2]]),
      enso_index = suppressWarnings(as.numeric(dat[[10]]))
    )
  } else {
    # Format: YR MON TOTAL ClimAdjust ANOM
    data_lines <- raw[!grepl("^\\s*(YR|#|$)", raw)]
    dat <- utils::read.table(
      text = paste(data_lines, collapse = "\n"),
      header = FALSE,
      col.names = c("year", "month", "total", "climadj", "anom"),
      stringsAsFactors = FALSE
    )
    out <- dplyr::tibble(
      year = as.integer(dat$year),
      month = as.integer(dat$month),
      enso_index = suppressWarnings(as.numeric(dat$anom))
    )
  }

  out <- out %>%
    dplyr::filter(!is.na(enso_index)) %>%
    dplyr::mutate(date = as.Date(sprintf("%04d-%02d-01", year, month))) %>%
    dplyr::arrange(date)

  out$enso_phase <- if (source %in% c("oni", "roni")) {
    .classify_enso_episodes(out$enso_index)
  } else {
    NA_character_
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  saveRDS(out, cache_file)
  yaml::write_yaml(
    list(
      source = source,
      url = url,
      downloaded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      nrow = nrow(out),
      year_range = c(min(out$year), max(out$year))
    ),
    meta_file
  )
  message("fetch_enso_index: cached ", nrow(out), " records to ", cache_file)
  out
}


#' fetch_enso_cfsv2_ensemble
#'
#' Fetch NOAA CPC CFSv2 seasonal-forecast Niño 3.4 ensemble member trajectories
#' and return them as a long tibble.  Mirrors `fetch_enso_index()`: download →
#' cache RDS + `_download_meta.yml`.  Provides the raw per-member forecast that
#' `sample_enso_cfsv2_trajectories()` stitches to the historical record.
#'
#' Access order (first hit wins):
#'   1. `<cache_dir>/enso_cfsv2_cache.rds`        — prior download cache
#'   2. `<cache_dir>/enso_cfsv2_companion.{rds,csv}` — pre-extracted members
#'      (NetCDF-free path for machines without libnetcdf / ncdf4)
#'   3. Live download of `url` (NetCDF; requires the `ncdf4` package)
#'
#' The CFSv2 file carries monthly Niño 3.4 SST *anomaly* (vs the 1991–2020 base)
#' over a degenerate 1×1×1 spatial box with dims (lon, lat, lev, `time`, `ens`):
#' a single forecast issuance where `time` = target months and `ens` = ensemble
#' members. `init_date` is set to the earliest target month (implicit issuance);
#' values outside the returned span are handled downstream (0-fill).
#'
#' @param url        CFSv2 Niño 3.4 monthly NetCDF URL.
#' @param cache_dir  Directory for cache + companion + meta files.
#' @param init       "latest" (default) or "YYYY-MM" to pin a start month.
#' @param refresh    If TRUE, ignore the cache and re-download.
#' @return tibble(init_date, target_date, lead, member, nino34_anomaly) filtered
#'         to the selected init.
fetch_enso_cfsv2_ensemble <- function(
  url = "https://www.cpc.ncep.noaa.gov/products/CFSv2/dataInd3/nino34Mon.nc",
  cache_dir = "generated_data/",
  init = "latest",
  refresh = FALSE
) {
  cache_file <- file.path(cache_dir, "enso_cfsv2_cache.rds")
  meta_file <- file.path(cache_dir, "enso_cfsv2_download_meta.yml")
  comp_rds <- file.path(cache_dir, "enso_cfsv2_companion.rds")
  comp_csv <- file.path(cache_dir, "enso_cfsv2_companion.csv")
  required <- c("init_date", "target_date", "lead", "member", "nino34_anomaly")

  # Filter a full (all-inits) table down to the requested init.
  select_init <- function(full) {
    full <- tibble::as_tibble(full)
    missing <- setdiff(required, names(full))
    if (length(missing) > 0) {
      stop(
        "fetch_enso_cfsv2_ensemble: CFSv2 data missing columns: ",
        paste(missing, collapse = ", ")
      )
    }
    full$init_date <- as.Date(full$init_date)
    full$target_date <- as.Date(full$target_date)
    avail <- sort(unique(full$init_date))
    keep <- if (identical(init, "latest")) {
      max(avail)
    } else {
      want <- as.Date(sprintf("%s-01", init))
      if (want %in% avail) want else avail[which.min(abs(avail - want))]
    }
    message(
      "fetch_enso_cfsv2_ensemble: using init ",
      keep,
      " (",
      length(unique(full$member[full$init_date == keep])),
      " members)"
    )
    dplyr::filter(full, init_date == keep)
  }

  # 1. Download cache
  if (!refresh && file.exists(cache_file)) {
    message("fetch_enso_cfsv2_ensemble: loading cached CFSv2 from ", cache_file)
    return(select_init(readRDS(cache_file)))
  }

  # 2. Companion file (NetCDF-free override)
  if (file.exists(comp_rds)) {
    message("fetch_enso_cfsv2_ensemble: reading companion RDS: ", comp_rds)
    return(select_init(readRDS(comp_rds)))
  }
  if (file.exists(comp_csv)) {
    message("fetch_enso_cfsv2_ensemble: reading companion CSV: ", comp_csv)
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("fetch_enso_cfsv2_ensemble: readr required to read companion CSV.")
    }
    return(select_init(readr::read_csv(comp_csv, show_col_types = FALSE)))
  }

  # 3. Live NetCDF download
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop(
      "fetch_enso_cfsv2_ensemble: `ncdf4` package required to read the CFSv2 ",
      "NetCDF; install.packages('ncdf4') (needs system libnetcdf) or drop a ",
      "companion file at ",
      comp_rds,
      " / ",
      comp_csv,
      "."
    )
  }

  message("fetch_enso_cfsv2_ensemble: downloading CFSv2 NetCDF from\n  ", url)
  tmp <- tempfile(fileext = ".nc")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch(
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
    error = function(e) {
      stop("fetch_enso_cfsv2_ensemble: download failed — ", conditionMessage(e))
    }
  )
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("fetch_enso_cfsv2_ensemble: downloaded file is empty: ", url)
  }

  nc <- ncdf4::nc_open(tmp)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # The CFSv2 Niño 3.4 file stores anomaly over a degenerate 1×1×1 spatial box
  # with dims (lon, lat, lev, time, ens): a SINGLE forecast issuance where
  # `time` = target months and `ens` = ensemble members. Pick the anomaly
  # variable (most dims); ncvar_get drops the length-1 spatial dims.
  vsel <- names(nc$var)[which.max(vapply(
    nc$var,
    function(v) v$ndims,
    integer(1)
  ))]
  vobj <- nc$var[[vsel]]
  arr <- ncdf4::ncvar_get(nc, vsel) # length-1 dims dropped

  all_names <- vapply(vobj$dim, function(d) d$name, character(1))
  all_lens <- vapply(vobj$dim, function(d) d$len, integer(1))
  all_vals <- lapply(vobj$dim, function(d) d$vals)
  keep <- which(all_lens > 1L) # non-degenerate dims, file order
  nd_names <- all_names[keep]
  if (length(nd_names) < 2L) {
    stop(
      "fetch_enso_cfsv2_ensemble: expected >=2 non-degenerate dims (time, ens); got: ",
      paste(sprintf("%s:%d", all_names, all_lens), collapse = ", ")
    )
  }

  find_nd <- function(patt) {
    h <- grep(patt, nd_names, ignore.case = TRUE)
    if (length(h)) h[1] else NA_integer_
  }
  ti <- find_nd("time|target|^t$")
  mi <- find_nd("ens|member|^m$")
  if (anyNA(c(ti, mi))) {
    stop(
      "fetch_enso_cfsv2_ensemble: could not resolve time/member dims from: ",
      paste(nd_names, collapse = ", ")
    )
  }

  # Decode target-time coordinate → month-start Date.
  t_units <- tryCatch(
    ncdf4::ncatt_get(nc, nd_names[ti], "units")$value,
    error = function(e) ""
  )
  decode_time <- function(vals, units) {
    if (!is.null(units) && nzchar(units)) {
      m <- regmatches(
        units,
        regexec(
          "(day|month|year)s?\\s+since\\s+([0-9]{1,4}-[0-9]{1,2}-[0-9]{1,2})",
          units
        )
      )[[1]]
      if (length(m) == 3L) {
        origin <- as.Date(m[3])
        d <- switch(
          m[2],
          day = origin + as.integer(round(vals)),
          month = origin %m+% months(as.integer(round(vals))),
          year = origin %m+% months(as.integer(round(vals * 12)))
        )
        return(lubridate::floor_date(d, "month"))
      }
    }
    v <- as.numeric(vals)
    if (all(v > 100000)) {
      return(as.Date(sprintf("%04d-%02d-01", v %/% 100, v %% 100)))
    }
    lubridate::floor_date(as.Date(v, origin = "1970-01-01"), "month")
  }
  target_dates <- decode_time(all_vals[[keep[ti]]], t_units)
  member_ids <- all_vals[[keep[mi]]]
  if (length(nd_names) != 2L) {
    stop(
      "fetch_enso_cfsv2_ensemble: expected exactly (time, member) after dropping ",
      "spatial dims; got: ",
      paste(nd_names, collapse = ", ")
    )
  }

  # Orient the anomaly matrix as [time × member].
  if (ti == 2L) {
    arr <- t(arr)
  }

  # The file is a lead-triangle: forecast members carry the FORECAST tail months
  # (present at the last target month), and a separate "observed/analysis" column
  # carries the recent lead-in months (present early, NA in the tail). Rebuild one
  # full trajectory per forecast member by concatenating the shared observed
  # lead-in with that member's forecast tail — the natural ensemble shape.
  last_row <- arr[which.max(target_dates), ]
  fc_cols <- which(!is.na(last_row)) # 40 forecast members
  obs_cols <- setdiff(seq_len(ncol(arr)), fc_cols) # observed/analysis col(s)
  obs_series <- if (length(obs_cols) > 0L) {
    arr[, obs_cols[1]]
  } else {
    rep(NA_real_, nrow(arr))
  }

  W <- arr[, fc_cols, drop = FALSE] # [n_time × n_fc]
  na_mask <- is.na(W)
  if (any(na_mask)) {
    obs_fill <- matrix(obs_series, nrow = nrow(W), ncol = ncol(W))
    W[na_mask] <- obs_fill[na_mask] # fill lead-in from observed
  }

  init_date <- min(target_dates, na.rm = TRUE) # single implicit issuance
  idx <- expand.grid(t = seq_len(nrow(W)), m = seq_len(ncol(W)))
  full <- dplyr::tibble(
    target_date = target_dates[idx$t],
    member = member_ids[fc_cols][idx$m],
    nino34_anomaly = as.vector(W)
  ) %>%
    dplyr::mutate(
      init_date = init_date,
      lead = as.integer(round(
        as.numeric(
          difftime(target_date, init_date, units = "days")
        ) /
          30.4
      ))
    ) %>%
    dplyr::filter(!is.na(nino34_anomaly)) %>%
    dplyr::select(init_date, target_date, lead, member, nino34_anomaly)

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  saveRDS(full, cache_file)
  yaml::write_yaml(
    list(
      url = url,
      downloaded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      variable = vsel,
      dims = as.list(stats::setNames(dim(arr), nd_names)),
      init = as.character(init_date),
      target_range = as.character(c(
        min(full$target_date),
        max(full$target_date)
      )),
      n_members = length(unique(full$member)),
      nrow = nrow(full)
    ),
    meta_file
  )
  message(
    "fetch_enso_cfsv2_ensemble: cached ",
    nrow(full),
    " records to ",
    cache_file
  )
  select_init(full)
}


#' fetch_cfsv2_hindcast_nino34
#'
#' Fetch NCEP CFSv2 Niño 3.4 hindcast (+ real-time) member trajectories from the
#' NMME archive and return them as a long tibble -- the training sample for the
#' hindcast calibration of the real-time CFSv2 forecast
#' (`03_prepare_enso_cfsv2_calibration.R`). Mirrors `fetch_enso_index()`: cache
#' RDS + `_download_meta.yml`.
#'
#' Source: IRI Data Library `NMME/.NCEP-CFSv2/.HINDCAST/.PENTAD_SAMPLES_FULL/.sst`
#' (starts Jan 1982 onward, hindcast and real-time appended in one series),
#' averaged server-side over the Niño 3.4 box so the whole record is ~1.2 MB.
#' Values are ABSOLUTE SST (degC), not anomalies; dims S (start, "months since
#' 1960-01-01", 360-day calendar), L (lead 0.5..9.5 months), M (28 member slots,
#' 24 filled for most starts).
#'
#' The IRI Data Library is being retired (announced shutdown by end of Oct
#' 2026), so the raw NetCDF is kept permanently at `archive_file` and is the
#' source of truth once present: it is re-downloaded only when `refresh = TRUE`.
#' Fallback if the IRIDL endpoint is gone: the CCSR OPeNDAP server
#' `https://forecast.ccsr.columbia.edu/data/NMME/NOAA-NCEP/CFSv2/pentad_samples/{hindcast,forecast}/tos`
#' (no server-side averaging: subset the Niño 3.4 box per start, ~628 KB each).
#'
#' Access order (first hit wins):
#'   1. `<cache_dir>/enso_cfsv2_hindcast_nino34.rds` -- parsed cache
#'   2. `archive_file`                               -- archived raw NetCDF
#'   3. Live download of `url` into `archive_file`   (requires `ncdf4` to parse)
#'
#' @param url          IRIDL expression URL returning the box-averaged `data.nc`.
#' @param cache_dir    Directory for the parsed cache + meta files.
#' @param archive_file Path of the archived raw NetCDF.
#' @param refresh      If TRUE, re-download the NetCDF and rebuild the cache.
#' @return tibble(init_date, lead, target_date, member, sst): one row per
#'   non-missing (start, lead, member); `lead` = 0 for the start month itself.
fetch_cfsv2_hindcast_nino34 <- function(
  url = paste0(
    "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NCEP-CFSv2/",
    ".HINDCAST/.PENTAD_SAMPLES_FULL/.sst/X/190.5/239.5/RANGE/Y/-4.5/4.5/RANGE/",
    "%5BX/Y%5Daverage/S/264/801/RANGE/data.nc"
  ),
  cache_dir = "generated_data/",
  archive_file = file.path(
    cache_dir,
    "enso_cfsv2_hindcast",
    "nmme_cfsv2_pentad_full_nino34.nc"
  ),
  refresh = FALSE
) {
  cache_file <- file.path(cache_dir, "enso_cfsv2_hindcast_nino34.rds")
  meta_file <- file.path(
    cache_dir,
    "enso_cfsv2_hindcast_nino34_download_meta.yml"
  )

  if (!refresh && file.exists(cache_file)) {
    message(
      "fetch_cfsv2_hindcast_nino34: loading cached hindcast from ",
      cache_file
    )
    return(readRDS(cache_file))
  }
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop(
      "fetch_cfsv2_hindcast_nino34: `ncdf4` package required to read the NMME NetCDF."
    )
  }

  downloaded <- FALSE
  if (refresh || !file.exists(archive_file)) {
    message(
      "fetch_cfsv2_hindcast_nino34: downloading NMME CFSv2 Nino3.4 from\n  ",
      url,
      "\n  (IRIDL is slow: expect ~1.5-3 min)"
    )
    dir.create(dirname(archive_file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(archive_file, ".part")
    on.exit(unlink(tmp), add = TRUE)
    old_to <- getOption("timeout")
    options(timeout = max(old_to, 1800))
    on.exit(options(timeout = old_to), add = TRUE)
    tryCatch(
      utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
      error = function(e) {
        stop(
          "fetch_cfsv2_hindcast_nino34: download failed -- ",
          conditionMessage(e)
        )
      }
    )
    if (!file.exists(tmp) || file.info(tmp)$size == 0) {
      stop("fetch_cfsv2_hindcast_nino34: downloaded file is empty: ", url)
    }
    file.rename(tmp, archive_file)
    downloaded <- TRUE
  } else {
    message(
      "fetch_cfsv2_hindcast_nino34: parsing archived NetCDF ",
      archive_file
    )
  }

  nc <- ncdf4::nc_open(archive_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  s_units <- ncdf4::ncatt_get(nc, "S", "units")$value
  if (!identical(s_units, "months since 1960-01-01")) {
    stop("fetch_cfsv2_hindcast_nino34: unexpected S units '", s_units, "'")
  }
  S <- as.integer(round(nc$dim$S$vals))
  L <- nc$dim$L$vals
  M <- as.integer(nc$dim$M$vals)
  arr <- ncdf4::ncvar_get(nc, "sst", collapse_degen = FALSE)
  dnames <- vapply(nc$var$sst$dim, function(d) d$name, character(1))
  arr <- aperm(arr, match(c("S", "L", "M"), dnames)) # -> [S x L x M]

  # 360-day "months since 1960-01-01": S counts whole months, so the start
  # month is exact; L = 0.5 is the start month itself (lead 0).
  init_date <- as.Date(sprintf(
    "%04d-%02d-01",
    1960L + S %/% 12L,
    S %% 12L + 1L
  ))
  lead <- as.integer(round(L - 0.5))
  idx <- expand.grid(s = seq_along(S), l = seq_along(L), m = seq_along(M))
  out <- dplyr::tibble(
    init_date = init_date[idx$s],
    lead = lead[idx$l],
    member = M[idx$m],
    sst = as.vector(arr)
  ) %>%
    dplyr::filter(is.finite(sst)) %>%
    dplyr::mutate(target_date = lubridate::`%m+%`(init_date, months(lead))) %>%
    dplyr::select(init_date, lead, target_date, member, sst) %>%
    dplyr::arrange(init_date, lead, member)

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  saveRDS(out, cache_file)
  n_mem <- dplyr::count(dplyr::distinct(out, init_date, member), init_date)$n
  yaml::write_yaml(
    list(
      url = url,
      archive_file = archive_file,
      downloaded = if (downloaded) {
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
      } else {
        format(file.info(archive_file)$mtime, "%Y-%m-%dT%H:%M:%S")
      },
      history = ncdf4::ncatt_get(nc, "sst", "history")$value,
      units = ncdf4::ncatt_get(nc, "sst", "units")$value,
      init_range = as.character(range(out$init_date)),
      n_inits = length(unique(out$init_date)),
      leads = range(out$lead),
      members_per_init = as.list(table(n_mem)),
      nrow = nrow(out)
    ),
    meta_file
  )
  message(
    "fetch_cfsv2_hindcast_nino34: cached ",
    nrow(out),
    " records to ",
    cache_file
  )
  out
}


# CFSv2 hindcast calibration (CPC "PDF + spread correction" replica) ----------
#
# CPC's documented method (products/CFSv2/cfsv2fcst/CFSv2SST8210.html): "the
# forecast value is mapped to its percentile within the 1991-2020 hindcasts and
# forecasts, and the corresponding percentile value from the 1991-2020
# observations is used as the corrected forecast", then "rescales the forecast
# anomalies by a factor chosen such that the ensemble spread matches the RMSE of
# the ensemble-mean forecast". Replicated here per (start month, lead) on the
# NMME CFSv2 archive (`fetch_cfsv2_hindcast_nino34()`), each MEMBER percentile-
# mapped (capped at the reference range, as CPC's output is). One deliberate
# deviation: the spread factor is applied to the RAW member deviations around
# the mapped ensemble mean, because a record event pushes every member past the
# 24-member reference range, the cap maps them all to one value, and a
# multiplicative factor cannot restore a zero spread (CPC's 120-member pool has
# a longer tail). Validated against CPC's own PDFC/PDFSPRDC files and plots in
# `03_prepare_enso_cfsv2_calibration.R`.

#' Percentile-map values through a reference model sample onto an observed
#' sample (empirical CDFs, linear interpolation, capped at the range ends).
#'
#' @param x            Numeric values to map (model anomalies).
#' @param model_sorted Sorted reference model anomalies (all members, all years).
#' @param obs_sorted   Sorted observed anomalies for the same (start, lead) cells.
#' @return Numeric vector, same length as `x`.
.cfsv2_quantile_map <- function(x, model_sorted, obs_sorted) {
  nm <- length(model_sorted)
  pm <- (seq_len(nm) - 0.5) / nm
  p <- stats::approx(model_sorted, pm, xout = x, rule = 2, ties = "mean")$y
  as.numeric(stats::quantile(obs_sorted, probs = p, type = 7, names = FALSE))
}

#' fit_cfsv2_calibration
#'
#' Fit the CPC-style PDF + spread correction of CFSv2 Niño 3.4 forecasts from
#' the NMME hindcast archive, for one or more observational targets.
#'
#' @param hindcast Tibble from `fetch_cfsv2_hindcast_nino34()`.
#' @param obs      Named list of observed monthly anomaly tibbles
#'                 (`fetch_enso_index()` output: date, enso_index), one per
#'                 target (e.g. `list(ersst = ..., oisst = ...)`).
#' @param window   Integer c(first, last) start years of the training sample
#'                 (CPC: 1991-2020, hindcasts and real-time forecasts pooled).
#' @return List: `clim` (sm, lead, clim: model climatology, degC), `targets`
#'   (per target: `ref` = sm, lead, model, obs list-columns of sorted samples;
#'   `spread` = sm, lead, g, rmse, raw_spread), `meta`.
fit_cfsv2_calibration <- function(hindcast, obs, window = c(1991L, 2020L)) {
  stopifnot(is.list(obs), length(obs) > 0L, !is.null(names(obs)))
  h <- hindcast %>%
    dplyr::mutate(
      sm = as.integer(format(init_date, "%m")),
      yr = as.integer(format(init_date, "%Y"))
    ) %>%
    dplyr::filter(yr >= window[1], yr <= window[2])
  if (!nrow(h)) {
    stop(
      "fit_cfsv2_calibration: no hindcast starts in window ",
      paste(window, collapse = "-")
    )
  }
  clim <- h %>%
    dplyr::group_by(sm, lead) %>%
    dplyr::summarise(clim = mean(sst), .groups = "drop")
  h <- h %>%
    dplyr::left_join(clim, by = c("sm", "lead")) %>%
    dplyr::mutate(x = sst - clim)

  targets <- lapply(names(obs), function(tg) {
    o <- dplyr::transmute(
      obs[[tg]],
      target_date = as.Date(date),
      obs = enso_index
    )
    d <- dplyr::inner_join(h, o, by = "target_date") %>%
      dplyr::filter(is.finite(obs))
    ref <- d %>%
      dplyr::group_by(sm, lead) %>%
      dplyr::summarise(
        model = list(sort(x)),
        obs = list(sort(obs[!duplicated(init_date)])),
        n_years = dplyr::n_distinct(init_date),
        .groups = "drop"
      )
    # Spread factor: RMSE of the PDF-corrected ensemble mean / mean raw spread.
    spread <- d %>%
      dplyr::group_by(sm, lead, init_date) %>%
      dplyr::summarise(x = list(x), o = obs[1], .groups = "drop") %>%
      dplyr::left_join(ref, by = c("sm", "lead")) %>%
      dplyr::mutate(
        mq = mapply(
          function(xx, ms, os) mean(.cfsv2_quantile_map(xx, ms, os)),
          x,
          model,
          obs
        ),
        sr = vapply(x, stats::sd, numeric(1))
      ) %>%
      dplyr::group_by(sm, lead) %>%
      dplyr::summarise(
        rmse = sqrt(mean((mq - o)^2)),
        raw_spread = sqrt(mean(sr^2)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(g = rmse / raw_spread)
    list(ref = ref, spread = spread)
  })
  names(targets) <- names(obs)

  list(
    clim = clim,
    targets = targets,
    meta = list(
      method = "CPC PDF + spread correction replica (member percentile map, capped; RMSE spread on raw deviations)",
      window = as.integer(window),
      n_starts = dplyr::n_distinct(h$init_date),
      leads = range(h$lead),
      targets = names(obs)
    )
  )
}

#' apply_cfsv2_calibration
#'
#' Apply a `fit_cfsv2_calibration()` correction to a real-time CFSv2 member
#' matrix (rows = target months in date order, columns = members), as built in
#' `sample_enso_cfsv2_trajectories()`. The first row whose members differ is the
#' first forecast month (lead 0) and its calendar month is the start month.
#' Forecast rows: each member is percentile-mapped and the ensemble is rebuilt
#' as `mean(mapped) + g * (member - mean(member))`. Earlier (analysis) rows,
#' identical across members, are returned unchanged: they are CFSv2's own
#' analysis of observed months, not forecasts, and a lead-0 map compresses them
#' (Jun 2026: 1.59 -> 1.23 vs observed ERSST 1.50). The sampler anchors them to
#' the latest published ONI season instead (`.cfsv2_anchor_analysis_rows()`).
#'
#' @param cfs_wide    Numeric matrix [target month x member], rownames = dates.
#' @param calibration Output of `fit_cfsv2_calibration()`.
#' @param target      Name of the observational target (e.g. "ersst").
#' @return The corrected matrix (same shape and dimnames).
apply_cfsv2_calibration <- function(cfs_wide, calibration, target = "ersst") {
  tg <- calibration$targets[[target]]
  if (is.null(tg)) {
    stop(
      "apply_cfsv2_calibration: target '",
      target,
      "' not in calibration (have: ",
      paste(names(calibration$targets), collapse = ", "),
      ")"
    )
  }
  dates <- as.Date(rownames(cfs_wide))
  if (is.unsorted(dates)) {
    stop("apply_cfsv2_calibration: rows must be in date order")
  }
  row_sd <- apply(cfs_wide, 1, stats::sd, na.rm = TRUE)
  f0 <- unname(which(is.finite(row_sd) & row_sd > 1e-8)[1])
  if (is.na(f0)) {
    stop("apply_cfsv2_calibration: no forecast rows (members never differ)")
  }
  sm0 <- as.integer(format(dates[f0], "%m"))
  max_lead <- max(tg$ref$lead)
  cell <- function(s, l) {
    i <- which(tg$ref$sm == s & tg$ref$lead == l)
    j <- which(tg$spread$sm == s & tg$spread$lead == l)
    if (length(i) != 1L || length(j) != 1L) {
      stop(
        "apply_cfsv2_calibration: no calibration cell for start month ",
        s,
        ", lead ",
        l
      )
    }
    list(model = tg$ref$model[[i]], obs = tg$ref$obs[[i]], g = tg$spread$g[j])
  }
  out <- cfs_wide
  for (i in seq.int(f0, nrow(cfs_wide))) {
    x <- cfs_wide[i, ]
    ok <- is.finite(x)
    if (!any(ok)) {
      next
    }
    cc <- cell(sm0, min(i - f0, max_lead))
    mq <- mean(.cfsv2_quantile_map(x[ok], cc$model, cc$obs))
    out[i, ok] <- mq + cc$g * (x[ok] - mean(x[ok]))
  }
  attr(out, "first_forecast_row") <- f0
  out
}

#' Anchor CFSv2 analysis rows to the latest published ONI season.
#'
#' The CFSv2 file's analysis rows (before the first forecast month) fill the raw
#' Niño 3.4 series wherever the observed monthly series has no value yet (the
#' detrended ERSST file lags ONI by ~3 months). They are on CFSv2's analysis
#' scale, which in a strong event can sit well above ERSST (Jul-Sep 2026: mean
#' 2.54 vs ONI JAS 2.16). This shifts those unobserved analysis rows by one
#' constant so that the 3-month mean of the latest ONI season they fall in
#' equals published ONI; forecast rows are left alone.
#'
#' @param cfs_wide  Matrix [target month x member] (rownames = dates).
#' @param f0        Index of the first forecast row.
#' @param hist_lkp  Named vector of observed monthly raw Niño 3.4 (names = dates).
#' @param oni_lkp   Named vector of published ONI (names = centre-month dates).
#' @return List(cfs_wide, delta, season_centre).
.cfsv2_anchor_analysis_rows <- function(cfs_wide, f0, hist_lkp, oni_lkp) {
  dates <- as.Date(rownames(cfs_wide))
  ana <- seq_len(f0 - 1L)
  unobs <- ana[!is.finite(as.numeric(hist_lkp[as.character(dates[ana])]))]
  none <- list(cfs_wide = cfs_wide, delta = 0, season_centre = NA)
  if (!length(unobs)) {
    return(none)
  }
  centres <- as.Date(names(oni_lkp))[is.finite(oni_lkp)]
  centres <- sort(
    centres[
      centres %in%
        (dates[unobs] %m-% months(1L)) |
        centres %in% dates[unobs] |
        centres %in% (dates[unobs] %m+% months(1L))
    ],
    decreasing = TRUE
  )
  for (cc in as.list(centres)) {
    win <- seq(cc %m-% months(1L), cc %m+% months(1L), by = "1 month")
    vals <- vapply(
      win,
      function(d) {
        h <- as.numeric(hist_lkp[as.character(d)])
        if (is.finite(h)) {
          return(h)
        }
        i <- match(d, dates)
        if (!is.na(i) && i < f0) mean(cfs_wide[i, ], na.rm = TRUE) else NA_real_
      },
      numeric(1)
    )
    if (anyNA(vals)) {
      next
    }
    shift_rows <- match(win, dates)
    shift_rows <- shift_rows[!is.na(shift_rows) & shift_rows %in% unobs]
    if (!length(shift_rows)) {
      next
    }
    delta <- (3 * as.numeric(oni_lkp[as.character(cc)]) - sum(vals)) /
      length(shift_rows)
    cfs_wide[unobs, ] <- cfs_wide[unobs, ] + delta
    return(list(cfs_wide = cfs_wide, delta = delta, season_centre = cc))
  }
  none
}


#' join_enso_lags
#'
#' Attach ENSO lag columns (enso_lag0 ... enso_lag{max_lag}) to a survival or
#' outbreak-size dataset by left-joining the ONI value at the calendar month
#' obtained by subtracting each lag from the reference date column.
#'
#' @param df        Data frame containing a Date column given by `date_col`.
#' @param enso_dat  Tibble from fetch_enso_index(): columns year, month, enso_index.
#' @param date_col  Name of the Date column in df to lag from. Default "t_start".
#' @param max_lag   Maximum lag in months to generate (inclusive). Default 12L.
#' @param max_lead  Maximum lead in months (inclusive): also writes
#'   enso_lead1 ... enso_lead{max_lead}, the ONI k months AFTER `date_col`
#'   (base.R::enso_lag_col_name()). Default 0L = no lead columns (prior behaviour).
#'   A lead month past the end of the ENSO series is an error, not an NA: the
#'   single-lag fit would silently read it as neutral (dlnm.R NA -> 0).
#' @return df with enso_lag0 ... enso_lag{max_lag} (and any enso_lead<k>)
#'   columns appended (dbl).
join_enso_lags <- function(
  df,
  enso_dat,
  date_col = "t_start",
  max_lag = 12L,
  max_lead = 0L
) {
  enso_wide <- enso_dat %>%
    dplyr::select(year, month, enso_index)

  max_lead <- as.integer(max_lead %||% 0L)
  result <- df
  for (lag in seq(-max_lead, max_lag)) {
    col_name <- enso_lag_col_name(lag)
    lagged <- result %>%
      dplyr::mutate(
        .lag_date = .data[[date_col]] %m+% months(-as.integer(lag)),
        .lag_year = lubridate::year(.lag_date),
        .lag_month = lubridate::month(.lag_date)
      ) %>%
      dplyr::left_join(
        enso_wide %>% dplyr::rename(!!col_name := enso_index),
        by = c(".lag_year" = "year", ".lag_month" = "month")
      ) %>%
      dplyr::select(-.lag_date, -.lag_year, -.lag_month)
    result <- lagged
  }
  if (max_lead > 0L) {
    has0 <- !is.na(result$enso_lag0)
    for (k in seq_len(max_lead)) {
      col <- enso_lag_col_name(-k)
      n_na <- sum(has0 & is.na(result[[col]]))
      if (n_na > 0L) {
        stop(
          sprintf(
            paste0(
              "join_enso_lags: %d row(s) have %s but no %s -- the lead ",
              "month is past the end of the ENSO series (last %d-%02d). ",
              "Lower --enso_max_lead or refresh the ENSO index."
            ),
            n_na,
            "enso_lag0",
            col,
            max(enso_wide$year),
            max(enso_wide$month[enso_wide$year == max(enso_wide$year)])
          ),
          call. = FALSE
        )
      }
    }
  }
  result
}

#' panel_max_enso_lead -- largest k with an enso_lead<k> column in `panel` (0 when
#' none). The forecast ENSO samplers are asked for every lead the panel carries
#' (sample_enso_*(max_lead =)), so the per-replicate override covers each lead
#' column a deployed fit can read (simulator_core.R::simulate_one_tu()).
panel_max_enso_lead <- function(panel) {
  k <- suppressWarnings(as.integer(sub(
    "^enso_lead",
    "",
    grep("^enso_lead[0-9]+$", names(panel), value = TRUE)
  )))
  if (length(k)) max(k) else 0L
}

#' add_enso_leads_from_lag0
#'
#' Write enso_lead1 ... enso_lead{max_lead} (base.R::enso_lag_col_name()) from the
#' panel's OWN enso_lag0 series: enso_lead<k> at month t = enso_lag0 at month
#' t + k (the same value for every gid, as ENSO is basin-wide). Used by the
#' projection panel, whose enso_lag0 is observed ONI up to the last ONI month and
#' the CFSv2 point estimate (+ continuation) after it, so a lead reads the same
#' best estimate the panel already uses for lag 0 (03_prepare_sim_panel_projection.R).
#' Months whose t + k lies past the panel's last month get NA.
#'
#' @param panel    data frame with `date_col` and enso_lag0.
#' @param max_lead integer >= 1.
#' @param date_col Date column [default "t_start"].
#' @param extra    optional further panel(s) (list of data frames with
#'   `date_col` + enso_lag0) supplying the ENSO series for months AFTER the
#'   last month of `panel` (e.g. the projection panel behind the extended
#'   smoothing panel); months already in `panel` keep `panel`'s value.
#' @return panel with the lead columns added (existing lead columns replaced).
add_enso_leads_from_lag0 <- function(
  panel,
  max_lead,
  date_col = "t_start",
  extra = list()
) {
  max_lead <- as.integer(max_lead)
  if (is.na(max_lead) || max_lead < 1L) {
    return(panel)
  }
  month_series <- function(d) {
    if (!"enso_lag0" %in% names(d)) {
      stop("add_enso_leads_from_lag0: panel has no enso_lag0 column")
    }
    mon <- unique(data.frame(t = as.Date(d[[date_col]]), v = d$enso_lag0))
    mon <- mon[!is.na(mon$v), , drop = FALSE]
    dup <- mon$t[duplicated(mon$t)]
    if (length(dup)) {
      stop(
        "add_enso_leads_from_lag0: enso_lag0 differs between gids in month(s) ",
        paste(head(format(unique(dup)), 5), collapse = ", ")
      )
    }
    mon
  }
  mon <- month_series(panel)
  for (d in extra) {
    e <- month_series(d)
    mon <- rbind(mon, e[e$t > max(mon$t), , drop = FALSE])
  }
  lkp <- stats::setNames(mon$v, format(mon$t))
  t0 <- as.Date(panel[[date_col]])
  for (k in seq_len(max_lead)) {
    panel[[enso_lag_col_name(
      -k
    )]] <- unname(lkp[format(lubridate::add_with_rollback(t0, months(k)))])
  }
  panel
}

# ---------------------------------------------------------------------------
# Forecast ENSO utilities — C3S plume ingestion and trajectory sampling
# ---------------------------------------------------------------------------

#' Read a C3S ENSO Niño 3.4 forecast plume from a companion data file or,
#' as a last resort, from a PNG via pixel digitization.
#'
#' Priority order:
#'   1. `<png_path_noext>.rds`  — companion RDS (tibble already in final shape)
#'   2. `<png_path_noext>.csv`  — companion CSV (read_csv with same schema)
#'   3. PNG pixel digitization  — crude fallback (requires `png` package)
#'
#' Returned tibble columns:
#'   date            <Date>    First day of the forecast month (YYYY-MM-01)
#'   q_level         <dbl>     Quantile level ∈ {0.10, 0.25, 0.50, 0.75, 0.90}
#'   nino34_anomaly  <dbl>     Niño 3.4 sea-surface temperature anomaly (°C)
#'
#' The PNG digitization branch assumes the C3S plume image uses the standard
#' ECMWF colour scheme: five quantile bands with distinct fills in the order
#' [10th, 25th, 50th, 75th, 90th].  It requires the `png` package and hard-
#' codes axis ranges; provide a companion .rds/.csv whenever possible.
#'
#' @param png_path  Path to the `.png` file (or its no-extension stem).
#' @return tibble with columns (date, q_level, nino34_anomaly).
read_enso_c3s_plume <- function(png_path) {
  stem <- tools::file_path_sans_ext(png_path)

  # 1. Companion RDS
  rds_path <- paste0(stem, ".rds")
  if (file.exists(rds_path)) {
    message("[read_enso_c3s_plume] reading companion RDS: ", rds_path)
    dat <- readRDS(rds_path)
    required <- c("date", "q_level", "nino34_anomaly")
    missing <- setdiff(required, names(dat))
    if (length(missing) == 0) {
      return(tibble::as_tibble(dat))
    }
    stop(
      "[read_enso_c3s_plume] companion RDS missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  # 2. Companion CSV
  csv_path <- paste0(stem, ".csv")
  if (file.exists(csv_path)) {
    message("[read_enso_c3s_plume] reading companion CSV: ", csv_path)
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("[read_enso_c3s_plume] readr required to read companion CSV.")
    }
    dat <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
      dplyr::mutate(date = as.Date(date))
    required <- c("date", "q_level", "nino34_anomaly")
    missing <- setdiff(required, names(dat))
    if (length(missing) == 0) {
      return(tibble::as_tibble(dat))
    }
    stop(
      "[read_enso_c3s_plume] companion CSV missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  # 3. PNG pixel digitization fallback
  if (!file.exists(png_path)) {
    stop(
      "[read_enso_c3s_plume] file not found: ",
      png_path,
      ". Provide a companion .rds or .csv alongside the PNG."
    )
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop(
      "[read_enso_c3s_plume] `png` package required for PNG digitization; ",
      "install.packages('png') or provide a companion .rds/.csv."
    )
  }

  warning(
    "[read_enso_c3s_plume] falling back to PNG pixel digitization. ",
    "Results are approximate — provide a companion .rds/.csv for ",
    "production use."
  )

  img <- png::readPNG(png_path) # [rows × cols × channels]; origin top-left
  nr <- nrow(img)
  nc <- ncol(img)

  # Hard-coded axis ranges for the ECMWF C3S Niño 3.4 plume (May 2026 image).
  # Forecast period: 2026-06 → 2027-05 (12 months).
  # Y axis: -3 to +3 °C anomaly.
  # Adjust these if the image axis ranges differ.
  date_start <- as.Date("2026-06-01")
  n_months <- 12L
  y_min <- -3.0
  y_max <- 3.0
  # Estimate plot bounding box by trimming near-white margins.
  # "Near white": all three channels > 0.90.
  is_bg <- apply(img[,, 1:3], c(1, 2), function(rgb) all(rgb > 0.90))
  row_has_content <- !apply(is_bg, 1, all)
  col_has_content <- !apply(is_bg, 2, all)
  r_lo <- which(row_has_content)[1]
  r_hi <- tail(which(row_has_content), 1)
  c_lo <- which(col_has_content)[1]
  c_hi <- tail(which(col_has_content), 1)

  plot_rows <- seq(r_lo, r_hi)
  plot_cols <- seq(c_lo, c_hi)
  H <- length(plot_rows)
  W <- length(plot_cols)

  # Map pixel → (date, Niño 3.4): linear scale.
  col_to_date <- function(col_px) {
    frac <- (col_px - c_lo) / (c_hi - c_lo)
    date_start + floor(frac * n_months) * 30L # approx 1-month grid
  }
  row_to_val <- function(row_px) {
    # Row 0 = top of image = y_max; row H-1 = bottom = y_min.
    frac <- (row_px - r_lo) / (r_hi - r_lo)
    y_max - frac * (y_max - y_min)
  }

  # C3S quantile band colours (approximate RGB centroids, tolerance 0.12).
  # Order: Q10, Q25, Q50, Q75, Q90.
  q_levels <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  q_colours <- list(
    c(0.40, 0.62, 0.80), # Q10: muted blue
    c(0.65, 0.80, 0.90), # Q25: light blue
    c(0.20, 0.40, 0.65), # Q50: medium blue (median line)
    c(0.85, 0.65, 0.40), # Q75: light orange
    c(0.80, 0.35, 0.20) # Q90: dark orange/red
  )
  tol <- 0.15

  rows_out <- list()
  for (col_px in plot_cols) {
    d <- col_to_date(col_px)
    r_col <- img[plot_rows, col_px, 1]
    g_col <- img[plot_rows, col_px, 2]
    b_col <- img[plot_rows, col_px, 3]
    for (qi in seq_along(q_levels)) {
      ref <- q_colours[[qi]]
      matches <- abs(r_col - ref[1]) < tol &
        abs(g_col - ref[2]) < tol &
        abs(b_col - ref[3]) < tol
      if (!any(matches)) {
        next
      }
      val <- row_to_val(plot_rows[which(matches)])
      rows_out[[length(rows_out) + 1L]] <- tibble::tibble(
        date = as.Date(format(d, "%Y-%m-01")),
        q_level = q_levels[qi],
        nino34_anomaly = mean(val, na.rm = TRUE)
      )
    }
  }

  if (length(rows_out) == 0) {
    stop(
      "[read_enso_c3s_plume] PNG digitization found no quantile-band pixels. ",
      "Provide a companion .rds/.csv."
    )
  }

  dplyr::bind_rows(rows_out) %>%
    dplyr::group_by(date, q_level) %>%
    dplyr::summarise(
      nino34_anomaly = mean(nino34_anomaly, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(date, q_level)
}


#' Sample per-replicate ENSO (ONI-equivalent) forecast trajectories.
#'
#' Builds, for each simulation month m and each ENSO lag l, an [M × n_rep]
#' matrix where M = number of simulation months.  Historical ONI months
#' (index date ≤ max(historical_oni$date)) are identical across all replicates;
#' forecast months beyond the ONI cache are sampled from the C3S plume.
#'
#' Niño 3.4 → ONI conversion (NOAA definition):
#'   ONI(M) = [Niño3.4(M-1) + Niño3.4(M) + Niño3.4(M+1)] / 3
#' The sampling therefore works on raw Niño 3.4, stitches historical (from
#' `fetch_enso_index(source="nino34")`) with the C3S plume, applies the 3-month
#' centered rolling mean, and then pre-lags to produce the final matrices.
#'
#' @param plume_data      Tibble from `read_enso_c3s_plume()`:
#'                        (date, q_level, nino34_anomaly).
#' @param n_rep           Number of forecast replicates to generate.
#' @param historical_oni  Tibble from `fetch_enso_index(source="oni")`:
#'                        must contain (year, month, enso_index).
#' @param historical_nino34 Tibble from `fetch_enso_index(source="nino34")`:
#'                        raw monthly Niño 3.4.  If NULL, falls back to
#'                        using historical_oni directly (less accurate).
#' @param sim_months      Date vector of simulation months (first day of each).
#'                        If NULL, derived from plume dates + lag range.
#' @param max_lag         Maximum ENSO lag to pre-compute (default 12L).
#' @param seed            Optional integer seed for reproducibility.
#' @return Named list of [M × n_rep] numeric matrices, one entry per lag name
#'         ("enso_lag0", "enso_lag1", …, "enso_lag{max_lag}").
sample_enso_forecast_trajectories <- function(
  plume_data,
  n_rep,
  historical_oni,
  historical_nino34 = NULL,
  sim_months = NULL,
  max_lag = 12L,
  seed = NULL,
  max_lead = 0L
) {
  max_lead <- as.integer(max_lead %||% 0L) # > 0: also enso_lead1..max_lead
  if (!is.null(seed)) {
    set.seed(seed)
  }
  stopifnot(
    all(c("date", "q_level", "nino34_anomaly") %in% names(plume_data)),
    all(c("year", "month", "enso_index") %in% names(historical_oni))
  )

  oni_max_date <- with(
    historical_oni,
    as.Date(sprintf(
      "%d-%02d-01",
      max(year),
      month[which.max(year * 12 + month)]
    ))
  )
  plume_dates <- sort(unique(plume_data$date))

  # If sim_months not given, infer from plume extent + enough room for all lags.
  if (is.null(sim_months)) {
    all_dates <- sort(unique(c(oni_max_date, plume_dates)))
    sim_months <- seq.Date(all_dates[1], tail(all_dates, 1), by = "1 month")
  }
  M <- length(sim_months)

  # Build a look-up for historical ONI: Date → scalar.
  oni_lkp <- stats::setNames(
    historical_oni$enso_index,
    as.Date(sprintf("%d-%02d-01", historical_oni$year, historical_oni$month))
  )

  # Build a look-up for historical raw Niño 3.4 (if available).
  nino34_lkp <- if (
    !is.null(historical_nino34) &&
      all(c("year", "month", "enso_index") %in% names(historical_nino34))
  ) {
    stats::setNames(
      historical_nino34$enso_index,
      as.Date(sprintf(
        "%d-%02d-01",
        historical_nino34$year,
        historical_nino34$month
      ))
    )
  } else {
    NULL
  }

  # Pivot plume to wide: row = date, column = q_level.
  plume_wide <- plume_data %>%
    dplyr::mutate(q_col = paste0("q", sprintf("%.0f", q_level * 100))) %>%
    tidyr::pivot_wider(
      id_cols = date,
      names_from = q_col,
      values_from = nino34_anomaly
    ) %>%
    dplyr::arrange(date)
  plume_q_cols <- setdiff(names(plume_wide), "date")
  plume_qs <- as.numeric(sub("^q", "", plume_q_cols)) / 100

  # For each replicate, build a complete monthly raw Niño 3.4 series from
  # historical start → forecast end, then apply the 3-month rolling mean
  # to get the ONI-equivalent series.  The series must cover enough months
  # on each side of sim_months to handle the max_lag requirement AND the
  # rolling mean edge trim (±1 month).
  series_start <- min(sim_months) %m-% months(max_lag + 1L)
  series_end <- max(sim_months) %m+% months(1L + max_lead)
  series_dates <- seq.Date(series_start, series_end, by = "1 month")
  L_series <- length(series_dates)

  # Pre-sample Niño 3.4 for forecast months per replicate.
  # Historical months: use nino34_lkp if available, else oni_lkp as proxy.
  hist_lkp_for_nino34 <- if (!is.null(nino34_lkp)) nino34_lkp else oni_lkp
  base_series <- as.numeric(hist_lkp_for_nino34[as.character(series_dates)])
  # Months beyond ONI cache: indices needing sampling.
  forecast_idx <- which(series_dates > oni_max_date)

  # Quantile-based inverse-CDF sampling for each forecast month.
  # We use the closest plume date for each series date.
  sample_from_plume <- function(target_date, u) {
    if (nrow(plume_wide) == 0) {
      return(rep(0, n_rep))
    }
    # Find closest plume date.
    diffs <- abs(as.numeric(plume_wide$date - target_date))
    row_i <- which.min(diffs)
    q_vals <- as.numeric(plume_wide[row_i, plume_q_cols])
    finite_mask <- is.finite(q_vals) & is.finite(plume_qs)
    if (sum(finite_mask) < 2L) {
      return(rep(0, n_rep))
    } # fallback: 0 anomaly
    # Interpolate quantile function at the replicate-specific quantile levels.
    stats::approx(
      plume_qs[finite_mask],
      q_vals[finite_mask],
      xout = u,
      rule = 2
    )$y
  }

  # Build [L_series × n_rep] Niño 3.4 matrix; historical rows are identical
  # across columns; forecast rows vary.
  # One quantile draw per replicate, shared across all forecast months, so that
  # individual trajectories are temporally coherent (same plume quantile
  # throughout the forecast horizon).
  u_rep <- stats::runif(n_rep)
  nino34_mat <- matrix(base_series, nrow = L_series, ncol = n_rep)
  for (fi in forecast_idx) {
    nino34_mat[fi, ] <- sample_from_plume(series_dates[fi], u_rep)
  }

  # Apply 3-month centered rolling mean → ONI-equivalent per replicate.
  # zoo::rollmean(k=3, align="center") drops first and last row.
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop(
      "sample_enso_forecast_trajectories: `zoo` package required for rollmean."
    )
  }
  oni_mat <- apply(nino34_mat, 2, function(col_vec) {
    rm <- as.numeric(zoo::rollmean(
      col_vec,
      k = 3L,
      align = "center",
      fill = NA
    ))
    rm
  }) # [L_series × n_rep]; first and last row are NA
  # Trim to the valid interior (drop first and last month from rolling mean).
  valid_idx <- which(!is.na(oni_mat[, 1]))
  oni_dates <- series_dates[valid_idx]
  oni_mat_trim <- oni_mat[valid_idx, , drop = FALSE]

  # Look up ONI per (sim_month - lag) for each lag; store in named list.
  result <- vector("list", max_lag + 1L)
  names(result) <- paste0("enso_lag", seq(0L, max_lag))

  for (lag in c(seq(0L, max_lag), -seq_len(max_lead))) {
    # leads (negative) after the lags
    lag_dates <- sim_months %m-% months(lag)
    lag_mat <- matrix(NA_real_, nrow = M, ncol = n_rep)
    for (mi in seq_len(M)) {
      d_key <- as.character(lag_dates[mi])
      # Check if this date is in the trimmed ONI matrix.
      idx <- match(d_key, as.character(oni_dates))
      if (!is.na(idx)) {
        lag_mat[mi, ] <- oni_mat_trim[idx, ]
      } else {
        # Fallback to scalar historical ONI (pre-rolling-mean available dates).
        scalar <- oni_lkp[d_key]
        if (is.na(scalar) && lag < 0L) {
          stop(
            sprintf(
              "sample_enso_forecast_trajectories: no ENSO for lead %d at %s.",
              -lag,
              d_key
            ),
            call. = FALSE
          )
        }
        lag_mat[mi, ] <- if (!is.na(scalar)) scalar else 0.0
      }
    }
    result[[enso_lag_col_name(lag)]] <- lag_mat
  }

  result
}


#' Sample per-replicate ENSO (ONI-equivalent) forecast trajectories from CFSv2.
#'
#' CFSv2 analogue of `sample_enso_forecast_trajectories()`.  Instead of inverse-CDF
#' sampling a C3S quantile plume, each replicate is assigned one CFSv2 ensemble
#' member (sampled with replacement) whose raw Niño 3.4 anomaly trajectory is
#' rebased to the historical seam, stitched onto the historical record, and
#' passed through the 3-month centered rolling mean → ONI-equivalent.  Returns
#' the identical named-list `[M × n_rep]` per-lag structure so downstream
#' consumers (utils_simulations loaders, simulate_one_tu) are unchanged.
#'
#' Member → replicate: `sample.int(n_members, n_rep, replace=TRUE)` under `seed`;
#' one member per replicate across all leads (temporally coherent, matching the
#' one-draw-per-replicate C3S behaviour).  Scales the ~40 members to any n_rep.
#'
#' Rebase (`rebase` arg): CFSv2 anomalies use the 1991–2020 base; the historical
#' raw series (`fetch_enso_index("nino34")`) uses a different base.  "seam" applies
#' an additive shift so the CFSv2 ensemble mean matches the historical raw value
#' at the latest overlapping month; "none" leaves values unshifted; a numeric
#' value (or a numeric string, e.g. the CLI's "0.3") is used as an explicit
#' offset. `rebase` is IGNORED when `calibration` is set: the calibrated path
#' anchors the analysis rows to observed ONI itself and forces `rebase` to
#' "none" (cfsv2_correction_provenance() records the effective value).
#'
#' @param cfsv2_data       Tibble from `fetch_enso_cfsv2_ensemble()`:
#'                         (init_date, target_date, lead, member, nino34_anomaly).
#' @param n_rep            Number of forecast replicates to generate.
#' @param historical_oni   Tibble from `fetch_enso_index(source="oni")`.
#' @param historical_nino34 Tibble from `fetch_enso_index(source="nino34")`; if
#'                         NULL, historical_oni is used as a Niño 3.4 proxy.
#' @param sim_months       Date vector of simulation months (first day of each).
#' @param max_lag          Maximum ENSO lag to pre-compute (default 12L).
#' @param rebase           "seam" (default) | "none" | numeric offset (a
#'                         numeric string is parsed). Ignored (forced to
#'                         "none") when `calibration` is set.
#' @param sd_correction    If TRUE, crudely rescale the raw CFSv2 anomaly by
#'                         k = sd(observed ONI) / sd(raw CFSv2 members) before
#'                         the seam rebase. Mutually
#'                         exclusive with `target_peak`.
#' @param target_peak      If set (numeric ONI-equivalent, e.g. 2.5), rescale
#'                         the raw CFSv2 anomaly by k = target_peak / raw_peak,
#'                         where raw_peak is the max of the 3-month-rolling-mean
#'                         ONI-equivalent of the ENSEMBLE-MEAN raw forecast over
#'                         the forecast horizon. Applied at
#'                         the same point as `sd_correction` (before the seam
#'                         rebase), and mutually exclusive with it — this is the
#'                         preferred correction, calibrated to an external
#'                         outlook rather than an incidental sd-ratio artefact.
#' @param calibration      Optional `fit_cfsv2_calibration()` object (or a path to
#'                         its RDS, `03_prepare_enso_cfsv2_calibration.R`):
#'                         applies the CPC-style hindcast PDF + spread correction
#'                         (`apply_cfsv2_calibration()`) before the seam rebase.
#'                         Mutually exclusive with `sd_correction` / `target_peak`.
#' @param calibration_target Observational target of the calibration: "ersst"
#'                         (ONI-consistent, default) or "oisst" (CPC's own).
#' @param seed             Optional integer seed for reproducibility.
#' @return Named list of [M × n_rep] numeric matrices, one per lag
#'         ("enso_lag0", …, "enso_lag{max_lag}").
# ---------------------------------------------------------------------------
# ENSO continuation past the CFSv2 span
# ---------------------------------------------------------------------------

#' fit_oni_ar1
#'
#' Lag-1 autocorrelation of a monthly Nino3.4 anomaly series (named numeric,
#' names "YYYY-MM-DD", e.g. the sampler's raw Nino3.4 lookup) over
#' `fit_window` ("YYYY_YYYY", inclusive years). Returns list(phi, n, window).
fit_oni_ar1 <- function(series, fit_window = "1991_2020") {
  d <- as.Date(names(series))
  yr <- as.integer(strsplit(fit_window, "_")[[1L]])
  keep <- is.finite(series) &
    !is.na(d) &
    as.integer(format(d, "%Y")) >= yr[1L] &
    as.integer(format(d, "%Y")) <= yr[2L]
  x <- as.numeric(series[keep])[order(d[keep])]
  if (length(x) < 24L) {
    stop(
      "fit_oni_ar1: fewer than 24 monthly values in ",
      fit_window,
      call. = FALSE
    )
  }
  x <- x - mean(x)
  phi <- sum(x[-1L] * x[-length(x)]) / sum(x^2)
  list(phi = as.numeric(phi), n = length(x), window = fit_window)
}

#' continue_enso_ar1
#'
#' Deterministic AR(1) decay of each replicate's seam value towards 0:
#' row k = phi^k * seam. Returns an [n_months x n_rep] matrix.
continue_enso_ar1 <- function(seam_values_by_rep, n_months, phi) {
  stopifnot(is.numeric(phi), length(phi) == 1L, abs(phi) < 1)
  v <- as.numeric(seam_values_by_rep)
  v[!is.finite(v)] <- 0
  outer(phi^seq_len(n_months), v)
}

# ---------------------------------------------------------------------------
# Analog continuation (skill scoring: utils/enso_continuation_skill.R). Not to be
# confused with the counterfactual `enso.counterfactual.scenario: analog`
# (sample_enso_analog_scenarios()), which replays an observed period.
# ---------------------------------------------------------------------------

#' enso_continuation_features
#'
#' Per-month analog features of one ONI-scale series: v = ONI(i),
#' tend = ONI(i) - ONI(i-3), phase = ONI(i-6), and `ok` = a usable analog start
#' (all features defined and at least `horizon` months of future).
enso_continuation_features <- function(x, mon, model, horizon) {
  n <- length(x)
  i <- seq_len(n)
  lagk <- function(k) c(rep(NA_real_, k), x[seq_len(n - k)])[seq_len(n)]
  data.frame(
    model = model,
    i = i,
    mon = as.integer(mon),
    v = x,
    tend = x - lagk(3L),
    phase = lagk(6L),
    ok = i > 6L & i + horizon <= n,
    stringsAsFactors = FALSE
  )
}

#' enso_continuation_library_digest
#'
#' xxhash64 of the library content (series + features + feature sds), recorded
#' in its meta and copied into run provenance (partition_meta.yml, manifests).
enso_continuation_library_digest <- function(lib) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    return(NA_character_)
  }
  digest::digest(
    list(lib$series, lib$features, unname(lib$feature_sd)),
    algo = "xxhash64"
  )
}

#' load_enso_continuation_library
#'
#' Read the continuation library (03_prepare_enso_continuation_library.R) and
#' index its candidate starts by calendar month. Stops with the build command
#' when the file is missing or malformed.
load_enso_continuation_library <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(
      sprintf(
        paste0(
          "ENSO analog continuation: library not found: '%s'. Build it with ",
          "`Rscript analysis/03_prepare_enso_continuation_library.R`."
        ),
        path %||% ""
      ),
      call. = FALSE
    )
  }
  lib <- qs2::qs_read(path)
  need <- c("series", "features", "feature_sd", "meta")
  if (!all(need %in% names(lib)) || is.null(lib$meta$horizon_months)) {
    stop(
      "ENSO analog continuation: '",
      path,
      "' is not a continuation library (missing ",
      paste(setdiff(need, names(lib)), collapse = ", "),
      "); rebuild it with ",
      "`Rscript analysis/03_prepare_enso_continuation_library.R`.",
      call. = FALSE
    )
  }
  cand <- lib$features[
    lib$features$ok,
    c("model", "i", "mon", "v", "tend", "phase")
  ]
  lib$by_month <- split(cand, factor(cand$mon, levels = 1:12))
  lib$path <- path
  lib
}

#' enso_continuation_provenance
#'
#' Run-identity record of a simulator.enso.continuation block: list(method) for
#' neutral / ar1 (plus phi when set), and for analog also the library path, its
#' content digest and k / blend_months / min_sep. Reads the library (stops with
#' the build command when it is missing), so a run or a manifest naming an
#' analog continuation cannot point at a library that does not exist.
#' The numeric CFSv2 rebase offset of `rebase` (a number, or a string that
#' parses as one, e.g. "0.3"); NULL for "seam" / "none" / anything else.
.rebase_offset <- function(rebase) {
  if (length(rebase) != 1L) {
    return(NULL)
  }
  v <- if (is.numeric(rebase)) {
    as.numeric(rebase)
  } else {
    suppressWarnings(as.numeric(as.character(rebase)))
  }
  if (is.finite(v)) v else NULL
}

#' cfsv2_correction_provenance
#'
#' The CFSv2 correction of a projection panel / forecast: the resolved
#' correction mode and the EFFECTIVE rebase ("none" under calibration, which
#' ignores `rebase`; a parsed numeric offset as a number). Returns an empty
#' list at the default (correction "sd", rebase "seam"), so a default stamp
#' gains no field (same pattern as enso_continuation_provenance()).
#' @param correction resolved mode: "none" | "sd" | "target_peak" | "calibration"
#'   (NA when no CFSv2 was used)
#' @param rebase the requested rebase ("seam" | "none" | numeric / numeric string)
cfsv2_correction_provenance <- function(correction, rebase) {
  if (is.null(correction) || is.na(correction)) {
    return(list())
  }
  eff <- if (identical(correction, "calibration")) {
    "none"
  } else {
    (.rebase_offset(rebase) %||% as.character(rebase))
  }
  if (identical(correction, "sd") && identical(eff, "seam")) {
    return(list())
  }
  list(cfsv2_correction = correction, cfsv2_rebase_effective = eff)
}

enso_continuation_provenance <- function(cont) {
  cont <- cont %||% list()
  method <- cont[["method"]] %||% "neutral"
  # The analog seam: recorded only when not the default last_cfsv2,
  # so a default run's provenance is unchanged.
  seam <- as.character(cont[["seam"]] %||% "last_cfsv2")
  add_seam <- function(out) {
    if (!identical(seam, "last_cfsv2")) {
      out$seam <- seam
    }
    out
  }
  if (!identical(method, "analog")) {
    out <- list(method = method)
    if (!is.null(cont[["phi"]])) {
      out$phi <- cont[["phi"]]
    }
    return(add_seam(out))
  }
  path <- cont[["library"]] %||% "generated_data/enso_continuation_library.qs"
  lib <- load_enso_continuation_library(path)
  add_seam(list(
    method = method,
    library = path,
    library_digest = lib$meta$digest %||% enso_continuation_library_digest(lib),
    k = as.integer(cont[["k"]] %||% 30L),
    blend_months = as.numeric(cont[["blend_months"]] %||% 3),
    min_sep = as.integer(cont[["min_sep"]] %||% 24L),
    max_per_model = {
      mpm <- as.numeric(cont[["max_per_model"]] %||% Inf)
      if (is.finite(mpm)) as.integer(mpm) else "none"
    },
    exclude_models = as.list(as.character(unlist(cont[["exclude_models"]])))
  ))
}

#' enso_analog_neighbours
#'
#' The k nearest library starts to one seam state (list(mon, v, tend, phase)):
#' same calendar month, squared distance over the features standardised by the
#' library's feature sds, at most one start per `min_sep` months within a model
#' and at most `max_per_model` starts per model (greedy in distance order);
#' models in `exclude_models` are skipped. Returns a data.frame(model, i).
#'
#' `max_per_model` weights the models more equally: without it a model counts
#' in proportion to its piControl length (CM2.1-NMME is 30 % of the library and
#' supplied 12 of 30 analogs after the 2026/27 El Nino, with the strongest
#' post-El-Nino cold bias of the 11 models).
#' `horizon` (NULL = the unrestricted candidates): keep only starts with at least
#' `horizon` months of future in their model's series (i + horizon <= length),
#' for a counterfactual pool whose trajectories outrun the library's own
#' --horizon_months. The deployed metric (`feature_sd`) is kept; contrast
#' enso_continuation_skill.R::enso_continuation_library_rehorizon(), which
#' recomputes it for the skill scoring.
enso_analog_neighbours <- function(
  lib,
  state,
  k = 30L,
  min_sep = 24L,
  feats = c("v", "tend", "phase"),
  max_per_model = Inf,
  exclude_models = character(0),
  horizon = NULL
) {
  cand <- lib$by_month[[as.character(state$mon)]]
  if (length(exclude_models)) {
    unknown <- setdiff(exclude_models, names(lib$series))
    if (length(unknown)) {
      stop(
        "ENSO analog continuation: exclude_models not in the library: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    cand <- cand[!cand$model %in% exclude_models, , drop = FALSE]
  }
  if (!is.null(horizon) && !is.null(cand)) {
    n_of <- vapply(lib$series, length, integer(1))
    cand <- cand[
      cand$i + as.integer(horizon) <= n_of[cand$model],
      ,
      drop = FALSE
    ]
  }
  if (is.null(cand) || !nrow(cand)) {
    stop(
      "ENSO analog continuation: no library start in calendar month ",
      state$mon,
      call. = FALSE
    )
  }
  max_per_model <- as.numeric(max_per_model %||% Inf)
  if (is.na(max_per_model) || max_per_model < 1) {
    stop(
      "ENSO analog continuation: max_per_model must be >= 1 (or Inf)",
      call. = FALSE
    )
  }
  d <- .enso_analog_distance(lib, cand, state, feats)
  cand <- cand[order(d), c("model", "i")]
  sel <- integer(0)
  for (r in seq_len(nrow(cand))) {
    if (length(sel) >= k) {
      break
    }
    same <- cand$model[sel] == cand$model[r]
    if (
      sum(same) < max_per_model &&
        (!any(same) || all(abs(cand$i[sel][same] - cand$i[r]) >= min_sep))
    ) {
      sel <- c(sel, r)
    }
  }
  cand[sel, , drop = FALSE]
}

#' Squared standardised distance of library candidates `cand` to a seam state
#' over `feats` (the deployed metric: lib$feature_sd).
.enso_analog_distance <- function(
  lib,
  cand,
  state,
  feats = c("v", "tend", "phase")
) {
  d <- 0
  for (f in feats) {
    d <- d + ((cand[[f]] - state[[f]]) / lib$feature_sd[[f]])^2
  }
  d
}

#' enso_analog_paths
#'
#' The faded-join paths of analogs `nn` (data.frame(model, i)) from one seam
#' value: row r, column j = analog(i_r + j) + (seam_v - analog(i_r)) * rho^j,
#' j = 1..n_months, rho = exp(-1 / blend_months) (Inf: constant offset).
#' path(0) would be seam_v itself, so column 1 is the month after the seam.
#' The one implementation of the join, shared by continue_enso_analog() and the
#' counterfactual pool (build_lenssen_counterfactual_pool()).
#' @return [nrow(nn) x n_months] numeric matrix
enso_analog_paths <- function(lib, nn, seam_v, n_months, blend_months = 3) {
  fade <- if (is.infinite(blend_months)) {
    rep(1, n_months)
  } else {
    exp(-seq_len(n_months) / blend_months)
  }
  out <- matrix(NA_real_, nrow = nrow(nn), ncol = n_months)
  for (r in seq_len(nrow(nn))) {
    s <- lib$series[[nn$model[r]]]
    i <- nn$i[r]
    out[r, ] <- s[i + seq_len(n_months)] + (seam_v - s[i]) * fade
  }
  out
}

#' enso_seam_state
#'
#' The analog seam state of month `seam_date` read from PUBLISHED ONI
#' (`oni_lkp`: named numeric keyed "YYYY-MM-DD"): list(mon, v = ONI(m),
#' tend = ONI(m) - ONI(m-3), phase = ONI(m-6)), the features of
#' enso_continuation_features(). Stops when a month is missing.
enso_seam_state <- function(oni_lkp, seam_date) {
  m <- as.Date(format(as.Date(seam_date), "%Y-%m-01"))
  at <- function(k) {
    d <- format(seq(m, by = "-1 month", length.out = k + 1L)[k + 1L])
    v <- as.numeric(oni_lkp[d])
    if (!length(v) || !is.finite(v)) {
      stop(
        "enso_seam_state: no published ONI for ",
        d,
        " (seam ",
        format(m),
        ").",
        call. = FALSE
      )
    }
    v
  }
  v <- at(0L)
  list(
    mon = as.integer(format(m, "%m")),
    v = v,
    tend = v - at(3L),
    phase = at(6L)
  )
}

#' continue_enso_analog
#'
#' Continue each replicate's ONI past its seam month: replicate r draws one of
#' the k nearest library analogs of its seam state (v, tend, phase in
#' `seam_*[r]`, calendar month `seam_month`) and follows it, joined with a
#' fading offset: path(j) = analog(i + j) + (seam_v - analog(i)) * rho^j,
#' rho = exp(-1 / blend_months) (blend_months = Inf: constant offset).
#' Neighbour sets are computed once per distinct seam state (replicates sharing
#' a CFSv2 member share it). Returns an [n_months x n_rep] matrix.
continue_enso_analog <- function(
  lib,
  seam_v,
  seam_tend,
  seam_phase,
  seam_month,
  n_months,
  k = 30L,
  blend_months = 3,
  min_sep = 24L,
  max_per_model = Inf,
  exclude_models = character(0)
) {
  n_rep <- length(seam_v)
  stopifnot(
    length(seam_tend) == n_rep,
    length(seam_phase) == n_rep,
    is.numeric(blend_months),
    blend_months > 0
  )
  if (n_months > lib$meta$horizon_months) {
    stop(
      sprintf(
        paste0(
          "ENSO analog continuation: %d month(s) needed past the CFSv2 span but the ",
          "library holds %d; rebuild it with a longer --horizon_months."
        ),
        n_months,
        lib$meta$horizon_months
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(c(seam_v, seam_tend, seam_phase)))) {
    stop("ENSO analog continuation: non-finite seam state", call. = FALSE)
  }
  key <- sprintf("%.10g|%.10g|%.10g", seam_v, seam_tend, seam_phase)
  nn_cache <- list()
  out <- matrix(NA_real_, nrow = n_months, ncol = n_rep)
  for (r in seq_len(n_rep)) {
    nn <- nn_cache[[key[r]]]
    if (is.null(nn)) {
      nn <- enso_analog_neighbours(
        lib,
        list(
          mon = seam_month,
          v = seam_v[r],
          tend = seam_tend[r],
          phase = seam_phase[r]
        ),
        k = k,
        min_sep = min_sep,
        max_per_model = max_per_model,
        exclude_models = exclude_models
      )
      nn_cache[[key[r]]] <- nn
    }
    p <- sample.int(nrow(nn), 1L)
    out[, r] <- enso_analog_paths(
      lib,
      nn[p, , drop = FALSE],
      seam_v[r],
      n_months,
      blend_months
    )[1L, ]
  }
  out
}

sample_enso_cfsv2_trajectories <- function(
  cfsv2_data,
  n_rep,
  historical_oni,
  historical_nino34 = NULL,
  sim_months = NULL,
  max_lag = 12L,
  rebase = "seam",
  sd_correction = FALSE,
  target_peak = NULL,
  calibration = NULL,
  calibration_target = "ersst",
  seed = NULL,
  beyond_cfsv2 = c("neutral", "ar1", "analog"),
  ar1_phi = NULL,
  ar1_fit_window = "1991_2020",
  legacy_seam = FALSE,
  max_lead = 0L,
  analog_library = NULL,
  analog_k = 30L,
  analog_blend_months = 3,
  analog_min_sep = 24L,
  analog_max_per_model = Inf,
  analog_exclude_models = character(0),
  analog_seam = ENSO_SEAMS
) {
  beyond_cfsv2 <- match.arg(beyond_cfsv2)
  # analog_seam: "last_cfsv2" (default) seams the
  # analog continuation at the last CFSv2 month, whose ONI is a 2-month mean;
  # "last_full_oni" seams it one month earlier, at the last full 3-month CFSv2
  # mean (enso.R::enso_seam_month()), and the continuation fills the last CFSv2
  # month onward. Only the analog continuation has a seam state.
  analog_seam <- match.arg(analog_seam)
  if (
    !identical(analog_seam, "last_cfsv2") && !identical(beyond_cfsv2, "analog")
  ) {
    stop(
      "sample_enso_cfsv2_trajectories: analog_seam '",
      analog_seam,
      "' needs beyond_cfsv2 = 'analog'.",
      call. = FALSE
    )
  }
  # beyond_cfsv2 = "analog": k-NN analog continuation from the piControl
  # library (continue_enso_analog()); `analog_library` is a path or a library
  # already read by load_enso_continuation_library().
  if (identical(beyond_cfsv2, "analog")) {
    if (isTRUE(legacy_seam)) {
      stop(
        "sample_enso_cfsv2_trajectories: the analog continuation needs legacy_seam = FALSE.",
        call. = FALSE
      )
    }
    if (!is.list(analog_library)) {
      analog_library <- load_enso_continuation_library(analog_library)
    }
  }
  # max_lead > 0 (deployed ENSO-lead fits): also return enso_lead1..max_lead,
  # the ONI k months AFTER each sim month, from the same per-replicate series
  # (published ONI, then CFSv2, then the beyond-CFSv2 continuation). The series
  # is extended by max_lead months so a lead never falls off its end; with
  # max_lead = 0 the output is unchanged.
  max_lead <- as.integer(max_lead %||% 0L)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  stopifnot(
    all(c("target_date", "member", "nino34_anomaly") %in% names(cfsv2_data)),
    all(c("year", "month", "enso_index") %in% names(historical_oni))
  )
  n_modes <- sum(
    isTRUE(sd_correction),
    !is.null(target_peak),
    !is.null(calibration)
  )
  if (n_modes > 1L) {
    stop(
      "sample_enso_cfsv2_trajectories: cfsv2_sd_correction, cfsv2_target_peak and ",
      "cfsv2_calibration are mutually exclusive; set only one."
    )
  }
  if (is.character(calibration)) {
    calibration <- readRDS(calibration)
  }
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("sample_enso_cfsv2_trajectories: `zoo` package required for rollmean.")
  }

  oni_max_date <- with(
    historical_oni,
    as.Date(sprintf(
      "%d-%02d-01",
      max(year),
      month[which.max(year * 12 + month)]
    ))
  )
  cfsv2_data$target_date <- as.Date(cfsv2_data$target_date)
  cfs_target_dates <- sort(unique(cfsv2_data$target_date))

  # If sim_months not given, infer from CFSv2 extent + ONI end.
  if (is.null(sim_months)) {
    all_dates <- sort(unique(c(oni_max_date, cfs_target_dates)))
    sim_months <- seq.Date(all_dates[1], tail(all_dates, 1), by = "1 month")
  }
  M <- length(sim_months)

  # Look-ups for historical ONI and raw Niño 3.4.
  oni_lkp <- stats::setNames(
    historical_oni$enso_index,
    as.Date(sprintf("%d-%02d-01", historical_oni$year, historical_oni$month))
  )
  nino34_lkp <- if (
    !is.null(historical_nino34) &&
      all(c("year", "month", "enso_index") %in% names(historical_nino34))
  ) {
    stats::setNames(
      historical_nino34$enso_index,
      as.Date(sprintf(
        "%d-%02d-01",
        historical_nino34$year,
        historical_nino34$month
      ))
    )
  } else {
    NULL
  }

  # Widen CFSv2 to [target_date × member] of raw Niño 3.4 anomaly.
  cfs_members <- sort(unique(cfsv2_data$member))
  n_members <- length(cfs_members)
  cfs_wide <- matrix(
    NA_real_,
    nrow = length(cfs_target_dates),
    ncol = n_members,
    dimnames = list(as.character(cfs_target_dates), as.character(cfs_members))
  )
  cfs_wide[cbind(
    match(as.character(cfsv2_data$target_date), rownames(cfs_wide)),
    match(as.character(cfsv2_data$member), colnames(cfs_wide))
  )] <- cfsv2_data$nino34_anomaly

  # Historical raw Niño 3.4 lookup (fallback to ONI proxy).
  hist_lkp_for_nino34 <- if (!is.null(nino34_lkp)) nino34_lkp else oni_lkp

  # Crude spread correction:
  # rescale the raw CFSv2 anomaly by k = sd(observed ONI) / sd(raw CFSv2
  # members), applied BEFORE the seam rebase below so the seam-anchor delta is
  # computed against the (already-scaled) ensemble and continuity at the seam
  # is preserved. Without this, a strong-El-Nino CFSv2 issue can push the
  # forecast far outside the fitted crossbasis boundary knots (unbounded
  # cubic extrapolation -> mu blow-up in the disease models).
  #
  # Peak-target correction:
  # the preferred alternative — rescale the raw CFSv2 anomaly by
  # k = target_peak / raw_peak, where raw_peak is the max of the 3-month
  # rolling-mean ONI-equivalent of the ENSEMBLE-MEAN raw forecast over the
  # forecast horizon. Same application point (before the seam rebase) so
  # seam-anchor continuity is preserved. Calibrated to an external outlook
  # (e.g. CPC/IRI) rather than an incidental sd-ratio artefact.
  if (isTRUE(sd_correction)) {
    sigma_obs <- stats::sd(historical_oni$enso_index, na.rm = TRUE)
    sigma_raw <- stats::sd(cfsv2_data$nino34_anomaly, na.rm = TRUE)
    k <- if (is.finite(sigma_obs) && is.finite(sigma_raw) && sigma_raw > 0) {
      sigma_obs / sigma_raw
    } else {
      1.0
    }
    message(sprintf(
      "[sample_enso_cfsv2_trajectories] sd_correction: k = sd(ONI)=%.3f / sd(raw CFSv2)=%.3f = %.3f",
      sigma_obs,
      sigma_raw,
      k
    ))
    cfs_wide <- cfs_wide * k
  } else if (!is.null(target_peak)) {
    # Locate the same seam-anchor date/values the seam-rebase step below will
    # use, so the scale factor below anticipates that additive shift (making
    # the seam-rebase step a ~0 no-op) instead of being computed independently
    # of it and then thrown off by it.
    hist_avail <- names(hist_lkp_for_nino34)[is.finite(hist_lkp_for_nino34)]
    ov <- intersect(rownames(cfs_wide), hist_avail)
    if (length(ov) > 0L) {
      b_key <- ov[which.max(as.Date(ov))]
      cfs_mean_anchor <- mean(cfs_wide[b_key, ], na.rm = TRUE)
      hist_val_anchor <- as.numeric(hist_lkp_for_nino34[b_key])
    } else {
      cfs_mean_anchor <- 0.0
      hist_val_anchor <- 0.0
    }

    # Raw (unscaled) ensemble-mean ONI-equivalent trajectory over the forecast
    # horizon, used to locate the peak.
    fc_mean_anom <- vapply(
      cfs_target_dates,
      function(d) {
        mean(cfs_wide[as.character(d), ], na.rm = TRUE)
      },
      numeric(1)
    )
    peak_series_dates <- seq.Date(
      min(cfs_target_dates) %m-% months(2L),
      max(cfs_target_dates),
      by = "1 month"
    )
    peak_base <- as.numeric(hist_lkp_for_nino34[as.character(
      peak_series_dates
    )])
    is_peak_fc <- peak_series_dates %in% cfs_target_dates
    peak_base[is_peak_fc] <- fc_mean_anom[match(
      peak_series_dates[is_peak_fc],
      cfs_target_dates
    )]
    peak_oni <- as.numeric(zoo::rollmean(
      peak_base,
      k = 3L,
      align = "center",
      fill = NA
    ))
    raw_oni_peak <- suppressWarnings(max(peak_oni[is_peak_fc], na.rm = TRUE))

    # Affine transform anchored at the seam: v' = hist_val_anchor + k * (v -
    # cfs_mean_anchor). This keeps the ensemble mean exactly at hist_val_anchor
    # at the seam date (so the seam-rebase delta below comes out ~0) while
    # scaling the rest of the trajectory (incl. its peak) around that anchor.
    denom <- raw_oni_peak - cfs_mean_anchor
    k <- if (is.finite(raw_oni_peak) && is.finite(denom) && abs(denom) > 1e-8) {
      (target_peak - hist_val_anchor) / denom
    } else {
      1.0
    }
    message(sprintf(
      "[sample_enso_cfsv2_trajectories] target_peak: k = (target(%.3f) - seam(%.3f)) / (raw_oni_peak(%.3f) - seam_raw(%.3f)) = %.3f",
      target_peak,
      hist_val_anchor,
      raw_oni_peak,
      cfs_mean_anchor,
      k
    ))
    cfs_wide <- hist_val_anchor + k * (cfs_wide - cfs_mean_anchor)
  } else if (!is.null(calibration)) {
    # Hindcast-calibrated CPC-style PDF + spread correction (see
    # apply_cfsv2_calibration()); rows must be in date order, which
    # cfs_target_dates (sorted) guarantees.
    message(sprintf(
      "[sample_enso_cfsv2_trajectories] calibration: CPC-style PDF + spread correction (target %s, window %s)",
      calibration_target,
      paste(calibration$meta$window, collapse = "-")
    ))
    cfs_wide <- apply_cfsv2_calibration(
      cfs_wide,
      calibration,
      target = calibration_target
    )
    # The calibrated forecast rows are already on the target's observational
    # scale, so the additive seam rebase below is replaced by anchoring the
    # unobserved analysis rows to the latest published ONI season.
    anc <- .cfsv2_anchor_analysis_rows(
      cfs_wide,
      attr(cfs_wide, "first_forecast_row"),
      hist_lkp_for_nino34,
      oni_lkp
    )
    cfs_wide <- anc$cfs_wide
    message(
      if (is.na(anc$season_centre)) {
        "[sample_enso_cfsv2_trajectories] calibration: observed Nino3.4 covers the analysis rows (no anchoring); seam rebase skipped"
      } else {
        sprintf(
          "[sample_enso_cfsv2_trajectories] calibration: analysis rows anchored to ONI %s (shift %+.3f); seam rebase skipped",
          format(anc$season_centre),
          anc$delta
        )
      }
    )
    rebase <- "none"
  }

  # Seam-anchor rebase: shift the ensemble so its mean matches historical raw
  # Niño 3.4 at the latest month present in both series.
  delta <- 0.0
  if (identical(rebase, "seam")) {
    hist_avail <- names(hist_lkp_for_nino34)[is.finite(hist_lkp_for_nino34)]
    ov <- intersect(rownames(cfs_wide), hist_avail)
    if (length(ov) > 0L) {
      b_key <- ov[which.max(as.Date(ov))]
      cfs_mean <- mean(cfs_wide[b_key, ], na.rm = TRUE)
      hist_val <- as.numeric(hist_lkp_for_nino34[b_key])
      if (is.finite(cfs_mean) && is.finite(hist_val)) {
        delta <- hist_val - cfs_mean
      }
    }
  } else if (!is.null(.rebase_offset(rebase))) {
    # A numeric offset, also as a string ("0.3" from the CLI).
    delta <- .rebase_offset(rebase)
  }
  cfs_wide <- cfs_wide + delta

  # Build the full monthly raw Niño 3.4 series (historical + forecast), then
  # apply the 3-month rolling mean to get the ONI-equivalent per replicate.
  series_start <- min(sim_months) %m-% months(max_lag + 1L)
  series_end <- max(sim_months) %m+% months(1L)
  # Past the CFSv2 span: one more month so the
  # centred mean at the last simulated month is defined from the continuation
  # itself rather than from the series edge. The legacy layout keeps the edge.
  if (!isTRUE(legacy_seam)) {
    series_end <- series_end %m+% months(1L)
  }
  if (max_lead > 0L) {
    series_end <- series_end %m+% months(max_lead)
  }
  series_dates <- seq.Date(series_start, series_end, by = "1 month")
  L_series <- length(series_dates)

  base_series <- as.numeric(hist_lkp_for_nino34[as.character(series_dates)])
  # Seam on the RAW series: every month without a finite raw Niño 3.4 value is
  # filled from CFSv2 (or the neutral fallback below). Not the ONI end: CPC
  # publishes the detrended Niño 3.4 file later than ONI, and an ONI-end seam
  # would leave the months in between NA, which the lag loop would fill with a
  # one-month neutral notch in the forecast ENSO input.
  forecast_idx <- which(!is.finite(base_series))

  # Assign one member per replicate (with replacement), coherent across leads.
  rep_member <- if (n_members > 0L) {
    sample.int(n_members, size = n_rep, replace = TRUE)
  } else {
    integer(0)
  }

  nino34_mat <- matrix(base_series, nrow = L_series, ncol = n_rep)
  for (fi in forecast_idx) {
    d_key <- as.character(series_dates[fi])
    if (d_key %in% rownames(cfs_wide) && n_members > 0L) {
      filled <- as.numeric(cfs_wide[d_key, ])[rep_member]
      filled[!is.finite(filled)] <- 0.0
      nino34_mat[fi, ] <- filled
    } else {
      nino34_mat[fi, ] <- 0.0 # beyond CFSv2 span → neutral fallback
    }
  }

  # Months AFTER the last CFSv2 target (the span end; earlier gaps keep the 0.0
  # fallback above): `after_idx`, and the last span row `last_row`.
  last_cfs <- if (length(cfs_target_dates)) {
    max(cfs_target_dates)
  } else {
    as.Date(NA)
  }
  after_idx <- if (is.na(last_cfs)) {
    integer(0)
  } else {
    intersect(forecast_idx, which(series_dates > last_cfs))
  }
  last_row <- if (length(after_idx)) {
    match(as.character(last_cfs), as.character(series_dates))
  } else {
    NA_integer_
  }
  if (
    !isTRUE(legacy_seam) &&
      length(after_idx) &&
      !is.na(last_row) &&
      identical(beyond_cfsv2, "ar1")
  ) {
    # Each replicate's member decays from its last CFSv2 value towards the
    # climatological 0: x_{t+k} = phi^k x_t. Deterministic (no innovation),
    # so the member spread at the seam is carried, not inflated.
    phi <- ar1_phi %||% fit_oni_ar1(hist_lkp_for_nino34, ar1_fit_window)$phi
    nino34_mat[after_idx, ] <- continue_enso_ar1(
      nino34_mat[last_row, ],
      n_months = length(after_idx),
      phi = phi
    )
    message(sprintf(
      "[sample_enso_cfsv2_trajectories] beyond CFSv2 (%s): AR(1) continuation, phi = %.3f, %d month(s)",
      format(last_cfs),
      phi,
      length(after_idx)
    ))
  }

  # Apply 3-month centered rolling mean → ONI-equivalent per replicate.
  oni_mat <- apply(nino34_mat, 2, function(col_vec) {
    as.numeric(zoo::rollmean(col_vec, k = 3L, align = "center", fill = NA))
  })
  if (
    !isTRUE(legacy_seam) &&
      length(after_idx) &&
      !is.na(last_row) &&
      identical(beyond_cfsv2, "neutral")
  ) {
    # Neutral past the span means ONI = 0 there -- not a centred mean that mixes
    # the last CFSv2 values with zeros (that would pull the last span month one
    # third of the way to neutral: exactly the month each group's min-lag track
    # reads at its horizon end; legacy_seam keeps that behaviour). The last span month's
    # ONI uses the months that exist (m-1, m).
    oni_mat[last_row, ] <- colMeans(nino34_mat[
      c(last_row - 1L, last_row),
      ,
      drop = FALSE
    ])
    oni_mat[after_idx, ] <- 0.0
  }
  if (
    length(after_idx) && !is.na(last_row) && identical(beyond_cfsv2, "analog")
  ) {
    # Analog continuation in ONI space (the library is ONI-scale). The last span
    # month's ONI as for neutral (months m-1, m), then each replicate's seam
    # state -- ONI(m), ONI(m) - ONI(m-3), ONI(m-6), each read as the lag loop
    # below reads it (published ONI first) -- selects its analogs. Drawn after
    # the member draw and only here, so neutral / ar1 keep their RNG stream.
    # analog_seam "last_full_oni": the seam is the row before (a full 3-month
    # CFSv2 mean) and the continuation also replaces the last CFSv2 month.
    full_seam <- identical(analog_seam, "last_full_oni")
    seam_row <- if (full_seam) last_row - 1L else last_row
    cont_idx <- if (full_seam) c(last_row, after_idx) else after_idx
    if (!full_seam) {
      oni_mat[last_row, ] <- colMeans(nino34_mat[
        c(last_row - 1L, last_row),
        ,
        drop = FALSE
      ])
    }
    if (seam_row <= 7L) {
      stop(
        "sample_enso_cfsv2_trajectories: analog continuation needs 6 months before the last ",
        "CFSv2 month inside the series.",
        call. = FALSE
      )
    }
    oni_at <- function(row) {
      sc <- oni_lkp[as.character(series_dates[row])]
      if (!is.na(sc)) rep(as.numeric(sc), n_rep) else oni_mat[row, ]
    }
    seam_v <- oni_at(seam_row)
    oni_mat[cont_idx, ] <- continue_enso_analog(
      analog_library,
      seam_v = seam_v,
      seam_tend = seam_v - oni_at(seam_row - 3L),
      seam_phase = oni_at(seam_row - 6L),
      seam_month = as.integer(format(series_dates[seam_row], "%m")),
      n_months = length(cont_idx),
      k = as.integer(analog_k),
      blend_months = as.numeric(analog_blend_months),
      min_sep = as.integer(analog_min_sep),
      max_per_model = as.numeric(analog_max_per_model %||% Inf),
      exclude_models = as.character(unlist(
        analog_exclude_models %||% character(0)
      ))
    )
    message(sprintf(
      "[sample_enso_cfsv2_trajectories] beyond CFSv2 (%s): analog continuation, k = %d, blend %g mo, %d month(s), library %s",
      format(last_cfs),
      as.integer(analog_k),
      as.numeric(analog_blend_months),
      length(after_idx),
      analog_library$meta$digest %||% "?"
    ))
  }
  valid_idx <- which(!is.na(oni_mat[, 1]))
  oni_dates <- series_dates[valid_idx]
  oni_mat_trim <- oni_mat[valid_idx, , drop = FALSE]

  # Pre-lag: look up ONI-equivalent per (sim_month - lag).
  result <- vector("list", max_lag + 1L)
  names(result) <- paste0("enso_lag", seq(0L, max_lag))

  for (lag in seq(0L, max_lag)) {
    lag_dates <- sim_months %m-% months(lag)
    lag_mat <- matrix(NA_real_, nrow = M, ncol = n_rep)
    for (mi in seq_len(M)) {
      d_key <- as.character(lag_dates[mi])
      # Published ONI first wherever it exists (it is what the models were
      # fit on); the per-replicate rolling mean of raw/CFSv2 Niño 3.4 for
      # every later month; 0.0 only before the record starts.
      scalar <- oni_lkp[d_key]
      idx <- match(d_key, as.character(oni_dates))
      if (!is.na(scalar)) {
        lag_mat[mi, ] <- scalar
      } else if (!is.na(idx)) {
        lag_mat[mi, ] <- oni_mat_trim[idx, ]
      } else {
        lag_mat[mi, ] <- 0.0
      }
    }
    result[[paste0("enso_lag", lag)]] <- lag_mat
  }
  # Leads (max_lead > 0): the same lookup at sim_month + k. Past the series end
  # (cannot happen after the extension above) is an error, not a silent 0.0.
  for (k in seq_len(max_lead)) {
    lead_dates <- sim_months %m+% months(k)
    lead_mat <- matrix(NA_real_, nrow = M, ncol = n_rep)
    for (mi in seq_len(M)) {
      d_key <- as.character(lead_dates[mi])
      scalar <- oni_lkp[d_key]
      idx <- match(d_key, as.character(oni_dates))
      if (!is.na(scalar)) {
        lead_mat[mi, ] <- scalar
      } else if (!is.na(idx)) {
        lead_mat[mi, ] <- oni_mat_trim[idx, ]
      } else {
        stop(
          sprintf(
            "sample_enso_cfsv2_trajectories: no ONI for lead %d at %s (series end %s).",
            k,
            d_key,
            format(series_end)
          ),
          call. = FALSE
        )
      }
    }
    result[[enso_lag_col_name(-k)]] <- lead_mat
  }

  result
}


#' Sample deterministic ENSO trajectories from a historical analog period.
#'
#' Replaces forecast-period ENSO (from the first element of \code{sim_months}
#' onward) with values cyclically drawn from a specified past ONI window.  All
#' replicates receive identical values — the analog scenario is deterministic.
#' Pre-horizon lag dates (i.e. \code{sim_month - lag < min(sim_months)}) are
#' filled from the historical ONI cache unchanged.
#'
#' @param historical_oni  Tibble from \code{fetch_enso_index()}: columns
#'   \code{year}, \code{month}, \code{enso_index}.  Must cover all
#'   pre-horizon lag dates needed by the largest lag.
#' @param sim_months      Date vector (first-of-month) of simulation months.
#'   The minimum element is treated as the horizon start; dates before it
#'   arise only from lag look-ups and are served from \code{historical_oni}.
#' @param analog_start    Date (or "YYYY-MM-DD" string): first month of the
#'   analog period to replay (inclusive).
#' @param analog_end      Date (or "YYYY-MM-DD" string): last month of the
#'   analog period (inclusive).
#' @param max_lag         Maximum ENSO lag to pre-compute (default 12L).
#' @param n_rep           Number of replicates.  All columns of the output
#'   matrices are identical (default 1L).
#' @return Named list of [M \eqn{\times} n_rep] numeric matrices, one per
#'   lag name ("enso_lag0", …, "enso_lag{max_lag}").
sample_enso_analog_scenarios <- function(
  historical_oni,
  sim_months,
  analog_start,
  analog_end,
  max_lag = 12L,
  n_rep = 1L,
  max_lead = 0L
) {
  max_lead <- as.integer(max_lead %||% 0L) # > 0: also enso_lead1..max_lead (cyclic analog months)
  stopifnot(all(c("year", "month", "enso_index") %in% names(historical_oni)))
  analog_start <- as.Date(analog_start)
  analog_end <- as.Date(analog_end)
  sim_months <- sort(as.Date(sim_months))
  horizon_start <- sim_months[1L]

  # Historical ONI look-up: Date string → scalar.
  oni_lkp <- stats::setNames(
    historical_oni$enso_index,
    as.character(as.Date(sprintf(
      "%d-%02d-01",
      historical_oni$year,
      historical_oni$month
    )))
  )

  analog_months_vec <- seq.Date(analog_start, analog_end, by = "1 month")
  analog_length <- length(analog_months_vec)
  if (analog_length < 1L) {
    stop("sample_enso_analog_scenarios: analog_end must be >= analog_start.")
  }

  # Per-date ONI look-up.
  # d < horizon_start → historical ONI from cache (unchanged).
  # d >= horizon_start → cyclically mapped onto the analog period.
  #   Offset k = (year(d) - year(horizon_start)) * 12 + month(d) - month(horizon_start)
  #   analog index = (k %% analog_length) + 1
  lookup_oni <- function(d) {
    d <- as.Date(d)
    if (d < horizon_start) {
      val <- oni_lkp[as.character(d)]
      return(if (!is.na(val)) val else 0.0)
    }
    yr <- as.integer(format(d, "%Y"))
    mo <- as.integer(format(d, "%m"))
    yr0 <- as.integer(format(horizon_start, "%Y"))
    mo0 <- as.integer(format(horizon_start, "%m"))
    k <- (yr - yr0) * 12L + (mo - mo0)
    analog_idx <- (k %% analog_length) + 1L
    val <- oni_lkp[as.character(analog_months_vec[analog_idx])]
    if (!is.na(val)) val else 0.0
  }

  M <- length(sim_months)
  result <- vector("list", max_lag + 1L)
  names(result) <- paste0("enso_lag", seq(0L, max_lag))

  for (lag in c(seq(0L, max_lag), -seq_len(max_lead))) {
    # leads (negative) after the lags
    lag_dates <- sim_months %m-% months(lag)
    oni_vals <- vapply(lag_dates, lookup_oni, numeric(1))
    # All replicates are identical (deterministic counterfactual).
    result[[enso_lag_col_name(lag)]] <- matrix(
      rep(oni_vals, n_rep),
      nrow = M,
      ncol = n_rep
    )
  }
  result
}


## ---------------------------------------------------------------------------
## ENSO analog sets (Lenssen counterfactual pool subsetting)
## ---------------------------------------------------------------------------

#' Valid Lenssen analog-set names.
#'
#' Single source of truth for validating `simulator.enso.counterfactual.analog_set`.
#' `"all"` is the full pool (default, and byte-identical to the pre-set behaviour);
#' the other three are ENSO-state subsets defined by [classify_enso_analog_sets()].
ENSO_ANALOG_SETS <- c("all", "nina", "neutral", "elnino")


#' Classify Lenssen analog trajectories into ENSO-state sets.
#'
#' Labels each row of an analog trajectory matrix as \code{"nina"},
#' \code{"neutral"} or \code{"elnino"} based on that trajectory's mean ONI over a
#' target peak season (default DJF).
#'
#' \strong{Why not lead-0?} The analog pool is built by matching every trajectory's
#' lead-0 state to the \emph{same} observed ONI at the handoff month, so the lead-0
#' spread across the pool is tiny (sd ~0.11) and carries no discriminating signal.
#' The sets must therefore be defined on how each trajectory \emph{evolves} — the
#' first boreal-winter peak inside the horizon is when ENSO actually peaks and when
#' the factual (CFSv2) signal is strongest, making it the natural contrast point.
#'
#' \strong{Why CPC thresholds and not empirical terciles?} With the operational CPC
#' cut-offs (\eqn{\pm}0.5 °C) the labels mean what "La Niña"/"El Niño" mean in the
#' literature and are directly comparable to the factual arm's ONI values.  Terciles
#' would be pool-relative and would silently redefine themselves on every rebuild.
#'
#' Base R only (no dplyr/lubridate) so it can be sourced and verified standalone.
#'
#' @param traj_mat         Numeric matrix [n_analogs x n_months] of monthly ONI
#'   values; column 1 is \code{horizon_start}.
#' @param horizon_start    Date (or coercible) of the trajectory's first column.
#' @param peak_months      Integer calendar months defining the peak season, in
#'   order.  Default \code{c(12L, 1L, 2L)} (DJF).
#' @param peak_season_year Optional integer: calendar year of the \emph{first}
#'   month of the target season (e.g. 2026 for DJF 2026/27).  \code{NULL} (default)
#'   takes the first complete occurrence inside the horizon.
#' @param thresholds       Named numeric of length 2: \code{nina} (upper bound,
#'   inclusive) and \code{elnino} (lower bound, inclusive).  Default CPC
#'   \code{c(nina = -0.5, elnino = 0.5)}.
#' @param na_action        \code{"error"} (default) to fail if any selected cell is
#'   NA; \code{"na"} to propagate NA labels.
#' @return List with \code{labels}, \code{stat}, \code{stat_name}, \code{columns},
#'   \code{window}, \code{thresholds}, \code{counts}, \code{definition}.
classify_enso_analog_sets <- function(
  traj_mat,
  horizon_start,
  peak_months = c(12L, 1L, 2L),
  peak_season_year = NULL,
  thresholds = ENSO_PHASE_THRESHOLDS,
  na_action = c("error", "na")
) {
  na_action <- match.arg(na_action)

  if (!is.matrix(traj_mat) || !is.numeric(traj_mat)) {
    stop(
      "classify_enso_analog_sets: `traj_mat` must be a numeric matrix [n_analogs x n_months]."
    )
  }
  if (nrow(traj_mat) < 1L || ncol(traj_mat) < 1L) {
    stop(
      "classify_enso_analog_sets: `traj_mat` must have at least one row and one column."
    )
  }

  peak_months <- as.integer(peak_months)
  if (
    length(peak_months) < 1L ||
      any(is.na(peak_months)) ||
      any(peak_months < 1L) ||
      any(peak_months > 12L)
  ) {
    stop("classify_enso_analog_sets: `peak_months` must be integers in 1..12.")
  }

  if (!all(c("nina", "elnino") %in% names(thresholds))) {
    stop(
      "classify_enso_analog_sets: `thresholds` must be a named numeric with ",
      "elements 'nina' and 'elnino'."
    )
  }
  if (!is.finite(thresholds[["nina"]]) || !is.finite(thresholds[["elnino"]])) {
    stop("classify_enso_analog_sets: `thresholds` must be finite.")
  }
  if (thresholds[["nina"]] >= thresholds[["elnino"]]) {
    stop(sprintf(
      paste0(
        "classify_enso_analog_sets: thresholds must satisfy nina < elnino, got ",
        "nina=%s, elnino=%s."
      ),
      format(thresholds[["nina"]]),
      format(thresholds[["elnino"]])
    ))
  }

  horizon_start <- as.Date(horizon_start)
  if (length(horizon_start) != 1L || is.na(horizon_start)) {
    stop(
      "classify_enso_analog_sets: `horizon_start` must be a single non-NA date."
    )
  }

  # Month grid for the trajectory columns.  Derived from horizon_start rather than
  # from column names so the classifier does not depend on cache dimnames.
  grid <- seq(horizon_start, by = "1 month", length.out = ncol(traj_mat))
  grid_month <- as.integer(format(grid, "%m"))
  grid_year <- as.integer(format(grid, "%Y"))

  # Locate the first (or requested) complete, contiguous occurrence of peak_months.
  L <- length(peak_months)
  n_cols <- ncol(traj_mat)
  starts <- if (n_cols >= L) seq_len(n_cols - L + 1L) else integer(0)
  is_match <- vapply(
    starts,
    function(i) {
      identical(grid_month[i:(i + L - 1L)], peak_months)
    },
    logical(1)
  )
  cand <- starts[is_match]

  if (!is.null(peak_season_year)) {
    peak_season_year <- as.integer(peak_season_year)
    cand <- cand[grid_year[cand] == peak_season_year]
  }

  season_lab <- paste(substr(month.abb[peak_months], 1L, 1L), collapse = "")
  if (length(cand) < 1L) {
    stop(sprintf(
      paste0(
        "classify_enso_analog_sets: peak season %s%s is absent or truncated in the ",
        "trajectory horizon (%s .. %s, %d months). Refusing to classify on a partial ",
        "season -- that would produce plausible-looking but wrong sets. Either rebuild ",
        "the cache with a longer --lead_months, or pass `peak_months` for a season that ",
        "fits inside the horizon."
      ),
      season_lab,
      if (is.null(peak_season_year)) "" else paste0(" ", peak_season_year),
      format(grid[1L]),
      format(grid[n_cols]),
      n_cols
    ))
  }

  col_start <- cand[1L]
  columns <- col_start:(col_start + L - 1L)
  window <- grid[columns]

  sub <- traj_mat[, columns, drop = FALSE]
  if (na_action == "error" && anyNA(sub)) {
    stop(sprintf(
      paste0(
        "classify_enso_analog_sets: %d NA value(s) in the peak-season columns (%s). ",
        "Pass na_action = 'na' to propagate NA labels instead of failing."
      ),
      sum(is.na(sub)),
      paste(columns, collapse = ",")
    ))
  }

  stat <- rowMeans(sub, na.rm = FALSE)
  stat_name <- paste0("peak_", tolower(season_lab))

  # Threshold phase of the statistic "djf_mean_oni": the mean of the
  # trajectory's monthly 3-month ONI values over the season, i.e. a 5-month
  # triangular mean of Nino 3.4 -- not ONI(DJF).
  labels <- classify_enso_phase(stat, "threshold", thresholds = thresholds)

  counts <- vapply(
    setdiff(ENSO_ANALOG_SETS, "all"),
    function(s) sum(labels == s, na.rm = TRUE),
    integer(1)
  )

  list(
    labels = labels,
    stat = stat,
    stat_name = stat_name,
    columns = as.integer(columns),
    window = window,
    thresholds = thresholds,
    counts = counts,
    definition = list(
      stat_name = stat_name,
      peak_months = peak_months,
      peak_season = season_lab,
      peak_season_year = as.integer(grid_year[col_start]),
      horizon_start = format(horizon_start),
      window_start = format(window[1L]),
      window_end = format(window[L]),
      columns = as.integer(columns),
      thresholds = as.list(thresholds),
      n_total = length(stat),
      counts = as.list(counts),
      n_na = sum(is.na(labels))
    )
  )
}


#' Resolve an analog-set name to trajectory row indices in a Lenssen cache.
#'
#' @param cache      List read from the Lenssen analog cache RDS; must contain
#'   \code{oni_trajectories} and, for any set other than \code{"all"},
#'   \code{set_labels}.
#' @param analog_set One of [ENSO_ANALOG_SETS].
#' @return Integer vector of row indices into \code{cache$oni_trajectories}.
enso_analog_set_rows <- function(cache, analog_set = "all") {
  analog_set <- as.character(analog_set)
  if (length(analog_set) != 1L || !analog_set %in% ENSO_ANALOG_SETS) {
    stop(sprintf(
      "enso_analog_set_rows: unknown analog_set '%s'; must be one of: %s.",
      paste(analog_set, collapse = ","),
      paste(ENSO_ANALOG_SETS, collapse = ", ")
    ))
  }

  if (is.null(cache$oni_trajectories)) {
    stop("enso_analog_set_rows: cache has no `oni_trajectories`.")
  }
  n_rows <- nrow(cache$oni_trajectories)

  if (identical(analog_set, "all")) {
    return(seq_len(n_rows))
  }

  labels <- cache$set_labels
  if (is.null(labels)) {
    stop(sprintf(
      paste0(
        "enso_analog_set_rows: analog_set='%s' requested but the cache has no ",
        "`set_labels` -- rebuild with `Rscript analysis/03_prepare_lenssen_analogs.R ",
        "--horizon_start <cache handoff>`."
      ),
      analog_set
    ))
  }
  if (length(labels) != n_rows) {
    stop(sprintf(
      paste0(
        "enso_analog_set_rows: cache is inconsistent -- %d set_labels for %d ",
        "trajectories. Rebuild the analog cache."
      ),
      length(labels),
      n_rows
    ))
  }

  rows <- which(labels == analog_set)
  if (length(rows) < 1L) {
    stop(sprintf(
      paste0(
        "enso_analog_set_rows: analog_set='%s' selects 0 of %d trajectories. ",
        "Available sets: %s."
      ),
      analog_set,
      n_rows,
      paste(
        sprintf("%s=%d", names(table(labels)), as.integer(table(labels))),
        collapse = ", "
      )
    ))
  }
  rows
}


#' Sample per-replicate ENSO trajectories from the Lenssen piControl analog cache.
#'
#' Implements the Scenario B (Lenssen) counterfactual: each simulation replicate
#' receives an independently sampled monthly ONI trajectory drawn from a pool of
#' pre-industrial-control (piControl) ENSO series assembled from 11 CMIP models
#' (Lenssen et al. 2024, \emph{GRL}, doi 10.1029/2023GL106988; data at
#' zenodo.org/records/10045687).  The piControl runs use 1850 forcing, providing
#' a physically defensible "no-warming world" counterfactual.
#'
#' Pre-horizon lag dates (i.e. \code{sim_month - lag < min(sim_months)}) are
#' filled from the historical ONI cache identical to the factual run; only
#' forecast-period dates draw from the piControl trajectories.
#'
#' \strong{Cache format} (\code{generated_data/lenssen_nino34_analogs.rds}):
#' \describe{
#'   \item{\code{oni_trajectories}}{Numeric matrix [n_analogs × n_months].  Each
#'     row is one trajectory segment of monthly ONI values (3-month centered
#'     rolling mean of Niño 3.4 already applied).  The trajectory must be long
#'     enough to cover all sim months (\code{n_months >= length(sim_months)}).}
#'   \item{\code{n_months}}{Integer: number of columns (months per trajectory).}
#'   \item{\code{n_analogs}}{Integer: number of available trajectories (rows).}
#'   \item{\code{set_labels}}{Character vector [n_analogs], row-aligned with
#'     \code{oni_trajectories}: the ENSO-state set each trajectory belongs to
#'     (see [classify_enso_analog_sets()]).  Required only when \code{analog_set}
#'     is not \code{"all"}.}
#'   \item{\code{meta}}{List with \code{source}, \code{models}, \code{created}.}
#' }
#' The cache is produced by the deferred Lenssen data-processing step
#' (download + analog matching + per-model variance rescaling); the cache file
#' must exist before running this scenario.
#'
#' @param cache_path      Path to the Lenssen analog cache RDS.  Default:
#'   \code{"generated_data/lenssen_nino34_analogs.rds"}.
#' @param sim_months      Date vector (first-of-month) of simulation months.
#'   \code{min(sim_months)} is treated as the horizon start.
#' @param n_rep           Number of simulation replicates (trajectories) to draw.
#' @param historical_oni  Tibble from \code{fetch_enso_index()}: columns
#'   \code{year}, \code{month}, \code{enso_index}.  Used for pre-horizon lags.
#' @param max_lag         Maximum ENSO lag to pre-compute (default 12L).
#' @param mean_offset     Additive shift applied to all trajectory values after
#'   variance scaling (default 0.0).  Use to correct the ~1°C cooler 1850
#'   mean state if needed.
#' @param log_var_scaling Log-scale variance adjustment: trajectories are
#'   multiplied by \code{exp(0.5 * log_var_scaling)} before the mean offset is
#'   added (default 0.0, i.e. no rescaling).
#' @param seed            Optional integer seed for trajectory sampling.
#' @param analog_set      One of [ENSO_ANALOG_SETS].  \code{"all"} (default) samples
#'   the whole pool and is byte-identical to the pre-set behaviour; the other names
#'   restrict sampling to that ENSO-state subset, yielding a counterfactual arm
#'   conditioned on an ENSO state rather than marginalised over it.  Appended last
#'   so existing positional call sites are unaffected.
#' @return Named list of [M \eqn{\times} n_rep] numeric matrices, one per
#'   lag name ("enso_lag0", …, "enso_lag{max_lag}").
sample_enso_lenssen_scenarios <- function(
  cache_path = "generated_data/lenssen_nino34_analogs.rds",
  sim_months,
  n_rep,
  historical_oni,
  max_lag = 12L,
  horizon_start = NULL,
  mean_offset = 0.0,
  log_var_scaling = 0.0,
  seed = NULL,
  analog_set = "all",
  max_lead = 0L
) {
  max_lead <- as.integer(max_lead %||% 0L) # > 0: also enso_lead1..max_lead
  if (!file.exists(cache_path)) {
    stop(sprintf(
      paste0(
        "sample_enso_lenssen_scenarios: Lenssen piControl analog cache not found: '%s'. ",
        "Run the Lenssen data-processing step (download zenodo/10045687, analog ",
        "matching, per-model variance rescaling) to produce this file before running ",
        "the Lenssen counterfactual scenario."
      ),
      cache_path
    ))
  }

  cache <- readRDS(cache_path)
  required_fields <- c("oni_trajectories", "n_months", "n_analogs")
  missing_fields <- setdiff(required_fields, names(cache))
  if (length(missing_fields) > 0L) {
    stop(
      "sample_enso_lenssen_scenarios: cache missing fields: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  traj_mat <- cache$oni_trajectories # [n_analogs × n_months]
  n_analogs <- cache$n_analogs
  n_months <- cache$n_months

  sim_months <- sort(as.Date(sim_months))
  M <- length(sim_months)
  horizon_start <- if (!is.null(horizon_start)) {
    as.Date(horizon_start)
  } else {
    sim_months[1L]
  }

  # Guard: the cache is built for a specific handoff (03_prepare_lenssen_analogs.R
  # selects same-calendar-month analogs whose lead-0 state matches the observed
  # ONI at the cache's horizon_start; column `traj_col` is then read as
  # horizon_start + (traj_col-1) months here).  If the run's horizon_start differs
  # from the cache's, every column is served under the wrong calendar month AND
  # anchored to the wrong observed ONI, producing a discontinuous jump at the
  # observed->scenario boundary.  Fail loudly rather than silently mis-align.
  cache_meta <- cache$meta
  if (!is.null(cache_meta) && !is.null(cache_meta$horizon_start)) {
    cache_horizon <- as.Date(cache_meta$horizon_start)
    if (!identical(cache_horizon, horizon_start)) {
      stop(sprintf(
        paste0(
          "sample_enso_lenssen_scenarios: cache handoff (%s) != requested horizon_start (%s). ",
          "The cache's column 1 is anchored to the observed ONI at the cache handoff ",
          "month, so consuming it at a different handoff mis-aligns the seasonal phase ",
          "and injects a level jump at the observed->scenario boundary. Rebuild with ",
          "`Rscript analysis/03_prepare_lenssen_analogs.R --horizon_start %s` or set ",
          "enso.counterfactual.lenssen_horizon_start to the cache handoff."
        ),
        format(cache_horizon),
        format(horizon_start),
        format(horizon_start)
      ))
    }
  } else {
    warning(
      "sample_enso_lenssen_scenarios: cache has no meta$horizon_start; cannot verify ",
      "the handoff month matches the cache anchor. Ensure the cache was built with ",
      "--horizon_start = ",
      format(horizon_start),
      "."
    )
  }

  if (n_months < M) {
    stop(sprintf(
      paste0(
        "sample_enso_lenssen_scenarios: each trajectory in the cache covers %d months ",
        "but %d simulation months are required.  Rebuild the cache with longer ",
        "trajectory segments."
      ),
      n_months,
      M
    ))
  }

  # Restrict the pool to one ENSO-state set.  Must happen BEFORE the affine
  # transform and, critically, before sample.int() below -- otherwise rows would be
  # drawn from the full pool and the set would have no effect.
  analog_set <- as.character(analog_set)
  if (
    length(analog_set) != 1L ||
      is.na(analog_set) ||
      !analog_set %in% ENSO_ANALOG_SETS
  ) {
    stop(sprintf(
      paste0(
        "sample_enso_lenssen_scenarios: unknown analog_set '%s'; must be one of: %s. ",
        "Set simulator.enso.counterfactual.analog_set in the run config."
      ),
      paste(analog_set, collapse = ","),
      paste(ENSO_ANALOG_SETS, collapse = ", ")
    ))
  }

  if (!identical(as.integer(nrow(traj_mat)), as.integer(n_analogs))) {
    stop(sprintf(
      paste0(
        "sample_enso_lenssen_scenarios: cache is inconsistent -- n_analogs=%d but ",
        "oni_trajectories has %d rows. Rebuild the analog cache."
      ),
      as.integer(n_analogs),
      nrow(traj_mat)
    ))
  }

  if (!identical(analog_set, "all")) {
    set_rows <- enso_analog_set_rows(cache, analog_set)
    traj_mat <- traj_mat[set_rows, , drop = FALSE]
    n_analogs <- nrow(traj_mat)
    message(sprintf(
      "   ENSO analog_set='%s': %d of %d analogs, n_rep=%d -> %.1fx trajectory reuse",
      analog_set,
      n_analogs,
      length(cache$set_labels),
      n_rep,
      n_rep / n_analogs
    ))
    # The effective Monte-Carlo sample size for ENSO-driven quantities is bounded by
    # n_analogs, not by n_rep: replicates that share a trajectory are not independent
    # draws of the ENSO state.  Report per-set intervals with that ceiling in mind.
  } else {
    message(sprintf(
      "   ENSO analog_set='all': full pool, %d analogs, n_rep=%d -> %.1fx trajectory reuse",
      n_analogs,
      n_rep,
      n_rep / n_analogs
    ))
  }

  # Apply mean-state offset and log-variance rescaling.
  if (mean_offset != 0.0 || log_var_scaling != 0.0) {
    traj_mat <- traj_mat * exp(0.5 * log_var_scaling) + mean_offset
  }

  # Sample n_rep trajectory row indices with replacement.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rep_idx <- sample.int(n_analogs, size = n_rep, replace = TRUE)

  # Historical ONI look-up (for pre-horizon lag dates; same for all replicates).
  oni_lkp <- stats::setNames(
    historical_oni$enso_index,
    as.character(as.Date(sprintf(
      "%d-%02d-01",
      historical_oni$year,
      historical_oni$month
    )))
  )

  result <- vector("list", max_lag + 1L)
  names(result) <- paste0("enso_lag", seq(0L, max_lag))

  for (lag in c(seq(0L, max_lag), -seq_len(max_lead))) {
    # leads (negative) after the lags
    lag_dates <- sim_months %m-% months(lag)
    lag_mat <- matrix(NA_real_, nrow = M, ncol = n_rep)
    for (mi in seq_len(M)) {
      d <- lag_dates[mi]
      if (d < horizon_start) {
        # Historical period: identical across all replicates.
        val <- oni_lkp[as.character(d)]
        lag_mat[mi, ] <- if (!is.na(val)) val else 0.0
      } else {
        # Forecast period: each replicate uses its sampled trajectory.
        yr <- as.integer(format(d, "%Y"))
        mo <- as.integer(format(d, "%m"))
        yr0 <- as.integer(format(horizon_start, "%Y"))
        mo0 <- as.integer(format(horizon_start, "%m"))
        traj_col <- (yr - yr0) * 12L + (mo - mo0) + 1L # 1-indexed
        if (traj_col < 1L || traj_col > n_months) {
          # A lead past the cached trajectory would silently read neutral.
          if (lag < 0L && traj_col > n_months) {
            stop(
              sprintf(
                paste0(
                  "sample_enso_lenssen_scenarios: lead %d at %s is past the cached trajectories ",
                  "(%d months from %s); rebuild the cache with --lead_months >= %d."
                ),
                -lag,
                format(d),
                n_months,
                format(horizon_start),
                traj_col
              ),
              call. = FALSE
            )
          }
          lag_mat[mi, ] <- 0.0
        } else {
          lag_mat[mi, ] <- traj_mat[rep_idx, traj_col]
        }
      }
    }
    result[[enso_lag_col_name(lag)]] <- lag_mat
  }

  result
}

#' Latest observed ONI month (first day of the month): the ENSO handoff the
#' forecast and counterfactual arms share. Read from the cached ONI series
#' only (`refresh = FALSE`), so it assumes the cache the simulator runs used is
#' still current.
latest_oni_month <- function(cache_dir = "generated_data/") {
  oni <- fetch_enso_index(
    source = "oni",
    cache_dir = cache_dir,
    refresh = FALSE
  )
  with(
    oni,
    as.Date(sprintf(
      "%d-%02d-01",
      max(year),
      month[which.max(year * 12 + month)]
    ))
  )
}

#' ENSO seams: where the forecast ENSO hands over to the continuation.
#' "last_cfsv2" (default, legacy) = the last CFSv2 month; "last_full_oni" = the
#' month before it, the last month whose ONI is a full 3-month CFSv2 mean (the
#' last CFSv2 month's ONI is a 2-month mean).
ENSO_SEAMS <- c("last_cfsv2", "last_full_oni")

#' enso_seam_month
#'
#' The last month whose ONI owes nothing to the continuation, under `seam`:
#' the last CFSv2 month ("last_cfsv2") or the month before it
#' ("last_full_oni"). The firm forecast horizon of a group is this month + its
#' minimum deployed lag (utils/make_forecast_configs.R), and the
#' CFSv2 sampler's analog seam sits here (sample_enso_cfsv2_trajectories(
#' analog_seam =)).
#' @param last_cfsv2_month Date (or coercible), first of month
#' @param seam one of ENSO_SEAMS
#' @return Date
enso_seam_month <- function(last_cfsv2_month, seam = ENSO_SEAMS) {
  seam <- match.arg(seam)
  if (identical(seam, "last_cfsv2")) {
    return(as.Date(last_cfsv2_month))
  }
  m <- as.Date(format(as.Date(last_cfsv2_month), "%Y-%m-01"))
  seq(m, by = "-1 month", length.out = 2L)[2L]
}

#' enso_event_launch
#'
#' Launch month of a past ENSO event under the event-attribution crossing rule
#' (analysis/configs/enso_events_v2.yml, launch_rule: crossing): walk back from
#' the peak month while ONI stays in the event's phase (>= threshold for El
#' Nino, <= -threshold for La Nina); the first month of that run is the
#' crossing, and the launch is the month before it.
#' @param oni   data frame with `date` (Date, first of month) and `enso_index`
#'   (fetch_enso_index()).
#' @param peak  peak month (Date or "YYYY-MM-DD").
#' @param phase "elnino" | "lanina".
#' @param threshold the ONI threshold (default the CPC 0.5).
#' @return the launch month (Date).
enso_event_launch <- function(
  oni,
  peak,
  phase = c("elnino", "lanina"),
  threshold = ENSO_PHASE_THRESHOLDS[["elnino"]]
) {
  phase <- match.arg(phase)
  peak <- as.Date(peak)
  o <- oni[order(oni$date), c("date", "enso_index")]
  # Threshold phase of the monthly ONI (statistic "event_run"); the registry's
  # "lanina" is classify_enso_phase()'s "nina".
  code <- if (identical(phase, "elnino")) "elnino" else "nina"
  in_phase <- classify_enso_phase(
    o$enso_index,
    "threshold",
    thresholds = c(nina = -threshold, elnino = threshold)
  ) ==
    code
  i <- match(peak, as.Date(o$date))
  if (is.na(i)) {
    stop(
      "enso_event_launch: peak month ",
      format(peak),
      " is not in the ONI series.",
      call. = FALSE
    )
  }
  if (!isTRUE(in_phase[i])) {
    stop(
      sprintf(
        "enso_event_launch: ONI at the peak %s (%.2f) is not in the %s phase.",
        format(peak),
        o$enso_index[i],
        phase
      ),
      call. = FALSE
    )
  }
  while (i > 1L && isTRUE(in_phase[i - 1L])) {
    i <- i - 1L
  }
  if (i == 1L) {
    stop(
      "enso_event_launch: the event's run starts at the first ONI month; no launch before it.",
      call. = FALSE
    )
  }
  as.Date(o$date[i - 1L])
}

# ---------------------------------------------------------------------------
# Lenssen piControl counterfactual cache
#
# The legacy builder's helpers for analysis/03_prepare_lenssen_analogs.R
# (tests/test_enso_counterfactual_match.R pins them against frozen copies), and the model
# list / Rda extraction shared with 03_prepare_enso_continuation_library.R.
# ---------------------------------------------------------------------------

#' The 11 Lenssen et al. (2024) piControl models (zip names of Zenodo record
#' 10045687; not all CMIP6: several are NMME / CMIP5-era runs).
LENSSEN_MODELS <- c(
  "CanESM5",
  "CCSM4-NMME",
  "CESM1-1-CAM5-CMIP5",
  "CESM1-NMME",
  "CESM2",
  "CM2.1-NMME",
  "CM2.5-NMME",
  "CM4",
  "ESM4",
  "GISS-E21G",
  "MIROC6"
)

#' lenssen_extract_rda
#'
#' Extract step01_InitVerifInfo.Rda from a model zip into
#' `<extract_base>/<model>/` (the legacy builder's extract_rda(); the
#' continuation library builder uses extract_base = <zenodo_dir>/extracted).
#' Reuses an already extracted file unless `refresh`.
lenssen_extract_rda <- function(
  model,
  zenodo_dir,
  extract_base = file.path(zenodo_dir, "extracted"),
  refresh = FALSE
) {
  zip_path <- file.path(zenodo_dir, paste0(model, ".zip"))
  model_dir <- file.path(extract_base, model)
  out_path <- file.path(model_dir, "step01_InitVerifInfo.Rda")

  if (file.exists(out_path) && !refresh) {
    return(out_path)
  }

  if (!file.exists(zip_path)) {
    stop(
      "Zip not found: ",
      zip_path,
      "\nCopy or symlink the Lenssen Zenodo archive files to ",
      zenodo_dir
    )
  }

  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  zip_entry <- paste0(
    model,
    "/experiments/tos_zos_NEW/step01_InitVerifInfo.Rda"
  )
  unzip(zip_path, files = zip_entry, exdir = model_dir, junkpaths = TRUE)

  if (!file.exists(out_path)) {
    stop("Extraction failed for ", model, ": expected ", out_path)
  }

  message("   extracted: ", out_path)
  out_path
}

#' Legacy matcher (03_prepare_lenssen_analogs.R select_analogs(), moved
#' verbatim): best-N analog init indices (piControl month positions, 1-based)
#' by absolute distance of the RAW model ONI to `target_oni` in the same
#' calendar month, greedy within-model minimum separation.
.lenssen_select_analogs_legacy <- function(
  nino34,
  timeMap,
  target_oni,
  target_month,
  n_analogs,
  min_sep,
  lead_months
) {
  # Same calendar-month filter (seasonal consistency)
  same_season_idx <- which(timeMap[, 2L] == target_month)

  # Exclude positions within lead_months of the end (forward trajectory would overflow)
  last_valid_idx <- length(nino34) - lead_months
  same_season_idx <- same_season_idx[same_season_idx <= last_valid_idx]

  if (length(same_season_idx) == 0L) {
    warning(
      "No valid analog months for model (target_month=",
      target_month,
      ")"
    )
    return(integer(0))
  }

  # Rank by absolute ONI distance
  dist_vec <- abs(nino34[same_season_idx] - target_oni)
  ranked <- same_season_idx[order(dist_vec)]

  # Greedy minimum-separation selection
  # Note: for same-season (yearly) matching, min_sep=12 allows adjacent years
  # (consecutive Mays are exactly 12 months apart, so abs diff >= 12 is TRUE).
  # Increase min_sep to 13 to require at least one gap year.
  selected <- integer(0)
  for (idx in ranked) {
    if (length(selected) >= n_analogs) {
      break
    }
    if (length(selected) == 0L || all(abs(idx - selected) >= min_sep)) {
      selected <- c(selected, idx)
    }
  }
  selected
}

#' Legacy forward Nino3.4 trajectories from init indices (0..lead_months);
#' The legacy builder's extract_trajectories().
.lenssen_extract_trajectories_legacy <- function(
  nino34,
  init_indices,
  lead_months
) {
  purrr::map(init_indices, function(idx) {
    nino34[idx + seq(0L, lead_months)]
  })
}

#' Legacy rescale: trajectories x obs_sd / sd(model series), AFTER matching;
#' The legacy builder's rescale_trajectories().
.lenssen_rescale_legacy <- function(trajectories, model_nino34, obs_sd) {
  model_sd <- sd(model_nino34, na.rm = TRUE)
  scale_factor <- obs_sd / model_sd
  list(
    trajectories_adj = lapply(trajectories, function(traj) traj * scale_factor),
    scale_factor = scale_factor
  )
}

# ---------------------------------------------------------------------------
# Library matcher for the counterfactual pool
#
# One matcher, two uses: the counterfactual analog pool is drawn from the SAME
# continuation library as the factual arm's analog continuation (per-model
# series rescaled to the observed ONI sd over the library's window, 1991-2020),
# with the same features (level v, 3-month tendency, ONI 6 months earlier),
# standardisation, model exclusion and faded join. The seam is the last observed
# month (launch - 1), so the counterfactual is "the analog continuation from the
# last observed state"; column 1 of the cache is still the launch month, so the
# sampler's cache contract is unchanged. Sets are classified on every DJF the
# reported months read (winters = all_read); a row whose other winters are not
# neutral is labelled "mixed" (it stays in "all").
# ---------------------------------------------------------------------------

#' Defaults of deploy.enso_counterfactual (the legacy builder). Read
#' in-function.
#' `min_sep_months`: minimum months between two analog starts within one model
#' (the continuation uses 24 for deploy.enso_continuation.min_sep; the
#' counterfactual 12, which allows adjacent same-season years).
#' `seam_at` (library mode only; the legacy builder pins lead 0 to the observed
#' launch ONI and ignores it): where the faded join meets the observed ONI.
#' "launch" (default): every analog path is joined
#' onto the observed ONI of the launch month itself (cache column 1 = observed)
#' and the offset fades from there; "launch_minus_1": joined at the month
#' before the launch (column 1 = path(1); the first library builder, whose
#' launch month departs from observed ONI by the analogs' own first step).
#' Selection matches the pre-launch (launch - 1) state in both.
ENSO_COUNTERFACTUAL_DEFAULTS <- list(
  match = "legacy",
  n_analogs = 50L,
  min_sep_months = 12L,
  blend_months = 3,
  winters = "primary",
  read_span = "direct",
  other_winters = "neutral",
  min_set_size = 30L,
  seam_at = "launch"
)

#' The values of enso_counterfactual.seam_at (ENSO_COUNTERFACTUAL_DEFAULTS).
ENSO_COUNTERFACTUAL_SEAMS <- c("launch", "launch_minus_1")

#' enso_counterfactual_cfg
#'
#' deploy.enso_counterfactual (or NULL) merged over ENSO_COUNTERFACTUAL_DEFAULTS
#' and validated. Optional keys are read with [[ ]].
enso_counterfactual_cfg <- function(x = NULL) {
  out <- ENSO_COUNTERFACTUAL_DEFAULTS
  for (k in names(out)) {
    if (!is.null(x[[k]])) out[[k]] <- x[[k]]
  }
  out$match <- as.character(out$match)
  if (!out$match %in% c("legacy", "library")) {
    stop(
      "enso_counterfactual.match must be 'legacy' or 'library'; got '",
      out$match,
      "'.",
      call. = FALSE
    )
  }
  if (!out$winters %in% c("primary", "all_read")) {
    stop(
      "enso_counterfactual.winters must be 'primary' or 'all_read'.",
      call. = FALSE
    )
  }
  if (!out$read_span %in% c("direct", "unit")) {
    stop(
      "enso_counterfactual.read_span must be 'direct' or 'unit'.",
      call. = FALSE
    )
  }
  if (!out$other_winters %in% c("neutral", "any")) {
    stop(
      "enso_counterfactual.other_winters must be 'neutral' or 'any'.",
      call. = FALSE
    )
  }
  out$seam_at <- as.character(out$seam_at)
  if (
    length(out$seam_at) != 1L || !out$seam_at %in% ENSO_COUNTERFACTUAL_SEAMS
  ) {
    stop(
      "enso_counterfactual.seam_at must be 'launch' or 'launch_minus_1'; got '",
      paste(out$seam_at, collapse = ","),
      "'.",
      call. = FALSE
    )
  }
  out$n_analogs <- as.integer(out$n_analogs)
  out$min_sep_months <- as.integer(out$min_sep_months)
  out$blend_months <- as.numeric(out$blend_months)
  out$min_set_size <- as.integer(out$min_set_size)
  out
}

#' enso_read_winters
#'
#' Peak-season years (the calendar year of the season's first month: 2026 for
#' DJF 2026/27) whose season overlaps the ONI months [launch, read_end] and
#' lies wholly at or after `launch` (inside a trajectory starting there).
#' Seasons that begin before the launch but overlap it are returned in the
#' attribute "partial" (observed in part, so not classifiable on the pool).
enso_read_winters <- function(launch, read_end, peak_months = c(12L, 1L, 2L)) {
  launch <- as.Date(format(as.Date(launch), "%Y-%m-01"))
  read_end <- as.Date(format(as.Date(read_end), "%Y-%m-01"))
  peak_months <- as.integer(peak_months)
  if (is.na(read_end) || read_end < launch) {
    return(structure(integer(0), partial = integer(0)))
  }
  y0 <- as.integer(format(launch, "%Y")) - 1L
  y1 <- as.integer(format(read_end, "%Y"))
  full <- integer(0)
  partial <- integer(0)
  for (y in y0:y1) {
    first <- as.Date(sprintf("%d-%02d-01", y, peak_months[1L]))
    ms <- seq(first, by = "1 month", length.out = length(peak_months))
    if (!any(ms >= launch & ms <= read_end)) {
      next
    }
    if (first >= launch) full <- c(full, y) else partial <- c(partial, y)
  }
  structure(full, partial = partial)
}

#' enso_read_end
#'
#' The last ONI month the reported months read: max over the windows of
#' (window end - minimum deployed lag), since a month t with lag L reads ONI at
#' t - L. `windows`: data.frame(group, unit, end); `lags`: named list group ->
#' integer lags. span "direct": each group's own window end and lags (the ENSO
#' its reported months read); "unit": each unit's end (max over its groups)
#' with the minimum lag over the unit's groups (adds the months read only by
#' members censored at reporting, which still feed the unit's contagion).
enso_read_end <- function(windows, lags, span = c("direct", "unit")) {
  span <- match.arg(span)
  windows$end <- as.Date(windows$end)
  minlag <- vapply(
    as.character(windows$group),
    function(g) min(as.integer(unlist(lags[[g]]))),
    integer(1)
  )
  if (identical(span, "direct")) {
    ends <- windows$end
    ml <- minlag
  } else {
    u <- as.character(windows$unit)
    ends <- as.Date(
      tapply(as.numeric(windows$end), u, max),
      origin = "1970-01-01"
    )
    ml <- tapply(minlag, u, min)
  }
  ends <- as.Date(ends)
  max(.enso_shift_months(ends, -as.integer(ml)))
}

#' First-of-month Dates `d` shifted by `k` months (vectorised, k may be negative).
.enso_shift_months <- function(d, k) {
  mi <- as.integer(format(as.Date(d), "%Y")) *
    12L +
    as.integer(format(as.Date(d), "%m")) -
    1L +
    as.integer(k)
  as.Date(sprintf("%d-%02d-01", mi %/% 12L, mi %% 12L + 1L))
}

#' build_lenssen_counterfactual_pool
#'
#' The library matcher's counterfactual pool: up to `n_per_model` analogs per
#' model nearest to the seam `state` (enso_analog_neighbours(k = Inf,
#' max_per_model = n_per_model, horizon = n_months), within-model separation
#' `min_sep` months, `exclude_models` dropped), joined to the seam value by the
#' faded offset (enso_analog_paths()). Deterministic (no RNG).
#' `pin` (seam_at = "launch"): the observed ONI of the launch month. NULL (the
#' default here; seam_at = "launch_minus_1") joins each analog at the seam
#' month (column j = path(j) from analog(i)); a value joins it at the launch
#' month instead: column 1 = `pin`, column j = analog(i + j) + (pin -
#' analog(i + 1)) * rho^(j - 1) (.enso_pinned_paths()). The selection (the analog
#' starts i, matched on the launch - 1 state) is the same in both.
#' @return list(traj [n x n_months], model, init_index, distance)
build_lenssen_counterfactual_pool <- function(
  lib,
  state,
  n_months,
  n_per_model = 50L,
  min_sep = 12L,
  exclude_models = character(0),
  blend_months = 3,
  pin = NULL
) {
  nn <- enso_analog_neighbours(
    lib,
    state,
    k = Inf,
    min_sep = as.integer(min_sep),
    max_per_model = as.numeric(n_per_model),
    exclude_models = as.character(unlist(exclude_models)),
    horizon = as.integer(n_months)
  )
  cand <- lib$by_month[[as.character(state$mon)]]
  key <- paste(cand$model, cand$i)
  m <- match(paste(nn$model, nn$i), key)
  dist <- .enso_analog_distance(lib, cand[m, , drop = FALSE], state)
  # Model blocks in LENSSEN_MODELS / library order, nearest first within a model
  # (the legacy cache's row layout).
  mo <- match(nn$model, names(lib$series))
  ord <- order(mo, dist)
  nn <- nn[ord, , drop = FALSE]
  dist <- dist[ord]
  traj <- if (is.null(pin)) {
    enso_analog_paths(lib, nn, state$v, as.integer(n_months), blend_months)
  } else {
    .enso_pinned_paths(
      lib,
      nn,
      as.numeric(pin),
      as.integer(n_months),
      blend_months
    )
  }
  list(
    traj = traj,
    model = as.character(nn$model),
    init_index = as.integer(nn$i),
    distance = as.numeric(dist)
  )
}

#' Faded-join paths pinned at the launch month (seam_at = "launch"): column 1 =
#' `pin`, columns 2..n_months = enso_analog_paths() from analog(i + 1) with seam
#' value `pin` (path(0) = pin, so the offset pin - analog(i + 1) fades over
#' blend_months from the launch). Reads analog months i + 1 .. i + n_months,
#' the same span as the unpinned paths.
.enso_pinned_paths <- function(lib, nn, pin, n_months, blend_months = 3) {
  if (length(pin) != 1L || !is.finite(pin)) {
    stop(
      "counterfactual pool: the launch-month pin must be one finite ONI value.",
      call. = FALSE
    )
  }
  out <- matrix(pin, nrow = nrow(nn), ncol = n_months)
  if (n_months > 1L) {
    nn1 <- nn
    nn1$i <- nn1$i + 1L
    out[, -1L] <- enso_analog_paths(lib, nn1, pin, n_months - 1L, blend_months)
  }
  out
}

#' enso_launch_oni
#'
#' The observed ONI of the launch month that seam_at = "launch" pins cache
#' column 1 to (NULL under "launch_minus_1", which needs none). Stops when the
#' launch month's ONI is not published.
enso_launch_oni <- function(oni_lkp, launch, seam_at = "launch") {
  if (identical(seam_at, "launch_minus_1")) {
    return(NULL)
  }
  d <- format(as.Date(format(as.Date(launch), "%Y-%m-01")))
  v <- suppressWarnings(as.numeric(oni_lkp[d]))
  if (!length(v) || !is.finite(v)) {
    stop(
      "enso_counterfactual.seam_at 'launch': no published ONI for the launch month ",
      d,
      "; use seam_at 'launch_minus_1' or a launch month with published ONI.",
      call. = FALSE
    )
  }
  v
}

#' classify_enso_analog_sets_multi
#'
#' classify_enso_analog_sets() over several DJF winters. Labels come from the
#' primary winter(s) (several primaries must agree, else "mixed"); under
#' other_winters = "neutral" a row whose other winters are not all neutral is
#' "mixed". Returns classify_enso_analog_sets()'s list plus stat_matrix /
#' winter_labels (rows x winters) and `winters`; with no other winter and one
#' primary it IS classify_enso_analog_sets()'s result.
classify_enso_analog_sets_multi <- function(
  traj_mat,
  horizon_start,
  peak_months = c(12L, 1L, 2L),
  primary_years = NULL,
  other_years = integer(0),
  thresholds = ENSO_PHASE_THRESHOLDS,
  other_winters = c("neutral", "any")
) {
  other_winters <- match.arg(other_winters)
  prim <- if (length(primary_years)) as.integer(primary_years) else NULL
  base <- classify_enso_analog_sets(
    traj_mat,
    horizon_start,
    peak_months,
    peak_season_year = prim[1L],
    thresholds = thresholds
  )
  prim <- prim %||% base$definition$peak_season_year
  other <- setdiff(as.integer(other_years), prim)
  if (length(prim) == 1L && !length(other)) {
    return(base)
  }
  yrs <- c(prim, other)
  cls <- lapply(yrs, function(y) {
    if (y == prim[1L]) {
      base
    } else {
      classify_enso_analog_sets(
        traj_mat,
        horizon_start,
        peak_months,
        peak_season_year = y,
        thresholds = thresholds
      )
    }
  })
  lab_mat <- do.call(cbind, lapply(cls, `[[`, "labels"))
  stat_mat <- do.call(cbind, lapply(cls, `[[`, "stat"))
  colnames(lab_mat) <- colnames(stat_mat) <- as.character(yrs)
  labels <- lab_mat[, 1L]
  if (length(prim) > 1L) {
    agree <- apply(
      lab_mat[, as.character(prim), drop = FALSE],
      1L,
      function(r) length(unique(r)) == 1L
    )
    labels[!agree & !is.na(labels)] <- "mixed"
  }
  if (length(other) && identical(other_winters, "neutral")) {
    bad <- apply(lab_mat[, as.character(other), drop = FALSE], 1L, function(r) {
      any(is.na(r) | r != "neutral")
    })
    labels[bad & !is.na(labels)] <- "mixed"
  }
  sets <- c(setdiff(ENSO_ANALOG_SETS, "all"), "mixed")
  counts <- vapply(sets, function(s) sum(labels == s, na.rm = TRUE), integer(1))
  out <- base
  out$labels <- labels
  out$counts <- counts
  out$stat_matrix <- stat_mat
  out$winter_labels <- lab_mat
  out$winters <- list(
    primary = as.integer(prim),
    other = as.integer(other),
    other_winters = other_winters
  )
  out$definition$counts <- as.list(counts)
  out$definition$winters <- list(
    primary = as.list(as.integer(prim)),
    other = as.list(as.integer(other)),
    other_winters = other_winters,
    window_by_winter = lapply(stats::setNames(cls, yrs), function(cl) {
      list(start = cl$definition$window_start, end = cl$definition$window_end)
    })
  )
  out
}

#' lenssen_n_lead_months
#'
#' Months the counterfactual sampler draws ENSO for: every panel month from the
#' launch to the panel end (simulations.R derives sim_months from the panel),
#' plus `max_lead` lead months; the cache needs n_months = lead_months + 1 at
#' least that many. (Moved from make_event_attribution_configs.R.)
lenssen_n_lead_months <- function(launch, panel_end, max_lead) {
  m <- (as.integer(format(as.Date(panel_end), "%Y")) -
    as.integer(format(as.Date(launch), "%Y"))) *
    12L +
    as.integer(format(as.Date(panel_end), "%m")) -
    as.integer(format(as.Date(launch), "%m")) +
    1L
  as.integer(m + as.integer(max_lead))
}

#' lenssen_build_command
#'
#' The 03_prepare_lenssen_analogs.R command building one cache. The legacy spec
#' gives the event generator's former command string exactly; library mode adds
#' the matcher flags.
lenssen_build_command <- function(
  launch,
  lead_months,
  peak_season_year,
  output_dir,
  cf = NULL,
  read_end = NULL,
  library = NULL,
  exclude_models = character(0)
) {
  cmd <- sprintf(
    paste(
      "Rscript analysis/03_prepare_lenssen_analogs.R --horizon_start %s --lead_months %d",
      "--set_peak_season_year %d --output_dir %s --oni_cache_dir generated_data/"
    ),
    format(as.Date(launch)),
    as.integer(lead_months),
    as.integer(peak_season_year),
    output_dir
  )
  cf <- enso_counterfactual_cfg(cf)
  if (identical(cf$match, "legacy")) {
    return(cmd)
  }
  paste0(
    cmd,
    sprintf(
      paste(
        " --match library --library %s --n_analogs %d --min_sep %d --blend_months %s",
        "--exclude_models '%s' --set_winters %s --set_other_winters %s --min_set_size %d",
        "--seam_at %s%s"
      ),
      library %||% "generated_data/enso_continuation_library.qs",
      cf$n_analogs,
      cf$min_sep_months,
      format(cf$blend_months),
      paste(as.character(unlist(exclude_models)), collapse = ","),
      cf$winters,
      cf$other_winters,
      cf$min_set_size,
      cf$seam_at,
      if (!is.null(read_end)) {
        paste0(" --set_read_end ", format(as.Date(read_end)))
      } else {
        ""
      }
    )
  )
}

#' lenssen_cache_spec
#'
#' The identity of a library-mode counterfactual cache: everything it is built
#' from. lenssen_cache_status() compares it field by field; its content_digest()
#' is stamped in the cache meta (`spec_digest`).
#' seam_at = "launch" appends `seam_at` and `launch_oni` (the pinned observed
#' ONI, enso_launch_oni(); required) so the pin is part of the identity and the
#' digest; "launch_minus_1" gives the spec of the first library builder
#' unchanged (its caches keep their digests).
lenssen_cache_spec <- function(
  launch,
  n_months,
  cf,
  lib,
  seam_state,
  peak_months = c(12L, 1L, 2L),
  thresholds = ENSO_PHASE_THRESHOLDS,
  primary_years = NULL,
  other_years = integer(0),
  read_end = NULL,
  exclude_models = character(0),
  launch_oni = NULL
) {
  cf <- enso_counterfactual_cfg(cf)
  launch <- as.Date(launch)
  spec <- .lenssen_cache_spec_base(
    launch,
    n_months,
    cf,
    lib,
    seam_state,
    peak_months,
    thresholds,
    primary_years,
    other_years,
    read_end,
    exclude_models
  )
  if (identical(cf$seam_at, "launch")) {
    if (
      is.null(launch_oni) ||
        length(launch_oni) != 1L ||
        !is.finite(as.numeric(launch_oni))
    ) {
      stop(
        "lenssen_cache_spec: seam_at 'launch' needs the launch month's observed ONI (enso_launch_oni()).",
        call. = FALSE
      )
    }
    spec$seam_at <- "launch"
    spec$launch_oni <- round(as.numeric(launch_oni), 6)
  }
  spec
}
.lenssen_cache_spec_base <- function(
  launch,
  n_months,
  cf,
  lib,
  seam_state,
  peak_months,
  thresholds,
  primary_years,
  other_years,
  read_end,
  exclude_models
) {
  list(
    match = "library",
    horizon_start = format(launch),
    seam_month = format(seq(launch, by = "-1 month", length.out = 2L)[2L]),
    seam_state = lapply(seam_state[c("v", "tend", "phase")], function(v) {
      round(as.numeric(v), 6)
    }),
    n_months = as.integer(n_months),
    n_analogs = cf$n_analogs,
    min_sep_months = cf$min_sep_months,
    blend_months = cf$blend_months,
    exclude_models = as.list(sort(as.character(unlist(exclude_models)))),
    library = lib$path %||% "",
    library_digest = lib$meta$digest %||% enso_continuation_library_digest(lib),
    sd_window = lib$meta$rescale$window %||% "",
    obs_sd = lib$meta$rescale$obs_sd %||% NA_real_,
    peak_months = as.integer(peak_months),
    thresholds = as.list(thresholds),
    winters = cf$winters,
    read_span = cf$read_span,
    other_winters = cf$other_winters,
    read_end = if (is.null(read_end)) "" else format(as.Date(read_end)),
    primary_years = as.list(as.integer(primary_years %||% integer(0))),
    other_years = as.list(as.integer(other_years))
  )
}

#' Content digest of a counterfactual cache: trajectories, labels, models,
#' init indices and the spec (base.R::content_digest()).
lenssen_cache_digest <- function(cache) {
  content_digest(list(
    traj = unname(cache$oni_trajectories),
    labels = cache$set_labels,
    model = cache$model,
    init_index = cache$init_index,
    spec = cache$meta$spec %||% list()
  ))
}

#' build_lenssen_cache
#'
#' Build a library-mode counterfactual cache from `spec` (lenssen_cache_spec()),
#' the observed ONI (fetch_enso_index() tibble) and the continuation library;
#' write it (RDS + `lenssen_nino34_analogs_meta.yml`) when `out_path` is given.
#' Deterministic. The cache keeps the legacy contract (oni_trajectories,
#' n_months, n_analogs, set_labels, set_stat, model, init_index, meta) plus
#' `distance`, `stat_matrix` and `meta$spec` / `meta$digest`.
build_lenssen_cache <- function(spec, oni, lib, out_path = NULL) {
  oni_lkp <- stats::setNames(
    as.numeric(oni$enso_index),
    format(as.Date(oni$date))
  )
  launch <- as.Date(spec$horizon_start)
  state <- enso_seam_state(oni_lkp, spec$seam_month)
  pool <- build_lenssen_counterfactual_pool(
    lib,
    state,
    n_months = spec$n_months,
    n_per_model = spec$n_analogs,
    min_sep = spec$min_sep_months,
    exclude_models = unlist(spec$exclude_models),
    blend_months = spec$blend_months,
    pin = if (identical(spec$seam_at, "launch")) spec$launch_oni
  )
  thr <- unlist(spec$thresholds)
  cl <- classify_enso_analog_sets_multi(
    pool$traj,
    launch,
    spec$peak_months,
    primary_years = unlist(spec$primary_years),
    other_years = if (identical(spec$winters, "all_read")) {
      unlist(spec$other_years)
    } else {
      integer(0)
    },
    thresholds = thr,
    other_winters = spec$other_winters
  )
  obs_launch <- suppressWarnings(as.numeric(oni_lkp[format(launch)]))
  sets <- sort(unique(stats::na.omit(cl$labels)))
  stat_mat <- cl$stat_matrix %||%
    matrix(
      cl$stat,
      ncol = 1L,
      dimnames = list(NULL, as.character(cl$definition$peak_season_year))
    )
  set_summary <- lapply(stats::setNames(nm = c("all", sets)), function(s) {
    rows <- if (identical(s, "all")) {
      seq_along(cl$labels)
    } else {
      which(cl$labels == s)
    }
    list(
      n = length(rows),
      model_share = as.list(round(prop.table(table(pool$model[rows])), 4)),
      seam_error = if (length(rows) && is.finite(obs_launch)) {
        round(mean(pool$traj[rows, 1L]) - obs_launch, 4)
      } else {
        NA_real_
      },
      winter_means = if (length(rows)) {
        as.list(round(colMeans(stat_mat[rows, , drop = FALSE]), 4))
      } else {
        list()
      }
    )
  })
  meta <- list(
    source = "Lenssen et al. 2024 GRL piControl via the ENSO continuation library (review item 8)",
    match = "library",
    state_matching = "library_v_tend_phase",
    models = as.list(unique(pool$model)),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    horizon_start = format(launch),
    seam_month = spec$seam_month,
    target_oni = round(state$v, 4),
    target_month = as.integer(state$mon),
    n_analogs_per_model = spec$n_analogs,
    min_sep_months = spec$min_sep_months,
    blend_months = spec$blend_months,
    lead_months = as.integer(spec$n_months) - 1L,
    obs_sd_oni = spec$obs_sd,
    sd_window = spec$sd_window,
    n_models = length(unique(pool$model)),
    n_total_analogs = nrow(pool$traj),
    exclude_models = spec$exclude_models,
    set_definition = c(
      cl$definition,
      list(per_set_model_share = lapply(set_summary[sets], `[[`, "model_share"))
    ),
    sets = set_summary,
    spec = spec,
    spec_digest = content_digest(spec)
  )
  # The launch-month pin (seam_at = "launch"); absent under launch_minus_1, whose
  # meta stays as before.
  if (identical(spec$seam_at, "launch")) {
    meta <- c(meta, list(seam_at = "launch", launch_oni = spec$launch_oni))
  }
  cache <- list(
    oni_trajectories = pool$traj,
    n_months = as.integer(spec$n_months),
    n_analogs = nrow(pool$traj),
    set_labels = cl$labels,
    set_stat = cl$stat,
    model = pool$model,
    init_index = pool$init_index,
    distance = pool$distance,
    stat_matrix = cl$stat_matrix,
    meta = meta
  )
  cache$meta$digest <- lenssen_cache_digest(cache)
  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(cache, out_path)
    yaml::write_yaml(
      cache$meta,
      file.path(dirname(out_path), "lenssen_nino34_analogs_meta.yml")
    )
  }
  cache
}

#' lenssen_cache_status
#'
#' Whether the cache at `path` is the one `spec` describes: built in library
#' mode with the same spec (field by field), its stored content digest intact,
#' and the `sets` at least `min_set_size` strong. list(ok, reason, path,
#' counts, small, digest).
lenssen_cache_status <- function(path, spec, sets = "all", min_set_size = 30L) {
  res <- function(
    ok,
    reason,
    counts = NULL,
    small = character(0),
    digest = NA_character_
  ) {
    list(
      ok = ok,
      reason = reason,
      path = path,
      counts = counts,
      small = small,
      digest = digest
    )
  }
  if (!file.exists(path)) {
    return(res(FALSE, "cache absent"))
  }
  cache <- readRDS(path)
  have <- cache$meta$spec
  if (is.null(have)) {
    return(res(FALSE, "a legacy cache (no library spec)"))
  }
  # Over the union of the fields: a field only one side has (seam_at /
  # launch_oni of a pinned cache vs a launch_minus_1 spec) is a difference.
  keys <- union(names(spec), names(have))
  diff <- keys[
    !vapply(keys, function(k) identical(have[[k]], spec[[k]]), logical(1))
  ]
  if (length(diff)) {
    return(res(FALSE, paste("spec differs:", paste(diff, collapse = ", "))))
  }
  if (!identical(cache$meta$digest, lenssen_cache_digest(cache))) {
    return(res(FALSE, "content digest does not match the stored one"))
  }
  counts <- vapply(
    stats::setNames(nm = sets),
    function(s) {
      if (identical(s, "all")) {
        nrow(cache$oni_trajectories)
      } else {
        sum(cache$set_labels == s, na.rm = TRUE)
      }
    },
    integer(1)
  )
  small <- names(counts)[counts < as.integer(min_set_size)]
  res(
    TRUE,
    if (length(small)) {
      sprintf(
        "sets under %d: %s",
        as.integer(min_set_size),
        paste(sprintf("%s=%d", small, counts[small]), collapse = ", ")
      )
    } else {
      "ok"
    },
    counts,
    small,
    cache$meta$digest
  )
}

#' lenssen_counterfactual_args
#'
#' The counterfactual block's Lenssen arguments, resolved once for every
#' consumer: cache path, the
#' analog anchor (`lenssen_horizon_start`, NULL / "" = the arm's own
#' `horizon_start`; an invalid date stops), mean offset, variance scaling and
#' analog set.
#' @param cf simulator.enso.counterfactual (list or NULL)
#' @param horizon_start the arm's horizon.start (Date)
#' @return list(cache_path, horizon_start, mean_offset, log_var_scaling, analog_set)
lenssen_counterfactual_args <- function(cf, horizon_start) {
  cf <- cf %||% list()
  lhs <- cf[["lenssen_horizon_start"]]
  # Guard NULL / "" explicitly: as.Date(NULL) is a length-0 Date and
  # `if (is.na(<length 0>))` would throw (the former handoff-jump bug).
  enso_horizon <- if (is.null(lhs) || !nzchar(as.character(lhs))) {
    as.Date(horizon_start)
  } else {
    as.Date(lhs)
  }
  if (is.na(enso_horizon)) {
    stop(sprintf(
      "enso.counterfactual.lenssen_horizon_start = '%s' is not a valid date.",
      lhs
    ))
  }
  list(
    cache_path = cf[["cache_path"]] %||%
      "generated_data/lenssen_nino34_analogs.rds",
    horizon_start = enso_horizon,
    mean_offset = as.numeric(cf[["mean_offset"]] %||% 0.0),
    log_var_scaling = as.numeric(cf[["log_var_scaling"]] %||% 0.0),
    analog_set = as.character(cf[["analog_set"]] %||% "all")
  )
}

#' lenssen_cache_provenance
#'
#' The counterfactual cache's identity for manifests and partition_meta.yml:
#' path, content digest, spec digest, match, set counts and winters. NULL for a
#' legacy cache (no spec), so legacy runs record nothing new.
lenssen_cache_provenance <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  cache <- readRDS(path)
  m <- cache$meta
  if (is.null(m$spec)) {
    return(NULL)
  }
  out <- list(
    cache_path = path,
    digest = m$digest,
    spec_digest = m$spec_digest,
    match = m$match,
    library_digest = m$spec$library_digest,
    seam_month = m$seam_month,
    winters = m$set_definition$winters %||%
      list(primary = list(m$set_definition$peak_season_year)),
    counts = m$set_definition$counts
  )
  # The launch-month pin (seam_at = "launch"); no element under launch_minus_1,
  # so that provenance is unchanged.
  if (!is.null(m$spec$seam_at)) {
    out$seam_at <- m$spec$seam_at
  }
  out
}

#' The shared, unversioned legacy counterfactual cache (the skeleton default).
#' Library mode refuses it, so a library build can never replace the input of
#' a legacy run.
LENSSEN_SHARED_CACHE <- "generated_data/lenssen_nino34_analogs.rds"

#' ensure_lenssen_counterfactual_cache
#'
#' Library mode (deploy.enso_counterfactual.match = "library"): build the
#' counterfactual cache at `path` in-process, or verify the one there, for a
#' run launched at `launch` whose reported months read ONI up to `read_end`.
#' Deterministic and reading only generated_data/ (library + ONI cache), so the
#' generators run it on the cluster. Legacy mode returns NULL and does nothing.
#'
#' @param path        cache RDS (a per-run path; the shared legacy path is refused)
#' @param launch      the counterfactual's analog anchor (cache column 1)
#' @param n_months    trajectory length the sampler needs (lead months + 1)
#' @param cf          deploy.enso_counterfactual (enso_counterfactual_cfg())
#' @param cont        deploy.enso_continuation (library, exclude_models)
#' @param read_end    last ONI month the reported months read (enso_read_end())
#' @param primary_year peak-season year the sets are labelled on (NULL = the
#'   first winter read)
#' @param sets        analog sets the run uses; `primary_set` must reach
#'   min_set_size (else stop), the others are reported
#' @param oni_dir     directory of the ONI cache (fetch_enso_index())
#' @return NULL (legacy) or list(path, status, spec, provenance, winters)
ensure_lenssen_counterfactual_cache <- function(
  path,
  launch,
  n_months,
  cf,
  cont = NULL,
  read_end,
  primary_year = NULL,
  sets = "all",
  primary_set = "all",
  oni_dir = "generated_data/",
  label = "[counterfactual]"
) {
  cf <- enso_counterfactual_cfg(cf)
  if (identical(cf$match, "legacy")) {
    return(NULL)
  }
  if (
    normalizePath(path, mustWork = FALSE) ==
      normalizePath(LENSSEN_SHARED_CACHE, mustWork = FALSE)
  ) {
    stop(
      label,
      " enso_counterfactual.match 'library' refuses the shared legacy cache path '",
      path,
      "': set simulator.enso.counterfactual.cache_path to a per-run path ",
      "(e.g. generated_data/lenssen_prod_run/lenssen_nino34_analogs.rds).",
      call. = FALSE
    )
  }
  lib <- load_enso_continuation_library(
    (cont %||% list())[["library"]] %||%
      "generated_data/enso_continuation_library.qs"
  )
  excl <- as.character(unlist(
    (cont %||% list())[["exclude_models"]] %||% character(0)
  ))
  oni <- fetch_enso_index(source = "oni", cache_dir = oni_dir, refresh = FALSE)
  oni_lkp <- stats::setNames(
    as.numeric(oni$enso_index),
    format(as.Date(oni$date))
  )
  launch <- as.Date(launch)
  seam <- format(seq(launch, by = "-1 month", length.out = 2L)[2L])
  state <- enso_seam_state(oni_lkp, seam)
  rw <- enso_read_winters(launch, read_end)
  # primary_year may hold several winters (a multi-winter event of
  # enso_events_v2.yml): a set label then needs every primary winter to agree.
  prim <- if (!is.null(primary_year)) as.integer(primary_year) else rw[1L]
  if (!length(prim) || anyNA(prim)) {
    # No DJF is read inside the reported window (e.g. firm windows that end at
    # the last CFSv2 month): label the sets on the first complete DJF after the
    # launch, as the legacy builder does (classify_enso_analog_sets() "auto").
    prim <- as.integer(format(launch, "%Y"))
    message(sprintf(
      "%s no DJF is read in [%s, %s]: sets labelled on DJF %d/%02d (the first after the launch).",
      label,
      format(launch),
      format(as.Date(read_end)),
      prim,
      (prim + 1L) %% 100L
    ))
  }
  other <- if (identical(cf$winters, "all_read")) {
    setdiff(as.integer(rw), prim)
  } else {
    integer(0)
  }
  spec <- lenssen_cache_spec(
    launch,
    n_months,
    cf,
    lib,
    state,
    primary_years = prim,
    other_years = other,
    read_end = read_end,
    exclude_models = excl,
    launch_oni = enso_launch_oni(oni_lkp, launch, cf$seam_at)
  )
  st <- lenssen_cache_status(
    path,
    spec,
    sets = sets,
    min_set_size = cf$min_set_size
  )
  if (!isTRUE(st$ok)) {
    message(sprintf(
      "%s building the library-mode Lenssen cache %s (%s).",
      label,
      path,
      st$reason
    ))
    build_lenssen_cache(spec, oni, lib, out_path = path)
    st <- lenssen_cache_status(
      path,
      spec,
      sets = sets,
      min_set_size = cf$min_set_size
    )
    if (!isTRUE(st$ok)) {
      stop(
        label,
        " cache ",
        path,
        " still not valid after the build: ",
        st$reason,
        call. = FALSE
      )
    }
  } else {
    message(sprintf(
      "%s library-mode Lenssen cache %s verified (digest %s).",
      label,
      path,
      st$digest
    ))
  }
  if (primary_set %in% st$small) {
    stop(
      sprintf(
        paste0(
          "%s analog set '%s' has %d trajectories, under enso_counterfactual.min_set_size %d ",
          "(winters %s). Raise n_analogs or relax the winter rule."
        ),
        label,
        primary_set,
        st$counts[[primary_set]],
        cf$min_set_size,
        paste(c(prim, other), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(st$small)) {
    warning(sprintf("%s %s", label, st$reason), call. = FALSE)
  }
  list(
    path = path,
    status = st,
    spec = spec,
    provenance = lenssen_cache_provenance(path),
    winters = list(primary = prim, other = other, partial = attr(rw, "partial"))
  )
}
