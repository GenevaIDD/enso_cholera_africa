# clim_compile.R — read raw climate files per GADM unit and compile monthly
#   stats and extremes (CSV or parquet source).
#
# Provides:
#   get_clim_files_gid(), read_clim_csv()       raw climate files of a gid
#   get_precip_vars(), get_temp_vars()          variable names
#   compute_baseline_stats()                    baseline climatology
#   compile_all_stats(), compile_all_extremes() monthly stats / extremes
#   compile_all_stats_parquet(), compile_all_extremes_parquet()
#
# Details:
#   Used by 02_preprocess_climate.R (stats into climate_stats.sqlite) and the
#   03 / 04 data builders.

# Climate data wrangling --------------------------------------------------

#' get_clim_files_gid
#' get all climate data files corresponding to a given GADM id
#'
#' @param gid
#' @param clim_dir
#'
get_clim_files_gid <- function(gid, clim_dir = "data/climate_data/") {
  country <- stringr::str_extract(gid, "[A-Z]{3}")
  gid_num <- stringr::str_remove(gid, stringr::str_c(country, "_")) %>%
    stringr::str_remove(country)

  this_dir <- dir(clim_dir, pattern = country, full.names = T, recursive = F)

  dir(this_dir, pattern = gid_num, full.names = T, recursive = T) %>%
    stringr::str_subset(stringr::str_c("v410_", gid_num, "_[A-Z]")) %>%
    stringr::str_subset("txt", negate = T) %>%
    stringr::str_subset("1980|1990|2000|2010|2020")
}

#' parse_file_info
#' Parse the information in a climate file name
#' @param str
#'
parse_file_info <- function(str) {
  country <- stringr::str_extract(str, "[A-Z]{3}")
  for (l in 2:1) {
    unit <- stringr::str_extract(
      str,
      stringr::str_c(
        "(?<=410_)",
        stringr::str_c(rep("[0-9]+", l + 1), collapse = "_")
      )
    )
    admin_level <- l
    if (!is.na(unit)) {
      break()
    }
  }

  if (is.na(unit)) {
    unit <- "country"
    admin_level <- 0
  }

  period <- stringr::str_extract(str, "[0-9]{4}_[0-9]{4}")
  variable <- stringr::str_extract(
    str,
    glue::glue("(?<={unit}_)(.)*(?=_{period})")
  )
  tibble(
    country = country,
    admin_level = admin_level,
    unit_num = unit,
    period = period,
    variable = variable,
    file = str,
    str = stringr::str_c(country, "_ADM", admin_level, "_", unit),
  )
}


# Climate data stats ------------------------------------------------------

read_clim_csv <- function(x) {
  # Namespaced: furrr's multisession workers do not attach readr.
  res <- readr::read_csv(x, na = c("", "NA", "-9999"), show_col_types = FALSE)

  if (nrow(res) == 0) {
    warning("File ", x, " is empty.")
  }

  res
}

#' compute_stats_precip
#' Compute the statistics of monthly precipitation
#'
#' @param x csv file to process
#'
compute_stats_precip <- function(x) {
  res <- read_clim_csv(x)

  if (nrow(res) == 0) {
    return(NULL)
  }

  res %>%
    select(Date, contains("2010")) %>%
    janitor::clean_names() %>%
    mutate(
      month = format(date, "%Y-%m") %>% stringr::str_c("-01") %>% as.Date(),
      pl1 = lag(chirps_pop_2010, 1),
      pl2 = lag(chirps_pop_2010, 2),
      pl3 = lag(chirps_pop_2010, 3),
      pl4 = lag(chirps_pop_2010, 4),
      pl5 = lag(chirps_pop_2010, 5),
      pl6 = lag(chirps_pop_2010, 6),
      pcum7 = chirps_pop_2010 + pl1 + pl2 + pl3 + pl4 + pl5 + pl6 + pl6
    ) %>%
    group_by(month) %>%
    summarise(
      precip = sum(chirps_pop_2010),
      max_precip = max(chirps_pop_2010),
      mean_cum7 = mean(pcum7, na.rm = T),
      max_cum7 = max(pcum7, na.rm = T)
    )
}

#' compute_stats_temp
#' Compute statistics of monthly temperature
#'
#' @param x csv file to process
#'
compute_stats_temp <- function(x) {
  stat <- stringr::str_extract(x, "daymin|daymax")

  if (is.na(stat)) {
    stat <- "daymean"
  }

  dict <- c("value")
  names(dict) <- stringr::str_c("temp_", stat)

  res <- read_clim_csv(x)

  if (nrow(res) == 0) {
    return(NULL)
  }

  res %>%
    select(Date, contains("2010")) %>%
    janitor::clean_names() %>%
    mutate(
      month = format(date, "%Y-%m") %>% stringr::str_c("-01") %>% as.Date()
    ) %>%
    group_by(month) %>%
    summarise(value = mean(era5_land_pop_2010)) %>%
    rename(dict)
}


#' get_precip_vars
#' Returns columns names of precipitation statistics
#'
get_precip_vars <- function() {
  c("precip", "max_precip", "mean_cum7", "max_cum7")
}

#' get_temp_vars
#' Returns columns names of precipitation statistics
#'
get_temp_vars <- function() {
  c("temp_daymean", "temp_daymax", "temp_daymin")
}


#' compute_baseline_stats
#' Compute the monthly statistics for the period 1980-2009
#'
#' @param dat the dataframe of statistics to preocess
#' @param vars the variables to compute statistics for
#'
compute_baseline_stats <- function(
  gid,
  stats = NULL,
  var = "precip",
  redo = FALSE,
  baseline
) {
  gid <- flatten_gid(gid)
  out_file <- make_clim_baseline_stats_file(
    var = var,
    gid = gid,
    baseline = baseline
  )

  if (!file.exists(out_file) | redo) {
    # Get climate stats if not provided
    if (is.null(stats)) {
      stats <- compile_stats(gid = gid, var = var, redo = redo)
    }

    cat("-- Computing baseline stats for", var, gid, baseline, "\n")

    # Get variables of interest
    clim_vars <- get_clim_vars(var)
    # Make the left and right time bounds
    time_bounds <- make_baseline_time_bounds(baseline)

    baseline_stats <- stats %>%
      filter(
        lubridate::year(month) >= time_bounds$TL &
          lubridate::year(month) <= time_bounds$TR
      ) %>%
      mutate(month_num = lubridate::month(month)) %>%
      tidyr::pivot_longer(any_of(clim_vars)) %>%
      group_by(month_num, name) %>%
      summarise(
        q10 = quantile(value, 0.1, na.rm = T),
        q90 = quantile(value, 0.9, na.rm = T),
        mean = mean(value, na.rm = T),
        sd = sd(value, na.rm = T)
      )

    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # a fresh checkout has no cache dirs
    saveRDS(baseline_stats, file = out_file)
  } else {
    cat("-- Reading pre-computed baseline stats for", var, gid, "\n")
    baseline_stats <- readRDS(out_file)
  }

  baseline_stats
}


make_baseline_time_bounds <- function(baseline) {
  time_bounds <- list(
    TL = as.numeric(stringr::str_extract(baseline, "[0-9]{4}")),
    TR = as.numeric(stringr::str_extract(baseline, "(?<=_)[0-9]{4}"))
  )

  if (any(is.na(time_bounds))) {
    stop(
      "Failed in extracting time bounds from ",
      baseline,
      ". Provide as YYYY_YYYY"
    )
  }

  time_bounds
}


#' get_clim_vars
#' Define variables of interest
#' @param var variable, one of precip or temp
#'
get_clim_vars <- function(var) {
  if (var == "temp") {
    clim_vars <- get_temp_vars()
  } else if (var == "precip") {
    clim_vars <- get_precip_vars()
  } else if (var == "all_parquet") {
    # The parquet source puts every variable in one file per admin unit, so it
    # has no temp/precip split to dispatch on. Branch added so the untouched
    # compute_baseline_stats() / compute_extremes() can be reused verbatim.
    clim_vars <- get_clim_vars_parquet()
  } else {
    stop("Variable ", var, " not found.")
  }
  clim_vars
}

#' get_stats_fun
#' Define function to compute statistics
#' @param var
#'
get_stats_fun <- function(var) {
  # Define variables of interest and stat functions
  if (var == "temp") {
    statsFun <- compute_stats_temp
  } else if (var == "precip") {
    statsFun <- compute_stats_precip
  } else {
    stop("Variable ", var, " not found.")
  }
  statsFun
}

#' get_match_pattern
#' Define pattern to look for in climate outputs
#' @param var
#'
get_match_pattern <- function(var) {
  # Define variables of interest and stat functions
  if (var == "temp") {
    match_pattern <- "ERA5_Land_2m_temperature"
  } else if (var == "precip") {
    match_pattern <- "CHIRPS_total_precipitation"
  } else {
    stop("Variable ", var, " not found.")
  }
  match_pattern
}


#' compile_stats
#' Computes and compiles the stats data for a given gid
#'
#' @param gid the gid for which to compile the data
#' @param var the variable to compile data for, one of `temp` or `precip`
#' @param redo whether to redo the compilation
#'
compile_stats <- function(gid, var = "temp", redo = TRUE) {
  gid <- flatten_gid(gid)
  out_file <- make_clim_stats_file(var = var, gid = gid)

  if (!file.exists(out_file) | redo) {
    cat("-- Computing stats for", var, gid, "\n")

    clim_files <- get_clim_files_gid(gid)

    if (length(clim_files) == 0) {
      warning("No files found for", gid)
      return(NULL)
    }

    # Define variables, function and pattern
    statsFun <- get_stats_fun(var)
    match_pattern <- get_match_pattern(var)

    # Get files corresponding to variable of interest
    var_files <- stringr::str_subset(clim_files, match_pattern)

    if (length(var_files) == 0) {
      warning("No files found for ", gid, " ", var)
      return(NULL)
    }

    # Parse the file information to group by unique variable
    var_files_info <- purrr::map_df(var_files, parse_file_info)

    # Compute statistics for each variable
    res_list <- var_files_info %>%
      group_by(variable) %>%
      group_map(function(x, y) {
        purrr::map_df(x$file, ~ statsFun(.))
      })

    has_res <- purrr::map_lgl(res_list, ~ nrow(.) > 0)

    if (!any(has_res)) {
      warning("All files empty for ", gid, " ", var)
      return(NULL)
    }

    res <- res_list[has_res] %>%
      purrr::reduce(full_join) %>%
      arrange(month)

    # Save
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # a fresh checkout has no cache dirs
    saveRDS(res, file = out_file)
  } else {
    cat("-- Reading pre-computed stats for", var, gid, "\n")
    res <- readRDS(out_file)
  }

  res %>%
    mutate(gid = gid)
}


#' compile_extremes
#' Computes and compiles the extremes data for a given gid
#'
#' @param gid the gid for which to compile the data
#' @param var the variable to compile data for, one of `temp` or `precip`
#' @param redo whether to redo the compilation
#'
compile_extremes <- function(
  gid,
  var = "temp",
  redo = TRUE,
  baseline = "1980-2009"
) {
  gid <- flatten_gid(gid)
  out_file <- make_clim_extremes_file(var = var, gid = gid, baseline = baseline)

  if (!file.exists(out_file) | redo) {
    cat("-- Computing extrems for", var, gid, baseline, "\n")

    stats <- compile_stats(gid = gid, var = var, redo = redo)

    if (is.null(stats)) {
      warning(
        "Extremes computation failed for ",
        gid,
        " ",
        var,
        " because no data was found."
      )
      return(NULL)
    }

    # Compute baseline statistics for 30yr period as defined by user
    baseline_stats <- compute_baseline_stats(
      gid = gid,
      stats = stats,
      var = var,
      redo = redo,
      baseline = baseline
    )

    # Define climate variables of interest
    clim_vars <- get_clim_vars(var)

    # Compute extremes for each variable
    res <- compute_extremes(
      dat = stats,
      baseline_stats = baseline_stats,
      vars = clim_vars
    ) %>%
      mutate(gid = gid)

    # Save
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # a fresh checkout has no cache dirs
    saveRDS(res, file = out_file)
  } else {
    cat("-- Reading pre-computed extrems for", var, gid, "\n")
    res <- readRDS(out_file)
  }

  res
}

#' compile_all_stats
#' Wrapper to compute all climate statistics for a given gid
#'
#' @param gid the gid for which to compile the data
#' @param var the variable to compile data for, one of `temp` or `precip`
#' @param redo whether to redo the compilation
#'
compile_all_stats <- function(gid, vars = c("temp", "precip"), redo = FALSE) {
  res <- purrr::map(vars, ~ compile_stats(gid = gid, var = ., redo = redo))

  if (all(purrr::map_lgl(res, is.null))) {
    return(NULL)
  } else {
    res <- res %>%
      purrr::reduce(full_join) %>%
      arrange(month)
  }

  cat(
    "---------------------------------
      -- Done stats computations for ",
    gid,
    "
      ---------------------------------\n"
  )

  res
}

#' compile_all_extremes
#' Wrapper to compute all climate extremes for a given gid
#'
#' @param gid the gid for which to compile the data
#' @param var the variable to compile data for, one of `temp` or `precip`
#' @param redo whether to redo the compilation
#'
compile_all_extremes <- function(
  gid,
  vars = c("temp", "precip"),
  redo = FALSE,
  baseline
) {
  res <- purrr::map_df(
    vars,
    ~ compile_extremes(gid = gid, var = ., redo = redo, baseline = baseline)
  )

  cat(
    "----------------------------------------------- \n",
    "-- Done extreme computations for ",
    gid,
    "\n",
    "----------------------------------------------- \n"
  )

  res
}


# Climate data stats — parquet source -------------------------------------
#
# Reads the parquet release (Jahn et al. 2026, doi 10.21203/rs.3.rs-9508034/v1):
# the same CHIRPS / ERA5-Land admin-level product as data/climate_data/, one
# file per admin unit, all variables, 1981-01..2024-12. These functions sit
# beside the CSV readers above, which stay unchanged (`--climate_source csv`
# is bit-identical). One file per gid also removes compile_stats()'s
# full_join, which multiplied rows when a gid substring-matched several files
# (malformed GADM Ghana ids).

#' make_clim_parquet_core
#' The region-independent tail of a unit's parquet filename, derived from its gid.
#'
#' The publisher's rule, verified against COM / GHA / LSO / STP / BEN / GNB / NER
#' / SLE: a filename is `<REGION>_<core>`, where `<core>` is built by splitting
#' the gid at its FIRST dot — head before `_v410`, dot-flattened tail after.
#'
#'   COM.1_1     -> COM_v410_1_1_...      (ADM1)
#'   STP.1.1_1   -> STP_v410_1_1_1_...    (ADM2)
#'   COM         -> COM_v410_...          (ADM0, no tail)
#'   GHA14_2     -> GHA14_2_v410_...      (GADM's malformed Ghana ADM1 id has no
#'                                         dot at all, so the whole id is head)
#'   GHA14.11_2  -> GHA14_v410_11_2_...   (malformed Ghana ADM2)
#'
#' The Ghana head-only form does not match a `v410_<num>_[A-Z]` pattern, so a
#' reader relying on one would substring-match those units onto their ADM2
#' siblings.
#'
#' `<REGION>` is deliberately NOT reconstructed: it is the publisher's own
#' grouping and is not derivable from the ISO3. Most African countries are
#' `Africa_`, but Niger is `MENA_` (and presumably the other Maghreb/Sahel
#' countries too). Hardcoding `Africa_` silently lost all 16 overlapping NER
#' gids in testing.
#'
#' @param gid a GADM id in the raw form carried by the GADM join / crosswalk
#'   (dotted, e.g. `AGO.4_1`, or Ghana-malformed, e.g. `GHA14.11_2`). NOT the
#'   flattened `AGO_4_1` form used as the database key — flattening is lossy
#'   here because it erases the dot the split depends on.
#'
make_clim_parquet_core <- function(gid) {
  parts <- stringr::str_split(gid, stringr::fixed("."), n = 2)[[1]]
  stem <- stringr::str_c(parts[1], "_v410")
  if (length(parts) > 1 && nzchar(parts[2])) {
    stem <- stringr::str_c(stem, "_", flatten_gid(parts[2]))
  }
  stringr::str_c(stem, "_CHIRPSv3_ERA5Land_1981_2024_observation.parquet.gzip")
}


#' get_clim_parquet_gid
#' Resolve a GADM id to its single parquet file.
#'
#' Matches on the full `_<core>` suffix (an exact string, not a substring
#' pattern), so the only freedom left is the region prefix. Substring matching is
#' what let `GHA5_2` pull in 23 sibling files; an exact suffix cannot, because
#' `..._4_1_1_CHIRPSv3_...` does not end with `..._4_1_CHIRPSv3_...`.
#'
#' @return the file path, or NULL if the country was never downloaded or the unit
#'   is not published. Multiple matches are fatal — that would mean one unit is
#'   published under two region prefixes, and picking one silently is exactly the
#'   failure mode this rewrite exists to remove.
#'
get_clim_parquet_gid <- function(gid, clim_dir = "data/climate_parquet/") {
  country <- stringr::str_extract(gid, "[A-Z]{3}")
  if (is.na(country)) {
    warning("Could not extract an ISO3 country code from gid '", gid, "'")
    return(NULL)
  }

  country_dir <- file.path(clim_dir, country)
  if (!dir.exists(country_dir)) {
    return(NULL)
  }

  suffix <- stringr::str_c("_", make_clim_parquet_core(gid))
  files <- list.files(country_dir)
  hits <- files[endsWith(files, suffix)]

  if (length(hits) == 0) {
    return(NULL)
  }
  if (length(hits) > 1) {
    stop(
      "get_clim_parquet_gid(): gid '",
      gid,
      "' matched ",
      length(hits),
      " parquet files in ",
      country_dir,
      ": ",
      paste(hits, collapse = ", "),
      "\n  A unit must resolve to exactly one file. Refusing to guess."
    )
  }

  file.path(country_dir, hits)
}


#' get_clim_vars_parquet
#' Monthly statistic columns produced by compute_stats_parquet().
#'
#' The first seven are get_precip_vars() + get_temp_vars(), which downstream
#' readers expect; the rest are extra columns the parquet release provides.
#'
get_clim_vars_parquet <- function() {
  c(
    get_precip_vars(), # precip, max_precip, mean_cum7, max_cum7
    get_temp_vars(), # temp_daymean, temp_daymax, temp_daymin
    "temp_absmax",
    "temp_absmin", # monthly extreme of the daily extreme
    "hurs_mean",
    "huss_mean"
  ) # relative / specific humidity
}


#' read_clim_parquet
#' Read one admin-unit parquet and mask its missing-value sentinel.
#'
#' @param weighting spatial aggregation across the unit's pixels: `simple`
#'   (unweighted), `static` (population 2000), `dynamic` (population, refreshed
#'   every 10 years). `dynamic` is the closest analogue to the CSV era's
#'   `*_pop_2010` columns.
#'
read_clim_parquet <- function(
  path,
  weighting = c("dynamic", "static", "simple")
) {
  weighting <- match.arg(weighting)

  vars <- c("pr", "tas", "tasmax", "tasmin", "hurs", "huss")
  cols <- stringr::str_c(vars, "_", weighting)

  res <- arrow::read_parquet(path, col_select = c("Date", dplyr::all_of(cols)))

  if (nrow(res) == 0) {
    warning("File ", path, " is empty.")
    return(NULL)
  }

  # Re-select rather than trust arrow's column ordering, then rename positionally
  # so downstream code sees bare `pr` / `tas` / ... regardless of `weighting`.
  res <- dplyr::select(res, Date, dplyr::all_of(cols))
  names(res) <- c("Date", vars)

  # CRITICAL: missing values are the numeric sentinel -9999, NOT NA. Verified on
  # a COM ADM1 unit whose metadata reports `pixels: 0` — all six ERA5 columns are
  # -9999 and is.na() finds nothing. Unmasked, a monthly mean is off by ~10^4.
  # (read_clim_csv() gets this for free via read_csv(na = ...); here it is manual.)
  res %>%
    mutate(across(dplyr::all_of(vars), ~ dplyr::na_if(., -9999))) %>%
    arrange(Date)
}


#' compute_stats_parquet
#' Monthly statistics for one admin unit, from its parquet file.
#'
#' Replaces compute_stats_precip() + compute_stats_temp() in a single pass — the
#' parquet holds precipitation, temperature and humidity together, so there is
#' nothing left to join.
#'
compute_stats_parquet <- function(x, weighting = "dynamic") {
  res <- read_clim_parquet(x, weighting = weighting)

  if (is.null(res)) {
    return(NULL)
  }

  # mean_cum7 / max_cum7 are trailing 7-day cumulative precipitation. The CSV
  # version (compute_stats_precip) wrote
  #   pcum7 = value + pl1 + pl2 + pl3 + pl4 + pl5 + pl6 + pl6
  # -- pl6 twice, eight terms -- so it was neither a 7- nor a 6-day window and
  # double-counted the oldest day. Corrected here to a true 7-day trailing sum,
  # which means new *_cum7 values are NOT comparable to the ones already in
  # climate_stats.sqlite. Recorded as cum7_definition in the provenance sidecar.
  res <- res %>%
    mutate(
      pcum7 = rowSums(do.call(
        cbind,
        lapply(0:6, function(k) dplyr::lag(pr, k))
      ))
    )

  res %>%
    mutate(month = as.Date(format(Date, "%Y-%m-01"))) %>%
    group_by(month) %>%
    summarise(
      # na.rm = FALSE, as in the CSV path: a month containing masked days yields
      # NA rather than a silently partial total.
      precip = sum(pr),
      max_precip = max(pr),
      # na.rm = TRUE, as in the CSV path -- the first six days of the series have
      # no 7-day window. .safe_agg_stat keeps an all-NA month as NA instead of
      # letting max() return -Inf.
      mean_cum7 = .safe_agg_stat(mean, pcum7),
      max_cum7 = .safe_agg_stat(max, pcum7),
      # Monthly MEAN of the daily statistic, matching compute_stats_temp().
      temp_daymean = mean(tas),
      temp_daymax = mean(tasmax),
      temp_daymin = mean(tasmin),
      # New: the monthly extreme of the daily extreme.
      temp_absmax = max(tasmax),
      temp_absmin = min(tasmin),
      hurs_mean = mean(hurs),
      huss_mean = mean(huss),
      .groups = "drop"
    ) %>%
    arrange(month)
}


#' .safe_agg_stat
#' Apply an aggregation with na.rm = TRUE, but return NA (not NaN / -Inf) when
#' every value is missing.
#'
.safe_agg_stat <- function(f, x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(f(x, na.rm = TRUE))
}


#' compile_all_stats_parquet
#' Drop-in replacement for compile_all_stats() reading the parquet source.
#'
#' No join anywhere in this path — that is the whole point.
#'
#' @param gid raw (dotted / Ghana-malformed) gid. The cached file and the
#'   returned `gid` column both use the flattened form, matching the CSV path
#'   and the climate_stats key.
#'
compile_all_stats_parquet <- function(
  gid,
  weighting = "dynamic",
  redo = FALSE,
  clim_dir = "data/climate_parquet/"
) {
  gid_flat <- flatten_gid(gid)
  out_file <- make_clim_stats_file(var = "all_parquet", gid = gid_flat)

  if (!file.exists(out_file) | redo) {
    # The raw gid is required here: flattening erases the dot that
    # makeClimParquetName() splits on.
    clim_file <- get_clim_parquet_gid(gid, clim_dir = clim_dir)

    if (is.null(clim_file)) {
      warning("No parquet file found for ", gid)
      return(NULL)
    }

    cat("-- Computing parquet stats for", gid, "\n")

    res <- compute_stats_parquet(clim_file, weighting = weighting)

    if (is.null(res) || nrow(res) == 0) {
      warning("Parquet file empty for ", gid)
      return(NULL)
    }

    if (anyDuplicated(res$month)) {
      stop(
        "compile_all_stats_parquet(): duplicated months for ",
        gid,
        ". One parquet file must yield one row per month; a duplicate means ",
        "the file itself is malformed. Refusing to return it -- this is the ",
        "class of defect that inflated climate_stats.sqlite to 8.5 GB."
      )
    }

    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # a fresh checkout has no cache dirs
    saveRDS(res, file = out_file)
  } else {
    cat("-- Reading pre-computed parquet stats for", gid_flat, "\n")
    res <- readRDS(out_file)
  }

  cat(
    "---------------------------------
      -- Done parquet stats for ",
    gid_flat,
    "
      ---------------------------------\n"
  )

  res %>%
    mutate(gid = gid_flat)
}


#' compile_all_extremes_parquet
#' Drop-in replacement for compile_all_extremes() reading the parquet source.
#'
#' compute_baseline_stats(), make_baseline_time_bounds() and compute_extremes() are
#' reused UNCHANGED — the extremes maths is source-agnostic, and keeping it
#' shared is what makes an old-vs-new comparison meaningful.
#'
compile_all_extremes_parquet <- function(
  gid,
  redo = FALSE,
  baseline,
  weighting = "dynamic",
  clim_dir = "data/climate_parquet/"
) {
  gid_flat <- flatten_gid(gid)
  out_file <- make_clim_extremes_file(
    var = "all_parquet",
    gid = gid_flat,
    baseline = baseline
  )

  if (!file.exists(out_file) | redo) {
    stats <- compile_all_stats_parquet(
      gid = gid,
      weighting = weighting,
      redo = redo,
      clim_dir = clim_dir
    )

    if (is.null(stats)) {
      warning(
        "Extremes computation failed for ",
        gid,
        " because no parquet data was found."
      )
      return(NULL)
    }

    cat("-- Computing parquet extremes for", gid_flat, baseline, "\n")

    baseline_stats <- compute_baseline_stats(
      gid = gid_flat,
      stats = stats,
      var = "all_parquet",
      redo = redo,
      baseline = baseline
    )

    res <- compute_extremes(
      dat = stats,
      baseline_stats = baseline_stats,
      vars = get_clim_vars_parquet()
    ) %>%
      mutate(gid = gid_flat)

    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE) # a fresh checkout has no cache dirs
    saveRDS(res, file = out_file)
  } else {
    cat("-- Reading pre-computed parquet extremes for", gid_flat, "\n")
    res <- readRDS(out_file)
  }

  cat(
    "----------------------------------------------- \n",
    "-- Done parquet extremes for ",
    gid_flat,
    "\n",
    "----------------------------------------------- \n"
  )

  res
}


#' compute_extremes
#' Computes the extremes using baseline statistics using the sd and the quantiles
#'
#' @param dat the stats data on which to compute extremes
#' @param baseline_stats the baseline statistics to use for definition of extremes
#' @param vars the variables for which to compute extremes
#'
compute_extremes <- function(dat, baseline_stats, vars) {
  dat %>%
    tidyr::pivot_longer(any_of(vars)) %>%
    mutate(month_num = lubridate::month(month)) %>%
    inner_join(baseline_stats, by = c("name", "month_num")) %>%
    mutate(
      hi_sd = value > (mean + 2 * sd),
      hi_q = value > q90,
      lo_sd = value < (mean - 2 * sd),
      lo_q = value < q10,
      hilo_sd = case_when(hi_sd ~ "high", lo_sd ~ "low", T ~ "none"),
      hilo_q = case_when(hi_q ~ "high", lo_q ~ "low", T ~ "none")
    )
}
