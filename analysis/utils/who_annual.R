# who_annual.R — WHO annual cholera case reports: country-year totals and
#   the exact-zero / below-threshold mask
#
# Provides:
#   build_who_annual_zero_mask()   country-year WHO cases + is_zero mask
#   build_who_annual_cases()       country-year WHO annual cases
#   who_window_totals()            WHO totals over a window
#
# Details:
#   Reads data/who_annual_report_africa/ (annual report + dashboard annual
#   totals). Used by 01_preprocess_national_zero_months.R (via
#   national_zero.R) and the WHO-vs-model results.

#' build_who_annual_zero_mask
#'
#' Country-year mask of whether annual suspected cholera cases are exactly
#' zero, combining two WHO sources in `data/who_annual_report_africa/`:
#'
#' * the WHO annual report (`who_annual_report_africa.rds`, `sCh`, 1970-2024);
#' * the WHO dashboard annual totals (`who_dashboard_annual_africa.rds`,
#'   `case_total` by `report_year`, 2023 onwards; the latest year is partial).
#'
#' Within the annual report, only `spatial_scale == "country"` rows are kept
#' (drops the Mayotte admin1 row) and each country-year takes the maximum
#' `sCh` over its rows. A few country-years carry more than one row: exact
#' duplicates (EGY 2022) and overlapping partial reports of the same period
#' where one row is the sum of the others (TZA 2022: 106, 75, 181), so the
#' maximum avoids double counting. The year is the calendar year of `TL`
#' (no annual-report row spans two calendar years).
#'
#' Where a country-year is present in both sources, `suspected_cases` is the
#' mean of the two values and `suspected_cases_max` their maximum (the
#' conservative value for a below-threshold test: a year is under the
#' threshold only if every WHO source says so).
#'
#' @param annual_report_path path to `who_annual_report_africa.rds`.
#' @param dashboard_path path to `who_dashboard_annual_africa.rds`.
#' @return tibble with one row per country-year: `country` (ISO3), `year`
#'   (integer), `suspected_cases` (numeric), `is_zero` (logical,
#'   `suspected_cases == 0`), `suspected_cases_max` (numeric).
build_who_annual_zero_mask <- function(
  annual_report_path = "data/who_annual_report_africa/who_annual_report_africa.rds",
  dashboard_path = "data/who_annual_report_africa/who_dashboard_annual_africa.rds"
) {
  report <- readRDS(annual_report_path) |>
    dplyr::filter(spatial_scale == "country") |>
    dplyr::transmute(
      country = iso3,
      year = as.integer(format(TL, "%Y")),
      cases = as.numeric(sCh)
    ) |>
    dplyr::summarise(cases = max(cases), .by = c(country, year))

  dashboard <- readRDS(dashboard_path) |>
    dplyr::transmute(
      country = iso_3_code,
      year = as.integer(report_year),
      cases = as.numeric(case_total)
    )

  dplyr::bind_rows(report, dashboard) |>
    dplyr::summarise(
      suspected_cases = mean(cases),
      suspected_cases_max = max(cases),
      .by = c(country, year)
    ) |>
    dplyr::mutate(is_zero = suspected_cases == 0) |>
    dplyr::relocate(suspected_cases_max, .after = is_zero) |>
    dplyr::arrange(country, year)
}

#' build_who_annual_cases
#'
#' Annual suspected cholera cases per country-year for a comparison of years
#' (09 report_who_enso_ratio). Unlike
#' build_who_annual_zero_mask(), an absent country-year is NA, never zero: a
#' country missing from the WHO dashboard is not reporting zero cases.
#'
#' * WHO annual report: `spatial_scale == "country"` rows, the maximum `sCh`
#'   per country-year (as build_who_annual_zero_mask()), year of `TL`.
#' * WHO dashboard weekly series: summed by its epidemiological year (`epiyr`),
#'   so the year starts with epiweek 1 (the dashboard's own annual file starts
#'   2026 on 2025-12-29, a week the weekly `epiyr` files under 2025).
#' * Both present: `cases` is their mean (as build_who_annual_zero_mask()), both
#'   values kept.
#' * `target_year` is partial (the weekly series to its last week, `cutoff`).
#'   A country whose `target_year` weekly series stops more than `censor_weeks`
#'   weeks before `cutoff` is censored (WHO stopped receiving its data): its
#'   `target_year` value is NA and kept as `cases_censored`.
#'
#' @param report    the WHO annual report table (who_annual_report_africa.rds)
#' @param weekly    the WHO dashboard weekly table (who_dashboard_weekly_africa.rds)
#' @param countries ISO3 codes to keep
#' @param years     integer years to tabulate (absent country-years are NA rows)
#' @param target_year the partial year (its weekly series defines `cutoff`)
#' @param censor_weeks weeks before `cutoff` beyond which a stopped series is censored
#' @return list(annual = tibble(country, year, cases_report, cases_dashboard,
#'   n_weeks, last_week, source, cases, partial_year, censored, cases_censored),
#'   cutoff = the last weekly week (Monday) of `target_year`, censored = ISO3)
build_who_annual_cases <- function(
  report,
  weekly,
  countries,
  years,
  target_year,
  censor_weeks = 8L
) {
  countries <- sort(unique(as.character(countries)))
  years <- sort(unique(as.integer(years)))
  target_year <- as.integer(target_year)
  rep_tab <- tibble::as_tibble(report) |>
    dplyr::filter(
      as.character(spatial_scale) == "country",
      iso3 %in% countries
    ) |>
    dplyr::transmute(
      country = as.character(iso3),
      year = as.integer(format(as.Date(TL), "%Y")),
      cases = as.numeric(sCh)
    ) |>
    dplyr::summarise(cases_report = max(cases), .by = c(country, year))
  wk <- tibble::as_tibble(weekly) |>
    dplyr::filter(iso_3_code %in% countries) |>
    dplyr::transmute(
      country = as.character(iso_3_code),
      date_wk = as.Date(date_wk),
      cases = as.numeric(cases),
      year = as.integer(epiyr)
    )
  dup <- dplyr::filter(dplyr::count(wk, country, date_wk), n > 1L)
  if (nrow(dup)) {
    stop(
      "build_who_annual_cases: ",
      nrow(dup),
      " duplicate (country, week) rows in the weekly series, e.g. ",
      dup$country[1],
      " ",
      format(dup$date_wk[1]),
      call. = FALSE
    )
  }
  if (anyNA(wk$cases)) {
    stop(
      "build_who_annual_cases: NA cases in the weekly series.",
      call. = FALSE
    )
  }
  dash <- wk |>
    dplyr::summarise(
      cases_dashboard = sum(cases),
      n_weeks = dplyr::n(),
      last_week = max(date_wk),
      .by = c(country, year)
    )
  tw <- wk$date_wk[wk$year == target_year]
  if (!length(tw)) {
    stop(
      "build_who_annual_cases: no weekly rows in ",
      target_year,
      call. = FALSE
    )
  }
  cutoff <- max(tw)
  cens <- dash$country[
    dash$year == target_year &
      dash$last_week < cutoff - 7L * as.integer(censor_weeks)
  ]
  dash_years <- sort(unique(wk$year))
  annual <- tidyr::expand_grid(country = countries, year = years) |>
    dplyr::left_join(rep_tab, by = c("country", "year")) |>
    dplyr::left_join(dash, by = c("country", "year")) |>
    dplyr::mutate(
      source = dplyr::case_when(
        !is.na(cases_report) & !is.na(cases_dashboard) ~ "report+dashboard",
        !is.na(cases_report) ~ "annual_report",
        !is.na(cases_dashboard) ~ "dashboard_weekly",
        year %in% dash_years ~ "not_on_dashboard",
        TRUE ~ "not_in_report"
      ),
      cases = dplyr::case_when(
        source == "report+dashboard" ~ (cases_report + cases_dashboard) / 2,
        source == "annual_report" ~ cases_report,
        source == "dashboard_weekly" ~ cases_dashboard,
        TRUE ~ NA_real_
      ),
      partial_year = year == target_year,
      censored = partial_year & country %in% cens,
      cases_censored = ifelse(censored, cases, NA_real_),
      cases = ifelse(censored, NA_real_, cases)
    ) |>
    dplyr::arrange(country, year)
  list(annual = annual, cutoff = cutoff, censored = sort(cens))
}

#' who_window_totals
#'
#' Cases of each country-year over a set of calendar months (e.g. May-August)
#' from the WHO monthly table of 03_prepare_who_monthly.R, kept only where every
#' one of those months is complete (`partial` FALSE): a year with a missing or
#' partial month has no window total.
#' @param monthly   the WHO monthly table (iso3, t_start, cases, partial)
#' @param countries ISO3 codes to keep
#' @param months    integer calendar months of the window
#' @return tibble(country, year, cases)
who_window_totals <- function(monthly, countries, months) {
  months <- sort(unique(as.integer(months)))
  tibble::as_tibble(monthly) |>
    dplyr::filter(
      iso3 %in% countries,
      as.integer(format(as.Date(t_start), "%m")) %in% months
    ) |>
    dplyr::mutate(
      country = as.character(iso3),
      year = as.integer(format(as.Date(t_start), "%Y"))
    ) |>
    dplyr::summarise(
      cases = sum(as.numeric(cases)),
      n_complete = sum(!as.logical(partial)),
      .by = c(country, year)
    ) |>
    dplyr::filter(n_complete == length(months)) |>
    dplyr::select(country, year, cases) |>
    dplyr::arrange(country, year)
}
