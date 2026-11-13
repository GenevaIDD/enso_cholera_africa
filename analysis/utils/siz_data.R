
makeSizeData <- function(x, baseline = "1990_2019") {

  # Get clim data
  clim_stats <- getStatsGid(gid = x$gid[1])
  clim_extrms <- getExtremesGid(gid = x$gid[1], baseline = baseline)
  
  clim_extrms_cat <- clim_extrms %>% 
    select(-hilo_sd, -mean, -sd) %>% 
    mutate(name = stringr::str_c("hlq_", name)) %>% 
    tidyr::pivot_wider(names_from = "name",
                values_from = "hilo_q")
  
  clim_dev <- clim_stats %>% 
    tidyr::pivot_longer(cols = c(getPrecipVars(), getTempVars())) %>% 
    inner_join(clim_extrms %>% select(month, name, mean, sd), by = c("month", "name")) %>%
    mutate(dev = value - mean,
           dev_sd = (value - mean)/sd) %>% 
    select(gid, month, name, mean, dev, dev_sd) %>% 
    tidyr::pivot_wider(names_from = "name",
                values_from = c("mean", "dev", "dev_sd"))
  
  x %>% 
    mutate(TL_month = floor_date(TL, unit = "month")) %>%   # onset month (as ENSO / dur / mon), not the nearest
    inner_join(clim_stats, by = c("gid" = "gid", "TL_month" = "month")) %>% 
    inner_join(clim_extrms_cat, by = c("gid" = "gid", "TL_month" = "month")) %>% 
    inner_join(clim_dev, by = c("gid" = "gid", "TL_month" = "month"))
}

#' .makeSizeClimStats (private helper for makeClimateLags mode = "event")
#'
#' Fetches all available monthly climate statistics for a vector of gids from
#' both SQLite databases, computes baseline deviations (same logic as
#' makeSizeData()), and returns a wide-format tibble (one row per gid × month)
#' with mean_* and dev_* columns in **unscaled** units.  readSizeData() applies
#' the /100 / /10 column-name-pattern scaling at load time, which automatically
#' covers any <var>_lagk columns stored in the same RDS.
#'
#' @param gids     character vector of GADM gids
#' @param baseline character, e.g. "1990_2019"
#' @return tibble (gid, month [Date], mean_precip, dev_precip, ...)
.makeSizeClimStats <- function(gids, baseline) {
  clim_stats  <- getStatsGid(gids)
  clim_extrms <- getExtremesGid(gids, baseline = baseline)

  clim_stats %>%
    tidyr::pivot_longer(cols       = c(getPrecipVars(), getTempVars()),
                        names_to   = "name",
                        values_to  = "value") %>%
    dplyr::inner_join(
      clim_extrms %>% dplyr::select(gid, month, name, mean, sd),
      by = c("gid", "month", "name")
    ) %>%
    dplyr::mutate(dev = value - mean) %>%
    dplyr::select(gid, month, name, mean, dev) %>%
    tidyr::pivot_wider(names_from  = "name",
                       values_from = c("mean", "dev"))
}

#' .addLogPopCentred
#'
#' P4. `log_pop_c = log(pop) - mean(log(pop))`, centred WITHIN `admin_level`.
#'
#' The centring is what makes `beta = 1` reproduce the offset-only model
#' exactly: uncentred, the deviation term would move the intercept too, and the
#' `fixed` and `free` paths would no longer be nested. Grouping by level (rather
#' than centring globally) is what makes each level's coefficient read as that
#' level's own departure from proportionality under `free_by_level`, and is
#' harmless under `free`.
#'
#' Centres are returned as an attribute because the SIMULATOR needs them: it
#' must centre frame units on the TRAINING centres, not on the frame's own mean,
#' or the elasticity is applied around the wrong pivot and the extrapolation is
#' silently rescaled.
.addLogPopCentred <- function(d, tag = "P4") {
  if (!"pop" %in% names(d)) stop(".addLogPopCentred: no `pop` column.", call. = FALSE)
  lvl <- if ("admin_level" %in% names(d)) as.character(d$admin_level)
         else rep("ALL", nrow(d))
  lp  <- log(as.numeric(d$pop))
  ok  <- is.finite(lp)
  ctr <- tapply(lp[ok], lvl[ok], mean)
  d$log_pop_c <- lp - unname(ctr[lvl])
  attr(d, "log_pop_centres") <- ctr
  message(sprintf("[%s] log_pop_c centres: %s", tag,
                  paste(sprintf("%s=%.4f", names(ctr), ctr), collapse = " ")))
  d
}

#' readSizeData
#'
#' Load + filter + scale outbreak-size analysis data. Ports the data-prep
#' logic from legacy 07_run_outbreak_size_model.R:26-63.
readSizeData <- function(opt, path = makeSizeDatFile(opt)) {
  if (!file.exists(path)) stop("readSizeData: missing input ", path)
  d <- readRDS(path)
  ep_rec <- assertEpisodeKey(opt, path, "size")

  # P0: resolve and validate the estimand ONCE, up front, so a mis-declared
  # config fails before any data work rather than producing a plausible-looking
  # fit on mixed units.
  estimand <- resolveSizeEstimand(opt)

  # Categorical time-to-last-outbreak (legacy bins).
  d <- d %>%
    dplyr::mutate(
      time_to_last_outbreak = cut(
        (pmax(0, time_to_last_outbreak) + 1) / 365,
        breaks         = c(0, 1, 2, 3, 4, 5, 6, Inf),
        include.lowest = TRUE, right = FALSE
      )
    ) %>%
    tidyr::replace_na(list(time_to_last_outbreak = "[6,Inf]"))

  # Column-class-based scaling (matches legacy regex behaviour).
  cols      <- colnames(d)
  precip_cols <- setdiff(grep("precip|cum7", cols, value = TRUE), grep("hlq", cols, value = TRUE))
  temp_all    <- setdiff(grep("temp", cols, value = TRUE),         grep("hlq", cols, value = TRUE))
  temp_cols     <- setdiff(temp_all, grep("dev", temp_all, value = TRUE))
  temp_dev_cols <- intersect(temp_all, grep("dev", temp_all, value = TRUE))
  hlq_cols      <- grep("hlq", cols, value = TRUE)

  if (length(precip_cols))   d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(precip_cols),   ~ . / 100))
  if (length(temp_cols))     d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(temp_cols),     ~ . / 10))
  if (length(temp_dev_cols)) d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(temp_dev_cols), ~ . / 1))
  if (length(hlq_cols))      d <- dplyr::mutate(d, dplyr::across(dplyr::all_of(hlq_cols),      ~ . == "high"))

  d <- d %>%
    dplyr::group_by(uid) %>%
    dplyr::mutate(suspected_cases = round(suspected_cases)) %>%
    dplyr::ungroup()

  # Country filter ("all" = no filter)
  ctry <- opt$country %||% "all"
  if (!identical(ctry, "all")) {
    keep <- trimws(strsplit(ctry, ",")[[1]])
    d <- dplyr::filter(d, country %in% keep)
  }

  # Minimum-population guard: drop locations whose WorldPop population is
  # implausibly small (broken raster extraction, e.g. pop ~ 5e-6 for COD LP
  # 14306). Such rows give meaningless incidence and would dominate mean_rate.
  # Locations with an *unknown* (NA) population are kept — they are not "below"
  # the floor and the size response (raw suspected_cases) does not use pop.
  # Threshold configurable via opt$size$min_pop (default 10 people).
  min_pop <- opt$size$min_pop %||% 10

  # Maximum-attack-rate guard: drop rows whose implied attack rate
  # (suspected_cases / pop) exceeds an epidemiologically plausible ceiling.
  # Such rows signal a broken/collapsed WorldPop population (e.g. pop ~ 892 for
  # several Katanga health districts pooled into COD_3_6_1) rather than real
  # incidence: the size fit would learn a huge compensating ADM2 random effect
  # on the tiny training pop, which then blows up mu = pop * mean_rate * exp(lp)
  # at simulation time (the per-gid frame pop is the sane average). Rows with an
  # *unknown* (NA) or non-positive population are kept (handled by min_pop / the
  # mean_rate pop > 0 guard). Threshold configurable via
  # opt$size$max_attack_rate (default 0.10, i.e. 10% of the catchment).
  max_ar <- opt$size$max_attack_rate %||% 0.10

  # P1: point `pop` at the declared denominator BEFORE the population guards,
  # so min_pop / max_attack_rate are evaluated on the population the model will
  # actually use. Guarding one population and fitting another would admit rows
  # that are implausible under the fitted denominator.
  d <- resolveSizePopColumn(d, estimand)

  # P3: re-key `admin_level` to the canonical GADM level BEFORE
  # applyAdminLevelFilter(). Order is load-bearing and this is the ONLY correct
  # position for it:
  #   - after the filter, the filter would have partitioned rows on corpus
  #     labels while mean_rate pooled them on canonical ones -- two keyings on
  #     one dataset, which is the defect being fixed;
  #   - before it, the ADM1/ADM2 overlap logic and the mean_rate cells agree by
  #     construction, which is also what closes the P1' filter-divergence
  #     residue (the fit's filterOverlappingADM1() reads the size data's own
  #     composition while dropOverlappingADM1FromFrame() reads the survival
  #     data's -- once both sides are canonically keyed they cannot disagree).
  # Gated on the estimand block being opted in, so pre-P0 configs are untouched.
  d <- .applyP3(d, opt, "size", estimand, tag = "readSizeData/P3")

  # Apply admin-level row filter; drop missing TL/TR; drop sub-floor populations;
  # drop implausible attack rates.
  d <- d %>%
    applyAdminLevelFilter(admin_level_opt = opt$size$admin_level) %>%
    dplyr::filter(!is.na(TL), !is.na(TR)) %>%
    dplyr::filter(is.na(pop) | pop >= min_pop) %>%
    dplyr::filter(is.na(pop) | !(pop > 0) | suspected_cases <= max_ar * pop) %>%
    dplyr::ungroup()

  if (!nrow(d)) stop("readSizeData: empty data after filters (country=", ctry, ")")

  # Episode identity + duration, computed identically to the duration track
  # (analysis/utils/dur_data.R: makeDurationPanel() / .expandActiveMonths()) so
  # a `log(duration_months)` covariate here means the same thing as the
  # duration model's `duration_months`, and `episode_id` joins 1:1 onto the
  # duration panel for the Phase-5 size-duration dependence check
  # (analysis/05_fit_postprocess_dependence.R). `duration_months` is the count
  # of CALENDAR months intersected by [TL, TR] (n_months), not TR-TL in days /
  # ~30 -- the two differ for episodes that straddle a month boundary, and only
  # the calendar-month count matches the duration panel's row count per episode.
  d <- d %>%
    dplyr::mutate(
      episode_id      = sprintf("%s|%s|%s", uid, format(as.Date(TL)), format(as.Date(TR))),
      duration_days   = as.numeric(as.Date(TR) - as.Date(TL)),
      duration_months = (lubridate::year(lubridate::floor_date(as.Date(TR), "month")) -
                          lubridate::year(lubridate::floor_date(as.Date(TL), "month"))) * 12L +
                         (lubridate::month(lubridate::floor_date(as.Date(TR), "month")) -
                          lubridate::month(lubridate::floor_date(as.Date(TL), "month"))) + 1L)

  # mean_rate is a country x admin_level table (buildMeanRateTable()), not a
  # single global scalar.
  # The table is built on the FULL corpus (all countries/levels; buildMeanRateTable()
  # applies the same min_pop/max_attack_rate guards internally, ahead of any
  # opt$country subsetting done above), so a per-country fit still resolves against
  # the global cell/group/global fallback tiers rather than a shrunken sample of one.
  #
  # Built through the SHARED builder (utils/anchor.R) rather than the size
  # original `buildMeanRateTable()`: the original hardcodes the geometric
  # estimator, and size's default is `pooled` (see config_defaults.yml
  # size.mean_rate.estimator). Under
  # `estimator: geometric` the two are proven identical cell-for-cell, so the
  # original stays valid and side-by-side. `attachBaseRate()` writes the same
  # `mean_rate` / `mr_tier` columns and the same `attr(d, "mean_rate_table")`
  # that persistMeanRateTable() and tests/gates_estimand_coherence.R read.
  mrt <- buildBaseRateTable(opt, "siz", path = path)
  d   <- attachBaseRate(d, opt, "siz", brt = mrt)
  d$ADM2 <- d$gid
  # re_adm2_leaf_scope: adm2_only -> is_adm2 / gid_leaf (utils/paths.R).
  d <- .leafColumnsFor(d, opt$size, ep_rec, tag = "readSizeData")

  # P0: the mean_rate table and this data must declare the SAME estimand. The
  # table is built by an independent call that re-reads `path`, so this is a
  # genuine cross-check, not a tautology -- it is the seam where `Q x rbar^H`
  # would appear.
  assertSizeEstimand(sizeEstimandOf(mrt), sizeEstimandStamp(opt),
                     what = "readSizeData: mean_rate table vs fit data",
                     require_stamp = FALSE)

  # P4: the elasticity deviation column. Built AFTER all row filters so the
  # centres are the means of the rows actually fitted -- centring on the
  # unfiltered corpus would put the pivot outside the training set, and the
  # simulator (which reuses these centres) would extrapolate around a
  # population no row ever had.
  if (!identical(estimand$pop_elasticity, "fixed"))
    d <- .addLogPopCentred(d, tag = "readSizeData/P4")

  attr(d, "size_estimand") <- sizeEstimandStamp(opt, extra = list(
    n_rows        = nrow(d),
    cases_total   = sum(d$suspected_cases, na.rm = TRUE),
    admin_level   = opt$size$admin_level %||% NA_character_,
    log_pop_centres = attr(d, "log_pop_centres")
  ))
  attr(d, "episode_stamp") <- ep_rec
  d
}
