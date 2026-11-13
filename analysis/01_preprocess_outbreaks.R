# 01_preprocess_outbreaks.R — sub-weekly outbreak time series → outbreak-level table


# Preamble ----------------------------------------------------------------

library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(lubridate)
library(sf)
library(optparse)
library(furrr)

source("analysis/utils/base.R")
source("analysis/utils/gid.R")          # filterAfrica
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")

# User-supplied options
option_list <- list(
  make_option(c("-r", "--redo"),
              default = TRUE, action ="store", type = "logical", help = "redo data processing"),
  make_option(c("-f", "--make_figures"),
              default = TRUE, action ="store", type = "logical", help = "redo data processing"),
  make_option("--indir",
              default = defaultOutbreakDataDir(),   # utils/paths.R, the one data vintage
              type    = "character",
              help    = "Folder with legacy-format RDS and GPK files [default: %default]"),
  make_option("--drop_cases_le0",
              default = TRUE, action = "store", type = "logical",
              help    = "drop outbreaks with suspected_cases <= 0 (NA kept) [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

print(opt)


# Parallel group_modify -----------------------------------------------------
#
# taxdat::aggregate_single_lp() loops over the rows of one group (roughly
# quadratic in the group size) and every group is independent. This is
# group_modify() with the groups spread over workers: the same keys, the same
# group order (group_by's sorted keys) and the same rows within a group, so the
# result is identical to `group_by(...) %>% group_modify(.f, ...)`.
# `scheduling` is furrr's: Inf sends each group as its own future (few, uneven
# groups), a small number chunks many small groups.
parallelGroupModify <- function(d, groups, .f, ..., scheduling = 1) {
  g    <- dplyr::group_by(d, dplyr::across(dplyr::all_of(groups)))
  keys <- dplyr::group_keys(g)
  res  <- furrr::future_map2(
    dplyr::group_split(g, .keep = FALSE),
    split(keys, seq_len(nrow(keys))),
    .f, ...,
    .options = furrr::furrr_options(scheduling = scheduling, seed = NULL))
  dplyr::bind_rows(lapply(seq_along(res), function(i) {
    dplyr::bind_cols(keys[rep(i, nrow(res[[i]])), , drop = FALSE], res[[i]])
  }))
}

# availableCores() honours the Slurm/cgroup allocation.
plan(multisession, workers = future::availableCores())


# Africa ISO3 filter --------------------------------------------------------
#
# The upstream WHO-region filter (00_ingest_outbreakextractr.R:96-99, `who_region
# %in% c("AFR", "EMR")`) is a WHO *administrative* filter, not a continent filter:
# WHO's EMR region spans North Africa/the Horn AND the Middle East/South Asia, so
# non-African EMR countries (e.g. Afghanistan, Iraq, Lebanon, Pakistan, Syria,
# Yemen) pass it and reach this script. This project models African cholera
# outbreaks only, so we apply a genuine geographic Africa filter here (not in
# the ingest script) -- warning loudly rather than
# silently dropping, so a future re-ingest that pulls in more non-African EMR
# countries doesn't go unnoticed.
#
# Uses the continent lookup already relied on elsewhere in the pipeline
# (06_fit_figures.R's `.afro_sf`) rather than a hand-maintained ISO3 list --
# this correctly keeps WHO-EMR-*administered* but geographically African
# countries (Sudan, Somalia, Egypt, Morocco, Algeria, Libya, Tunisia, Djibouti)
# while dropping genuinely non-African ones.
africa_iso3 <- africaIso3()   # utils/gid.R: giscoR Africa + the synthetic countries


# Load time series --------------------------------------------------------

if (opt$redo) {

  ts_in <- readRDS(file.path(opt$indir, "time_series_preoutbreak_extraction.rds")) %>%
    mutate(location_period_id = as.character(location_period_id)) %>%
    # get_admin_level() is str_count(., "::") - 1: already vectorised.
    mutate(admin_level = taxdat::get_admin_level(location))
  ts <- filterAfrica(ts_in, africa_iso3, "01/coverage-ts")

  # Aggregate
  aggregated_ts <- ts %>%
    parallelGroupModify(c("country", "location", "location_period_id", "admin_level"),
                        .f = taxdat::aggregate_single_lp,
                        verbose = FALSE,
                        cases_column = "sCh",
                        scheduling = 4) %>%
    mutate(admin_level = str_c("ADM", taxdat::get_admin_level(location))) %>%
    # TZA::Mainland is a split-country unit that taxdat treats as ADM1; shift
    # all Mainland locations down by 1 so Mainland itself becomes ADM0 and its
    # sub-locations become ADM1, ADM2, etc.
    mutate(admin_level = if_else(
      str_detect(location, "TZA::Mainland"),
      str_c("ADM", as.integer(str_extract(admin_level, "[0-9]+")) - 1L),
      admin_level
    ))

  saveRDS(aggregated_ts, file = makeAggTsFile())
  
  # Aggregate data by admin level
  aggregated_admin_level_ts <- aggregated_ts %>%
    parallelGroupModify(c("country", "admin_level"),
                        .f = taxdat::aggregate_single_lp,
                        verbose = FALSE,
                        cases_column = "sCh",
                        scheduling = Inf)
  
  saveRDS(aggregated_admin_level_ts, file = makeAggAdmTsFile())

} else {
  aggregated_ts <- readRDS(makeAggTsFile())
  aggregated_admin_level_ts <- readRDS(makeAggAdmTsFile())
}

if (opt$make_figures) {
  p_coverage <- aggregated_ts %>% 
    ggplot(aes(y = location)) +
    geom_linerange(aes(xmin = TL, xmax = TR)) +
    facet_grid(country + admin_level ~ ., scales = "free") +
    theme_bw() +
    scale_x_date() +
    theme(axis.text.y = element_blank())
  
  ggsave(p_coverage, filename = "figures/data_coverage.pdf", width = 12, height = 50, limitsize = FALSE)
  
}

# Load outbreaks ----------------------------------------------------------

if (opt$redo) {
  ob_in <- readRDS(file.path(opt$indir, "time_series_outbreak_extraction.rds")) %>%
    # Create unique outbreak identifiters
    mutate(uid = str_c(location, "-", outbreak_number))
  ob_ts <- ob_in %>%
    filterAfrica(africa_iso3, "01/outbreak-metadata") %>%
    group_by(uid, location, location_period_id, country) %>%
    summarise(TL = min(TL),
              TR = max(TR),
              suspected_cases = sum(sCh),
              # Confirmed cases: recorded, not modelled. NA when no week of the
              # outbreak reports cCh; otherwise the sum over the weeks that do
              # (n_weeks_confirmed says how many).
              confirmed_cases = if (all(is.na(cCh))) NA_real_ else sum(as.numeric(cCh), na.rm = TRUE),
              n_weeks_confirmed = sum(!is.na(cCh)),
              pop = mean(pop),
              threshold = mean(threshold)) %>% 
    mutate(admin_level = str_c("ADM", taxdat::get_admin_level(location))) %>%
    mutate(admin_level = if_else(
      str_detect(location, "TZA::Mainland"),
      str_c("ADM", as.integer(str_extract(admin_level, "[0-9]+")) - 1L),
      admin_level
    )) %>%
    taxdat::get_AFRO_region(ctry_col = "country") %>%
    fillSyntheticAfroRegion() %>%   # utils/gid.R: the synthetic countries (no-op for real data)
    mutate(admin_level = factor(admin_level, str_c("ADM", 0:3))) %>% 
    mutate(year = lubridate::year(TL),
           month = lubridate::month(TL))

  # Non-positive outbreaks: the upstream
  # >= 50-case filter lets an outbreak with suspected_cases <= 0 through, and it
  # becomes a canonical episode of every track. Same rule as
  # dur_data.R::outbreakEpisodeGuardFlags()$cases_le0 (NA never dropped).
  cases_le0 <- !is.na(ob_ts$suspected_cases) & ob_ts$suspected_cases <= 0
  if (isTRUE(opt$drop_cases_le0) && any(cases_le0)) {
    warning(sprintf("[01/outbreak-metadata] Dropping %d outbreak(s) with suspected_cases <= 0: %s",
                    sum(cases_le0), paste(ob_ts$uid[cases_le0], collapse = ", ")), call. = FALSE)
    ob_ts <- ob_ts[!cases_le0, ]
  }

  conf_gt_susp <- !is.na(ob_ts$confirmed_cases) & ob_ts$confirmed_cases > ob_ts$suspected_cases
  if (any(conf_gt_susp))
    warning(sprintf("[01/outbreak-metadata] %d outbreak(s) report more confirmed than suspected cases (kept): %s",
                    sum(conf_gt_susp), paste(ob_ts$uid[conf_gt_susp], collapse = ", ")), call. = FALSE)

  saveRDS(ob_ts, file = makeOutbreakMetadataFile())

  # Processing cascade (table S1): before / after the Africa filter. The
  # coverage series counts are written only when --redo rebuilt them above.
  writeStepCounts("01_preprocess_outbreaks", c(
    if (exists("ts_in")) list(
      n_lp_series_in         = dplyr::n_distinct(ts_in$location_period_id),
      n_countries_series_in  = dplyr::n_distinct(ts_in$country),
      n_lp_series_africa     = dplyr::n_distinct(ts$location_period_id),
      n_countries_series_africa = dplyr::n_distinct(ts$country)),
    list(
      n_outbreaks_in         = dplyr::n_distinct(ob_in$uid),
      n_countries_outbreaks_in = dplyr::n_distinct(ob_in$country),
      n_outbreaks_cases_le0  = sum(cases_le0),
      drop_cases_le0         = as.integer(isTRUE(opt$drop_cases_le0)),
      n_outbreaks_africa     = nrow(ob_ts),
      n_outbreaks_confirmed_reported = sum(!is.na(ob_ts$confirmed_cases)),
      n_outbreaks_confirmed_gt_suspected = sum(conf_gt_susp),
      n_lp_outbreaks_africa  = dplyr::n_distinct(ob_ts$location_period_id),
      n_countries_outbreaks_africa = dplyr::n_distinct(ob_ts$country))),
    indir = opt$indir)
} else {
  ob_ts <- readRDS(makeOutbreakMetadataFile())
}


if (opt$redo) {
  # Load shapefiles
  # Only the outbreak LPs reach the gpkg and the maps (both inner-join ob_ts),
  # so repair those alone instead of every LP of the shapefile corpus.
  ob_shp <- readRDS(file.path(opt$indir, "outbreak_shapefiles.rds")) %>%
    filter(lctn_pr %in% ob_ts$location_period_id) %>%
    st_make_valid() %>%
    st_simplify()

  # write to geopackage
  ob_shp %>%
    inner_join(ob_ts, by = c("lctn_pr" = "location_period_id")) %>%
    st_write(file.path(opt$indir, "outbreak_shapefiles.gpkg"),
             append = FALSE)
} else {
  ob_shp <- st_read(file.path(opt$indir, "outbreak_shapefiles.gpkg"))
}


if (opt$make_figures) {
  
  # plot map
  p_shp <- ob_shp %>%
    select(lctn_pr) %>%
    inner_join(ob_ts, by = c("lctn_pr" = "location_period_id")) %>%
    mutate(decade = round(lubridate::year(TL)/5)*5) %>%
    ggplot() +
    geom_sf(alpha = 0) +
    # facet_grid(admin_level~decade) +
    taxdat::map_theme()
  
  
  ggsave(p_shp, filename = "figures/outbreak_map.png", width = 10, height = 8)
}


if (opt$make_figures) {
  
  africa_list <- c(
    "SEN", "MRT", "MLI", "BFA", "NER", "TCD",    NA, "SDN",    NA,    NA,
    "SLE", "LBR",  "TGO", "BEN", "NGA", "CMR", "CAF", "SSD", "ETH", "SOM",
    "GIB", "GIN", "CIV", "GHA", "GNB", "COG", "COD", "UGA", "KEN",    NA,
    NA,    NA,    NA,    NA,   "AGO", "ZMB", "RWA", "BDI", "TZA", "MWI",
    NA,    NA,    NA,    NA,   "NAM", "BWA", "ZWE", "ZAF", "MOZ",    NA
  )
  
  dat <-  ob_shp %>%
    select(lctn_pr) %>%
    inner_join(ob_ts, by = c("lctn_pr" = "location_period_id")) %>%
    mutate(decade = round(lubridate::year(TL)/5)*5) %>% 
    filter(country %in% africa_list)
  
  # plot map
  p_shp <- dat %>%
    ggplot(aes(fill = admin_level)) +
    geom_sf(data = filter(dat, admin_level == "ADM0")) +
    geom_sf(data = filter(dat, admin_level == "ADM1")) +
    geom_sf(data = filter(dat, admin_level == "ADM2")) +
    geom_sf(data = filter(dat, admin_level == "ADM3")) +
    geom_sf(data = filter(dat, admin_level == "ADM0"), alpha = 0, lwd = 1.2) +
    # facet_grid(admin_level~decade) +
    taxdat::map_theme()
  
  
  ggsave(p_shp, filename = "figures/outbreak_map_admin.png", width = 10, height = 8)
}


if (opt$make_figures) {
  p_outbreaks <- ob_ts %>% 
    ggplot(aes(y = country)) +
    geom_linerange(data = aggregated_admin_level_ts,
                   aes(xmin = TL, xmax = TR),
                   color = "darkgray") +
    geom_point(aes(x = TL), fill = "white", color = "darkgray", pch = 21, size = 2.3) +
    geom_point(aes(x = TL, color = month), pch = 4) +
    facet_grid(admin_level ~  ., scales = "free", space = "free") +
    theme_bw() +
    scale_x_date() + 
    scale_color_gradientn(colours =c("red", "orange", "blue", "purple", "red"))
  # scale_color_viridis_c(option = "H")
  
  
  ggsave(p_outbreaks, filename = "figures/outbreak_occurrence.pdf", width = 12, height = 20, limitsize = FALSE)
  
  
  ob_ts %>% 
    mutate(year = year(TL)) %>% 
    group_by(year, AFRO_region) %>% 
    ggplot(aes(x = year)) +
    geom_bar() +
    facet_grid(AFRO_region ~. )
  
  
  # Look at seasonality of outbreaks
  ob_ts %>% 
    ggplot(aes(x = month, color = AFRO_region)) +
    geom_density() +
    facet_grid(admin_level ~ .) +
    theme_bw() +
    scale_x_continuous(breaks = c(1:12))
  
  
  p_outbreak_seas_country <- ob_ts %>% 
    ggplot(aes(x = month, color = admin_level)) +
    geom_density(bw = 1.5, aes(lty = admin_level)) +
    facet_wrap(country ~ ., scales = "free_y") +
    theme_bw() +
    scale_x_continuous(breaks = c(1:12)) +
    theme(panel.grid.minor = element_blank()) +
    scale_color_viridis_d(option = "H")
  
  ggsave(p_outbreak_seas_country, 
         filename = "figures/outbreak_seasonality_country.png", width = 12, height = 7)
  
  africa_list <- c(
    "SEN", "MRT", "MLI", "BFA", "NER", "TCD",    NA, "SDN",    NA,    NA,
    "SLE", "LBR",  "TGO", "BEN", "NGA", "CMR", "CAF", "SSD", "ETH", "SOM",
    "GIB", "GIN", "CIV", "GHA", "GNB", "COG", "COD", "UGA", "KEN",    NA,
    NA,    NA,    NA,    NA,   "AGO", "ZMB", "RWA", "BDI", "TZA", "MWI",
    NA,    NA,    NA,    NA,   "NAM", "BWA", "ZWE", "ZAF", "MOZ",    NA
  )
  
  plot_list <- map(africa_list, function(cntry) {
    if (is.na(cntry)) {
      return(NULL)
    }
    
    dat <- ob_ts %>% 
      filter(country == cntry)
    
    if (nrow(dat) == 0) {
      dat <- ob_ts %>% 
        ungroup() %>% 
        slice(1:2) %>% 
        mutate(month = c(1, 12),
               country = cntry)
      draw <- FALSE
    } else {
      draw <- "TRUE"
    }
    
    
    p <- dat %>% 
      ggplot(aes(x = month, color = admin_level)) +
      facet_wrap(country ~ ., scales = "free_y") +
      theme_bw() +
      scale_x_continuous(breaks = c(1, 6, 12), limits = c(1, 12)) +
      theme(panel.grid.minor = element_blank()) +
      scale_color_viridis_d(option = "H", drop=FALSE)  +
      guides(color = "none", lty = "none")
    
    if (draw) {
      p <- p + geom_density(bw = 1.5) +
        guides(color = "none", lty = "none")
    }
    
    if (!(cntry %in% c("GIN", "CIV", "GHA", "NAM", "BWA", "ZAF", "ZWE", "MOZ", "SOM"))) {
      p <- p + theme(axis.text.x  = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank())
    }
    
    if (!(cntry %in% c("SLE", "GIB", "AGO", "NAM", "SDN"))) {
      p <- p  + 
        theme(axis.text.y  = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())
    } else {
      p <- p  + 
        theme(axis.text.y  = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    p  + theme(plot.margin = unit(c(0, -.3, 0, -.3), "cm"))
    
  })
  
  p_outbreak_seas <- cowplot::plot_grid(plotlist = plot_list,
                                        ncol = 10,
                                        align = "hv",
                                        axis = "tblr") +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  ggsave(p_outbreak_seas, filename = "figures/outbreak_seasonality_country_map.png",
         width = 12, height = 7)
}
