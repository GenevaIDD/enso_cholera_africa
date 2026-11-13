

# Postprocessing ----------------------------------------------------------

getCovarType <- function(x) {
  if (is.na(x)) return("other")
  if (stringr::str_detect(x, "months")) {
    "months"
  } else if (x %in% c(getPrecipVars(), "cum_precip") || (x %in% getTempVars())) {
    "climate_var"
  }  else if (stringr::str_detect(x, "hlq")) {
    "climate anomaly"
  } else if (stringr::str_detect(x, "mean_")) {
    "climate mean"
  } else if (stringr::str_detect(x, "dev_")) {
    "climate anomaly"
  } else if (stringr::str_detect(x, "time_to_last")) {
    "years to last outbreak"
  } else if (stringr::str_detect(x, "enso")) {
    "enso teleconnection"
  } else {
    "other"
  }
}


# Plots -------------------------------------------------------------------

get_covar_levels <- function(surv_dat_all) {

  all_levels <- c(stringr::str_c("months", 1:12),
                  getPrecipVars(), stringr::str_c("hlq_", getPrecipVars(), "TRUE"),
                  stringr::str_c("mean_", getPrecipVars()),
                  stringr::str_c("dev_sd_", getPrecipVars()),
                  stringr::str_c("dev_", getPrecipVars()),
                  getTempVars(),  stringr::str_c("hlq_", getTempVars(), "TRUE"),
                  stringr::str_c("mean_", getTempVars()),
                  stringr::str_c("dev_", getTempVars()),
                  stringr::str_c("dev_sd_", getTempVars()),
                  stringr::str_c("time_to_last_outbreak", levels(surv_dat_all$time_to_last_outbreak)),
                  stringr::str_c("enso_lag", 0:12))
  all_levels
}

get_covar_labels <- function(surv_dat_all) {

  all_labels <- c(month.abb,
                  getPrecipVars(), stringr::str_c("hlq_", getPrecipVars(), "TRUE"),
                  stringr::str_c("mean_", getPrecipVars()),
                  stringr::str_c("dev_sd_", getPrecipVars()),
                  stringr::str_c("dev_", getPrecipVars()),
                  getTempVars(),  stringr::str_c("hlq_", getTempVars(), "TRUE"),
                  stringr::str_c("mean_", getTempVars()),
                  stringr::str_c("dev_", getTempVars()),
                  stringr::str_c("dev_sd_", getTempVars()),
                  stringr::str_c("time_to_last_outbreak", levels(surv_dat_all$time_to_last_outbreak)),
                  stringr::str_c("ENSO (lag ", 0:12, "m)")) %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace("dev", "anomaly")  %>%
    stringr::str_replace("precip", "precip \n(100 mm)") %>%
    stringr::str_replace("mean temp daymean", "mean temp\n(10 degree)") %>%
    stringr::str_replace("daymean", "\n(1 degree)") %>%
    stringr::str_remove("time to last outbreak")
  
  all_labels
}

plot_surv_dat <- function(lp, surv_dat_all, all_ts_dat) {
  # Get gid
  gid <- surv_dat_all$gid[surv_dat_all$location_period_id == lp][1]
  
  # Get data
  this_surv_dat <- surv_dat_all %>% 
    filter(location_period_id == lp) 
  
  this_ts_dat <- all_ts_dat %>% 
    filter(location_period_id == lp) 
  
  this_dat <- getStatsGid(gids = gid) %>% 
    filter(month >= min(this_surv_dat$t_start),
           month <= max(this_surv_dat$t_end))
  
  p <- this_surv_dat %>% 
    ggplot(aes(x = t_start, y = as.numeric(outbreak))) +
    geom_linerange(data = this_ts_dat,
                   inherit.aes = FALSE,
                   aes(xmin = TL, xmax = TR, y = -0.02), col = "gray", lwd = 3) +
    geom_step(direction = "vh") +
    geom_rect(data = this_dat,
              inherit.aes = FALSE,
              aes(fill = precip, ymax = -0.1, ymin = -0.15, 
                  xmin = month, xmax = month %m+% months(1))) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    ggnewscale::new_scale_fill() +
    geom_rect(data = this_dat,
              inherit.aes = FALSE,
              aes(fill = temp_daymean, ymax = -0.16, ymin = -0.21, 
                  xmin = month, xmax = month %m+% months(1))) +
    scale_fill_gradient(low = "yellow", high = "darkred") +
    theme_bw() +
    labs(title = stringr::str_c(gid, lp, sep = " - "), 
         y = "outbreak")
  
  p
}

save_plot_surv_dat <- function(lp, surv_dat_all, all_ts_dat) {
  
  p <- plot_surv_dat(lp, surv_dat_all, all_ts_dat)
  gid <- surv_dat_all$gid[surv_dat_all$location_period_id == lp][1]
  
  cat("Done plot", lp, "\n")
  ggsave(p, filename = stringr::str_c("figures/surv_dat_plots/survdat_plot_", gid, "_", lp, ".png"),
         width = 9, height = 5)
}


getVarTypeDic <- function() {
  c(
    "max_cum7" = "precip",
    "max_precip" = "precip",
    "precip" = "precip",
    "mean_cum7" = "precip",
    "temp_daymax" = "temp",
    "temp_daymean" = "temp",
    "temp_daymin" = "temp"
  )
}

plot_clim_ext_data <- function(gid, 
                               baseline = "1990_2019",
                               fill_col = "hilo_sd") {
  
  ext_data <- getExtremesGid(gid, baseline = baseline) %>% 
    rename(var = name) %>% 
    mutate(var_type = getVarTypeDic()[var]) %>% 
    filter(month > "2009-01-01", month < "2024-01-01")
  
  ggplot(ext_data, aes(fill = !!rlang::sym(fill_col))) +
    geom_tile(aes(x = month, y = var)) +
    facet_grid(var_type ~ ., scales = "free", space = "free") +
    scale_fill_manual(values = c("red", "blue", "gray")) +
    theme_bw()
  
}


plot_clim_data <- function(gid, 
                           baseline = "1990_2019") {
  
  clim_stats <- getStatsGid(gid)
  clim_extrms <- getExtremesGid(gid = gid, baseline = baseline) 
  
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
                values_from = c("mean", "dev", "dev_sd")) %>% 
    filter(month > "2009-01-01", month < "2024-01-01")
  
  clim_dev %>%
    select(gid, month, contains("dev_sd")) %>%
    tidyr::pivot_longer(cols = contains("dev"),
                 names_to = "var") %>%
    mutate(var_type = getVarTypeDic()[stringr::str_remove(var, "dev_sd_")]) %>%
    ggplot(aes(fill = value)) +
    geom_tile(aes(x = month, y = var)) +
    facet_grid(var_type ~ ., scales = "free", space = "free") +
    theme_bw() +
    scale_fill_gradient2(limits = c(-3, 3))

}

# ---------------------------------------------------------------------------
# Figures helpers
# ---------------------------------------------------------------------------

#' Decorate a size-track re_stats tier with covar_type / var_type / var_type2.
#'
#' Applied to each list element returned by computeNestedSizeReStatsBrms()
#' (keys $region, $country, $adm1). Not needed for occurrence re_stats — those
#' columns are already produced by computeNestedSurvReStatsBrms().
decorate_covar <- function(df) {
  df %>%
    dplyr::mutate(
      covar_type  = vapply(var, getCovarType, character(1)),
      covar       = var,
      var_type    = dplyr::case_when(
        stringr::str_detect(var, "dev_")  ~ "anomaly",
        stringr::str_detect(var, "mean_") ~ "climatology",
        TRUE ~ stringr::str_replace_all(var, "_", " ")
      ),
      var_type2   = dplyr::if_else(var_type %in% c("anomaly", "climatology") |
                                   covar_type == "enso teleconnection",
                                   "climatic", "non-climatic"),
      AFRO_region = stringr::str_replace(AFRO_region, "_", " ")
    )
}
