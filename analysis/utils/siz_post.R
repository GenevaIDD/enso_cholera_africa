
#' Posterior-predictive scatter as a tibble (for plotting downstream).
#'
#' @param response name of the observed-response column in fit$data.
#'   Default "suspected_cases" (size track); the monthly track passes
#'   "cases_month".
makeSizePpCheckTbl <- function(fit, response = "suspected_cases") {
  pred <- predict(fit) %>% tibble::as_tibble()
  pred$obs <- fit$data[[response]]
  pred
}

#' Fixed-effect coefficient summary on the response scale.
summariseSizeCoefsBrms <- function(fit) {
  posterior::summarise_draws(
    brms::as_draws(fit, variable = "^b_", regex = TRUE)
  )
}

#' computeNestedSizeReStatsBrms
#'
#' Port of the inline RE-unpack logic from legacy 07 (lines 136-276).
#' Returns list(region, country, adm1) — same schema as the occurrence-side
#' adapter. Effects are on the response (multiplicative) scale.
computeNestedSizeReStatsBrms <- function(fit, opt = NULL, drop_region = FALSE) {
  all_re <- brms::ranef(fit, summary = FALSE)

  # Apply exactly n_seps sub("_",":") calls per group, where n_seps is the
  # number of colons in the group key name. This correctly converts the
  # brms underscore-joined level labels to colon-separated component labels
  # regardless of nesting depth (handles both AFRO_region-rooted and
  # drop_region=TRUE country-rooted group structures).
  all_re_df <- purrr::map_df(seq_along(all_re), function(x) {
    y        <- all_re[[x]]
    key_name <- names(all_re)[[x]]
    n_seps   <- stringr::str_count(key_name, ":")
    var_names <- attr(y, "dimnames")[[3]]
    n_draws   <- dim(y)[1]
    purrr::map_df(seq_len(dim(y)[2]), function(z) {
      grp_label <- attr(y, "dimnames")[[2]][z]
      for (i in seq_len(n_seps)) {
        grp_label <- sub("_", ":", grp_label, fixed = TRUE)
      }
      # Robust slice: force (draws x vars) matrix even when vars == 1.
      slice <- y[, z, , drop = FALSE]            # dim: (draws, 1, vars)
      mat   <- matrix(slice, nrow = n_draws,
                      ncol = length(var_names),
                      dimnames = list(NULL, var_names))
      tibble::as_tibble(mat) %>%
        dplyr::mutate(.draw = seq_len(n_draws)) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(var_names),
                            names_to = "variable") %>%
        dplyr::mutate(variable = paste(grp_label, variable, sep = "-"),
                      .key = key_name)
    })
  })
  rm(all_re)  # raw draw array no longer needed

  # Variable layout produced above:
  #   non-drop region-level:  "<region>-<var>"           (e.g., "Eastern Africa-Intercept")
  #   non-drop country-level: "<region>:<ISO3>-<var>"    (e.g., "Eastern Africa:ETH-Intercept")
  #   non-drop adm1-level:    "<region>:<ISO3>:<ISO3>_N-<var>"
  #   drop_region country-level: "<ISO3>-<var>"           (e.g., "ETH-Intercept")
  #   drop_region adm1-level:    "<ISO3>:<ISO3>_N-<var>"  (e.g., "ETH:ETH_1-Intercept")
  if (drop_region) {
    re_draws <- all_re_df %>%
      dplyr::mutate(
        grp_part = stringr::str_extract(variable, ".*(?=-)"),
        var      = sub("^.*-", "", variable),
        .parts   = stringr::str_split(grp_part, ":"),
        n_parts  = lengths(.parts),
        # Under drop_region: 1-part = country tier, 2-part = adm1 tier
        AFRO_region = NA_character_,
        country     = purrr::map_chr(.parts, ~ .x[1]),
        adm1_unit   = purrr::map2_chr(.parts, n_parts,
                                      ~ if (.y >= 2) .x[2] else NA_character_)
      ) %>%
      dplyr::select(-grp_part, -.parts, -n_parts, -.key)
  } else {
    re_draws <- all_re_df %>%
      dplyr::mutate(
        grp_part = stringr::str_extract(variable, ".*(?=-)"),
        var      = sub("^.*-", "", variable),
        .parts   = stringr::str_split(grp_part, ":"),
        n_parts  = lengths(.parts),
        AFRO_region = purrr::map_chr(.parts, ~ .x[1]),
        country     = purrr::map2_chr(.parts, n_parts,
                                      ~ if (.y >= 2) .x[2] else NA_character_),
        adm1_unit   = purrr::map2_chr(.parts, n_parts,
                                      ~ if (.y >= 3) .x[3] else NA_character_)
      ) %>%
      dplyr::select(-grp_part, -.parts, -n_parts, -.key)
  }

  if (drop_region) {
    re_region  <- NULL
    re_country <- re_draws %>% dplyr::filter(is.na(adm1_unit))
    re_adm1    <- re_draws %>% dplyr::filter(!is.na(adm1_unit))
  } else {
    re_region  <- re_draws %>% dplyr::filter(is.na(country))
    re_country <- re_draws %>% dplyr::filter(!is.na(country), is.na(adm1_unit))
    re_adm1    <- re_draws %>% dplyr::filter(!is.na(adm1_unit))
  }
  rm(all_re_df, re_draws)  # long pivot no longer needed; tier subsets are sufficient

  mean_draws <- brms::as_draws_df(fit, variable = "^b_", regex = TRUE) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = !dplyr::contains("."),
                        names_to = "var", values_to = "mean") %>%
    dplyr::mutate(var = sub("^b_", "", var)) %>%
    dplyr::select(.draw, var, mean)

  summarise_grp <- function(df, ...) {
    df %>%
      dplyr::group_by(...) %>%
      dplyr::summarise(
        mean = mean(exp(value), na.rm = TRUE),
        q025 = unname(stats::quantile(exp(value), 0.025, na.rm = TRUE)),
        q975 = unname(stats::quantile(exp(value), 0.975, na.rm = TRUE)),
        .groups = "drop"
      )
  }

  region_tbl <- if (!drop_region) {
    re_region %>%
      dplyr::inner_join(mean_draws, by = c(".draw", "var")) %>%
      dplyr::mutate(value = value + mean) %>%
      summarise_grp(variable, AFRO_region, var)
  } else {
    NULL
  }

  country_tbl <- if (drop_region) {
    # b + r_country (no region level). AFRO_region is kept in the grouping
    # columns (constant NA_character_, set above) rather than dropped --
    # mirrors occ_post.R's computeNestedSurvReStatsBrms(), which always keeps
    # an AFRO_region column (NA under drop_region) so downstream figure code
    # (decorate_covar() et al.) doesn't need a drop_region-specific branch.
    re_country %>%
      dplyr::inner_join(mean_draws, by = c(".draw", "var")) %>%
      dplyr::mutate(value = value + mean) %>%
      summarise_grp(variable, AFRO_region, country, var)
  } else {
    re_country %>%
      dplyr::left_join(
        re_region %>% dplyr::select(.draw, var, region_mean = value),
        by = c(".draw", "var")
      ) %>%
      dplyr::inner_join(mean_draws, by = c(".draw", "var")) %>%
      dplyr::mutate(value = value + tidyr::replace_na(region_mean, 0) + mean) %>%
      summarise_grp(variable, AFRO_region, country, var)
  }

  adm1_tbl <- if (drop_region) {
    # b + r_country + r_adm1 (no region level); join on country for correct
    # pairing. AFRO_region kept in the grouping columns (constant NA) for the
    # same reason as country_tbl above.
    re_adm1 %>%
      dplyr::left_join(
        re_country %>% dplyr::select(.draw, var, country, country_mean = value),
        by = c(".draw", "var", "country")
      ) %>%
      dplyr::inner_join(mean_draws, by = c(".draw", "var")) %>%
      dplyr::mutate(value = value +
                            tidyr::replace_na(country_mean, 0) + mean) %>%
      summarise_grp(variable, AFRO_region, country, adm1_unit, var)
  } else {
    re_adm1 %>%
      dplyr::left_join(
        re_country %>% dplyr::select(.draw, var, country_mean = value),
        by = c(".draw", "var")
      ) %>%
      dplyr::left_join(
        re_region %>% dplyr::select(.draw, var, region_mean = value),
        by = c(".draw", "var")
      ) %>%
      dplyr::inner_join(mean_draws, by = c(".draw", "var")) %>%
      dplyr::mutate(value = value +
                            tidyr::replace_na(country_mean, 0) +
                            tidyr::replace_na(region_mean, 0) + mean) %>%
      summarise_grp(variable, AFRO_region, country, adm1_unit, var)
  }
  rm(re_region, re_country, re_adm1, mean_draws); gc()

  # --- phi distributional RE (only when phi_re != "constant") ---
  phi_re_stats <- NULL
  phi_re <- if (!is.null(opt)) opt$size$phi_re %||% "constant" else "constant"
  if (phi_re != "constant") {
    phi_ranef <- brms::ranef(fit, summary = FALSE, dpar = "phi")
    # brms names nested groups with the original column names from data.
    # For (1 | AFRO_region/country/ADM1) the keys are:
    #   "AFRO_region", "AFRO_region:country", "AFRO_region:country:ADM1"
    phi_group_keys <- switch(phi_re,
      region              = "AFRO_region",
      region_country      = c("AFRO_region", "AFRO_region:country"),
      region_country_adm1 = c("AFRO_region", "AFRO_region:country", "AFRO_region:country:ADM1")
    )
    summarise_phi_ranef <- function(arr) {
      # arr: [draws x groups x 1]
      n_draws <- dim(arr)[1]
      n_groups <- dim(arr)[2]
      grp_names <- dimnames(arr)[[2]]
      mat <- matrix(arr[, , 1, drop = TRUE], nrow = n_draws,
                    dimnames = list(NULL, grp_names))
      tibble::tibble(
        group    = grp_names,
        mean_phi = colMeans(exp(mat)),
        q025_phi = apply(exp(mat), 2, stats::quantile, 0.025),
        q975_phi = apply(exp(mat), 2, stats::quantile, 0.975)
      )
    }
    phi_re_stats <- lapply(phi_group_keys, function(k) {
      if (!is.null(phi_ranef[[k]])) summarise_phi_ranef(phi_ranef[[k]]) else NULL
    })
    names(phi_re_stats) <- phi_group_keys
  }

  list(region = region_tbl, country = country_tbl, adm1 = adm1_tbl,
       phi_re_stats = phi_re_stats)
}
