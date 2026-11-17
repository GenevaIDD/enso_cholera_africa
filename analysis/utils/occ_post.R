# occ_post.R — postprocessing of brms GLMM fits: nested random-effect
#   summaries for the figures, and backend-agnostic term names
#
# Provides:
#   compute_nested_surv_re_stats_brms()   nested RE summaries of a brms fit
#   adm1_lineage_levels()   parent region / country levels of an ADM1 level
#   canonicalize_term()     fixed-effect term name -> shared vocabulary
#
# Details:
#   Used by the 04 workers, 05_fit_postprocess*.R and the 06 figures.

#' adm1_lineage_levels
#'
#' Given a single brms ADM1 RE level name from the
#' `AFRO_region:country:ADM1` group (encoded "<region>_<country>_<gid>",
#' e.g. "Central Africa_CMR_CMR_1"), recover the parent region- and
#' country-group level names so the nested slope deviations can be summed
#' (FE + r_region + r_country + r_adm1).  Mirrors the level-parsing in
#' compute_nested_surv_re_stats_brms() so the post-processing adapter and the
#' per-ADM1 figure agree on the encoding.
#'
#' @param adm1_level character scalar — ADM1 group level name
#' @return list(region, country, adm1) of level-name strings
adm1_lineage_levels <- function(adm1_level) {
  adm1_unit <- stringr::str_extract(adm1_level, "[A-Z]{3}_[0-9]+$")
  country_level <- stringr::str_remove(adm1_level, paste0("_", adm1_unit, "$"))
  region_level <- stringr::str_remove(country_level, "_[A-Z]{3}$")
  list(region = region_level, country = country_level, adm1 = adm1_level)
}

# Poisson / brms postprocessing ------------------------------------------
#
# Adapter: a brms Poisson-GLMM fit (04_fit_*_brms*.R) -> the nested list of RE
# summaries the figures consume. brms names: fixed effects `b_<term>`,
# deviations `r_<group>[<level>,<term>]`; nested `(x | A/B/C)` expands to groups
# A, A:B, A:B:C. Groups are rooted at country (grouping mode) or AFRO_region.
# Output: list(re_stats_region, re_stats_country, re_stats_adm1) with columns
# variable, AFRO_region (+country +adm1_unit), var, var_type, var_type0,
# var_type1, var_type2, mean, q025, q975, covar_type, covar.

#' Parse a brms r_* parameter name into group, level, and term.
#'
#' @param par character vector of parameter names.
#' @return tibble(par, group, level, term).
.parse_brms_r_par <- function(par) {
  tibble::tibble(par = par) %>%
    dplyr::mutate(
      group = stringr::str_extract(par, "(?<=^r_)[^\\[]+"),
      inside = stringr::str_extract(par, "(?<=\\[).*(?=\\])"),
      level = stringr::str_extract(inside, "^[^,]+"),
      term = stringr::str_extract(inside, "(?<=,).+$")
    ) %>%
    dplyr::select(-inside)
}

#' Compute the var_type / covar metadata used by the figures.
.attach_var_type_meta <- function(df, formula_vars, surv_dat_all) {
  df %>%
    dplyr::mutate(
      var = .data$term,
      var_type = stringr::str_extract(
        .data$var,
        stringr::str_c(formula_vars, collapse = "|")
      ),
      var_type0 = dplyr::case_when(
        stringr::str_detect(
          var_type,
          paste(get_precip_vars(), collapse = "|")
        ) ~ "precipitation",
        stringr::str_detect(
          var_type,
          paste(get_temp_vars(), collapse = "|")
        ) ~ "temperature",
        TRUE ~ stringr::str_replace_all(var_type, "_", " ")
      ),
      var_type1 = dplyr::case_when(
        stringr::str_detect(var_type, "dev") ~ "anomaly",
        stringr::str_detect(var_type, "mean") ~ "climatology",
        stringr::str_detect(var_type, "^ttlo_") ~ "time to last outbreak",
        TRUE ~ stringr::str_replace_all(var_type, "_", " ")
      ),
      var_type2 = dplyr::if_else(
        var_type1 %in% c("anomaly", "climatology", "enso"),
        "climatic",
        "non-climatic"
      )
    )
}

#' Nested random-effect summaries of a brms occurrence fit.
#'
#' @param brms_fit a fitted `brmsfit` from analysis/04b_run_poisson_models.R.
#' @param surv_dat_all the same data frame passed to the fit (used for covar
#'   level/label metadata).
#' @param formula_vars character vector of fixed-effect term names that should
#'   be tagged with var_type metadata (defaults to those in the cholera model).
#' @return list of three tibbles (region, country, ADM1 RE summaries).
compute_nested_surv_re_stats_brms <- function(
  brms_fit,
  surv_dat_all,
  formula_vars = NULL,
  drop_region = FALSE
) {
  # All draws in long form
  draws_df <- posterior::as_draws_df(brms_fit) %>%
    tibble::as_tibble()

  # Fixed-effect draws keyed by term
  b_cols <- stringr::str_subset(names(draws_df), "^b_")
  b_long <- draws_df %>%
    dplyr::select(dplyr::all_of(c(".draw", b_cols))) %>%
    tidyr::pivot_longer(
      dplyr::all_of(b_cols),
      names_to = "par",
      values_to = "b"
    ) %>%
    dplyr::mutate(term = stringr::str_remove(par, "^b_"))

  # Group-level deviation draws
  r_cols <- stringr::str_subset(names(draws_df), "^r_")
  if (length(r_cols) == 0) {
    stop("No r_* draws found; was the model fit with random effects?")
  }
  r_long <- draws_df %>%
    dplyr::select(dplyr::all_of(c(".draw", r_cols))) %>%
    tidyr::pivot_longer(
      dplyr::all_of(r_cols),
      names_to = "par",
      values_to = "r"
    ) %>%
    dplyr::left_join(.parse_brms_r_par(r_cols), by = "par")

  # Derive formula_vars from the actual RE terms present in the fit.
  # Normalise enso_lag{N} → "enso" so _attachVarTypeMeta can still classify
  # ENSO terms as climatic (var_type2 == "enso" path in the case_when).
  if (is.null(formula_vars)) {
    formula_vars <- sub("^enso_lag\\d+$", "enso", unique(r_long$term))
  }

  # Split by nesting depth via group name
  if (drop_region) {
    r_region <- NULL
    r_country <- r_long %>% dplyr::filter(group == "country")
    r_adm1 <- r_long %>% dplyr::filter(group == "country:ADM1")
    r_adm2 <- r_long %>% dplyr::filter(group == "country:ADM1:gid")
  } else {
    r_region <- r_long %>% dplyr::filter(group == "AFRO_region")
    r_country <- r_long %>% dplyr::filter(group == "AFRO_region:country")
    r_adm1 <- r_long %>% dplyr::filter(group == "AFRO_region:country:ADM1")
    r_adm2 <- r_long %>% dplyr::filter(group == "AFRO_region:country:ADM1:gid")
  }

  # Helper to summarize exp(value) draws into mean / 95% CrI
  summarize_er <- function(df) {
    df %>%
      dplyr::summarise(
        mean = mean(exp(value)),
        q025 = stats::quantile(exp(value), 0.025),
        q975 = stats::quantile(exp(value), 0.975),
        .groups = "drop"
      )
  }

  # Region-level: b + r_region (NULL under drop_region)
  re_stats_region <- if (!drop_region) {
    r_region %>%
      dplyr::inner_join(
        b_long %>% dplyr::select(.draw, term, b),
        by = c(".draw", "term")
      ) %>%
      dplyr::mutate(value = b + r, AFRO_region = level) %>%
      dplyr::group_by(AFRO_region, term) %>%
      dplyr::group_modify(
        ~ summarize_er(
          .x %>%
            dplyr::rename(value = value)
        )
      ) %>%
      dplyr::ungroup() %>%
      .attach_var_type_meta(formula_vars, surv_dat_all) %>%
      dplyr::mutate(
        country = stringr::str_c("AFRO_region:", AFRO_region),
        variable = stringr::str_c(
          "r_AFRO_region[",
          AFRO_region,
          ",",
          term,
          "]"
        ),
        covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
        covar = factor(
          var,
          levels = get_covar_levels(surv_dat_all),
          labels = get_covar_labels(surv_dat_all)
        )
      )
  } else {
    NULL
  }

  # Country-level
  # Under non-drop: group = "AFRO_region:country", level = "<region>_<ISO3>".
  # Under drop_region: group = "country", level = ISO3 code.
  if (drop_region) {
    country_levels <- r_country %>%
      dplyr::distinct(level) %>%
      dplyr::mutate(country = level, AFRO_region = NA_character_)
  } else {
    country_levels <- r_country %>%
      dplyr::distinct(level) %>%
      dplyr::mutate(
        country = stringr::str_extract(level, "[A-Z]{3}$"),
        AFRO_region = stringr::str_remove(level, "_[A-Z]{3}$")
      )
  }

  r_country2 <- r_country %>%
    dplyr::left_join(country_levels, by = "level")

  re_stats_country <- if (drop_region) {
    # b + r_country (no region level)
    r_country2 %>%
      dplyr::inner_join(
        b_long %>% dplyr::select(.draw, term, b),
        by = c(".draw", "term")
      ) %>%
      dplyr::mutate(value = b + r) %>%
      dplyr::group_by(country, term) %>%
      dplyr::group_modify(
        ~ summarize_er(.x %>% dplyr::rename(value = value))
      ) %>%
      dplyr::ungroup() %>%
      .attach_var_type_meta(formula_vars, surv_dat_all) %>%
      dplyr::mutate(
        AFRO_region = NA_character_,
        variable = stringr::str_c("r_country[", country, ",", term, "]"),
        covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
        covar = factor(
          var,
          levels = get_covar_levels(surv_dat_all),
          labels = get_covar_labels(surv_dat_all)
        )
      )
  } else {
    # b + r_region + r_country
    r_country2 %>%
      dplyr::inner_join(
        r_region %>%
          dplyr::select(.draw, term, AFRO_region = level, r_region = r),
        by = c(".draw", "term", "AFRO_region")
      ) %>%
      dplyr::inner_join(
        b_long %>% dplyr::select(.draw, term, b),
        by = c(".draw", "term")
      ) %>%
      dplyr::mutate(value = b + r_region + r) %>%
      dplyr::group_by(AFRO_region, country, term) %>%
      dplyr::group_modify(
        ~ summarize_er(.x %>% dplyr::rename(value = value))
      ) %>%
      dplyr::ungroup() %>%
      .attach_var_type_meta(formula_vars, surv_dat_all) %>%
      dplyr::mutate(
        variable = stringr::str_c(
          "r_AFRO_region:country[",
          AFRO_region,
          "_",
          country,
          ",",
          term,
          "]"
        ),
        covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
        covar = factor(
          var,
          levels = get_covar_levels(surv_dat_all),
          labels = get_covar_labels(surv_dat_all)
        )
      )
  }

  # ADM1-level
  # Level encoding under non-drop: "<region>_<country>_<adm1_unit>".
  # Level encoding under drop_region: "<country>_<adm1_unit>".
  # adm1_unit = trailing "[A-Z]{3}_[0-9]+" — same extraction in both cases.
  adm1_levels <- r_adm1 %>%
    dplyr::distinct(level) %>%
    dplyr::mutate(
      adm1_unit = stringr::str_extract(level, "[A-Z]{3}_[0-9]+$"),
      country = stringr::str_extract(adm1_unit, "^[A-Z]{3}"),
      AFRO_region = if (drop_region) {
        NA_character_
      } else {
        # Strip trailing "_<adm1_unit>" then "_<country>" to recover region.
        # e.g. "Central.Africa_CMR_CMR_1" -> "Central.Africa_CMR" -> "Central.Africa"
        stringr::str_remove(
          stringr::str_remove(level, paste0("_", adm1_unit, "$")),
          "_[A-Z]{3}$"
        )
      }
    )

  r_adm1_2 <- r_adm1 %>%
    dplyr::left_join(adm1_levels, by = "level")

  re_stats_adm1 <- if (drop_region) {
    # b + r_country + r_adm1 (no region level)
    r_adm1_2 %>%
      dplyr::inner_join(
        r_country2 %>%
          dplyr::select(.draw, term, country, r_country = r),
        by = c(".draw", "term", "country")
      ) %>%
      dplyr::inner_join(
        b_long %>% dplyr::select(.draw, term, b),
        by = c(".draw", "term")
      ) %>%
      dplyr::mutate(value = b + r_country + r) %>%
      dplyr::group_by(country, adm1_unit, term) %>%
      dplyr::group_modify(
        ~ summarize_er(.x %>% dplyr::rename(value = value))
      ) %>%
      dplyr::ungroup() %>%
      .attach_var_type_meta(formula_vars, surv_dat_all) %>%
      dplyr::mutate(
        AFRO_region = NA_character_,
        variable = stringr::str_c(
          "r_country:ADM1[",
          country,
          "_",
          adm1_unit,
          ",",
          term,
          "]"
        ),
        covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
        covar = factor(
          var,
          levels = get_covar_levels(surv_dat_all),
          labels = get_covar_labels(surv_dat_all)
        )
      )
  } else {
    # b + r_region + r_country + r_adm1
    r_adm1_2 %>%
      dplyr::inner_join(
        r_country2 %>%
          dplyr::select(.draw, term, AFRO_region, country, r_country = r),
        by = c(".draw", "term", "AFRO_region", "country")
      ) %>%
      dplyr::inner_join(
        r_region %>%
          dplyr::select(.draw, term, AFRO_region = level, r_region = r),
        by = c(".draw", "term", "AFRO_region")
      ) %>%
      dplyr::inner_join(
        b_long %>% dplyr::select(.draw, term, b),
        by = c(".draw", "term")
      ) %>%
      dplyr::mutate(value = b + r_region + r_country + r) %>%
      dplyr::group_by(AFRO_region, country, adm1_unit, term) %>%
      dplyr::group_modify(
        ~ summarize_er(.x %>% dplyr::rename(value = value))
      ) %>%
      dplyr::ungroup() %>%
      .attach_var_type_meta(formula_vars, surv_dat_all) %>%
      dplyr::mutate(
        variable = stringr::str_c(
          "r_AFRO_region:country:ADM1[",
          AFRO_region,
          "_",
          adm1_unit,
          ",",
          term,
          "]"
        ),
        covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
        covar = factor(
          var,
          levels = get_covar_levels(surv_dat_all),
          labels = get_covar_labels(surv_dat_all)
        )
      )
  }

  # ADM2-level
  # Level encoding under non-drop: "<region>_<country>_<adm1_gid>_<adm2_gid>"
  # Level encoding under drop_region: "<country>_<adm1_gid>_<adm2_gid>"
  # GADM ADM2 gids have 4 components: "<country>_<adm1n>_<adm2n>_<version>",
  # except Ghana's, where GADM 4.1 omits the first separator ("GHA14_13_2");
  # adm1_unit restores it as gid.R::gid_to_adm1() does ("GHA_14", the ADM1 key).
  re_stats_adm2 <- if (nrow(r_adm2) > 0) {
    adm2_levels <- r_adm2 %>%
      dplyr::distinct(level) %>%
      dplyr::mutate(
        gid = stringr::str_extract(level, "[A-Z]{3}_?[0-9]+_[0-9]+_[0-9]+$"),
        adm1_unit = stringr::str_extract(
          sub("^([A-Z]{3})([0-9])", "\\1_\\2", gid),
          "^[A-Z]{3}_[0-9]+"
        ),
        country = stringr::str_extract(adm1_unit, "^[A-Z]{3}"),
        AFRO_region = if (drop_region) {
          NA_character_
        } else {
          # Strip "_<gid>", "_<adm1_unit>", "_<country>" to recover region.
          stringr::str_remove(
            stringr::str_remove(
              stringr::str_remove(level, paste0("_", gid, "$")),
              paste0("_", adm1_unit, "$")
            ),
            "_[A-Z]{3}$"
          )
        }
      )

    r_adm2_2 <- r_adm2 %>%
      dplyr::left_join(adm2_levels, by = "level")

    if (drop_region) {
      # b + r_country + r_adm1 + r_adm2 (no region level)
      r_adm2_2 %>%
        dplyr::inner_join(
          r_adm1_2 %>%
            dplyr::select(.draw, term, country, adm1_unit, r_adm1 = r),
          by = c(".draw", "term", "country", "adm1_unit")
        ) %>%
        dplyr::inner_join(
          r_country2 %>%
            dplyr::select(.draw, term, country, r_country = r),
          by = c(".draw", "term", "country")
        ) %>%
        dplyr::inner_join(
          b_long %>% dplyr::select(.draw, term, b),
          by = c(".draw", "term")
        ) %>%
        dplyr::mutate(value = b + r_country + r_adm1 + r) %>%
        dplyr::group_by(country, adm1_unit, gid, term) %>%
        dplyr::group_modify(
          ~ summarize_er(.x %>% dplyr::rename(value = value))
        ) %>%
        dplyr::ungroup() %>%
        .attach_var_type_meta(formula_vars, surv_dat_all) %>%
        dplyr::mutate(
          AFRO_region = NA_character_,
          variable = stringr::str_c(
            "r_country:ADM1:gid[",
            country,
            "_",
            gid,
            ",",
            term,
            "]"
          ),
          covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
          covar = factor(
            var,
            levels = get_covar_levels(surv_dat_all),
            labels = get_covar_labels(surv_dat_all)
          )
        )
    } else {
      # b + r_region + r_country + r_adm1 + r_adm2
      r_adm2_2 %>%
        dplyr::inner_join(
          r_adm1_2 %>%
            dplyr::select(
              .draw,
              term,
              AFRO_region,
              country,
              adm1_unit,
              r_adm1 = r
            ),
          by = c(".draw", "term", "AFRO_region", "country", "adm1_unit")
        ) %>%
        dplyr::inner_join(
          r_country2 %>%
            dplyr::select(.draw, term, AFRO_region, country, r_country = r),
          by = c(".draw", "term", "AFRO_region", "country")
        ) %>%
        dplyr::inner_join(
          r_region %>%
            dplyr::select(.draw, term, AFRO_region = level, r_region = r),
          by = c(".draw", "term", "AFRO_region")
        ) %>%
        dplyr::inner_join(
          b_long %>% dplyr::select(.draw, term, b),
          by = c(".draw", "term")
        ) %>%
        dplyr::mutate(value = b + r_region + r_country + r_adm1 + r) %>%
        dplyr::group_by(AFRO_region, country, adm1_unit, gid, term) %>%
        dplyr::group_modify(
          ~ summarize_er(.x %>% dplyr::rename(value = value))
        ) %>%
        dplyr::ungroup() %>%
        .attach_var_type_meta(formula_vars, surv_dat_all) %>%
        dplyr::mutate(
          variable = stringr::str_c(
            "r_AFRO_region:country:ADM1:gid[",
            AFRO_region,
            "_",
            gid,
            ",",
            term,
            "]"
          ),
          covar_type = purrr::map_chr(var, ~ get_covar_type(.)),
          covar = factor(
            var,
            levels = get_covar_levels(surv_dat_all),
            labels = get_covar_labels(surv_dat_all)
          )
        )
    }
  } else {
    NULL
  }

  list(
    re_stats_region = re_stats_region,
    re_stats_country = re_stats_country,
    re_stats_adm1 = re_stats_adm1,
    re_stats_adm2 = re_stats_adm2
  )
}


# Term-name canonicalization ----------------------------------------------

#' Canonicalize a fixed-effect term name to a backend-agnostic string.
#'
#' Maps brms encodings (factor() prefix, ttlo level suffixes, TRUE suffix on
#' logicals) to a shared vocabulary.
#' Unrecognised terms are returned unchanged.
#'
#' @param x  Character vector of term names.
#' @return   Character vector, same length.
canonicalize_term <- function(x) {
  x <- as.character(x)

  # months: factor(months)2 -> months2, factormonths2 -> months2
  x <- sub("^factor\\(months\\)", "months", x)
  x <- sub("^factormonths", "months", x)

  # time_to_last_outbreak -> ttlo_LO_HI
  # brms: time_to_last_outbreak12, time_to_last_outbreak6Inf
  x <- sub("^time_to_last_outbreak([0-9]+)(Inf|[0-9]+)", "ttlo_\\1_\\2", x)
  # The "no recorded prior onset" level (occurrence.ttlo_no_prior_level):
  # brms time_to_last_outbreaknone.
  x <- sub("^time_to_last_outbreaknone$", "ttlo_none", x)

  # logical TRUE suffix (brms): ongoing_outbreak_lag1TRUE -> ongoing_outbreak_lag1
  x <- sub("TRUE$", "", x)

  x
}
