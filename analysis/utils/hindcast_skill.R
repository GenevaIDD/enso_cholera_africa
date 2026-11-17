# hindcast_skill.R — forecast-skill scoring of the spatPomp hindcast
#   experiment (occurrence, cases, monthly and season case skill).
#
# Provides:
#   hindcast_skill_bundle()     one pass over the runs, references, roll-ups
#   hindcast_report_tables()    phase-11 numbers and notebook summary tables
#   hindcast_run_index(), hindcast_load_run(), hindcast_obs_panel()
#   score_occurrence_model(), score_occurrence_references()   occurrence
#   score_cases_run(), score_case_references()                cases
#   hindcast_lead_class(), hindcast_boot_spec()   lead classes, bootstrap CIs
#   crps_discrete(), crps_binom(), brier(), interval_score(), ...   primitives
#
# Details:
#   09 derives every hindcast artefact once (results_derive.R derive_hindcast_*)
#   and the notebook / 10 only plot. Keyed on the ENSO `group`; labelling is the
#   caller's job (results_core.R's GROUP_* registries). Callers source base.R,
#   gid.R, obs_sim_compare.R, sim_io.R, sim_metrics.R; at call time also
#   population.R, config.R, paths.R and provenance.R.

# ===========================================================================
# 1. Lead classes
# ===========================================================================

#' hindcast_lead_class
#'
#' Months since launch -> lead-class label. breaks c(7, 12) gives "0-7",
#' "8-12", "13+". Lead 0 is the launch month.
hindcast_lead_class <- function(lead, breaks = c(7L, 12L)) {
  breaks <- sort(as.integer(breaks))
  lo <- c(0L, breaks + 1L)
  labs <- c(
    sprintf("%d-%d", lo[-length(lo)], breaks),
    sprintf("%d+", lo[length(lo)])
  )
  # evaluated once per distinct lead (the case frames repeat each lead per
  # gid x replicate)
  u <- unique(lead)
  as.character(cut(
    u,
    breaks = c(-Inf, breaks, Inf),
    labels = labs,
    right = TRUE
  ))[match(lead, u)]
}

# first-of-month Date -> integer month index (year * 12 + month - 1), evaluated
# once per distinct date
.month_index <- function(d) {
  u <- unique(as.Date(d))
  lt <- as.POSIXlt(u)
  ((lt$year + 1900L) * 12L + lt$mon)[match(as.Date(d), u)]
}

#' hindcast_lead_months -- months since launch (lead 0 = the launch month). With
#' `country_launch` (ISO3 -> launch; season-window runs) each row counts from
#' its own country's launch (`iso3`).
hindcast_lead_months <- function(
  month,
  launch,
  country_launch = NULL,
  iso3 = NULL
) {
  L <- if (is.null(country_launch)) {
    as.Date(launch)
  } else {
    as.Date(unname(country_launch[iso3]))
  }
  as.integer(.month_index(month) - .month_index(L))
}

#' hindcast_lead_levels -- "all" followed by the strata, in order.
hindcast_lead_levels <- function(breaks = c(7L, 12L)) {
  breaks <- sort(as.integer(breaks))
  c("all", hindcast_lead_class(c(0L, breaks + 1L), breaks))
}

# shift first-of-month Dates by k whole years (persistence references)
.shift_years <- function(m, k = 1L) {
  as.Date(sprintf("%04d-%02d-01", lubridate::year(m) + k, lubridate::month(m)))
}

# replicate rows of a frame once per lead class it belongs to ("all" + its own)
.with_all_lead <- function(df) {
  dplyr::bind_rows(dplyr::mutate(df, lead_class = "all"), df)
}

# ===========================================================================
# 2. Run index and per-run loader
# ===========================================================================

#' hindcast_run_index
#'
#' One row per (group, launch year) of a hindcast manifest, with every path 09
#' reads. Uses for_each_hindcast_run() (no partition aggregation: that is 07c's
#' job). Stops when a run is skipped (a silently thinner skill table is the
#' failure mode), when the runs disagree on episodes.key, or -- when
#' `episode_key` is given -- when they were simulated under another key. The
#' same holds for episodes.national_zero_fill (national zero-fill occurrence panel) and
#' `national_zero_fill`.
hindcast_run_index <- function(
  man,
  episode_key = NULL,
  tag = "[hindcast-skill]",
  national_zero_fill = NULL,
  national_zero_threshold = NULL,
  occ_panel = NULL
) {
  res <- for_each_hindcast_run(
    man,
    option_list = list(),
    tag = tag,
    announce = "after",
    ensure_summary = NULL,
    fn = function(run) {
      ek <- episode_key_of(run$opt)
      if (!is.null(episode_key) && !identical(ek, episode_key)) {
        stop(
          sprintf(
            paste0(
              "%s %s/%s was simulated under episodes.key = '%s' but ",
              "the results config resolves '%s'."
            ),
            tag,
            run$group,
            run$year,
            ek,
            episode_key
          ),
          call. = FALSE
        )
      }
      nzf <- national_zero_fill_of(run$opt)
      if (
        !is.null(national_zero_fill) &&
          !identical(nzf, isTRUE(national_zero_fill))
      ) {
        stop(
          sprintf(
            paste0(
              "%s %s/%s was simulated under episodes.national_zero_fill = %s ",
              "but the results config resolves %s."
            ),
            tag,
            run$group,
            run$year,
            tolower(nzf),
            tolower(isTRUE(national_zero_fill))
          ),
          call. = FALSE
        )
      }
      nzt <- national_zero_threshold_of(run$opt)
      if (
        !is.null(national_zero_threshold) &&
          !identical(nzt, as.integer(national_zero_threshold))
      ) {
        stop(
          sprintf(
            paste0(
              "%s %s/%s was simulated under episodes.national_zero_threshold = %d ",
              "but the results config resolves %d."
            ),
            tag,
            run$group,
            run$year,
            nzt,
            as.integer(national_zero_threshold)
          ),
          call. = FALSE
        )
      }
      # The occurrence-panel spec fields (occ_panel: the results config's spec; the
      # index keeps its columns, so the hindcast_runs artefact is unchanged).
      if (!is.null(occ_panel)) {
        rs <- occ_panel_spec_of(run$opt)
        for (f in c("national_zero_rule", "at_risk_start", "adm1_overlap")) {
          if (!identical(rs[[f]], occ_panel[[f]])) {
            stop(
              sprintf(
                paste0(
                  "%s %s/%s was simulated under episodes.%s = '%s' ",
                  "but the results config resolves '%s'."
                ),
                tag,
                run$group,
                run$year,
                f,
                rs[[f]],
                occ_panel[[f]]
              ),
              call. = FALSE
            )
          }
        }
      }
      vd <- make_run_file(run$opt, kind = "spatpomp_valid_dir")
      base <- tibble::tibble(
        group = run$group,
        launch_year = as.integer(run$launch_year %||% run$year),
        launch = as.Date(run$entry$launch),
        horizon_end = as.Date(run$entry$horizon_end),
        run_tag = run$entry$run_tag %||% NA_character_,
        config = run$config,
        merged_dir = run$merged_dir,
        valid_dir = vd,
        monthly_summary_path = file.path(
          run$merged_dir,
          paste0("monthly_summary", sim_io_ext(run$out_fmt))
        ),
        metrics_path = file.path(
          vd,
          paste0("metrics", sim_io_ext(sim_metrics_fmt(run$out_fmt)))
        ),
        out_fmt = run$out_fmt,
        episode_key = ek,
        national_zero_fill = nzf,
        national_zero_threshold = nzt,
        unit = run$group,
        countries = list(NULL),
        country_horizons = list(NULL)
      )
      if (is.null(run$firm_horizon_by_country)) {
        return(base)
      }
      # v2 (registry) manifest: a unit run may hold several fit groups. Score it
      # as one row per FIT group -- its member countries, its own firm window --
      # so every downstream rollup stays keyed by the fit group. Under
      # per-country season ends a group's countries end at different months
      # (horizon cohorts): each country is read up to its own firm horizon
      # (country_horizons, .hindcast_read_part()) and the row's horizon_end is
      # the group's latest, the same for every unit row of the group. When all
      # members share one horizon this is the former single-window row.
      firm <- run$firm_horizon_by_country
      cg <- manifest_country_groups(man)
      mem <- names(firm)[cg[names(firm)] %in% run$unit_groups]
      run_members <- if (length(man$units)) {
        unlist(man$units[[run$group]]$members)
      } else {
        unlist(man$groups[[run$group]]$iso3)
      }
      mem <- intersect(mem, run_members)
      # Season-window manifest: a run scores only the countries whose season
      # starts at its launch; a fit group's season year is pooled from every run
      # holding one of its countries, so the row carries the group's season-year
      # horizon (the same for all of them) and each country's launch.
      season <- !is.null(run$scored_countries)
      if (season) {
        mem <- intersect(mem, run$scored_countries)
      }
      dplyr::bind_rows(lapply(run$unit_groups, function(g) {
        cs <- mem[cg[mem] == g]
        if (!length(cs)) {
          return(NULL)
        }
        g_all <- names(firm)[cg[names(firm)] == g]
        h_g <- if (season) {
          as.Date(
            man$groups[[g]]$years[[as.character(base$launch_year)]]$horizon_end
          )
        } else {
          max(as.Date(unname(firm[g_all])))
        }
        out <- dplyr::mutate(
          base,
          group = g,
          horizon_end = h_g,
          countries = list(cs),
          country_horizons = list(firm[cs])
        )
        if (season) {
          out$country_launch <- list(stats::setNames(
            rep(as.character(run$entry$launch), length(cs)),
            cs
          ))
          # each country's season end and floor flag (manifest countries.<iso>)
          out$country_season_end <- list(stats::setNames(
            vapply(
              cs,
              function(iso) {
                as.character(
                  man$countries[[iso]]$years[[run$year]]$season_end %||%
                    NA_character_
                )
              },
              character(1)
            ),
            cs
          ))
          out$country_floored <- list(stats::setNames(
            vapply(
              cs,
              function(iso) {
                isTRUE(man$countries[[iso]]$floored)
              },
              logical(1)
            ),
            cs
          ))
        }
        out
      }))
    }
  )
  if (res$n_skip > 0) {
    stop(
      sprintf(
        "%s %d manifest run(s) skipped (see warnings); refusing a thinner skill table.",
        tag,
        res$n_skip
      ),
      call. = FALSE
    )
  }
  if (!res$n_ok) {
    stop(tag, " no hindcast runs in the manifest.", call. = FALSE)
  }
  idx <- dplyr::bind_rows(res$results)
  if (length(unique(idx$episode_key)) > 1L) {
    stop(
      tag,
      " the manifest's runs disagree on episodes.key (",
      paste(unique(idx$episode_key), collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  if (length(unique(idx$national_zero_fill)) > 1L) {
    stop(
      tag,
      " the manifest's runs disagree on episodes.national_zero_fill.",
      call. = FALSE
    )
  }
  if (length(unique(idx$national_zero_threshold %||% 0L)) > 1L) {
    stop(
      tag,
      " the manifest's runs disagree on episodes.national_zero_threshold.",
      call. = FALSE
    )
  }
  miss <- idx$metrics_path[!file.exists(idx$metrics_path)]
  if (length(miss)) {
    stop(
      tag,
      " missing validation/metrics for ",
      length(miss),
      " run(s), e.g. ",
      miss[1],
      ". Run 07c --hindcast first.",
      call. = FALSE
    )
  }
  idx
}

#' hindcast_panel_paths
#'
#' The observed survival and monthly panels the hindcast runs must be scored
#' against: the prep artefacts of the runs' own episodes.key (what 07_hindcast.sh
#' passes 07c as OBS_PANEL), checked for a shared episode definition. Scoring a
#' gid_season run against the lp_season panel would silently change the support.
hindcast_panel_paths <- function(
  idx,
  baseline = "1990_2019",
  tag = "[hindcast-skill]",
  occ_panel = NULL
) {
  ek <- unique(idx$episode_key)
  nzf <- isTRUE(unique(idx$national_zero_fill %||% FALSE))
  nzt <- as.integer(unique(idx$national_zero_threshold %||% 0L))
  # occ_panel (checked run by run in hindcast_run_index()) supplies the item
  # 11/12 fields the index does not carry; NULL = their defaults.
  op <- occ_panel %||% list()
  spec <- .occ_panel_spec(
    key = ek,
    national_zero_fill = nzf,
    national_zero_threshold = if (nzf) nzt else 0L,
    national_zero_rule = op[["national_zero_rule"]] %||% "",
    at_risk_start = op[["at_risk_start"]] %||% "season",
    adm1_overlap = op[["adm1_overlap"]] %||% "drop"
  )
  ep <- list(baseline = baseline, episodes = .episodes_block(spec))
  out <- list(
    episode_key = ek,
    surv_path = as.character(make_surv_dat_file(ep)),
    mon_path = as.character(make_monthly_dat_file(ep))
  )
  assert_occ_panel_variant(
    ep,
    out$surv_path,
    track = paste(tag, "observed survival panel")
  )
  sim_episode_check(
    preps = list(surv = out$surv_path, mon_data = out$mon_path),
    what = paste(tag, "observed panels"),
    episode_key = ek
  )
  out
}

#' hindcast_obs_panel
#'
#' The survival panel reduced to at-risk (gid, month, y) rows (y = any onset),
#' via load_observed_from_panel(). Loaded once and shared by every run.
hindcast_obs_panel <- function(surv_panel) {
  if (is.character(surv_panel)) {
    surv_panel <- readRDS(surv_panel)
  }
  load_observed_from_panel(surv_panel) |>
    dplyr::transmute(gid, month, y = as.integer(Y_obs > 0))
}

#' hindcast_load_run
#'
#' Read one run once and return everything the scorers need.
#'
#' @param row one row of hindcast_run_index(), or several rows of the same
#'   (group, launch_year): a fit group split across simulation units (v2
#'   manifest) is read from every unit run and pooled into one run
#' @param obs_panel hindcast_obs_panel() frame
#' @param mon_all build_mon_obs_frame() frame
#' @param gid_pop load_gid_population() frame
#' @param lead_breaks lead-class breaks
#' @return list(meta, cells, sim, bases, monthly_freq)
#'   cells: at-risk (gid, month) cells: lead, lead_class, p_hat, y_obs, pit
#'   sim:   monthly_summary rows (gid, month, replicate_id, Y, cases_month)
#'   bases: obs_sim_bases() list (cases track)
hindcast_load_run <- function(
  row,
  obs_panel,
  mon_all,
  gid_pop,
  lead_breaks = c(7L, 12L)
) {
  # Season-window index rows (country_launch): one fit group's season year is
  # pooled from runs launched at different months; each country is read from
  # its own launch to its own firm horizon and its leads count from its launch.
  season <- "country_launch" %in% names(row)
  if (
    length(unique(row$group)) != 1L ||
      length(unique(row$launch_year)) != 1L ||
      (!season && length(unique(row$launch)) != 1L) ||
      length(unique(row$horizon_end)) != 1L
  ) {
    stop(
      "[hindcast-skill] hindcast_load_run: rows must share group, launch year, launch ",
      "and horizon_end.",
      call. = FALSE
    )
  }
  parts <- lapply(seq_len(nrow(row)), function(k) .hindcast_read_part(row[k, ]))
  c_launch <- NULL
  if (season) {
    c_launch <- as.Date(unlist(unname(row$country_launch)))
    names(c_launch) <- names(unlist(unname(row$country_launch)))
    c_hz <- as.Date(unlist(unname(row$country_horizons)))
    names(c_hz) <- names(unlist(unname(row$country_horizons)))
    if (anyDuplicated(names(c_launch))) {
      stop(
        "[hindcast-skill] hindcast_load_run: a country is scored in two runs of ",
        row$group[1L],
        "/",
        row$launch_year[1L],
        ".",
        call. = FALSE
      )
    }
    # observed cases inside each country's own season window only
    iso_m <- gid_iso3(mon_all$gid)
    mon_all <- mon_all[
      iso_m %in%
        names(c_launch) &
        mon_all$month >= c_launch[iso_m] &
        mon_all$month <= c_hz[iso_m],
      ,
      drop = FALSE
    ]
    row$launch <- min(row$launch)
  }
  rows_all <- row
  row <- row[1L, ]
  row$run_tag <- paste(
    unique(unlist(lapply(parts, `[[`, "run_tag"))),
    collapse = ";"
  )
  tag <- sprintf("[hindcast-skill] %s/%d", row$group, row$launch_year)
  ms <- dplyr::bind_rows(lapply(parts, `[[`, "ms"))
  pit <- dplyr::bind_rows(lapply(parts, `[[`, "pit"))
  met <- list(
    crps = dplyr::bind_rows(lapply(parts, function(p) p$met$crps)),
    crps_total_cases = list(
      country = dplyr::bind_rows(lapply(parts, function(p) {
        p$met$crps_total_cases$country
      }))
    ),
    monthly_freq = dplyr::bind_rows(lapply(parts, function(p) {
      p$met$monthly_freq
    }))
  )
  if (!nrow(met$monthly_freq)) {
    met$monthly_freq <- NULL
  }

  # One support for both tracks: the survival panel's at-risk cells of the sim
  # gids inside the window. 07c's PIT table must be exactly that set.
  sim_cells <- dplyr::distinct(ms, gid, month)
  atrisk <- obs_panel |>
    dplyr::filter(month >= row$launch, month <= row$horizon_end) |>
    dplyr::semi_join(sim_cells, by = c("gid", "month"))
  k_pit <- paste(pit$gid, pit$month)
  k_pan <- paste(atrisk$gid, atrisk$month)
  if (!setequal(k_pit, k_pan)) {
    stop(
      sprintf(
        paste0(
          "%s metrics.qs$pit has %d cells but the survival panel gives %d ",
          "at-risk cells (%d differ). Re-run 07c --hindcast with the ",
          "deploy's OBS_PANEL."
        ),
        tag,
        length(k_pit),
        length(k_pan),
        length(setdiff(union(k_pit, k_pan), intersect(k_pit, k_pan)))
      ),
      call. = FALSE
    )
  }

  cells <- pit |>
    dplyr::select(gid, month, p_hat, y_obs, pit) |>
    dplyr::mutate(
      lead = hindcast_lead_months(month, row$launch, c_launch, gid_iso3(gid)),
      lead_class = hindcast_lead_class(lead, lead_breaks),
      year = lubridate::year(month)
    )

  b <- obs_sim_bases(ms, mon_all, gid_pop)
  # leads on the case-track frames too (from each country's launch under season
  # windows); the scorers that recompute lead_class overwrite it identically
  lead_on <- function(df) {
    if (is.null(df) || !nrow(df)) {
      return(df)
    }
    lead <- map_distinct(
      tibble::tibble(country = df$country, month = df$month),
      function(u) {
        hindcast_lead_months(u$month, row$launch, c_launch, u$country)
      }
    )
    dplyr::mutate(
      df,
      lead = lead,
      lead_class = hindcast_lead_class(lead, lead_breaks)
    )
  }
  b$sim_in <- lead_on(b$sim_in)
  b$obs_month <- lead_on(b$obs_month)
  n_rep <- length(unique(ms$replicate_id))
  ks <- if (nrow(pit)) {
    unname(suppressWarnings(stats::ks.test(pit$pit, "punif")$statistic))
  } else {
    NA_real_
  }
  meta <- tibble::tibble(
    group = row$group,
    launch_year = row$launch_year,
    launch = row$launch,
    horizon_end = row$horizon_end,
    run_tag = row$run_tag,
    episode_key = row$episode_key,
    n_rep = n_rep,
    n_gid = length(unique(ms$gid)),
    n_cells = nrow(cells),
    n_gid_atrisk = length(unique(cells$gid)),
    n_pit = nrow(pit),
    pit_ks = ks,
    cover_popwt = b$cover_popwt,
    n_sim_gids = length(b$sim_gids),
    n_adm1_touched = length(b$R_sim),
    mean_crps_count_legacy = if (!is.null(met$crps) && nrow(met$crps)) {
      mean(met$crps$crps_count, na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_crps_total_cases_legacy = {
      v <- met$crps_total_cases$country$crps_cases
      if (length(v)) mean(v, na.rm = TRUE) else NA_real_
    }
  )
  out <- list(
    meta = meta,
    cells = cells,
    sim = dplyr::select(
      ms,
      dplyr::any_of(c("gid", "month", "replicate_id", "Y", "cases_month"))
    ),
    bases = b,
    monthly_freq = met$monthly_freq
  )
  # season windows: each country's launch, and the observed cases cut to the
  # countries' windows (the case scorers read these instead of the bundle's)
  if (season) {
    out$country_launch <- c_launch
    out$country_horizon <- c_hz
    out$mon_obs <- mon_all
  }
  out$country_windows <- .hindcast_country_windows(rows_all, row, ms, c_launch)
  out
}

# Each scored country's window of a loaded run: tibble(iso3, launch, season_end,
# firm_horizon_end, floored). Season windows: from the index's per-country
# columns; otherwise every country shares the run's launch, its firm horizon is
# its own (per-country season ends) or the run's, season_end / floored are NA.
.hindcast_country_windows <- function(rows_all, row, ms, c_launch) {
  pick <- function(col) {
    if (!col %in% names(rows_all)) {
      return(NULL)
    }
    v <- unlist(unname(rows_all[[col]]))
    if (is.null(v)) NULL else v
  }
  hz <- pick("country_horizons")
  isos <- if (!is.null(c_launch)) {
    names(c_launch)
  } else if (!is.null(pick("countries"))) {
    unique(pick("countries"))
  } else {
    sort(unique(gid_iso3(unique(ms$gid))))
  }
  se <- pick("country_season_end")
  fl <- pick("country_floored")
  tibble::tibble(
    iso3 = isos,
    launch = if (!is.null(c_launch)) {
      unname(c_launch[isos])
    } else {
      rep(as.Date(row$launch), length(isos))
    },
    season_end = if (!is.null(se)) {
      as.Date(unname(se[isos]))
    } else {
      as.Date(rep(NA_character_, length(isos)))
    },
    firm_horizon_end = if (!is.null(hz)) {
      as.Date(unname(hz[isos]))
    } else {
      rep(as.Date(row$horizon_end), length(isos))
    },
    floored = if (!is.null(fl)) unname(fl[isos]) else rep(NA, length(isos))
  )
}

#' .hindcast_read_part -- one index row's monthly summary, metrics and PIT
#' cells. A fit-group row of a unit run (v2 manifest) keeps only that group's
#' countries, inside its own firm window (07c scored the PIT the same way).
.hindcast_read_part <- function(row) {
  floor_m <- function(m) lubridate::floor_date(as.Date(m), "month")
  ms <- sim_io_read(row$monthly_summary_path, row$out_fmt) |>
    dplyr::mutate(gid = as.character(gid), month = map_distinct(month, floor_m))
  met <- sim_io_read(row$metrics_path, sim_metrics_fmt(row$out_fmt))
  pit <- met$pit |>
    dplyr::mutate(
      gid = as.character(gid),
      month = lubridate::floor_date(as.Date(month), "month")
    )
  cs <- if ("countries" %in% names(row)) row$countries[[1L]] else NULL
  ch <- if ("country_horizons" %in% names(row)) {
    row$country_horizons[[1L]]
  } else {
    NULL
  }
  if (!is.null(cs)) {
    # each country inside its own firm window (its horizon cohort); without
    # per-country horizons, the row's window
    keep <- function(df) {
      iso <- gid_iso3(df$gid)
      h <- if (is.null(ch)) {
        rep(row$horizon_end, nrow(df))
      } else {
        as.Date(unname(ch[iso]))
      }
      iso %in% cs & !is.na(h) & df$month <= h & df$month <= row$horizon_end
    }
    ms <- ms[keep(ms), , drop = FALSE]
    pit <- pit[keep(pit), , drop = FALSE]
  }
  list(ms = ms, pit = pit, met = met, run_tag = row$run_tag)
}

# ===========================================================================
# 3. Scoring primitives (lifted from the notebook's chunks)
# ===========================================================================

#' crps_discrete -- exact CRPS of a pmf on support 0..K against integer y.
crps_discrete <- function(pmf, y) {
  Fk <- cumsum(pmf)
  k <- seq_along(pmf) - 1L
  sum((Fk - (y <= k))^2)
}

#' crps_binom -- exact CRPS of Binomial(m, p) against y (the LOO references).
crps_binom <- function(m, p, y) {
  if (is.na(p) || is.na(m)) {
    return(NA_real_)
  }
  crps_discrete(stats::dbinom(0:m, m, p), y)
}

#' brier / logscore / auc_mw -- the notebook's scores chunk, verbatim.
brier <- function(p, y) mean((p - y)^2)
logscore <- function(p, y, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}
auc_mw <- function(p, y) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) {
    return(NA_real_)
  }
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
#' bss_climatology -- the notebook's BSS (reference = the pooled base rate).
bss_climatology <- function(p, y) {
  b <- mean(y)
  1 - brier(p, y) / (b * (1 - b))
}

#' reliability_curve -- decile (equal-count) bins of p (notebook rel_curve).
reliability_curve <- function(p, y, nb = 10L) {
  ok <- !is.na(p)
  p <- p[ok]
  y <- y[ok]
  br <- unique(stats::quantile(
    p,
    probs = seq(0, 1, length.out = nb + 1),
    na.rm = TRUE
  ))
  if (length(br) < 3) {
    return(tibble::tibble())
  }
  tibble::tibble(
    p = p,
    y = y,
    bin = cut(p, breaks = br, include.lowest = TRUE)
  ) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      mean_p = mean(p),
      obs_freq = mean(y),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(bin = as.character(bin))
}

#' roc_points -- ROC at every distinct threshold (ties grouped, so the curve has
#' at most one point per distinct p; equals the notebook's per-cell curve at
#' its vertices).
roc_points <- function(p, y) {
  ok <- !is.na(p)
  p <- p[ok]
  y <- y[ok]
  d <- tibble::tibble(p = p, y = y) |>
    dplyr::group_by(p) |>
    dplyr::summarise(tp = sum(y), fp = sum(1 - y), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(p))
  tibble::tibble(
    threshold = c(Inf, d$p),
    fpr = c(0, cumsum(d$fp) / sum(d$fp)),
    tpr = c(0, cumsum(d$tp) / sum(d$tp))
  )
}

#' coverage_at_width -- central-interval coverage from PIT values.
coverage_at_width <- function(pit, widths) {
  vapply(
    widths,
    function(w) mean(pit >= (1 - w) / 2 & pit <= (1 + w) / 2, na.rm = TRUE),
    numeric(1)
  )
}

#' pit_hist_bins -- counts of PIT values in `bins` equal-width bins.
pit_hist_bins <- function(pit, bins = 20L) {
  br <- seq(0, 1, length.out = bins + 1)
  h <- graphics::hist(pmin(pmax(pit, 0), 1), breaks = br, plot = FALSE)
  tibble::tibble(
    lo = br[-length(br)],
    hi = br[-1],
    n = h$counts,
    density = h$counts / max(1, sum(h$counts)) * bins
  )
}

#' randomised_rank -- verification rank with ties broken uniformly (the count
#' analogue of the PIT; rank-hist chunk). Returns u = rank / n in [0, 1].
randomised_rank <- function(ens, obs) {
  below <- sum(ens < obs)
  eq <- sum(ens == obs)
  (below + floor(stats::runif(1) * (eq + 1))) / length(ens)
}

#' interval_score -- Winkler score of the central interval [l, u] of width w.
interval_score <- function(l, u, y, width) {
  a <- 1 - width
  (u - l) + (2 / a) * (l - y) * (y < l) + (2 / a) * (y - u) * (y > u)
}

#' skill_ratio -- 1 - sum(model)/sum(ref) over rows where both are finite
#' (the notebook's crpss_at rule).
skill_ratio <- function(model, ref) {
  ok <- is.finite(model) & is.finite(ref)
  if (!any(ok) || .ref_score_zero(sum(ref[ok]))) {
    return(NA_real_)
  }
  1 - sum(model[ok]) / sum(ref[ok])
}

# A summed reference score that is zero up to round-off: a
# reference that is exact on every scored row (e.g. a calendar-month
# climatology that is a point mass at the observed zero) leaves the skill
# score undefined (NA), not +-1e17 from a -3e-17 sum. Scores are >= 0 and in
# case or probability units, so 1e-9 is far below any real reference score.
.ref_score_zero <- function(s) is.finite(s) & abs(s) < 1e-9

# ---------------------------------------------------------------------------
# 3b. Vectorised ensemble kernel: one sort per unit serves every rule
# ---------------------------------------------------------------------------
# Every ensemble score of the module (occurrence units, case totals, monthly
# units, every reference) goes through these. A set of ensembles is a matrix,
# one row per unit; ragged sets (references) are NA-padded on the right
# (pad_rows()). log1p() and pmax(., t) are monotone, so the one row sort is
# reused by every rule. Equal to crps_fair() / crps_ens() / stats::quantile()
# up to summation order (test_hindcast_skill_primitives.R).

#' sort_rows -- each row sorted ascending, NA last (one radix order over
#' (row, value) instead of one sort() per row).
sort_rows <- function(M) {
  if (!length(M)) {
    return(M)
  }
  matrix(
    M[order(row(M), M, na.last = TRUE, method = "radix")],
    nrow = nrow(M),
    byrow = TRUE
  )
}

#' pad_rows -- list of numeric vectors -> NA-padded matrix, one row each.
pad_rows <- function(x) {
  n <- lengths(x)
  M <- matrix(NA_real_, length(x), max(c(n, 1L)))
  if (sum(n)) {
    M[cbind(rep.int(seq_along(x), n), sequence(n))] <- as.numeric(unlist(
      x,
      use.names = FALSE
    ))
  }
  M
}

#' crps_sorted
#'
#' Ensemble CRPS of every row of a row-sorted, NA-padded matrix S against y:
#' the sorted energy form of crps_ens() (fair = FALSE) or crps_fair() (fair =
#' TRUE: the unbiased pairwise term, Ferro 2014). n = members per row; a
#' 1-member row is |x - y|, an empty row NA.
crps_sorted <- function(S, y, n = rowSums(!is.na(S)), fair = TRUE) {
  if (!nrow(S)) {
    return(numeric(0))
  }
  mae <- rowSums(abs(S - y), na.rm = TRUE) / n
  den <- if (fair) n * (n - 1) else n^2
  out <- mae - rowSums((2 * col(S) - n - 1) * S, na.rm = TRUE) / den
  out[n == 1L] <- mae[n == 1L]
  out[n == 0L] <- NA_real_
  out
}

#' quantile7_sorted -- stats::quantile(type = 7) of every row of a row-sorted,
#' NA-padded matrix (same interpolation); a units x length(p) matrix.
quantile7_sorted <- function(S, p, n = rowSums(!is.na(S))) {
  out <- matrix(NA_real_, nrow(S), length(p))
  ok <- which(n > 0)
  for (k in seq_along(p)) {
    idx <- 1 + pmax(n[ok] - 1, 0) * p[k]
    lo <- floor(idx)
    hi <- ceiling(idx)
    xl <- S[cbind(ok, lo)]
    xh <- S[cbind(ok, hi)]
    i <- which(idx > lo & xh != xl)
    h <- (idx - lo)[i]
    xl[i] <- (1 - h) * xl[i] + h * xh[i]
    out[ok, k] <- xl
  }
  out
}

#' score_sorted
#'
#' Scores of every row of a row-sorted matrix S against y, one numeric vector
#' per rule:
#'   crps_raw     fair CRPS of x
#'   crps_log1p   fair CRPS of log1p(x) vs log1p(y) (Bosse et al. 2023)
#'   twcrps       fair threshold-weighted CRPS, weight 1{z > tw}: the CRPS of
#'                pmax(x, tw) vs pmax(y, tw) (Gneiting & Ranjan 2011; the
#'                chaining form of Allen et al. 2023)
#'   crps_energy  energy-form CRPS of x (the legacy crps / crps_emp)
#'   pinball      mean over `taus` of the type-7 quantile loss
score_sorted <- function(
  S,
  y,
  rules,
  n = rowSums(!is.na(S)),
  taus = NULL,
  tw = NULL
) {
  out <- list()
  for (ru in rules) {
    out[[ru]] <- switch(
      ru,
      crps_raw = crps_sorted(S, y, n),
      crps_log1p = crps_sorted(log1p(S), log1p(y), n),
      twcrps = crps_sorted(pmax(S, tw), pmax(y, tw), n),
      crps_energy = crps_sorted(S, y, n, fair = FALSE),
      pinball = rowMeans(pinball_loss(
        quantile7_sorted(S, taus, n),
        y,
        rep(taus, each = nrow(S))
      )),
      stop("score_sorted: unknown rule '", ru, "'")
    )
  }
  out
}

#' rows_by_unit -- long (unit, value) pairs -> NA-padded matrix with one row
#' per unit 1..n_units (members in input order).
rows_by_unit <- function(unit, value, n_units) {
  o <- order(unit, method = "radix")
  u <- unit[o]
  M <- matrix(NA_real_, n_units, max(c(tabulate(u, n_units), 1L)))
  if (length(u)) {
    M[cbind(u, sequence(tabulate(u, n_units)[unique(u)]))] <- value[o]
  }
  M
}

#' sim_cells
#'
#' A replicate frame (gid, month, replicate_id, `value`) as distinct (gid,
#' month) cells x replicates, built once per run: every per-replicate sum of
#' the scorers is then a rowsum() over cell rows of X, never a pass over the
#' replicate rows (one row per gid x month x replicate).
#' @return list(gid, day (as.integer(month)) per cell, cell = the cell of every
#'   input row, reps (sorted replicate ids), X = cells x replicates matrix of
#'   `value` (summed over duplicate (cell, replicate) rows, under na.rm; 0 where
#'   a cell has no row of that replicate), present = cells x replicates logical,
#'   NULL when every (cell, replicate) has a row)
sim_cells <- function(df, value, na.rm = FALSE) {
  gk <- unique(df$gid)
  day <- as.integer(df$month)
  dk <- sort(unique(day))
  nd <- length(dk)
  gd <- (match(df$gid, gk) - 1L) * nd + match(day, dk)
  used <- which(tabulate(gd, length(gk) * nd) > 0L)
  pos <- integer(length(gk) * nd)
  pos[used] <- seq_along(used)
  cell <- pos[gd]
  rm(gd, day)
  reps <- sort(unique(df$replicate_id))
  nC <- length(used)
  R <- length(reps)
  cr <- cell + (match(df$replicate_id, reps) - 1L) * nC
  cnt <- tabulate(cr, nC * R)
  X <- matrix(0, nC, R)
  if (all(cnt <= 1L)) {
    X[cr] <- df[[value]]
  } else {
    X[sort(unique(cr))] <- rowsum(df[[value]], cr, na.rm = na.rm)[, 1]
  }
  list(
    gid = gk[(used - 1L) %/% nd + 1L],
    day = dk[(used - 1L) %% nd + 1L],
    cell = cell,
    reps = reps,
    X = X,
    present = if (any(cnt == 0L)) matrix(cnt > 0L, nC, R)
  )
}

# ===========================================================================
# 4. Occurrence track
# ===========================================================================

# Occurrence variant: support = at-risk (gid, month) rows of the survival
# panel inside each run's window (identical to metrics.qs$pit -- asserted);
# target = onset count over the run's WINDOW per gid ("window"), the
# calendar-year cut kept as "calendar_year"; references "insample_gid",
# "loo_gid" (headline: per-gid base rate from at-risk months OUTSIDE the
# scored window), "loo_gid_fullhist", "regional", "zero", "persistence";
# stratified by lead class.
#' onset_window_counts
#'
#' Per-replicate and observed onset counts over a run's at-risk cells, at one
#' grain: "window" (one unit per run x gid -- the fair target) or
#' "calendar_year" (one unit per run x gid x calendar year -- the legacy
#' sibling). Every lead class gets its own unit set
#' plus "all".
#' @param sim_index sim_cells(run$sim, "Y"), built once per run (both grains)
#' @return list(ens = units x replicates matrix of onset sums sc, one row per
#'   obs row, NA where the unit has no sim row of that replicate; reps =
#'   replicate ids of its columns;
#'   obs = (gid, unit, lead_class, n_atrisk, y_obs, months))
onset_window_counts <- function(
  run,
  grain = c("window", "calendar_year"),
  sim_index = NULL
) {
  grain <- match.arg(grain)
  unit_of <- function(df) {
    if (grain == "window") rep(run$meta$launch_year, nrow(df)) else df$year
  }
  cells <- run$cells |> dplyr::mutate(unit = unit_of(run$cells))
  cl <- .with_all_lead(dplyr::select(
    cells,
    gid,
    month,
    unit,
    lead_class,
    y_obs
  ))
  obs <- cl |>
    dplyr::group_by(gid, unit, lead_class) |>
    dplyr::summarise(
      n_atrisk = dplyr::n(),
      y_obs = sum(y_obs),
      months = list(month),
      .groups = "drop"
    )
  # replicate onset sums per unit: every sim cell adds to its own lead class
  # and to "all" (one rowsum() of the cell rows per half, no join of the
  # summary)
  si <- sim_index %||% sim_cells(run$sim, "Y")
  uk <- paste(obs$gid, obs$unit, obs$lead_class, sep = "\r")
  U <- length(uk)
  R <- length(si$reps)
  ck <- match(
    paste(cl$gid, as.integer(as.Date(cl$month))),
    paste(si$gid, si$day)
  )
  cu <- match(paste(cl$gid, cl$unit, cl$lead_class, sep = "\r"), uk)
  E <- matrix(NA_real_, U, R)
  for (w in split(which(!is.na(ck)), cl$lead_class[!is.na(ck)] == "all")) {
    s <- rowsum(si$X[ck[w], , drop = FALSE], cu[w])
    u <- as.integer(rownames(s))
    if (!is.null(si$present)) {
      s[rowsum(si$present[ck[w], , drop = FALSE] + 0, cu[w]) == 0] <- NA
    }
    E[u, ] <- s
  }
  list(ens = E, reps = si$reps, obs = obs)
}

#' score_occurrence_model
#'
#' Model scores per unit: ensemble CRPS (energy form, the legacy crps_emp) and
#' fair CRPS, ensemble quantiles, mean, P(>=1), randomised rank.
#'
#' defer_rank = TRUE leaves rank_u NA and appends the tie-break inputs
#' (.rank_below, .rank_eq, .rank_n), so the one random draw per unit can be made
#' later by fill_randomised_ranks(), in this same row order. That is what lets
#' hindcast_skill_bundle() score runs on parallel workers and still consume the
#' random stream exactly as a serial loop does.
score_occurrence_model <- function(
  counts,
  probs = c(.025, .25, .5, .75, .975),
  defer_rank = FALSE
) {
  ob <- counts$obs
  # units with an ensemble member, in (gid, unit, lead_class) order, one
  # ensemble row each (vectorised kernel), NA-padded to the largest ensemble
  has <- which(rowSums(!is.na(counts$ens)) > 0L)
  has <- has[order(
    ob$gid[has],
    ob$unit[has],
    ob$lead_class[has],
    method = "radix"
  )]
  uk <- o <- ob[has, , drop = FALSE]
  S <- sort_rows(counts$ens[has, , drop = FALSE])
  n <- rowSums(!is.na(S))
  y <- o$y_obs
  w <- max(c(n, 1L))
  S <- cbind(S, matrix(NA_real_, nrow(S), max(w - ncol(S), 0L)))[,
    seq_len(w),
    drop = FALSE
  ]
  q <- quantile7_sorted(S, probs, n)
  colnames(q) <- sprintf("q%03d", round(1000 * probs))
  out <- dplyr::bind_cols(
    tibble::tibble(
      gid = uk$gid,
      unit = uk$unit,
      lead_class = uk$lead_class,
      n_atrisk = o$n_atrisk,
      y_obs = y,
      months = o$months,
      n_rep = as.integer(n),
      mean_sc = rowSums(S, na.rm = TRUE) / n,
      p_pos = rowSums(S > 0, na.rm = TRUE) / n
    ),
    tibble::as_tibble(q),
    tibble::tibble(
      rank_u = NA_real_,
      crps_emp = crps_sorted(S, y, n, fair = FALSE),
      crps_fair = crps_sorted(S, y, n),
      .rank_below = as.integer(rowSums(S < y, na.rm = TRUE)),
      .rank_eq = as.integer(rowSums(S == y, na.rm = TRUE)),
      .rank_n = as.integer(n)
    )
  )
  # the non-deferred rank draws one runif per unit in this row order, exactly
  # as randomised_rank() unit by unit
  if (defer_rank) out else fill_randomised_ranks(out)
}

#' fill_randomised_ranks -- the deferred tie-break draws of
#' score_occurrence_model(defer_rank = TRUE): one runif per row, in row order,
#' with randomised_rank()'s arithmetic; drops the .rank_* helper columns.
fill_randomised_ranks <- function(units) {
  units$rank_u <- (units$.rank_below +
    floor(stats::runif(nrow(units)) * (units$.rank_eq + 1))) /
    units$.rank_n
  dplyr::select(units, -c(.rank_below, .rank_eq, .rank_n))
}

#' occ_reference_rate
#'
#' Per-gid reference onset probability per at-risk month for each scored unit,
#' leaving out every month of the scored run's window (overlapping launches
#' share calendar months, so leaving out only the scored launch would leak).
#'   loo_gid           at-risk cells of the gid in the hindcast cell pool
#'   loo_gid_fullhist  at-risk months of the gid anywhere in the survival panel
#'   regional          pooled over the group's gids (hindcast cell pool)
#' @param pool distinct (group, gid, month, y) hindcast at-risk cells over all runs
#' @param obs_panel hindcast_obs_panel() frame (for fullhist)
#' @param windows (group, launch_year, launch, horizon_end) per run
#' @return (group, launch_year, gid, p_ref, m_ref) -- m_ref = months the rate rests on
occ_reference_rate <- function(reference, pool, obs_panel, windows, units) {
  keys <- dplyr::distinct(units, group, launch_year, gid) |>
    dplyr::inner_join(windows, by = c("group", "launch_year"))
  rate <- function(src, by) {
    keys |>
      dplyr::inner_join(src, by = by, relationship = "many-to-many") |>
      dplyr::filter(month < launch | month > horizon_end) |>
      dplyr::group_by(group, launch_year, gid) |>
      dplyr::summarise(
        p_ref = sum(y) / dplyr::n(),
        m_ref = dplyr::n(),
        .groups = "drop"
      )
  }
  out <- switch(
    reference,
    loo_gid = rate(dplyr::select(pool, gid, month, y), "gid"),
    loo_gid_fullhist = rate(dplyr::select(obs_panel, gid, month, y), "gid"),
    regional = {
      kk <- keys
      r <- kk |>
        dplyr::distinct(group, launch_year, launch, horizon_end) |>
        dplyr::inner_join(
          dplyr::select(pool, group, month, y),
          by = "group",
          relationship = "many-to-many"
        ) |>
        dplyr::filter(month < launch | month > horizon_end) |>
        dplyr::group_by(group, launch_year) |>
        dplyr::summarise(
          p_ref = sum(y) / dplyr::n(),
          m_ref = dplyr::n(),
          .groups = "drop"
        )
      dplyr::inner_join(
        dplyr::distinct(kk, group, launch_year, gid),
        r,
        by = c("group", "launch_year")
      )
    },
    stop("occ_reference_rate: no rate for reference '", reference, "'")
  )
  dplyr::left_join(
    dplyr::distinct(keys, group, launch_year, gid),
    out,
    by = c("group", "launch_year", "gid")
  )
}

#' score_occurrence_references
#'
#' Reference scores for every unit of `units` (rows of score_occurrence_model()
#' bound across runs, with group/launch_year/grain/months), per reference x rule.
#'   insample_gid  empirical ensemble of the gid's own y_obs over every unit of
#'                 the same grain and lead class (observation included; >= 2
#'                 units) -- legacy
#'   loo_gid / loo_gid_fullhist / regional  Binomial(n_atrisk, p_ref), exact CRPS
#'   zero          point mass at 0 (CRPS = y)
#'   persistence   observed onsets in the same months one year earlier (point
#'                 forecast; NA when none of those months was at risk)
#' Rules: "crps_emp" scores empirical references with crps_ens(), "crps_fair"
#' with crps_fair(); parametric and point references are exact under both.
score_occurrence_references <- function(
  units,
  references,
  rules,
  pool,
  obs_panel,
  windows
) {
  base <- dplyr::select(
    units,
    group,
    launch_year,
    grain,
    gid,
    unit,
    lead_class,
    n_atrisk,
    y_obs,
    months
  )
  out <- list()
  for (ref in references) {
    sc <- switch(
      ref,
      insample_gid = {
        base |>
          dplyr::group_by(grain, lead_class, gid) |>
          dplyr::filter(dplyr::n() >= 2) |>
          dplyr::mutate(
            s_emp = vapply(y_obs, function(y) crps_ens(y_obs, y), numeric(1)),
            s_fair = vapply(y_obs, function(y) crps_fair(y_obs, y), numeric(1)),
            p_ref = mean(y_obs / n_atrisk)
          ) |>
          dplyr::ungroup()
      },
      zero = dplyr::mutate(base, s_emp = y_obs, s_fair = y_obs, p_ref = 0),
      persistence = {
        prev <- obs_panel |>
          dplyr::mutate(month = .shift_years(month, 1L)) |>
          dplyr::rename(y_prev = y)
        pr <- base |>
          dplyr::select(
            group,
            launch_year,
            grain,
            gid,
            unit,
            lead_class,
            months
          ) |>
          tidyr::unnest(months) |>
          dplyr::rename(month = months) |>
          dplyr::left_join(prev, by = c("gid", "month")) |>
          dplyr::group_by(group, launch_year, grain, gid, unit, lead_class) |>
          dplyr::summarise(
            f = if (all(is.na(y_prev))) NA_real_ else sum(y_prev, na.rm = TRUE),
            .groups = "drop"
          )
        base |>
          dplyr::left_join(
            pr,
            by = c("group", "launch_year", "grain", "gid", "unit", "lead_class")
          ) |>
          dplyr::mutate(
            s_emp = abs(f - y_obs),
            s_fair = s_emp,
            p_ref = pmin(f / n_atrisk, 1)
          ) |>
          dplyr::select(-f)
      },
      {
        r <- occ_reference_rate(ref, pool, obs_panel, windows, base)
        base |>
          dplyr::left_join(r, by = c("group", "launch_year", "gid")) |>
          dplyr::mutate(
            s_emp = mapply(crps_binom, n_atrisk, p_ref, y_obs),
            s_fair = s_emp
          )
      }
    )
    for (ru in rules) {
      use_fair <- identical(ru, "crps_fair")
      out[[length(out) + 1L]] <- sc |>
        dplyr::transmute(
          group,
          launch_year,
          grain,
          gid,
          unit,
          lead_class,
          reference = ref,
          rule = ru,
          score_ref = if (use_fair) s_fair else s_emp,
          p_ref
        )
    }
  }
  dplyr::bind_rows(out)
}

#' occ_cell_reference_prob
#'
#' Per-cell (monthly) reference probability for the Brier / log-score family,
#' same references as score_occurrence_references() ("insample_gid" = the gid's
#' base rate over all its hindcast cells; "persistence" = onset 12 months
#' earlier, NA when not at risk).
occ_cell_reference_prob <- function(
  cells,
  reference,
  pool,
  obs_panel,
  windows
) {
  switch(
    reference,
    insample_gid = cells |>
      dplyr::group_by(gid) |>
      dplyr::mutate(p_ref = mean(y_obs)) |>
      dplyr::ungroup(),
    zero = dplyr::mutate(cells, p_ref = 0),
    persistence = cells |>
      dplyr::left_join(
        obs_panel |>
          dplyr::transmute(
            gid,
            month = .shift_years(month, 1L),
            p_ref = as.numeric(y)
          ),
        by = c("gid", "month")
      ),
    {
      r <- occ_reference_rate(
        reference,
        pool,
        obs_panel,
        windows,
        dplyr::distinct(cells, group, launch_year, gid)
      )
      dplyr::left_join(cells, r, by = c("group", "launch_year", "gid"))
    }
  )
}

#' score_cells
#'
#' Monthly-cell scores of one subset: Brier, the legacy pooled-base-rate BSS,
#' BSS vs the cell reference (when p_ref is given), log score (floored at eps),
#' AUC, sharpness, n_hard_zero (onsets on p = 0 cells).
score_cells <- function(p, y, p_ref = NULL, eps = 1e-6) {
  ok <- !is.na(p)
  out <- tibble::tibble(
    n = sum(ok),
    base_rate = mean(y[ok]),
    brier = brier(p[ok], y[ok]),
    bss_pooled = bss_climatology(p[ok], y[ok]),
    logscore = logscore(p[ok], y[ok], eps),
    auc = auc_mw(p[ok], y[ok]),
    sharpness = mean(p[ok]),
    n_hard_zero = sum(p[ok] == 0 & y[ok] == 1)
  )
  if (!is.null(p_ref)) {
    okr <- ok & !is.na(p_ref)
    out$n_ref <- sum(okr)
    out$brier_ref <- brier(p_ref[okr], y[okr])
    out$bss <- 1 - brier(p[okr], y[okr]) / out$brier_ref
    out$logscore_ref <- logscore(p_ref[okr], y[okr], eps)
    out$n_hard_zero_ref <- sum(p_ref[okr] == 0 & y[okr] == 1)
  }
  out
}

# ===========================================================================
# 5. Cases track
# ===========================================================================

# Cases variant: observed attribution kept (finest-level dedup + population
# share s); supports "obs_positive_keys", "exact_gid", "atrisk_gids"
# (headline: sim gids at risk in the window -- surveillance-, not
# outcome-defined), "atrisk_keys", "footprint"; references "loo_climatology",
# "zero", "persistence"; rules crps_raw / crps_log1p / pinball.
#' case_support
#'
#' Support gids + monthly observed cases for one support of one run.
#'   obs_positive_keys  old basis B: sim gids in region keys with observed cases;
#'                      obs = finest-deduped ADM1 cases x s(r)
#'   exact_gid          old basis C: sim gids with cases recorded at exactly them;
#'                      obs = those exact-gid cases
#'   atrisk_gids        sim gids with >= 1 at-risk (surveilled) month in the
#'                      window; obs = ADM1 cases x pop share of those gids
#'   atrisk_keys        every sim gid in a region key holding an at-risk gid
#'   footprint          every sim gid
#' For every support but exact_gid the observed side is the same population-
#' share attribution: region key K contributes obs_adm1(K) x min(pop(G n K) /
#' adm1_pop(K), 1), which is exactly obs_B for obs_positive_keys.
#' @return list(gids, obs_month = (region_key, country, month, cases_month),
#'   shares = (region_key, s_g) population shares, NULL for exact_gid)
case_support <- function(run, support, mon_all, gid_pop) {
  b <- run$bases
  if (support == "exact_gid") {
    om <- mon_all |>
      dplyr::filter(month >= b$HS, month <= b$HE, gid %in% b$G_obs) |>
      dplyr::transmute(
        region_key,
        country = stringr::str_sub(region_key, 1, 3),
        month,
        cases_month
      )
    return(list(gids = b$G_obs, obs_month = om, shares = NULL))
  }
  atrisk <- unique(run$cells$gid)
  gids <- switch(
    support,
    obs_positive_keys = b$simB_gids,
    atrisk_gids = intersect(b$sim_gids, atrisk),
    atrisk_keys = b$sim_gids[
      region_key_of(b$sim_gids) %in% region_key_of(atrisk)
    ],
    footprint = b$sim_gids,
    stop("case_support: unknown support '", support, "'")
  )
  avail_yr <- sort(unique(gid_pop$year))
  pop_y <- dplyr::filter(
    gid_pop,
    year == avail_yr[which.min(abs(avail_yr - b$yr))]
  )
  sh <- cover_shares(gids, pop_y) |> dplyr::select(region_key, s_g = s)
  keys <- unique(region_key_of(gids))
  om <- b$obs_month |>
    dplyr::filter(region_key %in% keys) |>
    dplyr::left_join(sh, by = "region_key") |>
    dplyr::mutate(cases_month = cases_month * dplyr::coalesce(s_g, 1)) |>
    dplyr::select(region_key, country, month, cases_month)
  list(gids = gids, obs_month = om, shares = sh)
}

#' case_sim_index -- the run's windowed simulated cases as sim_cells() (cells x
#' replicates, NA cases count 0) plus each cell's region key, built once per run
#' and shared by every support's case_unit_arrays(). Columns are the run's
#' replicates (run$sim), as the case scorers always counted them.
case_sim_index <- function(run) {
  s <- run$bases$sim_in
  si <- sim_cells(s, "cases_month", na.rm = TRUE)
  reps <- sort(unique(run$sim$replicate_id))
  if (!identical(si$reps, reps)) {
    # a replicate without a windowed row: a 0 column
    X <- matrix(0, nrow(si$X), length(reps))
    X[, match(si$reps, reps)] <- si$X
    si$X <- X
    si$reps <- reps
  }
  si$region_key <- s$region_key[match(seq_along(si$gid), si$cell)]
  si$cell <- NULL
  si$present <- NULL
  si
}

#' case_unit_arrays
#'
#' One (run, support): the support's simulated replicate sums and observed
#' cases per (ADM1 region key, month), the base of every case score. The legacy
#' window totals, the monthly views, the monthly units and the season totals
#' are all row sums of these (rowsum() on the key rows), never a second pass
#' over the monthly summary. Observed cases are the run's (cut to each
#' country's window under season windows).
#' @return list(keys = tibble(region_key, country, month, lead, lead_class,
#'   has_sim, has_obs), M = keys x replicates (0 where absent), obs (0 where
#'   absent), gids, shares (case_support()))
case_unit_arrays <- function(
  run,
  support,
  mon_all,
  gid_pop,
  lead_breaks,
  sim_index = NULL
) {
  cs <- case_support(run, support, run$mon_obs %||% mon_all, gid_pop)
  si <- sim_index %||% case_sim_index(run)
  w <- si$gid %in% cs$gids # sim cells of the support
  om <- cs$obs_month
  rk <- unique(c(si$region_key[w], om$region_key))
  rk <- rk[order(rk, method = "radix", na.last = TRUE)]
  key_of <- function(r, d) match(r, rk) * 1e5 + d
  ks <- key_of(si$region_key[w], si$day[w])
  ko <- key_of(om$region_key, as.integer(as.Date(om$month)))
  uk <- sort(unique(c(ks, ko)))
  K <- length(uk)
  R <- length(si$reps)
  fill <- function(n, g, x) {
    # sums of x by integer cell g, 0 elsewhere
    out <- numeric(n)
    if (length(g)) {
      out[sort(unique(g))] <- rowsum(x, g, na.rm = TRUE)[, 1]
    }
    out
  }
  # replicate sums per (region key, month): one rowsum() of the support's cells
  M <- matrix(0, K, R)
  if (any(w)) {
    s <- rowsum(
      if (all(w)) si$X else si$X[w, , drop = FALSE],
      match(ks, uk),
      na.rm = TRUE
    )
    M[as.integer(rownames(s)), ] <- s
  }
  keys <- tibble::tibble(
    region_key = rk[uk %/% 1e5],
    month = as.Date(uk %% 1e5, origin = "1970-01-01")
  )
  keys$country <- stringr::str_sub(keys$region_key, 1, 3)
  keys$lead <- hindcast_lead_months(
    keys$month,
    run$meta$launch,
    run$country_launch,
    keys$country
  )
  keys$lead_class <- hindcast_lead_class(keys$lead, lead_breaks)
  keys$has_sim <- uk %in% ks
  keys$has_obs <- uk %in% ko
  list(
    keys = keys,
    M = M,
    obs = fill(K, match(ko, uk), om$cases_month),
    gids = cs$gids,
    shares = cs$shares
  )
}

# (scope, lead class) groups of a support's window totals, in the legacy
# order: pairs with simulated rows (own class + "all"), sorted; then pairs only
# observed, sorted; then the region's "all" row if still absent.
.case_total_groups <- function(scope, lc, has_sim, has_obs, region_all = NULL) {
  pairs <- function(w) {
    p <- unique(data.frame(
      scope = c(scope[w], scope[w]),
      lead_class = c(lc[w], rep("all", sum(w))),
      stringsAsFactors = FALSE
    ))
    p[
      order(p$scope, p$lead_class, method = "radix", na.last = TRUE),
      ,
      drop = FALSE
    ]
  }
  ps <- pairs(has_sim)
  po <- pairs(has_obs)
  key <- function(p) paste(p$scope, p$lead_class, sep = "\r")
  g <- rbind(ps, po[!key(po) %in% key(ps), , drop = FALSE])
  if (
    !is.null(region_all) && !paste(region_all, "all", sep = "\r") %in% key(g)
  ) {
    g <- rbind(
      g,
      data.frame(
        scope = region_all,
        lead_class = "all",
        stringsAsFactors = FALSE
      )
    )
  }
  rownames(g) <- NULL
  g
}

# row sums of `X` (matrix or vector) over every key row's own lead class and
# "all", per group of `g` (a .case_total_groups() frame)
.case_sum_groups <- function(X, scope, lc, g) {
  X <- as.matrix(X)
  gk <- paste(g$scope, g$lead_class, sep = "\r")
  idx <- c(
    match(paste(scope, lc, sep = "\r"), gk),
    match(paste(scope, "all", sep = "\r"), gk)
  )
  out <- matrix(0, nrow(g), ncol(X))
  keep <- !is.na(idx)
  if (any(keep)) {
    s <- rowsum(rbind(X, X)[keep, , drop = FALSE], idx[keep])
    out[as.integer(rownames(s)), ] <- s
  }
  out
}

#' score_cases_run
#'
#' Window case totals per support x grain (region = the run's ENSO group,
#' country) x lead class: ensemble quantiles / mean / predictive percentile /
#' interval flags, and model scores under every rule. CRPS rules use the fair
#' estimator (the LOO reference has ~10 members); `crps` is the legacy energy-
#' form CRPS on raw totals (= the notebook's tcw$crps). Sums of
#' case_unit_arrays() (pass `arrays`, one per support, to reuse them).
#' @return tibble, one row per support x grain x scope x lead class
score_cases_run <- function(
  run,
  supports,
  mon_all,
  gid_pop,
  lead_breaks,
  taus,
  arrays = NULL
) {
  rows <- list()
  for (sp in supports) {
    A <- arrays[[sp]] %||%
      case_unit_arrays(run, sp, mon_all, gid_pop, lead_breaks)
    for (grain in c("region", "country")) {
      scope <- if (grain == "region") {
        rep(run$meta$group, nrow(A$keys))
      } else {
        A$keys$country
      }
      g <- .case_total_groups(
        scope,
        A$keys$lead_class,
        A$keys$has_sim,
        A$keys$has_obs,
        region_all = if (grain == "region") run$meta$group
      )
      S <- sort_rows(.case_sum_groups(A$M, scope, A$keys$lead_class, g))
      o <- .case_sum_groups(A$obs, scope, A$keys$lead_class, g)[, 1]
      n <- rep(ncol(S), nrow(S))
      q <- quantile7_sorted(S, c(.025, .25, .5, .75, .975), n)
      sc <- score_sorted(
        S,
        o,
        c("crps_energy", "crps_raw", "crps_log1p", "pinball"),
        n,
        taus = taus
      )
      rows[[length(rows) + 1L]] <- tibble::tibble(
        support = sp,
        grain = grain,
        scope = g$scope,
        lead_class = g$lead_class,
        obs = o,
        q025 = q[, 1],
        q25 = q[, 2],
        med = q[, 3],
        q75 = q[, 4],
        q975 = q[, 5],
        mean = rowSums(S) / n,
        pctile = 100 * rowSums(S <= o) / n,
        crps = sc$crps_energy,
        score_crps_raw = sc$crps_raw,
        score_crps_log1p = sc$crps_log1p,
        score_pinball = sc$pinball,
        n_support_gids = length(A$gids)
      )
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      group = run$meta$group,
      launch_year = run$meta$launch_year,
      cover_popwt = run$meta$cover_popwt,
      obs_naive = dplyr::if_else(
        grain == "region" &
          support == "obs_positive_keys" &
          lead_class == "all",
        run$bases$obs_naive,
        NA_real_
      ),
      in50 = obs >= q25 & obs <= q75,
      in95 = obs >= q025 & obs <= q975,
      pos = dplyr::case_when(
        obs < q025 ~ "over",
        obs > q975 ~ "under",
        TRUE ~ "within"
      ),
      .before = 1
    )
}

#' case_monthly_views
#'
#' Monthly ensemble median / 95 % band vs observed for one support, at region
#' and country grain (the notebook's detail_views monthly frames): row sums of
#' case_unit_arrays() over the months with simulated rows.
case_monthly_views <- function(
  run,
  support,
  mon_all,
  gid_pop,
  arrays = NULL,
  lead_breaks = c(7L, 12L)
) {
  A <- arrays %||% case_unit_arrays(run, support, mon_all, gid_pop, lead_breaks)
  k <- A$keys
  out <- list()
  for (grain in c("region", "country")) {
    gk <- if (grain == "region") {
      as.character(k$month)
    } else {
      paste(k$country, k$month, sep = "\r")
    }
    s <- k$has_sim
    ug <- unique(gk[s])
    ug <- ug[order(
      if (grain == "region") ug else k$country[s][match(ug, gk[s])],
      k$month[s][match(ug, gk[s])],
      method = "radix"
    )]
    S <- sort_rows(rowsum(A$M[s, , drop = FALSE], gk[s])[ug, , drop = FALSE])
    q <- quantile7_sorted(S, c(.5, .025, .975))
    ob <- if (any(k$has_obs)) {
      rowsum(A$obs[k$has_obs], gk[k$has_obs])[, 1]
    } else {
      numeric(0)
    }
    first <- match(ug, gk)
    v <- tibble::tibble(
      month = k$month[first],
      med = q[, 1],
      lo = q[, 2],
      hi = q[, 3],
      obs = unname(dplyr::coalesce(ob[ug], 0))
    )
    v <- if (grain == "region") {
      dplyr::mutate(v, scope = run$meta$group)
    } else {
      dplyr::mutate(v, scope = k$country[first], .before = 1)
    }
    out[[grain]] <- dplyr::mutate(v, grain = grain)
  }
  dplyr::bind_rows(out) |>
    dplyr::mutate(
      group = run$meta$group,
      launch_year = run$meta$launch_year,
      support = support,
      .before = 1
    )
}

#' score_case_references
#'
#' Reference scores for the per-run totals (bound across runs), per reference x
#' rule, under the same estimator as the model (fair CRPS / pinball):
#'   loo_climatology  other launch years' observed totals, same group / support
#'                    / grain / scope / lead class
#'   zero             point mass at 0
#'   persistence      previous launch year's observed total (NA if absent)
score_case_references <- function(totals, references, taus) {
  gk <- c("group", "support", "grain", "scope", "lead_class")
  t0 <- dplyr::select(
    totals,
    group,
    launch_year,
    support,
    grain,
    scope,
    lead_class,
    obs
  )
  t0 <- t0[
    do.call(order, c(unname(as.list(t0[gk])), method = "radix")),
    ,
    drop = FALSE
  ]
  g <- match(
    do.call(paste, c(unname(as.list(t0[gk])), sep = "\r")),
    unique(do.call(paste, c(unname(as.list(t0[gk])), sep = "\r")))
  )
  rows_of <- split(seq_len(nrow(t0)), g)
  out <- list()
  for (ref in references) {
    ens <- lapply(seq_len(nrow(t0)), function(i) {
      d <- rows_of[[g[i]]]
      switch(
        ref,
        loo_climatology = t0$obs[setdiff(d, i)],
        zero = 0,
        persistence = t0$obs[d[t0$launch_year[d] == t0$launch_year[i] - 1L]],
        stop("score_case_references: unknown reference '", ref, "'")
      )
    })
    sc <- score_sorted(
      sort_rows(pad_rows(ens)),
      t0$obs,
      c("crps_raw", "crps_log1p", "pinball"),
      taus = taus
    )
    out[[ref]] <- dplyr::bind_cols(
      dplyr::select(t0, dplyr::all_of(gk), launch_year),
      tibble::as_tibble(sc)
    ) |>
      tidyr::pivot_longer(
        c(crps_raw, crps_log1p, pinball),
        names_to = "rule",
        values_to = "score_ref"
      ) |>
      dplyr::mutate(reference = ref)
  }
  dplyr::bind_rows(out)
}

# ===========================================================================
# 5b. Monthly and season case skill (fixed-support references)
# ===========================================================================
# On when results.derivation.hindcast.case_skill_supports is non-empty.
# Units per run and support: month = (area, month) in the country window
# (grain adm1 or country); season = (country, lead class, "all") total
# ("region_sum" adds the group's countries). Under case_mask "coverage" every
# unit sums only surveilled (ADM1, month) cells, for model, observation and
# references alike. References use the observed history on the scored run's
# fixed support, scaled to the scored year's population; registry
# .CASE_SKILL_REFERENCES (clim_pooled, clim_launch_state, persistence_season,
# persistence_last). A season reference with fewer than case_season_min_ref
# values is NA.

#' case_obs_history
#'
#' The observed history the case references draw on, built once per 09 call
#' and shared by every run: monthly observed cases per ADM1 region key
#' (finest-level dedup, as obs_sim_bases()) and per gid (exact_gid), their
#' surveillance coverage (coverage mask cells, plus any month with an observed
#' case row) and populations per region key / gid and calendar year.
#' @param mon_all build_mon_obs_frame() frame: the whole history, not cut to a
#'   window
#' @param obs_panel hindcast_obs_panel() frame (only its months: the axis)
#' @param coverage tibble(gid, month) from build_coverage_mask(); required by
#'   case_mask "coverage"
#' @param wpp load_wpp_table() table; required by case_pop_scale "wpp"
case_obs_history <- function(
  mon_all,
  obs_panel,
  gid_pop,
  cfg,
  coverage = NULL,
  wpp = NULL,
  tag = "[hindcast-skill]"
) {
  mask <- as.character(cfg[["case_mask"]] %||% "none")
  scale <- as.character(cfg[["case_pop_scale"]] %||% "none")
  if (!mask %in% c("coverage", "none")) {
    stop(
      tag,
      " results.derivation.hindcast.case_mask must be 'coverage' or 'none'.",
      call. = FALSE
    )
  }
  if (!scale %in% c("wpp", "none")) {
    stop(
      tag,
      " results.derivation.hindcast.case_pop_scale must be 'wpp' or 'none'.",
      call. = FALSE
    )
  }
  if (mask == "coverage" && is.null(coverage)) {
    stop(
      tag,
      " case_mask 'coverage' needs the surveillance coverage mask (build_coverage_mask(); ",
      "the aggregated pre-outbreak series are missing).",
      call. = FALSE
    )
  }
  if (scale == "wpp" && is.null(wpp)) {
    stop(
      tag,
      " case_pop_scale 'wpp' needs the WPP2024 national totals (make_pop_run_file(opt, ",
      "'pop_wpp')); refusing to fall back to the nearest gid-table year.",
      call. = FALSE
    )
  }
  rng <- range(c(mon_all$month, obs_panel$month))
  m0 <- .month_index(rng[1])
  T <- .month_index(rng[2]) - m0 + 1L
  mi <- function(d) .month_index(d) - m0 + 1L
  rows <- list(
    rk = mon_all |>
      dplyr::filter(!is.na(region_key), !is.na(lvl)) |>
      dplyr::group_by(region_key, month) |>
      dplyr::filter(lvl == max(lvl)) |>
      dplyr::summarise(
        cases = sum(cases_month, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::rename(id = region_key),
    gid = mon_all |>
      dplyr::group_by(id = gid, month) |>
      dplyr::summarise(cases = sum(cases_month, na.rm = TRUE), .groups = "drop")
  )
  hist <- list(m0 = m0, T = T, mask = mask, scale = scale)
  for (lv in names(rows)) {
    r <- rows[[lv]]
    ids <- sort(unique(r$id))
    H <- matrix(0, length(ids), T, dimnames = list(ids, NULL))
    H[cbind(match(r$id, ids), mi(r$month))] <- r$cases
    hist[[lv]] <- list(ids = ids, H = H)
  }
  if (mask == "coverage") {
    cv <- dplyr::distinct(dplyr::transmute(
      coverage,
      gid = as.character(gid),
      month = as.Date(month)
    ))
    cv <- cv[cv$month >= rng[1] & cv$month <= rng[2], , drop = FALSE]
    cell <- function(ids, id_of) {
      # coverage per history id, plus observed rows
      C <- matrix(FALSE, length(ids), T, dimnames = list(ids, NULL))
      i <- match(id_of(cv$gid), ids)
      ok <- !is.na(i)
      C[cbind(i[ok], mi(cv$month[ok]))] <- TRUE
      C
    }
    rk_of <- function(g) gid_to_adm1(g, strict = FALSE)
    hist$rk$C <- cell(
      sort(unique(stats::na.omit(c(hist$rk$ids, rk_of(cv$gid))))),
      rk_of
    )
    hist$gid$C <- cell(sort(unique(c(hist$gid$ids, cv$gid))), identity)
    for (lv in c("rk", "gid")) {
      o <- hist[[lv]]$H > 0
      hist[[lv]]$C[rownames(o), ] <- hist[[lv]]$C[rownames(o), , drop = FALSE] |
        o
    }
  }
  if (scale == "wpp") {
    yrs <- seq(
      rng[1] |> format("%Y") |> as.integer(),
      max(as.integer(format(rng[2], "%Y")), 2030L)
    )
    # the gids needed plus every gid of their ADM1s (resolve_gid_population()'s
    # parent-residual tier reads the siblings)
    gp <- gid_pop[
      gid_to_adm1(gid_pop$gid, strict = FALSE) %in%
        c(hist$rk$ids, region_key_of(hist$gid$ids)),
      ,
      drop = FALSE
    ]
    pop_of <- function(g) {
      if (!length(g)) {
        return(matrix(numeric(0), 0, length(yrs)))
      }
      p <- resolve_gid_population(
        rep(g, each = length(yrs)),
        rep(yrs, length(g)),
        gp,
        wpp = wpp,
        on_missing = "na"
      )$pop
      matrix(p, length(g), length(yrs), byrow = TRUE, dimnames = list(g, NULL))
    }
    a1 <- unique(gp$gid[
      gp$gadm_level == "ADM1" &
        gid_to_adm1(gp$gid, strict = FALSE) %in% hist$rk$ids
    ])
    P1 <- pop_of(a1)
    hist$rk$P <- matrix(
      NA_real_,
      length(hist$rk$ids),
      length(yrs),
      dimnames = list(hist$rk$ids, NULL)
    )
    if (length(a1)) {
      s <- rowsum(P1, region_key_of(a1))
      hist$rk$P[rownames(s), ] <- s
    }
    hist$gid$P <- pop_of(hist$gid$ids)
    hist$yr0 <- yrs[1]
  }
  hist
}

# area-level history of one run x support: areas = the support's ADM1 region
# keys; H (cases on the fixed support: ADM1 cases x the scored year's population
# share, or the exact gids' own cases), C (surveilled; NULL = no mask), P
# (population for the scaling; NULL = none), act (observed cases > 0)
.case_area_history <- function(A, support, hist) {
  areas <- sort(unique(A$keys$region_key))
  T <- hist$T
  pick <- function(X, ids, fill) {
    out <- matrix(fill, length(ids), ncol(X))
    i <- match(ids, rownames(X))
    out[!is.na(i), ] <- X[i[!is.na(i)], , drop = FALSE]
    out
  }
  if (support == "exact_gid") {
    g <- intersect(A$gids, hist$gid$ids)
    agg <- function(X, fn) {
      out <- matrix(if (is.logical(X)) FALSE else 0, length(areas), ncol(X))
      if (length(g)) {
        s <- rowsum(X[g, , drop = FALSE] + 0, region_key_of(g))
        i <- match(rownames(s), areas)
        out[i[!is.na(i)], ] <- fn(s[!is.na(i), , drop = FALSE])
      }
      out
    }
    H <- agg(hist$gid$H, identity)
    C <- if (!is.null(hist$gid$C)) {
      gc <- intersect(A$gids, rownames(hist$gid$C))
      out <- matrix(FALSE, length(areas), T)
      if (length(gc)) {
        s <- rowsum(hist$gid$C[gc, , drop = FALSE] + 0, region_key_of(gc))
        i <- match(rownames(s), areas)
        out[i[!is.na(i)], ] <- s[!is.na(i), , drop = FALSE] > 0
      }
      out
    }
    P <- if (!is.null(hist$gid$P)) agg(hist$gid$P, identity)
  } else {
    s <- A$shares
    sh <- dplyr::coalesce(s$s_g[match(areas, s$region_key)], 1)
    H <- pick(hist$rk$H, areas, 0) * sh
    C <- if (!is.null(hist$rk$C)) pick(hist$rk$C, areas, FALSE)
    P <- if (!is.null(hist$rk$P)) pick(hist$rk$P, areas, NA_real_)
  }
  list(areas = areas, H = H, C = C, P = P, act = H > 0)
}

# Reference values of every unit: cu = one row per (unit, cell) with the cell's
# area index ai, month index mi (history axis), year index yi, leave-out bounds
# e_lo / e_hi. Either `offsets` (s, k: shifted copies of the whole cell set,
# valid only if every shifted cell is in range, surveilled and outside the
# leave-out) or `s_cell` (one shifted cell per row: a point forecast, no
# leave-out). Values are scaled to the scored year's population. Returns
# tibble(unit, k, value).
.case_ref_values <- function(
  cu,
  ah,
  hist,
  n_units,
  offsets = NULL,
  s_cell = NULL
) {
  nc <- nrow(cu)
  if (is.null(offsets)) {
    ci <- seq_len(nc)
    oi <- rep(1L, nc)
    sm <- cu$mi + s_cell
    ks <- 0L
  } else {
    no <- nrow(offsets)
    ci <- rep(seq_len(nc), times = no)
    oi <- rep(seq_len(no), each = nc)
    sm <- cu$mi[ci] + offsets$s[oi]
    ks <- offsets$k
  }
  ai <- cu$ai[ci]
  ok <- sm >= 1L & sm <= hist$T
  v <- numeric(length(sm))
  v[ok] <- ah$H[cbind(ai[ok], sm[ok])]
  if (!is.null(ah$C)) {
    ok[ok] <- ah$C[cbind(ai[ok], sm[ok])]
  }
  if (!is.null(offsets)) {
    ok <- ok & !(sm >= cu$e_lo[ci] & sm <= cu$e_hi[ci])
  }
  if (!is.null(ah$P)) {
    ys <- (hist$m0 + sm - 1L) %/% 12L - hist$yr0 + 1L
    ys[!ok | ys < 1L | ys > ncol(ah$P)] <- NA
    r <- ah$P[cbind(ai, cu$yi[ci])] / ah$P[cbind(ai, ys)]
    v <- ifelse(v == 0, 0, v * r)
    ok <- ok & !is.na(v)
  }
  g <- cu$unit[ci] + (oi - 1L) * n_units
  agg <- rowsum(cbind(ifelse(ok, v, 0), !ok + 0), g)
  gg <- as.numeric(rownames(agg))[agg[, 2] == 0]
  tibble::tibble(
    unit = as.integer((gg - 1) %% n_units + 1),
    k = ks[(gg - 1) %/% n_units + 1],
    value = agg[agg[, 2] == 0, 1]
  )
}

# activity of every unit in the month before launch shifted by k years (k = 0:
# the scored year): any of its cells' areas with observed cases; NA when out of
# the history. Returns a units x length(ks) logical matrix.
.case_unit_state <- function(cu, ah, hist, n_units, ks) {
  matrix(
    vapply(
      ks,
      function(k) {
        m <- cu$L - 1L + 12L * k
        inr <- m >= 1L & m <= hist$T
        a <- rep(NA, nrow(cu))
        a[inr] <- ah$act[cbind(cu$ai[inr], m[inr])]
        s <- rowsum(cbind(a %in% TRUE, is.na(a)) + 0, cu$unit)
        out <- rep(NA, n_units)
        u <- as.integer(rownames(s))
        out[u] <- ifelse(s[, 2] > 0, NA, s[, 1] > 0)
        out
      },
      logical(n_units)
    ),
    nrow = n_units
  )
}

# reference registry: name -> list(kind ("clim" | "point"), month, season)
.CASE_SKILL_REFERENCES <- function() {
  list(
    clim_pooled = list(
      kind = "clim",
      launch_state = FALSE,
      month = TRUE,
      season = TRUE
    ),
    clim_launch_state = list(
      kind = "clim",
      launch_state = TRUE,
      month = TRUE,
      season = TRUE
    ),
    persistence_season = list(
      kind = "point",
      source = "season",
      month = TRUE,
      season = TRUE
    ),
    persistence_last = list(
      kind = "point",
      source = "last",
      month = TRUE,
      season = FALSE
    )
  )
}

# Score one unit set (month or season units) of a run x support: model,
# observation, every reference. `cu` maps units to cells (rows of A); `uk` is
# the unit table (one row per unit, in unit order).
.case_score_units <- function(uk, cu, A, ah, hist, cfg, type) {
  U <- nrow(uk)
  if (!U) {
    return(NULL)
  }
  rules <- as.character(unlist(
    cfg[["case_skill_rules"]] %||% c("crps_raw", "crps_log1p", "twcrps")
  ))
  tw <- as.numeric(
    cfg[[
      if (type == "month") {
        "case_tw_threshold_month"
      } else {
        "case_tw_threshold_season"
      }
    ]] %||%
      50
  )
  S <- sort_rows(rowsum(A$M[cu$k, , drop = FALSE], cu$unit, reorder = TRUE))
  y <- unname(rowsum(A$obs[cu$k], cu$unit, reorder = TRUE)[, 1])
  n <- rep(ncol(S), U)
  q <- quantile7_sorted(S, c(.5, .025, .975), n)
  out <- dplyr::mutate(
    uk,
    n_cells = tabulate(cu$unit, U),
    obs = y,
    mean = rowSums(S) / n,
    med = q[, 1],
    q025 = q[, 2],
    q975 = q[, 3],
    pctile = 100 * rowSums(S <= y) / n
  )
  sm <- score_sorted(S, y, rules, n, tw = tw)
  for (ru in rules) {
    out[[paste0("model__", ru)]] <- sm[[ru]]
  }
  h <- as.integer(cfg[["case_clim_halfwidth"]] %||% 1L)
  Ky <- ceiling(hist$T / 12) + 1L
  kk <- setdiff(-Ky:Ky, 0L)
  offs <- if (type == "month") {
    tidyr::expand_grid(k = kk, d = -h:h)
  } else {
    tibble::tibble(k = kk, d = 0L)
  }
  offs$s <- 12L * offs$k + offs$d
  min_ref <- if (type == "season") {
    as.integer(cfg[["case_season_min_ref"]] %||% 5L)
  } else {
    1L
  }
  reg <- .CASE_SKILL_REFERENCES()
  refs <- as.character(unlist(cfg[["case_skill_references"]] %||% names(reg)))
  clim <- NULL
  state <- NULL
  for (ref in refs) {
    spec <- reg[[ref]]
    if (is.null(spec)) {
      stop("case skill: unknown reference '", ref, "'", call. = FALSE)
    }
    if (!isTRUE(spec[[type]])) {
      next
    }
    vals <- if (spec$kind == "clim") {
      if (is.null(clim)) {
        clim <- .case_ref_values(cu, ah, hist, U, offsets = offs)
      }
      v <- clim
      if (isTRUE(spec$launch_state)) {
        if (is.null(state)) {
          state <- .case_unit_state(cu, ah, hist, U, c(0L, kk))
        }
        st <- state[cbind(v$unit, match(v$k, c(0L, kk)))]
        v <- v[
          !is.na(st) & !is.na(state[v$unit, 1]) & st == state[v$unit, 1],
          ,
          drop = FALSE
        ]
      }
      v
    } else {
      s_cell <- if (spec$source == "last") {
        (cu$L - 1L) - cu$mi
      } else {
        -12L * ceiling((cu$mi - cu$L + 1L) / 12)
      }
      .case_ref_values(cu, ah, hist, U, s_cell = s_cell)
    }
    R <- rows_by_unit(vals$unit, vals$value, U)
    nr <- rowSums(!is.na(R))
    # the season minimum applies to climatologies, not to point forecasts
    if (spec$kind == "clim" && min_ref > 1L) {
      R[nr < min_ref, ] <- NA
      nr[nr < min_ref] <- 0L
    }
    sr <- score_sorted(sort_rows(R), y, rules, nr, tw = tw)
    for (ru in rules) {
      out[[paste0(ref, "__", ru)]] <- sr[[ru]]
    }
    out[[paste0(ref, "__n")]] <- as.integer(nr)
  }
  out
}

#' case_skill_units
#'
#' Monthly and season case units of one run, every support in `arrays`
#' (case_unit_arrays()), scored for the model and every reference
#' (section 5b header). Wide tables: model__<rule>, <reference>__<rule>,
#' <reference>__n.
#' @return list(month, season)
case_skill_units <- function(run, arrays, hist, cfg, lead_breaks) {
  cw <- run$country_windows
  grains <- as.character(unlist(
    cfg[["case_month_grains"]] %||% c("country", "adm1")
  ))
  mi_of <- function(d) .month_index(d) - hist$m0 + 1L
  out_m <- out_s <- list()
  for (sp in names(arrays)) {
    A <- arrays[[sp]]
    if (!nrow(A$keys)) {
      next
    }
    ah <- .case_area_history(A, sp, hist)
    kt <- A$keys
    w <- match(kt$country, cw$iso3)
    L <- mi_of(cw$launch[w])
    F <- mi_of(cw$firm_horizon_end[w])
    K <- (F - L) %/% 12L
    cells <- tibble::tibble(
      k = seq_len(nrow(kt)),
      ai = match(kt$region_key, ah$areas),
      mi = mi_of(kt$month),
      L = L,
      e_lo = L - 12L * K,
      e_hi = F + 12L * K,
      country = kt$country,
      region_key = kt$region_key,
      month = kt$month,
      lead = kt$lead,
      lead_class = kt$lead_class
    )
    cells <- cells[
      !is.na(w) & cells$mi >= cells$L & cells$mi <= F,
      ,
      drop = FALSE
    ]
    if (!is.null(ah$C)) {
      inr <- cells$mi >= 1L & cells$mi <= hist$T
      cov <- rep(FALSE, nrow(cells))
      cov[inr] <- ah$C[cbind(cells$ai[inr], cells$mi[inr])]
      cells <- cells[cov, , drop = FALSE]
    }
    if (!nrow(cells)) {
      next
    }
    cells$yi <- if (!is.null(ah$P)) {
      (hist$m0 + cells$mi - 1L) %/% 12L - hist$yr0 + 1L
    } else {
      NA_integer_
    }
    tag_units <- function(cu, keys) {
      # cu: cells with a `ukey`; keys: unit columns
      uk <- dplyr::distinct(cu[keys])
      uk <- uk[
        do.call(order, c(unname(as.list(uk)), method = "radix")),
        ,
        drop = FALSE
      ]
      cu$unit <- match(
        do.call(paste, c(unname(as.list(cu[keys])), sep = "\r")),
        do.call(paste, c(unname(as.list(uk)), sep = "\r"))
      )
      list(uk = uk, cu = cu)
    }
    key <- function(df) {
      dplyr::mutate(
        df,
        group = run$meta$group,
        launch_year = run$meta$launch_year,
        support = sp,
        .before = 1
      )
    }
    # month units
    for (gr in grains) {
      cu <- dplyr::mutate(
        cells,
        scope = if (gr == "adm1") region_key else country
      )
      tu <- tag_units(cu, c("scope", "country", "month"))
      tu$uk <- dplyr::left_join(
        tu$uk,
        dplyr::distinct(cells, country, month, lead, lead_class),
        by = c("country", "month")
      )
      out_m[[length(out_m) + 1L]] <- key(dplyr::mutate(
        .case_score_units(
          dplyr::mutate(tu$uk, grain = gr, .before = 1),
          tu$cu,
          A,
          ah,
          hist,
          cfg,
          "month"
        ),
        area = scope
      ))
    }
    # season units: each lead class and "all", per country and the group's sum
    cs <- .with_all_lead(cells)
    for (gr in c("country", "region_sum")) {
      cu <- dplyr::mutate(
        cs,
        scope = if (gr == "country") country else run$meta$group
      )
      tu <- tag_units(cu, c("scope", "lead_class"))
      out_s[[length(out_s) + 1L]] <- key(dplyr::mutate(
        .case_score_units(
          dplyr::mutate(tu$uk, grain = gr, .before = 1),
          tu$cu,
          A,
          ah,
          hist,
          cfg,
          "season"
        ),
        n_months = as.integer(tapply(tu$cu$mi, tu$cu$unit, function(x) {
          length(unique(x))
        })),
        country = if (gr == "country") scope else NA_character_,
        area = scope
      ))
    }
  }
  list(month = dplyr::bind_rows(out_m), season = dplyr::bind_rows(out_s))
}

# long (reference, rule, score_model, score_ref) form of a wide unit table
.case_skill_long <- function(units, refs, rules) {
  dplyr::bind_rows(lapply(refs, function(ref) {
    dplyr::bind_rows(lapply(rules, function(ru) {
      col <- paste0(ref, "__", ru)
      if (!col %in% names(units)) {
        return(NULL)
      }
      dplyr::mutate(
        units[, !grepl("__", names(units)), drop = FALSE],
        reference = ref,
        rule = ru,
        score_model = units[[paste0("model__", ru)]],
        score_ref = units[[col]]
      )
    }))
  }))
}

.hindcast_case_skill_rollup <- function(
  month_units,
  season_units,
  cfg,
  boot,
  lead_breaks
) {
  refs <- as.character(unlist(
    cfg[["case_skill_references"]] %||% names(.CASE_SKILL_REFERENCES())
  ))
  rules <- as.character(unlist(
    cfg[["case_skill_rules"]] %||% c("crps_raw", "crps_log1p", "twcrps")
  ))
  area <- if (isTRUE(cfg[["skill_area_summaries"]])) "area"
  by <- c(
    "scope_level",
    "scope",
    "support",
    "grain",
    "lead_class",
    "reference",
    "rule"
  )
  levels_of <- function(df, country = TRUE) {
    dplyr::bind_rows(
      dplyr::mutate(df, scope_level = "region", scope = group),
      dplyr::mutate(df, scope_level = "africa", scope = "africa"),
      if (country) {
        dplyr::mutate(
          df[!is.na(df$country), , drop = FALSE],
          scope_level = "country",
          scope = country
        )
      }
    )
  }
  roll <- function(long) {
    if (nrow(long)) {
      .skill_rollup(long, by, .skill_stats_units(), boot, area)
    }
  }
  ml <- .case_skill_long(month_units, refs, rules)
  sl <- .case_skill_long(season_units, refs, rules)
  list(
    month_units = month_units,
    month_skill = roll(levels_of(.with_all_lead(ml))),
    season_units = season_units,
    season_skill = roll(dplyr::bind_rows(
      levels_of(sl[sl$grain == "country", , drop = FALSE]),
      levels_of(sl[sl$grain == "region_sum", , drop = FALSE], country = FALSE)
    )),
    skill_meta = list(
      supports = as.character(unlist(cfg[["case_skill_supports"]])),
      month_grains = as.character(unlist(
        cfg[["case_month_grains"]] %||% c("country", "adm1")
      )),
      references = refs,
      rules = rules,
      lead_breaks = lead_breaks,
      lead_levels = hindcast_lead_levels(lead_breaks),
      clim_halfwidth = as.integer(cfg[["case_clim_halfwidth"]] %||% 1L),
      tw_threshold_month = as.numeric(cfg[["case_tw_threshold_month"]] %||% 50),
      tw_threshold_season = as.numeric(
        cfg[["case_tw_threshold_season"]] %||% 50
      ),
      season_min_ref = as.integer(cfg[["case_season_min_ref"]] %||% 5L),
      mask = as.character(cfg[["case_mask"]] %||% "none"),
      pop_scale = as.character(cfg[["case_pop_scale"]] %||% "none"),
      bootstrap = if (!is.null(boot)) {
        list(n = nrow(boot$counts), level = boot$level, years = boot$years)
      },
      estimator = paste(
        "fair CRPS (crps_raw, crps_log1p) and fair threshold-weighted CRPS",
        "(twcrps) for model and references; pooled CRPSS = 1 - sum(model) /",
        "sum(reference) over rows where both are finite"
      )
    )
  )
}

# ===========================================================================
# 6. Orchestration: one pass over the runs, then references and roll-ups
# ===========================================================================

#' hindcast_quantile_probs -- ensemble quantiles kept per occurrence unit: the
#' 50/95 % bands plus both ends of every interval width (interval score).
hindcast_quantile_probs <- function(widths) {
  sort(unique(round(
    c(.025, .25, .5, .75, .975, (1 - widths) / 2, 1 - (1 - widths) / 2),
    6
  )))
}

# region rows + an "africa" pool, for any frame with a `group` column
.by_scope <- function(df) {
  dplyr::bind_rows(
    dplyr::mutate(df, scope_level = "region", scope = group),
    dplyr::mutate(df, scope_level = "africa", scope = "africa")
  )
}

#' hindcast_skill_bundle
#'
#' Score every hindcast run on both tracks. Reads each run once (the 44
#' monthly summaries are too large to hold together), keeps only the small
#' per-unit tables, then scores the references and rolls up.
#'
#' @param idx hindcast_run_index() tibble
#' @param obs_panel hindcast_obs_panel() frame
#' @param mon_all build_mon_obs_frame() frame
#' @param gid_pop load_gid_population() frame
#' @param cfg results.derivation.hindcast list (lead_breaks, occ_grains,
#'   occ_references, occ_rules, case_supports, case_references, case_rules,
#'   pinball_taus, interval_widths, reliability_bins, zero_floor, pit_seed,
#'   rank_seed)
#' @param workers forked processes for the per-run scoring (map_workers());
#'   the result is identical for any value
#' @param coverage tibble(gid, month) of surveilled cells (build_coverage_mask());
#'   required when case_skill_supports is set and case_mask is "coverage"
#' @param wpp load_wpp_table() table; required when case_skill_supports is set and
#'   case_pop_scale is "wpp"
#' @return list(runs, occ, cases) -- the three 09 artefacts
hindcast_skill_bundle <- function(
  idx,
  obs_panel,
  mon_all,
  gid_pop,
  cfg,
  tag = "[hindcast-skill]",
  workers = 1L,
  coverage = NULL,
  wpp = NULL
) {
  lead_breaks <- as.integer(unlist(cfg$lead_breaks))
  widths <- as.numeric(unlist(cfg$interval_widths))
  taus <- as.numeric(unlist(cfg$pinball_taus))
  probs <- hindcast_quantile_probs(widths)
  # Monthly / season case skill (off when case_skill_supports is empty): the
  # observed history its references draw on is built once and shared by every
  # run (forked workers read it copy-on-write).
  skill_sups <- as.character(unlist(cfg[["case_skill_supports"]]))
  hist <- if (length(skill_sups)) {
    case_obs_history(
      mon_all,
      obs_panel,
      gid_pop,
      cfg,
      coverage = coverage,
      wpp = wpp,
      tag = tag
    )
  }
  set.seed(cfg$rank_seed %||% 42L)

  # Runs are scored independently, on `workers` forked processes (map_workers(),
  # base.R). The only random draw in a run's scoring is the rank tie-break; it
  # is deferred (score_occurrence_model(defer_rank = TRUE)) and made below, run
  # by run in index order, so the stream after set.seed() above is consumed
  # exactly as by a serial loop and every artefact is the same for any worker
  # count.
  # One scoring entry per (fit group, launch year). Under a v2 manifest a fit
  # group whose countries span several simulation units has one index row per
  # unit; those rows are pooled into one run, so every rollup keyed by
  # (group, launch_year) sees the whole group once. Entries keep the index order
  # of their first row, so a v1 index scores exactly as before.
  entry_key <- paste(idx$group, idx$launch_year)
  entries <- split(
    seq_len(nrow(idx)),
    factor(entry_key, levels = unique(entry_key))
  )
  per <- map_workers(
    seq_along(entries),
    function(i) {
      rows <- idx[entries[[i]], ]
      message(sprintf(
        "%s scoring %s/%d (%d of %d)%s",
        tag,
        rows$group[1L],
        rows$launch_year[1L],
        i,
        length(entries),
        if (nrow(rows) > 1L) {
          sprintf(" [%d unit runs pooled]", nrow(rows))
        } else {
          ""
        }
      ))
      run <- hindcast_load_run(rows, obs_panel, mon_all, gid_pop, lead_breaks)
      key <- function(df) {
        dplyr::mutate(
          df,
          group = run$meta$group,
          launch_year = run$meta$launch_year,
          .before = 1
        )
      }
      # the replicate rows are read once, into (cell x replicate) matrices; the
      # row frames are then dropped (every scorer below reads the matrices)
      oi <- sim_cells(run$sim, "Y")
      si <- case_sim_index(run)
      run$sim <- NULL
      run$bases$sim_in <- NULL
      units <- dplyr::bind_rows(lapply(unlist(cfg$occ_grains), function(g) {
        dplyr::mutate(
          score_occurrence_model(
            onset_window_counts(run, g, sim_index = oi),
            probs,
            defer_rank = TRUE
          ),
          grain = g
        )
      }))
      rm(oi)
      # one aggregation per support serves the legacy totals, the views and the
      # monthly / season units
      sups <- unique(c(unlist(cfg$case_supports), skill_sups))
      arrays <- stats::setNames(
        lapply(sups, function(sp) {
          case_unit_arrays(
            run,
            sp,
            mon_all,
            gid_pop,
            lead_breaks,
            sim_index = si
          )
        }),
        sups
      )
      list(
        meta = run$meta,
        cells = key(run$cells),
        units = key(units),
        cases = score_cases_run(
          run,
          unlist(cfg$case_supports),
          mon_all,
          gid_pop,
          lead_breaks,
          taus,
          arrays = arrays
        ),
        monthly = dplyr::bind_rows(lapply(
          unlist(cfg$case_supports),
          function(sp) {
            case_monthly_views(run, sp, mon_all, gid_pop, arrays = arrays[[sp]])
          }
        )),
        skill_units = if (length(skill_sups)) {
          case_skill_units(run, arrays[skill_sups], hist, cfg, lead_breaks)
        },
        freq = if (!is.null(run$monthly_freq)) key(run$monthly_freq) else NULL
      )
    },
    workers = workers
  )
  for (i in seq_along(per)) {
    per[[i]]$units <- fill_randomised_ranks(per[[i]]$units)
  }
  pick <- function(k) dplyr::bind_rows(lapply(per, `[[`, k))
  runs <- pick("meta")
  cells <- pick("cells")
  units <- pick("units")
  # Season-window index: a (group, launch year) spans runs with different
  # launches; its window is the pooled run's (earliest launch .. group horizon).
  windows <- if ("country_launch" %in% names(idx)) {
    dplyr::distinct(runs, group, launch_year, launch, horizon_end)
  } else {
    dplyr::distinct(idx, group, launch_year, launch, horizon_end)
  }
  pool <- dplyr::distinct(cells, group, gid, month, y = y_obs)
  # one season-year resample set for every CRPSS of both tracks (NULL = no CIs)
  boot <- hindcast_boot_spec(runs$launch_year, cfg)

  cases <- .hindcast_cases_rollup(
    pick("cases"),
    pick("monthly"),
    cfg,
    taus,
    lead_breaks,
    boot
  )
  if (length(skill_sups)) {
    su <- lapply(per, `[[`, "skill_units")
    cases <- c(
      cases,
      .hindcast_case_skill_rollup(
        dplyr::bind_rows(lapply(su, `[[`, "month")),
        dplyr::bind_rows(lapply(su, `[[`, "season")),
        cfg,
        boot,
        lead_breaks
      )
    )
  }
  list(
    runs = runs,
    occ = .hindcast_occ_rollup(
      units,
      cells,
      pick("freq"),
      pool,
      obs_panel,
      windows,
      runs,
      cfg,
      probs,
      widths,
      lead_breaks,
      boot
    ),
    cases = cases
  )
}

#' hindcast_boot_spec
#'
#' The season-year block bootstrap shared by every CRPSS (cfg bootstrap_n > 0):
#' B multinomial resample counts over the scored season years, one draw for all
#' strata and both tracks, so strata are resampled jointly. The draws use their
#' own seed (bootstrap_seed) and the global RNG state is restored, so they never
#' shift another draw. NULL when bootstrap_n is 0 or absent.
hindcast_boot_spec <- function(years, cfg) {
  B <- as.integer(cfg[["bootstrap_n"]] %||% 0L)
  if (is.na(B) || B <= 0L) {
    return(NULL)
  }
  years <- sort(unique(years))
  had <- exists(".Random.seed", envir = globalenv())
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else {
      rm(".Random.seed", envir = globalenv())
    }
  )
  set.seed(as.integer(cfg[["bootstrap_seed"]] %||% 20261101L))
  list(
    years = years,
    counts = t(stats::rmultinom(B, length(years), rep(1, length(years)))),
    level = as.numeric(cfg[["bootstrap_level"]] %||% 0.95)
  )
}

#' .skill_rollup
#'
#' The one skill summary of every hindcast score table (occurrence, legacy case
#' totals, monthly units, season totals). `stats` are named expressions over a
#' stratum's rows (rlang::quos(); score_model / score_ref), evaluated per `by`
#' stratum. Optional extras:
#'   boot  (hindcast_boot_spec()) season-year bootstrap of the pooled CRPSS
#'         1 - sum(model) / sum(ref) over the rows where both are finite, and of
#'         the mean model-minus-reference difference: crpss_lo / crpss_hi,
#'         diff_mean / diff_lo / diff_hi
#'   area  column naming the scored area (gid, country, ADM1): the median of
#'         the per-area CRPSS and the share of areas with positive skill
#'         (n_areas / median_area_crpss / share_area_skill)
.skill_rollup <- function(df, by, stats, boot = NULL, area = NULL) {
  out <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(!!!stats, .groups = "drop")
  if (!is.null(boot) && nrow(out)) {
    out <- dplyr::left_join(out, .skill_boot(df, by, boot), by = by)
  }
  if (!is.null(area) && nrow(out)) {
    out <- dplyr::left_join(out, .skill_area(df, by, area), by = by)
  }
  out
}

# the per-stratum statistics of the occurrence / monthly / season tables
# (n_runs = distinct runs) and of the legacy case totals (one row per run)
.skill_stats_units <- function() {
  rlang::quos(
    n_units = dplyr::n(),
    n_scored = sum(is.finite(score_ref)),
    n_runs = dplyr::n_distinct(group, launch_year),
    mean_model = mean(score_model),
    mean_ref = mean(score_ref, na.rm = TRUE),
    crpss = skill_ratio(score_model, score_ref)
  )
}
.skill_stats_runs <- function() {
  rlang::quos(
    n_runs = dplyr::n(),
    mean_model = mean(score_model),
    median_model = stats::median(score_model),
    mean_ref = mean(score_ref, na.rm = TRUE),
    crpss = skill_ratio(score_model, score_ref)
  )
}

.skill_boot <- function(df, by, boot) {
  d <- df[is.finite(df$score_model) & is.finite(df$score_ref), , drop = FALSE]
  if (!nrow(d)) {
    return(dplyr::distinct(df[0, by, drop = FALSE]))
  }
  gd <- dplyr::group_by(d, dplyr::across(dplyr::all_of(by)))
  st <- dplyr::group_indices(gd)
  keys <- dplyr::group_keys(gd)
  ns <- nrow(keys)
  ny <- length(boot$years)
  cell <- st + (match(d$launch_year, boot$years) - 1L) * ns
  agg <- rowsum(cbind(d$score_model, d$score_ref, 1), cell)
  Sm <- Sr <- N <- matrix(0, ns, ny)
  uc <- sort(unique(cell))
  Sm[uc] <- agg[, 1]
  Sr[uc] <- agg[, 2]
  N[uc] <- agg[, 3]
  Ct <- t(boot$counts)
  SmB <- Sm %*% Ct
  SrB <- Sr %*% Ct
  cr <- 1 - SmB / SrB
  cr[!is.finite(cr) | .ref_score_zero(SrB)] <- NA
  df_ <- (SmB - SrB) / (N %*% Ct)
  df_[!is.finite(df_)] <- NA
  a <- (1 - boot$level) / 2
  cq <- quantile7_sorted(sort_rows(cr), c(a, 1 - a))
  dq <- quantile7_sorted(sort_rows(df_), c(a, 1 - a))
  dplyr::bind_cols(
    keys,
    tibble::tibble(
      crpss_lo = cq[, 1],
      crpss_hi = cq[, 2],
      diff_mean = rowSums(Sm - Sr) / rowSums(N),
      diff_lo = dq[, 1],
      diff_hi = dq[, 2]
    )
  )
}

.skill_area <- function(df, by, area) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, area)))) |>
    dplyr::summarise(
      .a = skill_ratio(score_model, score_ref),
      .groups = "drop"
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_areas = sum(is.finite(.a)),
      median_area_crpss = if (any(is.finite(.a))) {
        stats::median(.a[is.finite(.a)])
      } else {
        NA_real_
      },
      share_area_skill = if (any(is.finite(.a))) {
        mean(.a[is.finite(.a)] > 0)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

.hindcast_occ_rollup <- function(
  units,
  cells,
  freq,
  pool,
  obs_panel,
  windows,
  runs,
  cfg,
  probs,
  widths,
  lead_breaks,
  boot = NULL
) {
  keys <- c("group", "launch_year", "grain", "gid", "unit", "lead_class")
  rules <- unlist(cfg$occ_rules)
  refs <- score_occurrence_references(
    units,
    unlist(cfg$occ_references),
    rules,
    pool,
    obs_panel,
    windows
  )
  model <- units |>
    dplyr::select(dplyr::all_of(keys), n_atrisk, y_obs, crps_emp, crps_fair) |>
    tidyr::pivot_longer(
      c(crps_emp, crps_fair),
      names_to = "rule",
      values_to = "score_model"
    ) |>
    dplyr::filter(rule %in% rules)
  scores <- dplyr::inner_join(
    model,
    refs,
    by = c(keys, "rule"),
    relationship = "one-to-many"
  )
  skill_by <- c(
    "scope_level",
    "scope",
    "grain",
    "lead_class",
    "reference",
    "rule"
  )
  area <- if (isTRUE(cfg[["skill_area_summaries"]])) "gid"
  skill <- dplyr::bind_rows(
    .skill_rollup(
      .by_scope(scores),
      skill_by,
      .skill_stats_units(),
      boot,
      area
    ),
    .skill_rollup(
      dplyr::mutate(scores, scope_level = "region_year", scope = group),
      c(skill_by, "launch_year"),
      .skill_stats_units()
    )
  )

  # monthly-cell scores per reference and lead class (Brier family)
  n_rep <- stats::median(runs$n_rep)
  eps <- if (isTRUE(cfg$zero_floor)) 1 / (2 * n_rep) else 1e-6
  cl <- .by_scope(.with_all_lead(cells))
  cell_scores <- dplyr::bind_rows(lapply(
    unlist(cfg$occ_references),
    function(ref) {
      cr <- occ_cell_reference_prob(cells, ref, pool, obs_panel, windows) |>
        dplyr::select(group, launch_year, gid, month, p_ref)
      cl |>
        dplyr::left_join(cr, by = c("group", "launch_year", "gid", "month")) |>
        dplyr::group_by(scope_level, scope, lead_class) |>
        dplyr::group_modify(~ score_cells(.x$p_hat, .x$y_obs, .x$p_ref, eps)) |>
        dplyr::ungroup() |>
        dplyr::mutate(reference = ref)
    }
  ))
  cl_all <- dplyr::filter(cl, lead_class == "all")
  per_scope <- function(f) {
    cl_all |>
      dplyr::group_by(scope_level, scope) |>
      dplyr::group_modify(~ f(.x)) |>
      dplyr::ungroup()
  }
  reliability <- per_scope(function(d) {
    reliability_curve(d$p_hat, d$y_obs, cfg$reliability_bins %||% 10L)
  })
  roc <- per_scope(function(d) roc_points(d$p_hat, d$y_obs))
  pit_hist <- per_scope(function(d) pit_hist_bins(d$pit))
  cov_widths <- c(.10, .20, .30, .40, .50, .60, .70, .80, .90, .95, .99)
  coverage_curve <- per_scope(function(d) {
    tibble::tibble(
      width = cov_widths,
      coverage = coverage_at_width(d$pit, cov_widths)
    )
  })

  # ensemble interval score + rank histogram per occurrence unit
  qn <- function(p) sprintf("q%03d", round(1000 * p))
  interval <- dplyr::bind_rows(lapply(widths, function(w) {
    lo <- qn((1 - w) / 2)
    hi <- qn(1 - (1 - w) / 2)
    dplyr::mutate(
      units,
      width = w,
      IS = interval_score(.data[[lo]], .data[[hi]], y_obs, w)
    )
  }))
  interval_score_tbl <- .by_scope(interval) |>
    dplyr::group_by(scope_level, scope, grain, lead_class, width) |>
    dplyr::summarise(IS = mean(IS), .groups = "drop")
  rank_hist <- .by_scope(units) |>
    dplyr::group_by(scope_level, scope, grain, lead_class) |>
    dplyr::group_modify(~ pit_hist_bins(.x$rank_u, 10L)) |>
    dplyr::ungroup()

  # seasonal cycle (07c's monthly_freq), as the notebook's monthly-rollup chunk
  season <- season_cov <- NULL
  if (!is.null(freq) && nrow(freq)) {
    season <- .by_scope(freq) |>
      dplyr::group_by(scope_level, scope, month_of_year) |>
      dplyr::summarise(
        freq_obs = mean(freq_obs, na.rm = TRUE),
        sim_mean = mean(freq_sim_mean, na.rm = TRUE),
        sim_lo = mean(freq_sim_q025, na.rm = TRUE),
        sim_hi = mean(freq_sim_q975, na.rm = TRUE),
        .groups = "drop"
      )
    inb <- function(o, lo, hi) o >= lo & o <= hi
    season_cov <- .by_scope(freq) |>
      dplyr::group_by(scope_level, scope) |>
      dplyr::summarise(
        cov50 = mean(inb(freq_obs, freq_sim_q25, freq_sim_q75), na.rm = TRUE),
        cov95 = mean(inb(freq_obs, freq_sim_q025, freq_sim_q975), na.rm = TRUE),
        .groups = "drop"
      )
  }

  list(
    gid_window = dplyr::select(units, -months),
    scores = scores,
    skill = skill,
    cell_scores = cell_scores,
    reliability = reliability,
    roc = roc,
    pit_hist = pit_hist,
    coverage_curve = coverage_curve,
    interval_score = interval_score_tbl,
    rank_hist = rank_hist,
    season_cycle = season,
    season_coverage = season_cov,
    meta = list(
      lead_breaks = lead_breaks,
      lead_levels = hindcast_lead_levels(lead_breaks),
      quantile_probs = probs,
      eps_logscore = eps,
      references = unlist(cfg$occ_references),
      rules = rules,
      grains = unlist(cfg$occ_grains)
    )
  )
}

.hindcast_cases_rollup <- function(
  totals,
  monthly,
  cfg,
  taus,
  lead_breaks,
  boot = NULL
) {
  rules <- unlist(cfg$case_rules)
  refs <- score_case_references(totals, unlist(cfg$case_references), taus)
  keys <- c("group", "launch_year", "support", "grain", "scope", "lead_class")
  model <- totals |>
    dplyr::select(
      dplyr::all_of(keys),
      obs,
      med,
      crps_raw = score_crps_raw,
      crps_log1p = score_crps_log1p,
      pinball = score_pinball
    ) |>
    tidyr::pivot_longer(
      c(crps_raw, crps_log1p, pinball),
      names_to = "rule",
      values_to = "score_model"
    ) |>
    dplyr::filter(rule %in% rules)
  scores <- dplyr::inner_join(
    model,
    refs,
    by = c(keys, "rule"),
    relationship = "one-to-many"
  )
  # region grain -> region + africa pools; country grain -> per country
  scoped <- function(df) {
    dplyr::bind_rows(
      df |>
        dplyr::filter(grain == "region") |>
        dplyr::mutate(scope_level = "region", scope_key = group),
      df |>
        dplyr::filter(grain == "region") |>
        dplyr::mutate(scope_level = "africa", scope_key = "africa"),
      df |>
        dplyr::filter(grain == "country") |>
        dplyr::mutate(scope_level = "country", scope_key = scope)
    )
  }
  skill <- .skill_rollup(
    dplyr::mutate(scoped(scores), area = scope, scope = scope_key),
    c("scope_level", "scope", "support", "lead_class", "reference", "rule"),
    .skill_stats_runs(),
    boot,
    if (isTRUE(cfg[["skill_area_summaries"]])) "area"
  )
  calibration <- scoped(totals) |>
    dplyr::group_by(scope_level, scope = scope_key, support, lead_class) |>
    dplyr::summarise(
      n_runs = dplyr::n(),
      obs_total = sum(obs),
      med_total = sum(med),
      sim_obs_ratio = sum(med) / sum(obs),
      cov50 = mean(in50),
      cov95 = mean(in95),
      pct_over = mean(pos == "over"),
      pct_under = mean(pos == "under"),
      median_crps = stats::median(crps),
      mean_cover = mean(cover_popwt),
      .groups = "drop"
    )
  list(
    run_totals = totals,
    scores = scores,
    skill = skill,
    calibration = calibration,
    monthly = monthly,
    meta = list(
      lead_breaks = lead_breaks,
      lead_levels = hindcast_lead_levels(lead_breaks),
      supports = unlist(cfg$case_supports),
      references = unlist(cfg$case_references),
      rules = rules,
      taus = taus,
      estimator = "fair CRPS for crps_raw / crps_log1p (model and reference)"
    )
  )
}

# ===========================================================================
# 7. Report tables (phase-11 numbers and the notebook's summary section)
# ===========================================================================

#' hindcast_report_tables
#'
#' Kable-ready summary of the three bundles under the headline variants and
#' their legacy siblings. Pure function of the bundles.
hindcast_report_tables <- function(
  runs,
  occ,
  cases,
  headline_occ_reference = "loo_gid",
  headline_case_support = "atrisk_gids",
  headline_case_skill_reference = "clim_pooled"
) {
  out <- .hindcast_report_legacy(
    runs,
    occ,
    cases,
    headline_occ_reference,
    headline_case_support
  )
  if (is.null(cases$month_skill)) {
    return(out)
  }
  sm <- cases$skill_meta
  if (
    !headline_case_skill_reference %in% sm$references ||
      !"country" %in% sm$month_grains ||
      !headline_case_support %in% sm$supports
  ) {
    stop(
      "[hindcast-skill] the case-skill headline (support '",
      headline_case_support,
      "', reference '",
      headline_case_skill_reference,
      "', grain 'country') is not among ",
      "the derived case_skill_supports / case_skill_references / case_month_grains.",
      call. = FALSE
    )
  }
  # monthly / season case skill (section 5b): headline support, country units
  # pooled per region and Africa, all leads, against the headline reference
  pick <- function(tbl, rule, nm) {
    d <- tbl |>
      dplyr::filter(
        scope_level %in% c("region", "africa"),
        support == headline_case_support,
        grain == "country",
        lead_class == "all",
        reference == headline_case_skill_reference,
        rule == !!rule
      )
    cols <- intersect(c("crpss", "crpss_lo", "crpss_hi"), names(d))
    dplyr::rename_with(
      dplyr::select(d, scope, dplyr::all_of(cols)),
      ~ sub("^crpss", nm, .x),
      dplyr::all_of(cols)
    )
  }
  for (p in list(
    list(cases$month_skill, "crps_raw", "case_crpss_month"),
    list(cases$month_skill, "crps_log1p", "case_crpss_month_log"),
    list(cases$season_skill, "crps_raw", "case_crpss_season"),
    list(cases$season_skill, "crps_log1p", "case_crpss_season_log")
  )) {
    out$skill_tbl <- dplyr::left_join(
      out$skill_tbl,
      pick(p[[1]], p[[2]], p[[3]]),
      by = "scope"
    )
  }
  by_lead <- function(tbl) {
    tbl |>
      dplyr::filter(
        scope_level %in% c("region", "africa"),
        support == headline_case_support,
        grain == "country"
      ) |>
      dplyr::select(dplyr::any_of(c(
        "scope",
        "lead_class",
        "reference",
        "rule",
        "crpss",
        "crpss_lo",
        "crpss_hi",
        "n_units",
        "n_scored"
      )))
  }
  out$month_by_lead <- by_lead(cases$month_skill)
  out$season_by_lead <- by_lead(cases$season_skill)
  out$headline$case_skill_reference <- headline_case_skill_reference
  out$caveats <- c(
    out$caveats,
    paste0(
      "Monthly and season case skill: references on each run's fixed support (its gids ",
      "and population shares), scaled to the scored year's population",
      if (identical(cases$skill_meta$pop_scale, "wpp")) {
        " (WPP2024 national growth outside the gid table's years)"
      },
      "; ",
      if (identical(cases$skill_meta$mask, "coverage")) {
        "only surveilled ADM1-months are scored, for the model and every reference"
      } else {
        "months without observed cases count as zero"
      },
      "."
    ),
    paste0(
      "Region and Africa monthly / season skill pool the country units; the countries of a ",
      "region come from separate runs under season windows, so cross-country correlation ",
      "of the region-summed ensemble is lost (secondary region_sum table)."
    )
  )
  out
}

.hindcast_report_legacy <- function(
  runs,
  occ,
  cases,
  headline_occ_reference,
  headline_case_support
) {
  occ_pick <- function(grain, reference, rule, lead) {
    occ$skill |>
      dplyr::filter(
        scope_level %in% c("region", "africa"),
        grain == !!grain,
        reference == !!reference,
        rule == !!rule,
        lead_class == !!lead
      ) |>
      dplyr::select(scope, crpss)
  }
  short_lead <- occ$meta$lead_levels[2]
  cell <- occ$cell_scores |>
    dplyr::filter(
      scope_level %in% c("region", "africa"),
      lead_class == "all",
      reference == headline_occ_reference
    ) |>
    dplyr::select(scope, bss, bss_pooled, auc, base_rate)
  case_pick <- function(support, reference, rule) {
    cases$skill |>
      dplyr::filter(
        scope_level %in% c("region", "africa"),
        support == !!support,
        reference == !!reference,
        rule == !!rule,
        lead_class == "all"
      ) |>
      dplyr::select(scope, crpss)
  }
  cal <- cases$calibration |>
    dplyr::filter(
      scope_level %in% c("region", "africa"),
      lead_class == "all",
      support == headline_case_support
    ) |>
    dplyr::select(scope, sim_obs_ratio, cov50, cov95)
  ks <- dplyr::bind_rows(
    runs |>
      dplyr::group_by(scope = group) |>
      dplyr::summarise(
        pit_ks = mean(pit_ks),
        n_runs = dplyr::n(),
        n_cells = sum(n_cells),
        .groups = "drop"
      ),
    runs |>
      dplyr::summarise(
        scope = "africa",
        pit_ks = mean(pit_ks),
        n_runs = dplyr::n(),
        n_cells = sum(n_cells)
      )
  )
  skill_tbl <- ks |>
    dplyr::left_join(
      dplyr::rename(
        occ_pick("window", headline_occ_reference, "crps_fair", "all"),
        occ_crpss = crpss
      ),
      by = "scope"
    ) |>
    dplyr::left_join(
      dplyr::rename(
        occ_pick("window", headline_occ_reference, "crps_fair", short_lead),
        occ_crpss_short = crpss
      ),
      by = "scope"
    ) |>
    dplyr::left_join(
      dplyr::rename(
        occ_pick("calendar_year", "insample_gid", "crps_emp", "all"),
        occ_crpss_legacy = crpss
      ),
      by = "scope"
    ) |>
    dplyr::left_join(cell, by = "scope") |>
    dplyr::left_join(
      dplyr::rename(
        case_pick(headline_case_support, "loo_climatology", "crps_raw"),
        case_crpss = crpss
      ),
      by = "scope"
    ) |>
    dplyr::left_join(
      dplyr::rename(
        case_pick(headline_case_support, "loo_climatology", "crps_log1p"),
        case_crpss_log = crpss
      ),
      by = "scope"
    ) |>
    dplyr::left_join(cal, by = "scope")
  list(
    skill_tbl = skill_tbl,
    occ_by_lead = occ$skill |>
      dplyr::filter(
        scope_level %in% c("region", "africa"),
        grain == "window",
        rule == "crps_fair"
      ) |>
      dplyr::select(scope, lead_class, reference, crpss, n_units),
    cases_by_support = cases$calibration |>
      dplyr::filter(
        scope_level %in% c("region", "africa"),
        lead_class == "all"
      ) |>
      dplyr::left_join(
        cases$skill |>
          dplyr::filter(
            reference == "loo_climatology",
            rule == "crps_raw",
            lead_class == "all"
          ) |>
          dplyr::select(scope_level, scope, support, crpss),
        by = c("scope_level", "scope", "support")
      ),
    counts = tibble::tibble(
      n_runs = nrow(runs),
      n_groups = dplyr::n_distinct(runs$group),
      n_cells = sum(runs$n_cells),
      n_gid_windows = sum(
        occ$gid_window$grain == "window" &
          occ$gid_window$lead_class == "all"
      ),
      launch_years = paste(range(runs$launch_year), collapse = "-")
    ),
    headline = list(
      occ_reference = headline_occ_reference,
      case_support = headline_case_support,
      short_lead = short_lead
    ),
    caveats = c(
      paste0(
        "The occurrence, duration and monthly fits were trained on 2010-2023, ",
        "which includes every hindcast year: the hindcasts are not out-of-sample ",
        "for the fixed or random effects."
      ),
      paste0(
        "Occurrence support: at-risk (gid, month) rows of the survival panel inside ",
        "each run's window, for both the model and every reference."
      ),
      paste0(
        "Occurrence references leave out every month of the scored window ",
        "(overlapping launches share calendar months)."
      ),
      paste0(
        "Case support '",
        headline_case_support,
        "' is defined by surveillance ",
        "(sim gids at risk in the window), not by the outcome; the old basis B ",
        "(obs_positive_keys) drops every false alarm and is kept as a sibling."
      )
    )
  )
}
