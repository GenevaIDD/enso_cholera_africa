# sim_io.R — simulator I/O: replicate files, the compact monthly summary,
#   task plans, partition merging and the forward-run manifest readers
#
# Provides:
#   sim_io_save(), sim_io_read()     replicate artefacts (qs / rds)
#   sim_monthly_compact(), sim_monthly_expand(), sim_monthly_save()
#                                    compact monthly summary (format 2)
#   build_task_plan(), write_task_yamls()   TU x replicate-chunk task plan
#   build_obs_panel()                the filter's observation panel
#   aggregate_partitions(), combine_group_runs()   merge simulator outputs
#   manifest_units(), manifest_runs(), manifest_country_horizons(),
#   manifest_arm_launches()          forecast / hindcast / event manifests
#   censor_sim_by_country()          per-country firm-horizon censoring
#   assert_arm_launch(), assert_panel_covers_horizon()   run guards
#
# Details:
#   Shared by 07 (tasks, kernels, postprocess, figures), 08, the config
#   generators and the results phase.

#' Build per-replicate warmup panels via importance-weighted trajectory selection.
#'
#' Scores each of the 100 elag10 phire simulation replicates against observed
#' outbreak dynamics (Bernoulli agreement) in the scoring window, samples
#' n_rep trajectories proportional to softmax(score / temperature), then
#' assembles warmup panels for the 12-month window before the horizon start.
#'
#' @param sim_monthly_path Path to `monthly_summary.qs` from the elag10 phire
#'   simulation (100 replicates, columns: replicate_id, gid, month, Y).
#' @param surv_panel       Survival panel from `make_surv_dat_file()` — used to
#'   extract observed outbreak indicators and lag columns for the warmup window.
#' @param warmup_start     Start of the warmup window (Date): horizon_start - warmup_months.
#' @param warmup_end       End of the warmup window (Date): horizon_start (exclusive).
#' @param n_rep            Number of forecast replicates to produce.
#' @param scoring_start    Start of the scoring window (Date); default 21 months before warmup_end.
#' @param scoring_end      End of the scoring window (Date); default warmup_end - 1 month.
#' @param score_temperature Softmax temperature τ (default 1.0): lower → closer to argmax.
#' @return List of length n_rep.  Each element is a mini-panel tibble covering
#'   [warmup_start, warmup_end) with columns:
#'   gid, t_start, ongoing_outbreak_lag1, ongoing_outbreak_country_lag1,
#'   ongoing_outbreak_tu_lag1, time_to_last_outbreak, time_to_last_outbreak_cat.
build_hybrid_warmup_panels <- function(
  sim_monthly_path,
  surv_panel,
  warmup_start,
  warmup_end,
  n_rep,
  scoring_start = NULL,
  scoring_end = NULL,
  score_temperature = 1.0,
  tu_map = NULL
) {
  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("build_hybrid_warmup_panels: qs2 package required.")
  }

  warmup_start <- as.Date(warmup_start)
  warmup_end <- as.Date(warmup_end)

  if (is.null(scoring_start)) {
    scoring_start <- warmup_end %m-% months(21L)
  }
  if (is.null(scoring_end)) {
    scoring_end <- warmup_end %m-% months(1L)
  }
  scoring_start <- as.Date(scoring_start)
  scoring_end <- as.Date(scoring_end)

  # 1. Load simulation monthly summary.
  message(
    "[build_hybrid_warmup_panels] loading sim monthly summary: ",
    sim_monthly_path
  )
  sim_ms <- sim_io_read(sim_monthly_path, "qs") # a compact (format-2) file comes back long

  # Validate required columns.
  required_sim <- c("replicate_id", "gid", "month", "Y")
  missing_sim <- setdiff(required_sim, names(sim_ms))
  if (length(missing_sim) > 0) {
    warning(
      "[build_hybrid_warmup_panels] sim monthly_summary missing columns: ",
      paste(missing_sim, collapse = ", "),
      " — falling back to uniform warmup weighting."
    )
    sim_ms_ok <- FALSE
  } else {
    sim_ms_ok <- TRUE
    sim_ms <- sim_ms %>%
      dplyr::mutate(month = as.Date(month), gid = as.character(gid))
  }

  # 2. Extract observed outbreak indicator from surv_panel for the scoring window.
  obs_score <- if (sim_ms_ok) {
    surv_panel %>%
      dplyr::select(gid, t_start, outbreak) %>%
      dplyr::mutate(gid = as.character(gid), t_start = as.Date(t_start)) %>%
      dplyr::filter(
        t_start >= scoring_start,
        t_start <= scoring_end,
        !is.na(outbreak)
      ) %>%
      dplyr::group_by(gid, t_start) %>%
      dplyr::summarise(
        obs_Y = as.integer(max(outbreak, na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    tibble::tibble(
      gid = character(0),
      t_start = as.Date(character(0)),
      obs_Y = integer(0)
    )
  }

  n_obs_pairs <- nrow(obs_score)
  message(sprintf(
    "[build_hybrid_warmup_panels] scoring window %s – %s: %d observed (gid, month) pairs",
    scoring_start,
    scoring_end,
    n_obs_pairs
  ))

  # 3. Score replicates against observations.
  sim_reps <- sort(unique(sim_ms$replicate_id))
  n_sim_reps <- length(sim_reps)

  weights <- rep(1.0, n_sim_reps) # default: uniform
  if (sim_ms_ok && n_obs_pairs > 0) {
    sim_score_dat <- sim_ms %>%
      dplyr::filter(
        month >= scoring_start,
        month <= scoring_end,
        gid %in% obs_score$gid
      ) %>%
      dplyr::rename(t_start = month) %>%
      dplyr::inner_join(obs_score, by = c("gid", "t_start"))

    if (nrow(sim_score_dat) > 0) {
      scores <- sim_score_dat %>%
        dplyr::group_by(replicate_id) %>%
        dplyr::summarise(
          score = mean(Y == obs_Y, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(match(replicate_id, sim_reps))

      # Softmax over agreement scores.
      raw_scores <- scores$score
      log_w <- raw_scores / score_temperature
      log_w <- log_w - max(log_w, na.rm = TRUE) # numeric stability
      weights <- exp(log_w)
      weights <- weights / sum(weights)
      message(sprintf(
        "[build_hybrid_warmup_panels] replicate scores: min=%.3f median=%.3f max=%.3f",
        min(raw_scores),
        stats::median(raw_scores),
        max(raw_scores)
      ))
    } else {
      message(
        "[build_hybrid_warmup_panels] no sim–obs pairs found in scoring window; using uniform weights."
      )
    }
  }

  # 4. Sample n_rep trajectory indices with replacement.
  warmup_rep_idx <- sample(
    seq_len(n_sim_reps),
    size = n_rep,
    replace = TRUE,
    prob = weights
  )
  selected_sim_reps <- sim_reps[warmup_rep_idx]
  message(sprintf(
    "[build_hybrid_warmup_panels] sampled %d trajectory indices (from %d sim reps).",
    n_rep,
    n_sim_reps
  ))

  # 5. Extract observed lag columns for the warmup window.
  lag_obs_cols <- c(
    "ongoing_outbreak_lag1",
    "ongoing_outbreak_country_lag1",
    "ongoing_outbreak_tu_lag1",
    "time_to_last_outbreak",
    "time_to_last_outbreak_cat"
  )
  obs_warmup <- surv_panel %>%
    dplyr::select(gid, t_start, dplyr::any_of(lag_obs_cols)) %>%
    dplyr::mutate(gid = as.character(gid), t_start = as.Date(t_start)) %>%
    dplyr::filter(t_start >= warmup_start, t_start < warmup_end) %>%
    dplyr::group_by(gid, t_start) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::any_of(c(
          "ongoing_outbreak_lag1",
          "ongoing_outbreak_country_lag1"
        )),
        ~ as.integer(max(dplyr::coalesce(., 0L))),
        .names = "{.col}"
      ),
      dplyr::across(
        dplyr::any_of("ongoing_outbreak_tu_lag1"),
        ~ {
          v <- .[!is.na(.)]
          if (length(v) == 0L) NA_integer_ else as.integer(max(v))
        },
        .names = "{.col}"
      ),
      dplyr::across(
        dplyr::any_of("time_to_last_outbreak"),
        ~ suppressWarnings(min(., na.rm = TRUE)),
        .names = "{.col}"
      ),
      dplyr::across(
        dplyr::any_of("time_to_last_outbreak_cat"),
        ~ if (length(.) == 0L) NA_character_ else .[1],
        .names = "{.col}"
      ),
      .groups = "drop"
    )
  obs_gids_warmup <- unique(obs_warmup$gid)

  # Simulated warmup slice (all unique sim reps that were selected, for efficiency).
  unique_selected <- unique(selected_sim_reps)
  sim_warmup_all <- if (sim_ms_ok) {
    sim_ms %>%
      dplyr::filter(
        replicate_id %in% unique_selected,
        month >= warmup_start,
        month < warmup_end
      ) %>%
      dplyr::rename(t_start = month, sim_Y = Y)
  } else {
    tibble::tibble(
      replicate_id = integer(0),
      gid = character(0),
      t_start = as.Date(character(0)),
      sim_Y = integer(0)
    )
  }

  # Pre-warmup TTLO seed: for each (replicate, gid), the last simulated
  # outbreak month *before* warmup_start.  Used to seed the warmup-window TTLO
  # scan so gids whose most recent outbreak predates the warmup window are
  # assigned the correct (lower) TTLO bin rather than defaulting to [6,Inf].
  ttlo_seed_all <- if (sim_ms_ok) {
    sim_ms %>%
      dplyr::filter(
        replicate_id %in% unique_selected,
        month < warmup_start,
        Y == 1L
      ) %>%
      dplyr::group_by(replicate_id, gid) %>%
      dplyr::summarise(last_ob_pre_warmup = max(month), .groups = "drop")
  } else {
    tibble::tibble(
      replicate_id = integer(0),
      gid = character(0),
      last_ob_pre_warmup = as.Date(character(0))
    )
  }

  all_gids <- sort(unique(as.character(surv_panel$gid)))
  warmup_months_seq <- seq.Date(
    warmup_start,
    warmup_end %m-% months(1L),
    by = "1 month"
  )

  # Precompute TU-level activity per (replicate_id, t_start, country) from
  # sim_warmup_all so ongoing_outbreak_tu_lag1 can be filled from simulated
  # trajectories.  Mirrors the fit-time logic in 03_prepare_occ.R:
  #   country_active = any gid in country had sim_Y == 1
  #   o_tu = any OTHER country in the same TU had country_active == 1
  # Countries absent from tu_map$tu (NA TU) are left as NA — matching fit-time
  # behaviour where ongoing_outbreak_tu_lag1 is NA for non-TU locations.
  # Falls back to NULL (→ NA_integer_ for all rows) when tu_map is not supplied.
  sim_tu_activity <- if (!is.null(tu_map) && sim_ms_ok) {
    gid_country_lu <- tu_map %>%
      dplyr::select(gid, country) %>%
      dplyr::distinct() %>%
      dplyr::mutate(gid = as.character(gid))
    country_to_tu_lu <- tu_map %>%
      dplyr::select(country, tu) %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(tu))
    sta <- sim_warmup_all %>%
      dplyr::left_join(gid_country_lu, by = "gid") %>%
      dplyr::filter(!is.na(country)) %>%
      dplyr::group_by(replicate_id, t_start, country) %>%
      dplyr::summarise(
        country_active = as.integer(max(dplyr::coalesce(sim_Y, 0L))),
        .groups = "drop"
      ) %>%
      dplyr::left_join(country_to_tu_lu, by = "country") %>%
      dplyr::group_by(replicate_id, t_start, tu) %>%
      dplyr::mutate(tu_sum = sum(country_active, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        o_tu_country = as.integer((tu_sum - country_active) >= 1L)
      ) %>%
      dplyr::select(replicate_id, t_start, country, o_tu_country)
    sta
  } else {
    NULL
  }

  # 6. Build each of the n_rep warmup mini-panels.
  panels_list <- vector("list", n_rep)
  for (r_forecast in seq_len(n_rep)) {
    sim_rep <- selected_sim_reps[r_forecast]

    # Base grid for all gids × warmup months.
    base_grid <- expand.grid(
      gid = all_gids,
      t_start = warmup_months_seq,
      stringsAsFactors = FALSE
    ) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(t_start = as.Date(t_start))

    # Simulated lag columns for this replicate.
    sim_r <- sim_warmup_all %>%
      dplyr::filter(replicate_id == sim_rep) %>%
      dplyr::select(gid, t_start, sim_Y)

    panel_r <- base_grid %>%
      dplyr::left_join(sim_r, by = c("gid", "t_start")) %>%
      dplyr::mutate(
        ongoing_outbreak_lag1 = as.integer(dplyr::coalesce(sim_Y, 0L))
      )

    # Reconstruct country and TU lags via rolling-max over same (country/TU, month).
    # We need a country/TU mapping: derive from surv_panel.
    if (!exists(".warmup_ct_map", inherits = FALSE)) {
      .warmup_ct_map <- surv_panel %>%
        dplyr::select(gid, country, AFRO_region) %>%
        dplyr::distinct() %>%
        dplyr::mutate(gid = as.character(gid))
    }
    panel_r <- panel_r %>%
      dplyr::left_join(.warmup_ct_map, by = "gid") %>%
      dplyr::group_by(country, t_start) %>%
      dplyr::mutate(
        ongoing_outbreak_country_lag1 = as.integer(
          (sum(dplyr::coalesce(ongoing_outbreak_lag1, 0L)) -
            dplyr::coalesce(ongoing_outbreak_lag1, 0L)) >=
            1L
        )
      ) %>%
      dplyr::ungroup()

    # Compute ongoing_outbreak_tu_lag1 from sim trajectories.  Fall back to
    # NA_integer_ (backward-compatible) when no tu_map was supplied.
    panel_r <- if (!is.null(sim_tu_activity)) {
      tu_r <- sim_tu_activity %>%
        dplyr::filter(replicate_id == sim_rep) %>%
        dplyr::select(t_start, country, o_tu_country)
      panel_r %>%
        dplyr::left_join(tu_r, by = c("t_start", "country")) %>%
        dplyr::mutate(ongoing_outbreak_tu_lag1 = o_tu_country) %>%
        dplyr::select(-o_tu_country)
    } else {
      panel_r %>% dplyr::mutate(ongoing_outbreak_tu_lag1 = NA_integer_)
    }

    # Pre-warmup TTLO seed for this replicate (days from last outbreak before
    # warmup_start to the first warmup month, per gid).
    seed_r <- ttlo_seed_all %>%
      dplyr::filter(replicate_id == sim_rep) %>%
      dplyr::select(gid, .last_ob_pre_warmup = last_ob_pre_warmup)

    # Reconstruct time_to_last_outbreak from simulated outbreak sequence.
    # Scan backward for each gid: find most recent month with sim_Y == 1.
    # Seed last_seen from pre-warmup history so gids whose most recent
    # outbreak predates the warmup window receive the correct (lower) TTLO bin
    # instead of defaulting to [6,Inf] (which would inflate hazard at the
    # first forecast step).
    # Unit: days (matching 03_prepare_occ.R and the expected input of
    # make_time_to_last_outbreak_cat — do NOT divide by 365 here).
    panel_r <- panel_r %>%
      dplyr::left_join(seed_r, by = "gid") %>%
      dplyr::group_by(gid) %>%
      dplyr::arrange(t_start, .by_group = TRUE) %>%
      dplyr::mutate(
        .last_ob_month = dplyr::if_else(
          ongoing_outbreak_lag1 == 1L,
          t_start,
          as.Date(NA)
        ),
        time_to_last_outbreak = {
          lobs <- .last_ob_month
          pre_date <- .last_ob_pre_warmup[1L]
          result_ttlo <- numeric(dplyr::n())
          # Seed from pre-warmup history; NA if gid truly never had an outbreak.
          last_seen <- if (!is.na(pre_date)) {
            as.numeric(t_start[1L] - pre_date) # days
          } else {
            NA_real_
          }
          for (ii in seq_len(dplyr::n())) {
            if (!is.na(lobs[ii])) {
              last_seen <- as.numeric(t_start[ii] - lobs[ii])
            }
            result_ttlo[ii] <- if (is.na(last_seen)) NA_real_ else last_seen # days
          }
          result_ttlo
        }
      ) %>%
      dplyr::select(-.last_ob_month, -.last_ob_pre_warmup) %>%
      dplyr::ungroup()

    # Override with observed values where available.
    if (nrow(obs_warmup) > 0) {
      panel_r <- panel_r %>%
        dplyr::rows_update(
          obs_warmup %>%
            dplyr::select(
              gid,
              t_start,
              dplyr::any_of(c(
                "ongoing_outbreak_lag1",
                "ongoing_outbreak_country_lag1",
                "ongoing_outbreak_tu_lag1",
                "time_to_last_outbreak"
              ))
            ),
          by = c("gid", "t_start"),
          unmatched = "ignore"
        )
    }

    # Apply make_time_to_last_outbreak_cat to derive the factor column.
    # make_time_to_last_outbreak_cat overwrites the numeric time_to_last_outbreak
    # column in-place; rename it to time_to_last_outbreak_cat so the kernel's
    # build_warmup_state() can find it under the expected column name.
    panel_r <- panel_r %>%
      make_time_to_last_outbreak_cat() %>%
      dplyr::rename(time_to_last_outbreak_cat = time_to_last_outbreak)

    panels_list[[r_forecast]] <- panel_r %>%
      dplyr::select(
        gid,
        t_start,
        dplyr::any_of(c(
          "ongoing_outbreak_lag1",
          "ongoing_outbreak_country_lag1",
          "ongoing_outbreak_tu_lag1",
          "time_to_last_outbreak",
          "time_to_last_outbreak_cat"
        ))
      )
  }

  panels_list
}

# Simulator IO helpers — dispatch save/read/extension on output format.
sim_io_ext <- function(fmt) {
  switch(fmt %||% "qs", rds = ".rds", parquet = ".parquet", ".qs")
}
sim_io_save <- function(x, path, fmt) {
  fmt <- fmt %||% "qs"
  if (identical(fmt, "rds")) {
    saveRDS(x, path)
  } else if (identical(fmt, "parquet")) {
    if (!requireNamespace("nanoparquet", quietly = TRUE)) {
      stop(
        "Package 'nanoparquet' is required for parquet output. ",
        "Install with: install.packages('nanoparquet')"
      )
    }
    nanoparquet::write_parquet(x, path)
  } else {
    qs2::qs_save(x, path)
  }
}
#' Read a simulator output. A compact monthly summary (sim_monthly_compact())
#' comes back as the long table it encodes (identical() to the one written),
#' unless `expand = FALSE`.
sim_io_read <- function(path, fmt, expand = TRUE) {
  fmt <- fmt %||% "qs"
  x <- if (identical(fmt, "rds")) {
    readRDS(path)
  } else if (identical(fmt, "parquet")) {
    if (!requireNamespace("nanoparquet", quietly = TRUE)) {
      stop(
        "Package 'nanoparquet' is required for parquet input. ",
        "Install with: install.packages('nanoparquet')"
      )
    }
    nanoparquet::read_parquet(path)
  } else {
    qs2::qs_read(path)
  }
  if (isTRUE(expand) && is_sim_monthly_compact(x)) sim_monthly_expand(x) else x
}

# ---------------------------------------------------------------------------
# Compact monthly summary (format 2)
# ---------------------------------------------------------------------------
# simulate_one_tu()'s monthly_summary is a dense replicate x gid x month grid
# of mostly inactive cells. Format 2 stores the three axes once and each value
# column sparsely: the cells that differ from the default (0 for integers, NA
# for doubles), or values on the cells of an earlier column ("ref"). Same file
# name; sim_io_read() expands it to the identical long table. Few-cell readers
# (results_core.R::load_sim_run_monthly()) use sim_monthly_cells(). A table
# that is not a dense grid, or fails the round trip, is written long. Merged
# files hold one block per partition. Gate: tests/test_sim_monthly_compact.R.

SIM_MONTHLY_FORMAT <- "sim_monthly_compact"
SIM_MONTHLY_VERSION <- 2L

is_sim_monthly_compact <- function(x) {
  is.list(x) &&
    !is.data.frame(x) &&
    identical(x[["format"]], SIM_MONTHLY_FORMAT)
}

# Sorted cell indices as their first value then the gaps (small integers
# compress far better); .undiff_idx() inverts it.
.diff_idx <- function(i) if (length(i)) c(i[1L], diff(i)) else i
.undiff_idx <- function(d) cumsum(d)

# A value column's default: NA for a mostly-NA double column, else 0; and the
# cells that differ from it (NA counts as a value against a 0 default).
.sim_col_default_value <- function(x) {
  if (is.double(x) && sum(is.na(x)) * 2 > length(x)) NA else 0
}
.sim_col_nondefault <- function(x, d = .sim_col_default_value(x)) {
  if (is.na(d)) !is.na(x) else is.na(x) | x != 0
}

#' One column of a dense grid, encoded against its default; `refs` holds the
#' non-default cells of the columns encoded before it.
.encode_sim_col <- function(x, refs) {
  t <- typeof(x)
  if (!t %in% c("integer", "double", "logical") || !is.null(attributes(x))) {
    return(list(enc = "dense", x = x))
  }
  d <- .sim_col_default_value(x)
  nz <- which(.sim_col_nondefault(x, d))
  # values on an earlier column's non-default cells, when those cover this one's
  for (r in names(refs)) {
    if (length(refs[[r]]) <= 2L * length(nz) && all(nz %in% refs[[r]])) {
      return(list(
        enc = "ref",
        type = t,
        default = d,
        ref = r,
        val = x[refs[[r]]],
        nz = nz
      ))
    }
  }
  v <- x[nz]
  ones <- t != "double" && length(v) && !anyNA(v) && all(v == 1)
  list(
    enc = "idx",
    type = t,
    default = d,
    idx = .diff_idx(nz),
    val = if (ones) NULL else v,
    nz = nz
  )
}

# A column of n default values; the value 1 of a column's type.
.sim_col_default <- function(col, n) {
  if (is.na(col$default)) {
    rep(switch(col$type, integer = NA_integer_, double = NA_real_, NA), n)
  } else {
    vector(col$type, n)
  }
}
.sim_col_one <- function(col) switch(col$type, integer = 1L, double = 1, TRUE)

.decode_sim_col <- function(col, n, refs) {
  if (identical(col$enc, "dense")) {
    return(col$x)
  }
  x <- .sim_col_default(col, n)
  if (identical(col$enc, "ref")) {
    x[refs[[col$ref]]] <- col$val
  } else {
    x[.undiff_idx(col$idx)] <- if (is.null(col$val)) {
      .sim_col_one(col)
    } else {
      col$val
    }
  }
  x
}

#' The non-default cells of each encoded column of a block (decoded indices).
.sim_block_nz <- function(b) {
  refs <- list()
  for (nm in names(b$cols)) {
    col <- b$cols[[nm]]
    refs[[nm]] <- switch(
      col$enc,
      dense = NULL,
      idx = .undiff_idx(col$idx),
      ref = {
        r <- refs[[col$ref]]
        v <- col$val
        r[if (is.na(col$default)) !is.na(v) else (is.na(v) | v != 0)]
      }
    )
  }
  refs
}

#' sim_monthly_compact -- a dense monthly_summary as a format-2 block, or NULL
#' when it is not one (any row order or column set the axes cannot rebuild).
sim_monthly_compact <- function(ms) {
  if (
    !tibble::is_tibble(ms) ||
      !nrow(ms) ||
      length(ms) < 4L ||
      !identical(names(ms)[1:3], c("replicate_id", "gid", "month"))
  ) {
    return(NULL)
  }
  rids <- unique(ms$replicate_id)
  gids <- unique(ms$gid)
  months <- unique(ms$month)
  R <- length(rids)
  N <- length(gids)
  M <- length(months)
  if (
    as.double(R) * N * M != nrow(ms) ||
      !identical(ms$replicate_id, rep(rids, each = N * M)) ||
      !identical(ms$gid, rep(rep(gids, each = M), times = R)) ||
      !identical(ms$month, rep(rep(months, times = N), times = R))
  ) {
    return(NULL)
  }
  cols <- list()
  refs <- list()
  for (nm in names(ms)[-(1:3)]) {
    e <- .encode_sim_col(ms[[nm]], refs)
    if (!is.null(e$nz)) {
      refs[[nm]] <- e$nz
    }
    e$nz <- NULL
    cols[[nm]] <- e
  }
  list(
    format = SIM_MONTHLY_FORMAT,
    version = SIM_MONTHLY_VERSION,
    bind = FALSE,
    blocks = list(list(
      replicate_ids = rids,
      gids = gids,
      months = months,
      cols = cols
    ))
  )
}

.sim_block_n <- function(b) {
  length(b$replicate_ids) * length(b$gids) * length(b$months)
}

.sim_block_expand <- function(b, cols = NULL) {
  R <- length(b$replicate_ids)
  N <- length(b$gids)
  M <- length(b$months)
  n <- R * N * M
  out <- tibble::tibble(
    replicate_id = rep(b$replicate_ids, each = N * M),
    gid = rep(rep(b$gids, each = M), times = R),
    month = rep(rep(b$months, times = N), times = R)
  )
  refs <- .sim_block_nz(b)
  for (nm in names(b$cols)) {
    if (is.null(cols) || nm %in% cols) {
      out[[nm]] <- .decode_sim_col(b$cols[[nm]], n, refs)
    }
  }
  out
}

#' sim_monthly_expand -- the long table a compact object encodes: one block as
#' written, several (a merged file) as dplyr::bind_rows() of the partitions.
#' `cols` (NULL = all): the value columns to build (the long table with the
#' others dropped, in its column order), so a reader that drops columns never
#' builds them.
sim_monthly_expand <- function(x, cols = NULL) {
  if (!isTRUE(x$bind) && length(x$blocks) == 1L) {
    return(.sim_block_expand(x$blocks[[1L]], cols))
  }
  dplyr::bind_rows(lapply(x$blocks, .sim_block_expand, cols = cols))
}

#' sim_monthly_bind -- the partitions' monthly summaries as one object for the
#' merged file: the blocks of compact parts concatenated (expanding to
#' bind_rows() of the parts); any part in the long format -> the long bind.
sim_monthly_bind <- function(parts) {
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) && all(vapply(parts, is_sim_monthly_compact, logical(1)))) {
    return(list(
      format = SIM_MONTHLY_FORMAT,
      version = SIM_MONTHLY_VERSION,
      bind = TRUE,
      blocks = unlist(lapply(parts, `[[`, "blocks"), recursive = FALSE)
    ))
  }
  dplyr::bind_rows(lapply(parts, function(p) {
    if (is_sim_monthly_compact(p)) sim_monthly_expand(p) else p
  }))
}

#' sim_monthly_save -- write a monthly summary: format 2 for qs / rds when the
#' table is a dense grid that round-trips identically (else, and for parquet,
#' the long table). A compact object is written as is (qs / rds) or expanded
#' (parquet). Long format forced with options(cholera.sim_monthly_long = TRUE).
sim_monthly_save <- function(ms, path, fmt) {
  fmt <- fmt %||% "qs"
  if (is_sim_monthly_compact(ms)) {
    if (identical(fmt, "parquet")) {
      ms <- sim_monthly_expand(ms)
    }
    return(sim_io_save(ms, path, fmt))
  }
  if (
    !identical(fmt, "parquet") && !isTRUE(getOption("cholera.sim_monthly_long"))
  ) {
    cx <- sim_monthly_compact(ms)
    if (!is.null(cx) && identical(sim_monthly_expand(cx), ms)) {
      return(sim_io_save(cx, path, fmt))
    }
  }
  sim_io_save(ms, path, fmt)
}

# Conversion sidecar. utils/compact_sim_outputs.R rewrites an
# existing long-format file in format 2 only after it reads back identical,
# keeping its mtime; the file still gets smaller, and the 09 stamps record each
# input's size as well as its mtime (results_core.R::.results_stamp_stale_reason()),
# so every artefact that listed it would read as stale. The tool writes
# `<file>.format2` (original size and mtime, new size); a stamp whose recorded
# size and mtime are the original ones is still fresh for that file as long as
# the file is exactly what the tool left (same mtime, the recorded new size).

sim_monthly_sidecar_path <- function(path) paste0(path, ".format2")

#' Record a conversion of `path` (now on disk) from `orig_size` bytes / `orig_mtime`.
sim_monthly_write_sidecar <- function(path, orig_size, orig_mtime) {
  write.dcf(
    data.frame(
      format = SIM_MONTHLY_FORMAT,
      version = SIM_MONTHLY_VERSION,
      original_size = format(orig_size, scientific = FALSE),
      original_mtime = format(
        as.numeric(orig_mtime),
        digits = 17,
        scientific = FALSE
      ),
      new_size = format(file.size(path), scientific = FALSE),
      converted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    ),
    sim_monthly_sidecar_path(path)
  )
}

#' Whether each file in `path` is the in-place conversion of the file a stamp
#' recorded (`rec_mtime`, `rec_size`): its sidecar names that original size and
#' mtime, the file still has that mtime and the size the conversion left.
#' FALSE without a readable sidecar.
sim_monthly_converted_from <- function(path, rec_mtime, rec_size) {
  now_mtime <- as.numeric(file.mtime(path))
  now_size <- file.size(path)
  vapply(
    seq_along(path),
    function(i) {
      sc <- sim_monthly_sidecar_path(path[i])
      if (!file.exists(sc)) {
        return(FALSE)
      }
      d <- tryCatch(read.dcf(sc)[1L, ], error = function(e) NULL)
      if (
        is.null(d) ||
          !all(c("original_size", "original_mtime", "new_size") %in% names(d))
      ) {
        return(FALSE)
      }
      isTRUE(as.numeric(d[["original_size"]]) == rec_size[i]) &&
        isTRUE(as.numeric(d[["new_size"]]) == now_size[i]) &&
        isTRUE(
          abs(as.numeric(d[["original_mtime"]]) - as.numeric(rec_mtime[i])) <
            1e-3
        ) &&
        isTRUE(abs(now_mtime[i] - as.numeric(rec_mtime[i])) < 1e-3)
    },
    logical(1)
  )
}

#' sim_monthly_cells
#'
#' The rows of a monthly summary (compact or long) a reader of sums needs: every
#' cell where one of `value_cols` differs from its default (Y != 0, cases_month
#' not NA, ...) and every cell of the gids in `keep_gids`, in the long table's
#' row order -- a row subset of the long table (identical columns), so sums,
#' maxima and the key combinations of any grouping that contains one of
#' `keep_gids` per group are those of the full table. `keep_gids` is a function
#' of the gid axis (the full gid vector, first-seen order) returning the gids to
#' keep. A long table is filtered the same way.
#' @return list(rows, gids, replicate_ids, months): the subset and the full
#'   table's distinct gids / replicate ids / months in first-seen order.
sim_monthly_cells <- function(
  x,
  value_cols,
  keep_gids = function(g) character(0)
) {
  if (!is_sim_monthly_compact(x)) {
    g <- unique(x$gid)
    keep <- x$gid %in% keep_gids(g)
    for (v in intersect(value_cols, names(x))) {
      keep <- keep | .sim_col_nondefault(x[[v]])
    }
    return(list(
      rows = x[keep, , drop = FALSE],
      gids = g,
      replicate_ids = unique(x$replicate_id),
      months = unique(x$month)
    ))
  }
  gids <- unique(unlist(lapply(x$blocks, `[[`, "gids"), use.names = FALSE))
  kg <- keep_gids(gids)
  one <- function(b) {
    R <- length(b$replicate_ids)
    N <- length(b$gids)
    M <- length(b$months)
    n <- R * N * M
    refs <- .sim_block_nz(b)
    cell <- sort(unique(c(
      unlist(refs[intersect(value_cols, names(refs))], use.names = FALSE),
      if (any(b$gids %in% kg)) {
        gi <- which(b$gids %in% kg)
        off <- (seq_len(R) - 1) * N * M
        if (n <= .Machine$integer.max) {
          off <- as.integer(off)
        }
        as.vector(outer(seq_len(M), outer((gi - 1L) * M, off, `+`), `+`))
      }
    )))
    if (n <= .Machine$integer.max) {
      cell <- as.integer(cell)
    } # integer %/% is several times faster
    r0 <- (cell - 1L) %/% (N * M)
    g0 <- ((cell - 1L) %/% M) %% N
    m0 <- (cell - 1L) %% M
    out <- tibble::tibble(
      replicate_id = b$replicate_ids[r0 + 1L],
      gid = b$gids[g0 + 1L],
      month = b$months[m0 + 1L]
    )
    for (nm in names(b$cols)) {
      col <- b$cols[[nm]]
      out[[nm]] <- if (identical(col$enc, "dense")) {
        col$x[cell]
      } else {
        v <- .sim_col_default(col, length(cell))
        nz <- refs[[nm]]
        vals <- if (identical(col$enc, "ref")) {
          col$val[match(nz, refs[[col$ref]])]
        } else if (is.null(col$val)) {
          rep(.sim_col_one(col), length(nz))
        } else {
          col$val
        }
        hit <- match(nz, cell)
        v[hit[!is.na(hit)]] <- vals[!is.na(hit)]
        v
      }
    }
    out
  }
  rows <- if (length(x$blocks) == 1L && !isTRUE(x$bind)) {
    one(x$blocks[[1L]])
  } else {
    dplyr::bind_rows(lapply(x$blocks, one))
  }
  list(
    rows = rows,
    gids = gids,
    replicate_ids = unique(unlist(
      lapply(x$blocks, `[[`, "replicate_ids"),
      use.names = FALSE
    )),
    months = unique(unlist(lapply(x$blocks, `[[`, "months"), use.names = FALSE))
  )
}
# Serialisation format for NON-tabular simulator outputs (the validation
# `metrics` list, hindcast skill tables): parquet cannot hold a nested list, so
# rds stays rds and everything else (qs, parquet, NULL) is written as qs.
# Single rule shared by 07c's run_validation()/run_hindcast_skill() and 07d's
# metrics reader.
sim_metrics_fmt <- function(out_fmt) {
  if (identical(out_fmt %||% "qs", "rds")) "rds" else "qs"
}

# ---------------------------------------------------------------------------
# Simulator task plan — TU x replicate-chunk assignments
# ---------------------------------------------------------------------------
# Canonical, single-source-of-truth planner for the v1 simulator (07). Shared
# by both entry points: the local furrr full-run path (07b_simulator_v1.R) and
# the SLURM task-materialisation path (07a_simulator_tasks.R). One definition
# keeps the seed striding (no colliding (seed + rep_id) ranges across tasks)
# and the per-track posterior indices identical on both paths. (The spatPomp engine
# (08) has its own separate fan-out /
# seed-striding scheme and does not use this planner.)

#' Build the (TU, replicate-chunk) plan from the resolved simulator config.
#' The TU list is taken straight from `scope.transmission_units`; "all" expands
#' to the real TUs derived from the fitted/location `frame` (or the three known
#' TUs hmmW / hmmCE / hmmS when no frame is supplied).
#'
#' @param opt   resolved options list (uses opt$simulator)
#' @param frame optional location frame with a `tu` column; when supplied,
#'   `transmission_units == "all"` expands to its distinct TUs (fit-free —
#'   callers may pass a frame built without loading fitted models).
#' @return list of per-task rows (task_id, tu, chunk_id, replicate_ids,
#'   posterior_idx[/_occ/_dur], seed, output_part)
build_task_plan <- function(opt, frame = NULL) {
  sim <- opt$simulator

  if (identical(sim$scope$transmission_units, "all")) {
    if (!is.null(frame)) {
      real_tus <- sort(unique(stats::na.omit(frame$tu)))
      has_na <- any(is.na(frame$tu))
      # NA-TU locations (countries absent from HMM cluster file) get their own
      # synthetic bucket so they are simulated with o_tu = NA throughout,
      # matching fit-time behaviour.  Only added when partial NA exists; the
      # all-NA case is already collapsed to `_ALL_` before this point.
      tus <- if (has_na) c(real_tus, "_NA_TU_") else real_tus
    } else {
      tus <- c("hmmW", "hmmCE", "hmmS")
    }
  } else {
    tus <- as.character(sim$scope$transmission_units)
  }

  R <- as.integer(sim$n_replicates)
  cs <- as.integer(sim$parallel$chunk_size)
  chunks <- split(seq_len(R), ceiling(seq_len(R) / cs))

  # Posterior-draw strategy: only 'paired' exists (every track's draw index =
  # the replicate id); the size track is retired from the simulator.
  strategy <- sim$posterior_strategy %||% "paired"
  if (!identical(strategy, "paired")) {
    stop(sprintf(
      paste0(
        "build_task_plan: simulator.posterior_strategy must be ",
        "'paired'; got '%s' (the 'independent' strategy permuted ",
        "the retired size track's draws)."
      ),
      strategy
    ))
  }

  rows <- list()
  i <- 0L
  for (tu in tus) {
    for (k in seq_along(chunks)) {
      i <- i + 1L
      ri <- chunks[[k]]
      pidx_occ <- ri
      # Duration draw index: paired with occurrence (= rep id).
      pidx_dur <- pidx_occ
      rows[[i]] <- list(
        task_id = sprintf("task_%04d", i),
        tu = tu,
        chunk_id = sprintf("chunk_%04d", k),
        replicate_ids = ri,
        posterior_idx = pidx_occ, # backward-compat alias (= occ)
        posterior_idx_occ = pidx_occ,
        posterior_idx_dur = pidx_dur,
        # Stride the base seed by R so set.seed(seed + rep_id) inside
        # simulate_one_tu() yields disjoint (seed + rep_id) ranges per task —
        # the previous `+ (i - 1L)` collided across (task, replicate) pairs.
        seed = as.integer(sim$seed) + (i - 1L) * R,
        output_part = file.path(
          make_run_file(opt, kind = "sim_parts_dir"),
          paste0("tu=", tu),
          sprintf("chunk_%04d", k)
        )
      )
    }
  }
  rows
}

#' Persist the per-task YAMLs under tasks/. Used by both the local furrr path
#' and the SLURM submit script. Returns the vector of written paths.
write_task_yamls <- function(plan, tasks_dir) {
  paths <- character(length(plan))
  for (i in seq_along(plan)) {
    paths[i] <- file.path(tasks_dir, paste0(plan[[i]]$task_id, ".yml"))
    yaml::write_yaml(plan[[i]], paths[i])
  }
  paths
}

# ---------------------------------------------------------------------------
# build_obs_panel — the filter's observation panel, one row per (gid, month)
# ---------------------------------------------------------------------------
# The one observed panel read by the bpfilter (08 smoothing, 07b forecast)
# and the observed-side aggregations (08c, 09):
#   obs_Y = max(outbreak) per (gid, t_start).
# With `mon_dat`, adds `cases_month`, aggregated with the filter's footprint
# rule (spatpomp_core.R::.aggregate_footprint_counts(), call-time dependency).
# `cases_gap` months stay NA. Episode months after onset have no survival
# row, so they carry obs_Y = NA and an observed `cases_month`.
# `zero_fill = TRUE` writes cases_month = 0 at every non-onset survival row
# with no monthly case row ("observed" = has a survival-panel row).
#
# @param surv_dat  survival panel (gid, t_start, outbreak).
# @param mon_dat   NULL, or the monthly case frame (gid, month, cases_month,
#   [admin_level], [cases_gap]).
# @param zero_fill logical; see above. Ignored without `mon_dat`.
# @return tibble(gid, t_start, obs_Y[, cases_month]).
build_obs_panel <- function(surv_dat, mon_dat = NULL, zero_fill = FALSE) {
  obs <- surv_dat %>%
    dplyr::select(gid, t_start, outbreak) %>%
    dplyr::mutate(
      gid = as.character(gid),
      t_start = as.Date(t_start),
      obs_Y = as.integer(outbreak)
    ) %>%
    dplyr::group_by(gid, t_start) %>%
    dplyr::summarise(obs_Y = max(obs_Y, na.rm = TRUE), .groups = "drop")
  if (is.null(mon_dat)) {
    return(obs)
  }

  cm_in <- mon_dat %>%
    dplyr::mutate(
      gid = as.character(gid),
      t_start = as.Date(month),
      cases_month = if ("cases_gap" %in% names(mon_dat)) {
        ifelse(cases_gap %in% TRUE, NA_real_, as.numeric(cases_month))
      } else {
        as.numeric(cases_month)
      }
    )
  cm_in <- cm_in[,
    intersect(c("gid", "t_start", "admin_level", "cases_month"), names(cm_in)),
    drop = FALSE
  ]
  cm <- .aggregate_footprint_counts(
    cm_in,
    "cases_month",
    "cases_month",
    keys = c("gid", "t_start"),
    keep_na = TRUE
  )
  out <- dplyr::full_join(obs, cm, by = c("gid", "t_start"))
  if (isTRUE(zero_fill)) {
    # Only non-onset panel rows: an observed onset month without a monthly case
    # row is "magnitude unobserved" (NA), never "no episode".
    no_onset <- !is.na(out$obs_Y) & out$obs_Y == 0L
    no_case <- !(paste(out$gid, out$t_start) %in% paste(cm$gid, cm$t_start))
    out$cases_month[no_onset & no_case] <- 0
  }
  dplyr::arrange(out, gid, t_start)
}

# ---------------------------------------------------------------------------
# build_coverage_mask lives in utils/coverage.R, the project's one
# surveillance-coverage definition; source
# coverage.R next to this module where it is called.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# aggregate_partitions — merge all successful partitions into merged/
# ---------------------------------------------------------------------------

#' Aggregate simulation partition outputs into a single merged directory.
#'
#' Recursively scans `parts_dir` for `_SUCCESS` sentinel files. For each
#' sentinel, loads the sibling `event_table.rds` and `monthly_summary.rds`.
#' Missing siblings produce a warning and the partition is skipped.
#' The merged outputs are written to `merged_dir/`.
#'
#' @param parts_dir root directory containing partition subdirectories.
#' @param merged_dir target directory for merged outputs (created if absent).
#' @return invisible list(event_table, monthly_summary).
aggregate_partitions <- function(
  parts_dir,
  merged_dir,
  tus = NULL,
  out_fmt = "qs",
  return_long = TRUE
) {
  # When tus is supplied (a character vector of TU labels), only scan the
  # corresponding `tu=<TU>` subdirectories. This prevents stale partitions
  # from a previous scope or parameter change from being silently included.
  # Pass tus=NULL to scan everything (backward-compatible; use for "all" scope).
  sentinels <- partition_sentinels(parts_dir, tus) # utils/paths.R
  if (length(sentinels) == 0L) {
    warning(sprintf(
      "aggregate_partitions: no _SUCCESS sentinels found under '%s'.",
      parts_dir
    ))
    dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
    return(invisible(list(
      event_table = tibble::tibble(),
      monthly_summary = tibble::tibble()
    )))
  }

  event_list <- vector("list", length(sentinels))
  monthly_list <- vector("list", length(sentinels))
  n_loaded <- 0L

  ext <- sim_io_ext(out_fmt)
  for (i in seq_along(sentinels)) {
    part_dir <- dirname(sentinels[[i]])
    et_path <- file.path(part_dir, paste0("event_table", ext))
    ms_path <- file.path(part_dir, paste0("monthly_summary", ext))

    ok <- TRUE
    if (!file.exists(et_path)) {
      warning(sprintf(
        "aggregate_partitions: missing event_table%s in '%s' — skipping.",
        ext,
        part_dir
      ))
      ok <- FALSE
    }
    if (!file.exists(ms_path)) {
      warning(sprintf(
        "aggregate_partitions: missing monthly_summary%s in '%s' — skipping.",
        ext,
        part_dir
      ))
      ok <- FALSE
    }
    if (!ok) {
      next
    }

    event_list[[i]] <- sim_io_read(et_path, out_fmt)
    # Compact (format-2) parts stay compact: the merged file holds their blocks.
    monthly_list[[i]] <- sim_io_read(ms_path, out_fmt, expand = FALSE)
    n_loaded <- n_loaded + 1L
  }

  event_table <- dplyr::bind_rows(event_list)
  monthly_summary <- sim_monthly_bind(monthly_list)
  rm(monthly_list)

  dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
  sim_io_save(
    event_table,
    file.path(merged_dir, paste0("event_table", ext)),
    out_fmt
  )
  sim_monthly_save(
    monthly_summary,
    file.path(merged_dir, paste0("monthly_summary", ext)),
    out_fmt
  )

  n_rows <- if (is_sim_monthly_compact(monthly_summary)) {
    sum(vapply(monthly_summary$blocks, .sim_block_n, numeric(1)))
  } else {
    nrow(monthly_summary)
  }
  message(sprintf(
    "[aggregate_partitions] %d/%d partitions loaded; %d events; %d monthly rows.",
    n_loaded,
    length(sentinels),
    nrow(event_table),
    n_rows
  ))

  # return_long = FALSE: the merged monthly summary as written (compact when
  # the parts are), not expanded -- for callers that only want the file.
  if (isTRUE(return_long) && is_sim_monthly_compact(monthly_summary)) {
    monthly_summary <- sim_monthly_expand(monthly_summary)
  }
  invisible(list(event_table = event_table, monthly_summary = monthly_summary))
}

# ---------------------------------------------------------------------------
# combine_group_runs
#
# Merge the per-group `merged/` outputs into one full-Africa dataset, cut to
# `cutoff_month` (africa_min_horizon: the last month every group covers).
#   group_merged_dirs : named list group -> merged/ directory
#   combined_dir      : output directory (event_table + monthly_summary)
#   cutoff_month      : last month kept (inclusive), Date or "YYYY-MM-DD"
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Per-country censoring + manifest v2 readers
# ---------------------------------------------------------------------------
# A forward run keyed by a partition unit may mix fit groups with different
# lag-defined ("firm") horizons; 07b simulates the unit to the MAX member
# horizon and writes it raw. Every reader censors each country at its own firm
# horizon through censor_sim_by_country(), so the rule lives in one place.
# Call-time dependency: gid.R (flatten_gid()).

#' ISO3 of each gid (its "AAA_" prefix); NA where there is none.
gid_iso3 <- function(gid) {
  stringr::str_extract(flatten_gid(as.character(gid)), "^[A-Z]{3}")
}

#' censor_sim_by_country
#'
#' @param df              monthly summary (month_col = "month") or event table
#'                        (kind = "events", month_col = "start_month").
#' @param horizon_by_iso3 named character/Date: ISO3 -> last month kept.
#' @param month_col       the month column compared with the horizon.
#' @param keep_beyond     FALSE: drop rows past the country's horizon. TRUE:
#'                        keep them, flagged `beyond_firm_horizon`.
#' @param kind            "monthly" | "events". Events also get
#'                        `right_censored` = started by the horizon but ending
#'                        after it (`end_month` is left untouched).
#' @param country_col     column holding ISO3 when the frame has no gid
#'                        (NULL: derived from `gid`).
#' @return df, filtered or flagged. Errors on an ISO3 with no horizon.
censor_sim_by_country <- function(
  df,
  horizon_by_iso3,
  month_col = "month",
  keep_beyond = FALSE,
  kind = c("monthly", "events"),
  country_col = NULL
) {
  kind <- match.arg(kind)
  if (!nrow(df)) {
    return(df)
  }
  # ISO3 once per unique gid (the regex per row was most of a hindcast run's censor time)
  iso <- if (!is.null(country_col)) {
    as.character(df[[country_col]])
  } else {
    local({
      g <- as.character(df$gid)
      u <- unique(g)
      gid_iso3(u)[match(g, u)]
    })
  }
  hz <- as.Date(unlist(horizon_by_iso3))
  names(hz) <- names(unlist(horizon_by_iso3))
  miss <- sort(unique(iso[is.na(iso) | !iso %in% names(hz)]))
  if (length(miss)) {
    stop(
      sprintf(
        "censor_sim_by_country: no horizon for %d country(ies): %s",
        length(miss),
        paste(head(miss, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  m <- df[[month_col]]
  m <- if (inherits(m, "Date")) {
    m
  } else {
    local({
      u <- unique(m)
      as.Date(u)[match(m, u)]
    })
  }
  h <- unname(hz[iso])
  beyond <- m > h
  if (identical(kind, "events") && "end_month" %in% names(df)) {
    em <- as.Date(df$end_month)
    df$right_censored <- !beyond & !is.na(em) & em > h
  }
  if (isTRUE(keep_beyond)) {
    df$beyond_firm_horizon <- beyond
    df
  } else {
    df[!beyond, , drop = FALSE]
  }
}

#' manifest_units
#'
#' The simulation units of a forecast/hindcast manifest. v2 (registry mode):
#' `man$units`; v1: one unit per fit group (members = iso3, groups = the
#' group, the group's run fields). Each element carries members, groups,
#' horizon_end (forecast) and the run fields as in the manifest.
manifest_units <- function(man) {
  if (length(man$units)) {
    return(lapply(man$units, function(u) {
      u$members <- as.character(unlist(u$members))
      u$groups <- as.character(unlist(u$groups))
      u
    }))
  }
  lapply(stats::setNames(nm = names(man$groups)), function(g) {
    x <- man$groups[[g]]
    x$members <- as.character(unlist(x$iso3))
    x$groups <- g
    x
  })
}

#' manifest_runs
#'
#' Per unit, the merged dir + config of one arm of a FORECAST manifest:
#' named list unit -> list(config, merged_dir). `arm` "factual" or
#' "counterfactual" (with `analog_set`, default the primary set; the flat
#' `counterfactual_*` aliases cover the primary set of older manifests).
manifest_runs <- function(
  man,
  arm = c("factual", "counterfactual"),
  analog_set = NULL
) {
  arm <- match.arg(arm)
  units <- manifest_units(man)
  out <- lapply(units, function(u) {
    if (identical(arm, "factual")) {
      return(list(config = u$factual_config, merged_dir = u$factual_merged_dir))
    }
    aset <- analog_set %||%
      (u$primary_analog_set %||% man$primary_analog_set %||% "all")
    cs <- u$counterfactual_sets[[aset]]
    if (!is.null(cs)) {
      return(list(config = cs$config, merged_dir = cs$merged_dir))
    }
    if (identical(aset, u$primary_analog_set %||% man$primary_analog_set)) {
      return(list(
        config = u$counterfactual_config,
        merged_dir = u$counterfactual_merged_dir
      ))
    }
    NULL
  })
  out[!vapply(out, is.null, logical(1))]
}

#' event_manifest_view
#'
#' One event of an event-attribution manifest
#' (utils/make_event_attribution_configs.R) shaped as a v2 FORECAST manifest:
#' groups (with each group's window end as `horizon_end`), units (factual +
#' counterfactual runs, both launching at the event launch) and countries, so
#' every manifest_*() reader and 09's arm loaders read it unchanged.
read_event_manifest <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(
      "read_event_manifest: event manifest not found: '",
      path %||% "",
      "'",
      call. = FALSE
    )
  }
  m <- yaml::read_yaml(path)
  if (!identical(m$kind, "event_attribution") || !length(m$events)) {
    stop(
      "read_event_manifest: ",
      path,
      " is not an event-attribution manifest with events.",
      call. = FALSE
    )
  }
  m
}
event_manifest_view <- function(em, event_id) {
  e <- em$events[[event_id]]
  if (is.null(e)) {
    stop(
      "event_manifest_view: no event '",
      event_id,
      "' in the manifest.",
      call. = FALSE
    )
  }
  list(
    schema = 2L,
    kind = "event_attribution_view",
    event = event_id,
    grouping = em$grouping,
    grouping_digest = em$grouping_digest,
    partition_by = em$partition_by,
    partition_source = em$partition_source,
    analog_sets = list(em$analog_set),
    primary_analog_set = em$analog_set,
    attribution_start = e$launch,
    africa_min_horizon = min(vapply(e$groups, `[[`, "", "horizon_end")),
    africa_max_horizon = e$africa_max_horizon,
    groups = e$groups,
    units = e$units,
    countries = e$countries
  )
}

#' event_sub_windows_table
#'
#' The per-winter sub-windows of one event of an event-attribution manifest
#' (schema 3; utils/make_event_attribution_configs.R): one row
#' per sub-window x key, at the country level when the manifest has per-country
#' season ends (`countries`), else per fit group. NULL for an event without
#' sub-windows (every schema-2 manifest).
#' @return NULL or tibble(sub_window, k, peak, peak_season_year, level
#'   ("country" | "group"), key, start, end (Dates))
event_sub_windows_table <- function(em, event_id) {
  e <- em$events[[event_id]]
  if (is.null(e)) {
    stop(
      "event_sub_windows_table: no event '",
      event_id,
      "' in the manifest.",
      call. = FALSE
    )
  }
  sw <- e[["sub_windows"]]
  if (!length(sw)) {
    return(NULL)
  }
  dplyr::bind_rows(lapply(seq_along(sw), function(k) {
    w <- sw[[k]]
    lev <- if (length(w[["countries"]])) "country" else "group"
    b <- if (identical(lev, "country")) w[["countries"]] else w[["groups"]]
    tibble::tibble(
      sub_window = as.character(w[["id"]]),
      k = k,
      peak = as.Date(w[["peak"]]),
      peak_season_year = as.integer(w[["peak_season_year"]]),
      level = lev,
      key = names(b),
      start = as.Date(unname(vapply(
        b,
        function(x) as.character(x[["start"]]),
        ""
      ))),
      end = as.Date(unname(vapply(b, function(x) as.character(x[["end"]]), "")))
    )
  }))
}

#' assert_counterfactual_cache_identity
#'
#' The counterfactual runs a results job reads must have drawn
#' their ENSO from the cache the manifest records. `expected` is the manifest's
#' record (forecast manifest `enso_counterfactual`, or an event's
#' list(cache_path, digest) from `lenssen_cache` / `lenssen_digest`); NULL =
#' a legacy manifest, nothing checked. Stops when the cache file
#' is missing, when its stored content digest differs, or when a run's
#' partition_meta.yml (next to its merged/ dir; 07b spatPomp) records another
#' digest or none. Runs without a partition_meta.yml here are reported, not
#' failed (outputs pulled without it).
#' @param runs list of manifest_runs() entries (config, merged_dir)
#' @return invisibly the number of runs checked
#' The partition_meta.yml of a run, from its merged/ dir: 07b spatPomp writes it
#' next to the run's parts (<sim root>/<run_tag>/), i.e. beside or one level
#' above merged/. NA when neither exists here.
partition_meta_path <- function(merged_dir) {
  md <- as.character(merged_dir %||% "")
  cand <- unique(c(
    file.path(dirname(md), "partition_meta.yml"),
    file.path(dirname(dirname(md)), "partition_meta.yml")
  ))
  cand[file.exists(cand)][1L]
}

assert_counterfactual_cache_identity <- function(
  expected,
  runs = list(),
  label = "[09]"
) {
  if (is.null(expected) || is.null(expected$digest)) {
    return(invisible(0L))
  }
  path <- as.character(expected$cache_path %||% "")
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      sprintf(
        "%s counterfactual cache %s (digest %s) is missing.",
        label,
        path,
        expected$digest
      ),
      call. = FALSE
    )
  }
  have <- readRDS(path)$meta$digest %||% NA_character_
  if (!identical(have, expected$digest)) {
    stop(
      sprintf(
        "%s counterfactual cache %s has digest %s, the manifest records %s: it was rebuilt after the runs.",
        label,
        path,
        have,
        expected$digest
      ),
      call. = FALSE
    )
  }
  n <- 0L
  absent <- character(0)
  for (r in runs) {
    md <- r$merged_dir %||% ""
    pm_path <- partition_meta_path(md)
    if (is.na(pm_path)) {
      absent <- c(absent, md)
      next
    }
    pm_d <- (yaml::read_yaml(pm_path)$enso_counterfactual %||% list())$digest
    if (!identical(pm_d, expected$digest)) {
      stop(
        sprintf(
          "%s counterfactual run %s drew its ENSO from cache digest %s, not the manifest's %s (%s).",
          label,
          dirname(pm_path),
          pm_d %||% "<none: a legacy cache>",
          expected$digest,
          path
        ),
        call. = FALSE
      )
    }
    n <- n + 1L
  }
  if (length(absent)) {
    message(sprintf(
      "%s counterfactual cache identity: %d run(s) without partition_meta.yml here, not checked.",
      label,
      length(absent)
    ))
  }
  invisible(n)
}

#' assert_panel_covers_horizon
#'
#' Config-time check of the closed-horizon convention (utils/base.R): the
#' climate panel a forward arm reads must contain the horizon's
#' last month, else 07b would drop that month (build_climate_panel_for_tu()
#' only warns at run time). `panel` is a panel data frame or an RDS path; a
#' path that does not exist here (the panel lives on the cluster) is reported
#' and not checked.
#' @return TRUE when checked and covered, NA when not checked (invisibly)
assert_panel_covers_horizon <- function(panel, end, label) {
  end_month <- as.Date(format(as.Date(end), "%Y-%m-01"))
  if (is.character(panel)) {
    if (!length(panel) || !nzchar(panel) || !file.exists(panel)) {
      message(sprintf(
        "%s: climate panel '%s' not found here; its coverage of %s is not checked.",
        label,
        paste(panel, collapse = ""),
        format(end_month)
      ))
      return(invisible(NA))
    }
    path <- panel
    panel <- readRDS(panel)
  } else {
    path <- "<panel>"
  }
  have <- max(as.Date(panel$t_start))
  if (have < end_month) {
    stop(
      sprintf(
        paste0(
          "%s: climate panel %s ends %s, before the closed horizon's last month %s ",
          "(07b would not simulate it). Rebuild the panel to %s or narrow the horizon."
        ),
        label,
        path,
        format(have),
        format(end_month),
        format(end_month)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' manifest_last_firm_oni
#'
#' The last ONI month a forecast manifest's firm horizons are tied to
#' (horizon_end_g = this month + the group's minimum deployed lag): every
#' continuation lead past a firm horizon is counted from it. A manifest with
#' `firm_oni_end` (the last full ONI month) uses it; older manifests
#' tie the horizons to `last_cfsv2_month`. The one accessor for it: never read
#' `last_cfsv2_month` directly for leads or horizons.
#' @return character "YYYY-MM-DD" (NULL when the manifest has neither)
manifest_last_firm_oni <- function(man) {
  man[["firm_oni_end"]] %||% man[["last_cfsv2_month"]]
}

#' manifest_country_groups
#'
#' ISO3 -> fit group, from `countries` (v2) or the groups' iso3 lists (v1).
manifest_country_groups <- function(man) {
  if (length(man$countries)) {
    return(vapply(
      man$countries,
      function(x) as.character(x$group),
      character(1)
    ))
  }
  unlist(lapply(names(man$groups), function(g) {
    stats::setNames(
      rep(g, length(unlist(man$groups[[g]]$iso3))),
      unlist(man$groups[[g]]$iso3)
    )
  }))
}

#' manifest_country_horizons
#'
#' ISO3 -> horizon (character "YYYY-MM-DD"). which = "firm": the country's fit
#' group's lag-defined horizon; "unit": the horizon its unit was simulated to
#' (equal to firm in v1). For a HINDCAST manifest pass `year` (launch year):
#' firm = groups.<g>.years.<year>.horizon_end, unit = units.<u>.years.<year>.
#' A country's own entry wins when the manifest has one (per-country season
#' ends): countries.<iso>.firm_horizon_end / unit_horizon_end, or for a
#' hindcast countries.<iso>.years.<year>.*; country entries without such a
#' field read the group's horizon.
manifest_country_horizons <- function(
  man,
  which = c("firm", "unit"),
  year = NULL
) {
  which <- match.arg(which)
  cg <- manifest_country_groups(man)
  yr <- if (is.null(year)) NULL else as.character(year)
  grp_h <- function(g) {
    x <- man$groups[[g]]
    as.character(if (is.null(yr)) x$horizon_end else x$years[[yr]]$horizon_end)
  }
  field <- paste0(which, "_horizon_end")
  # Season-window hindcast (make_hindcast_configs.R --season_window): `year` is
  # a run key ("<year>-<MM>") and only the countries scored in that run have an
  # entry under it -- no group fallback.
  if (identical(man$season_end, "season_window") && !is.null(yr)) {
    own <- lapply(man$countries, function(x) {
      ((x[["years"]] %||% list())[[yr]] %||% list())[[field]]
    })
    own <- Filter(Negate(is.null), own)
    return(vapply(own, as.character, character(1)))
  }
  own_h <- function(iso) {
    x <- man$countries[[iso]]
    if (is.null(x)) {
      return(NULL)
    }
    v <- if (is.null(yr)) {
      x[[field]]
    } else {
      ((x[["years"]] %||% list())[[yr]] %||% list())[[field]]
    }
    if (is.null(v)) NULL else as.character(v)
  }
  if (identical(which, "firm") || !length(man$units)) {
    return(vapply(
      names(cg),
      function(iso) own_h(iso) %||% grp_h(cg[[iso]]),
      character(1)
    ))
  }
  unit_of <- vapply(
    man$countries,
    function(x) as.character(x$unit),
    character(1)
  )
  vapply(
    names(cg),
    function(iso) {
      own <- own_h(iso)
      if (!is.null(own)) {
        return(own)
      }
      u <- man$units[[unit_of[[iso]]]]
      as.character(
        if (is.null(yr)) u$horizon_end else u$years[[yr]]$horizon_end
      )
    },
    character(1)
  )
}

#' manifest_has_country_horizons
#'
#' TRUE when a manifest's countries carry their own firm horizons (per-country
#' season ends) -- for a hindcast, under `year` -- so a reader must censor each
#' country even when the manifest has no units (v1 layout, unit = fit group).
manifest_has_country_horizons <- function(man, year = NULL) {
  if (!length(man$countries)) {
    return(FALSE)
  }
  yr <- if (is.null(year)) NULL else as.character(year)
  any(vapply(
    man$countries,
    function(x) {
      v <- if (is.null(yr)) {
        x[["firm_horizon_end"]]
      } else {
        ((x[["years"]] %||% list())[[yr]] %||% list())[["firm_horizon_end"]]
      }
      !is.null(v)
    },
    logical(1)
  ))
}

# ---------------------------------------------------------------------------
# Per-arm launch (handoff) months. The two forecast arms may hand
# off from the same smoothing cache at different months (the counterfactual at
# the Lenssen analog anchor, the factual at a later smoothed month); readers
# take each arm's launch from the manifest instead of assuming the factual
# horizon.start is shared.
# ---------------------------------------------------------------------------

#' assert_arm_launch
#'
#' The 07b invariant for one forward-run config (07b refuses a smoothing handoff
#' that differs from horizon.start): simulator.horizon.start equals
#' filter.forecast_window.start and, when set, filter.handoff_month. The arm's
#' launch IS its handoff, so one date describes it. Shared by
#' make_forecast_configs.R (generation time) and deploy_core.R's template check.
#' @param sim   the `simulator` block of a config or template.
#' @param label names the config in errors.
#' @return the launch month (Date).
assert_arm_launch <- function(sim, label) {
  hs <- (sim$horizon %||% list())$start
  if (is.null(hs) || !nzchar(as.character(hs))) {
    stop(
      sprintf("%s: simulator.horizon.start is not set.", label),
      call. = FALSE
    )
  }
  launch <- as.Date(as.character(hs))
  if (is.na(launch)) {
    stop(
      sprintf("%s: simulator.horizon.start '%s' is not a date.", label, hs),
      call. = FALSE
    )
  }
  flt <- sim$filter %||% list()
  fw <- (flt$forecast_window %||% list())$start
  if (!is.null(fw) && as.Date(as.character(fw)) != launch) {
    stop(
      sprintf(
        paste0(
          "%s: filter.forecast_window.start (%s) differs from horizon.start (%s); ",
          "07b reads the forecast from horizon.start, so they must be equal."
        ),
        label,
        fw,
        format(launch)
      ),
      call. = FALSE
    )
  }
  ho <- flt$handoff_month
  if (
    !is.null(ho) &&
      nzchar(as.character(ho)) &&
      as.Date(as.character(ho)) != launch
  ) {
    stop(
      sprintf(
        paste0(
          "%s: filter.handoff_month (%s) differs from horizon.start (%s); ",
          "07b refuses a handoff that is not the forecast start."
        ),
        label,
        ho,
        format(launch)
      ),
      call. = FALSE
    )
  }
  # Both ends name one CLOSED forward horizon (utils/base.R):
  # 07b simulates forecast_window, the readers censor at horizon.end.
  he <- (sim$horizon %||% list())$end
  fwe <- (flt$forecast_window %||% list())$end
  if (
    !is.null(he) &&
      !is.null(fwe) &&
      nzchar(as.character(he)) &&
      nzchar(as.character(fwe)) &&
      as.Date(as.character(fwe)) != as.Date(as.character(he))
  ) {
    stop(
      sprintf(
        paste0(
          "%s: filter.forecast_window.end (%s) differs from horizon.end (%s); ",
          "both are the closed forward horizon, so they must be equal."
        ),
        label,
        fwe,
        he
      ),
      call. = FALSE
    )
  }
  launch
}

#' manifest_arm_launches
#'
#' Every run of a FORECAST manifest with its launch month: one row per
#' (entry, arm, analog_set). `entry` is the fit group (v1) or the partition
#' unit (v2); `groups` (list column) the fit groups it holds. The launch comes
#' from the manifest's `factual_launch` / `counterfactual_sets.<set>.launch`
#' keys (`source = "manifest"`), or -- manifests without those keys --
#' from each config's simulator.horizon.start (`source = "config"`).
manifest_arm_launches <- function(man) {
  units <- manifest_units(man)
  rows <- list()
  cfg_launch <- function(path, label) {
    if (is.null(path) || !file.exists(path)) {
      stop(
        sprintf(
          "manifest_arm_launches: %s has no launch key and its config '%s' is missing.",
          label,
          path %||% "<none>"
        ),
        call. = FALSE
      )
    }
    assert_arm_launch(yaml::read_yaml(path)$simulator, path)
  }
  one <- function(entry, groups, arm, aset, launch_key, cfg) {
    src <- if (!is.null(launch_key) && nzchar(as.character(launch_key))) {
      "manifest"
    } else {
      "config"
    }
    l <- if (src == "manifest") {
      as.Date(as.character(launch_key))
    } else {
      cfg_launch(
        cfg,
        sprintf(
          "entry '%s' %s%s",
          entry,
          arm,
          if (is.na(aset)) "" else paste0("[", aset, "]")
        )
      )
    }
    tibble::tibble(
      entry = entry,
      groups = list(groups),
      arm = arm,
      analog_set = aset,
      launch = l,
      source = src
    )
  }
  for (en in names(units)) {
    u <- units[[en]]
    rows[[length(rows) + 1L]] <- one(
      en,
      u$groups,
      "factual",
      NA_character_,
      u$factual_launch,
      u$factual_config
    )
    primary <- as.character(
      u$primary_analog_set %||% man$primary_analog_set %||% "all"
    )
    sets <- names(u$counterfactual_sets)
    if (!length(sets) && !is.null(u$counterfactual_config)) {
      sets <- primary
    }
    for (a in sets) {
      cs <- u$counterfactual_sets[[a]] %||% list()
      cfg <- cs$config %||%
        (if (identical(a, primary)) u$counterfactual_config else NULL)
      key <- cs$launch %||%
        (if (identical(a, primary)) u$counterfactual_launch else NULL)
      rows[[length(rows) + 1L]] <- one(
        en,
        u$groups,
        "counterfactual",
        a,
        key,
        cfg
      )
    }
  }
  dplyr::bind_rows(rows)
}

#' manifest_arm_starts
#'
#' The uniform launch of each arm over the entries holding any of `groups`
#' (counterfactual restricted to `sets`). Every entry of one arm must launch at
#' the same month (the readers window all groups together); a counterfactual
#' launching after the factual arm is refused (the factual series can be
#' stitched from the smoothing back to an earlier counterfactual launch, never
#' the reverse).
#' @return list(factual = Date, counterfactual = Date)
manifest_arm_starts <- function(man, groups = NULL, sets = NULL) {
  al <- manifest_arm_launches(man)
  if (!is.null(groups)) {
    al <- al[
      vapply(al$groups, function(g) any(g %in% groups), logical(1)),
      ,
      drop = FALSE
    ]
  }
  if (!is.null(sets)) {
    al <- al[al$arm == "factual" | al$analog_set %in% sets, , drop = FALSE]
  }
  if (!nrow(al)) {
    stop(
      "manifest_arm_starts: no manifest run matches the requested groups/sets.",
      call. = FALSE
    )
  }
  out <- list()
  for (arm in c("factual", "counterfactual")) {
    x <- al[al$arm == arm, , drop = FALSE]
    if (!nrow(x)) {
      stop(
        sprintf("manifest_arm_starts: no %s run in the manifest.", arm),
        call. = FALSE
      )
    }
    if (length(unique(x$launch)) != 1L) {
      stop(
        sprintf(
          "manifest_arm_starts: the %s arm does not launch at one month: %s",
          arm,
          paste(
            sprintf(
              "%s%s=%s",
              x$entry,
              ifelse(is.na(x$analog_set), "", paste0("[", x$analog_set, "]")),
              format(x$launch)
            ),
            collapse = ", "
          )
        ),
        call. = FALSE
      )
    }
    out[[arm]] <- x$launch[1L]
  }
  if (out$counterfactual > out$factual) {
    stop(
      sprintf(
        paste0(
          "manifest_arm_starts: the counterfactual arm launches at %s, after the ",
          "factual arm (%s); only an earlier (or equal) counterfactual launch is ",
          "supported."
        ),
        format(out$counterfactual),
        format(out$factual)
      ),
      call. = FALSE
    )
  }
  out
}

#' manifest_attrib_start
#'
#' Start of the attribution window. `alignment = "stitch"`: the earlier launch
#' (the factual series is the smoothing output up to its own launch);
#' `"common_start"`: the later launch, the first month both arms simulate.
#' Equal launches give the same month either way.
manifest_attrib_start <- function(
  man,
  groups = NULL,
  sets = NULL,
  alignment = c("stitch", "common_start")
) {
  alignment <- match.arg(alignment)
  s <- manifest_arm_starts(man, groups, sets)
  if (identical(alignment, "stitch")) {
    min(s$factual, s$counterfactual)
  } else {
    max(s$factual, s$counterfactual)
  }
}

#' The factual arm's launch (its uniform horizon.start across `groups`).
manifest_factual_start <- function(man, groups = NULL, sets = NULL) {
  manifest_arm_starts(man, groups, sets)$factual
}

#' Label of the factual / counterfactual ratio: "ENSO" when both arms launch at
#' the same month, else "ENSO + post-handoff information" (the factual arm is
#' conditioned on smoothed observations the counterfactual never saw). Read by
#' 09 (stored on the artefacts), 10 and 07d; never re-derived at a plot site.
attribution_ratio_label <- function(factual_start, counterfactual_start) {
  if (as.Date(factual_start) == as.Date(counterfactual_start)) {
    "ENSO"
  } else {
    "ENSO + post-handoff information"
  }
}

#' The smoothing run's `parts/` dir that sits beside its `cache/` dir (the
#' forward configs' filter.smoothing_input). deploy_core.R::deploy_smoothing_dir()
#' builds both from the historical run tag; this recovers one from the other.
smoothing_parts_dir_of <- function(smoothing_input) {
  si <- sub("/+$", "", as.character(smoothing_input %||% ""))
  if (!nzchar(si) || basename(si) != "cache") {
    stop(
      sprintf(
        "smoothing_parts_dir_of: smoothing_input '%s' does not end in /cache.",
        smoothing_input %||% ""
      ),
      call. = FALSE
    )
  }
  file.path(dirname(si), "parts")
}

combine_group_runs <- function(
  group_merged_dirs,
  combined_dir,
  cutoff_month,
  out_fmt = "qs",
  horizon_by_country = NULL
) {
  ext <- sim_io_ext(out_fmt)
  cutoff <- as.Date(cutoff_month)
  et_list <- list()
  ms_list <- list()
  for (g in names(group_merged_dirs)) {
    md <- group_merged_dirs[[g]]
    et_path <- file.path(md, paste0("event_table", ext))
    ms_path <- file.path(md, paste0("monthly_summary", ext))
    if (!file.exists(et_path) || !file.exists(ms_path)) {
      warning(sprintf(
        "combine_group_runs: missing merged outputs for group '%s' under '%s' — skipping.",
        g,
        md
      ))
      next
    }
    et <- sim_io_read(et_path, out_fmt)
    ms <- sim_io_read(ms_path, out_fmt)
    et_list[[g]] <- dplyr::filter(et, as.Date(start_month) <= cutoff)
    ms_list[[g]] <- dplyr::filter(ms, as.Date(month) <= cutoff)
    rm(ms)
    # Units mixing fit groups: each country also stops at its own firm horizon
    # (a no-op when every horizon is >= the common cutoff).
    if (!is.null(horizon_by_country)) {
      et_list[[g]] <- censor_sim_by_country(
        et_list[[g]],
        horizon_by_country,
        month_col = "start_month",
        kind = "events"
      )
      ms_list[[g]] <- censor_sim_by_country(ms_list[[g]], horizon_by_country)
    }
  }
  if (!length(ms_list)) {
    stop(
      "combine_group_runs: no group merged outputs found under the supplied dirs."
    )
  }

  event_table <- dplyr::bind_rows(et_list)
  monthly_summary <- dplyr::bind_rows(ms_list)

  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  sim_io_save(
    event_table,
    file.path(combined_dir, paste0("event_table", ext)),
    out_fmt
  )
  # Format 2 when every group's table is still a dense grid (one block each;
  # expands to the same bind_rows()), else the long table.
  cx <- if (
    !identical(out_fmt, "parquet") &&
      !isTRUE(getOption("cholera.sim_monthly_long"))
  ) {
    lapply(ms_list, function(m) {
      c1 <- sim_monthly_compact(m)
      if (!is.null(c1) && identical(sim_monthly_expand(c1), m)) c1
    })
  }
  cx <- if (length(cx) && !any(vapply(cx, is.null, logical(1)))) {
    sim_monthly_bind(unname(cx))
  }
  if (!is.null(cx) && identical(sim_monthly_expand(cx), monthly_summary)) {
    sim_io_save(
      cx,
      file.path(combined_dir, paste0("monthly_summary", ext)),
      out_fmt
    )
  } else {
    sim_io_save(
      monthly_summary,
      file.path(combined_dir, paste0("monthly_summary", ext)),
      out_fmt
    )
  }
  rm(cx)

  message(sprintf(
    "[combine_group_runs] %d group(s) combined; cutoff <= %s; %d events, %d monthly rows -> %s",
    length(ms_list),
    format(cutoff),
    nrow(event_table),
    nrow(monthly_summary),
    combined_dir
  ))
  invisible(list(event_table = event_table, monthly_summary = monthly_summary))
}

# ---------------------------------------------------------------------------
# Hindcast manifest iteration: the group x launch-year loop shared by 07c's
# run_hindcast_skill() and 07d's run_hindcast_figures().
# Call-time dependencies: config.R (make_options_from_config()),
# paths.R (make_run_file()).
# ---------------------------------------------------------------------------

#' read_hindcast_manifest
#'
#' Validate the --manifest / --obs-panel pair of a --hindcast run and read the
#' manifest (`configs/hindcast_horizons.yml`, written by
#' make_hindcast_configs.R). Stops with the same messages both CLIs used.
read_hindcast_manifest <- function(manifest_path, obs_panel_path) {
  if (
    is.null(manifest_path) ||
      !nzchar(manifest_path) ||
      !file.exists(manifest_path)
  ) {
    stop(
      "--hindcast requires --manifest <hindcast_horizons.yml>; not found: ",
      manifest_path
    )
  }
  if (
    is.null(obs_panel_path) ||
      !nzchar(obs_panel_path) ||
      !file.exists(obs_panel_path)
  ) {
    stop("--hindcast requires --obs-panel <observed survival panel RDS>.")
  }
  man <- yaml::read_yaml(manifest_path)
  if (is.null(man$groups) || !length(man$groups)) {
    stop("Manifest has no `groups`: ", manifest_path)
  }
  man
}

#' hindcast_entry_size
#'
#' Bytes of each entry's merged monthly_summary.* (the manifest's `merged_dir`),
#' 0 when the manifest names none or it is not built yet. A proxy for the run's
#' scoring cost, used only to order parallel work.
hindcast_entry_size <- function(entries) {
  vapply(
    entries,
    function(e) {
      md <- e$entry$merged_dir
      if (is.null(md) || !dir.exists(md)) {
        return(0)
      }
      f <- list.files(md, pattern = "^monthly_summary\\.", full.names = TRUE)
      if (length(f)) sum(file.size(f), na.rm = TRUE) else 0
    },
    numeric(1)
  )
}

#' for_each_hindcast_run
#'
#' Iterate a hindcast manifest's (group, launch year) entries: skip an entry
#' whose config is missing, resolve its options from the config, resolve its
#' merged dir (`merged_dir` in the manifest, else the run's spatpomp merged
#' dir), optionally build a missing monthly summary (`ensure_summary`), skip
#' an entry that still has none, and call `fn(run)` for the rest.
#'
#' @param man manifest list (read_hindcast_manifest())
#' @param option_list optparse option list to resolve each config with
#' @param tag log tag, e.g. "[hindcast]" or "[hindcast-fig]"
#' @param fn function(run) called per scored entry; `run` is
#'   list(group, year, entry, config, opt, out_fmt, merged_dir, unit_groups,
#'   firm_horizon_by_country, launch_year, scored_countries). `group` is the run
#'   key (a fit group, or a partition unit in a v2 manifest);
#'   `firm_horizon_by_country` is NULL for v1. `year` is the manifest key (a
#'   launch year, or "<year>-<MM>" in a season-window manifest), `launch_year`
#'   the integer season year; `scored_countries` is NULL except in a
#'   season-window manifest.
#'   Its return values are collected (NULLs dropped).
#' @param announce "before": print the run banner (with a leading blank line)
#'   right after the config check, as 07c did; "after": print it once the
#'   monthly summary is confirmed, as 07d did.
#' @param ensure_summary optional function(opt, merged_dir, out_fmt) called
#'   when the monthly summary is missing (07c aggregates the partitions).
#' @param workers forked processes for the per-run work (base.R::map_workers());
#'   1 (default) is the serial loop.
#' @return list(results, n_ok, n_skip)
for_each_hindcast_run <- function(
  man,
  option_list,
  tag,
  fn,
  announce = c("before", "after"),
  ensure_summary = NULL,
  workers = 1L
) {
  announce <- match.arg(announce)
  banner <- function(g, yr, cfg) {
    message(sprintf(
      "%s%s === group=%s launch_year=%s (%s) ===",
      if (identical(announce, "before")) "\n" else "",
      tag,
      g,
      yr,
      basename(cfg)
    ))
  }
  # Manifest entries in group x year order; a missing config is a skip here.
  entries <- list()
  n_skip <- 0L
  # Runs live per fit group (v1) or per partition unit (v2, registry mode);
  # `group` is the run key either way, `unit_groups` the fit groups it holds.
  run_keys <- if (length(man$units)) man$units else man$groups
  for (g in names(run_keys)) {
    yrs <- run_keys[[g]]$years %||% list()
    for (yr in names(yrs)) {
      ye <- yrs[[yr]]
      cfg <- ye$config
      if (is.null(cfg) || !file.exists(cfg)) {
        warning(sprintf(
          "%s group=%s year=%s: config missing (%s) — skipping.",
          tag,
          g,
          yr,
          cfg %||% "<null>"
        ))
        n_skip <- n_skip + 1L
        next
      }
      entries[[length(entries) + 1L]] <- list(
        group = g,
        year = yr,
        entry = ye,
        config = cfg,
        unit_groups = if (length(man$units)) {
          as.character(unlist(run_keys[[g]]$groups))
        } else {
          g
        }
      )
    }
  }

  # One run: options, (optional) aggregation, `fn`. Runs are independent (each
  # reads and writes only its own run directory), so they go to `workers`
  # forked processes (base.R::map_workers(): input order kept, warnings
  # re-raised here, a failed run stops the loop). workers = 1 is a serial
  # loop.
  one_run <- function(e) {
    g <- e$group
    yr <- e$year
    ye <- e$entry
    cfg <- e$config
    if (identical(announce, "before")) {
      banner(g, yr, cfg)
    }
    opt_run <- make_options_from_config(
      option_list = option_list,
      enforce_options = character(0),
      cli_argv = c("--config", cfg),
      model = NULL
    )
    out_fmt_run <- opt_run$simulator$output$format %||% "qs"
    md <- ye$merged_dir %||%
      make_run_file(opt_run, kind = "spatpomp_merged_dir")
    ms_path <- file.path(md, paste0("monthly_summary", sim_io_ext(out_fmt_run)))

    if (!file.exists(ms_path) && !is.null(ensure_summary)) {
      ensure_summary(opt_run, md, out_fmt_run)
    }
    if (!file.exists(ms_path)) {
      warning(sprintf(
        "%s group=%s year=%s: no monthly_summary at %s — skipping.",
        tag,
        g,
        yr,
        md
      ))
      return(list(ok = FALSE, out = NULL))
    }
    if (identical(announce, "after")) {
      banner(g, yr, cfg)
    }

    # v2 manifests: each member country's FIRM window end for this launch year
    # (readers censor the unit's raw output there); NULL for v1 (unit = group,
    # firm = the run's own horizon) unless its countries carry their own
    # (per-country season ends, make_hindcast_configs.R --season_end
    # by_country_from_panel).
    firm <- if (
      length(man$units) || manifest_has_country_horizons(man, year = yr)
    ) {
      manifest_country_horizons(man, "firm", year = yr)
    } else {
      NULL
    }
    # Season-window manifests: the run key is "<year>-<MM>" (launch_year in the
    # entry) and only `scored` countries are read. The unit's other members get
    # the month before the launch as their horizon, so every reader's censor
    # (censor_sim_by_country()) drops them.
    scored <- NULL
    if (identical(man$season_end, "season_window")) {
      scored <- as.character(unlist(ye$scored))
      pad <- setdiff(as.character(unlist(run_keys[[g]]$members)), names(firm))
      if (length(pad)) {
        firm <- c(
          firm,
          stats::setNames(
            rep(
              format(seq(as.Date(ye$launch), by = "-1 month", length.out = 2L)[
                2L
              ]),
              length(pad)
            ),
            pad
          )
        )
      }
    }
    list(
      ok = TRUE,
      out = fn(list(
        group = g,
        year = yr,
        entry = ye,
        config = cfg,
        opt = opt_run,
        out_fmt = out_fmt_run,
        merged_dir = md,
        unit_groups = e$unit_groups,
        firm_horizon_by_country = firm,
        launch_year = as.integer(
          ye$launch_year %||%
            suppressWarnings(as.integer(yr))
        ),
        scored_countries = scored
      ))
    )
  }
  # Parallel: hand runs out largest first (manifest merged monthly_summary size;
  # map_workers() dispatches dynamically in input order) so a big run does not
  # start last and set the makespan alone; results go back to manifest order.
  ord <- seq_along(entries)
  if (length(entries) > 1L && isTRUE(workers > 1L)) {
    ord <- order(-hindcast_entry_size(entries), ord)
  }
  res <- map_workers(entries[ord], one_run, workers = workers)[order(ord)]
  ok <- vapply(res, `[[`, logical(1), "ok")
  results <- Filter(Negate(is.null), lapply(res[ok], `[[`, "out"))
  list(results = results, n_ok = sum(ok), n_skip = n_skip + sum(!ok))
}
