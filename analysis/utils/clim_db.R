# clim_db.R — read the climate stats / extremes SQLite databases per unit.
#
# Provides:
#   get_stats_db(), get_extrm_db()      open the stats / extremes database
#   get_stats_gid(), get_extremes_gid() rows of one or more gids
#   assert_unique_clim_rows()           refuse duplicated climate rows
#
# Details:
#   Databases are written by 02_preprocess_climate.R
#   (generated_data/climate_stats.sqlite by default).

# Database utils ----------------------------------------------------------

# `path` exists so a rebuild can be staged into a scratch database and diffed
# against the live one before 02_preprocess_climate.R's dbWriteTable(overwrite =
# TRUE) destroys it. Unset, both functions behave exactly as before.
get_stats_db <- function(path = "generated_data/climate_stats.sqlite") {
  stats_db <- DBI::dbConnect(RSQLite::SQLite(), path, extended_types = TRUE)
  stats_db
}

get_extrm_db <- function(path = "generated_data/climate_extremes.sqlite") {
  stats_db <- DBI::dbConnect(RSQLite::SQLite(), path, extended_types = TRUE)
  stats_db
}


get_stats_gid <- function(
  gids,
  cols = NULL,
  vars = c("precip", "temp"),
  path = NULL
) {
  # Connection to DB (`path` NULL = get_stats_db()'s default; set it to read a
  # staged or relocated database)
  .con <- if (is.null(path)) get_stats_db() else get_stats_db(path)

  if (!is.null(cols)) {
    sel_cols <- cols
  } else {
    sel_cols <- purrr::map(vars, ~ get_clim_vars(.)) %>% unlist()
  }

  # climate_stats stores FLATTENED gids ("AGO_4_1"); callers (make_size_data(),
  # make_surv_data()) may pass the raw dotted GADM form ("AGO.4_1") now that the
  # crosswalk gid universe is in play. Query by the flattened form so the
  # lookup actually matches, then map the result's `gid` back to whatever
  # spelling the caller used, so the caller's later `inner_join(by = "gid")`
  # still lines up. flatten_gid() is a no-op on already-flattened gids, so this
  # is safe either way -- same fix as the coverage check in
  # 02_preprocess_climate.R (utils.R flatten_gid() comment there).
  gids_flat <- flatten_gid(gids)
  flat_to_orig <- stats::setNames(gids, gids_flat)

  query <- glue::glue_sql(
    .con = .con,
    "SELECT gid, month, {`sel_cols`*}
      FROM climate_stats
      WHERE gid IN ({gids_flat*});"
  )

  res <- dbGetQuery(conn = .con, statement = query) %>%
    as_tibble() %>%
    mutate(gid = unname(flat_to_orig[gid]))

  DBI::dbDisconnect(.con)

  res
}

get_extremes_gid <- function(
  gids,
  sel_vars = NULL,
  vars = c("precip", "temp"),
  baseline,
  path = NULL
) {
  # Connection to DB (`path` as in get_stats_gid())
  .con <- if (is.null(path)) get_extrm_db() else get_extrm_db(path)

  if (is.null(sel_vars)) {
    sel_vars <- purrr::map(vars, ~ get_clim_vars(.)) %>% unlist()
  }
  tabname <- stringr::str_c("clim_extremes_", baseline)

  # Same flatten-then-map-back fix as get_stats_gid() -- clim_extremes_* stores
  # flattened gids too.
  gids_flat <- flatten_gid(gids)
  flat_to_orig <- stats::setNames(gids, gids_flat)

  query <- glue::glue_sql(
    .con = .con,
    "SELECT gid, month, name, mean, sd, hilo_sd, hilo_q
      FROM {`tabname`}
      WHERE gid IN ({gids_flat*})
      AND name IN ({sel_vars*});"
  )

  res <- dbGetQuery(conn = .con, statement = query) %>%
    as_tibble() %>%
    mutate(gid = unname(flat_to_orig[gid]))

  DBI::dbDisconnect(.con)

  res
}

#' assert_unique_clim_rows
#' Refuse to write a climate table whose key is not unique.
#'
#' This is the guard whose absence let three Ghana gids grow to 33.9M / 33.9M /
#' 19.9M rows against 528 distinct months each -- 99.7% of an 8.5 GB
#' climate_stats.sqlite -- via compile_stats()'s purrr::reduce(full_join)
#' (utils.R:493). Every downstream join then silently multiplied those rows.
#' (Moved verbatim from 02_preprocess_climate.R, which still calls it.)
assert_unique_clim_rows <- function(dat, key, what) {
  dup <- duplicated(dat[, key, drop = FALSE])
  if (!any(dup)) {
    cat(sprintf(
      "---- %s: %d rows, key (%s) unique\n",
      what,
      nrow(dat),
      paste(key, collapse = ", ")
    ))
    return(invisible(TRUE))
  }
  bad <- unique(dat$gid[dup])
  stop(
    what,
    ": ",
    sum(dup),
    " duplicated rows on key (",
    paste(key, collapse = ", "),
    ") across ",
    length(bad),
    " gid(s): ",
    paste(utils::head(bad, 10), collapse = ", "),
    "\n  Refusing to write -- a cartesian blow-up here silently inflates ",
    "every downstream model dataset. Inspect the raw inputs for those gids."
  )
}
