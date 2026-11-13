
#' filterAfrica
#'
#' Filter `dat` to rows whose `country` (ISO3) is geographically in Africa,
#' warning (not silently) about any non-African rows dropped.
#'
#' @param dat data frame/tibble with a `country` column (ISO3 codes).
#' @param africa_iso3 character vector of African ISO3 codes.
#' @param tag character prefix identifying the caller, for the warning message.
filterAfrica <- function(dat, africa_iso3, tag) {
  dropped <- dplyr::filter(dat, !country %in% africa_iso3)
  if (nrow(dropped) > 0) {
    dropped_countries <- sort(unique(dropped$country))
    warning(sprintf(
      paste0("[%s] Dropping %d row(s) from %d non-African country/countries ",
             "(%s) -- these pass the upstream AFR/EMR WHO-region filter but ",
             "are outside Africa."),
      tag, nrow(dropped), length(dropped_countries),
      paste(dropped_countries, collapse = ", ")
    ), call. = FALSE)
  }
  dplyr::filter(dat, country %in% africa_iso3)
}

#' syntheticCountries
#'
#' The fake countries of the synthetic world (00_ingest_synthetic.R), from
#' analysis/groupings/synthetic_countries.csv: ISO 3166 user-assigned codes
#' (XA?) that no real country carries, so adding them to a list of real codes
#' changes nothing for the real data (no real row matches an XA? code).
#'
#' @param path the synthetic country registry
#' @return data.frame(iso3, name, afro_region, ...); zero rows if the file is missing
syntheticCountries <- function(path = "analysis/groupings/synthetic_countries.csv") {
  if (!file.exists(path)) return(data.frame(iso3 = character(0), afro_region = character(0)))
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' africaIso3
#'
#' ISO3 codes of the countries geographically in Africa (giscoR's
#' gisco_countrycode, continent "Africa"; keeps the WHO-EMR-administered North
#' and Horn countries), plus the synthetic countries. The one list
#' 01_preprocess_outbreaks.R and 01b_preprocess_outbreaks_monthly.R filter on
#' (gate B of 01b compares their totals).
#'
#' @return character vector of ISO3 codes
africaIso3 <- function() {
  e <- new.env()
  utils::data("gisco_countrycode", package = "giscoR", envir = e)
  cc <- e$gisco_countrycode
  # which(): giscoR has countries with an NA continent, which `[` would keep as NA.
  unique(c(cc$ISO3_CODE[which(cc$continent == "Africa")], syntheticCountries()$iso3))
}

#' fillSyntheticAfroRegion
#'
#' taxdat::get_AFRO_region() knows only real ISO3 codes, so it leaves the
#' synthetic countries NA, and the smoothing drops NA-region rows. Fills those
#' from the synthetic registry's afro_region; a no-op for real data.
#'
#' @param dat data frame with an `AFRO_region` column
#' @param ctry_col its ISO3 column
fillSyntheticAfroRegion <- function(dat, ctry_col = "country") {
  syn <- syntheticCountries()
  if (!nrow(syn) || !"AFRO_region" %in% names(dat)) return(dat)
  hit <- is.na(dat$AFRO_region) & dat[[ctry_col]] %in% syn$iso3
  if (any(hit)) dat$AFRO_region[hit] <- syn$afro_region[match(dat[[ctry_col]][hit], syn$iso3)]
  dat
}

# Utils -------------------------------------------------------------------

# TZA mainland: the corpus reports a whole-mainland location period under a
# sub-national admin_level label. It is ADM0 by geometry, and is treated as ADM0
# EVERYWHERE, through this one helper (03_prepare_occ.R, 03_prepare_siz.R);
# under `--episode_key gid_season` it is ALSO applied before the section-B
# `admin_level != "ADM0"` filter, so the unit never reaches the canonical episode table.
TZA_MAINLAND_UID <- "AFR::TZA::Mainland-1"

markTzaMainlandAdm0 <- function(dat) {
  dplyr::mutate(dat, admin_level = dplyr::case_when(uid == TZA_MAINLAND_UID ~ "ADM0",
                                                    TRUE ~ admin_level))
}

flattenGid <- function(str) {
  stringr::str_replace_all(str, "\\.", "_")
}

#' gidAdm1
#'
#' ADM1 key ("ISO3_<n>", e.g. "AGO_1") of a GADM gid, flat or dotted. The one
#' parser every ADM1-from-gid site goes through.
#'
#' GADM 4.1 writes Ghana's ids without the "." after the ISO3 ("GHA5_2",
#' "GHA14.13_2"; every other country "AGO.1_1", "AGO.1.2_1"). The pipeline keeps
#' those ids as is (climate, population and crosswalk artefacts all key on them),
#' but `str_extract(flattenGid(gid), "[A-Z]{3}_[0-9]+")` returns NA on them, and
#' that NA ADM1 silently dropped Ghana from every brms fit (NA grouping rows are
#' excluded) and counted it as one NA "ADM1" in report 11. So the separator is
#' restored here, the result is "GHA_5" / "GHA_14", and a gid that still does not
#' parse is an error, not an NA.
#'
#' @param gid character gids (NA passes through as NA).
#' @param tag label for the error message.
#' @param strict FALSE returns NA for an unparseable id instead of stopping:
#'   only for keying a whole GADM layer (map joins), which also holds ids that
#'   are not study units (e.g. the disputed-area "Z01.1_1").
gidAdm1 <- function(gid, tag = "gidAdm1", strict = TRUE) {
  gid  <- as.character(gid)
  adm1 <- stringr::str_extract(
    sub("^([A-Z]{3})([0-9])", "\\1_\\2", flattenGid(gid)), "^[A-Z]{3}_[0-9]+")
  bad <- unique(gid[!is.na(gid) & is.na(adm1)])
  if (strict && length(bad))
    stop(sprintf("[%s] %d gid(s) with no parseable ADM1 (expected ISO3[._]<n>...): %s",
                 tag, length(bad), paste(utils::head(bad, 10), collapse = ", ")))
  adm1
}


# GADM ids ----------------------------------------------------------------

# Function to get next gid
nextGid <- function(gid, offset = 0) {
  gid <- flattenGid(gid)
  gid_nums <- stringr::str_split(gid, "_")[[1]]
  ind <- length(gid_nums) - 1
  
  gid_nums[ind] <- as.numeric(gid_nums[ind]) - offset + 1
  
  stringr::str_c(gid_nums, collapse = "_")
}

getNextGid <- function(gid, 
                       ref_gids, 
                       n_steps = 10) {
  old_gid <- gid
  
  for (i in 1:n_steps) {
    offsets <- 0:(n_steps-1)
    for (o in offsets) {
      new_gid <- nextGid(old_gid)
      
      if (new_gid %in% ref_gids) {
        found <- T
        break()
      } else {
        old_gid <- new_gid
        found <- F
      }
    }
    
    if (found) {
      break()
    }
  }
  
  if (!found) {
    NA_character_
  } else {
    new_gid
  }
}


# Other missing could be due to missmatch between run levels, for this convert all to admin 1
getADM1Gid <- function(gid) {
  gid <- flattenGid(gid)
  gid_nums <- stringr::str_split(gid, "_")[[1]]

  stringr::str_c(c(gid_nums[1:2], 1), collapse = "_")
}

#' getNearestSeasGid
#'
#' Spatial nearest-neighbour fallback for the seasonality-gid match cascade in
#' 03_prepare_occ.R / 03_prepare_siz.R Section A, inserted AFTER the exact match and
#' the two getNextGid()/getADM1Gid() numeric fallbacks fail. For each target gid,
#' returns the candidate gid (from ref_gids) whose GADM centroid is nearest,
#' RESTRICTED to candidates in the same country (ISO3 prefix) as the target -- a
#' cross-country nearest match is never geographically meaningful, so a target whose
#' country has zero candidates in ref_gids returns NA (falls through to the January
#' default applied downstream in Section B).
#'
#' Batches by country: loadGadmUnits() reads a whole country's ADM1+ADM2 layers in one
#' sf::st_read() call, so this loops over unique countries, not rows.
#'
#' Callers must source analysis/utils/population_core.R themselves for loadGadmUnits()
#' -- gid.R does not source it (modules never source() each other; see
#' analysis/utils/README.md's dependency-table convention). New call-time (not
#' source-time) dependency: gid.R -> population_core.R::loadGadmUnits().
#'
#' @param target_gids character vector of flattened (underscore) gids to resolve.
#' @param ref_gids character vector of flattened candidate gids to search among; may
#'   span multiple countries, filtered internally to each target's own country.
#' @param gadm_path path to the validated GADM geopackage (data/geodata/gadm_valid.gpkg).
#' @return character vector, same length/order as target_gids: nearest ref_gids
#'   member, or NA_character_ when the target's country has no candidate in ref_gids
#'   or geometry can't be loaded.
getNearestSeasGid <- function(target_gids, ref_gids,
                              gadm_path = "data/geodata/gadm_valid.gpkg") {
  target_gids <- flattenGid(target_gids)
  ref_gids    <- flattenGid(ref_gids)
  out <- rep(NA_character_, length(target_gids))
  iso3_of <- function(g) stringr::str_extract(g, "^[A-Z]{3}")
  target_iso3 <- iso3_of(target_gids)
  ref_iso3    <- iso3_of(ref_gids)

  for (iso3 in unique(stats::na.omit(target_iso3))) {
    is_t <- which(target_iso3 == iso3)
    cand <- unique(ref_gids[ref_iso3 == iso3])
    if (!length(cand)) next  # country entirely outside the seasonality panel

    units <- tryCatch(
      loadGadmUnits(gadm_path, iso3, levels = c("ADM1", "ADM2")),
      error = function(e) NULL
    )
    if (is.null(units)) next
    units$gid <- flattenGid(units$gid)

    t_gid <- unique(target_gids[is_t])
    t_sf  <- units[units$gid %in% t_gid, , drop = FALSE]
    c_sf  <- units[units$gid %in% cand,  , drop = FALSE]
    if (!nrow(t_sf) || !nrow(c_sf)) next

    t_cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(t_sf)))
    c_cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(c_sf)))
    nn_idx <- sf::st_nearest_feature(t_cent, c_cent)
    nn_map <- setNames(c_sf$gid[nn_idx], t_sf$gid)
    out[is_t] <- unname(nn_map[target_gids[is_t]])
  }
  out
}

# Admin-level row filtering -----------------------------------------------

#' .overlapAdm1Set
#'
#' The ADM1 regions that carry ADM2 rows: the one definition behind the
#' ADM1-ADM2 overlap drop of the fit (filterOverlappingADM1()), the simulator
#' frame (dropOverlappingADM1FromFrame()) and the occurrence panel's orphan
#' classification (adm1OverlapParents()).
#' @param dat rows with admin_level and ADM1.
#' @return character vector of ADM1 codes.
.overlapAdm1Set <- function(dat) {
  dat %>%
    dplyr::filter(as.character(admin_level) == "ADM2", !is.na(ADM1)) %>%
    dplyr::pull(ADM1) %>% as.character() %>% unique()
}

#' filterOverlappingADM1
#'
#' For ADM1-ADM2 mixed data: keep ADM2 rows as-is, and only retain ADM1 rows
#' for which NO ADM2 sibling exists in the same ADM1 region.
filterOverlappingADM1 <- function(dat) {
  adm1_with_adm2 <- .overlapAdm1Set(dat)
  dplyr::filter(dat, !(admin_level == "ADM1" & ADM1 %in% adm1_with_adm2))
}

#' adm1OverlapParents
#'
#' The ADM1-level panel gids that the fit drops at load
#' (filterOverlappingADM1(): their ADM1 also carries ADM2 rows) and, per
#' canonical episode of such a parent, whether a child episode covers it.
#' A parent episode is COVERED when some ADM2 gid of the same ADM1 has an
#' episode whose month span overlaps the parent's ([month(TL), month(TR)]),
#' and an ORPHAN otherwise: its onset is lost to the fit and nothing in the
#' children reports it. Under episodes.adm1_overlap = censor_orphan the
#' children's rows inside an orphan span leave the panel (maskOccPanelRows()).
#'
#' @param surv panel rows (gid, admin_level, country).
#' @param ep canonical episodes (gid, year_seas, TL, TR).
#' @return tibble(gid, ADM1, country, year_seas, TL, TR, n_child_gids,
#'   n_child_episodes, class = "orphan" | "covered"); one row per parent episode.
adm1OverlapParents <- function(surv, ep) {
  g <- surv %>%
    dplyr::distinct(gid = as.character(gid), admin_level = as.character(admin_level),
                    country = as.character(country)) %>%
    dplyr::mutate(ADM1 = gidAdm1(gid, "adm1OverlapParents"))
  ov <- .overlapAdm1Set(g)
  parents  <- g %>% dplyr::filter(admin_level == "ADM1", ADM1 %in% ov)
  children <- g %>% dplyr::filter(admin_level == "ADM2", ADM1 %in% parents$ADM1)
  empty <- tibble::tibble(gid = character(), ADM1 = character(), country = character(),
                          year_seas = integer(), TL = as.Date(character()), TR = as.Date(character()),
                          n_child_gids = integer(), n_child_episodes = integer(), class = character())
  if (!nrow(parents)) return(empty)
  m <- function(d) lubridate::floor_date(as.Date(d), "month")
  pe <- ep %>% dplyr::filter(as.character(gid) %in% parents$gid) %>%
    dplyr::transmute(gid = as.character(gid), year_seas = as.integer(year_seas),
                     TL = as.Date(TL), TR = as.Date(TR)) %>%
    dplyr::left_join(dplyr::select(parents, gid, ADM1, country), by = "gid")
  if (!nrow(pe)) return(empty)
  ce <- ep %>% dplyr::filter(as.character(gid) %in% children$gid) %>%
    dplyr::transmute(c_gid = as.character(gid), c_TL = as.Date(TL), c_TR = as.Date(TR)) %>%
    dplyr::left_join(dplyr::select(children, c_gid = gid, ADM1), by = "c_gid")
  n_kids <- children %>% dplyr::count(ADM1, name = "n_child_gids")
  hit <- pe %>%
    dplyr::mutate(.k = dplyr::row_number()) %>%
    dplyr::inner_join(ce, by = "ADM1", relationship = "many-to-many") %>%
    dplyr::filter(m(c_TL) <= m(TR), m(c_TR) >= m(TL)) %>%
    dplyr::count(.k, name = "n_child_episodes")
  pe %>%
    dplyr::mutate(.k = dplyr::row_number()) %>%
    dplyr::left_join(hit, by = ".k") %>%
    dplyr::left_join(n_kids, by = "ADM1") %>%
    dplyr::mutate(n_child_gids = as.integer(dplyr::coalesce(n_child_gids, 0L)),
                  n_child_episodes = as.integer(dplyr::coalesce(n_child_episodes, 0L)),
                  class = ifelse(n_child_episodes > 0L, "covered", "orphan")) %>%
    dplyr::select(gid, ADM1, country, year_seas, TL, TR, n_child_gids, n_child_episodes, class) %>%
    dplyr::arrange(gid, TL)
}

#' dropOverlappingADM1FromFrame
#'
#' Simulator-frame analogue of `filterOverlappingADM1()`. `build_location_frame()`
#' ingests every distinct panel gid, including coarse ADM1-level units that also
#' have nested ADM2 children (mixed-admin fits). Keeping both double-counts
#' geography AND injects whole-state ADM1 populations that inflate the size mean
#' (mu = pop x mean_rate x exp(lp)). Using the survival data's authoritative
#' `admin_level`/`ADM1` columns, drop the coarse ADM1 gid wherever an ADM2 sibling
#' exists, leaving ADM1-only units untouched.
#'
#' This is the single source of truth shared by the smoothing pipeline
#' (`smoothing_core.R`) and the forecast loader (`utils/simulations.R`) so the
#' two pipelines cannot drift on this drop. Returns `frame` unchanged (with a
#' warning) when `surv_dat` is NULL or lacks the required columns.
#'
#' @param frame tibble from build_location_frame() (must have a `gid` column).
#' @param surv_dat survival data with `admin_level`, `ADM1`, `gid` columns.
#' @param tag character prefix for log messages (e.g. "smoothing").
#' @return frame, possibly with coarse-ADM1 overlap rows removed.
dropOverlappingADM1FromFrame <- function(frame, surv_dat, tag = "sim") {
  if (is.null(surv_dat) ||
      !all(c("admin_level", "ADM1", "gid") %in% names(surv_dat))) {
    warning("[", tag, "] survival data NULL or lacks admin_level/ADM1/gid — ",
            "ADM1-ADM2 overlap filter skipped.")
    return(frame)
  }
  adm1_with_adm2 <- .overlapAdm1Set(surv_dat)
  overlap_gids <- surv_dat %>%
    dplyr::filter(as.character(admin_level) == "ADM1",
                  as.character(ADM1) %in% adm1_with_adm2) %>%
    dplyr::pull(gid) %>% as.character() %>% unique()
  n_drop <- sum(frame$gid %in% overlap_gids)
  if (n_drop > 0L) {
    frame <- frame[!(frame$gid %in% overlap_gids), , drop = FALSE]
    message(sprintf(paste0(
      "[%s] ADM1-ADM2 overlap filter: dropped %d coarse ADM1 unit(s) ",
      "with nested ADM2 children (mirrors filterOverlappingADM1)."), tag, n_drop))
  }
  frame
}

#' applyAdminLevelFilter
#'
#' Filter rows to the requested admin level.
#' Valid values: "ADM1", "ADM2", "ADM1-ADM2".
#' Defaults to "ADM2" when admin_level_opt is NULL.
applyAdminLevelFilter <- function(dat, admin_level_opt) {
  switch(admin_level_opt %||% "ADM2",
    ADM1        = dplyr::filter(dat, admin_level == "ADM1"),
    ADM2        = dplyr::filter(dat, admin_level == "ADM2"),
    `ADM1-ADM2` = filterOverlappingADM1(dat),
    {
      warning("Unknown admin_level: ", admin_level_opt, " — no row filter applied.")
      dat
    }
  )
}
