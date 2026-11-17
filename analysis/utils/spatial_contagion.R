# spatial_contagion.R — distance-weighted neighbour-contagion covariate of
#   the occurrence track
#
# Provides:
#   build_adm1_dist_km()        cached ADM1-centroid distance matrix (km)
#   kernel_from_dist(), build_nbr_kernel(), build_partition_nbr_kernel()
#   neighbour_pressure(), nbr_from_adm1_active()   the covariate itself
#   build_contagion_table()     the prep-time covariate table
#   resolve_nbr_spec(), resolve_use_nbr(), nbr_cfg_from_sim(),
#   nbr_kernel_for_frame()      config -> kernel, at fit and simulate time
#   contagion_coverage_shares() surveillance coverage of the neighbourhoods
#
# Details:
#   nbr_rho(i, t) = sum_{j != i} active_j(t-1) * exp(-d_ij / rho), with
#   active_j(t-1) the "ADM1 unit j had an ongoing outbreak in month t-1"
#   indicator of 03_prepare_occ.R and d_ij the great-circle distance between
#   ADM1 centroids. Callers source gid.R (flatten_gid()) and paths.R
#   (make_adm1_dist_file()); packages sf, yaml.

#' build_adm1_dist_km
#'
#' Great-circle distance matrix (km) between GADM ADM1 centroids, cached.
#'
#' Reads the ADM_1 layer ONCE and filters to `iso3_set`. Deliberately does not
#' reuse load_gadm_units() (analysis/utils/population_core.R), which st_read()s the
#' whole continental ADM_1 layer and filters by a single GID_0 afterwards --
#' calling it per country would re-read the full layer once per country. The
#' returned object nonetheless uses the same flattened-`gid` key that
#' load_gadm_units() emits, so the two stay interchangeable downstream.
#'
#' The diagonal is set to `Inf`, not 0, so that kernel_from_dist() yields exactly
#' zero self-weight and a unit's own `ongoing_outbreak_lag1` can never leak into
#' its own neighbour pressure. This is the same leave-one-out discipline as
#' `(country_count_lag1 - as.integer(ongoing_outbreak_lag1)) > 0` in
#' 03_prepare_occ.R.
#'
#' The cached universe is ALL GADM ADM1 units in `iso3_set`, not just the ones
#' with observed outbreaks: a neighbour that was never active contributes 0 to
#' the kernel anyway, and restricting to active units would make the cache key
#' depend on the outbreak data.
#'
#' @param iso3_set  character vector of ISO3 country codes to include.
#' @param gadm_path path to the validated GADM geopackage.
#' @param cache_dir directory for the cached matrix + metadata sidecar.
#' @param refresh   if TRUE, rebuild even when a cache file exists.
#' @return numeric matrix, symmetric apart from an `Inf` diagonal, with
#'   dimnames = flattened ADM1 gids (e.g. "NGA_1"), distances in km.
build_adm1_dist_km <- function(
  iso3_set,
  gadm_path = "data/geodata/gadm_valid.gpkg",
  cache_dir = "generated_data/",
  refresh = FALSE
) {
  iso3_set <- sort(unique(stats::na.omit(as.character(iso3_set))))
  cache_file <- make_adm1_dist_file(cache_dir = cache_dir)
  meta_file <- sub("\\.rds$", "_meta.yml", cache_file)

  if (!refresh && file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    # Cache is only reusable if it covers every requested country. A stale cache
    # built on a smaller iso3 set would silently give those units zero neighbour
    # pressure, which looks like a real (null) result rather than a missing input.
    #
    # "Covered" = present in the matrix OR requested when the cache was built
    # (the sidecar's `iso3`). A requested country can legitimately yield no rows.
    # Counting only matrix rows would make such a country uncoverable, so every
    # call that requested it would rebuild the matrix from the geopackage --
    # which needs `sf` (absent on some HPC nodes) and races concurrent array
    # tasks writing the same file.
    #
    # The sidecar is trusted only when it records `gid_adm1: gid_to_adm1`. Ghana's
    # GADM GID_1 ("GHA1_2", no "." after the ISO3) cannot be keyed by a bare
    # regex, only by gid_to_adm1(), so a cache built without it lists GHA yet has no
    # GHA rows, and trusting it would
    # give every Ghana ADM1 zero neighbour pressure. Such a cache falls back to
    # matrix-row coverage and is rebuilt.
    cached_iso3 <- sort(unique(stringr::str_extract(
      rownames(cached),
      "^[A-Z]{3}"
    )))
    meta <- if (file.exists(meta_file)) yaml::read_yaml(meta_file) else list()
    built_iso3 <- if (identical(meta$gid_adm1, "gid_to_adm1")) {
      as.character(unlist(meta$iso3))
    } else {
      character(0)
    }
    if (all(iso3_set %in% union(cached_iso3, built_iso3))) {
      message(
        "build_adm1_dist_km: loading cached ADM1 distance matrix from ",
        cache_file
      )
      return(cached)
    }
    message(
      "build_adm1_dist_km: cache at ",
      cache_file,
      " is missing ",
      sum(!iso3_set %in% cached_iso3),
      " requested country/countries -- rebuilding."
    )
  }

  if (!file.exists(gadm_path)) {
    stop("build_adm1_dist_km: GADM geopackage not found at ", gadm_path)
  }

  message("build_adm1_dist_km: reading ADM_1 layer from ", gadm_path)
  adm1 <- sf::st_read(gadm_path, layer = "ADM_1", quiet = TRUE)
  adm1 <- adm1[adm1$GID_0 %in% iso3_set, , drop = FALSE]
  if (!nrow(adm1)) {
    stop(
      "build_adm1_dist_km: no ADM_1 features matched iso3_set (",
      paste(utils::head(iso3_set, 5), collapse = ", "),
      ", ...)"
    )
  }

  # Flatten AND trim -- BOTH are required (gid_to_adm1() does both). The survival /
  # outbreak side keys ADM1 as `gid_to_adm1(gid)`, i.e. "AGO_1". Flattening alone
  # turns GADM's "AGO.1_1" into "AGO_1_1", which shares not one key with the survival side -- the join then matches nothing,
  # every column coalesces to 0, and the covariate silently reads as a constant
  # zero rather than erroring. (Observed: 0/255 ADM1s matched before the trim.)
  # This is the third flatten_gid() site; the other two are in 03_prepare_occ.R,
  # each carrying its own warning about the same class of blowup.
  # gid_to_adm1() also restores the ISO3 separator GADM 4.1 omits for Ghana
  # ("GHA1_2" -> "GHA_1"), which a bare regex would drop (no Ghana rows).
  gids <- gid_to_adm1(adm1$GID_1, "build_adm1_dist_km")
  keep <- !is.na(gids) & !duplicated(gids)
  adm1 <- adm1[keep, , drop = FALSE]
  gids <- gids[keep]
  ord <- order(gids)
  adm1 <- adm1[ord, , drop = FALSE]
  gids <- gids[ord]

  cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(adm1)))
  dm <- sf::st_distance(cent) # metres (geodesic on 4326)
  D <- matrix(
    as.numeric(dm),
    nrow = length(gids),
    dimnames = list(gids, gids)
  ) /
    1000 # km
  diag(D) <- Inf

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(D, cache_file)
  yaml::write_yaml(
    list(
      built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      gadm_path = gadm_path,
      gadm_mtime = format(file.mtime(gadm_path), "%Y-%m-%d %H:%M:%S"),
      n_units = length(gids),
      gid_adm1 = "gid_to_adm1",
      iso3 = as.list(iso3_set)
    ),
    meta_file
  )
  message(
    "build_adm1_dist_km: wrote ",
    length(gids),
    "x",
    length(gids),
    " distance matrix to ",
    cache_file
  )
  D
}


#' kernel_from_dist
#'
#' Exponential distance-decay kernel. The `Inf` diagonal set by
#' build_adm1_dist_km() yields exactly 0 self-weight.
#'
#' @param D       distance matrix in km, from build_adm1_dist_km().
#' @param rho_km  decay bandwidth in km.
#' @return matrix of weights, same dimnames as `D`, zero on the diagonal.
kernel_from_dist <- function(D, rho_km) {
  if (
    !is.numeric(rho_km) ||
      length(rho_km) != 1L ||
      !is.finite(rho_km) ||
      rho_km <= 0
  ) {
    stop(
      "kernel_from_dist: rho_km must be a single positive finite number, got: ",
      paste(utils::capture.output(str(rho_km)), collapse = " ")
    )
  }
  exp(-D / rho_km)
}


# Motivation: `ongoing_outbreak_country_lag1` / `ongoing_outbreak_tu_lag1` are
# outcome-leakage artifacts and are not used, which leaves the occurrence model
# with no cross-unit spatial coupling. A GLM probe found the rho = 500 km kernel
# positive in 6/6 West-group countries with >= 20 onsets, robust to dropping
# `own` and to adding an ADM1 random intercept; rho = 200 km is inconsistent.
#
# Resolution: the kernel is built at ADM1 grain. Survival rows can be
# ADM0/ADM1/ADM2 (occurrence.admin_level = "ADM1-ADM2" in production), and ADM1
# is the only grain at which kernel weights are commensurable across rows AND at
# which self-exclusion cleanly separates the neighbour term from the unit's own
# `ongoing_outbreak_lag1`. ADM2 rows inherit their parent ADM1's value, so the
# covariate is constant within an ADM1 across its ADM2 children.
#' neighbour_pressure
#'
#' Apply a distance kernel to an ADM1 x month activity grid.
#'
#' `active_long` must already hold the LAGGED indicator keyed at the current
#' month -- i.e. pass `ongoing_outbreak_lag1` (which 03_prepare_occ.R defines as
#' "active at t_start - 30") rather than `ongoing_outbreak`. No time shifting
#' happens here, exactly as `country_lag1` sums `ongoing_outbreak_lag1` within
#' `(country, t_start)` without shifting.
#'
#' Units present in `active_long` but absent from `K` (e.g. a crosswalk-derived
#' ADM1 with no GADM ADM_1 feature) contribute 0 to every other unit's pressure
#' and receive `NA`, which the caller coalesces to 0.
#'
#' @param active_long data frame with columns ADM1 (chr), t_start, and the
#'   activity column named by `value_col` (logical or 0/1).
#' @param K           kernel matrix from kernel_from_dist().
#' @param normalise   if TRUE, divide each row by rowSums(K), turning the count
#'   of active neighbours into a prevalence among neighbours. Changes the
#'   estimand (and inflates coastal/border units, whose denominator covers less
#'   land) -- kept as a robustness sibling, not the modelled default.
#' @param value_col   name of the activity column in `active_long`.
#' @return tibble with columns ADM1, t_start, nbr (numeric).
neighbour_pressure <- function(
  active_long,
  K,
  normalise = FALSE,
  value_col = "ongoing_outbreak_lag1"
) {
  units <- rownames(K)
  months <- sort(unique(active_long$t_start))

  A <- matrix(
    0,
    nrow = length(units),
    ncol = length(months),
    dimnames = list(units, as.character(months))
  )
  ri <- match(as.character(active_long$ADM1), units)
  ci <- match(as.character(active_long$t_start), colnames(A))
  ok <- !is.na(ri) & !is.na(ci)
  A[cbind(ri[ok], ci[ok])] <- as.numeric(active_long[[value_col]][ok] %in% TRUE)

  E <- K %*% A
  if (isTRUE(normalise)) {
    rs <- rowSums(K)
    rs[rs == 0] <- 1 # isolated unit: pressure is 0/1 = 0, not 0/0 = NaN
    E <- E / rs
  }

  tibble::tibble(
    ADM1 = rep(units, times = ncol(E)),
    t_start = rep(months, each = nrow(E)),
    nbr = as.numeric(E)
  )
}


#' build_nbr_kernel
#'
#' THE definition of the neighbour-contagion kernel, shared by the fit side
#' (build_contagion_table() -> 03_prepare_occ.R) and every simulator kernel
#' (build_partition_nbr_kernel() -> v1 / spatPomp R / spatPomp C), so the
#' covariate the occurrence model was fitted on and the one the simulator
#' recomputes agree by construction rather than by parallel code.
#'
#'   K_ij = exp(-d_ij / rho)                       (kernel_from_dist())
#'   partition: K_ij <- 0 unless group(i) == group(j); units whose country is in
#'              no group get an all-zero row and column (fit side: NA -> 0)
#'   rownorm:   K_ij <- K_ij / sum_j K_ij  (over the SAME masked universe), so
#'              K %*% a is the distance-weighted share of active neighbours,
#'              bounded in [0, 1]; isolated rows stay 0 (0/1, not 0/0)
#'
#' The normalising denominator is taken over every ADM1 in `D` (after the
#' partition mask), NOT over the simulator frame's ADM1 subset: ADM1s dropped
#' from a frame (e.g. `unseen_adm1: restrict`) were still in the fit-time
#' denominator. Numerators agree because such ADM1s were never active.
#'
#' @param D         distance matrix (km) from build_adm1_dist_km(), Inf diagonal.
#' @param rho_km    decay bandwidth in km.
#' @param scale     "raw" | "rownorm".
#' @param partition NULL, or data.frame(iso3, group) (load_grouping()).
#' @return matrix, same dimnames as `D`, zero diagonal.
build_nbr_kernel <- function(D, rho_km, scale = "raw", partition = NULL) {
  K <- kernel_from_dist(D, rho_km)
  if (!is.null(partition)) {
    grp <- stats::setNames(
      as.character(partition$group),
      as.character(partition$iso3)
    )
    g <- unname(grp[stringr::str_extract(rownames(K), "^[A-Z]{3}")])
    same <- outer(g, g, function(a, b) !is.na(a) & !is.na(b) & a == b)
    K[!same] <- 0
  }
  if (identical(scale, "rownorm")) {
    rs <- rowSums(K)
    rs[rs == 0] <- 1
    K <- K / rs
  }
  K
}


#' build_contagion_table
#'
#' Convenience wrapper producing all four contagion columns for one ADM1
#' activity grid, in the (ADM1, t_start)-keyed shape 03_prepare_occ.R joins on.
#'
#' All rho x scale combinations are always written, following the repo's own
#' "all 13 ENSO lag columns are always written so lag-sensitivity configs can
#' select any lag without re-running prep" convention -- a later rho- or
#' scale-sensitivity run then needs no prep re-run.
#'
#' With `partition` supplied, the partition-restricted siblings
#' (`nbr_lag1_r<rho>[_rn]_pg`) are written IN ADDITION, through the shared
#' build_nbr_kernel() so they match the simulator kernel exactly. The four
#' unpartitioned columns are computed exactly as before (byte-identical).
#'
#' @param active_long ADM1 x month activity grid (see neighbour_pressure()).
#' @param D           distance matrix from build_adm1_dist_km().
#' @param rhos        integer vector of bandwidths in km.
#' @param partition   NULL (default) or data.frame(iso3, group).
#' @param partitions  NULL or named list (registry name -> data.frame(iso3,
#'   group)) for the `_part<token>` siblings.
#' @return tibble: ADM1, t_start, nbr_lag1_r<rho>, nbr_lag1_r<rho>_rn
#'   [, nbr_lag1_r<rho>_pg, nbr_lag1_r<rho>_rn_pg]
#'   [, nbr_lag1_r<rho>_part<token>, nbr_lag1_r<rho>_rn_part<token> ...].
build_contagion_table <- function(
  active_long,
  D,
  rhos = c(200L, 500L),
  partition = NULL,
  partitions = NULL
) {
  out <- NULL
  add <- function(tb, col) {
    names(tb)[names(tb) == "nbr"] <- col
    if (is.null(out)) {
      tb
    } else {
      dplyr::left_join(out, tb, by = c("ADM1", "t_start"))
    }
  }
  for (rho in rhos) {
    K <- kernel_from_dist(D, rho)
    for (nrm in c(FALSE, TRUE)) {
      out <- add(
        neighbour_pressure(active_long, K, normalise = nrm),
        nbr_covariate_name(rho, if (nrm) "rownorm" else "raw")
      )
    }
    if (!is.null(partition)) {
      for (sc in c("raw", "rownorm")) {
        Kp <- build_nbr_kernel(D, rho, scale = sc, partition = partition)
        out <- add(
          neighbour_pressure(active_long, Kp, normalise = FALSE),
          nbr_covariate_name(rho, sc, "group")
        )
      }
    }
    # Registry-partitioned siblings nbr_lag1_r<rho>[_rn]_part<token>, one pair
    # per registry in `partitions` (named list: registry name -> data.frame(iso3,
    # group)), through the same build_nbr_kernel() as the simulator.
    for (src in names(partitions %||% list())) {
      for (sc in c("raw", "rownorm")) {
        Kp <- build_nbr_kernel(
          D,
          rho,
          scale = sc,
          partition = partitions[[src]]
        )
        out <- add(
          neighbour_pressure(active_long, Kp, normalise = FALSE),
          nbr_covariate_name(rho, sc, "registry", src)
        )
      }
    }
  }
  out
}


#' build_partition_nbr_kernel
#'
#' Dense ADM1 x ADM1 distance-decay kernel for a single simulator partition,
#' aligned to the partition's own ADM1 ordering.
#'
#' The forward simulator (07b_simulator_v1.R) and the spatPomp filter
#' (07b_simulator_spatpomp.R) both recompute neighbour pressure dynamically from
#' the simulated / filtered ADM1-active state, rather than folding an exogenous
#' observed value into `lp_static`. Both need the SAME kernel matrix `K`; they
#' differ only in WHEN they supply the active vector (R: per simulated step;
#' C: embedded-matrix matvec in the Csnippet). This helper is that shared build.
#'
#' The returned matrix is indexed in `unique(frame_tu$adm1_unit)` order -- the
#' exact order `make_membership_arrays_c()` (spatpomp_core_c.R) assigns to
#' `adm1_id`, so the C kernel can flatten it row-major into `WNBR[]` and index it
#' by the same 0-based ADM1 id. Under `partition_by: country`, the partition's
#' spatPomp object contains exactly that country's ADM1 units, so restricting `K`
#' to those gids IS the within-country kernel -- cross-country neighbours are
#' absent from the object (the deferred cross-country work).
#'
#' Partition ADM1 units absent from the GADM-derived distance matrix (e.g. a
#' crosswalk ADM1 with no ADM_1 feature) get zero rows and columns: they exert
#' and receive zero neighbour pressure, matching neighbour_pressure()'s
#' coalesce-missing-to-0 behaviour on the fit side.
#'
#' @param frame_tu  partition frame; must carry `adm1_unit` (flattened gid, e.g.
#'   "NGA_1").
#' @param rho_km    decay bandwidth in km.
#' @param gadm_path path to the validated GADM geopackage.
#' @param cache_dir directory for the cached distance matrix.
#' @param refresh   if TRUE, rebuild the distance cache.
#' @param scale     "raw" (default) | "rownorm" -- see build_nbr_kernel().
#' @param partition NULL (default) | data.frame(iso3, group) -- see
#'   build_nbr_kernel(). Both default to the pre-existing raw, unpartitioned
#'   kernel; resolve_nbr_spec() derives them from the fit's FE name.
#' @return dense numeric matrix, dim n_adm1 x n_adm1 in partition ADM1 order,
#'   zero diagonal (self-excluded), dimnames = partition ADM1 gids.
build_partition_nbr_kernel <- function(
  frame_tu,
  rho_km,
  gadm_path = "data/geodata/gadm_valid.gpkg",
  cache_dir = "generated_data/",
  refresh = FALSE,
  scale = "raw",
  partition = NULL
) {
  adm1_order <- unique(as.character(frame_tu$adm1_unit))
  iso3 <- unique(stringr::str_extract(adm1_order, "^[A-Z]{3}"))
  # A partition-restricted kernel normalises over the WHOLE group, so the
  # distance cache must also cover member countries absent from this frame.
  if (!is.null(partition)) {
    grp_here <- unique(partition$group[partition$iso3 %in% iso3])
    iso3 <- union(iso3, partition$iso3[partition$group %in% grp_here])
  }

  D <- build_adm1_dist_km(
    iso3,
    gadm_path = gadm_path,
    cache_dir = cache_dir,
    refresh = refresh
  )
  # Raw + unpartitioned is kernel_from_dist() itself (the pre-existing path).
  K <- build_nbr_kernel(D, rho_km, scale = scale, partition = partition)

  n <- length(adm1_order)
  K_adm1 <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(adm1_order, adm1_order)
  )
  present <- adm1_order[adm1_order %in% rownames(K)]
  if (length(present)) {
    K_adm1[present, present] <- K[present, present]
  }
  diag(K_adm1) <- 0 # belt-and-braces; kernel_from_dist's Inf diagonal already 0
  K_adm1
}


#' nbr_from_adm1_active
#'
#' Apply a partition ADM1 kernel to a lag-1 ADM1-active 0/1 vector. Kept in one
#' place so the v1 forecast loop and the spatPomp R-function rprocess share the
#' exact matvec (`K %*% Aa1`), and so the C kernel's `WNBR` matvec has a single R
#' reference to match in the equivalence test.
#'
#' @param Aa1    numeric vector, length n_adm1, ADM1-active indicator (0/1) in
#'   the same order as `rownames(K_adm1)`.
#' @param K_adm1 kernel matrix from build_partition_nbr_kernel().
#' @return numeric vector, length n_adm1, neighbour pressure per ADM1.
nbr_from_adm1_active <- function(Aa1, K_adm1) {
  as.numeric(K_adm1 %*% Aa1)
}


#' resolve_nbr_spec
#'
#' Resolve EVERYTHING the simulator needs to rebuild the fit's neighbour
#' covariate, in one place, shared by the v1 driver, the spatPomp forecast
#' driver, the smoothing driver and (as a fallback) the three kernel builders:
#'
#'   use       -- config flag on AND an occurrence fit carries an nbr FE
#'   col       -- that FE's name (nbr_covariate_pattern(), paths.R)
#'   rho_km, scale, partition -- parsed from `col` (parse_nbr_covariate()), so
#'                the kernel can never silently disagree with the fitted
#'                covariate; a config rho_km that disagrees is an error
#'   registry  -- data.frame(iso3, group) when partition == "group", from
#'                load_grouping(cfg$grouping_source) (config.R, call-time)
#'   gadm_path, cache_dir -- passed through from `cfg`
#'
#' Under grouping, every group's occurrence fit that carries the term must
#' carry the SAME name (one kernel per partition); groups whose fit lacks it
#' still get beta_nbr = 0 downstream, exactly as before.
#'
#' @param cfg       the `simulator.nbr_contagion` list (enable, rho_km,
#'   gadm_path, cache_dir) plus optional `grouping_source`.
#' @param draws_occ ungrouped occurrence draw bundle (`$fe` S x P matrix).
#' @param grouping  grouped routing list with `$draws_occ_by_group`, or NULL.
#' @return list(use, col, rho_km, scale, partition, registry, gadm_path,
#'   cache_dir); `use = FALSE` (other fields default) when off.
resolve_nbr_spec <- function(cfg, draws_occ = NULL, grouping = NULL) {
  spec <- list(
    use = FALSE,
    col = NULL,
    rho_km = as.integer(cfg$rho_km %||% 500L),
    scale = "raw",
    partition = "none",
    registry = NULL,
    gadm_path = cfg$gadm_path %||% "data/geodata/gadm_valid.gpkg",
    cache_dir = cfg$cache_dir %||% "generated_data/"
  )
  if (!isTRUE(cfg$enable)) {
    return(spec)
  }
  fe_nbr <- function(d) {
    if (is.null(d) || is.null(d$fe)) {
      return(character(0))
    }
    grep(nbr_covariate_pattern(), colnames(d$fe), value = TRUE)
  }
  bundles <- if (!is.null(grouping) && !is.null(grouping$draws_occ_by_group)) {
    grouping$draws_occ_by_group
  } else {
    list(draws_occ)
  }
  cols <- unique(unlist(lapply(bundles, fe_nbr)))
  if (!length(cols)) {
    return(spec)
  }
  if (length(cols) > 1L) {
    stop(
      "resolve_nbr_spec: occurrence fits carry different neighbour-contagion ",
      "terms (",
      paste(cols, collapse = ", "),
      "); one kernel per run needs one ",
      "definition.",
      call. = FALSE
    )
  }
  p <- parse_nbr_covariate(cols)
  if (!is.null(cfg$rho_km) && as.integer(cfg$rho_km) != p$rho_km) {
    stop(
      sprintf(
        paste0(
          "resolve_nbr_spec: simulator.nbr_contagion.rho_km = %d but ",
          "the fit's term is %s (rho %d)."
        ),
        as.integer(cfg$rho_km),
        cols,
        p$rho_km
      ),
      call. = FALSE
    )
  }
  spec$use <- TRUE
  spec$col <- cols
  spec$rho_km <- p$rho_km
  spec$scale <- p$scale
  spec$partition <- p$partition
  if (identical(p$partition, "group")) {
    src <- cfg$grouping_source
    if (is.null(src) || !nzchar(src)) {
      stop(
        "resolve_nbr_spec: the fit's term ",
        cols,
        " is partition-restricted ",
        "but no grouping source was supplied (simulator.grouping.source).",
        call. = FALSE
      )
    }
    spec$registry <- load_grouping(src)
  } else if (identical(p$partition, "registry")) {
    # The registry is named by the fit's own term (`_part<token>`); a configured
    # simulator.nbr_contagion.partition_source must agree with it.
    src <- cfg$partition_source %||% ""
    if (nzchar(src)) {
      if (
        !identical(tolower(gsub("[^A-Za-z0-9]", "", src)), p$partition_token)
      ) {
        stop(
          sprintf(
            paste0(
              "resolve_nbr_spec: simulator.nbr_contagion.partition_source = ",
              "'%s' but the fit's term %s names registry token '%s'."
            ),
            src,
            cols,
            p$partition_token
          ),
          call. = FALSE
        )
      }
    } else {
      src <- grouping_source_from_token(p$partition_token)
    }
    spec$registry <- load_grouping(src)
    spec$partition_source <- src
  }
  spec
}


#' resolve_use_nbr
#'
#' Back-compat on/off guard: `resolve_nbr_spec(...)$use`. The term is active
#' iff BOTH the config flag is on (`simulator.nbr_contagion.enable: true`) AND
#' an occurrence fit carries an nbr fixed effect (any group, under grouping).
resolve_use_nbr <- function(cfg, draws_occ = NULL, grouping = NULL) {
  if (!isTRUE(cfg$enable)) {
    return(FALSE)
  }
  fe_has_nbr <- function(d) {
    if (is.null(d) || is.null(d$fe)) {
      return(FALSE)
    }
    any(grepl(nbr_covariate_pattern(), colnames(d$fe)))
  }
  if (!is.null(grouping) && !is.null(grouping$draws_occ_by_group)) {
    return(any(vapply(grouping$draws_occ_by_group, fe_has_nbr, logical(1))))
  }
  fe_has_nbr(draws_occ)
}


#' nbr_kernel_for_frame
#'
#' The one call every kernel builder makes: the partition kernel for
#' `frame_tu` under a resolved spec (resolve_nbr_spec()).
nbr_kernel_for_frame <- function(frame_tu, spec) {
  build_partition_nbr_kernel(
    frame_tu,
    spec$rho_km,
    gadm_path = spec$gadm_path,
    cache_dir = spec$cache_dir,
    scale = spec$scale,
    partition = spec$registry
  )
}


#' nbr_cfg_from_sim
#'
#' The `cfg` list resolve_nbr_spec() takes, from a merged `simulator:` block:
#' `nbr_contagion` plus the grouping source a partition-restricted term needs.
nbr_cfg_from_sim <- function(sim) {
  cfg <- sim$nbr_contagion %||% list()
  cfg$grouping_source <- sim$grouping$source
  # registry-partitioned term: the configured contagion registry, else the
  # partition registry (the contagion registry defaults to it).
  if (!nzchar(cfg$partition_source %||% "")) {
    cfg$partition_source <- sim$filter$partition_source %||% ""
  }
  cfg
}

# ---------------------------------------------------------------------------
# Contagion coverage report (STAMP ONLY)
# ---------------------------------------------------------------------------
# The neighbour-contagion covariate is not changed. Behind 03_prepare_occ.R
# --contagion_coverage_report these helpers measure how much of each unit's
# neighbourhood was under surveillance in the month the covariate reads (t-1),
# so the fit's implicit coverage can be reported (D33) next to the simulator's
# (every unit of the universe is simulated, coverage 1). Coverage definition
# "reporting": coverage.R::build_coverage_mask() over the
# panel's full LP member bridge with ADM1 windows, national zeros of the
# panel's own months file, UNION the recorded episode months -- never panel-row
# existence (masked occurrence-panel months would read as uncovered).

#' Gid x month surveillance coverage ("reporting"): build_coverage_mask() of
#' the member bridge (include_adm1 = TRUE) union the recorded episode months.
#' @param bridge LP member bridge (location_period_id, gid, country).
#' @param agg_ts,agg_adm_ts the aggregated series paths (make_agg_ts_file() /
#'   make_agg_adm_ts_file()).
#' @param episodes canonical episodes (gid, TL, TR).
#' @param nz_on,nz_path national zero-fill of the panel and its months file
#'   (passed EXPLICITLY; build_coverage_mask() has no default).
#' @return tibble(gid, month) of covered cells.
contagion_coverage_gid <- function(
  bridge,
  agg_ts,
  agg_adm_ts,
  episodes,
  nz_on,
  nz_path
) {
  cov <- build_coverage_mask(
    bridge,
    agg_ts,
    if (!is.null(agg_adm_ts) && file.exists(agg_adm_ts)) agg_adm_ts else NULL,
    include_adm1 = TRUE,
    national_zero = isTRUE(nz_on),
    national_zero_path = if (isTRUE(nz_on)) nz_path else NULL
  )
  ep <- episodes %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      gid = as.character(gid),
      m0 = lubridate::floor_date(as.Date(TL), "month"),
      m1 = lubridate::floor_date(as.Date(TR), "month")
    ) %>%
    dplyr::filter(!is.na(m0), !is.na(m1), m0 <= m1)
  epm <- if (nrow(ep)) {
    dplyr::bind_rows(lapply(seq_len(nrow(ep)), function(i) {
      tibble::tibble(
        gid = ep$gid[i],
        month = seq(ep$m0[i], ep$m1[i], by = "month")
      )
    }))
  } else {
    tibble::tibble(gid = character(), month = as.Date(character()))
  }
  dplyr::bind_rows(
    dplyr::mutate(cov, gid = as.character(gid), month = as.Date(month)),
    epm
  ) %>%
    dplyr::distinct(gid, month)
}

#' ADM1 x month coverage at lag `lag`: an ADM1 is covered at t_start when any
#' of its gids is covered in month(t_start) - lag.
#' @param cov_gid contagion_coverage_gid() output.
#' @param adm1_months tibble(ADM1, t_start) cells to evaluate.
#' @return `adm1_months` + covered (logical).
adm1_coverage_grid <- function(cov_gid, adm1_months, lag = 1L) {
  cv <- cov_gid %>%
    dplyr::mutate(ADM1 = gid_to_adm1(gid, "adm1_coverage_grid")) %>%
    dplyr::distinct(ADM1, month)
  adm1_months %>%
    dplyr::mutate(
      .m = lubridate::floor_date(as.Date(t_start), "month") %m-%
        months(as.integer(lag))
    ) %>%
    dplyr::left_join(
      dplyr::mutate(cv, covered = TRUE),
      by = c("ADM1", .m = "month")
    ) %>%
    dplyr::mutate(covered = dplyr::coalesce(covered, FALSE)) %>%
    dplyr::select(-.m)
}

#' Lag-1 activity of every ADM1 on a GLOBAL month grid (the panel
#' covariate's grid crosses each ADM1 only with its own country's panel months,
#' so a neighbour active at t-1 reads 0 when its country has no panel row at t).
#' @param act_spans tibble(ADM1, TL, TR) of recorded episodes.
#' @param adm1_set ADM1 codes; @param months t_start months.
#' @return tibble(ADM1, t_start, ongoing_outbreak_lag1).
adm1_activity_grid <- function(act_spans, adm1_set, months) {
  tidyr::expand_grid(
    ADM1 = as.character(adm1_set),
    t_start = sort(unique(as.Date(months)))
  ) %>%
    dplyr::left_join(
      dplyr::distinct(act_spans, ADM1, TL, TR),
      by = "ADM1",
      relationship = "many-to-many"
    ) %>%
    dplyr::group_by(ADM1, t_start) %>%
    dplyr::summarise(
      ongoing_outbreak_lag1 = active_last_month(TL, TR, t_start),
      .groups = "drop"
    )
}

#' Kernel-weighted neighbour coverage shares: share_i(t) = sum_j K_ij C_j(t) /
#' sum_j K_ij over the RAW kernel of each partition variant (build_nbr_kernel(),
#' the same masked universe as the covariate's rownorm denominator). One column
#' per rho x partition variant, named nbr_coverage_share_name().
#' @param cover_long tibble(ADM1, t_start, covered) (adm1_coverage_grid()).
#' @param D distance matrix (build_adm1_dist_km()); its rows are the universe.
#' @return tibble(ADM1, t_start, nbr_cov_lag1_r<rho>[...]).
contagion_coverage_shares <- function(
  cover_long,
  D,
  rhos,
  partition = NULL,
  partitions = NULL
) {
  variants <- list(list(part = NULL, kind = NULL, src = NULL))
  if (!is.null(partition)) {
    variants <- c(
      variants,
      list(list(part = partition, kind = "group", src = NULL))
    )
  }
  for (s in names(partitions %||% list())) {
    variants <- c(
      variants,
      list(list(part = partitions[[s]], kind = "registry", src = s))
    )
  }
  out <- NULL
  for (rho in rhos) {
    for (v in variants) {
      K <- build_nbr_kernel(D, rho, scale = "raw", partition = v$part)
      num <- neighbour_pressure(
        dplyr::mutate(cover_long, covered = covered %in% TRUE),
        K,
        normalise = TRUE,
        value_col = "covered"
      )
      nm <- nbr_coverage_share_name(rho, v$kind, v$src)
      names(num)[names(num) == "nbr"] <- nm
      out <- if (is.null(out)) {
        num
      } else {
        dplyr::left_join(out, num, by = c("ADM1", "t_start"))
      }
    }
  }
  out
}
