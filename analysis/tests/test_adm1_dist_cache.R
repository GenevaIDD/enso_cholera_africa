# test_adm1_dist_cache.R -- buildAdm1DistKm() cache coverage.
#
# HARD gate. The simulator's neighbour kernel (build_partition_nbr_kernel())
# requests the ADM1 distance matrix for every country of the partition's ENSO
# group. A country can legitimately yield no matrix rows (Ghana's GADM GID_1
# values, "GHA1_2", did until gidAdm1() keyed them). Coverage used to
# count only matrix rows, so a request including such a country rebuilt the
# matrix from the geopackage on every call -- which needs `sf` (absent on
# the cluster) and races concurrent array
# tasks writing the same file. Coverage now also trusts the sidecar's `iso3`
# (the countries the cache was BUILT for) -- but only a sidecar that records
# `gid_adm1: gidAdm1`: an older cache lists GHA with no GHA rows, and trusting it
# would give every Ghana ADM1 zero neighbour pressure.
#
# Pure filesystem: a fabricated cache + sidecar in a tempdir, and a geopackage
# path that does not exist, so any rebuild attempt fails loudly with
# "GADM geopackage not found" instead of needing `sf`.
#
# Run: Rscript analysis/tests/test_adm1_dist_cache.R   (from the project root)

suppressPackageStartupMessages({ library(stringr) })
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/gid.R")
source("analysis/utils/spatial_contagion.R")
source("analysis/tests/helpers.R")
PASS <- 0L; FAIL <- 0L

NO_GPKG <- file.path(tempdir(), "does_not_exist.gpkg")

# A 3-unit KEN/TZA matrix whose sidecar records it was built for KEN, TZA, GHA
# (GHA requested but, as in the real geopackage, contributing no rows).
mk_cache <- function(meta_iso3, gid_adm1 = "gidAdm1") {
  dir <- tempfile("adm1cache_"); dir.create(dir)
  g <- c("KEN_1", "KEN_2", "TZA_1")
  D <- matrix(c(Inf, 10, 50, 10, Inf, 40, 50, 40, Inf), 3, dimnames = list(g, g))
  f <- makeAdm1DistFile(cache_dir = dir)
  saveRDS(D, f)
  if (!is.null(meta_iso3))
    yaml::write_yaml(c(list(n_units = 3L), list(gid_adm1 = gid_adm1)[!is.null(gid_adm1)],
                       list(iso3 = as.list(meta_iso3))),
                     sub("\\.rds$", "_meta.yml", f))
  list(dir = dir, D = D)
}

rebuilds <- function(expr) {
  r <- try(expr, silent = TRUE)
  inherits(r, "try-error") && grepl("GADM geopackage not found", r)
}

cat("# [1] a country requested at build time but absent from the matrix is covered\n")
c1 <- mk_cache(c("GHA", "KEN", "TZA"))
got <- try(buildAdm1DistKm(c("KEN", "GHA"), gadm_path = NO_GPKG, cache_dir = c1$dir), silent = TRUE)
expect(!inherits(got, "try-error"), "KEN+GHA served from cache (no rebuild, no sf)")
expect(identical(got, c1$D), "returns the cached matrix unchanged")
got <- try(buildAdm1DistKm(c("TZA"), gadm_path = NO_GPKG, cache_dir = c1$dir), silent = TRUE)
expect(identical(got, c1$D), "a subset of the built set is served from cache")

cat("# [2] a country neither in the matrix nor in the sidecar still forces a rebuild\n")
expect(rebuilds(buildAdm1DistKm(c("KEN", "RWA"), gadm_path = NO_GPKG, cache_dir = c1$dir)),
       "RWA (never requested) -> rebuild attempted")

cat("# [3] no sidecar: coverage falls back to matrix rows (old behaviour)\n")
c3 <- mk_cache(NULL)
expect(identical(buildAdm1DistKm(c("KEN", "TZA"), gadm_path = NO_GPKG, cache_dir = c3$dir), c3$D),
       "matrix-row coverage still served from cache")
expect(rebuilds(buildAdm1DistKm(c("KEN", "GHA"), gadm_path = NO_GPKG, cache_dir = c3$dir)),
       "GHA without a sidecar -> rebuild attempted (no false coverage)")

cat("# [3b] a sidecar written before gidAdm1() is not trusted for coverage\n")
c3b <- mk_cache(c("GHA", "KEN", "TZA"), gid_adm1 = NULL)
expect(identical(buildAdm1DistKm(c("KEN", "TZA"), gadm_path = NO_GPKG, cache_dir = c3b$dir), c3b$D),
       "countries with matrix rows still served from cache")
expect(rebuilds(buildAdm1DistKm(c("KEN", "GHA"), gadm_path = NO_GPKG, cache_dir = c3b$dir)),
       "GHA listed only by a pre-gidAdm1 sidecar -> rebuild attempted")

cat("# [4] refresh = TRUE always rebuilds\n")
expect(rebuilds(buildAdm1DistKm(c("KEN"), gadm_path = NO_GPKG, cache_dir = c1$dir, refresh = TRUE)),
       "refresh bypasses the cache")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
