# test_gid_adm1.R -- ADM1 key from a GADM gid (gid.R::gidAdm1()).
#
# HARD gate. GADM 4.1 writes Ghana's ids without the "." after the ISO3
# ("GHA5_2", "GHA14.13_2"), so the old `str_extract(flattenGid(gid),
# "[A-Z]{3}_[0-9]+")` gave NA for the four Ghana gids of the survival panel.
# The NA ADM1 dropped Ghana from every brms fit (NA grouping rows excluded) and
# counted it as one NA "ADM1" (13.13 M) in report 11's covered population.
# Pins: the four Ghana gids map to "GHA_<n>", every other form is unchanged,
# and an unparseable gid is an error rather than an NA -- in gidAdm1(), in the
# inlined copies (build_location_frame(), occ_post.R's ADM2 level parser) and
# in build_adm1_obs_pop().
#
# Run: Rscript analysis/tests/test_gid_adm1.R   (from the project root)

suppressPackageStartupMessages({ library(dplyr); library(stringr) })
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/population.R")
source("analysis/utils/simulator_core.R")
source("analysis/utils/results_derive.R")
source("analysis/tests/helpers.R")
PASS <- 0L; FAIL <- 0L

errs <- function(expr, pattern = "no parseable ADM1") {
  r <- try(expr, silent = TRUE)
  inherits(r, "try-error") && grepl(pattern, r)
}

GHA <- c("GHA14_13_2", "GHA7_2", "GHA6_2", "GHA5_2")

cat("# [1] the four Ghana gids of the survival panel\n")
expect(identical(gidAdm1(GHA), c("GHA_14", "GHA_7", "GHA_6", "GHA_5")),
       "GHA14_13_2 / GHA7_2 / GHA6_2 / GHA5_2 -> GHA_14 / GHA_7 / GHA_6 / GHA_5")
expect(identical(gidAdm1(c("GHA14.13_2", "GHA1_2", "GHA1_1_2")),
                 c("GHA_14", "GHA_1", "GHA_1")),
       "dotted ADM2 and GADM GID_1 forms of Ghana parse the same way")

cat("# [2] canonical gids are unchanged\n")
expect(identical(gidAdm1(c("AGO.1_1", "AGO.1.2_1", "AGO_10_1", "AGO_10_3_1", "COD.12.5_1")),
                 c("AGO_1", "AGO_1", "AGO_10", "AGO_10", "COD_12")),
       "dotted and flat canonical gids -> same ADM1 as the old regex")
expect(identical(gidAdm1(c("KEN_1_1", NA)), c("KEN_1", NA)), "NA gid passes through as NA")

cat("# [3] an unparseable gid fails loudly\n")
expect(errs(gidAdm1(c("KEN_1_1", "garbage"))), "gidAdm1: 'garbage' -> error")
expect(errs(gidAdm1("gha5_2")), "gidAdm1: lower-case ISO3 -> error")
expect(errs(gidAdm1("GHA_x_2")), "gidAdm1: non-numeric ADM1 -> error")

expect(identical(gidAdm1(c("Z01.1_1", "GHA1_2"), strict = FALSE), c(NA, "GHA_1")),
       "strict = FALSE (GADM-layer map keys): NA instead of an error")

cat("# [4] build_location_frame() (inlined copy) agrees and fails loudly\n")
sd <- tibble::tibble(gid = c(GHA, "KEN_1_1"), country = c(rep("GHA", 4), "KEN"),
                     AFRO_region = c(rep("Western Africa", 4), "Eastern Africa"))
gp_file <- tempfile(fileext = ".qs")
qs2::qs_save(tibble::tibble(gid = sd$gid, gadm_level = c(2L, 1L, 1L, 1L, 1L),
                           iso3 = sd$country, year = 2020L, Q = 1e6), gp_file)
fr <- build_location_frame(sd, gid_pop_path = gp_file)
expect(identical(fr$adm1_unit[match(c(GHA, "KEN_1_1"), fr$gid)],
                 c(gidAdm1(GHA), "KEN_1")), "adm1_unit equals gidAdm1()")
expect(errs(build_location_frame(tibble::tibble(gid = "bogus", country = "KEN", AFRO_region = "E"))),
       "build_location_frame: unparseable gid -> error")

cat("# [5] occ_post.R ADM2 RE-level parser (inlined copy) on a Ghana level\n")
lvl <- c("Western Africa_GHA_GHA_14_GHA14_13_2", "Eastern Africa_KEN_KEN_1_KEN_1_2_1")
gid <- str_extract(lvl, "[A-Z]{3}_?[0-9]+_[0-9]+_[0-9]+$")
a1  <- str_extract(sub("^([A-Z]{3})([0-9])", "\\1_\\2", gid), "^[A-Z]{3}_[0-9]+")
expect(identical(gid, c("GHA14_13_2", "KEN_1_2_1")), "ADM2 gid extracted from both levels")
expect(identical(a1, c("GHA_14", "KEN_1")), "ADM1 of the ADM2 level matches gidAdm1()")

cat("# [6] build_adm1_obs_pop(): Ghana is its own ADM1s, never an NA unit\n")
surv <- tibble::tibble(
  gid         = c(GHA, "KEN_1_1"),
  admin_level = c("ADM2", "ADM1", "ADM1", "ADM1", "ADM1"),
  country     = c(rep("GHA", 4), "KEN"),
  AFRO_region = c(rep("Western Africa", 4), "Eastern Africa"))
surv$ADM1 <- gidAdm1(surv$gid)
gid_pop <- tibble::tibble(gid = surv$gid, year = 2025L, Q = c(1, 2, 3, 4, 5) * 1e6)
a <- build_adm1_obs_pop(list(POP_REF_YEAR = 2025L), surv %>% distinct(gid, country, AFRO_region),
                        surv, gid_pop, "test", 2025L, recover_dropped = TRUE)
expect(!anyNA(a$adm1_unit), "no NA adm1_unit")
expect(setequal(a$adm1_unit, c("GHA_14", "GHA_7", "GHA_6", "GHA_5", "KEN_1")),
       "one row per Ghana ADM1 + KEN_1")
bad <- surv; bad$gid[5] <- "XX9"
expect(errs(build_adm1_obs_pop(list(POP_REF_YEAR = 2025L), bad %>% distinct(gid, country, AFRO_region),
                               bad, gid_pop, "test", 2025L)),
       "build_adm1_obs_pop: unparseable gid -> error")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
