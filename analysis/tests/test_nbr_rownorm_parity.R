# test_nbr_rownorm_parity.R -- fit-side vs simulator-side neighbour covariate.
#
# HARD gate. The occurrence fit's contagion columns (buildContagionTable(),
# 03_prepare_occ.R) and the simulator kernels (build_partition_nbr_kernel() ->
# v1 / spatPomp R / spatPomp C) must yield the SAME value for the same ADM1
# activity. Both now go through buildNbrKernel(); this test pins that down for
# every scale x partition combination, including a frame that is a strict
# subset of its partition (unseen_adm1: restrict), where the row-normalising
# denominator must still be the fit-time one.
#
# Inputs: the cached ADM1 centroid distance matrix (generated_data/,
# makeAdm1DistFile()) and analysis/groupings/enso_grouping.csv. No MCMC, no fits.
#
# Run: Rscript analysis/tests/test_nbr_rownorm_parity.R   (from the project root)

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(stringr) })
source("analysis/utils/base.R")
source("analysis/utils/provenance.R")
source("analysis/utils/paths.R")
source("analysis/utils/config.R")   # load_grouping()
source("analysis/utils/gid.R")
source("analysis/utils/spatial_contagion.R")
source("analysis/tests/helpers.R")
PASS <- 0L; FAIL <- 0L
TOL  <- 1e-12

dist_file <- makeAdm1DistFile()
if (!file.exists(dist_file)) {
  cat("SKIP: ", dist_file, " not found (run 03_prepare_occ.R once).\n"); quit(status = 0L)
}
D   <- readRDS(dist_file)
reg <- load_grouping("enso_grouping")
reg <- reg[reg$iso3 %in% str_extract(rownames(D), "^[A-Z]{3}"), ]
grp <- "enso_south"
g_iso3 <- reg$iso3[reg$group == grp]
stopifnot(length(g_iso3) >= 2L)

# Synthetic lag-1 activity: every ADM1 in the matrix, 4 months, ~30% active.
set.seed(1L)
months <- seq(as.Date("2015-01-01"), by = "month", length.out = 4L)
active_long <- expand.grid(ADM1 = rownames(D), t_start = months,
                           stringsAsFactors = FALSE) %>%
  mutate(ongoing_outbreak_lag1 = runif(n()) < 0.3)

fit_tbl <- buildContagionTable(active_long, D, rhos = 500L, partition = reg)

cat("\n## 1. unpartitioned columns are unchanged by the partition argument\n")
fit_tbl0 <- buildContagionTable(active_long, D, rhos = 500L)
expect(identical(fit_tbl0, fit_tbl[, names(fit_tbl0)]),
       "raw/_rn columns byte-identical with and without `partition`")
expect(all(c("nbr_lag1_r500_pg", "nbr_lag1_r500_rn_pg") %in% names(fit_tbl)),
       "partitioned siblings written")
expect(all(fit_tbl$nbr_lag1_r500_rn_pg >= 0 & fit_tbl$nbr_lag1_r500_rn_pg <= 1 + TOL),
       "_rn_pg bounded in [0, 1]")

# Simulator frame: the group's ADM1s, minus a third of them (unseen_adm1: restrict).
g_adm1  <- rownames(D)[str_extract(rownames(D), "^[A-Z]{3}") %in% g_iso3]
frame_adm1 <- sort(sample(g_adm1, ceiling(2 * length(g_adm1) / 3)))
frame_tu <- tibble(adm1_unit = frame_adm1)
# Units outside the frame must be inactive for numerators to agree (they never
# carry simulated activity); zero them in the reference activity too.
sim_month <- months[2L]
A <- active_long %>% filter(t_start == sim_month)
A_vec <- setNames(as.numeric(A$ongoing_outbreak_lag1), A$ADM1)
A_vec[!(names(A_vec) %in% frame_adm1)] <- 0
ref_long <- active_long %>% filter(t_start == sim_month) %>%
  mutate(ongoing_outbreak_lag1 = ADM1 %in% names(A_vec)[A_vec > 0])
ref <- buildContagionTable(ref_long, D, rhos = 500L, partition = reg) %>%
  filter(ADM1 %in% frame_adm1)

cat("\n## 2. simulator kernel reproduces every fit-side column (frame subset of group)\n")
for (sc in c("raw", "rownorm")) for (pt in c("none", "group")) {
  col <- nbrCovariateName(500L, sc, pt)
  K   <- build_partition_nbr_kernel(frame_tu, 500L, scale = sc,
                                    partition = if (pt == "group") reg else NULL)
  sim <- setNames(nbr_from_adm1_active(A_vec[rownames(K)], K), rownames(K))
  want <- setNames(ref[[col]], ref$ADM1)[names(sim)]
  expect(max(abs(sim - want)) < TOL, sprintf("%s: max |sim - fit| < %g (got %.2e)",
                                             col, TOL, max(abs(sim - want))))
}

cat("\n## 3. raw unpartitioned kernel is the pre-existing one\n")
K_raw <- build_partition_nbr_kernel(frame_tu, 500L)
K_old <- kernelFromDist(D, 500L)[frame_adm1, frame_adm1]; diag(K_old) <- 0
expect(identical(unname(K_raw), unname(K_old)), "raw kernel identical to kernelFromDist subset")

cat("\n## 4. resolve_nbr_spec parses the fit's FE name\n")
mk <- function(nm) list(fe = matrix(0, 2, 2, dimnames = list(NULL, c("b_Intercept", nm))))
s0 <- resolve_nbr_spec(list(enable = FALSE), mk("nbr_lag1_r500"))
expect(!s0$use, "disabled flag -> use FALSE")
s1 <- resolve_nbr_spec(list(enable = TRUE, rho_km = 500L), mk("nbr_lag1_r500"))
expect(s1$use && s1$scale == "raw" && s1$partition == "none" && is.null(s1$registry),
       "raw term -> raw/none spec")
s2 <- resolve_nbr_spec(list(enable = TRUE, grouping_source = "enso_grouping"),
                       mk("nbr_lag1_r500_rn_pg"))
expect(s2$use && s2$scale == "rownorm" && s2$partition == "group" &&
       is.data.frame(s2$registry), "rn_pg term -> rownorm/group spec with registry")
expect(inherits(try(resolve_nbr_spec(list(enable = TRUE), mk("nbr_lag1_r500_rn_pg")),
                    silent = TRUE), "try-error"),
       "partitioned term without grouping source -> error")
expect(inherits(try(resolve_nbr_spec(list(enable = TRUE, rho_km = 200L), mk("nbr_lag1_r500")),
                    silent = TRUE), "try-error"),
       "config rho disagreeing with the fit's term -> error")
grp_arg <- list(draws_occ_by_group = list(a = mk("nbr_lag1_r500"), b = mk("nbr_lag1_r500_rn")))
expect(inherits(try(resolve_nbr_spec(list(enable = TRUE), NULL, grp_arg), silent = TRUE),
                "try-error"), "groups with different nbr terms -> error")
expect(identical(resolve_use_nbr(list(enable = TRUE), mk("nbr_lag1_r500_rn_pg")), TRUE),
       "resolve_use_nbr recognises the _pg term")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
