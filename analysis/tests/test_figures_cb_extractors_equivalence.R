# test_figures_cb_extractors_equivalence.R -- the ENSO crossbasis extractors lifted
# out of 09_postprocess_results.R (now 09b) into utils/figures_cb.R reproduce the derived
# manuscript artefact.
#
# HARD gate, outside --quick (it loads one deployed brms fit; no sampler runs).
# Recomputes extract_adm1_irr() / extract_country_erf() for one group of the
# occurrence track at the manuscript settings (the config's n_irr_draws,
# enso_conds and resolved best_lag_dir; seed 42) and
# all.equal()s them against the rows 09 stored in
# generated_data/results/manuscript/enso_relation_<RESULTS_SCHEMA>_track-occurrence.qs.
# SKIPs (exit 0) when the artefact or the deployed-fit archive is absent
# (laptop / CI), or when the artefact was derived with a different n_irr_draws.
#
# Run: Rscript analysis/tests/test_figures_cb_extractors_equivalence.R

suppressPackageStartupMessages({ library(dplyr) })
source("analysis/utils/results_core.R")   # + figures_cb, gid, enso, paths, provenance, config
source("analysis/utils/deploy_core.R")    # resolve_results_inputs_from_deploy(): best_lag_dir may come from the deploy

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# figures_cb extractors == 09's stored enso_relation artefact\n")

ART  <- sprintf("generated_data/results/manuscript/enso_relation_%s_track-occurrence.qs",
                RESULTS_SCHEMA[["enso_relation"]])
CFG  <- "analysis/configs/results/manuscript.yml"
# The deployed-fit archive 09 read: the config's resolved results.inputs.best_lag_dir
# (explicit, or derived from results.deploy_config), never a fixed path -- a fixed
# generated_data/enso_best_lag_fits compared the artefact with another run's fits
# once the manuscript moved to the _ep-gs(_nz) archives.
RES  <- make_options_from_config(option_list = list(), cli_argv = c("--config", CFG),
                                 model = NULL)$results
ARCH <- RES$inputs$best_lag_dir
# 09 extracts the ADM1 IRRs at results.derivation.enso_conds (results_context()'s
# ENSO_CONDS), not at figures_cb.R's fixed COMPOSITE_CONDS.
CONDS <- as.numeric(unlist(RES$derivation$enso_conds %||% COMPOSITE_CONDS))
GRP  <- "enso_west"

skip <- function(why) { cat("  SKIP  ", why, "\n"); cat("\n0 pass, 0 fail (skipped)\n"); quit(status = 0) }
if (!file.exists(ART))                       skip(paste("artefact absent:", ART))
if (!file.exists(file.path(ARCH, "MANIFEST.csv"))) skip(paste("deployed-fit archive absent:", ARCH))

n_draws <- as.integer(RES$derivation$n_irr_draws %||% 400L)
art <- qs2::qs_read(ART)
if (!GRP %in% art$adm1_irr$enso_group) skip(paste("group", GRP, "not in the artefact"))

mnf <- load_best_lag_manifest(ARCH)
ff  <- track_fit_files(mnf, "occurrence", GRP)
if (!file.exists(ff$fit) || !file.exists(ff$meta)) skip("deployed fit/meta file missing")

fit  <- readRDS(ff$fit); meta <- readRDS(ff$meta)
vm   <- meta[["enso"]]
grid <- build_onebasis_grid(vm, n_at = 200L)

new_adm1 <- extract_adm1_irr(fit, meta, enso_conds = CONDS,
                             group_label = GRP, n_draws = n_draws)
old_adm1 <- dplyr::filter(art$adm1_irr, enso_group == GRP)
expect(nrow(new_adm1) == nrow(old_adm1) && nrow(new_adm1) > 0L,
       sprintf("ADM1 IRR row count matches (%d)", nrow(old_adm1)))
expect(isTRUE(all.equal(as.data.frame(new_adm1), as.data.frame(old_adm1)[names(new_adm1)],
                        tolerance = 1e-9)),
       "extract_adm1_irr() == stored adm1_irr (all columns, tol 1e-9)")

new_cerf <- extract_country_erf(fit, meta, grid, group_label = GRP, n_draws = n_draws)
old_cerf <- dplyr::filter(art$country_erf, enso_group == GRP)
expect(nrow(new_cerf) == nrow(old_cerf) && nrow(new_cerf) > 0L,
       sprintf("country ERF row count matches (%d)", nrow(old_cerf)))
expect(isTRUE(all.equal(as.data.frame(new_cerf), as.data.frame(old_cerf)[names(new_cerf)],
                        tolerance = 1e-9)),
       "extract_country_erf() == stored country_erf (all columns, tol 1e-9)")

# the country-tier IRR extractor shares the scaffold: same countries, same draw subsample
new_cirr <- extract_country_irr(fit, meta, enso_conds = COMPOSITE_CONDS,
                                group_label = GRP, n_draws = n_draws)
expect(setequal(unique(new_cirr$country), unique(new_cerf$country)),
       "extract_country_irr() covers the same countries as extract_country_erf()")
expect(all(is.finite(new_cirr$IRR)) && all(new_cirr$lo <= new_cirr$hi),
       "extract_country_irr(): finite IRR, lo <= hi")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
