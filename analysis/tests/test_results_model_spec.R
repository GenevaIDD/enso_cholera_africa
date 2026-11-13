# test_results_model_spec.R -- the fit-level model specification of the
# deployed fits (report_model_spec, table S34;
# utils/results_derive.R: model_spec_of(), model_spec_check(),
# model_spec_tables()).
#
# HARD gate, quick-safe: no fit is loaded or sampled. Fit-like lists (the
# brmsfit fields model_spec_of() reads: formula$formula / $pforms, family,
# prior, coefficient names) are built from analysis/model_terms_deployed.yml
# itself, so a fit carrying exactly the documented terms must check clean.
#   [1] model_spec_of(): fixed / RE / offset terms, dispersion terms, RE groups
#       (a/b nesting expanded), family, link, response, coefficients, priors.
#   [2] model_spec_check() against the yml, per track: all ok for the
#       documented fits ({k} expanded, {nbr} matched to the fit's column); a
#       dropped term, an added climate term and a wrong family are flagged.
#   [3] model_spec_tables(): one summary row per fit, terms / priors / check
#       tables keyed by track, group and lag.
#   [4] registry and paths: the convergence step writes report_model_spec too;
#       its artefact path is track-keyed.
#
# Run: Rscript analysis/tests/test_results_model_spec.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
source("analysis/utils/results_core.R")
source("analysis/utils/results_derive.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")

yml <- yaml::read_yaml(MODEL_TERMS_DEPLOYED)
NBR <- "nbr_lag1_r500_rn_parthmmclusters"
ex <- function(entries) unlist(lapply(entries %||% list(), function(e)
  if (is.list(e)) vapply(seq_len(as.integer(e$count)), function(k) gsub("{k}", k, e$term, fixed = TRUE), character(1))
  else gsub("{nbr}", NBR, as.character(e), fixed = TRUE)))
fake_fit <- function(track, drop = character(0), add = character(0), family = NULL) {
  t <- yml$tracks[[track]]
  rhs <- setdiff(c(ex(t$baseline), ex(t$fixed_effects), ex(t$re_terms), add, ex(t$offsets)), drop)
  rhs <- gsub("{nbr}", NBR, rhs, fixed = TRUE)
  pforms <- if (length(t$shape_terms)) list(shape = stats::as.formula(paste("shape ~", paste(ex(t$shape_terms), collapse = " + "))))
  # yaml reads the bare `response: y` as the boolean TRUE (YAML 1.1).
  resp <- if (isTRUE(t$response)) "y" else as.character(t$response)
  list(formula = list(formula = stats::as.formula(paste(resp, "~", paste(rhs, collapse = " + "))), pforms = pforms),
       family = list(family = family %||% t$family, link = t$link),
       prior = data.frame(prior = c("normal(0, 1.5)", "normal(0, 0.5)", "normal(0, 2.5)"),
                          class = c("Intercept", "b", "sd"), coef = "", group = c("", "", "country"),
                          dpar = "", source = "user", stringsAsFactors = FALSE),
       variables = c("b_Intercept", paste0("b_", ex(t$fixed_effects)), "sd_country__Intercept", "lp__"))
}

cat("# [1] model_spec_of()\n")
so <- model_spec_of(fake_fit("occurrence"))
expect(so$response == "y" && so$family == "poisson" && so$link == "log", "response, family, link")
expect(setequal(so$terms$term[so$terms$kind == "offset"], c("offset(log_offset)", "offset(log(base_rate_occ))")),
       "offsets read from the formula")
expect(all(c("(1||country/ADM1/gid)", .norm_term(sprintf("(%s - 1 || country:ADM1)", NBR))) %in%
             so$terms$term[so$terms$kind == "re"]), "RE terms, normalised")
expect(setequal(so$re_groups, c("country", "country:ADM1", "country:ADM1:gid")), "RE groups, nesting expanded")
expect("s(t_start_seas,bs='cc',k=8)" %in% so$terms$term[so$terms$kind == "fixed"], "the cyclic baseline spline is a fixed term")
expect(NBR %in% so$fixef && !"Intercept" %in% so$fixef[0] && nrow(so$priors) == 3L, "coefficients (b_ prefix dropped) and priors")
sm <- model_spec_of(fake_fit("monthly"))
expect(identical(unique(sm$terms$component), c("mu", "shape")) &&
         "(1|country/ADM1)" %in% sm$terms$term[sm$terms$component == "shape"], "dispersion (shape) terms")

cat("# [2] model_spec_check() vs model_terms_deployed.yml\n")
for (tr in names(yml$tracks)) {
  ck <- model_spec_check(model_spec_of(fake_fit(tr)), .deployed_terms_expected(yml$tracks[[tr]]))
  expect(all(ck$ok), sprintf("%s: the documented fit checks clean (%s)", tr,
                             paste(ck$item[!ck$ok], collapse = ", ")))
}
ck <- model_spec_check(model_spec_of(fake_fit("occurrence", drop = "ongoing_outbreak_lag1", add = "mean_precip")),
                       .deployed_terms_expected(yml$tracks$occurrence))
fx <- filter(ck, item == "fixed")
expect(!fx$ok && fx$missing == "ongoing_outbreak_lag1" && fx$unexpected == "mean_precip",
       "a dropped term and an added climate term are flagged")
ckf <- model_spec_check(model_spec_of(fake_fit("duration", family = "negbinomial")), .deployed_terms_expected(yml$tracks$duration))
expect(!filter(ckf, item == "family")$ok && all(filter(ckf, item != "family")$ok), "a wrong family is flagged alone")

cat("# [3] model_spec_tables()\n")
specs <- list(enso_east = model_spec_of(fake_fit("occurrence")),
              enso_west = model_spec_of(fake_fit("occurrence", add = "mean_precip")))
tb <- model_spec_tables(specs, "occurrence", list(enso_east = 0L, enso_west = -6L), yml)
expect(nrow(tb$fits) == 2L && identical(tb$fits$matches_deployed_terms, c(TRUE, FALSE)) &&
         identical(tb$fits$lag, c(0L, -6L)), "one summary row per fit, the yml check per fit, the deployed lag")
expect(all(c("track", "enso_group", "lag") %in% names(tb$terms)) && all(tb$priors$track == "occurrence") &&
         nrow(tb$fixef) == 2L * length(specs$enso_east$fixef), "terms / priors / fixef keyed by track, group, lag")
expect(!nrow(model_spec_tables(specs, "occurrence", list(), NULL)$check), "no yml: no check table")

cat("# [4] registry and paths\n")
row <- Filter(function(r) identical(r$stem, "convergence"), .RESULTS_STEP_TABLE)[[1]]
expect(identical(row$kinds, c("convergence", "report_model_spec")) && !is.null(RESULTS_SCHEMA[["report_model_spec"]]),
       "the convergence step writes report_model_spec (schema registered)")
opt <- list(outdir = "generated_data/")
p <- makeResultsRunFile(opt, "results_report_model_spec", track = "monthly", schema = "v1",
                        config_path = "analysis/configs/results/manuscript.yml")
expect(grepl("report_model_spec_v1_track-monthly\\.qs$", p), "track-keyed artefact path")
expect(inherits(try(makeResultsRunFile(opt, "results_report_model_spec", schema = "v1",
                                       config_path = "analysis/configs/results/manuscript.yml"), silent = TRUE),
                "try-error"), "the track is required")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
