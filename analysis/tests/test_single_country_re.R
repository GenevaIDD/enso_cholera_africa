# test_single_country_re.R — HARD: grouped fits of one country keep their
#   country-level random effects under a tight sd prior.
#
# Usage:
#   Rscript analysis/tests/test_single_country_re.R
#
# Details:
#   No data, no fit. Why: trimming the country terms of a one-country group
#   renames the RE groups to `ADM1` / `ADM1:gid`, which the simulator /
#   smoothing extractor (simulator_core.R::.parse_brms_re_draws()) does not
#   map, so those REs would be silently zeroed at replay.

suppressPackageStartupMessages({
  library(dplyr)
})
for (m in c("base", "config", "formula", "priors", "paths", "simulator_core")) {
  source(sprintf("analysis/utils/%s.R", m))
}
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
cat("# single-country grouped fits keep their country REs\n")

one <- data.frame(
  country = "NGA",
  ADM1 = rep(c("NGA_1", "NGA_2", "NGA_3"), each = 2),
  gid = paste0("g", 1:6),
  x = rnorm(6)
)
two <- transform(one, country = rep(c("NGA", "NER"), each = 3))
re <- c(
  "(1 || country/ADM1/gid)",
  "(x - 1 || country)",
  "(x - 1 || country:ADM1)"
)

# [1] filter_re_terms_for_data(keep = "country") keeps the country level of
#     slash-nested intercepts and single-level `country` slope groups;
#     keep = character(0) trims them ((1 | ADM1[/gid]), slopes dropped).
cat("\n[1] filter_re_terms_for_data(keep=)\n")
old <- suppressMessages(filter_re_terms_for_data(re, one))
new <- suppressMessages(filter_re_terms_for_data(re, one, keep = "country"))
expect(
  identical(old, c("(1 || ADM1/gid)", "(x - 1 || country:ADM1)")),
  "keep = character(0): country trimmed / dropped (original behaviour)"
)
expect(
  identical(new, re),
  "keep = country: every term kept with its country level"
)
expect(
  identical(
    suppressMessages(filter_re_terms_for_data(re, two, keep = "country")),
    suppressMessages(filter_re_terms_for_data(re, two))
  ),
  "keep is inert when the country column has >= 2 levels"
)
one_adm <- transform(one, ADM1 = "NGA_1")
expect(
  identical(
    suppressMessages(filter_re_terms_for_data(
      "(1 | country/ADM1)",
      one_adm,
      keep = "country"
    )),
    "(1 | country)"
  ),
  "single ADM1 too: the ADM1 level is still trimmed, country kept"
)

# [2] .single_country_keep(): "country" only for grouped fits with exactly one
#     country; .flag_single_country() adds `single_country` only then.
cat("\n[2] .single_country_keep() / .flag_single_country()\n")
opt_g <- list(grouping = list(enable = TRUE))
opt_n <- list(grouping = list(enable = FALSE))
expect(
  identical(.single_country_keep(opt_g, one), "country"),
  "grouped, one country -> country"
)
expect(
  identical(.single_country_keep(opt_g, two), character(0)),
  "grouped, two countries -> none"
)
expect(
  identical(.single_country_keep(opt_n, one), character(0)),
  "not grouped -> none"
)
expect(
  identical(.single_country_keep(opt_g, NULL), character(0)),
  "no data -> none"
)
p <- list(re_terms = re)
expect(
  isTRUE(.flag_single_country(p, opt_g, one)$single_country),
  "flag set for one country"
)
expect(
  identical(.flag_single_country(p, opt_g, two), p),
  "list unchanged (no new element) otherwise"
)

# [3] .phi_re_term(keep=) keeps the country root of the shape RE.
cat("\n[3] .phi_re_term(keep=)\n")
expect(
  identical(
    suppressMessages(.phi_re_term(
      "region_country_adm1",
      TRUE,
      one,
      keep = "country"
    )),
    "(1 | country/ADM1)"
  ),
  "shape RE keeps its country root"
)
expect(
  identical(
    suppressMessages(.phi_re_term("region_country_adm1", TRUE, one)),
    "(1 | ADM1)"
  ),
  "without keep: trimmed to (1 | ADM1) (original behaviour)"
)

# [4] priors: the region/country sd tier becomes normal(0, 0.05) for flagged
#     parts (occurrence / duration / monthly), plus the shape-RE country sd;
#     unflagged parts are unchanged.
cat("\n[4] priors\n")
sd_of <- function(pr, grp, dpar = "") {
  pr <- as.data.frame(pr)
  pr$prior[pr$class == "sd" & pr$group == grp & pr$dpar == dpar]
}
dopt <- get_default_options()
dopt$grouping$enable <- TRUE
base_parts <- list(
  fixed_terms = "x",
  re_terms = c(
    "(1 | country/ADM1)",
    "(x - 1 || country)",
    "(x - 1 || country:ADM1)"
  ),
  off_terms = character(0),
  bh_terms = character(0),
  seas_terms = character(0)
)
flag_parts <- c(base_parts, list(single_country = TRUE))
for (trk in c("occ", "dur")) {
  f <- if (trk == "occ") build_occ_priors else build_duration_priors
  p0 <- suppressMessages(f(
    dopt,
    base_parts,
    list(
      prior_sd_region_country = "normal(0, 0.35)",
      prior_sd_adm1 = "normal(0, 0.35)"
    )
  ))
  p1 <- suppressMessages(f(
    dopt,
    flag_parts,
    list(
      prior_sd_region_country = "normal(0, 0.35)",
      prior_sd_adm1 = "normal(0, 0.35)"
    )
  ))
  expect(
    identical(sd_of(p0, "country"), "normal(0, 0.35)") &&
      identical(sd_of(p1, "country"), "normal(0, 0.05)"),
    sprintf("%s: country sd prior normal(0, 0.05) only when flagged", trk)
  )
  expect(
    identical(sd_of(p1, "country:ADM1"), sd_of(p0, "country:ADM1")),
    sprintf("%s: ADM1 tier unchanged", trk)
  )
}
mopt <- dopt
mopt$monthly$phi_re <- "region_country_adm1"
mparts <- function(flag) {
  x <- c(base_parts, list(phi_re = "region_country_adm1"))
  if (flag) {
    x$single_country <- TRUE
  }
  x
}
m0 <- suppressMessages(build_monthly_priors(mopt, list(), mparts(FALSE)))
m1 <- suppressMessages(build_monthly_priors(mopt, list(), mparts(TRUE)))
expect(
  identical(sd_of(m1, "country"), "normal(0, 0.05)") &&
    !identical(sd_of(m0, "country"), "normal(0, 0.05)"),
  "monthly: country sd prior normal(0, 0.05) only when flagged"
)
expect(
  identical(sd_of(m1, "country", "shape"), "normal(0, 0.05)") &&
    length(sd_of(m0, "country", "shape")) == 0L,
  "monthly: shape-RE country sd prior added only when flagged"
)

# [5] the extractor maps the kept group names (country:ADM1[:gid]) into its
#     adm1 / adm2 buckets, and does not map the old trimmed names.
cat("\n[5] simulator extractor maps the kept names\n")
S <- 4L
mk <- function(cols) {
  as.data.frame(
    setNames(replicate(length(cols), rnorm(S), simplify = FALSE), cols),
    check.names = FALSE
  )
}
kept <- .parse_brms_re_draws(mk(c(
  "b_Intercept",
  "r_country[NGA,Intercept]",
  "r_country:ADM1[NGA_NGA_1,Intercept]",
  "r_country:ADM1:gid[NGA_NGA_1_g1,Intercept]",
  "r_country:ADM1[NGA_NGA_1,enso_cb1]"
)))
expect(
  !is.null(kept$re_intercept$country) &&
    !is.null(kept$re_intercept$adm1) &&
    !is.null(kept$re_intercept$adm2),
  "country / country:ADM1 / country:ADM1:gid intercepts land in their buckets"
)
trim <- .parse_brms_re_draws(mk(c(
  "b_Intercept",
  "r_ADM1[NGA_1,Intercept]",
  "r_ADM1:gid[NGA_1_g1,Intercept]"
)))
expect(
  is.null(trim$re_intercept$adm1) && is.null(trim$re_intercept$adm2),
  "the old trimmed names (ADM1, ADM1:gid) are not mapped -- why they must not be produced"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1)
}
