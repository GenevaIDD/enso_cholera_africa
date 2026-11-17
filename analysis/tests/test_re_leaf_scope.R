# test_re_leaf_scope.R — HARD: `random_effects.re_adm2_leaf_scope` (size /
#   monthly / duration) and duration's ADM1-ADM2 episode set.
#
# Usage:
#   Rscript analysis/tests/test_re_leaf_scope.R
#
# Details:
#   The deepest intercept RE of the monthly and size models is keyed on
#   `gid`, which mixes ADM2 and ADM1-only units. Under `all`, an ADM1-only
#   unit gets a singleton leaf under its own ADM1 node (two REs for one
#   geography). Under `adm2_only` the nested intercept stops at ADM1 and the
#   leaf is a separate `(0 + is_adm2 | ...:ADM1:gid_leaf)` term that only
#   canonical-ADM2 rows load on; the simulator must read it back into its
#   adm2 bucket. Also duration_episode_rows(). Pure synthetic; no MCMC.

suppressPackageStartupMessages({
  library(dplyr)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/priors.R") # brms prior builders
source("analysis/utils/paths.R")
source("analysis/utils/base.R")
source("analysis/utils/dlnm.R")
source("analysis/utils/gid.R") # apply_admin_level_filter() for duration_episode_rows()

# Simulator RE parser, extracted by parsing (no spatPomp / brms needed).
local({
  ex <- parse("analysis/utils/simulator_core.R", keep.source = FALSE)
  want <- c(".parse_brms_re_draws", ".normalize_occ_term", "re_level_key")
  for (e in as.list(ex)) {
    if (
      is.call(e) &&
        identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% want
    ) {
      eval(e, envir = globalenv())
    }
  }
})

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
throws <- function(expr, pattern = NULL) {
  r <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  !is.null(r) && (is.null(pattern) || grepl(pattern, r))
}

mk_mon_opt <- function(scope = "all", grouping_enable = FALSE) {
  list(
    enso = list(lag_months = 4L, source = "oni", cache_dir = "generated_data/"),
    grouping = list(enable = grouping_enable),
    monthly = list(
      backend = "brms_cb",
      family = "negbinomial",
      link = "log",
      admin_level = "ADM1-ADM2",
      offset_pop = TRUE,
      seasonal_term = "harmonic2",
      phi_re = "constant",
      prev_cases = FALSE,
      covariates = list(),
      random_effects = list(
        re_region_country_adm1_intercept = TRUE,
        re_adm2_intercept = TRUE,
        re_adm2_leaf_scope = scope
      ),
      brms = list(re_uncorrelated = FALSE),
      dlnm = list(enable = FALSE, defaults = list(), variables = list())
    )
  )
}

n <- 32
dat <- data.frame(
  cases_month = rpois(n, 5),
  pop = 1e5,
  log_offset = log(30),
  AFRO_region = rep(c("East", "West"), each = n / 2),
  country = rep(c("KEN", "TZA", "NGA", "GHA"), each = n / 4),
  ADM1 = rep(paste0("A", 1:8), each = n / 8),
  gid = paste0("G", seq_len(n)),
  admin_level = rep(c("ADM1", "ADM2", "ADM2", "ADM2"), length.out = n),
  admin_level_canonical = rep(
    c("ADM1", "ADM2", "ADM2", "ADM2"),
    length.out = n
  ),
  month_num = rep(1:12, length.out = n),
  stringsAsFactors = FALSE
)

cat("# re_adm2_leaf_scope\n")

# size's `ADM2` column is a copy of gid (utils/siz_data.R). Under `all` the
# singleton leaves of ADM1-only units leave only the sum of the two REs
# identified, and the singletons still feed sd(leaf).
cat("\n[1] default 'all' is the historical nested leaf\n")
p_all <- build_monthly_formula_parts(mk_mon_opt("all"), dat)
expect(
  any(grepl(
    "(1 | AFRO_region/country/ADM1/gid)",
    p_all$re_terms,
    fixed = TRUE
  )),
  "monthly: (1 | AFRO_region/country/ADM1/gid)"
)
expect(!any(grepl("is_adm2", p_all$re_terms, fixed = TRUE)), "no is_adm2 term")
opt_nokey <- mk_mon_opt("all")
opt_nokey$monthly$random_effects$re_adm2_leaf_scope <- NULL
expect(
  identical(
    build_monthly_formula_parts(opt_nokey, dat)$re_terms,
    p_all$re_terms
  ),
  "absent key == 'all'"
)

cat(
  "\n[2] 'adm2_only': nested intercept stops at ADM1, separate ADM2-only leaf\n"
)
d2 <- .leaf_columns(dat, canonical_ok = TRUE)
expect(
  identical(d2$is_adm2, as.integer(dat$admin_level == "ADM2")) &&
    all(d2$gid_leaf[d2$is_adm2 == 0L] == "none") &&
    all(d2$gid_leaf[d2$is_adm2 == 1L] == d2$gid[d2$is_adm2 == 1L]),
  ".leaf_columns(): is_adm2 / gid_leaf"
)
p_adm2 <- build_monthly_formula_parts(mk_mon_opt("adm2_only"), d2)
expect(
  any(p_adm2$re_terms == "(1 | AFRO_region/country/ADM1)"),
  "monthly: (1 | AFRO_region/country/ADM1)"
)
expect(
  any(p_adm2$re_terms == "(0 + is_adm2 | AFRO_region:country:ADM1:gid_leaf)"),
  "monthly: (0 + is_adm2 | AFRO_region:country:ADM1:gid_leaf)"
)
expect(
  !any(grepl("/gid)", p_adm2$re_terms, fixed = TRUE)),
  "no nested /gid leaf"
)
p_grp <- build_monthly_formula_parts(
  mk_mon_opt("adm2_only", grouping_enable = TRUE),
  d2
)
expect(
  any(p_grp$re_terms == "(0 + is_adm2 | country:ADM1:gid_leaf)"),
  "grouped: (0 + is_adm2 | country:ADM1:gid_leaf)"
)
tiers <- .re_group_tiers(p_grp$re_terms)
expect(
  "country:ADM1:gid_leaf" %in% tiers$adm1,
  "the leaf takes the ADM1-tier sd prior"
)

cat("\n[3] the leaf needs the geometry-derived admin level\n")
expect(
  throws(.leaf_columns(dat, canonical_ok = FALSE), "canonical_admin_level"),
  ".leaf_columns() refuses corpus labels"
)
expect(
  throws(.leaf_scope_of(list(
    random_effects = list(re_adm2_leaf_scope = "gid")
  ))),
  "an unknown scope is an error"
)

cat("\n[4] size builder routes through the same helper\n")
opt_s <- list(
  enso = list(lag_months = 4L),
  grouping = list(enable = TRUE),
  size = list(
    backend = "brms",
    admin_level = "ADM1-ADM2",
    covariates = list(),
    phi_re = "constant",
    random_effects = list(
      re_region_country_adm1_intercept = TRUE,
      re_adm2_intercept = TRUE,
      re_adm2_leaf_scope = "adm2_only"
    ),
    brms = list(re_uncorrelated = FALSE)
  )
)
ps <- tryCatch(
  build_size_formula_parts(opt_s, d2 %>% mutate(ADM2 = gid)),
  error = function(e) e
)
expect(
  !inherits(ps, "error") &&
    any(ps$re_terms == "(1 | country/ADM1)") &&
    any(ps$re_terms == "(0 + is_adm2 | country:ADM1:gid_leaf)"),
  "size: same two terms under grouping"
)

cat("\n[5] fit filenames carry the scope (no file.exists() collision)\n")
o_all <- c(mk_mon_opt("all"), list(baseline = "1990_2019", country = "all"))
o_adm2 <- c(
  mk_mon_opt("adm2_only"),
  list(baseline = "1990_2019", country = "all")
)
f_all <- make_monthly_run_file(o_all, "mon_fit")
f_adm2 <- make_monthly_run_file(o_adm2, "mon_fit")
expect(
  !grepl("leaf-adm2", f_all) &&
    grepl("_leaf-adm2", f_adm2) &&
    identical(sub("_leaf-adm2", "", f_adm2), f_all),
  "only adm2_only adds _leaf-adm2"
)

cat("\n[6] the simulator reads the is_adm2 leaf back into its adm2 bucket\n")
draws_df <- data.frame(
  check.names = FALSE,
  `.draw` = 1:2,
  `b_Intercept` = c(1, 1),
  `r_country[KEN,Intercept]` = c(0.1, 0.1),
  `r_country:ADM1[KEN_KEN_1,Intercept]` = c(0.2, 0.2),
  `r_country:ADM1:gid_leaf[KEN_KEN_1_KEN_1_2_1,is_adm2]` = c(0.5, 0.7),
  `r_country:ADM1:gid_leaf[KEN_KEN_1_none,is_adm2]` = c(0.9, 0.9)
)
pr <- .parse_brms_re_draws(draws_df)
expect(
  !is.null(pr$re_intercept$adm2) &&
    "KEN_KEN_1_KEN_1_2_1" %in% colnames(pr$re_intercept$adm2),
  "is_adm2 leaf parsed as the adm2 intercept"
)
expect(length(pr$re_slope) == 0L, "and NOT as a slope")
frame <- data.frame(
  country = c("KEN", "KEN"),
  adm1_unit = c("KEN_1", "KEN_1"),
  gid = c("KEN_1_2_1", "KEN_1_1"),
  stringsAsFactors = FALSE
)
k <- re_level_key(frame, "adm2", drop_region = TRUE)
expect(
  k[1] %in%
    colnames(pr$re_intercept$adm2) &&
    !(k[2] %in% colnames(pr$re_intercept$adm2)),
  "an ADM2 unit matches its leaf; an ADM1-only unit matches none (RE 0)"
)

cat(
  "\n[7] duration: ADM1-ADM2 episode set + the same leaf, defaults unchanged\n"
)
mk_dur_opt <- function(level = "ADM1", adm2i = NULL, scope = NULL) {
  re <- list(re_region_country_adm1_intercept = TRUE)
  re$re_adm2_intercept <- adm2i
  re$re_adm2_leaf_scope <- scope
  list(
    enso = list(lag_months = 4L),
    grouping = list(enable = TRUE),
    duration = list(
      backend = "brms_cb",
      admin_level = level,
      basehaz_duration = "factor",
      basehaz_poisson = "none",
      covariates = list(),
      random_effects = re,
      brms_cb = list(re_uncorrelated = TRUE)
    )
  )
}
dd <- dat %>%
  mutate(y = 0L, mso_band = factor(rep(c("1", "2"), length.out = n)))
pd0 <- build_duration_formula_parts(mk_dur_opt(), dd)
expect(
  any(pd0$re_terms == "(1 | country/ADM1)") && !any(grepl("gid", pd0$re_terms)),
  "duration default: (1 | country/ADM1), no gid term"
)
expect(
  identical(
    build_duration_formula_parts(
      mk_dur_opt(adm2i = FALSE, scope = "all"),
      dd
    )$re_terms,
    pd0$re_terms
  ),
  "explicit defaults (re_adm2_intercept false, scope all) == absent keys"
)
dd2 <- .leaf_columns_for(
  dd,
  mk_dur_opt("ADM1-ADM2", TRUE, "adm2_only")$duration,
  list(key = "gid_season"),
  tag = "test"
)
pd2 <- build_duration_formula_parts(
  mk_dur_opt("ADM1-ADM2", TRUE, "adm2_only"),
  dd2
)
expect(
  any(pd2$re_terms == "(1 | country/ADM1)") &&
    any(pd2$re_terms == "(0 + is_adm2 | country:ADM1:gid_leaf)"),
  "duration ADM1-ADM2 + adm2_only: (1 | country/ADM1) + (0 + is_adm2 | country:ADM1:gid_leaf)"
)
expect(
  "country:ADM1:gid_leaf" %in% .re_group_tiers(pd2$re_terms)$adm1,
  "duration: the leaf takes the ADM1-tier sd prior"
)
expect(
  identical(
    .leaf_columns_for(dd, mk_dur_opt()$duration, list(key = "gid_season")),
    dd
  ),
  ".leaf_columns_for(): no-op under the default scope"
)

# duration_episode_rows() (dur_data.R), extracted by parsing: dur_data.R's own
# source closure is heavy and the function needs only gid.R.
local({
  ex <- parse("analysis/utils/dur_data.R", keep.source = FALSE)
  for (e in as.list(ex)) {
    if (
      is.call(e) &&
        identical(e[[1]], as.name("<-")) &&
        identical(as.character(e[[2]]), "duration_episode_rows")
    ) {
      eval(e, envir = globalenv())
    }
  }
})
ep <- data.frame(
  episode_id = 1:5,
  ADM1 = c("A1", "A1", "A2", "A3", "A3"),
  admin_level = c("ADM1", "ADM2", "ADM1", "ADM2", "ADM2"),
  stringsAsFactors = FALSE
)
expect(
  identical(
    duration_episode_rows(ep, list(admin_level = "ADM1-ADM2"))$episode_id,
    2:5
  ),
  "ADM1-ADM2 drops the ADM1 parent with ADM2 rows, keeps an ADM1-only unit"
)
expect(
  identical(duration_episode_rows(ep, list(admin_level = "ADM1")), ep) &&
    identical(duration_episode_rows(ep, list(admin_level = "ADM2")), ep) &&
    identical(duration_episode_rows(ep, NULL), ep),
  "ADM1 / ADM2 / absent: no row filter (historical duration semantics)"
)
expect(
  throws(duration_episode_rows(ep, list(admin_level = "ADM0")), "admin_level"),
  "an unknown duration admin_level is an error"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
