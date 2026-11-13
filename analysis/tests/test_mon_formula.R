# test_mon_formula.R — verify build_monthly_formula()/build_monthly_formula_parts()
# wiring for the fourth model track (monthly outbreak-size, response
# cases_month), plus the CB slope-RE assembly that 04_fit_mon_brms_cb.R
# performs inline (mirrors the size/duration _cb scripts' pattern).
#
# Hard-fail gate. Covers:
#   (a) LHS is cases_month; offset(log(pop)) + offset(log_offset) both
#       present (offset_pop toggle only controls the pop term); no mean_rate
#       offset anywhere (that's a size-track-only apparatus).
#   (b) intercept RE is (1 | AFRO_region/country/ADM1/gid) when
#       re_adm2_intercept: true (leaf col = gid, the GAM's s(gid, bs="re")),
#       vs (1 | AFRO_region/country/ADM1) when false.
#   (c) lf(shape ~ ...) present iff phi_re != "constant".
#   (d) seasonal_term = "cyclic" + a random slope on a seasonal column errors.
#   (e) CB slope RE terms (04_fit_mon_brms_cb.R:159-183 pattern) appear iff a
#       DLNM variable's re_adm1_slope is TRUE.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R")               # model-term assembly
source("analysis/utils/paths.R")   # nbrCovariatePattern()/nbrCovariateName() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/base.R")
source("analysis/utils/dlnm.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")  # expect/fail/expect_true/expect_equal

cat("# monthly formula builder wiring\n")

# --- Minimal opt + dat builders ----------------------------------------------

mk_opt <- function(phi_re = "constant", seasonal_term = "harmonic2",
                    re_adm2_intercept = TRUE, offset_pop = TRUE,
                    dlnm = list(enable = FALSE, defaults = list(), variables = list()),
                    grouping_enable = FALSE) {
  list(
    enso     = list(lag_months = 4L, source = "oni", cache_dir = "generated_data/"),
    grouping = list(enable = grouping_enable),
    monthly  = list(
      backend        = "brms_cb",
      family         = "negbinomial", link = "log",
      admin_level    = "ADM1-ADM2",
      offset_pop     = offset_pop,
      seasonal_term  = seasonal_term,
      phi_re         = phi_re,
      prev_cases     = FALSE,
      covariates     = list(),
      random_effects = list(re_region_country_adm1_intercept = TRUE,
                            re_adm2_intercept                = re_adm2_intercept),
      brms           = list(re_uncorrelated = FALSE),
      dlnm           = dlnm
    )
  )
}

# Synthetic monthly panel — ample distinct levels at every grouping level, plus
# the enso_lag*/mso lag columns build_cb_block() needs.
n <- 32
mon_dat <- data.frame(
  cases_month = rpois(n, 5),
  pop         = rep(1e5, n),
  log_offset  = log(30),
  AFRO_region = rep(c("East", "West"), each = n / 2),
  country     = rep(c("KEN", "TZA", "NGA", "GHA"), each = n / 4),
  ADM1        = rep(paste0("A", 1:8), each = n / 8),
  gid         = paste0("G", seq_len(n)),
  month_num   = rep(1:12, length.out = n),
  mso         = abs(rnorm(n, 50, 10)),
  stringsAsFactors = FALSE
)
enso_lags <- as.data.frame(matrix(rnorm(n * 13), nrow = n))
names(enso_lags) <- sprintf("enso_lag%d", 0:12)
mon_dat <- dplyr::bind_cols(mon_dat, enso_lags)

# Helper: extract the shape (dispersion) sub-formula text, or NULL. brms names
# this parameter "shape" for negbinomial (the lf() call in build_monthly_formula
# is literally `shape ~ 1 + ...`), not "phi" -- that key belongs to other
# families (e.g. Beta). Confirmed empirically via bf_obj$pforms.
phi_formula_text <- function(bf_obj) {
  pf <- bf_obj$pforms
  if (is.null(pf) || is.null(pf[["shape"]])) return(NULL)
  paste(deparse(pf[["shape"]]), collapse = " ")
}

# --- 1. LHS + offsets ---------------------------------------------------------
cat("\n[1] LHS is cases_month; offsets present; no mean_rate\n")

opt <- mk_opt()
bf  <- build_monthly_formula(opt, mon_dat)
form_txt <- paste(deparse(bf$formula), collapse = " ")
expect(grepl("^cases_month\\s*~", form_txt), "LHS is cases_month")
expect(grepl("offset(log(pop))", form_txt, fixed = TRUE), "offset(log(pop)) present")
expect(grepl("offset(log_offset)", form_txt, fixed = TRUE), "offset(log_offset) present (unconditional)")
expect(!grepl("mean_rate", form_txt, fixed = TRUE), "no mean_rate offset (size-track-only)")

opt_no_pop <- mk_opt(offset_pop = FALSE)
form_txt_no_pop <- paste(deparse(build_monthly_formula(opt_no_pop, mon_dat)$formula), collapse = " ")
expect(!grepl("offset(log(pop))", form_txt_no_pop, fixed = TRUE),
       "offset_pop: false drops offset(log(pop))")
expect(grepl("offset(log_offset)", form_txt_no_pop, fixed = TRUE),
       "offset(log_offset) stays present even when offset_pop is off")

# --- 2. Intercept RE nesting ---------------------------------------------------
cat("\n[2] intercept RE nesting depends on re_adm2_intercept\n")

parts_gid <- build_monthly_formula_parts(mk_opt(re_adm2_intercept = TRUE), mon_dat)
expect(any(grepl("AFRO_region/country/ADM1/gid", parts_gid$re_terms, fixed = TRUE)),
       "re_adm2_intercept: true nests down to gid")

parts_no_gid <- build_monthly_formula_parts(mk_opt(re_adm2_intercept = FALSE), mon_dat)
expect(!any(grepl("/gid", parts_no_gid$re_terms, fixed = TRUE)),
       "re_adm2_intercept: false stops at ADM1 (no /gid)")
expect(any(grepl("AFRO_region/country/ADM1", parts_no_gid$re_terms, fixed = TRUE)),
       "re_adm2_intercept: false still nests to ADM1")

# --- 3. phi sub-formula presence -----------------------------------------------
cat("\n[3] lf(shape ~ ...) present iff phi_re != constant\n")

bf_const <- build_monthly_formula(mk_opt(phi_re = "constant"), mon_dat)
expect(is.null(phi_formula_text(bf_const)), "phi_re = constant: no phi sub-formula")

for (pr in c("region", "region_country", "region_country_adm1")) {
  bf_pr <- build_monthly_formula(mk_opt(phi_re = pr), mon_dat)
  txt   <- phi_formula_text(bf_pr)
  expect(!is.null(txt), sprintf("phi_re = %s: phi sub-formula present", pr))
}

# --- 4. cyclic seasonal_term + seasonal RE slope errors ------------------------
cat("\n[4] cyclic seasonal_term + a seasonal RE slope errors\n")

opt_cyclic_bad <- mk_opt(seasonal_term = "cyclic")
opt_cyclic_bad$monthly$random_effects$slopes <- list(sin_month = list(adm1 = TRUE))
threw <- tryCatch({ build_monthly_formula_parts(opt_cyclic_bad, mon_dat); FALSE },
                   error = function(e) TRUE)
expect(isTRUE(threw), "cyclic + RE slope on a seasonal column raises an error")

opt_cyclic_ok <- mk_opt(seasonal_term = "cyclic")
threw_ok <- tryCatch({ build_monthly_formula_parts(opt_cyclic_ok, mon_dat); FALSE },
                      error = function(e) TRUE)
expect(!isTRUE(threw_ok), "cyclic without a seasonal RE slope is fine")

# --- 5. CB slope RE terms appear iff re_adm1_slope ------------------------------
cat("\n[5] CB slope RE terms appear iff a DLNM variable's re_adm1_slope is TRUE\n")

# Mirrors 04_fit_mon_brms_cb.R:130-183 -- the CB-term-append logic is inline in
# that fit script (not a shared function), so this reproduces the same call
# sequence against the real shared helpers (build_cb_block(), re_group_by_level(),
# nested_slope_re_terms(), filter_re_terms_for_data()) rather than re-deriving it.
.assemble_cb_parts <- function(opt, dat) {
  dlnm_cfg    <- opt$monthly$dlnm %||% list()
  drop_region <- isTRUE(opt$grouping$enable)
  cbb     <- build_cb_block(dat, dlnm_cfg, " (monthly model, test)")
  dat     <- cbb$dat
  cb_cols <- cbb$cb_cols
  cb_meta <- cbb$cb_meta

  parts <- build_monthly_formula_parts(opt, dat)
  if (length(cb_cols) == 0L) return(parts)

  dlnm_defaults <- dlnm_cfg$defaults %||% list()
  for (v in names(dlnm_cfg$variables)) {
    var_spec <- dlnm_cfg$variables[[v]] %||% list()
    klt <- var_spec$keep_linear_term %||% dlnm_defaults$keep_linear_term %||% FALSE
    if (!isTRUE(klt)) {
      parts$fixed_terms <- parts$fixed_terms[parts$fixed_terms != v]
      parts$re_terms    <- parts$re_terms[!grepl(sprintf("^\\(%s[[:space:]]", v), parts$re_terms)]
    }
  }
  parts$fixed_terms <- c(parts$fixed_terms, cb_cols)

  cb_groups <- re_group_by_level("gid", drop_region = drop_region)
  re_sep    <- "|"
  for (v in names(dlnm_cfg$variables)) {
    var_spec <- dlnm_cfg$variables[[v]] %||% list()
    re_slope <- var_spec$re_adm1_slope %||% (dlnm_defaults$re_adm1_slope %||% FALSE)
    if (!isTRUE(re_slope)) next
    v_cb_names <- cb_meta[[v]]$cb_col_names
    if (is.null(v_cb_names)) next
    re_levels <- as.character(unlist(
      var_spec$re_slope_levels %||% dlnm_defaults$re_slope_levels %||%
        c("region", "country", "adm1")))
    bad_levels <- setdiff(re_levels, names(cb_groups))
    if (length(bad_levels) > 0L)
      stop("monthly.dlnm variable '", v, "': unknown re_slope_levels entry/entries ",
           paste(bad_levels, collapse = ", "), "; allowed: ",
           paste(names(cb_groups), collapse = ", "))
    new_cb_re <- filter_re_terms_for_data(
      nested_slope_re_terms(v_cb_names, cb_groups, re_sep, levels = re_levels), dat
    )
    parts$re_terms <- c(parts$re_terms, new_cb_re)
  }
  parts
}

dlnm_off <- list(
  enable = TRUE,
  defaults = list(max_lag = 3L, argvar = list(fun = "bs", degree = 3L, df = 4L),
                  keep_linear_term = FALSE, re_adm1_slope = FALSE, single_lag = 0L),
  variables = list(mso = list(single_lag = 0L, re_adm1_slope = FALSE))
)
parts_off <- .assemble_cb_parts(mk_opt(dlnm = dlnm_off), mon_dat)
expect(!any(grepl("mso_cb", parts_off$re_terms)), "re_adm1_slope: FALSE -> no CB slope RE terms")
expect(any(grepl("^mso_cb", parts_off$fixed_terms)), "CB fixed-effect columns still appended")

dlnm_on <- list(
  enable = TRUE,
  defaults = list(max_lag = 3L, argvar = list(fun = "bs", degree = 3L, df = 4L),
                  keep_linear_term = FALSE, re_adm1_slope = FALSE, single_lag = 0L),
  variables = list(mso = list(single_lag = 0L, re_adm1_slope = TRUE))
)
parts_on <- .assemble_cb_parts(mk_opt(dlnm = dlnm_on), mon_dat)
expect(any(grepl("mso_cb", parts_on$re_terms)), "re_adm1_slope: TRUE -> CB slope RE terms appear")

# --- 6. re_slope_levels selects the CB slope nesting levels ---------------------
cat("\n[6] re_slope_levels drops the country:ADM1 CB slope level\n")

# Grouping mode (the production setting): re_group_by_level(drop_region = TRUE)
# gives country / country:ADM1 / country:ADM1:gid.
mk_dlnm_lv <- function(mso_spec) list(
  enable = TRUE,
  defaults = list(max_lag = 3L, argvar = list(fun = "bs", degree = 3L, df = 4L),
                  keep_linear_term = FALSE, re_adm1_slope = TRUE, single_lag = 0L),
  variables = list(mso = mso_spec, enso = list(single_lag = 4L))
)

mso_slopes <- function(parts) grep("^\\(mso_cb", parts$re_terms, value = TRUE)

parts_lv_default <- .assemble_cb_parts(
  mk_opt(dlnm = mk_dlnm_lv(list(single_lag = 0L)), grouping_enable = TRUE), mon_dat)
sl_default <- mso_slopes(parts_lv_default)
expect(any(grepl("| country)", sl_default, fixed = TRUE)),
       "default levels -> country-level mso slopes present")
expect(any(grepl("| country:ADM1)", sl_default, fixed = TRUE)),
       "default levels -> country:ADM1 mso slopes present (unchanged behaviour)")

parts_lv_country <- .assemble_cb_parts(
  mk_opt(dlnm = mk_dlnm_lv(list(single_lag = 0L,
                                re_slope_levels = c("region", "country"))),
         grouping_enable = TRUE), mon_dat)
sl_country <- mso_slopes(parts_lv_country)
expect(any(grepl("| country)", sl_country, fixed = TRUE)),
       "re_slope_levels [region, country] -> country-level mso slopes kept")
expect(!any(grepl("| country:ADM1)", sl_country, fixed = TRUE)),
       "re_slope_levels [region, country] -> no country:ADM1 mso slopes")
expect(any(grepl("| country:ADM1)", parts_lv_country$re_terms, fixed = TRUE) &
             grepl("^\\(enso_cb", parts_lv_country$re_terms)),
       "enso keeps its country:ADM1 slopes (per-variable, not global)")

threw_lv <- tryCatch({
  .assemble_cb_parts(
    mk_opt(dlnm = mk_dlnm_lv(list(single_lag = 0L, re_slope_levels = c("ADM1"))),
           grouping_enable = TRUE), mon_dat)
  FALSE
}, error = function(e) TRUE)
expect(isTRUE(threw_lv), "an unknown re_slope_levels entry errors instead of silently dropping slopes")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
