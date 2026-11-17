# test_cb_singlelag.R — SOFT: single-lag DLNM (onebasis) machinery.
#
# Usage:
#   Rscript analysis/tests/test_cb_singlelag.R
#
# Details:
#   Asserts: dlnm::onebasis on a synthetic ENSO column yields `df` basis
#   columns; extract_cb_surface's single-lag branch returns a finite ERF with
#   `single_lag` set (round-trip through a reconstructed onebasis stub);
#   make_run_file gives distinct paths for single_lag = 1 vs 2 (occurrence and
#   size tracks), so per-lag fits never overwrite each other.

suppressPackageStartupMessages({
  library(dlnm)
  library(splines)
  library(posterior)
  library(brms)
})

source("analysis/utils/config.R")
source("analysis/utils/formula.R") # model-term assembly
source("analysis/utils/base.R") # transitive dependency of dlnm.R, paths.R (utils/README.md)
source("analysis/utils/gid.R") # transitive dependency of clim_compile.R, clim_db.R (utils/README.md)
source("analysis/utils/paths.R") # make_run_file() -- used directly
source("analysis/utils/clim_compile.R") # transitive dependency of plots.R (utils/README.md)
source("analysis/utils/clim_db.R") # transitive dependency of plots.R (utils/README.md)
source("analysis/utils/provenance.R") # transitive dependency of paths.R (utils/README.md)
source("analysis/utils/occ_post.R") # adm1_lineage_levels() -- used directly
source("analysis/utils/dlnm.R") # extract_cb_surface() -- used directly
source("analysis/utils/plots.R") # transitive dependency of occ_post.R (utils/README.md)

source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

# 1) onebasis dimensionality -------------------------------------------------
set.seed(42)
x <- rnorm(500)
argvar <- list(fun = "bs", degree = 3L, df = 4L)
ob <- do.call(dlnm::onebasis, c(list(x), argvar))
if (ncol(ob) != 4L) {
  fail(sprintf("onebasis produced %d cols, expected df = 4", ncol(ob)))
}

# 2) extract_cb_surface single-lag round-trip ----------------------------------
# Build a cb_meta-style entry (as the run scripts store it) and a draws_matrix
# stand-in fit (extract_cb_surface only calls brms::as_draws_matrix(fit, ...)).
ob_attrs <- attributes(ob)
ob_attrs[c("dim", "dimnames", "class")] <- NULL
n_cols <- ncol(ob)
cb_names <- sprintf("enso_cb%d", seq_len(n_cols))
meta <- list(
  enso = list(
    single_lag = 3L,
    onebasis_attrs = ob_attrs,
    argvar = argvar,
    range = range(x),
    df = n_cols,
    cb_col_names = cb_names
  )
)
# Posterior coef draws on the basis columns (b_<col>), as a draws_matrix.
b_mat <- matrix(
  rnorm(200 * n_cols, 0, 0.3),
  nrow = 200,
  ncol = n_cols,
  dimnames = list(NULL, sprintf("b_%s", cb_names))
)
fit_stub <- posterior::as_draws_matrix(b_mat)

at_grid <- seq(range(x)[1], range(x)[2], length.out = 21L)
surfs <- tryCatch(
  extract_cb_surface(fit_stub, meta, at = at_grid),
  error = function(e) {
    fail(paste("extract_cb_surface errored:", conditionMessage(e)))
    NULL
  }
)
surf <- surfs$enso
if (is.null(surf$single_lag) || surf$single_lag != 3L) {
  fail("extract_cb_surface did not propagate single_lag = 3")
}
if (is.null(surf$matfit) || !all(is.finite(surf$matfit))) {
  fail("extract_cb_surface returned non-finite matfit for single-lag branch")
}
if (length(surf$at) != 21L) {
  fail(sprintf(
    "extract_cb_surface returned %d exposure points, expected 21",
    length(surf$at)
  ))
}

# 3) make_run_file distinct paths for single_lag 1 vs 2 ------------------------
mk_dlnm <- function(L) {
  list(
    enable = TRUE,
    defaults = list(
      max_lag = 11L,
      argvar = argvar,
      re_adm1_slope = TRUE,
      keep_linear_term = FALSE,
      single_lag = L
    ),
    variables = list(enso = list())
  )
}

# occurrence
occ1 <- make_options_from_config(
  list(),
  cli_argv = character(),
  model = "occurrence"
)
occ1$occurrence$backend <- "brms_cb"
occ2 <- occ1
occ1$occurrence$dlnm <- mk_dlnm(1L)
occ2$occurrence$dlnm <- mk_dlnm(2L)
if (
  basename(make_run_file(occ1, "fit")) == basename(make_run_file(occ2, "fit"))
) {
  fail("occurrence single_lag 1 vs 2 collide on fit basename")
}

# size
siz1 <- make_options_from_config(list(), cli_argv = character(), model = "size")
siz1$size$backend <- "brms_cb"
siz2 <- siz1
siz1$size$dlnm <- mk_dlnm(1L)
siz2$size$dlnm <- mk_dlnm(2L)
if (
  basename(make_run_file(siz1, "size_fit")) ==
    basename(make_run_file(siz2, "size_fit"))
) {
  fail("size single_lag 1 vs 2 collide on size_fit basename")
}

# 3b) ENSO leads: lag -3 and lag 3 never share a name -------------------------
if (
  !identical(
    enso_lag_col_name(c(-3L, 0L, 3L)),
    c("enso_lead3", "enso_lag0", "enso_lag3")
  )
) {
  fail("enso_lag_col_name(): -3/0/3 -> enso_lead3/enso_lag0/enso_lag3")
}
if (
  !identical(sl_token(c(-3L, 0L, 3L)), c("m3", "0", "3")) ||
    !identical(sl_tag(c(-3L, 3L)), c("slm3", "sl3"))
) {
  fail("sl_token()/sl_tag(): -3 -> m3/slm3, 3 -> 3/sl3")
}
occm3 <- occ1
occm3$occurrence$dlnm <- mk_dlnm(-3L)
occp3 <- occ1
occp3$occurrence$dlnm <- mk_dlnm(3L)
fm3 <- basename(make_run_file(occm3, "fit"))
fp3 <- basename(make_run_file(occp3, "fit"))
if (
  fm3 == fp3 ||
    !grepl("_slm3", fm3, fixed = TRUE) ||
    !grepl("_sl3", fp3, fixed = TRUE)
) {
  fail(
    "occurrence single_lag -3 vs 3: fit basenames must differ (_slm3 vs _sl3)"
  )
}
# the fit-time column lookup reads the lead column for a negative lag
dat_l <- tibble::tibble(enso_lag3 = rnorm(50), enso_lead3 = rnorm(50))
cb_l <- tryCatch(
  build_cb_block(
    dat_l,
    list(
      enable = TRUE,
      defaults = list(argvar = argvar, single_lag = -3L),
      variables = list(enso = list())
    )
  ),
  error = function(e) {
    fail(paste("build_cb_block(single_lag = -3) errored:", conditionMessage(e)))
    NULL
  }
)
if (
  !is.null(cb_l) &&
    (!identical(cb_l$cb_meta$enso$single_lag, -3L) ||
      !isTRUE(all.equal(cb_l$cb_meta$enso$range, range(dat_l$enso_lead3))))
) {
  fail("build_cb_block(single_lag = -3) must read enso_lead3")
}

# 4) Hierarchical CB slope plumbing (shared by 04b_cb + 04d_cb) --------------
# Both CB scripts now build crossbasis slope REs via nested_slope_re_terms()
# over re_group_by_level(), giving a region/country/ADM1 variance hierarchy
# rather than a single flat ADM1 term.
cb_cols <- sprintf("enso_cb%d", 1:2)

g_occ <- re_group_by_level("gid")
g_siz <- re_group_by_level("ADM2")
if (!identical(g_occ$adm1, "AFRO_region:country:ADM1")) {
  fail("re_group_by_level('gid') adm1 group is not fully nested")
}
if (!identical(g_siz$adm2, "AFRO_region:country:ADM1:ADM2")) {
  fail("re_group_by_level('ADM2') adm2 leaf col incorrect")
}

occ_terms <- nested_slope_re_terms(cb_cols, g_occ, "||")
# 2 columns × 3 levels = 6 terms
if (length(occ_terms) != 6L) {
  fail(sprintf("occ CB slopes: expected 6 terms, got %d", length(occ_terms)))
}
for (lvl_grp in c(
  "AFRO_region",
  "AFRO_region:country",
  "AFRO_region:country:ADM1"
)) {
  want <- sprintf("(enso_cb1 - 1 || %s)", lvl_grp)
  if (!want %in% occ_terms) {
    fail(sprintf("occ CB slopes: missing term %s", want))
  }
}

siz_terms <- nested_slope_re_terms(cb_cols, g_siz, "||")
if (length(siz_terms) != 6L) {
  fail(sprintf("size CB slopes: expected 6 terms, got %d", length(siz_terms)))
}
# occ and size share identical region/country/ADM1 group strings for slopes
if (!identical(occ_terms, siz_terms)) {
  fail(
    "occ and size CB slope terms diverge (should be identical at region/country/ADM1)"
  )
}

# 5) adm1_lineage_levels round-trip (figure <-> postprocess agreement) -------
lin <- adm1_lineage_levels("Eastern Africa_ETH_ETH_3")
if (
  !identical(lin$region, "Eastern Africa") ||
    !identical(lin$country, "Eastern Africa_ETH") ||
    !identical(lin$adm1, "Eastern Africa_ETH_ETH_3")
) {
  fail("adm1_lineage_levels did not recover region/country/adm1 level names")
}

message("PASS: test_cb_singlelag")
quit(status = 0)
