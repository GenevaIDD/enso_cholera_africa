# test_cb_block_extraction.R — HARD: build_cb_block() (analysis/utils/dlnm.R)
#   reproduces the crossbasis-construction loop of the brms_cb fit workers.
#
# Usage:
#   Rscript analysis/tests/test_cb_block_extraction.R
#
# Details:
#   Asserts byte-for-byte equality (bound columns, cb_cols, cb_meta) with a
#   frozen reference loop across a single-lag (onebasis) config, a
#   distributed-lag (crossbasis) config and a mixed config; the loop does not
#   depend on which track calls it. If this goes red, the extraction changed
#   behaviour: revert it.

suppressPackageStartupMessages({
  library(dplyr)
  library(dlnm)
})

source("analysis/utils/config.R")
source("analysis/utils/paths.R") # nbr_covariate_pattern()/nbr_covariate_name() -- call-time dependency of formula.R::build_formula_terms() and simulator_core.R::evaluate_lp() (utils/README.md dependency table)
source("analysis/utils/base.R")
source("analysis/utils/dlnm.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R") # expect/fail/expect_true/expect_equal

cat("# build_cb_block() extraction guard\n")

# --- Frozen pre-extraction reference loop ------------------------------------
# `.reference_cb_loop()` below is a frozen copy of the inline loop
# 04_fit_occ_brms_cb.R, 04_fit_siz_brms_cb.R and 04_fit_dur_brms_cb.R carried
# before build_cb_block(). The three configs stand in for "one occ, one siz,
# one dur config".
.reference_cb_loop <- function(dat, dlnm_cfg) {
  cb_cols <- character(0)
  cb_meta <- list()

  if (isTRUE(dlnm_cfg$enable) && length(dlnm_cfg$variables) > 0) {
    dlnm_defaults <- dlnm_cfg$defaults %||% list()

    for (v in names(dlnm_cfg$variables)) {
      var_spec <- dlnm_cfg$variables[[v]] %||% list()
      spec <- dlnm_defaults
      for (k in names(var_spec)) {
        spec[[k]] <- var_spec[[k]]
      }

      max_lag <- as.integer(spec$max_lag %||% 3L)
      argvar <- spec$argvar %||% list(fun = "bs", degree = 3L, df = 4L)
      arglag <- spec$arglag %||% list(fun = "ns", df = 3L)
      single_lag <- spec$single_lag

      if (!is.null(single_lag)) {
        L <- as.integer(single_lag)
        xcol <- if (identical(v, "enso")) {
          sprintf("enso_lag%d", L)
        } else if (L == 0L) {
          v
        } else {
          sprintf("%s_lag%d", v, L)
        }
        x <- as.numeric(dat[[xcol]])
        x[is.na(x)] <- 0

        ob <- do.call(dlnm::onebasis, c(list(x), argvar))
        cb_names <- sprintf("%s_cb%d", v, seq_len(ncol(ob)))
        colnames(ob) <- cb_names
        dat <- dplyr::bind_cols(
          dat,
          tibble::as_tibble(as.data.frame(unclass(ob)))
        )
        cb_cols <- c(cb_cols, cb_names)

        ob_attrs <- attributes(ob)
        ob_attrs[c("dim", "dimnames", "class")] <- NULL
        cb_meta[[v]] <- list(
          single_lag = L,
          onebasis_attrs = ob_attrs,
          argvar = argvar,
          range = range(x, na.rm = TRUE),
          df = ncol(ob),
          cb_col_names = cb_names
        )
      } else {
        Q <- make_lag_matrix(dat, v, max_lag)
        cb <- dlnm::crossbasis(
          Q,
          lag = c(0L, max_lag),
          argvar = argvar,
          arglag = arglag
        )
        cb_names <- sprintf("%s_cb%d", v, seq_len(ncol(cb)))
        colnames(cb) <- cb_names
        dat <- dplyr::bind_cols(dat, tibble::as_tibble(as.data.frame(cb)))
        cb_cols <- c(cb_cols, cb_names)

        meta_attrs <- attributes(cb)
        cb_meta[[v]] <- list(
          argvar = meta_attrs$argvar,
          arglag = meta_attrs$arglag,
          lag = meta_attrs$lag,
          df = meta_attrs$df,
          range = meta_attrs$range,
          cb_col_names = cb_names
        )
      }
    }
  }

  list(dat = dat, cb_cols = cb_cols, cb_meta = cb_meta)
}

# --- Synthetic data (enso_lag* + a climate-style var with lag columns) ------
set.seed(1)
n <- 300
base_dat <- data.frame(uid = rep(1:30, each = 10), t_start = rep(1:10, 30))
enso_lags <- as.data.frame(matrix(rnorm(n * 13), nrow = n))
names(enso_lags) <- sprintf("enso_lag%d", 0:12)
precip_lags <- as.data.frame(matrix(abs(rnorm(n * 4, 50, 10)), nrow = n))
names(precip_lags) <- c("mean_precip", sprintf("mean_precip_lag%d", 1:3))
dat <- dplyr::bind_cols(base_dat, enso_lags, precip_lags)

# --- Config 1: single-lag only (mimics a size single-lag ENSO config) -------
cfg_single <- list(
  enable = TRUE,
  defaults = list(argvar = list(fun = "bs", degree = 3L, df = 4L)),
  variables = list(enso = list(single_lag = 4L))
)

ref1 <- .reference_cb_loop(dat, cfg_single)
new1 <- build_cb_block(dat, cfg_single, " (size model)")
expect(identical(ref1$cb_cols, new1$cb_cols), "single-lag: cb_cols identical")
expect(identical(ref1$cb_meta, new1$cb_meta), "single-lag: cb_meta identical")
expect(
  identical(as.data.frame(ref1$dat), as.data.frame(new1$dat)),
  "single-lag: bound dat identical"
)

# --- Config 2: distributed-lag only (mimics an occurrence crossbasis config) -
cfg_dist <- list(
  enable = TRUE,
  defaults = list(
    max_lag = 3L,
    argvar = list(fun = "bs", degree = 3L, df = 4L),
    arglag = list(fun = "ns", df = 3L)
  ),
  variables = list(mean_precip = list())
)

ref2 <- .reference_cb_loop(dat, cfg_dist)
new2 <- build_cb_block(dat, cfg_dist, " (occurrence model)")
expect(
  identical(ref2$cb_cols, new2$cb_cols),
  "distributed-lag: cb_cols identical"
)
expect(
  identical(ref2$cb_meta, new2$cb_meta),
  "distributed-lag: cb_meta identical"
)
expect(
  identical(as.data.frame(ref2$dat), as.data.frame(new2$dat)),
  "distributed-lag: bound dat identical"
)

# --- Config 3: mixed single-lag + distributed-lag (mimics a duration config) -
cfg_mixed <- list(
  enable = TRUE,
  defaults = list(
    max_lag = 3L,
    argvar = list(fun = "bs", degree = 3L, df = 4L),
    arglag = list(fun = "ns", df = 3L)
  ),
  variables = list(
    enso = list(
      single_lag = 8L,
      argvar = list(fun = "bs", degree = 3L, df = 8L)
    ),
    mean_precip = list()
  )
)

ref3 <- .reference_cb_loop(dat, cfg_mixed)
new3 <- build_cb_block(dat, cfg_mixed, " (duration model)")
expect(identical(ref3$cb_cols, new3$cb_cols), "mixed: cb_cols identical")
expect(identical(ref3$cb_meta, new3$cb_meta), "mixed: cb_meta identical")
expect(
  identical(as.data.frame(ref3$dat), as.data.frame(new3$dat)),
  "mixed: bound dat identical"
)

# --- Config 4: disabled / no variables (message-only branches) --------------
new4 <- build_cb_block(dat, list(enable = FALSE), "")
expect(
  length(new4$cb_cols) == 0L && length(new4$cb_meta) == 0L,
  "disabled dlnm: empty cb_cols/cb_meta"
)
expect(
  identical(as.data.frame(new4$dat), as.data.frame(dat)),
  "disabled dlnm: dat passed through unchanged"
)

new5 <- build_cb_block(dat, list(enable = TRUE, variables = list()), "")
expect(
  length(new5$cb_cols) == 0L && length(new5$cb_meta) == 0L,
  "enabled but no variables: empty cb_cols/cb_meta"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
