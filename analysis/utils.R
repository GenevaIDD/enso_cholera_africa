# utils.R — back-compat facade that sources the analysis/utils/ modules
#
# Provides:
#   every function of the modules listed in `mods` below
#
# Details:
#   Kept so existing scripts and tests can source("analysis/utils.R"); new
#   code sources only the modules it needs. Must NOT be sourced with
#   local=/envir= (nested source() calls would escape to globalenv()).
#   Modules outside `mods` (e.g. figures_io.R) are deliberately not loaded.
#
# See: analysis/utils/README.md

local({
  d <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) {
      d <- dirname(normalizePath(of))
      break
    }
  }
  if (is.null(d)) {
    d <- "analysis"
  }
  mods <- file.path(
    d,
    "utils",
    paste0(
      c(
        "base",
        "gid",
        "paths",
        "loo_lags",
        "clim_compile",
        "clim_db",
        "occ_data",
        "siz_data",
        "enso",
        "population",
        "composite_lp",
        "provenance",
        "anchor",
        "occ_post",
        "siz_post",
        "dlnm",
        "plots",
        "sim_io",
        "coverage"
      ),
      ".R"
    )
  )
  miss <- mods[!file.exists(mods)]
  if (length(miss)) {
    stop(
      "utils.R facade: missing modules: ",
      paste(basename(miss), collapse = ", ")
    )
  }
  for (m in mods) {
    source(m)
  } # local=FALSE -> globalenv, intentional
})
