# install_occ_brms.R — install the R packages the brms occurrence pipeline
#   needs (01 / 02 preprocessing, 03_prepare_occ.R, 04 fits, 05 postprocess).
#
# Usage:
#   Rscript analysis/utils/install_occ_brms.R   # from the project root

install_if_missing <- function(
  pkgs,
  repos = c(CRAN = "https://cloud.r-project.org")
) {
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) == 0L) {
    message("All CRAN packages already installed.")
    return(invisible(character(0)))
  }
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = repos)
  invisible(missing)
}

# --------------------------------------------------------------------------
# 1. cmdstanr (r-universe, not on CRAN)
# --------------------------------------------------------------------------

if (!"cmdstanr" %in% rownames(installed.packages())) {
  message("Installing cmdstanr from r-universe...")
  install.packages(
    "cmdstanr",
    repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
  )
} else {
  message("cmdstanr already installed.")
}

# --------------------------------------------------------------------------
# 2. CRAN packages
# --------------------------------------------------------------------------

cran_pkgs <- c(
  # tidyverse sub-packages (decomposed from library(tidyverse))
  "dplyr",
  "ggplot2",
  "glue",
  "lubridate",
  "purrr",
  "readr",
  "stringr",
  "tibble",
  "tidyr",
  # model fitting
  "brms",
  "splines2",
  # distributed-lag non-linear models (crossbasis predictors)
  "dlnm",
  # postprocessing
  "loo",
  "posterior",
  # data / DB
  "DBI",
  "RSQLite",
  "sf",
  # parallel prep
  "furrr",
  # CFSv2 ENSO ensemble NetCDF ingestion (needs system libnetcdf)
  "ncdf4",
  # utilities
  "janitor",
  "nanoparquet",
  "optparse",
  "survival",
  "yaml"
)

install_if_missing(cran_pkgs)

# --------------------------------------------------------------------------
# 3. CmdStan (C++ backend used by brms via cmdstanr)
# --------------------------------------------------------------------------

cmdstan_ok <- tryCatch(
  {
    nzchar(cmdstanr::cmdstan_path())
  },
  error = function(e) FALSE
)

if (!cmdstan_ok) {
  message("Installing CmdStan...")
  cmdstanr::install_cmdstan()
} else {
  message("CmdStan already installed at: ", cmdstanr::cmdstan_path())
}

# --------------------------------------------------------------------------
# 4. taxdat (HopkinsIDD private package — used in 01_preprocess_outbreaks.R)
# --------------------------------------------------------------------------

if (!"taxdat" %in% rownames(installed.packages())) {
  message(
    "\ntaxdat is a private HopkinsIDD package and cannot be installed from CRAN.\n",
    "Install with:\n",
    "  remotes::install_github('HopkinsIDD/taxdat')\n"
  )
} else {
  message("taxdat already installed.")
}

message("\nDone.")
