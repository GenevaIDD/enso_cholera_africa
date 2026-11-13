# Project .Rprofile: opt-in renv activation (renv.lock, renv/activate.R).
#
# By default this file changes nothing: the pipeline runs on whatever library
# R finds, exactly as without this file. renv is activated only when BOTH
#   - the environment variable CHOLERA_USE_RENV is 1 / true / yes, and
#   - the renv package is installed (requireNamespace("renv")),
# e.g.  CHOLERA_USE_RENV=1 Rscript analysis/tests/run_all.R --quick
# First use in a fresh checkout: CHOLERA_USE_RENV=1 Rscript -e 'renv::restore()'.
#
# R reads only one user profile: a project .Rprofile shadows ~/.Rprofile, so the
# user's own profile (R_PROFILE_USER is not set when this file is read) is
# sourced first to keep the default behaviour identical.
local({
  home_prof <- path.expand("~/.Rprofile")
  if (file.exists(home_prof) &&
      normalizePath(home_prof, mustWork = FALSE) != normalizePath(".Rprofile", mustWork = FALSE))
    sys.source(home_prof, envir = globalenv())

  use_renv <- tolower(Sys.getenv("CHOLERA_USE_RENV", "")) %in% c("1", "true", "yes")
  if (use_renv) {
    if (file.exists("renv/activate.R") && requireNamespace("renv", quietly = TRUE)) {
      source("renv/activate.R")
    } else {
      message("CHOLERA_USE_RENV is set but renv or renv/activate.R is missing; ",
              "running on the default library.")
    }
  }
})
