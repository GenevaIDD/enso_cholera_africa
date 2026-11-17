# test_lp_shapefiles_vintage.R — HARD: population.lp_shapefiles reads the
#   shared outbreak-data vintage.
#
# Usage:
#   Rscript analysis/tests/test_lp_shapefiles_vintage.R
#
# Details:
#   Quick-safe (no data, no fit). utils/paths.R::default_outbreak_data_dir()
#   (constant OUTBREAK_DATA_VINTAGE) is the one default database vintage,
#   and the LP geometries 02a/02b read must follow it. config.R may not
#   source paths.R (test_utils_module_purity.R), so the default is a literal
#   there and in config_defaults.yml; this test pins both to the paths.R
#   constant. Each module is loaded into its own environment after base.R,
#   so neither sees the other.

suppressPackageStartupMessages(library(yaml))

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
# The defaults are pinned with the outbreak-data root at its default
# (CHOLERA_OUTBREAK_DATA_DIR unset; test_data_vintage_default.R [1b] covers the override).
Sys.unsetenv("CHOLERA_OUTBREAK_DATA_DIR")

env_paths <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env_paths)
sys.source("analysis/utils/paths.R", env_paths)

env_cfg <- new.env(parent = globalenv())
sys.source("analysis/utils/base.R", env_cfg)
sys.source("analysis/utils/config.R", env_cfg)

want_dir <- env_paths$default_outbreak_data_dir()
r_lp <- env_cfg$get_default_options()$population$lp_shapefiles
yml_lp <- yaml::read_yaml(
  "analysis/config_defaults.yml"
)$population$lp_shapefiles

cat("# population.lp_shapefiles vintage\n\n")
cat(sprintf(
  "  paths.R vintage dir : %s\n  get_default_options : %s\n  config_defaults.yml : %s\n\n",
  want_dir,
  r_lp,
  yml_lp
))

expect(
  identical(basename(dirname(r_lp)), env_paths$OUTBREAK_DATA_VINTAGE),
  "get_default_options()$population$lp_shapefiles sits in the OUTBREAK_DATA_VINTAGE directory"
)
expect(
  identical(dirname(r_lp), path.expand(want_dir)), # dirname() expands `~`
  "... and that directory is default_outbreak_data_dir()"
)
expect(
  identical(basename(r_lp), "outbreak_shapefiles.rds"),
  "... and the file is outbreak_shapefiles.rds"
)
expect(
  identical(yml_lp, r_lp),
  "config_defaults.yml population.lp_shapefiles equals get_default_options()"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
