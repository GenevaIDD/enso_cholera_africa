# preprocess_paths.R — KEY=value paths of the 00-02 data leg, for bash.
#
# Not a pipeline step (no number) and not a module: a CLI tool like
# deploy_prod_run.R. The 00-02 stage drivers, the preprocess campaign and the
# preprocess rsync tools read their paths here so that bash never builds an
# artefact name and never carries the outbreak vintage
# (tests/test_data_vintage_default.R): utils/paths.R stays the one authority.
#
# Prints, one per line on stdout:
#   OUTBREAK_VINTAGE=<paths.R::OUTBREAK_DATA_VINTAGE>
#   VINTAGE_DIR=<defaultOutbreakDataDir(), ~ expanded>
#   SECURE_ROOT=<SECURE_OUTBREAK_DIR, no trailing slash>
# and with --deploy (the deploy YAML's deploy.episodes block, which selects the
# national-zero months 01_preprocess_national_zero_months.R must build):
#   NZ_RULE=<national_zero_rule, "" = legacy>
#   NZ_THRESHOLD=<national_zero_threshold, 0 = exact zeros>
# plus what stages/03_prepare.sh passes to the 03 steps:
#   EPISODE_KEY=<episodes.key>
#   AT_RISK_START=<episodes.at_risk_start>   ADM1_OVERLAP=<episodes.adm1_overlap>
#   ENSO_CONT_ARGS=<03_prepare_sim_panel_projection.R flags of deploy.enso_continuation>
#   CFSV2_ARGS=<its CFSv2 correction flags, from the sim.forecast overlay's
#               simulator.enso.cfsv2_calibration[_target]; "" = the script's default>
# and with --outputs, the files the 01-02 steps write (to copy back from the
# cluster):
#   SEC=<path>   under SECURE_ROOT (outbreak-identifying); SECDIR= a folder
#   GEN=<path>   under the repo root (generated_data/); GENDIR= a folder
# The population artefacts follow -C <config> (its population: block); the
# national-zero months follow --deploy. The 02 crosswalk / climate outputs are
# the fixed generated_data/ names those scripts write.
#
# Usage:
#   Rscript analysis/utils/preprocess_paths.R
#   Rscript analysis/utils/preprocess_paths.R --deploy analysis/configs/deploy/prod_run<N>.yml
#   Rscript analysis/utils/preprocess_paths.R -C <config.yml> --deploy <D.yml> --outputs

suppressPackageStartupMessages(library(optparse))
source("analysis/utils/config.R")
source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/provenance.R")   # prepPolicyStampPath()

opt <- make_options_from_config(option_list = list(
  make_option("--deploy", default = NULL, type = "character",
              help = "deploy YAML whose deploy.episodes block names the national-zero months"),
  make_option("--outputs", action = "store_true", default = FALSE,
              help = "also list the files the 01-02 steps write (SEC= / GEN= lines)")
), model = NULL)

# `deploy` is also a config block that make_options_from_config() fills with its
# defaults (a list), so only a character value is a --deploy path.
deploy_yaml <- if (is.character(opt[["deploy"]])) opt[["deploy"]] else NULL
# The occurrence-panel spec of the deploy's episodes block, every absent key at
# its default (paths.R).
nz_spec <- NULL
dep <- NULL
if (!is.null(deploy_yaml)) {
  if (!file.exists(deploy_yaml)) stop("deploy YAML not found: ", deploy_yaml)
  dep <- yaml::read_yaml(deploy_yaml)[["deploy"]]
  ep <- dep[["episodes"]]
  if (is.null(ep)) stop(deploy_yaml, " has no deploy.episodes block")
  nz_spec <- .occPanelSpecFromEpisodes(ep)
}

# deploy.enso_continuation -> the projection panel's --enso_continuation* flags
# ("" when the block is absent: the script's defaults).
.ensoContArgs <- function(ec) {
  if (is.null(ec)) return("")
  fl <- c(method = "--enso_continuation", library = "--enso_continuation_library",
          k = "--enso_continuation_k", blend_months = "--enso_continuation_blend_months",
          min_sep = "--enso_continuation_min_sep", max_per_model = "--enso_continuation_max_per_model",
          exclude_models = "--enso_continuation_exclude_models", seam = "--enso_continuation_seam")
  out <- character(0)
  for (k in names(fl)) {
    v <- ec[[k]]
    if (is.null(v)) next
    v <- if (is.numeric(v) && is.infinite(v)) "Inf" else paste(unlist(v), collapse = ",")
    if (!nzchar(v)) next
    out <- c(out, fl[[k]], v)
  }
  paste(out, collapse = " ")
}
# sim.forecast overlay -> the projection panel's CFSv2 correction flags.
.cfsv2Args <- function(fc) {
  cal <- fc[["simulator.enso.cfsv2_calibration"]]
  if (is.null(cal) || !nzchar(cal)) return("")
  tg <- fc[["simulator.enso.cfsv2_calibration_target"]]   # absent = the script's default
  paste(c("--cfsv2_correction calibration --cfsv2_calibration", cal,
          if (!is.null(tg)) c("--cfsv2_calibration_target", tg)), collapse = " ")
}

kv <- function(k, v) cat(sprintf("%s=%s\n", k, v))
kv("OUTBREAK_VINTAGE", OUTBREAK_DATA_VINTAGE)
kv("VINTAGE_DIR", path.expand(defaultOutbreakDataDir()))
kv("SECURE_ROOT", sub("/+$", "", SECURE_OUTBREAK_DIR))

if (!is.null(nz_spec)) {
  if (!isTRUE(nz_spec$national_zero_fill))
    message("note: ", deploy_yaml, " has national_zero_fill off; the months are built anyway")
  kv("NZ_RULE", nz_spec$national_zero_rule)
  kv("NZ_THRESHOLD", nz_spec$national_zero_threshold)
  kv("EPISODE_KEY", nz_spec$key)
  kv("AT_RISK_START", nz_spec$at_risk_start)
  kv("ADM1_OVERLAP", nz_spec$adm1_overlap)
  kv("ENSO_CONT_ARGS", .ensoContArgs(dep[["enso_continuation"]]))
  kv("CFSV2_ARGS", .cfsv2Args(dep[["sim"]][["forecast"]]))
}

if (isTRUE(opt[["outputs"]])) {
  # 00 / 01 inside the vintage folder (time_series_*.rds, outbreak_shapefiles.*)
  kv("SECDIR", path.expand(defaultOutbreakDataDir()))
  # 01 / 01b / 01_nz at the secure root
  mon <- makeMonthlyCasesFile()
  for (f in c(makeOutbreakMetadataFile(), makeAggTsFile(), makeAggAdmTsFile(),
              mon, prepPolicyStampPath(mon)))
    kv("SEC", f)
  if (!is.null(nz_spec)) {
    nz <- makeNationalZeroMonthsFile(nz_spec$national_zero_threshold, nz_spec$national_zero_rule)
    for (f in c(nz, prepPolicyStampPath(nz), nationalZeroAuditPath(nz),
                nationalZeroConsistencyPath(nz_spec$national_zero_threshold)))
      kv("SEC", f)
  }
  # 02 population (02c), keyed by the config's population block
  for (k in c("pop_crosswalk", "pop_gid", "pop_wpp")) kv("GEN", makeRunFile(opt, kind = k))
  kv("GENDIR", makeRunFile(opt, kind = "pop_qc"))
  # 02 crosswalk / climate (fixed names in 02_preprocess_crosswalk.R / 02_preprocess_climate.R)
  for (f in c("generated_data/matched_ob_gids_crosswalk.rds", "generated_data/u_gids_crosswalk.rds",
              "generated_data/matched_ob_gids.rds", "generated_data/qc_climate_gids_missing.csv",
              "generated_data/climate_stats.sqlite", "generated_data/climate_extremes.sqlite",
              "generated_data/climate_db_provenance.yml"))
    kv("GEN", f)
  # per-gid climate caches (new gids of the vintage), and the step counts
  kv("GENDIR", dirname(makeClimStatsFile("x", "y")))
  kv("GENDIR", "generated_data/prep_counts")
}
