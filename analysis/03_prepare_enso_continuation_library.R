# 03_prepare_enso_continuation_library.R — build the seam-agnostic ENSO analog
# continuation library (Lenssen et al. 2024 piControl, 11 CMIP6 models).
#
# Read by utils/enso.R::loadEnsoContinuationLibrary() when a forward run sets
# simulator.enso.continuation.method: analog. Past the last CFSv2 target month,
# each replicate's ENSO is continued by a k-nearest-neighbour analog drawn from
# this library, matched on the replicate's own seam state (calendar month, ONI,
# 3-month tendency ONI(m) - ONI(m-3), ONI(m-6)) and joined with a fading offset
# (enso.R::continueEnsoAnalog()). Validated on 468 observed seams
# (1982-2020: best fair CRPS at every lead window vs neutral / ar1 / climatology).
#
# Unlike generated_data/lenssen_nino34_analogs.rds (03_prepare_lenssen_analogs.R,
# pre-matched to ONE start month and ONE ONI value for the counterfactual), this
# library stores whole series and per-month features, so any seam can be matched
# at run time. Rebuild: never needed for piControl (fixed archive); only when the
# rescale window or the horizon changes.
#
# Inputs:
#   <zenodo_dir>/extracted/<model>/step01_InitVerifInfo.Rda  (extracted by
#     03_prepare_lenssen_analogs.R, or here from <zenodo_dir>/<model>.zip);
#     `nino34` is the detrended 3-month-mean ONI-equivalent, `timeMap[, 2]` the
#     calendar month.
#   Observed ONI cache (fetchEnsoIndex("oni")): only for the rescale sd.
# Each model's series is rescaled to the observed ONI sd over --sd_window
# (03_prepare_lenssen_analogs.R uses the full record).
#
# Outputs:
#   generated_data/enso_continuation_library.qs        (paths.R::makeEnsoContinuationLibraryFile())
#   generated_data/enso_continuation_library_meta.yml
#
# --validate (no build; reads the library on disk and the ONI cache only):
# skill by lead of the analog continuation as deployed, on observed seams
# (default 1982-01 .. 2020-12, perfect seam, so optimistic), against the
# leave-out climatology and neutral references: fair CRPS, RMSE, spread, CRPSS
# vs climatology / neutral, mean over --validate_seeds with the seed range.
# The analog knobs (k, blend_months, min_sep, max_per_model, exclude_models,
# library) are those of analysis/config_defaults.yml's
# simulator.enso.continuation, overlaid by --deploy's deploy.enso_continuation.
# Logic: utils/enso_continuation_skill.R. 09
# reads the CSV for the continuation CRPSS of each group (S28).
#   generated_data/enso_continuation_skill.csv         (paths.R::makeEnsoContinuationSkillFile())
#   generated_data/enso_continuation_skill_meta.yml    (library digest, knobs, seeds, seams, leads)
#
# Usage:
#   Rscript analysis/03_prepare_enso_continuation_library.R
#   Rscript analysis/03_prepare_enso_continuation_library.R --horizon_months 60
#   Rscript analysis/03_prepare_enso_continuation_library.R --validate [--deploy analysis/configs/deploy/prod_run<N>.yml]

# Preamble ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(optparse)
})

source("analysis/utils/base.R")
source("analysis/utils/paths.R")
source("analysis/utils/enso.R")
source("analysis/utils/obs_sim_compare.R")          # crps_fair() (--validate)
source("analysis/utils/enso_continuation_skill.R")  # --validate

# Options -----------------------------------------------------------------

option_list <- list(
  make_option("--zenodo_dir", default = "data/lenssen/", type = "character",
              help = "Directory with the 11 model zips; extracted Rda under <dir>/extracted/ [default: %default]"),
  make_option("--refresh", default = FALSE, action = "store_true",
              help = "Re-extract the Rda files even if present [default: %default]"),
  make_option("--sd_window", default = "1991_2020", type = "character",
              help = "Observed ONI years (YYYY_YYYY) whose sd each model is rescaled to [default: %default]"),
  make_option("--horizon_months", default = 48L, type = "integer",
              help = "Months of future a library start must have to be a candidate [default: %default]"),
  make_option("--output_dir", default = "generated_data/", type = "character",
              help = "Directory for the library and its meta YAML (and the ONI cache) [default: %default]"),
  make_option("--validate", default = FALSE, action = "store_true",
              help = paste("Do not build: score the library on disk (skill by lead vs climatology / neutral)",
                           "and write enso_continuation_skill.csv + _meta.yml to --output_dir [default: %default]")),
  make_option("--deploy", default = "", type = "character",
              help = paste("--validate: deploy YAML whose deploy.enso_continuation overlays the defaults'",
                           "simulator.enso.continuation [default: none, the defaults only]")),
  make_option("--library", default = "", type = "character",
              help = "--validate: library to score [default: the continuation block's `library`]"),
  make_option("--validate_seeds", default = "20261018,20261019,20261020", type = "character",
              help = "--validate: comma-separated seeds; scores are seed means [default: %default]"),
  make_option("--validate_leads", default = 24L, type = "integer",
              help = "--validate: leads 1..N after the seam (<= the library's horizon) [default: %default]"),
  make_option("--validate_seams", default = "1982-01_2020-12", type = "character",
              help = "--validate: observed seam months YYYY-MM_YYYY-MM [default: %default]"),
  make_option("--validate_n_ens", default = 100L, type = "integer",
              help = "--validate: members per analog / climatology ensemble [default: %default]"),
  make_option("--validate_clim_from", default = "1950-01", type = "character",
              help = "--validate: first observed month of the climatology pool [default: %default]"),
  make_option("--validate_clim_exclude", default = 36L, type = "integer",
              help = "--validate: months either side of the seam left out of the climatology [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --validate: skill by lead of the deployed analog continuation ---------------

if (isTRUE(opt$validate)) {
  cont <- ensoContinuationSkillSettings("analysis/config_defaults.yml",
                                        if (nzchar(opt$deploy)) opt$deploy else NULL)
  lib_path <- if (nzchar(opt$library)) opt$library else cont$library
  lib <- loadEnsoContinuationLibrary(lib_path)
  digest_now <- ensoContinuationLibraryDigest(lib)
  digest_rec <- lib$meta$digest %||% digest_now
  if (!identical(digest_now, digest_rec))
    warning("library content digest ", digest_now, " differs from its recorded ", digest_rec,
            call. = FALSE)
  HV <- as.integer(opt$validate_leads)
  if (is.na(HV) || HV < 1L || HV > lib$meta$horizon_months)
    stop("--validate_leads must be in 1..", lib$meta$horizon_months, " (the library's horizon)")
  seeds <- as.integer(trimws(strsplit(opt$validate_seeds, ",")[[1L]]))
  if (!length(seeds) || anyNA(seeds)) stop("--validate_seeds must be comma-separated integers")
  sw <- strsplit(opt$validate_seams, "_")[[1L]]
  if (length(sw) != 2L) stop("--validate_seams must be YYYY-MM_YYYY-MM")
  seam_from <- as.Date(paste0(sw[1L], "-01")); seam_to <- as.Date(paste0(sw[2L], "-01"))
  clim_from <- as.Date(paste0(opt$validate_clim_from, "-01"))
  if (anyNA(c(seam_from, seam_to, clim_from))) stop("--validate_seams / --validate_clim_from: YYYY-MM")
  n_ens <- as.integer(opt$validate_n_ens)
  clim_excl <- as.integer(opt$validate_clim_exclude)

  message("-- ENSO continuation skill (--validate)")
  message("   library        : ", lib_path, " (digest ", digest_rec, ")")
  message("   settings from  : ", cont$source, " (configured method: ", cont$method, ")")
  if (!identical(cont$method, "analog"))
    message("   note           : the configured method is not analog; scoring its analog knobs anyway")
  message(sprintf("   analog knobs   : k %d, blend_months %s, min_sep %d, max_per_model %s, exclude %s",
                  cont$k, format(cont$blend_months), cont$min_sep, format(cont$max_per_model),
                  if (length(cont$exclude_models)) paste(cont$exclude_models, collapse = ",") else "none"))

  oni <- fetchEnsoIndex(source = "oni", cache_dir = opt$output_dir, refresh = FALSE)
  oni <- oni[order(oni$year, oni$month), ]
  obs   <- as.numeric(oni$enso_index)
  obs_d <- as.Date(sprintf("%d-%02d-01", oni$year, oni$month))
  obs_m <- as.integer(oni$month)
  seams <- ensoSkillSeams(obs_d, seam_from, seam_to, horizon = HV)
  if (!length(seams)) stop("--validate: no observed seam in ", opt$validate_seams, " with ", HV,
                           " months of observed future")
  message(sprintf("   seams          : %d (%s .. %s), leads 1..%d, %d member(s), seeds %s",
                  length(seams), format(obs_d[min(seams)], "%Y-%m"), format(obs_d[max(seams)], "%Y-%m"),
                  HV, n_ens, paste(seeds, collapse = ",")))

  methods <- list(
    neutral = ensoSkillMethodNeutral(),
    clim    = ensoSkillMethodClim(n_ens, pool_from = clim_from, exclude = clim_excl),
    analog  = ensoSkillMethodAnalog(lib, n_ens, k = cont$k, blend_months = cont$blend_months,
                                    min_sep = cont$min_sep, max_per_model = cont$max_per_model,
                                    exclude_models = cont$exclude_models))
  skill <- ensoContinuationSkill(obs, obs_d, obs_m, methods, seams, horizon = HV, seeds = seeds,
                                 references = c("clim", "neutral"))

  out_csv <- makeEnsoContinuationSkillFile(opt$output_dir)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(skill, out_csv, row.names = FALSE)
  git_head <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                       error = function(e) character(0), warning = function(w) character(0))
  meta <- list(
    created  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    script   = "analysis/03_prepare_enso_continuation_library.R --validate",
    git_head = if (length(git_head)) git_head[1L] else NA_character_,
    library  = list(path = lib_path, digest = digest_rec, digest_recomputed = digest_now,
                    horizon_months = lib$meta$horizon_months, models = as.list(names(lib$series))),
    continuation = list(
      settings_source = cont$source, configured_method = cont$method, scored_method = "analog",
      k = cont$k, blend_months = cont$blend_months, min_sep = cont$min_sep,
      max_per_model = if (is.finite(cont$max_per_model)) as.integer(cont$max_per_model) else "none",
      exclude_models = as.list(cont$exclude_models),
      features = list("v = ONI(s)", "tend = ONI(s) - ONI(s-3)", "phase = ONI(s-6)")),
    references = list(
      clim = sprintf(paste("same-calendar-month observed ONI from %s, %d month(s) either side of the",
                           "seam left out, %d draws per lead"),
                     format(clim_from, "%Y-%m"), clim_excl, n_ens),
      neutral = "ONI = 0"),
    seams = list(window = opt$validate_seams, n = length(seams),
                 first = format(obs_d[min(seams)], "%Y-%m"), last = format(obs_d[max(seams)], "%Y-%m"),
                 state = "observed ONI (perfect seam: optimistic vs a CFSv2-forecast seam)"),
    leads = HV,
    lead_definition = "lead k = k-th month after the seam (= continuation lead past the last CFSv2 month)",
    n_ens = n_ens, seeds = as.list(seeds),
    oni = list(cache_dir = opt$output_dir, last_month = format(max(obs_d), "%Y-%m")),
    scoring = paste("fair CRPS vs observed ONI per seam x lead, mean over seams; rmse of the ensemble",
                    "mean; spread = root mean ensemble variance; crpss_vs_<ref> = 1 - crps / crps_<ref>",
                    "of the seed-mean CRPS; _min / _max = range of the per-seed CRPSS"),
    columns = as.list(names(skill)))
  yaml::write_yaml(meta, sub("[.]csv$", "_meta.yml", out_csv))

  an <- skill[skill$method == "analog", ]
  message("   analog CRPSS vs climatology by lead [per-seed range], and vs neutral:")
  for (r in seq_len(nrow(an)))
    message(sprintf("     lead %2d: %6.3f  [%6.3f, %6.3f]   vs neutral %6.3f", an$lead[r],
                    an$crpss_vs_clim[r], an$crpss_vs_clim_min[r], an$crpss_vs_clim_max[r],
                    an$crpss_vs_neutral[r]))
  message("-- wrote ", out_csv, " (+ _meta.yml)")
  quit(save = "no", status = 0)
}
H <- as.integer(opt$horizon_months)
if (is.na(H) || H < 1L) stop("--horizon_months must be a positive integer")

model_names <- LENSSEN_MODELS   # utils/enso.R (shared with 03_prepare_lenssen_analogs.R)

message("-- ENSO continuation library (Lenssen piControl)")
message("   zenodo_dir     : ", opt$zenodo_dir)
message("   sd_window      : ", opt$sd_window)
message("   horizon_months : ", H)

# Observed ONI sd ----------------------------------------------------------

obs_oni <- fetchEnsoIndex(source = "oni", cache_dir = opt$output_dir, refresh = FALSE)
yr <- as.integer(strsplit(opt$sd_window, "_")[[1L]])
if (length(yr) != 2L || anyNA(yr)) stop("--sd_window must be YYYY_YYYY")
in_win <- obs_oni$year >= yr[1L] & obs_oni$year <= yr[2L] & is.finite(obs_oni$enso_index)
if (sum(in_win) < 120L) stop("fewer than 120 observed ONI months in --sd_window ", opt$sd_window)
sd_obs <- stats::sd(obs_oni$enso_index[in_win])
message(sprintf("   observed ONI sd (%s): %.4f", opt$sd_window, sd_obs))

# Extract: utils/enso.R::lenssenExtractRda() (shared with 03_prepare_lenssen_analogs.R)
# into <zenodo_dir>/extracted/<model>/.

# Build ---------------------------------------------------------------------

series <- list(); feats <- list(); scale <- list(); n_months <- list()
for (m in model_names) {
  e <- new.env()
  load(lenssenExtractRda(m, opt$zenodo_dir, file.path(opt$zenodo_dir, "extracted"),
                         isTRUE(opt$refresh)), envir = e)
  if (is.null(e$nino34) || is.null(e$timeMap))
    stop(m, ": step01_InitVerifInfo.Rda lacks nino34 / timeMap", call. = FALSE)
  keep <- is.finite(e$nino34)
  x    <- as.numeric(e$nino34[keep])
  mon  <- as.integer(e$timeMap[keep, 2L])
  if (any(!mon %in% 1:12)) stop(m, ": calendar months outside 1..12", call. = FALSE)
  k <- sd_obs / stats::sd(x)
  x <- x * k
  series[[m]]   <- x
  feats[[m]]    <- ensoContinuationFeatures(x, mon, model = m, horizon = H)
  scale[[m]]    <- round(k, 4)
  n_months[[m]] <- length(x)
  message(sprintf("   %-20s %6d months  scale %.3f  candidates %6d",
                  m, length(x), k, sum(feats[[m]]$ok)))
}
features <- dplyr::bind_rows(feats)
ok <- features$ok
feature_sd <- vapply(c("v", "tend", "phase"),
                     function(f) stats::sd(features[[f]][ok]), numeric(1))

lib <- list(
  series     = series,
  features   = features,
  feature_sd = feature_sd,
  meta = list(
    source         = "Lenssen et al. 2024 GRL piControl (Zenodo 10045687)",
    models         = as.list(model_names),
    n_months       = n_months,
    n_candidates   = sum(ok),
    candidates_by_month = as.list(stats::setNames(as.integer(table(factor(features$mon[ok], 1:12))),
                                                  month.abb)),
    rescale        = list(rule = "per-model sd to observed ONI sd", window = opt$sd_window,
                          obs_sd = round(sd_obs, 4), factors = scale),
    horizon_months = H,
    features       = list(v = "ONI(m)", tend = "ONI(m) - ONI(m-3)", phase = "ONI(m-6)"),
    feature_sd     = as.list(round(feature_sd, 4)),
    created        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
)
lib$meta$digest <- ensoContinuationLibraryDigest(lib)

# Write -----------------------------------------------------------------------

out <- makeEnsoContinuationLibraryFile(opt$output_dir)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
qs2::qs_save(lib, out)
yaml::write_yaml(lib$meta, sub("\\.qs$", "_meta.yml", out))
message(sprintf("-- wrote %s (%d candidates over %d models; digest %s)",
                out, sum(ok), length(series), lib$meta$digest))
