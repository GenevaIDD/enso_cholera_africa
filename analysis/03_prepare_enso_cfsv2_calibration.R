# 03_prepare_enso_cfsv2_calibration.R — fit the hindcast calibration of the
#   CFSv2 Niño 3.4 forecast (replica of CPC's PDF + spread correction).
#
# Usage:
#   Rscript analysis/03_prepare_enso_cfsv2_calibration.R
#   Rscript analysis/03_prepare_enso_cfsv2_calibration.R --window 1991,2020 \
#     --refresh_obs
#
# Details:
#   The forecast arm's factual ENSO input is the raw CPC CFSv2 ensemble; this
#   rebuilds CPC's correction from the NMME CFSv2 hindcast archive, for the
#   targets "ersst" (the ONI scale the models were fit on) and "oisst"
#   (CPC's own, to validate the replica). Inputs are fetched and cached by
#   utils/enso.R. Writes enso_cfsv2_calibration.rds (read by the sampler via
#   simulator.enso.cfsv2_calibration_target), its _meta.yml and the
#   replica-vs-CPC validation CSV.
#
# See: apply_cfsv2_calibration() in analysis/utils/enso.R

# 1. Preamble ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(optparse)
})

source("analysis/utils/base.R")
source("analysis/utils/enso.R")

# 2. Options -------------------------------------------------------------------

option_list <- list(
  make_option(
    "--window",
    default = "1991,2020",
    type = "character",
    help = "Training start years first,last (CPC: 1991,2020) [default: %default]"
  ),
  make_option(
    "--targets",
    default = "ersst,oisst",
    type = "character",
    help = "Observational targets to fit (ersst, oisst) [default: %default]"
  ),
  make_option(
    "--cache_dir",
    default = "generated_data/",
    type = "character",
    help = "Directory of the ENSO caches [default: %default]"
  ),
  make_option(
    "--cpc_dir",
    default = "generated_data/enso_cfsv2_hindcast/cpc_ftp",
    type = "character",
    help = "Archived CPC PDFC/PDFSPRDC files for the validation table [default: %default]"
  ),
  make_option(
    "--refresh_hindcast",
    default = FALSE,
    action = "store_true",
    help = "Re-download the NMME hindcast (IRIDL; ~1.5-3 min) [default: %default]"
  ),
  make_option(
    "--refresh_obs",
    default = FALSE,
    action = "store_true",
    help = "Re-download the ERSST / OISST Niño 3.4 caches [default: %default]"
  ),
  make_option(
    "--output",
    default = "generated_data/enso_cfsv2_calibration.rds",
    type = "character",
    help = "Calibration RDS [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
window <- as.integer(strsplit(opt$window, ",")[[1]])
targets <- trimws(strsplit(opt$targets, ",")[[1]])
stopifnot(length(window) == 2L, all(targets %in% c("ersst", "oisst")))

# 3. Data ----------------------------------------------------------------------
# Inputs (fetched + cached by utils/enso.R):
#   generated_data/enso_cfsv2_hindcast/nmme_cfsv2_pentad_full_nino34.nc
#     (IRIDL, archived)
#   generated_data/enso_nino34_cache.rds (ERSSTv6 Niño 3.4, CPC data/indices)
#   generated_data/enso_nino34_oisst_cache.rds  (OISST v2.1 Niño 3.4)
#   generated_data/enso_cfsv2_hindcast/cpc_ftp/ (optional: CPC PDFC/PDFSPRDC
#     files for the validation table)
# CPC publishes its PDF + spread-corrected forecast only as images (corrected
# data files stop in Jun 2026), hence the replica.

hind <- fetch_cfsv2_hindcast_nino34(
  cache_dir = opt$cache_dir,
  refresh = opt$refresh_hindcast
)
obs_src <- c(ersst = "nino34", oisst = "nino34_oisst")
obs <- lapply(targets, function(tg) {
  fetch_enso_index(
    obs_src[[tg]],
    cache_dir = opt$cache_dir,
    refresh = opt$refresh_obs
  )
})
names(obs) <- targets

# 4. Fit -----------------------------------------------------------------------
# CPC's documented method (products/CFSv2/cfsv2fcst/CFSv2SST8210.html):
#   1. per (start month, lead), each member anomaly is mapped to its
#      percentile among the 1991-2020 hindcast + real-time members and
#      replaced by the observed anomaly at that percentile (capped at the
#      reference range);
#   2. the ensemble spread is set to the RMSE of the corrected ensemble mean
#      over the same period (applied to the raw member deviations, see
#      `apply_cfsv2_calibration()` for why).
# Fitted for two observational targets: "oisst" (CPC's own; used to validate
# the replica against CPC's files) and "ersst" (CPC's current ERSSTv6
# detrended Niño 3.4, the series published ONI is built from, i.e. the scale
# the disease models were fit on). The sampler picks one via
# `simulator.enso.cfsv2_calibration_target`.

t0 <- Sys.time()
cal <- fit_cfsv2_calibration(hind, obs, window = window)
message(sprintf(
  "[03_prepare_enso_cfsv2_calibration] fitted %d start months x %d leads x %d targets in %.1f s",
  length(unique(cal$clim$sm)),
  length(unique(cal$clim$lead)),
  length(targets),
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))

# 5. Validation: replica vs CPC's own corrected files --------------------------
#
# For each archived CPC issue, the NMME real-time members of the same start are
# corrected with the "oisst" target (CPC's) and compared, per lead, with CPC's
# PDF-corrected ensemble mean and its PDF + spread-corrected member sd.

read_cpc_gr <- function(f) {
  ctl <- readLines(sub("\\.gr$", ".ctl", f))
  nz <- as.integer(sub(
    "zdef\\s+(\\d+).*",
    "\\1",
    grep("^zdef", ctl, value = TRUE)
  ))
  td <- strsplit(trimws(grep("^tdef", ctl, value = TRUE)), "\\s+")[[1]]
  nt <- as.integer(td[2])
  t0 <- as.Date(paste0("01", td[4]), "%d%b%Y")
  v <- readBin(f, "numeric", size = 4, n = nz * nt, endian = "big")
  v[v == -9999] <- NA
  a <- matrix(v, nrow = nz, ncol = nt)
  tibble(
    target_date = seq(t0, by = "month", length.out = nt),
    ens_mean = a[124, ],
    member_sd = apply(a[1:120, , drop = FALSE], 2, stats::sd)
  )
}

validation <- NULL
if ("oisst" %in% targets && dir.exists(opt$cpc_dir)) {
  issues <- sub(
    "CFSv2nino34Mon(\\d{6})PDFC\\.gr",
    "\\1",
    list.files(opt$cpc_dir, "^CFSv2nino34Mon\\d{6}PDFC\\.gr$")
  )
  validation <- bind_rows(lapply(issues, function(iss) {
    init <- as.Date(paste0(iss, "01"), "%Y%m%d")
    nm <- hind %>% filter(init_date == init)
    if (!nrow(nm)) {
      return(NULL)
    }
    W <- tidyr::pivot_wider(
      nm %>% select(target_date, member, sst),
      names_from = member,
      values_from = sst
    )
    clim <- cal$clim %>% filter(sm == month(init)) %>% arrange(lead)
    X <- as.matrix(W[, -1]) -
      clim$clim[match(nm$lead[match(W$target_date, nm$target_date)], clim$lead)]
    rownames(X) <- as.character(W$target_date)
    Xc <- apply_cfsv2_calibration(X, cal, target = "oisst")
    pc <- read_cpc_gr(file.path(
      opt$cpc_dir,
      sprintf("CFSv2nino34Mon%sPDFC.gr", iss)
    ))
    psc <- read_cpc_gr(file.path(
      opt$cpc_dir,
      sprintf("CFSv2nino34Mon%sPDFSPRDC.gr", iss)
    ))
    tibble(
      issue = iss,
      target_date = as.Date(rownames(Xc)),
      lead = as.integer(round(
        as.numeric(difftime(target_date, init, units = "days")) / 30.4
      )),
      nmme_raw_mean = rowMeans(X, na.rm = TRUE),
      replica_mean = rowMeans(Xc, na.rm = TRUE),
      replica_sd = apply(Xc, 1, stats::sd, na.rm = TRUE)
    ) %>%
      left_join(
        pc %>% select(target_date, cpc_pdfc_mean = ens_mean),
        by = "target_date"
      ) %>%
      left_join(
        psc %>% select(target_date, cpc_pdfsprdc_sd = member_sd),
        by = "target_date"
      )
  })) %>%
    filter(is.finite(cpc_pdfc_mean))
  if (nrow(validation)) {
    message(
      "[03_prepare_enso_cfsv2_calibration] replica vs CPC PDFC (",
      length(issues),
      " issues):"
    )
    print(
      validation %>%
        group_by(lead) %>%
        summarise(
          n = n(),
          rmse_mean = round(sqrt(mean((replica_mean - cpc_pdfc_mean)^2)), 3),
          cpc_correction = round(
            sqrt(mean((nmme_raw_mean - cpc_pdfc_mean)^2)),
            3
          ),
          sd_ratio = round(mean(replica_sd / cpc_pdfsprdc_sd), 2)
        ) %>%
        as.data.frame()
    )
  }
}

# 6. Write ---------------------------------------------------------------------

hind_meta <- file.path(
  opt$cache_dir,
  "enso_cfsv2_hindcast_nino34_download_meta.yml"
)
cal$meta$inputs <- list(
  hindcast_meta = if (file.exists(hind_meta)) {
    yaml::read_yaml(hind_meta)
  } else {
    NULL
  },
  obs = lapply(targets, function(tg) {
    f <- file.path(
      opt$cache_dir,
      sprintf("enso_%s_download_meta.yml", obs_src[[tg]])
    )
    if (file.exists(f)) yaml::read_yaml(f) else NULL
  }) |>
    stats::setNames(targets)
)
cal$meta$created <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
cal$validation <- validation

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
saveRDS(cal, opt$output)
spread_tbl <- bind_rows(lapply(names(cal$targets), function(tg) {
  cal$targets[[tg]]$spread %>% mutate(target = tg)
}))
yaml::write_yaml(
  c(
    cal$meta[c("method", "window", "n_starts", "leads", "targets", "created")],
    list(
      spread_factor_range = lapply(
        split(spread_tbl$g, spread_tbl$target),
        function(g) round(range(g), 3)
      ),
      validation_rmse_mean_by_lead = if (
        !is.null(validation) && nrow(validation)
      ) {
        as.list(round(
          tapply(
            (validation$replica_mean - validation$cpc_pdfc_mean)^2,
            validation$lead,
            function(z) sqrt(mean(z))
          ),
          3
        ))
      } else {
        NULL
      }
    )
  ),
  sub("\\.rds$", "_meta.yml", opt$output)
)
if (!is.null(validation) && nrow(validation)) {
  utils::write.csv(
    validation,
    sub("\\.rds$", "_validation.csv", opt$output),
    row.names = FALSE
  )
}
message("[03_prepare_enso_cfsv2_calibration] wrote ", opt$output)
