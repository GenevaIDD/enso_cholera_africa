# 03_prepare_rain_adm1.R — the ADM1 x month rainfall series that 09's
# rain_enso step (Figure 1c: rainfall anomaly-ONI correlation tiles and the
# rainy-season outlines) reads.
#
# climate_stats.sqlite (02_preprocess_climate.R) is 8.6 GB and exists only on
# the local machine; this step reduces it to the one table Figure 1d needs, a
# few MB, which is copied to the cluster so 09 can run there. Per ADM1: its own
# row of climate_stats where it has one, else the mean of its ADM2 rows
# (utils/rain_enso.R::adm1MonthlyPrecip()).
#
# Inputs:  generated_data/climate_stats.sqlite (table climate_stats: gid, month, precip)
# Outputs: generated_data/adm1_monthly_precip.qs  (adm1, month, precip, src, iso3, year, mo)
#          generated_data/adm1_monthly_precip_meta.yml  (source size / mtime, row and unit counts)
#
# Usage:
#   Rscript analysis/03_prepare_rain_adm1.R
#   Rscript analysis/03_prepare_rain_adm1.R --db generated_data/climate_stats.sqlite --outdir generated_data/

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})
for (m in c("base", "gid", "paths", "clim_db", "rain_enso"))
  source(sprintf("analysis/utils/%s.R", m))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--db", default = "generated_data/climate_stats.sqlite",
              help = "climate stats database [default %default]"),
  make_option("--outdir", default = "generated_data/",
              help = "output root [default %default]"))))

if (!file.exists(opt$db) || file.size(opt$db) == 0)
  stop("[03_prepare_rain_adm1] ", opt$db, " is missing or empty -- run 02_preprocess_climate.R ",
       "(locally: the database is not on the cluster)", call. = FALSE)

con <- getStatsDB(opt$db)
pr  <- DBI::dbGetQuery(con, "SELECT gid, month, precip FROM climate_stats")
DBI::dbDisconnect(con)
message(sprintf("[03_prepare_rain_adm1] %d gid-month rows read from %s", nrow(pr), opt$db))

adm1_m <- adm1MonthlyPrecip(pr)
src_tab <- adm1_m |> distinct(adm1, src) |> count(src)
message(sprintf("[03_prepare_rain_adm1] ADM1 units: %d (%s); months %s .. %s",
                n_distinct(adm1_m$adm1),
                paste(sprintf("%s=%d", src_tab$src, src_tab$n), collapse = ", "),
                min(adm1_m$month), max(adm1_m$month)))

out <- makeAdm1PrecipFile(opt$outdir)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
qs2::qs_save(adm1_m, out)
yaml::write_yaml(list(
  generator = "analysis/03_prepare_rain_adm1.R",
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  source = list(path = opt$db, size = format(file.size(opt$db), scientific = FALSE),
                mtime = format(file.mtime(opt$db), "%Y-%m-%d %H:%M:%S")),
  rows = nrow(adm1_m), n_adm1 = n_distinct(adm1_m$adm1),
  n_adm1_from_adm2 = sum(src_tab$n[src_tab$src == "adm2_mean"]),
  months = c(format(min(adm1_m$month)), format(max(adm1_m$month)))),
  sub("\\.qs$", "_meta.yml", out))
message("[03_prepare_rain_adm1] wrote ", out)
