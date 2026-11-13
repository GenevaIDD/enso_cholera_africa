# test_synthetic_world.R -- the synthetic world feeds the real 00 -> 03 steps.
#
# SOFT gate (see analysis/tests/run_all.R; skipped under --quick): builds a
# sandbox in a temp dir and runs the data phase of
# analysis/bash/tools/run_synthetic_pipeline.sh (00_ingest_synthetic.R, then the
# real 00 ingest, 01, 02 population / crosswalk / climate and every 03 panel;
# ~2 min, no fit, no network), then checks with the pipeline's own readers:
#   [1] no real country: every ISO3 in the database, GADM, WHO and seasonality
#       files is a synthetic XA? code, and no real country name appears
#   [2] the national-zero months carry the r3 rule digest 03 checks
#   [3] the occurrence, duration and monthly panels and the simulator panel are
#       gid_season panels of one episode definition (simEpisodeCheck)
#   [4] the occurrence panel: both synthetic groups and all four countries,
#       no NA AFRO_region (the synthetic fill), the deployed contagion column
#       of the synthetic partition without NA, onsets present
#   [5] the simulator panel has no NA climate in the fit window
#   [6] nothing was downloaded (ENSO caches are the synthetic ones)
#
# Run: Rscript analysis/tests/test_synthetic_world.R   (from the repository root)

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
for (p in c("taxdat", "sf", "arrow", "terra", "giscoR", "OutbreakExtractR"))
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("SKIP: package %s not installed\n", p)); quit(status = 0)
  }

repo <- getwd()
sandbox <- file.path(tempdir(), "synthetic_world_test")
unlink(sandbox, recursive = TRUE)
cat("# synthetic world: data phase in", sandbox, "\n")
st <- system2("bash", c("analysis/bash/tools/run_synthetic_pipeline.sh", "--sandbox", sandbox,
                        "--until", "data", "--cores", "2"), stdout = TRUE, stderr = TRUE)
cat(paste(st, collapse = "\n"), "\n")
expect(is.null(attr(st, "status")), "run_synthetic_pipeline.sh --until data exits 0")
if (!is.null(attr(st, "status"))) quit(status = 1)

setwd(sandbox)
Sys.setenv(CHOLERA_OUTBREAK_DATA_DIR = file.path(sandbox, "secure", "cholera_outbreak_data"))
suppressPackageStartupMessages({ library(dplyr) })
for (m in c("base", "gid", "paths", "provenance", "national_zero"))
  source(file.path("analysis/utils", paste0(m, ".R")))

# [1] no real country --------------------------------------------------------
cat("\n[1] only synthetic countries\n")
reg <- syntheticCountries()
e <- new.env(); utils::data("gisco_countrycode", package = "giscoR", envir = e)
real_iso3 <- e$gisco_countrycode$ISO3_CODE
real_names <- tolower(e$gisco_countrycode$cldr.short.en)
db <- read.csv(file.path(defaultOutbreakDataDir(), "combined_outbreaks_cholera.csv"), stringsAsFactors = FALSE)
iso_db <- unique(db$country_iso3)
gadm1 <- sf::st_drop_geometry(sf::st_read("data/geodata/gadm_valid.gpkg", layer = "ADM_1", quiet = TRUE))
who <- readRDS("data/who_annual_report_africa/who_annual_report_africa.rds")
seas <- read.csv("data/best_seasonality_estimates.csv")
all_iso <- unique(c(iso_db, gadm1$GID_0, who$iso3, seas$country))
expect(all(all_iso %in% reg$iso3) && !any(all_iso %in% real_iso3),
       sprintf("every ISO3 is synthetic (%s)", paste(all_iso, collapse = ", ")))
tokens <- unique(tolower(unlist(strsplit(c(db$location, gadm1$NAME_1, gadm1$COUNTRY), "[: ]+"))))
expect(!any(tokens %in% real_names[nchar(real_names) > 3]),
       "no real country name in locations or GADM names")

# [2] national-zero stamp ----------------------------------------------------
cat("\n[2] national-zero months\n")
nz <- makeNationalZeroMonthsFile(60L, "r3")
stamp <- readPrepPolicyStamp(nz)
expect(file.exists(nz) && identical(stamp$rule, "r3") &&
         identical(stamp$rule_params_digest, nationalZeroRuleDigest(nationalZeroRuleParams("r3"))),
       "national_zero_months_t60_r3 carries the r3 rule digest")

# [3] one episode definition -------------------------------------------------
cat("\n[3] episode stamps\n")
panel <- function(rx) { f <- list.files("generated_data", rx, full.names = TRUE); f[!grepl("policy|meta|csv$", f)][1] }
preps <- list(occ = panel("^outbreaks_survival_data_base-1990_2019_ep-gs_nz60-r3_ae_a1c\\.rds$"),
              dur = panel("^outbreaks_duration_data_base-1990_2019_ep-gs.*\\.rds$"),
              mon = panel("^outbreaks_monthly_data_base-1990_2019_ep-gs.*\\.rds$"),
              sim = panel("^sim_climate_panel_base-1990_2019_ep-gs\\.rds$"))
expect(all(vapply(preps, function(p) !is.na(p) && file.exists(p), TRUE)),
       sprintf("the four gid_season panels exist (%s)", paste(basename(unlist(preps)), collapse = ", ")))
ok <- tryCatch({ simEpisodeCheck(preps = preps, episode_key = "gid_season"); TRUE },
               error = function(e) { cat("  ", conditionMessage(e), "\n"); FALSE })
expect(ok, "simEpisodeCheck(): the panels share one gid_season episode definition")

# [4] occurrence panel -------------------------------------------------------
cat("\n[4] occurrence panel\n")
occ <- readRDS(preps$occ)
expect(setequal(unique(as.character(occ$country)), reg$iso3), "all four synthetic countries")
expect(!anyNA(occ$AFRO_region), "AFRO_region filled for every row")
nbr <- nbrCovariateName(500L, "rownorm", "registry", "hmm_clusters_synth")
expect(nbr %in% names(occ) && !anyNA(occ[[nbr]]), sprintf("%s present, no NA", nbr))
expect(sum(occ$outbreak, na.rm = TRUE) > 50, sprintf("onsets in the panel (%d)", as.integer(sum(occ$outbreak, na.rm = TRUE))))

# [5] simulator panel climate ------------------------------------------------
cat("\n[5] simulator panel climate\n")
sp <- readRDS(preps$sim)
win <- sp$t_start >= as.Date("2010-01-01") & sp$t_start <= as.Date("2024-12-01")
cl <- intersect(c("mean_precip", "dev_precip", "mean_temp_daymean", "dev_temp_daymean"), names(sp))
expect(length(cl) == 4L && !anyNA(sp[win, cl]), "no NA climate in 2010-2024")

# [6] no download ------------------------------------------------------------
cat("\n[6] offline\n")
meta <- yaml::read_yaml("generated_data/enso_oni_download_meta.yml")
expect(grepl("synthetic", meta$url %||% ""), "the ONI cache is still the synthetic one")

setwd(repo)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
