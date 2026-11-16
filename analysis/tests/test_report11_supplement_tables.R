# test_report11_supplement_tables.R -- the phase-11 helpers of the supplement
# update, on synthetic frames.
#
# HARD gate (see analysis/tests/run_all.R). Pins:
#   [1] afro_subregion_of() / with_afro(): the fixed ISO3 -> WHO-AFRO subregion
#       map (taxdat convention), the EMRO flag of SDN / SOM / DJI, the printed
#       column after the country column and the raw CSV columns
#   [2] report_data_id() / report_release_split(): D / T ids from CSV names,
#       copies into supp_data/ and github/, the manifest, an unmapped id listed
#       but not copied, an excluded slug withheld, a bad release value refused;
#       all-empty columns left out of the released copy, every other cell copied as text
#   [6] D35b / D35c: the WHO-AFRO subregion display order, the
#       per-block Kruskal-Wallis p joined onto every state row, printed cells
#
# Synthetic: no artefact, no data/.
# Run: Rscript analysis/tests/test_report11_supplement_tables.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
source("analysis/utils/results_core.R")
source("analysis/utils/results_report.R")

PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
throws <- function(expr, pattern = NULL) {
  r <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(r) && (is.null(pattern) || grepl(pattern, r))
}

cat("# [1] WHO-AFRO subregion\n")
MOD30 <- c(BEN = "Western Africa", CIV = "Western Africa", GHA = "Western Africa", GIN = "Western Africa",
           GNB = "Western Africa", LBR = "Western Africa", MLI = "Western Africa", NER = "Western Africa",
           NGA = "Western Africa", SLE = "Western Africa", TGO = "Western Africa",
           AGO = "Central Africa", CAF = "Central Africa", CMR = "Central Africa", COD = "Central Africa",
           COG = "Central Africa", TCD = "Central Africa",
           BDI = "Eastern Africa", ETH = "Eastern Africa", KEN = "Eastern Africa", SDN = "Eastern Africa",
           SOM = "Eastern Africa", SSD = "Eastern Africa", TZA = "Eastern Africa", UGA = "Eastern Africa",
           MOZ = "Southern Africa", MWI = "Southern Africa", NAM = "Southern Africa", ZMB = "Southern Africa",
           ZWE = "Southern Africa")
a <- afro_subregion_of(names(MOD30))
expect(identical(a$afro_subregion, unname(MOD30)),
       "the 30 modelled countries in the wishlist's subregions")
expect(identical(names(MOD30)[a$who_region == "EMRO"], c("SDN", "SOM")) && all(a$who_region %in% c("AFRO", "EMRO")),
       "SDN and SOM flagged EMRO, the others AFRO")
expect(all(is.na(afro_subregion_of(c("Africa", "enso_west", NA))$afro_subregion)), "non-country strings -> NA")
reg <- try(load_grouping("enso_grouping_prod_run_prior"), silent = TRUE)
if (!inherits(reg, "try-error"))
  expect(!anyNA(afro_subregion_of(reg$iso3 %||% reg[[1]])$afro_subregion), "every registry country has a subregion")
df <- tibble(g = c("enso_west", "enso_west", NA), Country = c("SOM", "BEN", "Africa"), x = 1:3)
w <- with_afro(df, "Country")
expect(identical(names(w), c("g", "Country", AFRO_COL, "x")) &&
         identical(w[[AFRO_COL]], c("Eastern Africa‡", "Western Africa", "")),
       "printed: 'AFRO subregion' right after the country column, '‡' on EMRO, empty for Africa")
r <- with_afro(df, raw = TRUE)
expect(identical(names(r), c("g", "Country", "afro_subregion", "who_region", "x")) &&
         identical(r$who_region, c("EMRO", "AFRO", NA)), "raw: afro_subregion + who_region (column auto-detected)")
expect(identical(with_afro(with_afro(df, "Country"), "Country"), w), "idempotent")
expect(identical(with_afro(df[0, ], "Country"), df[0, ]) && identical(with_afro(df["x"]), df["x"]),
       "empty frame or no country column: unchanged")

cat("# [2] data release split\n")
expect(identical(report_data_id(c("57_D12_enso_exposure_by_country.csv", "D30_who_scale.csv",
                                  "03_T3a_effects_occurrence.csv", "12_enso_lags.csv", "D4x.csv")),
                 c("D12", "D30", "T3", NA, NA)), "data ids from CSV names")
td <- tempfile("rel"); dir.create(td)
for (f in c("01_enso_lags.csv", "02_D2_episodes_by_group_country_level.csv", "03_D6_selected_lags.csv",
            "D30_who_scale.csv", "04_T3a_effects_occurrence.csv", "05_D99_new.csv"))
  write.csv(data.frame(a = 1), file.path(td, f), row.names = FALSE)
dir.create(file.path(td, "github")); write.csv(data.frame(a = 1), file.path(td, "github", "stale.csv"))
rel <- list(D2 = "supp", D6 = "github", D30 = "supp", T3 = "supp")
man <- report_release_split(td, rel, c(D6_selected_lags = "Selected lags", D30_who_scale = "WHO scale"))
expect(setequal(list.files(file.path(td, "supp_data")),
                c("D2_episodes_by_group_country_level.csv", "D30_who_scale.csv", "T3a_effects_occurrence.csv")) &&
         identical(list.files(file.path(td, "github")), "D6_selected_lags.csv"),
       "copies under the slug, by release; stale files cleared")
expect(file.exists(file.path(td, "data_release_manifest.csv")) &&
         identical(names(man), c("file", "source", "data_id", "release", "caption")),
       "manifest written with file, source, data_id, release, caption")
expect(identical(man$caption[man$data_id == "D6"], "Selected lags") &&
         identical(man$file[man$data_id == "D30"], "supp_data/D30_who_scale.csv"),
       "captions and release paths in the manifest")
expect(identical(man$release[man$data_id == "D99"], "") && identical(man$file[man$data_id == "D99"], "") &&
         !"01_enso_lags.csv" %in% man$source, "an unmapped id is listed, not copied; main-text tables ignored")
td2 <- tempfile("rel2_"); dir.create(td2)
src <- data.frame(scope = c("Africa", "GW"), lag = c("+6", "-4"), pr = c(1.25, 0.5), pr_lo = NA, ok = c(TRUE, FALSE),
                  note = c("", "a, b"), stringsAsFactors = FALSE)
report_write_csv(src, file.path(td2, "60_D14_event_attrib_pr_scope.csv"))
invisible(report_release_split(td2, list(D14 = "supp")))
rl <- readLines(file.path(td2, "supp_data", "D14_event_attrib_pr_scope.csv"))
back <- read.csv(file.path(td2, "supp_data", "D14_event_attrib_pr_scope.csv"), colClasses = "character", na.strings = character(0))
expect(identical(names(back), c("scope", "lag", "pr", "ok", "note")) &&
         identical(back$lag, c("+6", "-4")) && identical(back$pr, c("1.25", "0.5")) &&
         identical(back$note, c("", "a, b")) && identical(rl[2], '"Africa","+6",1.25,TRUE,""') &&
         "pr_lo" %in% names(read.csv(file.path(td2, "60_D14_event_attrib_pr_scope.csv"))),
       "release copy: all-empty column dropped, cells kept as text, numbered table untouched")
expect(throws(report_release_split(td, list(D2 = "zenodo")), "must be one of"), "a bad release value stops")
manx <- report_release_split(td, rel, exclude = "D6_selected_lags")
expect(!file.exists(file.path(td, "github", "D6_selected_lags.csv")) &&
         identical(manx$release[manx$data_id == "D6"], "withheld") && identical(manx$file[manx$data_id == "D6"], "") &&
         file.exists(file.path(td, "supp_data", "D30_who_scale.csv")),
       "an excluded slug is withheld: listed, not copied; the rest released")
dm <- get_default_options()$results$style$data_release
expect(all(unlist(dm) %in% c("supp", "github")) &&
         identical(sort(names(dm)[unlist(dm) == "supp"]),
                   sort(c("D2", paste0("D", 10:21), "D30", "T3", "T4", "T5"))),
       "default map: supp = D2, D10-D21, D30, T3-T5 (plan U9); the rest github")
expect(identical(get_default_options()$results$style$data_release_exclude, "D16_event_attrib_observed_adm1"),
       "default exclude: the D16 ADM1 observed counts")

cat("# [3] D34 run design, D35 analysis set by level\n")
rd <- tibble(arm = c("forecast", "counterfactual", "events", "smoothing"),
             n_replicates = c(5000L, 5000L, NA, 100L), max_draws = c(200L, 200L, NA, 200L),
             reps_per_draw = c(25, 25, NA, NA), draw_mode = c("cycle", "cycle", NA, "random"),
             handoff = c("shuffled", "shuffled", NA, NA), n_analogs_per_model = c(NA, 100L, NA, NA),
             n_models = c(NA, 10L, NA, NA), clamp_margin = c(0.25, 0.25, NA, 0),
             config = c("a.yml", "b.yml", NA, "c.yml"))
t34 <- report_run_design_table(rd)
expect(identical(t34$Arm, c("Forecast", "Counterfactual", "Smoothing")) &&
         identical(t34$Replicates, c("5,000", "5,000", "100")) && identical(t34$`Replicates per draw`, c("25", "25", "—")) &&
         identical(t34$`Analogs per model`, c("—", "100", "—")) && identical(t34$`Clamp margin`, c("0.25", "0.25", "0.00")) &&
         identical(t34$Handoff[3], "—"),
       "D34: arms without a config dropped, labels, dashes for what does not apply")
t1 <- tibble(admin_level = c("ADM0", "ADM1", "ADM2"), n_outbreaks = c(5L, 100L, 200L), n_gids = c(1L, 40L, 90L),
             n_years = 10L, total_cases = c(1e3, 5e4, 6e4), cases_q25 = 10, cases_med = 100, cases_q75 = 500,
             dur_wk_q25 = 2, dur_wk_med = 6, dur_wk_q75 = 12, pop_q25 = 1e5, pop_med = 5e5, pop_q75 = 1e6,
             ar_q25 = 0.1, ar_med = 0.5, ar_q75 = 1, n_episodes = c(4L, 90L, 180L),
             dur_mo_q25 = 1, dur_mo_med = 2, dur_mo_q75 = 4)
t35 <- report_analysis_set_table(t1)
expect(identical(t35$`Admin level`, c("ADM1", "ADM2")) && identical(t35$Outbreaks, c("100", "200")) &&
         identical(t35$`Cases, median [IQR]`[1], "100 [10–500]") && "Duration (months), median [IQR]" %in% names(t35),
       "D35: ADM0 dropped, Zheng columns, duration months when the artefact has them")

cat("# [4] 10 supplement maps (figures_supp.R): MU, GM, KC on synthetic frames\n")
source("analysis/utils/figures_supp.R")
# The unit classification is 09's (results_derive.R::modelled_units_tables(),
# kind modelled_units); sourced from a throwaway env with its call-time deps.
.md <- new.env()
for (m in c("base", "gid")) sys.source(sprintf("analysis/utils/%s.R", m), envir = .md)
.rd <- parse("analysis/utils/results_derive.R", keep.source = FALSE)
for (e in .rd) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                   as.character(e[[2]]) %in% c("gid_level", "modelled_units_tables", "rain_enso_lag_summary",
                                               "report_lags_loo_summary")) eval(e, .md)
rate <- tibble(gid = c("AGO_1_1", "AGO_2_1_1", "AGO_3_1", "ZZZ_1_1", "GHA5_2", "GHA14_13_2"),
               country = c("AGO", "AGO", "AGO", "ZZZ", "GHA", "GHA"),
               onsets = c(3L, 1L, 0L, 2L, 1L, 1L), enso_group = c("enso_south", "enso_south", "enso_south", NA,
                                                                   "enso_nga", "enso_nga"))
st <- .md$modelled_units_tables(rate, c("AGO.1_1", "AGO.2_1", "AGO.3_1", "AGO.4_1", "ZZZ.2_1", "NAM.1_1",
                                        "GHA5_2", "GHA14_2", "GHA1_2"))
expect(identical(st$units$level, c("ADM1", "ADM2", "ADM1", "ADM1", "ADM2")) &&
         identical(st$units$adm1, c("AGO_1", "AGO_2", "AGO_3", "GHA_5", "GHA_14")),
       "09 modelled_units: admin level and gidAdm1() key of flat gids (Ghana's GADM ids included)")
expect(identical(st$units$status, c("adm1", "adm2", "no_onset", "adm1", "adm2")) &&
         identical(st$gaps$adm1, c("AGO_4", "GHA_1")),
       "units by level, no-onset units, ADM1s of modelled countries without a unit (unfitted countries out)")
st$units <- st$units[1:3, ]
sq <- function(x, y) sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1), c(x, y + 1), c(x, y))))
afro <- sf::st_sf(ISO3_CODE = c("AGO", "NAM", "SOM"), geometry = sf::st_sfc(sq(0, 0), sq(1, 0), sq(0, 1)))
unit_sf <- sf::st_sf(gid = st$units$gid, status = st$units$status,
                     geometry = sf::st_sfc(sq(0, 0), sq(0.2, 0.2), sq(0.5, 0.5)))
gap_sf <- sf::st_sf(adm1 = "AGO_4", status = "no_episode", geometry = sf::st_sfc(sq(0.7, 0.7)))
bb <- sf::st_bbox(afro)
ok <- function(p) !inherits(try(ggplot2::ggplot_build(p), silent = TRUE), "try-error")
expect(ok(build_figS_modelled_units(unit_sf, gap_sf, afro, bb)), "figS_modelled_units builds")
gl <- c(enso_south = "South [GS]", enso_west = "West [GW]")
gm <- build_figS_group_map(afro, tibble(iso3 = c("AGO", "NAM", "SOM"), group = c("enso_south", "enso_south", "enso_west")),
                           c(enso_south = "#7570b3", enso_west = "#66a61e"), gl, c(enso_south = "GS", enso_west = "GW"),
                           bb, low_isos = "SOM")
gb <- ggplot2::ggplot_build(gm)
lab_layer <- Filter(function(d) "label" %in% names(d), gb$data)
expect(length(lab_layer) == 2L && setequal(lab_layer[[2]]$label, c("GS", "GW")) &&
         identical(lab_layer[[1]]$label, "\u2020"),
       "figS_group_map builds, one centroid label per group, a dagger on the low-probability country")
part <- tibble(id = rep(c("ns5_K2", "ns5_K3", "ns4_K2"), each = 2), iso3 = rep(c("AGO", "SOM"), 3),
               group = c("south", "south", "south", "klc3", "west", "west"), p = c(0.9, 0.6, 0.8, 0.5, 0.7, 0.7))
fits <- tibble(id = c("ns5_K2", "ns5_K3", "ns4_K2", "ns5_K2_mix"), basis = c("ns5", "ns5", "ns4", "ns5"),
               K = c(2L, 3L, 2L, 2L), mix = c(FALSE, FALSE, FALSE, TRUE))
d <- kl_classes_by_K_frame(part, fits, "ns5_K5")
expect(identical(sort(unique(d$id)), c("ns5_K2", "ns5_K3")) && all(d$K %in% 2:3),
       "kl classes by K: the pinned basis and prior only")
expect(throws(kl_classes_by_K_frame(part, fits, "bs4_K5"), "no sweep fit"), "a basis without fits stops")
kc <- build_figS_kl_classes_by_K(d, c(AGO = "enso_south", SOM = "enso_west"), c(enso_south = "#7570b3", enso_west = "#66a61e"),
                                 c(enso_west = "GW", enso_south = "GS"), ggplot2::theme_minimal())
kb <- ggplot2::ggplot_build(kc)
expect(ok(kc) && identical(levels(kb$plot$data$class_id), c("enso_west", "enso_south", "klc3")) &&
         identical(levels(kb$plot$data$iso3), c("AGO", "SOM")),
       "figS_kl_classes_by_K builds; deployed classes coloured, klc grey; rows by deployed group")

cat("# [5] wishlist T1-T3 builders\n")
ORD <- c("enso_west", "enso_nga", "enso_horn", "enso_east", "enso_south")
reg <- load_grouping("enso_grouping_prod_run_prior")
MECH <- list(path = "analysis/tables/enso_rain_mechanisms.yml",
             groups = yaml::read_yaml("analysis/tables/enso_rain_mechanisms.yml")$groups)   # 09 rain_enso$mechanisms
t1a <- report_t1a_mechanisms(MECH,
                             split(as.character(reg$iso3), as.character(reg$group)), ORD)
expect(identical(t1a$enso_group, ORD) && all(nzchar(t1a$mechanism)) && all(nzchar(t1a$references)) &&
         identical(t1a$members[t1a$enso_group == "enso_nga"], "GHA NGA"),
       "T1a: one curated row per group in display order, members from the registry")
expect(throws(report_t1a_mechanisms(MECH, list(), c(ORD, "enso_x")),
              "no entry"), "T1a: a group without an entry stops")
corr <- expand.grid(series = c("enso_east", "enso_south"), centre = 1:12, lag = -2:2, stringsAsFactors = FALSE)
corr$r <- 0.1; corr$p_r <- 0.5; corr$n <- 40L; corr$agree_r <- 0.5
corr$r[corr$series == "enso_east" & corr$centre == 11 & corr$lag == -2] <- 0.54
corr$r[corr$series == "enso_east" & corr$centre == 6 & corr$lag == 0] <- 0.9      # out of season: ignored
corr$r[corr$series == "enso_south" & corr$centre == 2 & corr$lag == 0] <- -0.55
mon <- expand.grid(group = c("enso_east", "enso_south"), mo = 1:12, stringsAsFactors = FALSE)
mon$weighting <- "pop"
mon$in_season <- (mon$group == "enso_east" & mon$mo %in% c(3:5, 11:12)) | (mon$group == "enso_south" & mon$mo %in% c(10:12, 1:4))
sea <- data.frame(weighting = "pop", group = c("enso_east", "enso_east", "enso_south"),
                  months = c("Mar-Apr-May", "Nov-Dec", "Oct-Nov-Dec-Jan-Feb-Mar-Apr"))
sel <- tibble(track = rep(c("occurrence", "duration", "monthly"), each = 2), enso_group = rep(c("enso_east", "enso_south"), 3),
              deployed_lag = c(0L, 4L, -1L, 11L, 2L, 5L))
rain <- list(corr = tibble::as_tibble(corr), month = mon, seasons = sea)
rain$group_lags <- .md$rain_enso_lag_summary(rain, c("enso_south", "enso_east"), "pop")   # 09, rain_enso v2
t1b <- report_t1b_rain_lags(rain, sel, c("enso_east", "enso_south"))
expect(identical(t1b$best_centre_month, c("Nov", "Feb")) && identical(t1b$best_lag, c(-2L, 0L)) &&
         isTRUE(all.equal(t1b$best_r, c(0.54, -0.55))) && identical(t1b$rainy_season[1], "Mar-Apr-May; Nov-Dec") &&
         identical(t1b$lag_dur, c(-1L, 11L)) && identical(t1b$lag_mon, c(2L, 5L)),
       "T1b: strongest in-season |r| (out-of-season cells ignored), seasons, deployed lags per track")
expect(identical(names(t1b), c("enso_group", "rainy_season", "best_centre_month", "best_lag", "best_r", "best_p",
                               "best_n", "best_agree", "season_mean_peak_lag", "season_mean_peak_r",
                               "lag_occ", "lag_dur", "lag_mon")),
       "T1b: CSV columns unchanged (09 rain_enso group_lags + the deployed lags joined in 11)")
expect(throws(report_t1b_rain_lags(rain[c("corr", "month", "seasons")], sel, "enso_east"), "group_lags"),
       "T1b: a rain_enso artefact without group_lags (v1) stops")
rk <- tibble(track = "occurrence", enso_group = "enso_east", single_lag = c(0L, 1L, -8L, 5L),
             rank = 1:4, elpd_loo = c(-910.1, -913.8, -930, -964.8), se_elpd = c(58.9, 59, 60, 61),
             elpd_diff = c(0, -3.7, -20, -54.7), se_diff = c(0, 3.3, 5, 11), p_loo = c(42.2, 41, 40, 39))
sel2 <- tibble(track = "occurrence", enso_group = "enso_east", deployed_lag = 1L)
t2 <- .md$report_lags_loo_summary(rk, sel2, 2)   # 09, report_lags v3 loo_summary
expect(identical(t2$rank1_lag, 0L) && identical(t2$deployed_lag, 1L) && identical(t2$deployed_elpd_diff, -3.7) &&
         identical(t2$rank2_lag, 1L) && identical(t2$lag0_elpd_diff, 0) && identical(t2$worst_lag, 5L) &&
         identical(t2$n_within_se, 2L) && identical(t2$elpd_loo_rank1, -910.1) && identical(t2$p_loo_rank1, 42.2),
       "T2: rank 1 / deployed / rank 2 / lag 0 / worst, lags within 2 SE (rank 1 included), rank-1 elpd and p_loo")
expect(identical(report_outside_range(c(-1.5, 1.0, 2.0, 1.0), c(-1.02, -1.02, -1.57, NA), c(1.44, 1.44, 2.59, NA)),
                 c(TRUE, FALSE, FALSE, FALSE)), "T3: outside the fitted range flagged; NA range never flagged")

cat("# [6] D35b / D35c by WHO-AFRO subregion\n")
expect(identical(AFRO_SUBREGION_ORDER, c("Western Africa", "Central Africa", "Eastern Africa", "Southern Africa")) &&
         setequal(AFRO_SUBREGION_ORDER, unique(AFRO_SUBREGION)), "subregion display order covers every subregion")
er <- tibble(region = rep(c("Africa", "Eastern Africa"), each = 2), enso_state = rep(c("La Niña", "All"), 2),
             share = c(0.4, 1, 0.25, 1), n_outbreaks = c(4L, 10L, 1L, 4L), n_years = c(2L, 5L, 1L, 3L),
             n_state_years = c(2L, 5L, 2L, 5L), per_year = c(2, 2, 0.5, 0.8),
             cases_q25 = 1, cases_med = 2, cases_q75 = 3, dur_wk_q25 = 1, dur_wk_med = 2, dur_wk_q75 = 3,
             ar_q25 = 0.1, ar_med = 0.2, ar_q75 = 0.3)
erp <- tibble(region = rep(c("Africa", "Eastern Africa"), each = 3),
              variable = rep(c("cases", "duration_weeks", "attack_rate_1k"), 2),
              p_value = c(0.0004, 0.2, NA, 0.5, 0.01, 0.03))
rw <- report_enso_region_table(er, erp, raw = TRUE)
expect(nrow(rw) == 4L && all(c("p_cases", "p_duration_weeks", "p_attack_rate_1k") %in% names(rw)) &&
         identical(rw$p_cases, c(0.0004, 0.0004, 0.5, 0.5)), "D35c raw: the block's p joined on every state row")
pr <- report_enso_region_table(er, erp)
expect(identical(pr$`Outbreaks (share)`[1], "4 (40.0%)") && identical(pr$`Outbreaks per season-year`[3], "0.5") &&
         identical(pr$`p cases / duration / attack rate`[1], "<0.001 / 0.200 / —"),
       "D35c printed: n (share %), p formatted, NA p as a dash")
expect(is.null(report_enso_region_table(er[0, ], erp)) && is.null(report_analysis_set_region_table(NULL)),
       "D35b / D35c: no rows -> NULL")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
