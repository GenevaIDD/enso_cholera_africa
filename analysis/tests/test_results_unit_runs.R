# test_results_unit_runs.R -- the results phase on a v2 (partition_by:
# registry) manifest.
#
#   [1] group_run_tag() / group_config_path() resolve a fit group to the unit
#       that holds it; results_run_configs() lists every unit config (v2) and
#       keeps the v1 order (per group: forecast, then each set).
#   [2] load_group_arm_runs() on a unit mixing two fit groups: one load per
#       unit run, rows relabelled with each country's FIT group, each country
#       censored at its own firm horizon; group_horizon = the groups' firm
#       horizons; gids carry the unit.
#   [3] results.scope.beyond_firm_horizon = "flag" keeps the later months,
#       flagged.
#   [4] a stamp predating beyond_firm_horizon is read as "drop" (no
#       re-derivation forced); a stamp with "flag" against "drop" differs.
#   [5] per-arm launches: manifest_arm_launches() from the
#       manifest keys or, without them, each config's horizon.start;
#       manifest_attrib_start() under stitch / common_start.
#   [6] load_sim_run_stitched(): the months before the factual launch come from
#       the smoothing events, paired through handoff_map.csv, zero-filled;
#       forecast months untouched; same columns; equal starts == plain load;
#       a run without the map, or naming a foreign smoothing replicate, is
#       refused; load_arm_run() reads a counterfactual plainly; the 09 loader
#       stitches under "stitch" and not under "common_start".
#   [7] a stamp predating derivation.arm_alignment reads as "stitch".
#   [8] attribution factual: a stamp predating
#       inputs.attrib_manifest reads as ""; results_attrib_ctx() /
#       results_attrib_manifest_files() / display_factual_start(); the post_fork
#       cell starts at the displayed factual's launch.
#   [9] a stamp holding a quoted "60" (older aligner output) reads as the config's
#       integer 60 (no decision diff); 50 differs.
#  [10] table S13's simulation environment: sim_software_versions()
#       (R, platform, SIM_SOFTWARE_PKGS; no CmdStan); report_software_sim_versions()
#       on synthetic partition_meta.yml files: one row per arm x unit x
#       software, distinct across sets; NA rows for a run without `software`
#       or without a partition_meta.yml; .report_software_sim_runs() picks the
#       displayed forecast, the attribution factual only when it is its own
#       run, and every set's counterfactual.
#
# Synthetic temp files only. Run: Rscript analysis/tests/test_results_unit_runs.R

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(qs2) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/results_core.R")
source("analysis/utils/results_derive.R")

td <- tempfile("resunit_"); dir.create(td)
cfg <- function(tag) { p <- file.path(td, paste0(tag, ".yml"))
  yaml::write_yaml(list(simulator = list(run_tag = tag, horizon = list(start = "2026-11-01"))), p); p }
months <- seq(as.Date("2026-11-01"), as.Date("2027-08-01"), by = "month")
part <- function(tag) {
  pd <- file.path(td, "sim", tag, "parts", "tu=eh"); dir.create(pd, recursive = TRUE)
  ms <- tidyr::expand_grid(replicate_id = 1:2, gid = c("KEN_1_1", "ETH_2_1"), month = format(months)) %>%
    mutate(Y = 1L, cases_month = 3)
  qs_save(ms, file.path(pd, "monthly_summary.qs")); file.create(file.path(pd, "_SUCCESS"))
}
part("fc_eh"); part("cf_eh_neutral")
man <- list(schema = 2L, primary_analog_set = "neutral",
  groups = list(enso_east = list(horizon_end = "2027-01-01", iso3 = list("KEN"), units = list("eh")),
                enso_horn = list(horizon_end = "2027-06-01", iso3 = list("ETH"), units = list("eh"))),
  units = list(eh = list(members = list("KEN", "ETH"), groups = list("enso_east", "enso_horn"),
                         horizon_end = "2027-06-01", factual_config = cfg("fc_eh"),
                         primary_analog_set = "neutral", counterfactual_config = cfg("cf_eh_neutral"),
                         counterfactual_sets = list(neutral = list(config = file.path(td, "cf_eh_neutral.yml"))))),
  countries = list(KEN = list(group = "enso_east", unit = "eh", firm_horizon_end = "2027-01-01"),
                   ETH = list(group = "enso_horn", unit = "eh", firm_horizon_end = "2027-06-01")))

cat("\n[1] run lookups\n")
expect(identical(group_run_tag(man, "enso_east", "forecast"), "fc_eh") &&
         identical(group_run_tag(man, "enso_horn", "counterfactual", "neutral"), "cf_eh_neutral"),
       "fit group -> its unit's run tags")
expect(identical(unname(results_run_configs(man, c("enso_east", "enso_horn"), "neutral")),
                 c(file.path(td, "fc_eh.yml"), file.path(td, "cf_eh_neutral.yml"))),
       "v2: one entry per unit (forecast, set)")
man1 <- list(groups = list(a = list(factual_config = "fa", primary_analog_set = "all", counterfactual_config = "ca",
                                     counterfactual_sets = list(nina = list(config = "cn"))),
                           b = list(factual_config = "fb", primary_analog_set = "all", counterfactual_config = "cb",
                                     counterfactual_sets = list(nina = list(config = "cnb")))))
expect(identical(unname(results_run_configs(man1, c("a", "b"), c("all", "nina"))), c("fa", "ca", "cn", "fb", "cb", "cnb")),
       "v1: per group, forecast then each set (order unchanged)")

cat("\n[2] load_group_arm_runs() on a mixed unit\n")
ctx <- new.env()
ctx$MAN <- man; ctx$SIM_ROOT <- file.path(td, "sim")
ctx$RES <- list(series = "latent", scope = list(beyond_firm_horizon = "drop"))
r <- suppressMessages(load_group_arm_runs(ctx, c("enso_east", "enso_horn"), "neutral", as.Date("2026-11-01")))
cm <- bind_rows(r$country_month)
expect(setequal(cm$enso_group, c("enso_east", "enso_horn")) &&
         all(cm$enso_group[cm$country == "KEN"] == "enso_east"), "rows relabelled with the country's fit group")
expect(max(cm$month[cm$country == "KEN"]) == as.Date("2027-01-01") &&
         max(cm$month[cm$country == "ETH"]) == as.Date("2027-06-01"), "each country stops at its firm horizon")
expect(setequal(unique(cm$run), c("forecast", "counterfactual")) && all(cm$unit == "eh"), "both arms, unit column")
expect(identical(r$group_horizon, c(enso_east = "2027-01-01", enso_horn = "2027-06-01")), "group_horizon = firm")
expect(all(as.character(cm$group_horizon_end[cm$enso_group == "enso_east"]) == "2027-01-01"), "group_horizon_end per row")
expect(identical(sort(r$gids$enso_group), c("enso_east", "enso_horn")) && all(r$gids$unit == "eh"), "gids carry group + unit")
am <- bind_rows(r$adm1_month)
expect(max(am$month[startsWith(am$adm1_unit, "KEN")]) == as.Date("2027-01-01"), "adm1 rows censored too")
r1 <- suppressMessages(load_group_arm_runs(ctx, "enso_horn", "neutral", as.Date("2026-11-01")))
expect(all(bind_rows(r1$country_month)$enso_group == "enso_horn"), "only the requested groups are returned")

cat("\n[2b] Africa window past a firm horizon (results.scope.africa, 2026-10-27)\n")
ctx$AFRICA_WIN <- results_africa_window(c(man, list(last_cfsv2_month = "2027-03-01")), c("enso_east", "enso_horn"),
                                        list(window_end = "2027-04-01", continuation_max_lead = "to_end"))
rx <- suppressMessages(load_group_arm_runs(ctx, c("enso_east", "enso_horn"), "neutral", as.Date("2026-11-01")))
cx2 <- bind_rows(rx$country_month)
expect(max(cx2$month[cx2$country == "KEN"]) == as.Date("2027-04-01") &&
         max(cx2$month[cx2$country == "ETH"]) == as.Date("2027-06-01"),
       "KEN read on to its end in the Africa window (2027-04), ETH still to its firm horizon")
expect(all(as.character(cx2$group_horizon_end[cx2$enso_group == "enso_east"]) == "2027-01-01") &&
         identical(rx$group_horizon, r$group_horizon), "group_horizon_end / group_horizon stay the firm horizons")
ctx$AFRICA_WIN <- NULL

cat("\n[3] beyond_firm_horizon = flag\n")
ctx$RES$scope$beyond_firm_horizon <- "flag"
rf <- suppressMessages(load_group_arm_runs(ctx, c("enso_east", "enso_horn"), "neutral", as.Date("2026-11-01")))
cf <- bind_rows(rf$country_month)
expect("beyond_firm_horizon" %in% names(cf) && max(cf$month[cf$country == "KEN"]) > as.Date("2027-01-01") &&
         all(cf$beyond_firm_horizon[cf$country == "KEN" & cf$month > as.Date("2027-01-01")]),
       "flag keeps the later months, flagged")

cat("\n[4] stamp back-compat\n")
res <- list(inputs = list(a = 1), scope = list(main_set = "neutral", beyond_firm_horizon = "drop"),
            derivation = list(), thresholds = list())
old_stamp <- list(results = list(inputs = list(a = 1), scope = list(main_set = "neutral"),
                                 derivation = list(), thresholds = list()))
expect(length(.results_decision_diff(old_stamp, res)) == 0L, "stamp without the key reads as 'drop'")
res_f <- res; res_f$scope$beyond_firm_horizon <- "flag"
expect(length(.results_decision_diff(old_stamp, res_f)) == 1L, "'flag' against an old stamp re-derives")

cat("\n[5] per-arm launches\n")
al <- manifest_arm_launches(man)
expect(nrow(al) == 2L && all(al$source == "config") && all(al$launch == as.Date("2026-11-01")),
       "v2 manifest without launch keys: launches from the configs' horizon.start")
man_k <- man
man_k$units$eh$factual_launch <- "2027-01-01"
man_k$units$eh$counterfactual_sets$neutral$launch <- "2026-11-01"
al_k <- manifest_arm_launches(man_k)
expect(all(al_k$source == "manifest") &&
         identical(format(al_k$launch[al_k$arm == "factual"]), "2027-01-01"),
       "launch keys read from the manifest")
expect(manifest_attrib_start(man_k) == as.Date("2026-11-01") &&
         manifest_attrib_start(man_k, alignment = "common_start") == as.Date("2027-01-01"),
       "stitch = earlier launch, common_start = later launch")
man_bad <- man_k; man_bad$units$eh$counterfactual_sets$neutral$launch <- "2027-03-01"
e <- tryCatch({ manifest_arm_starts(man_bad); NULL }, error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("after the factual arm", e, fixed = TRUE), "counterfactual launching later refused")
expect(identical(smoothing_parts_dir_of("generated_data/sim_spatpomp/h_run/cache/"),
                 "generated_data/sim_spatpomp/h_run/parts") &&
         !is.null(tryCatch({ smoothing_parts_dir_of("x/parts"); NULL }, error = function(e) 1)),
       "smoothing_parts_dir_of(): cache/ -> parts/, anything else refused")

cat("\n[6] factual series stitched from the smoothing\n")
sm_parts <- file.path(td, "sm", "parts")
sm_pd <- file.path(sm_parts, "tu=eh", "chunk_0001"); dir.create(sm_pd, recursive = TRUE)
ev <- tibble(rep_id = c(1L, 2L, 3L), t_start = as.Date(c("2026-10-01", "2026-11-01", "2026-12-01")),
             gid = c("KEN_1_1", "KEN_1_1", "ETH_2_1"), country = c("KEN", "KEN", "ETH"), Yout = 1L)
ev_mon <- tibble(rep_id = c(1L, 2L, 2L, 3L),
                 t_start = as.Date(c("2026-10-01", "2026-11-01", "2026-12-01", "2026-12-01")),
                 gid = c("KEN_1_1", "KEN_1_1", "KEN_1_1", "ETH_2_1"), country = c("KEN", "KEN", "KEN", "ETH"),
                 cases_month = c(99, 10, 5, 7))
qs_save(list(events = ev, events_mon = ev_mon), file.path(sm_pd, "part.qs")); invisible(file.create(file.path(sm_pd, "_SUCCESS")))
fmonths <- seq(as.Date("2027-01-01"), as.Date("2027-03-01"), by = "month")
fpart <- function(tag, map = TRUE) {
  pd <- file.path(td, "sim", tag, "parts", "tu=eh"); dir.create(pd, recursive = TRUE)
  ms <- tidyr::expand_grid(replicate_id = 1:4, gid = c("KEN_1_1", "ETH_2_1"), month = format(fmonths)) %>%
    mutate(Y = 0L, cases_month = 1)
  qs_save(ms, file.path(pd, "monthly_summary.qs")); file.create(file.path(pd, "_SUCCESS"))
  if (map) utils::write.csv(data.frame(tu = "eh", replicate_id = 1:4, smoothing_rep_id = c(2L, 3L, 2L, 1L),
                                       handoff_month = "2027-01-01", handoff_weighting = "uniform"),
                            file.path(pd, "handoff_map.csv"), row.names = FALSE)
}
fpart("fc_st"); fpart("fc_nomap", map = FALSE)
sev <- read_smoothing_parts_events(sm_parts, as.Date("2011-01-01"))
SIM <- file.path(td, "sim")
st <- suppressMessages(load_sim_run_stitched(SIM, "fc_st", "2026-11-01", "2027-03-01", "2027-01-01", sev))
pl <- load_sim_run_monthly(SIM, "fc_st", "2027-01-01", "2027-03-01")
cmst <- st$country_month
cell <- function(rid, ctry, m, col) cmst[[col]][cmst$replicate_id == rid & cmst$country == ctry & cmst$month == as.Date(m)]
expect(identical(cell(1L, "KEN", "2026-11-01", "n_cases"), 10) && identical(cell(1L, "KEN", "2026-11-01", "n_outbreaks"), 1L) &&
         identical(cell(1L, "KEN", "2026-12-01", "n_cases"), 5) && identical(cell(3L, "KEN", "2026-11-01", "n_cases"), 10) &&
         identical(cell(2L, "ETH", "2026-12-01", "n_cases"), 7) && identical(cell(2L, "ETH", "2026-12-01", "n_outbreaks"), 1L),
       "prefix = the paired smoothing replicate's events (one smoothing replicate may feed two)")
expect(sum(cmst$n_cases[cmst$replicate_id == 4L & cmst$month < as.Date("2027-01-01")]) == 0,
       "smoothing events before the window are not stitched in")
expect(sum(cmst$month < as.Date("2027-01-01")) == 4L * 2L * 2L,
       "prefix is a full replicate x country x month grid (zero-filled)")
fc_rows <- function(x) x %>% filter(month >= as.Date("2027-01-01")) %>% arrange(replicate_id, country, month)
expect(identical(fc_rows(cmst), fc_rows(pl$country_month)) &&
         identical(names(cmst), names(pl$country_month)) && identical(names(st$adm1_month), names(pl$adm1_month)),
       "forecast months untouched, same columns as load_sim_run_monthly()")
expect(min(st$adm1_month$month) == as.Date("2026-11-01") && identical(st$gids, pl$gids),
       "ADM1 table stitched over the window too; simulated gids unchanged")
expect(identical(st$stitch$months, as.Date(c("2026-11-01", "2026-12-01"))), "stitch element names the prefix months")
expect(identical(load_sim_run_stitched(SIM, "fc_st", "2027-01-01", "2027-03-01", "2027-01-01", sev), pl),
       "window starting at the factual launch: identical to load_sim_run_monthly()")
e <- tryCatch({ suppressMessages(load_sim_run_stitched(SIM, "fc_nomap", "2026-11-01", "2027-03-01", "2027-01-01", sev)); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("handoff_map.csv", e, fixed = TRUE), "a run without handoff_map.csv is refused")
sev_foreign <- sev; sev_foreign$rep_ids <- 1:2
e <- tryCatch({ suppressMessages(load_sim_run_stitched(SIM, "fc_st", "2026-11-01", "2027-03-01", "2027-01-01", sev_foreign)); NULL },
              error = function(e) conditionMessage(e))
expect(!is.null(e) && grepl("another smoothing run", e, fixed = TRUE), "a map naming a foreign smoothing replicate is refused")
expect(identical(load_arm_run(SIM, "fc_st", "counterfactual", "2026-11-01", "2027-03-01",
                              arm_starts = list(factual = as.Date("2027-01-01"), counterfactual = as.Date("2026-11-01")),
                              smoothing_events = sev),
                 load_sim_run_monthly(SIM, "fc_st", "2026-11-01", "2027-03-01")),
       "load_arm_run(): a counterfactual is read plainly")
# 09 path: v2 manifest, factual launched 2027-01, counterfactual 2026-11.
man_s <- man
man_s$units$eh$factual_config <- { p <- file.path(td, "fc_st.yml")
  yaml::write_yaml(list(simulator = list(run_tag = "fc_st", horizon = list(start = "2027-01-01"))), p); p }
ctx2 <- new.env()
ctx2$MAN <- man_s; ctx2$SIM_ROOT <- SIM; ctx2$SMOOTHING_PARTS_DIR <- sm_parts
ctx2$RES <- list(series = "latent", scope = list(beyond_firm_horizon = "drop"))
for (al_mode in c("stitch", "common_start")) {
  ctx2$ARM_ALIGNMENT <- al_mode
  w <- suppressMessages(results_arm_window(ctx2, c("enso_east", "enso_horn"), "neutral"))
  r2 <- suppressMessages(load_group_arm_runs(ctx2, c("enso_east", "enso_horn"), "neutral", w$attrib_start, win = w))
  f2 <- bind_rows(r2$country_month) %>% filter(run == "forecast")
  if (al_mode == "stitch")
    expect(w$attrib_start == as.Date("2026-11-01") && min(f2$month) == as.Date("2026-11-01") && length(r2$stitch) == 1L &&
             identical(w$ratio_label, "ENSO + post-handoff information"),
           "09 loader, stitch: window from the counterfactual launch, factual arm stitched")
  else
    expect(w$attrib_start == as.Date("2027-01-01") && min(f2$month) == as.Date("2027-01-01") && !length(r2$stitch) &&
             is.null(w$smoothing_events),
           "09 loader, common_start: window from the factual launch, no stitch, no smoothing read")
}

cat("\n[7] stamp back-compat: derivation.arm_alignment\n")
res7 <- list(inputs = list(a = 1), scope = list(main_set = "neutral"),
             derivation = list(smoothing_start = "2011-01-01", arm_alignment = "stitch"), thresholds = list())
old7 <- list(results = list(inputs = list(a = 1), scope = list(main_set = "neutral"),
                            derivation = list(smoothing_start = "2011-01-01"), thresholds = list()))
expect(length(.results_decision_diff(old7, res7)) == 0L, "stamp without the key reads as 'stitch'")
res7c <- res7; res7c$derivation$arm_alignment <- "common_start"
expect(length(.results_decision_diff(old7, res7c)) == 1L, "'common_start' against an old stamp re-derives")

cat("\n[8] attribution factual (2026-10-27)\n")
res8 <- list(inputs = list(a = 1, attrib_manifest = ""), scope = list(main_set = "neutral"),
             derivation = list(), thresholds = list())
old8 <- list(results = list(inputs = list(a = 1), scope = list(main_set = "neutral"),
                            derivation = list(), thresholds = list()))
expect(length(.results_decision_diff(old8, res8)) == 0L, "stamp without inputs.attrib_manifest reads as ''")
res8a <- res8; res8a$inputs$attrib_manifest <- "analysis/configs/forecast_attribution_x.yml"
expect(length(.results_decision_diff(old8, res8a)) == 1L, "an attribution manifest against an old stamp re-derives")
# The attribution view of man_k: the factual launches with the counterfactual.
man_a <- man_k; man_a$units$eh$factual_launch <- "2026-11-01"
man_a$units$eh$factual_config <- sub("\\.yml$", "_attrib.yml", man_k$units$eh$factual_config)
ctx8 <- list2env(list(EXPLOSIVE_CAP = 2e5, MAN = man_k, MAN_DISPLAY = man_k, MAN_ATTRIB = man_a, MANIFEST_PATH = "fh.yml",
                      MANIFEST_PATH_DISPLAY = "fh.yml", ATTRIB_MANIFEST_FILE = "fa.yml",
                      ATTRIB_FACTUAL = "attrib"), parent = emptyenv())
cx <- results_attrib_ctx(ctx8)
expect(identical(cx$MAN$units$eh$factual_launch, "2026-11-01") && identical(cx$MANIFEST_PATH, "fa.yml") &&
         identical(ctx8$MAN$units$eh$factual_launch, "2027-01-01"),
       "results_attrib_ctx(): a copy reading the attribution manifest; the context is untouched")
expect(identical(results_attrib_manifest_files(cx), c("fa.yml", "fh.yml")) &&
         identical(results_attrib_manifest_files(ctx8), "fh.yml"),
       "results_attrib_manifest_files(): both manifests under the attribution context, one otherwise")
expect(display_factual_start(cx, "enso_east", "neutral") == as.Date("2027-01-01") &&
         manifest_arm_starts(cx$MAN, "enso_east", "neutral")$factual == as.Date("2026-11-01"),
       "display_factual_start(): the displayed factual's launch, under the attribution context too")
ctx8f <- list2env(as.list(ctx8), parent = emptyenv()); ctx8f$ATTRIB_FACTUAL <- "forecast"
expect(identical(results_attrib_ctx(ctx8f), ctx8f), "results_attrib_ctx(): unchanged without an attribution factual")
cells <- .cell_starts(attrib_cells(cx, as.Date("2026-11-01"), display_factual_start(cx, "enso_east", "neutral")))
expect(cells[["attribution"]] == as.Date("2026-11-01") && cells[["post_fork"]] == as.Date("2027-01-01"),
       "post_fork starts at the displayed factual's launch, not the attribution factual's")

cat("\n[9] stamp normaliser: inputs.episode_national_zero_threshold\n")
res9 <- list(inputs = list(a = 1, episode_national_zero_threshold = 60L), scope = list(main_set = "neutral"),
             derivation = list(), thresholds = list())
old9 <- list(results = list(inputs = list(a = 1, episode_national_zero_threshold = "60"),
                            scope = list(main_set = "neutral"), derivation = list(), thresholds = list()))
expect(length(.results_decision_diff(old9, res9)) == 0L, "stamp \"60\" vs config 60L: no decision diff")
res9b <- res9; res9b$inputs$episode_national_zero_threshold <- 50L
d9 <- .results_decision_diff(old9, res9b)
expect(length(d9) == 1L && grepl("episode_national_zero_threshold", d9, fixed = TRUE),
       "stamp \"60\" vs config 50L: one diff, on the threshold")

cat("\n[10] simulation software (table S13)\n")
sv <- sim_software_versions()
expect(is.list(sv) && identical(names(sv)[1:2], c("R", "platform")) && !"CmdStan" %in% names(sv) &&
         all(c("spatPomp", "pomp") %in% names(sv)) && identical(names(sv)[-(1:2)], SIM_SOFTWARE_PKGS),
       "sim_software_versions(): R, platform and SIM_SOFTWARE_PKGS as a named list, no CmdStan")
expect(identical(software_versions()$software[1:3], c("R", "platform", "CmdStan")),
       "software_versions() default rows unchanged (R, platform, CmdStan first)")
mk_run <- function(tag, software = NULL, meta = TRUE) {
  md <- file.path(td, "s13", tag, "merged"); dir.create(md, recursive = TRUE)
  if (meta) yaml::write_yaml(c(list(run_tag = tag), if (!is.null(software)) list(software = software)),
                             file.path(dirname(md), "partition_meta.yml"))
  list(config = file.path(td, paste0(tag, ".yml")), merged_dir = md)
}
sw1 <- list(R = "4.3.2", platform = "x86_64", spatPomp = "0.36.0", pomp = "5.11")
runs10 <- list(forecast = list(eh = mk_run("fc10", sw1), wu = mk_run("fc10w", sw1)),
               counterfactual = c(list(eh = mk_run("cf10n", sw1)), list(eh = mk_run("cf10a", sw1)),
                                  list(wu = mk_run("cf10w"))),
               attribution = list(eh = mk_run("at10", meta = FALSE)))
sv10 <- report_software_sim_versions(runs10)
expect(identical(names(sv10), c("arm", "unit", "software", "version")),
       "columns arm, unit, software, version")
expect(nrow(dplyr::filter(sv10, arm == "forecast")) == 8L &&
         nrow(dplyr::filter(sv10, arm == "counterfactual", unit == "eh")) == 4L,
       "one row per arm x unit x software; the two sets' identical eh rows are one")
na10 <- dplyr::filter(sv10, is.na(software))
expect(nrow(na10) == 2L && setequal(paste(na10$arm, na10$unit), c("counterfactual wu", "attribution eh")) &&
         all(is.na(na10$version)),
       "NA row for a run without `software` (an older run) and one without partition_meta.yml")
expect(identical(sv10$version[sv10$arm == "forecast" & sv10$unit == "eh" & sv10$software == "spatPomp"], "0.36.0"),
       "the recorded version is read as is")
expect(nrow(report_software_sim_versions(list())) == 0L, "no runs: an empty table with the same columns")
ctx10 <- list2env(list(MAN = man, MAN_DISPLAY = man, MAN_ATTRIB = man, ATTRIB_FACTUAL = "forecast",
                       SETS = "neutral"), parent = emptyenv())
r10 <- .report_software_sim_runs(ctx10)
expect(identical(names(r10), c("forecast", "counterfactual")) && identical(names(r10$forecast), "eh") &&
         identical(r10$counterfactual$eh$config, file.path(td, "cf_eh_neutral.yml")),
       ".report_software_sim_runs(): forecast + each set's counterfactual; no attribution arm without its own factual")
ctx10$ATTRIB_FACTUAL <- "attrib"
expect(identical(names(.report_software_sim_runs(ctx10)), c("forecast", "attribution", "counterfactual")),
       "the attribution factual is its own arm when it is its own run")

unlink(td, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
