# test_results_ipcc_scope.R -- report_attrib_metrics v2's scope-level IPCC
# class table `ipcc_scope` (HARD gate, quick).
#
# Figure 4a's ratio bars are filled by the IPCC class of each scope (Africa,
# ENSO regions) at the headline m. 09 derives it in attrib_metrics_tables()
# from the scope's summed paired totals (scope_pair_vectors()).
#
#   [1] ipcc_class() == the case_when it replaces (attrib_exposure_classes()).
#   [2] ipcc_scope == a reference pd_scope_class() (below) on a synthetic
#       two-group, two-set run with an explosive replicate: Africa + regions,
#       both metrics, both sets, three magnitudes, explosive guard on.
#   [3] without ipcc_mags, attrib_metrics_tables() returns the v1 tables only
#       (the event cores' call), and with it every v1 table is unchanged.
#
# Synthetic inputs only. Run: Rscript analysis/tests/test_results_ipcc_scope.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
PASS <- 0L; FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/gid.R")
source("analysis/utils/paths.R")
source("analysis/utils/sim_io.R")
source("analysis/utils/results_core.R")
source("analysis/utils/results_derive.R")

cat("\n[1] ipcc_class() == case_when\n")
set.seed(3)
pe <- runif(500); ps <- runif(500); pe[1:20] <- 0.66; ps[1:20] <- 0.8
ref <- dplyr::case_when(pe >= .66 & ps >= .8 ~ "robust", pe >= .66 ~ "conflicting", TRUE ~ "none")
expect(identical(ipcc_class(pe, ps, .66, .8), ref), "500 share pairs, boundaries included")

# Reference pd_draw_pairs() / pd_scope_class(), without the draw join, on a
# draw-cache layout: window totals per (event, enso_group, unit,
# country, run, analog_set, replicate_id) and the explosive keys of each window.
pd_scope_class_ref <- function(dr, metric, aset, m, robust_exceed, robust_agree) {
  reps <- dplyr::filter(dr$reps, run == "forecast" | analog_set == aset)
  reps <- dplyr::anti_join(reps, dplyr::filter(dr$explosive, run == "forecast" | analog_set == aset),
                           by = c("event", "enso_group", "run", "analog_set", "replicate_id"))
  keys <- c("event", "enso_group", "unit", "country", "replicate_id")
  arm <- function(rn, col) reps %>% dplyr::filter(run == rn) %>%
    dplyr::select(dplyr::all_of(keys), !!col := dplyr::all_of(metric))
  pairs <- dplyr::inner_join(arm("forecast", "f"), arm("counterfactual", "cf"), by = keys)
  grp <- pairs %>% dplyr::group_by(enso_group, replicate_id) %>%
    dplyr::summarise(f = sum(f), cf = sum(cf), .groups = "drop")
  ids <- Reduce(intersect, split(grp$replicate_id, grp$enso_group))
  afr <- grp %>% dplyr::filter(replicate_id %in% ids) %>% dplyr::group_by(replicate_id) %>%
    dplyr::summarise(f = sum(f), cf = sum(cf), .groups = "drop") %>% dplyr::mutate(enso_group = "Africa")
  dplyr::bind_rows(grp, afr) %>% dplyr::group_by(key = enso_group) %>%
    dplyr::summarise(n_pairs = dplyr::n(), p_exceed = mean(f > (1 + m) * cf), p_sign = mean(f > cf),
                     .groups = "drop") %>%
    dplyr::mutate(class = dplyr::case_when(p_exceed >= robust_exceed & p_sign >= robust_agree ~ "robust",
                                           p_exceed >= robust_exceed ~ "conflicting", TRUE ~ "none"))
}

cat("\n[2] ipcc_scope == s82's pd_scope_class()\n")
set.seed(11)
adm1 <- tibble(adm1_unit = c("KEN_1", "KEN_2", "UGA_1", "UGA_2", "NGA_1", "NGA_2", "GHA_1", "GHA_2"),
               enso_group = rep(c("enso_east", "enso_nga"), each = 4))
months <- seq(as.Date("2026-05-01"), by = "month", length.out = 5)
arms <- tibble(run = c("forecast", "counterfactual", "counterfactual"), analog_set = c(NA, "all", "neutral"))
am <- tidyr::expand_grid(arms, adm1, replicate_id = 1:80, month = months) %>%
  mutate(unit = ifelse(enso_group == "enso_east", "hmmCE", ifelse(substr(adm1_unit, 1, 3) == "GHA", "hmmG", "hmmW")),
         cases = as.numeric(rnbinom(n(), mu = ifelse(run == "forecast", ifelse(enso_group == "enso_east", 70, 38), 35),
                                    size = .4)),
         outbreaks = rbinom(n(), 1, ifelse(run == "forecast", .3, .22)), group_horizon_end = as.Date("2026-09-01"))
am$cases[am$replicate_id == 3 & am$enso_group == "enso_east" & am$run == "forecast"] <- 1e5   # explosive
am$cases[am$replicate_id == 7 & am$enso_group == "enso_nga" & am$analog_set %in% "neutral"] <- 1e5
cm <- am %>% mutate(country = substr(adm1_unit, 1, 3)) %>%
  group_by(run, analog_set, enso_group, country, replicate_id, month, group_horizon_end) %>%
  summarise(n_cases = sum(cases), n_outbreaks = sum(outbreaks), .groups = "drop")
pop <- tibble(adm1_unit = adm1$adm1_unit, pop = c(1e5, 2e5, NA, 3e5, 4e5, 5e5, 6e5, 7e5))
gs <- c(enso_east = "East", enso_nga = "NGA+GHA")
st <- as.Date("2026-06-01"); mags <- c(0, 0.25, 0.5)
# attrib_metrics_windows() hands the core the frames cut to the window.
am_w <- filter(am, month >= st); cm_w <- filter(cm, month >= st)
o <- suppressMessages(attrib_metrics_tables(am_w, cm_w, pop, gs, c("all", "neutral"), st, 50, .05, c(1, 2),
                                            n_boot = 0L, cap = 2e5, ipcc_mags = mags, robust_exceed = .66, robust_agree = .8))
dr <- list(reps = am %>% filter(month >= st) %>% mutate(country = substr(adm1_unit, 1, 3)) %>%
             group_by(enso_group, unit, country, run, analog_set, replicate_id) %>%
             summarise(cases = sum(cases), onsets = sum(outbreaks), .groups = "drop") %>%
             mutate(event = "attribution", .before = 1L),
           explosive = mutate(explosive_replicate_keys(filter(cm, month >= st), st, 2e5), event = "attribution", .before = 1L))
ip <- o$ipcc_scope %>% filter(level %in% c("continent", "region")) %>%
  mutate(key = ifelse(level == "continent", "Africa", enso_group))
n_cmp <- 0L; bad <- 0L
for (metric in c("cases", "onsets")) for (aset in c("all", "neutral")) for (m in mags) {
  ref <- pd_scope_class_ref(dr, metric, aset, m, .66, .8)
  got <- ip %>% filter(metric == !!metric, analog_set == aset, magnitude == m) %>%
    select(key, n_pairs, p_exceed, p_sign, class) %>% arrange(key)
  ref <- ref %>% mutate(n_pairs = as.integer(n_pairs)) %>% arrange(key)
  n_cmp <- n_cmp + 1L
  if (!isTRUE(all.equal(as.data.frame(got), as.data.frame(ref), tolerance = 0))) {
    bad <- bad + 1L; print(got); print(ref)
  }
}
expect(bad == 0L, sprintf("%d (metric x set x m) cells: n_pairs, shares and class identical (%d differ)", n_cmp, bad))
expect(length(unique(ip$class)) >= 2L, sprintf("the synthetic run spans several classes (%s)",
                                              paste(sort(unique(ip$class)), collapse = ", ")))
expect(all(c("country", "region", "continent") %in% o$ipcc_scope$level), "country rows derived too")

cat("\n[3] the v1 tables\n")
o1 <- suppressMessages(attrib_metrics_tables(am_w, cm_w, pop, gs, c("all", "neutral"), st, 50, .05, c(1, 2), n_boot = 0L,
                                             cap = 2e5))
expect(is.null(o1$ipcc_scope), "no ipcc_mags: no ipcc_scope (event cores)")
expect(identical(o1, o[setdiff(names(o), "ipcc_scope")]), "with ipcc_mags: every v1 table unchanged")

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)
