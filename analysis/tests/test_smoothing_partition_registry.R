# test_smoothing_partition_registry.R — HARD: filter.partition_by: registry.
#
# Usage:
#   Rscript analysis/tests/test_smoothing_partition_registry.R
#
# Details:
#   No fit, no data. The partition (TU) registry is decoupled from the fit
#   grouping: smoke_grouping is the partition over units grouped by
#   enso_grouping, so TU grp_a = ETH (enso_horn) + KEN (enso_east) MIXES two
#   fit groups.

suppressPackageStartupMessages({
  library(qs2)
  library(lubridate)
})
source("analysis/utils/simulator_core.R")
source("analysis/utils/paths.R")
source("analysis/utils/spatpomp_core.R")
source("analysis/utils/smoothing_core.R")

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
err_of <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
warn_of <- function(expr) {
  w <- NULL
  withCallingHandlers(force(expr), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  w
}

# Four countries; fit groups from enso_grouping; partition from smoke_grouping.
frame4 <- function() {
  data.frame(
    tu = c(1L, 1L, 2L, 3L, 3L),
    country = c("ETH", "KEN", "KEN", "NGA", "CMR"),
    group = c("enso_horn", "enso_east", "enso_east", "enso_west", "enso_west"),
    adm1_unit = c("ETH.1", "KEN.1", "KEN.2", "NGA.1", "CMR.1"),
    stringsAsFactors = FALSE
  )
}

# --- [1] helpers -------------------------------------------------------------
# [1] helpers: partition_frame_column("registry"), load_partition_map(),
#     stamp_partition_column() (coverage error, unsafe labels).
cat("\n[1] registry partition helpers\n")
expect(
  identical(partition_frame_column("registry"), "part"),
  "partition_frame_column(registry) == 'part'"
)
pm <- load_partition_map("smoke_grouping")
expect(
  identical(
    unname(pm[c("ETH", "KEN", "NGA", "CMR")]),
    c("grp_a", "grp_a", "grp_b", "grp_b")
  ),
  "load_partition_map(smoke_grouping) maps ISO3 -> grp_a/grp_b"
)
expect(
  grepl("partition_source", err_of(load_partition_map("")), fixed = TRUE),
  "empty partition_source refused, naming the key"
)
fr <- stamp_partition_column(frame4(), "registry", "smoke_grouping")
expect(
  identical(fr$part, c("grp_a", "grp_a", "grp_a", "grp_b", "grp_b")),
  "stamp_partition_column writes part"
)
expect(
  identical(stamp_partition_column(frame4(), "enso_region"), frame4()),
  "stamp_partition_column is a no-op for other partition_by values"
)
bad <- rbind(
  frame4(),
  data.frame(
    tu = 4L,
    country = "ZMB",
    group = "enso_south",
    adm1_unit = "ZMB.1"
  )
)
e <- err_of(stamp_partition_column(bad, "registry", "smoke_grouping"))
expect(
  !is.null(e) && grepl("ZMB", e),
  "country absent from the partition registry refused, named"
)
td <- tempfile("grp")
dir.create(td)
writeLines(c("iso3,group", "ETH,a b"), file.path(td, "unsafe.csv"))
expect(
  !is.null(err_of(load_partition_map("unsafe", groupings_dir = td))),
  "path-unsafe TU label refused"
)
e <- err_of(assert_partition_column(data.frame(tu = 1L), "registry"))
expect(
  !is.null(e) && grepl("partition_source", e),
  "assert_partition_column(registry) hint names partition_source"
)

# --- [2] smoothing_apply_partition registry branch ------------------------------
# [2] smoothing_apply_partition() registry branch: part from the registry, tu
#     and the fit `group` left intact; partitions enumerate registry groups and
#     grp_a's frame holds both fit groups.
cat("\n[2] smoothing_apply_partition(registry)\n")
out <- smoothing_apply_partition(
  list(frame = frame4()),
  list(PARTITION_BY = "registry", PARTITION_SOURCE = "smoke_grouping")
)
expect(
  identical(out$frame$part, c("grp_a", "grp_a", "grp_a", "grp_b", "grp_b")),
  "part from the registry"
)
expect(identical(out$frame$tu, frame4()$tu), "tu left intact (no collapse)")
expect(identical(out$frame$group, frame4()$group), "fit group left intact")
expect(is.null(out$frame$tu_true), "no tu_true stamped")
expect(
  identical(smoothing_tu_list(out$frame), c("grp_a", "grp_b")),
  "partitions = registry groups"
)
fa <- out$frame[out$frame$part == "grp_a", ]
expect(
  setequal(fa$group, c("enso_horn", "enso_east")),
  "grp_a mixes two fit groups"
)
expect(setequal(fa$country, c("ETH", "KEN")), "grp_a holds ETH and KEN units")

# --- [3] exact cache resolver ---------------------------------------------------
# [3] exact cache resolver: TU grp_a never matches grp_a_x's cache; legacy
#     name without _Np= is a warned fallback.
cat("\n[3] exact cache resolver\n")
svtok <- sm_cache_schema_token()
dc <- tempfile("sm_")
dir.create(dc)
fn <- function(tu, np = "_Np=100") {
  sprintf("hist_filter_tu=%s%s_%s_nrep=2_20090101_20110101.qs", tu, np, svtok)
}
invisible(file.create(file.path(dc, fn("grp_a_x"))))
expect(
  length(.smoothing_cache_hits("grp_a", dc)) == 0L,
  "grp_a does not match grp_a_x's cache"
)
invisible(file.create(file.path(dc, fn("grp_a"))))
h <- .smoothing_cache_hits("grp_a", dc)
expect(
  length(h) == 1L && basename(h) == fn("grp_a"),
  "grp_a resolves to exactly its own cache"
)
dl <- tempfile("sm_")
dir.create(dl)
invisible(file.create(file.path(dl, "hist_filter_tu=grp_a_legacy.qs")))
w <- warn_of(h2 <- .smoothing_cache_hits("grp_a", dl))
expect(
  length(h2) == 1L && any(grepl("legacy", w)),
  "legacy name without _Np= found with a warning"
)

# --- [4] cache-dir identity -------------------------------------------------------
# [4] cache-dir identity (meta.yml): missing warns, mismatching partition
#     source/digest refused, per-country caches accepted as stitch source.
cat("\n[4] cache-dir meta.yml\n")
exp_reg <- smoothing_cache_meta("registry", "smoke_grouping", "enso_grouping")
expect(
  identical(exp_reg$partition_digest, registry_digest("smoke_grouping")),
  "meta carries the partition digest"
)
expect(
  identical(
    smoothing_cache_meta("country", "smoke_grouping")$partition_source,
    ""
  ),
  "partition_source blanked unless partition_by = registry"
)
dm <- tempfile("meta_")
dir.create(dm)
w <- warn_of(check_smoothing_cache_meta(dm, exp_reg))
expect(any(grepl("no meta.yml", w)), "missing meta.yml warns")
yaml::write_yaml(exp_reg, smoothing_cache_meta_path(dm))
expect(
  is.null(err_of(check_smoothing_cache_meta(dm, exp_reg))),
  "matching identity accepted"
)
other <- exp_reg
other$partition_digest <- "0000000000000000"
yaml::write_yaml(other, smoothing_cache_meta_path(dm))
e <- err_of(check_smoothing_cache_meta(dm, exp_reg))
expect(
  !is.null(e) && grepl("partition_digest", e),
  "registry edited since smoothing (digest) refused"
)
other <- exp_reg
other$partition_source <- "hmm_clusters"
yaml::write_yaml(other, smoothing_cache_meta_path(dm))
expect(
  !is.null(err_of(check_smoothing_cache_meta(dm, exp_reg))),
  "other partition registry refused"
)
yaml::write_yaml(smoothing_cache_meta("country"), smoothing_cache_meta_path(dm))
expect(
  is.null(err_of(check_smoothing_cache_meta(dm, exp_reg))),
  "per-country caches accepted (stitch source)"
)
yaml::write_yaml(
  smoothing_cache_meta("enso_region", grouping_source = "enso_grouping"),
  smoothing_cache_meta_path(dm)
)
expect(
  is.null(err_of(check_smoothing_cache_meta(dm, exp_reg))),
  "enso_region caches accepted as stitch source for a registry run (gid check guards content)"
)
yaml::write_yaml(
  smoothing_cache_meta("transmission_unit"),
  smoothing_cache_meta_path(dm)
)
expect(
  !is.null(err_of(check_smoothing_cache_meta(dm, exp_reg))),
  "transmission_unit caches refused for a registry run"
)
expect(
  is.null(err_of(check_smoothing_cache_meta(dm, NULL))),
  "NULL expectation = no check"
)

# --- [5] load_smoothing_distribution(registry) -------------------------------------
# [5] load_smoothing_distribution(partition_by = "registry"): native per-TU
#     cache loaded directly; otherwise assembled from member caches; members
#     holding different replicate sets refused.
cat("\n[5] load_smoothing_distribution(partition_by = registry)\n")
.sn <- .SP_UNIT_STATENAMES
fm <- as.Date(c("2015-02-01", "2015-03-01"))
mk_ft <- function(U, off, rep) {
  rn <- unlist(lapply(.sn, function(s) paste0(s, seq_len(U))))
  ft <- array(0, dim = c(length(rn), 1L, 2L), dimnames = list(rn, NULL, NULL))
  for (lp in seq_len(U)) {
    ft[paste0("o_tu", lp), 1L, ] <- .SP_OTU_SENTINEL
    ft[paste0("o", lp), 1L, 2L] <- 1000 * rep + off + lp
  }
  ft
}
mk_rec <- function(tu, gids, off, rep, ll) {
  list(
    tu = tu,
    rep_id = rep,
    draw_d = 3L,
    filter_traj = mk_ft(length(gids), off, rep),
    gids = gids,
    filt_months = fm,
    logLik = ll
  )
}
save_cache <- function(dir, tu, recs) qs2::qs_save(recs, file.path(dir, fn(tu)))
fw <- list(start = "2015-04-01", end = "2016-01-01")

# (a) stitched from member caches
ds <- tempfile("st_")
dir.create(ds)
save_cache(
  ds,
  "ETH",
  list(mk_rec("ETH", "ETH_1", 0, 1, -1), mk_rec("ETH", "ETH_1", 0, 2, -2))
)
save_cache(
  ds,
  "KEN",
  list(
    mk_rec("KEN", c("KEN_1", "KEN_2"), 10, 2, -20),
    mk_rec("KEN", c("KEN_1", "KEN_2"), 10, 1, -10)
  )
)
sm <- suppressWarnings(load_smoothing_distribution(
  "grp_a",
  ds,
  fw,
  partition_by = "registry",
  region_members = list(grp_a = c("ETH", "KEN"))
))
r1 <- sm$filter_results_by_tu$grp_a[[1L]]
expect(length(sm$filter_results_by_tu$grp_a) == 2L, "stitched: two replicates")
expect(
  identical(r1$gids, c("ETH_1", "KEN_1", "KEN_2")),
  "stitched: members' gids concatenated"
)
expect(
  isTRUE(all.equal(r1$logLik, -11)),
  "stitched: logLik summed over members (rep 1)"
)
expect(
  identical(sm$handoff_by_tu$grp_a, as.Date("2015-04-01")),
  "stitched: handoff = last filt month + 1"
)

# (b) members with different replicate sets refused
dx <- tempfile("st_")
dir.create(dx)
save_cache(
  dx,
  "ETH",
  list(mk_rec("ETH", "ETH_1", 0, 1, -1), mk_rec("ETH", "ETH_1", 0, 3, -3))
)
save_cache(
  dx,
  "KEN",
  list(mk_rec("KEN", "KEN_1", 10, 1, -10), mk_rec("KEN", "KEN_1", 10, 2, -20))
)
e <- err_of(suppressWarnings(load_smoothing_distribution(
  "grp_a",
  dx,
  fw,
  partition_by = "registry",
  region_members = list(grp_a = c("ETH", "KEN"))
)))
expect(
  !is.null(e) && grepl("different replicate sets", e),
  "different member replicate sets refused"
)

# (c) native per-TU cache preferred
dn <- tempfile("nat_")
dir.create(dn)
save_cache(dn, "grp_a", list(mk_rec("grp_a", c("ETH_1", "KEN_1"), 0, 1, -5)))
save_cache(dn, "grp_a_x", list(mk_rec("grp_a_x", "X_1", 0, 1, -99)))
yaml::write_yaml(exp_reg, smoothing_cache_meta_path(dn))
sm <- load_smoothing_distribution(
  "grp_a",
  dn,
  fw,
  partition_by = "registry",
  region_members = list(grp_a = c("ETH", "KEN")),
  expect_meta = exp_reg
)
expect(
  length(sm$filter_results_by_tu$grp_a) == 1L &&
    identical(sm$filter_results_by_tu$grp_a[[1L]]$gids, c("ETH_1", "KEN_1")),
  "native per-TU cache loaded directly (not grp_a_x's)"
)
bad_exp <- exp_reg
bad_exp$partition_source <- "hmm_clusters"
expect(
  !is.null(err_of(load_smoothing_distribution(
    "grp_a",
    dn,
    fw,
    partition_by = "registry",
    region_members = list(grp_a = c("ETH", "KEN")),
    expect_meta = bad_exp
  ))),
  "load refuses a cache dir of another partition identity"
)

# (d) unit = union of whole fit groups, assembled from fit-group caches
dg <- tempfile("grpc_")
dir.create(dg)
save_cache(dg, "enso_horn", list(mk_rec("enso_horn", "ETH_1", 0, 1, -1)))
save_cache(
  dg,
  "enso_east",
  list(mk_rec("enso_east", c("KEN_1", "KEN_2"), 10, 1, -10))
)
sm <- suppressWarnings(load_smoothing_distribution(
  "grp_a",
  dg,
  fw,
  partition_by = "registry",
  region_members = list(grp_a = c("ETH", "KEN")),
  region_pieces = list(grp_a = list(groups = c("enso_east", "enso_horn"))),
  tu_gids = list(grp_a = c("ETH_1", "KEN_1", "KEN_2"))
))
expect(
  setequal(
    sm$filter_results_by_tu$grp_a[[1L]]$gids,
    c("ETH_1", "KEN_1", "KEN_2")
  ),
  "unit assembled from its fit groups' caches when per-country caches are absent"
)

# (e) gid coverage check: records that do not tile the TU are refused
e <- err_of(suppressWarnings(load_smoothing_distribution(
  "grp_a",
  dn,
  fw,
  partition_by = "registry",
  region_members = list(grp_a = c("ETH", "KEN")),
  tu_gids = list(grp_a = c("ETH_1", "KEN_1", "KEN_2"))
)))
expect(
  !is.null(e) && grepl("another partition", e),
  "native cache missing a TU unit refused (gid check)"
)

# --- [6] pre-flight -----------------------------------------------------------------
# [6] check_partition_vs_grouping(): contagion blocks spanning TUs refused
#     under registry, warned under country; missing fit group refused.
cat("\n[6] check_partition_vs_grouping\n")
fr6 <- out$frame
fr6$tu <- fr6$part
nbr_enso <- list(use = TRUE, registry = load_grouping("enso_grouping"))
# enso blocks: enso_east = KEN only, enso_horn = ETH, enso_west = NGA+CMR -> all inside one TU
expect(
  is.null(err_of(check_partition_vs_grouping(fr6, "registry", NULL, nbr_enso))),
  "contagion blocks inside one TU each: accepted"
)
# a contagion registry that couples KEN (grp_a) with NGA (grp_b) spans TUs
span_reg <- data.frame(
  iso3 = c("ETH", "KEN", "NGA", "CMR"),
  group = c("b1", "b2", "b2", "b3")
)
e <- err_of(check_partition_vs_grouping(
  fr6,
  "registry",
  NULL,
  list(use = TRUE, registry = span_reg)
))
expect(
  !is.null(e) && grepl("b2", e) && grepl("span several TUs", e),
  "spanning contagion block refused under registry"
)
frc <- frame4()
frc$tu <- frc$country
w <- warn_of(check_partition_vs_grouping(frc, "country", NULL, nbr_enso))
expect(
  any(grepl("enso_west", w)),
  "spanning block only warned under partition_by = country"
)
expect(
  is.null(err_of(check_partition_vs_grouping(
    fr6,
    "registry",
    NULL,
    list(use = FALSE)
  ))),
  "contagion off: no refinement check"
)
frg <- fr6
frg$group[1L] <- NA
grp <- list(
  draws_occ_by_group = list(
    enso_horn = list(S = 100L),
    enso_east = list(S = 80L),
    enso_west = list(S = 90L)
  )
)
e <- err_of(check_partition_vs_grouping(frg, "registry", grp, NULL))
expect(!is.null(e) && grepl("ETH", e), "country without a fit group refused")
res <- suppressMessages(check_partition_vs_grouping(fr6, "registry", grp, NULL))
expect(
  identical(res$draw_ceiling[["grp_a"]], 80L),
  "mixed TU draw ceiling = min S over its groups"
)

unlink(c(td, dc, dl, dm, ds, dx, dn, dg), recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
