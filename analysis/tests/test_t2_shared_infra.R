# test_t2_shared_infra.R — HARD: shared infrastructure helpers leave existing
#   outputs unchanged.
#
# Usage:
#   Rscript analysis/tests/test_t2_shared_infra.R
#
# Details:
#   Quick-safe: no fit, no data, no network. Checks content_digest() /
#   registry_digest(), smoothing_cache_identity() and its optional-field rule,
#   manifest_add_provenance() and the .TTLO_* constants against frozen
#   reference copies.

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/spatpomp_core.R") # smoothing cache identity (sm_cache_schema_token())
source("analysis/utils/provenance.R") # manifest_add_provenance()
source("analysis/utils/occ_data.R") # make_time_to_last_outbreak_cat()
source("analysis/utils/simulator_core.R") # bin_ttlo(), ttlo_init_from_panel(), .TTLO_BREAKS

err_of <- function(expr) {
  tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

cat("# tier-2 T0 shared infrastructure\n")

# --- [1] content_digest ---------------------------------------------------------
# [1] base.R::content_digest(): canonical text, stability, sensitivity;
#     config.R::registry_digest() (delegating to it) keeps its value for every
#     shipped registry (pinned) and equals the file-bytes digest.
cat("\n[1] content_digest\n")
df <- data.frame(
  a = c(1L, NA),
  b = c("x\"y", NA),
  c = c(0.1, -0),
  d = as.Date(c("2020-01-01", NA)),
  f = factor(c("u", "v")),
  l = c(TRUE, NA),
  stringsAsFactors = FALSE
)
expect(
  identical(
    .content_canonical_text(df),
    paste0(
      "\"a\",\"b\",\"c\",\"d\",\"f\",\"l\"\n",
      "1,\"x\"\"y\",0.1,2020-01-01,\"u\",TRUE\n",
      "NA,NA,0,NA,\"v\",NA\n"
    )
  ),
  "data.frame -> fixed CSV text (quoted strings, NA bare, -0 as 0, Date ISO, factor labels)"
)
d1 <- content_digest(df)
expect(is.character(d1) && nchar(d1) == 16L, "16-hex xxhash64")
expect(
  identical(d1, "ca761b12586fa342"),
  "pinned value of the fixture data.frame"
)
expect(identical(d1, content_digest(df)), "stable across calls")
expect(
  identical(
    d1,
    digest::digest(
      charToRaw(.content_canonical_text(df)),
      algo = "xxhash64",
      serialize = FALSE
    )
  ),
  "= xxhash64 of the canonical text"
)
df_int <- df
df_int$a <- as.double(df_int$a)
expect(
  identical(content_digest(df_int), d1),
  "integer vs double storage: same digest"
)
expect(
  identical(content_digest(tibble::as_tibble(df)), d1),
  "tibble = data.frame"
)
df_attr <- df
attr(df_attr, "stamp") <- "x"
expect(identical(content_digest(df_attr), d1), "extra attributes ignored")
df2 <- df
df2$c[1] <- 0.2
expect(
  !identical(content_digest(df2), d1),
  "a changed value changes the digest"
)
expect(!identical(content_digest(df[2:1, ]), d1), "row order matters")
expect(
  !identical(content_digest(df[, rev(names(df))]), d1),
  "column order matters"
)
df3 <- df
names(df3)[1] <- "A"
expect(
  !identical(content_digest(df3), d1),
  "a renamed column changes the digest"
)
expect(
  !identical(content_digest(c("NA")), content_digest(NA_character_)),
  "\"NA\" vs NA differ"
)
expect(
  identical(content_digest(1 + 1e-17), content_digest(1)),
  "last-bit noise below 15 digits ignored"
)
expect(
  !identical(content_digest(1 + 1e-12), content_digest(1)),
  "a 1e-12 change is seen"
)
expect(
  identical(.content_canonical_text(c(b = 2, a = 1)), "\"b\",2\n\"a\",1\n"),
  "named vector"
)
expect(
  identical(.content_canonical_text(c("a", "b")), "\"a\"\n\"b\"\n"),
  "character vector"
)
expect(
  identical(.content_canonical_text(NULL), "") &&
    identical(.content_canonical_text(character(0)), ""),
  "NULL and zero-length -> empty text"
)
lst <- list(a = 1, 2, b = list(c = "z"), r = as.raw(1:3))
expect(
  identical(
    .content_canonical_text(lst),
    "[[\"a\"]] 2\n1\n[[2]] 2\n2\n[[\"b\"]] 14\n[[\"c\"]] 4\n\"z\"\n[[\"r\"]] 7\n010203\n"
  ),
  "list blocks: key + byte count + element text"
)
expect(
  !identical(
    content_digest(list(a = "x", b = "y")),
    content_digest(list(a = "x\ny"))
  ),
  "list block boundaries are unambiguous"
)
expect(
  identical(
    content_digest(as.raw(c(0x61, 0x62))),
    digest::digest(charToRaw("ab"), algo = "xxhash64", serialize = FALSE)
  ),
  "raw input hashed as its bytes"
)
expect(!is.null(err_of(content_digest(function() 1))), "a function is refused")
expect(
  !is.null(err_of(content_digest(data.frame(x = I(list(1, 2)))))),
  "a list column is refused"
)

# registry_digest: delegated, value unchanged. Pinned values
# (= digest::digest(file =, algo = "xxhash64"), the reference implementation). A
# registry whose CSV is edited later must update its pin here.
REGISTRY_PINS <- c(
  afro_subregion = "62befcaf6c26258d",
  country_adjacency = "7a21a85ae3ec10d8",
  enso_grouping_k2 = "7b03dee1b15652e5",
  enso_grouping_k3 = "6712b170492f7a29",
  enso_grouping_k4 = "5a1a8f14a0d32540",
  enso_grouping_nga = "52157a1207ea7f3b",
  enso_grouping = "e305c000e3c495ac",
  hmm_clusters = "8162e70ff102490e",
  smoke_grouping = "e88d890a62f5eb23"
)
for (s in names(REGISTRY_PINS)) {
  expect(
    identical(registry_digest(s), REGISTRY_PINS[[s]]),
    sprintf("registry_digest(%s) pinned", s)
  )
}
for (s in sub(
  "\\.csv$",
  "",
  list.files("analysis/groupings", pattern = "\\.csv$")
)) {
  expect(
    identical(
      registry_digest(s),
      digest::digest(file = registry_path(s), algo = "xxhash64")
    ),
    sprintf("registry_digest(%s) = file-bytes xxhash64", s)
  )
}
td <- tempfile("t0_reg_")
dir.create(td)
writeLines(c("iso3,group", "AAA,a", "BBB,b"), file.path(td, "fx.csv"))
expect(
  identical(
    registry_digest("fx", groupings_dir = td),
    digest::digest(file = file.path(td, "fx.csv"), algo = "xxhash64")
  ),
  "fixture registry = file-bytes xxhash64"
)
invisible(file.create(file.path(td, "empty.csv")))
expect(
  identical(
    registry_digest("empty", groupings_dir = td),
    digest::digest(file = file.path(td, "empty.csv"), algo = "xxhash64")
  ),
  "empty registry file = file-bytes xxhash64"
)
expect(
  !is.null(err_of(registry_digest("absent", groupings_dir = td))),
  "missing registry refused"
)
unlink(td, recursive = TRUE)

# --- [2] smoothing cache identity -------------------------------------------------
# [2] spatpomp_core.R::smoothing_cache_identity() = the reference 08c meta list
#     (frozen copy) and the same meta.yml bytes; the optional-field rule
#     (.cache_meta_add_optional() / .cache_meta_optional_check()); the 08c
#     construction appends `software`, which no reader checks: a cache built
#     under other versions loads.
cat("\n[2] smoothing_cache_identity / .cache_meta_optional_check\n")
# Frozen reference copy of the 08c construction (08c_smoothing_compile.R).
old_08c_meta <- function(S, by_tu) {
  meta <- smoothing_cache_meta(
    S$PARTITION_BY %||% "transmission_unit",
    S$PARTITION_SOURCE %||% "",
    grouping_source = if (isTRUE(S$GROUPED)) S$GRP_SOURCE else "",
    fits_dir = if (isTRUE(S$GROUPED)) S$GRP_FITS_DIR else "",
    Np = S$NP,
    nrep = S$N_REPS,
    schema = sm_cache_schema_token(S$COVARIATE_FEED %||% "latent")
  )
  meta$run_tag <- S$run_tag
  meta$tus <- names(by_tu)
  meta
}
S_fixtures <- list(
  default = list(
    N_REPS = 100L,
    NP = 25000L,
    COVARIATE_FEED = "latent",
    PARTITION_BY = "transmission_unit",
    PARTITION_SOURCE = "",
    run_tag = "s10_historical_smoothing",
    GROUPED = FALSE,
    GRP_SOURCE = NULL,
    GRP_FITS_DIR = NULL
  ),
  registry = list(
    N_REPS = 40L,
    NP = 5000L,
    COVARIATE_FEED = "detected",
    PARTITION_BY = "registry",
    PARTITION_SOURCE = "smoke_grouping",
    run_tag = "smoke_hist",
    GROUPED = TRUE,
    GRP_SOURCE = "enso_grouping",
    GRP_FITS_DIR = "generated_data/enso_best_lag_fits/x"
  ),
  sparse = list(
    PARTITION_BY = "country",
    GROUPED = TRUE,
    GRP_SOURCE = "hmm_clusters",
    GRP_FITS_DIR = "fits"
  )
)
by_tu_fx <- list(tu_b = 1, tu_a = 2)
for (nm in names(S_fixtures)) {
  S <- S_fixtures[[nm]]
  old <- old_08c_meta(S, by_tu_fx)
  new <- smoothing_cache_identity(S, tus = names(by_tu_fx))
  expect(
    identical(new, old),
    sprintf("%s: identity list identical to the pre-T0 construction", nm)
  )
  f_old <- tempfile(fileext = ".yml")
  f_new <- tempfile(fileext = ".yml")
  yaml::write_yaml(old, f_old)
  yaml::write_yaml(new, f_new)
  expect(
    identical(
      readBin(f_old, "raw", file.size(f_old)),
      readBin(f_new, "raw", file.size(f_new))
    ),
    sprintf("%s: meta.yml bytes identical", nm)
  )
  unlink(c(f_old, f_new))
}
# Registered optional fields (each with its legacy value; a default run writes none):
# pop_year.
expect(
  identical(.SM_META_OPTIONAL_FIELDS[["pop_year"]], 2020L),
  "registered optional fields carry their legacy value (pop_year = 2020L, item 9)"
)
expect(
  identical(
    smoothing_cache_identity(
      S_fixtures$default,
      tus = "x",
      optional = smoothing_cache_optional_fields(S_fixtures$default)
    ),
    smoothing_cache_identity(S_fixtures$default, tus = "x")
  ),
  "default smoothing: the optional fields add nothing (meta.yml unchanged)"
)
expect(
  identical(
    smoothing_cache_identity(S_fixtures$default, tus = "x", optional = list()),
    smoothing_cache_identity(S_fixtures$default, tus = "x")
  ),
  "empty optional = no change"
)
expect(
  !is.null(err_of(smoothing_cache_identity(
    S_fixtures$default,
    optional = list(not_registered = "x")
  ))),
  "an unregistered optional field is refused"
)
# Frozen copy of the 08c construction with software (08c_smoothing_compile.R:
# optional fields, then software = sim_software_versions()): the reference list,
# the non-legacy optional fields, `software` last.
t3b_08c_meta <- function(S, by_tu, software) {
  c(
    .cache_meta_add_optional(
      old_08c_meta(S, by_tu),
      smoothing_cache_optional_fields(S)
    ),
    list(software = software)
  )
}
sw_fx <- list(
  R = "4.3.2",
  platform = "x86_64-pc-linux-gnu",
  spatPomp = "0.36.0",
  pomp = "5.11"
)
for (nm in names(S_fixtures)) {
  S <- S_fixtures[[nm]]
  expect(
    identical(
      smoothing_cache_identity(
        S,
        tus = names(by_tu_fx),
        optional = smoothing_cache_optional_fields(S),
        software = sw_fx
      ),
      t3b_08c_meta(S, by_tu_fx, sw_fx)
    ),
    sprintf(
      "%s: the 08c identity list = the pre-T0 list + optional fields + software (last)",
      nm
    )
  )
}
expect(
  !"software" %in% c(names(.SM_META_OPTIONAL_FIELDS), .SM_META_RECORD_ONLY),
  "software is registered in neither .SM_META_OPTIONAL_FIELDS nor .SM_META_RECORD_ONLY"
)
td_sw <- tempfile("t3b_sw_")
dir.create(td_sw)
expect_meta <- smoothing_cache_identity(S_fixtures$registry, tus = "x")
yaml::write_yaml(
  smoothing_cache_identity(
    S_fixtures$registry,
    tus = "x",
    software = utils::modifyList(sw_fx, list(spatPomp = "9.9"))
  ),
  smoothing_cache_meta_path(td_sw)
)
expect(
  is.null(err_of(check_smoothing_cache_meta(td_sw, expect_meta))),
  "a cache whose meta.yml records other software versions still passes check_smoothing_cache_meta()"
)
unlink(td_sw, recursive = TRUE)
# The rule on a fixture registry (pop_year).
REG <- list(pop_year = "pop_ref_year")
m0 <- list(partition_by = "country")
expect(
  identical(
    .cache_meta_add_optional(m0, list(pop_year = "pop_ref_year"), REG),
    m0
  ),
  "writer: the legacy value is not written"
)
expect(
  identical(
    .cache_meta_add_optional(m0, list(pop_year = "season_year"), REG)$pop_year,
    "season_year"
  ),
  "writer: a non-legacy value is written"
)
expect(
  isTRUE(.cache_meta_optional_check(
    m0,
    "pop_year",
    "pop_ref_year",
    "pop_ref_year"
  )),
  "reader: absent field = legacy; legacy run accepted"
)
e <- err_of(.cache_meta_optional_check(
  m0,
  "pop_year",
  "season_year",
  "pop_ref_year"
))
expect(
  !is.null(e) &&
    grepl("pop_year", e) &&
    grepl("pop_ref_year", e) &&
    grepl("season_year", e),
  "reader: absent (legacy) cache vs non-legacy run refused, naming the field and both values"
)
m1 <- list(partition_by = "country", pop_year = "season_year")
expect(
  isTRUE(.cache_meta_optional_check(
    m1,
    "pop_year",
    "season_year",
    "pop_ref_year"
  )),
  "reader: match accepted"
)
e <- err_of(.cache_meta_optional_check(
  m1,
  "pop_year",
  "pop_ref_year",
  "pop_ref_year"
))
expect(
  !is.null(e) && grepl("'season_year' vs run 'pop_ref_year'", e, fixed = TRUE),
  "reader: mismatch refused, naming both values"
)
expect(
  isTRUE(.cache_meta_optional_check(list(n = 2), "n", 2L, 0L)),
  "reader: 2 (YAML double) = 2L"
)
# check_smoothing_cache_meta() is unchanged with no registered optional field.
dm <- tempfile("t0_cache_")
dir.create(dm)
exp_c <- smoothing_cache_meta("country")
yaml::write_yaml(
  c(exp_c, list(some_future_field = "x")),
  smoothing_cache_meta_path(dm)
)
expect(
  is.null(err_of(check_smoothing_cache_meta(dm, exp_c))),
  "an unregistered extra field in meta.yml is ignored (as before T0)"
)
unlink(dm, recursive = TRUE)

# --- [3] manifest_add_provenance --------------------------------------------------
# [3] provenance.R::manifest_add_provenance(): only a non-default value is
#     added, so a default manifest's YAML is byte-identical.
cat("\n[3] manifest_add_provenance\n")
man <- list(
  generated_at = "t",
  grouping = "enso_grouping",
  groups = list(a = 1)
)
expect(
  identical(
    manifest_add_provenance(man, "pop_year", "pop_ref_year", "pop_ref_year"),
    man
  ),
  "default value: manifest unchanged"
)
m2 <- manifest_add_provenance(man, "pop_year", "season_year", "pop_ref_year")
expect(
  identical(m2, c(man, list(pop_year = "season_year"))),
  "non-default value: appended last"
)
expect(
  identical(manifest_add_provenance(man, "x", NULL, "d"), man),
  "NULL value: nothing added"
)
expect(
  identical(manifest_add_provenance(man, "grouping", NULL, "d"), man),
  "NULL value: nothing removed"
)
expect(
  identical(
    manifest_add_provenance(
      man,
      "f",
      list(method = "neutral"),
      list(method = "neutral")
    ),
    man
  ),
  "list default: unchanged"
)
expect(
  identical(
    manifest_add_provenance(
      man,
      "f",
      list(method = "analog"),
      list(method = "neutral")
    )$f,
    list(method = "analog")
  ),
  "list non-default: added"
)
f_a <- tempfile(fileext = ".yml")
f_b <- tempfile(fileext = ".yml")
yaml::write_yaml(man, f_a)
yaml::write_yaml(
  manifest_add_provenance(man, "pop_year", "pop_ref_year", "pop_ref_year"),
  f_b
)
expect(
  identical(
    readBin(f_a, "raw", file.size(f_a)),
    readBin(f_b, "raw", file.size(f_b))
  ),
  "default value: manifest YAML bytes identical"
)
unlink(c(f_a, f_b))

# --- [4] TTLO constants -------------------------------------------------------------
# [4] base.R .TTLO_ORIGIN / .TTLO_LEVELS / .TTLO_NONE_DAYS, and the sites
#     rewired to them (occ_data.R::make_time_to_last_outbreak_cat(),
#     03_prepare_occ.R's origin, simulator_core.R bin_ttlo() /
#     ttlo_init_from_panel()) identical to frozen reference copies.
cat("\n[4] .TTLO_* constants\n")
expect(
  identical(.TTLO_ORIGIN, as.Date("2000-01-01")),
  ".TTLO_ORIGIN = as.Date(\"2000-01-01\")"
)
expect(
  identical(.TTLO_NONE_DAYS, 1e9),
  ".TTLO_NONE_DAYS = 1e9 (reserved for item 11)"
)
expect(
  identical(
    .TTLO_LEVELS,
    levels(cut(
      0,
      c(0, 1, 2, 3, 4, 5, 6, Inf),
      include.lowest = TRUE,
      right = FALSE
    ))
  ),
  ".TTLO_LEVELS = the labels cut() derives from the TTLO breaks"
)
expect(
  identical(.TTLO_BREAKS, c(0, 1, 2, 3, 4, 5, 6, Inf)),
  "simulator_core.R .TTLO_BREAKS unchanged"
)
# Frozen reference copies of the functions.
old_makeTimeToLastOutbreakCat <- function(df) {
  df %>%
    mutate(
      time_to_last_outbreak = cut(
        (pmax(0, time_to_last_outbreak) + 1) / 365,
        c(0, 1, 2, 3, 4, 5, 6, Inf),
        include.lowest = T,
        right = F
      )
    ) %>%
    tidyr::replace_na(list(time_to_last_outbreak = "[6,Inf]"))
}
old_bin_ttlo <- function(days) {
  years <- (pmax(0, days) + 1) / 365
  as.character(cut(
    years,
    c(0, 1, 2, 3, 4, 5, 6, Inf),
    include.lowest = TRUE,
    right = FALSE
  ))
}
old_ttlo_init_from_panel <- function(ttlo_cat_vec) {
  lo_year <- c(
    "[0,1)" = 0,
    "[1,2)" = 1,
    "[2,3)" = 2,
    "[3,4)" = 3,
    "[4,5)" = 4,
    "[5,6)" = 5,
    "[6,Inf]" = 6
  )
  chr <- as.character(ttlo_cat_vec)
  days <- lo_year[chr] * 365.25
  names(days) <- NULL
  days[is.na(days)] <- 6 * 365.25
  days
}
days_fx <- c(
  NA,
  -400,
  -1,
  0,
  1,
  363,
  364,
  365,
  729,
  730,
  1094,
  1095,
  1500,
  2000,
  2189,
  2190,
  2191,
  5000,
  1e9,
  Inf,
  seq(0, 3000, by = 0.5)
)
ttlo_df <- tibble::tibble(
  uid = seq_along(days_fx),
  time_to_last_outbreak = days_fx
)
new_cat <- make_time_to_last_outbreak_cat(ttlo_df)
expect(
  identical(new_cat, old_makeTimeToLastOutbreakCat(ttlo_df)),
  "make_time_to_last_outbreak_cat() identical to the pre-T0 function (factor, levels, NA level)"
)
expect(
  identical(levels(new_cat$time_to_last_outbreak), .TTLO_LEVELS),
  "panel factor levels = .TTLO_LEVELS"
)
expect(
  identical(bin_ttlo(days_fx), old_bin_ttlo(days_fx)),
  "bin_ttlo() identical to the pre-T0 function"
)
expect(
  identical(bin_ttlo(.TTLO_NONE_DAYS), "[6,Inf]"),
  "the reserved sentinel still bins to [6,Inf] at T0"
)
# "none" is not an unknown label (it maps to the sentinel), so the identity
# fixture omits it.
labs_fx <- c(.TTLO_LEVELS, NA, "[7,8)", rev(.TTLO_LEVELS))
expect(
  identical(ttlo_init_from_panel(labs_fx), old_ttlo_init_from_panel(labs_fx)),
  "ttlo_init_from_panel() identical to the pre-T0 function (characters)"
)
expect(
  identical(ttlo_init_from_panel("none"), .TTLO_NONE_DAYS),
  "item 11: the label none maps to the sentinel"
)
expect(
  identical(
    ttlo_init_from_panel(new_cat$time_to_last_outbreak),
    old_ttlo_init_from_panel(new_cat$time_to_last_outbreak)
  ),
  "ttlo_init_from_panel() identical to the pre-T0 function (panel factor)"
)
# 03_prepare_occ.R's TTLO block with the literal origin vs .TTLO_ORIGIN.
ob_fx <- tibble::tibble(
  location_period_id = rep(c("lp1", "lp2"), each = 5),
  t_start = rep(
    as.Date(c(
      "2010-01-01",
      "2010-02-01",
      "2011-03-01",
      "2012-01-01",
      "2015-06-01"
    )),
    2
  ),
  t_end = rep(
    as.Date(c(
      "2010-01-31",
      "2010-02-28",
      "2011-03-31",
      "2012-01-31",
      "2015-06-30"
    )),
    2
  ),
  outbreak = c(
    FALSE,
    TRUE,
    FALSE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE
  )
)
ttlo_block <- function(d, origin) {
  d %>%
    group_by(location_period_id) %>%
    mutate(
      rid = row_number(),
      last_outbreak_start = origin,
      last_outbreak_start = case_when(
        outbreak ~ t_end,
        TRUE ~ last_outbreak_start
      ),
      last_outbreak_start = cummax(as.numeric(last_outbreak_start)) %>%
        as.Date(origin = "1970-01-01") %>%
        lag(n = 1),
      time_to_last_outbreak = difftime(
        t_start,
        last_outbreak_start,
        units = "days"
      ) %>%
        as.numeric()
    ) %>%
    ungroup()
}
expect(
  identical(
    ttlo_block(ob_fx, .TTLO_ORIGIN),
    ttlo_block(ob_fx, as.Date("2000-01-01"))
  ),
  "03_prepare_occ.R TTLO block: .TTLO_ORIGIN = the literal origin"
)

cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
