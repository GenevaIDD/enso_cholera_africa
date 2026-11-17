# test_grouping_meta.R — HARD: grouping-registry sidecars, digests and label
#   helpers (utils/config.R).
#
# Usage:
#   Rscript analysis/tests/test_grouping_meta.R
#
# Details:
#   Checks the enso_grouping sidecar [1], the fallback without a sidecar [2],
#   unknown groups [3], registry_digest() [4], resolve_group_label() [5], the
#   label tokens [6] and assert_partition_labels() [7]. No data, no fit.

source("analysis/utils/base.R")
source("analysis/utils/config.R")
source("analysis/utils/enso.R")

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

cat("# grouping registry sidecars, digests, labels\n")

# --- [1] enso_grouping sidecar == the formerly hardcoded tables -------------
# The enso_grouping sidecar reproduces, value for value and name for name, the
# tables formerly hardcoded in utils/results_core.R (GROUP_LABELS /
# GROUP_SHORT / GROUP_ORDER / GROUP_SHORT_FIG) and utils/enso.R
# (enso_group_colours()). Phase-09 artefacts store GROUP_SHORT as data, so a
# drift here silently orphans stored scopes.
cat("\n[1] enso_grouping sidecar reproduces the hardcoded tables\n")
OLD_LABELS <- c(
  enso_east = "East Africa",
  enso_horn = "Horn of Africa",
  enso_west = "West Africa",
  enso_south = "Southern Africa",
  enso_central = "Central Africa"
)
OLD_SHORT <- c(
  enso_east = "East",
  enso_horn = "Horn",
  enso_west = "West",
  enso_south = "Southern",
  enso_central = "Central"
)
OLD_ORDER <- c(
  "enso_east",
  "enso_horn",
  "enso_west",
  "enso_south",
  "enso_central"
)
OLD_SHORT_FIG <- replace(OLD_SHORT, "enso_south", "South")
OLD_COLOURS <- c(
  enso_east = "#1b9e77",
  enso_horn = "#d95f02",
  enso_west = "#66a61e",
  enso_south = "#7570b3",
  enso_central = "#e7298a"
)
m <- load_grouping_meta("enso_grouping")
expect(identical(m$label, OLD_LABELS), "label == GROUP_LABELS")
expect(identical(m$short, OLD_SHORT), "short == GROUP_SHORT")
expect(identical(m$order, OLD_ORDER), "order == GROUP_ORDER")
expect(identical(m$short_fig, OLD_SHORT_FIG), "short_fig == GROUP_SHORT_FIG")
expect(identical(m$colour, OLD_COLOURS), "colour == old enso_group_colours()")
expect(
  identical(enso_group_colours(), OLD_COLOURS),
  "enso_group_colours() unchanged"
)
expect(identical(m$label_prefix, "enso_"), "label_prefix enso_")

# results_core.R globals (sourced in a child env so its side effects stay local)
rc <- new.env()
ok <- is.null(err_of(sys.source(
  "analysis/utils/results_core.R",
  envir = rc,
  chdir = FALSE
)))
expect(ok, "results_core.R sources")
if (ok) {
  expect(
    identical(rc$GROUP_LABELS, OLD_LABELS) &&
      identical(rc$GROUP_SHORT, OLD_SHORT) &&
      identical(rc$GROUP_ORDER, OLD_ORDER) &&
      identical(rc$GROUP_SHORT_FIG, OLD_SHORT_FIG) &&
      identical(rc$GROUP_COLOURS, OLD_COLOURS),
    "results_core GROUP_* globals identical to the old literals"
  )
  expect(
    identical(
      rc$results_group_meta(list(
        MAN = list(grouping = "smoke_grouping")
      ))$source,
      "smoke_grouping"
    ),
    "results_group_meta() follows MAN$grouping"
  )
  expect(
    identical(rc$results_group_meta(NULL)$source, "enso_grouping"),
    "results_group_meta() defaults to enso_grouping"
  )
}

# --- [2] fallback without sidecar ---------------------------------------------
# A registry without a sidecar falls back to its own labels.
cat("\n[2] fallback without sidecar\n")
td <- tempfile("grp")
dir.create(td)
writeLines(
  c("iso3,group", "AAA,zz", "BBB,yy", "CCC,zz"),
  file.path(td, "nosc.csv")
)
f <- load_grouping_meta("nosc", groupings_dir = td)
expect(identical(f$order, c("yy", "zz")), "order = sorted labels")
expect(
  identical(unname(f$label), c("yy", "zz")) &&
    identical(f$short, f$label) &&
    identical(f$short_fig, f$label),
  "label = short = short_fig = group string"
)
expect(
  length(f$colour) == 2L && all(grepl("^#", f$colour)),
  "fallback colours assigned"
)
expect(identical(f$label_prefix, ""), "empty label_prefix")

# --- [3] sidecar naming unknown groups is refused --------------------------------
# A sidecar naming a group absent from the registry is refused.
cat("\n[3] sidecar naming an unknown group\n")
writeLines(c("order: [yy, ww]"), file.path(td, "nosc.yml"))
e <- err_of(load_grouping_meta("nosc", groupings_dir = td))
expect(!is.null(e) && grepl("ww", e), "unknown sidecar group refused, named")
unlink(file.path(td, "nosc.yml"))

# --- [4] registry_digest ---------------------------------------------------------
# registry_digest() is stable, file-content-based and changes with content.
cat("\n[4] registry_digest\n")
d1 <- registry_digest("nosc", groupings_dir = td)
expect(identical(d1, registry_digest("nosc", groupings_dir = td)), "stable")
expect(is.character(d1) && nchar(d1) == 16L, "16-hex xxhash64")
writeLines(
  c("iso3,group", "AAA,zz", "BBB,yy", "CCC,yy"),
  file.path(td, "nosc.csv")
)
expect(
  !identical(d1, registry_digest("nosc", groupings_dir = td)),
  "changes with content"
)
expect(
  !is.null(err_of(registry_digest("absent", groupings_dir = td))),
  "missing registry refused"
)

# --- [5] resolve_group_label -------------------------------------------------------
# Short and full forms; unknown refused.
cat("\n[5] resolve_group_label\n")
expect(
  identical(
    resolve_group_label(m, c("east", "enso_horn")),
    c("enso_east", "enso_horn")
  ),
  "short and full forms resolve"
)
expect(
  !is.null(err_of(resolve_group_label(m, "nowhere"))),
  "unknown label refused"
)

# --- [6] tokens -----------------------------------------------------------------
# grouping_label_token() / grouping_source_from_token(): unique round trip for
# every shipped registry; ambiguous and unknown tokens refused.
cat("\n[6] grouping_label_token / grouping_source_from_token\n")
expect(
  identical(grouping_label_token("enso_east"), "ensoeast"),
  "token strips non-alnum"
)
srcs <- sub(
  "\\.csv$",
  "",
  list.files("analysis/groupings", pattern = "\\.csv$")
)
toks <- tolower(grouping_label_token(srcs))
expect(
  !anyDuplicated(toks),
  sprintf(
    "shipped registry tokens are unique (%s)",
    paste(toks, collapse = ", ")
  )
)
rt <- vapply(
  srcs,
  function(s) grouping_source_from_token(grouping_label_token(s)),
  character(1)
)
expect(
  identical(unname(rt), srcs),
  "token -> source round trip for every shipped registry"
)
writeLines(c("iso3,group", "AAA,a"), file.path(td, "a_b.csv"))
writeLines(c("iso3,group", "AAA,a"), file.path(td, "ab.csv"))
expect(
  !is.null(err_of(grouping_source_from_token("ab", groupings_dir = td))),
  "ambiguous token refused"
)
expect(
  !is.null(err_of(grouping_source_from_token("nope", groupings_dir = td))),
  "unknown token refused"
)

# --- [7] assert_partition_labels ------------------------------------------------
# Path-unsafe labels refused.
cat("\n[7] assert_partition_labels\n")
expect(
  is.null(err_of(assert_partition_labels(c(
    "enso_east",
    "grp_a",
    "hmmCE",
    "TU-1"
  )))),
  "safe labels accepted"
)
for (bad in c("a b", "_x", "a/b", "_NA_TU_", "x_Np=5", "y_set-neutral", "a*")) {
  expect(
    !is.null(err_of(assert_partition_labels(bad))),
    sprintf("'%s' refused", bad)
  )
}

unlink(td, recursive = TRUE)
cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
