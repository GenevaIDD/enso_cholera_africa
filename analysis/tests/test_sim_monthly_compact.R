# test_sim_monthly_compact.R — HARD: the compact (format-2) monthly summary
#   (utils/sim_io.R) against the long table it replaces on disk.
#
# Usage:
#   Rscript analysis/tests/test_sim_monthly_compact.R
#
# Details:
#   Quick-safe, synthetic temp files only. Every reader must see exactly what
#   it saw in the long format: round trip, long fallbacks, sim_monthly_cells(),
#   aggregate_partitions() / combine_group_runs(), the 09 loader,
#   utils/compact_sim_outputs.R and 09 stamp freshness across the conversion.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(qs2)
})
PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
for (m in c("base", "gid", "paths", "sim_io", "results_core")) {
  source(sprintf("analysis/utils/%s.R", m))
}

#' A monthly summary as simulate_one_tu() flattens its R x N x M arrays.
kernel_ms <- function(
  R = 7L,
  gids = c("KEN.1.1_1", "KEN.1.2_1", "KEN.2_1", "UGA.3.1_1"),
  M = 5L,
  cols = c("Y", "active", "mu_mon", "cases_month", "detected"),
  seed = 1L,
  draw = TRUE,
  rids = seq_len(R)
) {
  set.seed(seed)
  N <- length(gids)
  months <- seq(as.Date("2027-01-01"), by = "month", length.out = M)
  A <- array(rbinom(R * N * M, 1L, .3), c(R, N, M))
  Y <- A * array(rbinom(R * N * M, 1L, .4), c(R, N, M))
  D <- array(
    as.integer(A == 1L & runif(R * N * M) < .7 | runif(R * N * M) < .02),
    c(R, N, M)
  )
  MU <- array(NA_real_, c(R, N, M))
  C <- array(NA_integer_, c(R, N, M))
  act <- which(A == 1L)
  MU[act] <- rgamma(length(act), 2, .1)
  if (draw) {
    cs <- as.numeric(rnbinom(length(act), mu = MU[act], size = .7))
    cs[intersect(c(1L, 3L), seq_along(cs))] <- NA
    C[act] <- cs
  }
  ms <- tibble::tibble(
    replicate_id = rep(rids, each = N * M),
    gid = rep(rep(gids, each = M), times = R),
    month = rep(rep(as.character(months), times = N), times = R),
    Y = as.vector(aperm(Y, c(3L, 2L, 1L)))
  )
  fl <- function(a) as.vector(aperm(a, c(3L, 2L, 1L)))
  if ("active" %in% cols) {
    ms$active <- fl(A)
  }
  if ("mu_mon" %in% cols) {
    ms$mu_mon <- fl(MU)
    ms$cases_month <- fl(C)
  }
  if ("detected" %in% cols) {
    ms$detected <- fl(D)
  }
  ms
}
rt <- function(ms) {
  cx <- sim_monthly_compact(ms)
  !is.null(cx) && identical(sim_monthly_expand(cx), ms)
}
td <- tempfile("smc_")
dir.create(td)

# [1] round trip: simulate_one_tu()-shaped tables (every column set the kernel
#     writes, a never-drawn integer cases column, NA draws on active cells,
#     detected cells off the active ones, factor gids, double ids) expand
#     identical(); sim_io_read() expands; expand(cols =) = the long table's
#     columns; rds too.
cat("# [1] round trip\n")
expect(
  rt(kernel_ms()),
  "full column set (Y, active, mu_mon, cases_month, detected)"
)
expect(rt(kernel_ms(cols = "Y")), "legacy Y-only table")
expect(rt(kernel_ms(cols = c("Y", "active"))), "duration-only table")
expect(
  rt(kernel_ms(cols = c("Y", "active", "mu_mon"))),
  "mon track without detected"
)
ms0 <- kernel_ms(draw = FALSE)
expect(
  is.integer(ms0$cases_month) && rt(ms0),
  "never-drawn cases column stays integer NA"
)
expect(
  rt(kernel_ms(rids = c(5, 9, 11, 2, 3, 4, 8))),
  "double, unsorted replicate ids"
)
msf <- kernel_ms()
msf$gid <- factor(msf$gid)
expect(rt(msf), "factor gid")
expect(rt(kernel_ms(R = 1L, M = 1L)), "one replicate, one month")
cx <- sim_monthly_compact(kernel_ms())
expect(
  identical(cx$blocks[[1]]$cols$cases_month$enc, "ref"),
  "cases_month stored on the active cells"
)
f <- file.path(td, "ms.qs")
sim_monthly_save(kernel_ms(), f, "qs")
expect(
  is_sim_monthly_compact(qs_read(f)) &&
    identical(sim_io_read(f, "qs"), kernel_ms()),
  "sim_monthly_save() writes format 2; sim_io_read() returns the long table"
)
expect(
  is_sim_monthly_compact(sim_io_read(f, "qs", expand = FALSE)),
  "expand = FALSE keeps it compact"
)
fr <- file.path(td, "ms.rds")
sim_monthly_save(kernel_ms(), fr, "rds")
expect(
  is_sim_monthly_compact(readRDS(fr)) &&
    identical(sim_io_read(fr, "rds"), kernel_ms()),
  "rds"
)
keep <- c("replicate_id", "gid", "month", "Y", "active", "cases_month")
expect(
  identical(
    sim_monthly_expand(cx, cols = c("Y", "active", "cases_month")),
    kernel_ms()[, keep]
  ),
  "expand(cols =): the long table's columns, in its order"
)

# [2] a table that is not a dense grid (shuffled, censored) is written long;
#     parquet stays long; options(cholera.sim_monthly_long) forces long.
cat("# [2] long fallbacks\n")
sh <- kernel_ms()[sample(nrow(kernel_ms())), ]
expect(is.null(sim_monthly_compact(sh)), "shuffled rows: not a dense grid")
f2 <- file.path(td, "sh.qs")
sim_monthly_save(sh, f2, "qs")
expect(
  tibble::is_tibble(qs_read(f2)) && identical(qs_read(f2), sh),
  "written long"
)
cen <- kernel_ms() %>% filter(!(gid == "UGA.3.1_1" & month > "2027-03-01"))
expect(is.null(sim_monthly_compact(cen)), "censored grid: not dense")
options(cholera.sim_monthly_long = TRUE)
f3 <- file.path(td, "lg.qs")
sim_monthly_save(kernel_ms(), f3, "qs")
options(cholera.sim_monthly_long = NULL)
expect(
  tibble::is_tibble(qs_read(f3)),
  "options(cholera.sim_monthly_long = TRUE) writes long"
)
if (requireNamespace("nanoparquet", quietly = TRUE)) {
  f4 <- file.path(td, "ms.parquet")
  sim_monthly_save(cx, f4, "parquet")
  expect(
    identical(
      as.data.frame(sim_io_read(f4, "parquet")),
      as.data.frame(kernel_ms())
    ),
    "parquet: long"
  )
} else {
  cat("  SKIP  nanoparquet not installed\n")
}

# [3] sim_monthly_cells(): compact == long (rows, axes), a row subset of the
#     long table.
cat("# [3] sim_monthly_cells()\n")
kg <- function(g) g[!duplicated(sub("\\.[0-9]+_1$", "", g))]
for (nm in c("full", "Y")) {
  m <- if (nm == "full") kernel_ms(R = 20L) else kernel_ms(R = 20L, cols = "Y")
  a <- sim_monthly_cells(sim_monthly_compact(m), c("Y", "cases_month"), kg)
  b <- sim_monthly_cells(m, c("Y", "cases_month"), kg)
  expect(identical(a, b), sprintf("%s: compact == long", nm))
  rid <- paste(m$replicate_id, m$gid, m$month)
  expect(
    identical(
      m[match(paste(a$rows$replicate_id, a$rows$gid, a$rows$month), rid), ],
      a$rows
    ) &&
      !is.unsorted(match(
        paste(a$rows$replicate_id, a$rows$gid, a$rows$month),
        rid
      )),
    sprintf("%s: a row subset of the long table, in its order", nm)
  )
}

# [4] aggregate_partitions(): merged compact file expands to bind_rows() of the
#     parts; the return equals the long-format run's; mixed long and compact
#     parts; combine_group_runs() writes compact, reads back equal.
cat("# [4] aggregate_partitions() / combine_group_runs()\n")
mk_run <- function(root, compact, long_tu = NULL) {
  for (tu in c("ke", "ug")) {
    for (off in c(0L, 3L)) {
      pd <- file.path(
        root,
        "parts",
        paste0("tu=", tu),
        sprintf("repoff=%06d", off)
      )
      dir.create(pd, recursive = TRUE)
      g <- if (tu == "ke") {
        c("KEN.1.1_1", "KEN.2_1")
      } else {
        c("UGA.3.1_1", "UGA.4_1")
      }
      ms <- kernel_ms(
        R = 3L,
        gids = g,
        rids = off + 1:3,
        seed = off + nchar(tu)
      )
      if (compact && !identical(tu, long_tu)) {
        sim_monthly_save(ms, file.path(pd, "monthly_summary.qs"), "qs")
      } else {
        qs_save(ms, file.path(pd, "monthly_summary.qs"))
      }
      qs_save(
        tibble(replicate_id = off + 1L, gid = g[1], start_month = "2027-01-01"),
        file.path(pd, "event_table.qs")
      )
      file.create(file.path(pd, "_SUCCESS"))
    }
  }
}
mk_run(file.path(td, "L"), FALSE)
mk_run(file.path(td, "C"), TRUE)
mk_run(file.path(td, "X"), TRUE, long_tu = "ug")
aL <- suppressMessages(aggregate_partitions(
  file.path(td, "L", "parts"),
  file.path(td, "L", "merged")
))
aC <- suppressMessages(aggregate_partitions(
  file.path(td, "C", "parts"),
  file.path(td, "C", "merged")
))
aC2 <- suppressMessages(aggregate_partitions(
  file.path(td, "C", "parts"),
  file.path(td, "C", "merged2"),
  return_long = FALSE
))
aX <- suppressMessages(aggregate_partitions(
  file.path(td, "X", "parts"),
  file.path(td, "X", "merged")
))
mL <- qs_read(file.path(td, "L", "merged", "monthly_summary.qs"))
mC <- qs_read(file.path(td, "C", "merged", "monthly_summary.qs"))
expect(
  is_sim_monthly_compact(mC) && length(mC$blocks) == 4L,
  "merged file: one block per partition"
)
expect(
  identical(sim_monthly_expand(mC), mL),
  "merged compact expands to bind_rows() of the long parts"
)
expect(identical(aC, aL), "return value equals the long-format run's")
expect(
  is_sim_monthly_compact(aC2$monthly_summary) &&
    identical(sim_monthly_expand(aC2$monthly_summary), mL),
  "return_long = FALSE: the compact object as written"
)
expect(
  identical(aX, aL) &&
    tibble::is_tibble(qs_read(file.path(
      td,
      "X",
      "merged",
      "monthly_summary.qs"
    ))),
  "mixed long and compact parts: merged long, same table"
)
cL <- suppressMessages(combine_group_runs(
  list(a = file.path(td, "L", "merged")),
  file.path(td, "L", "comb"),
  "2027-03-01"
))
cC <- suppressMessages(combine_group_runs(
  list(a = file.path(td, "C", "merged")),
  file.path(td, "C", "comb"),
  "2027-03-01"
))
expect(identical(cC, cL), "combine_group_runs(): same return on compact input")
fc <- file.path(td, "C", "comb", "monthly_summary.qs")
expect(
  identical(
    sim_io_read(fc, "qs"),
    qs_read(file.path(td, "L", "comb", "monthly_summary.qs"))
  ),
  "combined file reads back as the long-format combined table"
)

# [5] the 09 loader: every load of tests/test_results_stitched_load_golden.R on
#     compact parts (and on detected-series parts) equals the golden / the
#     long-format loads.
cat("# [5] the 09 loader on compact parts\n")
counts <- c(PASS, FAIL)
source("analysis/tests/test_results_stitched_load_golden.R") # stitched_load_world(), stitched_loads(); resets the counters
PASS <- counts[1]
FAIL <- counts[2]
tw <- file.path(td, "stl")
dir.create(tw)
W <- stitched_load_world(tw)
for (p in Sys.glob(file.path(
  W$sim,
  "fc",
  "parts",
  "*",
  "monthly_summary.qs"
))) {
  sim_monthly_save(qs_read(p), p, "qs")
}
expect(
  all(vapply(
    Sys.glob(file.path(W$sim, "fc", "parts", "*", "monthly_summary.qs")),
    function(p) is_sim_monthly_compact(qs_read(p)),
    logical(1)
  )),
  "fixture parts rewritten in format 2"
)
gold <- readRDS(GOLDEN)
L <- stitched_loads(W)
for (nm in names(gold)) {
  expect(
    identical(gold[[nm]], L[[nm]]),
    sprintf("%s: compact parts == golden", nm)
  )
}
# detected series and an empty ADM1 window, long vs compact
dw <- file.path(td, "det")
for (fmt in c("long", "compact")) {
  pd <- file.path(dw, fmt, "r", "parts", "tu=k")
  dir.create(pd, recursive = TRUE)
  ms <- kernel_ms(
    R = 12L,
    gids = c("KEN.1.1_1", "KEN.1.2_1", "KEN.2_1", "UGA.3.1_1", "UGA3_1"),
    seed = 9L
  )
  if (fmt == "long") {
    qs_save(ms, file.path(pd, "monthly_summary.qs"))
  } else {
    sim_monthly_save(ms, file.path(pd, "monthly_summary.qs"), "qs")
  }
  file.create(file.path(pd, "_SUCCESS"))
}
ld <- function(fmt, ...) load_sim_run_monthly(file.path(dw, fmt), "r", ...)
expect(
  identical(
    ld("long", "2027-01-01", "2027-05-01", series = "detected"),
    ld("compact", "2027-01-01", "2027-05-01", series = "detected")
  ),
  "detected series"
)
expect(
  identical(
    ld(
      "long",
      "2027-02-01",
      "2027-03-01",
      horizon_by_country = list(KEN = "2027-03-01", UGA = "2027-05-01"),
      keep_beyond = TRUE
    ),
    ld(
      "compact",
      "2027-02-01",
      "2027-03-01",
      horizon_by_country = list(KEN = "2027-03-01", UGA = "2027-05-01"),
      keep_beyond = TRUE
    )
  ),
  "per-country horizons, flagged"
)
expect(
  identical(
    ld("long", "2030-01-01", "2030-05-01"),
    ld("compact", "2030-01-01", "2030-05-01")
  ),
  "ADM1 window outside the run (empty adm1_month)"
)

# [6] utils/compact_sim_outputs.R (CLI): --dry-run writes nothing; parts and
#     merged files of long-format runs are converted and read back identical
#     with their mtimes kept; a non-grid file stays long; a second pass finds
#     everything already compact.
cat("# [6] compact_sim_outputs.R\n")
root <- file.path(td, "conv")
dir.create(root)
for (r in c("runA", "runB")) {
  mk_run(file.path(root, r), FALSE)
  quiet <- suppressMessages(aggregate_partitions(
    file.path(root, r, "parts"),
    file.path(root, r, "merged"),
    return_long = FALSE
  ))
}
odd <- file.path(
  root,
  "runB",
  "parts",
  "tu=ke",
  "repoff=000000",
  "monthly_summary.qs"
)
qs_save(qs_read(odd)[sample(nrow(qs_read(odd))), ], odd) # not a dense grid any more
fs <- sort(list.files(
  root,
  "^monthly_summary\\.qs$",
  recursive = TRUE,
  full.names = TRUE
))
before <- lapply(stats::setNames(fs, fs), qs_read)
old_t <- as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
for (f in fs) {
  Sys.setFileTime(f, old_t)
}
cli <- function(...) {
  system2(
    "Rscript",
    c("analysis/utils/compact_sim_outputs.R", "--root", root, ...),
    stdout = TRUE,
    stderr = TRUE
  )
}
o <- cli("--pattern", "'run*'", "--dry-run")
expect(
  identical(lapply(stats::setNames(fs, fs), qs_read), before) &&
    all(file.mtime(fs) == old_t),
  "--dry-run reads and writes nothing"
)
lg <- file.path(td, "conv.csv")
o <- cli("--pattern", "'run*'", "--log", lg)
st <- utils::read.csv(lg, stringsAsFactors = FALSE)
expect(
  nrow(st) == length(fs) && is.null(attr(o, "status")),
  "every file processed, exit 0"
)
expect(
  all(vapply(
    fs,
    function(f) identical(sim_io_read(f, "qs"), before[[f]]),
    logical(1)
  )),
  "every file reads back as the long table it held"
)
expect(all(file.mtime(fs) == old_t), "modification times kept")
conv <- st$path[st$status == "compacted"]
expect(
  length(conv) == length(fs) - 2L &&
    all(vapply(
      conv,
      function(f) is_sim_monthly_compact(qs_read(f)),
      logical(1)
    )),
  "all but the shuffled part and its run's merged file compacted"
)
expect(
  identical(
    sort(st$path[st$status == "left_long"]),
    sort(c(odd, file.path(root, "runB", "merged", "monthly_summary.qs")))
  ),
  "the shuffled part and runB's merged file left long"
)
expect(
  is_sim_monthly_compact(qs_read(file.path(
    root,
    "runA",
    "merged",
    "monthly_summary.qs"
  ))) &&
    length(
      qs_read(file.path(root, "runA", "merged", "monthly_summary.qs"))$blocks
    ) ==
      4L,
  "merged file rebuilt from its compact parts (one block each)"
)
o <- cli("--pattern", "'runA'", "--log", lg)
expect(
  all(utils::read.csv(lg)$status == "already_compact"),
  "second pass: already compact"
)

# [7] 09 freshness: an artefact whose stamp recorded a long-format input stays
#     fresh after the tool converts that input (sidecar), and is stale without
#     the sidecar, with a sidecar naming another size, or after the file
#     changes again.
cat("# [7] 09 stamp freshness across the conversion\n")
fd <- file.path(td, "fresh", "run", "merged")
dir.create(fd, recursive = TRUE)
inp <- file.path(fd, "monthly_summary.qs")
qs_save(kernel_ms(R = 30L), inp)
Sys.setFileTime(inp, as.POSIXct("2026-11-04 14:20:00", tz = "UTC"))
cfg <- file.path(td, "fresh", "results.yml")
writeLines("results: {}", cfg)
res <- list(
  inputs = list(),
  scope = list(a = 1),
  derivation = list(),
  thresholds = list()
)
art <- file.path(td, "fresh", "artefact.qs")
invisible(suppressMessages(results_write(
  list(x = 1),
  art,
  res,
  cfg,
  "v1",
  inputs = inp
)))
stale <- function() {
  inherits(
    try(
      suppressMessages(suppressWarnings(results_require(art, cfg, res))),
      silent = TRUE
    ),
    "try-error"
  )
}
Sys.unsetenv("RESULTS_ALLOW_STALE_INPUTS")
expect(!stale(), "fresh before the conversion")
source("analysis/utils/compact_sim_outputs.R") # compact_one_file(); CLI guarded
cv <- compact_one_file(inp)
expect(
  identical(cv$status, "compacted") &&
    cv$bytes_after < cv$bytes_before &&
    file.exists(sim_monthly_sidecar_path(inp)),
  "converted (smaller), sidecar written"
)
expect(!stale(), "still fresh after the conversion (sidecar)")
sc <- sim_monthly_sidecar_path(inp)
keep <- readLines(sc)
invisible(file.rename(sc, paste0(sc, ".off")))
expect(stale(), "no sidecar: stale (size changed)")
invisible(file.rename(paste0(sc, ".off"), sc))
d <- read.dcf(sc)
d[1, "original_size"] <- "123"
write.dcf(d, sc)
expect(stale(), "sidecar naming another original size: stale")
writeLines(keep, sc)
expect(!stale(), "sidecar restored: fresh")
mt <- file.mtime(inp)
sim_monthly_save(kernel_ms(R = 31L), inp, "qs")
Sys.setFileTime(inp, mt)
expect(stale(), "file rewritten again (same mtime, other size): stale")
Sys.setFileTime(inp, mt + 60)
expect(stale(), "file touched later: stale")

unlink(td, recursive = TRUE)
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
