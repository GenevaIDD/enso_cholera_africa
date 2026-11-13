
`%||%` <- function(a, b) if (is.null(a)) b else a

# Weighted quantile (step-function empirical CDF, no external deps). Shared by
# utils/smoothing_core.R (smoothing ribbons) and utils/results_derive.R (phase-09
# smoothing ribbons); both reach the one copy through base.R, which every
# caller sources.
wquantile <- function(x, w, probs) {
  ord <- order(x, na.last = TRUE)
  xs  <- x[ord]; ws <- w[ord]
  ws  <- ws / sum(ws, na.rm = TRUE)
  cw  <- cumsum(ws)
  vapply(probs, function(p) {
    i <- which(cw >= p)[1L]
    if (is.na(i)) xs[length(xs)] else xs[i]
  }, numeric(1L))
}

# ENSO single-lag naming. A negative single_lag -k is a LEAD: the
# ONI k months AFTER the outcome month (joinEnsoLags() subtracts the lag from
# the date). Every writer and reader of a lag-keyed name goes through these, so
# lag -3 and lag 3 never share a column, config file or fit file; non-negative
# lags keep their historical names byte-for-byte.
#   ensoLagColName(-3) "enso_lead3"   ensoLagColName(3) "enso_lag3"
#   slToken(-3)        "m3"           slToken(3)        "3"     (config: _sl-<tok>_grp-)
#   slTag(-3)          "slm3"         slTag(3)          "sl3"   (fit-file tag)
ensoLagColName <- function(lag) {
  lag <- as.integer(lag)
  ifelse(lag < 0L, sprintf("enso_lead%d", -lag), sprintf("enso_lag%d", lag))
}
slToken <- function(lag) {
  lag <- as.integer(lag)
  ifelse(lag < 0L, sprintf("m%d", -lag), sprintf("%d", lag))
}
slTag <- function(lag) paste0("sl", slToken(lag))

# Time to last outbreak (TTLO) categories: one table
# for the panel builder, the occurrence loader and the R-side simulator bins.
#   .TTLO_ORIGIN    "last outbreak" date of the rows before an LP's first onset
#                   in 03_prepare_occ.R (so their TTLO is days since 2000-01-01).
#   .TTLO_LEVELS    the factor levels, in order, of time_to_last_outbreak:
#                   cut((pmax(0, days) + 1) / 365, c(0:6, Inf), include.lowest =
#                   TRUE, right = FALSE), NA -> the top level "[6,Inf]".
#                   "[0,1)" is the reference level. Used by occ_data.R::
#                   makeTimeToLastOutbreakCat() and simulator_core.R
#                   (bin_ttlo(), ttlo_init_from_panel(), evaluate_lp()).
#   .TTLO_NONE_DAYS the finite "no recorded prior onset" sentinel (days),
#                   used by the level "none".
.TTLO_ORIGIN    <- as.Date("2000-01-01")
.TTLO_LEVELS    <- c("[0,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,6)", "[6,Inf]")
.TTLO_NONE_DAYS <- 1e9

#' lapply() over `X`, on `workers` forked processes when workers > 1.
#'
#' Results come back in input order. A worker's warnings are collected and
#' re-raised here, in order, after the map (a forked process cannot raise them
#' in this session); its messages print as they happen. An error in any worker
#' stops with that worker's message. workers <= 1 (or one element) is a plain
#' lapply(), so serial callers are unchanged. Fork-based: no parallelism on
#' Windows, where workers > 1 falls back to serial with a message.
mapWorkers <- function(X, FUN, workers = 1L) {
  workers <- as.integer(workers %||% 1L)
  if (is.na(workers) || workers <= 1L || length(X) <= 1L) return(lapply(X, FUN))
  if (.Platform$OS.type == "windows") {
    message("mapWorkers: forked workers unavailable on Windows; running serially.")
    return(lapply(X, FUN))
  }
  wrapped <- function(x) {
    warns <- list()
    val <- withCallingHandlers(FUN(x), warning = function(w) {
      warns[[length(warns) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    list(value = val, warnings = warns)
  }
  out <- parallel::mclapply(X, wrapped, mc.cores = min(workers, length(X)),
                            mc.preschedule = FALSE)
  failed <- vapply(out, function(o) inherits(o, "try-error") || is.null(o), logical(1))
  if (any(failed)) {
    i <- which(failed)[1L]
    msg <- if (inherits(out[[i]], "try-error")) conditionMessage(attr(out[[i]], "condition"))
           else "the worker returned nothing (killed, e.g. out of memory?)"
    stop(sprintf("mapWorkers: element %d of %d failed: %s", i, length(X), msg), call. = FALSE)
  }
  for (o in out) for (w in o$warnings) warning(w)
  lapply(out, `[[`, "value")
}

# ---------------------------------------------------------------------------
# Content digest
# ---------------------------------------------------------------------------

#' content_digest
#'
#' xxhash64 (digest package, `serialize = FALSE`) of a canonical text
#' serialisation of `x`: the one content hash for identities stamped into a
#' cache meta, a manifest or a panel stamp. Unlike `digest::digest(x)` it does
#' not hash R's own serialisation, so the value does not move with the R
#' version, with attributes, or with integer-vs-double storage, and two tables
#' with the same content hash the same.
#'
#' Canonicalisation (`.content_canonical_text()`); all text UTF-8 and every
#' line ends in "\n":
#' - raw vector: its bytes as they are (no text form). `config.R::
#'   registry_digest()` hashes a registry file's bytes this way, so its value
#'   equals `digest::digest(file = path, algo = "xxhash64")`.
#' - NULL and zero-length vectors: the empty text.
#' - atomic vector: one value per line, in the order given. Named: one
#'   `"<name>",<value>` line per element.
#' - data.frame (tibbles included): CSV text; a header line of the quoted
#'   column names, then one line per row. Columns and rows are kept in the
#'   order given (sort them first when their order should not matter). A list
#'   column is an error.
#' - matrix: as `as.data.frame(x)`.
#' - list: one block per element, in the order given: a header line
#'   `[[<quoted name, or the index when unnamed>]] <n bytes>` followed by the
#'   element's own canonical text of exactly n bytes (a raw element written as
#'   lowercase hex), so block boundaries are unambiguous.
#' Values: character quoted with `"`, inner quotes doubled (so "NA" and NA
#' differ); NA written `NA` unquoted; logical `TRUE` / `FALSE`; a factor as its
#' labels; Date `YYYY-MM-DD`; POSIXct ISO 8601 in UTC with microseconds;
#' integer and double both through `sprintf("%.<digits>g")` (R runs it in the C
#' numeric locale; 1L and 1.0 hash the same), -0 written as 0, and `NaN`,
#' `Inf`, `-Inf` as R prints them. `digits` (default 15 significant digits)
#' keeps last-bit floating-point noise out of the digest. Attributes other than
#' names, class and dim are ignored. Functions, environments, complex values
#' and other non-data objects are an error.
#'
#' Users: `config.R::registry_digest()` (raw file bytes), the Lenssen
#' counterfactual cache digest, the WHO-scale table digest and the contagion
#' coverage digests.
#'
#' @return a 16-character hex string; NA when the optional `digest` package is
#'   absent (as `registry_digest()` and `provenance.R::digest_or_null()`: a
#'   missing package weakens a check, never breaks the pipeline).
content_digest <- function(x, digits = 15L) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  bytes <- if (is.raw(x)) x else charToRaw(enc2utf8(.content_canonical_text(x, digits)))
  digest::digest(bytes, algo = "xxhash64", serialize = FALSE)
}

# Canonical text of the values of one atomic vector, one string per element
# (rules in content_digest()).
.content_canonical_values <- function(v, digits = 15L) {
  if (is.factor(v)) v <- as.character(v)
  if (!length(v)) return(character(0))
  if (inherits(v, "Date")) {
    out <- format(v, "%Y-%m-%d")
  } else if (inherits(v, "POSIXct")) {
    out <- format(v, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
  } else if (is.character(v)) {
    out <- paste0("\"", gsub("\"", "\"\"", enc2utf8(v), fixed = TRUE), "\"")
  } else if (is.logical(v)) {
    out <- ifelse(v, "TRUE", "FALSE")
  } else if (is.integer(v) || is.double(v)) {
    v <- as.double(v)
    v[!is.na(v) & v == 0] <- 0
    return(sprintf(paste0("%.", as.integer(digits), "g"), v))   # NA / NaN / Inf as printed
  } else {
    stop(sprintf("content_digest: unsupported atomic type '%s'", typeof(v)), call. = FALSE)
  }
  out[is.na(v)] <- "NA"
  out
}

.content_lines <- function(lines) if (length(lines)) paste0(lines, "\n", collapse = "") else ""

# The canonical text content_digest() hashes (a raw `x` is hashed as bytes and
# never reaches this function at the top level).
.content_canonical_text <- function(x, digits = 15L) {
  if (is.null(x)) return("")
  if (inherits(x, "POSIXlt")) x <- as.POSIXct(x)
  if (is.matrix(x) && !is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (any(vapply(x, is.list, logical(1))))
      stop("content_digest: a data.frame list column is not supported", call. = FALSE)
    header <- paste(.content_canonical_values(names(x), digits), collapse = ",")
    rows <- if (nrow(x) > 0L && ncol(x) > 0L)
      do.call(paste, c(unname(lapply(x, .content_canonical_values, digits = digits)), sep = ","))
    else character(0)
    return(.content_lines(c(header, rows)))
  }
  if (is.list(x)) {
    nm <- names(x)
    blocks <- vapply(seq_along(x), function(i) {
      el  <- x[[i]]
      txt <- if (is.raw(el)) .content_lines(paste(as.character(el), collapse = ""))
             else .content_canonical_text(el, digits)
      key <- if (!is.null(nm) && !is.na(nm[i]) && nzchar(nm[i]))
        .content_canonical_values(nm[i]) else as.character(i)
      paste0("[[", key, "]] ", nchar(txt, type = "bytes"), "\n", txt)
    }, character(1))
    return(paste(blocks, collapse = ""))
  }
  if (is.raw(x)) return(.content_lines(paste(as.character(x), collapse = "")))
  if (is.atomic(x)) {
    vals <- .content_canonical_values(x, digits)
    if (!is.null(names(x)))
      vals <- paste(.content_canonical_values(names(x), digits), vals, sep = ",")
    return(.content_lines(vals))
  }
  stop(sprintf("content_digest: cannot digest an object of type '%s'", typeof(x)), call. = FALSE)
}

# ---------------------------------------------------------------------------
# Month windows: one convention
#
# FORWARD HORIZONS ARE CLOSED [start, end]: simulator.horizon,
# filter.forecast_window, the manifest's horizon_end / firm_horizon_end /
# unit_horizon_end / africa_min_horizon, hindcast and event windows. Every
# month from start to end inclusive is simulated and read.
# FILTER WINDOWS ARE HALF-OPEN [start, end) with end = the handoff month:
# filter.warmup_window and the smoothing's FILTER_START / FILTER_END. The end
# month itself belongs to the forward arm that launches there.
# ---------------------------------------------------------------------------

#' horizonMonths
#'
#' Months of a CLOSED forward horizon [start, end] (first-of-month Dates);
#' empty when end < start.
horizonMonths <- function(start, end) {
  s <- as.Date(format(as.Date(start), "%Y-%m-01"))
  e <- as.Date(format(as.Date(end), "%Y-%m-01"))
  if (is.na(s) || is.na(e) || e < s) return(as.Date(character(0)))
  seq(s, e, by = "1 month")
}

#' filterWindowMonths
#'
#' Months of a HALF-OPEN filter window [start, end) (first-of-month Dates; the
#' end, the handoff month, is excluded); empty when end <= start.
filterWindowMonths <- function(start, end) {
  m <- horizonMonths(start, end)
  m[m < as.Date(format(as.Date(end), "%Y-%m-01"))]
}
