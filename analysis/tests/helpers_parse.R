# analysis/tests/helpers_parse.R
#
# Shared parse-tree helpers for the static (no-run) tests: who calls a
# function, which strings a file holds, ... read off utils::getParseData()
# rather than off a regex over the text, so comments, strings and names that
# merely contain the target do not match.
#
# Usage: source("analysis/tests/helpers_parse.R") near the top of a test file.
# Base R only (utils::getParseData).
#
# Holds parse_data_of(), callers_of() and the `%||%` literal-site helpers
# (nullor_*, below; used by test_hardcoded_default_parity.R [6] and [7]).

#' parse_data_of -- the parse data of an R file.
#'
#' @param f path to an R file.
#' @return the data frame of utils::getParseData() (line1, col1, token, text,
#'   parent, terminal, ...), or NULL when the file does not parse.
parse_data_of <- function(f) {
  tryCatch(utils::getParseData(parse(f, keep.source = TRUE)),
           error = function(e) NULL)
}

#' callers_of -- every use of a function name in a set of R files.
#'
#' One row per terminal token naming `fn`:
#'   kind "call"       `fn(...)` or `pkg::fn(...)` (SYMBOL_FUNCTION_CALL)
#'   kind "reference"  the bare symbol passed as a value, e.g. lapply(x, fn)
#'   kind "definition" the symbol on the left of `<-` / `<<-` / `=`
#'   kind "string"     the quoted name, e.g. do.call("fn", ...), get("fn")
#' plus one row of kind "unparseable" (line NA) per file that does not parse,
#' so a caller cannot hide in a broken file.
#'
#' @param fn function name, character(1).
#' @param paths R files and/or directories (searched recursively for *.R).
#' @return data.frame(file, line, kind), possibly with zero rows.
callers_of <- function(fn, paths) {
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
    else p
  }), use.names = FALSE)
  rows <- lapply(unique(files), function(f) {
    pd <- parse_data_of(f)
    if (is.null(pd))
      return(data.frame(file = f, line = NA_integer_, kind = "unparseable",
                        stringsAsFactors = FALSE))
    term <- pd[pd$terminal, , drop = FALSE]
    term <- term[order(term$line1, term$col1), , drop = FALSE]
    nxt  <- c(term$token[-1L], "")
    quoted <- c(sprintf('"%s"', fn), sprintf("'%s'", fn))
    kind <- ifelse(term$token == "SYMBOL_FUNCTION_CALL" & term$text == fn, "call",
            ifelse(term$token == "SYMBOL" & term$text == fn,
                   ifelse(nxt %in% c("LEFT_ASSIGN", "EQ_ASSIGN"), "definition", "reference"),
            ifelse(term$token == "STR_CONST" & term$text %in% quoted, "string", NA_character_)))
    hit <- !is.na(kind)
    data.frame(file = rep(f, sum(hit)), line = term$line1[hit], kind = kind[hit],
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame(file = character(0), line = integer(0),
                                      kind = character(0), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# --- `%||%` literal sites ----------------------------------------------------
# A "site" is one operand of a `%||%` chain that ends in a literal:
# `a$k %||% b$k %||% "x"` gives two sites, (a$k, "x") and (b$k, "x"); the
# chain is read off the parse tree, so a literal on the line after its `%||%`
# is seen and comments / strings never match. Each operand is reduced to its
# key chain: `x$k`, `x[["k"]]` and `(x %||% list())$k` all read as x$k.

#' nullor_literal_value -- is a parsed right-hand side a literal?
#' Strings, numbers (incl. `L`, a unary minus), TRUE / FALSE, NA*, Inf.
#' @return list(ok = TRUE, value) or list(ok = FALSE).
nullor_literal_value <- function(e) {
  if (is.symbol(e) && identical(as.character(e), "Inf")) return(list(ok = TRUE, value = Inf))
  if (is.call(e) && identical(e[[1L]], as.name("-")) && length(e) == 2L &&
      is.numeric(e[[2L]]) && length(e[[2L]]) == 1L)
    return(list(ok = TRUE, value = -e[[2L]]))
  if (is.atomic(e) && length(e) == 1L && (is.character(e) || is.numeric(e) || is.logical(e)))
    return(list(ok = TRUE, value = e))
  list(ok = FALSE)
}

#' nullor_chain -- the key chain of an operand, or NULL when it is not one.
#' `x` -> "x"; `x$a$b` / `x[["a"]]$b` / `(x %||% list())$a$b` -> c("x", "a", "b").
nullor_chain <- function(e) {
  if (is.symbol(e)) return(as.character(e))
  if (!is.call(e)) return(NULL)
  f <- e[[1L]]
  if (identical(f, as.name("(")) && length(e) == 2L) return(nullor_chain(e[[2L]]))
  if (identical(f, as.name("%||%")) && length(e) == 3L && is.call(e[[3L]]) &&
      identical(e[[3L]][[1L]], as.name("list")) && length(e[[3L]]) == 1L)
    return(nullor_chain(e[[2L]]))
  if (identical(f, as.name("$")) && length(e) == 3L) {
    base <- nullor_chain(e[[2L]])
    k <- e[[3L]]
    if (is.null(base) || !(is.symbol(k) || is.character(k))) return(NULL)
    return(c(base, as.character(k)))
  }
  if (identical(f, as.name("[[")) && length(e) == 3L && is.character(e[[3L]]) &&
      length(e[[3L]]) == 1L) {
    base <- nullor_chain(e[[2L]])
    if (is.null(base)) return(NULL)
    return(c(base, e[[3L]]))
  }
  NULL
}

#' nullor_split_operands -- the operands of a left-associative `%||%` chain:
#' `a %||% b %||% c` (parsed `(a %||% b) %||% c`) -> list(a, b, c).
nullor_split_operands <- function(e) {
  if (is.call(e) && identical(e[[1L]], as.name("%||%")) && length(e) == 3L)
    return(c(nullor_split_operands(e[[2L]]), list(e[[3L]])))
  list(e)
}

#' nullor_literal_sites -- every `<operand> %||% ... %||% <literal>` site of a file.
#'
#' @param f path to an R file.
#' @param dedup TRUE (default): one row per (expr, literal), `line` = the
#'   first occurrence (the same fallback at several call sites is one fact).
#' @return data.frame(file (basename), line, expr (the operand's chain as
#'   `a$b$c`, or its deparse when it is no chain), literal (source text),
#'   chain (logical: the operand is a key chain), pos (operand position),
#'   n_operands). Operands that are themselves literals are dropped. Zero rows
#'   when the file has no site; stop() when it does not parse.
nullor_literal_sites <- function(f, dedup = TRUE) {
  pd <- parse_data_of(f)
  if (is.null(pd)) stop("nullor_literal_sites(): ", f, " does not parse.", call. = FALSE)
  ops <- pd[pd$token == "SPECIAL" & pd$text == "%||%", , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(ops))) {
    pid <- ops$parent[i]
    kids <- pd[pd$parent == pid & pd$token == "expr", , drop = FALSE]
    kids <- kids[order(kids$line1, kids$col1), , drop = FALSE]
    if (nrow(kids) != 2L) next
    rhs_txt <- utils::getParseText(pd, kids$id[2L])
    rhs <- tryCatch(str2lang(rhs_txt), error = function(e) NULL)
    if (is.null(rhs) || !isTRUE(nullor_literal_value(rhs)$ok)) next
    lhs <- str2lang(utils::getParseText(pd, kids$id[1L]))
    operands <- nullor_split_operands(lhs)
    for (k in seq_along(operands)) {
      o <- operands[[k]]
      if (isTRUE(nullor_literal_value(o)$ok)) next
      ch <- nullor_chain(o)
      rows[[length(rows) + 1L]] <- data.frame(
        file = basename(f), line = pd$line1[pd$id == pid],
        expr = if (is.null(ch)) paste(deparse(o, width.cutoff = 500L), collapse = " ")
               else paste(ch, collapse = "$"),
        literal = rhs_txt, chain = !is.null(ch), pos = k, n_operands = length(operands),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows))
    return(data.frame(file = character(0), line = integer(0), expr = character(0),
                      literal = character(0), chain = logical(0), pos = integer(0),
                      n_operands = integer(0), stringsAsFactors = FALSE))
  df <- do.call(rbind, rows)
  df <- df[order(df$line), , drop = FALSE]
  if (dedup) df <- df[!duplicated(df[c("expr", "literal")]), , drop = FALSE]
  rownames(df) <- NULL
  df
}

#' nullor_resolve -- the default key path(s) an operand reads.
#'
#' @param file basename of the site's file.
#' @param expr the operand's chain text (`a$b$c`).
#' @param root_map list of list(file = basename or "*", root = chain prefix
#'   text, prefix = list of key-path prefixes); the longest matching root
#'   wins, a file-specific entry over "*".
#' @param expr_map named list: exact expr -> key path (checked first).
#' @return list(kp = list of key paths or NULL, root, bare = TRUE when the
#'   operand IS a mapped root with no key after it).
nullor_resolve <- function(file, expr, root_map, expr_map = list()) {
  if (!is.null(expr_map[[expr]])) return(list(kp = list(expr_map[[expr]]), root = expr))
  tok <- strsplit(expr, "$", fixed = TRUE)[[1]]
  for (i in rev(seq_along(tok))) {
    root <- paste(tok[seq_len(i)], collapse = "$")
    hit <- Filter(function(m) identical(m$root, root) && m$file %in% c(file, "*"), root_map)
    if (length(hit)) {
      hit <- hit[order(vapply(hit, function(m) m$file == "*", logical(1)))][[1]]
      rest <- tok[-seq_len(i)]
      if (!length(rest)) return(list(kp = NULL, root = root, bare = TRUE))
      return(list(kp = lapply(hit$prefix, function(p) c(p, rest)), root = root))
    }
  }
  list(kp = NULL, root = tok[1])
}

#' nullor_norm_value -- a value in comparable form (numbers as double).
nullor_norm_value <- function(v) {
  if (is.null(v)) return(NULL)
  u <- unlist(v, use.names = FALSE)
  if (is.numeric(u)) as.numeric(u) else u
}

#' Literal fallbacks that are sentinels ("not set", "auto"), never a default:
#' exempt from the comparison by this fixed list (case-insensitive) and NA*.
NULLOR_SENTINELS <- c("", "(none)", "(auto)", "<none>", "<null>", "<absent>", "unknown",
                      "?", "(top-level)")

#' nullor_classify -- the status of every site against the defaults.
#'
#' @param sites nullor_literal_sites() rows.
#' @param defaults get_default_options().
#' @param root_map,expr_map as nullor_resolve().
#' @param ignored list of list(file = basename or "*", root = chain prefix
#'   text, reason): roots that are not config lists (locals, artefacts,
#'   function arguments, CLI options), each with the reason.
#' @param effective function(kp) -> the effective default at key path kp
#'   (default: the leaf of `defaults`, read through config_get_leaf()).
#' @return `sites` plus status (SENTINEL | NONCHAIN | IGNORED | BARE |
#'   UNCLASSIFIED | NOKEY | UNSET | MATCH | MISMATCH) and detail (the key
#'   path(s) and default(s) for a MISMATCH).
nullor_classify <- function(sites, defaults, root_map, expr_map = list(),
                            ignored = list(),
                            effective = function(kp) config_get_leaf(defaults, paste(kp, collapse = "."))) {
  has_path <- function(kp) {
    d <- defaults
    for (p in kp) { if (!is.list(d) || !(p %in% names(d))) return(FALSE); d <- d[[p]] }
    TRUE
  }
  status <- character(nrow(sites)); detail <- character(nrow(sites))
  for (i in seq_len(nrow(sites))) {
    s <- sites[i, ]
    val <- nullor_literal_value(str2lang(s$literal))$value
    if ((is.character(val) && tolower(val) %in% NULLOR_SENTINELS) || (length(val) == 1L && is.na(val))) {
      status[i] <- "SENTINEL"; next
    }
    if (!isTRUE(s$chain)) { status[i] <- "NONCHAIN"; next }
    res <- nullor_resolve(s$file, s$expr, root_map, expr_map)
    if (is.null(res$kp)) {
      tok <- strsplit(s$expr, "$", fixed = TRUE)[[1]]
      pre <- vapply(seq_along(tok), function(j) paste(tok[seq_len(j)], collapse = "$"), "")
      ign <- any(vapply(ignored, function(g) g$root %in% pre && g$file %in% c(s$file, "*"), logical(1)))
      status[i] <- if (isTRUE(res$bare)) "BARE" else if (ign) "IGNORED" else "UNCLASSIFIED"
      next
    }
    any_def <- FALSE; any_cmp <- FALSE; bad <- character(0)
    for (kp in res$kp) {
      if (!has_path(kp)) next
      any_def <- TRUE
      eff <- effective(kp)
      if (is.null(eff)) next
      any_cmp <- TRUE
      if (!identical(nullor_norm_value(eff), nullor_norm_value(val)))
        bad <- c(bad, sprintf("%s = %s", paste(kp, collapse = "."),
                              paste(format(unlist(eff)), collapse = ",")))
    }
    status[i] <- if (!any_def) "NOKEY" else if (!any_cmp) "UNSET" else if (length(bad)) "MISMATCH" else "MATCH"
    detail[i] <- paste(bad, collapse = "; ")
  }
  sites$status <- status; sites$detail <- detail
  sites
}
