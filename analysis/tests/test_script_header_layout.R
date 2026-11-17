# test_script_header_layout.R — HARD: every R file's header follows the
#   template of README.md "Code style".
#
# Usage:
#   Rscript analysis/tests/test_script_header_layout.R
#   Rscript analysis/tests/test_script_header_layout.R --all   # every offender
#
# Details:
#   Header = the comment lines that open the file (after a shebang). Checks: the
#   first line is "# <file>.R — ...", at most one more description line, then
#   blocks separated by a bare "#": Usage: or Provides: (required), Details:
#   (at most 8 lines), See: (one line), in that order; no rule lines, no #',
#   lines up to 80 characters; a test's description starts with its run_all
#   gate class (HARD: / SOFT:); numbered scripts use numbered sections
#   starting with "1. Preamble". scratchpads/: the first line only.

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
show_all <- "--all" %in% commandArgs(TRUE)

cat("# R script headers follow the README.md template\n")

files <- system("git ls-files '*.R'", intern = TRUE)
files <- files[!grepl("^renv/", files)]

# Gate class of every test registered in run_all.R (its run_one() calls)
gate <- character(0)
walk_calls <- function(e) {
  if (!is.call(e)) {
    return(invisible())
  }
  if (identical(e[[1]], as.name("run_one")) && is.character(e[[2]])) {
    hf <- e[["hard_fail"]]
    if (is.logical(hf)) {
      gate[[basename(e[[2]])]] <<- if (isTRUE(hf)) "HARD" else "SOFT"
    }
  }
  for (k in seq_along(e)[-1]) {
    walk_calls(e[[k]])
  }
}
for (e in parse("analysis/tests/run_all.R", keep.source = FALSE)) {
  walk_calls(e)
}

header_of <- function(x) {
  if (length(x) && grepl("^#!", x[1])) {
    x <- x[-1]
  }
  n <- 0L
  while (n < length(x) && grepl("^#", x[n + 1L])) {
    n <- n + 1L
  }
  x[seq_len(n)]
}

# Numbered scripts: column-0 section lines "# <n>. <Title> ----" padded to 80
# characters (longer titles: at least 4 dashes), the first one "1. Preamble",
# numbers increasing (2, 3, 3a, 3b, 4, ...); no other banner style left.
check_sections <- function(p) {
  x <- readLines(p, warn = FALSE)
  sec_rx <- "^# ([0-9]+)([a-z]?)\\. (\\S.*?) -{4,}$"
  is_sec <- grepl(sec_rx, x, perl = TRUE)
  bad <- character(0)
  old <- grepl("^# ={5,}|^# *Section [0-9]|^# -{3,} .* -{3,}$", x) |
    (grepl("^# .* -{4,}$", x) & !is_sec)
  if (any(old)) {
    bad <- c(bad, "a section banner not in '# <n>. Title ----' form")
  }
  if (!any(is_sec)) {
    return(c(bad, "no numbered sections"))
  }
  s <- x[is_sec]
  if (!grepl("^# 1\\. Preamble -{4,}$", s[1])) {
    bad <- c(bad, "section 1 is not Preamble")
  }
  if (any(nchar(s) != 80L & !grepl(" -{4}$", s))) {
    bad <- c(bad, "section line not padded to 80")
  }
  if (any(nchar(s) > 80L & grepl("-{5,}$", s))) {
    bad <- c(bad, "section line over 80 characters")
  }
  key <- as.integer(sub(sec_rx, "\\1", s, perl = TRUE)) *
    100L +
    match(sub(sec_rx, "\\2", s, perl = TRUE), c("", letters), 0L)
  if (is.unsorted(key, strictly = TRUE)) {
    bad <- c(bad, "section numbers not increasing")
  }
  bad
}

check_file <- function(p) {
  h <- header_of(readLines(p, warn = FALSE))
  own <- basename(p)
  if (!length(h)) {
    return("no header")
  }
  if (!startsWith(h[1], paste0("# ", own, " — "))) {
    return("first line is not '# <file>.R — ...'")
  }
  if (grepl("^analysis/scratchpads/", p)) {
    return(character(0))
  }
  bad <- character(0)
  if (any(grepl("^#'", h))) {
    bad <- c(bad, "roxygen #' in the header")
  }
  if (any(grepl("^#\\s*[=-]{4,}", h))) {
    bad <- c(bad, "rule line in the header")
  }
  if (any(nchar(h) > 80L)) {
    bad <- c(bad, "header line over 80 characters")
  }
  label <- "^# (Usage|Provides|Details):$|^# See: \\S"
  brk <- which(h == "#")
  desc_end <- if (length(brk)) brk[1] - 1L else length(h)
  if (desc_end > 2L || any(grepl(label, h[seq_len(desc_end)]))) {
    bad <- c(bad, "description over 2 lines")
  }
  # blocks: runs between bare "#" lines after the description
  starts <- brk + 1L
  ends <- c(brk[-1] - 1L, length(h))
  blocks <- Map(
    function(a, b) if (a <= b) h[a:b] else character(0),
    starts,
    ends
  )
  labs <- vapply(
    blocks,
    function(b) {
      if (length(b)) {
        sub("^# (Usage|Provides|Details|See):.*$", "\\1", b[1])
      } else {
        ""
      }
    },
    ""
  )
  known <- labs %in% c("Usage", "Provides", "Details", "See")
  if (any(!known)) {
    bad <- c(
      bad,
      "a block without a Usage: / Provides: / Details: / See: label"
    )
  }
  labs <- labs[known]
  blocks <- blocks[known]
  rank <- c(Usage = 1L, Provides = 1L, Details = 2L, See = 3L)
  if (!any(labs %in% c("Usage", "Provides"))) {
    bad <- c(bad, "no Usage: or Provides: block")
  }
  if (any(duplicated(rank[labs])) || is.unsorted(rank[labs])) {
    bad <- c(bad, "blocks out of order or repeated")
  }
  if (grepl("^analysis/tests/test_", p) && !"Usage" %in% labs) {
    bad <- c(bad, "a test needs Usage:")
  }
  if (grepl("^analysis/[0-9]", p) && !"Usage" %in% labs) {
    bad <- c(bad, "a numbered script needs Usage:")
  }
  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    if (labs[k] == "See") {
      if (length(b) != 1L) {
        bad <- c(bad, "See: is one line")
      }
      next
    }
    body <- b[-1]
    if (!length(body)) {
      bad <- c(bad, sprintf("empty %s: block", labs[k]))
    }
    if (any(!grepl("^#   \\S|^#     ", body))) {
      bad <- c(bad, sprintf("%s: lines not indented '#   '", labs[k]))
    }
    if (labs[k] == "Details" && length(body) > 8L) {
      bad <- c(bad, "Details: over 8 lines")
    }
  }
  if (grepl("^analysis/[0-9]", p)) {
    bad <- c(bad, check_sections(p))
  }
  g <- gate[own]
  if (!is.na(g) && !grepl(paste0("^# ", own, " — ", g, ": "), h[1])) {
    bad <- c(bad, sprintf("description does not start with %s:", g))
  }
  unique(bad)
}

res <- lapply(files, check_file)
names(res) <- files
off <- res[lengths(res) > 0L]
grp <- function(p) {
  sub(
    "^analysis/(scratchpads|tests|utils)/.*",
    "\\1",
    sub("^analysis/[^/]+\\.R$", "numbered", p)
  )
}
for (g in c("numbered", "utils", "tests", "scratchpads")) {
  inr <- files[grp(files) == g]
  o <- intersect(names(off), inr)
  expect(
    length(o) == 0L,
    sprintf(
      "%s: %d / %d files follow the template",
      g,
      length(inr) - length(o),
      length(inr)
    )
  )
}
rest <- setdiff(
  names(off),
  files[grp(files) %in% c("numbered", "utils", "tests", "scratchpads")]
)
expect(
  length(rest) == 0L,
  sprintf("other R files follow the template (%d offenders)", length(rest))
)

if (length(off)) {
  why <- table(unlist(off))
  cat("\nby rule:\n")
  for (w in names(sort(why, decreasing = TRUE))) {
    cat(sprintf("  %4d  %s\n", why[[w]], w))
  }
  shown <- if (show_all) names(off) else head(names(off), 15L)
  cat("\noffenders", if (!show_all) "(first 15; --all for every one)", ":\n")
  for (p in shown) {
    cat(sprintf("  %s: %s\n", p, paste(off[[p]], collapse = "; ")))
  }
}
cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
quit(status = as.integer(FAIL > 0L))
