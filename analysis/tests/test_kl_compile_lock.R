# test_kl_compile_lock.R — HARD: dlnm_kernel.R::kl_with_compile_lock()
#   serialises one kl program's cmdstanr compile across processes.
#
# Usage:
#   Rscript analysis/tests/test_kl_compile_lock.R
#
# Details:
#   Quick-safe: no Stan, no data. Array tasks share stan_dir, so the compile
#   is locked and kl_program_file() writes the .stan atomically; forked
#   processes stand in for array tasks.

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")
source("analysis/utils/base.R")
for (m in c("dlnm", "dlnm_kernel")) {
  source(sprintf("analysis/utils/%s.R", m))
}

cat("# kl compile lock\n")
td <- tempfile("kllock_")
dir.create(td)
prog <- file.path(td, "occ_kl_latent_x_thr.stan")
writeLines("model {}", prog)
log <- file.path(td, "log.txt")

# [1] two processes holding the lock never overlap (no interleaved markers)
cat("\n[1] two processes never overlap inside the lock\n")
work <- function(tag) {
  kl_with_compile_lock(
    prog,
    {
      cat("start", tag, "\n", file = log, append = TRUE)
      Sys.sleep(1.5)
      cat("end", tag, "\n", file = log, append = TRUE)
      tag
    },
    poll_sec = 0.2
  )
}
jobs <- lapply(c("a", "b"), function(tg) parallel::mcparallel(work(tg)))
res <- parallel::mccollect(jobs)
ev <- trimws(readLines(log))
expect(
  length(ev) == 4L &&
    all(sub(" .*", "", ev) == c("start", "end", "start", "end")) &&
    sub(".* ", "", ev[1]) == sub(".* ", "", ev[2]),
  sprintf(
    "markers are start/end pairs, never interleaved (%s)",
    paste(ev, collapse = " | ")
  )
)
expect(
  setequal(unlist(res), c("a", "b")),
  "both processes returned their value"
)

# [2] the lock directory is removed on exit, also when expr errors
cat("\n[2] the lock is released, also on error\n")
lockdir <- paste0(prog, ".lock")
expect(!dir.exists(lockdir), "no lock left after two clean runs")
e <- try(kl_with_compile_lock(prog, stop("compile failed")), silent = TRUE)
expect(
  inherits(e, "try-error") && !dir.exists(lockdir),
  "an error inside the lock still removes it"
)

# [3] a stale lock (older than stale_min) is broken; a live one times out
cat("\n[3] stale and live locks\n")
dir.create(lockdir)
Sys.setFileTime(lockdir, Sys.time() - 3 * 3600)
expect(
  identical(
    kl_with_compile_lock(prog, "ok", stale_min = 60, poll_sec = 0.1),
    "ok"
  ),
  "a lock older than stale_min is broken"
)
dir.create(lockdir)
e2 <- try(
  kl_with_compile_lock(
    prog,
    "never",
    wait_min = 0.02,
    stale_min = 60,
    poll_sec = 0.2
  ),
  silent = TRUE
)
expect(
  inherits(e2, "try-error") && grepl("waited", e2),
  "a live lock times out after wait_min"
)
unlink(lockdir, recursive = TRUE)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
