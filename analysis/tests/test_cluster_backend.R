# test_cluster_backend.R — HARD: analysis/bash/lib/cluster.sh and the local
#   backend.
#
# Usage:
#   Rscript analysis/tests/test_cluster_backend.R
#
# Details:
#   Every driver and payload reaches the cluster through lib/cluster.sh; a
#   laptop run goes through lib/local/sbatch. Checks cholera_sbatch_cmd [1],
#   cholera_env [2], CHOLERA_BACKEND=local [3] and the sbatch shim [4]. Temp
#   profiles and payloads only: no R model, no cluster.

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# lib/cluster.sh + lib/local/sbatch\n")

LIB <- normalizePath("analysis/bash/lib/cluster.sh")
SHIM <- normalizePath("analysis/bash/lib/local/sbatch")
TMP <- tempfile("cluster_backend_")
dir.create(TMP)

sh <- function(code, env = character()) {
  out <- suppressWarnings(system2(
    "bash",
    c("-c", shQuote(code)),
    stdout = TRUE,
    stderr = TRUE,
    env = env
  ))
  st <- attr(out, "status")
  list(out = out, status = if (is.null(st)) 0L else st)
}
# The assembled sbatch command of `args` (words of a cholera_sbatch_cmd call).
cmd_of <- function(args, profile = NULL, env = character()) {
  pe <- if (is.null(profile)) {
    sprintf("CHOLERA_CLUSTER_PROFILE=%s", file.path(TMP, "empty.profile"))
  } else {
    sprintf("CHOLERA_CLUSTER_PROFILE=%s", profile)
  }
  r <- sh(
    sprintf(
      'cd %s && source %s && cholera_sbatch_cmd %s && echo "CMD: ${CHOLERA_SBATCH_CMD[*]}"',
      shQuote(TMP),
      shQuote(LIB),
      args
    ),
    env = c(pe, env)
  )
  list(
    status = r$status,
    cmd = sub("^CMD: ", "", grep("^CMD: ", r$out, value = TRUE)),
    out = r$out
  )
}
words <- function(cmd) strsplit(cmd, " ", fixed = TRUE)[[1]]
has_w <- function(cmd, w) w %in% words(cmd)

writeLines("", file.path(TMP, "empty.profile"))
PROF <- file.path(TMP, "test.profile")
writeLines(
  c(
    'CHOLERA_PARTITION_SHORT="short1,short2"',
    'CHOLERA_PARTITION_LONG="long"',
    'CHOLERA_PARTITION_BIGMEM="big"',
    'CHOLERA_ACCOUNT="acct"',
    'CHOLERA_QOS_LONG="slowq"',
    'CHOLERA_MEM_BIGMEM="500G"',
    'CHOLERA_RES_LOO_COMPARE="--partition=only --mem=7G"',
    'CHOLERA_RES_RELAY="--mem=8G"',
    'CHOLERA_ENV_PURGE="echo PURGED"',
    'CHOLERA_ENV_SETUP="echo BASE_STACK"',
    'CHOLERA_ENV_SETUP_GEO="echo GEO_STACK"',
    'CHOLERA_LOG_DIR="mylogs"'
  ),
  PROF
)

# --- [1] cholera_sbatch_cmd ------------------------------------------------------
# cholera_sbatch_cmd under a temp profile: the class partition, account and
# qos; CHOLERA_RES_<KEY> above --defaults above the class; options the caller
# passes win over all three; --output in CHOLERA_LOG_DIR (%A_%a for arrays);
# the key of a payload; a bad class
# is refused; no profile adds no partition.
cat("\n[1] cholera_sbatch_cmd flag assembly\n")
a <- cmd_of("SHORT --parsable --job-name=x p/submit_results.sh A B", PROF)
expect(a$status == 0L && length(a$cmd) == 1L, "SHORT assembles one command")
expect(
  has_w(a$cmd, "--partition=short1,short2") &&
    has_w(a$cmd, "--account=acct") &&
    !any(grepl("^--qos", words(a$cmd))),
  "SHORT: the class partition, the global account, no qos (only LONG sets one)"
)
expect(
  has_w(a$cmd, "--output=mylogs/%x_%j.log"),
  "a single job logs to CHOLERA_LOG_DIR/%x_%j.log"
)
expect(
  grepl("p/submit_results.sh A B$", a$cmd),
  "the payload and its arguments come last, unchanged"
)
l <- cmd_of("LONG --array=0-3%2 p/submit_fit.sh b", PROF)
expect(
  has_w(l$cmd, "--partition=long") &&
    has_w(l$cmd, "--qos=slowq") &&
    has_w(l$cmd, "--output=mylogs/%x_%A_%a.log"),
  "LONG: its partition and qos; an array logs to %x_%A_%a"
)
b <- cmd_of("BIGMEM p/submit_smoothing.sh c", PROF)
expect(
  has_w(b$cmd, "--mem=500G") && has_w(b$cmd, "--partition=big"),
  "BIGMEM: the class memory and partition"
)
r <- cmd_of("SHORT --array=0-1 p/submit_loo_compare.sh j", PROF)
expect(
  has_w(r$cmd, "--partition=only") &&
    has_w(r$cmd, "--mem=7G") &&
    !has_w(r$cmd, "--partition=short1,short2"),
  "CHOLERA_RES_LOO_COMPARE replaces the class partition and adds its memory"
)
e <- cmd_of(
  "SHORT --partition=mine --mem 3G -c 2 --output=o.log p/submit_loo_compare.sh j",
  PROF
)
expect(
  has_w(e$cmd, "--partition=mine") &&
    !has_w(e$cmd, "--partition=only") &&
    !has_w(e$cmd, "--mem=7G") &&
    !any(grepl("^--output=mylogs", words(e$cmd))),
  "options the caller passes (--partition=, --mem X, -c, --output) win over the profile"
)
d <- cmd_of(
  "SHORT --res RELAY --defaults '--mem=1G --cpus-per-task=1 --time=00:10:00' --wrap 'echo hi'",
  PROF
)
expect(
  has_w(d$cmd, "--mem=8G") &&
    !has_w(d$cmd, "--mem=1G") &&
    has_w(d$cmd, "--cpus-per-task=1") &&
    has_w(d$cmd, "--time=00:10:00"),
  "--res RELAY: CHOLERA_RES_RELAY above --defaults, the other defaults kept"
)
d2 <- cmd_of(
  "SHORT --res RELAY --defaults '--mem=1G' --time=01:00:00 --wrap 'echo hi'",
  PROF
)
expect(
  !any(grepl("^--time=00", words(d2$cmd))) && has_w(d2$cmd, "--time=01:00:00"),
  "the caller's --time wins over --defaults"
)
n <- cmd_of("SHORT p/submit_results.sh")
expect(
  n$status == 0L &&
    !any(grepl("^--(partition|account|qos|mem)", words(n$cmd))) &&
    has_w(n$cmd, "--output=logs/%x_%j.log"),
  "an empty profile adds no partition / account / qos / memory, logs to logs/"
)
expect(
  cmd_of("HUGE p/submit_results.sh", PROF)$status != 0L,
  "an unknown class is refused"
)
miss <- sh(
  sprintf("source %s", shQuote(LIB)),
  env = "CHOLERA_CLUSTER_PROFILE=/no/such/profile"
)
expect(
  miss$status != 0L,
  "a CHOLERA_CLUSTER_PROFILE that does not exist is refused"
)

# --- [2] cholera_env ----------------------------------------------------------------
# cholera_env: the named stack, the fallback to CHOLERA_ENV_SETUP,
# CHOLERA_ENV_PURGE first with --purge.
cat("\n[2] cholera_env\n")
ev <- function(args) {
  sh(
    sprintf("source %s && cholera_env %s", shQuote(LIB), args),
    env = sprintf("CHOLERA_CLUSTER_PROFILE=%s", PROF)
  )$out
}
expect(identical(ev(""), "BASE_STACK"), "no stack: CHOLERA_ENV_SETUP")
expect(identical(ev("GEO"), "GEO_STACK"), "GEO: CHOLERA_ENV_SETUP_GEO")
expect(
  identical(ev("REPORTS"), "BASE_STACK"),
  "a stack the profile leaves unset falls back to CHOLERA_ENV_SETUP"
)
expect(
  identical(ev("--purge GEO"), c("PURGED", "GEO_STACK")),
  "--purge runs CHOLERA_ENV_PURGE first"
)

# --- [3] backend ----------------------------------------------------------------------
# CHOLERA_BACKEND=local puts lib/local first on PATH.
cat("\n[3] CHOLERA_BACKEND=local\n")
pth <- sh(
  sprintf("source %s && command -v sbatch", shQuote(LIB)),
  env = c("CHOLERA_BACKEND=local", sprintf("CHOLERA_CLUSTER_PROFILE=%s", PROF))
)
expect(
  pth$status == 0L && identical(pth$out, SHIM),
  "local: lib/local/sbatch is the sbatch on PATH"
)
expect(
  sh(sprintf("source %s", shQuote(LIB)), env = "CHOLERA_BACKEND=pbs")$status !=
    0L,
  "an unknown backend is refused"
)

# --- [4] the shim ------------------------------------------------------------------------
# lib/local/sbatch: a single job (the payload's #SBATCH job name and cpus,
# capped at CHOLERA_CORES), --parsable, "Submitted batch job N", --wrap,
# --export, --dependency ignored; arrays (ranges, lists, :step, %cap) run every
# task once with SLURM_ARRAY_TASK_ID / SLURM_CPUS_PER_TASK; a failed task makes
# sbatch exit non-zero; job ids increase.
cat("\n[4] lib/local/sbatch\n")
writeLines(
  c(
    "#!/bin/bash",
    "#SBATCH --job-name=demo",
    "#SBATCH -c 8",
    "echo \"task=${SLURM_ARRAY_TASK_ID:-none} cpus=${SLURM_CPUS_PER_TASK} job=${SLURM_JOB_NAME} extra=${EXTRA_VAR:-} args=$*\"",
    "echo \"${SLURM_ARRAY_TASK_ID:-none}\" >> ran.txt",
    "[[ \"${SLURM_ARRAY_TASK_ID:-}\" == 13 ]] && exit 3",
    "exit 0"
  ),
  file.path(TMP, "payload.sh")
)
run_shim <- function(args, cores = 4L) {
  sh(
    sprintf("cd %s && %s %s", shQuote(TMP), shQuote(SHIM), args),
    env = c(sprintf("CHOLERA_CORES=%d", cores), "CHOLERA_LOG_DIR=jl")
  )
}
logs <- function() {
  list.files(file.path(TMP, "jl"), pattern = "\\.log$", full.names = TRUE)
}
ran <- function() {
  f <- file.path(TMP, "ran.txt")
  x <- if (file.exists(f)) readLines(f) else character()
  unlink(f)
  x
}

s1 <- run_shim(
  "--parsable --dependency=afterok:999 --export=ALL,EXTRA_VAR=yes payload.sh a b"
)
id1 <- suppressWarnings(as.integer(s1$out))
expect(
  s1$status == 0L && length(s1$out) == 1L && !is.na(id1),
  "--parsable prints the bare job id; --dependency is ignored"
)
lg1 <- grep(sprintf("demo_%d\\.log$", id1), logs(), value = TRUE)
l1 <- if (length(lg1)) readLines(lg1) else ""
expect(
  length(lg1) == 1L &&
    identical(l1, "task=none cpus=4 job=demo extra=yes args=a b"),
  "the payload's #SBATCH job name; its -c 8 capped at CHOLERA_CORES 4; --export reaches the job; log <name>_<id>.log"
)
expect(identical(ran(), "none"), "a single job runs once")
s2 <- run_shim("payload.sh")
id2 <- suppressWarnings(as.integer(sub("^Submitted batch job ", "", s2$out)))
expect(
  s2$status == 0L &&
    grepl("^Submitted batch job [0-9]+$", s2$out) &&
    isTRUE(id2 > id1),
  "without --parsable: 'Submitted batch job N', ids increase"
)
invisible(ran())
s3 <- run_shim(
  "--parsable --array=0-4:2,7,9-10%2 --job-name=arr -c 1 payload.sh x"
)
expect(s3$status == 0L, "an array whose tasks succeed exits 0")
expect(
  setequal(ran(), c("0", "2", "4", "7", "9", "10")),
  "ranges, :step and lists expand to tasks 0 2 4 7 9 10, each run once"
)
al <- grep("/arr_[0-9]+_[0-9]+\\.log$", logs(), value = TRUE)
expect(
  length(al) == 6L &&
    all(vapply(
      al,
      function(f) grepl("cpus=1 job=arr", readLines(f)[1]),
      logical(1)
    )),
  "one log per task, --job-name and -c 1 honoured"
)
s4 <- run_shim("--array=0-3 payload.sh", cores = 2L)
c4 <- vapply(
  grep("/demo_[0-9]+_[0-9]+\\.log$", logs(), value = TRUE),
  function(f) sub(".*cpus=([0-9]+).*", "\\1", readLines(f)[1]),
  character(1)
)
expect(
  s4$status == 0L && all(c4 == "1"),
  "2 cores, 4 tasks of -c 8: two at a time, SLURM_CPUS_PER_TASK 1"
)
invisible(ran())
s5 <- run_shim("--array=12-14 payload.sh")
expect(
  s5$status != 0L && any(grepl("task 13 FAILED", s5$out, fixed = TRUE)),
  "a failed array task: non-zero exit, the task named on stderr"
)
expect(setequal(ran(), c("12", "13", "14")), "the other tasks still ran")
s6 <- run_shim("--parsable --wrap 'echo wrapped; touch wrap_ran'")
expect(
  s6$status == 0L && file.exists(file.path(TMP, "wrap_ran")),
  "--wrap runs its command"
)
expect(
  run_shim("--wrap 'exit 4'")$status == 4L,
  "a failed single job returns its exit status"
)

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0L) {
  quit(status = 1L)
}
