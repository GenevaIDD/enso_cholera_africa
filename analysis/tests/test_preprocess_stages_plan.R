# test_preprocess_stages_plan.R — HARD: the 00-02 data leg's Slurm drivers
#   plan the right chain and submit nothing under --plan.
#
# Usage:
#   Rscript analysis/tests/test_preprocess_stages_plan.R
#
# Details:
#   campaigns/preprocess.sh is the one command that rebuilds a new outbreak
#   vintage on the cluster; a wrong dependency edge lets a step read the
#   previous vintage's artefact. Checks the campaign chain, slices, each
#   stage driver's --plan, utils/preprocess_paths.R and the 03 stage plan.
#   Committed deploy / population configs are only read. No data, no
#   cluster: a stub sbatch first on PATH records any call and fails.

suppressPackageStartupMessages(library(yaml))

PASS <- 0L
FAIL <- 0L
source("analysis/tests/helpers.R")

cat("# 00-02 data leg on SLURM: campaign and stage plans\n")

CAMPAIGN <- "analysis/bash/campaigns/preprocess.sh"
STAGES <- "analysis/bash/stages"
DEP <- "analysis/configs/deploy/prod_run.yml"
POPCFG <- "analysis/configs/pop_extract/pop_extract.yml"

STUB_DIR <- tempfile("stub_bin_")
dir.create(STUB_DIR)
MARKER <- file.path(STUB_DIR, "sbatch_called")
writeLines(
  c(
    "#!/bin/bash",
    sprintf("echo \"$*\" >> '%s'", MARKER),
    "echo 'STUB sbatch: a --plan must never submit' >&2",
    "exit 97"
  ),
  file.path(STUB_DIR, "sbatch")
)
Sys.chmod(file.path(STUB_DIR, "sbatch"), "0755")

run_bash <- function(script, args) {
  out <- suppressWarnings(system2(
    "bash",
    c(script, args),
    stdout = TRUE,
    stderr = TRUE,
    env = sprintf("PATH=%s:%s", STUB_DIR, Sys.getenv("PATH"))
  ))
  st <- attr(out, "status")
  list(out = out, status = if (is.null(st)) 0L else st)
}
sb <- function(r, job) {
  grep(sprintf("--job-name=%s", job), r$out, fixed = TRUE, value = TRUE)
}
has <- function(r, pattern) any(grepl(pattern, r$out, fixed = TRUE))
dep_of <- function(line) {
  m <- regmatches(line, regexpr("--dependency=afterok:[^ ]+", line))
  if (!length(m)) "" else sub("--dependency=afterok:", "", m, fixed = TRUE)
}

# --- [1] the whole chain --------------------------------------------------------
# [1] the campaign's --plan chains 00 -> 01 (+01b, 01_nz) and 00 -> 02
#     population (array -> 02c), 02 crosswalk afterok 01 + 02c, 02 climate
#     afterok the crosswalk; nothing reaches sbatch (a stub first on PATH
#     records any call and fails).
cat("\n[1] campaign --plan: the full 00 -> 02 chain, nothing submitted\n")
r <- run_bash(CAMPAIGN, c("--config", POPCFG, "--deploy", DEP, "--plan"))
expect(r$status == 0L, "campaign --plan exits 0")
expect(!file.exists(MARKER), "campaign --plan never calls sbatch")

l00 <- sb(r, "prep_00_ingest")
l01 <- sb(r, "prep_01_outbreaks")
l1b <- sb(r, "prep_01b_monthly")
lnz <- sb(r, "prep_01_national_zero")
lxw <- sb(r, "prep_02_crosswalk")
lcl <- sb(r, "prep_02_climate_all")
lpa <- grep("submit_population_country.sh", r$out, fixed = TRUE, value = TRUE)
lpc <- grep(
  "02c_preprocess_population_compile.R",
  r$out,
  fixed = TRUE,
  value = TRUE
)
expect(
  all(lengths(list(l00, l01, l1b, lnz, lxw, lcl, lpa, lpc)) == 1L),
  "one sbatch line per job: 00, 01, 01b, 01_nz, population array, 02c, crosswalk, climate"
)
expect(dep_of(l00) == "", "00 waits on nothing")
expect(dep_of(l01) == "PLAN_00_ingest", "01 waits on 00")
expect(
  dep_of(l1b) == "PLAN_prep_01_outbreaks" &&
    dep_of(lnz) == "PLAN_prep_01_outbreaks",
  "01b and 01_nz wait on 01 (they read its outbreak metadata)"
)
expect(
  dep_of(lpa) == "PLAN_00_ingest",
  "the population array waits on 00 (LP shapefiles), not on 01"
)
expect(dep_of(lpc) == "PLAN_pop_array", "02c waits on the population array")
expect(
  setequal(
    strsplit(dep_of(lxw), ":")[[1]],
    c("PLAN_prep_01_outbreaks", "PLAN_pop_compile")
  ),
  "the crosswalk waits on 01 (gpkg) and 02c (population crosswalk)"
)
expect(
  "PLAN_02_crosswalk" %in% strsplit(dep_of(lcl), ":")[[1]],
  "climate waits on the crosswalk (u_gids_crosswalk.rds)"
)
expect(
  grepl("-f FALSE", l01, fixed = TRUE),
  "01 runs without figures on the cluster"
)
expect(
  grepl(
    "analysis/02_preprocess_climate.R --stage all -s TRUE",
    lcl,
    fixed = TRUE
  ),
  "climate runs --stage all with -s TRUE"
)
expect(
  all(grepl(
    "analysis/bash/slurm/submit_preprocess.sh",
    c(l00, l01, l1b, lnz, lxw, lcl, lpc),
    fixed = TRUE
  )),
  "single jobs go through slurm/submit_preprocess.sh"
)

# --- [2] the deploy's national-zero rule ----------------------------------------
# [2] 01_nz builds the deploy's deploy.episodes rule and threshold.
cat("\n[2] 01_nz builds the deploy's national-zero months\n")
ep <- yaml::read_yaml(DEP)$deploy$episodes
expect(
  grepl(
    sprintf(
      "--rule %s --zero_threshold %d",
      ep$national_zero_rule,
      as.integer(ep$national_zero_threshold)
    ),
    lnz,
    fixed = TRUE
  ),
  sprintf(
    "01_nz gets --rule %s --zero_threshold %s from %s",
    ep$national_zero_rule,
    ep$national_zero_threshold,
    basename(DEP)
  )
)
rx <- run_bash(
  file.path(STAGES, "01_preprocess.sh"),
  c("--nz-rule", "r9", "--nz-threshold", "7", "--plan")
)
expect(
  rx$status == 0L &&
    any(grepl("--rule r9 --zero_threshold 7", rx$out, fixed = TRUE)),
  "explicit --nz-rule / --nz-threshold reach 01_nz"
)
rn <- run_bash(file.path(STAGES, "01_preprocess.sh"), "--plan")
expect(
  rn$status != 0L,
  "01_preprocess.sh without --deploy or --nz-* is refused"
)

# --- [3] slices -----------------------------------------------------------------
# [3] --from / --until cut the chain: a step outside the slice is neither
#     planned nor waited on; bad slices and a missing --config are refused.
cat("\n[3] --from / --until\n")
rf <- run_bash(CAMPAIGN, c("--config", POPCFG, "--from", "02xw", "--plan"))
expect(
  rf$status == 0L,
  "--from 02xw --plan exits 0 (no --deploy needed past 01)"
)
expect(
  !length(sb(rf, "prep_00_ingest")) &&
    !length(sb(rf, "prep_01_outbreaks")) &&
    !any(
      grepl("population", rf$out, fixed = TRUE) &
        grepl("--array", rf$out, fixed = TRUE)
    ),
  "--from 02xw plans neither 00, 01 nor the population array"
)
xw <- sb(rf, "prep_02_crosswalk")
expect(
  length(xw) == 1L && dep_of(xw) == "",
  "--from 02xw: the crosswalk waits on nothing"
)
ru <- run_bash(CAMPAIGN, c("--deploy", DEP, "--until", "01", "--plan"))
expect(
  ru$status == 0L &&
    length(sb(ru, "prep_01_outbreaks")) == 1L &&
    !length(sb(ru, "prep_02_crosswalk")),
  "--until 01 needs no --config and plans nothing of 02"
)
expect(
  run_bash(
    CAMPAIGN,
    c(
      "--config",
      POPCFG,
      "--deploy",
      DEP,
      "--from",
      "02clim",
      "--until",
      "00",
      "--plan"
    )
  )$status !=
    0L,
  "--from after --until is refused"
)
expect(
  run_bash(
    CAMPAIGN,
    c("--config", POPCFG, "--deploy", DEP, "--from", "03", "--plan")
  )$status !=
    0L,
  "an unknown step is refused"
)
expect(
  run_bash(CAMPAIGN, c("--deploy", DEP, "--plan"))$status != 0L,
  "02 population / crosswalk without --config is refused"
)

# --- [4] stage drivers alone --------------------------------------------------
# [4] every stage driver's --plan exits 0 and prints its *_JOB_ID= line;
#     02_population.sh picks the payload of --task-unit and refuses --years
#     with the per-country unit.
cat("\n[4] each stage driver's --plan\n")
for (st in list(
  list("00_ingest.sh", character(0), "INGEST_JOB_ID="),
  list("01_preprocess.sh", c("--deploy", DEP), "NZ_JOB_ID="),
  list("02_population.sh", c("--config", POPCFG), "POP_COMPILE_JOB_ID="),
  list("02_crosswalk.sh", c("--config", POPCFG), "CROSSWALK_JOB_ID="),
  list("02_climate.sh", c("--stage", "match"), "CLIMATE_JOB_ID=")
)) {
  rs <- run_bash(file.path(STAGES, st[[1]]), c(st[[2]], "--plan"))
  expect(
    rs$status == 0L && has(rs, st[[3]]),
    sprintf("%s --plan exits 0 and prints %s", st[[1]], st[[3]])
  )
}
rcy <- run_bash(
  file.path(STAGES, "02_population.sh"),
  c("--config", POPCFG, "--task-unit", "country_year", "--plan")
)
expect(
  has(rcy, "slurm/submit_population.sh") &&
    has(rcy, "--task_unit country_year"),
  "--task-unit country_year plans the per-year payload and 02a --task_unit country_year"
)
expect(
  run_bash(
    file.path(STAGES, "02_population.sh"),
    c("--config", POPCFG, "--years", "2020", "--plan")
  )$status !=
    0L,
  "--years with the per-country unit is refused"
)
expect(
  run_bash(
    file.path(STAGES, "02_climate.sh"),
    c("--stage", "nope", "--plan")
  )$status !=
    0L,
  "02_climate.sh refuses an unknown --stage"
)
expect(!file.exists(MARKER), "no stage --plan called sbatch")

# --- [5] paths from R ---------------------------------------------------------------
# [5] utils/preprocess_paths.R names the vintage folder from paths.R and the
#     national-zero files of the deploy; no driver carries a vintage literal.
cat("\n[5] utils/preprocess_paths.R\n")
source("analysis/utils/paths.R")
kv <- suppressWarnings(system2(
  "Rscript",
  c(
    "analysis/utils/preprocess_paths.R",
    "-C",
    POPCFG,
    "--deploy",
    DEP,
    "--outputs"
  ),
  stdout = TRUE,
  stderr = FALSE
))
val <- function(k) {
  sub(paste0("^", k, "="), "", grep(paste0("^", k, "="), kv, value = TRUE))
}
expect(
  identical(val("VINTAGE_DIR"), path.expand(default_outbreak_data_dir())),
  "VINTAGE_DIR is paths.R::default_outbreak_data_dir()"
)
bare <- suppressWarnings(system2(
  "Rscript",
  "analysis/utils/preprocess_paths.R",
  stdout = TRUE,
  stderr = FALSE
))
expect(
  is.null(attr(bare, "status")) &&
    any(grepl("^VINTAGE_DIR=", bare)) &&
    !any(grepl("^NZ_", bare)),
  "without --deploy (the rsync-inputs call) it runs and prints no NZ_ keys"
)
expect(
  identical(val("NZ_RULE"), ep$national_zero_rule) &&
    identical(val("NZ_THRESHOLD"), as.character(ep$national_zero_threshold)),
  "NZ_RULE / NZ_THRESHOLD follow the deploy"
)
nz <- as.character(make_national_zero_months_file(
  as.integer(ep$national_zero_threshold),
  ep$national_zero_rule
))
expect(
  nz %in% val("SEC") && make_outbreak_metadata_file() %in% val("SEC"),
  "--outputs lists the deploy's national-zero months and the outbreak metadata"
)
new_bash <- c(
  CAMPAIGN,
  file.path(
    STAGES,
    c(
      "00_ingest.sh",
      "01_preprocess.sh",
      "02_crosswalk.sh",
      "02_climate.sh",
      "02_population.sh"
    )
  ),
  "analysis/bash/slurm/submit_preprocess.sh",
  "analysis/bash/slurm/submit_population_country.sh",
  # the data-leg copy tools, where the checkout has them (a release has none)
  list.files(
    "analysis/bash/tools",
    pattern = "^rsync_preprocess_.*\\.sh$",
    full.names = TRUE
  )
)
lit <- vapply(
  new_bash,
  function(f) {
    any(grepl(
      "outbreak_data_[A-Z][a-z]{2}_[0-9]{4}",
      readLines(f, warn = FALSE)
    ))
  },
  logical(1)
)
expect(!any(lit), "no data-leg driver carries an outbreak-vintage literal")

# --- [6] the 03 stage ------------------------------------------------------------
# [6] stages/03_prepare.sh --plan: the 03 chain after --after (default occ
#     panel -> the deploy's panel -> siz -> dur -> mon -> sim panels), WHO
#     monthly and ADM1 rainfall on --after alone; the deploy's panel spec and
#     projection flags reach the scripts; every panel step has --force.
cat("\n[6] stages/03_prepare.sh --plan\n")
r3 <- run_bash(
  file.path(STAGES, "03_prepare.sh"),
  c(
    "--deploy",
    DEP,
    "--after",
    "11:22",
    "--extended-end",
    "2026-08-01",
    "--plan"
  )
)
expect(r3$status == 0L, "03_prepare.sh --plan exits 0")
j3 <- c(
  "prep_03_occ_default",
  "prep_03_occ_deploy",
  "prep_03_siz",
  "prep_03_dur",
  "prep_03_mon",
  "prep_03_sim_panel",
  "prep_03_sim_projection",
  "prep_03_sim_extended",
  "prep_03_who_monthly",
  "prep_03_rain_adm1"
)
l3 <- lapply(j3, function(j) sb(r3, paste0(j, " ")))
expect(all(lengths(l3) == 1L), "one sbatch line per 03 step (10)")
names(l3) <- j3
chain <- c(
  "prep_03_occ_default",
  "prep_03_occ_deploy",
  "prep_03_siz",
  "prep_03_dur",
  "prep_03_mon",
  "prep_03_sim_panel",
  "prep_03_sim_projection",
  "prep_03_sim_extended"
)
ok_chain <- dep_of(l3[[chain[1]]]) == "11:22" &&
  all(vapply(
    2:length(chain),
    function(i) dep_of(l3[[chain[i]]]) == paste0("PLAN_", chain[i - 1L]),
    logical(1)
  ))
expect(ok_chain, "the panel steps run one after another, the first on --after")
expect(
  dep_of(l3$prep_03_who_monthly) == "11:22" &&
    dep_of(l3$prep_03_rain_adm1) == "11:22",
  "WHO monthly and ADM1 rainfall wait on --after only"
)
ep_dep <- yaml::read_yaml(DEP)$deploy$episodes
lo <- l3$prep_03_occ_deploy
expect(
  grepl(
    sprintf(
      "--national_zero_threshold %s --national_zero_rule %s",
      ep_dep$national_zero_threshold,
      ep_dep$national_zero_rule
    ),
    lo,
    fixed = TRUE
  ) &&
    grepl(
      sprintf(
        "--at_risk_start %s --adm1_overlap %s",
        ep_dep$at_risk_start,
        ep_dep$adm1_overlap
      ),
      lo,
      fixed = TRUE
    ) &&
    grepl("--contagion_coverage_report", lo, fixed = TRUE),
  "the deploy's occurrence panel carries deploy.episodes and the coverage report"
)
expect(
  !grepl("national_zero", l3$prep_03_occ_default, fixed = TRUE),
  "the default occurrence panel has no national zeros"
)
lp <- l3$prep_03_sim_projection
ec <- yaml::read_yaml(DEP)$deploy$enso_continuation
expect(
  grepl(sprintf("--enso_continuation %s", ec$method), lp, fixed = TRUE) &&
    grepl(sprintf("--enso_continuation_seam %s", ec$seam), lp, fixed = TRUE) &&
    grepl("--cfsv2_correction calibration", lp, fixed = TRUE),
  "the projection panel follows deploy.enso_continuation and the forecast's CFSv2 calibration"
)
expect(
  grepl("--end 2026-08-01", l3$prep_03_sim_extended, fixed = TRUE),
  "--extended-end reaches the extended panel"
)
expect(
  all(grepl(
    "--force",
    unlist(l3[setdiff(j3, "prep_03_rain_adm1")]),
    fixed = TRUE
  )),
  "every guarded 03 step runs with --force"
)
expect(
  has(r3, "PREPARE_JOB_IDS=PLAN_prep_03_occ_default:"),
  "PREPARE_JOB_IDS lists the 03 jobs"
)
r3n <- run_bash(
  file.path(STAGES, "03_prepare.sh"),
  c("--deploy", DEP, "--plan")
)
expect(
  r3n$status == 0L &&
    !length(sb(r3n, "prep_03_sim_extended ")) &&
    dep_of(sb(r3n, "prep_03_occ_default ")) == "",
  "without --extended-end the extended panel is skipped; without --after the first step waits on nothing"
)
r3b <- run_bash(
  file.path(STAGES, "03_prepare.sh"),
  c("--deploy", DEP, "--extended-end", "2026-08", "--plan")
)
expect(r3b$status != 0L, "a malformed --extended-end is refused")
expect(!file.exists(MARKER), "03_prepare.sh --plan never calls sbatch")

cat(sprintf("\n%d pass, %d fail\n", PASS, FAIL))
if (FAIL > 0) {
  quit(status = 1)
}
