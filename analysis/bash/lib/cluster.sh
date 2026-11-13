# lib/cluster.sh — the one place that knows the cluster. Sourced (never run)
# by every stage, campaign, tool and slurm payload under analysis/bash/.
#
#   source analysis/bash/lib/cluster.sh      # payloads (cwd = repo root: SLURM copies the script)
#   source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/cluster.sh"   # drivers
#
# Reads the cluster profile ${CHOLERA_CLUSTER_PROFILE:-analysis/bash/cluster.profile}
# (a bash file of variables; analysis/bash/cluster.profile.example lists
# them). No profile = no environment
# setup and no partition: sbatch uses the cluster's default partition and each
# payload's #SBATCH mem/time/cpus.
#
# Exposes:
#   cholera_env [--purge] [STACK]
#                         software environment of a stack: evals $CHOLERA_ENV_PURGE
#                         first with --purge (e.g. module purge), then
#                         $CHOLERA_ENV_SETUP_<STACK> (else $CHOLERA_ENV_SETUP;
#                         default nothing). Stacks in use: (none) = fits and
#                         simulator R, GEO = sf / GDAL R (00-03, 09, population),
#                         MANUSCRIPT = GEO + quarto, REPORTS = GEO + quarto + LibreOffice.
#   cholera_sbatch CLASS [--res KEY] [--defaults "OPTS"] [sbatch args] [payload [payload args]]
#                         sbatch with the profile's resources added. CLASS is
#                         SHORT (<= 12 h), LONG (days) or BIGMEM (large memory).
#                         Adds, unless the caller already passed the option:
#                           --partition   $CHOLERA_PARTITION_<CLASS>
#                           --account     $CHOLERA_ACCOUNT_<CLASS> | $CHOLERA_ACCOUNT
#                           --qos         $CHOLERA_QOS_<CLASS> | $CHOLERA_QOS
#                           --mem / --time / --cpus-per-task
#                                         $CHOLERA_{MEM,TIME,CPUS}_<CLASS> (unset by
#                                         default: the payload's #SBATCH lines apply)
#                           --output      $CHOLERA_LOG_DIR/%x_%j.log (%x_%A_%a.log for arrays)
#                         Per-payload overrides: $CHOLERA_RES_<KEY>, sbatch options
#                         (word-split) that win over the class defaults, e.g.
#                         CHOLERA_RES_LOO_COMPARE="--partition=short". KEY is
#                         --res KEY, else the payload's basename without the
#                         submit_ prefix and .sh, upper-cased
#                         (slurm/submit_loo_compare.sh -> LOO_COMPARE).
#                         --defaults: the driver's own sizing of this job (sbatch
#                         options as --opt=value words in one string, e.g. "--mem=32G
#                         --time=04:00:00"), below CHOLERA_RES_<KEY>.
#                         Precedence: caller's args > CHOLERA_RES_<KEY> > --defaults
#                         > class > #SBATCH header.
#   cholera_sbatch_cmd CLASS ...   same, but only fills the array CHOLERA_SBATCH_CMD
#                         (for --plan / --dry-run printing).
#
# Backends: CHOLERA_BACKEND=slurm (default) or local. local puts
# analysis/bash/lib/local/ first on PATH: its `sbatch` runs every job on this
# machine, synchronously ($CHOLERA_CORES array tasks at a time).

if [[ -z "${CHOLERA_CLUSTER_LIB_LOADED:-}" ]]; then
CHOLERA_CLUSTER_LIB_LOADED=1

CHOLERA_BASH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

_cholera_profile="${CHOLERA_CLUSTER_PROFILE:-$CHOLERA_BASH_DIR/cluster.profile}"
if [[ -n "${CHOLERA_CLUSTER_PROFILE:-}" && ! -f "$_cholera_profile" ]]; then
  echo "Error: CHOLERA_CLUSTER_PROFILE=$CHOLERA_CLUSTER_PROFILE not found." >&2
  return 1 2>/dev/null || exit 1
fi
# shellcheck disable=SC1090
[[ -f "$_cholera_profile" ]] && source "$_cholera_profile"
unset _cholera_profile

: "${CHOLERA_BACKEND:=slurm}"
: "${CHOLERA_LOG_DIR:=logs}"
export CHOLERA_BACKEND CHOLERA_LOG_DIR
case "$CHOLERA_BACKEND" in
  slurm) ;;
  local)
    : "${CHOLERA_CORES:=$(nproc 2>/dev/null || echo 1)}"
    export CHOLERA_CORES
    case ":$PATH:" in
      *":$CHOLERA_BASH_DIR/lib/local:"*) ;;
      *) export PATH="$CHOLERA_BASH_DIR/lib/local:$PATH" ;;
    esac ;;
  *) echo "Error: CHOLERA_BACKEND must be slurm or local (got '$CHOLERA_BACKEND')." >&2
     return 1 2>/dev/null || exit 1 ;;
esac

cholera_env() {
  if [[ "${1:-}" == "--purge" ]]; then
    shift
    if [[ -n "${CHOLERA_ENV_PURGE:-}" ]]; then eval "$CHOLERA_ENV_PURGE"; fi
  fi
  local stack="${1:-}" var setup
  var="CHOLERA_ENV_SETUP${stack:+_${stack^^}}"
  setup="${!var-${CHOLERA_ENV_SETUP:-}}"
  [[ -z "$setup" ]] && return 0
  eval "$setup"
}

# sbatch options that take a value as the next word when written without '='.
_cholera_opt_takes_value() {
  case "$1" in
    -[AaBbCcDdeiJLMmNnopqSstwx]) return 0 ;;
    --parsable|--hold|--exclusive|--requeue|--no-requeue|--wait|--test-only|--quiet|--verbose|--use-min-nodes|--overcommit|--contiguous|--spread-job|--oversubscribe|--ignore-pbs|--no-kill|--help|--usage|--version) return 1 ;;
    --*=*) return 1 ;;
    --*) return 0 ;;
    *) return 1 ;;
  esac
}

cholera_sbatch_cmd() {
  local class="${1:-}"; shift || true
  case "$class" in
    SHORT|LONG|BIGMEM) ;;
    *) echo "Error: cholera_sbatch: class must be SHORT, LONG or BIGMEM (got '$class')." >&2; return 1 ;;
  esac
  local key="" defaults=""
  while [[ "${1:-}" == "--res" || "${1:-}" == "--defaults" ]]; do
    if [[ "$1" == "--res" ]]; then key="$2"; else defaults="$2"; fi
    shift 2
  done
  local -a args=("$@")
  # Which options did the caller pass, and which word is the payload?
  local have=" " i=0 a payload="" array=false
  while (( i < ${#args[@]} )); do
    a="${args[$i]}"
    if [[ "$a" != -* ]]; then payload="$a"; break; fi
    case "$a" in
      --partition*|-p)        have+="partition " ;;
      --account*|-A)          have+="account " ;;
      --qos*|-q)              have+="qos " ;;
      --mem=*|--mem)          have+="mem " ;;
      --time*|-t)             have+="time " ;;
      --cpus-per-task*|-c)    have+="cpus " ;;
      --output*|-o)           have+="output " ;;
      --array*|-a)            array=true ;;
    esac
    if [[ "$a" == -?* && "$a" != *=* ]] && _cholera_opt_takes_value "$a"; then i=$((i + 2)); else i=$((i + 1)); fi
  done
  if [[ -z "$key" ]]; then
    if [[ -n "$payload" ]]; then
      key=$(basename "$payload" .sh); key=${key#submit_}
    else
      key=WRAP
    fi
  fi
  key=${key^^}

  local -a extra=()
  local res_var="CHOLERA_RES_${key}"
  local -a res=()
  # shellcheck disable=SC2206
  [[ -n "${!res_var:-}" ]] && res=(${!res_var})
  # shellcheck disable=SC2206
  [[ -n "$defaults" ]] && res+=($defaults)
  local r
  for r in ${res[@]+"${res[@]}"}; do
    case "$r" in
      --partition=*|-p*)      [[ "$have" == *" partition "* ]] && continue; have+="partition " ;;
      --account=*)            [[ "$have" == *" account "* ]] && continue; have+="account " ;;
      --qos=*)                [[ "$have" == *" qos "* ]] && continue; have+="qos " ;;
      --mem=*)                [[ "$have" == *" mem "* ]] && continue; have+="mem " ;;
      --time=*)               [[ "$have" == *" time "* ]] && continue; have+="time " ;;
      --cpus-per-task=*|-c*)  [[ "$have" == *" cpus "* ]] && continue; have+="cpus " ;;
      --output=*)             [[ "$have" == *" output "* ]] && continue; have+="output " ;;
    esac
    extra+=("$r")
  done

  local v
  v="CHOLERA_PARTITION_${class}"; [[ "$have" != *" partition "* && -n "${!v:-}" ]] && extra+=("--partition=${!v}")
  v="CHOLERA_ACCOUNT_${class}"; local acct="${!v:-${CHOLERA_ACCOUNT:-}}"
  [[ "$have" != *" account "* && -n "$acct" ]] && extra+=("--account=$acct")
  v="CHOLERA_QOS_${class}"; local qos="${!v:-${CHOLERA_QOS:-}}"
  [[ "$have" != *" qos "* && -n "$qos" ]] && extra+=("--qos=$qos")
  v="CHOLERA_MEM_${class}";  [[ "$have" != *" mem "*  && -n "${!v:-}" ]] && extra+=("--mem=${!v}")
  v="CHOLERA_TIME_${class}"; [[ "$have" != *" time "* && -n "${!v:-}" ]] && extra+=("--time=${!v}")
  v="CHOLERA_CPUS_${class}"; [[ "$have" != *" cpus "* && -n "${!v:-}" ]] && extra+=("--cpus-per-task=${!v}")
  if [[ "$have" != *" output "* ]]; then
    if [[ "$array" == true ]]; then extra+=("--output=${CHOLERA_LOG_DIR}/%x_%A_%a.log")
    else extra+=("--output=${CHOLERA_LOG_DIR}/%x_%j.log"); fi
  fi
  CHOLERA_SBATCH_CMD=(sbatch ${extra[@]+"${extra[@]}"} ${args[@]+"${args[@]}"})
}

cholera_sbatch() {
  cholera_sbatch_cmd "$@" || return 1
  mkdir -p "$CHOLERA_LOG_DIR"
  "${CHOLERA_SBATCH_CMD[@]}"
}

fi
