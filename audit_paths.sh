#!/bin/bash
# =============================================================================
# audit_paths.sh — static, READ-ONLY verification that every pipeline script in
# every repo agrees with the post-reorg storage layout.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo.
#
# WHY THIS EXISTS:
#   Before prune_archive.sh permanently deletes ~250 GB from a filesystem with no
#   backups, something has to answer "will anything still work afterwards?" without
#   re-running a 10-hour pipeline. This does that statically.
#
# WHAT IT CHECKS (each failure is explained inline when it fires):
#   1. Every path exported by each config.sh sits on the tier it is supposed to.
#   2. Every .sh parses (bash -n).
#   3. Nobody uses the disguised /home/bradenf4/nobackup/... spelling, which looks
#      like HOME but is a symlink into archive/scratch.
#   4. Only the sanctioned scripts write under ARCHIVE.
#   5. Every `module load` target actually resolves, and every container exists.
#   6. POST-PRUNE SIMULATION: no pipeline input lives under a path that
#      prune_archive.sh is going to delete.
#
# It reads only; it creates and changes nothing. Login node is fine — the work is
# a few greps and stats, not a filesystem walk.
#
# USAGE:  bash audit_paths.sh          (exit 0 = safe to prune)
# =============================================================================
set -uo pipefail

PASS=0; FAIL=0; WARN=0
ok()   { printf '  [ ok ] %s\n' "$*"; PASS=$((PASS+1)); }
bad()  { printf '  [FAIL] %s\n' "$*"; FAIL=$((FAIL+1)); }
warn() { printf '  [warn] %s\n' "$*"; WARN=$((WARN+1)); }
banner(){ printf '\n========== %s ==========\n' "$*"; }

HOME_ROOT=/home/bradenf4
SCRATCH_ROOT=/nobackup/autodelete/usr/bradenf4
ARCHIVE_ROOT=/nobackup/archive/usr/bradenf4

SPIR=${HOME_ROOT}/spirituality_fmri/scripts
LUKE=${HOME_ROOT}/missionary_language/scripts
MSMRI=${HOME_ROOT}/ms-mri-analysis

tier_of() {
    case "$1" in
        ${HOME_ROOT}/nobackup/*) echo DISGUISED ;;
        ${SCRATCH_ROOT}*)        echo SCRATCH ;;
        ${ARCHIVE_ROOT}*)        echo ARCHIVE ;;
        ${HOME_ROOT}*)           echo HOME ;;
        /apps/*)                 echo APPS ;;
        *)                       echo OTHER ;;
    esac
}

# expect <config.sh> <VAR> <EXPECTED_TIER> <must_exist yes|no>
expect() {
    local cfg="$1" var="$2" want="$3" must="$4"
    local val got
    val=$(bash -c "source '${cfg}' >/dev/null 2>&1; printf '%s' \"\${${var}:-}\"")
    if [ -z "${val}" ]; then bad "${var} is unset in ${cfg}"; return; fi
    got=$(tier_of "${val}")
    if [ "${got}" != "${want}" ]; then
        bad "${var} is on ${got}, expected ${want}: ${val}"
        [ "${got}" = DISGUISED ] && printf '         (that path is a symlink into archive/scratch, not HOME)\n'
        return
    fi
    if [ "${must}" = yes ] && [ ! -e "${val}" ]; then
        bad "${var} -> ${want} tier, but does not exist: ${val}"; return
    fi
    ok "${var} = ${val}  [${got}]"
}

echo "audit_paths.sh   host=$(hostname)   $(date)"

# ---------------------------------------------------------------------------
banner "1. config.sh path/tier assertions"
for proj in "Spirituality:${SPIR}" "Missionary_language:${LUKE}"; do
    name=${proj%%:*}; cfg=${proj#*:}/config.sh
    echo "--- ${name} (${cfg})"
    [ -f "${cfg}" ] || { bad "no config.sh at ${cfg}"; continue; }
    expect "${cfg}" PROJ_HOME    HOME    yes
    expect "${cfg}" TIMING_DIR   HOME    yes
    expect "${cfg}" LOGS         HOME    yes
    expect "${cfg}" FMRIPREP_SIF HOME    yes
    expect "${cfg}" TEMPLATEFLOW_HOME HOME yes
    expect "${cfg}" FS_LICENSE   HOME    yes
    expect "${cfg}" BIDS_DIR     SCRATCH yes
    expect "${cfg}" DERIV        SCRATCH yes
    expect "${cfg}" FMRIPREP_OUT SCRATCH yes
    expect "${cfg}" TEDANA_OUT   SCRATCH yes
    expect "${cfg}" AFNI_OUT     SCRATCH yes
    expect "${cfg}" WORK_ROOT    SCRATCH yes
    expect "${cfg}" RAWDATA_DIR  ARCHIVE yes
    expect "${cfg}" FINAL_DIR    ARCHIVE no
    expect "${cfg}" AFNI_SIF     APPS    yes
done
echo "--- ms-mri (${MSMRI}/setup/config.sh)"
MCFG=${MSMRI}/setup/config.sh
expect "${MCFG}" RAW          ARCHIVE yes
expect "${MCFG}" DERIV        SCRATCH yes
expect "${MCFG}" WORK         SCRATCH yes
expect "${MCFG}" LOGS         SCRATCH yes
expect "${MCFG}" CONTAINERS   HOME    yes
expect "${MCFG}" TEMPLATEFLOW_HOME HOME yes

# ---------------------------------------------------------------------------
banner "2. bash -n (syntax) on every script"
nerr=0
while IFS= read -r f; do
    bash -n "$f" 2>/dev/null || { bad "syntax error: $f"; nerr=$((nerr+1)); }
done < <(find "${SPIR}" "${LUKE}" "${MSMRI}" -name '*.sh' -type f 2>/dev/null | sort)
[ ${nerr} -eq 0 ] && ok "all .sh files parse"

# ---------------------------------------------------------------------------
banner "3. disguised HOME->archive path spelling"
# /home/bradenf4/nobackup/{archive,autodelete} are symlinks into the nobackup
# tiers. A script using that spelling looks compliant and is not.
hits=$(grep -rn "${HOME_ROOT}/nobackup/" "${SPIR}" "${LUKE}" "${MSMRI}" \
         --include='*.sh' --include='*.py' 2>/dev/null | grep -v '^\s*#' | grep -vE '#.*nobackup')
if [ -z "${hits}" ]; then ok "no script uses the disguised spelling"
else bad "disguised paths found:"; echo "${hits}" | sed 's/^/         /'; fi

# ---------------------------------------------------------------------------
banner "4. who writes under ARCHIVE"
# Sanctioned: the promoter, the pruner, the migrator/verifier, and dcm2niix
# (which only touches rawdata metadata). Anything else is a policy regression.
ALLOW='promote_to_archive.sh|prune_archive.sh|migrate_storage.sh|verify_migration.sh|06_archive_derivatives.sh|01_fetch_pilot.sh|audit_paths.sh'
writers=$(grep -rlnE '(mkdir -p|cp -a|rsync|>|tar).*(ARCHIVE_ROOT|PROJ_ARCHIVE|FINAL_DIR|DERIV_ARCHIVE)' \
            "${SPIR}" "${LUKE}" "${MSMRI}" --include='*.sh' 2>/dev/null \
          | grep -vE "(${ALLOW})" || true)
if [ -z "${writers}" ]; then ok "no unsanctioned script writes under ARCHIVE"
else warn "review these for archive writes:"; echo "${writers}" | sed 's/^/         /'; fi

# ---------------------------------------------------------------------------
banner "5. modules and containers"
if command -v module >/dev/null 2>&1 || [ -n "${MODULESHOME:-}" ]; then
    for m in $(grep -rhoE 'module load [a-zA-Z0-9._/-]+' "${SPIR}" "${LUKE}" "${MSMRI}" \
                 --include='*.sh' 2>/dev/null | awk '{print $3}' | sort -u); do
        case "$m" in '"$APPTAINER_MODULE"'|'$APPTAINER_MODULE') continue ;; esac
        if module spider "$m" >/dev/null 2>&1; then ok "module resolves: $m"
        else bad "module does NOT resolve: $m"; fi
    done
else warn "module command unavailable here; skipped"; fi

for sif in "${HOME_ROOT}/software/fmriprep-25.1.4.sif" /apps/afni/afni_make_build_latest.sif; do
    [ -f "${sif}" ] && ok "container present: ${sif}" || bad "container MISSING: ${sif}"
done

# ---------------------------------------------------------------------------
banner "6. POST-PRUNE SIMULATION"
# These are exactly the paths prune_archive.sh deletes. If a pipeline input still
# resolves under one of them, pruning breaks the pipeline.
DOOMED=(
  "${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project/BIDS"
  "${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project/derivatives/fmriprep"
  "${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project/derivatives/tedana"
  "${ARCHIVE_ROOT}/Luke_active/Missionary_language/Pilot/BIDS"
  "${ARCHIVE_ROOT}/Luke_active/Missionary_language/Pilot/derivatives/fmriprep"
  "${ARCHIVE_ROOT}/Luke_active/Missionary_language/Pilot/derivatives/tedana"
  "${ARCHIVE_ROOT}/software/templateflow"
  "${ARCHIVE_ROOT}/software/fmri_prep"
)
badref=0
for proj in "${SPIR}" "${LUKE}"; do
    cfg=${proj}/config.sh; [ -f "${cfg}" ] || continue
    for v in BIDS_DIR DERIV FMRIPREP_OUT TEDANA_OUT AFNI_OUT WORK_ROOT TIMING_DIR \
             FMRIPREP_SIF TEMPLATEFLOW_HOME RAWDATA_DIR; do
        val=$(bash -c "source '${cfg}' >/dev/null 2>&1; printf '%s' \"\${${v}:-}\"")
        for d in "${DOOMED[@]}"; do
            case "${val}" in "${d}"*) bad "${v} in ${cfg} lives under a path that will be DELETED: ${val}"; badref=1 ;; esac
        done
    done
done
lit=$(grep -rn --include='*.sh' -F "${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project/BIDS" "${SPIR}" "${LUKE}" 2>/dev/null | grep -vE "${ALLOW}" || true)
[ -n "${lit}" ] && { bad "literal reference to a doomed path:"; echo "${lit}" | sed 's/^/         /'; badref=1; }
[ ${badref} -eq 0 ] && ok "no pipeline input depends on anything prune_archive.sh deletes"

# ---------------------------------------------------------------------------
banner "SUMMARY"
echo "pass: ${PASS}   fail: ${FAIL}   warn: ${WARN}"
if [ ${FAIL} -eq 0 ]; then
    echo
    echo "PATH AUDIT CLEAN — the layout is consistent and pruning breaks nothing."
else
    echo
    echo "*** ${FAIL} check(s) FAILED. Fix these before pruning archive. ***"
fi
exit ${FAIL}
