#!/bin/bash
# =============================================================================
# Promote final products from SCRATCH to ARCHIVE.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo
# (upstream: /apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md).
#
# WHY THIS EXISTS:
#   All compute now happens on SCRATCH (/nobackup/autodelete), which deletes files
#   UNUSED for 12 weeks. The trigger is disuse, not age — so a project you set down in
#   September loses its derivatives in December, including runs that cost 10-12 h of
#   fMRIPrep per subject. Without a deliberate promote step, that purge is the de facto
#   retention policy for your results.
#
#   *** Run this before parking a project, and after any stage you consider final. ***
#
# WHAT COUNTS AS FINAL (decided 2026-08-19):
#   - group / second-level maps          (afni_group)
#   - test-retest reproducibility output (afni_reproducibility)
#   - per-run first-level _stats_REML buckets + small provenance
#
#   NOT promoted, because they are regenerable and enormous:
#     _errts / _errts_REML  (~1.0 GB per run)   residuals, recomputed by 3dREMLfit
#     _blurred, _blurred_scaled (~570 MB/run)   intermediates, recomputed in STEP 1-2
#     _stats+tlrc (non-REML)                    superseded by the REML bucket
#     fmriprep / tedana                         regenerable from BIDS (10-12 h/subject)
#     BIDS                                      regenerable from the DICOMs in archive
#
#   Per task dir this is ~15 MB kept out of ~1.7 GB, so 24 runs promote in well under
#   1 GB instead of 47 GB.
#
# WHY TAR:
#   Archive has a 1 MILLION FILE limit (146k used as of 2026-08-19) and is migrating to
#   a much slower backing store. One sequential tarball is far kinder to both than a
#   loose tree of thousands of small files.
#
# USAGE:
#   DRYRUN=1 bash promote_to_archive.sh    # list what would be promoted, write nothing
#   bash promote_to_archive.sh             # do it
#
# Output is small (<1 GB), so this is fine on a login node. If you ever promote
# something large, sbatch it instead — ORC asks that sustained I/O go through Slurm.
# =============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

DRYRUN=${DRYRUN:-0}
STAMP=$(date +%Y%m%d)

# Files worth keeping from each first-level task dir. Everything else is regenerable.
KEEP_PATTERNS=(
    '*_stats_REML+tlrc.BRIK'      # THE result: betas, t-stats, rating_slope
    '*_stats_REML+tlrc.HEAD'
    '*_xmat.1D'                   # design matrix — provenance, tiny
    '*_xmat.jpg'
    '*_stats.REML_cmd'            # exact command used
    '*_motion_censor.1D'          # which TRs were censored
    '*_motion_enorm.1D'
    '*_motion_demean.txt'         # the Friston-24 regressors actually fed to -ortvec
)

echo "promote_to_archive.sh   DRYRUN=${DRYRUN}   $(date)"
echo "  from SCRATCH: ${DERIV}"
echo "  to   ARCHIVE: ${FINAL_DIR}"
echo

if [ ! -d "${DERIV}" ]; then
    echo "ERROR: ${DERIV} does not exist. Has migrate_storage.sh run yet?"
    exit 1
fi

if [ "${DRYRUN}" != "1" ]; then
    mkdir -p "${FINAL_DIR}"
fi

FAILED=0

# --- 1. First-level buckets (filtered) ---------------------------------------
echo "========== first-level _stats_REML buckets =========="
FILELIST=$(mktemp)
FIND_ARGS=()
for p in "${KEEP_PATTERNS[@]}"; do
    FIND_ARGS+=(-name "$p" -o)
done
unset 'FIND_ARGS[${#FIND_ARGS[@]}-1]'          # drop the trailing -o

if [ -d "${AFNI_OUT}" ]; then
    ( cd "${AFNI_OUT}" && find . -type f \( "${FIND_ARGS[@]}" \) -print ) | sort > "${FILELIST}"
    N=$(wc -l < "${FILELIST}")
    echo "matched ${N} files"
    if [ "${N}" -eq 0 ]; then
        echo "WARNING: nothing matched — check that first-levels have actually run."
    elif [ "${DRYRUN}" = "1" ]; then
        echo "[dry run] would tar these into ${FINAL_DIR}/afni_firstlvl_stats_${STAMP}.tar.gz"
        head -20 "${FILELIST}"; [ "${N}" -gt 20 ] && echo "  ... and $((N - 20)) more"
    else
        TARBALL=${FINAL_DIR}/afni_firstlvl_stats_${STAMP}.tar.gz
        if tar -czf "${TARBALL}" -C "${AFNI_OUT}" -T "${FILELIST}"; then
            echo "wrote ${TARBALL} ($(du -h "${TARBALL}" | cut -f1))"
            # Verify the archive is readable and has the expected member count.
            M=$(tar -tzf "${TARBALL}" | grep -c . )
            echo "verify: ${M} members readable (expected ${N})"
            [ "${M}" -eq "${N}" ] || { echo "*** MEMBER COUNT MISMATCH ***"; FAILED=$((FAILED+1)); }
        else
            echo "FAILED to write ${TARBALL}"; FAILED=$((FAILED+1))
        fi
    fi
else
    echo "SKIP — ${AFNI_OUT} not present"
fi
rm -f "${FILELIST}"
echo

# --- 2. Whole-tree components (already small) --------------------------------
promote_tree() {
    local label="$1" src="$2" name="$3"
    echo "========== ${label} =========="
    if [ ! -d "${src}" ]; then
        echo "SKIP — ${src} not present"; echo; return 0
    fi
    echo "src: ${src} ($(du -sh "${src}" | cut -f1))"
    if [ "${DRYRUN}" = "1" ]; then
        echo "[dry run] would tar into ${FINAL_DIR}/${name}_${STAMP}.tar.gz"
    else
        local tb=${FINAL_DIR}/${name}_${STAMP}.tar.gz
        if tar -czf "${tb}" -C "$(dirname "${src}")" "$(basename "${src}")"; then
            echo "wrote ${tb} ($(du -h "${tb}" | cut -f1))"
            tar -tzf "${tb}" >/dev/null && echo "verify: archive reads cleanly" \
                || { echo "*** TARBALL UNREADABLE ***"; FAILED=$((FAILED+1)); }
        else
            echo "FAILED to write ${tb}"; FAILED=$((FAILED+1))
        fi
    fi
    echo
}

promote_tree "group / second-level maps" "${GROUP_OUT}" "afni_group"
promote_tree "test-retest reproducibility" "${REPRO_OUT}" "afni_reproducibility"

# --- 3. Timing files ----------------------------------------------------------
# These already live on /home (backed up daily, never purged), so archive is a third
# copy rather than a rescue. Cheap insurance: they cannot be regenerated on-cluster.
echo "========== timing files (extra copy) =========="
if [ -d "${TIMING_DIR}" ] && [ "${DRYRUN}" != "1" ]; then
    tar -czf "${FINAL_DIR}/timing_${STAMP}.tar.gz" -C "$(dirname "${TIMING_DIR}")" "$(basename "${TIMING_DIR}")" \
        && echo "wrote timing_${STAMP}.tar.gz ($(ls -1 "${TIMING_DIR}" | wc -l) files)"
elif [ "${DRYRUN}" = "1" ]; then
    echo "[dry run] would tar ${TIMING_DIR} ($(ls -1 "${TIMING_DIR}" 2>/dev/null | wc -l) files)"
fi
echo

# --- Summary ------------------------------------------------------------------
echo "========== SUMMARY =========="
if [ "${DRYRUN}" = "1" ]; then
    echo "Dry run — nothing written."
else
    echo "Promoted to ${FINAL_DIR}:"
    ls -lh "${FINAL_DIR}" 2>/dev/null | tail -n +2
    echo
    echo "REMINDER: archive is NOT backed up. Keep the offsite copy current too."
fi
[ "${FAILED}" -eq 0 ] || echo "*** ${FAILED} failure(s) — see above ***"
exit "${FAILED}"
