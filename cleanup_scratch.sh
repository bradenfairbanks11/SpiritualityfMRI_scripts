#!/bin/bash
#SBATCH --job-name=cleanup_scratch
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
# =============================================================================
# Clear the pre-reorg leftovers sitting loose at the root of SCRATCH.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo.
#
# CONTEXT: before the 2026-08-19 reorg, jobs wrote working directories straight to
# the scratch root — ~335 GB of them. Everything now lives under
# $SCRATCH_ROOT/{Spirituality,Missionary_language}/ and $SCRATCH_ROOT/work/.
# This removes the strays so the scratch root shows only the current layout.
#
# THIS IS MUCH LOWER STAKES THAN prune_archive.sh:
#   These are fMRIPrep/nipype WORKING directories — scratch intermediates, not
#   results. Their only value is resuming a crashed run. Every derivative they
#   produced was already copied and verified.
#
# THE ONE REAL RISK, AND THE GUARD AGAINST IT:
#   Deleting the work dir of a run that never finished means that run restarts from
#   zero. So a subject's work dir is only removed once all three completion markers
#   exist in the derivatives: the fMRIPrep .html report, tedana output, and
#   afni_firstlvl output. Anything short of that is kept and reported.
#
# WHAT IT PRESERVES BEFORE DELETING:
#   Loose run logs and run_test.sh from the superseded group dirs are copied to
#   ~/spirituality_fmri/logs/legacy_group/ (HOME — backed up) first. They are
#   provenance for how the old group analysis was invoked, and they are kilobytes.
#
# WHAT IT DELIBERATELY DOES NOT TOUCH:
#   ms_dataset/  — 48 GB, the parked MS project. Its derivatives exist nowhere else,
#                  and scratch purges after 12 weeks UNUSED. Deleting or keeping that
#                  is a real decision about a real project, not tidying. See SUMMARY.
#   work/        — the current WORK_ROOT.
#   Spirituality/, Missionary_language/ — the new layout.
#
# USAGE:
#   bash cleanup_scratch.sh                            # dry run, safe
#   sbatch --export=ALL,DRYRUN=0 cleanup_scratch.sh    # for real
# =============================================================================
set -uo pipefail

DRYRUN=${DRYRUN:-1}

HOME_ROOT=/home/bradenf4
SCRATCH_ROOT=/nobackup/autodelete/usr/bradenf4
SPIR_SCRATCH=${SCRATCH_ROOT}/Spirituality
DERIV=${SPIR_SCRATCH}/derivatives
LEGACY=${HOME_ROOT}/spirituality_fmri/logs/legacy_group

REMOVED=0
KEPT=0

banner() { printf '\n========== %s ==========\n' "$*"; }

echo "cleanup_scratch.sh   DRYRUN=${DRYRUN}   host=$(hostname)   started $(date)"
echo "job: ${SLURM_JOB_ID:-<none, running outside Slurm>}"
if [ "${DRYRUN}" != "1" ] && [ -z "${SLURM_JOB_ID:-}" ]; then
    echo
    echo "REFUSING TO RUN: deleting hundreds of GB is heavy metadata I/O and belongs"
    echo "in Slurm, not on a login node (ORC rules 2 and 8)."
    echo "Use:  sbatch --export=ALL,DRYRUN=0 cleanup_scratch.sh"
    exit 1
fi

cd "${SCRATCH_ROOT}" || { echo "cannot cd to ${SCRATCH_ROOT}"; exit 1; }

# --- 1. Preserve the small provenance bits from the superseded group dirs ------
banner "PRESERVE: legacy group run logs -> HOME"
if [ "${DRYRUN}" = "1" ]; then
    echo "[dry run] would create ${LEGACY} and copy:"
    ls spirituality_group/run_full*.log spirituality_group_test/run_test.sh 2>/dev/null
else
    mkdir -p "${LEGACY}"
    for f in spirituality_group/run_full*.log spirituality_group_test/run_test.sh; do
        [ -f "${f}" ] || continue
        cp -p "${f}" "${LEGACY}/$(echo "${f}" | tr '/' '_')" && echo "kept: ${f}"
    done
    echo "preserved into ${LEGACY}"
fi

# --- 2. Per-subject working directories ---------------------------------------
# Only removed once the subject's results are all present. A work dir whose run
# never finished is worth more than the disk it occupies.
banner "WORK DIRS: per-subject fMRIPrep/tedana intermediates"
for d in sub-*; do
    [ -d "${d}" ] || continue
    # Strip any _ses-N suffix to get the subject label the derivatives use.
    subj=${d%%_*}
    size=$(du -sh "${d}" 2>/dev/null | cut -f1)

    have_html=0; have_tedana=0; have_afni=0
    [ -f "${DERIV}/fmriprep/${subj}.html" ]  && have_html=1
    [ -d "${DERIV}/tedana/${subj}" ]         && have_tedana=1
    [ -d "${DERIV}/afni_firstlvl/${subj}" ]  && have_afni=1

    echo
    echo "${d}  (${size})  -> subject ${subj}"
    echo "   fmriprep report: ${have_html}   tedana: ${have_tedana}   afni_firstlvl: ${have_afni}"

    if [ "${have_html}" = "1" ] && [ "${have_tedana}" = "1" ] && [ "${have_afni}" = "1" ]; then
        echo "   COMPLETE — work dir is a disposable intermediate."
        if [ "${DRYRUN}" = "1" ]; then
            echo "   [dry run] would delete ${d}"
        else
            rm -rf "${d}" && echo "   deleted."
        fi
        REMOVED=$((REMOVED + 1))
    elif [ "$(find "${d}" -type f 2>/dev/null | wc -l)" -lt 20 ]; then
        # sub-02 and friends: a handful of bids_db/config.toml files from runs that
        # were abandoned before producing anything. Nothing to resume, nothing to lose.
        echo "   ABANDONED STUB — only $(find "${d}" -type f 2>/dev/null | wc -l) files, no derivatives. Removing."
        if [ "${DRYRUN}" = "1" ]; then
            echo "   [dry run] would delete ${d}"
        else
            rm -rf "${d}" && echo "   deleted."
        fi
        REMOVED=$((REMOVED + 1))
    else
        echo "   *** KEEPING — results are incomplete. This run may still be resumable."
        KEPT=$((KEPT + 1))
    fi
done

# --- 3. Superseded group output directories -----------------------------------
# Every result file in these is byte-identical to derivatives/afni_group. The only
# non-duplicate content is staged inputs/ (re-stageable from afni_firstlvl by
# second_level_afni.sh) and the run logs preserved in step 1.
banner "SUPERSEDED: old group output dirs"
for d in spirituality_group spirituality_group_test; do
    [ -d "${d}" ] || { echo "SKIP — already gone: ${d}"; continue; }
    echo "${d}  ($(du -sh "${d}" 2>/dev/null | cut -f1))  superseded by derivatives/afni_group"
    if [ "${DRYRUN}" = "1" ]; then
        echo "   [dry run] would delete ${d}"
    else
        rm -rf "${d}" && echo "   deleted."
    fi
    REMOVED=$((REMOVED + 1))
done

# --- 4. Empty and near-empty root-level run dirs ------------------------------
banner "STRAYS: root-level timestamped run dirs"
n_stray=0
for d in 2025[0-9]*-[0-9]*_* 2026[0-9]*-[0-9]*_* fmriprep_25_1_wf; do
    [ -d "${d}" ] || continue
    n_stray=$((n_stray + 1))
    if [ "${DRYRUN}" != "1" ]; then rm -rf "${d}"; fi
done
if [ "${DRYRUN}" = "1" ]; then
    echo "[dry run] would delete ${n_stray} stray run dir(s) (empty or bids_db/config.toml only)"
else
    echo "deleted ${n_stray} stray run dir(s)"
fi
REMOVED=$((REMOVED + n_stray))

# --- Summary -------------------------------------------------------------------
banner "SUMMARY"
echo "finished $(date)"
if [ "${DRYRUN}" = "1" ]; then
    echo "DRY RUN — nothing deleted. Would remove ${REMOVED} item(s); ${KEPT} kept as incomplete."
else
    echo "Removed ${REMOVED} item(s). ${KEPT} kept because their results are incomplete."
fi
echo
if [ "${DRYRUN}" = "1" ]; then
    echo "Scratch root as it stands RIGHT NOW (unchanged — this was a dry run):"
else
    echo "Scratch root now:"
fi
ls -d "${SCRATCH_ROOT}"/* 2>/dev/null | sed 's|^|  |'
echo
echo "STILL NEEDS A DECISION FROM YOU: ms_dataset/ (48 GB, the parked MS project)."
echo "  Its containers are safe in ~/software/containers, but its derivatives exist"
echo "  ONLY there — and scratch purges after 12 weeks unused. Either promote what"
echo "  matters to archive, or accept losing it. Untouched by this script either way."
exit 0
