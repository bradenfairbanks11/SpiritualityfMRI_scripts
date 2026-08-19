#!/bin/bash
#SBATCH --job-name=migrate_storage
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
# =============================================================================
# One-time storage migration: get both fMRI pipelines off /nobackup/archive.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo
# (upstream: /apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md).
#
# WHY THIS IS AN SBATCH JOB, NOT A LOGIN-NODE COMMAND:
#   ORC policy — "Use Slurm for sustained, computationally intensive, or high-I/O work.
#   Do not circumvent login-node resource limits." This copies hundreds of GB off a
#   filesystem that is mid-migration to a much slower backing store.
#
# WHAT IT DOES:
#   1. Copies the shared containers + TemplateFlow from ARCHIVE -> HOME  (~5.6 GB)
#   2. Copies BIDS + derivatives from ARCHIVE -> SCRATCH for both projects
#   3. Reports du -sh per component so future runs can be sized properly
#
# WHAT IT DOES NOT DO:
#   *** IT DELETES NOTHING. *** Every source tree is left exactly as it was.
#   Verification and deletion are a separate, deliberate, manual step. Archive is
#   NOT backed up — there is no undo there.
#
# USAGE:
#   DRYRUN=1 bash migrate_storage.sh     # login node, safe: sizes everything, copies nothing
#   sbatch migrate_storage.sh            # the real run
#
# Run the dry run FIRST. It prints the transfer sizes, which tell you whether the
# 12 h wall time above is generous or tight.
# =============================================================================
set -uo pipefail

DRYRUN=${DRYRUN:-0}

HOME_ROOT=/home/bradenf4
SCRATCH_ROOT=/nobackup/autodelete/usr/bradenf4
ARCHIVE_ROOT=/nobackup/archive/usr/bradenf4

SPIR_ARCHIVE=${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project
SPIR_SCRATCH=${SCRATCH_ROOT}/Spirituality
LUKE_ARCHIVE=${ARCHIVE_ROOT}/Luke_active/Missionary_language/Pilot
LUKE_SCRATCH=${SCRATCH_ROOT}/Missionary_language

RSYNC_OPTS=(-a --stats --human-readable)
[ "${DRYRUN}" = "1" ] && RSYNC_OPTS+=(--dry-run)

# The timing/ dir lives INSIDE the first-level output tree. Its authoritative copy is
# now ~/<project>/timing (backed up, never purged) because the generator runs on
# Braden's local machine and cannot be re-run on the cluster. Excluding it here keeps
# exactly one authoritative copy instead of two that can silently diverge.
RSYNC_OPTS+=(--exclude 'timing/')

FAILED=0

banner() { printf '\n========== %s ==========\n' "$*"; }

copy_component() {
    local label="$1" src="$2" dst="$3"
    banner "${label}"
    if [ ! -e "${src}" ]; then
        echo "SKIP — source does not exist: ${src}"
        return 0
    fi
    echo "src: ${src}"
    echo "dst: ${dst}"
    echo "size: $(du -sh "${src}" 2>/dev/null | cut -f1)"
    if [ "${DRYRUN}" = "1" ]; then
        echo "[dry run] would create ${dst}"
    else
        mkdir -p "${dst}"
    fi
    if rsync "${RSYNC_OPTS[@]}" "${src}/" "${dst}/"; then
        echo "OK: ${label}"
    else
        echo "FAILED: ${label} (rsync exit $?)"
        FAILED=$((FAILED + 1))
    fi
}

echo "migrate_storage.sh   DRYRUN=${DRYRUN}   host=$(hostname)   started $(date)"
echo "job: ${SLURM_JOB_ID:-<none, running outside Slurm>}"
if [ "${DRYRUN}" != "1" ] && [ -z "${SLURM_JOB_ID:-}" ]; then
    echo
    echo "REFUSING TO RUN: a real migration must go through Slurm, not a login node."
    echo "Use:  sbatch migrate_storage.sh     (or DRYRUN=1 bash migrate_storage.sh to size it)"
    exit 1
fi

# --- 1. Shared software: ARCHIVE -> HOME -------------------------------------
# Every job start was pulling a 2.3 GB container plus TemplateFlow off archive.
copy_component "SOFTWARE: TemplateFlow"  "${ARCHIVE_ROOT}/software/templateflow" \
                                         "${HOME_ROOT}/software/templateflow"

banner "SOFTWARE: fMRIPrep container"
SIF_SRC=${ARCHIVE_ROOT}/software/fmri_prep/my_images/fmriprep-25.1.4.sif
SIF_DST=${HOME_ROOT}/software/fmriprep-25.1.4.sif
echo "src: ${SIF_SRC}  ($(du -sh "${SIF_SRC}" 2>/dev/null | cut -f1))"
echo "dst: ${SIF_DST}"
if [ "${DRYRUN}" = "1" ]; then
    echo "[dry run] would copy the .sif"
else
    mkdir -p "${HOME_ROOT}/software"
    if rsync -a --stats --human-readable "${SIF_SRC}" "${SIF_DST}"; then
        echo "OK: fMRIPrep container"
    else
        echo "FAILED: fMRIPrep container"; FAILED=$((FAILED + 1))
    fi
fi
# NOTE: the FreeSurfer license is NOT copied — ~/.freesurfer_license.txt already exists
# and is byte-identical to the archive copy (verified 2026-08-19). config.sh points there.

# --- 1b. ms-mri containers: SCRATCH -> HOME (rescue from the purge clock) -----
# These sat under /nobackup/autodelete/.../ms_dataset/containers, i.e. the 12-week
# UNUSED purge tier, last touched 2026-06-30. The ms-mri project is parked, so they
# were on track to be deleted around 2026-09-22 — and compute nodes have no internet,
# so a purged container cannot be re-pulled from inside a job.
# ms-mri-analysis/setup/config.sh now points $CONTAINERS at ~/software/containers.
copy_component "SOFTWARE: ms-mri containers (20 GB)" \
               "${SCRATCH_ROOT}/ms_dataset/containers" \
               "${HOME_ROOT}/software/containers"

# ms-mri kept its own TemplateFlow cache. Merged into the shared one (rsync without
# --delete takes the union, so templates unique to either side survive).
copy_component "SOFTWARE: ms-mri TemplateFlow (merge into shared)" \
               "${SCRATCH_ROOT}/ms_dataset/templateflow" \
               "${HOME_ROOT}/software/templateflow"

# --- 2. Spirituality: ARCHIVE -> SCRATCH -------------------------------------
copy_component "SPIRITUALITY: BIDS"        "${SPIR_ARCHIVE}/BIDS"                        "${SPIR_SCRATCH}/BIDS"
copy_component "SPIRITUALITY: fmriprep"    "${SPIR_ARCHIVE}/derivatives/fmriprep"        "${SPIR_SCRATCH}/derivatives/fmriprep"
copy_component "SPIRITUALITY: tedana"      "${SPIR_ARCHIVE}/derivatives/tedana"          "${SPIR_SCRATCH}/derivatives/tedana"
copy_component "SPIRITUALITY: afni_firstlvl" "${SPIR_ARCHIVE}/derivatives/afni_firstlvl" "${SPIR_SCRATCH}/derivatives/afni_firstlvl"
copy_component "SPIRITUALITY: afni_group"  "${SPIR_ARCHIVE}/derivatives/afni_group"      "${SPIR_SCRATCH}/derivatives/afni_group"
copy_component "SPIRITUALITY: afni_reproducibility" \
               "${SPIR_ARCHIVE}/derivatives/afni_reproducibility" \
               "${SPIR_SCRATCH}/derivatives/afni_reproducibility"
# NOT copied: derivatives/afni_firstlvl_precap200 — that is the preserved pre-cap-200
# lineage kept deliberately for old-vs-new comparison. It is a historical artifact,
# not a pipeline input, so it stays in archive.

# --- 3. Missionary_language: ARCHIVE -> SCRATCH ------------------------------
copy_component "MISSIONARY_LANGUAGE: BIDS"     "${LUKE_ARCHIVE}/BIDS"                     "${LUKE_SCRATCH}/BIDS"
copy_component "MISSIONARY_LANGUAGE: fmriprep" "${LUKE_ARCHIVE}/derivatives/fmriprep"     "${LUKE_SCRATCH}/derivatives/fmriprep"
copy_component "MISSIONARY_LANGUAGE: tedana"   "${LUKE_ARCHIVE}/derivatives/tedana"       "${LUKE_SCRATCH}/derivatives/tedana"
copy_component "MISSIONARY_LANGUAGE: afni_firstlvl" \
               "${LUKE_ARCHIVE}/derivatives/afni_firstlvl" \
               "${LUKE_SCRATCH}/derivatives/afni_firstlvl"

# --- 4. Summary ---------------------------------------------------------------
banner "SUMMARY"
echo "finished $(date)"
if [ "${FAILED}" -eq 0 ]; then
    echo "All components completed with no rsync failures."
else
    echo "*** ${FAILED} component(s) FAILED — see above. Do NOT delete anything. ***"
fi
echo
echo "NOTHING WAS DELETED. The archive copies are untouched."
echo "Next: verify the copies, run one subject end-to-end on scratch, and only then"
echo "consider retiring the archive originals. Archive has no backups."
exit "${FAILED}"
