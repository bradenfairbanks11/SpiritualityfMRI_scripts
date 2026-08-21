#!/bin/bash
#SBATCH --job-name=verify_migration
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
# Verify the 2026-08-19 storage migration BEFORE anything is deleted from archive.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo.
#
# WHY THIS IS AN SBATCH JOB:
#   ORC policy — "Use Slurm for sustained, computationally intensive, or high-I/O work."
#   MODE=full reads every byte of ~250 GB on BOTH sides. That is emphatically not a
#   login-node activity, especially on an archive tier mid-migration to slower disks.
#
# WHAT IT DOES:
#   For each component that migrate_storage.sh copied, asks one question:
#       "Is every source file present at the destination, byte-for-byte?"
#   It uses `rsync --dry-run --itemize-changes`, which lists exactly the files rsync
#   WOULD still need to transfer. An empty list means the destination is complete.
#
#   Extra files at the destination are FINE and are not flagged. That matters for
#   ~/software/templateflow, which is deliberately the union of two source trees.
#
# MODES:
#   MODE=fast   (default) compare size + modification time. Minutes. Catches missing,
#               truncated, and half-copied files — the realistic rsync failure modes.
#   MODE=full   compare full file checksums. Hours. Catches silent bit corruption too.
#               *** Use MODE=full before deleting anything from archive. ***
#               Archive has no backups, so "probably fine" is not good enough there.
#
# OUTPUT:
#   On a clean run it writes ../logs/migration_verified.txt listing each PASSing
#   component. prune_archive.sh REFUSES to delete a component not listed in that file.
#   That interlock is the entire point of this script — do not hand-edit that file.
#
# USAGE:
#   sbatch verify_migration.sh                        # fast check
#   sbatch --export=ALL,MODE=full verify_migration.sh # the one that gates deletion
# =============================================================================
set -uo pipefail

MODE=${MODE:-fast}

HOME_ROOT=/home/bradenf4
SCRATCH_ROOT=/nobackup/autodelete/usr/bradenf4
ARCHIVE_ROOT=/nobackup/archive/usr/bradenf4

SPIR_ARCHIVE=${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project
SPIR_SCRATCH=${SCRATCH_ROOT}/Spirituality
LUKE_ARCHIVE=${ARCHIVE_ROOT}/Luke_active/Missionary_language/Pilot
LUKE_SCRATCH=${SCRATCH_ROOT}/Missionary_language

LOGS=${HOME_ROOT}/spirituality_fmri/logs
STAMP=${LOGS}/migration_verified.txt

# Same exclude migrate_storage.sh used: timing/ was deliberately NOT copied, so its
# absence at the destination is correct, not a failure.
RSYNC_OPTS=(-a --dry-run --itemize-changes --exclude 'timing/')
case "${MODE}" in
    # --size-only, NOT rsync's default size+mtime comparison. Reason: ~/software/
    # templateflow is deliberately the MERGE of two source trees, and the second
    # rsync stamped its own mtimes over files the first one wrote. Those files are
    # byte-identical but their clocks differ, so a default fast check reports all
    # 2495 of them as "different" — pure noise that buries a real failure.
    # --size-only still catches the failures fast mode exists to catch: missing
    # files and truncated files (truncation changes size). It cannot see a same-size
    # content difference, but neither can any mtime-based check, and that is exactly
    # what MODE=full is for.
    fast) RSYNC_OPTS+=(--size-only) ;;
    # --checksum reads every byte and compares digests. An mtime-only difference then
    # itemizes as '.f..t.....' (leading dot = nothing to transfer), which the parser
    # below ignores — so the merged templateflow passes here without special-casing.
    full) RSYNC_OPTS+=(--checksum) ;;
    *)    echo "MODE must be 'fast' or 'full', got '${MODE}'"; exit 2 ;;
esac

PASSED=0
FAILED=0
PASS_LIST=()

banner() { printf '\n========== %s ==========\n' "$*"; }

# Count files and bytes without a second recursive walk per side.
tree_stats() {
    local p="$1"
    if [ -f "${p}" ]; then
        echo "1 file, $(du -sh "${p}" 2>/dev/null | cut -f1)"
    elif [ -d "${p}" ]; then
        echo "$(find "${p}" -type f 2>/dev/null | wc -l) files, $(du -sh "${p}" 2>/dev/null | cut -f1)"
    else
        echo "MISSING"
    fi
}

verify_component() {
    local label="$1" src="$2" dst="$3"
    banner "${label}"

    if [ ! -e "${src}" ]; then
        echo "SKIP — source no longer exists: ${src}"
        echo "      (nothing to verify; not counted as a failure)"
        return 0
    fi
    if [ ! -e "${dst}" ]; then
        echo "src: ${src}   [$(tree_stats "${src}")]"
        echo "dst: ${dst}"
        echo "FAIL: destination does not exist."
        FAILED=$((FAILED + 1))
        return 1
    fi

    echo "src: ${src}   [$(tree_stats "${src}")]"
    echo "dst: ${dst}   [$(tree_stats "${dst}")]"
    echo "mode: ${MODE}"
    # A lower destination file count is EXPECTED for the afni_firstlvl components:
    # timing/ lives inside that tree and was deliberately excluded from the copy,
    # because its authoritative home is now ~/<project>/timing. Not a shortfall.

    # Trailing slash on src = "contents of", which is how migrate_storage.sh copied it.
    local src_arg="${src}"
    [ -d "${src}" ] && src_arg="${src}/"
    local dst_arg="${dst}"
    [ -d "${src}" ] && dst_arg="${dst}/"

    local out
    out=$(rsync "${RSYNC_OPTS[@]}" "${src_arg}" "${dst_arg}" 2>&1)
    local rc=$?

    if [ ${rc} -ne 0 ]; then
        echo "FAIL: rsync itself errored (exit ${rc}):"
        echo "${out}" | head -20
        FAILED=$((FAILED + 1))
        return 1
    fi

    # Itemized lines look like >f+++++++++ path. Column 2 == 'f' means a regular file
    # rsync would still send, i.e. missing or different at the destination. Directory
    # lines (.d..t....) are timestamp cosmetics on dirs we created with mkdir -p and
    # are deliberately ignored — they say nothing about file contents.
    local diffs
    diffs=$(echo "${out}" | awk '/^[>c<][f]/ {print}')
    local n
    n=$(echo -n "${diffs}" | grep -c . )

    if [ "${n}" -eq 0 ]; then
        echo "PASS: every source file is present and identical at the destination."
        PASSED=$((PASSED + 1))
        PASS_LIST+=("${label}")
    else
        echo "FAIL: ${n} file(s) missing or different at the destination."
        echo "--- first 25 ---"
        echo "${diffs}" | head -25
        FAILED=$((FAILED + 1))
    fi
}

echo "verify_migration.sh   MODE=${MODE}   host=$(hostname)   started $(date)"
echo "job: ${SLURM_JOB_ID:-<none, running outside Slurm>}"
if [ "${MODE}" = "full" ] && [ -z "${SLURM_JOB_ID:-}" ]; then
    echo
    echo "REFUSING TO RUN: MODE=full reads hundreds of GB and must go through Slurm."
    echo "Use:  sbatch --export=ALL,MODE=full verify_migration.sh"
    exit 1
fi

# --- Shared software ----------------------------------------------------------
# templateflow is the union of the archive copy and the ms-mri copy, so it is checked
# once per source. Extra destination files are not flagged, which makes that work.
verify_component "SOFTWARE: TemplateFlow (from archive)" \
                 "${ARCHIVE_ROOT}/software/templateflow" "${HOME_ROOT}/software/templateflow"
verify_component "SOFTWARE: TemplateFlow (from ms-mri)" \
                 "${SCRATCH_ROOT}/ms_dataset/templateflow" "${HOME_ROOT}/software/templateflow"
verify_component "SOFTWARE: fMRIPrep container" \
                 "${ARCHIVE_ROOT}/software/fmri_prep/my_images/fmriprep-25.1.4.sif" \
                 "${HOME_ROOT}/software/fmriprep-25.1.4.sif"
verify_component "SOFTWARE: ms-mri containers" \
                 "${SCRATCH_ROOT}/ms_dataset/containers" "${HOME_ROOT}/software/containers"
verify_component "SOFTWARE: Rlib" \
                 "${SPIR_SCRATCH}/derivatives/afni_group/Rlib" "${HOME_ROOT}/software/Rlib"

# --- Spirituality: archive -> scratch ----------------------------------------
verify_component "SPIRITUALITY: BIDS"     "${SPIR_ARCHIVE}/BIDS"                 "${SPIR_SCRATCH}/BIDS"
verify_component "SPIRITUALITY: fmriprep" "${SPIR_ARCHIVE}/derivatives/fmriprep" "${SPIR_SCRATCH}/derivatives/fmriprep"
verify_component "SPIRITUALITY: tedana"   "${SPIR_ARCHIVE}/derivatives/tedana"   "${SPIR_SCRATCH}/derivatives/tedana"
verify_component "SPIRITUALITY: afni_firstlvl" \
                 "${SPIR_ARCHIVE}/derivatives/afni_firstlvl" "${SPIR_SCRATCH}/derivatives/afni_firstlvl"
verify_component "SPIRITUALITY: afni_group" \
                 "${SPIR_ARCHIVE}/derivatives/afni_group" "${SPIR_SCRATCH}/derivatives/afni_group"
verify_component "SPIRITUALITY: afni_reproducibility" \
                 "${SPIR_ARCHIVE}/derivatives/afni_reproducibility" \
                 "${SPIR_SCRATCH}/derivatives/afni_reproducibility"

# --- Missionary_language: archive -> scratch ---------------------------------
verify_component "MISSIONARY_LANGUAGE: BIDS"     "${LUKE_ARCHIVE}/BIDS"                 "${LUKE_SCRATCH}/BIDS"
verify_component "MISSIONARY_LANGUAGE: fmriprep" "${LUKE_ARCHIVE}/derivatives/fmriprep" "${LUKE_SCRATCH}/derivatives/fmriprep"
verify_component "MISSIONARY_LANGUAGE: tedana"   "${LUKE_ARCHIVE}/derivatives/tedana"   "${LUKE_SCRATCH}/derivatives/tedana"
verify_component "MISSIONARY_LANGUAGE: afni_firstlvl" \
                 "${LUKE_ARCHIVE}/derivatives/afni_firstlvl" "${LUKE_SCRATCH}/derivatives/afni_firstlvl"

# --- Summary ------------------------------------------------------------------
banner "SUMMARY"
echo "finished $(date)"
echo "mode: ${MODE}    passed: ${PASSED}    failed: ${FAILED}"

if [ "${FAILED}" -eq 0 ] && [ "${PASSED}" -gt 0 ]; then
    {
        echo "# migration verification stamp — written by verify_migration.sh"
        echo "# DO NOT HAND-EDIT. prune_archive.sh trusts this file to decide what is"
        echo "# safe to delete from an un-backed-up filesystem."
        echo "mode=${MODE}"
        echo "date=$(date -Iseconds)"
        echo "job=${SLURM_JOB_ID:-none}"
        for c in "${PASS_LIST[@]}"; do echo "PASS|${c}"; done
    } > "${STAMP}"
    echo
    echo "Wrote verification stamp: ${STAMP}"
    if [ "${MODE}" = "full" ]; then
        echo "MODE=full passed. Deleting the verified archive copies is now safe."
    else
        echo "NOTE: this was the FAST check (size + mtime). Re-run with MODE=full"
        echo "      before deleting anything from archive."
    fi
else
    echo
    echo "*** VERIFICATION DID NOT PASS CLEANLY. DELETE NOTHING. ***"
    echo "No stamp written. Re-copy the failing components with migrate_storage.sh"
    echo "(rsync is restartable — it only sends what is missing) and verify again."
fi
exit "${FAILED}"
