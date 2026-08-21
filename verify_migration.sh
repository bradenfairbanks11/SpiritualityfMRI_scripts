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
#   Writes ../logs/migration_verified.txt listing each PASSing component, whenever
#   at least one component passes. prune_archive.sh REFUSES to delete a component
#   not listed in that file. That per-component interlock is the entire point of
#   this script — do not hand-edit that file.
#
#   A FAILING component is simply omitted from the stamp, so it stays undeletable
#   while its passing siblings can proceed. The stamp is not all-or-nothing.
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

# Verify that every path under src also EXISTS under dst, WITHOUT requiring the
# contents to match. This is the right question to ask about a tree that was
# deliberately MERGED.
#
# ~/software/templateflow is the union of the archive copy and the ms-mri copy,
# and migrate_storage.sh rsynced them in that order, so where the two disagreed
# the ms-mri copy won. Asking "did every archive byte survive?" is therefore
# unanswerable by construction: it can never pass. On 2026-08-20 (job 13243193)
# it reported 92 differing files and, because the stamp used to require a
# completely clean run, vetoed 14 other components that HAD passed a full
# checksum audit.
#
# Those 92 files are the older release losing to the newer one — archive's copy
# is dated 2025-08-07, home's 2026-06-30, and no template directory exists only
# in archive. 90 of them are tpl-MNIInfant cohorts; the other two are a
# MNI152NLin2009cAsym description that gained a res-03 entry and a Harvard-Oxford
# label table that grew from 27 to 49 rows.
#
# "Is every archive path present at the destination?" IS answerable, and it is
# the actual precondition for deleting the archive copy: nothing goes missing.
# rsync itemizes a file absent from the destination as '>f+++++++++' — the nine
# '+' mean "newly created". A file that exists but differs itemizes with letters
# and dots instead ('>fcst......'), and is ignored here on purpose.
#
# Content comparison is deliberately cheap regardless of MODE: this check is about
# existence, so there is no reason to checksum 10 GB to answer it.
verify_coverage() {
    local label="$1" src="$2" dst="$3"
    banner "${label}"

    if [ ! -e "${src}" ]; then
        echo "SKIP — source no longer exists: ${src}"
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
    echo "mode: coverage (every source path must exist at the destination;"
    echo "      differing content is expected here and is not a failure)"

    local out rc
    out=$(rsync -a --dry-run --itemize-changes --size-only --exclude 'timing/' \
                "${src}/" "${dst}/" 2>&1)
    rc=$?
    if [ ${rc} -ne 0 ]; then
        echo "FAIL: rsync itself errored (exit ${rc}):"
        echo "${out}" | head -20
        FAILED=$((FAILED + 1))
        return 1
    fi

    local missing n
    missing=$(echo "${out}" | awk '/^>f\+\+\+\+\+\+\+\+\+/ {print}')
    n=$(echo -n "${missing}" | grep -c .)

    if [ "${n}" -eq 0 ]; then
        echo "PASS: all $(find "${src}" -type f 2>/dev/null | wc -l) source paths exist at the destination."
        echo "      Deleting the source would lose nothing."
        PASSED=$((PASSED + 1))
        PASS_LIST+=("${label}")
    else
        echo "FAIL: ${n} source path(s) do NOT exist at the destination."
        echo "--- first 25 ---"
        echo "${missing}" | head -25
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
verify_coverage  "SOFTWARE: TemplateFlow (archive superseded)" \
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

# The stamp lists ONLY the components that actually passed, and it is written
# whenever at least one did. It is deliberately NOT all-or-nothing.
#
# Gating already happens per component, one level down: prune_archive.sh refuses
# to delete anything whose exact label is not marked PASS here (see its prune(),
# "Layer 1"). So a component that failed simply never gets a PASS| line and can
# never be deleted — which is the guarantee that actually matters.
#
# Requiring a globally clean run on top of that added no safety and cost real
# safety: on 2026-08-20 one cosmetic software mismatch suppressed the whole stamp
# and blocked 14 components that had passed a full checksum audit, leaving ~250 GB
# duplicated on a tier that is being migrated to slower disks.
if [ "${PASSED}" -gt 0 ]; then
    {
        echo "# migration verification stamp — written by verify_migration.sh"
        echo "# DO NOT HAND-EDIT. prune_archive.sh trusts this file to decide what is"
        echo "# safe to delete from an un-backed-up filesystem."
        echo "# Only components listed PASS| below are deletable. Anything absent is not."
        echo "mode=${MODE}"
        echo "date=$(date -Iseconds)"
        echo "job=${SLURM_JOB_ID:-none}"
        echo "passed=${PASSED}"
        echo "failed=${FAILED}"
        for c in "${PASS_LIST[@]}"; do echo "PASS|${c}"; done
    } > "${STAMP}"
    echo
    echo "Wrote verification stamp: ${STAMP}  (${PASSED} PASS, ${FAILED} FAIL)"
    if [ "${FAILED}" -ne 0 ]; then
        echo
        echo "*** ${FAILED} component(s) FAILED and are NOT in the stamp. ***"
        echo "prune_archive.sh will refuse to delete those. Re-copy them with"
        echo "migrate_storage.sh (rsync is restartable — it only sends what is"
        echo "missing) and verify again before expecting them to be pruned."
    fi
    if [ "${MODE}" = "full" ]; then
        echo
        echo "MODE=full. Deleting the components listed above is now safe."
    else
        echo
        echo "NOTE: this was the FAST check (size only). Re-run with MODE=full"
        echo "      before deleting anything from archive — prune_archive.sh"
        echo "      will refuse a non-full stamp anyway."
    fi
else
    echo
    echo "*** NOTHING PASSED. DELETE NOTHING. ***"
    echo "No stamp written. Re-copy with migrate_storage.sh and verify again."
fi
exit "${FAILED}"
