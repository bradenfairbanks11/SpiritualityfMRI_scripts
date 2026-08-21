#!/bin/bash
#SBATCH --job-name=prune_archive
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
# Reduce /nobackup/archive to its intended contents: RAW DATA + AFNI DERIVATIVES.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo.
#
# *** THIS IS THE ONLY SCRIPT IN THIS REPO THAT DELETES DATA. ***
# *** ARCHIVE HAS NO BACKUPS AND NO SNAPSHOTS. THERE IS NO UNDO. ***
#
# THE INTENDED END STATE:
#   ARCHIVE keeps  rawdata/   — the DICOMs. Unregenerable. The only true originals.
#                  derivatives/afni_*  — the analysis products worth keeping forever.
#   ARCHIVE drops  BIDS/      — regenerable from rawdata via dcm2niix.sh
#                  derivatives/fmriprep, derivatives/tedana
#                             — regenerable from BIDS, and already copied to scratch
#
# WHY THIS IS SAFE, IN FOUR LAYERS:
#   1. It refuses to delete anything without a MODE=full verification stamp from
#      verify_migration.sh proving the data exists elsewhere, byte-for-byte.
#   2. It RE-CHECKS each copy immediately before deleting it. The stamp could be days
#      old, and the destination is SCRATCH — which purges after 12 weeks unused. A
#      stale stamp must never authorize deleting the last copy of something.
#   3. DRYRUN=1 is the default. Deleting requires typing DRYRUN=0 on purpose.
#   4. A real delete must go through Slurm. `rm -rf` over hundreds of thousands of
#      files is heavy metadata I/O (ORC rules 2 and 8), not a login-node activity.
#
# WHY THIS IS AN SBATCH JOB, NOT `rm -rf` BY HAND:
#   Because a typo in an interactive `rm -rf` on this filesystem is unrecoverable, and
#   because every path here is derived, checked, and logged before it is touched.
#
# USAGE:
#   bash prune_archive.sh                   # dry run: shows exactly what would go
#   sbatch --export=ALL,DRYRUN=0 prune_archive.sh   # the real, irreversible run
#
# OPTIONAL EXTRAS (both default OFF — flip only if you mean it):
#   PRUNE_SCRIPTS=1  also delete the stale scripts/ snapshots in archive. The live
#                    copies are git repos under ~/, which IS backed up. Off by default
#                    because it is provenance, it is tiny, and it costs nothing to keep.
#   PRUNE_SOFTWARE=1 also delete archive software/ (containers, templateflow). Copies
#                    now live in ~/software. Off by default: it is outside the "rawdata
#                    + AFNI" rule you asked for and deserves its own decision.
# =============================================================================
set -uo pipefail

DRYRUN=${DRYRUN:-1}
PRUNE_SCRIPTS=${PRUNE_SCRIPTS:-0}
PRUNE_SOFTWARE=${PRUNE_SOFTWARE:-0}

HOME_ROOT=/home/bradenf4
SCRATCH_ROOT=/nobackup/autodelete/usr/bradenf4
ARCHIVE_ROOT=/nobackup/archive/usr/bradenf4

SPIR_ARCHIVE=${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project
SPIR_SCRATCH=${SCRATCH_ROOT}/Spirituality
LUKE_ARCHIVE=${ARCHIVE_ROOT}/Luke_active/Missionary_language/Pilot
LUKE_SCRATCH=${SCRATCH_ROOT}/Missionary_language

LOGS=${HOME_ROOT}/spirituality_fmri/logs
STAMP=${LOGS}/migration_verified.txt

DELETED=0
SKIPPED=0
RECLAIM_REPORT=()

banner() { printf '\n========== %s ==========\n' "$*"; }

# --- Guard rails --------------------------------------------------------------
echo "prune_archive.sh   DRYRUN=${DRYRUN}   host=$(hostname)   started $(date)"
echo "job: ${SLURM_JOB_ID:-<none, running outside Slurm>}"

if [ "${DRYRUN}" != "1" ] && [ -z "${SLURM_JOB_ID:-}" ]; then
    echo
    echo "REFUSING TO RUN: a real deletion must go through Slurm, not a login node."
    echo "Use:  sbatch --export=ALL,DRYRUN=0 prune_archive.sh"
    exit 1
fi

if [ ! -f "${STAMP}" ]; then
    echo
    echo "REFUSING TO RUN: no verification stamp at ${STAMP}"
    echo "Run this first, and let it pass:"
    echo "    sbatch --export=ALL,MODE=full verify_migration.sh"
    exit 1
fi

STAMP_MODE=$(grep -E '^mode=' "${STAMP}" | head -1 | cut -d= -f2)
STAMP_DATE=$(grep -E '^date=' "${STAMP}" | head -1 | cut -d= -f2-)
echo "stamp: mode=${STAMP_MODE}  written=${STAMP_DATE}"

if [ "${DRYRUN}" != "1" ] && [ "${STAMP_MODE}" != "full" ]; then
    echo
    echo "REFUSING TO DELETE: the stamp is from a '${STAMP_MODE}' verification."
    echo "A fast check compares size and mtime only — it cannot see silent corruption."
    echo "Archive has no backups, so deletion requires the full checksum pass:"
    echo "    sbatch --export=ALL,MODE=full verify_migration.sh"
    exit 1
fi

# --- The one function that removes anything -----------------------------------
# victim      : the archive path to delete
# proof_label : the verify_migration.sh component that proves it exists elsewhere
# proof_copy  : where that surviving copy lives (re-checked here, right now)
prune() {
    local label="$1" victim="$2" proof_label="$3" proof_copy="$4"
    banner "${label}"

    if [ ! -e "${victim}" ]; then
        echo "SKIP — already gone: ${victim}"
        return 0
    fi
    echo "delete : ${victim}"
    echo "size   : $(du -sh "${victim}" 2>/dev/null | cut -f1)"
    echo "proof  : ${proof_label}"
    echo "copy at: ${proof_copy}"

    # Layer 1 — the stamp must vouch for this exact component.
    if ! grep -qF "PASS|${proof_label}" "${STAMP}"; then
        echo "SKIP — '${proof_label}' is not marked PASS in the stamp. Not deleting."
        SKIPPED=$((SKIPPED + 1))
        return 1
    fi

    # Layer 2 — the surviving copy must still be there RIGHT NOW. Scratch purges.
    if [ ! -e "${proof_copy}" ]; then
        echo "SKIP — the surviving copy has VANISHED from ${proof_copy}."
        echo "       Scratch purges after 12 weeks unused. Deleting now would destroy"
        echo "       the last copy. Re-run migrate_storage.sh, then verify, then retry."
        SKIPPED=$((SKIPPED + 1))
        return 1
    fi

    # Layer 3 — and it must still be complete. Cheap size+mtime pass; the stamp
    # already covered checksums, this is about drift since the stamp was written.
    local src_arg="${victim}" dst_arg="${proof_copy}"
    if [ -d "${victim}" ]; then src_arg="${victim}/"; dst_arg="${proof_copy}/"; fi
    local n
    n=$(rsync -a --dry-run --itemize-changes --exclude 'timing/' "${src_arg}" "${dst_arg}" 2>/dev/null \
        | awk '/^[>c<][f]/' | grep -c .)
    if [ "${n}" -ne 0 ]; then
        echo "SKIP — the copy is no longer complete: ${n} file(s) missing or changed."
        echo "       Not deleting. Re-copy and re-verify first."
        SKIPPED=$((SKIPPED + 1))
        return 1
    fi
    echo "recheck: OK — copy is complete as of right now."

    if [ "${DRYRUN}" = "1" ]; then
        echo "[dry run] would delete ${victim}"
    else
        echo "DELETING ${victim} ..."
        if rm -rf "${victim}"; then
            echo "deleted."
        else
            echo "ERROR: rm failed for ${victim}"
            SKIPPED=$((SKIPPED + 1))
            return 1
        fi
    fi
    DELETED=$((DELETED + 1))
    return 0
}

# --- Spirituality --------------------------------------------------------------
prune "SPIRITUALITY: BIDS (regenerable from rawdata)" \
      "${SPIR_ARCHIVE}/BIDS" "SPIRITUALITY: BIDS" "${SPIR_SCRATCH}/BIDS"

prune "SPIRITUALITY: derivatives/fmriprep" \
      "${SPIR_ARCHIVE}/derivatives/fmriprep" "SPIRITUALITY: fmriprep" \
      "${SPIR_SCRATCH}/derivatives/fmriprep"

prune "SPIRITUALITY: derivatives/tedana" \
      "${SPIR_ARCHIVE}/derivatives/tedana" "SPIRITUALITY: tedana" \
      "${SPIR_SCRATCH}/derivatives/tedana"

# bids_filters is empty on both sides — tiny, generated, and needs no proof.
banner "SPIRITUALITY: derivatives/bids_filters (empty, generated)"
if [ -d "${SPIR_ARCHIVE}/derivatives/bids_filters" ]; then
    if [ -z "$(ls -A "${SPIR_ARCHIVE}/derivatives/bids_filters" 2>/dev/null)" ]; then
        if [ "${DRYRUN}" = "1" ]; then
            echo "[dry run] would remove empty dir"
        else
            rmdir "${SPIR_ARCHIVE}/derivatives/bids_filters" && echo "removed empty dir"
        fi
    else
        echo "SKIP — not empty after all; leaving it alone."
    fi
else
    echo "SKIP — already gone."
fi

# --- Missionary_language -------------------------------------------------------
prune "MISSIONARY_LANGUAGE: BIDS (regenerable from rawdata)" \
      "${LUKE_ARCHIVE}/BIDS" "MISSIONARY_LANGUAGE: BIDS" "${LUKE_SCRATCH}/BIDS"

prune "MISSIONARY_LANGUAGE: derivatives/fmriprep" \
      "${LUKE_ARCHIVE}/derivatives/fmriprep" "MISSIONARY_LANGUAGE: fmriprep" \
      "${LUKE_SCRATCH}/derivatives/fmriprep"

prune "MISSIONARY_LANGUAGE: derivatives/tedana" \
      "${LUKE_ARCHIVE}/derivatives/tedana" "MISSIONARY_LANGUAGE: tedana" \
      "${LUKE_SCRATCH}/derivatives/tedana"

# --- Opt-in extras -------------------------------------------------------------
if [ "${PRUNE_SCRIPTS}" = "1" ]; then
    banner "EXTRA: stale scripts/ snapshots in archive"
    for d in "${SPIR_ARCHIVE}/scripts" "${LUKE_ARCHIVE}/scripts"; do
        [ -e "${d}" ] || { echo "SKIP — already gone: ${d}"; continue; }
        echo "delete: ${d}  ($(du -sh "${d}" 2>/dev/null | cut -f1))"
        if [ "${DRYRUN}" = "1" ]; then echo "[dry run]"; else rm -rf "${d}" && echo "deleted."; fi
    done
else
    banner "EXTRA: scripts/ — NOT touched (PRUNE_SCRIPTS=0)"
    echo "Stale snapshots kept as provenance. Live copies are git repos under ~/."
fi

if [ "${PRUNE_SOFTWARE}" = "1" ]; then
    banner "EXTRA: archive software/"
    echo "delete: ${ARCHIVE_ROOT}/software  ($(du -sh "${ARCHIVE_ROOT}/software" 2>/dev/null | cut -f1))"
    if ! grep -qF "PASS|SOFTWARE: fMRIPrep container" "${STAMP}" \
       || ! grep -qF "PASS|SOFTWARE: TemplateFlow (from archive)" "${STAMP}"; then
        echo "SKIP — software components are not both PASS in the stamp."
        SKIPPED=$((SKIPPED + 1))
    elif [ "${DRYRUN}" = "1" ]; then
        echo "[dry run] would delete archive software/"
    else
        rm -rf "${ARCHIVE_ROOT}/software" && { echo "deleted."; DELETED=$((DELETED + 1)); }
    fi
else
    banner "EXTRA: archive software/ — NOT touched (PRUNE_SOFTWARE=0)"
fi

# --- Summary --------------------------------------------------------------------
banner "SUMMARY"
echo "finished $(date)"
if [ "${DRYRUN}" = "1" ]; then
    echo "DRY RUN — nothing was deleted. ${DELETED} component(s) would be removed,"
    echo "${SKIPPED} refused. Re-run with DRYRUN=0 under sbatch to do it for real."
else
    echo "${DELETED} component(s) deleted, ${SKIPPED} refused."
fi
echo
echo "What archive should hold now:"
echo "  rawdata/                    the DICOMs — unregenerable originals"
echo "  derivatives/afni_*          the analysis products worth keeping"
echo
echo "What is now regenerable-only (scratch copies, which PURGE after 12 weeks unused):"
echo "  BIDS, fmriprep, tedana"
echo "  -> If you park either project, run promote_to_archive.sh first."
exit "${SKIPPED}"
