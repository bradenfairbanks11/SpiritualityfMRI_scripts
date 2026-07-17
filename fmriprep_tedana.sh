#!/bin/bash

# Full multi-echo fMRI preprocessing pipeline: fMRIPrep + Tedana + MNI warp
# Submits two chained SLURM jobs per participant:
#   Job 1 — fMRIPrep (preprocessing + per-echo output)
#   Job 2 — Tedana ICA denoising + ANTs warp to MNI (runs only after Job 1 succeeds)

# ============================================================
# CONFIGURATION — edit these before running
# ============================================================

BIDS_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/BIDS
FMRIPREP_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/fmriprep
TEDANA_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/tedana

CONDA_ENV=tedenv

# ============================================================

mkdir -p logs

# Participants to process. Pass labels as args (e.g. "sub-04 sub-05" or "04 05");
# with NO args, falls back to ALL sub-* in BIDS (the original behavior).
#
# IDEMPOTENCY: participants whose final MNI tedana outputs already exist are
# skipped. Set FORCE=1 (e.g. `FORCE=1 ./fmriprep_tedana.sh 01`) to reprocess
# them anyway.
#
# SINGLE SESSION: set SES=<label> to process only one session of each listed
# subject (e.g. `SES=2 ./fmriprep_tedana.sh 05` -> only ses-2). fMRIPrep 25.x
# has no --session flag, so this writes a --bids-filter-file restricting the run
# to that session; the completeness guard then counts only that session, and the
# tedana job skips any task that is already done. Leave SES unset for all sessions.
if [ "$#" -gt 0 ]; then
    PARTICIPANT_DIRS=()
    for arg in "$@"; do
        PARTICIPANT_DIRS+=("${BIDS_DIR}/sub-${arg#sub-}")
    done
else
    echo "WARNING: no participant labels given -> defaulting to ALL sub-* in ${BIDS_DIR}"
    PARTICIPANT_DIRS=("${BIDS_DIR}"/sub-*)
fi

for PARTICIPANT_DIR in "${PARTICIPANT_DIRS[@]}"; do
    PARTICIPANT_ID=$(basename "${PARTICIPANT_DIR}")
    PARTICIPANT_LABEL=${PARTICIPANT_ID#sub-}

    # Skip anything that isn't actually a BIDS subject dir (typo'd label, etc.)
    if [ ! -d "${PARTICIPANT_DIR}" ]; then
        echo "ERROR: ${PARTICIPANT_DIR} not found in BIDS -- skipping."
        continue
    fi

    # ----------------------------------------------------------
    # OPTIONAL SESSION RESTRICTION (SES=<label>)
    # fMRIPrep 25.x has no --session flag, so restrict via a PyBIDS
    # --bids-filter-file. Also narrows the completeness guard below to
    # this one session. Leave SES unset to process all sessions.
    # ----------------------------------------------------------
    BIDS_FILTER_ARG=""
    GUARD_SES_GLOB="ses-*"
    WORK_SUBDIR="${PARTICIPANT_ID}"
    if [ -n "${SES:-}" ]; then
        SES_LABEL=${SES#ses-}
        if [ ! -d "${PARTICIPANT_DIR}/ses-${SES_LABEL}" ]; then
            echo "ERROR: ${PARTICIPANT_ID} has no ses-${SES_LABEL} in BIDS -- skipping."
            continue
        fi
        GUARD_SES_GLOB="ses-${SES_LABEL}"
        WORK_SUBDIR="${PARTICIPANT_ID}_ses-${SES_LABEL}"
        FILTER_FILE=${FMRIPREP_OUT}/../bids_filters/${PARTICIPANT_ID}_ses-${SES_LABEL}_filter.json
        mkdir -p "$(dirname "${FILTER_FILE}")"
        cat > "${FILTER_FILE}" <<JSON
{
    "t1w":   {"session": "${SES_LABEL}"},
    "t2w":   {"session": "${SES_LABEL}"},
    "bold":  {"session": "${SES_LABEL}"},
    "sbref": {"session": "${SES_LABEL}"},
    "fmap":  {"session": "${SES_LABEL}"}
}
JSON
        BIDS_FILTER_ARG="--bids-filter-file ${FILTER_FILE}"
        echo "SESSION FILTER: ${PARTICIPANT_ID} restricted to ses-${SES_LABEL} (${FILTER_FILE})"
    fi

    # ----------------------------------------------------------
    # IDEMPOTENCY GUARD
    # Count the session x task units this subject SHOULD produce (one per
    # echo-1 bold in BIDS) and how many already have their final MNI tedana
    # output. If every expected unit is present, skip the subject entirely.
    # Honours SES (counts only that session). Set FORCE=1 to bypass.
    # ----------------------------------------------------------
    if [ "${FORCE:-0}" != "1" ]; then
        EXPECTED_COUNT=0
        DONE_COUNT=0
        for BIDS_ECHO1 in "${PARTICIPANT_DIR}"/${GUARD_SES_GLOB}/func/${PARTICIPANT_ID}_ses-*_task-*_echo-1_bold.nii.gz; do
            [ -f "${BIDS_ECHO1}" ] || continue
            EXPECTED_COUNT=$((EXPECTED_COUNT + 1))
            FN=$(basename "${BIDS_ECHO1}")
            SES=$(echo "${FN}" | grep -oP 'ses-\K[^_]+')
            TASK=$(echo "${FN}" | grep -oP 'task-\K[^_]+')
            MNI_OUT=${TEDANA_OUT}/${PARTICIPANT_ID}/ses-${SES}/func/${PARTICIPANT_ID}_ses-${SES}_task-${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
            [ -f "${MNI_OUT}" ] && DONE_COUNT=$((DONE_COUNT + 1))
        done

        if [ "${EXPECTED_COUNT}" -gt 0 ] && [ "${DONE_COUNT}" -eq "${EXPECTED_COUNT}" ]; then
            echo "SKIP: ${PARTICIPANT_ID} already complete (${DONE_COUNT}/${EXPECTED_COUNT} MNI tedana outputs present). Use FORCE=1 to reprocess."
            continue
        fi
        if [ "${DONE_COUNT}" -gt 0 ]; then
            echo "NOTE: ${PARTICIPANT_ID} partially done (${DONE_COUNT}/${EXPECTED_COUNT}); resubmitting -- finished tasks will be skipped inside the tedana job."
        fi
    fi

    # ----------------------------------------------------------
    # JOB 1: fMRIPrep
    # --parsable makes sbatch print only the job ID number so
    # we can pass it to Job 2 as a dependency.
    # ----------------------------------------------------------
    FMRIPREP_JOB=$(sbatch --parsable <<EOT
#!/bin/bash
#SBATCH --job-name=${PARTICIPANT_ID}_fmriprep
#SBATCH --output=logs/${PARTICIPANT_ID}_fmriprep.out
#SBATCH --error=logs/${PARTICIPANT_ID}_fmriprep.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com

set -euo pipefail

module load apptainer/1.3.6-qycanb2

export TEMPLATEFLOW_HOME=/nobackup/archive/usr/bradenf4/software/templateflow

apptainer run --cleanenv \
    --env TEMPLATEFLOW_HOME=/nobackup/archive/usr/bradenf4/software/templateflow \
    --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4,/nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4,/nobackup/archive/usr/bradenf4/software/templateflow:/nobackup/archive/usr/bradenf4/software/templateflow \
    /nobackup/archive/usr/bradenf4/software/fmri_prep/my_images/fmriprep-25.1.4.sif \
    ${BIDS_DIR} \
    ${FMRIPREP_OUT} \
    participant \
    --participant-label ${PARTICIPANT_LABEL} \
    ${BIDS_FILTER_ARG} \
    --fs-license-file /nobackup/archive/usr/bradenf4/software/fmri_prep/preprocessing/license.txt \
    --work-dir /nobackup/autodelete/usr/bradenf4/${WORK_SUBDIR} \
    --nthreads 8 \
    --mem 64G \
    --me-output-echos \
    --me-t2s-fit-method curvefit
EOT
    )

    echo "Submitted fMRIPrep job ${FMRIPREP_JOB} for ${PARTICIPANT_ID}"

    # ----------------------------------------------------------
    # JOB 2: Tedana + MNI warp
    # --dependency=afterok means this job only starts if Job 1
    # exits with code 0 (success). If fMRIPrep fails, this is
    # automatically cancelled.
    # ----------------------------------------------------------
    sbatch --dependency=afterok:${FMRIPREP_JOB} <<EOT
#!/bin/bash
#SBATCH --job-name=${PARTICIPANT_ID}_tedana
#SBATCH --output=logs/${PARTICIPANT_ID}_tedana.out
#SBATCH --error=logs/${PARTICIPANT_ID}_tedana.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com

set -euo pipefail

# -------------------------------------------------------
# STEP 0: Load tools
# Tedana is in the conda environment; ANTs is a cluster module.
# -------------------------------------------------------
source activate ${CONDA_ENV}
# Ignore ~/.local site-packages: a stray user-site numpy 2.5.0 (installed 2026-06-30)
# shadows the env's numpy 2.4.6 and removed np.row_stack, which mapca still calls -> tedana crash.
export PYTHONNOUSERSITE=1
module load ants/2.5.2-nzivpdx
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

# -------------------------------------------------------
# STEP 1: Loop over sessions found in the fMRIPrep output
# -------------------------------------------------------
for SES_DIR in ${FMRIPREP_OUT}/${PARTICIPANT_ID}/ses-*/; do
    [ -d "\${SES_DIR}" ] || continue
    SES=\$(basename "\${SES_DIR}")
    FMRIPREP_FUNC=\${SES_DIR}func

    # -------------------------------------------------------
    # STEP 2: Loop over tasks by detecting echo-1 output files.
    # This automatically handles whichever tasks fMRIPrep produced
    # without needing a hardcoded task list.
    # -------------------------------------------------------
    for ECHO1_FILE in \$(ls \${FMRIPREP_FUNC}/${PARTICIPANT_ID}_\${SES}_task-*_echo-1_desc-preproc_bold.nii.gz 2>/dev/null | sort); do
        TASK=\$(basename "\${ECHO1_FILE}" | grep -oP 'task-\K[^_]+')

        # -------------------------------------------------------
        # IDEMPOTENCY: skip this session/task if its final MNI tedana
        # output already exists (FORCE=1 bypasses this at submit time).
        # -------------------------------------------------------
        EXISTING_MNI=${TEDANA_OUT}/${PARTICIPANT_ID}/\${SES}/func/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
        if [ "${FORCE:-0}" != "1" ] && [ -f "\${EXISTING_MNI}" ]; then
            echo "SKIP: ${PARTICIPANT_ID} \${SES} task-\${TASK} already has MNI tedana output."
            continue
        fi

        # -------------------------------------------------------
        # STEP 3: Collect all echo files for this session/task
        # Sorted to ensure ascending echo-time order for tedana.
        # -------------------------------------------------------
        ECHO_FILES=(\$(ls \${FMRIPREP_FUNC}/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_echo-*_desc-preproc_bold.nii.gz 2>/dev/null | sort))

        echo "Found \${#ECHO_FILES[@]} echo files for ${PARTICIPANT_ID} \${SES} task-\${TASK}"

        # -------------------------------------------------------
        # STEP 4: Extract echo times from BIDS JSON sidecar files.
        # BIDS stores EchoTime in seconds; tedana expects seconds.
        # -------------------------------------------------------
        ECHO_TIMES=()

        for ECHO_FILE in "\${ECHO_FILES[@]}"; do
            ECHO_NUM=\$(basename "\${ECHO_FILE}" | grep -oP 'echo-\K[0-9]+')
            JSON_FILE=${BIDS_DIR}/${PARTICIPANT_ID}/\${SES}/func/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_echo-\${ECHO_NUM}_bold.json

            if [ ! -f "\${JSON_FILE}" ]; then
                echo "ERROR: No JSON sidecar found for ${PARTICIPANT_ID} \${SES} task-\${TASK} echo-\${ECHO_NUM}"
                exit 1
            fi

            TE=\$(python3 -c "import json; d=json.load(open('\${JSON_FILE}')); print(round(d['EchoTime'], 4))")
            ECHO_TIMES+=(\$TE)
        done

        echo "Echo times (s): \${ECHO_TIMES[@]}"

        # -------------------------------------------------------
        # STEP 5: Run tedana ICA denoising
        # Input:  per-echo files + echo times
        # Output: desc-denoised_bold.nii.gz in native space
        # -------------------------------------------------------
        TEDANA_FUNC=${TEDANA_OUT}/${PARTICIPANT_ID}/\${SES}/func
        mkdir -p \${TEDANA_FUNC}

        BRAIN_MASK=\${FMRIPREP_FUNC}/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_desc-brain_mask.nii.gz

        tedana \
            -d "\${ECHO_FILES[@]}" \
            -e "\${ECHO_TIMES[@]}" \
            --fittype curvefit \
            --out-dir \${TEDANA_FUNC} \
            --prefix ${PARTICIPANT_ID}_\${SES}_task-\${TASK} \
            --mask \${BRAIN_MASK} \
            --n-threads 8 \
            --seed 42 \
            --overwrite

        # -------------------------------------------------------
        # STEP 6: Warp denoised output to MNI standard space
        # -------------------------------------------------------
        DENOISED_NATIVE=\${TEDANA_FUNC}/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_desc-denoised_bold.nii.gz
        if [ ! -f "\${DENOISED_NATIVE}" ]; then
            echo "ERROR: Tedana output not found: \${DENOISED_NATIVE}"
            exit 1
        fi

        MNI_REF=${FMRIPREP_OUT}/${PARTICIPANT_ID}/\${SES}/func/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_space-MNI152NLin2009cAsym_boldref.nii.gz
        COREG_XFM=${FMRIPREP_OUT}/${PARTICIPANT_ID}/\${SES}/func/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt
        # Resolve the T1w->MNI warp. fMRIPrep writes it to the session-agnostic
        # sub-XX/anat/ for MULTI-session subjects (with a run entity, e.g.
        # sub-01_run-1_from-T1w_to-MNI...) but nests it under sub-XX/ses-*/anat/ for
        # SINGLE-session subjects (e.g. sub-04/05). Glob both layouts; the from-T1w_to-MNI
        # pattern tolerates the run label and excludes the inverse (from-MNI_to-T1w) warp.
        NORM_XFM=\$(ls ${FMRIPREP_OUT}/${PARTICIPANT_ID}/anat/${PARTICIPANT_ID}_*from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 \
                       ${FMRIPREP_OUT}/${PARTICIPANT_ID}/ses-*/anat/${PARTICIPANT_ID}_*from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 \
                       2>/dev/null | head -n1 || true)
        DENOISED_MNI=\${TEDANA_FUNC}/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz

        # Verify all transform files exist before calling ANTs
        for f in "\${MNI_REF}" "\${COREG_XFM}" "\${NORM_XFM}"; do
            [ -f "\$f" ] || { echo "ERROR: Required ANTs file not found: \$f"; exit 1; }
        done

        antsApplyTransforms \
            -d 3 \
            -e 3 \
            -i \${DENOISED_NATIVE} \
            -r \${MNI_REF} \
            -t \${COREG_XFM} \
            -t \${NORM_XFM} \
            -o \${DENOISED_MNI} \
            -n LanczosWindowedSinc

        echo "Done: ${PARTICIPANT_ID} \${SES} task-\${TASK} — MNI-space file at \${DENOISED_MNI}"
    done
done
EOT

done
