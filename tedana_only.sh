#!/bin/bash

# Re-run ONLY the tedana + MNI-warp step (NO fMRIPrep) for participants whose
# fMRIPrep derivatives already exist. Submits one SLURM job per participant.
#
#   Usage:  bash tedana_only.sh sub-04 sub-05    (or:  bash tedana_only.sh 04 05)
#
# The job body below is kept identical to the tedana block in fmriprep_tedana.sh
# (same PYTHONNOUSERSITE guard and same broadened T1w->MNI glob). Use this to
# recover subjects whose fMRIPrep succeeded but whose tedana step failed, without
# paying for a 10-12 h fMRIPrep re-run.

# ============================================================
# CONFIGURATION
# ============================================================
BIDS_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/BIDS
FMRIPREP_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/fmriprep
TEDANA_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/tedana

CONDA_ENV=tedenv
# ============================================================

mkdir -p logs

if [ "$#" -eq 0 ]; then
    echo "ERROR: pass participant labels, e.g. 'bash tedana_only.sh sub-04 sub-05'"
    exit 1
fi

for arg in "$@"; do
    PARTICIPANT_ID=sub-${arg#sub-}
    PARTICIPANT_DIR=${FMRIPREP_OUT}/${PARTICIPANT_ID}

    # Skip anything without existing fMRIPrep output (typo'd label, not yet run, etc.)
    if [ ! -d "${PARTICIPANT_DIR}" ]; then
        echo "ERROR: no fMRIPrep output at ${PARTICIPANT_DIR} -- skipping ${PARTICIPANT_ID}."
        continue
    fi

    sbatch <<EOT
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

# STEP 0: tools. Tedana is in the conda env; ANTs is a cluster module.
source activate ${CONDA_ENV}
# Ignore ~/.local site-packages: a stray user-site numpy 2.5.0 (installed 2026-06-30)
# shadows the env's numpy 2.4.6 and removed np.row_stack, which mapca still calls -> tedana crash.
export PYTHONNOUSERSITE=1
module load ants/2.5.2-nzivpdx
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

# STEP 1: loop over sessions found in the fMRIPrep output
for SES_DIR in ${FMRIPREP_OUT}/${PARTICIPANT_ID}/ses-*/; do
    [ -d "\${SES_DIR}" ] || continue
    SES=\$(basename "\${SES_DIR}")
    FMRIPREP_FUNC=\${SES_DIR}func

    # STEP 2: loop over tasks by detecting echo-1 output files
    for ECHO1_FILE in \$(ls \${FMRIPREP_FUNC}/${PARTICIPANT_ID}_\${SES}_task-*_echo-1_desc-preproc_bold.nii.gz 2>/dev/null | sort); do
        TASK=\$(basename "\${ECHO1_FILE}" | grep -oP 'task-\K[^_]+')

        # STEP 3: collect all echo files for this session/task (ascending echo order)
        ECHO_FILES=(\$(ls \${FMRIPREP_FUNC}/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_echo-*_desc-preproc_bold.nii.gz 2>/dev/null | sort))
        echo "Found \${#ECHO_FILES[@]} echo files for ${PARTICIPANT_ID} \${SES} task-\${TASK}"

        # STEP 4: extract echo times (seconds) from BIDS JSON sidecars
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

        # STEP 5: run tedana ICA denoising (--overwrite so partial prior output is replaced)
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

        # STEP 6: warp denoised output to MNI standard space
        DENOISED_NATIVE=\${TEDANA_FUNC}/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_desc-denoised_bold.nii.gz
        if [ ! -f "\${DENOISED_NATIVE}" ]; then
            echo "ERROR: Tedana output not found: \${DENOISED_NATIVE}"
            exit 1
        fi

        MNI_REF=${FMRIPREP_OUT}/${PARTICIPANT_ID}/\${SES}/func/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_space-MNI152NLin2009cAsym_boldref.nii.gz
        COREG_XFM=${FMRIPREP_OUT}/${PARTICIPANT_ID}/\${SES}/func/${PARTICIPANT_ID}_\${SES}_task-\${TASK}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt
        # T1w->MNI warp: sub-XX/anat/ for multi-session subjects, sub-XX/ses-*/anat/ for
        # single-session subjects (e.g. sub-04/05). Glob both; the from-T1w_to-MNI pattern
        # tolerates the run label and excludes the inverse (from-MNI_to-T1w) warp.
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

    echo "Submitted tedana-only job for ${PARTICIPANT_ID}"

done
