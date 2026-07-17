#!/bin/bash
#SBATCH --job-name=spirituality_firstlvl
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com

# =============================================================================
# Spirituality fMRI — First-Level GLM with PARAMETRIC MODULATION (AFNI)
#
# On every trial the participant rates 1-4 (button box) how strongly they "felt
# the Spirit" (1 = lowest, 4 = highest). This script models each stimulus event
# amplitude-modulated by that rating (AFNI -stim_times_AM2), so the GLM estimates
# how the BOLD response scales with rating.
#
# One GLM per participant / session / task. Each GLM has:
#   6 motion nuisance regressors (-stim_base)
#   1 amplitude-modulated stimulus regressor (AM2 -> two betas):
#       rating[0] = mean response (overall activation)
#       rating[1] = parametric SLOPE (BOLD vs rating) <-- result of interest
#
# Pipeline per run: smooth -> scale -> motion -> 3dDeconvolve -> 3dREMLfit
#
# Timing files (made by generate_timing_files.py, AM format "onset*rating"):
#   ${TIMING_DIR}/byv-<BYV>_<ses>_task-<TASK>_ratingAM.1D
# Runs with no usable ratings have NO timing file and are skipped automatically
# (e.g. 0AK1RT ses-1 task-FHS).
#
# Participant mapping: BIDS uses sub-NN (assigned by scan-date order in
# dcm2niix), but timing files are keyed by the BYV alphanumeric ID. The script
# resolves sub-NN -> BYV via participants.tsv (see PARTICIPANTS_TSV).
#
# Run with:  sbatch first_level_afni.sh
# =============================================================================

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# --- Paths already known from the preprocessing pipeline ---
BIDS_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/BIDS
RAWDATA_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata
FMRIPREP_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/fmriprep
TEDANA_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/tedana

# --- Paths still needed ---
AFNI_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/afni_firstlvl                       
TIMING_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/afni_firstlvl/timing                   
AFNI_SIF=/apps/afni/afni_make_build_latest.sif  
PARTICIPANTS_TSV=${RAWDATA_DIR}/participants.tsv  
BIND=/nobackup/archive/usr/bradenf4      

# --- Acquisition parameters (known) ---
TR=1.388
SMOOTH_FWHM=4.0
MOTION_THRESH=0.3

# --- Tasks (BIDS labels: F=FHS, S=scripture, A=architecture per dcm2niix) ---
TASKS=(scripture FHS architecture)

# --- Fixed BLOCK stimulus duration per task (seconds) ---
# From generate_timing_files.py duration QC (median stimulus on-screen time):
#   scripture median 19.3s | FHS median 44.3s | architecture median 14.2s
block_dur_for() {
    case "$1" in
        scripture)    echo 19 ;;
        FHS)          echo 44 ;;
        architecture) echo 14 ;;
        *)            echo "ERROR: unknown task $1" >&2; exit 1 ;;
    esac
}

# ==============================================================================

mkdir -p logs "${AFNI_OUT}"

module load apptainer/1.3.6-qycanb2

AFNI="apptainer exec --bind ${BIND}:${BIND} ${AFNI_SIF} bash -c"

if [ ! -f "${PARTICIPANTS_TSV}" ]; then
    echo "ERROR: participants.tsv not found at ${PARTICIPANTS_TSV}"
    echo "       It maps BIDS sub-NN to the BYV id used in timing filenames."
    echo "       Generate it during dcm2niix conversion (see BIDS/dcm2niix.txt)."
    exit 1
fi

# Resolve a BIDS sub-NN label to its BYV alphanumeric id via participants.tsv.
# Expects a header row containing columns 'participant_id' and 'byv_id'
# (participant_id values like sub-01; byv_id values like 0AK1RT).
byv_for_sub() {
    local sub="$1"
    # NOTE: the awk variable is 'want', not 'sub' -- sub() is a built-in awk
    # function and reusing the name is a syntax error on some awk variants.
    awk -F'\t' -v want="${sub}" '
        NR==1 { for (i=1;i<=NF;i++){ if($i=="participant_id")p=i; if($i=="byv_id")b=i}; next }
        $p==want { print $b; exit }
    ' "${PARTICIPANTS_TSV}"
}

# ==============================================================================
# LOOP: participants -> sessions -> tasks
# ==============================================================================

# Participants to process. Pass labels as args (e.g. "sub-04 sub-05" or "04 05");
# with NO args, falls back to ALL sub-* in the tedana derivatives (original behavior).
# Under SLURM, pass them after the script name:  sbatch first_level_afni.sh sub-04 sub-05
if [ "$#" -gt 0 ]; then
    PARTICIPANT_DIRS=()
    for arg in "$@"; do
        PARTICIPANT_DIRS+=("${TEDANA_OUT}/sub-${arg#sub-}")
    done
else
    echo "WARNING: no participant labels given -> defaulting to ALL sub-* in ${TEDANA_OUT}"
    PARTICIPANT_DIRS=("${TEDANA_OUT}"/sub-*)
fi

for PARTICIPANT_DIR in "${PARTICIPANT_DIRS[@]}"; do
    if [ ! -d "${PARTICIPANT_DIR}" ]; then
        echo "ERROR: ${PARTICIPANT_DIR} not found in tedana output -- skipping."
        continue
    fi
    PARTICIPANT_ID=$(basename "${PARTICIPANT_DIR}")

    BYV=$(byv_for_sub "${PARTICIPANT_ID}")
    if [ -z "${BYV}" ]; then
        echo "Skipping ${PARTICIPANT_ID}: no byv_id in participants.tsv"
        continue
    fi

    for SES_DIR in "${PARTICIPANT_DIR}"/ses-*/; do
        [ -d "${SES_DIR}" ] || continue
        SES=$(basename "${SES_DIR}")
        TEDANA_FUNC=${SES_DIR}func

        for TASK in "${TASKS[@]}"; do

            BOLD_IN=${TEDANA_FUNC}/${PARTICIPANT_ID}_${SES}_task-${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
            CONFOUNDS=${FMRIPREP_OUT}/${PARTICIPANT_ID}/${SES}/func/${PARTICIPANT_ID}_${SES}_task-${TASK}_desc-confounds_timeseries.tsv
            TIMING=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_ratingAM.1D

            if [ ! -f "${BOLD_IN}" ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: BOLD not found"
                continue
            fi
            if [ ! -f "${TIMING}" ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: no timing file (dead/absent run) ${TIMING}"
                continue
            fi

            OUT_DIR=${AFNI_OUT}/${PARTICIPANT_ID}/${SES}/task-${TASK}
            mkdir -p "${OUT_DIR}"
            PREFIX=${PARTICIPANT_ID}_${SES}_task-${TASK}
            BLOCK_DUR=$(block_dur_for "${TASK}")

            echo "=== ${PARTICIPANT_ID} ${SES} task-${TASK} (byv-${BYV}, BLOCK ${BLOCK_DUR}s) ==="

            # ------------------------------------------------------------------
            # STEP 1: Spatial smoothing (4mm FWHM)
            # ------------------------------------------------------------------
            ${AFNI} "
                3dmerge \
                    -1blur_fwhm ${SMOOTH_FWHM} \
                    -doall \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred.nii.gz \
                    ${BOLD_IN}
            "

            # ------------------------------------------------------------------
            # STEP 2: Scale to percent signal change (cap at 200)
            # ------------------------------------------------------------------
            ${AFNI} "
                3dTstat \
                    -mean \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred_mean.nii.gz \
                    ${OUT_DIR}/${PREFIX}_blurred.nii.gz

                3dcalc \
                     e- ${OUT_DIR}/${PREFIX}_blurred.nii.gz \
                    -b ${OUT_DIR}/${PREFIX}_blurred_mean.nii.gz \
                    -expr 'min(200, a/b*100)*step(a)*step(b)' \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz
            "

            # ------------------------------------------------------------------
            # STEP 3: Extract 6 motion parameters from fMRIPrep confounds TSV.
            # Python parser (not 1dcat) so 'n/a' cells -> 0 (fMRIPrep writes n/a
            # in the first row of motion derivatives).
            # ------------------------------------------------------------------
            python3 -c "
import csv
with open('${CONFOUNDS}') as f:
    reader = csv.DictReader(f, delimiter='\t')
    cols = ['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']
    for row in reader:
        vals = ['0' if row[c] in ('n/a','') else row[c] for c in cols]
        print(' '.join(vals))
" > ${OUT_DIR}/${PREFIX}_motion.txt

            # ------------------------------------------------------------------
            # STEP 4: Demean motion + build censor file (FD > 0.3mm excluded)
            # ------------------------------------------------------------------
            ${AFNI} "
                1d_tool.py \
                    -infile ${OUT_DIR}/${PREFIX}_motion.txt \
                    -demean \
                    -write ${OUT_DIR}/${PREFIX}_motion_demean.txt

                1d_tool.py \
                    -infile ${OUT_DIR}/${PREFIX}_motion_demean.txt \
                    -show_censor_count \
                    -censor_prev_TR \
                    -censor_motion ${MOTION_THRESH} \
                    ${OUT_DIR}/${PREFIX}_motion
            "

            # ------------------------------------------------------------------
            # STEP 5: First-level GLM with parametric modulation (3dDeconvolve)
            #
            # -stim_times_AM2 builds TWO regressors from the AM timing file:
            #   rating[0] = mean response   (constant amplitude)
            #   rating[1] = parametric slope (mean-centered rating; the result)
            # BLOCK(${BLOCK_DUR},1) = fixed boxcar HRF, peak normalized to 1 so
            # betas are in percent signal change.
            # ------------------------------------------------------------------
            ${AFNI} "
                cd ${OUT_DIR}

                3dDeconvolve \
                    -input ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz \
                    -polort A \
                    -TR_times ${TR} \
                    -censor ${OUT_DIR}/${PREFIX}_motion_censor.1D \
                    \
                    -num_stimts 7 \
                    \
                    -stim_file 1 ${OUT_DIR}/${PREFIX}_motion_demean.txt'[0]' -stim_base 1 -stim_label 1 trans_x \
                    -stim_file 2 ${OUT_DIR}/${PREFIX}_motion_demean.txt'[1]' -stim_base 2 -stim_label 2 trans_y \
                    -stim_file 3 ${OUT_DIR}/${PREFIX}_motion_demean.txt'[2]' -stim_base 3 -stim_label 3 trans_z \
                    -stim_file 4 ${OUT_DIR}/${PREFIX}_motion_demean.txt'[3]' -stim_base 4 -stim_label 4 rot_x \
                    -stim_file 5 ${OUT_DIR}/${PREFIX}_motion_demean.txt'[4]' -stim_base 5 -stim_label 5 rot_y \
                    -stim_file 6 ${OUT_DIR}/${PREFIX}_motion_demean.txt'[5]' -stim_base 6 -stim_label 6 rot_z \
                    \
                    -stim_times_AM2 7 ${TIMING} 'BLOCK(${BLOCK_DUR},1)' \
                        -stim_label 7 rating \
                    \
                    -num_glt 2 \
                    -gltsym 'SYM: rating[0]' -glt_label 1 mean_response \
                    -gltsym 'SYM: rating[1]' -glt_label 2 rating_slope \
                    \
                    -tout -fout \
                    -x1D    ${OUT_DIR}/${PREFIX}_xmat.1D \
                    -xjpeg  ${OUT_DIR}/${PREFIX}_xmat.jpg \
                    -errts  ${OUT_DIR}/${PREFIX}_errts \
                    -bucket ${OUT_DIR}/${PREFIX}_stats \
                    -jobs 8
            "

            # ------------------------------------------------------------------
            # STEP 6: 3dREMLfit (accurate autocorrelation; use for group level)
            # ------------------------------------------------------------------
            ${AFNI} "
                3dREMLfit \
                    -matrix ${OUT_DIR}/${PREFIX}_xmat.1D \
                    -input  ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz \
                    -fout -tout \
                    -Rbuck  ${OUT_DIR}/${PREFIX}_stats_REML \
                    -Rvar   ${OUT_DIR}/${PREFIX}_REMLvar \
                    -Rerrts ${OUT_DIR}/${PREFIX}_errts_REML
            "

            echo "Done: ${PARTICIPANT_ID} ${SES} task-${TASK}"
            echo "  Parametric slope sub-brick: rating_slope (in ${PREFIX}_stats_REML)"

        done
    done
done

echo "All first-level jobs complete."
