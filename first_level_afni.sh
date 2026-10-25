#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   8 cores / 32 GB / 12 h : one job loops over every participant, session and
#                            task in turn (3dDeconvolve uses the 8 cores).
#   --chdir                : run inside scratch, so any stray relative-path file
#                            lands there instead of in the home directory.
# ==============================================================================
#SBATCH --job-name=spirituality_firstlvl
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# ==============================================================================
# first_level_afni.sh — FIRST-LEVEL GLM WITH PARAMETRIC MODULATION BY RATING
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   On every trial the participant rates, with a button box, how strongly they
#   "felt the Spirit" (1 = lowest, 4 = highest). For each participant, session
#   and task, this script fits one general linear model (GLM) that estimates:
#     - mean_response : the average BOLD response to the task's trials
#     - rating_slope  : how much the BOLD response goes up per rating point
#                       (the main result of the study)
#
# WHERE IT SITS IN THE PIPELINE
#   Before: fmriprep_tedana.sh  (preprocessing, multi-echo denoising, MNI warp)
#   After : second_level_afni.sh, reproducibility_afni.sh, first_level_IM_afni.sh
#
# STEPS, PER RUN (one run = one participant x session x task)
#   1. Smooth the tedana-denoised BOLD (4 mm FWHM) and label it as MNI space
#   2. Build the analysis mask (fMRIPrep brain mask AND voxels with data at
#      every time point), then convert to percent signal change inside it
#   3. Pull the 24 motion regressors and framewise displacement out of
#      fMRIPrep's confounds table
#   4. Build the motion censor file (which volumes to drop, FD > 0.3 mm)
#   5. Build the GLM design matrix (3dDeconvolve)
#   6. Fit the GLM with autocorrelation-corrected least squares (3dREMLfit)
#
# INPUTS (per run)
#   tedana BOLD   : <TEDANA_OUT>/sub-NN/ses-N/func/
#                     sub-NN_ses-N_task-T_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
#   confounds     : <FMRIPREP_OUT>/sub-NN/ses-N/func/sub-NN_ses-N_task-T_desc-confounds_timeseries.tsv
#   brain mask    : <FMRIPREP_OUT>/sub-NN/ses-N/func/
#                     sub-NN_ses-N_task-T_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz
#   rated trials  : <TIMING_DIR>/byv-<BYV>_ses-N_task-T_ratingAM.1D   ("onset*rating")
#   unrated trials: <TIMING_DIR>/byv-<BYV>_ses-N_task-T_unrated.1D    (only if any)
#   (Timing files are made on Braden's Mac by generate_timing_files.py.)
#
# OUTPUTS (per run), in <AFNI_OUT>/sub-NN/ses-N/task-T/
#   <PREFIX>_stats_REML+tlrc   THE RESULT: betas + t/F stats. Sub-bricks by label:
#                                mean_response#0_Coef / _Tstat
#                                rating_slope#0_Coef  / _Tstat
#   <PREFIX>_errts_REML+tlrc   residuals (used for smoothness estimates later)
#   <PREFIX>_blurred_scaled    the smoothed, scaled data that was fit
#   <PREFIX>_analysis_mask     the voxels the GLM was fit in (STEP 2)
#   <PREFIX>_motion_demean.txt, <PREFIX>_fd.1D, <PREFIX>_motion_censor.1D,
#   <PREFIX>_xmat.1D
#   (<PREFIX> = sub-NN_ses-N_task-T)
#
# HOW TO RUN
#   sbatch first_level_afni.sh                 # every participant in TEDANA_OUT
#   sbatch first_level_afni.sh sub-04 sub-05   # only these participants
#   Re-running overwrites that run's outputs in place.
#
# METHODS SUMMARY (the exact values used below)
#   Smoothing      : 4 mm FWHM Gaussian (3dmerge), after tedana, before scaling
#   Scaling        : percent of each voxel's temporal mean, capped at 200
#   Mask           : fMRIPrep's run-level brain mask (MNI152NLin2009cAsym),
#                    intersected with an extents mask (voxels whose smoothed
#                    data are > 0 at every time point, as in afni_proc.py)
#   Stimulus model : fixed-duration boxcar, BLOCK(d,1) (peak = 1, so betas are
#                    in percent signal change); d = 19 s scripture, 44 s FHS,
#                    14 s architecture
#   Rating effect  : -stim_times_AM2, which gives a constant-amplitude regressor
#                    (mean_response) and a mean-centred rating regressor
#                    (rating_slope)
#   Unrated trials : a separate nuisance regressor, same boxcar (left out of a
#                    run's model if all its trials fall on censored volumes)
#   Drift          : Legendre polynomials, order 1 + int(run length / 150 s)
#                    (-polort A)
#   Motion         : 24 regressors (6 rigid-body parameters, their derivatives,
#                    and the squares of both), demeaned
#   Censoring      : a volume AND the one before it are dropped when
#                    framewise displacement (FD; Power et al. 2012, as computed
#                    by fMRIPrep) is > 0.3 mm
#   Estimation     : 3dREMLfit (ARMA(1,1) noise model, voxelwise)
# ==============================================================================

set -e    # stop at the first command that fails


# ==============================================================================
# SECTION 1. SETTINGS
# ------------------------------------------------------------------------------
# Every path and parameter this script uses is defined here. Settings marked
# [shared] also appear in other scripts and must match them; after editing one,
# run  bash check_settings.sh  to confirm every copy still agrees.
# ==============================================================================

# --- 1a. Where the data live --------------------------------------------------
# Three storage areas on the BYU supercomputer:
#   /home/bradenf4                     backed up, never purged: code, timing, logs
#   /nobackup/autodelete/usr/bradenf4  scratch, NOT backed up, purged after 12
#                                      weeks unused: everything this job reads/writes
#   /nobackup/archive/usr/bradenf4     raw data (read only) and archived results
TEDANA_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/tedana          # [shared] input BOLD
FMRIPREP_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/fmriprep      # [shared] confounds + masks
AFNI_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_firstlvl     # [shared] output
TIMING_DIR=/home/bradenf4/spirituality_fmri/timing                                     # [shared] timing files
PARTICIPANTS_TSV=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata/participants.tsv  # [shared]
EXCLUSIONS_TSV=/home/bradenf4/spirituality_fmri/scripts/exclusions.tsv                 # [shared] runs to skip
LOGS=/home/bradenf4/spirituality_fmri/logs                                             # [shared]

# --- 1b. Software -------------------------------------------------------------
# AFNI runs inside an Apptainer image. The image only sees folders that are
# explicitly "bound" into it, so all three storage areas are bound.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]

# Ignore any Python packages in ~/.local. A stray numpy there once broke tedana;
# this script's python3 call (STEP 3) only needs the standard library.
export PYTHONNOUSERSITE=1

# --- 1c. Analysis parameters --------------------------------------------------
TR=1.388              # [shared] seconds between volumes
SMOOTH_FWHM=4.0       # smoothing kernel, mm FWHM
FD_THRESH=0.3         # censor volumes with framewise displacement > 0.3 mm (STEP 4)

# Percent signal change ceiling (STEP 2). 200 is the AFNI default
# (afni_proc.py -scale_max_val): a voxel at 200 is 100% above its own mean, far
# outside real BOLD changes (~0.1-5%), so values that high are edge, ventricle or
# near-zero-mean junk.
PSC_CAP=200

# Tasks, by their BIDS labels.
TASKS=(scripture FHS architecture)                                                     # [shared]

# Fixed stimulus duration per task, in seconds. These are the median on-screen
# times measured by generate_timing_files.py:
#   scripture 19.3 s | FHS 44.3 s | architecture 14.2 s
block_dur_for() {                                                                      # [shared]
    case "$1" in
        scripture)    echo 19 ;;
        FHS)          echo 44 ;;
        architecture) echo 14 ;;
        *)            echo "ERROR: unknown task $1" >&2; exit 1 ;;
    esac
}


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# byv_for_sub sub-NN  ->  prints that participant's BYV id (e.g. 0AK1RT).
# BIDS labels participants sub-NN (numbered in scan-date order), but the timing
# files are named by the lab's BYV id. participants.tsv maps one to the other via
# its 'participant_id' and 'byv_id' columns.
# (The awk variable is 'want', not 'sub': sub() is a built-in awk function.)
byv_for_sub() {
    awk -F'\t' -v want="$1" '
        NR==1 { for (i=1;i<=NF;i++){ if($i=="participant_id")p=i; if($i=="byv_id")b=i}; next }
        $p==want { print $b; exit }
    ' "${PARTICIPANTS_TSV}"
}

# xmat_signal_kept CENSORED.xmat.1D UNCENSORED.xmat.1D
#   ->  one line per design-matrix column:  label <TAB> signal_kept <TAB> in_run
# Both files are 3dDeconvolve design matrices for the same model: the first
# with censored volumes removed (-x1D), the second with every volume kept
# (-x1D_uncensored). For each column, signal_kept = its sum of squares in the
# censored matrix divided by its sum of squares in the uncensored one, i.e. the
# fraction of that regressor that survives censoring (0 if the column is empty
# even before censoring). in_run = "no" when the column is empty before
# censoring (its response lies past the end of the run), else "yes".
# Column names come from the matrix's "ColumnLabels" header line.
xmat_signal_kept() {                                                                   # [shared]
    awk '
        FNR == 1 { f++ }
        /ColumnLabels/ {
            s = $0; sub(/.*ColumnLabels = "/, "", s); sub(/".*/, "", s)
            n = split(s, lab, ";")
            for (i = 1; i <= n; i++) gsub(/^ +| +$/, "", lab[i])
            next
        }
        /^#/ || NF == 0 { next }
        { for (i = 1; i <= NF; i++) ss[f, i] += $i * $i }
        END {
            for (i = 1; i <= n; i++)
                printf "%s\t%.6f\t%s\n", lab[i], (ss[2, i] > 0 ? ss[1, i] / ss[2, i] : 0), (ss[2, i] > 0 ? "yes" : "no")
        }
    ' "$1" "$2"
}

# is_excluded sub-NN ses-N task  ->  succeeds (and prints the reason) if that run
# is listed in exclusions.tsv. That table has one row per excluded run:
#   byv_id  session  task  reason  date_added      ("*" = every session / task)
# Currently: 1VK4L2 ses-1, all tasks (asleep for most of the session).
is_excluded() {                                                                        # [shared]
    local byv
    [ -f "${EXCLUSIONS_TSV}" ] || return 1
    case "$1" in sub-*) byv=$(byv_for_sub "$1") ;; *) byv=$1 ;; esac
    [ -n "${byv}" ] || return 1
    awk -F'\t' -v b="${byv}" -v s="$2" -v t="$3" '
        NR==1 || /^#/ { next }
        $1==b && ($2==s || $2=="*") && ($3==t || $3=="*") { print $4; found=1; exit }
        END { exit(found ? 0 : 1) }' "${EXCLUSIONS_TSV}"
}


# ==============================================================================
# SECTION 3. SETUP AND CHECKS
# ==============================================================================

mkdir -p "${LOGS}" "${AFNI_OUT}"

module load "${APPTAINER_MODULE}"

# AFNI refuses to overwrite existing files by default, which (with set -e) would
# abort any re-run. Apptainer passes APPTAINERENV_* variables into the image, so
# this one setting reaches every AFNI command below.
export APPTAINERENV_AFNI_DECONFLICT=OVERWRITE

# Every AFNI command below runs as:  ${AFNI} "<commands>"
AFNI="apptainer exec ${BIND_ARGS} ${AFNI_SIF} bash -c"

if [ ! -f "${PARTICIPANTS_TSV}" ]; then
    echo "ERROR: participants.tsv not found at ${PARTICIPANTS_TSV}"
    echo "       It maps BIDS sub-NN to the BYV id used in timing filenames."
    exit 1
fi


# ==============================================================================
# SECTION 4. MAIN LOOP: participants -> sessions -> tasks
# ==============================================================================

# Participants: the ones named on the command line, or everyone in TEDANA_OUT.
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

            # --- Skip runs listed in exclusions.tsv ---------------------------
            if why=$(is_excluded "${PARTICIPANT_ID}" "${SES}" "${TASK}"); then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: EXCLUDED (exclusions.tsv) -- ${why}"
                continue
            fi

            # --- This run's input files ---------------------------------------
            BOLD_IN=${TEDANA_FUNC}/${PARTICIPANT_ID}_${SES}_task-${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
            CONFOUNDS=${FMRIPREP_OUT}/${PARTICIPANT_ID}/${SES}/func/${PARTICIPANT_ID}_${SES}_task-${TASK}_desc-confounds_timeseries.tsv
            TIMING=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_ratingAM.1D
            UNRATED=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_unrated.1D
            MASK=${FMRIPREP_OUT}/${PARTICIPANT_ID}/${SES}/func/${PARTICIPANT_ID}_${SES}_task-${TASK}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz

            # A run with no usable ratings has no timing file (e.g. 0AK1RT ses-1
            # FHS) and is skipped here, as is any run missing its BOLD or mask.
            if [ ! -f "${BOLD_IN}" ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: BOLD not found"
                continue
            fi
            if [ ! -f "${TIMING}" ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: no timing file (dead/absent run) ${TIMING}"
                continue
            fi
            if [ ! -f "${MASK}" ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: no fMRIPrep brain mask ${MASK}"
                continue
            fi

            OUT_DIR=${AFNI_OUT}/${PARTICIPANT_ID}/${SES}/task-${TASK}
            mkdir -p "${OUT_DIR}"
            PREFIX=${PARTICIPANT_ID}_${SES}_task-${TASK}
            BLOCK_DUR=$(block_dur_for "${TASK}")

            echo "=== ${PARTICIPANT_ID} ${SES} task-${TASK} (byv-${BYV}, BLOCK ${BLOCK_DUR}s) ==="

            # ------------------------------------------------------------------
            # STEP 1. Smooth (4 mm FWHM), then label the data as MNI space.
            #
            # The data really are in MNI152NLin2009cAsym, but tedana writes the
            # NIfTI header with a space code meaning "scanner space", so AFNI
            # would read them as native (+orig). AFNI copies a dataset's space
            # label to every output made from it, so fixing it once here makes
            # every later file correctly MNI (+tlrc). The fix is applied to our
            # smoothed copy, not to tedana's file, which is left untouched.
            # ------------------------------------------------------------------
            ${AFNI} "
                3dmerge \
                    -1blur_fwhm ${SMOOTH_FWHM} \
                    -doall \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred.nii.gz \
                    ${BOLD_IN}

                3drefit -space MNI ${OUT_DIR}/${PREFIX}_blurred.nii.gz
            "

            # ------------------------------------------------------------------
            # STEP 2. Build the analysis mask, then convert to percent signal
            #         change inside it.
            #
            # ANALYSIS MASK = fMRIPrep's brain mask AND the extents mask.
            #   Extents mask: voxels whose smoothed data are above 0 at EVERY
            #   time point (3dTstat -min > 0). Near the edge of the data, a
            #   voxel can have valid values in most volumes but 0 (or a small
            #   negative value, from the interpolation when tedana's output was
            #   warped to MNI) in a few. Those time points would become 0 below,
            #   and if their timing happens to line up with a regressor they can
            #   look like activation, often in a whole cluster of edge voxels
            #   (afni_proc.py -help, MASKING NOTE). AFNI's own pipeline removes
            #   such voxels by default with the same kind of mask: "every voxel
            #   in the extents mask has data at every time point."
            #   Saved as <PREFIX>_analysis_mask.nii.gz; first_level_IM_afni.sh
            #   uses the same file, so both models fit exactly the same voxels.
            #
            # SCALING:
            #   a = smoothed time series, b = its mean over time,
            #   c = the analysis mask
            #   a/b*100          each voxel's mean becomes 100, so values read as
            #                    percent of that voxel's mean
            #   min(PSC_CAP, ..) clip non-physiological values (see PSC_CAP)
            #   step(a)*step(b)  zero any time point that is <= 0 (inside the
            #                    analysis mask there are none, by construction)
            #   c*               zero everything outside the analysis mask, so
            #                    the GLM only fits those voxels
            # Order (smooth, then mask, then scale) matches afni_proc.py.
            # ------------------------------------------------------------------
            ANALYSIS_MASK=${OUT_DIR}/${PREFIX}_analysis_mask.nii.gz
            ${AFNI} "
                3dTstat \
                    -mean \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred_mean.nii.gz \
                    ${OUT_DIR}/${PREFIX}_blurred.nii.gz

                3dTstat \
                    -min \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred_min.nii.gz \
                    ${OUT_DIR}/${PREFIX}_blurred.nii.gz

                3dcalc \
                    -a ${MASK} \
                    -b ${OUT_DIR}/${PREFIX}_blurred_min.nii.gz \
                    -expr 'step(a)*step(b)' \
                    -datum byte \
                    -prefix ${ANALYSIS_MASK}

                3dcalc \
                    -a ${OUT_DIR}/${PREFIX}_blurred.nii.gz \
                    -b ${OUT_DIR}/${PREFIX}_blurred_mean.nii.gz \
                    -c ${ANALYSIS_MASK} \
                    -expr 'c*min(${PSC_CAP}, a/b*100)*step(a)*step(b)' \
                    -prefix ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz
            "
            NMASK=$(${AFNI} "3dBrickStat -count -non-zero ${MASK}" 2>/dev/null | tail -n1 | awk '{print $1}')
            NANAL=$(${AFNI} "3dBrickStat -count -non-zero ${ANALYSIS_MASK}" 2>/dev/null | tail -n1 | awk '{print $1}')
            echo "  analysis mask: ${NANAL} voxels (fMRIPrep brain mask ${NMASK}; extents mask removed $((NMASK - NANAL)))"

            # ------------------------------------------------------------------
            # STEP 3. Extract motion information from fMRIPrep's confounds table.
            #
            # (a) <PREFIX>_motion.txt: the 24 motion regressors. For each of the
            #     6 rigid-body parameters (translation x/y/z in mm, rotation
            #     x/y/z in radians): the parameter, its derivative, and the
            #     squares of both (the "Friston-24" set).
            # (b) <PREFIX>_fd.1D: framewise displacement, one value per volume.
            #     fMRIPrep computes it as in Power et al. (2012): the sum of the
            #     absolute volume-to-volume changes in the 6 parameters, with
            #     rotations converted to mm of movement on a 50 mm sphere.
            # fMRIPrep writes "n/a" in the first row of every derivative and of
            # FD (there is no previous volume); those become 0.
            # ------------------------------------------------------------------
            python3 -c "
import csv
with open('${CONFOUNDS}') as f:
    reader = csv.DictReader(f, delimiter='\t')
    cols = ['trans_x','trans_x_derivative1','trans_x_derivative1_power2','trans_x_power2',
            'trans_y','trans_y_derivative1','trans_y_derivative1_power2','trans_y_power2',
            'trans_z','trans_z_derivative1','trans_z_derivative1_power2','trans_z_power2',
            'rot_x','rot_x_derivative1','rot_x_power2','rot_x_derivative1_power2',
            'rot_y','rot_y_derivative1','rot_y_power2','rot_y_derivative1_power2',
            'rot_z','rot_z_derivative1','rot_z_power2','rot_z_derivative1_power2']
    with open('${OUT_DIR}/${PREFIX}_motion.txt', 'w') as mot, \
         open('${OUT_DIR}/${PREFIX}_fd.1D', 'w') as fd:
        for row in reader:
            vals = ['0' if row[c] in ('n/a','') else row[c] for c in cols]
            print(' '.join(vals), file=mot)
            v = row['framewise_displacement']
            print('0' if v in ('n/a','') else v, file=fd)
"

            # ------------------------------------------------------------------
            # STEP 4. Demean the motion regressors and build the censor file.
            #
            # Demeaning: each of the 24 motion columns is shifted to mean 0
            # (1d_tool.py -demean). -overwrite: 1d_tool.py does not obey
            # AFNI_DECONFLICT (SECTION 3), so without it a re-run would stop here
            # because the file already exists.
            #
            # Censoring: a volume is dropped when its FD is > FD_THRESH (0.3 mm).
            # FD at volume t measures the movement between volumes t-1 and t,
            # so volume t-1 is dropped too (the same rule as afni_proc.py's
            # -censor_prev_TR). Output: <PREFIX>_motion_censor.1D, one value per
            # volume, 1 = keep, 0 = drop.
            # ------------------------------------------------------------------
            ${AFNI} "
                1d_tool.py \
                    -infile ${OUT_DIR}/${PREFIX}_motion.txt \
                    -demean \
                    -overwrite \
                    -write ${OUT_DIR}/${PREFIX}_motion_demean.txt
            "

            awk -v thr="${FD_THRESH}" '
                { fd[NR] = $1 + 0 }
                END {
                    for (t = 1; t <= NR; t++) {
                        drop = (fd[t] > thr) || (t < NR && fd[t + 1] > thr)
                        print (drop ? 0 : 1)
                    }
                }' "${OUT_DIR}/${PREFIX}_fd.1D" > "${OUT_DIR}/${PREFIX}_motion_censor.1D"
            NCENS=$(grep -c '^0$' "${OUT_DIR}/${PREFIX}_motion_censor.1D" || true)
            NVOL=$(wc -l < "${OUT_DIR}/${PREFIX}_motion_censor.1D" | tr -d ' ')
            echo "  censored ${NCENS} of ${NVOL} volumes (FD > ${FD_THRESH} mm, plus the volume before each)"

            # ------------------------------------------------------------------
            # STEP 5. Build the GLM design matrix (3dDeconvolve).
            #
            # Regressors:
            #   rating (-stim_times_AM2): from "onset*rating" pairs, AFNI makes
            #     two regressors with a fixed boxcar of the task's duration:
            #       rating[0] = same height every trial  -> mean_response
            #       rating[1] = height = rating minus the run's mean rating
            #                   -> rating_slope (BOLD change per rating point)
            #     BLOCK(d,1) scales the boxcar response to peak at 1, so the
            #     betas are in percent signal change.
            #   unrated (only if the run has unrated trials): trials with no
            #     button press. The participant still saw them, so they get
            #     their own nuisance regressor; otherwise their response would
            #     fall into the baseline and pull mean_response toward 0.
            #     (Added 2026-10-20.) Their betas are not a result.
            #   motion (-ortvec): the 24 demeaned motion regressors.
            #   drift (-polort A): polynomials, order 1 + int(run length/150 s).
            # Volumes marked 0 in the censor file are left out (-censor).
            # The two -gltsym lines label rating[0] and rating[1] as
            # mean_response and rating_slope in the output.
            # ------------------------------------------------------------------
            NSTIM=1
            UNRATED_OPTS=""
            if [ -f "${UNRATED}" ]; then
                NSTIM=2
                UNRATED_OPTS="-stim_times 2 ${UNRATED} 'BLOCK(${BLOCK_DUR},1)' -stim_label 2 unrated"
                echo "  + unrated nuisance regressor: $(wc -w < "${UNRATED}" | tr -d ' ') trial(s)"
            fi

            # --- Check that every stimulus regressor still has data ------------
            # A regressor whose whole response falls on censored volumes is all
            # zeros, and 3dREMLfit (STEP 6) stops with "Cannot continue with all
            # zero column", which would end the job for every remaining run. The
            # unrated regressor is the one at risk: many runs have only 1-4
            # unrated trials, so one burst of movement can censor all of them.
            # A quick design-only 3dDeconvolve (-nodata: no data are read; same
            # number of volumes, TR, drift and censoring as the real model)
            # writes the design with and without censoring, and xmat_signal_kept
            # gives the fraction of each regressor left after censoring.
            #   unrated has nothing left      -> drop the unrated regressor
            #   rating[0] or rating[1] empty  -> skip this run with a warning
            #     (all rated trials censored, or every trial got the same
            #     rating so there is no slope to estimate)
            NVOL=$(wc -l < "${OUT_DIR}/${PREFIX}_motion_censor.1D" | tr -d ' ')
            DATA_TR=$(${AFNI} "3dinfo -tr ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz" 2>/dev/null | tail -n1)
            ${AFNI} "
                cd ${OUT_DIR}
                3dDeconvolve \
                    -nodata ${NVOL} ${DATA_TR} \
                    -polort A \
                    -censor ${OUT_DIR}/${PREFIX}_motion_censor.1D \
                    -num_stimts ${NSTIM} \
                    -stim_times_AM2 1 ${TIMING} 'BLOCK(${BLOCK_DUR},1)' -stim_label 1 rating \
                    ${UNRATED_OPTS} \
                    -x1D ${OUT_DIR}/${PREFIX}_designcheck.1D \
                    -x1D_uncensored ${OUT_DIR}/${PREFIX}_designcheck_uncensored.1D \
                    -x1D_stop
            " > /dev/null 2>&1
            KEPT_TABLE=$(xmat_signal_kept \
                "${OUT_DIR}/${PREFIX}_designcheck.1D" \
                "${OUT_DIR}/${PREFIX}_designcheck_uncensored.1D")
            kept_of() { printf '%s\n' "${KEPT_TABLE}" | awk -F'\t' -v c="$1" '$1 == c { print $2 }'; }
            is_empty() { awk -v v="${1:-0}" 'BEGIN { exit !(v + 0 <= 0) }'; }
            if is_empty "$(kept_of 'rating#0')" || is_empty "$(kept_of 'rating#1')"; then
                echo "  WARNING: skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: a rating regressor is empty (all rated trials censored, or every trial has the same rating)"
                echo "           (rating#0 kept $(kept_of 'rating#0'), rating#1 kept $(kept_of 'rating#1'); see ${PREFIX}_designcheck.1D)"
                continue
            fi
            if [ "${NSTIM}" = 2 ] && is_empty "$(kept_of 'unrated#0')"; then
                echo "  WARNING: unrated regressor dropped: all its trials are censored"
                NSTIM=1
                UNRATED_OPTS=""
            fi
            ${AFNI} "
                cd ${OUT_DIR}

                3dDeconvolve \
                    -input ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz \
                    -mask ${ANALYSIS_MASK} \
                    -polort A \
                    -TR_times ${TR} \
                    -censor ${OUT_DIR}/${PREFIX}_motion_censor.1D \
                    \
                    -num_stimts ${NSTIM} \
                    \
                    -ortvec ${OUT_DIR}/${PREFIX}_motion_demean.txt motion \
                    \
                    -stim_times_AM2 1 ${TIMING} 'BLOCK(${BLOCK_DUR},1)' \
                        -stim_label 1 rating \
                    ${UNRATED_OPTS} \
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
            # STEP 6. Fit the GLM with 3dREMLfit.
            #
            # Same design matrix, but the noise in each voxel is modelled as
            # temporally autocorrelated (ARMA(1,1)), which gives more accurate
            # t-statistics than ordinary least squares. The REML outputs are the
            # ones every later script uses.
            # ------------------------------------------------------------------
            ${AFNI} "
                3dREMLfit \
                    -matrix ${OUT_DIR}/${PREFIX}_xmat.1D \
                    -input  ${OUT_DIR}/${PREFIX}_blurred_scaled.nii.gz \
                    -mask   ${ANALYSIS_MASK} \
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
