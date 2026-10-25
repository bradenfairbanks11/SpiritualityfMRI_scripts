#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   8 cores / 32 GB / 12 h : one job loops over every participant, session and
#                            task in turn. (3dREMLfit, the slowest step, can
#                            split its work across cores with OpenMP.)
#   --chdir                : run inside scratch, so any stray relative-path file
#                            lands there instead of in the home directory.
# ==============================================================================
#SBATCH --job-name=spirituality_firstlvl_IM
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
# first_level_IM_afni.sh — SINGLE-TRIAL FIRST-LEVEL GLM (one beta per trial)
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   first_level_afni.sh estimates ONE average response per run. This script
#   uses the same data and the same nuisance model, but gives every rated trial
#   its OWN regressor (AFNI's -stim_times_IM; IM = "individually modulated").
#   The result is one beta map PER TRIAL. Each trial is also labelled "low"
#   (rating 1-2) or "high" (rating 3-4).
#   These trial-level betas are the input to AFNI's trial-level test-retest
#   reliability tools, 3dLMEr -TRR (whole brain) and TRR (region level):
#     Chen et al. (2021) NeuroImage 245:118647, "Trial and error: A hierarchical
#     modeling approach to test-retest reliability".
#
#   It does NOT redo preprocessing. It reuses first_level_afni.sh's files from
#   the same run folder, so smoothing, scaling, motion regressors and motion
#   censoring are identical to the main analysis by construction.
#
# WHERE IT SITS IN THE PIPELINE
#   Before: first_level_afni.sh  (must have run for the same run first)
#   After : trr_export_afni.sh   (splits the trial betas into one file per
#                                 trial and builds the TRR data tables)
#   QC    : check_trial_censoring.sh  (how much of each trial was censored)
#
# STEPS, PER RUN (one run = one participant x session x task)
#   1. Split the rating timing file into an onset-only timing file and a
#      per-trial table (trial number, onset, rating, low/high)
#   2. Build the design matrix with one regressor per trial (3dDeconvolve)
#   3. Drop trials with no data left after motion censoring (their regressor
#      would be all zeros, which 3dREMLfit cannot fit)
#   4. Rebuild the design matrix without them, if any were dropped
#   5. Fit it with autocorrelation-corrected least squares (3dREMLfit) and keep
#      only the trial betas
#   6. Check that the output has exactly one sub-brick per kept trial, in order
#
# INPUTS (per run)
#   from first_level_afni.sh, in <AFNI_OUT>/sub-NN/ses-N/task-T/:
#     <PREFIX>_blurred_scaled.nii.gz   smoothed, scaled BOLD          (its STEP 2)
#     <PREFIX>_motion_demean.txt       24 demeaned motion regressors  (its STEP 4)
#     <PREFIX>_motion_censor.1D        motion censor file             (its STEP 4)
#     <PREFIX>_analysis_mask.nii.gz    the voxels its GLM was fit in  (its STEP 2):
#                                      fMRIPrep brain mask AND voxels with data
#                                      at every time point
#   rated trials  : <TIMING_DIR>/byv-<BYV>_ses-N_task-T_ratingAM.1D   ("onset*rating")
#   unrated trials: <TIMING_DIR>/byv-<BYV>_ses-N_task-T_unrated.1D    (only if any)
#   (<PREFIX> = sub-NN_ses-N_task-T.) A run missing any of the first five files
#   is skipped with a message.
#
#   Which trials get a beta: exactly the trials in the *_ratingAM.1D file, in
#   that file's order, so the trial set is the same as in the main (AM) model.
#   Trials with no button press have no rating, so they are not in that file
#   and get no trial beta (there is nothing to call low or high). They are still
#   modelled: if the *_unrated.1D file exists, they enter as one nuisance
#   regressor, exactly as in first_level_afni.sh, so their response does not
#   fall into the baseline. That regressor's beta is dropped from the output
#   (STEP 3).
#
# OUTPUTS (per run), in <AFNI_OUT>/sub-NN/ses-N/task-T/single_trial/
#   <PREFIX>_IM_betas+tlrc     THE RESULT: one sub-brick per rated trial, labelled
#                                trial#0, trial#1, ... (REML betas, in percent
#                                signal change). Sub-brick k = trial k of the table.
#   <PREFIX>_IM_trials.tsv     one row per KEPT trial. Columns:
#                                trial  (0-based; = sub-brick index)
#                                onset  (s)   rating  (1-4)
#                                cond   (low = rating 1-2, high = rating 3-4)
#                                orig_trial  (its number among ALL rated trials)
#                                signal_kept (fraction of its response left
#                                             after censoring, 0-1)
#   <PREFIX>_IM_trial_censoring.tsv  every rated trial, with signal_kept and
#                                status (kept / dropped, and why)
#   <PREFIX>_IM_alltrials.tsv, <PREFIX>_IM_alltimes.1D  every rated trial
#   <PREFIX>_IM_times.1D       onsets of the kept trials (given to -stim_times_IM)
#   <PREFIX>_IM_xmat.1D, <PREFIX>_IM_xmat.jpg   the final design matrix and a
#                                picture of it (_IM_xmat_alltrials* = before
#                                dropping trials)
#
# HOW TO RUN
#   sbatch first_level_IM_afni.sh                 # every participant in AFNI_OUT
#   sbatch first_level_IM_afni.sh sub-04 sub-05   # only these participants
#   Re-running overwrites that run's outputs in place.
#
# METHODS SUMMARY (the exact values used below)
#   Data           : the smoothed (4 mm FWHM), percent-signal-change-scaled
#                    data from first_level_afni.sh, fit inside the same
#                    analysis mask
#   Stimulus model : one regressor per rated trial (-stim_times_IM), each a
#                    fixed-duration boxcar, BLOCK(d,1) (peak = 1, so betas are in
#                    percent signal change); d = 19 s scripture, 44 s FHS,
#                    14 s architecture
#   Unrated trials : one separate nuisance regressor, same boxcar
#   Drift          : Legendre polynomials, order 1 + int(run length / 150 s)
#                    (-polort A)
#   Motion         : the same 24 demeaned motion regressors as the main model
#   Censoring      : the same censor file as the main model (a volume AND the
#                    one before it are dropped when fMRIPrep's framewise
#                    displacement is > 0.3 mm)
#   Estimation     : 3dREMLfit (ARMA(1,1) noise model, voxelwise)
#   Rating split   : high = rating >= 3, low = rating 1-2 (for 3dLMEr -TRR,
#                    which needs exactly 2 conditions)
#   Dropped trials : a trial whose whole modelled response falls on censored
#                    volumes (or whose onset is after the end of the run) gets
#                    no beta and is left out of the model
#   Runs with fewer than 2 (kept) rated trials are skipped.
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

# --- 1c. Analysis parameters --------------------------------------------------
# TR is passed to 3dDeconvolve as -TR_times (STEP 2). Per AFNI's help, that only
# sets the time step of -iresp/-sresp outputs, which this script does not make;
# the TR used to build the regressors is read from the input dataset's header.
TR=1.388              # [shared] seconds between volumes

# Rating split for the 2-condition 3dLMEr -TRR model (STEP 1):
# a trial is "high" if its 1-4 rating is >= HIGH_MIN, otherwise "low".
HIGH_MIN=3

# Trials with no data left after motion censoring (STEP 3). signal_kept = the
# fraction of a trial's modelled response that falls on volumes that survive
# censoring (1 = untouched, 0 = all censored). A trial is dropped unless its
# signal_kept is ABOVE this value. 0 = drop only trials with nothing left, which
# is the minimum needed for 3dREMLfit to run. (Values are rounded to 6 decimals,
# so a trial with less than 0.0000005 of its response left also counts as 0.)
# check_trial_censoring.sh and <PREFIX>_IM_trial_censoring.tsv show how much of
# each trial was censored, if you want a stricter cut-off.
MIN_SIGNAL_KEPT=0                                                                      # [shared]

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
# its 'participant_id' and 'byv_id' columns. (Same lookup as first_level_afni.sh.)
byv_for_sub() {
    local sub="$1"
    awk -F'\t' -v want="${sub}" '
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

mkdir -p "${LOGS}"

module load "${APPTAINER_MODULE}"

# AFNI refuses to overwrite existing files by default, which (with set -e) would
# abort any re-run. Apptainer passes APPTAINERENV_* variables into the image, so
# this one setting reaches every AFNI command below.
export APPTAINERENV_AFNI_DECONFLICT=OVERWRITE

# Every AFNI command below runs as:  ${AFNI} "<commands>"
AFNI="apptainer exec ${BIND_ARGS} ${AFNI_SIF} bash -c"

if [ ! -f "${PARTICIPANTS_TSV}" ]; then
    echo "ERROR: participants.tsv not found at ${PARTICIPANTS_TSV}"
    exit 1
fi


# ==============================================================================
# SECTION 4. MAIN LOOP: participants -> sessions -> tasks
# ==============================================================================

# Participants: the ones named on the command line, or everyone in AFNI_OUT
# (i.e. everyone first_level_afni.sh has produced output for).
if [ "$#" -gt 0 ]; then
    PARTICIPANT_DIRS=()
    for arg in "$@"; do
        PARTICIPANT_DIRS+=("${AFNI_OUT}/sub-${arg#sub-}")
    done
else
    echo "WARNING: no participant labels given -> defaulting to ALL sub-* in ${AFNI_OUT}"
    PARTICIPANT_DIRS=("${AFNI_OUT}"/sub-*)
fi

for PARTICIPANT_DIR in "${PARTICIPANT_DIRS[@]}"; do
    if [ ! -d "${PARTICIPANT_DIR}" ]; then
        echo "ERROR: ${PARTICIPANT_DIR} not found in first-level output -- skipping."
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

        for TASK in "${TASKS[@]}"; do

            # --- Skip runs listed in exclusions.tsv ---------------------------
            if why=$(is_excluded "${PARTICIPANT_ID}" "${SES}" "${TASK}"); then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: EXCLUDED (exclusions.tsv) -- ${why}"
                continue
            fi

            # --- This run's input files ---------------------------------------
            RUN_DIR=${AFNI_OUT}/${PARTICIPANT_ID}/${SES}/task-${TASK}
            PREFIX=${PARTICIPANT_ID}_${SES}_task-${TASK}
            SCALED=${RUN_DIR}/${PREFIX}_blurred_scaled.nii.gz
            MOTION=${RUN_DIR}/${PREFIX}_motion_demean.txt
            CENSOR=${RUN_DIR}/${PREFIX}_motion_censor.1D
            TIMING=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_ratingAM.1D
            UNRATED=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_unrated.1D
            # The analysis mask first_level_afni.sh fit in (fMRIPrep brain mask AND
            # voxels with data at every time point), so both models use the same voxels.
            MASK=${RUN_DIR}/${PREFIX}_analysis_mask.nii.gz

            # Skip the run if any required input is missing (the unrated file is
            # optional and not checked here).
            missing=""
            for f in "${SCALED}" "${MOTION}" "${CENSOR}" "${TIMING}" "${MASK}"; do
                [ -f "${f}" ] || missing="${missing} ${f}"
            done
            if [ -n "${missing}" ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: missing${missing}"
                echo "  (run first_level_afni.sh for this run first)"
                continue
            fi

            OUT_DIR=${RUN_DIR}/single_trial
            mkdir -p "${OUT_DIR}"
            BLOCK_DUR=$(block_dur_for "${TASK}")
            ALLTIMES=${OUT_DIR}/${PREFIX}_IM_alltimes.1D        # every rated trial
            ALLTRIALS=${OUT_DIR}/${PREFIX}_IM_alltrials.tsv     # every rated trial
            CENSCHECK=${OUT_DIR}/${PREFIX}_IM_trial_censoring.tsv
            TIMES=${OUT_DIR}/${PREFIX}_IM_times.1D              # trials that get a beta
            TRIALS=${OUT_DIR}/${PREFIX}_IM_trials.tsv           # trials that get a beta

            # ------------------------------------------------------------------
            # STEP 1. Split the rating timing file into an onset-only timing
            #         file and a per-trial table (every rated trial).
            #
            # The *_ratingAM.1D file holds one "onset*rating" token per rated
            # trial. Each token is put on its own line (tr), then awk:
            #   - numbers the trials 0, 1, 2, ... in file order
            #   - labels each "high" if rating >= HIGH_MIN (3), else "low"
            #   - writes one row per trial to <PREFIX>_IM_alltrials.tsv
            #   - writes all the onsets, on one line, to <PREFIX>_IM_alltimes.1D
            # Tokens without a "*" are ignored. A run with fewer than 2 rated
            # trials is skipped.
            # ------------------------------------------------------------------
            tr ' \t' '\n\n' < "${TIMING}" | awk -v hi="${HIGH_MIN}" -v times="${ALLTIMES}" '
                BEGIN { print "trial\tonset\trating\tcond"; n = 0; line = "" }
                $0 ~ /\*/ {
                    split($0, p, "*")
                    printf "%d\t%s\t%s\t%s\n", n, p[1], p[2], (p[2] + 0 >= hi ? "high" : "low")
                    line = line (n ? " " : "") p[1]
                    n++
                }
                END { print line > times }
            ' > "${ALLTRIALS}"
            NALL=$(($(wc -l < "${ALLTRIALS}") - 1))
            if [ "${NALL}" -lt 2 ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: only ${NALL} trial(s) in ${TIMING}"
                continue
            fi

            # ------------------------------------------------------------------
            # STEP 2. Build the design matrix with one regressor per trial
            #         (3dDeconvolve), for EVERY rated trial.
            #
            # Regressors:
            #   trial (-stim_times_IM): every onset gets its own regressor, a
            #     boxcar of the task's duration. BLOCK(d,1) scales the boxcar
            #     response to peak at 1, so the betas are in percent signal
            #     change.
            #   unrated (only if the run has unrated trials): no-response trials
            #     as one nuisance regressor, same boxcar, same as
            #     first_level_afni.sh.
            #   motion (-ortvec): the 24 demeaned motion regressors.
            #   drift (-polort A): polynomials, order 1 + int(run length/150 s).
            # Same nuisance model as first_level_afni.sh.
            # Volumes marked 0 in the censor file are left out (-censor).
            # -local_times: the onsets are seconds from the start of the run
            #   (it applies to the -stim_times options written after it).
            # Two copies of the design are written: with censoring
            #   (_IM_xmat_alltrials.1D, censored volumes removed) and without
            #   (_IM_xmat_alltrials_uncensored.1D). STEP 3 compares them.
            # -x1D_stop: write the design and stop; nothing is fitted yet.
            # ------------------------------------------------------------------
            stim_opts() {   # $1 = trial onset file, $2 = 1 to include unrated
                if [ "$2" = 1 ]; then
                    echo "-num_stimts 2 -local_times -stim_times_IM 1 $1 'BLOCK(${BLOCK_DUR},1)' -stim_label 1 trial -stim_times 2 ${UNRATED} 'BLOCK(${BLOCK_DUR},1)' -stim_label 2 unrated"
                else
                    echo "-num_stimts 1 -local_times -stim_times_IM 1 $1 'BLOCK(${BLOCK_DUR},1)' -stim_label 1 trial"
                fi
            }
            USE_UNRATED=0
            if [ -f "${UNRATED}" ]; then
                USE_UNRATED=1
                echo "  + unrated nuisance regressor: $(wc -w < "${UNRATED}" | tr -d ' ') trial(s)"
            fi
            ${AFNI} "
                cd ${OUT_DIR}

                3dDeconvolve \
                    -input ${SCALED} \
                    -mask ${MASK} \
                    -polort A \
                    -TR_times ${TR} \
                    -censor ${CENSOR} \
                    -ortvec ${MOTION} motion \
                    $(stim_opts "${ALLTIMES}" "${USE_UNRATED}") \
                    -x1D ${OUT_DIR}/${PREFIX}_IM_xmat_alltrials.1D \
                    -xjpeg ${OUT_DIR}/${PREFIX}_IM_xmat_alltrials.jpg \
                    -x1D_uncensored ${OUT_DIR}/${PREFIX}_IM_xmat_alltrials_uncensored.1D \
                    -x1D_stop
            "

            # ------------------------------------------------------------------
            # STEP 3. Drop trials that have no data left after censoring.
            #
            # WHY: a trial whose whole modelled response falls on censored
            # volumes gets a regressor that is all zeros. Its beta cannot be estimated, and 3dREMLfit stops with
            # "Cannot continue with all zero column". With 'set -e' that would
            # end the job for every remaining participant.
            #
            # HOW: for each regressor, signal_kept = (sum of squares of the
            # regressor over the volumes that survive censoring) / (sum of
            # squares over all volumes), i.e. the fraction of that trial's
            # modelled response that is still in the data (1 = untouched,
            # 0 = nothing left). A trial is DROPPED when signal_kept is not
            # above MIN_SIGNAL_KEPT (see SECTION 1c; 0 = drop only trials with
            # nothing left). Dropped trials are left out of the model entirely:
            # their response volumes are already gone, so there is nothing to
            # model. The unrated regressor is dropped by the same rule.
            # A trial whose onset is after the end of the run (e.g. a timing
            # file that does not match the run) has no column at all, because
            # AFNI ignores such onsets; it is listed as dropped for that reason.
            #
            # Written here:
            #   <PREFIX>_IM_trial_censoring.tsv  every rated trial, with
            #       signal_kept and status (kept / dropped)
            #   <PREFIX>_IM_trials.tsv, <PREFIX>_IM_times.1D  the KEPT trials,
            #       renumbered 0..K-1 so that row k = sub-brick k of the betas.
            #       Extra columns: orig_trial (its number in _IM_alltrials.tsv)
            #       and signal_kept.
            # A run with fewer than 2 kept trials is skipped.
            # ------------------------------------------------------------------
            KEPT_TABLE=$(xmat_signal_kept \
                "${OUT_DIR}/${PREFIX}_IM_xmat_alltrials.1D" \
                "${OUT_DIR}/${PREFIX}_IM_xmat_alltrials_uncensored.1D")
            awk -F'\t' -v min="${MIN_SIGNAL_KEPT}" -v times="${TIMES}" -v trials="${TRIALS}" '
                NR == FNR { if ($1 ~ /^trial#/) { k = substr($1, 7) + 0; kept[k] = $2; inrun[k] = $3 }; next }
                FNR == 1  { print $0 "\tsignal_kept\tstatus"
                            print "trial\tonset\trating\tcond\torig_trial\tsignal_kept" > trials
                            n = 0; line = ""; next }
                {
                    k = $1 + 0
                    if (!(k in kept)) {
                        # AFNI leaves out any onset after the end of the run, so
                        # that trial has no column at all. (Onsets are in time
                        # order, so this can only happen to the last trials and
                        # does not shift the other trial numbers.)
                        print $0 "\t\tdropped: onset is after the end of the run"
                        next
                    }
                    if (kept[k] + 0 > min) {
                        status = "kept"
                        printf "%d\t%s\t%s\t%s\t%d\t%s\n", n, $2, $3, $4, k, kept[k] > trials
                        line = line (n ? " " : "") $2
                        n++
                    } else if (inrun[k] == "no") {
                        status = "dropped: onset is after the end of the run"
                    } else {
                        status = "dropped: response fully censored"
                    }
                    print $0 "\t" kept[k] "\t" status
                }
                END { print line > times }
            ' <(printf '%s\n' "${KEPT_TABLE}") "${ALLTRIALS}" > "${CENSCHECK}"

            NTRIAL=$(($(wc -l < "${TRIALS}") - 1))
            NDROP=$((NALL - NTRIAL))
            if [ "${NDROP}" -gt 0 ]; then
                echo "  WARNING: ${NDROP} of ${NALL} trial(s) dropped (no data left after censoring):"
                awk -F'\t' 'NR > 1 && $6 != "kept" { printf "    trial %s (onset %s s): %s\n", $1, $2, $6 }' "${CENSCHECK}"
            fi
            CHANGED=${NDROP}
            if [ "${USE_UNRATED}" = 1 ]; then
                UNRATED_KEPT=$(printf '%s\n' "${KEPT_TABLE}" | awk -F'\t' '$1 == "unrated#0" { print $2 }')
                if ! awk -v v="${UNRATED_KEPT:-0}" -v min="${MIN_SIGNAL_KEPT}" 'BEGIN { exit !(v + 0 > min) }'; then
                    echo "  WARNING: unrated regressor dropped (no data left after censoring)"
                    USE_UNRATED=0
                    CHANGED=1
                fi
            fi
            if [ "${NTRIAL}" -lt 2 ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: only ${NTRIAL} trial(s) left after censoring"
                continue
            fi
            NHIGH=$(awk -F'\t' 'NR>1 && $4=="high"' "${TRIALS}" | wc -l | tr -d ' ')

            echo "=== ${PARTICIPANT_ID} ${SES} task-${TASK}: ${NTRIAL} trials (${NHIGH} high / $((NTRIAL - NHIGH)) low), BLOCK ${BLOCK_DUR}s ==="

            # ------------------------------------------------------------------
            # STEP 4. The final design matrix.
            #
            # If nothing was dropped, the STEP 2 design is the final one and is
            # just copied. Otherwise the design is rebuilt exactly as in STEP 2,
            # but with only the kept trials (and without the unrated regressor
            # if it was dropped).
            # ------------------------------------------------------------------
            if [ "${CHANGED}" = 0 ]; then
                cp -f "${OUT_DIR}/${PREFIX}_IM_xmat_alltrials.1D"  "${OUT_DIR}/${PREFIX}_IM_xmat.1D"
                cp -f "${OUT_DIR}/${PREFIX}_IM_xmat_alltrials.jpg" "${OUT_DIR}/${PREFIX}_IM_xmat.jpg"
            else
                ${AFNI} "
                    cd ${OUT_DIR}

                    3dDeconvolve \
                        -input ${SCALED} \
                        -mask ${MASK} \
                        -polort A \
                        -TR_times ${TR} \
                        -censor ${CENSOR} \
                        -ortvec ${MOTION} motion \
                        $(stim_opts "${TIMES}" "${USE_UNRATED}") \
                        -x1D ${OUT_DIR}/${PREFIX}_IM_xmat.1D \
                        -xjpeg ${OUT_DIR}/${PREFIX}_IM_xmat.jpg \
                        -x1D_stop
                "
            fi

            # ------------------------------------------------------------------
            # STEP 5. Fit the GLM with 3dREMLfit and keep only the trial betas.
            #
            # As in first_level_afni.sh, the noise in each voxel is modelled as
            # temporally autocorrelated (ARMA(1,1)).
            # -Rbeta: save the betas (no t-statistics are needed here).
            # -nobout: leave out the betas of the baseline columns (drift and
            #   motion), so <PREFIX>_IM_allbetas holds trial#0 .. trial#K-1,
            #   followed by unrated#0 when the unrated regressor is present
            #   (verified on synthetic data).
            # 3dbucket then copies only the first K sub-bricks (the trials) into
            # <PREFIX>_IM_betas, and the temporary _IM_allbetas is deleted.
            # Old copies of both are removed first. trr_export_afni.sh later
            # splits _IM_betas into one file per trial.
            # ------------------------------------------------------------------
            rm -f "${OUT_DIR:?}/${PREFIX:?}"_IM_betas+tlrc.* "${OUT_DIR:?}/${PREFIX:?}"_IM_allbetas+tlrc.*
            ${AFNI} "
                3dREMLfit \
                    -matrix ${OUT_DIR}/${PREFIX}_IM_xmat.1D \
                    -input  ${SCALED} \
                    -mask   ${MASK} \
                    -nobout \
                    -Rbeta  ${OUT_DIR}/${PREFIX}_IM_allbetas

                3dbucket -prefix ${OUT_DIR}/${PREFIX}_IM_betas \
                    ${OUT_DIR}/${PREFIX}_IM_allbetas+tlrc'[0..$((NTRIAL - 1))]'
            "
            rm -f "${OUT_DIR:?}/${PREFIX:?}"_IM_allbetas+tlrc.*

            # ------------------------------------------------------------------
            # STEP 6. Check that sub-brick k really is trial k.
            #
            # The beta dataset must have exactly K sub-bricks (3dinfo -nv), and
            # sub-brick k must be labelled "trial#k" (3dinfo -label). If either
            # check fails, the trial table would not line up with the betas, so
            # the whole job stops with an error.
            # ------------------------------------------------------------------
            NV=$(${AFNI} "3dinfo -nv ${OUT_DIR}/${PREFIX}_IM_betas+tlrc" 2>/dev/null | tail -n1)
            NLAB=$(${AFNI} "3dinfo -label ${OUT_DIR}/${PREFIX}_IM_betas+tlrc" 2>/dev/null | tail -n1 \
                | tr '|' '\n' | awk '$0 == "trial#" (NR - 1)' | wc -l | tr -d ' ')
            if [ "${NV}" != "${NTRIAL}" ] || [ "${NLAB}" != "${NTRIAL}" ]; then
                echo "ERROR: ${PREFIX}_IM_betas has ${NV} sub-bricks (${NLAB} labelled trial#0..) but ${NTRIAL} trials -- trial table would be misaligned." >&2
                exit 1
            fi

            echo "Done: ${PARTICIPANT_ID} ${SES} task-${TASK} -> ${OUT_DIR}/${PREFIX}_IM_betas+tlrc (${NV} trials)"
        done
    done
done

echo "All single-trial first-level jobs complete."
