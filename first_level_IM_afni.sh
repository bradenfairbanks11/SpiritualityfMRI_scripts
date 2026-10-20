#!/bin/bash
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
# Run the job in scratch, not the submit dir. #SBATCH is parsed by Slurm before
# the shell runs, so this must be a literal absolute path -- ${WORK_ROOT} would
# not expand. Keeps stray relative-path output off /home.
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# =============================================================================
# Spirituality fMRI — SINGLE-TRIAL First-Level GLM (AFNI -stim_times_IM)
#
# Companion to first_level_afni.sh. Same data, same nuisance model, but every
# trial gets its OWN regressor (IM = "individually modulated"), so the output is
# one beta map PER TRIAL instead of one averaged beta. Those trial-level betas
# are the input to AFNI's trial-level test-retest reliability tools:
#   3dLMEr -TRR (whole-brain)  and  TRR (region-level)
#   Chen et al. (2021) NeuroImage 245:118647, "Trial and error: A hierarchical
#   modeling approach to test-retest reliability".
#
# Does NOT redo preprocessing. It reuses first_level_afni.sh's intermediates
# from the same run directory, so smoothing / scaling / motion / censoring are
# identical by construction:
#   <PREFIX>_blurred_scaled.nii.gz   (STEP 2)
#   <PREFIX>_motion_demean.txt       (STEP 4)
#   <PREFIX>_motion_censor.1D        (STEP 4)
# Run first_level_afni.sh first; runs without those files are skipped.
#
# Trials and ratings come from the SAME *_ratingAM.1D timing file ("onset*rating",
# sorted by onset), so the trial set is exactly the AM model's trial set.
# Trials with no button press have no rating, so they are not in the AM file and
# get no trial beta (nothing to split low/high on). They are still modeled: when
# byv-<BYV>_<ses>_task-<TASK>_unrated.1D exists they enter as one nuisance
# regressor, exactly as in first_level_afni.sh, so their response does not fall
# into the baseline. Its beta is trimmed from the output (see STEP 3).
#
# Output, per participant / session / task, in <run dir>/single_trial/:
#   <PREFIX>_IM_times.1D     onsets only (the IM timing file)
#   <PREFIX>_IM_trials.tsv   trial (0-based = sub-brick index), onset, rating,
#                            cond (low = rating 1-2, high = rating 3-4)
#   <PREFIX>_IM_betas+tlrc   one sub-brick per trial (REML betas, % signal
#                            change; baseline + motion betas excluded via -nobout)
#
# Run with:  sbatch first_level_IM_afni.sh sub-04 sub-05   (or no args = all)
# =============================================================================

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# --- Locate config.sh (see first_level_afni.sh for why this resolver exists) ---
CONFIG=""
for _c in "${PIPELINE_CONFIG:-}" \
          "${SLURM_SUBMIT_DIR:-}/config.sh" \
          "${SLURM_SUBMIT_DIR:-}/../config.sh" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/config.sh" \
          "/home/bradenf4/spirituality_fmri/scripts/config.sh"; do
    case "${_c}" in ""|"/config.sh"|"/../config.sh") continue ;; esac
    [ -f "${_c}" ] && { CONFIG="${_c}"; break; }
done
[ -n "${CONFIG}" ] || { echo "FATAL: cannot locate config.sh (set PIPELINE_CONFIG)" >&2; exit 1; }
source "${CONFIG}"
# Provides: FMRIPREP_OUT TEDANA_OUT AFNI_OUT TIMING_DIR PARTICIPANTS_TSV AFNI_SIF LOGS BIND_ARGS

TR=1.388

# --- Rating split for the 2-condition 3dLMEr -TRR model ---
# A trial is "high" if its 1-4 rating is >= HIGH_MIN, otherwise "low".
HIGH_MIN=3

TASKS=(scripture FHS architecture)

# Same fixed BLOCK durations as first_level_afni.sh (keep the two in sync).
block_dur_for() {
    case "$1" in
        scripture)    echo 19 ;;
        FHS)          echo 44 ;;
        architecture) echo 14 ;;
        *)            echo "ERROR: unknown task $1" >&2; exit 1 ;;
    esac
}

# ==============================================================================

mkdir -p "${LOGS}"

module load apptainer/1.3.6-qycanb2
export APPTAINERENV_AFNI_DECONFLICT=OVERWRITE
AFNI="apptainer exec ${BIND_ARGS} ${AFNI_SIF} bash -c"

if [ ! -f "${PARTICIPANTS_TSV}" ]; then
    echo "ERROR: participants.tsv not found at ${PARTICIPANTS_TSV}"
    exit 1
fi

# BIDS sub-NN -> BYV id (identical to first_level_afni.sh).
byv_for_sub() {
    local sub="$1"
    awk -F'\t' -v want="${sub}" '
        NR==1 { for (i=1;i<=NF;i++){ if($i=="participant_id")p=i; if($i=="byv_id")b=i}; next }
        $p==want { print $b; exit }
    ' "${PARTICIPANTS_TSV}"
}

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
            if why=$(is_excluded "${PARTICIPANT_ID}" "${SES}" "${TASK}"); then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: EXCLUDED (exclusions.tsv) -- ${why}"
                continue
            fi
            RUN_DIR=${AFNI_OUT}/${PARTICIPANT_ID}/${SES}/task-${TASK}
            PREFIX=${PARTICIPANT_ID}_${SES}_task-${TASK}
            SCALED=${RUN_DIR}/${PREFIX}_blurred_scaled.nii.gz
            MOTION=${RUN_DIR}/${PREFIX}_motion_demean.txt
            CENSOR=${RUN_DIR}/${PREFIX}_motion_censor.1D
            TIMING=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_ratingAM.1D
            UNRATED=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_unrated.1D
            MASK=${FMRIPREP_OUT}/${PARTICIPANT_ID}/${SES}/func/${PARTICIPANT_ID}_${SES}_task-${TASK}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz

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
            TIMES=${OUT_DIR}/${PREFIX}_IM_times.1D
            TRIALS=${OUT_DIR}/${PREFIX}_IM_trials.tsv

            # ------------------------------------------------------------------
            # STEP 1: split the AM timing file ("onset*rating" tokens, one row)
            # into an onset-only IM timing file and a per-trial table. Trial k
            # here == sub-brick k of the IM beta dataset (AFNI keeps file order).
            # ------------------------------------------------------------------
            tr ' \t' '\n\n' < "${TIMING}" | awk -v hi="${HIGH_MIN}" -v times="${TIMES}" '
                BEGIN { print "trial\tonset\trating\tcond"; n = 0; line = "" }
                $0 ~ /\*/ {
                    split($0, p, "*")
                    printf "%d\t%s\t%s\t%s\n", n, p[1], p[2], (p[2] + 0 >= hi ? "high" : "low")
                    line = line (n ? " " : "") p[1]
                    n++
                }
                END { print line > times }
            ' > "${TRIALS}"
            NTRIAL=$(($(wc -l < "${TRIALS}") - 1))
            if [ "${NTRIAL}" -lt 2 ]; then
                echo "Skipping ${PARTICIPANT_ID} ${SES} task-${TASK}: only ${NTRIAL} trial(s) in ${TIMING}"
                continue
            fi
            NHIGH=$(awk -F'\t' 'NR>1 && $4=="high"' "${TRIALS}" | wc -l | tr -d ' ')

            echo "=== ${PARTICIPANT_ID} ${SES} task-${TASK}: ${NTRIAL} trials (${NHIGH} high / $((NTRIAL - NHIGH)) low), BLOCK ${BLOCK_DUR}s ==="

            # ------------------------------------------------------------------
            # STEP 2: design matrix with one regressor per trial. Same nuisance
            # model as first_level_afni.sh (polort A, Friston-24 ortvec, FD
            # censor). -x1D_stop: 3dREMLfit does the actual fit.
            # Stimulus 2 'unrated' = no-response trials as one nuisance regressor
            # (only if the run has any), same as first_level_afni.sh.
            # ------------------------------------------------------------------
            NSTIM=1
            UNRATED_OPTS=""
            if [ -f "${UNRATED}" ]; then
                NSTIM=2
                UNRATED_OPTS="-stim_times 2 ${UNRATED} 'BLOCK(${BLOCK_DUR},1)' -stim_label 2 unrated"
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
                    -num_stimts ${NSTIM} \
                    -ortvec ${MOTION} motion \
                    -local_times \
                    -stim_times_IM 1 ${TIMES} 'BLOCK(${BLOCK_DUR},1)' \
                        -stim_label 1 trial \
                    ${UNRATED_OPTS} \
                    -x1D ${OUT_DIR}/${PREFIX}_IM_xmat.1D \
                    -xjpeg ${OUT_DIR}/${PREFIX}_IM_xmat.jpg \
                    -x1D_stop
            "

            # ------------------------------------------------------------------
            # STEP 3: REML fit (ARMA(1,1) noise, like the AM GLM). -nobout drops
            # baseline + motion betas, so -Rbeta holds trial#0 .. trial#N-1,
            # followed by unrated#0 when the unrated regressor is present
            # (verified on synthetic data). Only the N trial sub-bricks are kept
            # in _IM_betas, and their labels are checked, so sub-brick k is
            # always trial k of the table. trr_export_afni.sh later splits them
            # into one file per trial.
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
