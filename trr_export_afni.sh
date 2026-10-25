#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   4 cores / 16 GB / 2 h : one job loops over the three tasks; the work is
#                           mostly file splitting and ROI averaging.
#   --chdir               : run inside scratch, so any stray relative-path file
#                           lands there instead of in the home directory.
# ==============================================================================
#SBATCH --job-name=spirituality_trr_export
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# ==============================================================================
# trr_export_afni.sh — EXPORT TRIAL-LEVEL DATA FOR TEST-RETEST RELIABILITY (TRR)
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   first_level_IM_afni.sh estimated one beta map PER TRIAL for each participant,
#   session and task. This script packages those trial betas so that AFNI's
#   trial-level test-retest reliability tools (3dLMEr -TRR for the whole brain,
#   TRR for regions) can ask: how consistent is each participant's response
#   from session 1 to session 2?
#   For every task it keeps the participants who have trial betas for BOTH
#   sessions, and writes:
#     - one 3D beta map per trial,
#     - a table listing every trial (participant, session, trial, rating,
#       low/high condition, file name),
#     - a brain mask for the task,
#     - regions of interest (ROIs) = the significant clusters of the group
#       mean_response map from second_level_afni.sh,
#     - a table of each trial's mean beta in each ROI.
#   It only splits and summarizes files (AFNI C programs, no R). The R steps
#   (3dLMEr -TRR, TRR) are run on the local Mac by local_r_steps.sh.
#
# WHERE IT SITS IN THE PIPELINE
#   Before: first_level_IM_afni.sh  (single-trial betas + trial tables)
#           second_level_afni.sh    (group mean_response 3dMEMA maps + the
#                                    3dClustSim cluster-size table)
#   After : copy the whole TRR_OUT folder to the local Mac -> local_r_steps.sh
#
# STEPS, PER TASK
#   (once, before the tasks) read the cluster-size threshold from the
#   second-level 3dClustSim table.
#   1. Find the participants with single-trial betas for BOTH sessions
#      (skipping excluded runs).
#   2. Split each run's trial-beta bucket into one file per trial, and write
#      the trial table.
#   3. Check that every file listed in the trial table was written.
#   4. Build the task's analysis mask (intersection of the runs' brain masks).
#   5. Make ROIs from the group mean_response map (cluster-corrected).
#   6. Average each trial's beta within each ROI and write the ROI table.
#
# INPUTS
#   trial betas  : <AFNI_OUT>/sub-NN/ses-N/task-T/single_trial/sub-NN_ses-N_task-T_IM_betas+tlrc
#                    (sub-brick k = trial k, 0-based)
#   trial table  : <AFNI_OUT>/sub-NN/ses-N/task-T/single_trial/sub-NN_ses-N_task-T_IM_trials.tsv
#                    columns: trial (0-based), onset, rating, cond (low/high)
#   brain masks  : <FMRIPREP_OUT>/sub-NN/ses-N/func/
#                    sub-NN_ses-N_task-T_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz
#   group map    : <GROUP_OUT>/mema/mema_T_mean+tlrc     (second_level_afni.sh STEP 3)
#   cluster table: <GROUP_OUT>/clustsim/mema_clustsim.NN1_bisided.1D  (STEP 6 there)
#   exclusions   : <EXCLUSIONS_TSV> (+ <PARTICIPANTS_TSV> to map sub-NN -> BYV id)
#
# OUTPUTS, all in <TRR_OUT>/ (every path inside is relative, so the folder can
# be copied whole to another machine)
#   inputs/<sub>_<ses>_task-<T>_trial-NNN.nii.gz   one 3D beta map per trial
#                                                  (NNN = 0-based trial number)
#   trials_task-<T>.tsv        Subj sess trial rating cond InputFile
#                              -> the 3dLMEr -TRR data table (whole brain)
#   mask_task-<T>.nii.gz       intersection of every included run's fMRIPrep mask
#   rois_task-<T>.nii.gz       ROI map: cluster 1 = largest, 2 = next, ...
#   rois_task-<T>_clusters.txt 3dClusterize report (size, peak xyz per ROI)
#   roi_task-<T>.tsv           Subj sess trial rating cond ROI Y
#                              -> the region-level TRR data table
#
# HOW TO RUN
#   sbatch trr_export_afni.sh
#   Re-running overwrites the outputs in place.
#
# METHODS SUMMARY (the exact values used below)
#   Participants : per task, those with single-trial betas for both ses-1 and
#                  ses-2 and neither run in exclusions.tsv; a task needs >= 2.
#   Trial data   : the REML single-trial betas from first_level_IM_afni.sh
#                  (percent signal change), one file per trial. Session labels
#                  are written as ses1 / ses2. cond = low (rating 1-2) or high
#                  (rating 3-4), as set in first_level_IM_afni.sh.
#   Mask         : per task, the voxels inside the fMRIPrep MNI152NLin2009cAsym
#                  brain mask of EVERY included run (3dmask_tool -frac 1.0).
#   ROIs         : clusters of the group 3dMEMA mean_response map for that task
#                  (NOT the rating slope, and NOT session-1-vs-2 agreement, so
#                  ROIs are not chosen using the reliability being measured).
#                  Voxelwise threshold p < 0.001 two-sided, with positive and
#                  negative clusters formed separately ("bisided"); clusters
#                  must touch face-to-face (NN1); minimum cluster size = the
#                  second-level 3dClustSim value for p = 0.001, alpha = 0.05
#                  (NN1, bisided), rounded up; searched inside the task mask.
#   ROI values   : mean beta over each ROI's voxels for every trial
#                  (3dROIstats).
# ==============================================================================

set -o pipefail    # a pipeline fails if ANY command in it fails
                   # (note: no "set -e" -- a failing command does not stop the script)


# ==============================================================================
# SECTION 1. SETTINGS
# ------------------------------------------------------------------------------
# Every path and parameter this script uses is defined here. Settings marked
# [shared] also appear in other scripts and must match them; after editing one,
# run  bash check_settings.sh  to confirm every copy still agrees.
# ==============================================================================

# --- 1a. Where the data live --------------------------------------------------
# Three storage areas on the BYU supercomputer:
#   /home/bradenf4                     backed up, never purged: code, logs
#   /nobackup/autodelete/usr/bradenf4  scratch, NOT backed up, purged after 12
#                                      weeks unused: everything this job reads/writes
#   /nobackup/archive/usr/bradenf4     raw data (read only) and archived results
AFNI_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_firstlvl     # [shared] input: single-trial betas
FMRIPREP_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/fmriprep      # [shared] input: brain masks
GROUP_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_group       # [shared] input: group maps
TRR_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_trr           # [shared] output
PARTICIPANTS_TSV=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata/participants.tsv  # [shared]
EXCLUSIONS_TSV=/home/bradenf4/spirituality_fmri/scripts/exclusions.tsv                 # [shared] runs to skip
LOGS=/home/bradenf4/spirituality_fmri/logs                                             # [shared]

# Files made by second_level_afni.sh:
MEMADIR=${GROUP_OUT}/mema                                         # 3dMEMA maps (STEP 3 there)
CLUST_TABLE=${GROUP_OUT}/clustsim/mema_clustsim.NN1_bisided.1D    # cluster table (STEP 6 there)

# --- 1b. Software -------------------------------------------------------------
# AFNI runs inside an Apptainer image. The image only sees folders that are
# explicitly "bound" into it, so all three storage areas are bound.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]

# --- 1c. Analysis parameters --------------------------------------------------
# Tasks, by their BIDS labels.
TASKS=(scripture FHS architecture)                                                     # [shared]

# The two sessions compared for test-retest reliability.
SES1=ses-1
SES2=ses-2

# Voxelwise p-value for the ROI clusters. Must match CLUST_PTHR in
# second_level_afni.sh, which built the cluster-size table for this p-value.
CLUST_PTHR=0.001


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# byv_for_sub sub-NN  ->  prints that participant's BYV id (e.g. 0AK1RT).
# BIDS labels participants sub-NN (numbered in scan-date order), but the lab's
# records (and exclusions.tsv) use the BYV id. participants.tsv maps one to the
# other via its 'participant_id' and 'byv_id' columns.
# (The awk variable is 'want', not 'sub': sub() is a built-in awk function.)
byv_for_sub() {
    awk -F'\t' -v want="$1" '
        NR==1 { for (i=1;i<=NF;i++){ if($i=="participant_id")p=i; if($i=="byv_id")b=i}; next }
        $p==want { print $b; exit }
    ' "${PARTICIPANTS_TSV}"
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

mkdir -p "${LOGS}" "${TRR_OUT}/inputs"

module load "${APPTAINER_MODULE}"

# Every AFNI command below runs as:  ${AFNI} "<commands>"
# --env AFNI_DECONFLICT=OVERWRITE lets AFNI programs overwrite existing files,
# so a re-run does not fail on outputs from an earlier run.
AFNI="apptainer exec --env AFNI_DECONFLICT=OVERWRITE ${BIND_ARGS} ${AFNI_SIF} bash -c"

# --- Minimum cluster size (CLUST_K, in voxels) --------------------------------
# Read from second_level_afni.sh's 3dClustSim table (NN1, bisided). The table
# has one row per voxelwise p-value (-pthr); column 1 is p and column 2 is the
# cluster size for alpha = 0.05 (the only -athr used). The row for CLUST_PTHR is
# taken and its value rounded UP to a whole voxel, since a cluster must be at
# least this big. If the table is missing, ROI tables are skipped (STEP 5).
CLUST_K=""
if [ -f "${CLUST_TABLE}" ]; then
    CLUST_K=$(awk -v p="${CLUST_PTHR}" '!/^#/ && NF>=2 && ($1+0)==(p+0) { k=$2; print (k==int(k)) ? k : int(k)+1; exit }' "${CLUST_TABLE}")
fi
if [ -z "${CLUST_K}" ]; then
    echo "WARNING: no cluster-size threshold for p=${CLUST_PTHR} in ${CLUST_TABLE}" >&2
    echo "         (run second_level_afni.sh first). ROI tables will be skipped." >&2
fi


# ==============================================================================
# SECTION 4. MAIN LOOP: one pass per task
# ==============================================================================

for TASK in "${TASKS[@]}"; do
    echo "=============================================================="
    echo "task-${TASK}"
    echo "=============================================================="

    # --------------------------------------------------------------------------
    # STEP 1. Find the participants with single-trial betas for BOTH sessions.
    #
    # A participant is kept only if, for ses-1 AND ses-2, the run is not in
    # exclusions.tsv and <PREFIX>_IM_betas+tlrc exists. Reliability needs both
    # sessions, so one missing session drops the participant from this task.
    # The task is skipped if fewer than 2 participants remain.
    # --------------------------------------------------------------------------
    SUBS=()
    for d in "${AFNI_OUT}"/sub-*; do
        [ -d "${d}" ] || continue
        sid=$(basename "${d}")
        ok=1
        for SES in "${SES1}" "${SES2}"; do
            if why=$(is_excluded "${sid}" "${SES}" "${TASK}"); then
                echo "  EXCLUDED (exclusions.tsv): ${sid} ${SES} task-${TASK} -- ${why}"
                ok=0
            fi
            b=${AFNI_OUT}/${sid}/${SES}/task-${TASK}/single_trial/${sid}_${SES}_task-${TASK}_IM_betas+tlrc.HEAD
            [ -f "${b}" ] || ok=0
        done
        [ "${ok}" -eq 1 ] && SUBS+=("${sid}")
    done
    n=${#SUBS[@]}
    if [ "${n}" -lt 2 ]; then
        echo "skip task-${TASK}: only ${n} subject(s) with single-trial betas for both sessions"
        continue
    fi
    echo "subjects (n=${n}): ${SUBS[*]}"

    # --------------------------------------------------------------------------
    # STEP 2. Split each run's trial-beta bucket into one file per trial, and
    #         write the trial table.
    #
    # Each run's _IM_betas bucket has one sub-brick per trial; sub-brick k is
    # trial k of _IM_trials.tsv (0-based). 3dbucket copies sub-brick k to
    # inputs/<PREFIX>_trial-<k, 3 digits>.nii.gz (old trial files of that run
    # are deleted first; 3dbucket messages are hidden, so a failure shows up
    # in STEP 3 instead).
    # Trial table, one row per trial:
    #   Subj   participant (sub-NN)
    #   sess   ses1 / ses2 (no hyphen, so it is a valid R factor level)
    #   trial  trial number (column 1 of _IM_trials.tsv)
    #   rating the 1-4 rating (column 3)
    #   cond   low / high   (column 4)
    #   InputFile  inputs/<PREFIX>_trial-NNN.nii.gz (relative path)
    # The run's fMRIPrep brain mask is also added to the list for STEP 4.
    # --------------------------------------------------------------------------
    TABLE=${TRR_OUT}/trials_task-${TASK}.tsv
    printf 'Subj\tsess\ttrial\trating\tcond\tInputFile\n' > "${TABLE}"
    MASKS=""
    for sid in "${SUBS[@]}"; do
        for SES in "${SES1}" "${SES2}"; do
            ST=${AFNI_OUT}/${sid}/${SES}/task-${TASK}/single_trial
            PREFIX=${sid}_${SES}_task-${TASK}
            TRIALS=${ST}/${PREFIX}_IM_trials.tsv
            NTRIAL=$(($(wc -l < "${TRIALS}") - 1))     # rows minus the header
            rm -f "${TRR_OUT:?}/inputs/${PREFIX:?}"_trial-*.nii.gz
            ${AFNI} "
                for k in \$(seq 0 $((NTRIAL - 1))); do
                    3dbucket -prefix ${TRR_OUT}/inputs/${PREFIX}_trial-\$(printf %03d \$k).nii.gz ${ST}/${PREFIX}_IM_betas+tlrc\"[\$k]\" >/dev/null 2>&1
                done
            "
            seslab=${SES//-/}   # ses-1 -> ses1 (a valid R factor level)
            awk -F'\t' -v s="${sid}" -v ses="${seslab}" -v pre="${PREFIX}" '
                NR>1 { printf "%s\t%s\t%d\t%s\t%s\tinputs/%s_trial-%03d.nii.gz\n", s, ses, $1, $3, $4, pre, $1 }
            ' "${TRIALS}" >> "${TABLE}"
            MASKS="${MASKS} ${FMRIPREP_OUT}/${sid}/${SES}/func/${PREFIX}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
        done
    done

    # --------------------------------------------------------------------------
    # STEP 3. Check that every file listed in the trial table was written.
    #
    # Each InputFile is opened for reading; any that cannot be read is counted.
    # If any are missing, the rest of this task is skipped, because the tables
    # would no longer line up with the files.
    # --------------------------------------------------------------------------
    nmiss=$(awk -F'\t' -v root="${TRR_OUT}" 'NR>1 { f=root "/" $6; if ((getline line < f) < 0) m++; close(f) } END { print m+0 }' "${TABLE}")
    if [ "${nmiss}" -gt 0 ]; then
        echo "ERROR task-${TASK}: ${nmiss} trial file(s) listed in ${TABLE} were not written." >&2
        continue
    fi
    echo "  trial table: ${TABLE} ($(($(wc -l < "${TABLE}") - 1)) trials)"

    # --------------------------------------------------------------------------
    # STEP 4. Build the task's analysis mask.
    #
    # 3dmask_tool -frac 1.0 keeps only voxels that are inside EVERY included
    # run's fMRIPrep brain mask (both sessions of every participant).
    # --------------------------------------------------------------------------
    TMASK=${TRR_OUT}/mask_task-${TASK}.nii.gz
    ${AFNI} "3dmask_tool -overwrite -input ${MASKS} -frac 1.0 -prefix ${TMASK}"

    # --------------------------------------------------------------------------
    # STEP 5. Make ROIs from the group mean_response map.
    #
    # The map is second_level_afni.sh's 3dMEMA one-sample result for this
    # task's mean_response: sub-brick [0] = group effect, [1] = its
    # t-statistic. 3dClusterize:
    #   -ithr 1 -idat 0       threshold on the t-statistic, keep the effect values
    #   -mask TMASK           search only inside the STEP 4 mask
    #   -NN 1                 voxels are neighbours only if they share a face
    #   -bisided p=CLUST_PTHR two-sided p < 0.001; positive and negative voxels
    #                         form separate clusters
    #   -clust_nvox CLUST_K   keep clusters of at least CLUST_K voxels
    #   -pref_map ROIS        write the ROI map (1 = largest cluster, 2 = next..)
    #   -outvol_if_no_clust   still write an (empty) map if nothing survives
    # The cluster report goes to rois_task-<T>_clusters.txt. The number of ROIs
    # is the largest value in the ROI map (3dBrickStat -max). The task gets no
    # ROI table if the map or threshold is missing or no cluster survives.
    # --------------------------------------------------------------------------
    MEMA_MEAN=${MEMADIR}/mema_${TASK}_mean+tlrc
    ROIS=${TRR_OUT}/rois_task-${TASK}.nii.gz
    ROITAB=${TRR_OUT}/roi_task-${TASK}.tsv
    rm -f "${ROIS}" "${ROITAB}"
    if [ -z "${CLUST_K}" ] || [ ! -f "${MEMA_MEAN}.HEAD" ]; then
        echo "  WARN: no MEMA mean map or cluster threshold -> no ROI table for task-${TASK}" >&2
        continue
    fi
    # MEMA output: [0] = effect, [1] = t-statistic.
    ${AFNI} "3dClusterize -inset ${MEMA_MEAN} -ithr 1 -idat 0 -mask ${TMASK} \
        -NN 1 -bisided p=${CLUST_PTHR} -clust_nvox ${CLUST_K} \
        -pref_map ${ROIS} -outvol_if_no_clust" > "${TRR_OUT}/rois_task-${TASK}_clusters.txt" 2>&1
    NROI=$(${AFNI} "3dBrickStat -max ${ROIS}" 2>/dev/null | tail -n1 | awk '{print int($1)}')
    if [ -z "${NROI}" ] || [ "${NROI}" -lt 1 ]; then
        echo "  WARN task-${TASK}: no group cluster survived (p<${CLUST_PTHR}, k>=${CLUST_K}) -> no ROI table" >&2
        continue
    fi
    echo "  ROIs: ${NROI} cluster(s) -> ${ROIS}"

    # --------------------------------------------------------------------------
    # STEP 6. Average each trial's beta within each ROI; write the ROI table.
    #
    # 3dROIstats -mask ROIS -quiet prints one line per input file (in input
    # order), with one mean per ROI (ROI values 1..NROI). The trial files are
    # given in trial-table order, so line i belongs to table row i. If the line
    # count does not match the row count, no ROI table is written.
    # Output: one row per trial x ROI, with the trial's table columns plus
    #   ROI = ROI01, ROI02, ...   and   Y = that trial's mean beta in that ROI.
    # Each line is split explicitly with split(..., " ") (which ignores leading
    # blanks): changing FS between files is unreliable on BSD/macOS awk, where
    # the first record of the second file is still split with the old FS.
    # --------------------------------------------------------------------------
    FILES=$(awk -F'\t' -v root="${TRR_OUT}" 'NR>1 { printf "%s/%s ", root, $6 }' "${TABLE}")
    ${AFNI} "3dROIstats -mask ${ROIS} -quiet ${FILES}" > "${TRR_OUT}/.roistats_${TASK}.txt" 2>/dev/null
    NLINES=$(wc -l < "${TRR_OUT}/.roistats_${TASK}.txt" | tr -d ' ')
    NROWS=$(($(wc -l < "${TABLE}") - 1))
    if [ "${NLINES}" -ne "${NROWS}" ]; then
        echo "ERROR task-${TASK}: 3dROIstats gave ${NLINES} lines for ${NROWS} trials -- not writing ${ROITAB}" >&2
        continue
    fi
    awk -F'\t' 'BEGIN { print "Subj\tsess\ttrial\trating\tcond\tROI\tY" }
         NR==FNR { if (FNR>1) meta[FNR-1] = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5; next }
         { nv = split($0, v, " "); for (r = 1; r <= nv; r++) printf "%s\tROI%02d\t%s\n", meta[FNR], r, v[r] }' \
        "${TABLE}" "${TRR_OUT}/.roistats_${TASK}.txt" > "${ROITAB}"
    rm -f "${TRR_OUT}/.roistats_${TASK}.txt"
    echo "  ROI table: ${ROITAB}"
done

echo "=============================================================="
echo "TRR export complete: ${TRR_OUT}"
echo "Copy the WHOLE folder to the local machine, e.g. from the local machine:"
echo "    rsync -av <cluster>:${TRR_OUT}/ <local>/afni_trr/"
echo "then run local_r_steps.sh there."
echo "=============================================================="
