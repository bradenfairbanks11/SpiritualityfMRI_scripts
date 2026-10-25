#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   8 cores / 32 GB / 8 h : one job runs the whole group analysis. The 8 cores
#                           are passed to 3dttest++ -Clustsim (see NTHR).
#   --chdir               : run inside scratch, so any stray relative-path file
#                           lands there instead of in the home directory.
# ==============================================================================
#SBATCH --job-name=spirituality_grouplvl
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# ==============================================================================
# second_level_afni.sh — SECOND-LEVEL (GROUP) ANALYSIS OF THE RATING EFFECT
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   first_level_afni.sh fit one GLM per participant, session and task. Each GLM
#   gives two effects (as a beta "Coef" map and a t-statistic "Tstat" map):
#     rating_slope  : how much the BOLD response goes up per "felt the Spirit"
#                     rating point (1-4). THE PRIMARY RESULT.
#     mean_response : the average BOLD response to the task's trials (a sanity
#                     check / localizer).
#   This script combines those maps across participants. It asks, voxel by
#   voxel: is the effect different from zero in the group (per task), and does
#   the rating slope differ between tasks? It also builds the tables needed to
#   correct the group maps for multiple comparisons (cluster-size thresholds).
#
# WHERE IT SITS IN THE PIPELINE
#   Before: first_level_afni.sh  (per-run GLM; its REML stats and residuals)
#   After : trr_export_afni.sh   (uses this script's mean_response 3dMEMA maps
#                                 and cluster table to define ROIs)
#           promote_to_archive.sh (copies the results to archive; this script
#                                 no longer does that, see STEP 7)
#
# STEPS
#   0. Find every first-level result on disk (skipping excluded runs).
#   1. Pull out the slope and mean Coef/Tstat maps. If a participant has more
#      than one session of a task, average the sessions, so each participant
#      gives ONE map per task and effect.
#   2. Build the group brain mask.
#   3. 3dMEMA one-sample test per task (slope = primary, mean = sanity check).
#   4. 3dMEMA between-task slope contrasts (task X minus task Y).
#   5. 3dttest++ versions of STEPs 3 and 4 (with -Clustsim when n is big enough).
#   6. Cluster-size table for the 3dMEMA maps (3dFWHMx -> 3dClustSim).
#   7. (Nothing is copied to archive; this step only prints a reminder.)
#   The participant list is NOT hard-coded: whatever first-level output exists
#   is used, so re-running after more participants/sessions finish grows the
#   group automatically.
#
# INPUTS
#   first-level stats    : <AFNI_OUT>/sub-NN/ses-N/task-T/sub-NN_ses-N_task-T_stats_REML+tlrc
#                            sub-bricks rating_slope#0_Coef / _Tstat,
#                                       mean_response#0_Coef / _Tstat
#   first-level residuals: <AFNI_OUT>/sub-NN/ses-N/task-T/sub-NN_ses-N_task-T_errts_REML+tlrc
#   exclusions           : <EXCLUSIONS_TSV> (+ <PARTICIPANTS_TSV> to map sub-NN -> BYV id)
#
# OUTPUTS, all in <GROUP_OUT>/
#   inputs/sub-NN_task-T_<slope|mean>_Coef.nii.gz   one per participant x task x effect
#   inputs/sub-NN_task-T_<slope|mean>_T.nii.gz
#   inputs/contrast_X_vs_Y/sub-NN_bDelta.nii.gz, sub-NN_tDelta.nii.gz
#   mask/group_mask.nii.gz                          the group mask
#   mask/sub-NN_ses-N_task-T_modelled.nii.gz        each run's modelled voxels
#   mema/mema_T_slope+tlrc, mema/mema_T_mean+tlrc   STEP 3 (T = task)
#   mema/mema_X_vs_Y_slope+tlrc                     STEP 4
#   ttest/ttest_T_<slope|mean>+tlrc, ttest/ttest_X_vs_Y_slope+tlrc   STEP 5
#   clustsim/mema_clustsim.NN{1,2,3}_{1sided,2sided,bisided}.1D      STEP 6
#   clustsim/acf_params.txt                         per-run ACF estimates (STEP 6)
#
# HOW TO RUN
#   sbatch second_level_afni.sh
#   Re-running overwrites the outputs in place.
#
# METHODS SUMMARY (the exact values used below)
#   Effects        : rating_slope (primary) and mean_response, from the 3dREMLfit
#                    first-level fits; Coef and Tstat sub-bricks of each.
#   Exclusions     : runs listed in exclusions.tsv are skipped everywhere.
#   Sessions       : a participant with several sessions of a task is reduced
#                    to one map per task. Beta = plain (unweighted) mean of the
#                    session betas. Each session's squared standard error is
#                    SE^2 = (beta / t)^2 (set to 1e6 where t = 0); the combined
#                    t = mean beta / sqrt(mean(SE^2) / n_sessions).
#   Group mask     : voxels with non-zero first-level betas (rating_slope AND
#                    mean_response) in EVERY run that enters the group
#                    analysis (every participant, session and task), i.e. the
#                    voxels the first-level GLM modelled with data in all of
#                    them (fMRIPrep brain mask AND voxels with data at every
#                    time point, per run). No group voxel is zero in any input,
#                    so no test treats a missing voxel as a real 0.
#   3dMEMA         : mixed-effects meta-analysis, which uses each participant's
#                    beta AND t (so less precise participants count less).
#                    One-sample per task and effect. Zero = missing data
#                    (-missing_data 0). Other options at 3dMEMA defaults
#                    (this includes the Hartung-Knapp adjustment).
#                    Run only when n >= 2.
#   Task contrasts : scripture - architecture, FHS - architecture,
#                    scripture - FHS, on the rating slope. Per participant
#                    (those with both tasks): difference beta = bX - bY;
#                    difference t = (bX - bY) / sqrt(SE_X^2 + SE_Y^2); then a
#                    one-sample 3dMEMA.
#   3dttest++      : same one-sample tests (on betas only) and paired
#                    between-task slope tests. -Clustsim (randomization-based
#                    cluster tables; makes the output z-scores) is added only
#                    when there are >= 14 datasets in total.
#   Cluster table  : spatial smoothness (ACF model a, b, c) estimated with
#                    3dFWHMx -detrend -acf on every included first-level REML
#                    residual time series (uncensored volumes only) inside the
#                    group mask, averaged over runs;
#                    3dClustSim with voxelwise p = 0.001 and alpha = 0.05.
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
#   /home/bradenf4                     backed up, never purged: code, logs
#   /nobackup/autodelete/usr/bradenf4  scratch, NOT backed up, purged after 12
#                                      weeks unused: everything this job reads/writes
#   /nobackup/archive/usr/bradenf4     raw data (read only) and archived results
AFNI_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_firstlvl     # [shared] input: first-level results
GROUP_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_group       # [shared] output
FINAL_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/final    # [shared] archive (not written here)
PARTICIPANTS_TSV=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata/participants.tsv  # [shared]
EXCLUSIONS_TSV=/home/bradenf4/spirituality_fmri/scripts/exclusions.tsv                 # [shared] runs to skip
LOGS=/home/bradenf4/spirituality_fmri/logs                                             # [shared]
PROJ_HOME=/home/bradenf4/spirituality_fmri     # project home (only used in the final reminder message)

AFNI_FIRSTLVL=${AFNI_OUT}                      # the name used for the first-level folder below

# Archive folder for the group maps. NOT USED any more: STEP 7 used to copy the
# results here; promote_to_archive.sh does that now. Kept for reference.
GROUP_FINAL=${FINAL_DIR}/afni_group

# Sub-folders of GROUP_OUT (created in SECTION 3).
INPUTS=${GROUP_OUT}/inputs                     # per-participant maps (STEPs 1, 4)
MASKDIR=${GROUP_OUT}/mask                      # masks (STEP 2)
MEMADIR=${GROUP_OUT}/mema                      # 3dMEMA results (STEPs 3, 4)
TTESTDIR=${GROUP_OUT}/ttest                    # 3dttest++ results (STEP 5)
CLUSTDIR=${GROUP_OUT}/clustsim                 # smoothness + cluster tables (STEP 6)
GROUPMASK=${MASKDIR}/group_mask.nii.gz         # the group mask (STEP 2)

# --- 1b. Software -------------------------------------------------------------
# AFNI runs inside an Apptainer image. The image only sees folders that are
# explicitly "bound" into it, so all three storage areas are bound.
# 3dMEMA is an R program. Its R packages (data.table, snow) are INSIDE this
# image; compute nodes have no internet, so nothing can be installed at run
# time. SECTION 3 checks that the packages are there.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]

# --- 1c. Analysis parameters --------------------------------------------------
# Tasks, by their BIDS labels.
TASKS=(scripture FHS architecture)                                                     # [shared]

# Between-task slope contrasts, written X:Y, meaning "X minus Y" (STEPs 4, 5).
BETWEEN_TASK_PAIRS=(scripture:architecture FHS:architecture scripture:FHS)

# Sub-brick labels inside each first-level <PREFIX>_stats_REML+tlrc
# (they are the GLT results named in first_level_afni.sh STEP 5).
SLOPE_COEF="rating_slope#0_Coef"      # primary: beta of the rating slope
SLOPE_T="rating_slope#0_Tstat"        #          its t-statistic
MEAN_COEF="mean_response#0_Coef"      # sanity:  beta of the mean response
MEAN_T="mean_response#0_Tstat"        #          its t-statistic

# Voxelwise (uncorrected) p-value used to build the 3dClustSim cluster-size
# table (STEP 6). trr_export_afni.sh uses the same value (0.001).
CLUST_PTHR=0.001

# Number of parallel jobs for 3dttest++ -Clustsim: the CPUs Slurm gave this job
# (8, from #SBATCH --cpus-per-task), or 4 if run outside Slurm.
NTHR=${SLURM_CPUS_PER_TASK:-4}


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

# stats_dset sub-NN ses-N task  ->  prints the path of that run's first-level
# REML stats dataset (without the .HEAD/.BRIK ending).
stats_dset() {  # <sub> <ses> <task>
    echo "${AFNI_FIRSTLVL}/$1/$2/task-$3/$1_$2_task-$3_stats_REML+tlrc"
}

# extract_brick <stats_dset> <label> <out.nii.gz>
# Copies ONE sub-brick, chosen by its label (e.g. rating_slope#0_Coef), out of a
# stats dataset into its own NIfTI file (3dbucket). Any old copy is deleted first.
extract_brick() {  # <stats_dset> <label> <out.nii.gz>
    local dset="$1" lab="$2" out="$3"
    rm -f "${out}"
    ${AFNI} "3dbucket -prefix ${out} ${dset}'[${lab}]'"
}

# subs_with <task> <slope|mean>  ->  prints the participants that have an
# extracted Coef map for that task and effect. A participant must ALSO have had
# sessions found in STEP 0, so a stale map left in inputs/ by an earlier run of
# a since-excluded run is never pooled. (Used by STEPs 3 and 5 one-sample tests.)
subs_with() {  # <task> <effect>
    local task="$1" eff="$2" s
    for s in "${SUBJECTS[@]}"; do
        # Require sessions discovered in STEP 0 too, so a stale input file left by
        # an earlier run of a since-excluded run is never pooled.
        [ -n "${SESS[${s} ${task}]}" ] || continue
        [ -f "${INPUTS}/${s}_task-${task}_${eff}_Coef.nii.gz" ] && echo "${s}"
    done
}


# ==============================================================================
# SECTION 3. SETUP AND CHECKS
# ==============================================================================

module load "${APPTAINER_MODULE}"

# Every AFNI command below runs as:  ${AFNI} "<commands>"
# --cleanenv: the host's environment (e.g. a stray R_LIBS from a login shell) is
# NOT passed into the image, so it cannot override the image's R packages. Every
# variable a command needs is filled in on the host before the command string
# reaches the image. (This script does not set AFNI's overwrite option; instead
# it deletes each old output with rm -f before remaking it.)
AFNI="apptainer exec --cleanenv ${BIND_ARGS} ${AFNI_SIF} bash -c"

# Check that the image has 3dMEMA's R packages. This fails early with a clear
# message instead of deep inside 3dMEMA with an unclear R error.
# (The ERROR text below still says "config.sh"; AFNI_SIF is now set in SECTION 1b.)
if ! ${AFNI} "Rscript -e 'q(status=length(setdiff(c(\"data.table\",\"snow\"), rownames(installed.packages()))))'"; then
    echo "ERROR: 3dMEMA's R packages (snow, data.table) are missing from ${AFNI_SIF}."
    echo "Check that AFNI_SIF in config.sh points at the R-enabled image (see container_request.md)."
    exit 1
fi

mkdir -p "${LOGS}" "${GROUP_OUT}"/{inputs,mask,mema,ttest,clustsim}


# ==============================================================================
# SECTION 4. MAIN WORK
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0. Find every first-level result on disk.
#
# Looks for <AFNI_OUT>/sub-*/ses-*/task-*/*_stats_REML+tlrc.HEAD and records:
#   SESS["sub-NN task"] = the sessions found for that participant and task
#                         (e.g. " ses-1 ses-2")
#   SUBJECTS            = every participant found, sorted
# Runs in exclusions.tsv are skipped: their outputs may still be on disk from
# before they were excluded.
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 0: discovering first-level stats under ${AFNI_FIRSTLVL}"
echo "=============================================================="

declare -A SESS
SUBJECTS=()

for HEAD in "${AFNI_FIRSTLVL}"/sub-*/ses-*/task-*/*_stats_REML+tlrc.HEAD; do
    [ -e "${HEAD}" ] || continue
    # .../sub-XX/ses-Y/task-TASK/sub-XX_ses-Y_task-TASK_stats_REML+tlrc.HEAD
    task_dir=$(dirname "${HEAD}")
    ses_dir=$(dirname "${task_dir}")
    sub_dir=$(dirname "${ses_dir}")
    SUB=$(basename "${sub_dir}")
    SES=$(basename "${ses_dir}")
    TASK=$(basename "${task_dir}"); TASK=${TASK#task-}

    # Outputs of excluded runs may still be on disk from an earlier run.
    if why=$(is_excluded "${SUB}" "${SES}" "${TASK}"); then
        echo "  EXCLUDED (exclusions.tsv): ${SUB} ${SES} task-${TASK} -- ${why}"
        continue
    fi

    key="${SUB} ${TASK}"
    SESS["${key}"]="${SESS[${key}]} ${SES}"
    if [[ ! " ${SUBJECTS[*]} " == *" ${SUB} "* ]]; then
        SUBJECTS+=("${SUB}")
    fi
done

IFS=$'\n' SUBJECTS=($(sort <<<"${SUBJECTS[*]}")); unset IFS

if [ "${#SUBJECTS[@]}" -eq 0 ]; then
    echo "ERROR: no *_stats_REML+tlrc datasets found under ${AFNI_FIRSTLVL}"
    exit 1
fi

echo "Discovered subjects: ${SUBJECTS[*]}"
for SUB in "${SUBJECTS[@]}"; do
    for TASK in "${TASKS[@]}"; do
        s="${SESS[${SUB} ${TASK}]}"
        [ -n "${s}" ] && echo "  ${SUB}  task-${TASK}  ->${s}"
    done
done

# ------------------------------------------------------------------------------
# STEP 1. Extract the maps, and average sessions -> one map per participant x
#         task x effect.
#
# For each participant, task and effect (slope, mean) this writes:
#   inputs/<sub>_task-<TASK>_<effect>_Coef.nii.gz   the beta
#   inputs/<sub>_task-<TASK>_<effect>_T.nii.gz      its t-statistic
#
# One session : the Coef and Tstat sub-bricks are copied out as they are.
# Two or more : the sessions are combined ("fixed-effects" average):
#   beta_avg = mean of the session betas (plain mean, not weighted)
#   SE_i^2   = (beta_i / t_i)^2   the squared standard error of session i,
#              recovered from beta and t (since t = beta / SE). Where t_i = 0
#              the SE cannot be recovered, so SE_i^2 is set to 1e6, which drives
#              the combined t toward 0.
#   SE2_mean = mean of the SE_i^2
#   t_avg    = beta_avg / sqrt(SE2_mean / n)   (n = number of sessions; this is
#              the SE of a mean of n independent estimates). Where SE2_mean is
#              0 or below, t_avg is set to 0.
#
# Space label: the first-level data are in MNI152NLin2009cAsym, and
# first_level_afni.sh labels them as MNI (3drefit -space MNI), so AFNI reads the
# buckets as +tlrc. The NIfTI copies keep the MNI grid; 3drefit -space MNI is
# applied again to each output so clustering and whereami treat them as MNI.
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 1: extracting + session-averaging per subject x task"
echo "=============================================================="

for SUB in "${SUBJECTS[@]}"; do
    for TASK in "${TASKS[@]}"; do
        sessions=(${SESS[${SUB} ${TASK}]})
        [ "${#sessions[@]}" -eq 0 ] && continue

        for eff in slope mean; do
            if [ "${eff}" = slope ]; then clab="${SLOPE_COEF}"; tlab="${SLOPE_T}"
            else                          clab="${MEAN_COEF}";  tlab="${MEAN_T}"; fi

            outC=${INPUTS}/${SUB}_task-${TASK}_${eff}_Coef.nii.gz
            outT=${INPUTS}/${SUB}_task-${TASK}_${eff}_T.nii.gz

            if [ "${#sessions[@]}" -eq 1 ]; then
                D=$(stats_dset "${SUB}" "${sessions[0]}" "${TASK}")
                extract_brick "${D}" "${clab}" "${outC}"
                extract_brick "${D}" "${tlab}" "${outT}"
            else
                # Fixed-effects average across sessions (formulas in the STEP 1 box).
                #   a = Coef, b = Tstat in the first 3dcalc:
                #   step(abs(b))*(a/(b+equals(b,0)))^2   -> (beta/t)^2 where t != 0
                #   equals(b,0)*1000000                  -> 1e6 where t == 0
                n="${#sessions[@]}"
                tmp=${INPUTS}/.tmp_${SUB}_${TASK}_${eff}
                rm -f "${tmp}"_*.nii.gz
                clist=(); slist=()
                i=0
                for SES in "${sessions[@]}"; do
                    D=$(stats_dset "${SUB}" "${SES}" "${TASK}")
                    ci=${tmp}_c${i}.nii.gz
                    ti=${tmp}_t${i}.nii.gz
                    si=${tmp}_s${i}.nii.gz
                    extract_brick "${D}" "${clab}" "${ci}"
                    extract_brick "${D}" "${tlab}" "${ti}"
                    ${AFNI} "3dcalc -a ${ci} -b ${ti} \
                        -expr 'step(abs(b))*(a/(b+equals(b,0)))^2 + equals(b,0)*1000000' \
                        -prefix ${si}"
                    clist+=("${ci}"); slist+=("${si}")
                    i=$((i+1))
                done
                betaAvg=${tmp}_betaAvg.nii.gz
                se2mean=${tmp}_se2mean.nii.gz
                rm -f "${betaAvg}" "${se2mean}"
                ${AFNI} "3dMean -prefix ${betaAvg} ${clist[*]}"
                ${AFNI} "3dMean -prefix ${se2mean} ${slist[*]}"
                rm -f "${outC}" "${outT}"
                ${AFNI} "3dcopy ${betaAvg} ${outC}"
                # a = beta_avg, b = SE2_mean:  t = a * sqrt(n / b) where b > 0, else 0
                ${AFNI} "3dcalc -a ${betaAvg} -b ${se2mean} \
                    -expr 'step(b)*a*sqrt(${n}/(b+1-step(b)))' \
                    -prefix ${outT}"
                rm -f "${tmp}"_*.nii.gz
            fi

            # Label the space so clustering / whereami work (grid is truly MNI).
            ${AFNI} "3drefit -space MNI ${outC}; 3drefit -space MNI ${outT}"
        done
    done
done

echo "Extracted inputs:"
ls -1 "${INPUTS}"/*_Coef.nii.gz 2>/dev/null | sed 's/^/  /'

# ------------------------------------------------------------------------------
# STEP 2. Build the group brain mask.
#
# The group mask = voxels the first-level GLM modelled in EVERY run that goes
# into the group analysis (every participant x session x task found in STEP 0).
# For each run, a voxel counts as modelled when both its rating_slope and its
# mean_response beta are non-zero (step(abs(a)) = 1 for any non-zero value):
# first_level_afni.sh sets every voxel outside its analysis mask (fMRIPrep
# brain mask AND data at every time point) to exactly 0. Each run's mask is
# resampled to the reference grid with nearest-neighbour interpolation
# (3dresample -rmode NN, which keeps it 0/1), then 3dmask_tool -frac 1.0 keeps
# only voxels present in 100% of the run masks.
#
# WHY per run and not per participant: STEP 1 AVERAGES sessions, so a voxel
# with data in only one session would come out non-zero (half its beta) in the
# average. Requiring every run keeps such voxels out.
# WHY not fMRIPrep's masks (used until 2026-10-24): voxels inside fMRIPrep's mask
# can still have no data (tedana coverage, extents mask). Those would sit inside
# the group mask as zeros, which 3dMEMA can skip (-missing_data 0) but 3dttest++
# would treat as real values (Taylor et al. 2018, bioRxiv 10.1101/308643).
# Reference grid: the first slope Coef map from STEP 1 (all maps share the
# first-level grid).
# 3dinfo -same_grid then prints 1 for each map that is on the same grid as the
# mask (a check only; its result does not stop the script).
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 2: building group mask"
echo "=============================================================="

# Reference grid = first extracted slope Coef map.
REF=$(ls -1 "${INPUTS}"/*_task-*_slope_Coef.nii.gz 2>/dev/null | head -n1)
if [ -z "${REF}" ]; then echo "ERROR: no extracted Coef maps for masking"; exit 1; fi

runmasks=()
for SUB in "${SUBJECTS[@]}"; do
    for TASK in "${TASKS[@]}"; do
        for SES in ${SESS[${SUB} ${TASK}]}; do
            D=$(stats_dset "${SUB}" "${SES}" "${TASK}")
            rmask=${MASKDIR}/${SUB}_${SES}_task-${TASK}_modelled.nii.gz
            rm -f "${rmask}" "${rmask%.nii.gz}_native.nii.gz"
            ${AFNI} "3dcalc -a ${D}'[${SLOPE_COEF}]' -b ${D}'[${MEAN_COEF}]' \
                         -expr 'step(abs(a))*step(abs(b))' -datum byte \
                         -prefix ${rmask%.nii.gz}_native.nii.gz
                     3dresample -master ${REF} -rmode NN -prefix ${rmask} -input ${rmask%.nii.gz}_native.nii.gz"
            rm -f "${rmask%.nii.gz}_native.nii.gz"
            [ -f "${rmask}" ] && runmasks+=("${rmask}")
        done
    done
done
if [ "${#runmasks[@]}" -eq 0 ]; then echo "ERROR: no run masks could be built for the group mask"; exit 1; fi

rm -f "${GROUPMASK}"
echo "  intersecting ${#runmasks[@]} run masks (non-zero betas in every run)"
${AFNI} "3dmask_tool -input ${runmasks[*]} -frac 1.0 -prefix ${GROUPMASK}"
echo "  group mask: $(${AFNI} "3dBrickStat -count -non-zero ${GROUPMASK}" 2>/dev/null | tail -n1 | awk '{print $1}') voxels"
${AFNI} "3drefit -space MNI ${GROUPMASK}"

echo "  grid consistency check (want 1s):"
${AFNI} "3dinfo -same_grid ${GROUPMASK} ${INPUTS}/*_slope_Coef.nii.gz" || true

# ------------------------------------------------------------------------------
# STEP 3. 3dMEMA one-sample test, per task (slope = primary, mean = sanity).
#
# 3dMEMA ("Mixed Effects Meta Analysis") tests whether the group effect is
# different from zero. Unlike a plain t-test it is given each participant's
# beta AND t-statistic, so it knows how precise each participant's estimate is
# and lets less precise participants count less.
#   -set NAME  sub beta t  sub beta t ...   one entry per participant
#   -mask                                   only test voxels in the group mask
#   -missing_data 0                         a value of 0 at a voxel means
#                                           "no data for this participant here"
# All other 3dMEMA options are at their defaults (in AFNI 26.1.02's 3dMEMA.R,
# that includes the Hartung-Knapp adjustment of the output t-statistic).
# Output sub-brick [0] = group effect, [1] = its t-statistic.
# Skipped when fewer than 2 participants have the map.
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 3: 3dMEMA one-sample per task"
echo "=============================================================="

for TASK in "${TASKS[@]}"; do
    for eff in slope mean; do
        setarg=""
        n=0
        for SUB in $(subs_with "${TASK}" "${eff}"); do
            c=${INPUTS}/${SUB}_task-${TASK}_${eff}_Coef.nii.gz
            t=${INPUTS}/${SUB}_task-${TASK}_${eff}_T.nii.gz
            setarg="${setarg} ${SUB} ${c} ${t}"
            n=$((n+1))
        done
        if [ "${n}" -lt 2 ]; then
            echo "  skip MEMA ${TASK} ${eff}: only ${n} subject(s)"
            continue
        fi
        out=${MEMADIR}/mema_${TASK}_${eff}
        echo "  MEMA ${TASK} ${eff}  (n=${n})"
        rm -f "${out}"+*.HEAD "${out}"+*.BRIK*
        ${AFNI} "3dMEMA -prefix ${out} -mask ${GROUPMASK} -missing_data 0 \
            -set ${TASK}_${eff} ${setarg}"
    done
done

# ------------------------------------------------------------------------------
# STEP 4. 3dMEMA between-task contrasts of the rating slope (X minus Y).
#
# For each participant who has a slope map for BOTH tasks:
#   bDelta = bX - bY                        difference of the slopes
#   SE^2   = (b / t)^2 for each task        (as in STEP 1)
#   tDelta = bDelta / sqrt(SE_X^2 + SE_Y^2) its t-statistic
#   (where a task's t = 0, its SE^2 term is replaced by 1)
# The difference maps then go into a one-sample 3dMEMA, as in STEP 3.
# Participants here are chosen by whether their two slope files exist in
# inputs/ (not through subs_with).
# Output: inputs/contrast_X_vs_Y/<sub>_bDelta / _tDelta, mema/mema_X_vs_Y_slope.
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 4: 3dMEMA between-task slope contrasts"
echo "=============================================================="

for pair in "${BETWEEN_TASK_PAIRS[@]}"; do
    X=${pair%%:*}; Y=${pair##*:}
    cdir=${INPUTS}/contrast_${X}_vs_${Y}
    mkdir -p "${cdir}"
    setarg=""; n=0
    for SUB in "${SUBJECTS[@]}"; do
        cX=${INPUTS}/${SUB}_task-${X}_slope_Coef.nii.gz
        tX=${INPUTS}/${SUB}_task-${X}_slope_T.nii.gz
        cY=${INPUTS}/${SUB}_task-${Y}_slope_Coef.nii.gz
        tY=${INPUTS}/${SUB}_task-${Y}_slope_T.nii.gz
        { [ -f "${cX}" ] && [ -f "${cY}" ]; } || continue
        bD=${cdir}/${SUB}_bDelta.nii.gz
        tD=${cdir}/${SUB}_tDelta.nii.gz
        rm -f "${bD}" "${tD}"
        # a = bX, b = tX, c = bY, d = tY.
        # (a/(b+equals(b,0)))^2*step(abs(b)) = seX^2 where tX != 0, else 0;
        # "+ equals(b,0) + equals(d,0)" adds 1 for each task whose t is 0.
        ${AFNI} "3dcalc \
            -a ${cX} -b ${tX} -c ${cY} -d ${tY} \
            -expr '(a-c)' -prefix ${bD}"
        ${AFNI} "3dcalc \
            -a ${cX} -b ${tX} -c ${cY} -d ${tY} \
            -expr '(a-c) / sqrt( (a/(b+equals(b,0)))^2*step(abs(b)) + (c/(d+equals(d,0)))^2*step(abs(d)) + equals(b,0) + equals(d,0) )' \
            -prefix ${tD}"
        setarg="${setarg} ${SUB} ${bD} ${tD}"
        n=$((n+1))
    done
    if [ "${n}" -lt 2 ]; then
        echo "  skip MEMA ${X}-vs-${Y}: only ${n} subject(s) with both tasks"
        continue
    fi
    out=${MEMADIR}/mema_${X}_vs_${Y}_slope
    echo "  MEMA ${X} vs ${Y} slope  (n=${n})"
    rm -f "${out}"+*.HEAD "${out}"+*.BRIK*
    ${AFNI} "3dMEMA -prefix ${out} -mask ${GROUPMASK} -missing_data 0 \
        -set ${X}_minus_${Y} ${setarg}"
done

# ------------------------------------------------------------------------------
# STEP 5. 3dttest++ versions of STEPs 3 and 4.
#
# 3dttest++ is an ordinary t-test on the betas only (it ignores each
# participant's precision). It is run as a cross-check of the 3dMEMA results.
#   (a) one-sample t-test per task and effect (-setA = one map per participant)
#   (b) paired t-test of the slope between tasks (-paired, -setA = task X,
#       -setB = task Y, same participants in the same order)
# -Clustsim NTHR: after the t-test, 3dttest++ builds its own cluster-size
#   tables by randomly sign-flipping the residuals (about 10000 times). The
#   AFNI help says this needs at least 14 datasets in total (setA + setB), so it
#   is only added when n >= 14 (one-sample) or 2n >= 14 (paired). With
#   -Clustsim on, the output statistics are z-scores instead of t-values.
#   Without it, use the STEP 6 table instead.
# The command first moves into ttest/ and uses a relative -prefix (the original
# author found -Clustsim needs a relative prefix).
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 5: 3dttest++ (with -Clustsim)"
echo "=============================================================="

# (a) One-sample per task
for TASK in "${TASKS[@]}"; do
    for eff in slope mean; do
        setarg=""; n=0
        for SUB in $(subs_with "${TASK}" "${eff}"); do
            setarg="${setarg} ${SUB} ${INPUTS}/${SUB}_task-${TASK}_${eff}_Coef.nii.gz"
            n=$((n+1))
        done
        if [ "${n}" -lt 2 ]; then
            echo "  skip ttest ${TASK} ${eff}: only ${n} subject(s)"
            continue
        fi
        out=${TTESTDIR}/ttest_${TASK}_${eff}
        echo "  ttest ${TASK} ${eff}  (n=${n})"
        rm -f "${out}"+*.HEAD "${out}"+*.BRIK*
        # 3dttest++ -Clustsim needs >=14 samples AND a RELATIVE -prefix. Below
        # that, run without it and use STEP 6's residual-based table instead.
        CS=""
        if [ "${n}" -ge 14 ]; then CS="-Clustsim ${NTHR} -prefix_clustsim cs_${TASK}_${eff}"
        else echo "    (n=${n}<14: skipping built-in -Clustsim; correct via STEP 6 table)"; fi
        ${AFNI} "cd ${TTESTDIR}; 3dttest++ -prefix ttest_${TASK}_${eff} -mask ${GROUPMASK} \
            -setA ${TASK}_${eff} ${setarg} ${CS}"
    done
done

# (b) Paired between-task (slope): participants with BOTH tasks, same order in
#     setA and setB. As in STEP 4, chosen by whether the two files exist.
for pair in "${BETWEEN_TASK_PAIRS[@]}"; do
    X=${pair%%:*}; Y=${pair##*:}
    setA=""; setB=""; n=0
    for SUB in "${SUBJECTS[@]}"; do
        cX=${INPUTS}/${SUB}_task-${X}_slope_Coef.nii.gz
        cY=${INPUTS}/${SUB}_task-${Y}_slope_Coef.nii.gz
        { [ -f "${cX}" ] && [ -f "${cY}" ]; } || continue
        setA="${setA} ${SUB} ${cX}"
        setB="${setB} ${SUB} ${cY}"
        n=$((n+1))
    done
    if [ "${n}" -lt 2 ]; then
        echo "  skip paired ttest ${X}-vs-${Y}: only ${n} subject(s) with both"
        continue
    fi
    out=${TTESTDIR}/ttest_${X}_vs_${Y}_slope
    echo "  paired ttest ${X} vs ${Y} slope  (n=${n})"
    rm -f "${out}"+*.HEAD "${out}"+*.BRIK*
    # -Clustsim needs setA+setB >= 14 values (i.e. n>=7 pairs) AND a relative prefix.
    CS=""
    if [ "$((2*n))" -ge 14 ]; then CS="-Clustsim ${NTHR} -prefix_clustsim cs_${X}_vs_${Y}"
    else echo "    (2n=$((2*n))<14: skipping built-in -Clustsim; correct via STEP 6 table)"; fi
    ${AFNI} "cd ${TTESTDIR}; 3dttest++ -prefix ttest_${X}_vs_${Y}_slope -mask ${GROUPMASK} -paired \
        -setA ${X} ${setA} \
        -setB ${Y} ${setB} ${CS}"
done

# ------------------------------------------------------------------------------
# STEP 6. Cluster-size table for the group maps (3dFWHMx -> 3dClustSim).
#
# 3dMEMA has no built-in cluster simulation, and 3dttest++ -Clustsim needs at
# least 14 datasets, so this table is the general-purpose correction. It
# depends only on the mask, the smoothness of the noise and CLUST_PTHR, not on
# the number of participants, so it can be used for the 3dMEMA maps and for
# small-n 3dttest++ maps.
#   (a) 3dFWHMx -detrend -acf: for every first-level REML residual time series
#       (<PREFIX>_errts_REML+tlrc) that is not excluded, estimate how the noise
#       is correlated across space, inside the group mask. Only the uncensored
#       volumes are used (from <PREFIX>_motion_censor.1D; censored volumes are 0
#       in errts), and each voxel's trend is removed first, as afni_proc.py
#       does. Without these, one run (sub-19 ses-2) gave FWHM 5.4 mm instead of
#       8.1 mm. The fit is
#       ACF(r) = a*exp(-r^2/(2*b^2)) + (1-a)*exp(-r/c); the last line of the
#       output holds a, b, c (and an effective FWHM). One line per run goes into
#       clustsim/acf_params.txt (a b c FWHM run). (acf_est.1D, the empirical curve, is
#       overwritten by each run.)
#   (b) a, b, c are averaged over runs.
#   (c) 3dClustSim simulates noise with that smoothness inside the group mask
#       and reports the minimum cluster size (in voxels) that occurs by chance
#       with probability alpha = 0.05 (-athr) at voxelwise p = CLUST_PTHR
#       (-pthr). Output: clustsim/mema_clustsim.NN{1,2,3}_{1sided,2sided,
#       bisided}.1D (NN = which neighbouring voxels count as touching).
# Apply a table with 3dClusterize (an example command is printed below).
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 6: 3dFWHMx -> 3dClustSim for MEMA cluster thresholds"
echo "=============================================================="

acf_file=${CLUSTDIR}/acf_params.txt
: > "${acf_file}"
for ERRTS in "${AFNI_FIRSTLVL}"/sub-*/ses-*/task-*/*_errts_REML+tlrc.HEAD; do
    [ -e "${ERRTS}" ] || continue
    _t=$(basename "$(dirname "${ERRTS}")"); _s=$(basename "$(dirname "$(dirname "${ERRTS}")")")
    _p=$(basename "$(dirname "$(dirname "$(dirname "${ERRTS}")")")")
    is_excluded "${_p}" "${_s}" "${_t#task-}" >/dev/null && continue
    dset=${ERRTS%.HEAD}
    # Uncensored volumes only (censored ones are 0 in errts), as afni_proc.py does.
    censor=${dset%_errts_REML+tlrc}_motion_censor.1D
    if [ -f "${censor}" ]; then
        keep=$(${AFNI} "1d_tool.py -infile ${censor} -show_trs_uncensored encoded" 2>/dev/null | tail -n1)
    else
        echo "  WARNING: no ${censor##*/}; using every volume of ${dset##*/}"
        keep=""
    fi
    sel=""; [ -n "${keep}" ] && sel="'[${keep}]'"
    line=$(${AFNI} "3dFWHMx -detrend -acf ${CLUSTDIR}/acf_est.1D -mask ${GROUPMASK} -input ${dset}${sel} 2>/dev/null" | tail -n1)
    echo "${line} ${dset##*/}" >> "${acf_file}"
done

if [ -s "${acf_file}" ]; then
    # NOTE: printf MUST emit a trailing \n, else `read` returns nonzero and
    # `set -e` aborts the script even though a/b/c were assigned correctly.
    read a b c < <(awk '{a+=$1;b+=$2;c+=$3;n++} END{if(n)printf "%.6f %.6f %.6f\n", a/n,b/n,c/n}' "${acf_file}")
    echo "  mean ACF params: a=${a} b=${b} c=${c}"
    ${AFNI} "cd ${CLUSTDIR}; 3dClustSim -mask ${GROUPMASK} -acf ${a} ${b} ${c} \
        -pthr ${CLUST_PTHR} -athr 0.05 -prefix mema_clustsim"
    echo "  MEMA cluster-size table -> ${CLUSTDIR}/mema_clustsim.* "
    echo "  Apply with e.g.:"
    echo "    3dClusterize -inset ${MEMADIR}/mema_scripture_slope+tlrc -ithr 1 -idat 0 \\"
    echo "      -mask ${GROUPMASK} -NN 1 -1sided RIGHT_TAIL p=${CLUST_PTHR} -clust_nvox <from table>"
else
    echo "  WARNING: no residual datasets found; skipped 3dClustSim for MEMA."
fi

# ------------------------------------------------------------------------------
# STEP 7. Copying to archive is NOT done here (changed 2026-08-21).
#
# This step used to copy (cp -a) the group maps into GROUP_FINAL on archive.
# Two problems with doing that from inside this job:
#   1. It is archive I/O from a batch job. BYU RC: "Archive storage should
#      generally NOT be used directly from batch jobs", and archive is migrating
#      to a much slower backing store.
#   2. It wrote an untarred tree onto a tier with a 1 M inode cap, every single
#      run, silently accumulating near-duplicate copies.
# promote_to_archive.sh already does this correctly: from a login node, filtered
# to the products worth keeping, and tarred. Run it when you park the project.
# This step now only prints that reminder.
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 7: skipped — retention now belongs to promote_to_archive.sh"
echo "=============================================================="
echo "Group results stay on SCRATCH: ${GROUP_OUT}"
echo "SCRATCH IS PURGED AFTER 12 WEEKS UNUSED. Before parking this project run:"
echo "    bash ${PROJ_HOME}/scripts/promote_to_archive.sh"

echo "=============================================================="
echo "Group analysis complete."
echo "  n subjects discovered: ${#SUBJECTS[@]}  (${SUBJECTS[*]})"
echo "  scratch (working) : ${GROUP_OUT}"
echo "  archive (kept)    : nothing yet - run promote_to_archive.sh to persist"
echo "  NOTE: with only a few subjects this is UNDERPOWERED - expect no"
echo "        surviving clusters. Re-run as sub-02/sub-04 and more sessions"
echo "        finish; new first-level output is picked up automatically."
echo "=============================================================="
