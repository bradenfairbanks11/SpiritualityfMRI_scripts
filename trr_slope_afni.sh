#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   4 cores / 16 GB / 12 h : one brms model per ROI, 4 chains run in parallel
#                            (one core each). Not yet timed on real data; about
#                            70 ROIs over the three tasks at a few minutes each
#                            fit comfortably. Submitted as an array (see HOW TO
#                            RUN), each task gets its own job.
#   --chdir                : run inside scratch, so any stray relative-path file
#                            lands there instead of in the home directory.
# ==============================================================================
#SBATCH --job-name=spirituality_trr_slope
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# ==============================================================================
# trr_slope_afni.sh — TEST-RETEST RELIABILITY OF THE RATING SLOPE (REGION LEVEL)
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   Asks whether each participant's BOLD-vs-rating SLOPE in an ROI (how much
#   the single-trial beta changes per step of the 1-4 "felt the Spirit"
#   rating) is the same in session 1 and session 2.
#
#   For each task and each ROI, trr_slope.R fits one Bayesian multilevel model
#   (brms) to the single-trial ROI means:
#     Y ~ 0 + sess + sess:rc + (0 + sess | Subj) + (0 + sess:rc | Subj)
#   (rc = rating - 2.5). Each participant gets a ses-1 and a ses-2 slope; the
#   model's correlation between them is the slope TRR. This is AFNI's TRR
#   model (Chen et al. 2021) with the two-level condition replaced by the
#   continuous rating. AFNI's own TRR / 3dLMEr -TRR take only two conditions,
#   so they cannot do this, and reliability_models_afni.sh's TRR-con uses
#   stimulus (yield) labels, not ratings.
#
#   A new, separate script: it changes no other script and writes only to its
#   own folder, results_slope/.
#
# WHERE IT SITS IN THE PIPELINE
#   Before: trr_export_afni.sh   (roi_task-<T>.tsv: one row per trial x ROI)
#   After : nothing; the results are end products. Compare with
#           reliability_models_afni.sh trr_roi (overall-response TRR) and
#           reproducibility_afni.sh 3dICC on rating_slope (first-level ICC).
#
# STEPS
#   1. For each task, run trr_slope.R on roi_task-<T>.tsv.
#   2. Join the per-task summaries into trr_slope_summary.tsv.
#
# INPUTS
#   <TRR_OUT>/roi_task-<T>.tsv   Subj sess trial rating cond stimulus yield faith ROI Y
#   <SCRIPTS>/trr_slope.R
#
# OUTPUTS, in <TRR_OUT>/results_slope/
#   trr_slope_summary.tsv                 every task's rows (columns below)
#   trr_slope_summary_task-<T>.tsv        one task's rows
#   rating_counts_task-<T>.tsv            trials per rating, per Subj x sess
#   trr_slope_task-<T>_dropped.txt        participants left out (< 2 distinct
#                                         ratings in a session)
#   TRR_slope_task-<T>_<ROI>.txt / .rds   model summary / fitted model
#   trr_slope_task-<T>.log                trr_slope.R's messages
#   Summary columns: task ROI n_subj n_trials, then for slope and intercept
#   TRR each as _cor (model correlation; primary) and _afni (TRR.R style, for
#   comparison with trr_roi) with _mode/_lo/_hi (posterior mode, 95% HDI),
#   then mean_slope (+ 95% interval), max_rhat, n_divergent.
#
# HOW TO RUN (from /home/bradenf4/spirituality_fmri/scripts; copy trr_slope.R
# there too)
#   One job per task (array task k runs TASKS[k]):
#     sbatch --array=0-2 trr_slope_afni.sh
#   then, when all three are done, join the summaries:
#     bash trr_slope_afni.sh join
#   Or every task in one job, joined at the end:
#     sbatch trr_slope_afni.sh
#   After a fresh trr_export_afni.sh job:  sbatch --dependency=afterok:<job> ...
#   Re-running overwrites that task's outputs.
#
# METHODS SUMMARY (the exact values used below)
#   Data          : single-trial REML betas (percent signal change) from
#                   first_level_IM_afni.sh, averaged within each ROI (clusters
#                   of the group 3dMEMA mean_response map; trr_export_afni.sh).
#                   Only rated trials; participants with both sessions.
#   Model         : brms, Gaussian,
#                   Y ~ 0 + sess + sess:rc + (0+sess|Subj) + (0+sess:rc|Subj),
#                   rc = rating - 2.5; brms default priors (LKJ(1) on
#                   correlations); 4 chains x 1000 iterations (500 warm-up),
#                   adapt_delta 0.99, max_treedepth 15, seed 1234 (as AFNI TRR).
#   Inclusion     : >= 2 distinct rating values in EACH session.
#   Slope TRR     : posterior of cor(ses-1 slope, ses-2 slope), reported as
#                   mode and 95% highest density interval.
# ==============================================================================

set -o pipefail    # a pipeline fails if ANY command in it fails
                   # (note: no "set -e" -- one failed task does not stop the rest)


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
TRR_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_trr           # [shared] output
LOGS=/home/bradenf4/spirituality_fmri/logs                                             # [shared]
# The scripts folder (for trr_slope.R). A fixed path, because under sbatch
# this script runs from a copy in Slurm's spool folder.
SCRIPTS=/home/bradenf4/spirituality_fmri/scripts

RESULTS=${TRR_OUT}/results_slope

# --- 1b. Software -------------------------------------------------------------
# R runs inside the AFNI Apptainer image, which carries brms and coda. The
# image only sees folders that are explicitly "bound" into it, so all three
# storage areas are bound. Compute nodes have no internet, so nothing can be
# installed from inside a job.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]
R_PKGS='c("brms","coda")'

# --- 1c. Analysis parameters --------------------------------------------------
# Tasks, by their BIDS labels. With --array, array task k runs TASKS[k].
TASKS=(scripture FHS architecture)                                                     # [shared]

# Sampler settings, the same as reliability_models_afni.sh's region TRR.
SLOPE_CHAINS=4
SLOPE_ITER=1000


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# join_summaries  ->  stack every trr_slope_summary_task-<T>.tsv (one header)
join_summaries() {
    local out=${RESULTS}/trr_slope_summary.tsv first=1 f
    : > "${out}"
    for f in "${RESULTS}"/trr_slope_summary_task-*.tsv; do
        [ -f "${f}" ] || continue
        if [ ${first} -eq 1 ]; then cat "${f}" >> "${out}"; first=0
        else tail -n +2 "${f}" >> "${out}"; fi
    done
    echo "  summary: ${out}"
}


# ==============================================================================
# SECTION 3. SETUP AND CHECKS
# ==============================================================================

mkdir -p "${LOGS}" "${RESULTS}"

# "join" only stacks the per-task summaries (no R, no image needed).
if [ "${1:-}" = "join" ]; then
    join_summaries
    exit 0
fi

module load "${APPTAINER_MODULE}"

# Every R command below runs as:  ${AFNI_R} "<commands>"
# --cleanenv starts the image with an empty environment, so a host setting
# (e.g. a stray R_LIBS) cannot override the packages inside the image.
AFNI_R="apptainer exec --cleanenv ${BIND_ARGS} ${AFNI_SIF} bash -c"

if [ ! -f "${SCRIPTS}/trr_slope.R" ]; then
    echo "ERROR: no ${SCRIPTS}/trr_slope.R (copy it next to this script)" >&2; exit 1
fi

# Check that the image has the R packages (R exits with the number missing).
if ! ${AFNI_R} "Rscript -e 'q(status=length(setdiff(${R_PKGS}, rownames(installed.packages()))))'"; then
    echo "ERROR: R packages missing from ${AFNI_SIF} (need: brms coda)." >&2; exit 1
fi

# With --array, run only one task: TASKS[SLURM_ARRAY_TASK_ID].
ARRAY=0
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    if [ -z "${TASKS[${SLURM_ARRAY_TASK_ID}]:-}" ]; then
        echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} has no task (use --array=0-$(( ${#TASKS[@]} - 1 )))" >&2
        exit 1
    fi
    TASKS=("${TASKS[${SLURM_ARRAY_TASK_ID}]}")
    ARRAY=1
    echo "Array task ${SLURM_ARRAY_TASK_ID}: task-${TASKS[0]} only"
fi


# ==============================================================================
# SECTION 4. THE STEPS
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1. Fit the slope-TRR model for every ROI of each task.
#
# trr_slope.R reads roi_task-<T>.tsv, drops participants without two distinct
# ratings in each session, fits the model once per ROI (compiling it once per
# task) and writes that task's summary. A task with no ROI table (no group
# cluster survived, or trr_export_afni.sh not run) is skipped.
# ------------------------------------------------------------------------------
echo "=============================================================="
echo "STEP 1: slope TRR (region level)"
echo "=============================================================="
for TASK in "${TASKS[@]}"; do
    RT=${TRR_OUT}/roi_task-${TASK}.tsv
    if [ ! -f "${RT}" ]; then
        echo "skip ${TASK}: no ${RT} (no group ROI survived, or trr_export_afni.sh not run)"
        continue
    fi
    rm -f "${RESULTS}"/TRR_slope_task-"${TASK}"_* "${RESULTS}/trr_slope_summary_task-${TASK}.tsv"
    if ! ${AFNI_R} "Rscript ${SCRIPTS}/trr_slope.R ${RT} ${TASK} ${RESULTS} ${SLOPE_CHAINS} ${SLOPE_ITER}" \
            2>&1 | tee "${RESULTS}/trr_slope_task-${TASK}.log"; then
        echo "  WARN: trr_slope.R failed for task-${TASK}; see ${RESULTS}/trr_slope_task-${TASK}.log" >&2
    fi
done

# ------------------------------------------------------------------------------
# STEP 2. Join the per-task summaries (not in an array job: the other tasks may
# still be running -- run  bash trr_slope_afni.sh join  when all are done).
# ------------------------------------------------------------------------------
if [ ${ARRAY} -eq 0 ]; then
    join_summaries
else
    echo "  array job: run  bash trr_slope_afni.sh join  after every task finishes"
fi

echo "=============================================================="
echo "slope TRR done: ${TASKS[*]}"
echo "=============================================================="
