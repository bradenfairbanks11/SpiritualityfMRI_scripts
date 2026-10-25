#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   16 cores / 64 GB / 24 h : 3dLMEr -TRR fits one mixed model per voxel and
#                             splits the voxels across all 16 cores (-jobs). It
#                             has not been timed on real data yet; 24 h is a
#                             generous guess for one model. Submitted as an
#                             array (see HOW TO RUN), each model gets its own
#                             job with these resources.
#   --chdir                 : run inside scratch, so any stray relative-path file
#                             lands there instead of in the home directory.
# ==============================================================================
#SBATCH --job-name=spirituality_reliability
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# ==============================================================================
# reliability_models_afni.sh — TRIAL-LEVEL TEST-RETEST RELIABILITY MODELS (R)
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   Runs the R-based reliability models on the files that trr_export_afni.sh
#   and reproducibility_afni.sh wrote. All of them run in R inside the AFNI
#   image, which carries R and every package they need.
#     trr_vox  Whole-brain trial-level test-retest reliability (TRR) with
#              3dLMEr -TRR (Chen et al. 2021, NeuroImage 245:118647). The model
#              needs exactly two conditions; they are defined by the STIMULUS
#              each trial showed (labels added by trr_export_afni.sh from
#              stimulus_conditions.tsv):
#                yield: high- vs low-yield stimuli    (all three tasks)
#                faith: LDS vs other-faith images     (architecture only)
#              One model per task x contrast (4 models). Each gives TRR-ave
#              (reliability of the average of the two conditions) and TRR-con
#              (reliability of the difference between them).
#     trr_roi  Region-level TRR with AFNI's TRR program (Bayesian multilevel
#              model, brms), one model per task x ROI, single condition: the
#              reliability of each participant's overall response in that ROI.
#     dice     Is each participant's ses-1 vs ses-2 map overlap (Dice) larger
#              than the overlap between different participants? Runs
#              dice_within_vs_between.R (Gorgolewski et al. 2013) on
#              reproducibility_afni.sh's dice_within_between.tsv.
#   (The voxelwise 3dICC map is NOT here: reproducibility_afni.sh runs it on
#   the cluster itself, in its STEP 11.)
#
#   Replaces local_r_steps.sh, which ran these steps on the Mac while the old
#   cluster AFNI image could not run R programs.
#
# WHERE IT SITS IN THE PIPELINE
#   Before: trr_export_afni.sh       (trial tables, trial files, masks, ROI tables)
#           reproducibility_afni.sh  (dice_within_between.tsv)
#   After : nothing; the results are the end products.
#
# STEPS
#   1. trr_vox: for each model, build its data table, run 3dLMEr -TRR, and
#      summarise the two reliability maps (medians; share of 0/+-1 voxels).
#   2. trr_roi: for each task and ROI, run TRR and record its mode and 95%
#      highest density interval.
#   3. dice:    run dice_within_vs_between.R.
#
# INPUTS
#   <TRR_OUT>/trials_task-<T>.tsv   Subj sess trial rating cond stimulus yield faith InputFile
#   <TRR_OUT>/inputs/*.nii.gz       one beta map per trial (InputFile, relative)
#   <TRR_OUT>/mask_task-<T>.nii.gz  voxels with a non-zero beta for every trial
#   <TRR_OUT>/rois_task-<T>.nii.gz, roi_task-<T>.tsv   ROIs and per-trial ROI means
#   <REPRO_OUT>/dice_within_between.tsv
#   <SCRIPTS>/dice_within_vs_between.R
#
# OUTPUTS, in <TRR_OUT>/results/ (dice: printed to the job log)
#   lmer_trr_task-<T>_con-<C>.tbl          the 3dLMEr data table of that model
#   lmer_trr_task-<T>_con-<C>_dropped.txt  participants left out (too few trials)
#   LMEr_TRR_task-<T>_con-<C>+tlrc         3dLMEr output (TRR-ave, TRR-con, ...)
#                                          (and its own _log.txt)
#   LMEr_TRR_task-<T>_con-<C>.log          3dLMEr's own messages
#   trr_vox_summary_task-<T>_con-<C>.tsv   that model's summary rows
#   trr_vox_summary.tsv                    all models' summary rows together
#   trr_task-<T>_<ROI>.tbl, TRR_task-<T>_<ROI>.{txt,RData,log}, *-ave.pdf
#   trr_roi_summary.tsv                    task ROI n_subj n_trials TRR_mode HDI95_lo HDI95_hi
#
# HOW TO RUN (from /home/bradenf4/spirituality_fmri/scripts)
#   The four 3dLMEr models are slow, so give each its own job (an array):
#     sbatch --array=0-3 reliability_models_afni.sh trr_vox
#   (array task k runs model k of TRR_VOX_RUNS, SECTION 1c). Then, in one job:
#     sbatch reliability_models_afni.sh trr_roi dice
#   With no step names, every step runs in one job, models one after another:
#     sbatch reliability_models_afni.sh
#   Re-running overwrites that step's outputs.
#
# METHODS SUMMARY (the exact values used below)
#   Whole-brain TRR : 3dLMEr -TRR, model
#                       0+sess+cond:sess+(0+sess|Subj)+(0+cond:sess|Subj)
#                     (the 3dLMEr -TRR help example), cond coded +0.5 / -0.5
#                     and given with -qVars cond; data = single-trial REML
#                     betas (percent signal change) from first_level_IM_afni.sh;
#                     mask = voxels with a non-zero beta for every trial of
#                     every included run. No -bounds (the help calls outlier
#                     bounds arbitrary).
#                     Contrasts: yield high (+0.5) vs low (-0.5), every task;
#                     faith LDS (+0.5) vs other (-0.5), architecture.
#                     A participant enters a model only with >= 2 trials of
#                     EACH condition in EACH session.
#   Region TRR      : AFNI TRR, single condition, 4 chains x 1000 iterations,
#                     default (Gaussian) likelihood; Y = mean single-trial beta
#                     in the ROI; reported as posterior mode and 95% HDI.
#   Dice test       : see dice_within_vs_between.R.
# ==============================================================================

set -o pipefail    # a pipeline fails if ANY command in it fails
                   # (note: no "set -e" -- one failed model does not stop the rest)


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
REPRO_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_reproducibility  # [shared] output
LOGS=/home/bradenf4/spirituality_fmri/logs                                             # [shared]
# The scripts folder (for dice_within_vs_between.R). A fixed path, because
# under sbatch this script runs from a copy in Slurm's spool folder.
SCRIPTS=/home/bradenf4/spirituality_fmri/scripts

RESULTS=${TRR_OUT}/results

# --- 1b. Software -------------------------------------------------------------
# AFNI runs inside an Apptainer image. The image only sees folders that are
# explicitly "bound" into it, so all three storage areas are bound. The image
# also carries R and the R packages used here (3dLMEr: lmerTest, phia, snow;
# TRR: brms, coda, ggplot2). Compute nodes have no internet, so nothing can be
# installed from inside a job.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]
R_PKGS='c("lmerTest","phia","snow","brms","coda","ggplot2")'

# --- 1c. Analysis parameters --------------------------------------------------
# Tasks, by their BIDS labels.
TASKS=(scripture FHS architecture)                                                     # [shared]

# Whole-brain 3dLMEr -TRR models, as task:column. The column is a column of
# trials_task-<T>.tsv; contrast_levels (SECTION 2) gives the level coded +0.5
# and the level coded -0.5. Per run there are about 12/12 (scripture), 7/5
# (FHS) and 14/16 (architecture) high/low-yield trials.
# With --array, array task k runs TRR_VOX_RUNS[k].
TRR_VOX_RUNS=(scripture:yield FHS:yield architecture:yield architecture:faith)

# A participant enters a 3dLMEr model only with at least this many trials of
# EACH condition in EACH session (the difference needs both conditions).
MIN_PER_COND=2

# Region-level TRR sampler settings. AFNI's TRR.R stores the posterior
# correlation draws in a vector pre-sized to 2000 and fills chains*iterations/2
# of them, so fewer than 2000 draws would leave zeros that drag the reported
# TRR toward 0. Keep CHAINS * ITERATIONS / 2 >= 2000 (checked in SECTION 3).
TRR_CHAINS=4
TRR_ITER=1000


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# contrast_levels column  ->  prints "plus_level minus_level"
contrast_levels() {
    case "$1" in
        yield) echo "high low" ;;
        faith) echo "lds other" ;;
        *)     echo "ERROR: unknown contrast column $1" >&2; return 1 ;;
    esac
}

# want step  ->  succeeds if that step was asked for on the command line
want() { local s; for s in "${STEPS[@]}"; do [ "${s}" = "$1" ] && return 0; done; return 1; }


# ==============================================================================
# SECTION 3. SETUP AND CHECKS
# ==============================================================================

# Steps to run: the ones named on the command line, or all of them.
STEPS=("$@")
[ ${#STEPS[@]} -eq 0 ] && STEPS=(trr_vox trr_roi dice)
for s in "${STEPS[@]}"; do
    case "${s}" in trr_vox|trr_roi|dice) ;; *) echo "ERROR: unknown step '${s}' (use trr_vox, trr_roi, dice)" >&2; exit 1 ;; esac
done

mkdir -p "${LOGS}" "${RESULTS}"

module load "${APPTAINER_MODULE}"

# Every command below runs as:  ${AFNI_R} "<commands>"
# --cleanenv starts the image with an empty environment, so a host setting
# (e.g. a stray R_LIBS) cannot override the packages inside the image. Every
# variable these commands need is filled in on the host before the command
# reaches the image; --env still applies under --cleanenv.
# AFNI_DECONFLICT=OVERWRITE lets a re-run replace its own files.
AFNI_R="apptainer exec --cleanenv --env AFNI_DECONFLICT=OVERWRITE --env AFNI_ENVIRON_WARNINGS=NO ${BIND_ARGS} ${AFNI_SIF} bash -c"

# Cores for 3dLMEr: everything Slurm gave this job.
JOBS=${SLURM_CPUS_PER_TASK:-4}

if [ $(( TRR_CHAINS * TRR_ITER / 2 )) -lt 2000 ]; then
    echo "ERROR: TRR_CHAINS*TRR_ITER/2 must be >= 2000 (see SECTION 1c)" >&2; exit 1
fi

# Check that the image has the R packages (R exits with the number missing).
if ! ${AFNI_R} "Rscript -e 'q(status=length(setdiff(${R_PKGS}, rownames(installed.packages()))))'"; then
    echo "ERROR: R packages missing from ${AFNI_SIF} (need: lmerTest phia snow brms coda ggplot2)." >&2
    echo "       Check AFNI_SIF (SECTION 1b) and container_request.md." >&2
    exit 1
fi

# With --array, run only one 3dLMEr model: TRR_VOX_RUNS[SLURM_ARRAY_TASK_ID].
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    if [ -z "${TRR_VOX_RUNS[${SLURM_ARRAY_TASK_ID}]:-}" ]; then
        echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} has no model (TRR_VOX_RUNS has ${#TRR_VOX_RUNS[@]}: use --array=0-$(( ${#TRR_VOX_RUNS[@]} - 1 )))" >&2
        exit 1
    fi
    TRR_VOX_RUNS=("${TRR_VOX_RUNS[${SLURM_ARRAY_TASK_ID}]}")
    echo "Array task ${SLURM_ARRAY_TASK_ID}: model ${TRR_VOX_RUNS[0]} only"
fi


# ==============================================================================
# SECTION 4. THE STEPS
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1. trr_vox -- 3dLMEr -TRR, whole brain, one model per task x contrast.
#
# For each model:
#   (a) Data table. From trials_task-<T>.tsv (columns found by header name):
#       Subj sess cond InputFile, with cond = +0.5 for the PLUS level of the
#       contrast column and -0.5 for the MINUS level (the 3dLMEr -TRR help's
#       coding, given with -qVars cond). A participant with fewer than
#       MIN_PER_COND trials of either level in either session is left out and
#       listed in _dropped.txt. A model with fewer than 3 participants is
#       skipped.
#   (b) 3dLMEr -TRR with the model from the 3dLMEr help (see METHODS SUMMARY),
#       inside the task's mask, on JOBS cores. Its messages go to the .log.
#   (c) Summary: for TRR-ave and TRR-con, the median inside the task mask and
#       inside the ROIs (if trr_export_afni.sh made any), plus the share of mask
#       voxels that are exactly 0, 1 or -1 -- the 3dLMEr help says numerical
#       failures show up as many such values. The sub-brick is found by its
#       label with 3dinfo -label2index, because AFNI may read the '-' in
#       'TRR-ave' inside a [label] selector as a range.
#       Each model writes its own summary file; trr_vox_summary.tsv is then
#       rebuilt from all of them (so parallel array tasks do not collide).
# ------------------------------------------------------------------------------
if want trr_vox; then
    echo "=============================================================="
    echo "STEP 1: 3dLMEr -TRR (whole brain; yield = high vs low, faith = LDS vs other)"
    echo "=============================================================="
    for RUNSPEC in "${TRR_VOX_RUNS[@]}"; do
        TASK=${RUNSPEC%%:*}; CON=${RUNSPEC#*:}
        LEVELS=$(contrast_levels "${CON}") || continue
        read -r PLUS MINUS <<< "${LEVELS}"
        TABLE=${TRR_OUT}/trials_task-${TASK}.tsv
        TMASK=${TRR_OUT}/mask_task-${TASK}.nii.gz
        if [ ! -f "${TABLE}" ] || [ ! -f "${TMASK}" ]; then
            echo "skip ${TASK} ${CON}: no ${TABLE} / mask (run trr_export_afni.sh first)"
            continue
        fi

        # (a) data table
        DTAB=${RESULTS}/lmer_trr_task-${TASK}_con-${CON}.tbl
        DROPPED=${RESULTS}/lmer_trr_task-${TASK}_con-${CON}_dropped.txt
        rm -f "${DROPPED}"
        awk -F'\t' -v root="${TRR_OUT}" -v minc="${MIN_PER_COND}" -v drop="${DROPPED}" \
            -v col="${CON}" -v plus="${PLUS}" -v minus="${MINUS}" '
            FNR==1  { for (i = 1; i <= NF; i++) h[$i] = i
                      if (!(col in h) || !("InputFile" in h)) { print "ERROR: no " col " / InputFile column" > "/dev/stderr"; exit 1 }
                      c = h[col]; f = h["InputFile"]
                      if (NR != FNR) print "Subj sess cond InputFile"
                      next }
            NR==FNR { cnt[$1 SUBSEP $2 SUBSEP $c]++; next }
            {
                ok = (cnt[$1,"ses1",plus] >= minc && cnt[$1,"ses1",minus] >= minc &&
                      cnt[$1,"ses2",plus] >= minc && cnt[$1,"ses2",minus] >= minc)
                if (!ok) { bad[$1] = 1; next }
                if ($c != plus && $c != minus) next
                printf "%s %s %s %s/%s\n", $1, $2, ($c == plus ? "0.5" : "-0.5"), root, $f
            }
            END { for (s in bad) print s > drop }
        ' "${TABLE}" "${TABLE}" > "${DTAB}" || { echo "  WARN: could not build ${DTAB}" >&2; continue; }
        NSUB=$(awk 'NR>1 {print $1}' "${DTAB}" | sort -u | wc -l | tr -d ' ')
        NTR=$(($(wc -l < "${DTAB}") - 1))
        if [ -s "${DROPPED}" ]; then
            echo "  ${TASK} ${CON}: left out (< ${MIN_PER_COND} trials of a condition in a session): $(sort "${DROPPED}" | tr '\n' ' ')"
        fi
        if [ "${NSUB}" -lt 3 ]; then
            echo "skip ${TASK} ${CON}: only ${NSUB} participant(s) with both conditions in both sessions"
            continue
        fi

        # (b) 3dLMEr -TRR
        OUTP=${RESULTS}/LMEr_TRR_task-${TASK}_con-${CON}
        rm -f "${OUTP:?}"+*.HEAD "${OUTP:?}"+*.BRIK*
        echo "  ${TASK} ${CON} (+0.5 = ${PLUS}, -0.5 = ${MINUS}): n=${NSUB} participants, ${NTR} trials, ${JOBS} cores -- started $(date '+%F %T')"
        if ! ${AFNI_R} "3dLMEr -prefix ${OUTP} -TRR -jobs ${JOBS} -mask ${TMASK} \
                -qVars cond \
                -model '0+sess+cond:sess+(0+sess|Subj)+(0+cond:sess|Subj)' \
                -dataTable @${DTAB}" > "${OUTP}.log" 2>&1; then
            echo "  WARN: 3dLMEr -TRR failed for ${TASK} ${CON}; see ${OUTP}.log" >&2
            continue
        fi
        # 3dLMEr writes the output in the view of its inputs (+tlrc for these
        # MNI-space betas); find it rather than assume the view.
        ODSET=$(ls -1 "${OUTP}"+*.HEAD 2>/dev/null | head -n1); ODSET=${ODSET%.HEAD}
        if [ -z "${ODSET}" ]; then
            echo "  WARN: 3dLMEr wrote no ${OUTP}+*.HEAD; see ${OUTP}.log" >&2
            continue
        fi
        echo "  finished $(date '+%F %T') -> ${ODSET}"

        # (c) summary
        MSUM=${RESULTS}/trr_vox_summary_task-${TASK}_con-${CON}.tsv
        : > "${MSUM}"
        ROIS=${TRR_OUT}/rois_task-${TASK}.nii.gz
        for LAB in TRR-ave TRR-con; do
            IX=$(${AFNI_R} "3dinfo -label2index ${LAB} ${ODSET}" 2>/dev/null | tail -n1)
            case "${IX}" in
                ''|*[!0-9]*) IX="" ;;
            esac
            if [ -z "${IX}" ]; then
                echo "  WARN: no ${LAB} sub-brick in ${ODSET} (3dinfo -label)" >&2
                continue
            fi
            for MS in "wholebrain:${TMASK}" "rois:${ROIS}"; do
                MT=${MS%%:*}; MF=${MS#*:}
                [ -f "${MF}" ] || continue
                MED=$(${AFNI_R} "3dBrickStat -median -mask ${MF} ${ODSET}'[${IX}]'" 2>/dev/null | tail -n1 | awk '{print $NF}')
                printf '%s\t%s\t%d\t%d\t%s\t%s\t%s\n' "${TASK}" "${CON}" "${NSUB}" "${NTR}" "${LAB}" "${MT}" "${MED:-NA}" >> "${MSUM}"
            done
            DEGEN=$(${AFNI_R} "3dmaskdump -noijk -mask ${TMASK} ${ODSET}'[${IX}]'" 2>/dev/null \
                | awk '{ n++; if ($1==0 || $1==1 || $1==-1) d++ } END { if (n) printf "%.3f", d/n }')
            printf '%s\t%s\t%d\t%d\t%s\t%s\t%s\n' "${TASK}" "${CON}" "${NSUB}" "${NTR}" "${LAB}_frac_0_or_pm1" wholebrain "${DEGEN:-NA}" >> "${MSUM}"
        done
        cat "${MSUM}"
    done

    # Rebuild the combined summary from every model's file (written to a
    # temporary name first, then renamed, so a reader never sees half a file).
    VOX_SUMMARY=${RESULTS}/trr_vox_summary.tsv
    { printf 'task\tcontrast\tn_subj\tn_trials\tmetric\tmask\tmedian\n'
      cat "${RESULTS}"/trr_vox_summary_task-*_con-*.tsv 2>/dev/null || true; } > "${VOX_SUMMARY}.tmp.$$"
    mv -f "${VOX_SUMMARY}.tmp.$$" "${VOX_SUMMARY}"
    echo "  summary (all models so far): ${VOX_SUMMARY}"
fi

# ------------------------------------------------------------------------------
# STEP 2. trr_roi -- region-level TRR, single condition, one model per task x ROI.
#
# From roi_task-<T>.tsv (columns found by header name), each ROI's rows become
# a table "Subj sess Y" (Y = that trial's mean beta in the ROI), and TRR fits
# it (-subject Subj -repetition sess -Y Y; TRR_CHAINS chains of TRR_ITER
# iterations). TRR writes <prefix>.txt; its mode and 95% HDI lines are copied
# into trr_roi_summary.tsv (NA if the fit failed).
# A task with no ROI table (no group cluster survived) is skipped.
# ------------------------------------------------------------------------------
if want trr_roi; then
    echo "=============================================================="
    echo "STEP 2: TRR (region level; single condition per task)"
    echo "=============================================================="
    ROI_SUMMARY=${RESULTS}/trr_roi_summary.tsv
    printf 'task\tROI\tn_subj\tn_trials\tTRR_mode\tHDI95_lo\tHDI95_hi\n' > "${ROI_SUMMARY}"
    for TASK in "${TASKS[@]}"; do
        RT=${TRR_OUT}/roi_task-${TASK}.tsv
        if [ ! -f "${RT}" ]; then
            echo "skip ${TASK}: no ${RT} (no group ROI survived, or trr_export_afni.sh not run)"
            continue
        fi
        for ROI in $(awk -F'\t' 'NR==1 { for (i = 1; i <= NF; i++) if ($i == "ROI") c = i; next } { print $c }' "${RT}" | sort -u); do
            TAB=${RESULTS}/trr_task-${TASK}_${ROI}.tbl
            awk -F'\t' -v r="${ROI}" 'NR==1 { for (i = 1; i <= NF; i++) h[$i] = i; print "Subj sess Y"; next }
                 $h["ROI"]==r && $h["Y"]!="" { print $h["Subj"], $h["sess"], $h["Y"] }' "${RT}" > "${TAB}"
            NSUB=$(awk 'NR>1 {print $1}' "${TAB}" | sort -u | wc -l | tr -d ' ')
            NTR=$(($(wc -l < "${TAB}") - 1))
            OUTP=${RESULTS}/TRR_task-${TASK}_${ROI}
            rm -f "${OUTP:?}".txt "${OUTP:?}".RData "${OUTP:?}"-ave.pdf
            echo "  ${TASK} ${ROI}: n=${NSUB}, ${NTR} trials"
            if ! ${AFNI_R} "cd ${RESULTS} && TRR -prefix ${OUTP} -chains ${TRR_CHAINS} -iterations ${TRR_ITER} \
                    -Y Y -subject Subj -repetition sess -dataTable ${TAB}" > "${OUTP}.log" 2>&1; then
                echo "  WARN: TRR failed for ${TASK}/${ROI}; see ${OUTP}.log" >&2
                printf '%s\t%s\t%d\t%d\tNA\tNA\tNA\n' "${TASK}" "${ROI}" "${NSUB}" "${NTR}" >> "${ROI_SUMMARY}"
                continue
            fi
            MODE=$(sed -n 's/^TRR across the two repetitions: mode - \([-0-9.eE]*\);.*/\1/p' "${OUTP}.txt" 2>/dev/null | tail -n1)
            HDI=$(sed -n 's/^TRR: 95% highest density interval - \[\(.*\)\]/\1/p' "${OUTP}.txt" 2>/dev/null | tail -n1 | tr -d ' ')
            printf '%s\t%s\t%d\t%d\t%s\t%s\t%s\n' "${TASK}" "${ROI}" "${NSUB}" "${NTR}" \
                "${MODE:-NA}" "${HDI%%,*}" "${HDI##*,}" >> "${ROI_SUMMARY}"
            echo "     TRR mode = ${MODE:-NA}  95% HDI [${HDI:-NA}]"
        done
    done
    echo "  summary: ${ROI_SUMMARY}  (ROI peak coordinates: rois_task-<T>_clusters.txt)"
fi

# ------------------------------------------------------------------------------
# STEP 3. dice -- within- vs between-participant Dice.
#
# Runs dice_within_vs_between.R (from SCRIPTS) on reproducibility_afni.sh's
# dice_within_between.tsv. The R script prints its results to this job's log.
# ------------------------------------------------------------------------------
if want dice; then
    echo "=============================================================="
    echo "STEP 3: within- vs between-participant Dice"
    echo "=============================================================="
    WB=${REPRO_OUT}/dice_within_between.tsv
    if [ ! -f "${WB}" ]; then
        echo "skip: no ${WB} (run reproducibility_afni.sh first)"
    elif [ ! -f "${SCRIPTS}/dice_within_vs_between.R" ]; then
        echo "skip: no ${SCRIPTS}/dice_within_vs_between.R" >&2
    else
        ${AFNI_R} "Rscript ${SCRIPTS}/dice_within_vs_between.R ${WB}" || echo "  WARN: dice_within_vs_between.R failed" >&2
    fi
fi

echo "=============================================================="
echo "reliability models done: ${STEPS[*]}"
echo "=============================================================="
