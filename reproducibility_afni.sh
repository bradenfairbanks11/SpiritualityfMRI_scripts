#!/bin/bash
#SBATCH --job-name=spirituality_repro
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com

# =============================================================================
# Spirituality fMRI — TEST-RETEST RELIABILITY (ses-1 vs ses-2), AFNI
#
# Quantifies reproducibility at TWO levels, and does NOT re-run any GLM (it
# consumes the _stats_REML buckets from first_level_afni.sh, so GLM settings are
# identical by construction).
#
# LEVEL 1 — WITHIN-SUBJECT, per subject/task/contrast (ses-1 vs ses-2):
#   * spatial_r  (3ddot -docor)  : Pearson r of the two sessions' maps -> pattern
#       CONSISTENCY. Blind to a systematic session offset/scaling. Computed on
#       both the beta (Coef) and the t-stat map.
#   * spatial_icc (irr::icc in R) : ICC(2,1) absolute agreement of the same
#       paired-voxel betas (voxels = targets, the 2 sessions = measurements).
#       Same granularity as spatial_r but ALSO penalizes session offset/scaling.
#       Computed with the established `irr` package (see spatial_icc.R), not
#       hand-rolled, for provenance.
#   * dice       (3ddot -dodice) : overlap of the two sessions' t-maps thresholded
#       at p<PTHR (two-sided) -> do the same regions survive in both sessions?
#   * scatter    (3dmaskdump)    : paired (ses-1, ses-2) betas dumped for plotting.
#
# LEVEL 2 — ACROSS-SUBJECT (group), per task/contrast:
#   * 3dICC ICC(2,1) two-way random effects -> a voxelwise ICC MAP, estimated by
#       pooling all subjects that have BOTH sessions (subjects = targets, the 2
#       sessions = measurements). See the GROUP section below.
#   NOTE: group ICC measures reliability of BETWEEN-SUBJECT differences. A region
#   where everyone activates strongly AND identically has low between-subject
#   variance -> LOW ICC even though activation is strong. So high ICC != "task
#   reliably activates" (that is the group activation map). At today's tiny n of
#   two-session subjects the group ICC is a wired-up placeholder that only becomes
#   meaningful as subjects accrue -- it is flagged loudly at run time.
#
# MASKING (fixes whole-brain inflation): every metric above is computed within a
# GROUP-LEVEL mean_response ACTIVATION mask (a functional ROI from the across-
# subject group map, independent of the ses1-vs-ses2 comparison -> avoids
# double-dipping; robust contrast only -> not gated on the noisy rating_slope),
# AND within the whole-brain mask as a reference. The summary TSV carries a `mask`
# column (activation | wholebrain) to distinguish them.
#
# Contrasts compared (sub-brick labels in the REML bucket):
#   mean_response  -> mean_response#0_Coef , mean_response#0_Tstat  (task activation)
#   rating_slope   -> rating_slope#0_Coef  , rating_slope#0_Tstat   (parametric slope)
#
# All subjects/sessions share the MNI152NLin2009cAsym grid (same warp template),
# so comparisons are voxelwise with no resampling (mask is NN-resampled to be safe).
#
# Run with:  sbatch reproducibility_afni.sh sub-05     (one or several labels)
#            sbatch reproducibility_afni.sh            (all sub-* with first-levels)
# The GROUP ICC always pools ALL sub-* that have both sessions, regardless of args.
# =============================================================================

# NOTE: intentionally NOT using `set -e`. A single AFNI/container failure on one
# subject/task must not abort the whole batch -- each task is guarded explicitly
# and skips forward on failure (see mask guard + threshold guard below). pipefail
# is on so a failing stage in a `cmd | tail | awk` pipeline is detectable.
set -o pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# --- Paths (match first_level_afni.sh / second_level_afni.sh) ---
PROJECT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project
TEDANA_OUT=${PROJECT}/derivatives/tedana
AFNI_OUT=${PROJECT}/derivatives/afni_firstlvl
AFNI_SIF=/apps/afni/afni_make_build_latest.sif
BIND=/nobackup/archive/usr/bradenf4

# --- New output root for reproducibility results ---
REPRO_OUT=${PROJECT}/derivatives/afni_reproducibility

# --- Group-level mean_response activation map (from second_level_afni.sh STEP 3).
#     3dMEMA one-sample output: [0]=effect, [1]=stat. Used to build the functional
#     ROI mask. If absent, the script falls back to whole-brain only. ---
GROUP_MEMA=${PROJECT}/derivatives/afni_group/mema

# --- Persistent R library shared with second_level_afni.sh (container R lacks the
#     packages; installed once into archive and exposed via R_LIBS_USER). ---
RLIB=${PROJECT}/derivatives/afni_group/Rlib

# --- Tasks (same labels as first_level_afni.sh) ---
TASKS=(scripture FHS architecture)

# --- Contrasts to compare (REML sub-brick label prefixes) ---
CONTRASTS=(mean_response rating_slope)

# --- Dice threshold: voxelwise p, two-sided (liberal on purpose, see header) ---
PTHR=0.05

# --- Activation-mask threshold: voxelwise p, two-sided, on the group mean map ---
ACT_PTHR=0.05

# --- The two sessions we compare (test-retest = ses-1 vs ses-2) ---
SES1=ses-1
SES2=ses-2

# --- Where this script (and spatial_icc.R) live, so we can find the R helper.
#     Under SLURM this is the submit dir; otherwise the script's own dir. ---
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ==============================================================================

mkdir -p logs "${REPRO_OUT}" "${RLIB}"

module load apptainer/1.3.6-qycanb2

# AFNI_DECONFLICT=OVERWRITE so re-runs overwrite intermediate files cleanly.
AFNI="apptainer exec --env AFNI_DECONFLICT=OVERWRITE --bind ${BIND}:${BIND} ${AFNI_SIF} bash -c"
# R-enabled invocation for the ICC steps (irr for spatial ICC, 3dICC for the group
# map). Exposes the persistent R library so those R programs find their packages.
AFNI_R="apptainer exec --env AFNI_DECONFLICT=OVERWRITE --env R_LIBS_USER=${RLIB} --bind ${BIND}:${BIND} ${AFNI_SIF} bash -c"

# --- Bootstrap: install the ICC R packages into RLIB on first run ---
#   spatial ICC (irr) + 3dICC (blme, lme4, metafor, snow).
ICC_PKGS='c("irr","lme4","blme","metafor","snow")'
if ! ${AFNI_R} "Rscript -e 'q(status=length(setdiff(${ICC_PKGS}, rownames(installed.packages()))))'"; then
    echo "Installing ICC R packages into ${RLIB} (irr, lme4, blme, metafor, snow) ..."
    ${AFNI_R} "Rscript -e 'options(repos=c(CRAN=\"https://cloud.r-project.org\")); p<-setdiff(${ICC_PKGS}, rownames(installed.packages())); if(length(p)) install.packages(p, lib=Sys.getenv(\"R_LIBS_USER\"))'"
fi

# --- Stage the spatial-ICC R helper somewhere the container can read (under BIND) ---
ICC_R=${REPRO_OUT}/spatial_icc.R
if [ -f "${SCRIPT_DIR}/spatial_icc.R" ]; then
    cp -f "${SCRIPT_DIR}/spatial_icc.R" "${ICC_R}"
else
    echo "WARNING: spatial_icc.R not found in ${SCRIPT_DIR}; spatial_icc will record NA" >&2
    ICC_R=""
fi

# Global summary table (one row per subject / task / contrast / map / metric / mask).
SUMMARY=${REPRO_OUT}/reproducibility_summary.tsv
if [ ! -f "${SUMMARY}" ]; then
    printf 'participant\ttask\tcontrast\tmap\tmetric\tmask\tvalue\n' > "${SUMMARY}"
fi

# Group ICC summary (one row per task / contrast).
ICC_SUMMARY=${REPRO_OUT}/reproducibility_icc_summary.tsv

add_metric() {  # args: participant task contrast map metric mask value
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$7" >> "${SUMMARY}"
}

# Run an AFNI command that is expected to print a single scalar on its last line
# and echo that value. If nothing comes back (AFNI/container error, unknown flag,
# missing sub-brick), warn loudly and echo "NA" so the gap is visible in the TSV
# instead of a silent blank cell.
afni_scalar() {  # args: description  afni_command_string
    local desc="$1" cmd="$2" val
    val=$(${AFNI} "${cmd}" 2>/dev/null | tail -n1 | awk '{print $1}')
    if [ -z "${val}" ]; then
        echo "  WARN: ${desc} -- AFNI produced no value; recording NA" >&2
        val=NA
    fi
    printf '%s' "${val}"
}

# Spatial ICC(2,1) of a masked pair of beta maps, via the irr package. Echoes the
# ICC value, or NA on any failure (missing helper, degenerate mask, R error).
spatial_icc() {  # args: description  mask  s1_coef_selector  s2_coef_selector
    local desc="$1" mask="$2" a="$3" b="$4" pairs val
    if [ -z "${ICC_R}" ]; then printf 'NA'; return; fi
    # Serial script -> a single reusable temp file is fine and portable (avoids
    # GNU mktemp's dislike of a suffix after the X's).
    pairs="${OUT}/.icc_pairs.txt"
    rm -f "${pairs}"
    ${AFNI} "3dmaskdump -noijk -mask ${mask} -o ${pairs} ${a} ${b}" >/dev/null 2>&1
    if [ ! -s "${pairs}" ]; then
        echo "  WARN: spatial_icc ${desc} -- empty maskdump; recording NA" >&2
        rm -f "${pairs}"; printf 'NA'; return
    fi
    val=$(${AFNI_R} "Rscript ${ICC_R} ${pairs}" 2>/dev/null | tail -n1 | awk '{print $1}')
    rm -f "${pairs}"
    [ -z "${val}" ] && val=NA
    printf '%s' "${val}"
}

# ==============================================================================
# LOOP: participants -> tasks   (ses-1 vs ses-2 compared inside each task)
# ==============================================================================

# Participants to process. Pass labels as args (e.g. "sub-05" or "05"); with NO
# args, defaults to ALL sub-* that have first-level output.
if [ "$#" -gt 0 ]; then
    PARTICIPANT_DIRS=()
    for arg in "$@"; do
        PARTICIPANT_DIRS+=("${AFNI_OUT}/sub-${arg#sub-}")
    done
else
    echo "WARNING: no participant labels given -> defaulting to ALL sub-* in ${AFNI_OUT}"
    PARTICIPANT_DIRS=("${AFNI_OUT}"/sub-*)
fi

# Count task-comparisons that actually ran. If NOTHING runs (e.g. the whole node
# has a broken container), we exit non-zero at the end so the FAIL email fires
# rather than reporting a silently "successful" empty run.
PROCESSED=0

for PARTICIPANT_DIR in "${PARTICIPANT_DIRS[@]}"; do
    if [ ! -d "${PARTICIPANT_DIR}" ]; then
        echo "ERROR: ${PARTICIPANT_DIR} not found in first-level output -- skipping."
        continue
    fi
    PARTICIPANT_ID=$(basename "${PARTICIPANT_DIR}")

    # Idempotent re-runs: drop any prior rows for this participant from the
    # summary so a repeated run refreshes rather than duplicates them.
    if [ -f "${SUMMARY}" ]; then
        awk -F'\t' -v s="${PARTICIPANT_ID}" 'NR==1 || $1!=s' "${SUMMARY}" > "${SUMMARY}.tmp" \
            && mv "${SUMMARY}.tmp" "${SUMMARY}"
    fi

    for TASK in "${TASKS[@]}"; do

        # First-level REML buckets for the two sessions (AFNI +tlrc HEAD/BRIK).
        S1=${AFNI_OUT}/${PARTICIPANT_ID}/${SES1}/task-${TASK}/${PARTICIPANT_ID}_${SES1}_task-${TASK}_stats_REML+tlrc
        S2=${AFNI_OUT}/${PARTICIPANT_ID}/${SES2}/task-${TASK}/${PARTICIPANT_ID}_${SES2}_task-${TASK}_stats_REML+tlrc

        if [ ! -f "${S1}.HEAD" ] || [ ! -f "${S2}.HEAD" ]; then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- missing a session's REML bucket (need both ${SES1} and ${SES2})"
            continue
        fi

        # tedana bolds (only needed to build the brain mask)
        BOLD1=${TEDANA_OUT}/${PARTICIPANT_ID}/${SES1}/func/${PARTICIPANT_ID}_${SES1}_task-${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
        BOLD2=${TEDANA_OUT}/${PARTICIPANT_ID}/${SES2}/func/${PARTICIPANT_ID}_${SES2}_task-${TASK}_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz
        if [ ! -f "${BOLD1}" ] || [ ! -f "${BOLD2}" ]; then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- missing a tedana bold for the mask"
            continue
        fi

        OUT=${REPRO_OUT}/${PARTICIPANT_ID}/task-${TASK}
        mkdir -p "${OUT}"
        echo "=== ${PARTICIPANT_ID} task-${TASK}: comparing ${SES1} vs ${SES2} ==="

        # ------------------------------------------------------------------
        # STEP A: whole-brain mask = intersection of each session's automask.
        # (3dDeconvolve now runs with -mask, so buckets are already zero outside
        # the brain; this automask intersection stays as the reference mask for
        #  voxel; restricting to brain keeps the metrics meaningful.)
        # ------------------------------------------------------------------
        MASK=${OUT}/mask_both.nii.gz
        ${AFNI} "
            3dAutomask -overwrite -prefix ${OUT}/mask_${SES1}.nii.gz ${BOLD1}
            3dAutomask -overwrite -prefix ${OUT}/mask_${SES2}.nii.gz ${BOLD2}
            3dmask_tool -overwrite -input ${OUT}/mask_${SES1}.nii.gz ${OUT}/mask_${SES2}.nii.gz -inter -prefix ${MASK}
        "

        # If the mask never got built, the container/AFNI failed on this node
        # (this is exactly what killed the 2026-07-20 run: the compute node had
        # no /dev/fuse). Skip this task instead of aborting the whole batch.
        if [ ! -f "${MASK}" ]; then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- mask build failed (AFNI/container error on this node?)"
            continue
        fi
        PROCESSED=$((PROCESSED + 1))

        # ------------------------------------------------------------------
        # STEP A2: group mean_response ACTIVATION mask (functional ROI).
        # From the across-subject group map (mema_<task>_mean, [1]=stat), so it is
        # independent of either session's own threshold (no double-dipping), and it
        # uses the ROBUST contrast only. Thresholded, NN-resampled to this grid, and
        # intersected with the brain mask. Falls back to whole-brain only if absent.
        # ------------------------------------------------------------------
        ACTMASK=${OUT}/mask_activation.nii.gz
        rm -f "${ACTMASK}"
        MEMA_MEAN=${GROUP_MEMA}/mema_${TASK}_mean+tlrc
        if [ -f "${MEMA_MEAN}.HEAD" ]; then
            ATHR=$(${AFNI} "p2dsetstat -quiet -inset ${MEMA_MEAN}'[1]' -pval ${ACT_PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            if [ -n "${ATHR}" ]; then
                ${AFNI} "
                    3dcalc -overwrite -a ${MEMA_MEAN}'[1]' -expr 'step(abs(a)-${ATHR})' -prefix ${OUT}/.act_thr.nii.gz
                    3dresample -overwrite -master ${MASK} -rmode NN -prefix ${OUT}/.act_thr_rs.nii.gz -input ${OUT}/.act_thr.nii.gz
                    3dmask_tool -overwrite -input ${MASK} ${OUT}/.act_thr_rs.nii.gz -inter -prefix ${ACTMASK}
                "
            fi
        fi

        # Build the mask list: activation (if a non-empty ROI was built) + wholebrain.
        MASKS=()
        NACT=""
        if [ -f "${ACTMASK}" ]; then
            NACT=$(${AFNI} "3dBrickStat -non-zero -count ${ACTMASK}" 2>/dev/null | tail -n1 | awk '{print $1}')
        fi
        if [ -n "${NACT}" ] && [ "${NACT}" -gt 0 ] 2>/dev/null; then
            MASKS+=("activation:${ACTMASK}")
        else
            echo "  WARN: ${TASK} -- no group activation ROI (mema_${TASK}_mean missing/empty); activation metrics -> whole-brain only" >&2
        fi
        MASKS+=("wholebrain:${MASK}")

        # ------------------------------------------------------------------
        # Per contrast: spatial_r (beta & t), spatial_icc (beta), Dice (t) within
        # each mask; scatter + diff once (whole-brain).
        # ------------------------------------------------------------------
        for CON in "${CONTRASTS[@]}"; do
            COEF="${CON}#0_Coef"
            TSTAT="${CON}#0_Tstat"

            # Session-thresholded binary t-maps for Dice (mask-independent; the mask
            # is applied by 3ddot -dodice). Built once per contrast.
            THR1=$(${AFNI} "p2dsetstat -quiet -inset ${S1}'[${TSTAT}]' -pval ${PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            THR2=$(${AFNI} "p2dsetstat -quiet -inset ${S2}'[${TSTAT}]' -pval ${PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            ACT1=${OUT}/act_${SES1}_${CON}.nii.gz
            ACT2=${OUT}/act_${SES2}_${CON}.nii.gz
            DICE_OK=0
            if [ -n "${THR1}" ] && [ -n "${THR2}" ]; then
                ${AFNI} "
                    3dcalc -overwrite -a ${S1}'[${TSTAT}]' -expr 'step(abs(a)-${THR1})' -prefix ${ACT1}
                    3dcalc -overwrite -a ${S2}'[${TSTAT}]' -expr 'step(abs(a)-${THR2})' -prefix ${ACT2}
                "
                DICE_OK=1
            else
                echo "  WARN: ${CON} -- empty p2dsetstat threshold (thr1='${THR1}' thr2='${THR2}'); dice=NA" >&2
            fi

            for MSPEC in "${MASKS[@]}"; do
                MTAG=${MSPEC%%:*}
                M=${MSPEC#*:}

                # spatial correlation (Pearson r) -- consistency
                R_BETA=$(afni_scalar  "r(beta) ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}"  "3ddot -docor -mask ${M} ${S1}'[${COEF}]'  ${S2}'[${COEF}]'")
                R_TSTAT=$(afni_scalar "r(t) ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}"     "3ddot -docor -mask ${M} ${S1}'[${TSTAT}]' ${S2}'[${TSTAT}]'")
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" beta  spatial_r "${MTAG}" "${R_BETA}"
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat spatial_r "${MTAG}" "${R_TSTAT}"

                # spatial ICC(2,1) absolute agreement (beta) -- via irr
                SPICC=$(spatial_icc "${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "${M}" "${S1}'[${COEF}]'" "${S2}'[${COEF}]'")
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" beta spatial_icc "${MTAG}" "${SPICC}"

                # Dice overlap of thresholded t-maps
                if [ "${DICE_OK}" -eq 1 ]; then
                    DICE=$(afni_scalar "dice ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "3ddot -dodice -mask ${M} ${ACT1} ${ACT2}")
                else
                    DICE=NA
                fi
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice "${MTAG}" "${DICE}"

                echo "  ${CON} [${MTAG}]: r(beta)=${R_BETA} r(t)=${R_TSTAT} sICC=${SPICC} dice=${DICE}"
            done

            # --- scatter data (paired beta values, whole-brain) for plotting ---
            # Columns: i j k ses1 ses2  -> ready for a scatter/regression plot.
            SCATTER=${OUT}/scatter_${CON}.txt
            rm -f "${SCATTER}"
            ${AFNI} "3dmaskdump -mask ${MASK} -o ${SCATTER} ${S1}'[${COEF}]' ${S2}'[${COEF}]'"

            # --- difference map (ses2 - ses1) of the beta, for eyeballing ---
            ${AFNI} "3dcalc -overwrite -a ${S2}'[${COEF}]' -b ${S1}'[${COEF}]' -expr 'a-b' -prefix ${OUT}/ses2_minus_ses1_${CON}_Coef.nii.gz"
        done
    done
done

if [ "${PROCESSED}" -eq 0 ]; then
    echo "ERROR: no subject/task comparison ran (all skipped or failed). Not a successful run." >&2
    exit 1
fi

echo "Within-subject reproducibility complete (${PROCESSED} task-comparisons). Summary: ${SUMMARY}"

# =============================================================================
# GROUP VOXELWISE ICC(2,1)  —  AFNI 3dICC, across-subject reliability map
# =============================================================================
# Pools ALL subjects that have BOTH sessions for a task (subjects = targets, the
# 2 sessions = repeated measurements) and fits, per voxel, the two-way random-
# effects model  '1+(1|session)+(1|Subj)'  = ICC(2,1). Output is an ICC MAP per
# task x contrast, plus a whole-mask median in reproducibility_icc_summary.tsv.
#
# Inputs are the beta (Coef) maps, extracted to MNI-labelled NIfTI (like
# second_level_afni.sh). The mask is the group mean_response activation ROI
# (falls back to input coverage). Runs regardless of the participant args above
# (the group is defined by whoever has both sessions on disk).
#
# REMEMBER: high ICC == subjects' DIFFERENCES replicate across sessions, NOT
# "the task reliably activates" (that is the group activation map). At small n
# this map is very unstable -- see the loud warning below.
# =============================================================================
echo "=============================================================="
echo "GROUP ICC(2,1): voxelwise 3dICC across two-session subjects"
echo "=============================================================="

GROUP_OUT=${REPRO_OUT}/group
mkdir -p "${GROUP_OUT}/inputs"

# Fresh group summary each run (recomputed over everyone with both sessions).
printf 'task\tcontrast\tn\tmetric\tvalue\n' > "${ICC_SUMMARY}"

for TASK in "${TASKS[@]}"; do
    # Discover subjects with BOTH ses-1 and ses-2 REML buckets for this task.
    ICC_SUBS=()
    for d in "${AFNI_OUT}"/sub-*; do
        [ -d "${d}" ] || continue
        sid=$(basename "${d}")
        h1=${AFNI_OUT}/${sid}/${SES1}/task-${TASK}/${sid}_${SES1}_task-${TASK}_stats_REML+tlrc.HEAD
        h2=${AFNI_OUT}/${sid}/${SES2}/task-${TASK}/${sid}_${SES2}_task-${TASK}_stats_REML+tlrc.HEAD
        { [ -f "${h1}" ] && [ -f "${h2}" ]; } && ICC_SUBS+=("${sid}")
    done
    n=${#ICC_SUBS[@]}
    if [ "${n}" -lt 2 ]; then
        echo "GROUP ICC ${TASK}: only ${n} subject(s) with both sessions -> skip (need >=2)."
        continue
    fi
    if [ "${n}" -lt 10 ]; then
        echo "GROUP ICC ${TASK}: *** WARNING n=${n} is very small -- the voxelwise ICC will be" >&2
        echo "                    highly unstable/uninterpretable until more two-session subjects" >&2
        echo "                    accrue. Wired up now; treat results as placeholder. ***" >&2
    fi
    echo "GROUP ICC ${TASK}: subjects with both sessions (n=${n}): ${ICC_SUBS[*]}"

    for CON in "${CONTRASTS[@]}"; do
        COEF="${CON}#0_Coef"

        # Extract each subject/session beta to MNI-labelled NIfTI; build data table.
        DT=""
        for sid in "${ICC_SUBS[@]}"; do
            for SES in "${SES1}" "${SES2}"; do
                SB=${AFNI_OUT}/${sid}/${SES}/task-${TASK}/${sid}_${SES}_task-${TASK}_stats_REML+tlrc
                NI=${GROUP_OUT}/inputs/${sid}_${SES}_task-${TASK}_${CON}_Coef.nii.gz
                rm -f "${NI}"
                ${AFNI} "3dbucket -prefix ${NI} ${SB}'[${COEF}]'; 3drefit -space MNI ${NI}"
                seslab=${SES//-/}          # ses-1 -> ses1 (a valid factor level)
                DT="${DT} ${sid} ${seslab} ${NI}"
            done
        done

        # Group mask: full coverage of the inputs, intersected with the group
        # mean_response activation ROI when available.
        GMASK=${GROUP_OUT}/groupmask_task-${TASK}.nii.gz
        rm -f "${GMASK}"
        ${AFNI} "3dmask_tool -overwrite -input ${GROUP_OUT}/inputs/*_task-${TASK}_${CON}_Coef.nii.gz -frac 1.0 -prefix ${GMASK}"
        MEMA_MEAN=${GROUP_MEMA}/mema_${TASK}_mean+tlrc
        if [ -f "${MEMA_MEAN}.HEAD" ] && [ -f "${GMASK}" ]; then
            GATHR=$(${AFNI} "p2dsetstat -quiet -inset ${MEMA_MEAN}'[1]' -pval ${ACT_PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            if [ -n "${GATHR}" ]; then
                ${AFNI} "
                    3dcalc -overwrite -a ${MEMA_MEAN}'[1]' -expr 'step(abs(a)-${GATHR})' -prefix ${GROUP_OUT}/.gact.nii.gz
                    3dresample -overwrite -master ${GMASK} -rmode NN -prefix ${GROUP_OUT}/.gact_rs.nii.gz -input ${GROUP_OUT}/.gact.nii.gz
                    3dmask_tool -overwrite -input ${GMASK} ${GROUP_OUT}/.gact_rs.nii.gz -inter -prefix ${GROUP_OUT}/groupmask_act_task-${TASK}.nii.gz
                "
                [ -f "${GROUP_OUT}/groupmask_act_task-${TASK}.nii.gz" ] && GMASK=${GROUP_OUT}/groupmask_act_task-${TASK}.nii.gz
            fi
        fi

        ICC_OUT=${GROUP_OUT}/ICC_task-${TASK}_${CON}
        rm -f "${ICC_OUT}"+*.HEAD "${ICC_OUT}"+*.BRIK*
        echo "  3dICC ${TASK} ${CON}: n=${n} subjects x 2 sessions, model 1+(1|session)+(1|Subj)"
        if ! ${AFNI_R} "3dICC -prefix ${ICC_OUT} -jobs ${SLURM_CPUS_PER_TASK:-4} -mask ${GMASK} \
                -model '1+(1|session)+(1|Subj)' \
                -dataTable Subj session InputFile ${DT}"; then
            echo "  WARN: 3dICC failed for ${TASK}/${CON} (R package / container / grid?); recording NA." >&2
            printf '%s\t%s\t%d\t%s\t%s\n' "${TASK}" "${CON}" "${n}" median_icc NA >> "${ICC_SUMMARY}"
            continue
        fi

        # Median ICC within the mask. 3dICC output sub-brick [0] is the ICC value
        # (sub-brick [1] is its F-stat). Verify labels once with: 3dinfo ICC_*+tlrc
        MED=$(${AFNI} "3dBrickStat -median -mask ${GMASK} ${ICC_OUT}+tlrc'[0]'" 2>/dev/null | tail -n1 | awk '{print $NF}')
        [ -z "${MED}" ] && MED=NA
        printf '%s\t%s\t%d\t%s\t%s\n' "${TASK}" "${CON}" "${n}" median_icc "${MED}" >> "${ICC_SUMMARY}"
        echo "  ICC map -> ${ICC_OUT}+tlrc   median ICC (in mask) = ${MED}"
    done
done

echo "=============================================================="
echo "Reliability run complete."
echo "  within-subject metrics : ${SUMMARY}   (mask = activation | wholebrain)"
echo "  group ICC maps         : ${GROUP_OUT}/ICC_task-*_*+tlrc"
echo "  group ICC medians      : ${ICC_SUMMARY}"
echo "  NOTE: group ICC needs many subjects to be meaningful; small-n results"
echo "        above are placeholders and will stabilize as the sample grows."
echo "=============================================================="

# =============================================================================
# ===== OPTIONAL WITHIN-SUBJECT ROLL-UP: mean +/- SD across subjects ==========
# =============================================================================
# Aggregates reproducibility_summary.tsv across processed subjects into a per-
# (task, contrast, map, metric, mask) table. Uncomment once several subjects
# exist. Value is column 7; keys are columns 2..6 (mask included).
#
# GROUP_SUMMARY=${REPRO_OUT}/reproducibility_group_summary.tsv
# awk -F'\t' '
#     NR==1 { next }                                        # skip header
#     $7=="NA" { next }                                     # ignore NA cells
#     {
#         key = $2 SUBSEP $3 SUBSEP $4 SUBSEP $5 SUBSEP $6  # task,contrast,map,metric,mask
#         n[key]++; sum[key]+=$7; sumsq[key]+=$7*$7
#     }
#     END {
#         print "task\tcontrast\tmap\tmetric\tmask\tn\tmean\tsd"
#         for (k in n) {
#             split(k, a, SUBSEP)
#             m  = sum[k]/n[k]
#             sd = (n[k]>1) ? sqrt((sumsq[k] - n[k]*m*m)/(n[k]-1)) : 0
#             printf "%s\t%s\t%s\t%s\t%s\t%d\t%.4f\t%.4f\n", a[1],a[2],a[3],a[4],a[5], n[k], m, sd
#         }
#     }
# ' "${SUMMARY}" | sort > "${GROUP_SUMMARY}"
# echo "Within-subject roll-up: ${GROUP_SUMMARY}"
