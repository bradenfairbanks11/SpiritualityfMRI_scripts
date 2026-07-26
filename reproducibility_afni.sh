#!/bin/bash
#SBATCH --job-name=spirituality_repro
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com

# =============================================================================
# Spirituality fMRI — TEST-RETEST REPRODUCIBILITY (ses-1 vs ses-2), AFNI
#
# For each subject/task, compares the FIRST session's first-level result to the
# SECOND session's, to ask: how reproducible is the activation across sessions?
#
# COMPARISON-ONLY: this does NOT re-run any GLM. It consumes the _stats_REML
# buckets already produced by first_level_afni.sh (so the GLM settings are
# identical by construction) and computes three metrics per task:
#
#   1. Spatial correlation (3ddot -docor) of ses-1 vs ses-2 maps within a brain
#      mask -> one Pearson r per map. Computed for BOTH the beta (coef) and the
#      t-stat map of each contrast.
#   2. Dice overlap (3ddot -dodice) of the two sessions' t-maps thresholded at
#      p < PTHR (two-sided) -> do the same regions survive in both sessions?
#   3. Scatter data (3dmaskdump) -> paired (ses-1, ses-2) voxel values dumped to
#      a table for a scatter/regression plot outside AFNI.
#
# Contrasts compared (sub-brick labels in the REML bucket):
#   mean_response  -> mean_response#0_Coef , mean_response#0_Tstat  (task activation)
#   rating_slope   -> rating_slope#0_Coef  , rating_slope#0_Tstat   (parametric slope)
#
# Both sessions are in the same MNI152NLin2009cAsym grid (same warp template),
# so the comparison is voxelwise with no resampling.
#
# DATA NOTE: only subjects with BOTH ses-1 and ses-2 first-levels are usable.
# Today that is sub-05 only; the script skips any subject/task missing a session.
#
# Run with:  sbatch reproducibility_afni.sh sub-05     (one or several labels)
#            sbatch reproducibility_afni.sh            (all sub-* with first-levels)
# =============================================================================

# NOTE: intentionally NOT using `set -e`. A single AFNI/container failure on one
# subject/task must not abort the whole batch -- each task is guarded explicitly
# and skips forward on failure (see mask guard + threshold guard below). pipefail
# is on so a failing stage in a `cmd | tail | awk` pipeline is detectable.
set -o pipefail

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# --- Paths (match first_level_afni.sh) ---
TEDANA_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/tedana
AFNI_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/afni_firstlvl
AFNI_SIF=/apps/afni/afni_make_build_latest.sif
BIND=/nobackup/archive/usr/bradenf4

# --- New output root for reproducibility results ---
REPRO_OUT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives/afni_reproducibility

# --- Tasks (same labels as first_level_afni.sh) ---
TASKS=(scripture FHS architecture)

# --- Contrasts to compare (REML sub-brick label prefixes) ---
CONTRASTS=(mean_response rating_slope)

# --- Dice threshold: voxelwise p, two-sided (liberal on purpose, see header) ---
PTHR=0.05

# --- The two sessions we compare (test-retest = ses-1 vs ses-2) ---
SES1=ses-1
SES2=ses-2

# ==============================================================================

mkdir -p logs "${REPRO_OUT}"

module load apptainer/1.3.6-qycanb2

# AFNI_DECONFLICT=OVERWRITE so re-runs overwrite intermediate files cleanly.
AFNI="apptainer exec --env AFNI_DECONFLICT=OVERWRITE --bind ${BIND}:${BIND} ${AFNI_SIF} bash -c"

# Global summary table (one row per subject / task / contrast / map / metric).
SUMMARY=${REPRO_OUT}/reproducibility_summary.tsv
if [ ! -f "${SUMMARY}" ]; then
    printf 'participant\ttask\tcontrast\tmap\tmetric\tvalue\n' > "${SUMMARY}"
fi

add_metric() {  # args: participant task contrast map metric value
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" >> "${SUMMARY}"
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

        # First-level REML buckets for the two sessions (AFNI +orig HEAD/BRIK).
        S1=${AFNI_OUT}/${PARTICIPANT_ID}/${SES1}/task-${TASK}/${PARTICIPANT_ID}_${SES1}_task-${TASK}_stats_REML+orig
        S2=${AFNI_OUT}/${PARTICIPANT_ID}/${SES2}/task-${TASK}/${PARTICIPANT_ID}_${SES2}_task-${TASK}_stats_REML+orig

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
        # STEP A: brain mask = intersection of each session's automask.
        # (3dDeconvolve ran without -mask, so the buckets hold values for every
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
        # Per contrast: spatial correlation (beta & t), Dice (t), scatter, diff.
        # ------------------------------------------------------------------
        for CON in "${CONTRASTS[@]}"; do
            COEF="${CON}#0_Coef"
            TSTAT="${CON}#0_Tstat"

            # --- STEP B: spatial correlation (Pearson r) within the mask ---
            R_BETA=$(afni_scalar  "r(beta) ${PARTICIPANT_ID}/${TASK}/${CON}"  "3ddot -docor -mask ${MASK} ${S1}'[${COEF}]'  ${S2}'[${COEF}]'")
            R_TSTAT=$(afni_scalar "r(t) ${PARTICIPANT_ID}/${TASK}/${CON}"     "3ddot -docor -mask ${MASK} ${S1}'[${TSTAT}]' ${S2}'[${TSTAT}]'")
            add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" beta  spatial_r "${R_BETA}"
            add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat spatial_r "${R_TSTAT}"

            # --- STEP C: Dice overlap of thresholded t-maps ---
            # Threshold each session at p<PTHR two-sided, using that run's own DoF
            # (p2dsetstat reads the degrees of freedom from the sub-brick header).
            THR1=$(${AFNI} "p2dsetstat -quiet -inset ${S1}'[${TSTAT}]' -pval ${PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            THR2=$(${AFNI} "p2dsetstat -quiet -inset ${S2}'[${TSTAT}]' -pval ${PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            if [ -z "${THR1}" ] || [ -z "${THR2}" ]; then
                # Empty threshold would make 3dcalc's expr 'step(abs(a)-)' a syntax
                # error. Record dice=NA and move on instead of crashing the task.
                echo "  WARN: ${CON} -- empty p2dsetstat threshold (thr1='${THR1}' thr2='${THR2}'); recording dice=NA" >&2
                DICE=NA
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice NA
            else
                ACT1=${OUT}/act_${SES1}_${CON}.nii.gz
                ACT2=${OUT}/act_${SES2}_${CON}.nii.gz
                ${AFNI} "
                    3dcalc -overwrite -a ${S1}'[${TSTAT}]' -expr 'step(abs(a)-${THR1})' -prefix ${ACT1}
                    3dcalc -overwrite -a ${S2}'[${TSTAT}]' -expr 'step(abs(a)-${THR2})' -prefix ${ACT2}
                "
                DICE=$(afni_scalar "dice ${PARTICIPANT_ID}/${TASK}/${CON}" "3ddot -dodice -mask ${MASK} ${ACT1} ${ACT2}")
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice "${DICE}"
            fi

            # --- STEP D: scatter data (paired beta values within the mask) ---
            # Columns: i j k ses1 ses2  -> ready for a scatter/regression plot.
            SCATTER=${OUT}/scatter_${CON}.txt
            rm -f "${SCATTER}"
            ${AFNI} "3dmaskdump -mask ${MASK} -o ${SCATTER} ${S1}'[${COEF}]' ${S2}'[${COEF}]'"

            # --- difference map (ses2 - ses1) of the beta, for eyeballing ---
            ${AFNI} "3dcalc -overwrite -a ${S2}'[${COEF}]' -b ${S1}'[${COEF}]' -expr 'a-b' -prefix ${OUT}/ses2_minus_ses1_${CON}_Coef.nii.gz"

            echo "  ${CON}: r(beta)=${R_BETA}  r(t)=${R_TSTAT}  dice=${DICE}"
        done
    done
done

if [ "${PROCESSED}" -eq 0 ]; then
    echo "ERROR: no subject/task comparison ran (all skipped or failed). Not a successful run." >&2
    exit 1
fi

echo "Reproducibility comparison complete (${PROCESSED} task-comparisons). Summary: ${SUMMARY}"

# =============================================================================
# ===== GROUP ROLL-UP: uncomment once >=2 two-session subjects are processed ===
# =============================================================================
# Aggregates reproducibility_summary.tsv across ALL processed subjects into a
# per-task group table (mean +/- SD of each metric). Meaningless at n=1 (only
# sub-05 has two sessions today), so it is left commented out until more
# two-session subjects exist. To activate: delete the leading "# " on each line.
#
# GROUP_SUMMARY=${REPRO_OUT}/reproducibility_group_summary.tsv
# awk -F'\t' '
#     NR==1 { next }                                   # skip header
#     {
#         key = $2 SUBSEP $3 SUBSEP $4 SUBSEP $5        # task, contrast, map, metric
#         n[key]++
#         sum[key]   += $6
#         sumsq[key] += $6*$6
#     }
#     END {
#         print "task\tcontrast\tmap\tmetric\tn\tmean\tsd"
#         for (k in n) {
#             split(k, a, SUBSEP)
#             m  = sum[k]/n[k]
#             sd = (n[k]>1) ? sqrt((sumsq[k] - n[k]*m*m)/(n[k]-1)) : 0
#             printf "%s\t%s\t%s\t%s\t%d\t%.4f\t%.4f\n", a[1], a[2], a[3], a[4], n[k], m, sd
#         }
#     }
# ' "${SUMMARY}" | sort > "${GROUP_SUMMARY}"
# echo "Group summary: ${GROUP_SUMMARY}"
