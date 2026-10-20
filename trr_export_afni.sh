#!/bin/bash
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

# =============================================================================
# Spirituality fMRI — export TRIAL-LEVEL data for test-retest reliability (TRR)
#
# Runs after first_level_IM_afni.sh and second_level_afni.sh. Pure AFNI C
# programs (no R), so it runs in the cluster container even though the R-based
# AFNI programs do not. The R side (3dLMEr -TRR, TRR) runs in local_r_steps.sh.
#
# For every task, pools the subjects with single-trial betas for BOTH sessions
# and writes a self-contained folder ${TRR_OUT} (copy it whole to the local
# machine; every path inside is relative):
#
#   inputs/<sub>_<ses>_task-<T>_trial-NNN.nii.gz   one 3D beta map per trial
#   trials_task-<T>.tsv     Subj sess trial rating cond InputFile
#                           -> becomes the 3dLMEr -TRR data table (whole brain)
#   mask_task-<T>.nii.gz    intersection of every included run's fMRIPrep mask
#   rois_task-<T>.nii.gz    ROIs = clusters of the group mean_response MEMA map
#                           (voxel p < CLUST_PTHR, cluster size from the
#                           second-level 3dClustSim table, NN1, bisided)
#   rois_task-<T>_clusters.txt   3dClusterize report (size, peak xyz per ROI)
#   roi_task-<T>.tsv        Subj sess trial rating cond ROI Y
#                           -> one TRR (region-level) data table per ROI
#
# ROI selection uses the group AVERAGE activation (sessions averaged in
# second_level_afni.sh STEP 1), never the session-1-vs-2 agreement, so ROIs are
# not selected on the reliability being measured.
#
# Run with:  sbatch trr_export_afni.sh
# =============================================================================

set -o pipefail

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
# Provides: FMRIPREP_OUT AFNI_OUT GROUP_OUT TRR_OUT AFNI_SIF LOGS BIND_ARGS

TASKS=(scripture FHS architecture)
SES1=ses-1
SES2=ses-2

# Must match second_level_afni.sh (CLUST_PTHR there builds the cluster table).
CLUST_PTHR=0.001
MEMADIR=${GROUP_OUT}/mema
CLUST_TABLE=${GROUP_OUT}/clustsim/mema_clustsim.NN1_bisided.1D

mkdir -p "${LOGS}" "${TRR_OUT}/inputs"

module load apptainer/1.3.6-qycanb2
AFNI="apptainer exec --env AFNI_DECONFLICT=OVERWRITE ${BIND_ARGS} ${AFNI_SIF} bash -c"

# Cluster-size threshold (voxels) for alpha = 0.05 at voxel p = CLUST_PTHR.
# The table has one row per -pthr; column 2 is the -athr 0.05 column. Rounded
# UP, since a cluster must be at least this big.
CLUST_K=""
if [ -f "${CLUST_TABLE}" ]; then
    CLUST_K=$(awk -v p="${CLUST_PTHR}" '!/^#/ && NF>=2 && ($1+0)==(p+0) { k=$2; print (k==int(k)) ? k : int(k)+1; exit }' "${CLUST_TABLE}")
fi
if [ -z "${CLUST_K}" ]; then
    echo "WARNING: no cluster-size threshold for p=${CLUST_PTHR} in ${CLUST_TABLE}" >&2
    echo "         (run second_level_afni.sh first). ROI tables will be skipped." >&2
fi

for TASK in "${TASKS[@]}"; do
    echo "=============================================================="
    echo "task-${TASK}"
    echo "=============================================================="

    # --- subjects with single-trial betas for BOTH sessions ---
    SUBS=()
    for d in "${AFNI_OUT}"/sub-*; do
        [ -d "${d}" ] || continue
        sid=$(basename "${d}")
        ok=1
        for SES in "${SES1}" "${SES2}"; do
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

    # --- split each run's IM bucket into one file per trial; build the table ---
    TABLE=${TRR_OUT}/trials_task-${TASK}.tsv
    printf 'Subj\tsess\ttrial\trating\tcond\tInputFile\n' > "${TABLE}"
    MASKS=""
    for sid in "${SUBS[@]}"; do
        for SES in "${SES1}" "${SES2}"; do
            ST=${AFNI_OUT}/${sid}/${SES}/task-${TASK}/single_trial
            PREFIX=${sid}_${SES}_task-${TASK}
            TRIALS=${ST}/${PREFIX}_IM_trials.tsv
            NTRIAL=$(($(wc -l < "${TRIALS}") - 1))
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

    # every exported file must exist, or the downstream tables silently misalign
    nmiss=$(awk -F'\t' -v root="${TRR_OUT}" 'NR>1 { f=root "/" $6; if ((getline line < f) < 0) m++; close(f) } END { print m+0 }' "${TABLE}")
    if [ "${nmiss}" -gt 0 ]; then
        echo "ERROR task-${TASK}: ${nmiss} trial file(s) listed in ${TABLE} were not written." >&2
        continue
    fi
    echo "  trial table: ${TABLE} ($(($(wc -l < "${TABLE}") - 1)) trials)"

    # --- common analysis mask ---
    TMASK=${TRR_OUT}/mask_task-${TASK}.nii.gz
    ${AFNI} "3dmask_tool -overwrite -input ${MASKS} -frac 1.0 -prefix ${TMASK}"

    # --- ROIs from the group mean_response map, then per-trial ROI means ---
    MEMA_MEAN=${MEMADIR}/mema_${TASK}_mean+tlrc
    ROIS=${TRR_OUT}/rois_task-${TASK}.nii.gz
    ROITAB=${TRR_OUT}/roi_task-${TASK}.tsv
    rm -f "${ROIS}" "${ROITAB}"
    if [ -z "${CLUST_K}" ] || [ ! -f "${MEMA_MEAN}.HEAD" ]; then
        echo "  WARN: no MEMA mean map or cluster threshold -> no ROI table for task-${TASK}" >&2
        continue
    fi
    # MEMA output: [0] = effect, [1] = Z/t stat.
    ${AFNI} "3dClusterize -inset ${MEMA_MEAN} -ithr 1 -idat 0 -mask ${TMASK} \
        -NN 1 -bisided p=${CLUST_PTHR} -clust_nvox ${CLUST_K} \
        -pref_map ${ROIS} -outvol_if_no_clust" > "${TRR_OUT}/rois_task-${TASK}_clusters.txt" 2>&1
    NROI=$(${AFNI} "3dBrickStat -max ${ROIS}" 2>/dev/null | tail -n1 | awk '{print int($1)}')
    if [ -z "${NROI}" ] || [ "${NROI}" -lt 1 ]; then
        echo "  WARN task-${TASK}: no group cluster survived (p<${CLUST_PTHR}, k>=${CLUST_K}) -> no ROI table" >&2
        continue
    fi
    echo "  ROIs: ${NROI} cluster(s) -> ${ROIS}"

    # 3dROIstats -quiet prints one line per input file, in input order, with one
    # mean per ROI (ROI values 1..NROI), tab-indented. Pair line i with table
    # row i. Each line is split explicitly with split(..., " ") (ignores leading
    # blanks): changing FS between files is unreliable on BSD/macOS awk, where
    # the first record of the second file is still split with the old FS.
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
