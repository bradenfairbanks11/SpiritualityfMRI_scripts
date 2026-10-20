#!/bin/bash
#SBATCH --job-name=spirituality_repro
#SBATCH --output=/home/bradenf4/spirituality_fmri/logs/%x_%j.out
#SBATCH --error=/home/bradenf4/spirituality_fmri/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com
# Run the job in scratch, not the submit dir. #SBATCH is parsed by Slurm before
# the shell runs, so this must be a literal absolute path -- ${WORK_ROOT} would
# not expand. Keeps stray relative-path output off /home.
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

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
#   * dice_clust (3ddot -dodice) : same overlap, but each session's t-map is
#       cluster-corrected first (voxel p<.001, per-session 3dClustSim, alpha .05).
#       Both-empty pairs score 0 for every Dice (Gorgolewski et al. 2013).
#   * scatter    (3dmaskdump)    : paired (ses-1, ses-2) betas dumped for plotting.
#
# GROUP DICE BASELINE: within- vs between-subject Dice for every task/contrast
# (dice_within_between.tsv), tested by dice_within_vs_between.R.
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
#   A third, DISABLED option -- a leave-one-subject-out ROI (mask=activation_loo)
#   -- lives in STEP A3 below. It removes the residual circularity of building the
#   ROI from a group map that includes the subject being measured. Only worth
#   enabling once the sample is large; see the notes on that block.
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

# --- All paths come from config.sh; edit there, not here ---
# Absolute form on purpose: under sbatch, $0 points at a spool copy, not this repo.
# --- Locate config.sh ---------------------------------------------------------
# Under sbatch, BOTH $0 and ${BASH_SOURCE[0]} point at Slurm's spool copy
# (/var/spool/slurmd/job<N>/slurm_script), NOT at this repo. Verified 2026-08-21
# with probe job 13292417. Deriving the path from BASH_SOURCE alone therefore
# fails inside every batch job — the script dies before doing any work.
#
# Resolved in order of decreasing reliability:
#   $PIPELINE_CONFIG   explicit override (used by audit_paths.sh fixture tests)
#   $SLURM_SUBMIT_DIR  where sbatch was invoked — correct for the normal workflow
#   dirname BASH_SOURCE  correct when run directly with bash, wrong under sbatch
#   the install path   last-resort absolute
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
# Provides: TEDANA_OUT AFNI_OUT REPRO_OUT GROUP_OUT RLIB AFNI_SIF LOGS BIND_ARGS

# --- Group-level mean_response activation map (from second_level_afni.sh STEP 3).
#     3dMEMA one-sample output: [0]=effect, [1]=stat. Used to build the functional
#     ROI mask. If absent, the script falls back to whole-brain only. ---
GROUP_MEMA=${GROUP_OUT}/mema

# NOTE: RLIB now lives on /home (see config.sh), not in the group output tree.
# Compute nodes have NO internet — that is why the original bootstrap install of
# `irr` failed from inside a job and had to be redone from a login node. Keeping the
# library on a backed-up, never-purged filesystem means that cannot recur.

# --- Tasks (same labels as first_level_afni.sh) ---
TASKS=(scripture FHS architecture)

# --- Contrasts to compare (REML sub-brick label prefixes) ---
CONTRASTS=(mean_response rating_slope)

# --- Dice threshold: voxelwise p, two-sided (liberal on purpose, see header) ---
PTHR=0.05

# --- Cluster-corrected Dice (metric dice_clust) ---
# Each session's t-map is thresholded at voxel p < CLUST_PTHR (bisided) and only
# clusters >= k voxels are kept, where k comes from 3dClustSim run on THAT
# session's own residual smoothness (3dFWHMx -acf on its errts_REML) for
# cluster-level alpha CLUST_ATHR. p < .001 is the recommended lower limit for the
# cluster-forming threshold (Woo, Krishnan & Wager 2014, NeuroImage 91:412-419;
# Eklund et al. 2016 PNAS showed p = .01 inflates false positives badly). The
# ACF model is AFNI's post-Eklund fix (Cox et al. 2017, Brain Connectivity 7:152).
CLUST_PTHR=0.001
CLUST_ATHR=0.05

# --- Activation-mask threshold: voxelwise p, two-sided, on the group mean map ---
ACT_PTHR=0.05

# --- The two sessions we compare (test-retest = ses-1 vs ses-2) ---
SES1=ses-1
SES2=ses-2

# --- Where this script (and spatial_icc.R) live, so we can find the R helper.
#     Under SLURM this is the submit dir; otherwise the script's own dir. ---
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ==============================================================================

mkdir -p "${LOGS}" "${REPRO_OUT}"

module load apptainer/1.3.6-qycanb2

# AFNI_DECONFLICT=OVERWRITE so re-runs overwrite intermediate files cleanly.
AFNI="apptainer exec --env AFNI_DECONFLICT=OVERWRITE ${BIND_ARGS} ${AFNI_SIF} bash -c"
# R-enabled invocation for the ICC steps (irr for spatial ICC, 3dICC for the group
# map). Exposes the persistent R library so those R programs find their packages.
AFNI_R="apptainer exec --env AFNI_DECONFLICT=OVERWRITE --env R_LIBS_USER=${RLIB} ${BIND_ARGS} ${AFNI_SIF} bash -c"

# --- Precondition: the ICC R packages must already be in RLIB ---
#   spatial ICC (irr) + 3dICC (blme, lme4, metafor, snow).
# This used to attempt install.packages() inline. Compute nodes have NO internet,
# so that could never reach CRAN from inside a job. Check and fail fast instead.
ICC_PKGS='c("irr","lme4","blme","metafor","snow")'
if ! ${AFNI_R} "Rscript -e 'q(status=length(setdiff(${ICC_PKGS}, rownames(installed.packages()))))'"; then
    echo "ERROR: the ICC R packages are missing from ${RLIB}."
    echo "       (need: irr, lme4, blme, metafor, snow)"
    echo
    echo "They cannot be installed from here — compute nodes have no internet."
    echo "Run this once on a LOGIN NODE, then resubmit:"
    echo "    bash ${PROJ_HOME}/scripts/setup_rlib.sh"
    exit 1
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

# Dice overlap of two binary maps within a mask, echoed as "dice n_a n_b"
# (n = suprathreshold voxels inside the mask). If BOTH maps are empty the Dice is
# 0, not undefined -- the convention of Gorgolewski et al. (2013, NeuroImage
# 69:231-243), "to penalise for lack of signal". Echoes "NA NA NA" on AFNI failure.
dice_pair() {  # args: description  mask  map_a  map_b
    local desc="$1" mask="$2" a="$3" b="$4" na nb d
    read -r na nb d <<< "$(${AFNI} "3dBrickStat -mask ${mask} -non-zero -count ${a}; \
        3dBrickStat -mask ${mask} -non-zero -count ${b}; \
        3ddot -dodice -mask ${mask} ${a} ${b}" 2>/dev/null | awk '{print $1}' | tr '\n' ' ')"
    if [ -z "${na}" ] || [ -z "${nb}" ]; then
        echo "  WARN: dice ${desc} -- AFNI produced no value; recording NA" >&2
        printf 'NA NA NA'; return
    fi
    if [ "${na}" -eq 0 ] && [ "${nb}" -eq 0 ]; then
        d=0
    fi
    [ -z "${d}" ] && d=NA
    printf '%s %s %s' "${d}" "${na}" "${nb}"
}

# Minimum cluster size (voxels) for one session: 3dFWHMx -acf on that session's
# REML residuals -> 3dClustSim. Cached next to `prefix` and recomputed only when
# the residuals are newer than the cached table. Echoes k (rounded up) or "".
session_clust_k() {  # args: errts_dataset(no .HEAD)  mask  prefix
    local errts="$1" mask="$2" pre="$3" acf a b c k
    if [ ! -f "${pre}.NN1_bisided.1D" ] || [ "${errts}.HEAD" -nt "${pre}.NN1_bisided.1D" ]; then
        acf=$(${AFNI} "3dFWHMx -acf ${pre}_acf.1D -mask ${mask} -input ${errts}" 2>/dev/null | tail -n1)
        read -r a b c _ <<< "${acf}"
        if [ -z "${c}" ]; then
            echo "  WARN: 3dFWHMx -acf failed on ${errts}; no cluster threshold" >&2
            return
        fi
        ${AFNI} "OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} 3dClustSim -mask ${mask} -acf ${a} ${b} ${c} \
            -pthr ${CLUST_PTHR} -athr ${CLUST_ATHR} -quiet -prefix ${pre}" >/dev/null 2>&1
    fi
    k=$(awk '!/^#/ && NF>=2 { k=$2; print (k==int(k)) ? k : int(k)+1; exit }' "${pre}.NN1_bisided.1D" 2>/dev/null)
    printf '%s' "${k}"
}

# Binary cluster-corrected map: voxel p < CLUST_PTHR (bisided, from the bucket's
# own t dof), clusters >= k voxels (NN1). An empty map is still written
# (-outvol_if_no_clust) so "no surviving clusters" is a value, not a missing file.
clust_map() {  # args: stats_bucket  tstat_label  mask  k  out.nii.gz
    local idx="${5%.nii.gz}_idx.nii.gz"
    rm -f "${idx}" "${5}"
    ${AFNI} "3dClusterize -inset ${1} -ithr '${2}' -mask ${3} -NN 1 -bisided p=${CLUST_PTHR} \
            -clust_nvox ${4} -pref_map ${idx} -outvol_if_no_clust \
        && 3dcalc -overwrite -a ${idx} -expr 'step(a)' -prefix ${5}" >/dev/null 2>&1
    rm -f "${idx}"
    [ -f "${5}" ]
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
        if why=$(is_excluded "${PARTICIPANT_ID}" "${SES1}" "${TASK}") \
           || why=$(is_excluded "${PARTICIPANT_ID}" "${SES2}" "${TASK}"); then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- EXCLUDED (exclusions.tsv): ${why}"
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

        # ------------------------------------------------------------------
        # STEP A3 (DISABLED) — LEAVE-ONE-SUBJECT-OUT activation ROI
        #
        # WHY: the STEP A2 ROI comes from a group MEMA that INCLUDES this subject,
        # so ~1/n of the evidence selecting the voxels IS the subject being
        # measured. That is a weak form of double dipping. STEP 1 of
        # second_level_afni.sh averages the two sessions before MEMA, so the ROI is
        # blind to the ses1-vs-ses2 DIFFERENCE -- but it still sees their SUM, and
        # selecting on a large |b1+b2| preferentially keeps voxels where the two
        # sessions happened to be large together, nudging spatial_r / spatial_icc
        # upward. LOO removes the dependence: sub-X's ROI is built from everyone
        # EXCEPT sub-X, so selection noise and measurement noise are independent.
        #
        # COST: the LOO MEMA has n-1 subjects, so dof drops and the threshold rises
        # (at n=4: t(2)=4.303 vs t(3)=3.182). Each subject also gets a DIFFERENT
        # ROI, so activation_loo values are not strictly comparable across subjects
        # -- keep the shared-ROI `activation` rows as the comparable set and read
        # the LOO rows as the bias-free companion.
        #
        # DO NOT ENABLE BELOW ~8-10 SUBJECTS. At small n the LOO ROI is so unstable
        # that you trade a known bias for unknown noise. LOO_MIN_N enforces this.
        #
        # PREREQUISITES before uncommenting:
        #   1. second_level_afni.sh must persist its per-subject session-averaged
        #      inputs (they currently exist only in autodelete scratch, which is
        #      purged). In its STEP 7 change
        #          for d in mask mema ttest clustsim; do
        #      to
        #          for d in mask mema ttest clustsim inputs; do
        #      and re-run it.
        #   2. 3dMEMA needs R packages data.table + snow. The bootstrap at the top
        #      of THIS script installs snow but not data.table -- either add
        #      "data.table" to ICC_PKGS, or rely on second_level_afni.sh having
        #      already installed it into the shared RLIB.
        #   3. rm -rf ${REPRO_OUT}/loo_mema whenever the subject set changes: the
        #      LOO maps are cached per (task, excluded subject) and go stale.
        # ------------------------------------------------------------------
        # LOO_MIN_N=8                                   # min OTHER subjects required
        # LOO_INPUTS=${PROJECT}/derivatives/afni_group/inputs
        # LOO_GROUPMASK=${PROJECT}/derivatives/afni_group/mask/group_mask.nii.gz
        # LOO_MEMA_DIR=${REPRO_OUT}/loo_mema
        # ACTMASK_LOO=${OUT}/mask_activation_loo.nii.gz
        # rm -f "${ACTMASK_LOO}"
        # mkdir -p "${LOO_MEMA_DIR}"
        #
        # if [ ! -d "${LOO_INPUTS}" ]; then
        #     echo "  WARN: LOO ROI -- ${LOO_INPUTS} not found (see prerequisite 1); skipping" >&2
        # else
        #     # Build the -set list from every OTHER subject with a session-averaged mean map.
        #     LOO_SET=""; NLOO=0
        #     for c in "${LOO_INPUTS}"/sub-*_task-${TASK}_mean_Coef.nii.gz; do
        #         [ -e "${c}" ] || continue
        #         s=$(basename "${c}"); s=${s%%_task-*}
        #         [ "${s}" = "${PARTICIPANT_ID}" ] && continue        # <-- leave this subject out
        #         t=${LOO_INPUTS}/${s}_task-${TASK}_mean_T.nii.gz
        #         [ -f "${t}" ] || continue
        #         LOO_SET="${LOO_SET} ${s} ${c} ${t}"
        #         NLOO=$((NLOO + 1))
        #     done
        #
        #     if [ "${NLOO}" -lt "${LOO_MIN_N}" ]; then
        #         echo "  LOO ROI ${PARTICIPANT_ID}/${TASK}: only ${NLOO} other subject(s), need >=${LOO_MIN_N} -- skipping" >&2
        #     else
        #         LOO_MEMA=${LOO_MEMA_DIR}/mema_${TASK}_mean_without-${PARTICIPANT_ID}
        #         if [ ! -f "${LOO_MEMA}+tlrc.HEAD" ]; then
        #             echo "  LOO MEMA ${TASK} excluding ${PARTICIPANT_ID} (n=${NLOO})"
        #             ${AFNI_R} "3dMEMA -prefix ${LOO_MEMA} -mask ${LOO_GROUPMASK} \
        #                 -missing_data 0 -set ${TASK}_mean ${LOO_SET}"
        #         fi
        #         if [ -f "${LOO_MEMA}+tlrc.HEAD" ]; then
        #             # Same threshold -> binarize -> NN-resample -> intersect as STEP A2.
        #             LTHR=$(${AFNI} "p2dsetstat -quiet -inset ${LOO_MEMA}'[1]' -pval ${ACT_PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
        #             if [ -n "${LTHR}" ]; then
        #                 ${AFNI} "
        #                     3dcalc -overwrite -a ${LOO_MEMA}'[1]' -expr 'step(abs(a)-${LTHR})' -prefix ${OUT}/.loo_thr.nii.gz
        #                     3dresample -overwrite -master ${MASK} -rmode NN -prefix ${OUT}/.loo_thr_rs.nii.gz -input ${OUT}/.loo_thr.nii.gz
        #                     3dmask_tool -overwrite -input ${MASK} ${OUT}/.loo_thr_rs.nii.gz -inter -prefix ${ACTMASK_LOO}
        #                 "
        #             else
        #                 echo "  WARN: LOO ROI ${PARTICIPANT_ID}/${TASK} -- empty p2dsetstat threshold" >&2
        #             fi
        #         else
        #             echo "  WARN: LOO MEMA failed for ${PARTICIPANT_ID}/${TASK} (R packages? see prerequisite 2)" >&2
        #         fi
        #     fi
        # fi

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
        # --- LOO ROI (STEP A3): uncomment together with the STEP A3 block above ---
        # NACT_LOO=""
        # if [ -f "${ACTMASK_LOO}" ]; then
        #     NACT_LOO=$(${AFNI} "3dBrickStat -non-zero -count ${ACTMASK_LOO}" 2>/dev/null | tail -n1 | awk '{print $1}')
        # fi
        # if [ -n "${NACT_LOO}" ] && [ "${NACT_LOO}" -gt 0 ] 2>/dev/null; then
        #     MASKS+=("activation_loo:${ACTMASK_LOO}")
        # fi
        MASKS+=("wholebrain:${MASK}")

        # ------------------------------------------------------------------
        # STEP B: per-session cluster-size thresholds for dice_clust. Each
        # session uses its own residual smoothness (sessions can differ in
        # motion/censoring), computed once per task and shared by contrasts.
        # ------------------------------------------------------------------
        E1=${AFNI_OUT}/${PARTICIPANT_ID}/${SES1}/task-${TASK}/${PARTICIPANT_ID}_${SES1}_task-${TASK}_errts_REML+tlrc
        E2=${AFNI_OUT}/${PARTICIPANT_ID}/${SES2}/task-${TASK}/${PARTICIPANT_ID}_${SES2}_task-${TASK}_errts_REML+tlrc
        K1=""; K2=""
        [ -f "${E1}.HEAD" ] && K1=$(session_clust_k "${E1}" "${MASK}" "${OUT}/clustsim_${SES1}")
        [ -f "${E2}.HEAD" ] && K2=$(session_clust_k "${E2}" "${MASK}" "${OUT}/clustsim_${SES2}")
        if [ -z "${K1}" ] || [ -z "${K2}" ]; then
            echo "  WARN: ${TASK} -- no cluster threshold (k1='${K1}' k2='${K2}'; missing errts_REML?); dice_clust=NA" >&2
        else
            echo "  cluster thresholds (p<${CLUST_PTHR}, alpha=${CLUST_ATHR}, NN1 bisided): ${SES1} k=${K1}  ${SES2} k=${K2}"
            add_metric "${PARTICIPANT_ID}" "${TASK}" all tstat clust_k_${SES1} wholebrain "${K1}"
            add_metric "${PARTICIPANT_ID}" "${TASK}" all tstat clust_k_${SES2} wholebrain "${K2}"
        fi

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

            # Cluster-corrected binary maps (dice_clust), whole-brain clustering;
            # the mask is applied afterwards by dice_pair, like the p<.05 maps.
            CL1=${OUT}/clust_${SES1}_${CON}.nii.gz
            CL2=${OUT}/clust_${SES2}_${CON}.nii.gz
            CLUST_OK=0
            if [ -n "${K1}" ] && [ -n "${K2}" ]; then
                if clust_map "${S1}" "${TSTAT}" "${MASK}" "${K1}" "${CL1}" \
                   && clust_map "${S2}" "${TSTAT}" "${MASK}" "${K2}" "${CL2}"; then
                    CLUST_OK=1
                else
                    echo "  WARN: ${CON} -- 3dClusterize failed; dice_clust=NA" >&2
                fi
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

                # Dice overlap of thresholded t-maps (both empty -> 0, see dice_pair)
                DICE=NA
                if [ "${DICE_OK}" -eq 1 ]; then
                    read -r DICE _ _ <<< "$(dice_pair "${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "${M}" "${ACT1}" "${ACT2}")"
                fi
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice "${MTAG}" "${DICE}"

                # Cluster-corrected Dice + suprathreshold voxel counts per session
                # (so empty maps -- no surviving cluster -- are visible in the TSV).
                DICEC=NA; NC1=NA; NC2=NA
                if [ "${CLUST_OK}" -eq 1 ]; then
                    read -r DICEC NC1 NC2 <<< "$(dice_pair "clust ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "${M}" "${CL1}" "${CL2}")"
                fi
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice_clust "${MTAG}" "${DICEC}"
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat clust_nvox_${SES1} "${MTAG}" "${NC1}"
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat clust_nvox_${SES2} "${MTAG}" "${NC2}"

                echo "  ${CON} [${MTAG}]: r(beta)=${R_BETA} r(t)=${R_TSTAT} sICC=${SPICC} dice=${DICE} dice_clust=${DICEC} (nvox ${NC1}/${NC2})"
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
        is_excluded "${sid}" "${SES1}" "${TASK}" >/dev/null && continue
        is_excluded "${sid}" "${SES2}" "${TASK}" >/dev/null && continue
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
        GMASK=${GROUP_OUT}/groupmask_task-${TASK}_${CON}.nii.gz
        rm -f "${GMASK}"
        ${AFNI} "3dmask_tool -overwrite -input ${GROUP_OUT}/inputs/*_task-${TASK}_${CON}_Coef.nii.gz -frac 1.0 -prefix ${GMASK}"
        MEMA_MEAN=${GROUP_MEMA}/mema_${TASK}_mean+tlrc
        if [ -f "${MEMA_MEAN}.HEAD" ] && [ -f "${GMASK}" ]; then
            GATHR=$(${AFNI} "p2dsetstat -quiet -inset ${MEMA_MEAN}'[1]' -pval ${ACT_PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
            if [ -n "${GATHR}" ]; then
                ${AFNI} "
                    3dcalc -overwrite -a ${MEMA_MEAN}'[1]' -expr 'step(abs(a)-${GATHR})' -prefix ${GROUP_OUT}/.gact.nii.gz
                    3dresample -overwrite -master ${GMASK} -rmode NN -prefix ${GROUP_OUT}/.gact_rs.nii.gz -input ${GROUP_OUT}/.gact.nii.gz
                    3dmask_tool -overwrite -input ${GMASK} ${GROUP_OUT}/.gact_rs.nii.gz -inter -prefix ${GROUP_OUT}/groupmask_act_task-${TASK}_${CON}.nii.gz
                "
                [ -f "${GROUP_OUT}/groupmask_act_task-${TASK}_${CON}.nii.gz" ] && GMASK=${GROUP_OUT}/groupmask_act_task-${TASK}_${CON}.nii.gz
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

# =============================================================================
# WITHIN- vs BETWEEN-SUBJECT DICE  (Gorgolewski et al. 2013, NeuroImage 69:231)
# =============================================================================
# A within-subject Dice only means something relative to how much ANY two people
# overlap (everyone activates visual cortex). Following Gorgolewski et al.:
#   within  = subject X ses-1 vs subject X ses-2
#   between = every pair of DIFFERENT subjects, within the same session
#             (ses-1 pairs and ses-2 pairs), later averaged per task
# A task is "reliable" at the individual level if within > between. The test
# (percentile bootstrap of the Harrell-Davis median) is in
# dice_within_vs_between.R, which reads the TSV written here.
#
# All comparisons for a task use ONE common mask (intersection of every
# subject's mask_both; optionally ∩ the group activation ROI), so within and
# between Dice are computed over the same voxels. That is why these within values
# can differ slightly from the per-subject `dice` rows above, which use each
# subject's own mask. Done for both binary-map variants: p05 (act_*) and
# clust (clust_*). Both-empty pairs score 0, as in dice_pair.
# =============================================================================
echo "=============================================================="
echo "WITHIN vs BETWEEN-SUBJECT DICE (Gorgolewski et al. 2013)"
echo "=============================================================="

DICE_WB=${REPRO_OUT}/dice_within_between.tsv
printf 'task\tcontrast\tvariant\tmask\tkind\tsession\tsubjA\tsubjB\tdice\n' > "${DICE_WB}"
WB_DIR=${REPRO_OUT}/group/dice_wb
mkdir -p "${WB_DIR}"

# Runs INSIDE the container: one AFNI process per comparison instead of one
# container start per comparison (thousands of pairs at n=10).
#   input lines : kind session subjA subjB fileA fileB
#   output lines: kind session subjA subjB dice
DICE_BATCH=${WB_DIR}/dice_batch.sh
cat > "${DICE_BATCH}" <<'EOF'
#!/bin/bash
mask=$1
while read -r kind ses sa sb fa fb; do
    if [ ! -f "${fa}" ] || [ ! -f "${fb}" ]; then
        echo "${kind} ${ses} ${sa} ${sb} NA"; continue
    fi
    na=$(3dBrickStat -mask "${mask}" -non-zero -count "${fa}" 2>/dev/null | awk '{print $1}')
    nb=$(3dBrickStat -mask "${mask}" -non-zero -count "${fb}" 2>/dev/null | awk '{print $1}')
    if [ "${na}" = 0 ] && [ "${nb}" = 0 ]; then
        d=0
    else
        d=$(3ddot -dodice -mask "${mask}" "${fa}" "${fb}" 2>/dev/null | awk '{print $1}')
    fi
    echo "${kind} ${ses} ${sa} ${sb} ${d:-NA}"
done < "$2"
EOF

for TASK in "${TASKS[@]}"; do
    WB_SUBS=()
    for d in "${REPRO_OUT}"/sub-*/task-${TASK}; do
        sid=$(basename "$(dirname "${d}")")
        is_excluded "${sid}" "${SES1}" "${TASK}" >/dev/null && continue
        is_excluded "${sid}" "${SES2}" "${TASK}" >/dev/null && continue
        [ -f "${d}/mask_both.nii.gz" ] && WB_SUBS+=("${sid}")
    done
    n=${#WB_SUBS[@]}
    if [ "${n}" -lt 2 ]; then
        echo "within/between ${TASK}: only ${n} subject(s) -> skip (need >=2 for a between pair)"
        continue
    fi
    echo "within/between ${TASK}: n=${n} (${WB_SUBS[*]})"

    # Common whole-brain mask, and its intersection with the group activation ROI.
    WB_MASK=${WB_DIR}/mask_wholebrain_task-${TASK}.nii.gz
    WB_INPUTS=""
    for sid in "${WB_SUBS[@]}"; do WB_INPUTS="${WB_INPUTS} ${REPRO_OUT}/${sid}/task-${TASK}/mask_both.nii.gz"; done
    ${AFNI} "3dmask_tool -overwrite -input ${WB_INPUTS} -frac 1.0 -prefix ${WB_MASK}" >/dev/null 2>&1
    WB_MASKS=("wholebrain:${WB_MASK}")
    MEMA_MEAN=${GROUP_MEMA}/mema_${TASK}_mean+tlrc
    if [ -f "${MEMA_MEAN}.HEAD" ]; then
        WATHR=$(${AFNI} "p2dsetstat -quiet -inset ${MEMA_MEAN}'[1]' -pval ${ACT_PTHR} -2sided" 2>/dev/null | tail -n1 | awk '{print $1}')
        if [ -n "${WATHR}" ]; then
            WA_MASK=${WB_DIR}/mask_activation_task-${TASK}.nii.gz
            ${AFNI} "
                3dcalc -overwrite -a ${MEMA_MEAN}'[1]' -expr 'step(abs(a)-${WATHR})' -prefix ${WB_DIR}/.wact.nii.gz
                3dresample -overwrite -master ${WB_MASK} -rmode NN -prefix ${WB_DIR}/.wact_rs.nii.gz -input ${WB_DIR}/.wact.nii.gz
                3dmask_tool -overwrite -input ${WB_MASK} ${WB_DIR}/.wact_rs.nii.gz -inter -prefix ${WA_MASK}
            " >/dev/null 2>&1
            [ -f "${WA_MASK}" ] && WB_MASKS=("activation:${WA_MASK}" "${WB_MASKS[@]}")
        fi
    fi

    for CON in "${CONTRASTS[@]}"; do
        for VAR in p05 clust; do
            case "${VAR}" in p05) MAPPRE=act ;; clust) MAPPRE=clust ;; esac
            LIST=${WB_DIR}/.pairs_${TASK}_${CON}_${VAR}.txt
            : > "${LIST}"
            for ((i = 0; i < n; i++)); do
                si=${WB_SUBS[i]}
                Di=${REPRO_OUT}/${si}/task-${TASK}
                echo "within ${SES1//-/}-${SES2//-/} ${si} ${si} ${Di}/${MAPPRE}_${SES1}_${CON}.nii.gz ${Di}/${MAPPRE}_${SES2}_${CON}.nii.gz" >> "${LIST}"
                for ((j = i + 1; j < n; j++)); do
                    sj=${WB_SUBS[j]}
                    Dj=${REPRO_OUT}/${sj}/task-${TASK}
                    for SES in "${SES1}" "${SES2}"; do
                        echo "between ${SES//-/} ${si} ${sj} ${Di}/${MAPPRE}_${SES}_${CON}.nii.gz ${Dj}/${MAPPRE}_${SES}_${CON}.nii.gz" >> "${LIST}"
                    done
                done
            done
            for MSPEC in "${WB_MASKS[@]}"; do
                MTAG=${MSPEC%%:*}
                M=${MSPEC#*:}
                ${AFNI} "bash ${DICE_BATCH} ${M} ${LIST}" 2>/dev/null \
                    | awk -v t="${TASK}" -v c="${CON}" -v v="${VAR}" -v m="${MTAG}" \
                        'NF==5 { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", t, c, v, m, $1, $2, $3, $4, $5 }' >> "${DICE_WB}"
            done
            rm -f "${LIST}"
        done
    done
done
echo "  within/between Dice -> ${DICE_WB}"
echo "  test it with: Rscript dice_within_vs_between.R ${DICE_WB}"

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
