#!/bin/bash
# ==============================================================================
# SLURM JOB SETTINGS
# ------------------------------------------------------------------------------
# Slurm reads these #SBATCH lines BEFORE the script runs, so they must be literal
# values (shell variables do not expand here).
#   8 cores / 32 GB / 4 h : one job loops over every participant and task in
#                           turn, then runs the group steps. The 8 cores are
#                           used by 3dClustSim (STEP 5) and 3dICC (STEP 11),
#                           which read the core count from SLURM_CPUS_PER_TASK.
#   --chdir               : run inside scratch, so any stray relative-path file
#                           lands there instead of in the home directory.
# ==============================================================================
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
#SBATCH --chdir=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# ==============================================================================
# reproducibility_afni.sh — TEST-RETEST RELIABILITY OF THE FIRST-LEVEL MAPS
#                           (session 1 vs session 2)
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   Every participant is scanned twice (ses-1 and ses-2) on the same three
#   tasks. This script asks: how similar are a person's brain maps from the two
#   sessions? It does NOT fit any new GLM. It reads the first-level results made
#   by first_level_afni.sh, so the GLM settings are identical for both sessions
#   by construction.
#
#   It measures reliability at two levels:
#
#   LEVEL 1 — WITHIN ONE PERSON (per participant x task x contrast):
#     spatial_r   Pearson correlation between the two sessions' maps, over all
#                 voxels in a mask (3ddot -docor). Done on the beta map and on
#                 the t-statistic map. It measures whether the PATTERN is the
#                 same; it ignores a constant shift or scaling between sessions.
#     spatial_icc Intraclass correlation ICC(2,1), absolute agreement, of the
#                 same paired voxel betas (voxels = the things measured, the 2
#                 sessions = the two measurements). Unlike spatial_r it is
#                 lowered by a shift or scaling between sessions. Computed with
#                 the R package irr (spatial_icc.R).
#     dice        Overlap of the two sessions' t-maps, each thresholded at
#                 p < 0.05 two-sided (3ddot -dodice). Do the same voxels pass
#                 threshold in both sessions?
#     dice_clust  The same overlap, but each session's t-map is first
#                 cluster-corrected (voxel p < 0.001, cluster size from that
#                 session's own 3dClustSim, cluster-level alpha 0.05).
#     For every Dice: if NEITHER map has any voxels, the Dice is set to 0
#     (Gorgolewski et al. 2013, NeuroImage 69:231-243).
#     Also written: the paired session betas (for scatter plots) and a
#     session-2-minus-session-1 difference map.
#
#   LEVEL 2 — ACROSS PEOPLE (per task x contrast):
#     group ICC   A voxelwise ICC(2,1) map from AFNI's 3dICC, pooling every
#                 participant who has both sessions (participants = the things
#                 measured, the 2 sessions = the two measurements).
#     NOTE: this ICC measures how reliably people's DIFFERENCES from each other
#     repeat across sessions. A region where everyone activates strongly and
#     identically has little between-person variance, and so a LOW ICC, even
#     though the activation is strong. High ICC does not mean "the task
#     reliably activates this region" (that is the group activation map). With
#     few two-session participants the map is unstable; the script prints a
#     loud warning when n < 10.
#
#   WITHIN- vs BETWEEN-PERSON DICE (per task x contrast):
#     A within-person Dice only means something compared with how much ANY two
#     people overlap. Following Gorgolewski et al. (2013), the script also
#     computes the Dice for every pair of DIFFERENT participants within the same
#     session, and writes all within and between values to one table. The test
#     (within > between) is done afterwards by dice_within_vs_between.R.
#
#   MASKS. Every Level-1 metric is computed twice, inside two masks:
#     activation  voxels significant (p < 0.05, two-sided) in the GROUP
#                 mean_response map from second_level_afni.sh, inside the
#                 whole-brain mask. Because the group map pools all participants
#                 and averages the two sessions, it does not depend on the
#                 ses-1 vs ses-2 comparison being measured. Only mean_response
#                 (the robust contrast) is used to define it, never the noisier
#                 rating_slope. If the group map is missing, this mask is
#                 skipped and only the whole-brain mask is used.
#     wholebrain  voxels the first-level GLM actually modelled in BOTH
#                 sessions: non-zero mean_response AND rating_slope betas in
#                 ses-1 and in ses-2 (see STEP 1).
#     The summary table's `mask` column says which one each row used.
#     A third, DISABLED mask (activation_loo, "leave one participant out") is
#     written out in STEP 3; see the notes there.
#
#   CONTRASTS compared (sub-brick labels in each first-level REML bucket):
#     mean_response   mean_response#0_Coef, mean_response#0_Tstat  (task activation)
#     rating_slope    rating_slope#0_Coef,  rating_slope#0_Tstat   (BOLD change per
#                                                                   rating point)
#
#   All participants and sessions are on the same MNI152NLin2009cAsym grid, so
#   maps are compared voxel by voxel with no resampling. (The group activation
#   map is still resampled to each mask's grid with nearest neighbour, to be
#   safe.)
#
# WHERE IT SITS IN THE PIPELINE
#   Before: first_level_afni.sh   (the per-session GLM results compared here)
#           second_level_afni.sh  (the group mean_response map used for the
#                                  activation mask; optional)
#   After : reliability_models_afni.sh dice  (runs dice_within_vs_between.R,
#                                             which tests within- vs
#                                             between-person Dice)
#
# STEPS
#   LEVEL 1, for each participant x task (SECTION 4):
#     1. Whole-brain mask: voxels that are brain in both sessions
#     2. Activation mask from the group mean_response map
#     3. (DISABLED) leave-one-participant-out activation mask
#     4. Make the list of masks to use (activation if built, then wholebrain)
#     5. Per-session minimum cluster size for dice_clust (3dFWHMx -> 3dClustSim)
#     For each contrast:
#     6. Make the binary maps: p < 0.05 t-maps and cluster-corrected t-maps
#     7. In each mask: spatial_r, spatial_icc, dice, dice_clust
#     8. Write the paired betas (scatter) and the ses-2 minus ses-1 map
#   LEVEL 2, for each task x contrast (SECTION 5):
#     9.  Find every participant with both sessions
#     10. Copy each participant's beta maps out to NIfTI
#     11. Build the group mask and run 3dICC; record the median ICC
#   WITHIN vs BETWEEN DICE, for each task x contrast (SECTION 6):
#     12. Write the small helper that computes many Dice values in one go
#     13. Find participants and build the common masks
#     14. List every within and between pair and compute their Dice
#
# INPUTS
#   first-level betas/t : <AFNI_OUT>/sub-NN/ses-N/task-T/sub-NN_ses-N_task-T_stats_REML+tlrc
#   first-level residuals: <AFNI_OUT>/sub-NN/ses-N/task-T/sub-NN_ses-N_task-T_errts_REML+tlrc
#                         (for dice_clust only)
#   group mean map       : <GROUP_OUT>/mema/mema_T_mean+tlrc   (3dMEMA, from
#                          second_level_afni.sh STEP 3; sub-brick [1] = its
#                          t-statistic; optional)
#   spatial_icc.R        : in the folder sbatch was run from (see SCRIPT_DIR)
#   participants.tsv, exclusions.tsv (runs to skip)
#   (<AFNI_OUT>, <GROUP_OUT>: see SECTION 1a)
#
# OUTPUTS, in <REPRO_OUT> = /nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_reproducibility
#   reproducibility_summary.tsv      THE LEVEL-1 RESULT. One row per
#                                    participant / task / contrast / map /
#                                    metric / mask. Columns:
#                                    participant task contrast map metric mask value
#   reproducibility_icc_summary.tsv  median group ICC per task / contrast
#   dice_within_between.tsv          every within- and between-person Dice
#   spatial_icc.R                    copy of the R helper the job actually ran
#   sub-NN/task-T/                   per-participant working files:
#     mask_ses-1.nii.gz, mask_ses-2.nii.gz, mask_both.nii.gz  (whole-brain masks)
#     mask_activation.nii.gz                                  (activation mask)
#     clustsim_ses-N_acf.1D, clustsim_ses-N.NN*_*.1D          (3dFWHMx / 3dClustSim)
#     act_ses-N_<contrast>.nii.gz      p < 0.05 binary t-maps
#     clust_ses-N_<contrast>.nii.gz    cluster-corrected binary t-maps
#     scatter_<contrast>.txt           paired betas: i j k ses1 ses2
#     ses2_minus_ses1_<contrast>_Coef.nii.gz
#   group/                           LEVEL-2 files:
#     inputs/sub-NN_ses-N_task-T_<contrast>_Coef.nii.gz
#     groupmask_task-T_<contrast>.nii.gz, groupmask_act_task-T_<contrast>.nii.gz
#     ICC_task-T_<contrast>+tlrc       THE LEVEL-2 RESULT (voxelwise ICC map)
#   group/dice_wb/                   within/between Dice files:
#     dice_batch.sh, mask_wholebrain_task-T.nii.gz, mask_activation_task-T.nii.gz
#
# HOW TO RUN (from /home/bradenf4/spirituality_fmri/scripts, so that
#             spatial_icc.R is found next to this script)
#   sbatch reproducibility_afni.sh                 # every sub-* in AFNI_OUT
#   sbatch reproducibility_afni.sh sub-04 sub-05   # only these participants
#   The arguments only limit LEVEL 1. The group ICC (SECTION 5) and the
#   within/between Dice (SECTION 6) always use every participant on disk.
#   Re-running replaces that participant's rows in reproducibility_summary.tsv.
#   Afterwards:  sbatch reliability_models_afni.sh dice
#                (runs dice_within_vs_between.R on dice_within_between.tsv)
#
# METHODS SUMMARY (the exact values used below)
#   Data           : first-level 3dREMLfit outputs (first_level_afni.sh) for
#                    ses-1 and ses-2; contrasts mean_response and rating_slope;
#                    runs listed in exclusions.tsv are skipped
#   Whole-brain    : voxels with non-zero first-level betas (mean_response and
#   mask             rating_slope) in BOTH sessions, i.e. exactly the voxels
#                    the GLM modelled with data in both (complete-case voxels;
#                    Sochat et al. 2015 Front Neurosci found this more accurate
#                    for image similarity than treating zeros as values)
#   Activation     : group one-sample 3dMEMA mean_response map
#   mask             (second_level_afni.sh), |t| above the p < 0.05 two-sided
#                    threshold (p2dsetstat), nearest-neighbour resampled, and
#                    intersected with the whole-brain mask. Both positive and
#                    negative effects are included.
#   spatial_r      : Pearson r of the two sessions' maps within the mask
#                    (3ddot -docor), on betas and on t-statistics
#   spatial_icc    : ICC(2,1), two-way random, absolute agreement, single
#                    measure (irr::icc), voxels as targets, on betas
#   dice           : each session's t-map thresholded at p < 0.05 two-sided
#                    (threshold from that map's own degrees of freedom);
#                    Dice = 2|A and B| / (|A| + |B|) within the mask
#                    (3ddot -dodice); both maps empty -> Dice = 0
#   dice_clust     : voxel threshold p < 0.001, "bisided" (positive and
#                    negative clusters formed separately), clusters of
#                    face-touching voxels (NN1) of at least k voxels kept
#                    (3dClusterize). k is per session: 3dFWHMx -acf on that
#                    session's REML residuals within the whole-brain mask gives
#                    the spatial autocorrelation (ACF) parameters, and
#                    3dClustSim -acf gives the smallest cluster size with
#                    cluster-level alpha < 0.05, rounded up. Clustering uses
#                    the whole-brain mask; Dice is then computed within each
#                    mask; both empty -> 0 (Gorgolewski et al. 2013)
#   Group ICC      : 3dICC, model 1+(1|session)+(1|Subj) = ICC(2,1), on beta
#                    maps of every participant with both sessions (needs n >= 2);
#                    mask = voxels non-zero in every input, intersected with the
#                    activation mask (p < 0.05) when available; the median ICC
#                    inside the mask is reported
#   Within vs      : within = same participant, ses-1 vs ses-2; between = every
#   between Dice     pair of different participants within ses-1 and within
#                    ses-2; all computed within one common mask per task
#                    (intersection of every participant's whole-brain mask, and
#                    that intersected with the activation mask); both for the
#                    p < 0.05 and the cluster-corrected maps. Tested in
#                    dice_within_vs_between.R (percentile bootstrap of the
#                    Harrell-Davis median of within minus mean between)
# ==============================================================================

# This script deliberately does NOT use `set -e` (stop at the first failure).
# One AFNI or container failure on one participant/task must not stop the whole
# batch: each task checks its own results (see the mask check in STEP 1 and the
# threshold checks in STEP 6) and skips ahead if something failed.
# `pipefail` makes a pipeline such as `cmd | tail | awk` count as failed when
# any part of it fails, not only the last part.
set -o pipefail


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
AFNI_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_firstlvl     # [shared] first-level results (input)
REPRO_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_reproducibility  # [shared] output
PARTICIPANTS_TSV=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata/participants.tsv  # [shared]
EXCLUSIONS_TSV=/home/bradenf4/spirituality_fmri/scripts/exclusions.tsv                 # [shared] runs to skip
LOGS=/home/bradenf4/spirituality_fmri/logs                                             # [shared]

# Group results from second_level_afni.sh. Only its 3dMEMA folder is read here
# (GROUP_MEMA): the one-sample mean_response map mema_<task>_mean+tlrc, whose
# sub-brick [0] is the group effect and [1] its t-statistic. It defines the
# activation mask. If it is missing, the script uses the whole-brain mask only.
# CAUTION: SECTION 5 later REUSES the name GROUP_OUT for a different folder,
# <REPRO_OUT>/group. GROUP_MEMA is set here, before that, so it is unaffected.
GROUP_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_group      # [shared]
GROUP_MEMA=${GROUP_OUT}/mema

# The two summary tables this script writes (see OUTPUTS in the header).
SUMMARY=${REPRO_OUT}/reproducibility_summary.tsv          # Level 1, one row per metric
ICC_SUMMARY=${REPRO_OUT}/reproducibility_icc_summary.tsv  # Level 2, one row per task/contrast

# Folder that holds this script and its R helper spatial_icc.R.
# Under sbatch, $0 and BASH_SOURCE point at Slurm's private copy of the script,
# not at the scripts folder, so the folder sbatch was run from
# (SLURM_SUBMIT_DIR) is used instead. Run directly with bash, it is the script's
# own folder. So submit the job from /home/bradenf4/spirituality_fmri/scripts.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# --- 1b. Software -------------------------------------------------------------
# AFNI runs inside an Apptainer image. The image only sees folders that are
# explicitly "bound" into it, so all three storage areas are bound.
# The image also carries R and every R package used here (irr for spatial_icc;
# lme4, blme, metafor, snow for 3dICC; data.table for 3dMEMA). Compute nodes
# have no internet, so nothing can be installed from inside a job.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]

# R packages that must be inside the image (checked in SECTION 3):
# spatial ICC needs irr; 3dICC needs lme4, blme, metafor and snow.
ICC_PKGS='c("irr","lme4","blme","metafor","snow")'

# --- 1c. Analysis parameters --------------------------------------------------
# Tasks, by their BIDS labels.
TASKS=(scripture FHS architecture)                                                     # [shared]

# Contrasts to compare: the label prefixes of the sub-bricks in each REML bucket
# (<label>#0_Coef = beta, <label>#0_Tstat = t-statistic).
CONTRASTS=(mean_response rating_slope)

# The two sessions compared (test = ses-1, retest = ses-2).
SES1=ses-1
SES2=ses-2

# Voxelwise threshold for the plain `dice` metric: p < 0.05, two-sided.
PTHR=0.05

# Cluster-corrected Dice (metric dice_clust). Each session's t-map is
# thresholded at voxel p < CLUST_PTHR (bisided) and only clusters of at least k
# voxels are kept. k comes from 3dClustSim run on THAT session's own residual
# smoothness (3dFWHMx -acf on its errts_REML), for cluster-level alpha
# CLUST_ATHR. p < .001 is the recommended upper limit for the cluster-forming
# p-value (Woo, Krishnan & Wager 2014, NeuroImage 91:412-419; Eklund et al. 2016
# PNAS showed p = .01 inflates false positives badly). The ACF model is AFNI's
# post-Eklund fix (Cox et al. 2017, Brain Connectivity 7:152).
CLUST_PTHR=0.001
CLUST_ATHR=0.05

# Threshold for the activation mask: p < 0.05, two-sided, on the group
# mean_response map.
ACT_PTHR=0.05


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
# If exclusions.tsv or the participant's BYV id cannot be found, the run is
# treated as NOT excluded.
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

# add_metric participant task contrast map metric mask value
#   Appends one row to the Level-1 summary table (SUMMARY).
add_metric() {  # args: participant task contrast map metric mask value
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$7" >> "${SUMMARY}"
}

# afni_scalar "description" "AFNI command"  ->  prints one number.
#   Runs an AFNI command that prints a single number on its last line, and
#   prints the first word of that line. If nothing comes back (AFNI/container
#   error, unknown option, missing sub-brick), it warns and prints "NA", so the
#   gap shows up in the table as NA instead of a silent blank cell.
afni_scalar() {  # args: description  afni_command_string
    local desc="$1" cmd="$2" val
    val=$(${AFNI} "${cmd}" 2>/dev/null | tail -n1 | awk '{print $1}')
    if [ -z "${val}" ]; then
        echo "  WARN: ${desc} -- AFNI produced no value; recording NA" >&2
        val=NA
    fi
    printf '%s' "${val}"
}

# spatial_icc "description" mask session1_beta session2_beta  ->  prints the ICC.
#   1. 3dmaskdump writes the two sessions' beta values, one row per voxel inside
#      the mask (-noijk: values only, no voxel coordinates), to a temp file.
#   2. spatial_icc.R (run with R inside the image) computes ICC(2,1), absolute
#      agreement, single measure, with irr::icc, treating voxels as targets.
#   Prints NA on any failure (R helper missing, empty mask, R error).
spatial_icc() {  # args: description  mask  s1_coef_selector  s2_coef_selector
    local desc="$1" mask="$2" a="$3" b="$4" pairs val
    if [ -z "${ICC_R}" ]; then printf 'NA'; return; fi
    # The script runs one comparison at a time, so one reusable temp file in the
    # current output folder (OUT) is enough. (This avoids mktemp, whose
    # templates differ between systems.)
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

# dice_pair "description" mask map_a map_b  ->  prints "dice n_a n_b".
#   map_a and map_b are binary (0/1) maps. n_a and n_b are the number of
#   non-zero voxels of each inside the mask (3dBrickStat -non-zero -count), and
#   dice is their Dice overlap inside the mask (3ddot -dodice).
#   If BOTH maps are empty the Dice is set to 0, not left undefined: the
#   convention of Gorgolewski et al. (2013, NeuroImage 69:231-243), "to penalise
#   for lack of signal". Prints "NA NA NA" if AFNI returned no counts.
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

# session_clust_k residuals mask prefix  ->  prints k, the minimum cluster size
# (in voxels) for one session.
#   1. 3dFWHMx -detrend -acf fits the spatial autocorrelation of that session's
#      REML residuals inside the mask, using only the uncensored volumes
#      (<PREFIX>_motion_censor.1D), as afni_proc.py does,
#      ACF(r) = a*exp(-r^2/(2b^2)) + (1-a)*exp(-r/c),
#      and prints the parameters a b c on its last line. (It also writes the
#      ACF curve to <prefix>_acf.1D.)
#   2. 3dClustSim -acf a b c simulates noise with that smoothness inside the
#      mask and writes tables of the cluster size needed for cluster-level
#      alpha CLUST_ATHR at voxel p CLUST_PTHR (<prefix>.NN1_bisided.1D, etc.).
#      It uses SLURM_CPUS_PER_TASK cores (4 if that is not set).
#   3. k is read from the NN1 bisided table (first data row, second column)
#      and rounded up to a whole voxel.
#   The tables are re-made only if they do not exist yet, the residuals are
#   newer than them, or they were made with an older smoothness method
#   (<prefix>.method). Prints nothing if 3dFWHMx fails.
session_clust_k() {  # args: errts_dataset(no .HEAD)  mask  prefix
    local errts="$1" mask="$2" pre="$3" acf a b c k censor keep sel
    local method="detrend+uncensored"
    if [ ! -f "${pre}.NN1_bisided.1D" ] || [ "${errts}.HEAD" -nt "${pre}.NN1_bisided.1D" ] \
       || [ "$(cat "${pre}.method" 2>/dev/null)" != "${method}" ]; then
        # Uncensored volumes only (censored ones are 0 in errts), as afni_proc.py does.
        censor=${errts%_errts_REML+tlrc}_motion_censor.1D
        keep=""
        [ -f "${censor}" ] && keep=$(${AFNI} "1d_tool.py -infile ${censor} -show_trs_uncensored encoded" 2>/dev/null | tail -n1)
        sel=""; [ -n "${keep}" ] && sel="'[${keep}]'"
        acf=$(${AFNI} "3dFWHMx -detrend -acf ${pre}_acf.1D -mask ${mask} -input ${errts}${sel}" 2>/dev/null | tail -n1)
        read -r a b c _ <<< "${acf}"
        if [ -z "${c}" ]; then
            echo "  WARN: 3dFWHMx -acf failed on ${errts}; no cluster threshold" >&2
            return
        fi
        ${AFNI} "OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} 3dClustSim -mask ${mask} -acf ${a} ${b} ${c} \
            -pthr ${CLUST_PTHR} -athr ${CLUST_ATHR} -quiet -prefix ${pre}" >/dev/null 2>&1
        printf '%s\n' "${method}" > "${pre}.method"
    fi
    k=$(awk '!/^#/ && NF>=2 { k=$2; print (k==int(k)) ? k : int(k)+1; exit }' "${pre}.NN1_bisided.1D" 2>/dev/null)
    printf '%s' "${k}"
}

# clust_map stats_bucket t_label mask k out.nii.gz
#   Writes a binary (0/1) cluster-corrected map. 3dClusterize thresholds the
#   t sub-brick at voxel p < CLUST_PTHR, "bisided" (positive and negative
#   voxels form separate clusters; the p-value is converted to a t threshold
#   using the bucket's own degrees of freedom), joins face-touching voxels into
#   clusters (-NN 1), and keeps clusters of at least k voxels. 3dcalc
#   step(a) turns the cluster-number map into 0/1.
#   -outvol_if_no_clust writes an all-zero map when no cluster survives, so
#   "no surviving clusters" is a value (Dice 0 or a count of 0), not a missing
#   file. Succeeds only if the output file exists.
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
# SECTION 3. SETUP AND CHECKS
# ==============================================================================

mkdir -p "${LOGS}" "${REPRO_OUT}"

module load "${APPTAINER_MODULE}"

# Every AFNI command below runs as:  ${AFNI} "<commands>"
# AFNI refuses to overwrite existing files by default; AFNI_DECONFLICT=OVERWRITE
# (passed into the image with --env) lets a re-run replace its own files.
AFNI="apptainer exec --env AFNI_DECONFLICT=OVERWRITE ${BIND_ARGS} ${AFNI_SIF} bash -c"
# The R-based steps (spatial_icc.R, 3dICC) run as:  ${AFNI_R} "<commands>"
# The R packages ship INSIDE the AFNI image. --cleanenv starts the image with an
# empty environment, so a host setting (e.g. a stray R_LIBS) cannot override
# those packages. Every variable these commands need is filled in on the host
# before the command reaches the image, and --env still applies under --cleanenv.
AFNI_R="apptainer exec --cleanenv --env AFNI_DECONFLICT=OVERWRITE ${BIND_ARGS} ${AFNI_SIF} bash -c"

# Check that the image has the R packages listed in ICC_PKGS (SECTION 1b).
# The R one-liner exits with the number of missing packages (0 = all present).
# Failing here is clearer than an obscure R error deep inside 3dICC.
if ! ${AFNI_R} "Rscript -e 'q(status=length(setdiff(${ICC_PKGS}, rownames(installed.packages()))))'"; then
    echo "ERROR: the ICC R packages are missing from ${AFNI_SIF}."
    echo "       (need: irr, lme4, blme, metafor, snow)"
    echo "Check that AFNI_SIF in config.sh points at the R-enabled image (see container_request.md)."
    exit 1
fi

# Copy the spatial-ICC R helper into REPRO_OUT, a bound folder the image can
# read; this also keeps a record of the exact helper the job ran. If it is not
# found in SCRIPT_DIR, every spatial_icc value is recorded as NA.
ICC_R=${REPRO_OUT}/spatial_icc.R
if [ -f "${SCRIPT_DIR}/spatial_icc.R" ]; then
    cp -f "${SCRIPT_DIR}/spatial_icc.R" "${ICC_R}"
else
    echo "WARNING: spatial_icc.R not found in ${SCRIPT_DIR}; spatial_icc will record NA" >&2
    ICC_R=""
fi

# Start the Level-1 summary table (header row) if it does not exist yet. An
# existing table is kept, so rows for participants not re-run this time stay.
if [ ! -f "${SUMMARY}" ]; then
    printf 'participant\ttask\tcontrast\tmap\tmetric\tmask\tvalue\n' > "${SUMMARY}"
fi


# ==============================================================================
# SECTION 4. LEVEL 1 — WITHIN-PERSON RELIABILITY
#            loop: participants -> tasks   (ses-1 vs ses-2 compared in each)
# ==============================================================================

# Participants: the ones named on the command line ("sub-05" or "05"), or every
# sub-* folder in the first-level output.
if [ "$#" -gt 0 ]; then
    PARTICIPANT_DIRS=()
    for arg in "$@"; do
        PARTICIPANT_DIRS+=("${AFNI_OUT}/sub-${arg#sub-}")
    done
else
    echo "WARNING: no participant labels given -> defaulting to ALL sub-* in ${AFNI_OUT}"
    PARTICIPANT_DIRS=("${AFNI_OUT}"/sub-*)
fi

# Count the participant x task comparisons that actually ran. If none ran
# (e.g. the container is broken on this node), the script exits with an error
# after this loop, so the FAIL email is sent instead of an empty "success".
PROCESSED=0

for PARTICIPANT_DIR in "${PARTICIPANT_DIRS[@]}"; do
    if [ ! -d "${PARTICIPANT_DIR}" ]; then
        echo "ERROR: ${PARTICIPANT_DIR} not found in first-level output -- skipping."
        continue
    fi
    PARTICIPANT_ID=$(basename "${PARTICIPANT_DIR}")

    # Re-runs replace, not duplicate: delete this participant's old rows from
    # the summary table (keeping the header) before new ones are added.
    if [ -f "${SUMMARY}" ]; then
        awk -F'\t' -v s="${PARTICIPANT_ID}" 'NR==1 || $1!=s' "${SUMMARY}" > "${SUMMARY}.tmp" \
            && mv "${SUMMARY}.tmp" "${SUMMARY}"
    fi

    for TASK in "${TASKS[@]}"; do

        # --- This task's input files, and checks ------------------------------
        # First-level REML buckets for the two sessions (AFNI +tlrc HEAD/BRIK).
        S1=${AFNI_OUT}/${PARTICIPANT_ID}/${SES1}/task-${TASK}/${PARTICIPANT_ID}_${SES1}_task-${TASK}_stats_REML+tlrc
        S2=${AFNI_OUT}/${PARTICIPANT_ID}/${SES2}/task-${TASK}/${PARTICIPANT_ID}_${SES2}_task-${TASK}_stats_REML+tlrc

        # Skip if either session's results are missing ...
        if [ ! -f "${S1}.HEAD" ] || [ ! -f "${S2}.HEAD" ]; then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- missing a session's REML bucket (need both ${SES1} and ${SES2})"
            continue
        fi
        # ... or if either session is listed in exclusions.tsv.
        if why=$(is_excluded "${PARTICIPANT_ID}" "${SES1}" "${TASK}") \
           || why=$(is_excluded "${PARTICIPANT_ID}" "${SES2}" "${TASK}"); then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- EXCLUDED (exclusions.tsv): ${why}"
            continue
        fi

        OUT=${REPRO_OUT}/${PARTICIPANT_ID}/task-${TASK}
        mkdir -p "${OUT}"
        echo "=== ${PARTICIPANT_ID} task-${TASK}: comparing ${SES1} vs ${SES2} ==="

        # ------------------------------------------------------------------
        # STEP 1. Whole-brain mask = voxels the GLM modelled in BOTH sessions.
        #
        # The first-level GLM only fits voxels that are inside fMRIPrep's
        # brain mask AND have data in the tedana output (first_level_afni.sh
        # STEP 2 sets everything else to 0), so every beta outside that set is
        # exactly 0. Those voxels must not enter the comparison: a 0 in both
        # sessions would count as perfect agreement in spatial_r and
        # spatial_icc, and would also enter the smoothness estimate and the
        # cluster-size simulation (STEP 5).
        #
        # So, for each session, mask_ses-N = voxels where BOTH the
        # mean_response and the rating_slope beta are non-zero
        # (step(abs(a)) = 1 for any non-zero value). mask_both = voxels in
        # both sessions' masks (3dmask_tool -inter). This is the "wholebrain"
        # mask for every metric. Because it is built from the betas
        # themselves, it can never disagree with what was modelled.
        # (Until 2026-10-24 this mask was 3dAutomask of each session's tedana
        # BOLD, which did not match the GLM's mask.)
        # ------------------------------------------------------------------
        MASK=${OUT}/mask_both.nii.gz
        rm -f "${OUT:?}"/mask_${SES1}.nii.gz "${OUT:?}"/mask_${SES2}.nii.gz "${MASK:?}"
        ${AFNI} "
            3dcalc -overwrite -a ${S1}'[mean_response#0_Coef]' -b ${S1}'[rating_slope#0_Coef]' \
                -expr 'step(abs(a))*step(abs(b))' -datum byte -prefix ${OUT}/mask_${SES1}.nii.gz
            3dcalc -overwrite -a ${S2}'[mean_response#0_Coef]' -b ${S2}'[rating_slope#0_Coef]' \
                -expr 'step(abs(a))*step(abs(b))' -datum byte -prefix ${OUT}/mask_${SES2}.nii.gz
            3dmask_tool -overwrite -input ${OUT}/mask_${SES1}.nii.gz ${OUT}/mask_${SES2}.nii.gz -inter -prefix ${MASK}
        "
        echo "  whole-brain mask: $(${AFNI} "3dBrickStat -count -non-zero ${MASK}" 2>/dev/null | tail -n1 | awk '{print $1}') voxels (non-zero betas in both sessions)"

        # If the mask was not made, AFNI or the container failed on this node
        # (this is what killed the 2026-07-20 run: the compute node had no
        # /dev/fuse). Skip this task instead of stopping the whole batch.
        if [ ! -f "${MASK}" ]; then
            echo "skip: ${PARTICIPANT_ID} task-${TASK} -- mask build failed (AFNI/container error on this node?)"
            continue
        fi
        PROCESSED=$((PROCESSED + 1))

        # ------------------------------------------------------------------
        # STEP 2. Activation mask from the GROUP mean_response map.
        #
        # Source: mema_<task>_mean+tlrc from second_level_afni.sh, sub-brick
        # [1] = group t-statistic. Because it pools all participants, it does
        # not depend on either session's own threshold (avoids "double
        # dipping"), and it uses only the robust mean_response contrast.
        #   p2dsetstat     converts p < ACT_PTHR (two-sided) into a t value,
        #                  using the map's own degrees of freedom
        #   3dcalc         step(abs(a)-thr): 1 where |t| > thr (positive OR
        #                  negative effects), else 0
        #   3dresample     put it on this mask's grid, nearest neighbour
        #   3dmask_tool    keep only voxels also in the whole-brain mask
        # If the group map is missing, no activation mask is made (STEP 4 then
        # uses the whole-brain mask only).
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
        # STEP 3 (DISABLED). LEAVE-ONE-PARTICIPANT-OUT activation mask
        #
        # WHY: the STEP 2 mask comes from a group map that INCLUDES this
        # participant, so about 1/n of the evidence that selects the voxels
        # comes from the person being measured. That is a weak form of double
        # dipping. second_level_afni.sh (its STEP 1) averages the two sessions
        # before 3dMEMA, so the mask cannot see the ses-1 vs ses-2 DIFFERENCE,
        # but it does see their SUM. Selecting voxels with a large |b1+b2|
        # favours voxels where both sessions happened to be large together,
        # which nudges spatial_r / spatial_icc upward. Leaving the participant
        # out removes this: sub-X's mask is built from everyone EXCEPT sub-X, so
        # the noise that selects voxels and the noise being measured are
        # independent.
        #
        # COST: the leave-one-out map has n-1 participants, so its degrees of
        # freedom drop and the threshold rises (at n=4: t(2)=4.303 vs
        # t(3)=3.182). Each participant also gets a DIFFERENT mask, so
        # activation_loo values are not strictly comparable across people --
        # keep the shared-mask `activation` rows as the comparable set and read
        # the leave-one-out rows as the bias-free companion.
        #
        # DO NOT ENABLE BELOW ~8-10 PARTICIPANTS. With few participants the
        # leave-one-out mask is so unstable that you trade a known bias for
        # unknown noise. LOO_MIN_N enforces this.
        #
        # PREREQUISITES before uncommenting:
        #   1. second_level_afni.sh must keep its per-participant
        #      session-averaged inputs (they currently exist only in autodelete
        #      scratch, which is purged). In its STEP 7 change
        #          for d in mask mema ttest clustsim; do
        #      to
        #          for d in mask mema ttest clustsim inputs; do
        #      and re-run it.
        #   2. 3dMEMA needs R packages data.table + snow (both in the AFNI image)
        #      and must be called through ${AFNI_R}, not ${AFNI}.
        #   3. rm -rf ${REPRO_OUT}/loo_mema whenever the set of participants
        #      changes: the leave-one-out maps are cached per (task, left-out
        #      participant) and go stale.
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
        #             # Same threshold -> binarize -> NN-resample -> intersect as STEP 2.
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

        # ------------------------------------------------------------------
        # STEP 4. The list of masks to compute every metric in.
        #
        # Each entry is "tag:file". The activation mask is used only if STEP 2
        # made it AND it has at least one voxel (3dBrickStat -non-zero -count);
        # the whole-brain mask is always used.
        # ------------------------------------------------------------------
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
        # --- Leave-one-out mask (STEP 3): uncomment together with STEP 3 above ---
        # NACT_LOO=""
        # if [ -f "${ACTMASK_LOO}" ]; then
        #     NACT_LOO=$(${AFNI} "3dBrickStat -non-zero -count ${ACTMASK_LOO}" 2>/dev/null | tail -n1 | awk '{print $1}')
        # fi
        # if [ -n "${NACT_LOO}" ] && [ "${NACT_LOO}" -gt 0 ] 2>/dev/null; then
        #     MASKS+=("activation_loo:${ACTMASK_LOO}")
        # fi
        MASKS+=("wholebrain:${MASK}")

        # ------------------------------------------------------------------
        # STEP 5. Minimum cluster size k for each session (for dice_clust).
        #
        # Each session uses its own residual smoothness, because sessions can
        # differ in motion and censoring (see session_clust_k in SECTION 2).
        # Computed once per task and shared by both contrasts. Both k values
        # are written to the summary table as metric clust_k_ses-N
        # (contrast "all", map "tstat", mask "wholebrain").
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

        # --- Per contrast: STEPS 6-8 ------------------------------------------
        for CON in "${CONTRASTS[@]}"; do
            COEF="${CON}#0_Coef"
            TSTAT="${CON}#0_Tstat"

            # ------------------------------------------------------------------
            # STEP 6. The binary maps used for the two Dice metrics.
            #
            # (a) p < PTHR maps, for `dice` (act_ses-N_<contrast>.nii.gz):
            #     p2dsetstat turns p < 0.05 two-sided into a t threshold, using
            #     each session's own degrees of freedom; 3dcalc then marks
            #     voxels with |t| above it as 1. These maps are not masked; the
            #     mask is applied later by 3ddot -dodice. Made once per contrast.
            # (b) cluster-corrected maps, for `dice_clust`
            #     (clust_ses-N_<contrast>.nii.gz): see clust_map in SECTION 2.
            #     Clustering is done in the whole-brain mask; each mask from
            #     STEP 4 is applied afterwards by dice_pair, as for (a).
            # ------------------------------------------------------------------
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

            # ------------------------------------------------------------------
            # STEP 7. The reliability metrics, inside each mask from STEP 4.
            #
            #   spatial_r   (beta and t)  3ddot -docor, Pearson r -> consistency
            #   spatial_icc (beta)        ICC(2,1) absolute agreement, via irr
            #   dice        (t)           Dice of the p < 0.05 maps
            #   dice_clust  (t)           Dice of the cluster-corrected maps,
            #                             plus each session's voxel count inside
            #                             the mask (clust_nvox_ses-N), so empty
            #                             maps (no surviving cluster) are
            #                             visible in the table
            # Both empty -> Dice 0 (see dice_pair). A metric that could not be
            # computed is recorded as NA.
            # ------------------------------------------------------------------
            for MSPEC in "${MASKS[@]}"; do
                MTAG=${MSPEC%%:*}
                M=${MSPEC#*:}

                # spatial correlation (Pearson r)
                R_BETA=$(afni_scalar  "r(beta) ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}"  "3ddot -docor -mask ${M} ${S1}'[${COEF}]'  ${S2}'[${COEF}]'")
                R_TSTAT=$(afni_scalar "r(t) ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}"     "3ddot -docor -mask ${M} ${S1}'[${TSTAT}]' ${S2}'[${TSTAT}]'")
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" beta  spatial_r "${MTAG}" "${R_BETA}"
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat spatial_r "${MTAG}" "${R_TSTAT}"

                # spatial ICC(2,1), absolute agreement (beta)
                SPICC=$(spatial_icc "${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "${M}" "${S1}'[${COEF}]'" "${S2}'[${COEF}]'")
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" beta spatial_icc "${MTAG}" "${SPICC}"

                # Dice of the p < 0.05 t-maps
                DICE=NA
                if [ "${DICE_OK}" -eq 1 ]; then
                    read -r DICE _ _ <<< "$(dice_pair "${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "${M}" "${ACT1}" "${ACT2}")"
                fi
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice "${MTAG}" "${DICE}"

                # Dice of the cluster-corrected t-maps + voxel counts per session
                DICEC=NA; NC1=NA; NC2=NA
                if [ "${CLUST_OK}" -eq 1 ]; then
                    read -r DICEC NC1 NC2 <<< "$(dice_pair "clust ${PARTICIPANT_ID}/${TASK}/${CON}/${MTAG}" "${M}" "${CL1}" "${CL2}")"
                fi
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat dice_clust "${MTAG}" "${DICEC}"
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat clust_nvox_${SES1} "${MTAG}" "${NC1}"
                add_metric "${PARTICIPANT_ID}" "${TASK}" "${CON}" tstat clust_nvox_${SES2} "${MTAG}" "${NC2}"

                echo "  ${CON} [${MTAG}]: r(beta)=${R_BETA} r(t)=${R_TSTAT} sICC=${SPICC} dice=${DICE} dice_clust=${DICEC} (nvox ${NC1}/${NC2})"
            done

            # ------------------------------------------------------------------
            # STEP 8. Files for looking at the data (whole-brain mask only).
            #
            # scatter_<contrast>.txt: 3dmaskdump writes one row per voxel in
            #   the whole-brain mask with columns  i j k ses1 ses2  (voxel
            #   indices, then the two sessions' betas), ready for a scatter plot.
            # ses2_minus_ses1_<contrast>_Coef.nii.gz: the beta difference map
            #   (ses-2 minus ses-1), for visual inspection.
            # ------------------------------------------------------------------
            SCATTER=${OUT}/scatter_${CON}.txt
            rm -f "${SCATTER}"
            ${AFNI} "3dmaskdump -mask ${MASK} -o ${SCATTER} ${S1}'[${COEF}]' ${S2}'[${COEF}]'"

            ${AFNI} "3dcalc -overwrite -a ${S2}'[${COEF}]' -b ${S1}'[${COEF}]' -expr 'a-b' -prefix ${OUT}/ses2_minus_ses1_${CON}_Coef.nii.gz"
        done
    done
done

# Stop with an error if no comparison ran at all (see PROCESSED above). Note
# that this also skips SECTIONS 5 and 6.
if [ "${PROCESSED}" -eq 0 ]; then
    echo "ERROR: no subject/task comparison ran (all skipped or failed). Not a successful run." >&2
    exit 1
fi

echo "Within-subject reproducibility complete (${PROCESSED} task-comparisons). Summary: ${SUMMARY}"


# ==============================================================================
# SECTION 5. LEVEL 2 — GROUP VOXELWISE ICC(2,1) WITH 3dICC
# ------------------------------------------------------------------------------
# Pools every participant who has BOTH sessions for a task (participants = the
# things measured, the 2 sessions = the two measurements) and fits, in every
# voxel, the two-way random-effects model  1+(1|session)+(1|Subj),  which is
# ICC(2,1). Output: one ICC map per task x contrast, plus the median ICC inside
# the mask in reproducibility_icc_summary.tsv.
#
# Inputs are the beta (Coef) maps, copied out to NIfTI and labelled as MNI space
# (as second_level_afni.sh does). This section ignores the participant
# arguments: the group is everyone on disk with both sessions.
#
# REMEMBER: high ICC = people's DIFFERENCES from each other repeat across
# sessions, NOT "the task reliably activates" (that is the group activation
# map). With few participants this map is very unstable (see the warning below).
# ==============================================================================
echo "=============================================================="
echo "GROUP ICC(2,1): voxelwise 3dICC across two-session subjects"
echo "=============================================================="

# From here on GROUP_OUT means THIS script's group folder, <REPRO_OUT>/group
# (not second_level_afni.sh's afni_group folder from SECTION 1a).
GROUP_OUT=${REPRO_OUT}/group
mkdir -p "${GROUP_OUT}/inputs"

# Start a fresh group summary table every run (it is recomputed from everyone
# with both sessions).
printf 'task\tcontrast\tn\tmetric\tvalue\n' > "${ICC_SUMMARY}"

for TASK in "${TASKS[@]}"; do

    # ------------------------------------------------------------------
    # STEP 9. Find every participant with both sessions for this task.
    #
    # A participant is included if both sessions' REML buckets exist and
    # neither session is listed in exclusions.tsv. At least 2 are needed;
    # with fewer than 10 a warning is printed.
    # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # STEP 10. Copy each participant x session beta map out to NIfTI.
        #
        # 3dbucket copies the <contrast>#0_Coef sub-brick to its own file;
        # 3drefit -space MNI labels it as MNI space. Each file is added to the
        # 3dICC data table (DT) as:  participant  session  file
        # (the session label is written ses1 / ses2, without the dash).
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # STEP 11. Group mask, 3dICC, and the median ICC.
        #
        # Group mask: voxels that are non-zero in EVERY input file in
        # group/inputs for this task and contrast (3dmask_tool -frac 1.0),
        # then intersected with the group activation mask (same p < 0.05
        # recipe as STEP 2) when the group mean_response map exists.
        # 3dICC: model 1+(1|session)+(1|Subj) = ICC(2,1), using
        # SLURM_CPUS_PER_TASK cores (4 if not set). If it fails, NA is recorded.
        # Median: 3dBrickStat -median of sub-brick [0] inside the mask.
        # ------------------------------------------------------------------
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

        # Median ICC inside the mask. The code assumes 3dICC's sub-brick [0] is
        # the ICC and [1] its F-statistic; confirm once with: 3dinfo -label ICC_*+tlrc
        MED=$(${AFNI} "3dBrickStat -median -mask ${GMASK} ${ICC_OUT}+tlrc'[0]'" 2>/dev/null | tail -n1 | awk '{print $NF}')
        [ -z "${MED}" ] && MED=NA
        printf '%s\t%s\t%d\t%s\t%s\n' "${TASK}" "${CON}" "${n}" median_icc "${MED}" >> "${ICC_SUMMARY}"
        echo "  ICC map -> ${ICC_OUT}+tlrc   median ICC (in mask) = ${MED}"
    done
done


# ==============================================================================
# SECTION 6. WITHIN- vs BETWEEN-PERSON DICE  (Gorgolewski et al. 2013,
#            NeuroImage 69:231)
# ------------------------------------------------------------------------------
# A within-person Dice only means something compared with how much ANY two
# people overlap (everyone activates visual cortex, for example). Following
# Gorgolewski et al.:
#   within  = participant X ses-1 vs participant X ses-2
#   between = every pair of DIFFERENT participants, within the same session
#             (ses-1 pairs and ses-2 pairs)
# A task is reliable at the individual level if within > between. The test
# (dice_within_vs_between.R: between values averaged per task/contrast/variant/
# mask, then a percentile bootstrap of the Harrell-Davis median of within minus
# between) reads the table written here.
#
# All comparisons for a task use ONE common mask (the intersection of every
# participant's mask_both, and optionally that intersected with the group
# activation mask), so within and between Dice cover the same voxels. That is
# why these within values can differ slightly from the per-participant `dice`
# rows in SECTION 4, which use each participant's own mask. Done for both kinds
# of binary map from STEP 6: p05 (act_*) and clust (clust_*). Both-empty pairs
# score 0, as in dice_pair. This section uses every participant with SECTION 4
# output on disk, not only the ones named on the command line.
# ==============================================================================
echo "=============================================================="
echo "WITHIN vs BETWEEN-SUBJECT DICE (Gorgolewski et al. 2013)"
echo "=============================================================="

DICE_WB=${REPRO_OUT}/dice_within_between.tsv
printf 'task\tcontrast\tvariant\tmask\tkind\tsession\tsubjA\tsubjB\tdice\n' > "${DICE_WB}"
WB_DIR=${REPRO_OUT}/group/dice_wb
mkdir -p "${WB_DIR}"

# ------------------------------------------------------------------
# STEP 12. Write dice_batch.sh, a small helper that runs INSIDE the image.
#
# There are thousands of pairs at n=10. Starting the container once per pair
# would be slow, so the whole list of pairs is handed to this helper, which
# runs once per list. For each pair it counts the voxels of each map inside
# the mask (3dBrickStat), sets Dice = 0 if both are empty, and otherwise
# computes the Dice (3ddot -dodice); a missing map gives NA.
#   arguments   : mask  pair_list
#   input lines : kind session subjA subjB fileA fileB
#   output lines: kind session subjA subjB dice
# ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # STEP 13. Participants and common masks for this task.
    #
    # Participants: every <REPRO_OUT>/sub-*/task-<task> folder that has a
    # mask_both.nii.gz (i.e. STEP 1 ran for it, in this run or an earlier
    # one), minus any run listed in exclusions.tsv. At least 2 are needed.
    # Masks: "wholebrain" = voxels in EVERY participant's mask_both
    # (3dmask_tool -frac 1.0); "activation" = that, intersected with the
    # group activation mask (same p < 0.05 recipe as STEP 2), when the group
    # mean_response map exists.
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # STEP 14. List every pair and compute its Dice.
    #
    # For each contrast and each map variant (p05 = act_*, clust = clust_*):
    #   - one "within" line per participant (ses-1 vs ses-2), session label
    #     "ses1-ses2"
    #   - one "between" line per pair of different participants, for ses-1
    #     and again for ses-2
    # The list goes to dice_batch.sh once per mask, and each result line is
    # appended to dice_within_between.tsv with the task, contrast, variant
    # and mask in front.
    # ------------------------------------------------------------------
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


# ==============================================================================
# SECTION 7. FINAL REPORT
# ==============================================================================
echo "=============================================================="
echo "Reliability run complete."
echo "  within-subject metrics : ${SUMMARY}   (mask = activation | wholebrain)"
echo "  group ICC maps         : ${GROUP_OUT}/ICC_task-*_*+tlrc"
echo "  group ICC medians      : ${ICC_SUMMARY}"
echo "  NOTE: group ICC needs many subjects to be meaningful; small-n results"
echo "        above are placeholders and will stabilize as the sample grows."
echo "=============================================================="

# ------------------------------------------------------------------------------
# OPTIONAL (DISABLED): LEVEL-1 ROLL-UP, mean +/- SD across participants
# ------------------------------------------------------------------------------
# Averages reproducibility_summary.tsv across participants into one row per
# (task, contrast, map, metric, mask). Uncomment once several participants
# exist. The value is column 7; the grouping keys are columns 2-6 (mask
# included). NA cells are skipped.
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
