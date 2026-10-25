#!/bin/bash
# ==============================================================================
# check_trial_censoring.sh — HOW MUCH OF EACH TRIAL WAS CENSORED FOR MOTION?
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   A quality-control (QC) check to run BEFORE the single-trial model
#   (first_level_IM_afni.sh), which gives every trial its own regressor.
#   Motion censoring removes volumes. If it removes most of the volumes that
#   cover a trial's response, that trial's beta rests on very little data (and
#   if it removes all of them, the trial cannot be fitted at all).
#
#   For every rated trial this script computes signal_kept, the exact number
#   first_level_IM_afni.sh uses to decide which trials to drop:
#     signal_kept = the fraction of the trial's modelled response that falls on
#                   volumes that survive censoring (1 = untouched, 0 = nothing
#                   left). "Modelled response" = the trial's regressor, i.e. a
#                   boxcar of the task's duration convolved with AFNI's
#                   hemodynamic response (BLOCK). It is measured as the
#                   regressor's sum of squares over the kept volumes divided by
#                   its sum of squares over all volumes, so volumes near the
#                   peak of the response count more than volumes on its edges.
#   It then shows which trials the single-trial model WILL drop (signal_kept not
#   above MIN_SIGNAL_KEPT) and which are low (below FLAG_BELOW), so a cut-off can
#   be chosen before that model is run.
#
#   It changes no data and excludes nothing. It writes one table, plus
#   temporary design files that it deletes.
#
# WHERE IT SITS IN THE PIPELINE
#   Before: first_level_afni.sh  (makes the motion censor file this reads)
#   Then  : first_level_IM_afni.sh  (the single-trial model; same calculation)
#   Not a pipeline step: nothing downstream reads its output.
#
# STEPS, PER RUN (one run = one participant x session x task)
#   1. Skip the run if it is listed in exclusions.tsv, or if one of its input
#      files is missing
#   2. Write the onsets of every rated trial to a temporary timing file
#   3. Build the single-trial design matrix WITHOUT fitting anything
#      (3dDeconvolve -nodata), with and without censoring
#   4. Compute signal_kept for every trial (and the unrated regressor), write
#      one row per trial to the QC table, and print a one-line run summary
#   Then, after all runs: print how many trials will be dropped or are low,
#   and list them.
#
# INPUTS (per run)
#   censor file  : <AFNI_OUT>/sub-NN/ses-N/task-T/sub-NN_ses-N_task-T_motion_censor.1D
#                    one value per volume, 1 = keep, 0 = censored. Made by
#                    first_level_afni.sh STEP 4: a volume is censored when its
#                    framewise displacement (FD, from fMRIPrep) is > 0.3 mm,
#                    and the volume before it is censored too. Same file
#                    first_level_IM_afni.sh uses. (Censor files made before
#                    2026-10-24 used an older, more lenient rule: re-run
#                    first_level_afni.sh first.)
#   scaled BOLD  : <AFNI_OUT>/sub-NN/ses-N/task-T/sub-NN_ses-N_task-T_blurred_scaled.nii.gz
#                    only its header is read, for the TR (the single-trial
#                    model also takes its TR from this header)
#   rated trials : <TIMING_DIR>/byv-<BYV>_ses-N_task-T_ratingAM.1D   ("onset*rating")
#   unrated      : <TIMING_DIR>/byv-<BYV>_ses-N_task-T_unrated.1D    (only if any)
#
# OUTPUTS
#   <AFNI_OUT>/qc_trial_censoring.tsv   one row per rated trial. Columns:
#       sub ses task trial onset rating signal_kept im_model low
#     trial      = 0-based number among ALL rated trials of the run, in timing-
#                  file order (= orig_trial in first_level_IM_afni.sh's
#                  <PREFIX>_IM_trials.tsv)
#     im_model   = "kept" or "dropped: ..." -- what first_level_IM_afni.sh
#                  will do with this trial at its current MIN_SIGNAL_KEPT
#     low        = "LOW" when a kept trial's signal_kept < FLAG_BELOW, else empty
#   On screen: one summary line per run, then a list of every dropped or low
#   trial. The table is rewritten from scratch every time.
#
# HOW TO RUN (login node: one quick AFNI call per run, a few minutes in total)
#   bash check_trial_censoring.sh                  # every participant in AFNI_OUT
#   bash check_trial_censoring.sh sub-04 sub-05    # only these participants
#   FLAG_BELOW=0.7 bash check_trial_censoring.sh   # mark trials below 0.7 as LOW
#
# METHODS SUMMARY (the exact values used below)
#   Censor rule (from first_level_afni.sh): FD > 0.3 mm, plus the volume before
#   Trial regressor: BLOCK(d,1), d = 19 s scripture, 44 s FHS, 14 s
#                    architecture; -stim_times_IM, onsets from the start of the
#                    run (-local_times); drift -polort A. Identical to
#                    first_level_IM_afni.sh STEP 2.
#   signal_kept    : sum of squares of the trial's regressor over uncensored
#                    volumes / over all volumes
#   Dropped        : signal_kept not above MIN_SIGNAL_KEPT (0.3), or onset after
#                    the end of the run -- the same rule as
#                    first_level_IM_afni.sh STEP 3
# ==============================================================================

set -uo pipefail   # -u: stop on an unset variable; pipefail: a pipe fails if
                   # any part of it fails. (No -e: one bad run should not stop
                   # the check for the others.)


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
#                                      weeks unused: first-level outputs
#   /nobackup/archive/usr/bradenf4     raw data (read only) and archived results
AFNI_OUT=/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/afni_firstlvl     # [shared] censor files in; QC table out
TIMING_DIR=/home/bradenf4/spirituality_fmri/timing                                     # [shared] timing files
PARTICIPANTS_TSV=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata/participants.tsv  # [shared]
EXCLUSIONS_TSV=/home/bradenf4/spirituality_fmri/scripts/exclusions.tsv                 # [shared] runs to skip

# Temporary design files go in a fresh folder here and are deleted at the end.
# (A short path on purpose: AFNI silently fails to write design matrices whose
# path is longer than about 250 characters.)
TMP_PARENT=/nobackup/autodelete/usr/bradenf4/work/Spirituality

# --- 1b. Software -------------------------------------------------------------
# AFNI runs inside an Apptainer image. The image only sees folders that are
# explicitly "bound" into it, so all three storage areas are bound.
AFNI_SIF=/apps/afni/26.3.04.sif                                                        # [shared] AFNI 26.3.04
BIND_ARGS="--bind /home/bradenf4:/home/bradenf4 --bind /nobackup/autodelete/usr/bradenf4:/nobackup/autodelete/usr/bradenf4 --bind /nobackup/archive/usr/bradenf4:/nobackup/archive/usr/bradenf4"  # [shared]
APPTAINER_MODULE=apptainer/1.3.6-qycanb2                                               # [shared]

# --- 1c. Analysis parameters --------------------------------------------------
# The drop rule of first_level_IM_afni.sh: a trial is dropped unless its
# signal_kept is ABOVE this value (0.3 since 2026-10-24).
# Must match first_level_IM_afni.sh, so that "im_model" is what it will do.
MIN_SIGNAL_KEPT=0.3                                                                    # [shared]

# Mark a trial "LOW" when its signal_kept is below this. For looking only; it
# does not change what first_level_IM_afni.sh does. Can be changed for one run
# from the command line (see HOW TO RUN).
FLAG_BELOW=${FLAG_BELOW:-0.5}

# Tasks, by their BIDS labels.
TASKS=(scripture FHS architecture)                                                     # [shared]

# Fixed stimulus duration per task, in seconds. These are the median on-screen
# times measured by generate_timing_files.py:
#   scripture 19.3 s | FHS 44.3 s | architecture 14.2 s
block_dur_for() {                                                                      # [shared]
    case "$1" in
        scripture)    echo 19 ;;
        FHS)          echo 44 ;;
        architecture) echo 14 ;;
        *)            echo "ERROR: unknown task $1" >&2; exit 1 ;;
    esac
}


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# byv_for_sub sub-NN  ->  prints that participant's BYV id (e.g. 0AK1RT).
# BIDS labels participants sub-NN (numbered in scan-date order), but the timing
# files are named by the lab's BYV id. participants.tsv maps one to the other via
# its 'participant_id' and 'byv_id' columns.
# (The awk variable is 'want', not 'sub': sub() is a built-in awk function.)
byv_for_sub() {
    awk -F'\t' -v want="$1" '
        NR==1 { for (i=1;i<=NF;i++){ if($i=="participant_id")p=i; if($i=="byv_id")b=i}; next }
        $p==want { print $b; exit }
    ' "${PARTICIPANTS_TSV}"
}

# xmat_signal_kept CENSORED.xmat.1D UNCENSORED.xmat.1D
#   ->  one line per design-matrix column:  label <TAB> signal_kept <TAB> in_run
# Both files are 3dDeconvolve design matrices for the same model: the first
# with censored volumes removed (-x1D), the second with every volume kept
# (-x1D_uncensored). For each column, signal_kept = its sum of squares in the
# censored matrix divided by its sum of squares in the uncensored one, i.e. the
# fraction of that regressor that survives censoring (0 if the column is empty
# even before censoring). in_run = "no" when the column is empty before
# censoring (its response lies past the end of the run), else "yes".
# Column names come from the matrix's "ColumnLabels" header line.
xmat_signal_kept() {                                                                   # [shared]
    awk '
        FNR == 1 { f++ }
        /ColumnLabels/ {
            s = $0; sub(/.*ColumnLabels = "/, "", s); sub(/".*/, "", s)
            n = split(s, lab, ";")
            for (i = 1; i <= n; i++) gsub(/^ +| +$/, "", lab[i])
            next
        }
        /^#/ || NF == 0 { next }
        { for (i = 1; i <= NF; i++) ss[f, i] += $i * $i }
        END {
            for (i = 1; i <= n; i++)
                printf "%s\t%.6f\t%s\n", lab[i], (ss[2, i] > 0 ? ss[1, i] / ss[2, i] : 0), (ss[2, i] > 0 ? "yes" : "no")
        }
    ' "$1" "$2"
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

module load "${APPTAINER_MODULE}"

# Every AFNI command below runs as:  ${AFNI} "<commands>"
AFNI="apptainer exec ${BIND_ARGS} ${AFNI_SIF} bash -c"

if [ ! -f "${PARTICIPANTS_TSV}" ]; then
    echo "ERROR: participants.tsv not found at ${PARTICIPANTS_TSV}"
    exit 1
fi

# Temporary folder for the design files; removed when the script ends.
mkdir -p "${TMP_PARENT}"
TMP=$(mktemp -d "${TMP_PARENT}/trialcens.XXXXXX")
trap 'rm -rf "${TMP:?}"' EXIT

# Participants: the ones named on the command line, or everyone in AFNI_OUT.
if [ "$#" -gt 0 ]; then
    PARTICIPANT_DIRS=()
    for arg in "$@"; do PARTICIPANT_DIRS+=("${AFNI_OUT}/sub-${arg#sub-}"); done
else
    PARTICIPANT_DIRS=("${AFNI_OUT}"/sub-*)
fi

OUT=${AFNI_OUT}/qc_trial_censoring.tsv
printf 'sub\tses\ttask\ttrial\tonset\trating\tsignal_kept\tim_model\tlow\n' > "${OUT}"

echo "Trial censoring check: signal_kept per trial (dropped if not above ${MIN_SIGNAL_KEPT}; LOW if below ${FLAG_BELOW})"
printf '%-8s %-6s %-13s %6s %9s %7s %4s %8s\n' sub ses task trials run_cens% dropped LOW unrated


# ==============================================================================
# SECTION 4. MAIN LOOP: participants -> sessions -> tasks
# ==============================================================================

for PDIR in "${PARTICIPANT_DIRS[@]}"; do
    [ -d "${PDIR}" ] || { echo "  ${PDIR} not found -- skipping"; continue; }
    SUB=$(basename "${PDIR}")
    BYV=$(byv_for_sub "${SUB}")
    [ -n "${BYV}" ] || { echo "  ${SUB}: no byv_id in participants.tsv -- skipping"; continue; }

    for SES_DIR in "${PDIR}"/ses-*/; do
        [ -d "${SES_DIR}" ] || continue
        SES=$(basename "${SES_DIR}")
        for TASK in "${TASKS[@]}"; do

            # ------------------------------------------------------------------
            # STEP 1. Skip excluded runs and runs with missing inputs.
            # ------------------------------------------------------------------
            if why=$(is_excluded "${SUB}" "${SES}" "${TASK}"); then
                printf '%-8s %-6s %-13s  EXCLUDED (%s)\n' "${SUB}" "${SES}" "${TASK}" "${why:0:40}"
                continue
            fi
            PREFIX=${SUB}_${SES}_task-${TASK}
            RUN_DIR=${AFNI_OUT}/${SUB}/${SES}/task-${TASK}
            CENSOR=${RUN_DIR}/${PREFIX}_motion_censor.1D
            SCALED=${RUN_DIR}/${PREFIX}_blurred_scaled.nii.gz
            TIMING=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_ratingAM.1D
            UNRATED=${TIMING_DIR}/byv-${BYV}_${SES}_task-${TASK}_unrated.1D
            [ -f "${CENSOR}" ] && [ -f "${SCALED}" ] && [ -f "${TIMING}" ] || continue
            BLOCK_DUR=$(block_dur_for "${TASK}")

            # ------------------------------------------------------------------
            # STEP 2. Onsets of every rated trial, in timing-file order, on one
            # line (the "onset*rating" tokens with the rating removed). Same
            # trial order and numbering as first_level_IM_afni.sh STEP 1.
            # ------------------------------------------------------------------
            ALLTIMES=${TMP}/${PREFIX}_alltimes.1D
            tr ' \t' '\n\n' < "${TIMING}" | awk -F'*' '/\*/ { printf "%s%s", (n++ ? " " : ""), $1 } END { print "" }' > "${ALLTIMES}"

            # ------------------------------------------------------------------
            # STEP 3. Build the design matrix only (no data are fitted).
            #
            # 3dDeconvolve -nodata NVOL TR builds the design for a run of NVOL
            # volumes with that TR without reading any data. NVOL = the number
            # of lines in the censor file; TR = the scaled BOLD's header
            # (3dinfo -tr). The stimulus part is identical to
            # first_level_IM_afni.sh STEP 2; motion regressors are left out
            # because they do not change a trial regressor's own values.
            # -x1D           = design with censored volumes removed
            # -x1D_uncensored = the same design with every volume kept
            # ------------------------------------------------------------------
            NVOL=$(wc -l < "${CENSOR}" | tr -d ' ')
            UNRATED_OPTS="-num_stimts 1"
            if [ -f "${UNRATED}" ]; then
                UNRATED_OPTS="-num_stimts 2 -stim_times 2 ${UNRATED} 'BLOCK(${BLOCK_DUR},1)' -stim_label 2 unrated"
            fi
            ${AFNI} "
                cd ${TMP}
                3dDeconvolve \
                    -nodata ${NVOL} \$(3dinfo -tr ${SCALED}) \
                    -polort A \
                    -censor ${CENSOR} \
                    -local_times \
                    ${UNRATED_OPTS} \
                    -stim_times_IM 1 ${ALLTIMES} 'BLOCK(${BLOCK_DUR},1)' -stim_label 1 trial \
                    -x1D ${TMP}/${PREFIX}_x.1D \
                    -x1D_uncensored ${TMP}/${PREFIX}_xu.1D \
                    -x1D_stop
            " > /dev/null 2>&1
            if [ ! -f "${TMP}/${PREFIX}_x.1D" ] || [ ! -f "${TMP}/${PREFIX}_xu.1D" ]; then
                echo "  ${SUB} ${SES} ${TASK}: could not build the design matrix -- skipping"
                continue
            fi

            # ------------------------------------------------------------------
            # STEP 4. signal_kept per trial, the table rows, and the run summary.
            #
            # Same drop rule as first_level_IM_afni.sh STEP 3:
            #   no column at all (onset after the end of the run) -> dropped
            #   signal_kept not above MIN_SIGNAL_KEPT              -> dropped
            # ------------------------------------------------------------------
            KEPT_TABLE=$(xmat_signal_kept "${TMP}/${PREFIX}_x.1D" "${TMP}/${PREFIX}_xu.1D")
            UNRATED_KEPT=$(printf '%s\n' "${KEPT_TABLE}" | awk -F'\t' '$1 == "unrated#0" { print $2 }')
            RUN_CENS=$(awk '$1 + 0 == 0 { c++ } END { printf "%.1f", (NR ? 100 * c / NR : 0) }' "${CENSOR}")

            printf '%s\n' "${KEPT_TABLE}" | awk -F'\t' \
                -v min="${MIN_SIGNAL_KEPT}" -v low="${FLAG_BELOW}" -v out="${OUT}" \
                -v vsub="${SUB}" -v vses="${SES}" -v vtask="${TASK}" \
                -v cens="${RUN_CENS}" -v unrated="${UNRATED_KEPT:--}" '
                NR == FNR { if ($1 ~ /^trial#/) { k = substr($1, 7) + 0; kept[k] = $2; inrun[k] = $3 }; next }
                $0 ~ /\*/ {
                    split($0, p, "*"); k = n++
                    # With -nodata, AFNI keeps an empty column for an onset after
                    # the end of the run (in_run = no); the real model leaves
                    # that trial out. Both count as "after the end of the run".
                    if (!(k in kept) || inrun[k] == "no") { v = ""; im = "dropped: onset is after the end of the run" }
                    else { v = kept[k]; im = (v + 0 > min ? "kept" : (v + 0 > 0 ? "dropped: signal_kept not above " min : "dropped: response fully censored")) }
                    flag = (im == "kept" && v + 0 < low ? "LOW" : "")
                    if (im != "kept") ndrop++
                    if (flag == "LOW") nlow++
                    printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n", vsub, vses, vtask, k, p[1], p[2], v, im, flag >> out
                }
                END { printf "%-8s %-6s %-13s %6d %8s%% %7d %4d %8s\n", vsub, vses, vtask, n, cens, ndrop, nlow, unrated }
            ' - <(tr ' \t' '\n\n' < "${TIMING}")
        done
    done
done


# ==============================================================================
# SECTION 5. SUMMARY
# ==============================================================================

NALL=$(($(wc -l < "${OUT}") - 1))
NDROP=$(awk -F'\t' 'NR>1 && $8 != "kept"' "${OUT}" | wc -l | tr -d ' ')
NLOW=$(awk -F'\t' 'NR>1 && $9 == "LOW"' "${OUT}" | wc -l | tr -d ' ')
echo
echo "${NALL} rated trials checked."
echo "  ${NDROP} will be DROPPED by first_level_IM_afni.sh (MIN_SIGNAL_KEPT = ${MIN_SIGNAL_KEPT})"
echo "  ${NLOW} have signal_kept below ${FLAG_BELOW} (LOW)"
if [ "$((NDROP + NLOW))" -gt 0 ]; then
    echo "Dropped or LOW trials:"
    awk -F'\t' 'NR>1 && ($8 != "kept" || $9 == "LOW") {
        printf "  %s %s %-12s trial %-3s onset %7.1fs  signal_kept %-8s %s\n", $1, $2, $3, $4, $5, ($7 == "" ? "-" : $7), $8 }' "${OUT}"
fi
echo "Column 'unrated' above = signal_kept of the unrated regressor (- = run has none)."
echo "Full table: ${OUT}"
