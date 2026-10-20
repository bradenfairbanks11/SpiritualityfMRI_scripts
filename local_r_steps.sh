#!/bin/bash
# =============================================================================
# Spirituality fMRI — the R-based AFNI steps, run on the LOCAL machine
#
# The cluster's AFNI container cannot run 3dICC (and the group steps here are
# the R-based ones), so the cluster exports small, self-contained folders and
# this script does the R-dependent modelling locally. Not a Slurm job.
#
# Steps (run all, or name some:  bash local_r_steps.sh icc trr_roi):
#   icc      3dICC ICC(2,1) group maps, from the inputs reproducibility_afni.sh
#            left in afni_reproducibility/group/ (same model and mask choice as
#            the cluster script would have used)
#   trr_vox  3dLMEr -TRR whole-brain trial-level test-retest reliability
#            (Chen et al. 2021, NeuroImage 245:118647). 3dLMEr -TRR requires
#            exactly 2 conditions: trials are split by rating, low (1-2) vs
#            high (3-4). Outputs TRR-ave (reliability of the overall response)
#            and TRR-con (reliability of the high-minus-low difference).
#   trr_roi  TRR (region-level, Bayesian, brms) for a SINGLE condition per
#            task: reliability of the overall response in each group ROI
#   dice     within- vs between-subject Dice test (dice_within_vs_between.R)
#
# Expected layout (copy these from the cluster with rsync):
#   ${LOCAL_ROOT}/afni_trr/                 <- ${TRR_OUT} (trr_export_afni.sh)
#   ${LOCAL_ROOT}/afni_reproducibility/     <- ${REPRO_OUT}: group/ and
#                                              dice_within_between.tsv
# LOCAL_ROOT defaults to ~/spirituality_local; override with
#   SPIRITUALITY_LOCAL=/path bash local_r_steps.sh
#
# PRECONDITION for icc / trr_vox: the R-based AFNI programs must start. On
# 2026-10-20 they segfaulted on this Mac because AFNI's R_io.so was linked
# against R 4.5 while R 4.6 is installed; the preflight below catches that.
# =============================================================================

set -o pipefail

LOCAL_ROOT=${SPIRITUALITY_LOCAL:-${HOME}/spirituality_local}
TRR_DIR=${LOCAL_ROOT}/afni_trr
REPRO_DIR=${LOCAL_ROOT}/afni_reproducibility
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS=${TRR_DIR}/results

# Leave two cores free for the machine.
NCPU=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)
JOBS=${JOBS:-$(( NCPU > 3 ? NCPU - 2 : 1 ))}

ICC_TASKS=(scripture FHS architecture)
ICC_CONTRASTS=(mean_response rating_slope)
# FHS is left out of the voxelwise rating split: at most 12 trials per run and
# 6 of its runs have < 3 high-rated trials (timing files, 2026-10-20).
TRR_VOX_TASKS=(scripture architecture)
TRR_ROI_TASKS=(scripture FHS architecture)
# A subject enters the voxelwise model only with >= MIN_PER_COND trials of
# EACH condition in EACH session (the high-minus-low difference needs both).
MIN_PER_COND=2

# TRR (region-level) sampler settings. AFNI's TRR.R stores the posterior
# correlation draws in a vector pre-sized to 2000 and fills chains*iterations/2
# of them, so fewer than 2000 draws would leave zeros that drag the reported
# TRR toward 0. Keep CHAINS * ITERATIONS / 2 >= 2000.
TRR_CHAINS=4
TRR_ITER=1000

export AFNI_DECONFLICT=OVERWRITE
export AFNI_ENVIRON_WARNINGS=NO

STEPS=("$@")
[ ${#STEPS[@]} -eq 0 ] && STEPS=(icc trr_vox trr_roi dice)
want() { local s; for s in "${STEPS[@]}"; do [ "${s}" = "$1" ] && return 0; done; return 1; }

mkdir -p "${RESULTS}"

if [ $(( TRR_CHAINS * TRR_ITER / 2 )) -lt 2000 ]; then
    echo "FATAL: TRR_CHAINS*TRR_ITER/2 must be >= 2000 (see note above)" >&2; exit 1
fi

# --- Preflight: R-based AFNI programs must start (catches the R_io.so mismatch) ---
if want icc || want trr_vox; then
    for prog in 3dICC 3dLMEr TRR Rscript; do
        if ! command -v "${prog}" >/dev/null 2>&1; then
            echo "FATAL: ${prog} is not on PATH (is ~/abin on PATH in this shell?)" >&2; exit 1
        fi
    done
    if ! 3dLMEr -help >/dev/null 2>&1 || ! 3dICC -help >/dev/null 2>&1; then
        echo "FATAL: AFNI's R programs (3dICC / 3dLMEr) crash on start-up." >&2
        echo "  AFNI R bridge links against: $(otool -L "$(dirname "$(command -v afni)")/R_io.so" 2>/dev/null | awk '/libR.dylib/ {print $1}')" >&2
        echo "  R on PATH                  : $(Rscript --version 2>&1 | head -n1)" >&2
        echo "  Those R versions must match (install the matching R side by side, or" >&2
        echo "  an AFNI build linked to this R). Steps trr_roi and dice do not need it:" >&2
        echo "      bash local_r_steps.sh trr_roi dice" >&2
        exit 1
    fi
fi

# =============================================================================
# icc — 3dICC ICC(2,1), across-subject voxelwise reliability map
# =============================================================================
if want icc; then
    echo "=============================================================="
    echo "3dICC (group ICC(2,1))"
    echo "=============================================================="
    G=${REPRO_DIR}/group
    ICC_SUMMARY=${G}/reproducibility_icc_summary.tsv
    printf 'task\tcontrast\tn\tmetric\tvalue\n' > "${ICC_SUMMARY}"
    for TASK in "${ICC_TASKS[@]}"; do
        for CON in "${ICC_CONTRASTS[@]}"; do
            # Subjects with both sessions' extracted betas (same rule as the cluster).
            DT=""; n=0
            for f1 in "${G}"/inputs/sub-*_ses-1_task-${TASK}_${CON}_Coef.nii.gz; do
                [ -f "${f1}" ] || continue
                sid=$(basename "${f1}"); sid=${sid%%_*}
                f2=${G}/inputs/${sid}_ses-2_task-${TASK}_${CON}_Coef.nii.gz
                [ -f "${f2}" ] || continue
                DT="${DT} ${sid} ses1 ${f1} ${sid} ses2 ${f2}"
                n=$((n + 1))
            done
            if [ "${n}" -lt 2 ]; then
                echo "skip ${TASK}/${CON}: ${n} subject(s) with both sessions"
                continue
            fi
            # Same mask the cluster script would have used: activation ROI if built.
            GMASK=${G}/groupmask_act_task-${TASK}_${CON}.nii.gz
            [ -f "${GMASK}" ] || GMASK=${G}/groupmask_task-${TASK}_${CON}.nii.gz
            if [ ! -f "${GMASK}" ]; then
                echo "skip ${TASK}/${CON}: no group mask in ${G} (re-run reproducibility_afni.sh)" >&2
                continue
            fi
            [ "${n}" -lt 10 ] && echo "  NOTE ${TASK}/${CON}: n=${n} -- voxelwise ICC is unstable at this n" >&2
            ICC_OUT=${G}/ICC_task-${TASK}_${CON}
            rm -f "${ICC_OUT:?}"+*.HEAD "${ICC_OUT:?}"+*.BRIK*
            echo "  3dICC ${TASK} ${CON}: n=${n}, mask $(basename "${GMASK}")"
            if ! 3dICC -prefix "${ICC_OUT}" -jobs "${JOBS}" -mask "${GMASK}" \
                    -model '1+(1|session)+(1|Subj)' \
                    -dataTable Subj session InputFile ${DT} > "${ICC_OUT}.log" 2>&1; then
                echo "  WARN: 3dICC failed for ${TASK}/${CON}; see ${ICC_OUT}.log" >&2
                printf '%s\t%s\t%d\t%s\t%s\n' "${TASK}" "${CON}" "${n}" median_icc NA >> "${ICC_SUMMARY}"
                continue
            fi
            MED=$(3dBrickStat -median -mask "${GMASK}" "${ICC_OUT}+tlrc[0]" 2>/dev/null | tail -n1 | awk '{print $NF}')
            printf '%s\t%s\t%d\t%s\t%s\n' "${TASK}" "${CON}" "${n}" median_icc "${MED:-NA}" >> "${ICC_SUMMARY}"
            echo "  -> ${ICC_OUT}+tlrc  median ICC = ${MED:-NA}"
        done
    done
    echo "  summary: ${ICC_SUMMARY}"
fi

# =============================================================================
# trr_vox — 3dLMEr -TRR, whole brain, low vs high rated trials
# =============================================================================
if want trr_vox; then
    echo "=============================================================="
    echo "3dLMEr -TRR (whole brain; cond = low (1-2) vs high (3-4) rating)"
    echo "=============================================================="
    VOX_SUMMARY=${RESULTS}/trr_vox_summary.tsv
    printf 'task\tn_subj\tn_trials\tmetric\tmask\tmedian\n' > "${VOX_SUMMARY}"
    for TASK in "${TRR_VOX_TASKS[@]}"; do
        TABLE=${TRR_DIR}/trials_task-${TASK}.tsv
        TMASK=${TRR_DIR}/mask_task-${TASK}.nii.gz
        if [ ! -f "${TABLE}" ] || [ ! -f "${TMASK}" ]; then
            echo "skip ${TASK}: no ${TABLE} / mask (run trr_export_afni.sh, copy afni_trr/)"
            continue
        fi
        # 3dLMEr data table. cond is dummy-coded -0.5 (low) / +0.5 (high) and
        # passed with -qVars, as in the 3dLMEr -TRR help example. Subjects
        # missing a condition in a session are dropped and listed.
        DTAB=${RESULTS}/lmer_trr_task-${TASK}.tbl
        awk -F'\t' -v root="${TRR_DIR}" -v minc="${MIN_PER_COND}" -v drop="${RESULTS}/lmer_trr_task-${TASK}_dropped.txt" '
            NR==FNR { if (FNR>1) cnt[$1 SUBSEP $2 SUBSEP $5]++; next }
            FNR==1  { print "Subj sess cond InputFile"; next }
            {
                ok = (cnt[$1,"ses1","low"] >= minc && cnt[$1,"ses1","high"] >= minc &&
                      cnt[$1,"ses2","low"] >= minc && cnt[$1,"ses2","high"] >= minc)
                if (!ok) { bad[$1] = 1; next }
                printf "%s %s %s %s/%s\n", $1, $2, ($5=="high" ? "0.5" : "-0.5"), root, $6
            }
            END { for (s in bad) print s > drop }
        ' "${TABLE}" "${TABLE}" > "${DTAB}"
        NSUB=$(awk 'NR>1 {print $1}' "${DTAB}" | sort -u | wc -l | tr -d ' ')
        NTR=$(($(wc -l < "${DTAB}") - 1))
        if [ -s "${RESULTS}/lmer_trr_task-${TASK}_dropped.txt" ]; then
            echo "  ${TASK}: dropped (< ${MIN_PER_COND} trials of a condition in a session): $(sort "${RESULTS}/lmer_trr_task-${TASK}_dropped.txt" | tr '\n' ' ')"
        fi
        if [ "${NSUB}" -lt 3 ]; then
            echo "skip ${TASK}: only ${NSUB} subject(s) with both conditions in both sessions"
            continue
        fi
        OUTP=${RESULTS}/LMEr_TRR_task-${TASK}
        rm -f "${OUTP:?}"+*.HEAD "${OUTP:?}"+*.BRIK*
        echo "  ${TASK}: n=${NSUB} subjects, ${NTR} trials, ${JOBS} jobs (expect hours)"
        if ! 3dLMEr -prefix "${OUTP}" -TRR -jobs "${JOBS}" -mask "${TMASK}" \
                -qVars cond \
                -model '0+sess+cond:sess+(0+sess|Subj)+(0+cond:sess|Subj)' \
                -dataTable @"${DTAB}" > "${OUTP}.log" 2>&1; then
            echo "  WARN: 3dLMEr -TRR failed for ${TASK}; see ${OUTP}.log" >&2
            continue
        fi
        # Medians of the two reliability maps in the whole mask and in the ROIs.
        # The help warns numerical failures show up as many exact 0 / 1 / -1
        # voxels, so their share is reported too.
        ROIS=${TRR_DIR}/rois_task-${TASK}.nii.gz
        for LAB in TRR-ave TRR-con; do
            # Index lookup, not a [label] selector: AFNI may read the '-' in
            # 'TRR-ave' as a sub-brick range.
            IX=$(3dinfo -label2index "${LAB}" "${OUTP}+tlrc" 2>/dev/null | tail -n1)
            if [ -z "${IX}" ]; then
                echo "  WARN: no ${LAB} sub-brick in ${OUTP}+tlrc (3dinfo -label)" >&2
                continue
            fi
            for MS in "wholebrain:${TMASK}" "rois:${ROIS}"; do
                MT=${MS%%:*}; MF=${MS#*:}
                [ -f "${MF}" ] || continue
                MED=$(3dBrickStat -median -mask "${MF}" "${OUTP}+tlrc[${IX}]" 2>/dev/null | tail -n1 | awk '{print $NF}')
                printf '%s\t%d\t%d\t%s\t%s\t%s\n' "${TASK}" "${NSUB}" "${NTR}" "${LAB}" "${MT}" "${MED:-NA}" >> "${VOX_SUMMARY}"
            done
            DEGEN=$(3dmaskdump -noijk -mask "${TMASK}" "${OUTP}+tlrc[${IX}]" 2>/dev/null \
                | awk '{ n++; if ($1==0 || $1==1 || $1==-1) d++ } END { if (n) printf "%.3f", d/n }')
            printf '%s\t%d\t%d\t%s\t%s\t%s\n' "${TASK}" "${NSUB}" "${NTR}" "${LAB}_frac_0_or_pm1" wholebrain "${DEGEN:-NA}" >> "${VOX_SUMMARY}"
        done
        echo "  -> ${OUTP}+tlrc"
    done
    echo "  summary: ${VOX_SUMMARY}"
fi

# =============================================================================
# trr_roi — TRR (region-level), single condition, one model per task x ROI
# =============================================================================
if want trr_roi; then
    echo "=============================================================="
    echo "TRR (region-level; single condition per task)"
    echo "=============================================================="
    ROI_SUMMARY=${RESULTS}/trr_roi_summary.tsv
    printf 'task\tROI\tn_subj\tn_trials\tTRR_mode\tHDI95_lo\tHDI95_hi\n' > "${ROI_SUMMARY}"
    for TASK in "${TRR_ROI_TASKS[@]}"; do
        RT=${TRR_DIR}/roi_task-${TASK}.tsv
        if [ ! -f "${RT}" ]; then
            echo "skip ${TASK}: no ${RT} (no group ROI survived, or export not copied)"
            continue
        fi
        for ROI in $(awk -F'\t' 'NR>1 {print $6}' "${RT}" | sort -u); do
            TAB=${RESULTS}/trr_task-${TASK}_${ROI}.tbl
            awk -F'\t' -v r="${ROI}" 'NR==1 { print "Subj sess Y"; next } $6==r && $7!="" { print $1, $2, $7 }' "${RT}" > "${TAB}"
            NSUB=$(awk 'NR>1 {print $1}' "${TAB}" | sort -u | wc -l | tr -d ' ')
            NTR=$(($(wc -l < "${TAB}") - 1))
            OUTP=${RESULTS}/TRR_task-${TASK}_${ROI}
            rm -f "${OUTP:?}".txt "${OUTP:?}".RData "${OUTP:?}"-ave.pdf
            echo "  ${TASK} ${ROI}: n=${NSUB}, ${NTR} trials"
            if ! TRR -prefix "${OUTP}" -chains "${TRR_CHAINS}" -iterations "${TRR_ITER}" \
                    -Y Y -subject Subj -repetition sess -dataTable "${TAB}" > "${OUTP}.log" 2>&1; then
                echo "  WARN: TRR failed for ${TASK}/${ROI}; see ${OUTP}.log" >&2
                printf '%s\t%s\t%d\t%d\tNA\tNA\tNA\n' "${TASK}" "${ROI}" "${NSUB}" "${NTR}" >> "${ROI_SUMMARY}"
                continue
            fi
            MODE=$(sed -n 's/^TRR across the two repetitions: mode - \([-0-9.eE]*\);.*/\1/p' "${OUTP}.txt" | tail -n1)
            HDI=$(sed -n 's/^TRR: 95% highest density interval - \[\(.*\)\]/\1/p' "${OUTP}.txt" | tail -n1 | tr -d ' ')
            printf '%s\t%s\t%d\t%d\t%s\t%s\t%s\n' "${TASK}" "${ROI}" "${NSUB}" "${NTR}" \
                "${MODE:-NA}" "${HDI%%,*}" "${HDI##*,}" >> "${ROI_SUMMARY}"
            echo "     TRR mode = ${MODE:-NA}  95% HDI [${HDI:-NA}]"
        done
    done
    echo "  summary: ${ROI_SUMMARY}  (ROI peak coordinates: rois_task-<T>_clusters.txt)"
fi

# =============================================================================
# dice — within- vs between-subject Dice (Gorgolewski et al. 2013)
# =============================================================================
if want dice; then
    echo "=============================================================="
    echo "within vs between-subject Dice"
    echo "=============================================================="
    WB=${REPRO_DIR}/dice_within_between.tsv
    if [ -f "${WB}" ]; then
        Rscript "${SCRIPT_DIR}/dice_within_vs_between.R" "${WB}"
    else
        echo "skip: no ${WB} (run reproducibility_afni.sh, copy afni_reproducibility/)"
    fi
fi

echo "=============================================================="
echo "local R steps done: ${STEPS[*]}"
echo "=============================================================="
