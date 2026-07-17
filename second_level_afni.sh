#!/bin/bash
#SBATCH --job-name=spirituality_grouplvl
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bradenfairbanks@gmail.com

# =============================================================================
# Spirituality fMRI — SECOND-LEVEL (GROUP) ANALYSIS (AFNI)
#
# Aggregates the first-level parametric-modulation results across subjects.
# The first level (first_level_afni.sh) fit one 3dDeconvolve/3dREMLfit GLM per
# subject/session/task with -stim_times_AM2, yielding per GLM:
#     rating_slope  = parametric SLOPE of BOLD vs the 1-4 "felt the Spirit"
#                     rating  (PRIMARY result of interest)
#     mean_response = overall activation (sanity / localizer)
#
# This script:
#   STEP 0  Discover which <sub>/<ses>/task-<TASK>_stats_REML+orig exist.
#   STEP 1  Extract slope + mean Coef/Tstat to NIfTI (MNI), and for subjects
#           with >1 session, fixed-effects-average their sessions so each
#           subject contributes ONE value per task.
#   STEP 2  Build a group brain mask.
#   STEP 3  3dMEMA one-sample per task  (slope primary, mean as sanity).
#   STEP 4  3dMEMA between-task slope contrasts.
#   STEP 5  3dttest++ (same one-sample + paired between-task) with -Clustsim.
#   STEP 6  Cluster correction for the 3dMEMA maps (3dFWHMx -> 3dClustSim).
#
# The subject list is NOT hard-coded: whatever first-level output exists is
# picked up automatically. Re-run after sub-02/sub-04 (and more sessions)
# finish and the group grows on its own.
#
# NOTE on space: the first-level BOLD was warped to MNI152NLin2009cAsym so all
# subjects share a grid, but AFNI mislabels the datasets +orig. Extracting to
# NIfTI writes the true MNI sform from the grid; 3drefit -space MNI152_2009
# then sets the label so clustering / whereami work.
#
# Run with:  sbatch second_level_afni.sh
# =============================================================================

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PROJECT=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project
AFNI_FIRSTLVL=${PROJECT}/derivatives/afni_firstlvl
FMRIPREP_OUT=${PROJECT}/derivatives/fmriprep

# Compute in fast autodelete scratch (archive is slow / discouraged for compute),
# then persist the small final maps back to archive at the end (STEP 7).
ARCHIVE_ROOT=/nobackup/archive/usr/bradenf4
AUTODELETE_ROOT=/nobackup/autodelete/usr/bradenf4
GROUP_OUT=${AUTODELETE_ROOT}/spirituality_group      # working dir (fast scratch)
GROUP_FINAL=${PROJECT}/derivatives/afni_group        # permanent home (archive)

# 3dMEMA is an R program; the AFNI container's R (3.4.4) lacks its packages
# (data.table, snow). We install them once into this persistent library on
# archive and expose it to the container via R_LIBS_USER (see bootstrap below).
RLIB=${GROUP_FINAL}/Rlib

AFNI_SIF=/apps/afni/afni_make_build_latest.sif

# --- Tasks ---
TASKS=(scripture FHS architecture)

# --- Between-task slope contrasts (X:Y  ->  X minus Y) ---
BETWEEN_TASK_PAIRS=(scripture:architecture FHS:architecture scripture:FHS)

# --- Sub-brick labels inside <PREFIX>_stats_REML (GLT results) ---
SLOPE_COEF="rating_slope#0_Coef"      # primary
SLOPE_T="rating_slope#0_Tstat"
MEAN_COEF="mean_response#0_Coef"      # sanity
MEAN_T="mean_response#0_Tstat"

# Uncorrected voxelwise p used to derive the ClustSim cluster-size table.
CLUST_PTHR=0.001

# CPUs for the -Clustsim / 3dClustSim randomization (SLURM alloc, else 4).
NTHR=${SLURM_CPUS_PER_TASK:-4}

# ==============================================================================

module load apptainer/1.3.6-qycanb2
mkdir -p "${RLIB}"
# Bind both archive (source data) and autodelete (fast scratch) into the
# container, and expose the persistent R library so 3dMEMA finds its packages.
AFNI="apptainer exec --bind ${ARCHIVE_ROOT}:${ARCHIVE_ROOT} --bind ${AUTODELETE_ROOT}:${AUTODELETE_ROOT} --env R_LIBS_USER=${RLIB} ${AFNI_SIF} bash -c"

# --- Bootstrap: install 3dMEMA's R packages into RLIB on first run ---
if ! ${AFNI} "Rscript -e 'q(status=length(setdiff(c(\"data.table\",\"snow\"), rownames(installed.packages()))))'"; then
    echo "Installing 3dMEMA R packages (snow, data.table) into ${RLIB} ..."
    ${AFNI} "Rscript -e 'options(repos=c(CRAN=\"https://cloud.r-project.org\")); install.packages(c(\"snow\",\"data.table\"), lib=Sys.getenv(\"R_LIBS_USER\"))'"
fi

mkdir -p logs "${GROUP_OUT}"/{inputs,mask,mema,ttest,clustsim}
INPUTS=${GROUP_OUT}/inputs
MASKDIR=${GROUP_OUT}/mask
MEMADIR=${GROUP_OUT}/mema
TTESTDIR=${GROUP_OUT}/ttest
CLUSTDIR=${GROUP_OUT}/clustsim
GROUPMASK=${MASKDIR}/group_mask.nii.gz

# ==============================================================================
# STEP 0 — Discover available first-level stats
# Build SESS["<sub> <task>"] = "ses-1 ses-2 ..."  and a subject list.
# ==============================================================================
echo "=============================================================="
echo "STEP 0: discovering first-level stats under ${AFNI_FIRSTLVL}"
echo "=============================================================="

declare -A SESS
SUBJECTS=()

for HEAD in "${AFNI_FIRSTLVL}"/sub-*/ses-*/task-*/*_stats_REML+orig.HEAD; do
    [ -e "${HEAD}" ] || continue
    # .../sub-XX/ses-Y/task-TASK/sub-XX_ses-Y_task-TASK_stats_REML+orig.HEAD
    task_dir=$(dirname "${HEAD}")
    ses_dir=$(dirname "${task_dir}")
    sub_dir=$(dirname "${ses_dir}")
    SUB=$(basename "${sub_dir}")
    SES=$(basename "${ses_dir}")
    TASK=$(basename "${task_dir}"); TASK=${TASK#task-}

    key="${SUB} ${TASK}"
    SESS["${key}"]="${SESS[${key}]} ${SES}"
    if [[ ! " ${SUBJECTS[*]} " == *" ${SUB} "* ]]; then
        SUBJECTS+=("${SUB}")
    fi
done

IFS=$'\n' SUBJECTS=($(sort <<<"${SUBJECTS[*]}")); unset IFS

if [ "${#SUBJECTS[@]}" -eq 0 ]; then
    echo "ERROR: no *_stats_REML+orig datasets found under ${AFNI_FIRSTLVL}"
    exit 1
fi

echo "Discovered subjects: ${SUBJECTS[*]}"
for SUB in "${SUBJECTS[@]}"; do
    for TASK in "${TASKS[@]}"; do
        s="${SESS[${SUB} ${TASK}]}"
        [ -n "${s}" ] && echo "  ${SUB}  task-${TASK}  ->${s}"
    done
done

# Path to a first-level stats dataset for a given sub/ses/task.
stats_dset() {  # <sub> <ses> <task>
    echo "${AFNI_FIRSTLVL}/$1/$2/task-$3/$1_$2_task-$3_stats_REML+orig"
}

# ==============================================================================
# STEP 1 — Extract + session-average -> one Coef/Tstat per subject x task
# For each effect writes:
#   inputs/<sub>_task-<TASK>_<effect>_Coef.nii.gz
#   inputs/<sub>_task-<TASK>_<effect>_T.nii.gz
# ==============================================================================
echo "=============================================================="
echo "STEP 1: extracting + session-averaging per subject x task"
echo "=============================================================="

# Extract one labelled sub-brick from a stats dset to a NIfTI file.
extract_brick() {  # <stats_dset> <label> <out.nii.gz>
    local dset="$1" lab="$2" out="$3"
    rm -f "${out}"
    ${AFNI} "3dbucket -prefix ${out} ${dset}'[${lab}]'"
}

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
                # Fixed-effects average across sessions.
                #   betaAvg  = mean(b_i)
                #   se2_i    = (b_i/t_i)^2      (t_i==0 -> 1e6, drives t->0)
                #   se2mean  = mean(se2_i)
                #   tAvg     = betaAvg / sqrt(se2mean/n)
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
                ${AFNI} "3dcalc -a ${betaAvg} -b ${se2mean} \
                    -expr 'step(b)*a*sqrt(${n}/(b+1-step(b)))' \
                    -prefix ${outT}"
                rm -f "${tmp}"_*.nii.gz
            fi

            # Label the space so clustering / whereami work (grid is truly MNI).
            ${AFNI} "3drefit -space MNI152_2009 ${outC}; 3drefit -space MNI152_2009 ${outT}"
        done
    done
done

echo "Extracted inputs:"
ls -1 "${INPUTS}"/*_Coef.nii.gz 2>/dev/null | sed 's/^/  /'

# ==============================================================================
# STEP 2 — Group brain mask
# Primary: intersection of fMRIPrep MNI brain masks (one per subject).
# Fallback: coverage of the extracted slope Coef maps.
# ==============================================================================
echo "=============================================================="
echo "STEP 2: building group mask"
echo "=============================================================="

# Reference grid = first extracted slope Coef map.
REF=$(ls -1 "${INPUTS}"/*_task-*_slope_Coef.nii.gz 2>/dev/null | head -n1)
if [ -z "${REF}" ]; then echo "ERROR: no extracted Coef maps for masking"; exit 1; fi

submasks=()
for SUB in "${SUBJECTS[@]}"; do
    bm=$(ls -1 "${FMRIPREP_OUT}/${SUB}"/ses-*/func/${SUB}_ses-*_task-*_space-MNI152NLin2009cAsym*desc-brain_mask.nii.gz 2>/dev/null | head -n1)
    [ -z "${bm}" ] && continue
    rm_out=${MASKDIR}/${SUB}_brainmask_rs.nii.gz
    rm -f "${rm_out}"
    ${AFNI} "3dresample -master ${REF} -rmode NN -prefix ${rm_out} -input ${bm}"
    submasks+=("${rm_out}")
done

rm -f "${GROUPMASK}"
if [ "${#submasks[@]}" -ge 1 ]; then
    echo "  intersecting ${#submasks[@]} fMRIPrep brain masks"
    ${AFNI} "3dmask_tool -input ${submasks[*]} -frac 1.0 -prefix ${GROUPMASK}"
else
    echo "  WARNING: no fMRIPrep brain masks found -> falling back to Coef coverage"
    ${AFNI} "3dmask_tool -input ${INPUTS}/*_slope_Coef.nii.gz -frac 1.0 -prefix ${GROUPMASK}"
fi
${AFNI} "3drefit -space MNI152_2009 ${GROUPMASK}"

echo "  grid consistency check (want 1s):"
${AFNI} "3dinfo -same_grid ${GROUPMASK} ${INPUTS}/*_slope_Coef.nii.gz" || true

# Helper: subjects that have an extracted map for a given task/effect.
subs_with() {  # <task> <effect>
    local task="$1" eff="$2" s
    for s in "${SUBJECTS[@]}"; do
        [ -f "${INPUTS}/${s}_task-${task}_${eff}_Coef.nii.gz" ] && echo "${s}"
    done
}

# ==============================================================================
# STEP 3 — 3dMEMA one-sample per task (slope primary, mean sanity)
# ==============================================================================
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

# ==============================================================================
# STEP 4 — 3dMEMA between-task slope contrasts
# Per subject (having both tasks): bDelta = bX - bY ;
#   seDelta = sqrt(seX^2 + seY^2) ; tDelta = bDelta/seDelta ; then 1-sample MEMA.
# ==============================================================================
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
        # seX^2=(cX/tX)^2, seY^2=(cY/tY)^2 (guard t==0); tDelta=(cX-cY)/sqrt(seX2+seY2)
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

# ==============================================================================
# STEP 5 — 3dttest++ (one-sample per task; paired between-task) with -Clustsim
# ==============================================================================
echo "=============================================================="
echo "STEP 5: 3dttest++ (with -Clustsim)"
echo "=============================================================="

# One-sample per task
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

# Paired between-task (slope), subjects with BOTH tasks, matched order
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

# ==============================================================================
# STEP 6 — Residual-based cluster correction (3dFWHMx -> 3dClustSim)
# 3dMEMA has no built-in ClustSim, and 3dttest++ -Clustsim needs >=14 samples,
# so this table is the general-purpose corrector: it depends only on the mask,
# the noise smoothness, and CLUST_PTHR -- NOT on subject count -- so it applies
# to BOTH the MEMA maps and the (small-n) 3dttest++ maps. Apply via 3dClusterize.
# ==============================================================================
echo "=============================================================="
echo "STEP 6: 3dFWHMx -> 3dClustSim for MEMA cluster thresholds"
echo "=============================================================="

acf_file=${CLUSTDIR}/acf_params.txt
: > "${acf_file}"
for ERRTS in "${AFNI_FIRSTLVL}"/sub-*/ses-*/task-*/*_errts_REML+orig.HEAD; do
    [ -e "${ERRTS}" ] || continue
    dset=${ERRTS%.HEAD}
    line=$(${AFNI} "3dFWHMx -acf ${CLUSTDIR}/acf_est.1D -mask ${GROUPMASK} -input ${dset} 2>/dev/null" | tail -n1)
    echo "${line}" >> "${acf_file}"
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

# ==============================================================================
# STEP 7 — Persist small final maps from scratch back to archive
# (autodelete is periodically purged; the group results are small.)
# ==============================================================================
echo "=============================================================="
echo "STEP 7: copying final maps to ${GROUP_FINAL}"
echo "=============================================================="
mkdir -p "${GROUP_FINAL}"
for d in mask mema ttest clustsim; do
    [ -d "${GROUP_OUT}/${d}" ] || continue
    mkdir -p "${GROUP_FINAL}/${d}"
    cp -a "${GROUP_OUT}/${d}"/. "${GROUP_FINAL}/${d}/" 2>/dev/null || true
done

echo "=============================================================="
echo "Group analysis complete."
echo "  n subjects discovered: ${#SUBJECTS[@]}  (${SUBJECTS[*]})"
echo "  scratch (working) : ${GROUP_OUT}"
echo "  archive (kept)    : ${GROUP_FINAL}  <- mema/ ttest/ mask/ clustsim/"
echo "  NOTE: with only a few subjects this is UNDERPOWERED - expect no"
echo "        surviving clusters. Re-run as sub-02/sub-04 and more sessions"
echo "        finish; new first-level output is picked up automatically."
echo "=============================================================="
