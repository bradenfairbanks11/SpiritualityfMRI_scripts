#!/bin/bash
# ==============================================================================
# dcm2niix.sh — CONVERT RAW DICOMS TO A BIDS DATASET
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   The scanner saves images as DICOM folders, one folder per series. This
#   script converts every session's DICOMs to compressed NIfTI files (plus a
#   JSON "sidecar" of scan settings for each) with dcm2niix, and names and
#   files them in BIDS layout (sub-NN/ses-N/anat, func, fmap).
#   Along the way it:
#     - numbers participants sub-01, sub-02, ... in scan-date order, and
#       numbers each participant's sessions ses-1, ses-2, ... in date order
#     - names each task run (scripture / FHS / architecture) from the
#       participant's counterbalanced run order in run_order_plain.tsv
#     - pairs each spin-echo fieldmap with its BOLD run (assign_fieldmaps.py)
#     - writes participants.tsv, the map from BIDS sub-NN to the lab's BYV id
#       that every later script uses
#   Sessions that are already fully converted are skipped (unless --force).
#
# WHERE IT SITS IN THE PIPELINE
#   Before: nothing (this is the first step; raw DICOMs come from the scanner)
#   After : fmriprep_tedana.sh  (preprocessing of the BIDS dataset)
#
# STEPS, PER SESSION (one DICOM folder = one participant x session)
#   1. Work out the participant number (sub-NN) and session number (ses-N)
#   2. Rename the inner series folder to "Nielsen_Active" (if needed)
#   3. Make the BIDS output folders; skip the session if already converted
#   4. Look up the session's task order in run_order_plain.tsv
#   5. Convert the multi-echo BOLD runs (4 echoes each), named by task
#   6. Convert the single-band reference (SBRef) images, named by task
#   7. Convert the T1w and T2w anatomical scans
#   8. Convert the AP and PA spin-echo fieldmaps
#   9. Pair fieldmaps with BOLD runs (assign_fieldmaps.py)
#   Then, after all sessions:
#  10. Write participants.tsv (sub-NN <-> BYV id, number of sessions)
#
# INPUTS
#   DICOM sessions : <RAWDATA_DIR>/*Nielsen*/   one folder per session. The
#                    folder name starts with the scan date (YYYY.MM.DD), and its
#                    3rd "_"-separated field is the participant's BYV id
#                    ("BYV" + alphanumeric id). Inside it, one folder (any name
#                    containing "Nielsen") holds the series folders:
#                      MR Minn_2.5mm_ME4_S4p2_TASK<n>          multi-echo BOLD
#                      MR Minn_2.5mm_ME4_S4p2_TASK<n>_SBRef    its SBRef
#                      MR T1w_MPR, MR T2w_SPC                  anatomicals
#                      MR CMRR_SEFM_2.5mm_S4p2_AP<n> / _PA<n>  fieldmaps
#   run order      : <RAWDATA_DIR>/run_order_plain.tsv
#                    column 1 = alphanumeric id, column 2 = session-1 order,
#                    column 3 = session-2 order; first row is a header.
#                    An order is 3 letters (F = FHS, S = scripture,
#                    A = architecture); letter k = the task of the k-th run.
#
# OUTPUTS
#   In <BIDS_DIR>/sub-NN/ses-N/   (each .nii.gz comes with a .json sidecar)
#     func/sub-NN_ses-N_task-T_echo-E_bold.nii.gz    E = 1..4
#     func/sub-NN_ses-N_task-T_echo-E_sbref.nii.gz
#     anat/sub-NN_ses-N_T1w.nii.gz, anat/sub-NN_ses-N_T2w.nii.gz
#     fmap/sub-NN_ses-N_dir-AP_run-K_epi.nii.gz, ..._dir-PA_run-K_epi.nii.gz
#     The fmap and func JSONs are then edited by assign_fieldmaps.py
#     (B0FieldIdentifier, B0FieldSource, IntendedFor).
#   <RAWDATA_DIR>/participants.tsv   columns participant_id, byv_id, n_sessions
#   In <RAWDATA_DIR>, each session's inner series folder is renamed to
#   "Nielsen_Active" (STEP 2).
#   NOTE: <RAWDATA_DIR> is on the archive storage area, so this script writes to
#   archive (participants.tsv and the folder renames).
#
# HOW TO RUN (from the scripts folder; there are no #SBATCH lines)
#   bash dcm2niix.sh            # convert new or partly converted sessions only
#   bash dcm2niix.sh --force    # re-convert every session
#   Needs the dcm2niix and python3 commands on the PATH (this script does not
#   load any module). Run it with bash, not sbatch: it finds
#   assign_fieldmaps.py next to itself (see script_dir in SECTION 1), which
#   does not work under sbatch.
#
# METHODS SUMMARY (the exact values used below)
#   Converter      : dcm2niix, called identically for every series:
#                      -b y   write a BIDS JSON sidecar
#                      -ba y  anonymize the sidecar
#                      -z y   gzip-compress the NIfTI (.nii.gz)
#                      -i y   ignore derived, localizer and 2D images
#                      -w 1   overwrite a file of the same name (so a re-run
#                             never makes duplicate "_a" copies)
#   BOLD           : multi-echo (4 echoes), one file per echo (echo-%e)
#   Task labels    : from each participant's counterbalanced order
#                    (run_order_plain.tsv), assigned to the BOLD runs in
#                    series-name order
#   Fieldmaps      : spin-echo EPI pairs, phase-encoding AP and PA; each pair is
#                    matched to the BOLD run acquired right after it
#                    (assign_fieldmaps.py)
# ==============================================================================


# ==============================================================================
# SECTION 1. SETTINGS
# ------------------------------------------------------------------------------
# Every path this script uses is defined here.
# ==============================================================================

# --- 1a. Where the data live --------------------------------------------------
# Three storage areas on the BYU supercomputer:
#   /home/bradenf4                     backed up, never purged: code, timing, logs
#   /nobackup/autodelete/usr/bradenf4  scratch, NOT backed up, purged after 12
#                                      weeks unused: the BIDS dataset
#   /nobackup/archive/usr/bradenf4     raw data and archived results
RAWDATA_DIR=/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/rawdata   # raw DICOMs (archive)
BIDS_DIR=/nobackup/autodelete/usr/bradenf4/Spirituality/BIDS                             # BIDS output (scratch)

# The names the rest of this script uses for those two folders.
home_dir=${RAWDATA_DIR}
bids_dir=${BIDS_DIR}

# --- 1b. Software -------------------------------------------------------------
# dcm2niix and python3 are whatever is on the PATH (no module is loaded here).
#
# Folder this script lives in, so assign_fieldmaps.py (STEP 9) can be found
# next to it. This works when the script is run with bash. Under sbatch,
# BASH_SOURCE points at Slurm's spool copy of the script instead, so the
# assign_fieldmaps.py call would fail. (The scripts live in
# /home/bradenf4/spirituality_fmri/scripts on the cluster.)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ignore any Python packages in ~/.local. A stray numpy there once broke tedana;
# this script's python3 call (STEP 9) only needs the standard library.
export PYTHONNOUSERSITE=1

# --- 1c. Conversion parameters ------------------------------------------------
# The dcm2niix options and the series-name patterns are written out in each
# STEP below (see METHODS SUMMARY for what every dcm2niix option means).


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# session_is_complete DICOM_DIR ANAT_DIR FUNC_DIR FMAP_DIR
#   Succeeds (returns 0) if the session's BIDS folders already hold everything
#   its DICOMs can produce:
#     - 4 *_bold.nii.gz files (one per echo) per TASK series   (func)
#     - one *_epi.nii.gz per AP or PA fieldmap series          (fmap)
#     - a *_T1w.nii.gz if there is a "MR T1w_MPR" series       (anat)
#     - a *_T2w.nii.gz if there is a "MR T2w_SPC" series       (anat)
#   ("Has at least" that many files counts as complete. SBRef files are not
#   counted.) It uses the SAME series-name patterns as the conversion steps, so
#   "expected" always matches "what this script would actually write". An
#   anatomy-only session (functional series purged) therefore counts as
#   complete and is left alone, while a session whose BOLD has not been
#   converted yet counts as incomplete and is (re)built.
session_is_complete() {
    local dicom_dir="$1" anat_dir="$2" func_dir="$3" fmap_dir="$4"
    local src_task src_ap src_pa src_t1 src_t2 exp_bold exp_epi
    src_task=$(find "${dicom_dir}" -maxdepth 1 -type d -name "*Minn_2.5mm_ME4_S4p2_TASK*" ! -name "*SBRef*" 2>/dev/null | wc -l)
    src_ap=$(find "${dicom_dir}" -maxdepth 1 -type d -name "*CMRR_SEFM_2.5mm_S4p2_AP*" 2>/dev/null | wc -l)
    src_pa=$(find "${dicom_dir}" -maxdepth 1 -type d -name "*CMRR_SEFM_2.5mm_S4p2_PA*" 2>/dev/null | wc -l)
    [ -d "${dicom_dir}/MR T1w_MPR" ] && src_t1=1 || src_t1=0
    [ -d "${dicom_dir}/MR T2w_SPC" ] && src_t2=1 || src_t2=0
    exp_bold=$(( src_task * 4 ))   # ME4 protocol = 4 echoes per task run
    exp_epi=$(( src_ap + src_pa ))

    local have_bold have_epi have_t1 have_t2
    have_bold=$(find "${func_dir}" -maxdepth 1 -name "*_bold.nii.gz" 2>/dev/null | wc -l)
    have_epi=$(find "${fmap_dir}"  -maxdepth 1 -name "*_epi.nii.gz"  2>/dev/null | wc -l)
    have_t1=$(find "${anat_dir}"   -maxdepth 1 -name "*_T1w.nii.gz"  2>/dev/null | wc -l)
    have_t2=$(find "${anat_dir}"   -maxdepth 1 -name "*_T2w.nii.gz"  2>/dev/null | wc -l)

    [ "${have_bold}" -ge "${exp_bold}" ] \
        && [ "${have_epi}" -ge "${exp_epi}" ] \
        && [ "${have_t1}"  -ge "${src_t1}" ] \
        && [ "${have_t2}"  -ge "${src_t2}" ]
}

# convert_dicom SERIES_PATTERN OUTPUT_DIR FILENAME_PREFIX SUFFIX   (STEP 8)
#   Converts every series folder in the CURRENT session's ${dicom_dir} whose
#   name contains SERIES_PATTERN, numbering them run-1, run-2, ... in the
#   (alphabetical) order bash lists them. Output names:
#     <FILENAME_PREFIX>_run-<K>_<SUFFIX>.nii.gz
#   The leading * in the match tolerates the "MR " series prefix
#   (e.g. "MR CMRR_SEFM_..._AP1").
convert_dicom() {
    local sequence_type=$1
    local output_dir=$2
    local filename_prefix=$3
    local suffix=$4

    count=0
    for i in "${dicom_dir}/"*"${sequence_type}"*; do
        if [ -d "$i" ]; then
            count=$((count+1))
            dcm2niix -b y -ba y -z y -i y -w 1 -o "${output_dir}" \
                -f "${filename_prefix}_run-${count}_${suffix}" "${i}/"
        fi
    done
}


# ==============================================================================
# SECTION 3. SETUP: command-line options and the list of sessions
# ==============================================================================

# --- Command-line option ------------------------------------------------------
# By default, any session whose BIDS output already contains everything its
# source DICOMs can produce is SKIPPED (see session_is_complete), so a re-run
# only touches new or partly converted sessions and never re-writes a finished
# one. --force (or -f) re-converts everything. Re-conversion is always safe:
# every dcm2niix call uses -w 1 (overwrite in place), so a redo overwrites the
# matching files instead of making duplicate "_a" copies (dcm2niix's default
# -w 2), which would otherwise confuse assign_fieldmaps.py / fMRIPrep.
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force|-f) FORCE=1 ;;
        *) echo "WARNING: unknown argument '$arg' (only --force is recognized)" ;;
    esac
done

# --- Participant and session counters -----------------------------------------
# Built up as the sessions are read in date order (STEP 1):
declare -A byv_to_sub    # maps BYV ID -> zero-padded sub number
declare -A byv_ses_count # maps BYV ID -> sessions seen so far
sub_counter=0

# --- Every session folder, oldest first ---------------------------------------
# Folder names start with the scan date (YYYY.MM.DD), so an alphabetical sort is
# also a chronological sort.
mapfile -t sorted_folders < <(find "${home_dir}" -maxdepth 1 -type d -name "*Nielsen*" | sort)


# ==============================================================================
# SECTION 4. MAIN LOOP: one session folder at a time, in date order
# ==============================================================================

for ses in "${sorted_folders[@]}"; do

    # --------------------------------------------------------------------------
    # STEP 1. Work out sub-NN and ses-N.
    #
    # The BYV id is the 3rd "_"-separated field of the folder name. The first
    # time a BYV id is seen it gets the next participant number (01, 02, ...);
    # each time it is seen, its session count goes up by 1. Because folders are
    # read oldest first, participants are numbered in order of their first scan
    # and sessions in date order.
    # --------------------------------------------------------------------------
    folder_name=$(basename "$ses")
    byv_id=$(echo "$folder_name" | cut -d'_' -f3)

    # Assign a new sub number the first time this BYV ID is seen
    if [[ -z "${byv_to_sub[$byv_id]}" ]]; then
        sub_counter=$((sub_counter + 1))
        byv_to_sub[$byv_id]=$(printf "%02d" $sub_counter)
        byv_ses_count[$byv_id]=0
    fi

    # Increment session count and set working variables
    byv_ses_count[$byv_id]=$((${byv_ses_count[$byv_id]} + 1))
    sub=${byv_to_sub[$byv_id]}
    number_ses=${byv_ses_count[$byv_id]}

    echo "Processing subject ${sub} in session ${ses}"

    # --------------------------------------------------------------------------
    # STEP 2. Rename the inner series folder to "Nielsen_Active".
    #
    # The folder holding the series has different names in different sessions
    # (e.g. "Nielsen 2 Active"). The first sub-folder (as listed by find) whose name contains
    # "nielsen" (any case) is renamed to Nielsen_Active. Nothing happens if
    # Nielsen_Active already exists, so the step is safe to repeat and never
    # moves a folder *into* an existing Nielsen_Active.
    # (This renames a folder inside RAWDATA_DIR, on archive.)
    # --------------------------------------------------------------------------
    if [ ! -d "${ses}/Nielsen_Active" ]; then
        inner_src=$(find "${ses}" -maxdepth 1 -mindepth 1 -type d -iname "*Nielsen*" | head -n1)
        [ -n "${inner_src}" ] && mv "${inner_src}" "${ses}/Nielsen_Active"
    fi

    dicom_dir="${ses}/Nielsen_Active"
    out_dir=${bids_dir}/sub-${sub}/ses-${number_ses}

    # --------------------------------------------------------------------------
    # STEP 3. Make the BIDS output folders, and skip the session if it is
    #         already fully converted (unless --force).
    #
    # The participant map (byv_to_sub / byv_ses_count) was already updated in
    # STEP 1, so participants.tsv stays correct even for skipped sessions.
    # (The output folders are made before the check, so they exist even for a
    # skipped session.)
    # --------------------------------------------------------------------------
    anat_dir=${out_dir}/anat
    func_dir=${out_dir}/func
    fmap_dir=${out_dir}/fmap
    mkdir -p ${anat_dir} ${func_dir} ${fmap_dir}

    # Start of every output file name: sub-NN_ses-N
    name_string=sub-${sub}_ses-${number_ses}

    if [ "${FORCE}" != "1" ] && session_is_complete "${dicom_dir}" "${anat_dir}" "${func_dir}" "${fmap_dir}"; then
        echo "  Skipping sub-${sub} ses-${number_ses}: already fully converted (use --force to redo)"
        continue
    fi

    # --------------------------------------------------------------------------
    # STEP 4. Look up this session's task order.
    #
    # The alphanumeric id is the BYV id without its "BYV" prefix. In
    # run_order_plain.tsv (header row skipped), the row whose column 1 is that
    # id gives the order: column 2 for session 1, column 3 for session 2.
    # The order is 3 letters, F = FHS, S = scripture, A = architecture, and
    # letter k is the task of the k-th task run (slot k).
    # If no order is found, the session is skipped with an ERROR message.
    # --------------------------------------------------------------------------
    alnum_id="${byv_id#BYV}"
    order_col=$((number_ses + 1))  # ses 1 -> col 2, ses 2 -> col 3
    order_code=$(awk -v id="$alnum_id" -v col="$order_col" \
        'NR>1 && $1==id {print $col}' "${home_dir}/run_order_plain.tsv")

    if [[ -z "$order_code" ]]; then
        echo "ERROR: ${alnum_id} not found in run_order_plain.tsv or no order for session ${number_ses}. Skipping ${ses}."
        continue
    fi

    declare -A task_map
    task_map[F]="FHS"
    task_map[S]="scripture"
    task_map[A]="architecture"

    slot1_task="${task_map[${order_code:0:1}]}"
    slot2_task="${task_map[${order_code:1:1}]}"
    slot3_task="${task_map[${order_code:2:1}]}"
    task_names=("$slot1_task" "$slot2_task" "$slot3_task")

    # --------------------------------------------------------------------------
    # STEP 5. Convert the multi-echo BOLD runs.
    #
    # Every "...Minn_2.5mm_ME4_S4p2_TASK..." series folder except the SBRef
    # ones, sorted by name (TASK1, TASK2, TASK3 = acquisition order). The k-th
    # folder gets the k-th task from STEP 4. dcm2niix writes one file per echo:
    # %e in the -f name is replaced by the echo number.
    # The leading/trailing * tolerate the "MR " series prefix and the numeric
    # suffix (e.g. "MR Minn_2.5mm_ME4_S4p2_TASK1").
    # --------------------------------------------------------------------------
    mapfile -t ep2d_dirs < <(find "${dicom_dir}" -maxdepth 1 -type d \
        -name "*Minn_2.5mm_ME4_S4p2_TASK*" ! -name "*SBRef*" | sort)

    for i in "${!ep2d_dirs[@]}"; do
        task="${task_names[$i]}"
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${func_dir}" \
            -f "${name_string}_task-${task}_echo-%e_bold" "${ep2d_dirs[$i]}/"
    done

    # --------------------------------------------------------------------------
    # STEP 6. Convert the single-band reference (SBRef) images.
    #
    # One SBRef series per BOLD run ("...TASK...SBRef..."), sorted by name and
    # given task names in the same way as STEP 5.
    # --------------------------------------------------------------------------
    mapfile -t sbref_dirs < <(find "${dicom_dir}" -maxdepth 1 -type d \
        -name "*Minn_2.5mm_ME4_S4p2_TASK*SBRef*" | sort)

    for i in "${!sbref_dirs[@]}"; do
        task="${task_names[$i]}"
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${func_dir}" \
            -f "${name_string}_task-${task}_echo-%e_sbref" "${sbref_dirs[$i]}/"
    done

    # --------------------------------------------------------------------------
    # STEP 7. Convert the anatomical scans: T1w and T2w.
    #
    # Only the folders named exactly "MR T1w_MPR" and "MR T2w_SPC" are used
    # (the standard images, not the ND variants). Each is converted only if it
    # exists.
    # --------------------------------------------------------------------------
    if [ -d "${dicom_dir}/MR T1w_MPR" ]; then
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${anat_dir}" \
            -f "${name_string}_T1w" "${dicom_dir}/MR T1w_MPR/"
    fi

    if [ -d "${dicom_dir}/MR T2w_SPC" ]; then
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${anat_dir}" \
            -f "${name_string}_T2w" "${dicom_dir}/MR T2w_SPC/"
    fi

    # --------------------------------------------------------------------------
    # STEP 8. Convert the spin-echo fieldmaps: AP and PA.
    #
    # Each phase-encoding direction is converted separately with convert_dicom
    # (SECTION 2), so the AP series become dir-AP_run-1, run-2, ... and the PA
    # series become dir-PA_run-1, run-2, ...
    # --------------------------------------------------------------------------
    convert_dicom "CMRR_SEFM_2.5mm_S4p2_AP" "${fmap_dir}" "${name_string}_dir-AP" "epi"
    convert_dicom "CMRR_SEFM_2.5mm_S4p2_PA" "${fmap_dir}" "${name_string}_dir-PA" "epi"

    # --------------------------------------------------------------------------
    # STEP 9. Pair each fieldmap with its BOLD run (assign_fieldmaps.py).
    #
    # Each AP/PA pair is acquired right before its task run. assign_fieldmaps.py
    # matches every fieldmap to the run acquired right after it (by the
    # AcquisitionTime in the JSON sidecars, or by SeriesNumber if any file in
    # the session lacks one) and writes, per run:
    # B0FieldIdentifier on the fieldmap JSONs, B0FieldSource on the func JSONs
    # (same string), and an IntendedFor list on the fieldmaps. fMRIPrep uses
    # these to apply the right distortion correction to each run. This replaced
    # an older version of this script that pointed every fieldmap at every run.
    # Anatomy-only sessions are handled: the script reports "nothing to pair".
    # --------------------------------------------------------------------------
    python3 "${script_dir}/assign_fieldmaps.py" "${func_dir}" "${fmap_dir}" "ses-${number_ses}"

    echo "Finished processing subject ${sub}"
done


# ==============================================================================
# SECTION 5. WRITE participants.tsv
# ------------------------------------------------------------------------------
# STEP 10. Write the BYV-id <-> sub-NN map built in STEP 1, so later scripts can
# use it: the first-level scripts loop over BIDS sub-NN folders but the timing
# files are named by BYV id, so they read this file. Columns:
#   participant_id  (sub-NN)
#   byv_id          (alphanumeric, "BYV" prefix removed -- matches run_order_plain.tsv)
#   n_sessions      (number of session folders found for this participant,
#                    including skipped ones)
# Rows are sorted by participant_id. The file is written to RAWDATA_DIR so it
# lives alongside the source data (it is rewritten from scratch every run).
# ==============================================================================
participants_tsv="${home_dir}/participants.tsv"
printf "participant_id\tbyv_id\tn_sessions\n" > "${participants_tsv}"
for byv_id in "${!byv_to_sub[@]}"; do
    printf "sub-%s\t%s\t%s\n" \
        "${byv_to_sub[$byv_id]}" "${byv_id#BYV}" "${byv_ses_count[$byv_id]}"
done | sort -t$'\t' -k1,1V >> "${participants_tsv}"

echo "Wrote participant map: ${participants_tsv}"
