#!/bin/bash

# Define home directory
home_dir=~/nobackup/archive/Nielsen_active/Spirituality/Project/rawdata
bids_dir=~/nobackup/archive/Nielsen_active/Spirituality/Project/BIDS

# Directory this script lives in (so we can call assign_fieldmaps.py alongside it)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Options ------------------------------------------------------------------
# By default, any session whose BIDS output already contains everything its
# source DICOMs can produce is SKIPPED, so a re-run only touches new or
# partially-converted sessions and never re-writes a finished one. Pass --force
# to re-convert everything regardless. Re-conversion is always safe: every
# dcm2niix call below uses -w 1 (overwrite in place), so a redo overwrites the
# matching files instead of spawning duplicate "_a" copies (dcm2niix's default
# -w 2 behavior), which would otherwise confuse assign_fieldmaps.py / fMRIPrep.
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force|-f) FORCE=1 ;;
        *) echo "WARNING: unknown argument '$arg' (only --force is recognized)" ;;
    esac
done

# Return 0 (true) if a session's BIDS output already contains everything the
# source DICOMs can produce: T1w/T2w, 4 echoes per TASK run (ME4), and the AP+PA
# fieldmaps. Deliberately uses the SAME series-name patterns as the conversion
# steps below, so "expected" always matches "what this script would actually
# write" -- anat-only sessions (func purged) count as complete and are left
# untouched, while a session missing convertible func (e.g. anat present but
# BOLD not yet converted) counts as incomplete and gets (re)built.
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

# Collect participant folders sorted by date (YYYY.MM.DD prefix sorts chronologically)
declare -A byv_to_sub    # maps BYV ID -> zero-padded sub number
declare -A byv_ses_count # maps BYV ID -> sessions seen so far
sub_counter=0

mapfile -t sorted_folders < <(find "${home_dir}" -maxdepth 1 -type d -name "*Nielsen*" | sort)

# Loop through all subject folders in chronological order
for ses in "${sorted_folders[@]}"; do
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

    # Normalize the inner series folder to "Nielsen_Active". Handles names with
    # spaces (e.g. "Nielsen 2 Active"), is idempotent (skips if already done), and
    # never moves a folder *into* an existing Nielsen_Active.
    if [ ! -d "${ses}/Nielsen_Active" ]; then
        inner_src=$(find "${ses}" -maxdepth 1 -mindepth 1 -type d -iname "*Nielsen*" | head -n1)
        [ -n "${inner_src}" ] && mv "${inner_src}" "${ses}/Nielsen_Active"
    fi

    dicom_dir="${ses}/Nielsen_Active"
    out_dir=${bids_dir}/sub-${sub}/ses-${number_ses}

    # Create output directories
    anat_dir=${out_dir}/anat
    func_dir=${out_dir}/func
    fmap_dir=${out_dir}/fmap
    mkdir -p ${anat_dir} ${func_dir} ${fmap_dir}

    # Set output string
    name_string=sub-${sub}_ses-${number_ses}

    # ----------------------------------------------------------
    # SKIP-IF-DONE: leave already-complete sessions untouched (unless --force),
    # so a full re-run is idempotent and only builds new/incomplete sessions.
    # The participant map (byv_to_sub / byv_ses_count) is already updated above,
    # so participants.tsv stays correct even for skipped sessions.
    # ----------------------------------------------------------
    if [ "${FORCE}" != "1" ] && session_is_complete "${dicom_dir}" "${anat_dir}" "${func_dir}" "${fmap_dir}"; then
        echo "  Skipping sub-${sub} ses-${number_ses}: already fully converted (use --force to redo)"
        continue
    fi

    # ----------------------------------------------------------
    # TASK NAMING: Look up this participant's task order from
    # run_order_plain.tsv using their alphanumeric ID (the part
    # of the folder name that follows "BYV"). The order code is
    # a 3-letter string where F=FHS, S=scripture, A=architecture,
    # and position (1st/2nd/3rd letter) = slot (1st/2nd/3rd run).
    # Column 2 = session 1 order, column 3 = session 2 order.
    # ----------------------------------------------------------
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

    # ----------------------------------------------------------
    # FUNCTIONAL BOLD: Convert each TASK folder (exclude SBRef)
    # in acquisition order, assigning task names by slot position.
    # ----------------------------------------------------------
    # Leading/trailing * tolerate the "MR " series prefix and the numeric suffix
    # (e.g. "MR Minn_2.5mm_ME4_S4p2_TASK1"); exclude the SBRef companions.
    mapfile -t ep2d_dirs < <(find "${dicom_dir}" -maxdepth 1 -type d \
        -name "*Minn_2.5mm_ME4_S4p2_TASK*" ! -name "*SBRef*" | sort)

    for i in "${!ep2d_dirs[@]}"; do
        task="${task_names[$i]}"
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${func_dir}" \
            -f "${name_string}_task-${task}_echo-%e_bold" "${ep2d_dirs[$i]}/"
    done

    # ----------------------------------------------------------
    # SINGLE-BAND REFERENCE: One SBRef per BOLD run
    # ----------------------------------------------------------
    mapfile -t sbref_dirs < <(find "${dicom_dir}" -maxdepth 1 -type d \
        -name "*Minn_2.5mm_ME4_S4p2_TASK*SBRef*" | sort)

    for i in "${!sbref_dirs[@]}"; do
        task="${task_names[$i]}"
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${func_dir}" \
            -f "${name_string}_task-${task}_echo-%e_sbref" "${sbref_dirs[$i]}/"
    done

    # ----------------------------------------------------------
    # ANATOMICAL: T1w and T2w (standard variants only, not ND)
    # ----------------------------------------------------------
    if [ -d "${dicom_dir}/MR T1w_MPR" ]; then
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${anat_dir}" \
            -f "${name_string}_T1w" "${dicom_dir}/MR T1w_MPR/"
    fi

    if [ -d "${dicom_dir}/MR T2w_SPC" ]; then
        dcm2niix -b y -ba y -z y -i y -w 1 -o "${anat_dir}" \
            -f "${name_string}_T2w" "${dicom_dir}/MR T2w_SPC/"
    fi

    # ----------------------------------------------------------
    # FIELDMAPS: AP and PA spin-echo EPI pairs
    # ----------------------------------------------------------
    convert_dicom() {
        local sequence_type=$1
        local output_dir=$2
        local filename_prefix=$3
        local suffix=$4

        count=0
        # Leading * tolerates the "MR " series prefix (e.g. "MR CMRR_SEFM_..._AP1")
        for i in "${dicom_dir}/"*"${sequence_type}"*; do
            if [ -d "$i" ]; then
                count=$((count+1))
                dcm2niix -b y -ba y -z y -i y -w 1 -o "${output_dir}" \
                    -f "${filename_prefix}_run-${count}_${suffix}" "${i}/"
            fi
        done
    }

    convert_dicom "CMRR_SEFM_2.5mm_S4p2_AP" "${fmap_dir}" "${name_string}_dir-AP" "epi"
    convert_dicom "CMRR_SEFM_2.5mm_S4p2_PA" "${fmap_dir}" "${name_string}_dir-PA" "epi"

    # ----------------------------------------------------------
    # FIELDMAP PAIRING: refined per-run B0FieldIdentifier/B0FieldSource +
    # IntendedFor via assign_fieldmaps.py. Replaces the old coarse
    # all-fieldmaps->all-runs IntendedFor, which fMRIPrep-era re-runs
    # would clobber and which gave every run the same blanket pairing.
    # (Handles anat-only sessions gracefully: "nothing to pair".)
    # ----------------------------------------------------------
    python3 "${script_dir}/assign_fieldmaps.py" "${func_dir}" "${fmap_dir}" "ses-${number_ses}"

    echo "Finished processing subject ${sub}"
done

# ----------------------------------------------------------
# PARTICIPANTS MAP: write the BYV-id <-> sub-NN mapping that was built in memory
# above so downstream scripts can resolve it. The first-level analysis
# (first_level_afni.sh) keys timing files by the BYV alphanumeric id but loops
# over BIDS sub-NN dirs, so it reads this file. Columns:
#   participant_id  (sub-NN)
#   byv_id          (alphanumeric, "BYV" prefix stripped -- matches run_order_plain.tsv)
#   n_sessions      (number of sessions converted for this participant)
# Written to rawdata so it lives alongside the source data.
# ----------------------------------------------------------
participants_tsv="${home_dir}/participants.tsv"
printf "participant_id\tbyv_id\tn_sessions\n" > "${participants_tsv}"
for byv_id in "${!byv_to_sub[@]}"; do
    printf "sub-%s\t%s\t%s\n" \
        "${byv_to_sub[$byv_id]}" "${byv_id#BYV}" "${byv_ses_count[$byv_id]}"
done | sort -t$'\t' -k1,1V >> "${participants_tsv}"

echo "Wrote participant map: ${participants_tsv}"
