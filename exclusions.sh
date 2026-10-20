# =============================================================================
# exclusions.sh — the ONE place pipeline scripts ask "is this run excluded?"
#
# Sourced by config.sh, so every script that sources config.sh has is_excluded.
# The list itself is exclusions.tsv next to this file (tab-separated):
#
#   byv_id  session  task  reason  date_added
#   1VK4L2  ses-1    *     Asleep for most of the session ...  2026-10-20
#
# session / task may be "*" (= all). Keyed by the BYV id because that is what
# the lab records use; BIDS sub-NN labels are mapped through participants.tsv.
#
# WHY A LIST AND NOT JUST "don't make timing files": every downstream script
# (second level, reproducibility, TRR export) discovers subjects from whatever
# is ON DISK, so first-level outputs from an earlier run of an excluded session
# would otherwise be picked up again. Each discovery loop calls is_excluded.
#
# To exclude a run: add a row here with the reason, commit, and re-run from the
# first level (outputs already on disk are ignored, not deleted).
# =============================================================================

EXCLUSIONS_TSV=${EXCLUSIONS_TSV:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/exclusions.tsv"}
export EXCLUSIONS_TSV

# is_excluded <sub-NN | BYV id> <ses-N> <task>
#   exit 0 + prints the reason  -> excluded
#   exit 1                      -> not excluded (also when the list is absent)
# Use as:  if why=$(is_excluded "${SUB}" "${SES}" "${TASK}"); then ...; fi
is_excluded() {
    local who="$1" ses="$2" task="$3" byv="$1"
    [ -f "${EXCLUSIONS_TSV}" ] || return 1
    case "${who}" in
        sub-*)
            byv=$(awk -F'\t' -v w="${who}" '
                NR==1 { for (i=1;i<=NF;i++) { if ($i=="participant_id") p=i; if ($i=="byv_id") b=i }; next }
                $p==w { print $b; exit }' "${PARTICIPANTS_TSV}" 2>/dev/null)
            ;;
    esac
    [ -n "${byv}" ] || return 1
    awk -F'\t' -v b="${byv}" -v s="${ses}" -v t="${task}" '
        NR==1 || /^#/ { next }
        $1==b && ($2==s || $2=="*") && ($3==t || $3=="*") { print $4; found=1; exit }
        END { exit(found ? 0 : 1) }' "${EXCLUSIONS_TSV}"
}
