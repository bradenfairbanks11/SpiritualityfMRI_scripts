#!/bin/bash
# ==============================================================================
# check_settings.sh — CONFIRM THAT SHARED SETTINGS MATCH ACROSS SCRIPTS
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   Each pipeline script carries its own copy of the settings it needs (paths,
#   the AFNI image, TR, task durations, ...), so it can be read on its own. A
#   setting that appears in more than one script is marked  # [shared]  at the
#   end of its line. This script collects every [shared] setting from every
#   script and reports any that have different values in different scripts.
#
#   It reads files only and changes nothing. It is not part of the pipeline;
#   run it after editing any [shared] setting.
#
# WHAT IS COMPARED
#   - variable settings:   NAME=value        # [shared]
#   - function settings:   name() {          # [shared]   (the whole function
#                          ...                              body is compared)
#                          }
#
# HOW TO RUN (from the scripts folder, on the cluster or a Mac)
#   bash check_settings.sh
#   Exit status 0 = everything matches, 1 = at least one mismatch.
# ==============================================================================

set -uo pipefail

# The scripts to check: every .sh file in the same folder as this one.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# --- Collect every [shared] setting as:  name <TAB> script <TAB> value --------
# Variables: the text after "NAME=", up to the [shared] marker.
# Functions: every line from "name() {" to the closing "}" at line start.
collect() {
    for f in "${SCRIPT_DIR}"/*.sh; do
        [ "$(basename "$f")" = "check_settings.sh" ] && continue
        awk -v file="$(basename "$f")" '
            # Strip the trailing "# [shared] ..." comment and surrounding spaces.
            function clean(s) { sub(/[[:space:]]*#[[:space:]]*\[shared\].*$/, "", s); sub(/^[[:space:]]+/, "", s); return s }
            infn {
                body = body "\\n" $0
                if ($0 ~ /^}/) { printf "%s()\t%s\t%s\n", fname, file, body; infn = 0 }
                next
            }
            /#[[:space:]]*\[shared\]/ {
                line = clean($0)
                if (line ~ /^[A-Za-z_][A-Za-z_0-9]*\(\)[[:space:]]*\{/) {
                    fname = line; sub(/\(\).*/, "", fname); body = ""; infn = 1
                } else if (match(line, /^(export[[:space:]]+)?[A-Za-z_][A-Za-z_0-9]*=/)) {
                    name = substr(line, 1, RLENGTH - 1); sub(/^export[[:space:]]+/, "", name)
                    printf "%s\t%s\t%s\n", name, file, substr(line, RLENGTH + 1)
                }
            }' "$f"
    done
}

ROWS=$(collect)
if [ -z "${ROWS}" ]; then
    echo "No [shared] settings found in ${SCRIPT_DIR}/*.sh"
    exit 0
fi

# --- Compare: for each setting name, list the distinct values -----------------
NBAD=0
NAMES=$(printf '%s\n' "${ROWS}" | cut -f1 | sort -u)
for name in ${NAMES}; do
    vals=$(printf '%s\n' "${ROWS}" | awk -F'\t' -v n="${name}" '$1==n { print $3 }' | sort -u)
    nvals=$(printf '%s\n' "${vals}" | wc -l | tr -d ' ')
    nfiles=$(printf '%s\n' "${ROWS}" | awk -F'\t' -v n="${name}" '$1==n' | wc -l | tr -d ' ')
    if [ "${nvals}" -eq 1 ]; then
        printf '  ok        %-22s (%s scripts)\n' "${name}" "${nfiles}"
    else
        NBAD=$((NBAD + 1))
        printf '  MISMATCH  %-22s has %s different values:\n' "${name}" "${nvals}"
        printf '%s\n' "${ROWS}" | awk -F'\t' -v n="${name}" '$1==n { v = $3; if (length(v) > 90) v = substr(v, 1, 87) "..."; printf "              %-28s %s\n", $2, v }'
    fi
done

echo
if [ "${NBAD}" -eq 0 ]; then
    echo "All shared settings match."
    exit 0
else
    echo "${NBAD} shared setting(s) differ between scripts. Make every copy match, then re-run."
    exit 1
fi
