#!/bin/bash
# =============================================================================
# One-time (and top-up) install of the R packages that 3dMEMA, 3dICC, and the
# spatial-ICC helper need, into ${RLIB} on HOME.
#
# BYU HPC agent/operating instructions: see ./BYU_ORC_AGENTS.md in this repo.
#
# *** RUN THIS ON A LOGIN NODE. IT WILL NOT WORK ANYWHERE ELSE. ***
#   Compute nodes have NO internet access, so install.packages() cannot reach
#   CRAN from inside a job. second_level_afni.sh and reproducibility_afni.sh used
#   to attempt this install inline, in the job body — which could only ever have
#   succeeded on the rare occasion the packages were already present. They now
#   check and fail fast with a pointer here instead.
#
# WHY ${RLIB} IS ON HOME: /home is backed up and never purged. Scratch is purged
# after 12 weeks unused, and a purged R library cannot be rebuilt from inside a
# job for the same no-internet reason.
#
# This is small, serial, network-bound work — a login node is the right place for
# it and Slurm would be the wrong one.
#
# USAGE:
#   bash setup_rlib.sh
# =============================================================================
set -uo pipefail

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

# 3dMEMA needs snow + data.table; 3dICC needs lme4/blme/metafor/snow; the spatial
# ICC helper needs irr. Union, installed in one pass.
PKGS='c("snow","data.table","irr","lme4","blme","metafor")'

echo "setup_rlib.sh   host=$(hostname)   started $(date)"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo
    echo "REFUSING TO RUN: this is inside Slurm job ${SLURM_JOB_ID}."
    echo "Compute nodes have no internet. Run it from a login node:"
    echo "    bash ${PROJ_HOME}/scripts/setup_rlib.sh"
    exit 1
fi

mkdir -p "${RLIB}"
module load apptainer/1.3.6-qycanb2

R_IN_AFNI="apptainer exec ${BIND_ARGS} --env R_LIBS_USER=${RLIB} ${AFNI_SIF} bash -c"

echo "library : ${RLIB}"
echo "packages: ${PKGS}"
echo

${R_IN_AFNI} "Rscript -e '
  options(repos = c(CRAN = \"https://cloud.r-project.org\"))
  want    <- ${PKGS}
  lib     <- Sys.getenv(\"R_LIBS_USER\")
  missing <- setdiff(want, rownames(installed.packages()))
  if (!length(missing)) {
      cat(\"All\", length(want), \"packages already present in\", lib, \"\n\")
  } else {
      cat(\"Installing:\", paste(missing, collapse = \", \"), \"\n\")
      install.packages(missing, lib = lib)
  }
'"

echo
echo "--- verifying ---"
${R_IN_AFNI} "Rscript -e '
  want    <- ${PKGS}
  missing <- setdiff(want, rownames(installed.packages()))
  if (length(missing)) {
      cat(\"STILL MISSING:\", paste(missing, collapse = \", \"), \"\n\"); quit(status = 1)
  }
  cat(\"OK - all\", length(want), \"packages installed.\n\")
'"
rc=$?

echo
if [ ${rc} -eq 0 ]; then
    echo "Done. second_level_afni.sh and reproducibility_afni.sh can now run."
else
    echo "*** Install did not complete. Check the CRAN errors above. ***"
    echo "If R itself is too old for a package, see the 3dICC notes — the AFNI"
    echo "container ships an old R and that is a known constraint."
fi
exit ${rc}
