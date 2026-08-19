#!/bin/bash
# ============================================================================
# Central path config for the Spirituality (Nielsen lab) multi-echo fMRI pipeline.
#
# Source this at the top of every script in this directory:
#     source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
#
# Use the absolute form above rather than "$(dirname "$0")" — under sbatch, $0 points
# at a copy in the Slurm spool directory, not at this repo.
#
# ---------------------------------------------------------------------------
# Storage strategy (decided 2026-08-19). See ./BYU_ORC_AGENTS.md and
# https://rc.byu.edu/wiki/?id=Storage
#
#   HOME    /home/bradenf4              2 TB. Backed up daily + snapshots. Never purged.
#                                       -> code, unregenerable inputs, containers, R libs
#   SCRATCH /nobackup/autodelete/...    20 TB. NOT backed up. PURGED after 12 weeks UNUSED.
#                                       -> all job I/O: BIDS, fmriprep, tedana, AFNI, work dirs
#   ARCHIVE /nobackup/archive/...       20 TB. NOT backed up. Never purged. 1 M file limit.
#                                       -> raw DICOMs (read-only) + promoted final products
#
# BYU RC storage docs: "Archive storage should generally NOT be used directly from
# batch jobs." The archive is also mid-migration to a much slower backing store.
# Nothing in the compute path may read or write ARCHIVE except RAWDATA_DIR (read-only).
#
# The 12-week SCRATCH purge is triggered by files being UNUSED, not by age. A parked
# project loses its derivatives. Use ./promote_to_archive.sh before setting work down.
# ============================================================================

# ---- Tier roots (env-overridable, so the pipeline can be render-tested against a fixture) ----
export HOME_ROOT=${HOME_ROOT:-/home/bradenf4}
export SCRATCH_ROOT=${SCRATCH_ROOT:-/nobackup/autodelete/usr/bradenf4}
export ARCHIVE_ROOT=${ARCHIVE_ROOT:-/nobackup/archive/usr/bradenf4}

# ---- Project roots ----
export PROJECT_NAME=Spirituality
export PROJ_HOME=${HOME_ROOT}/spirituality_fmri
export PROJ_SCRATCH=${SCRATCH_ROOT}/${PROJECT_NAME}
export PROJ_ARCHIVE=${ARCHIVE_ROOT}/Nielsen_active/Spirituality/Project

# ---- Source data: ARCHIVE, READ-ONLY. Never write here. ----
export RAWDATA_DIR=${PROJ_ARCHIVE}/rawdata
export PARTICIPANTS_TSV=${RAWDATA_DIR}/participants.tsv
export RUN_ORDER_TSV=${RAWDATA_DIR}/run_order_plain.tsv

# ---- Unregenerable inputs: HOME (backed up) ----
# The *_ratingAM.1D timing files are produced by generate_timing_files.py on Braden's
# LOCAL machine — that generator is not in this repo and does not run on the cluster.
# They must never live on a purging filesystem. Upload new ones straight to here.
export TIMING_DIR=${PROJ_HOME}/timing

# ---- All compute I/O: SCRATCH ----
export BIDS_DIR=${PROJ_SCRATCH}/BIDS
export DERIV=${PROJ_SCRATCH}/derivatives
export FMRIPREP_OUT=${DERIV}/fmriprep
export TEDANA_OUT=${DERIV}/tedana
export AFNI_OUT=${DERIV}/afni_firstlvl
export REPRO_OUT=${DERIV}/afni_reproducibility
export GROUP_OUT=${DERIV}/afni_group
export BIDS_FILTER_DIR=${DERIV}/bids_filters
export WORK_ROOT=${SCRATCH_ROOT}/work/${PROJECT_NAME}

# ---- Final products: ARCHIVE (+ offsite copy) ----
# Written only by promote_to_archive.sh, never by a pipeline job.
export FINAL_DIR=${PROJ_ARCHIVE}/final

# ---- Logs: HOME (small, backed up) ----
export LOGS=${PROJ_HOME}/logs

# ---- Software: HOME ----
# Containers and TemplateFlow moved off archive 2026-08-19 — every job start was pulling
# 2.4 GB off the tier that is about to get much slower.
export FMRIPREP_SIF=${HOME_ROOT}/software/fmriprep-25.1.4.sif
export TEMPLATEFLOW_HOME=${HOME_ROOT}/software/templateflow
export FS_LICENSE=${HOME_ROOT}/.freesurfer_license.txt
export AFNI_SIF=/apps/afni/afni_make_build_latest.sif
export CONDA_ENV=tedenv
# R library for 3dICC / spatial ICC (irr, lme4, blme, metafor, snow).
# On HOME, not scratch: compute nodes have NO internet, so a purged Rlib cannot be
# rebuilt from inside a job — the install has to be redone from a login node.
export RLIB=${HOME_ROOT}/software/Rlib

# ---- Apptainer binds ----
# Must span EVERY tier a container touches. Inputs are on SCRATCH, raw data on ARCHIVE,
# containers and timing files on HOME — a single-root bind silently fails inside the
# container rather than erroring loudly.
export BIND_ARGS="--bind ${HOME_ROOT}:${HOME_ROOT} --bind ${SCRATCH_ROOT}:${SCRATCH_ROOT} --bind ${ARCHIVE_ROOT}:${ARCHIVE_ROOT}"

# ---- Guard against user-site package poisoning ----
# A stray numpy in ~/.local shadows the tedenv conda env for anything not run in a
# container. This is what broke tedana in July 2026 (np.row_stack removed in numpy 2.5).
export PYTHONNOUSERSITE=1

# ---- Create output dirs on demand (SCRATCH and HOME only — never ARCHIVE) ----
mkdir -p "${DERIV}" "${WORK_ROOT}" "${BIDS_FILTER_DIR}" "${LOGS}"
