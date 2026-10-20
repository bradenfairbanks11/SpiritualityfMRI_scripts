# CLAUDE.md — Spirituality fMRI scripts (Nielsen lab)

## BYU HPC agent instructions (read first)

These scripts run on BYU Office of Research Computing HPC systems. Before doing significant
work here, read the local copy of the ORC agent instructions:

    ./BYU_ORC_AGENTS.md

Upstream version: 2026-08-11, synced 2026-08-21 (verified byte-identical to upstream) from
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md`
(mirror: https://rc.byu.edu/documentation/BYU_ORC_AGENTS.md).

Check the copy's mtime; if it is more than 7 days old, refresh it from `/apps` (preferred)
or the URL, then continue. If neither is reachable, continue with the existing copy.

Those instructions outrank this file. The stated priority order is:
ORC policy and administrator instructions > `BYU_ORC_AGENTS.md` > this file > the user's
request > general agent defaults.

## Storage tiers — read `config.sh` before touching any path

Every path in this pipeline comes from `config.sh`. Do not hardcode paths in scripts; add
them there. The tiers:

| Tier | Purpose | Policy |
|---|---|---|
| `/home/bradenf4` | code, timing files, containers, R libs | backed up daily, snapshots, never purged |
| `/nobackup/autodelete/usr/bradenf4` | **all** job I/O | not backed up, **purged after 12 weeks unused** |
| `/nobackup/archive/usr/bradenf4` | raw DICOMs (read-only) + promoted final products | not backed up, never purged, 1 M file limit |

**Do not compute out of archive.** BYU RC storage docs: *"Archive storage should generally
NOT be used directly from batch jobs."* The archive is migrating to a much slower store.
The only archive paths in the compute path are `RAWDATA_DIR` (read-only) and `FINAL_DIR`
(written solely by `promote_to_archive.sh`).

**The scratch purge trigger is disuse, not age.** Run `promote_to_archive.sh` before parking
this project, or the derivatives — 10–12 h of fMRIPrep per subject — are deleted.

## Traps that have already bitten this repo

- **`#SBATCH` directives are parsed by Slurm before the shell runs.** Shell variables do not
  expand there. Log paths in directly-submitted scripts must be literal absolutes.
- **Under `sbatch`, BOTH `$0` and `${BASH_SOURCE[0]}` point at Slurm's spool copy**
  (`/var/spool/slurmd/job<N>/slurm_script`) — not at this repo. Verified 2026-08-21
  with probe job 13292417. Earlier guidance here said `${BASH_SOURCE[0]}` was the fix
  for `$0`; it is not, and every directly-submitted job script died on its
  `source config.sh` line until this was corrected. Use the `CONFIG=""` resolver
  block the scripts now carry: `$PIPELINE_CONFIG`, then `$SLURM_SUBMIT_DIR`, then
  `dirname ${BASH_SOURCE[0]}` (right only when run with `bash`), then the install path.
- **Apptainer `--bind` must span all three tiers** (`${BIND_ARGS}`). Inputs are on scratch,
  raw data on archive, containers and timing files on home. A single-root bind fails
  *silently* inside the container.
- **Timing files cannot be regenerated on the cluster.** `generate_timing_files.py` runs on
  Braden's local machine. New `*_ratingAM.1D` **and `*_unrated.1D`** files go straight to
  `~/spirituality_fmri/timing/` (upload the whole folder; a missing `_unrated.1D` silently
  puts that run's no-response trials back into the baseline).
  `first_level_afni.sh` silently skips runs with no timing file, so a missing one looks like
  a clean run that did nothing.
- **Compute nodes have no internet.** Any pip/R install must be done from a login node. This
  is why `RLIB` lives on `/home`.
- **`~/.local` poisons every conda env.** A stray numpy there shadowed `tedenv` and broke
  tedana (`np.row_stack` removed in numpy 2.5). `config.sh` sets `PYTHONNOUSERSITE=1`.
- **Always re-run `assign_fieldmaps.py` after any dcm2niix re-conversion** — dcm2niix writes
  a coarse blanket `IntendedFor` that clobbers per-run fieldmap pairing.

## Slurm conventions

Every job specifies nodes, cores, memory, and time. Do **not** add `--partition`, `--qos`, or
`--constraint` unless something genuinely requires it — each shrinks the eligible node pool
and lengthens the queue. Wait ≥60 s between `squeue`/`sacct` checks; prefer `--dependency`
(as the fMRIPrep→tedana chain already does) over polling loops.

## Filesystem hygiene

`derivatives/` is BIDS + fMRIPrep + tedana + AFNI — exactly the many-small-files shape the
shared filesystems are sensitive to. Do not run `find`, `du -a`, or `ls -lR` over `rawdata/`
or `derivatives/`; scope listings to a specific `sub-XX/ses-Y/` path.

## CUI

Human-subjects neuroimaging. Confirmed with Braden on 2026-08-19 that this project is not
CUI or export-controlled. That confirmation does not extend to new projects.

## Project notes

Multi-echo fMRI, parametric "felt the Spirit" rating GLM (AFNI `-stim_times_AM2`); the
`rating_slope` sub-brick is the scientific target. Pipeline: `dcm2niix.sh` →
`assign_fieldmaps.py` → `fmriprep_tedana.sh` → `first_level_afni.sh` →
`second_level_afni.sh` / `reproducibility_afni.sh`, then `promote_to_archive.sh`.

Trial-level reliability branch (added 2026-10-20): `first_level_IM_afni.sh` (one beta per
trial via `-stim_times_IM`, reusing first_level_afni.sh's scaled data + motion files) →
`trr_export_afni.sh` (per-trial NIfTIs, data tables, group-cluster ROIs into `${TRR_OUT}`) →
copy `afni_trr/` + `afni_reproducibility/` to the local Mac → `local_r_steps.sh`.

- **R-based AFNI programs do not run on the cluster** (3dICC fails in the container), so
  3dICC, `3dLMEr -TRR` and `TRR` run locally via `local_r_steps.sh`. Cluster scripts must
  stay pure-AFNI-C for anything new; `reproducibility_afni.sh` still *attempts* 3dICC and
  records NA, leaving its inputs in `afni_reproducibility/group/` for the local run.
- **`3dLMEr -TRR` hard-codes exactly 2 conditions × 2 sessions** (`runTRR` in 3dLMEr.R).
  Trials are split by rating: low = 1–2, high = 3–4 (`HIGH_MIN` in first_level_IM_afni.sh).
  FHS is excluded from the voxelwise run (too few high-rated trials); region-level `TRR`
  handles a single condition and covers all three tasks.
- **AFNI's TRR.R pre-sizes its posterior vector to 2000 draws.** Keep
  chains × iterations / 2 ≥ 2000 or the reported TRR is biased toward 0.
- `dice_clust` = cluster-corrected Dice (voxel p<.001, per-session ACF 3dClustSim, α=.05).
  Both-empty map pairs score Dice 0 and within- vs between-subject Dice is tested with
  `dice_within_vs_between.R` — both per Gorgolewski et al. 2013 (NeuroImage 69:231).
- Not yet in `promote_to_archive.sh`: `afni_trr/` and `afni_firstlvl/*/single_trial/`.

Commit script edits to git rather than leaving `.bak` copies.
