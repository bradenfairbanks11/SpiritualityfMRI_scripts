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

## Settings live INSIDE each script (since 2026-10-24) — and must agree

Braden's rule: every analysis script must be readable on its own (for writing the methods
paper), so dcm2niix.sh, first_level_afni.sh, first_level_IM_afni.sh, second_level_afni.sh,
trr_export_afni.sh, reproducibility_afni.sh, reliability_models_afni.sh and check_trial_censoring.sh define every path
and parameter in their own SECTION 1 and do NOT source config.sh. Settings used by more than
one script carry a trailing `# [shared]` and must be byte-identical across copies: after
changing one, edit every copy and run `bash check_settings.sh` (exits 1 on any mismatch).
Layout every script follows (first_level_afni.sh is the template): SLURM box → header
(what / pipeline position / steps / inputs / outputs / how to run / METHODS SUMMARY) →
SECTION 1 settings → 2 helpers → 3 setup/checks → 4 main steps, each step boxed and explained
in plain language. Keep new scripts and edits in that style.
`config.sh` is still sourced by fmriprep_tedana.sh, tedana_only.sh and promote_to_archive.sh
(Braden chose to leave those three as they are). The storage tiers:

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
  `source config.sh` line until this was corrected. Scripts that still source config.sh
  use the `CONFIG=""` resolver block; self-contained scripts use literal absolute paths
  (e.g. exclusions.tsv at /home/bradenf4/spirituality_fmri/scripts/exclusions.tsv).
- **Apptainer `--bind` must span all three tiers** (`${BIND_ARGS}`). Inputs are on scratch,
  raw data on archive, containers and timing files on home. A single-root bind fails
  *silently* inside the container.
- **Timing files cannot be regenerated on the cluster.** `generate_timing_files.py` runs on
  Braden's local machine. New `*_ratingAM.1D`, `*_unrated.1D` **and `*_stimuli.tsv`** files go
  straight to `~/spirituality_fmri/timing/` (upload the whole folder; a missing `_unrated.1D`
  silently puts that run's no-response trials back into the baseline, and a missing
  `_stimuli.tsv` makes trr_export_afni.sh skip the task).
  `first_level_afni.sh` silently skips runs with no timing file, so a missing one looks like
  a clean run that did nothing.
- **Compute nodes have no internet.** Any pip/R install must be done from a login node.
  R packages now live INSIDE the AFNI image (below), so there is nothing to install for R.
- **`~/.local` poisons every conda env.** A stray numpy there shadowed `tedenv` and broke
  tedana (`np.row_stack` removed in numpy 2.5). Every script that runs host Python sets
  `export PYTHONNOUSERSITE=1`.
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
`reliability_models_afni.sh` (cluster, inside the AFNI image: 3dLMEr -TRR as
`sbatch --array=0-3 ... trr_vox`, one model per array task; then `... trr_roi dice`).
It replaced the Mac-only `local_r_steps.sh` on 2026-10-24, once the new image could run R.

- **TRR conditions are stimulus labels (since 2026-10-24), not the rating split.**
  `stimulus_conditions.tsv` (repo; copy to `~/spirituality_fmri/scripts/`) labels every stimulus
  yield = high/low (high-yield list from a separate pre-study survey; source of record is the
  RTF in Box `data/Inter-session_behavioral_ICC_results/high_yield_stimuli/`) and, for
  architecture, faith = lds/other (file-name prefix). trr_export_afni.sh joins each trial to
  its stimulus via `_IM_trials.tsv` orig_trial = row of `<timing>/byv-*_stimuli.tsv` (onsets
  must match, else the task stops). reliability_models_afni.sh `TRR_VOX_RUNS` runs 3dLMEr -TRR for
  scripture/FHS/architecture yield + architecture faith. Scanner check (2026-10-24): yield
  separates in-scanner ratings strongly for scripture/FHS (~1.4-1.9 points, every subject) but
  weakly for architecture (+0.25-0.6), where only the LDS images follow the survey.

- **AFNI image: `/apps/afni/26.3.04.sif` (since 2026-10-24).** ORC built it to the spec in
  `container_request.md`: Ubuntu 24.04, R >= 4.4, every R package inside the image,
  `R_LIBS_USER=NULL`. R-based AFNI calls use `apptainer exec --cleanenv` so host env cannot
  override the image's packages. The old image (`afni_make_build_latest.sif`: Ubuntu 18.04,
  R 3.4.4, AFNI 24.3.01) could not run 3dICC/3dLMEr/TRR; `RLIB`/`setup_rlib.sh` were its
  workaround and are now unused. Do not mix outputs from the two AFNI versions in one analysis.
- **Motion censoring = fMRIPrep `framewise_displacement` > 0.3 mm, plus the volume before**
  (first_level_afni.sh STEP 4, since 2026-10-24; Braden's decision). Before that the code
  censored at 1d_tool.py euclidean-norm >= 1 on mm + RADIANS (rotations barely counted) while
  its comment claimed "FD > 0.3mm". Outputs made before 2026-10-24 used the old rule.
- **All-zero regressors crash 3dREMLfit** ("Cannot continue with all zero column"), and with
  `set -e` that ends the job for every remaining run. Guarded since 2026-10-24:
  first_level_IM_afni.sh drops trials with signal_kept <= MIN_SIGNAL_KEPT (= 0.3 since 2026-10-24;
  signal_kept = censored / uncensored sum of squares, from -x1D vs -x1D_uncensored; the unrated
  regressor is dropped only at 0) and renumbers
  `_IM_trials.tsv` so row k = sub-brick k (orig_trial keeps the original number);
  first_level_afni.sh does a -nodata design check, drops a fully censored unrated regressor and
  skips a run whose rating[0]/rating[1] column is empty. AFNI design-matrix paths must stay
  under ~250 characters or 3dDeconvolve silently fails to write them.
- **First-level analysis mask = fMRIPrep brain mask AND extents mask** (smoothed data > 0 at every
  time point; afni_proc.py applies the same kind of mask by default), saved as
  `<PREFIX>_analysis_mask.nii.gz` and reused by first_level_IM_afni.sh (since 2026-10-24).
  `1d_tool.py` ignores AFNI_DECONFLICT, so its output calls need `-overwrite` or re-runs crash.
- **second_level_afni.sh group mask = non-zero rating_slope AND mean_response betas in EVERY run**
  that enters the group analysis (per run, not per session-average: STEP 1 averages sessions, so
  a voxel modelled in one session only would look non-zero). Was one fMRIPrep mask per subject
  until 2026-10-24, which left data-less zero voxels inside the mask for 3dttest++ (no -zskip).
- **reproducibility_afni.sh whole-brain mask = non-zero mean_response AND rating_slope betas in
  both sessions** (since 2026-10-24; was 3dAutomask of tedana BOLD, which did not match the GLM
  mask and let zero-beta voxels count as agreement in spatial_r/spatial_icc and enter 3dFWHMx /
  3dClustSim). Group masks elsewhere are already "non-zero in every input" or intersections.
- **trr_export_afni.sh mask_task-T = non-zero single-trial beta for every trial of every included
  run** (3dmask_tool -frac 1.0 on the _IM_betas buckets; since 2026-10-24, was the intersection of
  fMRIPrep brain masks). Used as the 3dLMEr -TRR mask and the 3dClusterize search mask for ROIs.
- 3dICC runs on the cluster inside `reproducibility_afni.sh`; `3dLMEr -TRR`, `TRR` and the
  Dice test run in `reliability_models_afni.sh`. None of these R steps has run on real data in
  the new image yet.
- **`3dLMEr -TRR` hard-codes exactly 2 conditions × 2 sessions** (`runTRR` in 3dLMEr.R).
  The conditions are stimulus labels (yield, faith; see above). The old rating split
  (`HIGH_MIN`, low = 1–2, high = 3–4) is still written to `_IM_trials.tsv` for reference only.
  Region-level `TRR` handles a single condition and covers all three tasks.
- **AFNI's TRR.R pre-sizes its posterior vector to 2000 draws.** Keep
  chains × iterations / 2 ≥ 2000 or the reported TRR is biased toward 0.
- `dice_clust` = cluster-corrected Dice (voxel p<.001, per-session ACF 3dClustSim, α=.05).
  Both-empty map pairs score Dice 0 and within- vs between-subject Dice is tested with
  `dice_within_vs_between.R` — both per Gorgolewski et al. 2013 (NeuroImage 69:231).
- Not yet in `promote_to_archive.sh`: `afni_trr/` and `afni_firstlvl/*/single_trial/`.

## Excluded runs — `exclusions.tsv` (added 2026-10-20)

Runs that must not enter ANY analysis are listed in `exclusions.tsv` (byv_id, session,
task or `*`, reason, date). Each self-contained script carries its own copy of
`is_excluded <sub-NN|BYV> <ses> <task>` (marked [shared]; config.sh still sources
exclusions.sh for the three scripts that use it); every first-level loop and every downstream
discovery loop (second level, reproducibility, TRR export) calls it. Old outputs of an
excluded run stay on disk but are ignored. Add new exclusions there with a reason — never
by deleting timing files (downstream scripts discover from disk, so that does not exclude).
Current: 1VK4L2 ses-1, all tasks (asleep).

Commit script edits to git rather than leaving `.bak` copies.
