#!/usr/bin/env python3
"""
==============================================================================
qc_carpets.py -- BEFORE/AFTER-TEDANA CARPET PLOTS AND MOTION-COUPLING QC
==============================================================================

WHAT THIS SCRIPT DOES
  tedana (multi-echo ICA denoising) removes noise components from each run.
  This script shows, for every run, how much the data changed, by comparing:
      BEFORE = tedana's desc-optcom_bold   (the "optimal combination" of the
                                            4 echoes, NOT yet denoised)
      AFTER  = tedana's desc-denoised_bold (the same data after the noise
                                            components were removed)
  Both are in the same native (boldref) space and come out of the same tedana
  run, so the only difference between them is what tedana removed.

  Two quality measures are used:
    - FD (framewise displacement, from fMRIPrep's confounds table): how far the
      head moved from one volume to the next, in mm.
    - DVARS: how much the whole-brain signal changes from one volume to the
      next. Here it is "raw" DVARS: the root-mean-square, over all brain
      voxels, of the volume-to-volume signal difference (in the data's own
      units, not standardized).
  If motion is corrupting the signal, DVARS rises whenever FD rises, so the
  correlation between FD and DVARS ("motion-artifact coupling") is high. Good
  denoising should LOWER that correlation and lower DVARS.

  A "carpet plot" (Power 2017 style) is an image with one row per voxel and
  one column per volume; each row is z-scored, so motion shows up as vertical
  bands across many voxels at once.

WHO CALLS IT / WHERE IT SITS IN THE PIPELINE
  No other script calls it. It is run by hand, as an optional QC check, after
  fmriprep_tedana.sh has finished (it needs both the fMRIPrep and the native-
  space tedana outputs). It is not a Slurm job script (no #SBATCH lines).
  It needs a Python with numpy, nibabel, pandas and matplotlib (for example the
  tedenv conda environment used for tedana).

STEPS
  1. Find every run: one per fMRIPrep echo-1 preprocessed BOLD file.
  2. For each run (skipped if any of its 4 input files is missing):
     2a. Load the brain mask and the BEFORE and AFTER data inside the mask.
     2b. Read FD from fMRIPrep's confounds table.
     2c. Compute DVARS for BEFORE and AFTER.
     2d. Correlate FD with DVARS (BEFORE and AFTER) and store summary numbers.
     2e. Draw a 4-panel figure: FD, DVARS before/after, carpet BEFORE,
         carpet AFTER.
  3. Write the summary table (one row per processed run).
  4. Draw the dataset-wide summary bar figure.
  5. Print the summary table to the screen.

INPUTS (per run; <pre> = sub-NN_ses-N_task-T)
  <DERIV>/fmriprep/sub-NN/ses-N/func/<pre>_echo-1_desc-preproc_bold.nii.gz
                                       (only used to find the runs)
  <DERIV>/fmriprep/sub-NN/ses-N/func/<pre>_desc-brain_mask.nii.gz
  <DERIV>/fmriprep/sub-NN/ses-N/func/<pre>_desc-confounds_timeseries.tsv
                                       (column "framewise_displacement")
  <DERIV>/tedana/sub-NN/ses-N/func/<pre>_desc-optcom_bold.nii.gz     (BEFORE)
  <DERIV>/tedana/sub-NN/ses-N/func/<pre>_desc-denoised_bold.nii.gz   (AFTER)
  DERIV defaults to /nobackup/autodelete/usr/bradenf4/Spirituality/derivatives

OUTPUTS, in <DERIV>/tedana/qc_carpets/  (the folder is created if needed)
  <pre>_carpet_beforeafter.png     one 4-panel figure per run
  tedana_carpet_qc_summary.tsv     one row per processed run, columns:
      run, sub, mean_FD, meanDVARS_before, meanDVARS_after,
      DVARS_reduction_pct, FDcoupling_before, FDcoupling_after, coupling_drop
  tedana_carpet_qc_summary.png     bar charts of FD-DVARS coupling and % DVARS
                                   reduction, one bar (pair) per run
  Existing files with these names are overwritten.

HOW TO RUN
  python3 qc_carpets.py
  To read from / write to a different derivatives folder, set DERIV first:
  DERIV=/some/other/derivatives python3 qc_carpets.py
  (DERIV is the only path setting read from the environment.)

METHODS SUMMARY (the exact values used below)
  Data compared  : tedana desc-optcom_bold (before) vs desc-denoised_bold
                   (after), native boldref space, inside fMRIPrep's run-level
                   brain mask (mask voxels > 0)
  FD             : fMRIPrep "framewise_displacement" (first volume is n/a ->
                   NaN), cut to the number of BOLD volumes
  DVARS          : raw DVARS = sqrt(mean over mask voxels of (x[t] - x[t-1])^2),
                   no intensity normalization; first volume = NaN
  Coupling       : Pearson r between FD and DVARS over volumes where FD and
                   both DVARS series are all defined (numpy corrcoef)
  DVARS reduction: 100 x (1 - mean DVARS after / mean DVARS before)
  coupling_drop  : r(before) - r(after)
  Carpet display : each voxel's time series z-scored (a voxel with zero
                   standard deviation is divided by 1 instead); 3000 voxels
                   drawn at random without replacement (numpy default_rng,
                   seed 42; the same voxels for before and after), kept in
                   mask order; grey scale from z = -2.5 to +2.5
  FD reference   : dashed line drawn at 0.5 mm (display only)
  Figures        : per-run 13 x 9 in at 110 dpi; summary 14 x 5 in at 120 dpi
==============================================================================
"""

# ==============================================================================
# LIBRARIES
# ------------------------------------------------------------------------------
# numpy (arrays), nibabel (reads NIfTI images), pandas (tables),
# matplotlib (figures). The "Agg" backend draws figures straight to files, so
# no screen/display is needed (important on the supercomputer).
# All Python warnings are silenced so the screen output stays readable.
# ==============================================================================
import os, glob, warnings
import numpy as np
import nibabel as nib
import pandas as pd
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore")


# ==============================================================================
# SECTION 1. SETTINGS
# ------------------------------------------------------------------------------
# DERIV is the derivatives folder on scratch (/nobackup/autodelete, purged after
# 12 weeks unused). Its default is written below. If an environment variable
# called DERIV is set when the script starts, that value is used instead.
# The other three folders are always built from DERIV, so with the default
# they are:
#   FMRIPREP = /nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/fmriprep
#   TEDANA   = /nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/tedana
#   OUTDIR   = /nobackup/autodelete/usr/bradenf4/Spirituality/derivatives/tedana/qc_carpets
# The output folder is created here, as soon as the script starts.
#
# RNG: the random-number generator used to pick which voxels are drawn in the
# carpet plots. Fixed seed 42, so re-running on the same set of runs picks the
# same voxels. (One generator is shared by all runs in turn, so each run's
# voxel choice also depends on how many runs came before it.)
# N_CARPET_ROWS: how many voxels (rows) to show in each carpet plot. Showing
# every brain voxel would not fit on the image, so 3000 are sampled.
# This affects only the pictures, not the DVARS/FD numbers.
# ==============================================================================
DERIV = os.environ.get(
    "DERIV",
    "/nobackup/autodelete/usr/bradenf4/Spirituality/derivatives",
)
FMRIPREP = f"{DERIV}/fmriprep"
TEDANA = f"{DERIV}/tedana"
OUTDIR = f"{DERIV}/tedana/qc_carpets"
os.makedirs(OUTDIR, exist_ok=True)
RNG = np.random.default_rng(42)
N_CARPET_ROWS = 3000  # subsample voxels for display


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

def load_masked(path, mask):
    """Load a 4-D BOLD image and keep only the voxels inside the brain mask.

    Returns a 2-D array: one row per in-mask voxel, one column per volume
    (time point), as 32-bit floats. The whole image is read into memory first.
    """
    img = nib.load(path)
    data = np.asarray(img.dataobj, dtype=np.float32)  # whole 4-D image read into memory
    return data[mask]


def dvars(ts):
    """Raw DVARS for every volume.

    ts = voxels x volumes. For each pair of neighbouring volumes, take the
    signal difference in every voxel, square it, average over voxels, and take
    the square root. There is no difference for the first volume, so it is set
    to NaN (missing); the result has one value per volume.
    """
    d = np.diff(ts, axis=1)
    dv = np.sqrt(np.mean(d ** 2, axis=0))
    return np.concatenate([[np.nan], dv])  # align to T, first=nan


def zscore_rows(ts):
    """Z-score each voxel's time series (subtract its mean, divide by its SD).

    Used only for the carpet pictures, so every voxel is on the same scale.
    A voxel whose SD is 0 (flat signal) is divided by 1 instead, to avoid
    dividing by zero.
    """
    mu = ts.mean(axis=1, keepdims=True)
    sd = ts.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (ts - mu) / sd


def find_runs():
    """List every run and the paths of its input files.

    A run is found from fMRIPrep's echo-1 preprocessed BOLD file
    (sub-*/ses-*/func/*_echo-1_desc-preproc_bold.nii.gz). From its name the
    prefix sub-NN_ses-N_task-T is taken, and the BEFORE, AFTER, mask and
    confounds paths are built from it. The files are NOT checked here; that
    happens in main(). Runs come out in alphabetical order.
    """
    runs = []
    for e1 in sorted(glob.glob(f"{FMRIPREP}/sub-*/ses-*/func/*_echo-1_desc-preproc_bold.nii.gz")):
        base = os.path.basename(e1)
        pre = base.split("_echo-1_")[0]           # sub-XX_ses-Y_task-Z
        sub = pre.split("_")[0]
        ses = pre.split("_")[1]
        func_ted = f"{TEDANA}/{sub}/{ses}/func"
        func_fp = f"{FMRIPREP}/{sub}/{ses}/func"
        runs.append(dict(
            key=pre, sub=sub, ses=ses,
            before=f"{func_ted}/{pre}_desc-optcom_bold.nii.gz",
            after=f"{func_ted}/{pre}_desc-denoised_bold.nii.gz",
            mask=f"{func_fp}/{pre}_desc-brain_mask.nii.gz",
            confounds=f"{func_fp}/{pre}_desc-confounds_timeseries.tsv",
        ))
    return runs


# ==============================================================================
# SECTION 3. MAIN PROCEDURE
# ==============================================================================

def main():
    """Run the QC for every run, then write the summary table and figure."""

    # --------------------------------------------------------------------------
    # STEP 1. FIND THE RUNS
    # --------------------------------------------------------------------------
    runs = find_runs()
    print(f"Found {len(runs)} runs")
    summary = []
    for r in runs:
        # ----------------------------------------------------------------------
        # STEP 2. ONE RUN AT A TIME
        # If any of the 4 input files is missing, print SKIP and go to the next
        # run (the "for ... else" below: the else part runs only when no file
        # was missing).
        # ----------------------------------------------------------------------
        for p in (r["before"], r["after"], r["mask"], r["confounds"]):
            if not os.path.exists(p):
                print(f"  SKIP {r['key']}: missing {os.path.basename(p)}")
                break
        else:
            print(f"  processing {r['key']} ...", flush=True)

            # ------------------------------------------------------------------
            # STEP 2a. LOAD THE MASK AND THE BEFORE/AFTER DATA
            # Brain voxels = mask value > 0. T = number of volumes.
            # ------------------------------------------------------------------
            mask = np.asarray(nib.load(r["mask"]).dataobj) > 0
            before = load_masked(r["before"], mask)
            after = load_masked(r["after"], mask)
            T = before.shape[1]

            # ------------------------------------------------------------------
            # STEP 2b. READ FRAMEWISE DISPLACEMENT
            # fMRIPrep writes "n/a" for the first volume; that (and anything
            # else non-numeric) becomes NaN. Cut to T values.
            # ------------------------------------------------------------------
            fd = pd.read_csv(r["confounds"], sep="\t")["framewise_displacement"]
            fd = pd.to_numeric(fd, errors="coerce").to_numpy()[:T]

            # ------------------------------------------------------------------
            # STEP 2c. DVARS BEFORE AND AFTER DENOISING
            # ------------------------------------------------------------------
            dv_b = dvars(before)
            dv_a = dvars(after)

            # ------------------------------------------------------------------
            # STEP 2d. MOTION-ARTIFACT COUPLING AND SUMMARY NUMBERS
            # Pearson r between FD and DVARS, using only volumes where FD and
            # both DVARS values exist (so the first volume is always dropped).
            # rb = r before denoising, ra = r after. Lower = less motion
            # artifact left in the signal.
            # Stored per run: mean FD, mean DVARS before/after, % reduction in
            # mean DVARS, r before/after, and the drop in r (before - after).
            # Means ignore NaN values.
            # ------------------------------------------------------------------
            valid = ~(np.isnan(fd) | np.isnan(dv_b) | np.isnan(dv_a))
            rb = np.corrcoef(fd[valid], dv_b[valid])[0, 1]
            ra = np.corrcoef(fd[valid], dv_a[valid])[0, 1]

            summary.append(dict(
                run=r["key"], sub=r["sub"], mean_FD=np.nanmean(fd),
                meanDVARS_before=np.nanmean(dv_b), meanDVARS_after=np.nanmean(dv_a),
                DVARS_reduction_pct=100 * (1 - np.nanmean(dv_a) / np.nanmean(dv_b)),
                FDcoupling_before=rb, FDcoupling_after=ra,
                coupling_drop=rb - ra,
            ))

            # ------------------------------------------------------------------
            # STEP 2e. THE PER-RUN FIGURE (4 panels stacked, shared time axis)
            #   Panel 1: FD trace (red), dashed reference line at 0.5 mm.
            #            Title: run name, r before -> after, mean DVARS
            #            before -> after and % reduction.
            #   Panel 2: DVARS before (orange) and after (teal).
            #   Panel 3: carpet BEFORE (optcom).  Panel 4: carpet AFTER (denoised).
            # Carpets: z-score every voxel, then pick the same 3000 random
            # voxels (sorted back into mask order) for both carpets; grey
            # scale from z = -2.5 (black) to +2.5 (white).
            # ------------------------------------------------------------------
            cb = zscore_rows(before)
            ca = zscore_rows(after)
            idx = np.sort(RNG.choice(cb.shape[0], min(N_CARPET_ROWS, cb.shape[0]), replace=False))
            cb, ca = cb[idx], ca[idx]

            fig = plt.figure(figsize=(13, 9))
            gs = GridSpec(4, 1, height_ratios=[1.1, 1.1, 3, 3], hspace=0.08)
            t = np.arange(T)

            ax0 = fig.add_subplot(gs[0])
            ax0.plot(t, fd, color="crimson", lw=0.8)
            ax0.axhline(0.5, color="gray", ls="--", lw=0.6)
            ax0.set_ylabel("FD (mm)")
            ax0.set_xlim(0, T); ax0.set_xticklabels([])
            ax0.set_title(f"{r['key']}   |   FD-DVARS coupling: "
                          f"{rb:.2f} → {ra:.2f}   |   mean DVARS "
                          f"{np.nanmean(dv_b):.1f} → {np.nanmean(dv_a):.1f} "
                          f"({summary[-1]['DVARS_reduction_pct']:.0f}% ↓)",
                          fontsize=11, loc="left")

            ax1 = fig.add_subplot(gs[1])
            ax1.plot(t, dv_b, color="darkorange", lw=0.8, label="before")
            ax1.plot(t, dv_a, color="teal", lw=0.8, label="after")
            ax1.set_ylabel("DVARS"); ax1.set_xlim(0, T); ax1.set_xticklabels([])
            ax1.legend(loc="upper right", fontsize=8, ncol=2)

            vmax = 2.5
            ax2 = fig.add_subplot(gs[2])
            ax2.imshow(cb, aspect="auto", cmap="gray", vmin=-vmax, vmax=vmax,
                       interpolation="nearest", extent=[0, T, 0, cb.shape[0]])
            ax2.set_ylabel("BEFORE\n(optcom)"); ax2.set_xticklabels([]); ax2.set_yticks([])

            ax3 = fig.add_subplot(gs[3])
            ax3.imshow(ca, aspect="auto", cmap="gray", vmin=-vmax, vmax=vmax,
                       interpolation="nearest", extent=[0, T, 0, ca.shape[0]])
            ax3.set_ylabel("AFTER\n(denoised)"); ax3.set_yticks([])
            ax3.set_xlabel("volume (TR)")

            out = f"{OUTDIR}/{r['key']}_carpet_beforeafter.png"
            fig.savefig(out, dpi=110, bbox_inches="tight")
            plt.close(fig)
            del before, after, cb, ca      # free memory before the next run

    # --------------------------------------------------------------------------
    # STEP 3. WRITE THE SUMMARY TABLE (tab-separated, one row per processed run)
    # --------------------------------------------------------------------------
    df = pd.DataFrame(summary)
    df.to_csv(f"{OUTDIR}/tedana_carpet_qc_summary.tsv", sep="\t", index=False)

    # --------------------------------------------------------------------------
    # STEP 4. THE DATASET-WIDE SUMMARY FIGURE
    #   Left : FD-DVARS correlation per run, before (orange) vs after (teal).
    #   Right: % reduction in mean DVARS per run.
    # --------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    x = np.arange(len(df)); w = 0.38
    axes[0].bar(x - w/2, df.FDcoupling_before, w, label="before", color="darkorange")
    axes[0].bar(x + w/2, df.FDcoupling_after, w, label="after", color="teal")
    axes[0].set_xticks(x); axes[0].set_xticklabels(df.run, rotation=90, fontsize=7)
    axes[0].set_ylabel("corr(FD, DVARS)")
    axes[0].set_title("Motion-artifact coupling (lower = cleaner)"); axes[0].legend()

    axes[1].bar(x, df.DVARS_reduction_pct, color="seagreen")
    axes[1].set_xticks(x); axes[1].set_xticklabels(df.run, rotation=90, fontsize=7)
    axes[1].set_ylabel("% reduction in mean DVARS")
    axes[1].set_title("Signal-fluctuation reduction after tedana")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/tedana_carpet_qc_summary.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    # --------------------------------------------------------------------------
    # STEP 5. PRINT THE SUMMARY TABLE (rounded to 3 decimals) AND OUTPUT FOLDER
    # --------------------------------------------------------------------------
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print("\n=== SUMMARY ===")
    print(df.round(3).to_string(index=False))
    print(f"\nOutputs -> {OUTDIR}")


if __name__ == "__main__":
    main()
