#!/usr/bin/env python3
"""
Before/after tedana carpet-plot QC for the Spirituality multi-echo dataset.

For each run compares:
    BEFORE = tedana desc-optcom_bold  (optimal combination, pre-denoise)
    AFTER  = tedana desc-denoised_bold (post ICA denoise)
both in the SAME native boldref space, so the only difference is the
noise components tedana removed.

Outputs, per run, a PNG with:
    - framewise displacement (FD) trace
    - DVARS before vs after
    - carpet plot before
    - carpet plot after
Plus a dataset-wide summary table (TSV) and bar figure quantifying the
reduction in motion-artifact coupling (corr between FD and DVARS).
"""
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

DERIV = "/nobackup/archive/usr/bradenf4/Nielsen_active/Spirituality/Project/derivatives"
FMRIPREP = f"{DERIV}/fmriprep"
TEDANA = f"{DERIV}/tedana"
OUTDIR = f"{DERIV}/tedana/qc_carpets"
os.makedirs(OUTDIR, exist_ok=True)
RNG = np.random.default_rng(42)
N_CARPET_ROWS = 3000  # subsample voxels for display


def load_masked(path, mask):
    """Return (nvox_in_mask, T) float32 array of in-mask voxel time series."""
    img = nib.load(path)
    data = np.asarray(img.dataobj, dtype=np.float32)  # 698G free -> safe
    return data[mask]


def dvars(ts):
    """Raw DVARS: RMS over voxels of the temporal difference. ts=(nvox,T)."""
    d = np.diff(ts, axis=1)
    dv = np.sqrt(np.mean(d ** 2, axis=0))
    return np.concatenate([[np.nan], dv])  # align to T, first=nan


def zscore_rows(ts):
    mu = ts.mean(axis=1, keepdims=True)
    sd = ts.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (ts - mu) / sd


def find_runs():
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


def main():
    runs = find_runs()
    print(f"Found {len(runs)} runs")
    summary = []
    for r in runs:
        for p in (r["before"], r["after"], r["mask"], r["confounds"]):
            if not os.path.exists(p):
                print(f"  SKIP {r['key']}: missing {os.path.basename(p)}")
                break
        else:
            print(f"  processing {r['key']} ...", flush=True)
            mask = np.asarray(nib.load(r["mask"]).dataobj) > 0
            before = load_masked(r["before"], mask)
            after = load_masked(r["after"], mask)
            T = before.shape[1]

            fd = pd.read_csv(r["confounds"], sep="\t")["framewise_displacement"]
            fd = pd.to_numeric(fd, errors="coerce").to_numpy()[:T]

            dv_b = dvars(before)
            dv_a = dvars(after)

            # motion-artifact coupling: corr(FD, DVARS) over valid frames
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

            # ---- figure ----
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
            del before, after, cb, ca

    df = pd.DataFrame(summary)
    df.to_csv(f"{OUTDIR}/tedana_carpet_qc_summary.tsv", sep="\t", index=False)

    # summary bar figure: FD-DVARS coupling before vs after
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

    pd.set_option("display.width", 200, "display.max_columns", 20)
    print("\n=== SUMMARY ===")
    print(df.round(3).to_string(index=False))
    print(f"\nOutputs -> {OUTDIR}")


if __name__ == "__main__":
    main()
