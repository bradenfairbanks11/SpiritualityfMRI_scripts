# QA checklist — Spirituality multi-echo pipeline

How to tell if preprocessing went wrong, stage by stage. Baselines in **bold**
come from the known-good sub-03/ses-2 run. Apply this to every new subject
(next up: sub-04, sub-05).

Pipeline: `dcm2niix` → `assign_fieldmaps.py` → **fMRIPrep 25.1.4** (per-echo) →
**tedana** ICA denoise → **antsApplyTransforms** warp to MNI → AFNI first-level.

---

## 0. Did the jobs even run? (the most basic failure)
```bash
squeue -u bradenf4                                   # still running / pending
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
tail -50 scripts/logs/sub-05_fmriprep.err            # fMRIPrep crash trace
tail -50 scripts/logs/sub-05_tedana.err
```
- `State=COMPLETED ExitCode=0:0` = good. `FAILED`/`OUT_OF_MEMORY`/`TIMEOUT` = bad.
- A tedana job stuck in `PENDING (DependencyNeverSatisfied)` means fMRIPrep failed.
- Classic past failure: **"No BOLD images found"** → the BIDS `func/` was empty
  (exactly what silently happened to sub-04 before re-conversion).

## 1. Conversion / BIDS (most important for FRESH subjects)
A mislabeled or missing file here silently corrupts everything downstream.
```bash
# echoes present per task (expect 1 2 3 4), fmap count (expect 6)
ls BIDS/sub-05/ses-1/func/*task-scripture_echo-*_bold.nii.gz
ls BIDS/sub-05/ses-1/fmap/*_epi.nii.gz | wc -l
```
- **4 echoes** per task (TE 0.0134/0.0313/0.0492/0.0671 s), 3 tasks, 6 fieldmaps.
- **Task labels must match the intended run order** (sub-04 = AFS, sub-05 = SFA
  from `run_order_plain.tsv`). A wrong label misaligns the GLM timing → confirm
  the per-run task names in the fMRIPrep report match what you expect.
- Fieldmaps must be an **opposing PE pair** (AP=`j-`, PA=`j`) and
  `assign_fieldmaps.py` must report **"all runs paired cleanly"**.

## 2. fMRIPrep report — `derivatives/fmriprep/sub-XX.html`
Open in a browser (scp it + its `figures/` dir locally, or X-forward firefox).
Flicker/overlay figures are the point — look, don't just read.

| Section | What "good" looks like | Failure signature |
|---|---|---|
| Anatomical: brain mask + segmentation | red outline hugs brain; GM/WM/CSF labels follow tissue | mask clips brain or grabs skull/dura; segmentation bleeds |
| Spatial normalization (T1↔MNI flicker) | gyri/ventricles line up across the flicker | gross shift/rotation, warped anatomy |
| FreeSurfer surfaces (white/pial) | pial hugs cortex, doesn't cross into skull/sinus | pial leaks into dura; white cuts through GM |
| **SDC** (fieldmap, before/after flicker) | "after" pulls frontal/temporal signal back onto the T1 | no change (fieldmap not applied) or over-warp |
| BOLD→T1 registration (BBR) | tissue boundaries align in the flicker | boundaries offset by mm |
| **Carpet / BOLD summary** | flat carpet; FD/DVARS traces low & stable | dark/bright bands, step changes, spikes tracking FD |
| Errors (bottom) | **"No errors to report!"** | any listed node crash |

**Motion (confounds):** sub-03 baseline **meanFD 0.11–0.15 mm, maxFD ≤0.85,
1–3% of volumes over 0.3**. Worry if a run has **meanFD > 0.3** or **>20–30%
censored** — the first-level GLM loses too many volumes.
```bash
# per-run mean FD + % censored (FD>0.3)
python3 - <<'PY'
import csv,statistics,glob
for f in sorted(glob.glob('derivatives/fmriprep/sub-05/ses-1/func/*confounds_timeseries.tsv')):
    r=list(csv.DictReader(open(f),delimiter='\t'))
    fd=[float(x['framewise_displacement']) for x in r if x['framewise_displacement'] not in ('','n/a')]
    c=sum(v>0.3 for v in fd); print(f.split('task-')[1][:12], f"meanFD={statistics.mean(fd):.3f} cens={100*c/len(r):.0f}%")
PY
```

## 3. tedana report — `derivatives/tedana/sub-XX/ses-Y/func/*_tedana_report.html`
Interactive kappa/rho scatter + component carpet + T2*/S0 maps.
- **T2\* and S0 maps** anatomically sensible (T2* shows GM/WM contrast).
- **Component classification** — sub-03 baseline: **47–57 total components,
  ~22–25 accepted / ~25–33 rejected per run** (≈ half accepted is healthy).
- kappa = BOLD-like (T2*-dependent, keep); rho = noise (S0-dependent, reject).
- Failure signatures: **0 / near-0 accepted** (denoised series gutted), very few
  total components, a large **`desc-decayFitFailures_mask`** or tiny
  **`desc-adaptiveGoodSignal_mask`** (T2* fit failed over much of the brain —
  usually bad/too-motiony data). `*_report.txt` has a plain-language summary.
```bash
# accepted vs rejected per run
python3 - <<'PY'
import csv,glob
from collections import Counter
for f in sorted(glob.glob('derivatives/tedana/sub-05/ses-1/func/*_desc-tedana_metrics.tsv')):
    r=list(csv.DictReader(open(f),delimiter='\t'))
    print(f.split('task-')[1][:12], dict(Counter(x['classification'] for x in r)))
PY
```

## 4. Warp to MNI (custom antsApplyTransforms step)
Final file: `*_space-MNI152NLin2009cAsym_desc-tedana_bold.nii.gz`.
- Overlay on the MNI template (`fsleyes`/`afni`) → brain fills the template,
  no shift/flip. A wrong transform pick shows as gross misalignment.

## 5. First-level AFNI (later, needs timing files)
- `*_xmat.jpg`: regressors sane, not collinear; motion regressors present.
- `_stats_REML` exists; the `rating_slope` (rating[1]) sub-brick is the target.
- `3dREMLfit.err` FDR "no voxels survived" warning is **benign**, not a crash.

---
### Watch-list specifically for sub-04 / sub-05
- Brand-new conversions with the newer `MR …TASKn` export naming → **verify task
  labels vs run order** in the report (sub-04 AFS, sub-05 SFA).
- First pass through FreeSurfer for both → check anatomical + surface panels.
- sub-05 is a new subject → eyeball T1 quality.
