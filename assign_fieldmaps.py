#!/usr/bin/env python3
"""
assign_fieldmaps.py -- pair spin-echo (PEPOLAR) fieldmaps with their BOLD runs.

Each CMRR SEFM AP/PA pair is acquired right before its task run. This script
matches every fieldmap to the run acquired immediately *after* it, using the
DICOM AcquisitionTime stored in the JSON sidecars (falling back to SeriesNumber
if times are unavailable).

The functional data are multi-echo (4 echoes) and may also include SBRef
images, so a "run" is the *group* of all func files sharing the same task-<label>
(every echo of the BOLD and its SBRef). Each fieldmap is matched to a run and
the script writes:

  * B0FieldIdentifier on each fieldmap JSON       (e.g. "ses-2_sefm_run1")
  * B0FieldSource     on every func JSON of the run (same string)
  * a per-run IntendedFor on each fieldmap (list of all that run's func NIfTIs;
    back-compat with fMRIPrep <20.2)

so fMRIPrep / SDCFlows applies the correct, run-specific distortion correction.

Sessions with no functional data (anatomy-only) are reported and skipped.

Usage: assign_fieldmaps.py FUNC_DIR FMAP_DIR SES_LABEL
"""
import glob
import json
import os
import re
import sys

# Warn if the time between a fieldmap and the run it is matched to exceeds this.
MAX_REASONABLE_GAP_S = 600


def load(path):
    with open(path) as fh:
        return json.load(fh)


def save(path, data):
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)


def acq_seconds(meta):
    """Seconds-since-midnight from AcquisitionTime, or None if unavailable."""
    raw = meta.get("AcquisitionTime")
    if raw is None:
        return None
    digits = str(raw).replace(":", "")
    try:
        return int(digits[0:2]) * 3600 + int(digits[2:4]) * 60 + float(digits[4:] or 0)
    except (ValueError, IndexError):
        return None


def task_of(path):
    """The task-<label> entity from a BIDS filename, or the filename itself."""
    m = re.search(r"_task-([A-Za-z0-9]+)", os.path.basename(path))
    return m.group(1) if m else os.path.basename(path)


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: assign_fieldmaps.py FUNC_DIR FMAP_DIR SES_LABEL")
    func_dir, fmap_dir, ses = sys.argv[1:4]

    func = [[jf, load(jf)] for jf in
            glob.glob(os.path.join(func_dir, "*_bold.json")) +
            glob.glob(os.path.join(func_dir, "*_sbref.json"))]
    fmaps = [[jf, load(jf)] for jf in glob.glob(os.path.join(fmap_dir, "*_epi.json"))]

    if not func:
        print("  [fieldmap] no functional runs in this session (anatomy-only); nothing to pair")
        return
    if not fmaps:
        print("  [fieldmap] WARNING: %d functional run file(s) but no *_epi.json "
              "fieldmaps in %s" % (len(func), fmap_dir))
        return

    # Choose an ordering basis valid for every file in the session.
    if all(acq_seconds(m) is not None for _, m in func + fmaps):
        keyfn = lambda item: acq_seconds(item[1])
        basis = "AcquisitionTime"
        have_time = True
    else:
        keyfn = lambda item: float(item[1].get("SeriesNumber", 0))
        basis = "SeriesNumber (AcquisitionTime missing)"
        have_time = False
    print("  [fieldmap] pairing by %s" % basis)

    # Group all func files (every echo of bold + sbref) into runs by task label.
    runs = {}
    for item in func:
        runs.setdefault(task_of(item[0]), []).append(item)

    # A run's time is the earliest of its files; order runs by acquisition.
    run_list = sorted(
        ((min(keyfn(f) for f in files), task, files) for task, files in runs.items()),
        key=lambda r: r[0],
    )
    run_keys = [r[0] for r in run_list]

    # Assign each fieldmap to the first run acquired at or after it.
    slots = {i: [] for i in range(len(run_list))}
    orphans = []
    for item in fmaps:
        k = keyfn(item)
        target = next((i for i, rk in enumerate(run_keys) if rk >= k), None)
        if target is None:
            orphans.append(item)
        else:
            slots[target].append(item)

    print("  [fieldmap] functional runs in acquisition order:")
    for i, (_, task, files) in enumerate(run_list):
        print("      run %d: task-%s (%d func files)" % (i + 1, task, len(files)))

    problems = 0
    for i, (run_time, task, files) in enumerate(run_list):
        ident = "%s_sefm_run%d" % (ses, i + 1)
        pair = slots[i]

        # Func side: tag every echo (bold + sbref) of this run.
        niftis = []
        for fjf, fmeta in files:
            fmeta["B0FieldSource"] = ident
            save(fjf, fmeta)
            niftis.append(ses + "/func/" + os.path.basename(fjf).replace(".json", ".nii.gz"))

        # Fieldmap side: identify the group and (for old fMRIPrep) its targets.
        pedirs = []
        for gjf, gmeta in pair:
            gmeta["B0FieldIdentifier"] = ident
            gmeta["IntendedFor"] = niftis
            save(gjf, gmeta)
            pedirs.append(gmeta.get("PhaseEncodingDirection", "?"))

        # --- QC for this run ---
        warns = []
        if not pair:
            warns.append("no fieldmap matched this run")
        elif len(set(pedirs)) < 2:
            warns.append("no opposing phase-encode pair (got %s); "
                         "PEPOLAR/TOPUP needs AP+PA" % ",".join(pedirs))
        gap_str = ""
        if pair and have_time:
            gap = run_time - max(keyfn(f) for f in pair)
            gap_str = "  gap=%.0fs" % gap
            if gap < 0:
                warns.append("fieldmap acquired AFTER its run (check ordering)")
            elif gap > MAX_REASONABLE_GAP_S:
                warns.append("unusually large gap (%.0fs) to run" % gap)

        marker = "  " if not warns else "!!"
        print("  [fieldmap] %s run %d (task-%s) -> %s  | %d fmap [%s], %d func img%s"
              % (marker, i + 1, task, ident, len(pair),
                 ",".join(pedirs) or "none", len(niftis), gap_str))
        for w in warns:
            print("  [fieldmap]        WARNING: %s" % w)
        problems += len(warns)

    if orphans:
        problems += len(orphans)
        print("  [fieldmap] WARNING: %d fieldmap(s) acquired after the last run "
              "and left unassigned:" % len(orphans))
        for jf, _ in orphans:
            print("        ", os.path.basename(jf))

    n_runs_with_fmap = sum(1 for i in slots if slots[i])
    if n_runs_with_fmap != len(run_list):
        problems += 1
        print("  [fieldmap] WARNING: %d run(s) but only %d had a fieldmap"
              % (len(run_list), n_runs_with_fmap))

    print("  [fieldmap] %s"
          % ("all runs paired cleanly" if problems == 0
             else "%d issue(s) above -- review before running fMRIPrep" % problems))


if __name__ == "__main__":
    main()
