#!/usr/bin/env python3
"""
==============================================================================
assign_fieldmaps.py -- PAIR EACH SPIN-ECHO FIELDMAP WITH THE BOLD RUN IT CORRECTS
==============================================================================

WHAT THIS SCRIPT DOES
  Echo-planar (EPI) BOLD images are stretched and squashed near air/tissue
  boundaries (for example the sinuses and ear canals). This is called
  susceptibility distortion. To undo it, a short "spin-echo fieldmap" is
  collected: two quick EPI images with opposite phase-encoding directions
  (AP = anterior-to-posterior, PA = posterior-to-anterior). The two images are
  distorted in opposite directions, so together they let fMRIPrep (through its
  SDCFlows module, "PEPOLAR"/TOPUP-style correction) estimate and remove the
  distortion. In this study one AP/PA pair (CMRR SEFM sequence) is acquired
  right BEFORE each task run.

  fMRIPrep only knows which fieldmap belongs to which run if the BIDS JSON
  "sidecar" files say so. This script works that out from the scan times and
  writes the answer into the JSON files:

    * B0FieldIdentifier  written into each fieldmap JSON   (e.g. "ses-2_sefm_run1")
    * B0FieldSource      written into every functional JSON of the matched run
                         (the same string, so the two are linked)
    * IntendedFor        written into each fieldmap JSON: the list of that run's
                         functional NIfTI files. This is the older way of
                         linking, kept so fMRIPrep versions before 20.2 also work.

  A "run" here means ALL functional files that share the same task-<label>:
  every echo of the multi-echo BOLD (this study has 4 echoes) and every SBRef
  (single-band reference) image. All of them get the same B0FieldSource.

  Each fieldmap is assigned to the FIRST run whose acquisition time is the same
  as or later than the fieldmap's own time. Times come from the DICOM
  AcquisitionTime stored in the JSON sidecars. If ANY file in the session lacks
  a usable AcquisitionTime, the whole session is ordered by SeriesNumber
  instead.

  The script then prints a QC report (one line per run, plus WARNING lines).
  It does not stop on a problem; it reports it. The last line is either
  "all runs paired cleanly" or "<N> issue(s) above -- review before running
  fMRIPrep".

WHO CALLS IT / WHERE IT SITS IN THE PIPELINE
  Called automatically by dcm2niix.sh, once per participant x session, right
  after the DICOMs of that session are converted to BIDS NIfTI + JSON:
      python3 "${script_dir}/assign_fieldmaps.py" "${func_dir}" "${fmap_dir}" "ses-${number_ses}"
  Before: dcm2niix.sh (DICOM -> BIDS conversion)
  After : fmriprep_tedana.sh (fMRIPrep reads B0FieldIdentifier/B0FieldSource)
  Re-run it after ANY dcm2niix re-conversion: dcm2niix writes a blanket
  IntendedFor that would otherwise replace the per-run pairing.

STEPS
  1. Read the three command-line arguments (func folder, fmap folder, session).
  2. Load every functional JSON (*_bold.json, *_sbref.json) and every fieldmap
     JSON (*_epi.json). Stop early if the session has no functional data
     (anatomy-only) or no fieldmaps.
  3. Choose how to put the scans in time order: AcquisitionTime if every file
     has one, otherwise SeriesNumber.
  4. Group the functional files into runs by their task-<label>, and sort the
     runs by time (a run's time = the earliest of its files).
  5. Give each fieldmap to the first run acquired at or after it. Fieldmaps
     acquired after the last run are left unassigned ("orphans").
  6. For each run: write B0FieldSource into its functional JSONs, write
     B0FieldIdentifier + IntendedFor into its fieldmap JSONs, and check it
     (has a fieldmap? has both AP and PA? fieldmap before the run? gap not
     too long?).
  7. Report orphan fieldmaps, runs without a fieldmap, and a final verdict.

INPUTS
  FUNC_DIR  : the session's BIDS func folder, e.g. .../sub-01/ses-2/func
              read: *_bold.json, *_sbref.json
  FMAP_DIR  : the session's BIDS fmap folder, e.g. .../sub-01/ses-2/fmap
              read: *_epi.json
  SES_LABEL : the session label used in the identifier and in IntendedFor
              paths, e.g. ses-2

OUTPUTS
  No new files. The JSON sidecars above are REWRITTEN IN PLACE (the whole file
  is re-saved with 2-space indentation):
    - each functional JSON gains/overwrites "B0FieldSource"
    - each assigned fieldmap JSON gains/overwrites "B0FieldIdentifier" and
      "IntendedFor"
  (Orphan fieldmaps are not changed.)
  A QC report is printed to the screen (lines start with "  [fieldmap]").

HOW TO RUN
  python3 assign_fieldmaps.py FUNC_DIR FMAP_DIR SES_LABEL
  e.g.
  python3 assign_fieldmaps.py BIDS/sub-01/ses-2/func BIDS/sub-01/ses-2/fmap ses-2
  Re-running is safe: the same inputs give the same tags.

METHODS SUMMARY (the exact rules used below)
  Fieldmap type   : spin-echo EPI pairs with opposite phase encoding (AP/PA),
                    BIDS suffix _epi, used by fMRIPrep/SDCFlows (PEPOLAR).
  Ordering basis  : DICOM AcquisitionTime (HHMMSS.ffffff, converted to seconds
                    after midnight); if missing in any file of the session,
                    SeriesNumber for the whole session.
  Run definition  : all *_bold and *_sbref files sharing one task-<label>
                    (all echoes); run time = earliest file in the run.
  Assignment rule : each fieldmap -> the first run whose time >= the
                    fieldmap's time.
  Linking fields  : B0FieldIdentifier = B0FieldSource = "<ses>_sefm_run<k>",
                    k = run's position in acquisition order (1, 2, 3, ...);
                    IntendedFor = "<ses>/func/<file>.nii.gz" for every file of
                    the run.
  QC warnings     : no fieldmap for a run; fieldmap(s) without two different
                    PhaseEncodingDirection values; fieldmap acquired after its
                    run; gap from the latest fieldmap to the run start > 600 s;
                    fieldmaps after the last run.
==============================================================================
"""
import glob
import json
import os
import re
import sys


# ==============================================================================
# SECTION 1. SETTINGS
# ------------------------------------------------------------------------------
# The only tunable value. Fieldmaps are acquired right before their run, so a
# gap of more than 10 minutes (600 s) between a fieldmap and "its" run is
# suspicious and produces a WARNING (it does not change the pairing).
# ==============================================================================
MAX_REASONABLE_GAP_S = 600


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

def load(path):
    """Read one JSON sidecar file and return its contents as a dictionary."""
    with open(path) as fh:
        return json.load(fh)


def save(path, data):
    """Write a dictionary back to a JSON sidecar file (overwrites it, 2-space indent)."""
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)


def acq_seconds(meta):
    """Return the scan's AcquisitionTime as seconds after midnight, or None.

    AcquisitionTime looks like "14:03:22.512500" (or "140322.5125"). The colons
    are removed, then hours (first 2 digits) x 3600 + minutes (next 2) x 60 +
    seconds (the rest, may include a decimal). If the field is missing or cannot
    be read as a time, None is returned (which makes the whole session fall back
    to SeriesNumber ordering, see STEP 3).
    """
    raw = meta.get("AcquisitionTime")
    if raw is None:
        return None
    digits = str(raw).replace(":", "")
    try:
        return int(digits[0:2]) * 3600 + int(digits[2:4]) * 60 + float(digits[4:] or 0)
    except (ValueError, IndexError):
        return None


def task_of(path):
    """Return the task label from a BIDS file name (the part after "_task-").

    e.g. "sub-01_ses-2_task-scripture_echo-1_bold.json" -> "scripture".
    Only letters and digits are kept. If the name has no task- part, the whole
    file name is returned (so that file becomes a "run" of its own).
    """
    m = re.search(r"_task-([A-Za-z0-9]+)", os.path.basename(path))
    return m.group(1) if m else os.path.basename(path)


# ==============================================================================
# SECTION 3. MAIN PROCEDURE
# ==============================================================================

def main():
    """Pair the fieldmaps of one session with its runs and tag the JSON files."""

    # --------------------------------------------------------------------------
    # STEP 1. READ THE COMMAND-LINE ARGUMENTS
    # Exactly three are required: FUNC_DIR FMAP_DIR SES_LABEL. Anything else
    # prints the usage line and exits with an error.
    # --------------------------------------------------------------------------
    if len(sys.argv) != 4:
        sys.exit("usage: assign_fieldmaps.py FUNC_DIR FMAP_DIR SES_LABEL")
    func_dir, fmap_dir, ses = sys.argv[1:4]

    # --------------------------------------------------------------------------
    # STEP 2. LOAD THE JSON SIDECARS
    # func  = every BOLD echo and every SBRef of the session.
    # fmaps = every spin-echo fieldmap image (AP and PA) of the session.
    # Each entry is [file path, contents as a dictionary].
    # An anatomy-only session (no func files) is normal: report it and stop.
    # Func files but no fieldmaps is a problem: warn and stop (nothing written).
    # --------------------------------------------------------------------------
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

    # --------------------------------------------------------------------------
    # STEP 3. CHOOSE HOW TO PUT SCANS IN TIME ORDER
    # The same yardstick must be used for every file, or the comparison in
    # STEP 5 would be meaningless. So: if EVERY func and fmap file has a
    # readable AcquisitionTime, sort by it (in seconds). If even one lacks it,
    # sort the whole session by SeriesNumber (the scanner's running count of
    # series; a missing SeriesNumber counts as 0).
    # keyfn(item) returns the sort value for one [path, contents] entry.
    # have_time records whether real times (and therefore gaps in seconds) exist.
    # --------------------------------------------------------------------------
    if all(acq_seconds(m) is not None for _, m in func + fmaps):
        keyfn = lambda item: acq_seconds(item[1])
        basis = "AcquisitionTime"
        have_time = True
    else:
        keyfn = lambda item: float(item[1].get("SeriesNumber", 0))
        basis = "SeriesNumber (AcquisitionTime missing)"
        have_time = False
    print("  [fieldmap] pairing by %s" % basis)

    # --------------------------------------------------------------------------
    # STEP 4. GROUP FUNCTIONAL FILES INTO RUNS AND SORT THE RUNS IN TIME
    # All files with the same task label (all echoes of the BOLD + the SBRefs)
    # form one run. Note: this assumes each task is acquired once per session;
    # two runs of the same task would be merged into one "run" here.
    # A run's time is the EARLIEST time among its files.
    # run_list = [(run time, task label, [files]), ...] sorted earliest first.
    # run_keys = just the run times, in the same order.
    # --------------------------------------------------------------------------
    runs = {}
    for item in func:
        runs.setdefault(task_of(item[0]), []).append(item)

    run_list = sorted(
        ((min(keyfn(f) for f in files), task, files) for task, files in runs.items()),
        key=lambda r: r[0],
    )
    run_keys = [r[0] for r in run_list]

    # --------------------------------------------------------------------------
    # STEP 5. ASSIGN EACH FIELDMAP TO A RUN
    # Rule: a fieldmap belongs to the FIRST run (in time order) whose time is
    # the same as or later than the fieldmap's time, i.e. the next run after it.
    # slots[i]  = the fieldmaps assigned to run i (0-based).
    # orphans   = fieldmaps acquired after the last run (no run follows them).
    # --------------------------------------------------------------------------
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

    # --------------------------------------------------------------------------
    # STEP 6. WRITE THE TAGS AND CHECK EACH RUN
    # For run number k (1, 2, 3 ... in time order) the shared tag is
    # "<SES_LABEL>_sefm_run<k>", e.g. "ses-2_sefm_run1".
    #   6a. Functional side: put B0FieldSource = tag in every JSON of the run,
    #       and collect the run's NIfTI paths for IntendedFor. Paths are written
    #       as "<ses>/func/<name>.nii.gz", i.e. relative to the participant
    #       folder (the BIDS IntendedFor convention).
    #   6b. Fieldmap side: put B0FieldIdentifier = tag and IntendedFor = that
    #       list in every fieldmap JSON assigned to the run.
    #   6c. QC: warn if the run has no fieldmap; if its fieldmaps do not include
    #       two different PhaseEncodingDirection values (PEPOLAR needs AP + PA);
    #       and, only when real times are available, if the gap between the
    #       LATEST assigned fieldmap and the run start is negative (fieldmap
    #       after its run) or longer than MAX_REASONABLE_GAP_S.
    #       Each warning adds 1 to the "problems" count.
    # The printed line starts with "!!" when the run has any warning.
    # --------------------------------------------------------------------------
    problems = 0
    for i, (run_time, task, files) in enumerate(run_list):
        ident = "%s_sefm_run%d" % (ses, i + 1)
        pair = slots[i]

        # 6a. Functional side: tag every echo (bold + sbref) of this run.
        niftis = []
        for fjf, fmeta in files:
            fmeta["B0FieldSource"] = ident
            save(fjf, fmeta)
            niftis.append(ses + "/func/" + os.path.basename(fjf).replace(".json", ".nii.gz"))

        # 6b. Fieldmap side: identify the group and (for old fMRIPrep) its targets.
        pedirs = []
        for gjf, gmeta in pair:
            gmeta["B0FieldIdentifier"] = ident
            gmeta["IntendedFor"] = niftis
            save(gjf, gmeta)
            pedirs.append(gmeta.get("PhaseEncodingDirection", "?"))

        # 6c. QC for this run.
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

    # --------------------------------------------------------------------------
    # STEP 7. SESSION-LEVEL CHECKS AND FINAL VERDICT
    #   - Orphan fieldmaps (after the last run): each one adds 1 problem. Their
    #     JSON files were not changed.
    #   - If fewer runs got a fieldmap than there are runs: add 1 problem.
    #     (A run with no fieldmap is therefore counted twice: once in STEP 6c
    #     and once here.)
    #   - Print "all runs paired cleanly" or the number of issues.
    # --------------------------------------------------------------------------
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
