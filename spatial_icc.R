#!/usr/bin/env Rscript
# ==============================================================================
# spatial_icc.R -- SPATIAL ICC(2,1) BETWEEN TWO SESSIONS' BETA MAPS (ONE SUBJECT)
# ==============================================================================
#
# WHAT THIS SCRIPT DOES
#   Answers: "for one participant, task and contrast, how well does the beta
#   map from session 1 agree with the beta map from session 2?"
#   It reads a two-column list of numbers (one row per brain voxel: the
#   session-1 beta and the session-2 beta of that voxel) and prints ONE number,
#   the intraclass correlation coefficient (ICC) between the two columns.
#
#   What the ICC means here
#     An ICC asks how much of the total variation in a set of measurements is
#     due to real differences between the things being measured ("targets"),
#     as opposed to disagreement between repeated measurements of the same
#     target. Here the targets are the VOXELS and the repeated measurements are
#     the two SESSIONS. ICC = 1 means both sessions give every voxel the same
#     value; ICC near 0 means the sessions do not agree.
#
#     ICC(2,1), "absolute agreement" (Shrout & Fleiss 1979 notation; called
#     ICC(A,1) by McGraw & Wong 1996):
#       "2"     = two-way model: both the targets (voxels) and the measurement
#                 occasions (sessions) are treated as sources of variance.
#       "1"     = reliability of a SINGLE measurement (one session), not of the
#                 average of the two.
#       "absolute agreement" = a systematic difference between sessions (e.g.
#                 every voxel 0.2 higher in session 2) counts as disagreement
#                 and lowers the ICC. (A "consistency" ICC, like ICC(3,1),
#                 would ignore such a constant shift.)
#     This is why it is reported next to the Pearson spatial correlation
#     (spatial_r): Pearson r only asks whether the two maps rise and fall
#     together, so a map that is uniformly shifted or scaled in session 2 can
#     still have r = 1. ICC(2,1) is lowered by a shift or a scale change.
#
#   The value comes from the established `irr` R package, not a hand-written
#   formula, so the exact ICC form is named and citable:
#     irr::icc(ratings, model = "twoway", type = "agreement", unit = "single")
#   Only the ICC value is used. (irr::icc also returns an F-test and a
#   confidence interval; they are not used. Neighbouring voxels are not
#   independent, so those would not be valid anyway.)
#
# WHO CALLS IT / WHERE IT SITS IN THE PIPELINE
#   Called by reproducibility_afni.sh (its spatial_icc function), once per
#   participant x task x contrast x mask. reproducibility_afni.sh:
#     1. copies this file to <REPRO_OUT>/spatial_icc.R (a folder the AFNI
#        container can see),
#     2. dumps the paired in-mask betas with
#          3dmaskdump -noijk -mask <mask> -o <pairs file> <ses-1 Coef> <ses-2 Coef>
#     3. runs, inside the AFNI 26.3.04 Apptainer image with --cleanenv
#        (so the image's own R and irr package are used):
#          Rscript <REPRO_OUT>/spatial_icc.R <pairs file>
#     4. keeps the last line printed (the ICC or NA) and records it as the
#        metric "spatial_icc" in reproducibility_summary.tsv.
#   Before: first_level_afni.sh (makes the beta maps)
#   After : nothing; the number goes straight into the summary table.
#
# STEPS
#   1. Check the command line and that the input file exists.
#   2. Load the irr package.
#   3. Read the file and keep its LAST TWO columns (session 1, session 2).
#   4. Drop voxels where either value is not a finite number; check that at
#      least 3 voxels remain and that neither session's values are all equal.
#   5. Compute ICC(2,1), absolute agreement, with irr::icc.
#   6. Print the value to 6 decimal places.
#   If any step fails, the script prints "NA" instead, writes the reason to
#   the error stream (stderr), and still exits with status 0 (success), so one
#   bad map never stops the batch job.
#
# INPUTS
#   <paired_values.txt> : plain text, whitespace-separated, no header, at least
#     two numeric columns. The LAST TWO columns are used as
#     (session-1 beta, session-2 beta). This fits both
#       3dmaskdump -noijk mask s1 s2   (2 columns: s1 s2)        <- what the
#                                                                  pipeline uses
#       3dmaskdump mask s1 s2          (5 columns: i j k s1 s2)
#
# OUTPUTS
#   One line on standard output: the ICC rounded to 6 decimals (e.g. 0.734512),
#   or "NA". No files are written.
#
# HOW TO RUN
#   Rscript spatial_icc.R <paired_values.txt>
#   (Requires the R package irr; the AFNI 26.3.04 image contains it.)
#
# METHODS SUMMARY (the exact values used below)
#   Measure    : intraclass correlation, ICC(2,1) absolute agreement
#                (Shrout & Fleiss 1979; = ICC(A,1), McGraw & Wong 1996)
#   Software   : R package irr, irr::icc(m, model = "twoway",
#                type = "agreement", unit = "single"); defaults r0 = 0,
#                conf.level = 0.95 (CI/F-test not used)
#   Targets    : voxels inside the mask given to 3dmaskdump (one row each)
#   Raters     : the two sessions (ses-1, ses-2) -> 2 columns
#   Values     : first-level GLM beta (Coef) of the contrast, per voxel
#   Exclusions : rows with a non-finite value in either session are dropped
#   Returns NA : no/missing input file; irr not installed; file unreadable or
#                < 2 columns or < 3 rows; < 3 finite voxel pairs; a session
#                with zero standard deviation; irr::icc error or non-finite
#                result
# ==============================================================================


# ==============================================================================
# SECTION 1. SETTINGS
# ------------------------------------------------------------------------------
# None. The only input is the file name on the command line; every analysis
# choice is written directly in the irr::icc call in SECTION 4.
# ==============================================================================


# ==============================================================================
# SECTION 2. HELPER FUNCTIONS
# ==============================================================================

# fail_NA(msg): give up gracefully.
#   Writes "spatial_icc.R: <msg>" to stderr (only if msg is not empty), prints
#   "NA" on standard output (the caller reads this as the result) and stops the
#   script with exit status 0, so the calling batch job keeps going.
fail_NA <- function(msg) {
    if (nzchar(msg)) message("spatial_icc.R: ", msg)
    cat("NA\n")
    quit(save = "no", status = 0)   # status 0 so the batch never aborts on this
}


# ==============================================================================
# SECTION 3. INPUT CHECKS AND READING THE DATA
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1. COMMAND LINE: one argument, the paired-values file, which must exist.
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) fail_NA("no input file given")
f <- args[1]
if (!file.exists(f)) fail_NA(sprintf("input file not found: %s", f))

# ------------------------------------------------------------------------------
# STEP 2. LOAD THE irr PACKAGE (quietly). If it is not installed, return NA.
# (The message mentions R_LIBS_USER, but R looks in all of its library
# folders; in the AFNI image irr is installed inside the image itself.)
# ------------------------------------------------------------------------------
ok <- suppressWarnings(suppressMessages(require(irr, quietly = TRUE)))
if (!ok) fail_NA("R package 'irr' not available in R_LIBS_USER")

# ------------------------------------------------------------------------------
# STEP 3. READ THE FILE
# Whitespace-separated numbers, no header line. comment.char = "" means no
# character is treated as a comment marker; fill = TRUE pads short rows with NA
# instead of failing. Any read error -> NULL -> NA. At least 2 columns and 3
# rows are required.
# ------------------------------------------------------------------------------
d <- tryCatch(
    utils::read.table(f, header = FALSE, colClasses = "numeric",
                      comment.char = "", fill = TRUE),
    error = function(e) NULL)
if (is.null(d) || ncol(d) < 2 || nrow(d) < 3)
    fail_NA("could not read >=2 columns / >=3 rows")


# ==============================================================================
# SECTION 4. COMPUTE THE SPATIAL ICC
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 4. BUILD THE VOXEL x SESSION MATRIX AND CHECK IT
# m = the last two columns: column 1 = session 1, column 2 = session 2; one row
# per voxel (this is the "n targets x m raters" layout irr::icc expects).
# Rows where either value is NA, NaN or infinite are removed.
# At least 3 voxel pairs must remain, and each session must vary across voxels
# (a standard deviation of 0 makes the ICC undefined).
# ------------------------------------------------------------------------------
# Last two columns = the two sessions.
m <- as.matrix(d[, (ncol(d) - 1):ncol(d)])
# Keep only finite paired rows.
m <- m[is.finite(m[, 1]) & is.finite(m[, 2]), , drop = FALSE]
if (nrow(m) < 3) fail_NA("fewer than 3 finite paired voxels")
if (stats::sd(m[, 1]) == 0 || stats::sd(m[, 2]) == 0)
    fail_NA("zero variance in a session (degenerate ICC)")

# ------------------------------------------------------------------------------
# STEP 5. ICC(2,1), ABSOLUTE AGREEMENT
#   model = "twoway"    : voxels and sessions both modelled (the "2" in ICC(2,1))
#   type  = "agreement" : a constant shift between sessions lowers the ICC
#   unit  = "single"    : reliability of one session (the ",1")
# If irr::icc errors, or its value is missing or not finite, return NA.
# ------------------------------------------------------------------------------
res <- tryCatch(
    irr::icc(m, model = "twoway", type = "agreement", unit = "single"),
    error = function(e) NULL)
if (is.null(res) || is.null(res$value) || !is.finite(res$value))
    fail_NA("irr::icc returned no finite value")

# ------------------------------------------------------------------------------
# STEP 6. PRINT THE RESULT (6 decimal places) ON STANDARD OUTPUT.
# ------------------------------------------------------------------------------
cat(sprintf("%.6f\n", res$value))
