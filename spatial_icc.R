#!/usr/bin/env Rscript
# =============================================================================
# spatial_icc.R  --  per-subject spatial ICC(2,1), absolute agreement
#
# Companion to reproducibility_afni.sh. Computes a WITHIN-SUBJECT reproducibility
# number for one subject/task/contrast: the intraclass correlation between the
# two sessions' voxelwise beta values, treating VOXELS as the targets and the
# two SESSIONS as the repeated measurements.
#
# This is the ICC analogue of the Pearson spatial correlation: same paired-voxel
# scatter, but ICC(2,1) absolute agreement additionally penalizes a systematic
# between-session offset/scaling (the "slope != 1" blind spot of Pearson r).
#
# Uses the established `irr` package rather than a hand-rolled formula, so the
# exact Shrout & Fleiss form is named explicitly and is reproducible/citable:
#   irr::icc(model="twoway", type="agreement", unit="single")  ==  ICC(2,1), A,1
#
# Input : a whitespace-delimited file with (at least) two numeric columns; the
#         LAST TWO columns are taken as (session-1 beta, session-2 beta). This
#         matches `3dmaskdump -noijk mask s1 s2` (2 cols) and also a full
#         `3dmaskdump` dump (i j k s1 s2 -> last two used).
# Output: a single line with the ICC value to 6 dp, or "NA" if it cannot be
#         computed (too few/again-degenerate voxels, zero variance, load error).
#
# Usage : Rscript spatial_icc.R <paired_values.txt>
# =============================================================================

fail_NA <- function(msg) {
    if (nzchar(msg)) message("spatial_icc.R: ", msg)
    cat("NA\n")
    quit(save = "no", status = 0)   # status 0 so the batch never aborts on this
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) fail_NA("no input file given")
f <- args[1]
if (!file.exists(f)) fail_NA(sprintf("input file not found: %s", f))

ok <- suppressWarnings(suppressMessages(require(irr, quietly = TRUE)))
if (!ok) fail_NA("R package 'irr' not available in R_LIBS_USER")

d <- tryCatch(
    utils::read.table(f, header = FALSE, colClasses = "numeric",
                      comment.char = "", fill = TRUE),
    error = function(e) NULL)
if (is.null(d) || ncol(d) < 2 || nrow(d) < 3)
    fail_NA("could not read >=2 columns / >=3 rows")

# Last two columns = the two sessions.
m <- as.matrix(d[, (ncol(d) - 1):ncol(d)])
# Keep only finite paired rows.
m <- m[is.finite(m[, 1]) & is.finite(m[, 2]), , drop = FALSE]
if (nrow(m) < 3) fail_NA("fewer than 3 finite paired voxels")
if (stats::sd(m[, 1]) == 0 || stats::sd(m[, 2]) == 0)
    fail_NA("zero variance in a session (degenerate ICC)")

res <- tryCatch(
    irr::icc(m, model = "twoway", type = "agreement", unit = "single"),
    error = function(e) NULL)
if (is.null(res) || is.null(res$value) || !is.finite(res$value))
    fail_NA("irr::icc returned no finite value")

cat(sprintf("%.6f\n", res$value))
