#!/usr/bin/env Rscript
# =============================================================================
# dice_within_vs_between.R  --  is a task reliable at the INDIVIDUAL level?
#
# Companion to reproducibility_afni.sh, which writes dice_within_between.tsv:
#   within  rows: subject X ses-1 vs subject X ses-2
#   between rows: subjects X vs Y (X != Y), same session (ses-1 and ses-2 pairs)
#
# Follows Gorgolewski et al. (2013), "Single subject fMRI test-retest
# reliability metrics and confounding factors", NeuroImage 69:231-243:
#   * the between-subject Dice is averaged over ALL between pairs (both
#     sessions) for a task -> one baseline value B per task/contrast/variant/mask
#   * each subject's within-subject Dice w_i is compared with B
#   * significance: "a percentile bootstrap test of the Harrell-Davis (HD)
#     median" of the within-minus-between difference
#
# Implementation choices the paper does not pin down (stated so they can be
# cited or changed):
#   * the bootstrap resamples SUBJECTS (their d_i = w_i - B); B is held fixed,
#     as the paper treats the averaged between-subject overlap as the baseline
#   * the test is one-sided (H1: within > between), which is the paper's
#     hypothesis; a two-sided p is also reported
#   * NA Dice rows (missing maps) are dropped and counted
#   * fewer than MIN_N subjects -> no p-value (a bootstrap of 3 values is
#     meaningless); descriptives are still reported
#
# The HD median is computed with the Harrell & Davis (1982) weights directly
# (identical to Hmisc::hdquantile(x, 0.5)); if Hmisc is installed the two are
# cross-checked once at start-up.
#
# Usage : Rscript dice_within_vs_between.R <dice_within_between.tsv> [n_boot]
# Output: <same dir>/dice_within_vs_between_summary.tsv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript dice_within_vs_between.R <dice_within_between.tsv> [n_boot]")
infile <- args[1]
n_boot <- if (length(args) >= 2) as.integer(args[2]) else 5000L
MIN_N  <- 5L
set.seed(20131113)   # fixed: results reproducible run to run

# Harrell, F.E. & Davis, C.E. (1982). A new distribution-free quantile
# estimator. Biometrika 69:635-640. Weighted sum of order statistics with
# Beta((n+1)q, (n+1)(1-q)) weights.
hd_quantile <- function(x, q = 0.5) {
    x <- sort(x[is.finite(x)])
    n <- length(x)
    if (n == 0) return(NA_real_)
    if (n == 1) return(x)
    a <- (n + 1) * q
    b <- (n + 1) * (1 - q)
    w <- diff(pbeta((0:n) / n, a, b))
    sum(w * x)
}

if (requireNamespace("Hmisc", quietly = TRUE)) {
    z <- c(0.31, 0.12, 0.55, 0.47, 0.09, 0.73, 0.28)
    stopifnot(isTRUE(all.equal(hd_quantile(z), unname(Hmisc::hdquantile(z, 0.5)))))
}

d <- utils::read.delim(infile, stringsAsFactors = FALSE,
                       na.strings = c("NA", "nan", "NaN", ""))
need <- c("task", "contrast", "variant", "mask", "kind", "session", "subjA", "subjB", "dice")
miss <- setdiff(need, names(d))
if (length(miss)) stop("missing column(s): ", paste(miss, collapse = ", "))
d$dice <- suppressWarnings(as.numeric(d$dice))

keys <- unique(d[, c("task", "contrast", "variant", "mask")])
out <- vector("list", nrow(keys))

for (r in seq_len(nrow(keys))) {
    k <- keys[r, ]
    g <- d[d$task == k$task & d$contrast == k$contrast &
           d$variant == k$variant & d$mask == k$mask, ]
    wi <- g[g$kind == "within", ]
    bw <- g[g$kind == "between", ]
    n_na <- sum(is.na(g$dice))
    wi <- wi[!is.na(wi$dice), ]
    bw <- bw[!is.na(bw$dice), ]

    w <- wi$dice
    B <- if (nrow(bw)) mean(bw$dice) else NA_real_
    n <- length(w)
    dd <- w - B

    est <- ci_lo <- ci_hi <- p1 <- p2 <- NA_real_
    if (n >= 1 && is.finite(B)) est <- hd_quantile(dd)
    if (n >= MIN_N && is.finite(B)) {
        boot <- replicate(n_boot, hd_quantile(sample(dd, n, replace = TRUE)))
        ci    <- stats::quantile(boot, c(0.025, 0.975), names = FALSE)
        ci_lo <- ci[1]; ci_hi <- ci[2]
        # (count + 1) / (n_boot + 1): a bootstrap p is never exactly 0, only
        # "smaller than the resolution 1/n_boot allows".
        p1    <- (sum(boot <= 0) + 1) / (n_boot + 1)  # H1: within > between
        p2    <- min(1, 2 * (min(sum(boot <= 0), sum(boot >= 0)) + 1) / (n_boot + 1))
    }

    out[[r]] <- data.frame(
        task = k$task, contrast = k$contrast, variant = k$variant, mask = k$mask,
        n_subj = n, n_between_pairs = nrow(bw), n_na_dropped = n_na,
        within_mean = if (n) mean(w) else NA_real_,
        within_hd_median = if (n) hd_quantile(w) else NA_real_,
        between_mean = B,
        diff_hd_median = est, diff_ci95_lo = ci_lo, diff_ci95_hi = ci_hi,
        p_one_sided = p1, p_two_sided = p2,
        n_boot = if (is.finite(p1)) n_boot else NA_integer_,
        stringsAsFactors = FALSE)
}

res <- do.call(rbind, out)
res <- res[order(res$task, res$contrast, res$variant, res$mask), ]
outfile <- file.path(dirname(infile), "dice_within_vs_between_summary.tsv")
num <- vapply(res, is.double, logical(1))
res[num] <- lapply(res[num], function(v) round(v, 4))
utils::write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat(sprintf("within vs between Dice: %d group(s) -> %s\n", nrow(res), outfile))
if (any(res$n_subj < MIN_N))
    cat(sprintf("NOTE: groups with < %d subjects get descriptives only (no bootstrap p).\n", MIN_N))
print(res[, c("task", "contrast", "variant", "mask", "n_subj",
              "within_mean", "between_mean", "diff_hd_median", "p_one_sided")],
      row.names = FALSE)
