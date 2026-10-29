#!/usr/bin/env Rscript
# =============================================================================
# trr_slope.R  --  test-retest reliability of the RATING SLOPE, one task,
#                  every ROI (region level, Bayesian multilevel model)
#
# Companion to trr_slope_afni.sh, which calls it once per task on the ROI
# table that trr_export_afni.sh wrote (roi_task-<T>.tsv; one row per trial x
# ROI: Subj sess trial rating cond stimulus yield faith ROI Y).
#
# Question: is each participant's BOLD-vs-rating slope (how much the ROI's
# single-trial beta changes per step of the 1-4 "felt the Spirit" rating) the
# same in session 1 and session 2?
#
# Model, per ROI (brms; the structure of AFNI's TRR program with -cond, with
# the two-level condition replaced by the continuous rating):
#     Y ~ 0 + sess + sess:rc + (0 + sess | Subj) + (0 + sess:rc | Subj)
#   rc = rating - 2.5, so the session terms are the response at mid-scale
#   and the sess:rc terms are the slopes (beta change per rating step).
#   (0 + sess | Subj)      each participant's ses-1 and ses-2 intercepts,
#                          correlated  -> intercept TRR
#   (0 + sess:rc | Subj)   each participant's ses-1 and ses-2 slopes,
#                          correlated  -> SLOPE TRR (the main result)
#   Trial-to-trial noise is the residual, so the TRR is not attenuated by it
#   the way an ICC of first-level slope maps is (Chen et al. 2021,
#   NeuroImage 245:118647).
#   Sampler settings and priors as in AFNI's TRR.R: brms defaults (LKJ(1) on
#   the correlations), adapt_delta = 0.99, max_treedepth = 15, seed 1234.
#
# Two TRR summaries are reported for each correlation:
#   *_cor   the model's correlation parameter (as 3dLMEr -TRR reports it)
#           -> the primary estimate
#   *_afni  the correlation, across participants, of the sampled participant
#           effects, draw by draw (how AFNI's TRR.R computes its TRR) -> only
#           for comparison with trr_roi_summary.tsv
#   Each is summarised as the posterior mode (peak of the density, as TRR.R)
#   and the 95% highest density interval (coda::HPDinterval, as TRR.R).
#
# Participants: one needs at least MIN_DISTINCT ratings in EACH session (a
# slope needs at least two rating values); others are dropped and listed in
# <outdir>/trr_slope_task-<T>_dropped.txt.
#
# The Stan model is compiled once per task and re-used for every ROI.
#
# Usage : Rscript trr_slope.R <roi_task-T.tsv> <task> <outdir> [chains] [iterations]
# Output: <outdir>/trr_slope_summary_task-<T>.tsv   one row per ROI
#         <outdir>/rating_counts_task-<T>.tsv       trials per rating, per Subj x sess
#         <outdir>/trr_slope_task-<T>_dropped.txt   participants left out
#         <outdir>/TRR_slope_task-<T>_<ROI>.txt     model summary + TRR lines
#         <outdir>/TRR_slope_task-<T>_<ROI>.rds     the fitted brms model
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript trr_slope.R <roi_task-T.tsv> <task> <outdir> [chains] [iterations]")
infile  <- args[1]
task    <- args[2]
outdir  <- args[3]
CHAINS  <- if (length(args) >= 4) as.integer(args[4]) else 4L
ITER    <- if (length(args) >= 5) as.integer(args[5]) else 1000L
MIN_DISTINCT <- 2L     # distinct rating values needed in EACH session
CENTER       <- 2.5    # middle of the 1-4 rating scale

suppressPackageStartupMessages({ library(brms); library(coda) })
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- data --------------------------------------------------------------------
d <- read.delim(infile, stringsAsFactors = FALSE)
need <- c("Subj", "sess", "rating", "ROI", "Y")
if (length(miss <- setdiff(need, names(d)))) stop("missing column(s) in ", infile, ": ", paste(miss, collapse = " "))
d <- d[is.finite(d$Y) & d$rating %in% 1:4, need]
d$sess <- factor(d$sess, levels = c("ses1", "ses2"))
d <- d[!is.na(d$sess), ]
d$rc <- d$rating - CENTER

# Trials per rating, per participant x session (the same for every ROI, so
# taken from the first ROI).
one <- d[d$ROI == sort(unique(d$ROI))[1], ]
counts <- as.data.frame.matrix(table(paste(one$Subj, one$sess, sep = "\t"), factor(one$rating, levels = 1:4)))
counts <- data.frame(Subj = sub("\t.*", "", rownames(counts)), sess = sub(".*\t", "", rownames(counts)),
                     r1 = counts[, "1"], r2 = counts[, "2"], r3 = counts[, "3"], r4 = counts[, "4"])
counts$n_distinct <- rowSums(counts[, c("r1", "r2", "r3", "r4")] > 0)
write.table(counts, file.path(outdir, sprintf("rating_counts_task-%s.tsv", task)),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Drop participants without MIN_DISTINCT ratings in both sessions.
ok <- tapply(counts$n_distinct >= MIN_DISTINCT, counts$Subj, function(x) length(x) == 2 && all(x))
drop <- names(ok)[!ok]
writeLines(if (length(drop)) drop else character(0),
           file.path(outdir, sprintf("trr_slope_task-%s_dropped.txt", task)))
d <- d[!d$Subj %in% drop, ]
cat(sprintf("task-%s: %d participants kept, %d dropped (< %d distinct ratings in a session)%s\n",
            task, length(unique(d$Subj)), length(drop), MIN_DISTINCT,
            if (length(drop)) paste0(": ", paste(drop, collapse = " ")) else ""))
if (length(unique(d$Subj)) < 3) stop("fewer than 3 participants left for task-", task)

# --- helpers (mode and HDI as in AFNI's TRR.R) --------------------------------
post_mode <- function(x) { dd <- density(x); dd$x[which.max(dd$y)] }
summ <- function(x) { h <- HPDinterval(as.mcmc(x)); c(mode = post_mode(x), lo = h[1], hi = h[2]) }
# Draw-by-draw correlation, across participants, of two columns of a
# ranef(summary = FALSE) array (draws x participants x terms).
afni_cor <- function(re, a, b) vapply(seq_len(dim(re)[1]), function(i) cor(re[i, , a], re[i, , b]), numeric(1))

form <- bf(Y ~ 0 + sess + sess:rc + (0 + sess | Subj) + (0 + sess:rc | Subj))
fit0 <- NULL
rows <- list()

for (roi in sort(unique(d$ROI))) {
    dr <- d[d$ROI == roi, ]
    prefix <- file.path(outdir, sprintf("TRR_slope_task-%s_%s", task, roi))
    cat(sprintf("  %s %s: n=%d, %d trials\n", task, roi, length(unique(dr$Subj)), nrow(dr)))
    fit <- tryCatch({
        if (is.null(fit0)) {
            fit0 <- brm(form, data = dr, chains = CHAINS, cores = CHAINS, iter = ITER, seed = 1234, refresh = 0,
                        control = list(adapt_delta = 0.99, max_treedepth = 15))
            fit0
        } else {
            update(fit0, newdata = dr, chains = CHAINS, cores = CHAINS, iter = ITER, seed = 1234, refresh = 0,
                   control = list(adapt_delta = 0.99, max_treedepth = 15), recompile = FALSE)
        }
    }, error = function(e) { cat("  WARN: fit failed for ", roi, ": ", conditionMessage(e), "\n", sep = ""); NULL })
    if (is.null(fit)) {
        rows[[roi]] <- data.frame(task = task, ROI = roi, n_subj = length(unique(dr$Subj)), n_trials = nrow(dr))
        next
    }

    dr_post <- as_draws_df(fit)
    cor_icpt  <- dr_post[["cor_Subj__sessses1__sessses2"]]
    cor_slope <- dr_post[["cor_Subj__sessses1:rc__sessses2:rc"]]
    if (is.null(cor_icpt) || is.null(cor_slope))   # brms named them differently: say what it did name
        stop("correlation draws not found; cor_ draws are: ",
             paste(grep("^cor_", names(dr_post), value = TRUE), collapse = " "))
    re <- ranef(fit, summary = FALSE)$Subj
    afni_icpt  <- afni_cor(re, "sessses1", "sessses2")
    afni_slope <- afni_cor(re, "sessses1:rc", "sessses2:rc")
    fe <- fixef(fit, summary = FALSE)
    mean_slope <- (fe[, "sessses1:rc"] + fe[, "sessses2:rc"]) / 2

    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    max_rhat <- max(rhat(fit), na.rm = TRUE)

    s_cs <- summ(cor_slope); s_as <- summ(afni_slope)
    s_ci <- summ(cor_icpt);  s_ai <- summ(afni_icpt)
    ms <- c(mean(mean_slope), quantile(mean_slope, c(.025, .975)))
    rows[[roi]] <- data.frame(task = task, ROI = roi, n_subj = length(unique(dr$Subj)), n_trials = nrow(dr),
        slopeTRR_cor_mode = s_cs[1], slopeTRR_cor_lo = s_cs[2], slopeTRR_cor_hi = s_cs[3],
        slopeTRR_afni_mode = s_as[1], slopeTRR_afni_lo = s_as[2], slopeTRR_afni_hi = s_as[3],
        icptTRR_cor_mode = s_ci[1], icptTRR_cor_lo = s_ci[2], icptTRR_cor_hi = s_ci[3],
        icptTRR_afni_mode = s_ai[1], icptTRR_afni_lo = s_ai[2], icptTRR_afni_hi = s_ai[3],
        mean_slope = ms[1], mean_slope_lo = ms[2], mean_slope_hi = ms[3],
        max_rhat = max_rhat, n_divergent = n_div)

    saveRDS(fit, paste0(prefix, ".rds"))
    sink(paste0(prefix, ".txt"))
    print(summary(fit))
    cat(sprintf("\nSlope TRR (model correlation): mode %.3f; 95%% HDI [%.3f, %.3f]\n", s_cs[1], s_cs[2], s_cs[3]))
    cat(sprintf("Slope TRR (TRR.R style):       mode %.3f; 95%% HDI [%.3f, %.3f]\n", s_as[1], s_as[2], s_as[3]))
    cat(sprintf("Intercept TRR (model corr.):   mode %.3f; 95%% HDI [%.3f, %.3f]\n", s_ci[1], s_ci[2], s_ci[3]))
    cat(sprintf("Intercept TRR (TRR.R style):   mode %.3f; 95%% HDI [%.3f, %.3f]\n", s_ai[1], s_ai[2], s_ai[3]))
    cat(sprintf("Mean slope (avg of sessions):  %.4f [%.4f, %.4f]\n", ms[1], ms[2], ms[3]))
    cat(sprintf("Max Rhat %.3f; divergent transitions %d\n", max_rhat, n_div))
    sink()
    cat(sprintf("     slope TRR = %.2f [%.2f, %.2f]  (Rhat %.3f, %d divergent)\n",
                s_cs[1], s_cs[2], s_cs[3], max_rhat, n_div))
}

# A failed fit has only the first four columns; give it NA for the rest.
full <- names(rows[[which.max(vapply(rows, ncol, 1L))]])
out <- do.call(rbind, lapply(rows, function(r) { r[setdiff(full, names(r))] <- NA; r[full] }))
write.table(out, file.path(outdir, sprintf("trr_slope_summary_task-%s.tsv", task)),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("  summary: ", file.path(outdir, sprintf("trr_slope_summary_task-%s.tsv", task)), "\n", sep = "")
