#!/usr/bin/env Rscript
# Label-free operator gates of docs/plans/tcaclr_confirm_plan.md section 6, gate 3, on one frozen expression matrix.
# NO truth is read and nothing is scored: the CLR matrix and the operator matrices are computed only to evaluate
#   (i) the P1 implementation gate: tcaCLR at t = 1 (analysis/explore_stack/stack_lib.R::sk_a3r on the lazy-walk P^1) equals the
#       Stouffer combination of the within-candidate row z-scores of CLR, tolerance 1e-10 (docs/plans/fable_theory_memo.md 3.4 P1, 5.4 R6);
#   (ii) the reachability diagnostic at t = 9: median row SD of log P^9 over reachable off-diagonal entries >= 1 nat and
#       unreachable share of off-diagonal pairs <= 5% (otherwise the arm is "not evaluable");
#   plus the package-operator equivalence check and the P2/P3 descriptors (row-SD quantiles and CV, row skewness).
# Settings fixed by the plan: B-spline MI (bins "hg", no transform), row z clipped at 0, Stouffer; softmax over the top 50, keff 10;
# lazy walk 0.5 I + 0.5 A; t = 9. Score matrices are NOT saved; only the gate CSV and a log are written.
# Usage: Rscript analysis/tcaclr_confirm/label_free_gates.R --expr <tsv genes x profiles> --out <dir> [--threads 8] [--arm A|B]
args <- commandArgs(trailingOnly = TRUE)
opt <- list(expr = NULL, out = NULL, threads = 8L, arm = "?"); i <- 1L
while (i <= length(args)) { k <- sub("^--", "", args[i]); stopifnot(k %in% names(opt)); v <- args[i + 1L]
  opt[[k]] <- if (k == "threads") as.integer(v) else v; i <- i + 2L }
stopifnot(!is.null(opt$expr), !is.null(opt$out)); dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "label_free_gates.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
vmhwm <- function() { l <- readLines("/proc/self/status"); v <- sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", grep("^VmHWM:", l, value = TRUE)); as.numeric(v) / 1048576 }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/robust_z/robust_z_lib.R"); source("analysis/explore_stack/stack_lib.R")
t0 <- Sys.time()
ex <- utils::read.delim(opt$expr, check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; rm(ex); storage.mode(X) <- "double"
G <- nrow(X); say(sprintf("arm %s: %d genes x %d profiles, non-finite %d, duplicated names %d", opt$arm, G, ncol(X), sum(!is.finite(X)), sum(duplicated(rownames(X)))))
stopifnot(all(is.finite(X)))
fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
fit$calibrate(method = "normal", combine = "stouffer"); C <- fit$clr_scores
fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator; fit$release(); rm(X, fit); invisible(gc(verbose = FALSE))
A <- rz_operator(C); eq_op <- max(abs(A - A_pkg)); rm(A_pkg)
say(sprintf("MI + CLR + operator done in %.0f s; max|rz_operator(CLR) - package operator| = %.3g; VmHWM %.1f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), eq_op, vmhwm()))
P <- rz_lazy_P(A)
## ---- P1 implementation gate (t = 1) --------------------------------------------------------------------------------------
P1 <- as.matrix(P); S1 <- sk_a3r(P1); rm(P1)
cand <- A > 0; diag(cand) <- FALSE; ncand <- rowSums(cand)
Zi <- matrix(0, G, G)
for (i in which(ncand >= 2L)) { j <- which(cand[i, ]); v <- C[i, j]; sdv <- stats::sd(v)
  z <- if (is.finite(sdv) && sdv > 0) (v - mean(v)) / sdv else rep(NA_real_, length(j)); z[is.na(z) | z < 0] <- 0; Zi[i, j] <- z }
Sind <- (Zi + t(Zi)) / sqrt(2); diag(Sind) <- 0; rm(Zi)
big <- ncand > 10L; keep <- outer(big, big, "&")
d_all <- max(abs(S1 - Sind)); d_big <- max(abs((S1 - Sind)[keep]))
say(sprintf("P1: rows %d, rows with <= 10 candidates %d (min %d, median %g); max|diff| all entries %.3g; rows > 10 candidates %.3g", G, sum(!big), min(ncand), stats::median(ncand), d_all, d_big))
rm(S1, Sind, keep, cand); invisible(gc(verbose = FALSE))
## ---- reachability diagnostic (t = 9) --------------------------------------------------------------------------------------
Pt <- as.matrix(P); for (s in 2:9) { Pt <- as.matrix(Pt %*% P); if (s %% 4L == 0L) invisible(gc(verbose = FALSE)) }
rch <- sk_reach(Pt)
L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA; rm(Pt); invisible(gc(verbose = FALSE))
n <- rowSums(!is.na(L)); mu <- rowMeans(L, na.rm = TRUE); D <- L - mu
sdr <- sqrt(rowSums(D^2, na.rm = TRUE) / (n - 1)); sdr[n < 2] <- NA
skew <- rowSums(D^3, na.rm = TRUE) / n / (rowSums(D^2, na.rm = TRUE) / n)^1.5
Pl <- exp(L); Pl <- Pl - rowMeans(Pl, na.rm = TRUE); skew_lin <- rowSums(Pl^3, na.rm = TRUE) / n / (rowSums(Pl^2, na.rm = TRUE) / n)^1.5
rm(L, D, Pl); invisible(gc(verbose = FALSE))
med_sd <- stats::median(sdr, na.rm = TRUE); unreach <- unname(rch["unreach_share"])
say(sprintf("reachability: unreachable share %.5f; rows with < 2 reachable entries %d; median row SD of log P^9 %.3f nats (IQR %.3f-%.3f; CV %.3f); median row skewness log %.3f, linear %.3f; VmHWM %.1f GB",
  unreach, rch[["rows_lt2_reachable"]], med_sd, stats::quantile(sdr, 0.25, na.rm = TRUE), stats::quantile(sdr, 0.75, na.rm = TRUE),
  stats::sd(sdr, na.rm = TRUE) / mean(sdr, na.rm = TRUE), stats::median(skew, na.rm = TRUE), stats::median(skew_lin, na.rm = TRUE), vmhwm()))
res <- data.frame(arm = opt$arm, genes = G, check = c("operator_equals_package_max_abs", "P1_t1_identity_max_abs_rows_gt10_candidates", "P1_t1_identity_max_abs_all_rows",
    "rows_with_le10_candidates", "median_candidates_per_row", "reach_unreachable_share", "reach_rows_lt2_reachable", "median_row_sd_logP9_nats",
    "row_sd_q25", "row_sd_q75", "row_sd_cv", "median_row_skew_logP9", "median_row_skew_linearP9"),
  value = c(eq_op, d_big, d_all, sum(!big), stats::median(ncand), unreach, rch[["rows_lt2_reachable"]], med_sd, stats::quantile(sdr, 0.25, na.rm = TRUE),
    stats::quantile(sdr, 0.75, na.rm = TRUE), stats::sd(sdr, na.rm = TRUE) / mean(sdr, na.rm = TRUE), stats::median(skew, na.rm = TRUE), stats::median(skew_lin, na.rm = TRUE)),
  threshold = c("1e-12", "1e-10", "report", "report", "report", "<= 0.05", "report", ">= 1", "report", "report", "report", "report", "report"), stringsAsFactors = FALSE)
res$pass <- c(eq_op < 1e-12, d_big <= 1e-10, NA, NA, NA, unreach <= 0.05, NA, med_sd >= 1, NA, NA, NA, NA, NA)
res$peak_gb <- vmhwm(); res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
utils::write.csv(res, file.path(opt$out, "label_free_gates.csv"), row.names = FALSE); print(res, row.names = FALSE)
gate <- all(res$pass, na.rm = TRUE); say(if (gate) "GATES PASSED (P1 implementation, reachability)" else "A GATE FAILED: the arm is not evaluable until explained")
say(sprintf("peak memory (VmHWM) %.1f GB; wall %.0f s", vmhwm(), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
