#!/usr/bin/env Rscript
# Reachability-confinement check on the real K562 matrix (docs/plans/fable_sim_reading.md section 2; PI approval 17 Nov 2026).
# STRUCTURE ONLY: for member rows of each CORUM complex (the scoreboard truth: corum_sets on the truth-gene K562 matrix) report the
# number of reachable off-diagonal entries of log P^9 and the share of them inside the complex. No method is scored against any truth;
# no AUPR, enrichment or per-pair score is computed. Spent data (K562 truth-gene matrix, CORUM). Usage (repo root): Rscript analysis/sim_planted/real_reach_k562.R
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/robust_z/robust_z_lib.R"); source("analysis/human_truths.R")
step_walk <- function(Pt, P) as.matrix(Pt %*% P)
X <- as.matrix(read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", row.names = 1, check.names = FALSE)); G <- nrow(X)
cat(sprintf("matrix %d genes x %d profiles\n", G, ncol(X)))
fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = 4L)$calibrate(method = "normal", combine = "stouffer")
fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A <- fit$operator; fit$release()
P <- rz_lazy_P(A, 0.5); Pt <- as.matrix(P); for (t in 2:9) Pt <- step_walk(Pt, P)
L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA; rm(Pt)
reach <- !is.na(L); nreach <- rowSums(reach); rm(L)
cat(sprintf("unreachable share of all pairs %.4f; rows with fewer than 2 reachable: %d; median reachable per row %d of %d\n",
            1 - sum(nreach) / (G * (G - 1)), sum(nreach < 2), as.integer(median(nreach)), G - 1))
cs <- corum_sets(rownames(X)); nm <- names(cs$members)
rows <- lapply(seq_along(cs$members), function(k) {
  m <- cs$members[[k]]; n <- length(m)
  inm <- rowSums(reach[m, m, drop = FALSE])
  rp <- mean(grepl("^RP[LS][0-9]+[A-Z]?$", rownames(X)[m]))
  data.frame(complex = nm[k], size = n, ribosomal = rp >= 0.3, reachable_offdiag = mean(nreach[m]), in_module_reachable = mean(inm),
             module_size_minus1 = n - 1, inmodule_share = mean(inm / nreach[m], na.rm = TRUE))
})
d <- do.call(rbind, rows)
d$class <- ifelse(d$ribosomal, "ribosomal", ifelse(d$size <= 6, "3-6", ifelse(d$size <= 19, "7-19", "20+")))
s <- do.call(rbind, lapply(split(d, d$class), function(x) data.frame(class = x$class[1], complexes = nrow(x), median_size = median(x$size),
  median_reachable = median(x$reachable_offdiag), median_inmodule_share = median(x$inmodule_share), min_inmodule_share = min(x$inmodule_share),
  complexes_closed = sum(x$reachable_offdiag < 1.5 * x$size), complexes_inmodule_over_half = sum(x$inmodule_share > 0.5))))
print(format(s, digits = 3), row.names = FALSE)
dir.create("results/sim_planted/addendum", showWarnings = FALSE, recursive = TRUE)
write.csv(d, "results/sim_planted/addendum/real_reach_k562_by_complex.csv", row.names = FALSE)
write.csv(s, "results/sim_planted/addendum/real_reach_k562_by_class.csv", row.names = FALSE)
