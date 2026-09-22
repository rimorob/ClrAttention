# Design-A diffusion depth sweep (FIXED MI / FIXED operator, vary t only).
#
# This script REPLACES the earlier version, which re-estimated MI at every
# depth (rejected: it confounds depth with graph re-estimation and risks
# rediscovering diffusion-injected similarity).
#
# Design A: estimate MI once from X, build one row-stochastic operator Ahat,
# hold MI/Ahat/P fixed, and sweep only the diffusion depth t in
# E(t) = P^t X with P = (1-alpha) I + alpha Ahat.
#
# Per depth t:
#   * fresh module finder on E(t): complete-linkage clustering on 1-|cor|,
#     cut at fixed height (PROVISIONAL procedure -- the exact per-depth
#     similarity graph and thresholding are not yet settled).
#   * hard-module recovery vs planted modules (Hungarian + independent-best).
#   * diffusion affiliation profiles m_ik(t) = sum_{j in C_k} (P^t)_{ij}
#     against cores C_k = t=0 connected components (size >= 2), fixed.
#
# Synthetic: 26 genes x 120 samples. Planted A(8)/B(7)/C(6) + 2 singleton
# bridge genes (g22: A-leaning with C loading; g23: B-leaning with C loading)
# + 3 noise genes.
# Portable setup: run from a repo checkout (loads the package in place via
# pkgload) or with the package installed (library(clr)). Outputs go next to
# this script, so the demo works on any machine.
demo_dir <- (function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
repo_root <- normalizePath(file.path(demo_dir, ".."), mustWork = FALSE)
if (requireNamespace("clr", quietly = TRUE)) {
  library(clr)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE)
} else {
  stop("Need either the installed 'clr' package or pkgload to load it from this repo checkout.")
}
outdir <- file.path(demo_dir, "depth_sweep")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(7)

## ---- 1. synthetic data with combinatorial genes ----
planted <- list(
  module_A = paste0("g", 1:8),
  module_B = paste0("g", 9:15),
  module_C = paste0("g", 16:21)
)
combo_A <- "g22"   # singleton bridge: leans A, partial C loading
combo_B <- "g23"   # singleton bridge: leans B, partial C loading
noise_genes <- paste0("g", 24:26)
genes <- c(unlist(planted, use.names = FALSE), combo_A, combo_B, noise_genes)
n_samples <- 120L
X <- matrix(rnorm(length(genes) * n_samples), nrow = length(genes),
            dimnames = list(genes, paste0("s", seq_len(n_samples))))
factors <- list(A = rnorm(n_samples), B = rnorm(n_samples), C = rnorm(n_samples))
for (g in planted$module_A) X[g, ] <- 0.85 * factors$A + sqrt(1 - 0.85^2) * X[g, ]
for (g in planted$module_B) X[g, ] <- 0.85 * factors$B + sqrt(1 - 0.85^2) * X[g, ]
for (g in planted$module_C) X[g, ] <- 0.85 * factors$C + sqrt(1 - 0.85^2) * X[g, ]
mk_bridge <- function(g, f1, f2, a1 = 0.62, a2 = 0.62) {
  resid_sd <- sqrt(max(1e-6, 1 - a1^2 - a2^2))
  X[g, ] <<- a1 * factors[[f1]] + a2 * factors[[f2]] + resid_sd * X[g, ]
}
mk_bridge(combo_A, "A", "C")
mk_bridge(combo_B, "B", "C")

## ---- 2. CLR pipeline ONCE (Design A: fixed MI, fixed operator) ----
att <- ClrAttention$new(X)
att$estimate_mi()$calibrate()$select_threshold(B = 100)
cat("HC threshold:", att$threshold, "\n")
# TOY-ONLY conductance knob. The HC threshold works as designed (kills weak
# cross-talk), but then this toy's operator is block-diagonal and diffusion
# depth does nothing. Real data has natural cross-talk; the toy needs help.
# We rebuild the operator once at a lower tau so the graph has modest
# conductance. This choice is a property of the toy, not of the method.
tau_demo <- 1.2
att$build_operator(tau = tau_demo)
Ahat <- att$operator
alpha <- att$params$operator$alpha
cat("demo tau:", tau_demo, " mean degree:", mean(rowSums(Ahat > 0)),
    " alpha:", alpha, "\n")
G <- nrow(Ahat)
P <- (1 - alpha) * diag(G) + alpha * Ahat
rownames(P) <- colnames(P) <- rownames(Ahat)

## ---- 3. cores: t=0 hard modules (size >= 2), fixed across t ----
# PROVISIONAL per-depth module finder (not settled): complete-linkage
# hierarchical clustering on 1 - |cor|, cut at fixed height h_c.
# Complete (not single) linkage resists chaining through bridge genes:
# connected components on this operator merge A+C through the bridges
# into a 17-gene mega-core, while complete linkage keeps them separate.
h_c <- 0.6
find_modules_cor <- function(Ct) {
  d <- stats::as.dist(1 - abs(Ct))
  cl <- stats::cutree(stats::hclust(d, method = "complete"), h = h_c)
  found <- split(rownames(Ct), cl)
  names(found) <- paste0("module_", seq_along(found))
  found[order(-lengths(found))]
}
t0_modules <- find_modules_cor(cor(t(X)))
cores <- t0_modules[lengths(t0_modules) >= 2L]
cat("cores:", length(cores),
    paste(names(cores), sprintf("(n=%d)", lengths(cores)), collapse = " "), "\n")
K <- length(cores)
Cmat <- matrix(0, G, K, dimnames = list(rownames(Ahat), names(cores)))
for (k in seq_len(K)) Cmat[cores[[k]], k] <- 1

## ---- 4. depth sweep ----
T_steps <- 30L
Pt <- diag(G); rownames(Pt) <- colnames(Pt) <- rownames(Ahat)
Et <- X

res <- data.frame(
  t = integer(), n_modules = integer(), largest = integer(),
  hungarian_wf1 = numeric(), best_mean_f1 = numeric(),
  combo_entropy = numeric(), unexplained = numeric(),
  within_cor = numeric(), between_cor = numeric()
)
combo_genes <- c(combo_A, combo_B)
affil_combo <- array(NA_real_, dim = c(length(combo_genes), K, T_steps + 1L),
                     dimnames = list(combo_genes, names(cores), NULL))
mt <- ModuleTester$new(planted)
inA <- planted$module_A; inB <- planted$module_B; inC <- planted$module_C

for (t in 0:T_steps) {
  if (t > 0) { Pt <- Pt %*% P; Et <- P %*% Et }

  # fresh module finder on E(t): complete-linkage |cor| clustering
  Ct <- cor(t(Et))
  found <- find_modules_cor(Ct)
  mt$set_found(found)
  sh <- mt$summary(method = "hungarian")
  sb <- mt$summary(method = "best")

  # affiliation profiles against fixed cores
  M <- Pt %*% Cmat
  rownames(M) <- rownames(Ahat)
  unexpl <- 1 - rowSums(M)
  Mc <- M[combo_genes, , drop = FALSE]
  affil_combo[, , t + 1L] <- Mc
  rs <- rowSums(Mc); rs[rs <= 0] <- 1
  Pc <- Mc / rs
  ent <- -rowSums(Pc * log(pmax(Pc, 1e-12))) / log(K)

  # collapse diagnostics on diffused profiles
  within <- mean(c(Ct[inA, inA][upper.tri(Ct[inA, inA])],
                   Ct[inB, inB][upper.tri(Ct[inB, inB])],
                   Ct[inC, inC][upper.tri(Ct[inC, inC])]))
  between <- mean(c(as.vector(Ct[inA, inB]), as.vector(Ct[inA, inC]),
                    as.vector(Ct[inB, inC])))

  res <- rbind(res, data.frame(
    t = t, n_modules = length(found),
    largest = max(lengths(found)),
    hungarian_wf1 = sh$weighted_f1, best_mean_f1 = sb$mean_f1,
    combo_entropy = mean(ent), unexplained = mean(pmax(unexpl, 0)[combo_genes]),
    within_cor = within, between_cor = between
  ))
  if (t %% 5 == 0)
    cat(sprintf(paste0("t=%2d modules=%d largest=%d hungW-F1=%.3f ",
                       "bestF1=%.3f comboH=%.3f\n"),
                t, length(found), max(lengths(found)),
                sh$weighted_f1, sb$mean_f1, mean(ent)))
}
write.csv(res, file.path(outdir, "depth_metrics.csv"), row.names = FALSE)

## ---- plots ----
# 1. recovery vs depth
png(file.path(outdir, "f1_vs_depth.png"), width = 800, height = 520, res = 110)
par(mar = c(5, 4, 4, 2))
plot(res$t, res$hungarian_wf1, type = "b", pch = 16, col = "#1b6ca8",
     ylim = c(0, 1.05), xlab = "diffusion depth t", ylab = "F1",
     main = "Module recovery vs diffusion depth (Design A: fixed MI/operator)")
lines(res$t, res$best_mean_f1, type = "b", pch = 17, col = "#e8833a")
legend("topright", bty = "n", pch = c(16, 17), col = c("#1b6ca8", "#e8833a"),
       legend = c("Hungarian weighted F1", "independent-best mean F1"))
dev.off()

# 2. module count / collapse
png(file.path(outdir, "collapse_vs_depth.png"), width = 800, height = 520, res = 110)
par(mar = c(5, 4, 4, 4))
plot(res$t, res$n_modules, type = "b", pch = 16, col = "#2ca02c",
     xlab = "diffusion depth t", ylab = "number of found modules",
     main = "Module count and profile homogenization vs depth")
par(new = TRUE)
plot(res$t, res$between_cor, type = "l", lty = 2, col = "#d62728",
     axes = FALSE, xlab = "", ylab = "", ylim = range(c(res$within_cor, res$between_cor)))
lines(res$t, res$within_cor, lty = 1, col = "#d62728")
axis(4); mtext("mean correlation", side = 4, line = 2.5)
legend("center", bty = "n", lty = c(1, 2), col = "#d62728",
       legend = c("within planted modules", "between planted modules"))
dev.off()

# 3. affiliation trajectories for combinatorial genes
png(file.path(outdir, "affiliation_trajectories.png"), width = 960, height = 640, res = 110)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
cols <- hcl.colors(K, "Dark 3")
for (g in combo_genes) {
  matplot(0:T_steps, t(affil_combo[g, , ]), type = "l", lty = 1, lwd = 2,
          col = cols, ylim = c(0, 1),
          xlab = "depth t", ylab = "affiliation mass m_ik(t)",
          main = paste0(g, " (planted combo)"))
  legend("topright", bty = "n", cex = 0.8, lty = 1, lwd = 2, col = cols,
         legend = names(cores))
}
dev.off()

# 4. combo entropy + unexplained mass vs depth
png(file.path(outdir, "affiliation_summary.png"), width = 800, height = 520, res = 110)
par(mar = c(5, 4, 4, 2))
plot(res$t, res$combo_entropy, type = "b", pch = 16, col = "#9467bd",
     ylim = c(0, 1.05), xlab = "diffusion depth t", ylab = "",
     main = "Combinatorial genes: affiliation entropy and unexplained mass")
lines(res$t, res$unexplained, type = "b", pch = 17, col = "#8c564b")
legend("topleft", bty = "n", pch = c(16, 17), col = c("#9467bd", "#8c564b"),
       legend = c("normalized affiliation entropy", "unexplained mass"))
dev.off()

cat("\nDone. Outputs in", outdir, "\n")
print(res[, c("t", "n_modules", "largest", "hungarian_wf1", "best_mean_f1",
              "combo_entropy")], row.names = FALSE)
