# End-to-end small synthetic test of the module-evaluation pipeline.
#
# 1. Simulate expression for 24 genes x 100 samples with 3 planted modules
#    (sizes 8/7/6, within-module correlation ~0.72) + 3 pure-noise genes.
# 2. Full CLR-attention pipeline: MI -> CLR scores -> permutation threshold
#    (higher criticism) -> sparse attention operator.
# 3. find_modules() on the operator (connected components).
# 4. ModuleTester: found modules vs planted modules (Hungarian + soft).
#
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
outdir <- demo_dir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

## ---- 1. synthetic data ----
planted <- list(
  module_A = paste0("g", 1:8),
  module_B = paste0("g", 9:15),
  module_C = paste0("g", 16:21)
)
noise_genes <- paste0("g", 22:24)
genes <- c(unlist(planted, use.names = FALSE), noise_genes)
n_samples <- 100L
X <- matrix(rnorm(length(genes) * n_samples), nrow = length(genes),
            dimnames = list(genes, paste0("s", seq_len(n_samples))))
for (m in planted) {
  f <- rnorm(n_samples)                       # shared latent factor
  for (g in m) X[g, ] <- 0.85 * f + sqrt(1 - 0.85^2) * X[g, ]
}

## ---- 2. CLR attention pipeline ----
att <- ClrAttention$new(X)
att$estimate_mi()$calibrate()$select_threshold(B = 100)$build_operator()
op <- att$operator
cat("threshold:", att$threshold, " mean degree:", mean(rowSums(op > 0)), "\n")

## ---- 3. module finder ----
found <- find_modules(op)
cat("found", length(found), "modules:",
    paste(names(found), sprintf("(n=%d)", lengths(found)), collapse = " "), "\n")

## ---- 4. evaluate ----
mt <- ModuleTester$new(planted)
mt$set_found(found)
h <- mt$align(method = "hungarian")
b <- mt$align(method = "best")
mt$align(method = "soft", epsilon = 0.2)
s_soft <- mt$soft_alignment
u <- mt$soft_uncertainty
P <- mt$soft_assignment

cat("\n=== Hungarian headline alignment ===\n"); print(h)
cat("\n=== Soft (eps=0.2) expected alignment ===\n"); print(s_soft)
cat("\n=== Soft uncertainty (dummy mass, entropy) ===\n"); print(u)
cat("\n=== Headline summary ===\n"); print(mt$summary())

write.csv(h, file.path(outdir, "hungarian_alignment.csv"), row.names = FALSE)
write.csv(u, file.path(outdir, "soft_uncertainty.csv"))

## ---- plots ----
# 1. expression heatmap, rows ordered by planted module
ord_genes <- c(unlist(planted, use.names = FALSE), noise_genes)
Xs <- t(scale(t(X[ord_genes, ])))
modcol <- rep(c("#1b6ca8", "#e8833a", "#2ca02c", "#bbbbbb"),
              times = c(8, 7, 6, 3))
png(file.path(outdir, "expression_heatmap.png"), width = 800, height = 640, res = 110)
par(mar = c(4, 8, 4, 2))
image(x = seq_len(ncol(Xs)), y = seq_len(nrow(Xs)), z = t(Xs[nrow(Xs):1, ]),
      col = hcl.colors(64, "RdBu", rev = TRUE), axes = FALSE,
      xlab = "samples", ylab = "",
      main = "Synthetic expression (rows ordered by planted module)")
axis(2, at = seq_len(nrow(Xs)), labels = rev(rownames(Xs)), las = 1, cex.axis = 0.7)
abline(h = c(8.5, 15.5, 21.5), col = "black", lwd = 2)
mtext(c("A (8)", "B (7)", "C (6)", "noise (3)"), side = 4, at = c(4.5, 12, 18.5, 23),
      las = 1, cex = 0.8, line = 0.5)
dev.off()

# 2. soft assignment heatmap (regulons x found modules + dummy)
Paug <- cbind(P, dummy = u$dummy_mass)
png(file.path(outdir, "soft_assignment_heatmap.png"), width = 900, height = 560, res = 110)
par(mar = c(6, 9, 4, 2))
image(x = seq_len(ncol(Paug)), y = seq_len(nrow(Paug)),
      z = t(Paug[nrow(Paug):1, , drop = FALSE]),
      col = hcl.colors(64, "YlOrRd"), axes = FALSE,
      xlab = "", ylab = "",
      main = "Soft assignment P(planted module -> found module), eps = 0.2")
axis(1, at = seq_len(ncol(Paug)), labels = colnames(Paug), las = 2)
axis(2, at = seq_len(nrow(Paug)), labels = rev(rownames(Paug)), las = 1)
for (i in seq_len(nrow(Paug))) for (j in seq_len(ncol(Paug))) {
  v <- Paug[nrow(Paug) - i + 1, j]
  if (v > 0.02) text(j, i, sprintf("%.2f", v),
                     col = if (v > 0.45) "white" else "black", cex = 0.9)
}
dev.off()

# 3. per-regulon F1: Hungarian vs soft-expected vs independent-best
f1mat <- rbind(hungarian = h$f1, soft = s_soft$f1, best = b$f1)
colnames(f1mat) <- h$reference
png(file.path(outdir, "f1_by_module.png"), width = 800, height = 520, res = 110)
par(mar = c(5, 4, 4, 2))
bp <- barplot(f1mat, beside = TRUE, ylim = c(0, 1.12),
              col = c("#1b6ca8", "#e8833a", "#9a9a9a"),
              main = "Per-module F1: Hungarian headline vs soft expectation vs independent-best",
              ylab = "F1", xlab = "planted module", legend.text = TRUE,
              args.legend = list(x = "topright", bty = "n"))
text(x = as.vector(bp), y = as.vector(f1mat) + 0.045,
     labels = sprintf("%.2f", as.vector(f1mat)), cex = 0.85)
dev.off()

cat("\nDone. Outputs in", outdir, "\n")
