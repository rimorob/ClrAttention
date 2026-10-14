#!/usr/bin/env Rscript
# Supplementary figure: B-spline mutual information against |Pearson r| for
# every gene pair, on M3D (466 experiments, one platform) and on EcoMAC
# (2,198 arrays, two platforms) without and with per-platform centring.
# The curve is the MI of a bivariate Gaussian, -1/2 log2(1 - r^2) bits: pairs
# far above it carry dependence that correlation does not capture (non-linear,
# or mixtures such as platform clusters). Reported per panel: the share of
# pairs whose MI exceeds the Gaussian value by more than `excess` bits, and
# the Spearman correlation between MI and |r|.
#
# Usage (repo root):  Rscript analysis/mi_vs_rho.R [--threads N] [--out results/figures]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(threads = NULL, out = "results/figures", excess = 0.1)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key == "threads") as.integer(args[i + 1L]) else
    if (key == "excess") as.numeric(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(clr))

m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
Xm <- as.matrix(raw[, -1]); Xm <- Xm[apply(Xm, 1, function(x) diff(range(x))) > 0, ]
E <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
Ec <- E
for (g in unique(plat)) { k <- plat == g
  mu <- rowMeans(E[, k]); s <- apply(E[, k], 1, stats::sd); s[s <= 0] <- 1
  Ec[, k] <- (E[, k] - mu) / s }

sets <- list(`M3D, 466 experiments (one platform)` = Xm,
             `EcoMAC, 2,198 arrays, uncorrected` = E,
             `EcoMAC, per-platform centred` = Ec)
pairs_of <- function(X) {
  mi <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$mi
  r <- abs(stats::cor(t(X)))
  ut <- upper.tri(r)
  data.frame(r = r[ut], mi = mi[ut])
}
dat <- lapply(sets, pairs_of)

ramp <- c("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
ink <- "#0b0b0b"; ink2 <- "#52514e"
f <- file.path(opt$out, "mi_vs_rho.png")
grDevices::png(f, width = 2400, height = 900, res = 200, bg = "#fcfcfb")
graphics::par(mfrow = c(1, 3), mar = c(3.6, 4, 3.4, 1), bg = "#fcfcfb", fg = ink2)
ymax <- max(vapply(dat, function(d) stats::quantile(d$mi, 0.9999), 0))
rr <- seq(0, 0.999, length.out = 400)
summ <- list()
for (nm in names(dat)) {
  d <- dat[[nm]]
  gauss <- -0.5 * log2(1 - pmin(d$r, 0.999)^2)
  ex <- mean(d$mi - gauss > opt$excess)
  rho <- stats::cor(d$mi, d$r, method = "spearman")
  summ[[nm]] <- data.frame(set = nm, n_pairs = nrow(d), share_excess = ex, spearman_mi_r = rho)
  # 2-D histogram, log counts, sequential blue
  bx <- seq(0, 1, length.out = 121); by <- seq(0, ymax, length.out = 121)
  h <- table(cut(d$r, bx, include.lowest = TRUE), cut(pmin(d$mi, ymax), by, include.lowest = TRUE))
  z <- log10(1 + matrix(h, length(bx) - 1))
  graphics::image(bx, by, z, col = grDevices::colorRampPalette(ramp)(64), xlab = "", ylab = "",
                  axes = FALSE, useRaster = TRUE)
  graphics::axis(1, col = ink2, col.axis = ink2, cex.axis = 0.8, lwd = 0.6)
  graphics::axis(2, col = ink2, col.axis = ink2, cex.axis = 0.8, lwd = 0.6, las = 1)
  graphics::lines(rr, -0.5 * log2(1 - rr^2), col = "#eb6834", lwd = 2)
  graphics::mtext("|Pearson r|", 1, 2.2, col = ink2, cex = 0.8)
  graphics::mtext("B-spline MI (bits)", 2, 2.6, col = ink2, cex = 0.8)
  graphics::mtext(nm, 3, 1.6, adj = 0, font = 2, col = ink, cex = 0.82)
  graphics::mtext(sprintf("%s pairs; %.1f%% exceed Gaussian MI by > %.1f bits; Spearman %.2f",
                          format(nrow(d), big.mark = ","), 100 * ex, opt$excess, rho),
                  3, 0.4, adj = 0, col = ink2, cex = 0.6)
}
graphics::legend("bottomright", legend = "Gaussian MI, -1/2 log2(1 - r^2)", col = "#eb6834", lwd = 2,
                 bty = "n", text.col = ink2, cex = 0.8)
grDevices::dev.off()
s <- do.call(rbind, summ)
utils::write.csv(s, file.path(opt$out, "mi_vs_rho_summary.csv"), row.names = FALSE)
print(s); cat("wrote", f, "\n")
