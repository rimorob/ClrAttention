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
# The pairwise MI is cached as 2-D bin counts in <cache>/mi_vs_rho_bins.rds
# (delete it to recompute); the figure is drawn from the cache.
#
# Usage (repo root):  Rscript analysis/mi_vs_rho.R [--threads N] [--out results/figures]
#                                                  [--cache results/mi_vs_rho]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(threads = NULL, out = "results/figures", excess = 0.1, cache = "results/mi_vs_rho")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key == "threads") as.integer(args[i + 1L]) else
    if (key == "excess") as.numeric(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(clr))

nbx <- 120L; nby <- 120L
cache_file <- file.path(opt$cache, "mi_vs_rho_bins.rds")
if (!file.exists(cache_file)) {
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
  summ <- list(); bins <- list()
  for (nm in names(sets)) {
    d <- pairs_of(sets[[nm]])
    gauss <- -0.5 * log2(1 - pmin(d$r, 0.999)^2)
    summ[[nm]] <- data.frame(set = nm, n_pairs = nrow(d), share_excess = mean(d$mi - gauss > opt$excess),
                             spearman_mi_r = stats::cor(d$mi, d$r, method = "spearman"),
                             mi_max = max(d$mi), share_below = mean(d$mi < gauss - opt$excess))
    bins[[nm]] <- d                                   # binned below, once the common y range is known
    cat(nm, "done\n")
  }
  ymax <- max(vapply(bins, function(d) max(d$mi), 0)) * 1.0001
  bx <- seq(0, 1, length.out = nbx + 1L); by <- seq(0, ymax, length.out = nby + 1L)
  bins <- lapply(bins, function(d) {
    h <- table(factor(findInterval(d$r, bx, rightmost.closed = TRUE), seq_len(nbx)),
               factor(findInterval(d$mi, by, rightmost.closed = TRUE), seq_len(nby)))
    matrix(as.integer(h), nbx)
  })
  dir.create(opt$cache, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(summary = do.call(rbind, summ), bins = bins, bx = bx, by = by), cache_file)
}
cache <- readRDS(cache_file)
s <- cache$summary
utils::write.csv(s, file.path(opt$out, "mi_vs_rho_summary.csv"), row.names = FALSE)
print(s)

## ---- figure ----------------------------------------------------------------------
source("analysis/plot_theme.R")
bx <- cache$bx; by <- cache$by
mids <- function(b) (b[-1] + b[-length(b)]) / 2
lev <- names(cache$bins)
tiles <- do.call(rbind, lapply(lev, function(nm) {
  h <- cache$bins[[nm]]
  g <- expand.grid(i = seq_len(nrow(h)), j = seq_len(ncol(h)))
  g$n <- h[cbind(g$i, g$j)]
  g <- g[g$n > 0, ]
  data.frame(set = nm, r = mids(bx)[g$i], mi = mids(by)[g$j], n = g$n)
}))
tiles$set <- factor(tiles$set, lev)
rr <- seq(0, 0.9995, length.out = 600)
curve <- data.frame(r = rr, mi = -0.5 * log2(1 - rr^2))
curve <- curve[curve$mi <= max(by), ]
corner <- data.frame(set = factor(lev, lev), pct = sprintf("%.1f%%", 100 * s$share_excess),
                     sub = sprintf("of %s pairs lie > %.1f bits above the curve\nSpearman(MI, |r|) = %.2f",
                                   format(s$n_pairs, big.mark = ","), opt$excess, s$spearman_mi_r))
curve_lab <- data.frame(set = factor(lev[1], lev), r = 0.5, mi = 0.43,
                        lab = "MI expected for a linear\n(Gaussian) relationship")
bias_lab <- data.frame(set = factor(lev[1], lev), r = 0.03, mi = 0.93,
                       lab = sprintf("Below the curve (%.1f%% of pairs by > %.1f bits):\ndownward estimator bias (7–9 joint bins),\nnot weaker dependence",
                                     100 * s$share_below[1], opt$excess))
floor_lab <- data.frame(set = factor(lev[1], lev), r = 0.03, mi = 0.29,
                        lab = "Floor near r = 0:\nfinite-sample bias")
p <- ggplot() +
  geom_tile(data = tiles, aes(r, mi, fill = n), width = diff(bx)[1], height = diff(by)[1]) +
  geom_line(data = curve, aes(r, mi), colour = ink[["primary"]], linewidth = 0.6, linetype = "22") +
  geom_text(data = curve_lab, aes(r, mi, label = lab), hjust = 1, vjust = 0, size = 2.6,
            colour = ink[["primary"]], lineheight = 0.9) +
  geom_segment(data = bias_lab, aes(x = 0.42, y = 0.8, xend = 0.74, yend = 0.3), colour = ink[["secondary"]],
               linewidth = 0.3, arrow = grid::arrow(length = unit(0.12, "cm"))) +
  geom_text(data = bias_lab, aes(r, mi, label = lab), hjust = 0, vjust = 1, size = 2.3,
            colour = ink[["secondary"]], lineheight = 0.9, fontface = "italic") +
  geom_text(data = floor_lab, aes(r, mi, label = lab), hjust = 0, vjust = 0, size = 2.3,
            colour = ink[["secondary"]], lineheight = 0.9, fontface = "italic") +
  geom_text(data = corner, aes(x = 0.02, y = max(by) * 0.985, label = pct), hjust = 0, vjust = 1,
            size = 6.5, fontface = "bold", colour = ink[["primary"]]) +
  geom_text(data = corner, aes(x = 0.02, y = max(by) * 0.86, label = sub), hjust = 0, vjust = 1,
            size = 2.5, colour = ink[["secondary"]], lineheight = 0.95) +
  facet_wrap(~set, nrow = 1) +
  scale_fill_gradientn(colours = seq_blue[-1], trans = "log10", labels = label_comma(),
                       name = "Gene pairs per bin (log scale)",
                       guide = guide_colourbar(barwidth = unit(6, "cm"), barheight = unit(0.3, "cm"),
                                               title.position = "top")) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(limits = c(0, max(by)), expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = sprintf("Excess MI over correlation comes from mixing platforms (%.1f%% of pairs vs %.1f–%.1f%%)",
                       100 * s$share_excess[2], 100 * s$share_excess[1], 100 * s$share_excess[3]),
       subtitle = "B-spline mutual information (shared histogram bins) against |Pearson r| for every gene pair; full range shown. Pairs far above the dashed curve carry\ndependence that correlation does not capture: in uncorrected EcoMAC these are platform mixtures, and per-platform centring removes nearly all of them.",
       x = "|Pearson r|", y = "Mutual information (bits)") +
  theme_clr() +
  theme(legend.position = "bottom", legend.justification = "right", panel.spacing.x = unit(24, "pt"),
        panel.grid.major = element_line(colour = grDevices::adjustcolor(ink[["grid"]], 0.6), linewidth = 0.3))
f <- file.path(opt$out, "mi_vs_rho.png")
save_fig(p, f, width = 11, height = 4.9)
ggsave(sub("png$", "pdf", f), p, width = 11, height = 4.9, device = grDevices::cairo_pdf, bg = ink[["surface"]])
cat("wrote", f, "\n")
