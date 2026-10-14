#!/usr/bin/env Rscript
# Gene-gene scatter plots across compendia, to show how platform and batch
# structure bends pairwise relationships: M3D E. coli v4 build 6 (907 chips,
# one Affymetrix platform) against EcoMAC (Carrera et al. 2014; 2,198 arrays
# on the E. coli Antisense (GPL199) and Genome 2.0 (GPL3154) platforms).
# Each panel reports Pearson r and B-spline MI (HG bins) for the pair.
#
# Platform per EcoMAC array: GEO series platform (data/ecomac/gse_platform.tsv),
# per-sample platform for series that mix platforms (gsm_platform.tsv); arrays
# without a GEO series (ASAP, M3D chip ids) are labelled "no GEO record".
#
# Usage (repo root):
#   Rscript analysis/ecomac_scatter.R [--pairs lexA:recA,...] [--out results/figures]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(pairs = "lexA:recA,lexA:sulA,crp:cyaA,fnr:narG", out = "results/figures",
            ecomac = "data/ecomac", m3d = "data/E_coli_v4_Build_6")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- args[i + 1L]; i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(clr))

## ---- EcoMAC ---------------------------------------------------------------------
E <- as.matrix(utils::read.delim(file.path(opt$ecomac, "EcoMAC.tsv"), header = FALSE))
gn <- utils::read.delim(file.path(opt$ecomac, "gene_names.tsv"), header = FALSE,
                        stringsAsFactors = FALSE)
stopifnot(nrow(E) == nrow(gn))
rownames(E) <- gn[[3]]
smp <- utils::read.delim(file.path(opt$ecomac, "samples.tsv"), stringsAsFactors = FALSE,
                         check.names = FALSE)
stopifnot(nrow(smp) == ncol(E))
gp <- utils::read.delim(file.path(opt$ecomac, "gse_platform.tsv"), header = FALSE,
                        stringsAsFactors = FALSE, col.names = c("gse", "gpl"))
sp <- utils::read.delim(file.path(opt$ecomac, "gsm_platform.tsv"), header = FALSE,
                        stringsAsFactors = FALSE, col.names = c("gsm", "gpl"))
plat <- gp$gpl[match(smp[["GEO accession"]], gp$gse)]
mixed <- !is.na(plat) & grepl(",", plat)
plat[mixed] <- sp$gpl[match(smp[["CEL file name"]][mixed], sp$gsm)]
plat[is.na(plat) | plat == ""] <- "none"
platform <- ifelse(plat == "GPL199", "Antisense (GPL199)",
            ifelse(plat == "GPL3154", "Genome 2.0 (GPL3154)", "no GEO record"))
utils::write.csv(data.frame(array = seq_len(ncol(E)), geo = smp[["GEO accession"]],
                            sample = smp[["CEL file name"]], author = smp[["Author"]],
                            platform = platform),
                 file.path(opt$ecomac, "array_platform.csv"), row.names = FALSE)
cat("EcoMAC arrays by platform:\n"); print(table(platform))

## ---- M3D build 6 ----------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
raw <- utils::read.delim(file.path(opt$m3d, paste0(base, "_chips907probes4297.tab")),
                         check.names = FALSE, stringsAsFactors = FALSE)
M <- as.matrix(raw[, -1])
pd <- utils::read.delim(file.path(opt$m3d, paste0(base, ".probe_set_descriptions")),
                        stringsAsFactors = FALSE)
nm <- pd[[3]][match(raw[[1]], pd[[1]])]                    # common name per probe set
rownames(M) <- ifelse(is.na(nm) | nm == "", raw[[1]], nm)

## ---- plotting --------------------------------------------------------------------
# Categorical slots 1-3 of the validated reference palette (all-pairs safe for scatter).
cols <- c("Antisense (GPL199)" = "#2a78d6", "Genome 2.0 (GPL3154)" = "#eb6834",
          "no GEO record" = "#1baf7a")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"
pair_mi <- function(x, y) {
  f <- ClrAttention$new(rbind(x, y))$estimate_mi(bins = "hg", transform = "none", threads = 1L)
  f$mi[1, 2]
}
panel <- function(x, y, col, main, xlab, ylab) {
  plot(x, y, type = "n", xlab = "", ylab = "", axes = FALSE, main = "")
  grid(col = grid_col, lty = 1)
  axis(1, col = ink2, col.axis = ink2, cex.axis = 0.8, lwd = 0.6)
  axis(2, col = ink2, col.axis = ink2, cex.axis = 0.8, lwd = 0.6, las = 1)
  mtext(xlab, 1, 2.2, col = ink2, cex = 0.8); mtext(ylab, 2, 2.6, col = ink2, cex = 0.8)
  points(x, y, pch = 16, cex = 0.45, col = grDevices::adjustcolor(col, 0.55))
  mtext(main, 3, 1.3, adj = 0, font = 2, col = ink, cex = 0.85)
  mtext(sprintf("n = %d   Pearson r = %.2f   MI = %.2f", length(x), stats::cor(x, y),
                pair_mi(x, y)), 3, 0.3, adj = 0, col = ink2, cex = 0.72)
}
for (pr in strsplit(opt$pairs, ",")[[1]]) {
  g <- strsplit(pr, ":")[[1]]
  if (!all(g %in% rownames(E)) || !all(g %in% rownames(M))) { cat("skip", pr, "\n"); next }
  f <- file.path(opt$out, sprintf("scatter_%s_%s.png", g[1], g[2]))
  grDevices::png(f, width = 2400, height = 1150, res = 200, bg = "#fcfcfb")
  graphics::par(mfrow = c(1, 3), mar = c(3.6, 4, 3.2, 1), oma = c(2.2, 0, 0, 0), bg = "#fcfcfb",
                fg = ink2)
  panel(M[g[1], ], M[g[2], ], cols[["Antisense (GPL199)"]],
        "M3D build 6 (one platform)", g[1], g[2])
  o <- sample(ncol(E))                                  # no platform drawn last on top
  panel(E[g[1], o], E[g[2], o], cols[platform[o]], "EcoMAC, all arrays", g[1], g[2])
  g2 <- platform == "Genome 2.0 (GPL3154)"
  panel(E[g[1], g2], E[g[2], g2], cols[["Genome 2.0 (GPL3154)"]],
        "EcoMAC, Genome 2.0 arrays only", g[1], g[2])
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", legend = sprintf("%s (%d)", names(cols), as.integer(table(factor(platform, names(cols))))),
         pch = 16, col = cols, horiz = TRUE, bty = "n", text.col = ink2, cex = 0.85, inset = 0.01)
  grDevices::dev.off()
  cat("wrote", f, "\n")
}
