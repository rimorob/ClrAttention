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
#   Rscript analysis/ecomac_scatter.R [--pairs lexA:recA,lexA:sulA] [--name scatter_lexA]
#                                     [--out results/figures]
# Writes one figure, <name>.png/.pdf: one row per pair, columns M3D | EcoMAC all
# arrays | EcoMAC Genome 2.0 only, shared axis ranges along each row (x shared
# throughout). Arrays treated with a DNA-damaging agent (quinolones, mitomycin C,
# UV, nalidixic acid; from the M3D and EcoMAC condition annotations) are
# outlined by a 90% normal ellipse: the SOS-induced cluster.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(pairs = "lexA:recA,lexA:sulA", name = "scatter_lexA", out = "results/figures",
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
ed <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_descriptions")),
                        stringsAsFactors = FALSE)
dmg_m3d <- grepl("(norfloxacin|ciprofloxacin|mitomycin|nalidixic)[^,]*:", 
                 ed$description[match(colnames(M), ed$chip_name)], ignore.case = TRUE)
tr <- smp[["Treatment (conditions)"]]
dmg_eco <- grepl("^(nor|norf|norfloxacin)( |$)|ciprofloxacin|nalidixic|^MMC$|^UV$", tr)
pd <- utils::read.delim(file.path(opt$m3d, paste0(base, ".probe_set_descriptions")),
                        stringsAsFactors = FALSE)
nm <- pd[[3]][match(raw[[1]], pd[[1]])]                    # common name per probe set
rownames(M) <- ifelse(is.na(nm) | nm == "", raw[[1]], nm)

## ---- plotting --------------------------------------------------------------------
# Platform hues lie outside the method colour key (plot_theme.R); shapes repeat
# the platform so the figure survives greyscale printing.
source("analysis/plot_theme.R")
plat_cols <- c(platform_cols, "M3D (one platform)" = ink[["secondary"]])
plat_shapes <- c("Antisense (GPL199)" = 16, "Genome 2.0 (GPL3154)" = 17, "no GEO record" = 15,
                 "M3D (one platform)" = 16)
pair_mi <- function(x, y) {      # B-spline MI, HG bins; clr reports MI in bits (log2)
  f <- ClrAttention$new(rbind(x, y))$estimate_mi(bins = "hg", transform = "none", threads = 1L)
  f$mi[1, 2]
}
cnt <- table(factor(platform, names(plat_cols)[1:3]))
cols_lab <- c(m3d = "M3D build 6: 907 arrays\n(= 466 experiments after averaging replicates)",
              all = sprintf("EcoMAC, all %s arrays\n(three platform groups pooled)", format(ncol(E), big.mark = ",")),
              g2 = sprintf("EcoMAC, Genome 2.0 arrays only\n(%d arrays, one platform)", cnt[["Genome 2.0 (GPL3154)"]]))
pairs <- strsplit(strsplit(opt$pairs, ",")[[1]], ":")
pts <- list(); st <- list()
set.seed(1); o <- sample(ncol(E))                         # no platform drawn last on top
g2 <- platform == "Genome 2.0 (GPL3154)"
for (g in pairs) {
  if (!all(g %in% rownames(E)) || !all(g %in% rownames(M))) { cat("skip", g, "\n"); next }
  rl <- sprintf("%s vs %s", g[2], g[1])
  sets <- list(m3d = list(x = M[g[1], ], y = M[g[2], ], p = rep("M3D (one platform)", ncol(M)), d = dmg_m3d),
               all = list(x = E[g[1], o], y = E[g[2], o], p = platform[o], d = dmg_eco[o]),
               g2  = list(x = E[g[1], g2], y = E[g[2], g2], p = platform[g2], d = dmg_eco[g2]))
  for (k in names(sets)) { s <- sets[[k]]
    pts[[paste(rl, k)]] <- data.frame(row = rl, col = cols_lab[[k]], x = s$x, y = s$y, platform = s$p, dmg = s$d)
    st[[paste(rl, k)]] <- data.frame(row = rl, col = cols_lab[[k]], n = length(s$x), r = stats::cor(s$x, s$y),
                                     mi = pair_mi(s$x, s$y), ndmg = sum(s$d))
  }
}
pts <- do.call(rbind, pts); st <- do.call(rbind, st)
rows <- unique(pts$row)
pts$row <- factor(pts$row, rows); st$row <- factor(st$row, rows)
pts$col <- factor(pts$col, cols_lab); st$col <- factor(st$col, cols_lab)
print(st)
utils::write.csv(st, file.path(opt$out, paste0(opt$name, "_stats.csv")), row.names = FALSE)
st$lab <- sprintf("n = %s arrays\nPearson r = %.2f\nMI = %.2f bits", trimws(format(st$n, big.mark = ",")), st$r, st$mi)
ell <- pts[pts$dmg & ave(pts$dmg, pts$row, pts$col, FUN = sum) >= 10, ]
# Label the SOS cluster in each panel that has one, to the upper left of its ellipse.
ell_lab <- do.call(rbind, lapply(split(ell, list(ell$row, ell$col), drop = TRUE), function(e)
  data.frame(row = e$row[1], col = e$col[1], x = stats::median(e$x) - 2.7 * stats::mad(e$x), y = stats::median(e$y) + stats::mad(e$y))))
ell_lab$lab <- "SOS-induced arrays\n(DNA-damage treatment)"
p <- ggplot(pts, aes(x, y)) +
  geom_point(aes(colour = platform, shape = platform), size = 0.75, alpha = 0.5, stroke = 0) +
  stat_ellipse(data = ell, level = 0.9, colour = ink[["primary"]], linetype = "22", linewidth = 0.45) +
  geom_text(data = ell_lab, aes(label = lab), size = 2.4, colour = ink[["primary"]], hjust = 1, vjust = 0,
            fontface = "italic",
            lineheight = 0.9) +
  geom_label(data = st, aes(x = -Inf, y = Inf, label = lab), hjust = -0.04, vjust = 1.08, size = 2.55,
             linewidth = 0, fill = grDevices::adjustcolor(ink[["surface"]], 0.85), colour = ink[["primary"]],
             lineheight = 0.95) +
  facet_grid(row ~ col, scales = "free_y", switch = "y") +
  scale_colour_manual(values = plat_cols, breaks = names(plat_cols)[c(4, 1:3)],
                      labels = c("M3D (one platform)",
                                 sprintf("EcoMAC %s: %d arrays", names(cnt), as.integer(cnt))), name = NULL) +
  scale_shape_manual(values = plat_shapes, breaks = names(plat_cols)[c(4, 1:3)],
                     labels = c("M3D (one platform)",
                                sprintf("EcoMAC %s: %d arrays", names(cnt), as.integer(cnt))), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1), nrow = 2)) +
  labs(title = "Mixing platforms creates the ‘Christmas tree’ and dilutes a known LexA–target relationship",
       subtitle = sprintf("lexA against two SOS genes it represses (rows). In M3D the relationship is monotonic and the SOS-induced arrays extend it. Pooled EcoMAC arrays\nsplit into platform clouds at different heights, so r and MI fall well below M3D; Genome 2.0 arrays alone (1 DNA-damage array of %d) show only a weak relation.\n'No GEO record': %d EcoMAC arrays (M3D/ASAP chip ids) without a GEO series, so their platform could not be looked up.", cnt[["Genome 2.0 (GPL3154)"]], cnt[["no GEO record"]]),
       x = "lexA expression (log2, RMA)", y = "Target expression (log2, RMA)",
       caption = "Dashed ellipse: 90% normal ellipse of arrays annotated with a DNA-damaging treatment (quinolone, mitomycin C, UV or nalidixic acid). MI: B-spline estimator with shared histogram bins.") +
  theme_clr() +
  theme(legend.position = "top", legend.justification = "left", legend.text = element_text(size = 9),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90, hjust = 0.5, face = "bold"),
        panel.spacing.x = unit(10, "pt"), panel.spacing.y = unit(12, "pt"))
f <- file.path(opt$out, paste0(opt$name, ".png"))
save_fig(p, f, width = 10, height = 7.6)
ggsave(sub("png$", "pdf", f), p, width = 10, height = 7.6, device = grDevices::cairo_pdf, bg = ink[["surface"]])
cat("wrote", f, "\n")
