#!/usr/bin/env Rscript
# Stage 2 summary: heatmaps of the pooled-AUPR curves (depth x keff; tcaCLR minus the registered cell, one panel per clip combination; aCLR minus the
# registered aCLR cell, one panel per clipS), the registered cell, the best cell (upper bound, descriptive only) and the grid cell nearest to the registered one.
#   Rscript analysis/tcaclr_depth/stage2_summary.R --dirs results/tcaclr_depth/stage2_armB_cyc2008,results/tcaclr_depth/stage2_k562_corum,... --out results/tcaclr_depth/stage2_summary
# Registered cell = the reference operator (keff 10, cap 50, clipped, t 9; clipO clipped for tcaCLR). Revision 2 additions not built here: regret, Spearman criterion
# validity against the stage 1 criteria, AUPR by complex-size stratum.
suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(TRUE); opt <- list(dirs = "", out = "results/tcaclr_depth/stage2_summary")
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1L]
stopifnot(nzchar(opt$dirs)); dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
D <- do.call(rbind, lapply(strsplit(opt$dirs, ",")[[1]], function(d) utils::read.csv(file.path(d, "curves.csv"), stringsAsFactors = FALSE)))
D$clipS <- as.logical(D$clipS); D$clipO <- as.logical(D$clipO)
grid_cap <- function(k) pmin(4L * k, 320L)
D$is_ref <- D$keff == 10 & D$cap == 50 & D$clipS & D$t == 9 & D$cap != grid_cap(D$keff)
lab <- function(clipS, clipO) ifelse(is.na(clipO), paste0("clipS ", ifelse(clipS, "clipped", "unclipped")), paste0("clipS ", ifelse(clipS, "clipped", "unclipped"), " / clipO ", ifelse(clipO, "clipped", "unclipped")))
tabs <- list(); idx <- 0L
for (tr in unique(D$truth)) for (sc in unique(D$scope[D$truth == tr])) {
  x <- D[D$truth == tr & D$scope == sc, ]; g <- x[!x$is_ref, ]
  for (me in c("tcaCLR", "aCLR")) {
    r <- x[x$method == me & x$is_ref & (is.na(x$clipO) | x$clipO), ]; if (!nrow(r)) next
    gm <- g[g$method == me, ]; b <- gm[which.max(gm$aupr), ]; n <- gm[gm$keff == 10 & gm$t == 9 & gm$clipS & (is.na(gm$clipO) | gm$clipO), ]
    idx <- idx + 1L
    tabs[[idx]] <- data.frame(truth = tr, scope = sc, method = me, n_pairs = r$n_pairs, n_pos = r$n_pos, base_rate = r$n_pos / r$n_pairs,
      registered_aupr = r$aupr, nearest_grid_cell_aupr_keff10_cap40_t9 = n$aupr, best_aupr = b$aupr, best_minus_registered = b$aupr - r$aupr,
      best_keff = b$keff, best_cap = b$cap, best_clipS = b$clipS, best_clipO = b$clipO, best_t = b$t,
      n_cells_above_registered = sum(gm$aupr > r$aupr), n_cells = nrow(gm))
    gm$d <- gm$aupr - r$aupr; gm$panel <- lab(gm$clipS, gm$clipO); gm$kf <- factor(gm$keff); gm$tf <- factor(gm$t)
    lim <- max(abs(gm$d), 1e-6)
    p <- ggplot(gm, aes(kf, tf, fill = d)) + geom_tile(colour = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%+.3f", d)), size = 1.9) +
      geom_point(data = gm[which.max(gm$aupr), ], shape = 0, size = 5, stroke = 0.8, colour = "black", show.legend = FALSE) +
      facet_wrap(~panel, nrow = 1) + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lim, lim), name = "AUPR minus\nregistered") +
      labs(x = "keff (cap = 4 x keff, at most 320)", y = "depth t", title = sprintf("%s, %s: %s pooled AUPR minus the registered cell", tr, sc, me),
           subtitle = sprintf("registered cell (keff 10, cap 50, t 9, clipped) AUPR %.4f, base rate %.4f; square = best cell (%.4f, descriptive upper bound)", r$aupr, r$n_pos / r$n_pairs, b$aupr)) +
      theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 7))
    ggsave(file.path(opt$out, sprintf("heatmap_%s_%s_%s.png", tr, sc, me)), p, width = if (me == "tcaCLR") 15 else 8, height = 5.5, dpi = 150)
  }
}
tab <- do.call(rbind, tabs); utils::write.csv(tab, file.path(opt$out, "stage2_registered_and_best.csv"), row.names = FALSE)
op <- options(width = 220, digits = 4); print(tab, row.names = FALSE); options(op)
