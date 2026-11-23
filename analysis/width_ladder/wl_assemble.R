#!/usr/bin/env Rscript
# Width ladder, assemble step (docs/plans/width_ladder_plan.md): response-surface tables (width by metric, per depth rule), best width per truth with edge flags, t table,
# ledger of every setting scored, descriptive heatmaps and curves (PNG). DESCRIPTIVE: no interval, no test; the best width is a selection on spent truths.
#   Rscript analysis/width_ladder/wl_assemble.R [--datasets m3d,ecomac,k562,yeast,kemmeren,hepg2] [--out results/width_ladder] [--tag ""]
# Reads <out>/<ds>/cells<tag>.csv (+ edge<tag>.tsv, gate.csv). Datasets whose files are missing are skipped with a message.
suppressPackageStartupMessages(source("analysis/width_ladder/wl_lib.R"))
O <- gl_options(commandArgs(TRUE), list(datasets = "m3d,ecomac,k562,yeast,kemmeren,hepg2", out = WL$out, tag = ""))
say <- gl_say(file.path(O$out, sprintf("assemble%s.log", O$tag))); dss <- strsplit(O$datasets, ",")[[1]]
W <- WL$widths; WID <- sprintf("W%g", W); MSN <- c("MS1", "MS2", "MS3"); CMP <- c("CLR", "tcaCLR", "tcaCLR2")
surf <- list(); tt <- list(); led <- list(); gates <- list()
for (ds in dss) {
  f <- file.path(O$out, ds, sprintf("cells%s.csv", O$tag)); if (!file.exists(f)) { say("skipped (no cells): ", ds); next }
  x <- utils::read.csv(f, stringsAsFactors = FALSE); fe <- file.path(O$out, ds, sprintf("edge%s.tsv", O$tag))
  if (file.exists(fe)) tt[[ds]] <- utils::read.delim(fe, stringsAsFactors = FALSE)
  fg <- file.path(O$out, ds, "gate.csv"); if (file.exists(fg)) { g <- utils::read.csv(fg, stringsAsFactors = FALSE); gates[[ds]] <- data.frame(dataset = ds, n_compared = sum(g$compared), max_abs_diff = max(g$abs_diff[g$compared]), n_fail = sum(!g$pass)) }
  key <- unique(x[, c("truth", "scope", "metric")])
  for (r in c("R9", "RE")) for (k in seq_len(nrow(key))) {
    s <- x[x$truth == key$truth[k] & x$scope == key$scope[k] & x$metric == key$metric[k], ]
    val <- function(nm, rule = r) { v <- s$value[s$score == nm & s$rule %in% c(rule, "both")]; if (length(v)) v[1] else NA_real_ }
    w <- vapply(WID, val, 0); row <- data.frame(dataset = ds, rule = r, truth = key$truth[k], scope = key$scope[k], metric = key$metric[k], stringsAsFactors = FALSE)
    for (nm in WID) row[[nm]] <- w[[nm]]; for (nm in MSN) row[[nm]] <- val(nm); for (nm in CMP) row[[nm]] <- val(nm, c("R9", "RE", "both"))
    if (any(is.finite(w))) { b <- which.max(w); row$best_keff <- W[b]; row$best_value <- w[[b]]
      row$best_at_edge <- if (b == 1L) "low edge (2)" else if (b == length(W)) "high edge (128)" else if (b == 2L) "next to low edge (3)" else if (b == length(W) - 1L) "next to high edge (96)" else "interior"
      row$best_over_tcaCLR <- row$best_value / row$tcaCLR; row$best_over_CLR <- row$best_value / row$CLR; row$MS1_over_best <- row$MS1 / row$best_value } 
    surf[[length(surf) + 1L]] <- row }
  ids <- unique(x[, c("rule", "score")]); for (k in seq_len(nrow(ids))) { n <- sum(x$rule == ids$rule[k] & x$score == ids$score[k])
    led[[length(led) + 1L]] <- data.frame(dataset = ds, rule = ids$rule[k], score = ids$score[k], kind = if (grepl("^W", ids$score[k])) "single width (setting)" else if (grepl("^MS", ids$score[k])) "multi-scale (setting)" else "comparator (not a setting)",
      n_metric_cells = n, stringsAsFactors = FALSE) }
}
if (!length(surf)) stop("no dataset has cells")
S <- do.call(rbind, surf); utils::write.csv(S, file.path(O$out, sprintf("surface_all%s.csv", O$tag)), row.names = FALSE)
for (ds in unique(S$dataset)) utils::write.csv(S[S$dataset == ds, ], file.path(O$out, ds, sprintf("surface%s.csv", O$tag)), row.names = FALSE)
B <- S[, c("dataset", "rule", "truth", "scope", "metric", "best_keff", "best_value", "best_at_edge", "tcaCLR", "tcaCLR2", "CLR", "MS1", "MS2", "MS3", "best_over_tcaCLR")]
utils::write.csv(B, file.path(O$out, sprintf("best_width%s.csv", O$tag)), row.names = FALSE)
L <- do.call(rbind, led); utils::write.csv(L, file.path(O$out, sprintf("ledger%s.csv", O$tag)), row.names = FALSE)
if (length(tt)) { T <- do.call(rbind, tt); utils::write.csv(T, file.path(O$out, sprintf("t_by_width%s.csv", O$tag)), row.names = FALSE) }
if (length(gates)) utils::write.csv(do.call(rbind, gates), file.path(O$out, sprintf("gate_summary%s.csv", O$tag)), row.names = FALSE)
say(sprintf("ledger: %d rows; settings %d, setting x metric cells %d", nrow(L), sum(L$kind != "comparator (not a setting)"), sum(L$n_metric_cells[L$kind != "comparator (not a setting)"])))
cat("best-width edge counts (single widths, per rule):\n"); print(table(S$rule, S$best_at_edge))

## ---- descriptive PNGs ---------------------------------------------------------------------------------------------------------
pal <- grDevices::hcl.colors(11, "Blue-Red 3", rev = FALSE)
cols <- c(WID, MSN, "CLR", "tcaCLR", "tcaCLR2"); lab <- c(as.character(W), MSN, CMP)
for (ds in unique(S$dataset)) for (r in c("R9", "RE")) {
  s <- S[S$dataset == ds & S$rule == r, ]; if (!nrow(s)) next
  s$row <- paste(s$truth, s$scope, s$metric); M <- as.matrix(s[, cols]); ref <- s$tcaCLR; R <- log2(M / ref)
  grDevices::png(file.path(O$out, ds, sprintf("heatmap_%s_%s%s.png", ds, r, O$tag)), width = 1500, height = 120 + 46 * nrow(s), res = 130)
  op <- graphics::par(mar = c(5, 17, 3, 1)); lim <- max(0.05, stats::quantile(abs(R[is.finite(R)]), 0.98)); Rc <- R; Rc[] <- pmax(-lim, pmin(lim, R))
  graphics::image(seq_along(cols), seq_len(nrow(s)), t(Rc[rev(seq_len(nrow(s))), , drop = FALSE]), col = pal, zlim = c(-lim, lim), axes = FALSE, xlab = "keff (ladder), multi-scale scores, comparators", ylab = "")
  graphics::axis(1, seq_along(cols), lab, las = 2, cex.axis = 0.8); graphics::axis(2, seq_len(nrow(s)), rev(s$row), las = 2, cex.axis = 0.65)
  graphics::abline(v = c(length(WID), length(WID) + 3) + 0.5, lwd = 2)
  for (k in seq_len(nrow(s))) if (is.finite(s$best_keff[k])) graphics::points(match(sprintf("W%g", s$best_keff[k]), cols), nrow(s) + 1 - k, pch = 19, cex = 0.9)
  graphics::title(sprintf("%s, depth rule %s: log2(score / tcaCLR keff 10, t 9); dot = best single width (descriptive, spent truths)", ds, r), cex.main = 0.8)
  graphics::par(op); grDevices::dev.off()
  pm <- unique(s[, c("truth", "scope", "metric")]); np <- nrow(pm); nc <- min(4L, np); nr <- ceiling(np / nc)
  grDevices::png(file.path(O$out, ds, sprintf("curves_%s_%s%s.png", ds, r, O$tag)), width = 420 * nc, height = 340 * nr, res = 110); op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
  for (k in seq_len(np)) { y <- s[s$truth == pm$truth[k] & s$scope == pm$scope[k] & s$metric == pm$metric[k], ]; v <- unlist(y[, WID]); rg <- range(c(v, unlist(y[, c(MSN, CMP)])), finite = TRUE)
    graphics::plot(log2(W), v, type = "b", pch = 19, ylim = rg, xaxt = "n", xlab = "keff", ylab = "value", main = paste(pm$truth[k], pm$scope[k], pm$metric[k], sep = "\n"), cex.main = 0.75); graphics::axis(1, log2(W), W, cex.axis = 0.7)
    graphics::abline(h = unlist(y[, c("CLR", "tcaCLR", "tcaCLR2")]), lty = c(3, 2, 4), col = c("grey40", "#1f77b4", "#d62728")); graphics::abline(h = unlist(y[, MSN]), lty = 1, col = c("#2ca02c", "#9467bd", "#ff7f0e"), lwd = 1) }
  graphics::par(op); grDevices::dev.off() }
say("PNGs written; legend of the curves in README: dotted grey CLR, dashed blue tcaCLR, dash-dot red tcaCLR2, solid green/purple/orange MS1/MS2/MS3")
