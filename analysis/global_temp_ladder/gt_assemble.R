#!/usr/bin/env Rscript
# Global-temperature ladder, assemble step (docs/plans/global_temperature_ladder_plan.md): response-surface tables (temperature by metric, per depth rule), best temperature per truth with edge
# flags, effective-partner table, t table, ledger of every setting scored, gate summary, descriptive heatmaps and curves (PNG). DESCRIPTIVE: no interval, no test; the best temperature is a
# selection on spent truths. Differences from analysis/width_ladder/wl_assemble.R (that file is NOT changed): (1) the heatmap device is sized from the number of rows AND the margins
# (the width-ladder script used height = 120 + 46 n px at 130 dpi with margins of 8 lines, so a one-row dataset (yeast, Kemmeren, HepG2) failed with "figure margins too large"), and every
# plot is wrapped so that one failing figure cannot stop the tables; (2) the axis is log2 tau with the V2 reference temperature drawn as a vertical line.
#   Rscript analysis/global_temp_ladder/gt_assemble.R [--datasets m3d,ecomac,k562,yeast,kemmeren,hepg2] [--out results/global_temp_ladder] [--tag ""]
suppressPackageStartupMessages(source("analysis/global_temp_ladder/gt_lib.R"))
O <- gl_options(commandArgs(TRUE), list(datasets = "m3d,ecomac,k562,yeast,kemmeren,hepg2", out = GTL$out, tag = ""))
dir.create(O$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(O$out, sprintf("assemble%s.log", O$tag))); dss <- strsplit(O$datasets, ",")[[1]]
MSN <- c("MS1", "MS2", "MS3"); CMP <- c("CLR", "tcaCLR", "tcaCLR2")
surf <- list(); tt <- list(); led <- list(); gates <- list(); eff <- list(); chk <- list()
for (ds in dss) {
  f <- file.path(O$out, ds, sprintf("cells%s.csv", O$tag)); if (!file.exists(f)) { say("skipped (no cells): ", ds); next }
  x <- utils::read.csv(f, stringsAsFactors = FALSE); fe <- file.path(O$out, ds, sprintf("edge%s.tsv", O$tag)); if (file.exists(fe)) tt[[ds]] <- utils::read.delim(fe, stringsAsFactors = FALSE)
  ft <- file.path(O$out, ds, sprintf("tau_table%s.csv", O$tag)); if (file.exists(ft)) eff[[ds]] <- utils::read.csv(ft, stringsAsFactors = FALSE)
  fg <- file.path(O$out, ds, "gate.csv"); if (file.exists(fg)) { g <- utils::read.csv(fg, stringsAsFactors = FALSE)
    gates[[ds]] <- data.frame(dataset = ds, n_compared_gate1 = sum(g$compared & g$gate == 1L), max_abs_diff_gate1 = max(g$abs_diff[g$compared & g$gate == 1L]), n_compared_gate2 = sum(g$compared & g$gate == 2L),
                              max_abs_diff_gate2 = if (any(g$gate == 2L & g$compared)) max(g$abs_diff[g$compared & g$gate == 2L]) else NA_real_, n_fail = sum(!g$pass)) }
  LAD <- unique(x$score[grepl("^T[+-]", x$score)]); LAD <- LAD[order(as.numeric(sub("^T", "", LAD)))]; L2 <- as.numeric(sub("^T", "", LAD))
  key <- unique(x[, c("truth", "scope", "metric")])
  for (r in c("R9", "RE")) for (k in seq_len(nrow(key))) {
    s <- x[x$truth == key$truth[k] & x$scope == key$scope[k] & x$metric == key$metric[k], ]
    val <- function(nm, rule = r) { v <- s$value[s$score == nm & s$rule %in% c(rule, "both")]; if (length(v)) v[1] else NA_real_ }
    w <- vapply(LAD, val, 0); row <- data.frame(dataset = ds, rule = r, truth = key$truth[k], scope = key$scope[k], metric = key$metric[k], stringsAsFactors = FALSE)
    for (j in seq_along(LAD)) row[[LAD[j]]] <- w[[j]]; row$V2 <- val("V2"); for (nm in MSN) row[[nm]] <- val(nm); for (nm in CMP) row[[nm]] <- val(nm, c("R9", "RE", "both"))
    if (any(is.finite(w))) { b <- which.max(w); n <- length(w); row$best_log2tau <- L2[b]; row$best_value <- w[[b]]
      row$best_at_edge <- if (b == 1L) sprintf("low edge (%g)", L2[1]) else if (b == n) sprintf("high edge (%g)", L2[n]) else if (b == 2L) sprintf("next to low edge (%g)", L2[2]) else if (b == n - 1L) sprintf("next to high edge (%g)", L2[n - 1L]) else "interior"
      row$best_over_tcaCLR <- row$best_value / row$tcaCLR; row$best_over_CLR <- row$best_value / row$CLR; row$best_over_V2 <- row$best_value / row$V2; row$MS1_over_best <- row$MS1 / row$best_value }
    surf[[length(surf) + 1L]] <- row }
  ids <- unique(x[, c("rule", "score")]); for (k in seq_len(nrow(ids))) { n <- sum(x$rule == ids$rule[k] & x$score == ids$score[k])
    led[[length(led) + 1L]] <- data.frame(dataset = ds, rule = ids$rule[k], score = ids$score[k], kind = if (grepl("^T[+-]", ids$score[k])) "single ladder temperature (setting)" else if (ids$score[k] == "V2") "V2 reference temperature (setting)" else if (grepl("^MS", ids$score[k])) "multi-scale (setting)" else "comparator (not a setting)",
      n_metric_cells = n, stringsAsFactors = FALSE) }
}
if (!length(surf)) stop("no dataset has cells")
rbind_fill <- function(l) { cn <- unique(unlist(lapply(l, names))); do.call(rbind, lapply(l, function(d) { for (m in setdiff(cn, names(d))) d[[m]] <- NA; d[, cn] })) }   # datasets differ in the ladder columns only in smoke subsets
S <- rbind_fill(surf)
utils::write.csv(S, file.path(O$out, sprintf("surface_all%s.csv", O$tag)), row.names = FALSE)
for (ds in unique(S$dataset)) utils::write.csv(S[S$dataset == ds, ], file.path(O$out, ds, sprintf("surface%s.csv", O$tag)), row.names = FALSE)
bc <- intersect(c("dataset", "rule", "truth", "scope", "metric", "best_log2tau", "best_value", "best_at_edge", "V2", "tcaCLR", "tcaCLR2", "CLR", "MS1", "MS2", "MS3", "best_over_tcaCLR", "best_over_V2"), names(S))
utils::write.csv(S[, bc], file.path(O$out, sprintf("best_tau%s.csv", O$tag)), row.names = FALSE)
L <- do.call(rbind, led); utils::write.csv(L, file.path(O$out, sprintf("ledger%s.csv", O$tag)), row.names = FALSE)
if (length(tt)) utils::write.csv(do.call(rbind, tt), file.path(O$out, sprintf("t_by_tau%s.csv", O$tag)), row.names = FALSE)
if (length(eff)) utils::write.csv(do.call(rbind, eff), file.path(O$out, sprintf("effective_partners%s.csv", O$tag)), row.names = FALSE)
if (length(gates)) utils::write.csv(do.call(rbind, gates), file.path(O$out, sprintf("gate_summary%s.csv", O$tag)), row.names = FALSE)
say(sprintf("ledger: %d rows; settings %d, setting x metric cells %d", nrow(L), sum(L$kind != "comparator (not a setting)"), sum(L$n_metric_cells[L$kind != "comparator (not a setting)"])))
cat("best-temperature edge counts (13 single ladder temperatures, per rule):\n"); print(table(S$rule, sub(" \\(.*", "", S$best_at_edge)))

## ---- descriptive PNGs ---------------------------------------------------------------------------------------------------------
pal <- grDevices::hcl.colors(11, "Blue-Red 3", rev = FALSE)
for (ds in unique(S$dataset)) for (r in c("R9", "RE")) {
  s <- S[S$dataset == ds & S$rule == r, ]; if (!nrow(s)) next
  LAD <- grep("^T[+-]", names(s), value = TRUE); LAD <- LAD[order(as.numeric(sub("^T", "", LAD)))]; L2 <- as.numeric(sub("^T", "", LAD)); cols <- c(LAD, "V2", MSN, CMP); lab <- c(sub("^T", "", LAD), "V2", MSN, CMP)
  s$row <- paste(s$truth, s$scope, s$metric); M <- as.matrix(s[, cols]); R <- log2(M / s$tcaCLR); n <- nrow(s)
  ok <- try({
    # device sized in inches from the margins: left 16 lines, bottom 6, top 3, right 1 lines (0.2 in each at cex 1) plus 0.32 in per row and 4 in of plot floor
    mar <- c(6, 16, 3, 1); h_in <- max(4.5, (mar[1] + mar[3]) * 0.2 + 0.32 * n + 1.5); w_in <- 11.5
    grDevices::png(file.path(O$out, ds, sprintf("heatmap_%s_%s%s.png", ds, r, O$tag)), width = w_in, height = h_in, units = "in", res = 130)
    op <- graphics::par(mar = mar); lim <- max(0.05, stats::quantile(abs(R[is.finite(R)]), 0.98)); Rc <- R; Rc[] <- pmax(-lim, pmin(lim, R))
    graphics::image(seq_along(cols), seq_len(n), t(Rc[rev(seq_len(n)), , drop = FALSE]), col = pal, zlim = c(-lim, lim), axes = FALSE, xlab = "log2 tau (ladder), V2 reference, multi-scale scores, comparators", ylab = "")
    graphics::axis(1, seq_along(cols), lab, las = 2, cex.axis = 0.8); graphics::axis(2, seq_len(n), rev(s$row), las = 2, cex.axis = 0.65)
    graphics::abline(v = c(length(LAD), length(LAD) + 1, length(LAD) + 4) + 0.5, lwd = 2)
    for (k in seq_len(n)) if (is.finite(s$best_log2tau[k])) graphics::points(match(sprintf("T%+.2f", s$best_log2tau[k]), cols), n + 1 - k, pch = 19, cex = 0.9)
    graphics::title(sprintf("%s, depth rule %s: log2(score / tcaCLR keff 10, t 9); dot = best single ladder temperature (descriptive, spent truths)", ds, r), cex.main = 0.8)
    graphics::par(op); grDevices::dev.off() }, silent = TRUE); if (inherits(ok, "try-error")) { while (grDevices::dev.cur() > 1L) grDevices::dev.off(); say("heatmap failed: ", ds, " ", r, ": ", conditionMessage(attr(ok, "condition"))) }
  ok <- try({
    pm <- unique(s[, c("truth", "scope", "metric")]); np <- nrow(pm); nc <- min(4L, np); nr <- ceiling(np / nc)
    grDevices::png(file.path(O$out, ds, sprintf("curves_%s_%s%s.png", ds, r, O$tag)), width = max(480, 420 * nc), height = max(420, 340 * nr), res = 110); op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
    tv2 <- if (!is.null(eff[[ds]])) eff[[ds]]$log2_tau[eff[[ds]]$op == "V2"] else NA_real_
    for (k in seq_len(np)) { y <- s[s$truth == pm$truth[k] & s$scope == pm$scope[k] & s$metric == pm$metric[k], ]; v <- unlist(y[, LAD]); rg <- range(c(v, unlist(y[, c(MSN, CMP, "V2")])), finite = TRUE)
      graphics::plot(L2, v, type = "b", pch = 19, ylim = rg, xlab = "log2 tau", ylab = "value", main = paste(pm$truth[k], pm$scope[k], pm$metric[k], sep = "\n"), cex.main = 0.75)
      graphics::abline(h = unlist(y[, c("CLR", "tcaCLR", "tcaCLR2")]), lty = c(3, 2, 4), col = c("grey40", "#1f77b4", "#d62728")); graphics::abline(h = unlist(y[, MSN]), lty = 1, col = c("#2ca02c", "#9467bd", "#ff7f0e"), lwd = 1)
      if (is.finite(tv2)) { graphics::abline(v = tv2, lty = 3, col = "#8c564b"); graphics::points(tv2, y$V2, pch = 17, col = "#8c564b", cex = 1.2) } }
    graphics::par(op); grDevices::dev.off() }, silent = TRUE); if (inherits(ok, "try-error")) { while (grDevices::dev.cur() > 1L) grDevices::dev.off(); say("curves failed: ", ds, " ", r, ": ", conditionMessage(attr(ok, "condition"))) }
}
say("PNGs written; legend of the curves in README: dotted grey CLR, dashed blue tcaCLR, dash-dot red tcaCLR2, solid green/purple/orange MS1/MS2/MS3, brown triangle and vertical line = V2 reference temperature of the dataset")
