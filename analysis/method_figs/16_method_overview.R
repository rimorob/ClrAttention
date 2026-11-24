# method_overview.png: the ladder CLR -> aCLR -> tcaCLR -> tcaCLR2 -> global-temperature variant -> multi-scale. One line per step (what it changes), the axis it touches,
# and the same 11-gene M3D neighbourhood read through each step (small matrices, each scaled to its own maximum; computed in 00_/01_compute_*.R).
# Layout inspired by a PaperBanana draft (numbered rows, axis tags, dashed exploratory rows); all content here is written and computed in code.
# Run:  nice -n 19 Rscript analysis/method_figs/16_method_overview.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); va <- readRDS("analysis/method_figs/cache/m3d_variants.rds")
H <- 6.7; RH <- 0.88
mats <- list(ex$S_sub, ex$aclr_sub, ex$S2_sub, va$S2_t2_sub, va$S2_v2_sub, ex$ms1)
steps <- list(
  list(n = "CLR", d = c("Scores each pair against each gene's own background: row z of the MI,", "clipped at 0, combined symmetrically (Stouffer)."), ax = list(c("normalisation", "#f6b595", ORT))),
  list(n = "aCLR", d = c("Turns the scores into attention: softmax over each gene's top 50 positive", "scores (about 10 effective partners), 9-step lazy walk, symmetric read-out."), ax = list(c("per-gene temperature", "#c9c2f4", VIT), c("walk depth", "#f5bfd3", MGT))),
  list(n = "tcaCLR", d = c("Adds a second CLR: row z of log P^9 over reachable pairs, clipped at 0,", "combined symmetrically; replaces the walk read-out."), ax = list(c("second CLR", "#f6b595", ORT))),
  list(n = "tcaCLR2", d = c("Replaces the fixed depth 9 by a depth from 10 per-gene permutations:", "t2 = max(1, round(ln 0.05 / ln edge)); no truth read."), ax = list(c("depth rule", "#f5bfd3", MGT))),
  list(n = "Global temperature", d = c("Replaces the per-gene temperature by one tau for all genes (mean", "effective partners 10); same walk and second CLR."), ax = list(c("temperature", "#c9c2f4", VIT)), expl = TRUE),
  list(n = "Multi-scale", d = c("Averages each pair's rank percentile over 13 temperatures", "(log2 tau -4 to +2); no tuning, no selection."), ax = list(c("temperature scale", "#c9c2f4", VIT)), expl = TRUE))
draw <- function() {
  txt("The method ladder: what each step changes", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  txt("axis touched", 4.75, Y(0.46), 7, ink[["secondary"]], just = "left", face = "italic"); txt("same 11 genes", 6.5, Y(0.46), 7, ink[["secondary"]], face = "italic")
  top <- 0.6
  for (k in seq_along(steps)) { st <- steps[[k]]; y0 <- top + (k - 1) * (RH + 0.06); ex_ <- isTRUE(st$expl)
    rect(0.1, Y(y0 + RH), 6.8, RH, fill = if (ex_) "#f1f0f6" else "#f7f6f3", col = if (ex_) ink[["muted"]] else NA, lty = if (ex_) 2 else 1, lwd = 0.8)
    grid.circle(0.38, Y(y0 + 0.2), 0.12, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA)); txt(k, 0.38, Y(y0 + 0.2) - 0.003, 8, "white", "bold")
    txt(st$n, 0.65, Y(y0 + 0.2), 8.5, face = "bold", just = "left")
    if (ex_) pill("exploratory", 0.65 + tw(st$n, 8.5, "bold") + 0.12, Y(y0 + 0.2), 0.17, 7, ink[["primary"]], "white", "bold", just = "left", padx = 0.14)
    lines_txt(st$d, 0.65, Y(y0 + 0.45), 0.15, 7, ink[["secondary"]])
    yy <- Y(y0 + 0.28); for (a in st$ax) { w <- pill(a[1], 4.75, yy, 0.2, 7, a[2], a[3], "plain", just = "left"); yy <- yy - 0.26 }
    M <- mats[[k]]; cell <- RH * 0.82 / nrow(M); x0 <- 6.15; yt <- Y(y0 + 0.07); diag(M) <- NA
    mat_draw(x0, yt, cell, replace(M, is.na(M), 0), rep("", nrow(M)), ramp = RAMP$grey, numbers = FALSE, top_labels = FALSE, left_labels = FALSE, gamma = 1.3, vmax = max(M, na.rm = TRUE))
  }
  for (k in 1:(length(steps) - 1)) { y0 <- top + (k - 1) * (RH + 0.06); harrow(0.38, Y(y0 + RH - 0.34), 0.38, Y(y0 + RH + 0.12), ink[["muted"]], 1) }
  yb <- top + length(steps) * (RH + 0.06) + 0.02
  lines_txt(c("Steps 1 to 4 are cumulative; steps 5 and 6 are exploratory variants (not preregistered). Small matrices: the 11 genes of the other figures",
              "(lexA, five top CLR partners, three genes ranked by the walk, two background genes), each matrix scaled to its own maximum."), 0.1, Y(yb + 0.1), 0.15, 7, ink[["secondary"]])
}
mf_save("overview", H, draw)
