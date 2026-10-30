# Fig. 1 (method): CLR as attention, one traced cell through a worked toy example.
# Layout logic follows Faith et al. 2007 (PLoS Biol, Fig. 1A): one matrix drawn large, gene i as a
# horizontal band and gene j as a vertical band crossing at the traced cell, the two background
# distributions drawn next to their band with the observed value marked in the tail, formulas beside
# the arrows, a big arrow to the output matrix.
#
# REAL TOY: every matrix is computed with the package's own functions on a small synthetic
# expression matrix with planted structure (a chain module, a one-factor module, noise genes):
#   D     bspline_mi() via ClrAttention$estimate_mi()        (diagonal zeroed as clr_calibrate does)
#   CLR   clr_calibrate(method = "normal", combine = "stouffer") via $calibrate()
#   A     .softmax_rows(S, keff, cap) via $build_operator(softmax_keff, softmax_cap)
#   P^t   P = (1-alpha) I + alpha A, repeated product, K_t = (P^t + P^t')/2: the exact lines of
#         analysis/baselines.R (lines "P <- Matrix::Diagonal(...)" ... "St <- (Pt + t(Pt)) / 2"),
#         which is how the walk is evaluated in the paper (ClrAttention$diffuse() applies the same
#         P to expression profiles instead).
# The z-score curves are recomputed from the same formulas as src/core/clr_core.cpp and ASSERTED
# equal to the clr_calibrate() output, so the picture cannot drift from the code.
# Toy parameters differ from the real method only where a 14-gene matrix forces it:
#   cap 5 (real 50), k_eff 3 (real 10; with <= k_eff candidates a row would be uniform), t = 1/3/9 shown.
#
# Colour key as in analysis/plot_theme.R: blue = dependence (MI), orange = CLR normalization,
# violet = softmax, magenta = depth / walk, grey = neutral. Gene i = deep blue, gene j = deep rose
# (the 2007 convention) are IDENTITY colours drawn as bands and frames, never as matrix fills.
#
# Output: docs/whitepapers/figures/fig1_method.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_method.R

suppressPackageStartupMessages({ library(grid); library(clr); library(Matrix) })
source("analysis/plot_theme.R")
out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"

# ---- toy data and the real pipeline --------------------------------------------------------------
SEED <- 7; G <- 14; N <- 60
CAP <- 5L; KEFF <- 3; ALPHA <- 0.5; TS <- c(1L, 3L, 9L)
set.seed(SEED)
X <- matrix(rnorm(G * N), G, N)
chain <- 4:9; modB <- 10:13                       # genes 1-3 and 14: noise
for (g in chain[-1]) X[g, ] <- 0.8 * X[g - 1, ] + sqrt(1 - 0.8^2) * X[g, ]   # AR(1) chain module
fB <- rnorm(N); for (g in modB) X[g, ] <- 0.7 * fB + sqrt(1 - 0.7^2) * X[g, ]
rownames(X) <- sprintf("g%d", seq_len(G))
I <- 4L; J <- 5L; K <- 9L                         # gene i, gene j, third gene k (chain ends)

fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", threads = 1)$
  calibrate(method = "normal", combine = "stouffer")
fit$build_operator(alpha = ALPHA, softmax_keff = KEFF, softmax_cap = CAP)
D <- fit$mi; diag(D) <- 0                         # clr_calibrate zeroes the diagonal first
S <- fit$clr_scores                               # CLR score matrix (paper: S_ij)
A <- fit$operator
P <- Matrix::Diagonal(G, 1 - ALPHA) + ALPHA * Matrix::Matrix(A, sparse = TRUE)
Kt <- lapply(TS, function(tt) { Pt <- as.matrix(P); for (s in seq_len(tt - 1L)) Pt <- as.matrix(Pt %*% P)
  (Pt + t(Pt)) / 2 })
names(Kt) <- TS

# z-scores exactly as clr_core.cpp (sample sd over all G entries incl. the zeroed diagonal)
zrow <- function(r) { m <- mean(D[r, ]); s <- sd(D[r, ]); list(mu = m, sd = s, z = (D[r, ] - m) / s) }
zi <- zrow(I); zj <- zrow(J)
zraw <- vapply(c(zi$z[J], zj$z[I]), identity, 0)
stopifnot(abs((pmax(zraw[1], 0) + pmax(zraw[2], 0)) / sqrt(2) - S[I, J]) < 1e-10)
Zfull <- t(apply(D, 1, function(r) (r - mean(r)) / sd(r))); Zfull[Zfull < 0] <- 0
stopifnot(max(abs((Zfull + t(Zfull)) / sqrt(2) - replace(S, cbind(1:G, 1:G), 0))) < 1e-10)
stopifnot(all(abs(rowSums(A) - 1) < 1e-9),
          abs(exp(-sum(A[I, A[I, ] > 0] * log(A[I, A[I, ] > 0]))) - KEFF) < 1e-4)   # k_eff = 3 holds
stopifnot(A[I, K] == 0, Kt[["1"]][I, K] < Kt[["3"]][I, K], Kt[["3"]][I, K] > 0)
cat(sprintf("toy: i=%d j=%d k=%d  D_ij=%.3f  z_ij=%.2f z_ji=%.2f  S_ij=%.2f  A_ij=%.2f  A_ik=%g  K_t[i,k]=%s\n",
            I, J, K, D[I, J], zi$z[J], zj$z[I], S[I, J], A[I, J], A[I, K],
            paste(sprintf("t%d:%.3f", TS, sapply(Kt, function(m) m[I, K])), collapse = " ")))

# ---- drawing helpers (canvas in inches; native units) --------------------------------------------
W <- 7.2; H <- 7.25
col_i <- "#1b4f9c"; col_j <- "#b23a72"            # identity colours (2007: blue / pink)
fill_i <- "#2a78d6"; fill_j <- "#e87ba4"
ramp_of <- function(hi, mid) colorRampPalette(c("#fcfcfb", mid, hi))(100)
ramps <- list(D = ramp_of("#184f95", "#9ec5f4"), S = ramp_of("#a8390f", "#f6b595"),
              A = ramp_of("#3a2aa0", "#c9c2f4"), K = ramp_of("#a8345f", "#f5bfd3"))
txt <- function(label, x, y, size = 7.5, col = ink[["primary"]], face = "plain", just = "centre",
                rot = 0)
  grid.text(label, x, y, default.units = "native", just = just, rot = rot,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face))
ptext <- function(expr, x, y, size = 7.5, col = ink[["primary"]], just = "centre")
  grid.text(expr, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1, arrow = NULL)
  grid.segments(x0, y0, x1, y1, default.units = "native", arrow = arrow,
                gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt", fill = col))
rect <- function(x, y, w, h, fill = NA, col = NA, lwd = 1, alpha = 1, just = c("left", "bottom"))
  grid.rect(x, y, w, h, default.units = "native", just = just,
            gp = gpar(fill = fill, col = col, lwd = lwd, alpha = alpha, linejoin = "mitre"))
head_arrow <- function(...) seg(..., arrow = arrow(length = unit(0.075, "inches"), type = "closed"))
big_arrow <- function(x0, x1, y, col, h = 0.17, hl = 0.17) {   # filled block arrow
  grid.polygon(c(x0, x1 - hl, x1 - hl, x1, x1 - hl, x1 - hl, x0),
               c(y + h / 3, y + h / 3, y + h / 1.4, y, y - h / 1.4, y - h / 3, y - h / 3),
               default.units = "native", gp = gpar(fill = col, col = NA))
}
# a matrix: top-left (x0, ytop), cell size c; values mapped through a ramp with gamma
heat <- function(x0, ytop, cell, M, ramp, vmax = max(M), gamma = 0.7) {
  n <- nrow(M); r <- rep(seq_len(n), times = n); cc <- rep(seq_len(n), each = n)
  v <- pmin(pmax(as.vector(M) / vmax, 0), 1)^gamma
  fl <- ramp[pmax(1, ceiling(v * 100))]
  fl[r == cc] <- "#ececea"
  rect(x0 + (cc - 1) * cell, ytop - r * cell, cell, cell, fill = fl, col = "white", lwd = 0.25)
}
# 2007-style bands: translucent row-i band and column-j band extending past the matrix, plus frames
bands <- function(x0, ytop, cell, n, i, j, ext = 0.13, i_on = TRUE, j_on = TRUE) {
  if (i_on) { rect(x0 - ext, ytop - i * cell, n * cell + 2 * ext, cell, fill = fill_i, alpha = 0.28)
    rect(x0, ytop - i * cell, n * cell, cell, col = col_i, lwd = 1.3) }
  if (j_on) { rect(x0 + (j - 1) * cell, ytop - n * cell - ext, cell, n * cell + 2 * ext, fill = fill_j, alpha = 0.28)
    rect(x0 + (j - 1) * cell, ytop - n * cell, cell, n * cell, col = col_j, lwd = 1.3) }
}
mark <- function(x0, ytop, cell, i, j, col = ink[["primary"]], lwd = 1.6)   # traced cell
  rect(x0 + (j - 1) * cell, ytop - i * cell, cell, cell, col = col, lwd = lwd)
gene_tags <- function(x0, ytop, cell, n, i, j, size = 7.5) {
  txt("gene i", x0 - 0.16, ytop - (i - 0.5) * cell, size, col_i, "bold", just = "right")
  txt("gene j", x0 + (j - 0.5) * cell, ytop + 0.12, size, col_j, "bold")
}
chip <- function(num, x, y) {
  grid.circle(x, y, 0.105, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
  txt(num, x, y - 0.002, 8.5, "white", "bold")
}
tag <- function(label, x, y, col) {     # transformer counterpart label, width from the text
  w <- convertWidth(grobWidth(textGrob(label, gp = gpar(fontfamily = FONT, fontsize = 7))),
                    "native", TRUE) + 0.18
  grid.roundrect(x, y, w, 0.2, default.units = "native", just = c("left", "centre"),
                 r = unit(0.04, "inches"), gp = gpar(fill = "#f1f0ed", col = NA))
  rect(x, y - 0.1, 0.035, 0.2, fill = col)
  txt(label, x + 0.09, y, 7, ink[["secondary"]], just = "left")
}
panel_head <- function(num, title, x, y, tagtxt, col) {
  chip(num, x + 0.105, y)
  txt(title, x + 0.27, y, 9.5, face = "bold", just = "left")
  tag(tagtxt, x, y - 0.27, col)
}

# distribution of one row's MI values, with its observed value in the tail (as 2007 Fig. 1A)
dist_panel <- function(x0, x1, ybase, h, zz, obs_idx, line, fill, who) {
  vals <- D[who, ]; mu <- zz$mu; s <- zz$sd; obs <- vals[obs_idx]
  lo <- 0; hi <- max(vals) * 1.08
  xs <- seq(lo, hi, length.out = 200); dens <- dnorm(xs, mu, s); ymax <- max(dens)
  X_ <- function(v) x0 + (v - lo) / (hi - lo) * (x1 - x0)
  Y_ <- function(d) ybase + d / ymax * h
  grid.polygon(c(X_(xs), X_(hi), X_(lo)), c(Y_(dens), ybase, ybase), default.units = "native",
               gp = gpar(fill = fill, col = NA, alpha = 0.28))
  tl <- xs >= obs
  grid.polygon(c(X_(xs[tl]), X_(hi), X_(obs)), c(Y_(dens[tl]), ybase, ybase),
               default.units = "native", gp = gpar(fill = fill, col = NA, alpha = 0.85))
  grid.lines(X_(xs), Y_(dens), default.units = "native", gp = gpar(col = line, lwd = 1.3))
  seg(x0, ybase, x1, ybase, ink[["secondary"]], 0.8)
  others <- setdiff(seq_along(vals), c(who, obs_idx))    # the row's other MI values as a rug
  seg(X_(vals[others]), ybase, X_(vals[others]), ybase - 0.05, ink[["secondary"]], 0.6)
  seg(X_(vals[who]), ybase, X_(vals[who]), ybase - 0.05, ink[["muted"]], 0.6)   # zeroed diagonal
  seg(X_(obs), ybase - 0.09, X_(obs), Y_(dnorm(obs, mu, s)) + 0.0, line, 1.6)
  grid.points(X_(obs), ybase, default.units = "native", pch = 19, size = unit(5, "pt"),
              gp = gpar(col = line))
  list(X = X_, Y = Y_, obs = obs, mu = mu)
}

# bar chart of one matrix row over the G genes (zeros drawn as a baseline tick only)
bars <- function(x0, ybase, w, h, v, col, vmax, hi_idx = integer(0), cols_other = NULL) {
  n <- length(v); bw <- w / n
  for (g in seq_len(n)) {
    bh <- v[g] / vmax * h
    if (bh > 0) rect(x0 + (g - 1) * bw + 0.012, ybase, bw - 0.024, bh, fill = col)
  }
  seg(x0, ybase, x0 + w, ybase, ink[["secondary"]], 0.8)
}

# ---- the figure ---------------------------------------------------------------------------------
draw_fig <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, W), yscale = c(0, H)))
  cellA <- 0.12
  short_tags <- function(x0, ytop, cell, n, i, j) {   # short i / j tags (long ones appear in panel 1)
    txt("i", x0 - 0.1, ytop - (i - 0.5) * cell, 8, col_i, "bold.italic", just = "right")
    txt("j", x0 + (j - 0.5) * cell, ytop + 0.09, 8, col_j, "bold.italic")
  }

  # ============ ROW 1 =====================================================================
  yh <- H - 0.17; ytop1 <- H - 0.98; mh <- G * cellA
  # ---- panel 1: D -----------------------------------------------------------------------
  x1 <- 0.5
  panel_head("1", "Dependence", 0.1, yh, "query–key score", pal[["blue"]])
  heat(x1, ytop1, cellA, D, ramps$D, vmax = max(D))
  bands(x1, ytop1, cellA, G, I, J, ext = 0.1)
  mark(x1, ytop1, cellA, I, J)
  gene_tags(x1, ytop1, cellA, G, I, J)
  ptext(expression(D[ij] == MI(i, j)), x1 + G * cellA / 2, ytop1 - mh - 0.27, 8)
  big_arrow(2.33, 2.62, ytop1 - mh / 2, pal[["orange"]], h = 0.22, hl = 0.15)

  # ---- panel 2: CLR ------------------------------------------------------------------------
  panel_head("2", "Context normalization (CLR)", 2.7, yh, "scaling the logits", pal[["orange"]])
  dx0 <- 2.72; dx1 <- 4.12
  y_a <- ytop1 - 0.72; y_b <- ytop1 - 1.92
  ptext(expression("row " * italic(i) * " of " * D), dx0, ytop1 + 0.04, 7.5, col_i, just = "left")
  di <- dist_panel(dx0, dx1, y_a, 0.5, zi, J, col_i, fill_i, I)
  ptext(expression("row " * italic(j) * " of " * D), dx0, y_b + 0.72, 7.5, col_j, just = "left")
  dj <- dist_panel(dx0, dx1, y_b, 0.5, zj, I, col_j, fill_j, J)
  txt("MI", dx1 + 0.03, y_a, 7, ink[["secondary"]], just = "left")
  txt("MI", dx1 + 0.03, y_b, 7, ink[["secondary"]], just = "left")
  ptext(expression(D[ij]), di$X(di$obs), y_a + 0.2, 7.5, col_i)
  ptext(expression(D[ij]), dj$X(dj$obs), y_b + 0.2, 7.5, col_j)
  ptext(substitute(z[ij] == (D[ij] - mu[i]) / sigma[i] ~ "=" ~ v, list(v = sprintf("%.1f", zi$z[J]))),
        dx0, y_a - 0.23, 7.5, col_i, just = "left")
  ptext(substitute(z[ji] == (D[ij] - mu[j]) / sigma[j] ~ "=" ~ v, list(v = sprintf("%.1f", zj$z[I]))),
        dx0, y_b - 0.23, 7.5, col_j, just = "left")
  txt("ticks: the row's 14 values", dx0, y_b - 0.42, 7, ink[["secondary"]], just = "left")

  # combination: two lines converge on one arrow, formula beside it
  cy <- ytop1 - mh / 2 - 0.02; cx <- 4.4
  seg(dx1 + 0.16, y_a + 0.1, cx, cy, col_i, 1.2)
  seg(dx1 + 0.16, y_b + 0.1, cx, cy, col_j, 1.2)
  x2 <- 5.42
  big_arrow(cx, x2 - 0.25, cy, pal[["orange"]], h = 0.2, hl = 0.14)
  ptext(expression(S[ij] == frac(z[ij]^"+" + z[ji]^"+", sqrt(2))), (cx + x2) / 2 - 0.12, cy - 0.42, 8)
  txt("z⁺ = max(z, 0)", (cx + x2) / 2 - 0.12, cy - 0.78, 7, ink[["secondary"]])
  heat(x2, ytop1, cellA, S, ramps$S, vmax = max(S))
  bands(x2, ytop1, cellA, G, I, J, ext = 0.1)
  mark(x2, ytop1, cellA, I, J)
  short_tags(x2, ytop1, cellA, G, I, J)
  ptext(substitute(S[ij] == v, list(v = sprintf("%.1f", S[I, J]))), x2 + G * cellA / 2,
        ytop1 - mh - 0.27, 8)

  # ============ ROW 2 =====================================================================
  sep <- H - 3.5
  seg(0.1, sep, W - 0.1, sep, ink[["grid"]], 1)
  yh2 <- sep - 0.2; ytop2 <- yh2 - 0.76; c3 <- 0.09; mh3 <- G * c3
  # ---- panel 3: softmax -------------------------------------------------------------------
  panel_head("3", "Softmax attention", 0.1, yh2, "softmax over the keys", pal[["violet"]])
  xs3 <- 0.5; xa3 <- 2.62; gw <- G * c3
  heat(xs3, ytop2, c3, S, ramps$S, vmax = max(S))
  bands(xs3, ytop2, c3, G, I, J, ext = 0.06)
  mark(xs3, ytop2, c3, I, J, lwd = 1.4)
  heat(xa3, ytop2, c3, A, ramps$A, vmax = max(A), gamma = 0.6)
  bands(xa3, ytop2, c3, G, I, J, ext = 0.06)
  mark(xa3, ytop2, c3, I, J, lwd = 1.4)
  short_tags(xs3, ytop2, c3, G, I, J); short_tags(xa3, ytop2, c3, G, I, J)
  # bar strips of row i, aligned to the matrix columns
  by <- 0.9; bh <- 0.5
  keep <- which(A[I, ] > 0)
  rect(xs3 + (J - 1) * c3, by - 0.02, c3, bh + 0.1, fill = fill_j, alpha = 0.28)
  rect(xa3 + (J - 1) * c3, by - 0.02, c3, bh + 0.1, fill = fill_j, alpha = 0.28)
  sv <- S[I, ]; vm <- max(sv)
  for (g in seq_len(G)) if (sv[g] > 0)
    rect(xs3 + (g - 1) * c3 + 0.012, by, c3 - 0.024, sv[g] / vm * bh,
         fill = if (g %in% keep) pal[["orange"]] else "#f6b595")
  seg(xs3, by, xs3 + gw, by, ink[["secondary"]], 0.8)
  for (g in seq_len(G)) if (A[I, g] > 0)
    rect(xa3 + (g - 1) * c3 + 0.012, by, c3 - 0.024, A[I, g] / max(A[I, ]) * bh, fill = pal[["violet"]])
  seg(xa3, by, xa3 + gw, by, ink[["secondary"]], 0.8)
  ptext(expression("row " * italic(i) * " of S: keep the top 5"), xs3 + gw / 2, by - 0.14, 7)
  ptext(expression("(dark bars; real method: 50)"), xs3 + gw / 2, by - 0.27, 7, ink[["secondary"]])
  ptext(expression("row " * italic(i) * " of A: sums to 1"), xa3 + gw / 2, by - 0.14, 7)
  ptext(substitute(A[ij] == v, list(v = sprintf("%.2f", A[I, J]))), xa3 + (J - 0.5) * c3 + 0.1,
        by + A[I, J] / max(A[I, ]) * bh - 0.05, 7.5, col_j, just = "left")
  ptext(substitute(S[ij] == v, list(v = sprintf("%.1f", S[I, J]))), xs3 + (J - 0.5) * c3 + 0.1,
        by + bh - 0.05, 7.5, col_j, just = "left")
  # arrow with the formula
  zx <- (xs3 + gw + xa3) / 2; ay <- ytop2 - mh3 / 2
  big_arrow(xs3 + gw + 0.1, xa3 - 0.1, ay, pal[["violet"]], h = 0.2, hl = 0.12)
  ptext(expression(A[i] == softmax(S[i] / tau[i])), zx, ay - 0.27, 7)
  txt("τ per gene so that", zx, ay - 0.5, 7, ink[["secondary"]])
  ptext(substitute(k[eff] == v ~ "(real: 10)", list(v = KEFF)), zx, ay - 0.64, 7)

  # ---- panel 4: walk -----------------------------------------------------------------------
  x4 <- 4.2
  big_arrow(xa3 + gw + 0.1, x4 - 0.1, ytop2 - mh3 / 2, pal[["magenta"]], h = 0.2, hl = 0.1)
  panel_head("4", "Walk: depth", x4, yh2, "stacked attention layers", pal[["magenta"]])
  ptext(expression(P == 0.5 * I + 0.5 * A), x4 + 1.9, yh2 - 0.02, 8, just = "left")
  c4 <- 0.062; mh4 <- G * c4; gw4 <- G * c4
  xs4 <- x4 + 0.12 + (0:2) * 0.97
  y4 <- ytop2 - 0.1
  offs <- lapply(Kt, function(m) { diag(m) <- 0; m }); kmax <- max(unlist(lapply(offs, function(m) m[I, ])))
  bby <- 1.05; bbh <- 0.45
  for (q in seq_along(TS)) {
    M <- Kt[[q]]; off <- offs[[q]]
    heat(xs4[q], y4, c4, M, ramps$K, vmax = max(off), gamma = 0.7)
    bands(xs4[q], y4, c4, G, I, J, ext = 0.04, j_on = FALSE)
    rect(xs4[q] + (K - 1) * c4, y4 - I * c4, c4, c4, col = ink[["primary"]], lwd = 1.5)
    txt(sprintf("t = %d", TS[q]), xs4[q] + gw4 / 2, y4 + 0.14, 8, face = "bold")
    txt("k", xs4[q] + (K - 0.5) * c4, y4 + 0.04, 7, ink[["primary"]], "italic")
    # row i of K_t as bars on the same gene axis (gene i itself omitted); gene k marked
    rect(xs4[q] + (K - 1) * c4, bby - 0.02, c4, bbh + 0.06, fill = ink[["muted"]], alpha = 0.3)
    for (g in seq_len(G)) if (off[I, g] > 0)
      rect(xs4[q] + (g - 1) * c4 + 0.008, bby, c4 - 0.016, off[I, g] / kmax * bbh, fill = pal[["magenta"]])
    seg(xs4[q], bby, xs4[q] + gw4, bby, ink[["secondary"]], 0.8)
    txt("k", xs4[q] + (K - 0.5) * c4, bby - 0.1, 7, ink[["primary"]], "italic")
    ptext(substitute(K[ik] == v, list(v = sprintf("%.3f", M[I, K]))), xs4[q] + gw4 / 2,
          bby - 0.25, 7.5)
    if (q < 3) head_arrow(xs4[q] + gw4 + 0.03, y4 - mh4 / 2, xs4[q + 1] - 0.03, y4 - mh4 / 2,
                          pal[["magenta"]], 1.2)
  }
  txt("i", xs4[1] - 0.06, y4 - (I - 0.5) * c4, 8, col_i, "bold.italic", just = "right")
  ptext(expression("row " * italic(i) * " of " * K[t]), xs4[1], bby + bbh + 0.1, 7, ink[["secondary"]], just = "left")
  ptext(expression(paste(K[t] == (P^t + (P^t)^T) / 2, ";  boxed: gene ", italic(k), ", not attended directly (", A[ik] == 0, ")")),
        x4 + 0.1, 0.58, 7, just = "left")

  # ============ footer ================================================================================
  txt("No learned parameters: every step is computed from the data.", 0.1, 0.3, 7.5,
      ink[["primary"]], "bold", just = "left")
  ptext(substitute(paste("Toy: ", g, " genes, ", n, " samples, planted modules; ", k[eff] == ke, ", top ", cp,
                         ", t = ", tt, " (real method: 10, 50, 9)."),
                   list(g = G, n = N, ke = KEFF, cp = CAP, tt = paste(TS, collapse = ", "))),
        0.1, 0.13, 7, ink[["secondary"]], just = "left")

  popViewport()
}

dev_png <- function() png(file.path(out_dir, "fig1_method.png"), width = W * 300, height = H * 300,
                          res = 300, type = "cairo", bg = "white")
dev_png(); draw_fig(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "fig1_method.pdf"), width = W, height = H, bg = "white")
draw_fig(); invisible(dev.off())
cat("wrote", file.path(out_dir, "fig1_method.{png,pdf}"), "\n")
