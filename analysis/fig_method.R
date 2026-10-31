# Fig. 1 (method): CLR as attention, one traced cell through a worked toy example.
# Layout (rev. 2): four numbered quadrants (1 -> 2, S connector down, 3 -> 4) after the best PaperBanana
# draft (layout only; its matrices were invented, all content here is computed).
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
stopifnot(sum(A[I, ] > 0) == CAP, A[I, K] == 0, Kt[["1"]][I, K] < Kt[["3"]][I, K], Kt[["3"]][I, K] > 0)
cat(sprintf("toy: i=%d j=%d k=%d  D_ij=%.3f  z_ij=%.2f z_ji=%.2f  S_ij=%.2f  A_ij=%.2f  A_ik=%g  K_t[i,k]=%s\n",
            I, J, K, D[I, J], zi$z[J], zj$z[I], S[I, J], A[I, J], A[I, K],
            paste(sprintf("t%d:%.3f", TS, sapply(Kt, function(m) m[I, K])), collapse = " ")))

# ---- drawing helpers (canvas in inches; native units) --------------------------------------------
W <- 7.2; H <- 7.85
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

# distribution of one row's MI values (normal fit, as CLR uses), observed value in the tail (as 2007 Fig. 1A)
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
  seg(X_(obs), ybase - 0.09, X_(obs), Y_(dnorm(obs, mu, s)), line, 1.6)
  grid.points(X_(obs), ybase, default.units = "native", pch = 19, size = unit(5.5, "pt"),
              gp = gpar(col = line))
  list(X = X_, Y = Y_, obs = obs, mu = mu)
}

# ---- the figure ---------------------------------------------------------------------------------
# win = c(x0, x1, y0, y1) in canvas inches crops the drawing to one window (used for per-panel exports)
draw_fig <- function(win = c(0, W, 0, H), solo = FALSE) {   # solo: omit inter-panel dividers, arrows, connector
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = win[1:2], yscale = win[3:4], clip = "on"))
  short_tags <- function(x0, ytop, cell, n, i, j) {   # short i / j tags (long ones appear in panel 1)
    txt("i", x0 - 0.09, ytop - (i - 0.5) * cell, 8, col_i, "bold.italic", just = "right")
    txt("j", x0 + (j - 0.5) * cell, ytop + 0.09, 8, col_j, "bold.italic")
  }
  XD <- 3.6                                          # vertical divider between the two columns
  sep <- H - 3.75                                    # horizontal divider between the two rows
  if (!solo) {
    seg(0.1, sep, W - 0.1, sep, ink[["grid"]], 1)
    seg(XD, H - 0.08, XD, sep + 0.0, ink[["grid"]], 1)
    seg(XD, sep, XD, 0.62, ink[["grid"]], 1)
  }

  # ============ ROW 1 =====================================================================
  yh <- H - 0.17
  # ---- panel 1: D -----------------------------------------------------------------------
  c1 <- 0.165; m1 <- G * c1; ytop1 <- H - 1.02; x1 <- 0.64
  panel_head("1", "Dependence", 0.1, yh, "query–key score", pal[["blue"]])
  heat(x1, ytop1, c1, D, ramps$D, vmax = max(D))
  bands(x1, ytop1, c1, G, I, J, ext = 0.1)
  mark(x1, ytop1, c1, I, J)
  gene_tags(x1, ytop1, c1, G, I, J)
  ptext(substitute(D[ij] == MI(i, j) ~ "=" ~ v, list(v = sprintf("%.2f", D[I, J]))), x1 + m1 / 2,
        ytop1 - m1 - 0.25, 8)
  if (!solo) big_arrow(3.0, 3.74, ytop1 - m1 / 2, pal[["orange"]], h = 0.24, hl = 0.17)

  # ---- panel 2: CLR ------------------------------------------------------------------------
  panel_head("2", "Context normalization (CLR)", XD + 0.1, yh, "scaling the logits", pal[["orange"]])
  dx0 <- 3.9; dx1 <- 5.25
  y_a <- H - 1.55; y_b <- H - 2.9
  ptext(expression("row " * italic(i) * " of " * D), dx0, y_a + 0.62, 7.5, col_i, just = "left")
  di <- dist_panel(dx0, dx1, y_a, 0.5, zi, J, col_i, fill_i, I)
  ptext(expression("row " * italic(j) * " of " * D), dx0, y_b + 0.62, 7.5, col_j, just = "left")
  dj <- dist_panel(dx0, dx1, y_b, 0.5, zj, I, col_j, fill_j, J)
  txt("MI", dx1, y_a - 0.17, 7, ink[["secondary"]], just = "right")
  txt("MI", dx1, y_b - 0.17, 7, ink[["secondary"]], just = "right")
  ptext(expression(D[ij]), di$X(di$obs), y_a + 0.2, 7.5, col_i)
  ptext(expression(D[ij]), dj$X(dj$obs), y_b + 0.2, 7.5, col_j)
  ptext(substitute(z[ij] == (D[ij] - mu[i]) / sigma[i] ~ "=" ~ v, list(v = sprintf("%.1f", zi$z[J]))),
        dx0, y_a - 0.34, 7.5, col_i, just = "left")
  ptext(substitute(z[ji] == (D[ij] - mu[j]) / sigma[j] ~ "=" ~ v, list(v = sprintf("%.1f", zj$z[I]))),
        dx0, y_b - 0.34, 7.5, col_j, just = "left")

  # S matrix (right column), the formula above it, two lines converging on the block arrow
  c2 <- 0.1; m2 <- G * c2; x2 <- 5.6; ytop2m <- H - 1.78
  ptext(expression(S[ij] == frac(z[ij]^"+" + z[ji]^"+", sqrt(2))), x2 + m2 / 2, H - 0.96, 8)
  txt("z⁺ = max(z, 0)", x2 + m2 / 2, H - 1.36, 7, ink[["secondary"]])
  cy <- ytop2m - m2 / 2
  seg(dx1 + 0.08, y_a + 0.0, 5.3, cy + 0.04, col_i, 1.2)
  seg(dx1 + 0.08, y_b + 0.0, 5.3, cy - 0.04, col_j, 1.2)
  big_arrow(5.3, x2 - 0.04, cy, pal[["orange"]], h = 0.2, hl = 0.12)
  heat(x2, ytop2m, c2, S, ramps$S, vmax = max(S))
  bands(x2, ytop2m, c2, G, I, J, ext = 0.07)
  mark(x2, ytop2m, c2, I, J)
  short_tags(x2, ytop2m, c2, G, I, J)
  ptext(substitute(S[ij] == v, list(v = sprintf("%.1f", S[I, J]))), x2 + m2 / 2, ytop2m - m2 - 0.2, 8)

  # ============ ROW 2 =====================================================================
  yh2 <- sep - 0.22; ytop3 <- yh2 - 0.85; c3 <- 0.08; gw <- G * c3
  # ---- panel 3: softmax -------------------------------------------------------------------
  panel_head("3", "Softmax attention", 0.1, yh2, "softmax over the keys", pal[["violet"]])
  xs3 <- 0.56; xa3 <- 2.28
  heat(xs3, ytop3, c3, S, ramps$S, vmax = max(S))
  bands(xs3, ytop3, c3, G, I, J, ext = 0.05)
  mark(xs3, ytop3, c3, I, J, lwd = 1.4)
  heat(xa3, ytop3, c3, A, ramps$A, vmax = max(A), gamma = 0.6)
  bands(xa3, ytop3, c3, G, I, J, ext = 0.05)
  mark(xa3, ytop3, c3, I, J, lwd = 1.4)
  short_tags(xs3, ytop3, c3, G, I, J); short_tags(xa3, ytop3, c3, G, I, J)
  txt("S", xs3 + gw / 2, ytop3 + 0.09, 8, ink[["primary"]], "bold")
  txt("A", xa3 + gw / 2, ytop3 + 0.09, 8, ink[["primary"]], "bold")
  # S (panel 2) -> S (panel 3): same matrix, elbow connector
  cx_ <- 2.12; ye <- ytop3 + 0.22; xe <- 1.62
  if (!solo) {
    grid.lines(c(x2 + m2 - 0.12, x2 + m2 - 0.12, cx_, cx_, xe, xe),
               c(ytop2m - m2 - 0.05, sep, sep, ye, ye, ytop3 + 0.03), default.units = "native",
               arrow = arrow(length = unit(0.075, "inches"), type = "closed"),
               gp = gpar(col = pal[["orange"]], lwd = 1.6, fill = pal[["orange"]], linejoin = "mitre"))
    txt("same S matrix", (x2 + m2 - 0.12 + cx_) / 2, sep + 0.1, 7, pal[["orange"]], "bold")
  }
  # arrow S -> A
  ay <- ytop3 - 0.5 * gw
  big_arrow(xs3 + gw + 0.08, xa3 - 0.08, ay, pal[["violet"]], h = 0.2, hl = 0.11)
  txt("softmax", (xs3 + gw + xa3) / 2, ay + 0.2, 7, pal[["violet"]], "bold")
  # row i of S and of A as bar strips under the matrices, aligned to the columns
  by <- 1.25; bh <- 0.45
  keep <- which(A[I, ] > 0)
  rect(xs3 + (J - 1) * c3, by - 0.02, c3, bh + 0.1, fill = fill_j, alpha = 0.28)
  rect(xa3 + (J - 1) * c3, by - 0.02, c3, bh + 0.1, fill = fill_j, alpha = 0.28)
  sv <- S[I, ]; vm <- max(sv)
  for (g in seq_len(G)) if (sv[g] > 0)
    rect(xs3 + (g - 1) * c3 + 0.01, by, c3 - 0.02, sv[g] / vm * bh,
         fill = if (g %in% keep) pal[["orange"]] else "#f6c9b2")
  seg(xs3, by, xs3 + gw, by, ink[["secondary"]], 0.8)
  for (g in seq_len(G)) if (A[I, g] > 0)
    rect(xa3 + (g - 1) * c3 + 0.01, by, c3 - 0.02, A[I, g] / max(A[I, ]) * bh, fill = pal[["violet"]])
  seg(xa3, by, xa3 + gw, by, ink[["secondary"]], 0.8)
  ptext(expression("row " * italic(i) * " of " * S), xs3 + gw / 2, by - 0.13, 7.5)
  ptext(substitute("dark: top " * cp * " kept (real: 50)", list(cp = CAP)), xs3 + gw / 2, by - 0.26, 7,
        ink[["secondary"]])
  ptext(expression("row " * italic(i) * " of " * A), xa3 + gw / 2, by - 0.13, 7.5)
  txt("sums to 1", xa3 + gw / 2, by - 0.26, 7, ink[["secondary"]])
  ptext(substitute(S[ij] == v, list(v = sprintf("%.1f", S[I, J]))), xs3 + (J - 0.5) * c3 + 0.07,
        by + bh - 0.06, 7.5, col_j, just = "left")
  ptext(substitute(A[ij] == v, list(v = sprintf("%.2f", A[I, J]))), xa3 + (J - 0.5) * c3 + 0.07,
        by + A[I, J] / max(A[I, ]) * bh - 0.06, 7.5, col_j, just = "left")
  ptext(expression(A[i] == softmax(S[i] / tau[i])), (xs3 + xa3 + gw) / 2, 0.82, 8)
  ptext(substitute(paste(tau[i], " set per gene so that ", k[eff] == v, " (real: 10)"), list(v = KEFF)),
        (xs3 + xa3 + gw) / 2, 0.66, 7.5)

  # ---- panel 4: walk -----------------------------------------------------------------------
  x4 <- XD + 0.25; c4 <- 0.072; gw4 <- G * c4; step4 <- 1.12
  if (!solo) big_arrow(xa3 + gw + 0.07, x4 + 0.0, ytop3 - gw / 2 + 0.0, pal[["magenta"]], h = 0.2, hl = 0.11)
  panel_head("4", "Walk: depth", XD + 0.1, yh2, "stacked attention layers", pal[["magenta"]])
  ptext(expression(P == 0.5 * I + 0.5 * A), XD + 3.5 - 0.1, yh2 - 0.02, 8, just = "right")
  xs4 <- x4 + (0:2) * step4
  offs <- lapply(Kt, function(m) { diag(m) <- 0; m }); kmax <- max(unlist(lapply(offs, function(m) m[I, ])))
  by4 <- 1.25; bbh <- 0.5
  for (q in seq_along(TS)) {
    M <- Kt[[q]]; off <- offs[[q]]
    heat(xs4[q], ytop3, c4, M, ramps$K, vmax = max(off), gamma = 0.7)
    bands(xs4[q], ytop3, c4, G, I, J, ext = 0.04, j_on = FALSE)
    rect(xs4[q] + (K - 1) * c4, ytop3 - I * c4, c4, c4, col = ink[["primary"]], lwd = 1.6)
    txt(sprintf("t = %d", TS[q]), xs4[q] + gw4 / 2, ytop3 + 0.14, 8.5, face = "bold")
    txt("k", xs4[q] + (K - 0.5) * c4, ytop3 + 0.03, 7, ink[["primary"]], "bold.italic")
    # row i of K_t as bars on the same gene axis (gene i itself omitted); gene k marked
    rect(xs4[q] + (K - 1) * c4, by4 - 0.02, c4, bbh + 0.06, fill = ink[["muted"]], alpha = 0.3)
    for (g in seq_len(G)) if (off[I, g] > 0)
      rect(xs4[q] + (g - 1) * c4 + 0.008, by4, c4 - 0.016, off[I, g] / kmax * bbh, fill = pal[["magenta"]])
    seg(xs4[q], by4, xs4[q] + gw4, by4, ink[["secondary"]], 0.8)
    if (off[I, K] == 0)   # a true zero: hairline across gene k's slot
      seg(xs4[q] + (K - 1) * c4, by4, xs4[q] + K * c4, by4, ink[["primary"]], 2.4)
    grid.points(xs4[q] + (K - 0.5) * c4, by4 - 0.09, default.units = "native", pch = 17,
                size = unit(4, "pt"), gp = gpar(col = ink[["primary"]]))
    ptext(substitute(K[ik] == v, list(v = sprintf("%.3f", M[I, K]))), xs4[q] + gw4 / 2,
          by4 - 0.26, 7.5)
    if (q < 3) head_arrow(xs4[q] + gw4 + 0.015, ytop3 - gw4 / 2, xs4[q + 1] - 0.015, ytop3 - gw4 / 2,
                          pal[["magenta"]], 1.2)
  }
  txt("i", xs4[1] - 0.06, ytop3 - (I - 0.5) * c4, 8, col_i, "bold.italic", just = "right")
  ptext(expression("row " * italic(i) * " of " * K[t] * " (gene " * italic(i) * " omitted)"),
        xs4[1], by4 + bbh + 0.12, 7, ink[["secondary"]], just = "left")
  ptext(expression(K[t] == (P^t + (P^t)^T) / 2), x4 + 1.5 * step4 - 0.05, 0.82, 8)
  ptext(expression(paste("box / triangle: gene ", italic(k), ", not attended directly (", A[ik] == 0, ")")),
        x4 - 0.05, 0.66, 7, just = "left")

  # ============ footer ================================================================================
  seg(0.1, 0.55, W - 0.1, 0.55, ink[["grid"]], 1)
  txt("No learned parameters: every step is computed from the data.", 0.1, 0.42, 7.5,
      ink[["primary"]], "bold", just = "left")
  ptext(substitute(paste("Toy: ", g, " genes, ", n, " samples, planted modules; ", k[eff] == ke, ", top ", cp,
                         ", t = ", tt, " (real method: 10, 50, 9)."),
                   list(g = G, n = N, ke = KEFF, cp = CAP, tt = paste(TS, collapse = ", "))),
        0.1, 0.27, 7, ink[["secondary"]], just = "left")
  txt("Curves: normal fits (mean, sd) of the row's 14 values (ticks). Light = low, dark = high; grey diagonal = self, excluded.",
      0.1, 0.12, 7, ink[["secondary"]], just = "left")
  popViewport()
}

# Optional per-panel export (reference images for the PaperBanana panel-by-panel track):
#   FIG_PANELS_DIR=<dir> Rscript analysis/fig_method.R   writes <dir>/panel{1..4}.png (the same
#   drawing cropped to each quadrant, 300 dpi) and exits without touching fig1_method.{png,pdf}.
pdir <- Sys.getenv("FIG_PANELS_DIR")
if (nzchar(pdir)) {
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  sepy <- H - 3.75; m <- 0.04
  wins <- list(c(m, 3.6 - m, sepy + m, H - m), c(3.6 + m, W - m, sepy + m, H - m),
               c(m, 3.6 - m, 0.6, sepy - m), c(3.6 + m, W - m, 0.6, sepy - m))
  for (q in 1:4) {
    w <- wins[[q]]
    png(file.path(pdir, sprintf("panel%d.png", q)), width = (w[2] - w[1]) * 300,
        height = (w[4] - w[3]) * 300, res = 300, type = "cairo", bg = "white")
    draw_fig(w, solo = TRUE); invisible(dev.off())
  }
  cat("wrote", file.path(pdir, "panel{1..4}.png"), "\n"); quit(save = "no")
}
dev_png <- function() png(file.path(out_dir, "fig1_method.png"), width = W * 300, height = H * 300,
                          res = 300, type = "cairo", bg = "white")
dev_png(); draw_fig(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "fig1_method.pdf"), width = W, height = H, bg = "white")
draw_fig(); invisible(dev.off())
cat("wrote", file.path(out_dir, "fig1_method.{png,pdf}"), "\n")
